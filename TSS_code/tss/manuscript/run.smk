# Homer and other helpers
from tss.config import HOMER_PATH
import tss.utils.sankey_helper as sh
import os
from tss.visualize.fig_utils import helper_save as hs
from tss.data.cumsum import get_cum_and_plot, create_cdf

# Visual libraries
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
mpl.use('Agg')
import seaborn as sns

annotation = config['annotation']
annos=list(annotation.keys())
ncrna = config['ncrna_keys']


def get_gff(wildcards):
    return annotation[wildcards.genome]

def get_ncrna(wildcards):
    return ",".join(ncrna[wildcards.genome])


def get_overlap(in_f, landmark_f, out_f):
    # Create temporary bed file and then convert to peak-like
    # Create bed
    cmd = f"{HOMER_PATH}/pos2bed.pl {in_f} > {in_f}.bed"

    os.system(cmd)
    #-v says get the non-intersecting peaks
    cmd = f"bedtools intersect -s -v -a {in_f}.bed -b {landmark_f} > {out_f}.bed"

    os.system(cmd)

    bed_df = pd.read_csv(f"{out_f}.bed", sep='\t', header=None)
    bed_df[3].to_csv(out_f, index=None, header=False) #Save just the index
    return


report: "workflow.rst"

rule all:
    input: #expand("data/{genome}/anno/peaks.distance.tsv", genome=list(config['annotation'].keys())),
           #expand("data/{genome}/anno/peaks.merge.tss.minimal.anno", genome=list(config['annotation'].keys())),
           expand("data/{genome}/genome/introns.bed", genome=annos),
           expand('data/peaksInfo/peaks_scoring.thresh{thresh}.tsv', thresh=config['cpm_thresh']),
           expand("data/{genome}/refType/peaks.refType.tsv", genome=annos),
           expand("data/peaksInfo/tissues_scoring.thresh{thresh}.tsv", thresh=config['cpm_thresh']),
           expand("data/{genome}/allTSS/thresh{thresh}_numSamp{numSamp}_peaksSankey.png",genome=annos, thresh=config['cpm_thresh'],numSamp=config['numsamples_thresh']),
           expand("data/{genome}/allTSS/SF4_pc_refPromoters_thresh{thresh}_numSamp{numSamp}.png",genome=annos, thresh=config['cpm_thresh'],numSamp=config['numsamples_thresh']),
           expand("data/{genome}/allTSS/F1_tissuesCaptured_thresh{thresh}_numSamp{numSamp}.png", genome=annos, thresh=config['cpm_thresh'],numSamp=config['numsamples_thresh'])



rule preprocess:
    """ Sets peaks as 1bp """
    #input: "peaksInfo/samples.merge"
    params: "data/peaksInfo/samples.merge"
    output: "data/peaksInfo/peaks.merge.tss", "data/peaksInfo/peaks.merge.tss.minimal"
    run:
        #print('samples merge', params[0])
        peaks_df = pd.read_csv(params[0], sep="\t")
        peaks_df["Start"] = (
                1.0 * (peaks_df["End"] + peaks_df["Start"]) / 2).astype(int)
        peaks_df["End"] = peaks_df["Start"]  # length 1
        peaks_df["actual_start"] = peaks_df["End"]
        peaks_df.to_csv(output[0], sep="\t", index=False)
        peaks_df[['Chr', 'Start', 'End', 'Strand', 'Stat']].to_csv(output[1], sep='\t')


rule homeranno:
    input: "data/peaksInfo/peaks.merge.tss"
    params:
        genome=config['genome_fasta'],
        anno=get_gff,
        HOMER_PATH=HOMER_PATH
    output: "data/{genome}/anno/peaks.merge.tss.anno"
    shell: '{params.HOMER_PATH}/annotatePeaks.pl {input} {params.genome} -gff {params.anno} > {output}'


rule combine_merge_and_anno:
    # def combine_merge_with_anno(merge_file, anno_file):
#     """ Function that adds the column Annotation from anno_file to the merge_file"""
    input:
        merged="data/peaksInfo/peaks.merge.tss.minimal",
        homeranno="data/{genome}/anno/peaks.merge.tss.anno"
    output: peaks_df="data/{genome}/anno/peaks.merge.tss.minimal.anno"
    run:
        merge_df = pd.read_csv(input.merged, index_col=0, sep='\t')
        anno_df = pd.read_csv(input.homeranno, sep='\t', na_filter=False)
        anno_df.index.name = "ID"
        anno_df = anno_df['Annotation']
        merge_df = pd.concat((merge_df, anno_df), axis=1)
       # print('merge_df')
        #print(merge_df.head())
        merge_df['Start'] = merge_df['Start'].fillna(-1).astype(int)
        merge_df['End'] = merge_df['End'].fillna(-1).astype(int)
        merge_df.to_csv(output.peaks_df, sep='\t')
        assert(len(merge_df)==len(anno_df))

# rule expression:
#     input:
#         merged="{genome}/merged/samples.merge",
#     output:
#         peak_expr='{genome}/expression/samples.expr.tsv',
#         out_cfg='{genome}/expression/.outcfg'
#     params:
#         raw_peaks=config['raw_peak'],
#         meta_f=config['meta_f']
#     shell: "python peak_sourceInfo.py run_peakexpr {input} {output.peak_expr} {params}"

def calc_scores(df, cpm_thresh, num_cols):
    new_cols = ["num_above_cpm", "isAbove_numsamples"]
    above_thresh = (df>cpm_thresh)
    df[new_cols[0]] = above_thresh.sum(axis=1)
    df[new_cols[1]] = (above_thresh.sum(axis=1))>num_cols#>params.sample_thresh]
    return df, new_cols, above_thresh


rule peak_scoring:
    #input:
    params:
        samples_expr="data/peaksInfo/samples.merge.peaksexpression",
        cpm_thresh=lambda wildcards: wildcards.thresh,
    output: 'data/peaksInfo/peaks_scoring.thresh{thresh}.tsv'
    run:
        df = pd.read_csv(params.samples_expr, sep='\t', index_col=0)
        print('df')
        print(df.head())
        #print('cpm thresh', )
        above_thresh = (df>int(params.cpm_thresh))
        df['num_above_cpm'] = above_thresh.sum(axis=1)
        df.index.name="ID"
        df[['num_above_cpm']].to_csv(output[0], sep='\t')


rule tissue_merging:
    params:
        tissues_expr="data/peaksInfo/tissues.merge.peaksexpression",
        cpm_thresh=lambda wildcards: wildcards.thresh,
        cho="CHO",
        collapse_bmdm=True
    output:
        tissue_merge="data/peaksInfo/tissuesMerge.thresh{thresh}.peaksexpression",
        tissue_score="data/peaksInfo/tissues_scoring.thresh{thresh}.tsv",
        tissue_thresh="data/peaksInfo/tissuesMerge.thresh{thresh}.bin.peaksexpression",
    run:
        tissues = pd.read_csv(params.tissues_expr, index_col=0, sep='\t')
        #merge bmdm
        if params.collapse_bmdm:
            bmdm_key = "BMDM"
            bmdm_cols = tissues.columns[tissues.columns.str.contains("BMDM", case=False)]
            tissues[bmdm_key] = tissues.loc[:, bmdm_cols].mean(axis=1)
            tissues=tissues.drop(bmdm_cols,axis=1)
        tissues.to_csv(output.tissue_merge, sep='\t')

        #tissues, new_cols, above_thresh = calc_scores(tissues, params.cpm_thresh, params.sample_thresh)
        above_thresh = (tissues>int(params.cpm_thresh))
        above_thresh.to_csv(output.tissue_thresh)
        tissues['num_above_cpm'] = above_thresh.sum(axis=1)
        out_df = tissues[['num_above_cpm']].copy()
        out_df.index.name="ID"

        # Add in CHO - Tissues
        cho_peaks = above_thresh.loc[above_thresh[params.cho]].index
        tissue_peaks = above_thresh.loc[above_thresh.drop(params.cho, axis=1).any(axis=1)].index
        all_peaks = above_thresh.loc[above_thresh.any(axis=1)].index

        out_df.loc[set(cho_peaks).intersection(set(tissue_peaks)), "isCHO"] = "Overlap"
        out_df.loc[set(cho_peaks)-(set(tissue_peaks)), "isCHO"] = "CHO"
        out_df.loc[set(tissue_peaks)-set(cho_peaks), "isCHO"] = "Tissues"
        out_df.loc[set(all_peaks)-(set(cho_peaks).union(set(tissue_peaks))), "isCHO"] = "None"
        out_df.to_csv(output.tissue_score, sep='\t')


rule generate_genome:
    params:
        gff=get_gff,
        pc="mRNA",
        ncrna=get_ncrna
    output: "data/{genome}/genome/refTSS.tsv"
    shell: "python peak_distance.py genome {output} {params}"


rule distance:
    #""" Distance to nearest TSS"""
    input:
        minimal_anno = "data/{genome}/anno/peaks.merge.tss.minimal.anno",
        refTSS = "data/{genome}/genome/refTSS.tsv"
    output:
        peaks_distance = "data/{genome}/anno/peaks.distance.tsv"
    params: cpu=32
    shell: "python peak_distance.py distance {input} {output} --cpu {params}"


rule get_introns:
    input: get_gff
    output: "data/{genome}/genome/introns.bed"
    script:
        "getIntrons.R"


rule filter_out_cds:
    input: "data/{genome}/anno/peaks.merge.tss.minimal.anno"
        #rules.distance.output,
    params: "data/{genome}/genome/CDS.bed"
    output: "data/{genome}/filters_genoebody/peaks.nocds.index"
    run: get_overlap(input[0], params[0], output[0])


rule filter_out_intron:
    input:
        "data/{genome}/anno/peaks.merge.tss.minimal.anno", #rules.distance.output,
        rules.get_introns.output,
    output: "data/{genome}/filters_genoebody/peaks.nointron.index"
    run: get_overlap(input[0], input[1], output[0])



# rule create_type:
#     input:
#         peaks_distance="data/{genome}/anno/peaks.distance.tsv"
#     output:
#         peaks_type="data/{genome}/tss_types/peaks_type.tsv"
#

rule create_refType:
    input:
        peaks_distance="data/{genome}/anno/peaks.distance.tsv",
        refTSS="data/{genome}/genome/refTSS.tsv",
        homeranno="data/{genome}/anno/peaks.merge.tss.minimal.anno"
    output: "data/{genome}/refType/peaks.refType.tsv"
    params:
        intergenic_dist=10000,
        bins="1000,1000",
        pc='mRNA',
        ncrna=get_ncrna
    shell: "python peak_refType.py {input} {output} {params}"



rule addGeneBody:
    input:
        rules.filter_out_cds.output,
        rules.filter_out_intron.output,
        "data/{genome}/refType/peaks.refType.tsv"
    output: "data/{genome}/refType/peaks.refType.geneBody.tsv"
    run:
        cds = pd.read_csv(input[0])['0'].values #the peaks that were NOT in cds
        intron = pd.read_csv(input[1])['0'].values #NOT in intron
        peaks = pd.read_csv(input[2],sep='\t', index_col=0)
        #print('peaks', peaks.head())

        peaks['isCDS'] = True
        peaks.loc[peaks.index.isin(cds), 'isCDS']=False
        peaks['isIntron'] = True
        peaks.loc[peaks.index.isin(intron), 'isIntron']=False
        peaks.to_csv(output[0], sep='\t')


rule create_totalTSS:
    input:
        scores='data/peaksInfo/peaks_scoring.thresh{thresh}.tsv',
        tissue_scores='data/peaksInfo/tissues_scoring.thresh{thresh}.tsv',
        refTypeAll="data/{genome}/refType/peaks.refType.geneBody.tsv",
    output:
        "data/{genome}/allTSS/thresh{thresh}_numSamp{numSamp}_peaksSankey.csv",
        report("data/{genome}/allTSS/thresh{thresh}_numSamp{numSamp}_peaksSankey.png"),
        allTSS="data/{genome}/allTSS/allTSS_thresh{thresh}_numSamp{numSamp}.tsv",
    params:
        keepCDS=False,
        keepIntron=False,
        #cpmThresh=config['cpm_thresh'],
        numSamplesThresh=lambda wildcards: wildcards.numSamp
    run:
        scores = pd.read_csv(input.scores,sep='\t')
        tissue_scores = pd.read_csv(input.tissue_scores,sep='\t')
        refType=pd.read_csv(input.refTypeAll,sep='\t')
        print(scores['num_above_cpm'].head())
        refType['isAbove_numsamples'] = scores['num_above_cpm']>int(params.numSamplesThresh)
        #refType['isAbove_numsamples'] = scores['isAbove_numsamples']
        refType['isCHO'] = tissue_scores['isCHO']

        # for f in filts:
        #     refType = refType.loc[refType[f]]
        ordered_cols = ["refType", "isCHO", "isAbove_numsamples", "isCDS", "isIntron"]
        sh.wrap_sankey(refType.fillna('Other'), ordered_cols, var_name=None, node_opacity=0.8,
                link_opacity=0.4, name="TSSs detected",
                add_colname=True, out_f=output[0].replace('.csv','.raw'))
        print('peaks shape', refType.shape)

        # Run with CDS and intron filter, vs run without
        ordered_cols = ["refType", "isCHO"]
        refType['isCHO'] = refType['isCHO'].fillna('None')

        refType = refType.loc[(refType['isAbove_numsamples'])]
        sh.wrap_sankey(refType, ordered_cols, var_name=None, node_opacity=0.8,
                link_opacity=0.4, name="TSSs detected",
                add_colname=True, out_f=output[0].replace('.csv','.keepBody'))

        refType = refType.loc[~(refType['isCDS'])]
        refType.to_csv(output.allTSS.replace('.tsv', '.rmCDS.tsv'))

        print('after cds', refType.shape)
        sh.wrap_sankey(refType, ordered_cols, var_name=None, node_opacity=0.8,
        link_opacity=0.4, name="TSSs detected",
        add_colname=True, out_f=output[0].replace('.csv','.rmCDS'))

        refType = refType.loc[~(refType['isIntron'])]
        print('after intron', refType.shape)
        print('after above numsamples', refType.shape)
        refType = refType.dropna(subset=['refType'],axis=0)
        print('after removing Nan from refType', refType.shape)
        sh.wrap_sankey(refType, ordered_cols, var_name=None, node_opacity=0.8,
                link_opacity=0.4, name="TSSs detected",
                add_colname=True, out_f=output[0].replace('.csv','.csv'))
        #refType.to_csv(output.allTSS_f, sep='\t')

        refType = refType["refType"].isin(['mRNA', 'Intergenic', 'ncRNA'],axis=0)
        sh.wrap_sankey(refType, ordered_cols, var_name=None, node_opacity=0.8,
            link_opacity=0.4, name="TSSs detected",
            add_colname=True, out_f=output[0].replace('.csv',''))


        refType.to_csv(output.allTSS_f)

    #shell: "python summarizePeaks.py {input} {output}"



def plot_genes_promoters_ref(promoters, refTSS, out_f, gene_key):
    promoters['Gene'] = promoters['Nearest TSS'].map(refTSS[gene_key].to_dict())
    genes_captured= set(promoters['Gene'])
    num_genes = len(genes_captured)
    prom_captured=set(promoters['Nearest TSS'])
    num_prom=len(prom_captured)
    f, ax = plt.subplots(nrows=1,ncols=2, dpi=300)
    plt.rc('figure', titlesize=8)
    sns.countplot(promoters.groupby('Gene').size(), ax=ax[0])
    ax[0].title.set_text(f"eTSS per RefSeq gene\nGenes captured={num_genes}")
    sns.countplot(promoters.groupby('Nearest TSS').size(), ax=ax[1])
    ax[1].title.set_text(f"eTSS per RefSeq promoter\nPromoters captured={num_prom}")
    plt.subplots_adjust()
    hs(out_f)
    return


rule refPromoters_captured:
    input:
        allTSS=rules.create_totalTSS.output['allTSS'],
        refTSS="data/{genome}/genome/refTSS.tsv"
    output:
        pc="data/{genome}/allTSS/SF4_pc_refPromoters_thresh{thresh}_numSamp{numSamp}.png",
        ncrna="data/{genome}/allTSS/SF4_ncrna_refPromoters_thresh{thresh}_numSamp{numSamp}.png",
        #pc_tiss="data/{genome}/F1/pc_genes_over_tissues.png",
    params:
        gene_key=lambda wildcards: config['gene_key'][wildcards.genome]
    run:
        df = pd.read_csv(input.allTSS, index_col=0)
        refTSS = pd.read_csv(input.refTSS, index_col=0, sep='\t')
        pc = (df.copy())[df['refType'] == 'mRNA']
        plot_genes_promoters_ref(pc, refTSS, output.pc)
        ncrna = df[df['refType'] == "ncRNA"]
        #print('num ncrna', len(ncrna))
        plot_genes_promoters_ref(ncrna, refTSS, output.ncrna, params.gene_key)


rule refPromoters_tissues_captured:
    input:
        allTSS=rules.create_totalTSS.output.allTSS,
        tissue_merge="data/peaksInfo/tissuesMerge.thresh{thresh}.peaksexpression",
        refTSS="data/{genome}/genome/refTSS.tsv"
    params:
        cpm_thresh = lambda wildcards: wildcards.thresh,
        gene_key=lambda wildcards: config['gene_key'][wildcards.genome],
        numSamp=lambda wildcards: wildcards.numSamp
    output:
        "data/{genome}/allTSS/F1_tissuesCaptured_thresh{thresh}_numSamp{numSamp}.png"
    run:
        refTSS = pd.read_csv(input.refTSS, index_col=0, sep='\t')
        #df = pd.read_csv(input.allTSS, sep='\t', index_col=0)
        df = pd.read_csv(input.allTSS, index_col=0)
        promoters = (df.copy())[df['refType'] == 'mRNA']
        promoters['Gene'] = promoters['Nearest TSS'].map(refTSS[params.gene_key].to_dict())
        genes_captured= set(promoters['Gene'])
        tissues = pd.read_csv(input.tissue_merge,index_col=0,sep='\t')
        # Filter for the promoter TSSs
        tissues = tissues.loc[promoters.index]
        pc_bin = tissues>int(params.cpm_thresh)
        #gene_bin = tissues>int(params.cpm_thresh)
        print('tissues shape', tissues.shape)
        gene_bin = pd.DataFrame(columns=pc_bin.columns, index=genes_captured)
        gene_to_peaks = {x:promoters[promoters['Gene']==x].index.values for x in genes_captured}
        def wrap_prom_to_gene(gene_t):
            def prom_to_gene_bin(tiss, gene):
                return pc_bin.loc[gene_to_peaks[gene], tiss].any()
            gene_t.loc[:] = gene_t.index
            return gene_t.apply(prom_to_gene_bin, args=(gene_t.name,))

        gene_bin = gene_bin.apply(wrap_prom_to_gene, axis=1)
        print('gene_bin')
        print(gene_bin.head())
        gene_bin.to_csv(output[0]+'.geneBin.tsv',sep='\t')
        cdf_number_gene,gene_tissues_order = create_cdf(gene_bin, start_tissues=("CHO",))
        prom_peaks = {}
        for t in gene_tissues_order:
            col = t
            prom_peaks[t] = set(gene_bin[~(gene_bin[col]==0)].index.values)

        prom_cum_peaks_n, prom_cum_peaks = get_cum_and_plot(prom_peaks, gene_tissues_order);
        mrna_peaks = {}
        for t in gene_tissues_order:
            col = t
            mrna_peaks[t] = set(gene_bin[~(gene_bin[col]==0)].index.values)
        mrna_cum_peaks_n, mrna_cum_peaks = get_cum_and_plot(mrna_peaks, gene_tissues_order,f_save=output[0])
        with open(output[0]+'.csv','w') as f:
            f.write('\t'.join(mrna_cum_peaks_n)+'\n'+'\t'.join(gene_tissues_order))


rule retrieve_peaks:
    input:
        gene_centric = "data/processed/{anno}/gene_centric_tss/files_used.json",
        merged_f = "data/processed/{anno}/merged/samples.merge",
        expr_peaks_f = "data/processed/{anno}/merged/samples.merge.peaksexpression.log10",
        peaks_dir = config["peaks_dir"],
    output: "data/processed/{anno}/retrieve_peaks/files_used.json"
    shell: "python tss/data/retrieve_tss_peaks.py {output} {input}"
    #log: 'data/processed/{anno}/retrieve_peaks/retrieve_peaks.log'


# rule revisedTSS:
#     input:
#         rules.create_totalTSS.output[-2]
#     output:
#         "data/{genome}/revTSS/pc_promoters.bed",
#         "data/{genome}/revTSS/pc_promoters.meta.tsv",
#         "data/{genome}/revTSS/ncrna_promoters.bed",
#         "data/{genome}/revTSS/ncrna_promoters.meta.tsv"
#     params:
#         outdir=lambda wildcards, output: os.path.dirname(output[0]),
#         center_on="" # minimumDist, eTSS, mostSignificant. eTSS uses CHO if found, otherwise mostSignificant
#
#     shell: "python {input} {outdir}"
#     input:
#     params:

#rule create_ncrna:

#     shell:
#
#

#
rule stableRNA_merge:
    input: rules.combine_merge_and_anno.output
    output:
        raw_merge="data/stable/stable.merged"
    params:
        stableRNA="data/peakInfo/csrna/stable.tissues",
        dist='given'
    shell:"mergePeaks -d {params.dist} -strand + {input} {params} -venn {output.raw_merge}.venn > {output.raw_merge}"

rule stableRNA_clean:
    input:
        rules.combine_merge_and_anno.output,
        rules.stableRNA_merge.output
    output: peaks_stable="data/stable/peaks_stable.tsv"
    run:
        ### For overlapping TSS
        minimal = pd.read_csv(input[0], sep='\t', index_col=0)
        merged_df = pd.read_csv(input[0], sep='\t', index_col=0)
        #merged_df = merged_df[~(merged_df.isnull().any(axis=1))]
        merged_df.to_csv(output[0], sep='\t')
        #Homer.merge_peaks(input_files, output_file, dist='given', type_merge=''):









