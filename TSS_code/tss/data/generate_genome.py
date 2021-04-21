import os
from os.path import dirname
from os.path import join
from tss.utils.Homer import *
from tqdm import tqdm
import click
import logging


def process_genome(ref_fa, annotation, genome_dir, with_descriptives=False):
    # genome_raw_dir = os.path.join(doc["supplemental"], "genome",
    #                               "ncbi_anno_103")
    # ref_fa = doc["ref_fa"]
    # annotation = doc["annotation"]
    # ref_fa = os.path.join(genome_raw_dir,
    #                       'GCF_003668045.1_CriGri-PICR_genomic.fna')
    # annotation = os.path.join(genome_raw_dir,
    #                           'alt_CriGri-PICR_top_level.gff3')

    ref_start_f = join(genome_dir, 'start_site_mRNA.tsv')
    mrna_peak = join(genome_dir, 'mRNA.peak')
    exon_peak = join(genome_dir, 'exon.peak')
    cds_peak = join(genome_dir, 'CDS.peak')

    cds_bed = join(genome_dir, 'CDS.bed')
    mrna_bed = join(genome_dir, 'mRNA.bed')
    mrna_peak_bed = join(genome_dir, 'mRNA.peak.bed')
    cds_gff = join(genome_dir, 'CDS.gff')
    mrna_gff = join(genome_dir, 'mRNA.gff')

    all_genes_f = join(genome_dir, "all_genes.tsv")
    #ref_start_homer_f = join(genome_dir, 'start_site_mRNA_homer.tsv')

    createPeakFileFromGFF(annotation, output_file=mrna_peak,
        anno_of_interest='mRNA', is_start=True)

    createPeakFileFromGFF(annotation, output_file=exon_peak,
        anno_of_interest='exon', is_start=True)

    createPeakFileFromGFF(annotation, output_file=cds_peak,
        anno_of_interest='CDS', is_start=True)


    # expand the keys into columns and save
    mrna_peak_df = pd.read_csv(mrna_peak, sep="\t", index_col=0)
    for ind in mrna_peak_df.index.values:
        curr = ind.split(';')
        for i in curr:
            curr_split = i.split('=')
            mrna_peak_df.set_value(ind, curr_split[0], curr_split[1])

    mrna_peak_df.set_index("transcript_id", inplace=True,drop=False)
    mrna_peak_df.index.drop_duplicates()  # There is one duplicate transcript
    mrna_peak_df['Length'] = mrna_peak_df['End'] - mrna_peak_df['Start'] + 1
    mrna_peak_df.to_csv(mrna_peak, sep="\t") #, index="transcript_id")


    ### Create CDS and mRNA gff and bed file
    cmd = "awk '{ if ($3==\"mRNA\") {print}  }' %s > %s" % (annotation, mrna_gff)
    print(cmd)
    os.system(cmd)

    cmd = "awk '{ if ($3==\"CDS\") {print}  }' %s > %s" % (annotation, cds_gff)
    print(cmd)
    os.system(cmd)

    cmd = "gff2bed < {cds_gff} > {cds_bed}".format(cds_gff=cds_gff,
                                                   cds_bed=cds_bed)
    print(cmd)
    os.system(cmd)

    cmd = "gff2bed < {mrna_gff} > {mrna_bed}".format(mrna_gff=mrna_gff,
                                                     mrna_bed=mrna_bed)
    print(cmd)
    os.system(cmd)

    cmd = f"pos2bed.pl {mrna_peak} > {mrna_peak_bed}"
    print(cmd)
    os.system(cmd)


    annotation_start_site = pd.DataFrame(
        columns=mrna_peak_df.columns.values)
    for group, start_site in tqdm(
            mrna_peak_df.groupby(['actual_start', 'Chr', 'gene'])):
        if len(start_site) == 1:
            annotation_start_site = annotation_start_site.append(
                start_site)
        else:
            # Pick the most confident txn and the longest transcript
            annotation_start_site = annotation_start_site.append(
                start_site.sort_values(['Length'],
                                       ascending=[False]).iloc[0])
    annotation_start_site.sort_values(['Chr', 'Start'], inplace=True)

    annotation_start_site['Start'] = annotation_start_site[
        'Start'].astype(int)
    annotation_start_site['End'] = annotation_start_site['End'].astype(
        int)

    # Drop any transcript_id duplicates. When I checked with GCF, there was 1
    print(f"Number of duplicate txn IDs: {(annotation_start_site.index.duplicated()).sum()}")
    annotation_start_site = annotation_start_site[~(annotation_start_site.index.duplicated())]

    #annotation_start_site["transcript_id"] = annotation_start_site.index
    annotation_start_site.to_csv(ref_start_f, sep='\t')

    # Index and then get the chromosomes from genome
    cmd = "samtools faidx {ref_fa}".format(ref_fa=ref_fa)
    os.system(cmd)
    print(cmd)
    cmd = "cut -f1,2 {ref_fa_ind} > {chrom}".format(
        ref_fa_ind=ref_fa + ".fai",
        chrom=join(genome_dir, "chrom.sizes"))
    print(cmd)
    os.system(cmd)

    with open(all_genes_f, 'w') as f:
        genes = list(pd.read_csv(mrna_peak, sep="\t").groupby(
            "gene").groups.keys())
        genes.sort()
        f.write("\n".join(genes))

    if with_descriptives:
        genome_ann_mrna = pd.read_csv(mrna_peak, sep='\t', index_col=0)
        genome_ann_exon = pd.read_csv(exon_peak, sep='\t', index_col=0)

        #### Number of unique gene_names
        print('Number of unique genes: ', len(set(genome_ann_mrna["gene"])))
        #### Number of unique start sites for mrna
        ss = set()
        for ind, val in genome_ann_mrna.iterrows():
            ss.add(
                str(val['Start']) + '_' + val['Strand'] + '_' + val['Chr'])
        print('Number of unique start sites for mrna', len(ss))

        # print('Number of transcripts per gene: ',1.0*len(genome_ann_mrna)/len(unique_genes))

        print('Number of transcripts', len(genome_ann_mrna))
        print('Number of exons', len(genome_ann_exon))
    return


def expand_anno_id(df, break_char="=",colname=8):
    df = df.copy()
    for ind, val in tqdm(df.iterrows()):
        curr = val[colname].split(';')

        for i in curr:

            v = i.strip().replace('"', "")
            if len(v) == 0:
                continue
            curr_split = v.split(break_char)

            df.at[ind, curr_split[0]] = curr_split[1]

    return df


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('ref_fa', type=click.Path(exists=True))
@click.argument('annotation', type=click.Path(exists=True))
@click.argument('mrna_peak', type=click.Path(exists=False))
def main(ref_fa, annotation, mrna_peak):
    genome_dir = dirname(mrna_peak)
    process_genome(ref_fa, annotation, genome_dir,
                   with_descriptives=False)
    return


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    main()

    # not used in this stub but often useful for finding various files

    # project_dir = Path(__file__).resolve().parents[2]
    # print(project_dir)
    #
    # #load_dotenv(find_dotenv())
    # main(['/data/isshamie/TSS/TSS/parameters/ncbi_anno_103.yaml'])
    # #main()

