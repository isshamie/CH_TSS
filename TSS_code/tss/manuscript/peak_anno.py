import click
import pandas as pd
from os.path import join
import numpy as np
from tss.data import generate_genome
import pickle
from tss.manuscript.utils import check_strand
from pandarallel import pandarallel

def proteincoding():
    return


def intergenic(peaks, peaks_expr, dist_to_tss_thresh, cpm_thresh):
    intergenic = peaks[(peaks["Annotation"].str.contains("Intergenic") & (np.abs(
                        peaks["Distance to TSS"].fillna(0).astype(
                            int)) > dist_to_tss_thresh))]
    return intergenic.index


def ncrna(annoPeaks, genome_anno_df, expr_df):
    anno_types = genome_anno_df[2].unique()
    not_include = ['mRNA', 'gene', 'region', 'exon', 'CDS',
                   'primary_transcript', 'guide_RNA', 'pseudogene',
                   'transcript', 'match', 'cDNA_match']
    anno_types_keep = set(anno_types) - set(not_include)

    print(anno_types_keep)
    anno_df_otherrna = genome_anno_df[genome_anno_df[2].isin(anno_types_keep)]
    print(anno_df_otherrna.shape)
    anno_df_otherrna = generate_genome.expand_anno_id(anno_df_otherrna,
                                                      colname=8)

    # Rename first keys:
    anno_df_otherrna = anno_df_otherrna.rename(
        {0: "Chr", 3: "Start", 4: "End", 6: "Strand"}, axis=1)
    anno_df_otherrna = anno_df_otherrna.set_index("ID", drop=False)

    ncrnas = annoPeaks[
        annoPeaks['Annotation'].fillna('').str.contains("ncRNA")]
    ncrnas_prom = ncrnas[
        (ncrnas["Annotation"].fillna('').str.contains("promoter"))]
    ncrnas_prom = ncrnas_prom.dropna(axis=1, how='all')
    ncrnas_prom = generate_genome.expand_anno_id(ncrnas_prom,
                                                 colname='Nearest PromoterID')

    ncrnas_prom["isSameStrand"] = ncrnas_prom.apply(check_strand, args=(
    anno_df_otherrna,), axis=1)
    # remove the non-same strands
    ncrnas_prom = ncrnas_prom[ncrnas_prom["isSameStrand"]]

    ncrnas_prom_expression = expr_df.loc[ncrnas_prom.index]
    print(f"Number of ncrnas: {ncrnas.shape[0]}")
    ncrnas_prom_expression = expr_df.loc[ncrnas_prom.index]
    gene_ID = dict()
    for biotype, df in anno_df_otherrna.groupby(2):
        print(biotype)
        gene_ID[biotype] = df['ID'].values

    peak_ID = dict()
    for i in gene_ID:
        print(i)
        if i == "V_gene_segment" or i == "C_gene_segment":
            continue
        print(f"Number of the type: {len(gene_ID[i])}")
        print('Parent')
        print(ncrnas_prom['Parent'].isin(gene_ID[i]).sum())

        peak_ID[i] = ncrnas_prom[
            ncrnas_prom['Parent'].isin(gene_ID[i])].index

    smallncRNAs = set()
    allncRNAs = set()
    for ind in peak_ID:
        allncRNAs.update(peak_ID[ind])
        if ind != 'lnc_RNA':
            smallncRNAs.update(peak_ID[ind])

    curr_ncrna_anno = annoPeaks.loc[allncRNAs]
    curr_ncrna_anno['Absolute'] = np.absolute(
        curr_ncrna_anno["Distance to TSS"])

    ncbi_ncrna = (curr_ncrna_anno[curr_ncrna_anno['Absolute'] ==
                                  curr_ncrna_anno.groupby(
                                      'Nearest PromoterID')[
                                      'Absolute'].transform('min')])

    allncRNAs_df = annoPeaks.loc[allncRNAs]
    allncRNAs_df['Absolute'] = np.absolute(annoPeaks["Distance to TSS"])
    print(len(allncRNAs_df["Nearest PromoterID"].values))
    print("Number of promoters",
          len(set(allncRNAs_df["Nearest PromoterID"].values)))
    return



def run(tissues_expression_f, annoPeaks_f, params):
    expr_df = pd.read_csv(tissues_expression_f, sep="\t", index_col=0)
    annoPeaks = pd.read_csv(annoPeaks_f, sep="\t", index_col=0)
    annoPeaks = annoPeaks.dropna(subset=['Annotation'])
    return



@click.command()
@click.argument('peak_expression_in', type=click.Path(exists=True))
@click.argument('config', type=click.Path(exists=True))
def main(peak_expression_in, config):
    return
