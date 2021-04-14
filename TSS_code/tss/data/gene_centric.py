from tss.data import annotation
import numpy as np
import pandas as pd
from collections import defaultdict
import os
import json
import pickle
import click
from os.path import dirname
from tss import pipeline


def add_tissue_info(gene_centric_f, peaks_expr_f, meta_f, df_tissues_f):
    # Add in tissues
    ## invert and create a peaks_to_tissue
    peaks_tissue = pd.read_csv(peaks_expr_f, index_col=0, sep="\t")
    meta_samples = pd.read_csv(meta_f, sep="\t", index_col=0)

    tissue_to_peaks = defaultdict(list)
    for val in peaks_tissue.columns.values:
        count = 0
        for t in np.unique(meta_samples["Tissue"].values):
            if t in val:
                count += 1
                tissue_to_peaks[t].append(val)
                if count > 1:
                    print(
                        "This sample has multiple tissues associated with it")
                    print("Sample: ", val)

    peaks_to_tissue = dict()
    for t in tissue_to_peaks:
        for p in tissue_to_peaks[t]:
            peaks_to_tissue[os.path.basename(p)] = t

    df = pickle.load(open(gene_centric_f, "rb"))
    df["Tissues"] = ""
    for ind, val in df.iterrows():
        if len(val["samples"]) != 0:
            tis = set()
            for t in val["samples"]:
                tis.add(peaks_to_tissue[t])
            df.at[ind, "Tissues"] = ",".join(tis)

    pickle.dump(df, open(df_tissues_f, "wb"))
    return


def wrap_gene_centric(peaks_f, peaks_expression_f, gene_df_f,tss_f,
                      txn_df_f, peak_bins, allow_intron,
                      txn_expression_f):
    print("gene-based")
    annotation.wrap_create_anno_centric(peaks_f,
                                        peaks_expression_f,
                                        tss_f,
                                        peak_bin=peak_bins,
                                        anno_col='Nearest gene',
                                        tss_df_col='gene',
                                        f_save=gene_df_f,
                                        allow_intron=False)
    print("transcript-based")
    annotation.wrap_create_anno_centric(peaks_f,
                                        peaks_expression_f,
                                        tss_f,
                                        peak_bin=peak_bins,
                                        anno_col='Nearest TSS',
                                        tss_df_col='transcript_id',
                                        f_save=txn_df_f,
                                        allow_intron=allow_intron)
    print("Converting to expression matrix")
    annotation.df_to_TSS_expression(txn_df_f, peaks_expression_f,
                                    f_out=txn_expression_f)

    return


def run(p, meta_f, peak_bins, allow_intron):
    # Set variable names
    #
    p_merged_f = p["merged"]
    p_stage_f = p["gene_centric"]
    p_reference = p["reference"]


    # Global and reference variables
    tss_f = p_reference["start_site_mRNA"]
    peaks_f = p_merged_f["no cds peak"]
    peaks_expression_f = p_merged_f["peak samples expression log2"]

    print('tss_f', tss_f)
    print('peaks_f', peaks_f)
    print('peaks_expression_f', peaks_expression_f)
    # Current filenames
    gene_df_f = p_stage_f["gene centric peaks from samples"]
    txn_df_f = p_stage_f["transcript centric peaks from samples"]
    txn_expression_f = p_stage_f["transcript centric peaks matrix"]

    # Run analysis
    wrap_gene_centric(peaks_f, peaks_expression_f, gene_df_f, tss_f,
                      txn_df_f, peak_bins, allow_intron,
                      txn_expression_f)


    gene_tissues_f = p_stage_f["gene centric peaks from samples plus " \
                          "tissues"]
    add_tissue_info(gene_df_f, peaks_expression_f, meta_f, gene_tissues_f)

    txn_tissues_f = p_stage_f["transcript centric peaks from samples plus tissues"]
    add_tissue_info(txn_df_f, peaks_expression_f, meta_f, txn_tissues_f)


    return



CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('txn_promoter_f', type=click.Path(exists=False))
@click.argument('merged_out', type=click.Path(exists=True))
@click.argument('start_site_mrna', type=click.Path(exists=True))
@click.argument('meta_f', type=click.Path(exists=True))
@click.argument('allow_intron', type=click.BOOL)
@click.argument('peak_bin_l', type=click.INT)
@click.argument('peak_bin_r', type=click.INT)
def main(txn_promoter_f, merged_out, start_site_mrna, meta_f, allow_intron, peak_bin_l, peak_bin_r):

    merged_dir = dirname(merged_out)
    gene_centric_dir = dirname(txn_promoter_f)
    ref_dir = dirname(start_site_mrna)

    p = pipeline.create_filenames_dict()
    p = pipeline.create_fullnames(p,'gene_centric', gene_centric_dir)
    p = pipeline.create_fullnames(p, 'reference', ref_dir)
    p = pipeline.create_fullnames(p, 'merged', merged_dir)

    run(p, meta_f, (peak_bin_l, peak_bin_r), allow_intron)

    # Save the parameters
    with open(p["gene_centric"]["files_used"], 'w') as fp:
        json.dump(p, fp)

    return

#
if __name__ == "__main__":
    main()