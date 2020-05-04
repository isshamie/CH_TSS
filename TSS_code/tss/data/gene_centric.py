from tss.data import annotation
import numpy as np
import pandas as pd
from collections import defaultdict
import os
import json
import pickle


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


def run(p):
    # Set variable names
    #
    p_merged_f = p["merged"]["filenames"]
    p_stage = p["gene_centric"]
    p_stage_p = p["gene_centric"]["params"]
    p_stage_f = p["gene_centric"]["filenames"]
    p_reference = p["reference"]
    p_global = p["global"]


    # Global and reference variables
    ref_fa = p_reference["GENOME_FA"]
    annotation = p_reference["GENOME_GFF3"]
    cds_f = p_reference["CDS_F"]
    tss_f = p_reference["TSS_F"]

    meta_f = p_global["META_FILE"]

    # Prior Filenames
    #peaks_f = p_merged_f["merged peak samples"]
    peaks_f = p_merged_f["no cds tsv"]

    peaks_expression_f = p_merged_f["peak samples expression"]
    #peaks_with_tss_distances_f = p_merged_f["peak samples distance
    # to TSS"]
    #peaks_with_tss_distances_size1_f = p_merged_f["peak samples
    # distance to TSS size 1nt"]
    #peaks_with_tss_distances_size1_noCDS_f = p_merged_f["no cds tsv"]


    # params
    peak_bins = p_stage_p["peak_bins"]
    allow_intron = p_stage_p["allow_intron"]

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


    # Save the parameters
    with open(os.path.join(p_stage["folder"],'params_used.json'), 'w') as fp:
        json.dump(p, fp)


    return
