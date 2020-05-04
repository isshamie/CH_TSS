import pandas as pd
from collections import defaultdict
import numpy as np
import os
import json


def merge_columns(df, mapping_dict):
    '''Function that merges columns by taking the mean of them.
    Merges based on if they have the same element in the meta_samples meta_column.
    Returns:
        new_df: Dataframe but with the columns of interest merged. Also column names are now based
                on the unique meta_samples[meta_column].
    '''

    vals = mapping_dict.keys()#np.unique(meta_samples[meta_column].values)
    new_df = pd.DataFrame(index=df.index, columns=vals)
    for i in vals:
        if not mapping_dict[i] == []:#(meta_samples[meta_column] == i).any():
            new_col = (df.loc[:, mapping_dict[i]])
            new_col = new_col.mean(axis=1)
            new_df.loc[:, i] = new_col
    new_df = new_df.loc[:, ~(new_df.isnull().all())]
    return new_df


def create_tissue_to_peaks(meta_f, peaks_expr_f, tissues_expr_f):
    meta_samples = pd.read_csv(meta_f, sep="\t",index_col=0)

    peaks_tissue = pd.read_csv(peaks_expr_f, index_col=0, sep="\t")
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
    peaks_tissue_merged = merge_columns(peaks_tissue,
                                        mapping_dict=tissue_to_peaks)
    peaks_tissue_merged.to_csv(tissues_expr_f, sep="\t")
    return


def run(p):
    # Filenames
    p_stage = p["tissues_collapse"]
    p_merged_f = p["merged"]["filenames"]
    meta_f = p["global"]["META_FILE"]
    peak_samples_expression_log10_f = p_merged_f["peak samples expression log10"]
    tissues_expr_log10_f = p_stage["filenames"]["peak tissues expression log10"]

    peak_samples_expression_f = p_merged_f["peak samples expression"]
    tissues_expr_f = p_stage["filenames"]["peak tissues expression"]


    # Run functions
    create_tissue_to_peaks(meta_f, peak_samples_expression_log10_f, tissues_expr_log10_f)
    create_tissue_to_peaks(meta_f, peak_samples_expression_f, tissues_expr_f)


    peaks_tissue = pd.read_csv(peak_samples_expression_log10_f,
        index_col=0, sep="\t")

    # Save the parameters
    with open(os.path.join(p_stage["folder"],'params_used.json'), 'w') as fp:
        json.dump(p, fp)
    return
