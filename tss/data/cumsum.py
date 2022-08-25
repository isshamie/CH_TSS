# Create Plots for Figure 1. Breaks things up by type of RNAs.
import pandas as pd
import tqdm
import numpy as np
import matplotlib.pyplot as plt
from mplh.fig_utils import helper_save, legend_from_color
from mplh.color_utils import get_colors


def to_collapse(curr, collapse_mapping):
    collapse = set()
    for i in curr:
        collapse.add(collapse_mapping[i])
    return list(collapse)


def create_binary(df, tissues, colname="tissues"):
    df = df[~(df[colname] == "")]
    bin_df = pd.DataFrame(index=df.index, columns=tissues, dtype=bool)
    bin_df[:] = False
    for ind, val in tqdm.tqdm(df.iterrows()):
        # print(val)
        if val[colname] == "":
            print(val)
        else:
            bin_df.loc[ind, val[colname].split(',')] = True

    return bin_df


## Create CDF of genes
def create_cdf(bin_df, start_tissues=()):
    ## Steps:
    # 1. Create dictionary of tissue combo to number of genes
    # 1. Double for loop and get the pairwise union overlap
    # 2. Pop them from list
    # 3. While not done with all columns
    # -Loop through and find the next maximum overlap
    # -Add to dictionary
    # -Pop from list
    all_cols = list(bin_df.columns.values)
    if len(start_tissues)>0: #Predefined initial list
        cdf_tissues = list(start_tissues)
        #cdf_number = dict(
        cdf_number = []
        overlap = set()
        for t in start_tissues:
            overlap = set((bin_df[bin_df[t]].index)).union(overlap)
            cdf_number.append(len(overlap))
            all_cols.remove(t)
        max_overlap = overlap
        max_num = cdf_number[-1]
    else:
        cdf_tissues = []
        cdf_number = []

        max_overlap = 0
        max_combo = ()
        max_num = 0
        for i in all_cols:
            for j in all_cols[1:]:
                overlap = set((bin_df[bin_df[i]].index)).union(
                    set((bin_df[bin_df[j]].index)))
                if len(overlap) > max_num:
                    max_overlap = overlap
                    max_num = len(overlap)
                    max_combo = (i, j)

        cdf_tissues.append(max_combo[np.argmax(
            [bin_df[max_combo[0]].sum(), bin_df[max_combo[1]].sum()])])
        cdf_tissues.append(max_combo[np.argmin(
            [bin_df[max_combo[0]].sum(), bin_df[max_combo[1]].sum()])])

        cdf_number.append(
            max(bin_df[max_combo[0]].sum(), bin_df[max_combo[1]].sum()))
        cdf_number.append(max_num)

        all_cols.remove(max_combo[0])
        all_cols.remove(max_combo[1])

    while len(all_cols) > 0:
        no_addition = True
        for i in all_cols:
            overlap = max_overlap.union(set((bin_df[bin_df[i]].index)))
            if len(overlap) > max_num:
                max_num = len(overlap)
                max_ind = i
                no_addition = False
        if no_addition:
            break

        max_overlap = max_overlap.union(
            set((bin_df[bin_df[max_ind]].index)))
        cdf_tissues.append(max_ind)
        cdf_number.append(max_num)
        all_cols.remove(max_ind)

    return cdf_number, cdf_tissues


# def create_cdf_with_list(bin_df, columns):
#     all_vals = set()
#     for i in columns:
#         overlap = max_overlap.union(set((bin_df[bin_df[i]].index)))
#         if len(overlap) > max_num:
#             max_num = len(overlap)
#             max_ind = i
#             no_addition = False
#
#     return cdf_number, cdf_tissues


def plot_cum(cum_peaks_n, tissues_order, f_save=None, ax=None, c=None):
    if ax is None:
        f, ax = plt.subplots()
    if c is not None:
        ax.scatter(np.arange(0, len(cum_peaks_n)), cum_peaks_n, c=c)
    else:
        ax.scatter(np.arange(0, len(cum_peaks_n)), cum_peaks_n)
    plt.xticks(np.arange(0, len(cum_peaks_n)))
    ax.set_xticklabels(tissues_order);
    ax.tick_params(axis="x", rotation=45)

    if f_save is not None:
        helper_save(f_save)
    return


def get_cum_and_plot(peaks_dict, tissues_order, f_save=None, to_plot=True):
    """
    peaks_dict: Dictionary where each k is an element in tissues_order and value is a set of peak IDs
    tissues_order: List thats ordered with the keys of peaks_dict
    """
    cum_peaks_n = []
    cum_peaks = set()
    for t in tissues_order:
        print(len(cum_peaks))
        cum_peaks = cum_peaks.union(peaks_dict[t])
        cum_peaks_n.append(len(cum_peaks))

    if to_plot:
        plot_cum(cum_peaks_n, tissues_order, f_save)
    return cum_peaks_n, cum_peaks


## Create binary of expression matrix
def create_binary_from_expr(expr_df, threshold=1):
    bin_df = (expr_df > threshold)
    return bin_df



def plot_cumsum(all_peaks_n, gene_tissues_order, cumulative_f):
    df = pd.DataFrame(all_peaks_n,
                      columns=["None"] + gene_tissues_order)
    colors, name_dict = get_colors(scheme='categorical', n_colors=3,
                                   names=df.index.values)
    colors
    f, ax = plt.subplots()
    for ind, val in df.iterrows():
        plot_cum(val.values, df.columns, f_save=None, ax=ax, c=colors[
            ind])  # plt.scatter(np.arange(len(val.values)val.values)
    legend_from_color(colors, curr_ax=ax)
    helper_save(cumulative_f)
    return

