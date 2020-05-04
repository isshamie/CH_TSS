import tss.data.stats as stats
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from tss.visualize.fig_utils import helper_save, determine_rows_cols
import os


########################################
def hist_plot(hist_out,include_norm=True, f=None, to_save=True,
              to_fwhm=True):
    """Visualize histograms created by hist."""

    if f is None:
        f = plt.figure()
    df = pd.read_csv(hist_out,sep='\t',header=0,names=['Distance from TSS','Coverage','+ Tags','- Tags'])
    if to_fwhm:
        max_val = np.max(df['+ Tags'])
        val = stats.fwhm(df["Distance from TSS"], df['+ Tags'], k=3)
        print(('Max value: {max_val}'.format(max_val=max_val)))
        print(("Full-width at half-maximum: {val} (nts)".format(val=val)))

    plt.plot(df['Distance from TSS'],df['+ Tags'], label='+ Tags')
    plt.plot(df['Distance from TSS'],df['- Tags'], label='- Tags')
    plt.xlim([-500,500])
    plt.xlabel('Distance from TSS')
    plt.ylabel('Reads per bp per TSS')
    plt.axvline(x=0,c='k')
    plt.legend(loc='upper right')

    if to_save:
        plt.savefig(hist_out+'.png')

    if include_norm:
        plt.figure()
        df = pd.read_csv(hist_out+'Norm',sep='\t',header=0,names=['Distance from TSS','Coverage','+ Tags','- Tags'])
        plt.plot(df['Distance from TSS'],df['+ Tags'],label='+ Tags')
        plt.plot(df['Distance from TSS'],df['- Tags'],label='- Tags')
        plt.xlim([-500,500])
        plt.xlabel('Distance from TSS (normalized per read)')
        plt.ylabel('Reads per bp per TSS')
        plt.axvline(x=0,c='k')
        plt.legend(loc='upper right')
        #plt.savefig(os.path.splitext(hist_out)[0]+'Norm.png')
        if to_save:
            plt.savefig(hist_out+'Norm.png')


def hist_plot_norm(hist_out, f=None,to_save=True,
              to_fwhm=True):
    df = pd.read_csv(hist_out + 'Norm', sep='\t', header=0,
                     names=['Distance from TSS', 'Coverage', '+ Tags',
                            '- Tags'])
    if to_fwhm:
        max_val = np.max(df['+ Tags'])
        val = stats.fwhm(df["Distance from TSS"], df['+ Tags'], k=3)
        print(('Max value: {max_val}'.format(max_val=max_val)))
        print(("Full-width at half-maximum: {val} (nts)".format(val=val)))

    if f is None:
        f = plt.figure()
    plt.plot(df['Distance from TSS'], df['+ Tags'], label='+ Tags')
    plt.plot(df['Distance from TSS'], df['- Tags'], label='- Tags')
    plt.xlim([-500, 500])
    plt.xlabel('Distance from TSS (normalized per read)')
    plt.ylabel('Reads per bp per TSS')
    plt.axvline(x=0, c='k')
    plt.legend(loc='upper right')
    # plt.savefig(os.path.splitext(hist_out)[0]+'Norm.png')
    if to_save:
        plt.savefig(hist_out + 'Norm.png')


########################################
def wrap_hist_plot(hist_outs, hist_save=None, names=None,
                   to_norm=False):
    """
    Takes multipled histogram files and plots in subplots, keeping
    the limits the same.

    :param
    hist_outs: List of histogram files to plot.
               The basenames will be used for the respective titles.
    hist_save: File name to save to. If None, will not save figure.
    :return: None
    """
    xlim = [np.infty, -np.infty]
    ylim = [np.infty, -np.infty]

    num_samples = len(hist_outs)
    nrows,ncols = determine_rows_cols(num_samples)
    f = plt.figure()
    axs = []
    for ind,fname in enumerate(hist_outs):
        axs.append(plt.subplot(nrows,ncols, ind+1))
        if to_norm:
            hist_plot_norm(fname, f=f, to_save=False)
        else:
            hist_plot(fname, include_norm=False, f=f, to_save=False)
        xlim[0] = min(axs[ind].get_xlim()[0], xlim[0])
        ylim[0] = min(axs[ind].get_ylim()[0], ylim[0])
        xlim[1] = max(axs[ind].get_xlim()[1], xlim[1])
        ylim[1] = max(axs[ind].get_ylim()[1], ylim[1])
        if names is None:
            curr_label = os.path.basename(fname)
        else:
            curr_label = names[ind]
        axs[ind].set_title(curr_label)

    [ax.set_xlim(xlim) for ax in axs]
    [ax.set_ylim(ylim) for ax in axs]
    #[ax.set_facecolor('white') for ax in axs]
    helper_save(hist_save)


def wrap_heat_plot(heat_outs,heat_save=None, names=None,
                   num_peaks=1000,is_norm=True):
    """
    Takes multipled heat-map files and plots in subplots, keeping
    the limits the same.

    :param
    heat_outs: List of heat matrices files to plot.
               The basenames will be used for the respective titles.
    hist_save: File name to save to. If None, will not save figure.
    :return: None
    """
    xlim = [np.infty, -np.infty]
    ylim = [np.infty, -np.infty]

    num_samples = len(heat_outs)
    nrows,ncols = determine_rows_cols(num_samples)
    f = plt.figure()
    axs = []
    for ind,fname in enumerate(heat_outs):
        axs.append(plt.subplot(nrows,ncols, ind+1))
        heat_plot(fname, save_f=heat_save,f=f,curr_ax=axs[ind],
                  num_peaks=num_peaks, is_norm=is_norm)
        xlim[0] = min(axs[ind].get_xlim()[0], xlim[0])
        ylim[0] = min(axs[ind].get_ylim()[0], ylim[0])
        xlim[1] = max(axs[ind].get_xlim()[1], xlim[1])
        ylim[1] = max(axs[ind].get_ylim()[1], ylim[1])
        if names is None:
            curr_label = os.path.basename(fname)
        else:
            curr_label = names[ind]
        axs[ind].set_title(curr_label)

    [ax.set_xlim(xlim) for ax in axs]
    [ax.set_ylim(ylim) for ax in axs]
    helper_save(heat_save)


########################################
def heat_plot(heat_file, sort_bins=(-1, 1), num_peaks=1000, \
              is_norm=True, save_f='', f=None, curr_ax=None):
    """
    Creates a heat plot with the input coming from annotatePeaks -ghist
    """
    if f is None:
        plt.figure()
    heat_df = pd.read_csv(heat_file, sep='\t', index_col=0)
    centr = int(np.floor(heat_df.shape[1] / 2))
    # print(np.sum(heat_df.iloc[:,centr-sort_bins[0]:centr+sort_bins[1]+1],axis=1))
    # print(heat_df.iloc[:,centr-sort_bins[0]:centr+sort_bins[1]+1])
    heat_df = heat_df.iloc[:min(num_peaks, heat_df.shape[0])]
    heat_df = heat_df.iloc[np.argsort(np.sum(
        heat_df.iloc[:, centr - sort_bins[0]:centr + sort_bins[1] + 1],
        axis=1))[::-1]]

    if is_norm:
        heat_df = heat_df.divide(np.sum(heat_df, axis=1), axis='index')

    if curr_ax is None:
        sns.heatmap(heat_df, robust=True, xticklabels=4, yticklabels=False)
    else:
        sns.heatmap(heat_df, robust=True, xticklabels=4,
                    yticklabels=False,ax=curr_ax)
    helper_save(save_f)

    return heat_df

