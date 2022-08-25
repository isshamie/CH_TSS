import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

matplotlib.use('Agg')


#################################################
### Plot helpers
#################################################
def helper_save(f_save, remove_bg=True):
    """ Function to save as png and svg if f_save is not None."""
    if remove_bg:
        f = plt.gcf()
        f.patch.set_facecolor("w")
        axs = f.get_axes()
        for ax in axs:
            ax.set_facecolor('white')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)

    if f_save is not None:
        #Remove suffixes
        f_save = f_save.replace('.csv', '')
        f_save = f_save.replace('.tsv', '')
        f_save = f_save.replace('.txt', '')

        if '.png' in f_save:
            plt.savefig(f_save, bbox_inches="tight", transparent=True)
            plt.savefig(f_save.split('.png')[0] + '.svg')
        else:
            plt.savefig(f_save + '.png',bbox_inches='tight', transparent=True)
            plt.savefig(f_save + '.svg')
            plt.savefig(f_save + '.pdf')
        return


def determine_rows_cols(num_samples):
    nrows = 2
    ncols = 2
    while nrows * ncols <= num_samples:
        ncols = ncols + 1
        if nrows * ncols <= num_samples:
            nrows = nrows + 1

    return nrows,ncols


#######################################
def wrap_plots(group,func, f_save, names=None, *args):

    xlim = [np.infty, -np.infty]
    ylim = [np.infty, -np.infty]

    num_samples = len(group)
    nrows,ncols = determine_rows_cols(num_samples)

    f = plt.figure()
    axs = []

    for ind, fname in enumerate(group):
        axs.append(plt.subplot(nrows, ncols, ind + 1))

        func(fname, fname + "_nucl", *args,f=f)

        # heat_plot(fname, save_f=heat_save, f=f, curr_ax=axs[ind],
        #           num_peaks=num_peaks, is_norm=is_norm)
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
    helper_save(f_save)

    return

#
#
# def stack_bars(x, ys, ylabels, f_save=None, use_text=True,num_dec=0):
#     # stack bars
#     f = plt.figure()
#     for ind, curr_y in enumerate(ys):
#         plt.bar(x,curr_y , bottom= ys[:ind].sum(axis=0), label=ylabels[ind])
#
#     # add text annotation corresponding to the percentage of each data.
#     if use_text:
#         for xpos, ypos, yval in zip(x, ys[0]/2, ys[0]):
#             plt.text(xpos, ypos, "%.1f"%(yval*100), ha="center", va="center")
#
#         for ind, val in enumerate(ys):
#             if ind == 0:
#                 continue
#             for xpos, ypos, yval in zip(x,ys[:ind].sum(axis=0)+val/2, val):
#                 plt.text(xpos, ypos, f"%.{num_dec}f"%(yval*100), ha="center", va="center")
#
#     #plt.ylim(0,110)
#     plt.title("Normalized stacked barplot of different RNA types")
#     plt.legend(bbox_to_anchor=(1.01,0.5), loc='center left')
#     if f_save is not None:
#         helper_save(f_save)
#     plt.xticks(rotation=90)
#     #plt.savefig('normalized_stacked_barplot_with_number.png', bbox_inches='tight', pad_inches=0.02)
#     return
#

#
# def show_values_on_bars(axs):
#     def _show_on_single_plot(ax):
#         for p in ax.patches:
#             _x = p.get_x() + p.get_width() / 2
#             _y = p.get_y() + p.get_height()
#             value = '{:.2f}'.format(p.get_height())
#             ax.text(_x, _y, value, ha="center")
#
#     if isinstance(axs, np.ndarray):
#         for idx, ax in np.ndenumerate(axs):
#             _show_on_single_plot(ax)
#     else:
#         _show_on_single_plot(axs)
