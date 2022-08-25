from tss.visualize.fig_utils import wrap_plots
#from tss.utils.Homer import homer_nucleotide
import matplotlib.pyplot as plt
from os.path import join, basename
import numpy as np
from tss.utils import helper
import sarge
import pandas as pd
from tss.utils import create_output
from mplh.fig_utils import helper_save
import os
from tss import pipeline
import logging
import json
import click
from os.path import dirname
#from tss.pipeline import create_fullnames, create_filenames_dict
from tss import pipeline


# #######################################
def wrap_plots(group,func, f_save, names=None, *args):
    xlim = [np.infty, -np.infty]
    ylim = [np.infty, -np.infty]

    num_samples = len(group)
    nrows,ncols = helper.determine_rows_cols(num_samples)
    if len(group) == 3: #Hardcode just 3 columns
        nrows, ncols = 1, 3

    f = plt.figure(dpi=300)
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
            curr_label = basename(fname)
        else:
            curr_label = names[ind]
        axs[ind].set_title(curr_label)
        if ind == len(group)-1:
            axs[ind].legend(bbox_to_anchor=(1.04,1),loc="upper left")

    [ax.set_xlim(xlim) for ax in axs]
    [ax.set_ylim(ylim) for ax in axs]
    helper_save(f_save)
    return


#######################################
def homer_nucleotide(input_file, output_file, ref_fa, size=1000,
                     desired_lims=None, only_plot=False, f=None):
    """
    Input
    * input_file: peak file
    * output_file where to save the nucleotide frequency plot
    * ref_fa: Reference annotation

    Will also plot the nuc freqs and save to {output_file}.png
    """

    # cmd = 'annotatePeaks.pl f06_annoPeaks/merge_bg_2.anno_promoter /data/genome/hamster/picr/picr.fa -size 1000 -hist 1 -di > nuc_freq.txt'
    if not only_plot:
        cmd = 'annotatePeaks.pl %s %s -size %s -hist 1 -di > %s' % (
            input_file, ref_fa, size, output_file)
        print(cmd)
        sarge.run(cmd)
        tmp = pd.read_csv(output_file, sep='\t', index_col=0)
    else:  # Alrady have the frequency plots, just want to plot
        tmp = pd.read_csv(input_file, sep='\t', index_col=0)
    if f is None:
        f = plt.figure(dpi=300)
        ax = f.add_subplot(111)
        to_return = False
        to_legend = True
    else:
        ax = f.gca()
        to_return = True
        to_legend = False
    ax.plot(tmp.index.values, tmp['C frequency'], linewidth=2)
    ax.plot(tmp.index.values, tmp['A frequency'], linewidth=2,alpha=0.8)
    ax.plot(tmp.index.values, tmp['T frequency'], linewidth=2,alpha=0.6)
    ax.plot(tmp.index.values, tmp['G frequency'], linewidth=2,alpha=0.4)
    if desired_lims is None or desired_lims == ():
        ax.vlines(0, ax.get_ylim()[0], ax.get_ylim()[1], linewidth=0.5)
    else:
        ax.vlines(0, desired_lims[0], desired_lims[1], linewidth=0.5)
    if to_legend:
        ax.legend(loc="upper left")
    ax.set_xlabel('bp from TSS')
    ax.set_ylabel('Frequency')

    helper_save(output_file)

    if to_return:
        return f
    else:
        return
    #plt.savefig(output_file + '.png')


def run(p, ref_fa):
    #
    # f1 = "Results/output/TSS1.exp.bed"
    # f2 = "Results/Figures/Figure3/alt/A.TSS1_mrna"
    #
    # f = plt.figure()
    # wrap_plots([f1,f2],homer_nucleotide,"Results/Figures/Figure3/C.combine",["Experimental","RefSeq"],ref_fa,1000,(0.1,0.65),False)
    #
    # f1 = "Results/output/TSS1.exp.bed"
    # f2 = "Results/Figures/Figure3/alt/A.TSS1_mrna"
    #
    # f = plt.figure()
    # wrap_plots([f1,f2],homer_nucleotide,"Results/Figures/Figure3/C_200.combine",["Experimental","RefSeq"],ref_fa,200,(0.1,0.65),False)

    eTSS_bed = p["output_peaks"]["eTSS exp bed"]
    refcenter_eTSS_bed = p["output_peaks"]["eTSS refseq_centered bed"]

    # meta_f = p["output_peaks"]["eTSS exp meta"]
    # mRNA_peak_file = p["reference"]["mRNA.peak"]
    #plt.figure()
    wrap_plots([eTSS_bed,refcenter_eTSS_bed],homer_nucleotide,
               p["nucleotide_out"]["compare_old_new"],["Experimental","RefSeq"],ref_fa,200,(0.1,0.65),False)


    # Refseq
    ## Create RefSeq peak file for this  # eTSS expand the nucleotides by 75bp both sides:
    #shift=75 #150
    # create_output.exp_bed_to_refseq(eTSS_bed,meta_f,refseq_f=mRNA_peak_file,save_f=refcenter_eTSS_bed,is_unique=True,
    #                                shift=shift)

    # ### Extract sequences- needed? It will with motifs
    # bed_df = read_bed_file(peak_f)
    # bed_df["Start"] -= shift
    # bed_df["End"] += shift
    # peak_f = join(nucleotide_out,f"TSS1_{shift}bp.exp.bed")
    # write_bed_file(bed_df,bed_f=peak_f)

    # ## eTSS
    # eTSS_seq_f = p["nucleotide_out"]["eTSS exp seq"]
    # cmd = "homerTools extract {peak_f} {ref_fa} -fa > {seq_f}".format(peak_f=eTSS_bed,
    #                                                                   seq_f=eTSS_seq_f,
    #                                                                   ref_fa=ref_fa)
    # os.system(cmd)
    # print(cmd)
    # ref_centered_seq_f = p["nucleotide_out"]["eTSS refseq_centered seq"] # mrna_filt + ".fa"
    # cmd = "homerTools extract {peak_f} {ref_fa} -fa > {seq_f}".format(peak_f=refcenter_eTSS_bed,
    #                                                                   seq_f=ref_centered_seq_f,ref_fa=ref_fa)
    # os.system(cmd)
    # print(cmd)


    # homer_nucleotide(input_file, output_file, ref_fa, size=1000,
#                      desired_lims=None, only_plot=False, f=None)


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('tss_bed_f', type=click.Path(exists=True))
@click.argument('nucleotides_out', type=click.Path(exists=False))
@click.argument('ref_fa', type=click.Path(exists=True))
def main(tss_bed_f, nucleotides_out, ref_fa):
    p = pipeline.create_filenames_dict()
    p = pipeline.create_fullnames(p, 'output_peaks', dirname(tss_bed_f))
    p = pipeline.create_fullnames(p, 'nucleotide_out', dirname(nucleotides_out))
    if p == {}:
        print("The directories were not setup accordingly for the pipeline")
    run(p, ref_fa=ref_fa)
    return


if __name__ == '__main__':
    main()

