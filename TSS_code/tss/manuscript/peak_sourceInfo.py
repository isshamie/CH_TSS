import tss.utils.Homer as Homer
from tss.config import HOMER_PATH as HP
from tss.data import annotation
import glob
import json
import os
import numpy as np
import pandas as pd
from collections import defaultdict
from os.path import dirname
import click
from tss.manuscript.utils import merge_columns

#all_peak_merge_files = glob.glob("f04_peaks/*")
#print(all_peak_merge_files)
from tss.data.data_io import read_bed_file

# def merge_peaks(input_dir, output_file):
#     print(input_dir, output_file)
#     return

def merge_sample_peaks(raw_peaks_dir, peak_f, anno_f, seq_f, dist,
                       ref_fa, ref_anno):
    """
    @param p: The parameters dictionary
    :param outdir:
    :param peaks_dir:
    :param dist:
    :return:
    """
    # Run merging
    print(raw_peaks_dir)
    all_peak_merge_files = glob.glob(raw_peaks_dir + "/*")
    Homer.merge_peaks(all_peak_merge_files, peak_f, dist=dist)
    # Need to clean up the merged file. Do this by removing the
    # folder names, so it is just the basename that is kept.
    with open(peak_f) as f:
        data = f.read()
    data = data.replace(raw_peaks_dir, "")
    with open(peak_f, "w") as f:
        f.write(data)

    Homer.annotate_peaks(peak_file=peak_f, output_file=anno_f,
                         ref_fa=ref_fa, ref_anno=ref_anno)


    #seq_f = anno_f + ".fa"
    cmd = 'homerTools extract {f_in} {gen} -fa > {f_out}'.format(
        f_in=peak_f, gen=ref_fa, f_out=seq_f)
    print(cmd)
    os.system(cmd)

    return



def add_stablerna():
    return


def extract_source_peaks():
    return


def create_tissues():
    return


def create_tissue_peaks_id():
    return


def numSamples():
    return


def numTissues():
    return




def create_tissue_to_peaks(meta_f, peaks_expr_f, tissues_expr_f, merge_bmdm=True):

    meta_samples = pd.read_csv(meta_f, sep="\t",index_col=0)
    #print(meta_samples.head())
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




@click.group()
@click.command()
# @click.argument('outdir', type=click.Path(exists=True))
# @click.pass_context
def cli():
    return


# def combine_merge_with_anno(merge_file, anno_file):
#     """ Function that adds the column Annotation from anno_file to the merge_file"""
#     merge_df = pd.read_csv(merge_file, index_col=0, sep='\t')
#     anno_df = pd.read_csv(anno_file, sep='\t', na_filter=False, index_col=0)
#     anno_df.index.name = "ID"
#     anno_df = anno_df['Annotation']
#     merge_df = pd.concat((merge_df, anno_df), axis=1)
#     merge_df.to_csv(merge_file+'.anno',sep='\t')
#     return merge_df


@cli.command()
@click.argument('homerAnno')
@click.argument('minimalAnno')
@click.argument('peaks')
def run_minimal(homerAnno, minimalAnno, peaks):
    # This should be in the merging part before really
    print("Combining homers merge and anno")
    #combine_merge_with_anno(minimalAnno, anno_file=homerAnno)
    print("Creating size-1 peak and converting to bed.")
    #Put peaks to size 1
    peaks_df = pd.read_csv(minimalAnno, sep="\t")
    peaks_df["Start"] = (
                1.0 * (peaks_df["End"] + peaks_df["Start"]) / 2).astype(
        int)
    peaks_df["End"] = peaks_df["Start"] #length 1
    peaks_df["actual_start"] = peaks_df["End"]
    peaks_df.to_csv(peaks, sep="\t", index=False)

    # Convert to bed
    cmd = f"pos2bed.pl {peaks} > {peaks.replace('.peak','.bed')}"
    print(cmd)
    os.system(cmd)
    return


def run_peaks_start(in_f, out_f):
    peaks_df = pd.read_csv(in_f, sep="\t")
    peaks_df["Start"] = (
            1.0 * (peaks_df["End"] + peaks_df["Start"]) / 2).astype(int)
    peaks_df["End"] = peaks_df["Start"]  # length 1
    peaks_df["actual_start"] = peaks_df["End"]
    peaks_df.to_csv(out_f, sep="\t", index=False)
    return

@cli.command()
@click.argument('merged_f')
@click.argument('raw_peak_dir')
@click.argument('out_peak_expr_f')
@click.argument('meta_f')
def run_peakexpr(merged_f, out_peak_expr_f, raw_peak_dir, meta_f):
    # # Convert into matrix
    print("Converting to expression matrix")
    annotation.convert_merged_vals_to_expression_matrix(
        merged_f, peak_folder=raw_peak_dir,
        output_file=out_peak_expr_f)
    print("Taking the log(x+1) of the expression")
    #Log
    peaks_df = pd.read_csv(out_peak_expr_f, sep="\t", index_col=0)
    peaks_df_log10 = np.log10(peaks_df + 1)
    peaks_df_log10.to_csv(out_peak_expr_f + '.log10', sep="\t")
    peaks_df_log2 = np.log2(peaks_df + 1)
    peaks_df_log2.to_csv(out_peak_expr_f + '.log2', sep="\t")

    print("Merging over replicates")
    create_tissue_to_peaks(meta_f, out_peak_expr_f, out_peak_expr_f+".tissues")
    create_tissue_to_peaks(meta_f, out_peak_expr_f+'.log2',
                           out_peak_expr_f + ".log2.tissues")
    create_tissue_to_peaks(meta_f, out_peak_expr_f+'.log10',
                           out_peak_expr_f + ".log10.tissues")
    return

# CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
# @click.command(context_settings=CONTEXT_SETTINGS)
# @click.argument('start_site_mrna', type=click.Path(exists=True))
# @click.argument('ref_anno', type=click.Path(exists=True))
# @click.argument('no_cds_peak', type=click.Path(exists=False))
# @click.argument('peaks_dir', type=click.Path(exists=True))
# @click.argument('ref_fa', type=click.Path(exists=True))
# @click.argument('mode', type=click.STRING)
# def main(start_site_mrna, ref_anno, no_cds_peak, peaks_dir, ref_fa, mode):
#
#     merged_dir = dirname(no_cds_peak)
#     genome_dir = dirname(start_site_mrna)
#
#     p = pipeline.create_filenames_dict()
#     p = pipeline.create_fullnames(p,'merged', merged_dir)
#     p = pipeline.create_fullnames(p, 'reference', genome_dir)
#     run(p, peaks_dir, ref_fa, ref_anno, mode)
#
#     # Save the parameters
#     with open(p["merged"]["files_used"], 'w') as fp:
#         json.dump(p, fp)
#
#     return
#
#
# if __name__ == "__main__":
#     main()