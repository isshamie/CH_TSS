import tss.utils.Homer as Homer
from tss.config import HOMER_PATH as HP
from tss.data import annotation
import glob
import json
import os
import numpy as np
import pandas as pd
#all_peak_merge_files = glob.glob("f04_peaks/*")
#print(all_peak_merge_files)


# def merge_peaks(input_dir, output_file):
#     print(input_dir, output_file)
#     return

def merge_sample_peaks(raw_peaks_dir, peak_f, anno_f, seq_f, dist,
                       ref_fa, annotation_f):
    """
    @param p: The parameters dictionary
    :param outdir:
    :param peaks_dir:
    :param dist:
    :return:
    """
    # Run merging
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
                         ref_fa=ref_fa, annotation=annotation_f)


    #seq_f = anno_f + ".fa"
    cmd = 'homerTools extract {f_in} {gen} -fa > {f_out}'.format(
        f_in=peak_f, gen=ref_fa, f_out=seq_f)
    print(cmd)
    os.system(cmd)

    return


def merge_tags(DATA_DIR, tags_cap_f, tags_input_f, TISSUES):
    ### Parameters setup
    """
    :param outdir:
    :param data_folder:
    :param tissues:
    """

    # p_global = p["global"]
    # p_merged = p["merged"]
    # DATA_DIR = p_global["DATA_DIR"]

    ## Cap Tags
    all_tags_start = []
    for t in TISSUES:
        for s in ['GROCap', 'START']:
            # curr = glob.glob(j + '*f03_tags/trim*.peak')
            curr_raw = glob.glob(DATA_DIR + t + '/*/')

            # print(curr_raw)
            for j in curr_raw:
                curr_type = j.split('/')[-2]
                if curr_type == 'GROCap':

                    new = glob.glob(os.path.join(DATA_DIR, s, j,
                                                 'f03_tags/*GROCap*'))
                elif curr_type == 'START':
                    new = np.array(glob.glob(
                        os.path.join(DATA_DIR, s, j, 'f03_tags/*')))
                    if not len(new) == 0:
                        new = (new[(np.array(
                            ['input' not in x for x in new]))])  # Remove the input
                else:
                    continue
                all_tags_start.extend(new)

    curr_out = tags_cap_f

    cmd = '{HP}/makeTagDirectory {out_f} -single -d {files}'.format(
           HP=HP,out_f=curr_out, files=' '.join(all_tags_start))
    print(cmd)
    os.system(cmd)

    ## Input Tags
    all_tags_input = []
    for t in TISSUES:
        for s in ['GROCap','START']:
            #curr = glob.glob(j + '*f03_tags/trim*.peak')
            j = glob.glob(os.path.join(DATA_DIR, t, s))
            if len(j) > 0:
                j = j[0]
            else:
                continue
            curr_type = s
            if curr_type == 'GROCap':
                #print(os.path.join(j,'f03_tags/*5GRO*'))
                new = np.array(glob.glob(os.path.join(DATA_DIR,s,j,'f03_tags/*GRO*')))
                if not len(new) == 0:
                    new = (new[(np.array(['Cap' not in x for x in new]))]) #Remove the input
            elif curr_type == 'START':
                new = np.array(glob.glob(os.path.join(DATA_DIR,s,j,'f03_tags/*input*')))
            else:
                continue
            all_tags_input.extend(new)

    curr_out = tags_input_f
    cmd = '{HP}makeTagDirectory {out_f} -single -d {files}'.format(
            HP=HP,out_f=curr_out, files=' '.join(all_tags_input))
    print(cmd)
    os.system(cmd)
    return


# def wrap_merging(p):
#     merge_sample_peaks(p)
#     merge_tags(p)
#     return


def wrap_peaks_breakdown(merged_f, anno_f, tss_f, minimal_f,
                         minimal_anno_f, peaks_with_tss_distances_f,
                         peaks_expression_f, raw_peak_dir, peaks_bed,
                         peaks_with_tss_distances_size1_f):
    """
    Input files
    merged_f: is the input peaks file,
    anno_f: TSS with Homer annotations
    tss_f: landmark file

    Output files
    peaks_expression_f: n_peaks-by-n_samples matrix
    minimal_f, minimal_anno_f: Condensed columns
    peaks_with_tss_distances_f: Cleaned up peaks with tss distances. This is done using the center of the peaks.
    Logged values of the peak expression matrices: As stated
    peaks_with_tss_distances_size1: 2 files,  bed format and the Homer peaks format that take the center of the peak value.
        peaks format
    """

    # # Convert into matrix
    annotation.convert_merged_vals_to_expression_matrix(
        merged_f, peak_folder=raw_peak_dir,
        output_file=peaks_expression_f)

    annotation.combine_merge_with_anno(minimal_f,
                                       anno_file=anno_f)

    annotation.wrap_distance_to_landmarks(minimal_anno_f,
                                          landmark_file=tss_f,
                                          main_landmark='transcript_id',
                                          landmark_cols=['gene'],
                                          output_f=peaks_with_tss_distances_f,
                                          is_parallel=True,
                                          num_par=36)

    #Put peaks to size 1
    peaks_df = pd.read_csv(peaks_with_tss_distances_f, sep="\t")
    peaks_df["Start"] = (
                1.0 * (peaks_df["End"] + peaks_df["Start"]) / 2).astype(
        int)
    peaks_df["End"] = peaks_df["Start"] + 1
    peaks_df.to_csv(peaks_with_tss_distances_size1_f, sep="\t")

    peaks_df[["Chr", "Start", "End", "ID", "Stat", "Strand"]].to_csv(
        peaks_bed, header=None, sep="\t", index=False)

    #Log
    peaks_expression_file = peaks_expression_f
    peaks_df = pd.read_csv(peaks_expression_file, sep="\t", index_col=0)
    peaks_df_log10 = np.log10(peaks_df + 1)
    peaks_df_log10.to_csv(peaks_expression_file + '.log10', sep="\t")

    peaks_df_log2 = np.log10(peaks_df + 1)
    peaks_df_log2.to_csv(peaks_expression_file + '.log2', sep="\t")

    return


def cds_overlap(peaks_bed, cds_f, peaks_with_tss_distances_size1_f,
                peaks_with_tss_distances_size1_noCDS_f, no_cds_tsv_f):
    """
    Function that uses bedtools intersect to not have any peaks that
    overlap with Coding sequence since that would be harder to call a start site (even if alternative start sites).
    Will make bed file of this and tsv file.
    :param peaks_bed:
    :param cds_f:
    :param peaks_with_tss_distances_size1_f:
    :param peaks_with_tss_distances_size1_noCDS_f:
    :param no_cds_tsv_f:
    :return:
    """
    cmd = "bedtools intersect -s -v -a {peaks_bed} -b {cds}  > {out_f}".format(
        peaks_bed=peaks_bed, cds=cds_f, out_f=peaks_with_tss_distances_size1_noCDS_f)
    print(cmd)
    os.system(cmd)

    peaks_df = pd.read_csv(peaks_with_tss_distances_size1_f,
                           sep="\t", index_col=0)

    peaks_df[["Chr", "Start", "End", "ID", "Stat", "Strand"]].to_csv(
        peaks_bed, header=None, sep="\t", index=False)
    out_f = peaks_with_tss_distances_size1_noCDS_f

    nocds = pd.read_csv(out_f, sep="\t", header=None)
    peaks_df[peaks_df["ID"].isin(nocds[3])].to_csv(
        no_cds_tsv_f, index=False, sep="\t")
    return


def run(p):
    # Set variable names
    #

    p_stage = p["merged"]
    p_stage_f = p_stage["filenames"]
    p_reference = p["reference"]
    p_global = p["global"]

    # Prior inputs
    raw_peaks = p["move_peaks"]["folder"]


    # Global variables
    DATA_DIR = p_global["DATA_DIR"]
    TISSUES = p_global["TISSUES"]
    ref_fa = p_reference["GENOME_FA"]
    annotation = p_reference["GENOME_GFF3"]
    cds_f = p_reference["CDS_F"]
    tss_f = p_reference["TSS_F"]

    # params
    dist = p_stage["params"]["dist"]

    # Filenames
    merged_f = p_stage_f["merged peak samples"]
    anno_f = p_stage_f["merged peak samples annotation"]
    seq_f = p_stage_f["merged peak samples sequences"]

    tags_cap_f = p_stage_f["merged tags cap"]
    tags_input_f = p_stage_f["merged tags input"]

    peaks_expression_f = p_stage_f["peak samples expression"]
    minimal_f = p_stage_f["peak samples minimal"]
    minimal_anno_f = p_stage_f["peak samples minimal anno"]

    peaks_with_tss_distances_f = p_stage_f["peak samples distance to TSS"]
    peaks_with_tss_distances_size1_f = p_stage_f["peak samples distance to TSS size 1nt"]
    peaks_with_tss_distances_size1_noCDS_f = p_stage_f["peaks samples " \
                                               "distance to TSS size 1nt no CDS"]
    no_cds_tsv_f = p_stage_f["no cds tsv"]

    peaks_bed = p_stage_f["peaks_bed"]


    # Run Functions
    merge_sample_peaks(raw_peaks, merged_f, anno_f, seq_f, dist, ref_fa,
                       annotation)
    #merge_tags(DATA_DIR, tags_cap_f, tags_input_f, TISSUES)

    wrap_peaks_breakdown(merged_f, anno_f, tss_f, minimal_f,
                        minimal_anno_f, peaks_with_tss_distances_f,
                        peaks_expression_f, raw_peak_dir=raw_peaks,
                        peaks_bed=peaks_bed,
                        peaks_with_tss_distances_size1_f=peaks_with_tss_distances_size1_f)

    cds_overlap(peaks_bed, cds_f, peaks_with_tss_distances_size1_f,
                peaks_with_tss_distances_size1_noCDS_f, no_cds_tsv_f)

    # Save the parameters
    with open(os.path.join(p_stage["folder"],'params_used.json'), 'w') as fp:
        json.dump(p, fp)
    return

