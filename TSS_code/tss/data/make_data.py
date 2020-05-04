## Functions to process data
import pandas as pd
import os
import numpy as np
import tss.utils.Homer as Homer
import sarge
import sys


def make_tag_directory(in_bam,tag_dir,ref_fa):
    '''make tag directory which extracts mapping position into tsv file
    '''
    cmd = ('makeTagDirectory {o_dir} -genome {g} -checkGC \
            -single {bam}').format(o_dir=tag_dir,g=ref_fa,bam=in_bam)
    print(cmd);sys.stdout.flush()
    sarge.run(cmd)


########################################
def make_bedgraph_file(in_dir, out_files):
    ''' Create bedgraph file for experiment. - strand will have negative values to be used'''
    cmd = (
        'makeUCSCfile {in_dir} -o {out_file} -strand + -fragLength 1').format(
        in_dir=in_dir, out_file=out_files[0])
    print(cmd);
    sys.stdout.flush()
    sarge.run(cmd)

    cmd = (
        'makeUCSCfile {in_dir} -o {out_file} -strand - -fragLength 1').format(
        in_dir=in_dir, out_file=out_files[1])
    print(cmd);
    sys.stdout.flush()
    sarge.run(cmd)
    print('here')
    print(out_files)
    cmd = "gunzip -c %s.gz | awk '{$4=-$4; print}' | gzip > %s_tmp.gz" % (
    out_files[1], out_files[
        1])  # .format(out_f_neg=out_files[1]) #Change - strand to have neg values
    print(cmd)
    sarge.run(cmd)
    out = out_files[0].replace('_pos',
                               '')  # trim_CHO--mSTART-JHS823_S21_R1_001_pos.bedgraph.gz

    cmd = ("cat {out_f_0}.gz {out_f_1}_tmp.gz > {final}.gz").format(
        out_f_0=out_files[0], out_f_1=out_files[1],
        final=out)  # Concatenate the 2 into new merged file
    print(cmd)
    sarge.run(cmd)
    cmd = ("rm {out_f_0}.gz {out_f_1}_tmp.gz {out_f_1}.gz").format(
        out_f_0=out_files[0],
        out_f_1=out_files[1])  # replace pos and neg
    print(cmd)
    sarge.run(cmd)
    return


#@requires(f03_tags)
#@requires(f04_peaks)
def convert_merged_vals_to_expression_matrix(merged_file, peak_folder='.', output_file=None):
    """Function that takes in a merged peak file and converts into a peak-by-sample expression matrix.
    Really just removes some of the meta-info associated with the merged peaks."""

    peaks_merged = pd.read_csv(merged_file, sep='\t', index_col=0,
                               low_memory=False)

    ## Drop duplicates
    dups = peaks_merged[peaks_merged.index.duplicated(keep=False)]
    print('Number of duplicate indices (having same index name, happens bc both +/- strand on exact same bp): ', len(
        dups) / 2)
    peaks_merged = peaks_merged[~peaks_merged.index.isin(dups.index)]

    peak_tissue_matrix = pd.DataFrame(
        np.zeros([peaks_merged.shape[0], len(peaks_merged.columns.values[7:])], dtype='int64'),
        columns=peaks_merged.columns.values[7:])
    peak_tissue_matrix.set_index(peaks_merged.index.values, inplace=True)
    for col in peak_tissue_matrix.columns.values:
        print((os.path.join(peak_folder, col)))
        try:
            orig_peaks_df = Homer.read_peak_file(os.path.join(peak_folder, col))
            # curr_peaks = pd.read_csv(os.path.join(peak_folder,col), sep='\t', index_col=0)
        except IOError:
            print('File couldnt be loaded!')
            continue

        # Load in the column file and get the peak values
        filt = peaks_merged.loc[~peaks_merged.loc[:, col].isnull()]  # Get rows where the column had a peak
        curr_new_inds = filt.index.values # Get the indices for all-indices
        orig_id = np.array(filt[col]) # Get the names of peaks from the original file

        # There are cases where the peaks from original file got merged into the same peak when merging.
        # In this case, add artifical peaks with the id as in the merge file, and take the average of the
        # Normalized Tag Count as the value
        for i in orig_id[[',' in x for x in orig_id]]:
            split_peaks = i.split(',')
            orig_peaks_df.loc[i, 'Normalized Tag Count'] = np.mean(
                orig_peaks_df.loc[split_peaks, 'Normalized Tag Count'])

        # Filter orig_peaks_df to only with keys.
        if len(orig_peaks_df.index.intersection(orig_id)) != len(orig_id):
            print('The labels from the merge is different from the original. Not adding this file to the matrix')
            continue

        orig_peaks_df = orig_peaks_df.reindex(orig_id)
        vals = np.array(orig_peaks_df.loc[:, 'Normalized Tag Count'])
        ##
        # Normalized Tag Count instead of p-value here.
        # vals = np.array((curr_peaks.loc[curr_keys, 'p-value vs Control'] +
        #                  curr_peaks.loc[curr_keys, 'p-value vs Local'])/2)
        ##
        vals = np.nan_to_num(vals)  # Turn nan values into 0
        peak_tissue_matrix.loc[curr_new_inds, col] = vals
        # peak_tissue_matrix.loc[curr_AI,col] = (curr_peaks.loc[curr_keys,'Stat'])

    # In case any peaks have no tissues with the expression, which shouldn't really happen
    peak_tissue_matrix = peak_tissue_matrix.loc[:, (peak_tissue_matrix > 0).sum() != 0]
    peak_tissue_matrix = peak_tissue_matrix.loc[np.sum(peak_tissue_matrix > 0, axis=1) != 0, :]

    if output_file is not None:
        peak_tissue_matrix.to_csv(output_file,sep='\t')
        peaks_merged = peaks_merged[['Chr', 'Start', 'End', 'Strand', 'Stat']]
        peaks_merged.to_csv(merged_file + '.minimal', sep='\t')

    return peak_tissue_matrix


def combine_merge_with_anno(merge_file,anno_file):
    """ Function that adds the column Annotation
    from anno_file to the merge_file"""
    merge_df = pd.read_csv(merge_file, index_col=0, sep='\t')
    anno_df = pd.read_csv(anno_file, sep='\t', na_filter=False, index_col=0)
    anno_df.index.name = "ID"
    anno_df = anno_df['Annotation']
    merge_df = pd.concat((merge_df, anno_df), axis=1)
    merge_df.to_csv(merge_file+'.anno',sep='\t')
    return merge_df

