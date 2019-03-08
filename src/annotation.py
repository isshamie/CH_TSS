import yaml
import sys
import sys
import pandas as pd
import matplotlib
import seaborn as sns
import pickle
import numpy as np
from matplotlib import pyplot as plt
from itertools import product
import glob
import re
import inspect
import tqdm
import os
import Homer
from cycler import cycler
sys.path.append("/home/isshamie/software/homebrew/parallel_functions/")
import parallel_functions



def convert_merged_vals_to_expression_matrix(merged_file, peak_folder='.', output_file=None):
    """Function that takes in a merged peak file and converts into a peak-by-sample expression matrix.
    Really just removes some of the meta-info associated with the merged peaks."""

    peaks_merged = pd.read_csv(merged_file, sep='\t', index_col=0,low_memory=False)

    ## Drop duplicates
    dups = peaks_merged[peaks_merged.index.duplicated(keep=False)]
    print 'Number of duplicate indices (having same index name, happens bc both +/- strand on exact same bp): ', len(
        dups) / 2
    peaks_merged = peaks_merged[~peaks_merged.index.isin(dups.index)]

    peak_tissue_matrix = pd.DataFrame(
        np.zeros([peaks_merged.shape[0], len(peaks_merged.columns.values[7:])], dtype='int64'),
        columns=peaks_merged.columns.values[7:])
    peak_tissue_matrix.set_index(peaks_merged.index.values, inplace=True)
    for col in peak_tissue_matrix.columns.values:
        print(os.path.join(peak_folder, col))
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
        for i in orig_id[map(lambda x: ',' in x, orig_id)]:
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
    """ Function that adds the column Annotation from anno_file to the merge_file"""
    merge_df = pd.read_csv(merge_file, index_col=0, sep='\t')
    anno_df = pd.read_csv(anno_file, sep='\t', na_filter=False, index_col=0)
    anno_df.index.name = "ID"
    anno_df = anno_df['Annotation']
    merge_df = pd.concat((merge_df, anno_df), axis=1)
    merge_df.to_csv(merge_file+'.anno',sep='\t')
    return merge_df


###
# Get distances to regions of interest for the peaks
###
def distance_to_landmarks(anno_peaks, landmark_df, main_landmark='transcript_id',
                          landmark_cols=('gene','gene_id')):
    """ Funcntion that adds distance to the nearest landmark from the annotation file, regardless if its on the same
        strand or not. It will also note if it is on the same strand by adding that column"""
    for i in landmark_cols:
        anno_peaks['Nearest ' + i] = ''

    for ind, val in (anno_peaks.iterrows()):
        filt = landmark_df[val['Chr'] == landmark_df['Chr']]
        peak_start = (val['Start'] + val['End']) / 2
        if len(filt) != 0:
            filt2 = np.abs(filt['actual_start'] - peak_start).argmin()
            anno_peaks.set_value(ind, 'Nearest TSS', landmark_df.loc[filt2, main_landmark])

            for land in landmark_cols:
                anno_peaks.set_value(ind, 'Nearest ' + land, landmark_df.loc[filt2, land])
            #anno_peaks.set_value(ind, 'Nearest gene', landmark_df.loc[filt2, 'gene'])
            #anno_peaks.set_value(ind, 'Nearest gene_id', landmark_df.loc[filt2, 'gene_id'])

            ## Get distance to nearest gene. If on - strand, tss 'End' is the beginning
            if landmark_df.loc[filt2, 'Strand'] == '+':
                anno_peaks.set_value(ind, 'Distance to TSS', peak_start - landmark_df.loc[filt2, 'actual_start'])
            elif landmark_df.loc[filt2, 'Strand'] == '-':
                anno_peaks.set_value(ind, 'Distance to TSS', landmark_df.loc[filt2, 'actual_start'] - peak_start)
            anno_peaks.set_value(ind, 'isSameStrand', val['Strand'] == landmark_df.loc[filt2, 'Strand'])
    return anno_peaks


def wrap_distance_to_landmarks(peaks_file, landmark_file='/data/isshamie/genome/start_site_mRNA_updated_final_sort.tsv',
                          main_landmark = 'transcript_id',
                               landmark_cols=('gene','gene_id'),
                               is_bed=False,
                               output_f=None,is_parallel=False, num_par=4):
    """ Wrapper for distance_to_landmarks"""
    landmark_df = pd.read_csv(landmark_file, sep="\t", index_col=0)
    if is_bed:
        anno_peaks = Homer.read_bed_file(peaks_file)
    else:
        anno_peaks = pd.read_csv(peaks_file,sep="\t",index_col=0)

    # Only keep ID, Chr, Start, End, Strand, Stat, Annotation
    to_keep = ['Chr', 'Start', 'End', 'Strand', 'Stat']
    if 'Annotation' in anno_peaks.columns.values:
        to_keep.append('Annotation')
    anno_peaks = anno_peaks[to_keep]
    anno_peaks['Nearest TSS'] = ''
    anno_peaks['Distance to TSS'] = np.infty
    anno_peaks['isSameStrand'] = False

    if is_parallel:
        anno_peaks = parallel_functions.parallel_df(anno_peaks,distance_to_landmarks,
                                                    func_args = (landmark_df, main_landmark, landmark_cols),
                                                    num_processes=num_par)
    else:
        anno_peaks = distance_to_landmarks(anno_peaks, landmark_df, main_landmark,landmark_cols)
    if not output_f is None:
        anno_peaks.to_csv(output_f,sep='\t')
    return anno_peaks


###############################################
# Helper functions for creating annotation
def get_peaks_within_distance(peaks, distance=(-1000, +100)):
    filt = peaks[(peaks['Distance to TSS'] > distance[0]) & (peaks['Distance to TSS'] < distance[1])]
    return filt


def get_genes_peaks(peaks, gene, distance=(-1000, +100)):
    filt = peaks[(peaks['Nearest gene'] == gene) &
                 ((peaks['Distance to TSS'] > distance[0]) & (peaks['Distance to TSS'] < distance[1]))]
    return filt


def get_gene_ids_peaks(peaks, gene, distance=(-1000, +100)):
    filt = peaks[(peaks['Nearest gene_id'] == gene) &
                 ((peaks['Distance to TSS'] > distance[0]) & (peaks['Distance to TSS'] < distance[1]))]
    return filt


def get_nearest_tss(peaks_gene):
    # peaks of one gene
    peaks_gene['abs'] = np.abs(peaks_gene['Distance to TSS'])
    idx = peaks_gene['abs'].idxmin()
    return idx, peaks_gene.loc[idx, 'Distance to TSS']


def tss_per_gene(peaks_gene):
    # peaks of one gene
    return len(set(peaks_gene['Nearest TSS']))


def tissues_gene_expressed(peaks_gene, peaks_tissue_matrix):
    # peaks of one gene
    curr_gene_vals = peaks_tissue_matrix[peaks_tissue_matrix.isin(peaks_gene.idnex)]
    tiss = (curr_gene_vals > 1).any(axis=1)
    tiss = tiss.columns[tiss]
    return tiss


###############################################


###
# Create an updated annotation
###
def create_anno_centric_df(peaks_df, tss_df, peak_tissue_matrix, peak_bin=(-1000, 100), anno_col='Nearest TSS',
                           tss_df_col='transcript_id',
                           f_save=None, allow_intron = False,
                           has_gene_id=False,
                           allow_downstream_exons=True, allow_cds=True):
    gene_list = tss_df[tss_df_col].unique()

    prom_col1 = 'bin_%d_%d' % (peak_bin[0], peak_bin[1])  # This column is for + and - strands within the bin range
    prom_col2 = 'sameStrand_bin_%d_%d' % (peak_bin[0], peak_bin[1])  # This column is only for the same strand

    df = pd.DataFrame(index=gene_list, columns=['peaks', 'closest_to_X_numberOfPeaks', prom_col1, prom_col2,
                                                'samples', 'minDistance', 'minDistancePeak',
                                                'maxValue', 'maxPeakId', 'maxSample',
                                                'Chr', 'Start', 'End', 'Strand',
                                                 'gene', 'transcript'])

    if has_gene_id:
        df["gene_id"] = np.nan
    df['closest_to_X_numberOfPeaks'] = 0
    df['peaks'] = [[]] * len(df)
    df['minDistance'] = np.infty
    df['samples'] = [[]] * len(df)
    df['maxSamples'] = [[]] * len(df)
    df[prom_col1] = [[]] * len(df)
    df[prom_col2] = [[]] * len(df)
    for i in tqdm.tqdm_notebook(peaks_df.groupby([anno_col])):
        curr_gene = i[0]
        if curr_gene == '':
            continue
        df.loc[curr_gene, 'closest_to_X_numberOfPeaks'] += len(i[1])
        curr_peaks = i[1]
        df.at[curr_gene, 'peaks'] = list(curr_peaks.index)

        ## Filter to only within the distance of peak_bin
        curr_peaks = get_peaks_within_distance(curr_peaks, distance=peak_bin)

        #These contain peaks nearby, but not necessarily TSS (opp strand, intron,..)
        df.at[curr_gene, prom_col1] = list(curr_peaks.index)


        ## Filter to only peaks with isSameStrand
        curr_peaks = curr_peaks[curr_peaks['isSameStrand']]
        
        
        ## Filter ones with annotation of intron if flagged
        if not allow_intron:
            curr_peaks = curr_peaks[~curr_peaks['Annotation'].str.contains('intron')]
            curr_peaks = curr_peaks[
                ~curr_peaks['Annotation'].str.contains('TTS')]
        ###NOT DONE YET
        if not allow_downstream_exons:
            print("NOT implemented yet. Please set to True")
            curr_peaks = curr_peaks[curr_peaks['Annotation'].str.contains('exon [2,9]')]
            return
        if not allow_cds:
            print("NOT implemented yet. Please set to True")
            #curr_peaks = curr_peaks[curr]
        ###TILL HERE
        
        ## This is where the actual TSS peaks are placed
        df.at[curr_gene, prom_col2] = list(curr_peaks.index)
        
        
        ## Check if we still have any peaks left
        if len(curr_peaks) == 0:
            continue 
        if isinstance(curr_peaks, pd.Series):
            print('is a series', curr_peaks)

        ## IF we are here, then we have a TSS! Now to handle potentially having multiple peaks

        ## Find the closest peak
        minPeak, minPeakValue = get_nearest_tss(curr_peaks)
        df.at[curr_gene, 'minDistance'] = minPeakValue
        df.at[curr_gene, 'minDistancePeak'] = minPeak

        ## Store which tissues that have a peak
        ## Try statement b/c few peaks were not in the expression matrix for some reason
        try:
            # Get the samples that had the peak
            curr_tissues = peak_tissue_matrix.columns[(peak_tissue_matrix.loc[curr_peaks.index] > 0).any()].values
            curr_tissues = list(map(lambda x: os.path.basename(x),
                                    curr_tissues))
            df.at[curr_gene, 'samples'] = curr_tissues
            curr_tss_df_peaks = peak_tissue_matrix[peak_tissue_matrix.index.isin(curr_peaks.index)]
            # # If
            # if 'CHO' in curr_tissues:
            #     df.at[curr_gene, 'is_in_cho_promoter'] = True

            ## Get the peak that had the highest Stat (Here Counts per Million)
            max_peak = curr_tss_df_peaks.max(axis=1).idxmax()
            max_value = curr_tss_df_peaks.max(axis=1).max()
            df.at[curr_gene, 'maxValue'] = max_value
            df.at[curr_gene, 'maxPeakId'] = max_peak
            ## Samples that had the maxPeak
            max_tissues = peak_tissue_matrix.columns[(peak_tissue_matrix.loc[max_peak] > 0)].values
            max_tissues = list(map(lambda x: os.path.basename(x),
                                   max_tissues))
            df.at[curr_gene, 'maxSamples'] = max_tissues

            ## Get peak info for max peak
            start = curr_peaks.loc[max_peak, 'Start']
            end = curr_peaks.loc[max_peak, 'End']
            chrom = curr_peaks.loc[max_peak, 'Chr']
            strand = curr_peaks.loc[max_peak, 'Strand']
            gene = curr_peaks.loc[max_peak, 'Nearest gene']
            if has_gene_id:
                gene_id = curr_peaks.loc[max_peak, 'Nearest gene_id']
            mrna = curr_peaks.loc[max_peak, 'Nearest TSS']

            ## Assigning values
            df.at[curr_gene, 'transcript'] = mrna
            df.at[curr_gene, 'gene'] = gene
            if has_gene_id:
                df.at[curr_gene, 'gene_id'] = gene_id
            df.at[curr_gene, 'Start'] = start
            df.at[curr_gene, 'End'] = end
            df.at[curr_gene, 'Chr'] = chrom
            df.at[curr_gene, 'Strand'] = strand

        except KeyError:
            print(curr_peaks.index + ' not in peak expression file. Maybe it was dropped from being a duplicate? Or no peak had a value in it.')

    df['hasGene'] = df.fillna('')[prom_col2].apply(lambda x: len(x) > 0)
    df['Number of SS'] = df[prom_col2].fillna('').apply(lambda x: len(x))

    print('Total Number: ', len(df))
    print('Total Number with start sites: ', np.sum(df['hasGene']))
    print('Fraction with start sites: ', 1.0 * np.sum(df['hasGene']) / len(df))
    return df


def wrap_create_anno_centric(peaks_file, peaks_expression_file, tss_file,
                             peak_bin=(-1000, 100), anno_col='Nearest TSS',
                             tss_df_col='transcript_id',
                             f_save=None, allow_intron=False):

    tss_df = pd.read_csv(tss_file, sep='\t', index_col=0)
    peaks_df = pd.read_csv(peaks_file,sep='\t',index_col=0,na_filter=False)
    peak_tissue_matrix = pd.read_csv(peaks_expression_file,sep='\t',index_col=0)
    df = create_anno_centric_df(peaks_df, tss_df, peak_tissue_matrix,
                           peak_bin=peak_bin, anno_col=anno_col,
                           tss_df_col=tss_df_col,
                           f_save=f_save, allow_intron=allow_intron)

    if f_save is not None:
        df.to_csv(f_save + '.tsv',sep='\t')
        df.to_pickle(f_save + '.p')
        meta_f = f_save + '.params'
        meta = {'peak_bin':peak_bin,'tss_df_col':tss_df_col,'allow_intron':allow_intron,
                'peaks_file':peaks_file, 'peaks_expression_file':peaks_expression_file, 'tss_file':tss_file}
        pickle.dump(meta,open(meta_f,'wb'))
    return df


###
# Create a landmark(gene)-by-expression matrix.
###
def df_to_TSS_expression(f_in, f_mat_in, col_name=None, f_out=None):
    """Create geneXpeakValue matrix where rows are genes and cols are tissues.
    For each gene, take all the peaks that match that, and for each tissue, take its max expression for each tissue
    """
    df = pickle.load(open(f_in,'rb'))
    peak_tissue_matrix = pd.read_csv(f_mat_in,sep='\t',index_col=0)

    if col_name is None:
        col_name = df.columns.values[df.columns.str.contains('sameStrand_bin_')][0] # Default tss name
    gene_tissue_matrix = pd.DataFrame(index=df.index, columns=peak_tissue_matrix.columns, dtype='float32').fillna(0)

    for ind,val in tqdm.tqdm_notebook(df.iterrows()):
        if type(val[col_name]) == list:
            curr_peaks = peak_tissue_matrix[peak_tissue_matrix.index.isin(val[col_name])]
            gene_tissue_matrix.loc[ind] = curr_peaks.max()
    gene_tissue_matrix.fillna(0,inplace=True)
    if f_out is not None:
        gene_tissue_matrix.to_csv(f_out,sep='\t')
    return gene_tissue_matrix


###
# Takes the resulting annotation df and converts into a bed file. The name includes the samples that had the maxValue
###
def convert_anno_to_bed(txn_df,tss_annotation,f_save):
    tss = pd.read_csv(tss_annotation,sep='\t',index_col=0)
    txn_hasTSS = txn_df[txn_df['hasGene']]
    txn_hasTSS = txn_hasTSS[~(txn_hasTSS.isnull().any(axis=1))]

    # Putting into bed format
    tss_peakCenter_bed = pd.DataFrame()
    for ind, val in txn_hasTSS.iterrows():
        cs = tss[tss['transcript_id'] == ind]['cs'].values[0]

        # Drop the experiment type and just keep the tissue
        tissue_names = ','.join(set(list(map(lambda x: x.split('_')[0], val['tissues']))))

        title = '%s;%s;%s;%s;cs=%d' % (val['gene'], val['gene_id'], ind, tissue_names, cs)
        strand = tss[tss['transcript_id'] == ind]['Strand'].values[0]

        tss_peakCenter_bed = tss_peakCenter_bed.append(
            {
                'chrom': val['Chr'],
                'chromStart': val['Start'],
                'chromEnd': val['End'],
                'name': title,
                'score': val['maxValue'],
                'strand': strand
            },
            ignore_index=True)

    tss_peakCenter_bed = tss_peakCenter_bed[[
        'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'
    ]]
    tss_peakCenter_bed['chromStart'] = tss_peakCenter_bed['chromStart'].astype(int)
    tss_peakCenter_bed['chromEnd'] = tss_peakCenter_bed['chromEnd'].astype(int)
    tss_peakCenter_bed.to_csv(f_save + '.bed', sep='\t', header=None, index=None)
    return


###
# Gets all the peaks that had a TSS in it
###
def retrieve_all_peaks_from_anno(anno_f, merged_f, output_f,
                                 col_name=None, is_max=False):
    """Will filter the peaks to the ones annotated as TSSs.

        Parameters
        ----------
        anno_f : The tss file
        merged_f: The original merged peaks file which links the new
        peak ids to the original merged file peak ids.
        output_f: The file to save the output as
        col_name: The column of the list of peaks for each
        gene/transcript
        is_max: If True, will only take the peaks with the highest
        avg count (this is in the column maxPeakId)
        """
    anno_df = pickle.load(open(anno_f, 'rb'))
    peaks_df = pd.read_csv(merged_f, sep='\t', index_col=0)
    if col_name is None:
        col_name = anno_df.columns.values[
            anno_df.columns.str.contains('sameStrand_bin_')][0]
    if is_max:
        all_peaks = anno_df["maxPeakId"].values
    else:
        all_peaks = [item for sublist in list(anno_df[col_name].values) for
                     item in sublist]
    peaks_df_filt = peaks_df[peaks_df.index.isin(all_peaks)]
    peaks_df_filt.to_csv(output_f, sep='\t')
    return


def retrieve_sample_peaks_from_anno(anno_f, merged_f, output_f,
                                    sample_name, sample_f=None,
                                    col_name=None,
                                    use_sample_peaks=False):
    """Will filter peaks annotated as TSSs and found in the
    sample. Can retrieve the merged peaks or the original samples peaks.

        Parameters
        ----------
        anno_f : The tss file
        merged_f: The original merged peaks file which links the new
        peak ids to the original merged file peak ids.
        output_f: The file to save the output as
        sample_name: the name of the sample in the merged file
        sample_f: the sample peak file to extract the information
            from. Optional.
        use_sample_peaks: To use the sample peaks or merged. Optional
        col_name: The column of the list of peaks for each
            gene/transcript
        """

    anno_df = pickle.load(open(anno_f, 'rb'))
    peaks_df = pd.read_csv(merged_f, sep='\t', index_col=0)
    if col_name is None:
        col_name = anno_df.columns.values[
            anno_df.columns.str.contains('sameStrand_bin_')][0]
    all_peaks = [item for sublist in list(anno_df[col_name].values) for
                 item in sublist]
    peaks_df_filt = peaks_df[peaks_df.index.isin(all_peaks)]

    # if use_sample_peaks: peaks from original peak file or
    #  else take from the merged file
    peaks_to_keep = []
    for ind, val in peaks_df_filt.iterrows():
        if sample_name in val["Parent files"].split("|"):
            if use_sample_peaks:
                peaks_to_keep.append(peaks_df_filt.loc[ind, sample_name])
            else:
                peaks_to_keep.append(ind)
    if use_sample_peaks:
        sample_df = Homer.read_peak_file(
            sample_f)  # original peak file
        sample_df_filt = sample_df[sample_df.index.isin(peaks_to_keep)]
    else:
        sample_df_filt = peaks_df_filt.loc[peaks_to_keep]

    sample_df_filt.to_csv(output_f, sep='\t')
    return sample_df_filt


def add_max_info(merged_df, expr_peaks):
    """ Add the maximum peak information"""
    merged_df["max_file"] = ""
    merged_df["max_file_peak_name"] = ""
    merged_df["max_Chr"] = ""
    merged_df["max_Start"] = -1
    merged_df["max_End"] = -1
    merged_df["max_Stat"] = 0.0
    merged_df["cho_file"] = ""
    merged_df["cho_file_peak_name"] = ""
    merged_df["cho_Chr"] = ""
    merged_df["cho_Start"] = -1
    merged_df["cho_End"] = -1
    merged_df["cho_Stat"] = 0.0

    cho_cols = merged_df.columns[merged_df.columns.str.contains("CHO")]

    for ind, val in tqdm.tqdm(merged_df.iterrows()):
        max_f = expr_peaks.loc[ind].idxmax()
        merged_df.at[ind, "max_file"] = max_f
        merged_df.at[ind, "max_file_peak_name"] = merged_df.at[
            ind, max_f]

        # CHO related
        max_cho = expr_peaks.loc[ind, cho_cols].idxmax()
        if expr_peaks.loc[ind, max_cho] > 0:
            merged_df.at[ind, "cho_file"] = max_cho
            merged_df.at[ind, "cho_file_peak_name"] = merged_df.at[
                ind, max_cho]

    for sample_f, df in merged_df.groupby("max_file"):
        # pd.read_csv(f,sep="\t")
        print(sample_f)
        sample_df = Homer.read_peak_file(sample_f)

        # Loop through the merged_df in which the max peak was found in sample_f
        for ind, val in df.iterrows():
            # For each peak, find the Chr, Start, End from the sample_df.
            if "," in val[
                "max_file_peak_name"]:  # Multiple peaks present when merged
                get_p = \
                sample_df.loc[val["max_file_peak_name"].split(",")][
                    "Normalized Tag Count"].idxmax()
                get_p = sample_df.loc[get_p]
            else:
                get_p = sample_df.loc[val["max_file_peak_name"]]

            merged_df.at[ind, "max_Chr"] = get_p["Chr"]
            merged_df.at[ind, "max_Start"] = get_p["Start"]
            merged_df.at[ind, "max_End"] = get_p["End"]
            merged_df.at[ind, "max_Stat"] = get_p[
                "Normalized Tag Count"]

    # CHO
    for sample_f, df in merged_df.groupby("cho_file"):
        print(sample_f)
        if sample_f == "":
            continue
        sample_df = Homer.read_peak_file(sample_f)

        # Loop through the merged_df in which the max peak was found in sample_f
        for ind, val in df.iterrows():
            # For each peak, find the Chr, Start, End from the sample_df.
            if "," in val[
                "cho_file_peak_name"]:  # Multiple peaks present when merged
                get_p = \
                sample_df.loc[val["cho_file_peak_name"].split(",")][
                    "Normalized Tag Count"].idxmax()
                get_p = sample_df.loc[get_p]
            else:
                get_p = sample_df.loc[val["cho_file_peak_name"]]

            merged_df.at[ind, "cho_Chr"] = get_p["Chr"]
            merged_df.at[ind, "cho_Start"] = get_p["Start"]
            merged_df.at[ind, "cho_End"] = get_p["End"]
            merged_df.at[ind, "cho_Stat"] = get_p[
                "Normalized Tag Count"]

    return merged_df


def wrap_add_max_info(merged_f, expr_peaks_f,out_f=None):
    #merged_f = "Results/merged/samples.merge"
    merged_df = pd.read_csv(merged_f, sep="\t", index_col=0)
    #peak_f = "Results/merged/samples.merge.peaksexpression.log10"
    expr_peaks = pd.read_csv(expr_peaks_f, index_col=0, sep="\t")
    merged_df = add_max_info(merged_df, expr_peaks)
    if out_f is not None:
        merged_df.to_csv(out_f,sep="\t")
    return merged_df


##
# Peak filtering
##
def filter_only_with_divergent(df, divergent_threshold=500):
    """ Filter if there are divergent transcripts"""
    inds_to_keep = []
    for ind, val in df[df['isSameStrand']].iterrows():
        curr = df[
            (val['Chr'] == df['Chr']) & ~(df['isSameStrand'])]
        if val['Strand'] == '-':
            curr = curr[
                (val['End'] - curr['Start'] < divergent_threshold) & (
                            val['End'] - curr['Start'] > 0)]
        else:
            curr = curr[
                (val['Start'] - curr['End'] < divergent_threshold) & (
                            val['Start'] - curr['End'] > 0)]

        # If there are divergent peaks, add this peak
        if len(curr) > 0:
            inds_to_keep.append(ind)

    df = df.loc[inds_to_keep]
    return df


##
# Filtering peaks where the peak index is the name and Distance to
# TSS is a column.
##
def filt_peaks(peaks, distance=(-1000, +100), thresh=0, thresh_col=None,
               divergent_thr=0, save_f=None):
    """ If only one distance measure is given, assumes that the absolute
        value of the 'Distance to TSS' column in peaks has to be
        greater than the distance.
    """

    # The first two statements are for looking for a distance greater
    # than the TSS
    if type(distance) is int:
        peaks = peaks[np.abs(peaks['Distance to TSS']) > distance]
    elif len(distance) == 1:
        peaks = peaks[np.abs(peaks['Distance to TSS']) > distance[0]]

    # Looks for between the first and second elements gievn in distance
    else:
        peaks = peaks[(peaks['Distance to TSS'] > distance[0]) & (
                peaks['Distance to TSS'] < distance[1])]

    # If there's a score threshold
    if thresh_col is not None:
        peaks = peaks[peaks[thresh_col] > thresh]

    # If there's a divergent threshold
    if divergent_thr != 0:
        peaks = filter_only_with_divergent(peaks,
                                           divergent_threshold=divergent_thr)

    # Save the filtered peaks and the params used
    if save_f is not None:
        peaks.to_csv(save_f, sep="\t")
        pickle.dump([distance,thresh,thresh_col,divergent_thr],
                    open(save_f + ".params.p","wb"))
    return peaks
