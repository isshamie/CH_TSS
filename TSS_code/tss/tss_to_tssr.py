import pandas as pd
import tqdm


def enhancer_tss_to_tssr(df, dist_thresh=1000, dist_to_gene=1000,
                         value='mean'):
    """ Collapses nearby peaks and returns new region IDs and info about their peaks

    :param df: peak df with Chr, Start, End, Strand keys.
    :param dist_thresh: how wide to merge
    :param dist_to_gene: If set, 'Distance to TSS' needs to be set to use as a filter
    :param value: NOT IMPLEMENTED. The way to integrate expression values. e.g. 'mean'
    :return tssr_id: {new id: [old ids associated with new]}
    :return tssr_info: pd.DataFrame where index is new_id, and columns are information about number of peaks associated.
    """

    # Calculate the distance of adjacent peaks.
    df = df.copy().sort_values(["Chr", "Start"])
    dist_to_df = pd.DataFrame(
        df.iloc[1:]["End"].values - df.iloc[:-1]["End"].values,
        index=df.iloc[1:].index, columns=["Distance"])

    # initial TSS
    # current TSS attributes
    curr_chr = df.iloc[0]["Chr"]
    curr_len = df.iloc[0]["End"] - df.iloc[0]["Start"]
    curr_strand = df.iloc[0]["Strand"]

    # Setup enhancer dicts
    curr_enh_id = 0
    tssr_id = {curr_enh_id: [
        0]}  # [curr_enh_id] = df.index.values[curr_enh_id]
    tssr_info = {
        curr_enh_id: dict(num_tss=1, length=curr_len, chrom=curr_chr,
                          num_pos_strand=int(curr_strand == "+"),
                          num_neg_strand=int(curr_strand != "+"))}

    # Loop through peaks, adding to the current new peaks id until
    # distance threshold is passed or different strand.
    for ind, val in tqdm.tqdm_notebook(dist_to_df.iterrows()):
        curr_strand = df.loc[ind, 'Strand']
        if tssr_info[curr_enh_id]["length"] + val[
            'Distance'] > dist_thresh or curr_chr != df.loc[ind, 'Chr']:
            # Reset indices and add to dict
            curr_enh_id += 1
            tssr_id[curr_enh_id] = [ind]
            curr_chr = df.loc[ind, 'Chr']
            curr_len = df.loc[ind, 'End'] - df.loc[ind, 'Start']
            tssr_info[curr_enh_id] = dict(num_tss=1, length=curr_len,
                                          chrom=curr_chr,
                                          num_pos_strand=int(
                                              curr_strand == "+"),
                                          num_neg_strand=int(
                                              curr_strand == "-"))
        else:
            # Update ID
            tssr_id[curr_enh_id].append(ind)
            # Update metrics
            tssr_info[curr_enh_id]["num_tss"] += 1
            tssr_info[curr_enh_id]["num_pos_strand"] += int(
                curr_strand == "+")
            tssr_info[curr_enh_id]["num_neg_strand"] += int(
                curr_strand == "-")
            tssr_info[curr_enh_id]["length"] += val['Distance']

    tssr_info = pd.DataFrame(tssr_info).transpose()
    return tssr_id, tssr_info


def filter_distal_peaks(tssr_id, distal_peaks):
    distal_peaks_enh = {}
    for t in distal_peaks:
        print('t', t)
        distal_peaks_enh[t] = distal_peaks[t].copy()
        for curr_enh_id in tssr_id:
            if len(tssr_id[curr_enh_id]) == 1:
                continue
            else: # just take one of the peaks
                distal_peaks_enh[t] = distal_peaks_enh[t] - set(
                    tssr_id[curr_enh_id][1:])
    return distal_peaks_enh


def get_tissue_peaks(distal_peaks, tssr_id):
    """
    Get the first peak index for each enhancer_id region
    :param distal_peaks: dict where each key is a sample, and the value is a set of peak ids.
    :param tssr_id: {enh_id:[peak_ids]}. Only the first peak_id will be kept in the distal_peaks for each enh_id
    :return: distal_peaks_enh #updated to regions, not peaks.
    """
    distal_peaks_enh = {}
    for t in distal_peaks:
        print('t', t)
        distal_peaks_enh[t] = distal_peaks[t].copy()
        for curr_enh_id in tssr_id:
            if len(tssr_id[curr_enh_id]) == 1:
                continue
            distal_peaks_enh[t] = distal_peaks_enh[t] - set(
                tssr_id[curr_enh_id][
                1:])  # distal_peaks_enh[t] = distal_peaks[t]-set(tssr_id)
    return distal_peaks_enh


