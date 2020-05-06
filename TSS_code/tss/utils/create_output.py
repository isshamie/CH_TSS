""" The goal of these functions is to go from the gene-centric TSS
dataframes to a final reported output for people to use. The output
is going to be modeled after the FANTOM project, where there can be
multiple peaks per isoform start site.
Each TSS peak is assigned an ID:
 ID: pN@GENE_transcriptID; where GENE indicates gene symbol or
 transcript name
 and N indicates the rank in the ranked list of promoter activities for
 that gene e.g. p1@SPI1. If the peak is not experimental, use p0@GENE.

Chr:
Start:
End:

CS:
Binary score, where the most significant bit is the most important
predictor.

1. Filter if Max <1CPM
2. ATAC-peak nearby
3. Only 1 tissue and No sample replicates
4. TATAA Box detected
5. Inr motif
e.g. 11001:
1: Max greater than 1 CPM, 1: ATAC-seq peak seen, 0: only 1 sample
0: TATAA not seen 1: Inr motif seen


Tissues: Comma separated string

CDS Upstream
The ranking & confidence score could be varied, along with additional
filters.

May also put:
Gene Upstream Position
Coding Region Start	Gene
Downstream Position


Also, create ATAC-seq bed file with regions of open-chromatin region
in CHO with associated peak ID.
"""

from tss.utils.Homer import *
import tqdm
#import parallel_functions as pf


def construct_nonExperimental_peakID(txn_f, tss_f):
    tss = pd.read_csv(tss_f, sep="\t", index_col=0)
    txn = pickle.load(open(txn_f, "rb"))
    tss = tss[tss.index.isin(txn[~txn["hasGene"]].index)]
    bed_df = pd.DataFrame(
        columns=["Chr", "Start", "End", "Strand", "Stat", "ID"])
    non_peaks_map = dict()
    meta_df = pd.DataFrame(
        columns=["Tissues", "cs", "CHO ATAC Region", "ID", "Gene",
                 "Gene ID", "Transcript"])
    for ind, val in tqdm.tqdm(tss.iterrows()):
        # name = "%s_%s_%s" % (
        #     val["gene"], val["gene_id"], ind)
        name = "%s_%s" % (val["gene"], ind)
        curr_peak = "p0@%s" % name  # p1@SPI1
        bed_df.loc[curr_peak] = [val["Chr"], val["Start"], val["End"],
                                 val["Strand"], 0, curr_peak]

        meta_df.at[curr_peak, "Gene"] = val["gene"]
        #meta_df.at[curr_peak, "Gene ID"] = val["geneID"]
        meta_df.at[curr_peak, "Transcript"] = ind
        meta_df.at[curr_peak, "Is Experimental"] = 0
        non_peaks_map[name] = curr_peak
    return bed_df,meta_df


def construct_peakID(txn_anno_f, expr_f, anno_f, rank_func,
                     col_name=None):
    txn_anno = pickle.load(open(txn_anno_f, "rb"))
    expr_peaks = pd.read_csv(expr_f, sep="\t", index_col=0)
    anno_df = pd.read_csv(anno_f, sep="\t", index_col=0)
    #print(anno_df.head())
    #print(len(anno_df))
    anno_df = anno_df[
        ["Chr", "Start", "End", "Strand", "Stat", "cho_Chr",
         "cho_Start", "cho_End", "cho_Stat", "max_Chr", "max_Start",
         "max_End", "max_Stat"]]
    peaks_map = dict()
    cs = dict()
    bed_df = pd.DataFrame(
        columns=["Chr", "Start", "End", "Strand", "Stat", "ID"])
    meta_df = pd.DataFrame(
        columns=["Tissues", "cs", "CHO ATAC Region", "ID", "Gene",
                 "Gene ID", "Transcript"])

    if col_name is None:
        col_name = txn_anno.columns[
            txn_anno.columns.str.contains("sameStrand")][0]

    # Loop through transcripts with peaks
    for ind, val in (txn_anno[txn_anno[
        "hasGene"]].iterrows()):
        curr_peaks = val[col_name]
        # Create the IDs
        ordered_list = rank_func(curr_peaks, expr_peaks)
        for ind2, p in enumerate(ordered_list):
            # curr_peak = "p%d@%s_%s_%s" % (
            #     ind2 + 1, val["gene"], val["gene_id"],
            #     val["transcript"])  # p1@SPI1
            curr_peak = "p%d@%s_%s" % (
                ind2 + 1, val["gene"], val["transcript"])  # p1@SPI1
            peaks_map[p] = curr_peak

            # Get their genome coordinates
            p_chr, p_start, p_end, p_stat = get_peak_info(p, anno_df,
                                                          is_cho=True,
                                                          is_orig=True)
            bed_df.loc[curr_peak] = [p_chr, p_start, p_end,
                                     val["Strand"], p_stat, curr_peak]

            _, cs_bin = cs_score(p, expr_peaks)
            cs[p] = cs_bin

            meta_df.at[curr_peak, "Tissues"] = val["Tissues"]
            meta_df.at[curr_peak, "cs"] = cs_bin
            meta_df.at[curr_peak, "ID"] = p
            meta_df.at[curr_peak, "Gene"] = val["gene"]
            #meta_df.at[curr_peak, "Gene ID"] = val["gene_id"]
            meta_df.at[curr_peak, "Transcript"] = val["transcript"]
            meta_df.at[curr_peak, "Is Experimental"] = 1


            # meta_df.at["CHO ATAC Region"] = atac_region(atac_f,
            #                                         bed_df.loc[curr_peak])

    ## This will make the file 0-based index inclusive
    bed_df["Start"] = bed_df["Start"] - 1
    bed_df["End"] = bed_df["End"]

    return bed_df, meta_df


def construct_all_peaks(txn_f, expr_f, anno_f, rank_func, tss_f, atac_f,
                        all_atac=None, out_f=None):
    """ Wrapper for constructing the experimental TSS and the genes
    that didnt have experimental TSSs. Assumes the input txn_f and tss_f
    are 1-based index, and this will convert to bed file, which is 0-based index.
    This does this by doing Start-1 at the end"""
    
    print("Constructing TSS not observed experimentally")
    peak2_df, meta2_df = construct_nonExperimental_peakID(txn_f, tss_f)
    bed2_df = peak2_df.copy()
    bed2_df["Start"] -= 1
    if out_f is not None:
        bed2_df[["Chr", "Start", "End", "ID", "Stat",
                "Strand"]].to_csv(out_f + ".ref.bed", sep="\t",
                                  header=None,
                                  index=False)
        peak2_df[["Chr", "Start", "End", "ID", "Stat",
                "Strand"]].to_csv(out_f + ".ref.tsv", sep="\t",
                                  header=None,
                                  index=False)
        meta2_df.to_csv(out_f + ".ref.meta.tsv", sep="\t")

    print("Constructing TSS observed experimentally")
    peak_df, meta_df = construct_peakID(txn_f, expr_f, anno_f, rank_func)
    bed_df = peak_df.copy()
    bed_df["Start"] -= 1
    if out_f is not None:
        bed_df[["Chr", "Start", "End",  "ID", "Stat",
                "Strand"]].to_csv(out_f + ".exp.bed", sep="\t",
                                  header=None,
                                  index=False)
        peak_df[["Chr", "Start", "End",  "ID", "Stat",
                "Strand"]].to_csv(out_f + ".exp.tsv", sep="\t",
                                  header=None,
                                  index=False)
        meta_df.to_csv(out_f + ".exp.meta.tsv", sep="\t")

    bed_df = pd.concat((bed_df, bed2_df), axis=0)
    bed_df = bed_df.sort_values(["Chr", "Start", "End", "Strand"])
    bed_df = bed_df[["Chr", "Start", "End",  "ID", "Stat", "Strand"]]

    meta_df = pd.concat((meta_df, meta2_df), axis=0)
    meta_df = meta_df.loc[bed_df.index]

    print("Including ATAC information")
    meta_df = atac_region(meta_df, bed_df, atac_bed=atac_f, all_atac=all_atac, thresh=250)
    meta_df["Is Experimental"] = meta_df["Is Experimental"].astype(int)

    bed_df.fillna("", inplace=True)
    meta_df.fillna("", inplace=True)
    if out_f is not None:
        bed_df[["Chr", "Start", "End", "ID", "Stat",
                "Strand"]].to_csv(out_f + ".bed", sep="\t", header=None,
                                  index=False)
        peak_df[["Chr", "Start", "End", "ID", "Stat",
                "Strand"]].to_csv(out_f + ".tsv", sep="\t", header=None,
                                  index=False)
        meta_df.to_csv(out_f + ".meta.tsv", sep="\t")
    return bed_df, meta_df


def cs_score(peak, peaks, extra=None, on_thresh=0.3):
    """
    Function that constructs confidence score for the peak.
    :param peak: The peak id of interest. This is is related to the
    peaks
    :param peaks: peaks expression dataframe
    :param on_thresh: Threshold for considering peak expressed.
    Default is log10(2)
    :param extra: This would be additional function that adds to the
    dictionary if not None
    :return:
    cs_dict: Dictionary of the mapping
    cs_score: The binnary variable that is the code
    """

    cs_dict = dict()
    cs_dict["Max over thresh "] = (peaks.loc[peak] > on_thresh).any()
    cs_dict["Multiple samples"] = (peaks.loc[
                                       peak] > on_thresh).sum() >= 2
    # cs_dict["has ATAC"] = peaks_meta.loc[peak]["has ATAC"]

    cs_bin = "0b%d%d" % (
        cs_dict["Max over thresh "], cs_dict["Multiple samples"])
    # + str(cs_dict["has ATAC"])

    if extra is not None:
        cs_dict, cs_bin = extra(cs_dict, cs_bin)
    return cs_dict, cs_bin


def rank_median_samples(peaks_of_interest, peaks, on_thresh=0):
    """A ranking function of peaks.
       Ranks peaks by taking the number of samples it was found in
       and multiplying it by the median value for all samples in
       which it passed the threshold"""
    vals = dict()
    for p in peaks_of_interest:
        curr_p = peaks.loc[p]
        num_samples = (curr_p > on_thresh).sum()
        vals[p] = num_samples * np.median(curr_p[curr_p > on_thresh])

    ordered_list = sorted(vals, key=vals.get, reverse=True)
    return ordered_list


def rank_sum_samples(peaks_of_interest, peaks, on_thresh=0):
    """A ranking function of peaks.
       Ranks peaks by taking the sum of peak values (log CPM)
       of samples that passed the threshold"""
    vals = dict()
    for p in peaks_of_interest:
        curr_p = peaks.loc[p]
        vals[p] = (curr_p[curr_p > on_thresh]).sum()

    ordered_list = sorted(vals, key=vals.get, reverse=True)
    return ordered_list


# def rank_harmonic_samples(peaks_of_interest, peaks, on_thresh=0):
#     """A ranking function of peaks.
#        Ranks peaks by taking the sum of peak values (log CPM)
#        of samples that passed the threshold"""
#     vals = dict()
#     for p in peaks_of_interest:
#         curr_p = peaks.loc[p]
#         vals[p] = (curr_p[curr_p > on_thresh]).sum()
#
#     ordered_list = sorted(vals, key=vals.get, reverse=True)
#     return ordered_list

def get_peak_info(p, peaks_df, is_cho=True, is_orig=True, to_log=True):
    # print("Before")
    # print(peaks_df.loc[p])
    # print(peaks_df.loc[p].dtype)
    if to_log:
        peaks_df.at[p, "Stat"] = np.log10(
            1.0 * peaks_df.at[p, "Stat"] + 1)
        peaks_df.at[p, "cho_Stat"] = np.log10(
            peaks_df.at[p, "cho_Stat"] + 1)
        peaks_df.at[p, "max_Stat"] = np.log10(
            peaks_df.at[p, "max_Stat"] + 1)
    # print("After")
    # print(peaks_df.loc[p])

    if not is_orig:
        return peaks_df.at[p, "Chr"], peaks_df.at[p, "Start"], \
               peaks_df.at[p, "End"], peaks_df.at[p, "Stat"]

    else:
        # See if CHO has a peak
        if is_cho:
            if peaks_df.at[p, "cho_Start"] != -1:  # Not -1
                return peaks_df.at[p, "cho_Chr"], peaks_df.at[
                    p, "cho_Start"], peaks_df.at[p, "cho_End"], \
                       peaks_df.at[p, "cho_Stat"]
        return peaks_df.at[p, "max_Chr"], peaks_df.at[p, "max_Start"], \
               peaks_df.at[p, "max_End"], peaks_df.at[p, "max_Stat"]


def atac_region(meta_df, bed_df, atac_bed=None, all_atac=None,
                colname="CHO ATAC Region", thresh=250):
    if atac_bed is not None:
        atac_df = pd.read_csv(atac_bed, sep="\t")
        cols = atac_df.columns.values.astype(str)
        cols[:6] = ["Chr", "Start", "End", "ID", "Stat", "Strand"]
        atac_df.columns = cols
        atac_df = atac_df.drop_duplicates(
            subset=["Chr", "Start", "End"])

        meta_df[colname] = ""
        for i, val in tqdm.tqdm(bed_df.iterrows()):
            ind = val["ID"]
            curr = atac_df[atac_df["Chr"] == val["Chr"]]
            curr = curr[np.abs(val["Start"] - curr["Start"]) < thresh]
            if len(curr) > 0:
                for ind2, val2 in curr.iterrows():
                    if meta_df.at[ind, colname] == "":
                        meta_df.at[ind, colname] = "%d:%d" % (
                            val2["Start"], val2["End"])
                    else:
                        meta_df.at[ind, colname] = "%s,%d:%d" % (
                            meta_df.at[ind, colname], val2["Start"],
                            val2["End"])

    if all_atac is not None:
        all_df = pd.read_csv(all_atac, sep="\t")
        cols = all_df.columns.values.astype(str)
        cols[:6] = ["Chr", "Start", "End", "ID", "Stat", "Strand"]
        all_df.columns = cols
        all_df = all_df.drop_duplicates(
            subset=["Chr", "Start", "End"])

        meta_df["has ATAC"] = 0
        for i, val in tqdm.tqdm(bed_df.iterrows()):
            ind = val["ID"]
            curr = all_df[all_df["Chr"] == val["Chr"]]
            curr = curr[np.abs(val["Start"] - curr["Start"]) < thresh]
            if len(curr) > 0:
                meta_df.at[ind, "has ATAC"] = 1
    return meta_df


def par_atac_region(all_df, atac_bed=None, all_atac=None,
                    colname="CHO ATAC Region", thresh=250):
    if atac_bed is not None:
        atac_df = pd.read_csv(atac_bed, sep="\t")
        cols = atac_df.columns.values.astype(str)
        cols[:6] = ["Chr", "Start", "End", "ID", "Stat", "Strand"]
        atac_df.columns = cols
        atac_df = atac_df.drop_duplicates(
            subset=["Chr", "Start", "End"])

        all_df[colname] = ""
        for i, val in tqdm.tqdm(all_df.iterrows()):
            ind = val["ID"]
            curr = atac_df[atac_df["Chr"] == val["Chr"]]
            curr = curr[np.abs(val["Start"] - curr["Start"]) < thresh]
            if len(curr) > 0:
                for ind2, val2 in curr.iterrows():
                    if all_df.at[ind, colname] == "":
                        all_df.at[ind, colname] = "%d:%d" % (
                            val2["Start"], val2["End"])
                    else:
                        all_df.at[ind, colname] = "%s,%d:%d" % (
                            all_df.at[ind, colname], val2["Start"],
                            val2["End"])

    if all_atac is not None:
        all_df = pd.read_csv(all_atac, sep="\t")
        cols = all_df.columns.values.astype(str)
        cols[:6] = ["Chr", "Start", "End", "ID", "Stat", "Strand"]
        all_df.columns = cols
        all_df = all_df.drop_duplicates(
            subset=["Chr", "Start", "End"])

        all_df["has ATAC"] = 0
        for i, val in tqdm.tqdm(all_df.iterrows()):
            ind = val["ID"]
            curr = all_df[all_df["Chr"] == val["Chr"]]
            curr = curr[np.abs(val["Start"] - curr["Start"]) < thresh]
            if len(curr) > 0:
                all_df.at[ind, "has ATAC"] = 1
    return all_df


#######################################################################
def exp_bed_to_refseq(bed_f, meta_f, refseq_f, save_f="",
                      is_unique=False, shift=0,is_bed=True):
    """Goes back to the actual locations of the transcripts in the
    output_bed file and creates a bed file of the original locations
    bed_f: The output bed file
    meta_f: The output meta file
    save_f: Where to save the file. If empty string, will return the
    file.
    refseq_f: Centered annotation peak
    unique_only: Boolean to remove duplicate locations
    sites
    file for where the original annotation TSSs are. e.g
    "/data/isshamie/genome/picr_final/mRNA_final
    .peak"
    is_bed: If true, will -1 for the tss, since that is 1-based index.
    """

    tss = pd.read_csv(refseq_f, sep="\t", index_col=0)
    tss = tss[~tss.index.duplicated()]
    meta_df = pd.read_csv(meta_f, sep="\t", index_col=0)
    bed_df = pd.read_csv(bed_f, sep="\t", header=None)

    out_meta = meta_df[meta_df["Transcript"].isin(tss.index)]
    out_bed = bed_df[bed_df[3].isin(out_meta.index)]

    out_bed = out_bed.copy()
    out_bed = out_bed.set_index(3)

    out_bed[0] = np.array(tss.loc[out_meta["Transcript"], "Chr"])
    out_bed[1] = np.array(tss.loc[out_meta["Transcript"]]["actual_start"]-shift)
    out_bed[2] = np.array(tss.loc[out_meta["Transcript"]]["actual_start"]+shift)

    if is_bed: #0-based index
        out_bed[1] -= 1

    out_bed[3] = out_bed.index
    out_bed = out_bed[[0, 1, 2, 3, 4, 5]]

    if is_unique:
        out_bed.drop_duplicates(subset=(0, 1, 2),keep='first',
                                inplace=True)
    if save_f == "":
        save_f = bed_f.strip(".bed") + "_refseq_centered.bed"
    out_bed.to_csv(save_f, sep="\t", header=None, index=False)

    return out_bed





#
# txn_f = "/mnt/jabba/data/isshamie/TSS/Analysis/10_15_Results/Results/tss_annotation/txn_df_02_tissues.p"
# #expr_f = "Results/merged/samples.merge.peaksexpression.log10"
# expr_f = "/mnt/jabba/data/isshamie/TSS/Analysis/10_15_Results/Results/merged/tissues.merge.peaksexpression.log10"
# anno_f = "/mnt/jabba/data/isshamie/TSS/Analysis/10_15_Results/Results/tss_annotation_peaks/all_peaks_txn_df_02_maxinfo.tsv"
# rank_func = rank_median_samples
# atac_f = "/mnt/jabba/data/isshamie/TSS/Analysis/10_15_Results/Results/ATAC_results/CHO/ATAC_ppr.naive_overlap.narrowPeak.sort"
# tss_f = "/mnt/jabba/data/isshamie/genome/picr_final/mRNA_final.peak" \
#         ""#tss_annotation
# all_atac = "/mnt/jabba/data/isshamie/TSS/Analysis/10_15_Results/Results/ATAC_results/ATAC_merge.bed"
# bed_df, meta_df = construct_all_peaks(txn_f, expr_f, anno_f, rank_func, tss_f, atac_f,
#                         all_atac=all_atac, out_f="/mnt/jabba/data/isshamie/TSS/Analysis/10_15_Results/Results/output/TSS1")
