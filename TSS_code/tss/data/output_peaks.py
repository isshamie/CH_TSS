from tss.utils import create_output

def run(p):
    #create_output()

    tss_annotation = p['reference']['TSS_F']
    mRNA_peak_file =  p['reference']["mRNA_peak_file"]

    txn_f = p["gene centric"]["filenames"]["transcript centric peaks matrix"]#"Results/tss_annotation/txn_df_02_tissues.p"
    # expr_f = "Results/merged/samples.merge.peaksexpression.log10"

    expr_f = p["gene centric"]["filenames"]["peak tissues expression log10"]
    #"Results/merged/tissues.merge.peaksexpression.log10"
    anno_f = "Results/tss_annotation_peaks/all_peaks_txn_df_02_maxinfo.tsv"
    rank_func = create_output.rank_median_samples
    atac_f = "Results/ATAC_results/CHO/ATAC_ppr.naive_overlap.narrowPeak.sort"
    tss_f = mRNA_peak_file
    all_atac = "Results/ATAC_results/ATAC_merge.bed"


    out_f = p["output_peaks"]["filenames"]["final peaks"]

    bed_df, meta_df = create_output.construct_all_peaks(txn_f,
                                                        expr_f, anno_f,
                                          rank_func, tss_f, atac_f,
                                          all_atac=all_atac,
                                          out_f=out_f)

    return

