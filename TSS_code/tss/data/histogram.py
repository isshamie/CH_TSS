from tss.utils import create_output
from tss.utils.histograms import hist_plot, heat_plot, \
    wrap_hist_plot, wrap_heat_plot
from tss.utils.Homer import hist


def run(p):
    p_stage = p["histogram"]
    p_stage_f = p_stage["filenames"]
    p_reference = p["reference"]
    p_global = p["global"]

    # Prior inputs
    raw_peaks = p["move_peaks"]["folder"]
    mRNA_peak_file = p_global["mRNA_peak_file"]

    ## Create RefSeq peak file for this
    peak_f = "Results/output/TSS1.exp.bed"
    meta_f = "Results/output/TSS1.exp.meta"
    save_f = "Results/Figures/Figure3/A.TSS1_mrna"
    tag_f = "Results/merged/tags_TSS_merged"

    # mrna_filt = "Results/Figures/Figure3/A.mrna"
    create_output.exp_bed_to_refseq(peak_f, meta_f,
                                    refseq_f=mRNA_peak_file,
                                    save_f=save_f, is_unique=True)

    # Experimental
    output_file_exp = "Results/Figures/Figure3/A.exp.hist"
    print('output histogram exp', output_file_exp)
    hist(tag_f, output_file_exp, ref_fa, anno_gff, mode='peak',
         peak=peak_f, region=4000, res=25, pc=3)
    hist_plot(output_file_exp)
    heat_df = heat_plot(output_file_exp + 'MatS',
                        save_f=output_file_exp + '.heat.png')

    # Refseq
    output_file_mrna = "Results/Figures/Figure3/A.TSS1_mrna.hist"
    input_file_mrna = save_f
    print('output histogram mrna', output_file_mrna)
    hist(tag_f, output_file_mrna, ref_fa, anno_gff, mode='peak',
         peak=save_f, region=4000, res=25, pc=3)
    hist_plot(output_file_mrna)
    heat_df = heat_plot(output_file_mrna + 'MatS',
                        save_f=output_file_mrna + '.heat.png')

    output_file_exp = "Results/Figures/Figure3/A.exp.hist"
    output_file_mrna = "Results/Figures/Figure3/A.TSS1_mrna.hist"

    wrap_hist_plot([output_file_exp, output_file_mrna],
                   hist_save="Results/Figures/Figure3/A.combine.png",
                   names=["Experimental", "RefSeq"])

    output_file_exp = "Results/Figures/Figure3/A.exp.hist"
    output_file_mrna = "Results/Figures/Figure3/A.TSS1_mrna.hist"

    wrap_hist_plot([output_file_exp, output_file_mrna],
                   hist_save="Results/Figures/Figure3/A.combine.png",
                   names=["Experimental", "RefSeq"], to_norm=True)

    ## Run without normalization
    heat_df = heat_plot(output_file_exp + 'MatS', is_norm=False,
                        save_f=output_file_exp + '.noNorm.heat.png')
    heat_df = heat_plot(output_file_mrna + 'MatS', is_norm=False,
                        save_f=output_file_mrna + '.noNorm.heat.png')

    wrap_heat_plot(
        [output_file_exp + 'MatS', output_file_mrna + 'MatS'],
        heat_save="Results/Figures/Figure3/B.combine.png",
        names=["Experimental", "RefSeq"])
    return

