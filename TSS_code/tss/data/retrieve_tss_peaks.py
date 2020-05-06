import glob
from tss.data import annotation
import os
from os.path import basename, join, dirname
import click
from tss import pipeline
import logging
import json


def run(p, outdir, peaks_dir, merged_f,expr_peaks_f):
    """
    :param outdir:
    :param peaks_dir:
    :param anno_gene_f:
    :param merged_f:
    :param expr_peaks_f:
    :return:
    """
    #outdir = join(outdir, "tss_annotation_peaks/")
    #anno_f = join(outdir, "tss_annotation/gene_df_02.p")
    #merged_f = join(outdir, "merged/samples.merge")
    anno_gene_f  =  p["gene_centric"]["gene centric peaks from samples"]

    output_f = join(outdir,"all_peaks_gene_df.tsv")
    annotation.retrieve_all_peaks_from_anno(anno_gene_f, merged_f, output_f)

    output_f = join(outdir, "gene_df_maxVal_all_peaks.tsv")
    annotation.retrieve_all_peaks_from_anno(anno_gene_f, merged_f, output_f,
                                            is_max=True)

    for sample_name in glob.glob(peaks_dir + "/*"):
        # Merged sample file
        output_f = os.path.join(outdir, "merged_" + os.path.basename(
            sample_name) + '.tsv')
        annotation.retrieve_sample_peaks_from_anno(anno_gene_f, merged_f,
                                                   output_f,
                                                   sample_name,
                                                   sample_f=None,
                                                   col_name=None,
                                                   use_sample_peaks=False)
        # Original sample file
        output_f = os.path.join(outdir, "sample_" + os.path.basename(
            sample_name) + '.tsv')
        annotation.retrieve_sample_peaks_from_anno(anno_gene_f, merged_f,
                                                   output_f,
                                                   sample_name,
                                                   sample_f=sample_name,
                                                   col_name=None,
                                                   use_sample_peaks=True,
                                                   peaks_dir=peaks_dir)

    anno_txn_f = p["gene_centric"]["transcript centric peaks from samples"]
    all_peaks_txn_f = join(outdir,"all_peaks_txn_df.tsv")
    annotation.retrieve_all_peaks_from_anno(anno_txn_f, merged_f, all_peaks_txn_f,)

    #expr_peaks_f = "Results/merged/samples.merge.peaksexpression.log10"
    out_f = join(outdir,"all_peaks_txn_df_maxinfo.tsv")
    annotation.wrap_add_max_info(all_peaks_txn_f, expr_peaks_f, out_f,
                                 peaks_dir=peaks_dir)

    return



CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('retrieve_peaks_log_f', type=click.Path(exists=False))
@click.argument('gene_centric_out_f', type=click.Path(exists=True))
@click.argument('merged_f', type=click.Path(exists=True))
@click.argument('expr_peaks_f', type=click.Path(exists=True))
@click.argument('peaks_dir', type=click.Path(exists=True))
def main(retrieve_peaks_log_f, merged_f,gene_centric_out_f, expr_peaks_f, peaks_dir):

    p = pipeline.create_filenames_dict()
    p = pipeline.create_fullnames(p,'gene_centric', dirname(gene_centric_out_f))

    outdir = dirname(retrieve_peaks_log_f)
    p = pipeline.create_fullnames(p, 'gene_retrieve_peaks',outdir)


    run(p,outdir, peaks_dir, merged_f,expr_peaks_f)
    # Save the parameters
    with open(p["gene_retrieve_peaks"]["files_used"], 'w') as fp:
        json.dump(p, fp)

    return


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    main()