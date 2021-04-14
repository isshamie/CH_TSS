# import click
# import logging
# from pathlib import Path
# from tss.utils.config import read_config_file, check_required
# import tss.data.merging as merging
# import tss.data.tissues_collapse as tissues_collapse
# import tss.data.gene_centric as gene_centric
# import tss.data.output_peaks as output_peaks
# import tss.data.atac as atac
import numpy as np
import shutil
from os.path import join, exists, basename
import os
import glob
# merge_samples

# Pipeline:
# 1. merge_samples
# 	-merge_peaks
# 	-annotate_peaks
# 	-bash: homerTools extract
# 	-bash: makeTagDirectory
# 	-make_bedgraph_file
# 2. peaks_break_into_files
#
# 3. peaks_to_annotation.ipynb
#
# 4. Merge_into_tissues
#
# 5. retrieve_tss_peaks
#
# 6. Rerun Figures2
#
# 7. ATAC: Create species
#
# 8. ATAC: run pipeline
#
# 9. ATAC: Copy folders
#
# 10. Separate peaks based on expression  (needed?)
#
# 11. create_output
#
# 12. Rerun Figure3
#
# 13. Tag Histograms Compare Experimental to RefSeq
#
# 14. Nucleotide Frequency
#
# 15. Rerun  CDF offi TSS with Samples.ipynb
#
# 16. Compare RNA_seq
#
# 17. Igv snapshot


def create_fullnames(p, stage, outdir):
    """Takes the stage of interest, and for each of the filenames
    to be created, make it an absolute value with the outdir"""
    if stage not in p:
        print(f"Stage {stage} not here")
        return {}
    for f in p[stage]:
        p[stage][f] = join(outdir, p[stage][f])
    return p


def create_filenames_dict():
    params = dict()
    params["reference"] = dict()
    params["reference"]["mRNA_peak_f"] = "mRNA.peak"
    params["reference"]["start_site_mRNA"] = "start_site_mRNA.tsv"
    params["reference"]["CDS_peak_f"] = "CDS.peak"
    params["reference"]["mRNA bed"] = "mRNA.bed"
    params["reference"]["CDS bed"] = "CDS.bed"
    # Merged
    params["merged"] = dict()
    params["merged"]["merged peak samples"] = "samples.merge"
    params["merged"]["merged peak samples annotation"] = "samples.merge.anno"
    params["merged"]["merged peak samples sequences"] = "samples.merge.anno.fa"

    params["merged"]["merged tags cap"] = "tags.cap.merged"
    params["merged"]["merged tags input"] = "tags.input.merged"
    params["merged"]["merged tags cap bedgraph"] = "tags.cap.merged.bedgraph"
    params["merged"]["merged tags input bedgraph"] = "tags.cap.merged.bedgraph"

    params["merged"]["peak samples expression"] = "samples.merge.peaksexpression"
    params["merged"]["peak samples expression log10"] = "samples.merge.peaksexpression.log10"
    params["merged"]["peak samples expression log2"] = "samples.merge.peaksexpression.log2"
    params["merged"]["peak samples minimal"] = "samples.merge.minimal"
    params["merged"]["peak samples minimal anno"] = "samples.merge.minimal.anno"

    params["merged"]["peak samples distance to TSS"] = "peaks_with_tss_distances.tsv"
    params["merged"]["peaks_size1 peak"] = "peaks_with_tss_distances_size1.peak"
    params["merged"]["peaks_size1 bed"] = "peaks_with_tss_distances_size1.bed"
    params["merged"]["no cds bed"] =  "peaks_with_tss_distances_size1.noCDS.bed"
    params["merged"]["no cds peak"] = "peaks_with_tss_distances_size1.noCDS.peak"
    params["merged"]["files_used"] = 'files_used.json'

    #Flip to make gene table
    params['gene_centric'] = dict()
    params['gene_centric']["gene centric peaks from samples"] = "gene_df.p"
    params['gene_centric']["gene centric peaks from samples plus tissues"] = "gene_df_tissues.p"
    params['gene_centric']["gene centric peaks from samples max"] = "gene_df.p"
    params['gene_centric']["transcript centric peaks from samples"] = "txn_df.p"
    params['gene_centric']["transcript centric peaks from samples plus tissues"] = "txn_df_tissues.p"
    params['gene_centric']["transcript centric peaks matrix"] = "txn_df_matrix"
    params['gene_centric']["transcript centric peaks from samples max"] = "max_all_peaks_txn_df.p"
    params["gene_centric"]["files_used"] = 'files_used.json'


    # Collaps annotation onto tissues
    params["tissues_collapse"] = dict()
    params["tissues_collapse"]["peak tissues expression"] = "tissues.merge.peaksexpression"
    params["tissues_collapse"]["peak tissues expression log10"] = "tissues.merge.peaksexpression.log10"
    params["tissues_collapse"]["peak tissues expression log2"] = "tissues.merge.peaksexpression.log2"

    # elif stage == "retrieve_tss_peaks":
    params['gene_retrieve_peaks'] = dict()
    params['gene_retrieve_peaks']["gene_all_peaks"] = "all_peaks_gene_df.tsv"
    params['gene_retrieve_peaks']["gene_all_peaks_max"] = "gene_df_maxVal_all_peaks.tsv"
    params['gene_retrieve_peaks']["txn_all_peaks"] = "all_peaks_txn_df.tsv"
    params["gene_retrieve_peaks"]["txn_maxInfo.tsv"] = "all_peaks_txn_df_maxinfo.tsv"
    params["gene_retrieve_peaks"]["files_used"] = 'files_used.json'


    # Create the output eTSS with names
    params["output_peaks"] = dict()
    params["output_peaks"]["eTSS bed"] = "TSS.bed"
    params["output_peaks"]["eTSS meta"] = "TSS.meta.tsv"
    params["output_peaks"]["eTSS exp bed"] = "TSS.exp.bed"
    params["output_peaks"]["eTSS exp meta"] = "TSS.exp.meta.tsv"
    params["output_peaks"]["eTSS refseq_centered bed"] = "refseq_centered.TSS.bed"
    params["output_peaks"]["files_used"] = 'files_used.json'

    # Get the number of TSS captured
    params["tss_results"] = dict()
    params["tss_results"]["gene fraction unique"] = "gene_df_tissues_frac_unique"
    params["tss_results"]["txn fraction unique"] = "txn_df_tissues_frac_unique"
    params["tss_results"]["gene venn"] = "gene_df_tissues_cho_venn"
    params["tss_results"]["txn venn"] = "txn_df_tissues_cho_venn"

    params["barplots"] = dict()
    params["barplots"]["gene"] = "gene_tissues_expressed.png"
    params["barplots"]["transcript"] = "transcript_tissues_expressed.png"
    params["barplots"]["files_used"] = 'files_used.json'

    # Calculate nucleotide plots before vs after
    params["nucleotide_out"] = dict()
    params["nucleotide_out"]["eTSS exp seq"] = "TSS.exp.fa"
    params["nucleotide_out"]["eTSS refseq_centered seq"] = "TSS.refseq_centered.fa"
    params["nucleotide_out"]["compare_old_new"] = "eTSS_nucleotide_compare.png"
    # Calculate the histogram of the peaks on the TSS
    params["histograms"] = dict()

    return params

#
# def list_stage_filenames(stage):
#     params = dict()
#     if stage == "reference":
#         params["mRNA_peak_file"] = "mRNA.peak"
#         params["start_site_mRNA"] = "start_site_mRNA.tsv"
#
#     elif stage == 'merged':
#         params["merged peak samples"] = "samples.merge"
#         params["merged peak samples annotation"] = "samples.merge.anno"
#         params["merged peak samples sequences"] = "samples.merge.anno.fa"
#
#         params["merged tags cap"] = "tags.cap.merged"
#         params["merged tags input"] = "tags.input.merged"
#         params["merged tags cap bedgraph"] = "tags.cap.merged.bedgraph"
#         params["merged tags input bedgraph"] = "tags.cap.merged.bedgraph"
#
#         params["peak samples expression"] = "samples.merge.peaksexpression"
#         params["peak samples expression log10"] = "samples.merge.peaksexpression.log10"
#         params["peak samples expression log2"] = "samples.merge.peaksexpression.log2"
#         params["peak samples minimal"] = "samples.merge.minimal"
#         params["peak samples minimal anno"] = "samples.merge.minimal.anno"
#
#         params["peak samples distance to TSS"] = "peaks_with_tss_distances.tsv"
#         params["peaks_size1 peak"] = "peaks_with_tss_distances_size1.peak"
#         params["peaks_size1 bed"] = "peaks_with_tss_distances_size1.bed"
#         params["no cds bed"] =  "peaks_with_tss_distances_size1.noCDS.bed"
#         params["no cds peak"] = "peaks_with_tss_distances_size1.noCDS.peak"
#
#     elif stage == 'gene_centric':
#         params["gene centric peaks from samples"] = \
#             "all_peaks_gene_df.p"
#         params["gene centric peaks from samples plus tissues"] = \
#             "all_peaks_gene_df_tissues.p"
#         params["gene centric peaks from samples max"] = \
#             "max_all_peaks_gene_df.p"
#         params["transcript centric peaks from samples"] \
#             ="all_peaks_txn_df.p"
#         params["transcript centric peaks from samples plus tissues"] \
#             = "all_peaks_txn_df_tissues.p"
#         params["transcript centric peaks matrix"] = "txn_df_matrix"
#         params["transcript centric peaks from samples max"] = \
#             "max_all_peaks_txn_df.p"
#     #elif stage == "retrieve_tss_peaks":
#
#     elif stage == "tissues_collapse":
#         params["peak tissues expression"] = "tissues.merge.peaksexpression"
#         params["peak tissues expression log10"] = "tissues.merge.peaksexpression.log10"
#     elif stage == "output_peaks": #Create Output peaks
#         params["final peaks"] = "TSS"
#
#     elif stage == "tss_results":
#         params["gene fraction unique"] = "gene_df_tissues_frac_unique"
#         params["txn fraction unique"] = "txn_df_tissues_frac_unique"
#
#         params["gene venn"] = "gene_df_tissues_cho_venn"
#         params["txn venn"] = "txn_df_tissues_cho_venn"
#     # elif stage == "histogram":
#     #     params["histograms"]
#     return params
#
#
# def set_absolute(config, stage):
#     """
#     Sets absolute file name for stage and the filenames
#     :param config:
#     :param stage:
#     :return:
#     """
#
#     if not config[stage]["prior"]:
#         config[stage]["absolute"] = join(config["global"]["RESULTS"],
#                                          config["global"]["PREFIX"], stage, config[stage]["folder"])
#     else:
#         config[stage]["absolute"] = config[stage]["folder"]
#
#     filenames_dict = list_stage_filenames(stage)
#     config[stage]["filenames"] = dict()
#     for f in filenames_dict:
#         config[stage]["filenames"][f] = join(config[stage][
#                                                  "folder"],
#                                              filenames_dict[f])
#     return config[stage]
#
#
# def check_dir_filenames(folder, filenames):
#     for f in filenames:
#         if not exists(join(folder, f)):
#             return False
#     return True
#
#
# def check_dir_params(curr_params, old_params_folder):
#     old_params_f = join(old_params_folder, 'params_used.json')
#     if not os.path.exists(old_params_f):
#         return False
#     old_params = read_config_file(old_params_f)
#     for i in old_params:
#         if i in curr_params and not (old_params[i] == curr_params[i]):
#             return False
#     for i in curr_params:
#         if i in old_params and (not old_params[i] == curr_params[i]):
#             return False
#     return True
#
#
# def set_stage_config(config, stage):
#     """
#     Sets up the config for that stage. Checks to see if there are any
#     filenames, and if there
#     are,
#     see if all exist. If they dont, return False. Otherwise, return
#     True.
#     :param stage:
#     :param config:
#     :return:
#     """
#     stage_config = config[stage]
#     filenames = list_stage_filenames(stage)
#     folder = stage_config["folder"]
#     make_new = True
#     # Check if folder is overwrittern and use that if it is.
#     if folder is not None:
#         if check_dir_filenames(folder, filenames):
#             stage_config["folder"] = folder
#             make_new = False
#     else:
#         # Look through for numbered dirs in current folder, and check if
#         # the params are the same. If they are then use that
#         folder = join(config["global"]["RESULTS"], config[
#             "global"]["PREFIX"], stage)
#         if not os.path.exists(folder):
#             os.mkdir(folder)
#         else:
#             subfolders = [f.path for f in os.scandir(folder) if f.is_dir()]
#             for i in subfolders:
#                 if check_dir_params(stage_config["params"], i):
#                     stage_config["folder"] = i
#                     make_new = False
#     # Otherwise, create a new number with updated filenames
#     if make_new:
#         subfolders = np.sort([f.path for f in os.scandir(folder) if \
#                 f.is_dir()])
#         previous_folders = list(map(lambda x: str(basename(x)),
#                                     subfolders))
#         max_num = 0
#         for p in subfolders:
#             curr_num = basename(p)
#             if curr_num.isnumeric():
#                 if not check_dir_filenames(folder, filenames) and \
#                         config["global"]["clean"]:
#                     shutil.rmtree(p)
#                     max_num = int(curr_num)-1
#                     break
#                 max_num = max(max_num, int(curr_num))
#         stage_config["folder"] = os.path.join(config["global"]["RESULTS"], config[
#             "global"]["PREFIX"], stage, str(max_num+1))
#         os.makedirs(stage_config["folder"])
#         stage_config["prior"] = False
#     else:
#         stage_config["prior"] = True
#     return stage_config
#
#
# def create_params(config):
#     """
#     Sets up all the files and directories to be created for the
#     pipeline. Returns a dictionary of all the file names.
#     If the file is in the params dictionary initially that overwrites the default location.
#     :param config:
#     :return:
#     """
#     print(config["stages"])
#     for stage in config["stages"]:
#
#         #if os.path.exists()
#         config[stage] = set_stage_config(config, stage)
#         config[stage] = set_absolute(config, stage)
#
#     # merged_dir = join(out, config["merged"])
#     # tss_annotation = join(out, config["tss annotation"])
#     # tss_annotation_peaks = join(out, config["tss annotation peaks"])
#     # final_peaks = join(out, config["final peaks"])
#     return config
#
#
# def run(config, stage):
#     if stage == "merged":
#         merging.run(config)
#     elif stage == "gene_centric":
#         gene_centric.run(config)
#     elif stage == "tissues_collapse":
#         tissues_collapse.run(config)
#     elif stage == "output_peaks":
#         output_peaks.run(config)
#     elif stage == "atac":
#         atac.run(config)
#     return
#
#
# CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
# # @click.command(context_settings=CONTEXT_SETTINGS)
# # # run settings
# # @click.option("-P", "--prefix", default=None, help="Job prefix")
# # @click.option("-N", "--cores", default=None, help="Number of cores for batch jobs", type=int)
# # @click.option(
# #     "-y", "--yolo", default=False, is_flag=True, help="Disable overwrite protection"
# # )
# @click.command(context_settings=CONTEXT_SETTINGS)
# @click.argument('parameter')#, type=click.Path(exists=True))
# #@click.argument('output_prefix', type=click.Path(), default='test')
# def main(parameter):
#     """ Runs data processing scripts to turn raw data from (../raw) into
#         cleaned data ready to be analyzed (saved in ../processed).
#     """
#     logger = logging.getLogger(__name__)
#     logger.info('making final data set from raw data')
#     print('parameter', parameter)
#     config = read_config_file(parameter)
#     print(config)
#     check_required(config["global"], ["PREFIX", "DATA_DIR", "RESULTS"])
#
#
#     # Create the main folder if not made
#     main_folder = os.path.join(config["global"]["RESULTS"], config[
#         "global"]["PREFIX"])
#     if not os.path.exists(main_folder):
#         os.mkdir(main_folder)
#
#     # Create the folders needed
#     config = create_params(config)
#     print(config)
#
#     for stage in config["stages"]:
#         print("Running", stage)
#         if config[stage]["prior"] is None or not config[stage]["prior"]:
#
#             # Save the config and run
#             run(config, stage)
#         elif "overwrite" in config[stage] and config[stage]["overwrite"]:
#             # Save the config and run
#             # Need to overwrite as well
#             run(config, stage)
#
#
# if __name__ == '__main__':
#     log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
#     logging.basicConfig(level=logging.INFO, format=log_fmt)
#     # not used in this stub but often useful for finding various files
#     project_dir = Path(__file__).resolve().parents[2]
#     print(project_dir)
#
#     #load_dotenv(find_dotenv())
#     main(['/data/isshamie/TSS/TSS/parameters/ncbi_anno_103.yaml'])
#     #main()
#
