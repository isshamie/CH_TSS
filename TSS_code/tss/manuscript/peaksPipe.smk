import os
import glob
from tss.pipeline import create_filenames_dict

#PARAMS = "parameters/{sample}.txt"
#configfile: "config.yaml"

configfile: "parameters/params.yaml"
print("config",config)
results = config['results']
#os.chdir(results)

print(config["peaks_dir"])
filenames = create_filenames_dict()
print(config["META_SAMPLES"])

annotation = {'synapse': '/data/isshamie/genome/hamster/syn20999279_picr.gff',
              'GCF': "/data/isshamie/genome/hamster/ncbi_anno_103/GCF_003668045.1_CriGri-PICR/GCF_003668045.1_CriGri-PICR_genomic.gff"
              }
  # "alt":"/data/isshamie/genome/hamster/ncbi_anno_103/alt_CriGri-PICR_top_level.gff3"}
peak_bins = (-1000, 1000)

rule all:
    input:
        expand("data/processed/{anno}/eTSS/TSS.bed", anno=annotation.keys()),
        expand("figures/{anno}/barplots_eTSS/{merge}/transcript_tissues_expressed.png", anno=list(annotation.keys()), merge=config["MERGE_TISSUES"]),
        expand("figures/{anno}/nucleotides/eTSS_nucleotide_compare.png", anno=list(annotation.keys()))


def get_input(wildcards):
    return annotation[wildcards.anno]


rule process_genome:
    input:
        anno =  get_input,
        genome =  config["GENOME_FA"]
        #is_gff = config["is_gff"]
    output:
        mRNA_peak = "data/processed/{anno}/genome/mRNA.peak",
        start_mRNA = "data/processed/{anno}/genome/start_site_mRNA.tsv"
    params: filenames
    shell: "python tss/data/generate_genome.py {input.genome} {input.anno} {output.mRNA_peak}"


rule merge_peaks:
    input:
        start_mRNA = "data/processed/{anno}/genome/start_site_mRNA.tsv",
        ref_anno = get_input
    output:
        merged_f = "data/processed/{anno}/merged/samples.merge"
    log: 'data/processed/{anno}/merged/merged.log'
    params:
        peaks_dir = config["peaks_dir"],
        genome = config["GENOME_FA"]
    shell: "python tss/data/merging.py {input} {output.merged_f} {params} to_merge &> {log}"


rule peak_source:
    input:
    output:


rule meta_score:

rule meta_anno:
    input:
    output:


rule run_filters:
    input:
    output:


rule output_and_visualize:
    input:
    output:



rule peak_expression_and_location:
    input:
        start_mRNA = "data/processed/{anno}/genome/start_site_mRNA.tsv",
        ref_anno = get_input,
        merged_f = "data/processed/{anno}/merged/samples.merge"
    output:
        peak_f = "data/processed/{anno}/merged/peaks_with_tss_distances_size1.noCDS.peak",
        out_f = "data/processed/{anno}/merged/files_used.json",
        expr_peaks_f = "data/processed/{anno}/merged/samples.merge.peaksexpression.log10",
    log: 'data/processed/{anno}/merged/anno_and_exp.log'
    params:
        peaks_dir = config["peaks_dir"],
        genome = config["GENOME_FA"]
    shell: "python tss/data/merging.py {input.start_mRNA} {input.ref_anno} {output.peak_f} {params} to_anno_and_exp &> {log}"
        #out = os.path.join(config["DATA_PROCESSED_DIR"],output[0])


rule peak_merge_and_preprocess:
    pass


rule peak_anno:
    input: rules.merge_peaks.input
    output: "data/processed/{anno}/merged/peaks_anno",
    shell: ""

rule gene_centric_tss:
    input:
        merged_out = "data/processed/{anno}/merged/files_used.json",
        start_mRNA = "data/processed/{anno}/genome/start_site_mRNA.tsv",
        meta_f = config["META_SAMPLES"]
    output:
        txn_f = "data/processed/{anno}/gene_centric_tss/txn_df.p",
        out_f = "data/processed/{anno}/gene_centric_tss/files_used.json",
        txn_tissues = "data/processed/{anno}/gene_centric_tss/txn_df_tissues.p"
    log: 'data/processed/{anno}/gene_centric_tss/gene_centric.log'
    params:
        allow_intron = False,
        peak_bin_l = -1000,
        peak_bin_r = 1000
    shell: 'python tss/data/gene_centric.py {output.txn_f} {input} -- {params} &> {log}'


rule tissues_collapse:
    input:
        merged_out = "data/processed/{anno}/merged/samples.merge.peaksexpression.log10",
    params:
        meta = config["META_SAMPLES"]
    output:
        tissues_expr_f = "data/processed/{anno}/tissues_collapse/tissues.merge.peaksexpression.log10"
    shell: "python tss/data/tissues_collapse.py {output} {params.meta} {input}"


rule retrieve_peaks:
    input:
        gene_centric = "data/processed/{anno}/gene_centric_tss/files_used.json",
        merged_f = "data/processed/{anno}/merged/samples.merge",
        expr_peaks_f = "data/processed/{anno}/merged/samples.merge.peaksexpression.log10",
        peaks_dir = config["peaks_dir"],
    output: "data/processed/{anno}/retrieve_peaks/files_used.json"
    shell: "python tss/data/retrieve_tss_peaks.py {output} {input}"
    #log: 'data/processed/{anno}/retrieve_peaks/retrieve_peaks.log'



rule create_output:
    input:
        gene_retrieve_peaks = "data/processed/{anno}/retrieve_peaks/files_used.json",
        gene_centric_out_f = "data/processed/{anno}/gene_centric_tss/txn_df_tissues.p",
        reference_out = "data/processed/{anno}/genome/mRNA.peak",
        cho_atac_f =  config["atac_f"],
        all_atac = config["all_atac"],
        tissues_expr_f = "data/processed/{anno}/tissues_collapse/tissues.merge.peaksexpression.log10"
    output:
          tss_bed = "data/processed/{anno}/eTSS/TSS.bed",
          refseq_center =  "data/processed/{anno}/eTSS/refseq_centered.TSS.bed",
          eTSS_exp_bed = "data/processed/{anno}/eTSS/TSS.exp.bed"
    log: 'data/processed/{anno}/eTSS/eTSS.log'
    shell:
        "python tss/data/wrap_create_eTSS.py {output.tss_bed} {input}  &> {log}"


rule barplot:
    input:
        bed_eTSS =  "data/processed/{anno}/eTSS/TSS.exp.bed",
        gene_centric_out = "data/processed/{anno}/gene_centric_tss/files_used.json"
    output:
          "figures/{anno}/barplots_eTSS/{merge}/transcript_tissues_expressed.png"
    log: "figures/{anno}/barplots_eTSS/{merge}/barplots.log"
    params: meta_samples = config['META_SAMPLES']
    shell: "python tss/data/tss_barplots.py {output} {input} {params} {wildcards.merge} &> {log}"


rule nucleotides:
    input:
        bed_eTSS =  "data/processed/{anno}/eTSS/TSS.exp.bed"
    output:
        nuc_plot = "figures/{anno}/nucleotides/eTSS_nucleotide_compare.png"
    log: "figures/{anno}/nucleotides/nucleotides.log"
    shell: "python tss/data/nucleotide_plot.py {input} {output} {config[GENOME_FA]} &> {log}"


rule distance_to_anno:
    input:
        "data/processed/{anno}/eTSS/TSS.bed",
         "data/processed/{anno}/genome/start_site_mRNA.tsv"
    output: "figures/{anno}/distance_to_anno/distance_to_tss.png"
    log: "figures/{anno}/distance_to_anno/distance_to_anno.log"
    shell: "python tss/data/distance_to_anno.py {input} {output} &> {log}"

# rule histograms:
#     input:
#         tag_dir =
#     output:
#         histograms_dir =
#     log: "data/processed/{anno}/histograms/compare_tags_eTSS"
#     shell: "python tss/data/tags_on_annotation {input} {output} &> {log}"
#

rule motifs:
    input:
        meta_eTSS = "data/processed/{anno}/eTSS/TSS.meta.tsv",
        bed_eTSS = "data/processed/{anno}/eTSS/TSS.bed",
        refseq_center = "data/processed/{anno}/eTSS/refseq_centered.TSS.bed"
    output:
        "figures/motifs/eTSS_inr_tata.png"
