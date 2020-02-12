

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
# 15. Rerun  CDF of TSS with Samples.ipynb
#
# 16. Compare RNA_seq
#
# 17. Igv snapshot
#
