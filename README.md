Companion code for the paper:  
'A Chinese hamster transcription start site atlas that enables targeted editing of CHO cells.'  Isaac Shamie, Sascha H Duttke, Karen J la Cour Karottki, Claudia Z Han, Anders H Hansen, Hooman Hefzi, Kai Xiong, Shangzhong Li, Samuel J Roth, Jenhan Tao, Gyun Min Lee, Christopher K Glass, Helene Faustrup Kildegaard, Christopher Benner, Nathan E Lewis, NAR Genomics and Bioinformatics, Volume 3, Issue 3, September 2021, lqab061, https://doi.org/10.1093/nargab/lqab061
 

#### Data
All sequencing data are submitted to the Gene Expression Omnibus (GEO) with GEO ID GSE159044. You can run the full pipeline using the raw sequencing data, 
The [Supplementary Data](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nargab/3/3/10.1093_nargab_lqab061/1/lqab061_supplemental_files.zip?Expires=1663299647&Signature=erKuj1xbKppkpzaeEhjz4fJfqX3yJKr62OQ-ZhXprLW4ucz0~ceooNWqfEfM84pNTcGkSgh6EQ2Qly03NfBAGenDgvoVWOIwbvayNRCIwoKtCO7Gsgr2RuV~pcM~q4txQsnVqOl3ZbyEipNIkbgG2GmBifbKWFlPcJXfkXvrWLn~f5SHqGszzF8WKyBniQPx2y3Vcm8Nbn~0Dz7OVFp5hWYPa~i-5DcAY8hmF9IMib50XFPiIj56dSKu~~CS5cEUYjF6EWa0vXdP1hsbt9vlTUmKDe-20eOMpQTHXry1bO4050N9QUShHoRWONZ~r4gJLFnFRjEBsnDTI1uqVXHsBw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA) provided in the manuscript is also uploaded to Synapse (synapse.org), with ID syn22969187. This includes our revised protein-coding promoter TSS annotation, in which each of TSS has an associated RefSeq transcript and gene association. This is done for both NCBI RefSeq ([Supplementary Data S2](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nargab/3/3/10.1093_nargab_lqab061/1/lqab061_supplemental_files.zip?Expires=1663299647&Signature=erKuj1xbKppkpzaeEhjz4fJfqX3yJKr62OQ-ZhXprLW4ucz0~ceooNWqfEfM84pNTcGkSgh6EQ2Qly03NfBAGenDgvoVWOIwbvayNRCIwoKtCO7Gsgr2RuV~pcM~q4txQsnVqOl3ZbyEipNIkbgG2GmBifbKWFlPcJXfkXvrWLn~f5SHqGszzF8WKyBniQPx2y3Vcm8Nbn~0Dz7OVFp5hWYPa~i-5DcAY8hmF9IMib50XFPiIj56dSKu~~CS5cEUYjF6EWa0vXdP1hsbt9vlTUmKDe-20eOMpQTHXry1bO4050N9QUShHoRWONZ~r4gJLFnFRjEBsnDTI1uqVXHsBw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)) and with RefSeq in conjunction with the proteogenomics annotation reported in (42) ([Supplementary Data S3](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nargab/3/3/10.1093_nargab_lqab061/1/lqab061_supplemental_files.zip?Expires=1663299647&Signature=erKuj1xbKppkpzaeEhjz4fJfqX3yJKr62OQ-ZhXprLW4ucz0~ceooNWqfEfM84pNTcGkSgh6EQ2Qly03NfBAGenDgvoVWOIwbvayNRCIwoKtCO7Gsgr2RuV~pcM~q4txQsnVqOl3ZbyEipNIkbgG2GmBifbKWFlPcJXfkXvrWLn~f5SHqGszzF8WKyBniQPx2y3Vcm8Nbn~0Dz7OVFp5hWYPa~i-5DcAY8hmF9IMib50XFPiIj56dSKu~~CS5cEUYjF6EWa0vXdP1hsbt9vlTUmKDe-20eOMpQTHXry1bO4050N9QUShHoRWONZ~r4gJLFnFRjEBsnDTI1uqVXHsBw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)). 
Open-chromatin regions merged across samples are provided on synapse as a bed file as well. The genome is taken from NCBI.

#### Required Software
- Homer (http://homer.ucsd.edu/homer/index.html) 
- Snakemake (min version 6)
- numpanpar - (https://github.com/isshamie/parallel_helper) parallel utility package for parallelizing over pandas dataframes and numpy arrays
- Lewis' Lab repo https://github.com/LewisLabUCSD/NGS-Pipeline (branch isaac)

#### Steps to run
1. Create new conda environment (https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). Install required software. 
2. Download repo and install repo package by running
```
git clone git@github.com:isshamie/CH_TSS.git
cd CH_TSS 
pip install -e . #local & editable installation. Can drop the -e`
``` 
** Can run the sequencing pipeline or download the processed data from [this Google Drive folder](https://drive.google.com/drive/folders/1k_ppajLUyKSC1pXJDMDGHN3SLT4haYaU?usp=sharing).  
If you want to run the full pipeline, run steps 3-6. Otherwise downnload data and continue to 7.  
  
<ins>Sequence alignment and peak detection:</ins> 
3\. Download data from GEO Accession GSE159044 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159044)  
4\. Install both the NGS-Pipeline (https://github.com/LewisLabUCSD/NGS-Pipeline) and  https://github.com/kundajelab/atac_dnase_pipelines for ATAC-seq   
5\. Update run_pipeline.sh and run_atac.sh to have the full path to the data.   
6\. Run both `./run_pipeline.sh` and `./run_atac.sh` to run the TSS sequencing and ATAC-seq, respectively.   
 
<ins>Downstream analysis (snakemake) </ins>
7\. The parameter file is parameters/params.yaml. Update parameters/params.yaml to have proper paths.  
8\. Download required software (see above) and Python packages (pandas, snakemake>6.0, matplotlib-venn, scipy)  
9\. Create the new TSS annotation and figures using the snakefile by running
`snakemake -s snakefile`  
10\. Run total-RNA pipeline:
First, use NGS-pipeline to get bam files for the total-RNA data.
Then, get counts and process by runninng the jupyter notebooks 
found in notebooks/csRNA_pipeline in the numbered order. 
Step can be skipped, but won't be able to run a few notebooks below
(F1d and F3b, below)  

11\. Run additional notebooks to finish additional figures.  
<ins>Notebooks:</ins>  
- F1d_sankey.ipynb: Generates the Sankey diagram
- F1e_gene_count.ipynb: Generates Figure 1E, total genes captured
- F2ab_histograms_tss: Read histograms around annotation, using Homer
- F2c_refTSS_Nuc: Generates nucleotide plots around TSSs for F2c.  
- F3a_barplots.01 and F3a_barplots.02: Run in succession. These merge
the bone-marrow WT and 1hKLA and then gets the barplots for number
of genes expressed in X number of tissues.
- F3b_RNASeq_Gene_TPM_CDF: Gets the CDF for the total-RNA-seq genes. 
Requires to run the total-RNA pipeline (#10 above)
- F3ef_homer_motifs: Wrapper to run Homer motif detection
- F4a_silenced_glycosyltransferases: Detection of TSSs for important gene class
- SF2_compare_experiments: Compares the different TSS experiments 
- SF4a_ATAC.ipynb: Open-chromatin around CHO TSSs.  
- SF5_RNA_CHO: RNA-seq expression from CHO public data in genes grouped by TSS status 
- SF6_promoter_usage: Promoter usage in all annotated genes and conserved genes  
*Note that some notebooks may have some hard-coded paths. This should
be minimal, but there may be some re-running while fixing the location
paths. Additionally, extra figures are made which involve varied parameters.
