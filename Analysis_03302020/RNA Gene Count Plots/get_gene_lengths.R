#!/usr/bin/Rscript
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)


gene_lengths <- function(GTFfile, FASTAfile,id){
  
  #Load the annotation and reduce it
  GTF <- import.gff(GTFfile, format="gtf")#, genome="CriGri", asRangedData=F, feature.type="exon")
  grl <- reduce(split(GTF, elementMetadata(GTF)$id))
  reducedGTF <- unlist(grl, use.names=T)
  elementMetadata(reducedGTF)$id <- rep(names(grl), elementNROWS(grl))
  
  #Open the fasta file
  FASTA <- FaFile(FASTAfile)
  open(FASTA)
  
  #Add the GC numbers
  elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
  elementMetadata(reducedGTF)$widths <- width(reducedGTF)
  
  #Create a list of the ensembl_id/GC/length
  calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
  }
  output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$id), calc_GC_length))
  colnames(output) <- c("Length", "GC")

  
  return(output)
}