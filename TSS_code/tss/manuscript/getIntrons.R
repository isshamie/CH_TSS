gff <- snakemake@input[[1]]
out_f <- snakemake@output[[1]]

library('GenomicFeatures')
library(rtracklayer)
#gff <- 'syn20999279_picr.gff'

txdb <- makeTxDbFromGFF(gff, format='gff')
# get intron info
all.introns <- intronicParts(txdb)
export.bed(object=all.introns,con=out_f)