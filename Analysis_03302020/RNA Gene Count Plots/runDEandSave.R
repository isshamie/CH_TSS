library("EnhancedVolcano")
library("dplyr")

runDEandSave <- function(dds, fsave,contrast_v){
  print("Running DESeq")
  print(fsave)
  
  # Run DESeq
  if (file.exists(paste0("quantification/dds_",fsave,".rds"))){
    print("Already run de, loading")
    dds <- readRDS(paste0("quantification/dds_",fsave,".rds"))
    
  }else{
    dds <- DESeq(dds,parallel = T)  
    # Save DESeq object
    saveRDS(object = dds,file = paste0("quantification/dds_",fsave,".rds"))
  }
  
  # Save counts 
  counts(dds, normalized=TRUE) %>% as.data.frame() %>%
    write.csv(paste0("quantification/",fsave,"counts_normalized.csv"))
  
  # Save VSD
  vsd <- varianceStabilizingTransformation(dds)
  assay(vsd) %>% as.data.frame() %>% 
    write.csv(paste0("quantification/",fsave,"_counts_variance_controlled.csv"))
  
  
  ####
  # Get contrast results and print summary
  
  #########################################
  # First see if there are multiple contrasts
  for (row in 1:nrow(contrast_v)){ 
    print(contrast_v[row,])
    res <- results(dds, contrast = contrast_v[row,])
    summary(res)
    
    # Save results as csv and tsv
    res %>% as.data.frame %>% na.omit() %>%
      as.data.table(keep.rownames = "Gene") %>% 
      fwrite(file=paste0("DE_out/res_",fsave,"_",contrast_v[row,2],"__",contrast_v[row,3],".csv"))
    
    res %>% as.data.frame %>% na.omit %>%
      as.data.table(keep.rownames = "Gene") %>% 
      fwrite(file=paste0("DE_out/res_",fsave,"_",contrast_v[row,2],"__",contrast_v[row,3],".tsv"),sep="\t")
    
    
    # Remove the na p-adjusted numbers and save de genes (ones less than 0.1 adjusted)
    res %>% as.data.frame %>% na.omit %>%  as.data.table(keep.rownames = "Gene") %>% filter(padj < 0.1) %>% arrange(padj) %>%
      fwrite(file=paste0("quantification/genes_de_",fsave,"_",contrast_v[row,2],"__",contrast_v[row,3],".csv"))
    
    # Save volcano
    # png(paste0("Figures/",fsave,"_",contrast_v[row,2],"__",contrast_v[row,3],"_volcano.png"))
    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    xlim = c(-2.5, 2.5),
                    ylim = c(0,5))
    ggsave(paste0("Figures/",fsave,"_",contrast_v[row,2],"__",contrast_v[row,3],"_volcano.png"))
    
    #dev.off()
    
    # Save MA
    png(paste0("Figures/",fsave,"_",contrast_v[row,2],"__",contrast_v[row,3],"MAplot.png"))
    plotMA(res, ylim=c(-2,2))
    dev.off()
    
  }
  
}
  

  
  