# Biosignature analysis
# The R script is created by Jyotirmoy Das
# The package is licensed under GNL 3.0, free to distribute
# Copyright (c) Jyotirmoy Das
## For Illumina 450K data analysis
DiffMeth450K <- function(diffmeth){
  {
    message("[===========================]")
    message("[<<<< Analysis START for 450K data>>>>>]")
  }
  
  #if (!requireNamespace("BiocManager", quietly = TRUE))
  #    install.packages("BiocManager")
  
  #BiocManager::install("ChAMP")
  #BiocManager::install(c("minfi","ChAMPdata","Illumina450ProbeVariants.db","sva","IlluminaHumanMethylation450kmanifest","limma","RPMM","DNAcopy","preprocessCore","impute","marray","wateRmelon","goseq","plyr","GenomicRanges","RefFreeEWAS","qvalue","isva","doParallel","bumphunter","quadprog","shiny","shinythemes","plotly","RColorBrewer","DMRcate","dendextend","IlluminaHumanMethylationEPICmanifest","FEM","matrixStats","missMethyl","combinat"))
  {
    message("[===========================]")
    message("[<<<< Loading ChAMP package >>>>>]")
  }
  suppressMessages(library("ChAMP"))
  require(tcltk)
  choose_directory = function(caption = 'Select data directory') {
    if (exists('utils::choose.dir')) {
      choose.dir(caption = caption)
    } else {
      tk_choose.dir(caption = caption)
    }
  }
  setwd(choose_directory())
  myLoad <- champ.load(arraytype="450K")
  #champ.QC()
  # Normalize the file with BMIQ
  myNorm <- champ.norm(beta=myLoad$beta, rgSet=myLoad$rgSet,
                       mset=myLoad$mset,
                       resultsDir="./CHAMP_Normalization/", method="BMIQ",
                       plotBMIQ=FALSE,
                       arraytype="450K",
                       cores=24)
  # Save the normalized file
  write.table(myNorm, "norm_excludeExp13_HLADR_Exp-Control.txt",
            sep = "\t",
            quote = F,
            row.names = T)
  # Calculate the Differential Methylation CpG sites
  myDMC <- champ.DMP(beta = myNorm,
                     pheno = myLoad$pd$Sample_Group,
                     compare.group = NULL,
                     adjPVal = 0.05,
                     adjust.method = "BH",
                     arraytype = "450K")
  
  write.table(myDMC[[1]], "DMCs_excludeExp13_HLADR_Exp-Control.txt",
              sep = "\t",
              quote = F,
              row.names = T)
  
  # Calculate the Differentially Methylated Regions
#  myDMR <- champ.DMR(beta=myNorm, pheno=myLoad$pd$Sample_Group,
#                     compare.group=NULL, arraytype="450K",
#                     method = "Bumphunter",
#                     minProbes=7,
#                     adjPvalDmr=0.05,
#                     cores=6)
#  write.table(myDMR[[1]], "DMRs_450K.txt",
#              sep = "\t",
#              quote = F,
#              row.names = F)
}
