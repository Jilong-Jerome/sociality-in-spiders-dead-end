library(GENESPACE)
library(data.table)
args <- commandArgs(T)
ste <- args[1]
runwd <- file.path(paste("/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/steps/genespace/",ste,sep=""))
gids <- c("LIN","MIM","BI","SARA","TENT","DUM")
#gids <- c("human","chimp","rhesus")
gpar <- init_genespace(
  genomeIDs = gids, 
  speciesIDs = gids, 
  versionIDs = gids, 
  ploidy = rep(1,6),
  #ploidy = rep(1,3),
  wd = runwd, 
  gffString = "gff", 
  pepString = "pep",
  path2orthofinder = "orthofinder", 
  path2mcscanx = "/home/jilong/software/MCScanX",
  rawGenomeDir = file.path(runwd, "rawGenomes"))
parse_annotations(
  gsParam = gpar, 
  gffEntryType = "gene", 
  gffIdColumn = "locus",
  gffStripText = "locus=", 
  headerEntryIndex = 1,
  headerSep = " ", 
  headerStripText = "locus=",
  troubleshoot = TRUE)
gpar <- run_orthofinder(gsParam = gpar)
gpar <- synteny(gsParam = gpar)
g_vec <- c("LIN","SARA","DUM","MIM","TENT","TENT","DUM","BI","TENT","DUM","SARA","TENT","DUM","BI","SARA","LIN","BI","DUM","LIN","SARA","BI","BI","MIM","LIN","SARA","SARA","BI","MIM","LIN","SARA","MIM","BI","DUM","BI","DUM","DUM","SARA","LIN","LIN","LIN","LIN","LIN","LIN")
c_vec <- c("LIN.HiC_1","SARA.HiC_9","DUM.HiC_4","MIM.HiC_6","TENT.HiC_9","TENT.HiC_6","DUM.HiC_5","BI.HiC_11","TENT.HiC_12","DUM.HiC_13","SARA.HiC_3","TENT.HiC_5","DUM.HiC_8","BI.HiC_12","SARA.HiC_14","LIN.HiC_9","BI.HiC_8","DUM.HiC_9","LIN.HiC_10","SARA.HiC_5","BI.HiC_5","BI.HiC_15","MIM.HiC_17","LIN.HiC_11","SARA.HiC_1","SARA.HiC_16","BI.HiC_17","MIM.HiC_16","LIN.HiC_13","SARA.HiC_15","MIM.HiC_5","BI.HiC_4","DUM.HiC_7","BI.HiC_9","DUM.HiC_12","DUM.HiC_14","SARA.HiC_7","LIN.HiC_3","LIN.HiC_5","LIN.HiC_12","LIN.HiC_7","LIN.HiC_21","LIN.HiC_20")
#g_vec <- c("MIM")
#c_vec <- c("MIM.HiC_6")
inversion_tab <- data.frame(genome=g_vec, chr=c_vec)
ripdat <- plot_riparianHits(gsParam = gpar,refGenome = "DUM",chrLabCex = 0.3, invertTheseChrs = inversion_tab,gapProp=0.02,blackBg = FALSE)
