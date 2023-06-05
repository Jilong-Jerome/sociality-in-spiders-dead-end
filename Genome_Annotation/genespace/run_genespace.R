library(GENESPACE)
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
#ripdat <- plot_riparianHits(gpar)
pg <- pangenome(gpar)
