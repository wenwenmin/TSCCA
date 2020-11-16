setwd("C:/A盘-项目/project-4-TensorSCCA/My_Data/miRTarBase")

miRTarBase.dat = read.csv("C:/A盘-项目/project-4-TensorSCCA/My_Data/miRTarBase/hsa_MTI.csv", header = TRUE,sep = ",")

colnames(miRTarBase.dat) = 
  c("miRTarBase.ID",
    "miRNA",
    "Species.miRNA",
    "TargetGene",
    "TargetGene.EntrezGeneID",
    "Species.TargetGene",
    "Experiments",
    "SupportType",
    "References.PMID")

miRTarBase.dat  = data.frame(apply(miRTarBase.dat, 2, as.character),stringsAsFactors=FALSE)
miRTarBase.dat$TargetGene.EntrezGeneID = gsub(" ", "", miRTarBase.dat$TargetGene.EntrezGeneID)#去掉空格

save(miRTarBase.dat, file = "miRTarBaseData.RData") 