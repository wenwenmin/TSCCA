setwd("C:/AÅÌ-ÏîÄ¿/project-4-TensorSCCA/My_Data/Cancer_marks")

allOnco_May2018 = read.csv("allOnco_May2018.tsv",sep="\t", header=T, stringsAsFactors=F)
miRCancerOctober2017 = read.csv("miRCancerOctober2017.txt",sep="\t", header=T, stringsAsFactors=F)

save(allOnco_May2018,miRCancerOctober2017, file = "Cancer_marks.data.RData")
