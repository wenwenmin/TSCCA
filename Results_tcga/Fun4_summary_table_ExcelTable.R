load("../My_Data/1_tcga_data/TCGA_Tensor_Data.RData")
load("../My_Data/3_miRTarBase/miRTarBaseData.RData")
load("Fun1_Result_ku10_kv100_kw20_J50.RData")

# 癌症的miRNA和基因
load("../My_Data/geneEdge_CancerGene.RData")
#################################################################################
#################################################################################
# 获取miR-gene模块
module.num = length(out$d)
modules = list();modules.size = c()
for(i in 1:module.num){
  modules[[i]] = list(miRs = which(out$U[,i]!=0),genes = which(out$V[,i]!=0),cancers = which(out$W[,i]!=0),SingularValue_d = out$d[i])
  modules.size = rbind(modules.size,c(length(which(out$U[,i]!=0)), length(which(out$V[,i]!=0)), length(which(out$W[,i]!=0))))
}

# 获取miRNA名称
miR.nameList1 = miR.nameList2 = rep("NA",length(DiffExp_Tensor.attribute$miRNames))
for(i in 1:length(DiffExp_Tensor.attribute$miRNames)){
  temp = unlist(strsplit(DiffExp_Tensor.attribute$miRNames[i], "[|]"))
  miR.nameList1[i] = temp[1]
  miR.nameList2[i] = temp[2]
}
DiffExp_Tensor.attribute$miRNames_C1 = miR.nameList1
DiffExp_Tensor.attribute$miRNames_C2 = miR.nameList2

# 获取gene名称
gene.nameList1 = gene.nameList2 = rep("NA",length(DiffExp_Tensor.attribute$miRNames))
for(i in 1:length(DiffExp_Tensor.attribute$geneNames)){
  temp = unlist(strsplit(DiffExp_Tensor.attribute$geneNames[i], "[|]"))
  gene.nameList1[i] = temp[1]
  gene.nameList2[i] = temp[2]
}
DiffExp_Tensor.attribute$geneNames_C1 = gene.nameList1
DiffExp_Tensor.attribute$geneNames_C2 = gene.nameList2
#################################################################################
#################################################################################
## 一、得到每一个模块中valided_miRGene_pairs
miRTarBase.dat.filter = miRTarBase.dat[miRTarBase.dat$TargetGene.EntrezGeneID%in%DiffExp_Tensor.attribute$geneNames_C2,]
temp0 = unique(miRTarBase.dat.filter[c("miRNA", "TargetGene.EntrezGeneID")])
# 175658 edges
# ----------------------------------------------------------
## 获取每一个模块中的miRNA-gene边，但是基因用的是EntrezGeneID
valid_miRGene_pairs = pvalues = rep(0,module.num)
modules.valid_miRGene_pairs = list()
for(i in 1:module.num){
  temp1 = temp0[temp0$miRNA%in%DiffExp_Tensor.attribute$miRNames_C1[modules[[i]]$miRs],]
  temp2 = temp1[temp1$TargetGene.EntrezGeneID%in%DiffExp_Tensor.attribute$geneNames_C2[modules[[i]]$genes],]
  valid_miRGene_pairs[i] = dim(temp2)[1]  
  modules.valid_miRGene_pairs[[i]] = temp2
}
# ----------------------------------------------------------
## 将基因名称EntrezGeneID格式换成gene symbol格式
modules.miRgenePair = list()
for(j in 1:length(modules.valid_miRGene_pairs)){
  row.miRgenePair = c()
  temp = modules.valid_miRGene_pairs[[j]]
  for(i in 1:dim(temp)[1]){
    temp1 = temp[i,1]
    #TargetGene.EntrezGeneID into Official symbol
    temp2 = gene.nameList1[which(DiffExp_Tensor.attribute$geneNames_C2==temp[i,2])]
    row.miRgenePair = c(row.miRgenePair,paste(temp1,temp2,sep = "->"))
  }
  modules.miRgenePair[[j]] = row.miRgenePair
}
#################################################################################
#################################################################################
## 判断模块中的miRNA和基因是否是癌症Mark
# get cancer genes and miRNAs
cancerGeneList = gene.nameList1[which(gene.nameList1%in%allOnco_May2018$symbol==TRUE)] # 1164
cancermiRList = miR.nameList1[which(tolower(miR.nameList1)%in%miRCancerOctober2017$mirId==TRUE)] #259

# 得到模块基因和小RNA成员
getModuleCancerMarkBinaryVector = function(module){
  miRs = miR.nameList1[module$miRs]
  genes = gene.nameList1[module$genes]
  #使用sapply对miRs向量中的每一个元素判断是否是癌症基因
  CancerMiR = sapply(miRs, function(x) ifelse(x%in%cancermiRList, "Yes", "No"))
  CancerGene = sapply(genes, function(x) ifelse(x%in%cancerGeneList, "Yes", "No"))
  return(list(miRs=miRs,CancerMiR=CancerMiR,genes=genes,CancerGene=CancerGene))
}
# 使用lapply对modules列表中的每一个元素作用getModuleCancerMarkBinaryVector
Res_CancerAnalysis = lapply(modules, getModuleCancerMarkBinaryVector)
#################################################################################
#################################################################################
# 我们用一个Tables存储我们想要的模块特征，一个模块一个小表,最后合并50个小表
Tables = list()
for(i in 1:length(modules)){
  # 我们用一个subTable存储我们想要的模块特征
  cat(paste("module",i,"\n"))
  subTable = data.frame(
    Gene=rep("",100),
    nonzero_v=rep("",100),	
    CancerGene=rep("",100),	
    MiR=rep("",100),	
    nonzero_u=rep("",100),	
    CancerMiR=rep("",100),	
    MGEdge=rep("",100), 
    Cancer=rep("",100),
    nonzero_w=rep("",100),stringsAsFactors=F)
  
  # 得到模块miRNA成员及权重
  miRs = miR.nameList1[modules[[i]]$miRs]
  miRs.weight = out$U[out$U[,i]!=0,i]  
  subTable$MiR[1:length(miRs)] = miRs
  subTable$nonzero_u[1:length(miRs)] = miRs.weight
  
  # 得到模块基因成员及权重
  genes = gene.nameList1[modules[[i]]$genes]
  genes.weight = out$V[out$V[,i]!=0,i]
  subTable$Gene[1:length(genes)] = genes
  subTable$nonzero_v[1:length(genes)] = genes.weight
  
  # 得到模块癌症成员及权重
  cancers = DiffExp_Tensor.attribute$cancerNames[modules[[i]]$cancers]
  cancers.weight = out$W[out$W[,i]!=0,i]
  subTable$Cancer[1:length(cancers)] = cancers
  subTable$nonzero_w[1:length(cancers)] = cancers.weight
  
  # 得到模块验证的边
  M_miRTarget = modules.miRgenePair[[i]]
  subTable$MGEdge[1:length(M_miRTarget)] = M_miRTarget
  
  # 得到模块miRNA和基因是否为癌症基因的二元向量
  MCancerMiR = Res_CancerAnalysis[[i]]$CancerMiR
  MCancerGene = Res_CancerAnalysis[[i]]$CancerGene
  subTable$CancerMiR[1:length(MCancerMiR)] = MCancerMiR
  subTable$CancerGene[1:length(MCancerGene)] = MCancerGene
  
  Tables[[i]] = subTable
}
# *******************************************************************************
# 小面我们希望把这50个表组合成一个大表
# Step1:先让大表等于其中的一个小表
BigTable = Tables[[1]]
flags = rep(paste("M",1,sep=""),ncol(BigTable))
BigTable = rbind(flags,BigTable)
# Step2:组合后面49个模块表
for(i in 2:length(Tables)){
  subTable = Tables[[i]]
  flags = rep(paste("M",i,sep=""),ncol(subTable))
  subTable = rbind(flags,subTable)
  BigTable = cbind(BigTable,subTable)
}
#################################################################################
#################################################################################
# 输出表格
library("openxlsx")
fileName1 = "Summary of Modules"
wb <- createWorkbook()
addWorksheet(wb, fileName1)
writeData(wb, sheet = fileName1,rowNames = F, BigTable)
saveWorkbook(wb,"Res_Fun4_Table/BigTable.xlsx",overwrite=F)


write.table(BigTable, file = "Res_Fun4_Table/BigTable[csvTable].csv", append = FALSE, quote = F, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = F, col.names = TRUE, qmethod = "double")
#################################################################################
#################################################################################