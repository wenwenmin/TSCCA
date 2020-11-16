load("Fun1_Result_ku10_kv100_kw20_J50.RData")
load("../My_Data/1_tcga_data/TCGA_Tensor_Data.RData")
load("../My_Data/geneEdge_CancerGene.RData")

# modules.valid_miRGene_pairs
# load("Result.FigTable/Fun6_Result_modules.valid_miRGene_pairs.RData")

#######################################################################################################
#######################################################################################################
# 50 miRNA-gene modules
ModuleNum = length(out$d)
modules = list()
for(i in 1:ModuleNum){
  modules[[i]] = list(miRs = which(out$U[,i]!=0),genes = which(out$V[,i]!=0),cancers = which(out$W[,i]!=0),SingularValue_d = out$d[i])
}

# miRNAs
miR.nameList1 = miR.nameList2 = rep("NA",length(DiffExp_Tensor.attribute$miRNames))
for(i in 1:length(DiffExp_Tensor.attribute$miRNames)){
  temp = unlist(strsplit(DiffExp_Tensor.attribute$miRNames[i], "[|]"))
  miR.nameList1[i] = temp[1] # hsa-let-7a-5p....
  miR.nameList2[i] = temp[2] # MIMAT0000062....
}

# genes
gene.nameList1 = gene.nameList2 = rep("NA",length(DiffExp_Tensor.attribute$miRNames))
for(i in 1:length(DiffExp_Tensor.attribute$geneNames)){
  temp = unlist(strsplit(DiffExp_Tensor.attribute$geneNames[i], "[|]"))
  gene.nameList1[i] = temp[1] # A1BG....
  gene.nameList2[i] = temp[2] # 53947...
}

#######################################################################################################
#######################################################################################################
# **get cancer genes and miRNAs
cancerGeneList = gene.nameList1[which(gene.nameList1%in%allOnco_May2018$symbol==TRUE)] # 1164
cancermiRList = miR.nameList1[which(tolower(miR.nameList1)%in%miRCancerOctober2017$mirId==TRUE)] #259


#######################################################################################################
#######################################################################################################
# **PPI edges
# raw.PathwayCommons,"Gene1_symbol","InteractionType","Gene2_symbol"
temp = raw.PathwayCommons[raw.PathwayCommons$Gene1_symbol%in%gene.nameList1,]
temp1 = temp[temp$Gene2_symbol%in%gene.nameList1,]

PC.edges.temp = unique(temp1[c("Gene1_symbol", "Gene2_symbol")]) # 229435 edges
# -------------------------------------  
# removing.duplicate.genePair = function(PC.edges){
#   delRows = NULL # the rows to be removed
#   for(i in 1:nrow(PC.edges)){
#     #  j = which(tab$V1 == tab$V2[i] & tab$V2 == tab$V1[i])
#     j = which(PC.edges$Gene1_symbol == PC.edges$Gene2_symbol[i] & PC.edges$Gene2_symbol == PC.edges$Gene1_symbol[i])
#     j = j [j > i]
#     if (length(j) > 0){
#       delRows = c(delRows, j)
#     }
#   }
#   PC.edges = PC.edges[-delRows,]
# }
# -------------------------------------
# PC.edges = removing.duplicate.genePair(PC.edges.temp) # 214129
PC.edges.temp2 <- data.frame(t(apply(PC.edges.temp,1,sort)), stringsAsFactors=FALSE)
PC.edges <- PC.edges.temp2[!duplicated(PC.edges.temp2),] # 214129
colnames(PC.edges) <- c("Gene1_symbol", "Gene2_symbol")

#######################################################################################################
#######################################################################################################
# ***Hypergeometric test
# define a data.frmae to save the result
modules.geneEdge_cancerMark = data.frame(
  num.cancermiRs = rep(0,ModuleNum), 
  num.cancermiRs.p = rep(0,ModuleNum), 
  num.cancermiRs.FoldChange=rep(0,ModuleNum),
  
  num.cancergenes = rep(0,ModuleNum),
  num.cancergenes.p = rep(0,ModuleNum), 
  num.cancergenes.FoldChange=rep(0,ModuleNum),
  
  num.edges = rep(0,ModuleNum),
  num.edges.p = rep(0,ModuleNum), 
  num.edges.FoldChange=rep(0,ModuleNum))
modules.cancermiR = modules.cancergene = modules.geneEdge = list()

for(i in 1:ModuleNum){
  miRs = miR.nameList1[modules[[i]]$miRs]
  genes = gene.nameList1[modules[[i]]$genes]
  miRs.cancer = miRs[miRs%in%cancermiRList]
  genes.cancer = genes[genes%in%cancerGeneList]
  
  edges = PC.edges[(PC.edges$Gene1_symbol%in%genes)&(PC.edges$Gene2_symbol%in%genes),]
  
  modules.cancermiR[[i]] = miRs.cancer
  modules.cancergene[[i]] = genes.cancer
  modules.geneEdge[[i]] = edges
  
  modules.geneEdge_cancerMark$num.cancermiRs[i] = length(miRs.cancer)
  modules.geneEdge_cancerMark$num.cancermiRs.FoldChange[i] = (length(miRs.cancer)/length(miRs))/(length(cancermiRList)/length(miR.nameList1))
  
  modules.geneEdge_cancerMark$num.cancergenes[i] = length(genes.cancer)
  modules.geneEdge_cancerMark$num.cancergenes.FoldChange[i] = (length(genes.cancer)/length(genes))/(length(cancerGeneList)/length(gene.nameList1))
  
  modules.geneEdge_cancerMark$num.edges[i] = dim(edges)[1]
  modules.geneEdge_cancerMark$num.edges.FoldChange[i] = (dim(edges)[1]/choose(length(genes),2))/(dim(PC.edges)[1]/choose(length(gene.nameList1),2))
  
  # cancer miRs
  modules.geneEdge_cancerMark$num.cancermiRs.p[i] = 
    phyper(length(miRs.cancer)-1,
           length(cancermiRList), 
           length(miR.nameList1)-length(cancermiRList), 
           length(miRs),lower.tail = F)
  # cancer genes
  # phyper(CG.Num1-1,all.CG.num, all.gene.num-all.CG.num, length(which(V1[,1]!=0)),lower.tail = F)
  modules.geneEdge_cancerMark$num.cancergenes.p[i] = 
    phyper(length(genes.cancer)-1,
           length(cancerGeneList), 
           length(gene.nameList1)-length(cancerGeneList), 
           length(genes),lower.tail = F)
  
  # gene edges
  modules.geneEdge_cancerMark$num.edges.p[i] = 
    phyper(dim(edges)[1]-1,
           dim(PC.edges)[1], 
           length(gene.nameList1)*length(gene.nameList1)/2-dim(PC.edges)[1], 
           length(genes)*length(genes)/2,lower.tail = F)
}
ResTable = modules.geneEdge_cancerMark
row.names(ResTable) = paste("Module",1:nrow(ResTable))

# save excel
fileName1 = "Result of module analysis"
wb <- createWorkbook()
addWorksheet(wb, fileName1)
writeData(wb, rowNames = T, sheet = fileName1, ResTable)
saveWorkbook(wb,"Res_Fun5_Table/Table_geneEdge_CancerGene_enrichment_analysis.xlsx",overwrite=F)

save(modules.geneEdge_cancerMark, modules.cancermiR, modules.cancergene,
     modules.geneEdge, file = "Fun5_Result.RData")


tscca_ResList = list(ResTable   = modules.geneEdge_cancerMark, 
                        cancermiR  = modules.cancermiR, 
                        cancergene = modules.cancergene, 
                        geneEdge   = modules.geneEdge)
save(tscca_ResList, file = "Fun5_tscca_ResList.RData")

#######################################################################################################
#######################################################################################################
# ******module overlapping analysis
# modules.venn = list()
# for(i in 1:20){
#   modules.venn[[i]] = list(
#     miRgenelist = c(miR.nameList1[which(out$U[,i]!=0)],gene.nameList1[which(out$V[,i]!=0)]),
#     miRlist = c(miR.nameList1[which(out$U[,i]!=0)]),
#     genelist = c(gene.nameList1[which(out$V[,i]!=0)]))
# }
# 
# overlap.mat = matrix(0,20,20)
# for(i in 1:20){
#   for(j in 1:20){
#     overlap.mat[i,j] = length(intersect(modules.venn[[i]]$miRgenelist,modules.venn[[j]]$miRgenelist))/
#       length(union(modules.venn[[i]]$miRgenelist,modules.venn[[i]]$miRgenelist))
#   }
# }
# 
# diag(overlap.mat) = 0
# write.csv(overlap.mat,file="Result.FigTable/Fun7_overlap.mat.csv")
#######################################################################################################
#######################################################################################################