module_analysis = function(modules, 
                           DiffExp_Tensor, 
                           DiffExp_Tensor.attribute,
                           allOnco_May2018, 
                           miRCancerOctober2017,
                           miRTarBase.dat, 
                           raw.PathwayCommons){
  #-----------------------------------------------------------------------------
  # Input: out is the output of scca or tscca
  # Except "out", other data can be obtained from below
  # load("Fun2_scca_ACCData_Result.RData")
  # load("../My_Data/1_tcga_data/TCGA_Tensor_Data.RData")
  # load("../My_Data/geneEdge_CancerGene.RData")
  # load("../My_Data/3_miRTarBase/miRTarBaseData.RData")
  
  #-----------------------------------------------------------------------------
  ## Step1: extract 50 miRNA-gene modules
  # ModuleNum = length(out$d)
  # modules = list()
  # for(i in 1:ModuleNum){
  #   modules[[i]] = list(miRs = which(out$U[,i]!=0),genes = which(out$V[,i]!=0),SingularValue_d = out$D[i])
  # }
  
  ModuleNum = length(modules)
  
  # miRNAs
  miR.nameList1 = miR.nameList2 = rep("NA",length(DiffExp_Tensor.attribute$miRNames))
  for(i in 1:length(DiffExp_Tensor.attribute$miRNames)){
    temp = unlist(strsplit(DiffExp_Tensor.attribute$miRNames[i], "[|]"))
    miR.nameList1[i] = temp[1] # hsa-let-7a-5p....
    miR.nameList2[i] = temp[2] # MIMAT0000062....
  }
  DiffExp_Tensor.attribute$miRNames_C1 = miR.nameList1
  DiffExp_Tensor.attribute$miRNames_C2 = miR.nameList2
  # genes
  gene.nameList1 = gene.nameList2 = rep("NA",length(DiffExp_Tensor.attribute$miRNames))
  for(i in 1:length(DiffExp_Tensor.attribute$geneNames)){
    temp = unlist(strsplit(DiffExp_Tensor.attribute$geneNames[i], "[|]"))
    gene.nameList1[i] = temp[1] # A1BG....
    gene.nameList2[i] = temp[2] # 53947...
  }
  DiffExp_Tensor.attribute$geneNames_C1 = gene.nameList1
  DiffExp_Tensor.attribute$geneNames_C2 = gene.nameList2
  
  #-----------------------------------------------------------------------------
  ## Step3: get cancer-genes and miRNAs, and gene edges
  cancerGeneList = gene.nameList1[which(gene.nameList1%in%allOnco_May2018$symbol==TRUE)] # 1164
  cancermiRList = miR.nameList1[which(tolower(miR.nameList1)%in%miRCancerOctober2017$mirId==TRUE)] #259
  
  # PPI edges
  # raw.PathwayCommons,"Gene1_symbol","InteractionType","Gene2_symbol"
  temp = raw.PathwayCommons[raw.PathwayCommons$Gene1_symbol%in%gene.nameList1,]
  temp1 = temp[temp$Gene2_symbol%in%gene.nameList1,]
  
  PC.edges.temp = unique(temp1[c("Gene1_symbol", "Gene2_symbol")]) # 229435 edges
  
  # PC.edges = removing.duplicate.genePair(PC.edges.temp) # 214129
  PC.edges.temp2 <- data.frame(t(apply(PC.edges.temp,1,sort)), stringsAsFactors=FALSE)
  PC.edges <- PC.edges.temp2[!duplicated(PC.edges.temp2),] # 214129
  colnames(PC.edges) <- c("Gene1_symbol", "Gene2_symbol")
  
  
  #-----------------------------------------------------------------------------
  ## Step3: define a data.frmae to save the result
  sumTable = data.frame(
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
  
  #-----------------------------------------------------------------------------
  ## Step4: analysis for each module
  for(i in 1:ModuleNum){
    miRs = miR.nameList1[modules[[i]]$miRs]
    genes = gene.nameList1[modules[[i]]$genes]
    miRs.cancer = miRs[miRs%in%cancermiRList]
    genes.cancer = genes[genes%in%cancerGeneList]
    
    edges = PC.edges[(PC.edges$Gene1_symbol%in%genes)&(PC.edges$Gene2_symbol%in%genes),]
    
    modules.cancermiR[[i]] = miRs.cancer
    modules.cancergene[[i]] = genes.cancer
    modules.geneEdge[[i]] = edges
    
    sumTable$num.cancermiRs[i] = length(miRs.cancer)
    sumTable$num.cancermiRs.FoldChange[i] = (length(miRs.cancer)/length(miRs))/(length(cancermiRList)/length(miR.nameList1))
    
    sumTable$num.cancergenes[i] = length(genes.cancer)
    sumTable$num.cancergenes.FoldChange[i] = (length(genes.cancer)/length(genes))/(length(cancerGeneList)/length(gene.nameList1))
    
    sumTable$num.edges[i] = dim(edges)[1]
    sumTable$num.edges.FoldChange[i] = (dim(edges)[1]/choose(length(genes),2))/(dim(PC.edges)[1]/choose(length(gene.nameList1),2))
    
    ## Hypergeometric test
    # cancer miRs
    sumTable$num.cancermiRs.p[i] = 
      phyper(length(miRs.cancer)-1,
             length(cancermiRList), 
             length(miR.nameList1)-length(cancermiRList), 
             length(miRs),lower.tail = F)
    
    # cancer genes
    # phyper(CG.Num1-1,all.CG.num, all.gene.num-all.CG.num, length(which(V1[,1]!=0)),lower.tail = F)
    sumTable$num.cancergenes.p[i] = 
      phyper(length(genes.cancer)-1,
             length(cancerGeneList), 
             length(gene.nameList1)-length(cancerGeneList), 
             length(genes),lower.tail = F)
    
    # gene edges
    sumTable$num.edges.p[i] = 
      phyper(dim(edges)[1]-1,
             dim(PC.edges)[1], 
             length(gene.nameList1)*length(gene.nameList1)/2-dim(PC.edges)[1], 
             length(genes)*length(genes)/2,lower.tail = F)
  }
  
  #-----------------------------------------------------------------------------
  ## Step5: analysis for each module on miRTarBaseData 
  mg_num_edge = get_miRGene_edge(modules, miRTarBase.dat, DiffExp_Tensor.attribute)
  
  
  #-----------------------------------------------------------------------------
  ## step6: save result
  ResTable = cbind(sumTable,mg_num_edge)
  row.names(ResTable) = paste("Module",1:nrow(ResTable))
  
  
  resList = list(ResTable   = ResTable, 
                 cancermiR  = modules.cancermiR, 
                 cancergene = modules.cancergene, 
                 geneEdge   = modules.geneEdge)
  
  return(resList)
}

get_miRGene_edge = function(modules, miRTarBase.dat, DiffExp_Tensor.attribute){
  # ---------------------------------------------------------------------------------
  ## miRTarBaseData analysis
  # load("../My_Data/3_miRTarBase/miRTarBaseData.RData")
  # mg_num_edge = get_miRGene_edge(modules, miRTarBase.dat, DiffExp_Tensor.attribute)
  # ---------------------------------------------------------------------------------
  filter1 = miRTarBase.dat[miRTarBase.dat$miRNA%in%DiffExp_Tensor.attribute$miRNames_C1,]
  filter2 = filter1[filter1$TargetGene.EntrezGeneID%in%DiffExp_Tensor.attribute$geneNames_C2,]
  filter3 = unique(filter2[c("miRNA", "TargetGene.EntrezGeneID")]) # there are 54321 edges
  
  ModuleNum = length(modules)
  
  mg_edge_mat = matrix(0,ModuleNum,2)
  colnames(mg_edge_mat) = c("mg_edge","p_val")
  
  module_mg_edge = list()
  for(i in 1:ModuleNum){
    # The miRNAs of module i belong to miRTarBase table
    temp1 = filter3[filter3$miRNA%in%DiffExp_Tensor.attribute$miRNames_C1[modules[[i]]$miRs],]
    
    # The genes of the module i belongs to miRTarBase table
    temp2 = temp1[temp1$TargetGene.EntrezGeneID%in%DiffExp_Tensor.attribute$geneNames_C2[modules[[i]]$genes],]
    
    # miRNA-gene edges for module i 
    module_mg_edge[[i]] = temp2
    
    # Hypergeometric test
    q = nrow(temp2)   # the number of white balls
    m = nrow(filter3) # (m = 54321) the number of white balls in the urn
    n = 7889*523 - m  # the number of black balls in the urn
    k = 10*100        # the number of balls drawn from the urn
    
    mg_edge_mat[i,1] = q
    mg_edge_mat[i,2] = phyper(q-1, m, n, k, lower.tail = F)
  }
  return(mg_edge_mat = mg_edge_mat)
}

# ------------------------------------------------------------------------------
get_tscca_modularity = function(out, DiffExp_Tensor){
  # TSCCA 
  # out = diff_out[[1]]
  ModuleNum = length(out$d)
  modularity = c()
  for(i in 1:ModuleNum){
    modularity[i] = mean(abs(DiffExp_Tensor[which(out$U[,i]!=0), which(out$V[,i]!=0), which(out$W[,i]!=0)]))
  }
  return(modularity)
}
# ------------------------------------------------------------------------------
get_scca_modularity = function(out, cancerID, DiffExp_Tensor){
  # SCCA 
  # out = diff_out[[2]]
  ModuleNum = length(out$d)
  modularity = sc_modularity = c()
  for(i in 1:ModuleNum){
    modularity[i] = mean(abs(DiffExp_Tensor[which(out$U[,i]!=0), which(out$V[,i]!=0), ]))
    # single cancer modularity
    sc_modularity[i] = mean(abs(DiffExp_Tensor[which(out$U[,i]!=0), which(out$V[,i]!=0), cancerID]))
  }
  return(list(modularity=modularity, sc_modularity=sc_modularity))
} 
# ------------------------------------------------------------------------------