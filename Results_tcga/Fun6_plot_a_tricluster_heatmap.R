# ------------------------------------------------------------------------------
# 画热图的R包
library('ComplexHeatmap')
# 调整热图颜色的包
library("circlize")
library("ggplot2")

# ------------------------------------------------------------------------------
load("../My_Data/1_tcga_data/TCGA_Tensor_Data.RData")

# 获取33种癌症的名称
dims = dim(DiffExp_Tensor)
cancerNames = DiffExp_Tensor.attribute$cancerNames
cancerNames[7] = "COAD2" # 由于这个名称太长，设置小一点替代

# ------------------------------------------------------------------------------
# function 1: a heatmap for a tricluster
tricluster_heatmap = function(module,rMdoule,DiffExp_Tensor,cancerNames,fileName){
  # get CorMat
  ModuleRes = getJointCorMat(module,DiffExp_Tensor,cancerNames)
  rModuleRes = getJointCorMat(rMdoule,DiffExp_Tensor,cancerNames)
  
  # 得到考虑模块和随机模块组合的矩阵
  Mat = rbind(ModuleRes$corMat, rModuleRes$corMat)
  # mean(abs(ModuleRes$corMat))
  
  #split.col = factor(ModuleRes$Split1)
  split.col = factor(ModuleRes$Split2, levels=ModuleRes$factorLevels)
  split.row = factor(c(rep("module genes",100),rep("random genes",100)))
  
  # 画热图的R包
  ht = Heatmap(Mat, name=paste("Mat","\nCor",sep=""), 
               cluster_rows=F, cluster_columns=F,
               #根据离散的变量将热图划分成几块模式
               row_split=split.row,column_split=split.col,
               #调节热图x和y轴的标签的字体大小
               row_title_gp=gpar(font=1,fontsize=11),
               column_title_gp=gpar(font=1,fontsize=11))
  
  #png(paste("Result_Fun4_heatmap/", method,"_module_heatmap.png",sep=""), width=1600, height=400, res = 140)
  png(fileName, width=1600, height=400, res = 140)
  print(ht)
  dev.off()
  return(ht)
}

# ------------------------------------------------------------------------------
# function 2
getJointCorMat = function(module, DiffExp_Tensor, cancerNames){
  cancers=module$cancers
  corMat = Split1 = Split2 = Split3=c()
  for(k in 1:length(module$cancers)){
    corMat = cbind(corMat, t(DiffExp_Tensor[module$miRs,module$genes,cancers[k]]))
    Split2 = c(Split2,rep(cancerNames[cancers[k]],10))
  }
  factorLevels = cancerNames[cancers]
  return(list(corMat=corMat,Split2=Split2,factorLevels=factorLevels))
}

# ------------------------------------------------------------------------------
# # function 3
get_TSCCA_heatmap = function(out, DiffExp_Tensor, cancerNames){
  # Input data
  tscca_heatmap_list = list()
  # get a random module
  set.seed(0)
  rMdoule = list(miRs=sample(1:dims[1], 10, replace=F), 
                 genes=sample(1:dims[2], 100, replace=F), 
                 cancers=sample(1:dims[3], 20, replace=F))
  # For each module and get a heat map
  for(i in 1:50){
    print(i)
    miRs    = which(out$U[,i]!=0)
    genes   = which(out$V[,i]!=0)
    cancers = which(out$W[,i]!=0)
    module  = list(miRs=miRs, genes=genes,cancers=cancers)
    
    fileName = paste("Res_Fun2_fig/TSCCA_Module_", i,"_module_heatmap.png",sep="")
    module_fig = tricluster_heatmap(module, rMdoule, DiffExp_Tensor, cancerNames, fileName=fileName)
    tscca_heatmap_list[[i]] = module_fig
  }
  return(tscca_heatmap_list)
}

# ------------------------------------------------------------------------------
tscca_heatmap_list = get_TSCCA_heatmap(tscca.out, DiffExp_Tensor, cancerNames)




