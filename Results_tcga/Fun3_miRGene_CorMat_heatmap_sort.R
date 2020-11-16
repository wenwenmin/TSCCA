# ------------------------------------------------------------------------------
# Function 1: 得到模块在不同癌症组合的miRNA-gene相关矩阵
getJointCorMat = function(module,cancerOrder,cancerNames,topNum=10){
  corMat = Split1 = Split2 = Split3=c()
  for(k in 1:topNum){
    corMat = cbind(corMat, t(DiffExp_Tensor[module$miRs,module$genes,cancerOrder[k]]))
    # 确保热图中的块排序是安装癌症权重值大小排序的
    Split1 = c(Split1,rep(paste(letters[k],cancerNames[cancerOrder[k]],sep=""),10))
    Split2 = c(Split2,rep(cancerNames[cancerOrder[k]],10))
    Split3 = c(Split3,cancerNames[cancerOrder[k]])
  }
  factorLevels = cancerNames[cancerOrder][1:topNum]
  return(list(corMat=corMat,Split1=Split1,Split2=Split2,Split3=Split3,factorLevels=factorLevels))
}

# ------------------------------------------------------------------------------
# Function 2: get 50 heat map figures
get_TSCCA_heatmap = function(out, DiffExp_Tensor, folderName){
  # Input data
  # out = TSCCA(DiffExp_Tensor, ku, kv, kw, J)
  
  # Step 1 获取TSCCA识别的50个miRNA-gene模块
  modules = list()
  for(i in 1:length(out$d)){
    modules[[i]] = list(miRs=which(out$U[,i]!=0), genes=which(out$V[,i]!=0), cancer.weight=out$W[,i])
  }
  order.modules = modules
  for(i in 1:length(out$d)){
    order.modules[[i]]$miRs=
      order.modules[[i]]$miRs[order(-abs(out$U[order.modules[[i]]$miRs,i]))]
    order.modules[[i]]$genes=
      order.modules[[i]]$genes[order(-abs(out$V[order.modules[[i]]$genes,i]))]
  }
  
  # Step 2 获取33种癌症的名称
  cancerNames = DiffExp_Tensor.attribute$cancerNames
  cancerNames[7] = "COAD2" # 由于这个名称太长，设置小一点替代
  
  # Step 3 get heat map
  tscca_heatmap_list = list()
  for(i in 1:50){
    print(i)
    cancerOrder = order(-abs(order.modules[[i]]$cancer.weight))
    ModuleRes   = getJointCorMat(order.modules[[i]],cancerOrder,cancerNames,20)
    
    # 产生一个随机模块
    dims = dim(DiffExp_Tensor)
    set.seed(i)
    miRs    = sample(1:dims[1], 10,  replace=F)
    genes   = sample(1:dims[2], 100, replace=F)
    rMdoule = list(miRs=miRs,genes=genes)
    
    # 得到随机模块的相关性矩阵
    rModuleRes = getJointCorMat(rMdoule, cancerOrder, cancerNames, 20)
    
    # 得到考虑模块和随机模块组合的矩阵
    Mat = rbind(ModuleRes$corMat, rModuleRes$corMat)
    
    split.col = factor(ModuleRes$Split2,levels=ModuleRes$factorLevels)
    split.row = factor(c(rep("module genes",100),rep("random genes",100)))
    
    # 画热图的R包
    ht = Heatmap(Mat, name=paste("M",i,"\nCor",sep=""), cluster_rows=F, cluster_columns=F,
                 #根据离散的变量将热图划分成几块模式
                 row_split=split.row,column_split=split.col,
                 #调节热图x和y轴的标签的字体大小
                 row_title_gp=gpar(font=1,fontsize=11),column_title_gp=gpar(font=1,fontsize=11))
    tscca_heatmap_list[[i]] = ht
    
    png(paste(folderName,i,"_Heatmap.png",sep=""), width=1600,height=400,res = 140)
    print(ht)
    dev.off()
  }
  return(tscca_heatmap_list)
}

# ------------------------------------------------------------------------------
load("../My_Data/1_tcga_data/TCGA_Tensor_Data.RData")
load("Fun1_Result_ku10_kv100_kw20_J50.RData")

# 画热图的R包
library('ComplexHeatmap')
# 调整热图颜色的包
library("circlize")

folderName = "Res_Fun3_fig/TSCCA_Module_"
tscca_heatmap_list = get_TSCCA_heatmap(tscca.out, DiffExp_Tensor, folderName = folderName)



