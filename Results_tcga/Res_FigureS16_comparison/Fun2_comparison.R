source('Functions/SCCA_rank1.R')
source('Functions/TSCCA_rank1.R')
source('Functions/SA_modularity.R')
source('Functions/SCP_rank1.R')
source('Functions/SCCA_rank1_Final.R')
source('Functions/SNCP_rank1.R')

ptm0 = proc.time()

load("../My_Data/1_tcga_data/TCGA_Tensor_Data.RData")
# ------------------------------------------------------------------------------
load("../My_Data/1_tcga_data/TCGA_JointCor_Data.RData")
corMatSum = JointCor
# ------------------------------------------------------------------------------

permutation_modularity_scores = rep(0,50)
for(j in 1:length(permutation_modularity_scores)){
  set.seed(j)
  permut.miRs  = sample(1:523,  10,  replace = FALSE)
  permut.genes = sample(1:7889, 100, replace = FALSE)
  permut.cancers = sample(1:33, 20,  replace = FALSE)
  permutation_modularity_scores[j] = mean(abs(DiffExp_Tensor[permut.miRs,permut.genes,permut.cancers]))
}

# parameters
ku = 10; kv = 100; kw = 20; J = 1; parameters=list(ku=ku,kv=kv,kw=kw,J=J)
### comparison in different methods
repeatNum = 20
# save result matrix
time_tab = modu_tab = matrix(0, nrow=repeatNum, ncol=7)
d5 = d6 = d7 = c()
Res_comparison = list()
for(i in 1:repeatNum){
  print(i)
  
  # method 1
  out1 = SCCA_rank1(corMatSum, ku, kv, seed0=i)
  modu_tab[i,1] = mean(abs(DiffExp_Tensor[which(out1$u!=0), which(out1$v!=0), ]))
  time_tab[i,1] = out1$time
  
  # method 2 ACC data
  out2 = SCCA_rank1(DiffExp_Tensor[,,1], ku, kv, seed0=i)
  modu_tab[i,2] = mean(abs(DiffExp_Tensor[which(out2$u!=0), which(out2$v!=0), ]))
  time_tab[i,2] = out2$time
  
  # method 3 BRCA data
  out3 = SCCA_rank1(DiffExp_Tensor[,,3], ku, kv, seed0=i)
  modu_tab[i,3] = mean(abs(DiffExp_Tensor[which(out3$u!=0), which(out3$v!=0), ]))
  time_tab[i,3] = out3$time
  
  # method 4
  out4 = SA.modularity(DiffExp_Tensor, ku, kv, kw, niter=50000, flag=1, seed0=i)
  modu_tab[i,4] = out4$opt_obj
  time_tab[i,4] = out4$time
  
  # method 5
  out5 = TSCCA_rank1(DiffExp_Tensor, ku, kv, kw, Init_flag=2, niter=20, seed0=i)
  modu_tab[i,5] = out5$modularity
  time_tab[i,5] = out5$time
  d5 = c(d5,out5$d)
  
  # method 6
  out6 = SCP_rank1(DiffExp_Tensor, ku, kv, kw, Init_flag=2, niter=20, seed0=i)
  modu_tab[i,6] = out6$modularity
  time_tab[i,6] = out6$time
  d6 = c(d6,out6$d)
  
  # method 7
  out7 = SNCP_rank1(DiffExp_Tensor, ku, kv, kw, Init_flag=2, niter=20, seed0=i)
  modu_tab[i,7] = out7$modularity
  time_tab[i,7] = out7$time
  d7 = c(d7,out7$d)
  
  Res_comparison[[i]] = list(out1=out1,out2=out2,out3=out3,out4=out4,out5=out5,out6=out6,out7=out7)
}

# colnames(modu_tab) = colnames(time_tab) = c("SCCA_sumDat", "SCCA_ACCDat","SCCA_BRCAData","SA","TSCCA","SCP","SNCP")
colnames(modu_tab) = colnames(time_tab) = c("JointData",
                                            "ACC_Data",
                                            "BRCA_Data",
                                            "Modularity_SA",
                                            "TSCCA",
                                            "SCP",
                                            "SNCP")
print(apply(modu_tab,2,mean))
print(apply(time_tab,2,mean))
print(proc.time() - ptm0)
save(modu_tab,time_tab,Res_comparison,parameters,file = "Fun2_comparison_result.RData")

################################################################################
################################################################################

library("ggplot2")
library("ggpubr")

get_boxPlot = function(modu_tab,time_tab){
  # input modu_tab and time_tab
  factorLevels = c("ACC_Data","BRCA_Data","JointData","Modularity_SA","SNCP","SCP","TSCCA")
  modu_dat <- stack(as.data.frame(modu_tab))
  time_dat <- stack(as.data.frame(time_tab))
  
  # matrix -->ggplot data_frame
  modu_dat$ind = factor(modu_dat$ind, levels=factorLevels)
  time_dat$ind = factor(time_dat$ind, levels=factorLevels)
  
  # ggplot2 box plot
  p1 = ggplot(modu_dat, aes(x=ind, y=values, color=ind)) + 
    geom_boxplot(size=1.2) + 
    theme(axis.line=element_line(size=rel(2)),
          axis.ticks=element_line(size=rel(2)),
          axis.title = element_text(size = rel(1.3)),
          axis.text = element_text(colour = "black",size=rel(1.2)),
          axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "none",
          legend.title =element_text(size=1, colour="white"), 
          legend.text=element_text(size=rel(1.3), colour="black")) +
    ylab("Modularity score") + xlab(NULL)
  ggsave("Fun6_boxplotA.png", width = 4.2, height = 5, units = "in")
  
  p2 = ggplot(time_dat, aes(x=ind, y=values, color=ind)) + 
    geom_boxplot(size=1.2) +
    theme(axis.line=element_line(size=rel(2)),
          axis.ticks=element_line(size=rel(2)),
          axis.title = element_text(size = rel(1.3)),
          axis.text = element_text(colour = "black",size=rel(1.2)),
          axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "none",
          legend.title =element_text(size=1, colour="white"), 
          legend.text=element_text(size=rel(1.3), colour="black")) +
    ylab("Time (s)") + xlab(NULL)
  
  return(list(p1=p1,p2=p2))
}
figs = get_boxPlot(modu_tab,time_tab)
ggarrange(figs$p1, figs$p2, labels = c("A", "B"), ncol = 2, nrow = 1)
ggsave("Fun2_comparison_figure.png", width = 8.7, height = 5.4, units = "in")

