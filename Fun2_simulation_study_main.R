############################################################################
# source('Functions/Fun1_setDefaultPath.R')
source('Functions/TSCCA.R') # TSCCA algorithm 
source('Functions/SCCA.R')
source('Functions/SA_modularity.R')
source('Functions/SA_modularity_TSCCAInit.R')
source('Functions/Evaluation.R')
############################################################################
library("irlba")
library("ggplot2")
library("cowplot")
library("ComplexHeatmap")
require("reshape2")
library("circlize")
############################################################################
############################################################################
get_tscca_bicluster = function(tscca.out){
  est_bicluster = list()
  est_bicluster[[1]] = list(u=which(tscca.out$U[,1]!=0), v=which(tscca.out$V[,1]!=0))
  est_bicluster[[2]] = list(u=which(tscca.out$U[,2]!=0), v=which(tscca.out$V[,2]!=0))
  est_bicluster[[3]] = list(u=which(tscca.out$U[,3]!=0), v=which(tscca.out$V[,3]!=0))
  est_bicluster
}

get_SA_bicluster = function(out.SA){
  est_bicluster = list()
  est_bicluster[[1]] = list(u=out.SA[[1]]$u, v=out.SA[[1]]$v)
  est_bicluster[[2]] = list(u=out.SA[[2]]$u, v=out.SA[[2]]$v)
  est_bicluster[[3]] = list(u=out.SA[[3]]$u, v=out.SA[[3]]$v)
  est_bicluster
}
############################################################################
############################################################################
n_repeat_times = 50
res_cs_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
res_ce_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
res_relevance_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
res_recovery_dat = matrix(0, nrow = n_repeat_times, ncol = 7)

algorithm_output_res = list()
for(kk in 1:n_repeat_times){
  print(kk)
  ##########################################################################
  # 1 generate data
  set.seed(kk) # 1
  p = 300; q = 30; cancerTypes = 4
  Tensor = array(0, c(p,q,cancerTypes))
  sd0 = 0.2
  Tensor = array(rnorm(p*q*cancerTypes, mean = 0, sd = sd0), c(p,q,cancerTypes))
  # Homogeneous block-structures 
  Tensor[  1:100,  1:10, 1] <- matrix(rnorm(100*10, mean =  0.5, sd = sd0), nrow=100) 
  Tensor[101:200, 11:20, 1] <- matrix(rnorm(100*10, mean = -0.5, sd = sd0), nrow=100) 
  
  Tensor[  1:100,  1:10, 2] <- matrix(rnorm(100*10, mean = -0.5, sd = sd0), nrow=100) 
  Tensor[201:300, 21:30, 2] <- matrix(rnorm(100*10, mean =  0.5, sd = sd0), nrow=100) 
  
  Tensor[101:200, 11:20, 3] <- matrix(rnorm(100*10, mean =  0.5, sd = sd0), nrow=100) 
  Tensor[201:300, 21:30, 3] <- matrix(rnorm(100*10, mean = -0.5, sd = sd0), nrow=100) 
  ##############################################################################
  ##############################################################################
  # 每一个pattern由两部分组成
  true_bicluster = list()
  true_bicluster[[1]] = list(u=  1:100, v=1:10)
  true_bicluster[[2]] = list(u=101:200, v=11:20)
  true_bicluster[[3]] = list(u=201:300, v=21:30)
  ##############################################################################
  ##############################################################################
  # sparse level
  ku = 100; # extract top 100 genes
  kv = 10;  # extract top 10  miRNAs
  kw = 3;   # extract top 3   cancers 
  # Runing SCCA
  # Test SCCA in data 1
  out.cancer1 = SCCA(Tensor[,,1], ku, kv, J=3)
  # Test SCCA in data 2
  out.cancer2 = SCCA(Tensor[,,2], ku, kv, J=3)
  # Test SCCA in data 3
  out.cancer3 = SCCA(Tensor[,,3], ku, kv, J=3)
  # Test SCCA in data 4
  out.cancer4 = SCCA(Tensor[,,4], ku, kv, J=3)
  # Test SCCA in sum data
  corMatSum = Tensor[,,1] + Tensor[,,2] + Tensor[,,3] + Tensor[,,4]
  out.scca = SCCA(corMatSum, ku, kv, J=3)
  # Test TSCCA in tensor data
  tscca.out = TSCCA(Tensor, ku, kv, kw, J=3)
  # Test SA_modulaity in tensor data
  Res = SA_modularity_TSCCAInit(tscca.out, Tensor, ku, kv, kw, J=3, SA_niter=8000)
  out.SA = Res$SA_res
  # Example
  # out.SA.rank1 = SA.modularity(Tensor, ku, kv, kw, niter=4000)
  # plot(out.SA.rank1$obj)
  # -----------------------------------------------------------------
  algs_res = list(alg1_biclusters = get_tscca_bicluster(out.cancer1),
                  alg2_biclusters = get_tscca_bicluster(out.cancer2),
                  alg3_biclusters = get_tscca_bicluster(out.cancer3),
                  alg4_biclusters = get_tscca_bicluster(out.cancer4),
                  alg5_biclusters = get_tscca_bicluster(out.scca),
                  alg6_biclusters = get_SA_bicluster(out.SA),
                  alg7_biclusters = get_tscca_bicluster(tscca.out))
  algorithm_output_res[[kk]] = algs_res
  # -----------------------------------------------------------------
  res_cs_dat[kk,1] = get_consensus_score(true_bicluster,algs_res[[1]])
  res_cs_dat[kk,2] = get_consensus_score(true_bicluster,algs_res[[2]])
  res_cs_dat[kk,3] = get_consensus_score(true_bicluster,algs_res[[3]])
  res_cs_dat[kk,4] = get_consensus_score(true_bicluster,algs_res[[4]])
  res_cs_dat[kk,5] = get_consensus_score(true_bicluster,algs_res[[5]])
  res_cs_dat[kk,6] = get_consensus_score(true_bicluster,algs_res[[6]])
  res_cs_dat[kk,7] = get_consensus_score(true_bicluster,algs_res[[7]])
  # -----------------------------------------------------------------
  res_ce_dat[kk,1] = get_clustering_error_score(true_bicluster,algs_res[[1]])
  res_ce_dat[kk,2] = get_clustering_error_score(true_bicluster,algs_res[[2]])
  res_ce_dat[kk,3] = get_clustering_error_score(true_bicluster,algs_res[[3]])
  res_ce_dat[kk,4] = get_clustering_error_score(true_bicluster,algs_res[[4]])
  res_ce_dat[kk,5] = get_clustering_error_score(true_bicluster,algs_res[[5]])
  res_ce_dat[kk,6] = get_clustering_error_score(true_bicluster,algs_res[[6]])
  res_ce_dat[kk,7] = get_clustering_error_score(true_bicluster,algs_res[[7]])
  # -----------------------------------------------------------------
  res_recovery_dat[kk,1] = get_recovery_score(true_bicluster,algs_res[[1]])
  res_recovery_dat[kk,2] = get_recovery_score(true_bicluster,algs_res[[2]])
  res_recovery_dat[kk,3] = get_recovery_score(true_bicluster,algs_res[[3]])
  res_recovery_dat[kk,4] = get_recovery_score(true_bicluster,algs_res[[4]])
  res_recovery_dat[kk,5] = get_recovery_score(true_bicluster,algs_res[[5]])
  res_recovery_dat[kk,6] = get_recovery_score(true_bicluster,algs_res[[6]])
  res_recovery_dat[kk,7] = get_recovery_score(true_bicluster,algs_res[[7]])
  # -----------------------------------------------------------------
  res_relevance_dat[kk,1] = get_relevance_score(true_bicluster,algs_res[[1]])
  res_relevance_dat[kk,2] = get_relevance_score(true_bicluster,algs_res[[2]])
  res_relevance_dat[kk,3] = get_relevance_score(true_bicluster,algs_res[[3]])
  res_relevance_dat[kk,4] = get_relevance_score(true_bicluster,algs_res[[4]])
  res_relevance_dat[kk,5] = get_relevance_score(true_bicluster,algs_res[[5]])
  res_relevance_dat[kk,6] = get_relevance_score(true_bicluster,algs_res[[6]])
  res_relevance_dat[kk,7] = get_relevance_score(true_bicluster,algs_res[[7]])
  # -----------------------------------------------------------------
  # print("-------------------------")
  # print("get_consensus_score")
  # print(res_cs_dat[kk,])
  # print("-------------------------")
  # print("get_clustering_error_score")
  # print(res_ce_dat[kk,])
  # print("-------------------------")
  # print("get_recovery_score")
  # print(res_recovery_dat[kk,])
  # print("-------------------------")
  # print("get_relevance_score")
  # print(res_relevance_dat[kk,])
}
df_names = c("cancer1","cancer2","cancer3","cancer4","Sum","Modularity_SA","TSCCA")
df.cs = data.frame(method = factor(x=df_names, levels=df_names),
                   avg=apply(res_cs_dat,2,mean), sd =apply(res_cs_dat,2,sd))

df.ce = data.frame(method = factor(x=df_names,levels=df_names),
                   avg=apply(res_ce_dat,2,mean),sd =apply(res_ce_dat,2,sd))

df.recovery = data.frame(method = factor(x=df_names,levels=df_names),
                    avg=apply(res_recovery_dat,2,mean), sd =apply(res_recovery_dat,2,sd))

df.relevance = data.frame(method = factor(x=df_names,levels=df_names),
                    avg=apply(res_relevance_dat,2,mean), sd =apply(res_relevance_dat,2,sd))

save.image("ResultData_simulationStudy_Final.RData")
############################################################################
############################################################################
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
fig = ggplot(df.recovery, aes(x=method, y=avg)) + theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  geom_bar(stat="identity", color="black", fill="#329637", size=0.62, width=0.6) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2, size=0.62) +
  theme(axis.text = element_text(colour = "black", size = rel(1.5))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "black"))+
  # X和Y轴的线设
  theme(axis.line  = element_line(size = rel(1.5), colour = "black")) +
  # 坐标轴上的小线设
  theme(axis.ticks = element_line(size = rel(1.5), colour = "black")) +
  theme(axis.title.y = element_text(size = rel(1.5))) +
  labs(title=NULL, x =NULL, y = "recovery")
#ggsave("Fig_result_recovery+.png", width = 7.5, height = 4.5, units = "in")
ggsave("Fig_result_recovery.png", width = 5, height = 4.5, units = "in")

fig = ggplot(df.ce, aes(x=method, y=avg)) + theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  geom_bar(stat="identity", color="black", fill="#F29700", size=0.62, width=0.6) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2, size=0.62) +
  theme(axis.text = element_text(colour = "black", size = rel(1.5))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "black"))+
  # X和Y轴的线设
  theme(axis.line  = element_line(size = rel(1.5), colour = "black")) +
  # 坐标轴上的小线设
  theme(axis.ticks = element_line(size = rel(1.5), colour = "black")) +
  theme(axis.title.y = element_text(size = rel(1.5))) +
  labs(title=NULL, x =NULL, y = "CE")
ggsave("Fig_result_CE_score.png", width = 5, height = 4.5, units = "in")
################################################################################
################################################################################