source('Functions/TSCCA.R') # TSCCA algorithm 
source('Functions/SCCA.R')
source('Functions/SA_modularity.R')
source('Functions/SA_modularity_TSCCAInit.R')
source('Functions/Evaluation.R')
# ------------------------------------------------------------------------------
library("irlba")
library("ggplot2")
library("cowplot")
library("ComplexHeatmap")
require("reshape2")
library("circlize")
# ------------------------------------------------------------------------------
get_scca_cluster = function(scca.out){
  est_cluster = list()
  est_cluster[[1]] = list(u=which(scca.out$U[,1]!=0), v=which(scca.out$V[,1]!=0),w=c(1,2,3,4))
  est_cluster[[2]] = list(u=which(scca.out$U[,2]!=0), v=which(scca.out$V[,2]!=0),w=c(1,2,3,4))
  est_cluster[[3]] = list(u=which(scca.out$U[,3]!=0), v=which(scca.out$V[,3]!=0),w=c(1,2,3,4))
  return(est_cluster)
}

get_tscca_cluster = function(tscca.out){
  est_cluster = list()
  est_cluster[[1]] = list(u=which(tscca.out$U[,1]!=0), v=which(tscca.out$V[,1]!=0), w=which(tscca.out$W[,1]!=0))
  est_cluster[[2]] = list(u=which(tscca.out$U[,2]!=0), v=which(tscca.out$V[,2]!=0), w=which(tscca.out$W[,2]!=0))
  est_cluster[[3]] = list(u=which(tscca.out$U[,3]!=0), v=which(tscca.out$V[,3]!=0), w=which(tscca.out$W[,3]!=0))
  return(est_cluster) 
}

get_SA_cluster = function(out.SA){
  est_cluster = list()
  est_cluster[[1]] = list(u=out.SA[[1]]$u, v=out.SA[[1]]$v, w=out.SA[[1]]$w)
  est_cluster[[2]] = list(u=out.SA[[2]]$u, v=out.SA[[2]]$v, w=out.SA[[2]]$w)
  est_cluster[[3]] = list(u=out.SA[[3]]$u, v=out.SA[[3]]$v, w=out.SA[[3]]$w)
  est_cluster
}
# ------------------------------------------------------------------------------
# Generate simulation data
n_repeat_times = 50
res_cs_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
res_ce_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
res_relevance_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
res_recovery_dat = matrix(0, nrow = n_repeat_times, ncol = 7)

algorithm_output_res = list()
for(kk in 1:n_repeat_times){
  print(kk)
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
  true_cluster = list()
  true_cluster[[1]] = list(u=  1:100, v=1:10)
  true_cluster[[2]] = list(u=101:200, v=11:20)
  true_cluster[[3]] = list(u=201:300, v=21:30)
  ##############################################################################
  ##############################################################################
  # sparse level
  ku = 100; # extract top 100 genes
  kv = 10;  # extract top 10  miRNAs
  kw = 2;   # extract top 3   cancers 
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
  out.scca  = SCCA(corMatSum, ku, kv, J=3)
  # Test TSCCA in tensor data
  tscca.out = TSCCA(Tensor, ku, kv, kw, J=3)
  # Test SA_modulaity in tensor data
  Res       = SA_modularity_TSCCAInit(tscca.out, Tensor, ku, kv, kw, J=3, SA_niter=8000)
  out.SA    = Res$SA_res
  
  # -----------------------------------------------------------------
  algs_res = list(alg1_clusters = get_scca_cluster(out.cancer1),
                  alg2_clusters = get_scca_cluster(out.cancer2),
                  alg3_clusters = get_scca_cluster(out.cancer3),
                  alg4_clusters = get_scca_cluster(out.cancer4),
                  alg5_clusters = get_scca_cluster(out.scca),
                  alg6_clusters = get_SA_cluster(out.SA),
                  alg7_clusters = get_tscca_cluster(tscca.out))
  algorithm_output_res[[kk]] = algs_res
  # -----------------------------------------------------------------
  res_cs_dat[kk,1] = get_consensus_score(true_cluster,algs_res[[1]])
  res_cs_dat[kk,2] = get_consensus_score(true_cluster,algs_res[[2]])
  res_cs_dat[kk,3] = get_consensus_score(true_cluster,algs_res[[3]])
  res_cs_dat[kk,4] = get_consensus_score(true_cluster,algs_res[[4]])
  res_cs_dat[kk,5] = get_consensus_score(true_cluster,algs_res[[5]])
  res_cs_dat[kk,6] = get_consensus_score(true_cluster,algs_res[[6]])
  res_cs_dat[kk,7] = get_consensus_score(true_cluster,algs_res[[7]])
  # -----------------------------------------------------------------
  res_ce_dat[kk,1] = get_clustering_error_score(true_cluster,algs_res[[1]])
  res_ce_dat[kk,2] = get_clustering_error_score(true_cluster,algs_res[[2]])
  res_ce_dat[kk,3] = get_clustering_error_score(true_cluster,algs_res[[3]])
  res_ce_dat[kk,4] = get_clustering_error_score(true_cluster,algs_res[[4]])
  res_ce_dat[kk,5] = get_clustering_error_score(true_cluster,algs_res[[5]])
  res_ce_dat[kk,6] = get_clustering_error_score(true_cluster,algs_res[[6]])
  res_ce_dat[kk,7] = get_clustering_error_score(true_cluster,algs_res[[7]])
  # -----------------------------------------------------------------
  res_recovery_dat[kk,1] = get_recovery_score(true_cluster,algs_res[[1]])
  res_recovery_dat[kk,2] = get_recovery_score(true_cluster,algs_res[[2]])
  res_recovery_dat[kk,3] = get_recovery_score(true_cluster,algs_res[[3]])
  res_recovery_dat[kk,4] = get_recovery_score(true_cluster,algs_res[[4]])
  res_recovery_dat[kk,5] = get_recovery_score(true_cluster,algs_res[[5]])
  res_recovery_dat[kk,6] = get_recovery_score(true_cluster,algs_res[[6]])
  res_recovery_dat[kk,7] = get_recovery_score(true_cluster,algs_res[[7]])
  # -----------------------------------------------------------------
  res_relevance_dat[kk,1] = get_relevance_score(true_cluster,algs_res[[1]])
  res_relevance_dat[kk,2] = get_relevance_score(true_cluster,algs_res[[2]])
  res_relevance_dat[kk,3] = get_relevance_score(true_cluster,algs_res[[3]])
  res_relevance_dat[kk,4] = get_relevance_score(true_cluster,algs_res[[4]])
  res_relevance_dat[kk,5] = get_relevance_score(true_cluster,algs_res[[5]])
  res_relevance_dat[kk,6] = get_relevance_score(true_cluster,algs_res[[6]])
  res_relevance_dat[kk,7] = get_relevance_score(true_cluster,algs_res[[7]])
}
# ------------------------------------------------------------------------------
df_names = c("cancer1","cancer2","cancer3","cancer4","Sum","Modularity_SA","TSCCA")
df.cs = data.frame(method = factor(x=df_names, levels=df_names),
                   avg=apply(res_cs_dat,2,mean), sd =apply(res_cs_dat,2,sd))

df.ce = data.frame(method = factor(x=df_names,levels=df_names),
                   avg=apply(res_ce_dat,2,mean),sd =apply(res_ce_dat,2,sd))

df.recovery = data.frame(method = factor(x=df_names,levels=df_names),
                         avg=apply(res_recovery_dat,2,mean), sd =apply(res_recovery_dat,2,sd))

df.relevance = data.frame(method = factor(x=df_names,levels=df_names),
                          avg=apply(res_relevance_dat,2,mean), sd =apply(res_relevance_dat,2,sd))
# ------------------------------------------------------------------------------
# td = Sys.time()
# save.image("Result_Fun3.RData")
# ------------------------------------------------------------------------------
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
fig = ggplot(df.recovery, aes(x=method, y=avg)) + theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  geom_bar(stat="identity", color="black", fill="#329637", size=0.62, width=0.6) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2, size=0.62) +
  theme(axis.text = element_text(colour = "black", size = rel(1.5))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "black"))+
  theme(axis.line  = element_line(size = rel(1.5), colour = "black")) +
  theme(axis.ticks = element_line(size = rel(1.5), colour = "black")) +
  theme(axis.title.y = element_text(size = rel(1.5))) +
  labs(title=NULL, x =NULL, y = "recovery")
ggsave("Fig1_result_recovery.png", width = 5, height = 4.5, units = "in")

fig = ggplot(df.ce, aes(x=method, y=avg)) + theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  geom_bar(stat="identity", color="black", fill="#F29700", size=0.62, width=0.6) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2, size=0.62) +
  theme(axis.text = element_text(colour = "black", size = rel(1.5))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "black"))+
  theme(axis.line  = element_line(size = rel(1.5), colour = "black")) +
  theme(axis.ticks = element_line(size = rel(1.5), colour = "black")) +
  theme(axis.title.y = element_text(size = rel(1.5))) +
  labs(title=NULL, x =NULL, y = "CE")
ggsave("Fig2_result_CE_score.png", width = 5, height = 4.5, units = "in")
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
get_1D_comparison_result = function(true_biclusters, algorithm_output_res, dim=1){
  # Input
  n_repeat_times = 50
  res_ce_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
  res_recovery_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
  
  for(kk in 1:n_repeat_times){
    # print(kk)
    algs_res =  algorithm_output_res[[kk]]
    # --------------------------------------------------------------------------
    res_ce_dat[kk,1] = get_1D_CE(true_biclusters,algs_res[[1]], dim)
    res_ce_dat[kk,2] = get_1D_CE(true_biclusters,algs_res[[2]], dim)
    res_ce_dat[kk,3] = get_1D_CE(true_biclusters,algs_res[[3]], dim)
    res_ce_dat[kk,4] = get_1D_CE(true_biclusters,algs_res[[4]], dim)
    res_ce_dat[kk,5] = get_1D_CE(true_biclusters,algs_res[[5]], dim)
    res_ce_dat[kk,6] = get_1D_CE(true_biclusters,algs_res[[6]], dim)
    res_ce_dat[kk,7] = get_1D_CE(true_biclusters,algs_res[[7]], dim)
    # --------------------------------------------------------------------------
    res_recovery_dat[kk,1] = get_1D_recovery(true_biclusters,algs_res[[1]], dim)
    res_recovery_dat[kk,2] = get_1D_recovery(true_biclusters,algs_res[[2]], dim)
    res_recovery_dat[kk,3] = get_1D_recovery(true_biclusters,algs_res[[3]], dim)
    res_recovery_dat[kk,4] = get_1D_recovery(true_biclusters,algs_res[[4]], dim)
    res_recovery_dat[kk,5] = get_1D_recovery(true_biclusters,algs_res[[5]], dim)
    res_recovery_dat[kk,6] = get_1D_recovery(true_biclusters,algs_res[[6]], dim)
    res_recovery_dat[kk,7] = get_1D_recovery(true_biclusters,algs_res[[7]], dim)
    # --------------------------------------------------------------------------
  }
  resMat = data.frame(ce_avg = apply(res_ce_dat,2,mean),
                      ce_sd  = apply(res_ce_dat,2,sd),
                      recovery_avg = apply(res_recovery_dat,2,mean),
                      recovery_sd  = apply(res_recovery_dat,2,sd))
  
  row.names(resMat) = c("cancer1","cancer2","cancer3","cancer4","Sum","Modularity_SA","TSCCA")
  return(resMat)
}
# ------------------------------------------------------------------------------
get_2D_comparison_result = function(true_cluster,algorithm_output_res){
  # Input
  n_repeat_times = 50
  res_ce_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
  res_recovery_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
  
  for(kk in 1:n_repeat_times){
    # print(kk)
    algs_res =  algorithm_output_res[[kk]]
    # ----------------------------------------------------------------------------
    res_ce_dat[kk,1] = get_clustering_error_score(true_cluster,algs_res[[1]])
    res_ce_dat[kk,2] = get_clustering_error_score(true_cluster,algs_res[[2]])
    res_ce_dat[kk,3] = get_clustering_error_score(true_cluster,algs_res[[3]])
    res_ce_dat[kk,4] = get_clustering_error_score(true_cluster,algs_res[[4]])
    res_ce_dat[kk,5] = get_clustering_error_score(true_cluster,algs_res[[5]])
    res_ce_dat[kk,6] = get_clustering_error_score(true_cluster,algs_res[[6]])
    res_ce_dat[kk,7] = get_clustering_error_score(true_cluster,algs_res[[7]])
    # ----------------------------------------------------------------------------
    res_recovery_dat[kk,1] = get_recovery_score(true_cluster,algs_res[[1]])
    res_recovery_dat[kk,2] = get_recovery_score(true_cluster,algs_res[[2]])
    res_recovery_dat[kk,3] = get_recovery_score(true_cluster,algs_res[[3]])
    res_recovery_dat[kk,4] = get_recovery_score(true_cluster,algs_res[[4]])
    res_recovery_dat[kk,5] = get_recovery_score(true_cluster,algs_res[[5]])
    res_recovery_dat[kk,6] = get_recovery_score(true_cluster,algs_res[[6]])
    res_recovery_dat[kk,7] = get_recovery_score(true_cluster,algs_res[[7]])
    # ----------------------------------------------------------------------------
  }
  resMat = data.frame(ce_avg = apply(res_ce_dat,2,mean),
                      ce_sd  = apply(res_ce_dat,2,sd),
                      recovery_avg = apply(res_recovery_dat,2,mean),
                      recovery_sd  = apply(res_recovery_dat,2,sd))
  
  row.names(resMat) = c("cancer1","cancer2","cancer3","cancer4","Sum","Modularity_SA","TSCCA")
  
  return(resMat)
}
# ------------------------------------------------------------------------------
get_3D_comparison_result = function(true_biclusters,algorithm_output_res){
  # Input
  n_repeat_times = 50
  res_ce_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
  res_recovery_dat = matrix(0, nrow = n_repeat_times, ncol = 7)
  
  for(kk in 1:n_repeat_times){
    # print(kk)
    algs_res =  algorithm_output_res[[kk]]
    # --------------------------------------------------------------------------
    res_ce_dat[kk,1] = get_3D_CE(true_biclusters,algs_res[[1]])
    res_ce_dat[kk,2] = get_3D_CE(true_biclusters,algs_res[[2]])
    res_ce_dat[kk,3] = get_3D_CE(true_biclusters,algs_res[[3]])
    res_ce_dat[kk,4] = get_3D_CE(true_biclusters,algs_res[[4]])
    res_ce_dat[kk,5] = get_3D_CE(true_biclusters,algs_res[[5]])
    res_ce_dat[kk,6] = get_3D_CE(true_biclusters,algs_res[[6]])
    res_ce_dat[kk,7] = get_3D_CE(true_biclusters,algs_res[[7]])
    # --------------------------------------------------------------------------
    res_recovery_dat[kk,1] = get_3D_recovery(true_biclusters,algs_res[[1]])
    res_recovery_dat[kk,2] = get_3D_recovery(true_biclusters,algs_res[[2]])
    res_recovery_dat[kk,3] = get_3D_recovery(true_biclusters,algs_res[[3]])
    res_recovery_dat[kk,4] = get_3D_recovery(true_biclusters,algs_res[[4]])
    res_recovery_dat[kk,5] = get_3D_recovery(true_biclusters,algs_res[[5]])
    res_recovery_dat[kk,6] = get_3D_recovery(true_biclusters,algs_res[[6]])
    res_recovery_dat[kk,7] = get_3D_recovery(true_biclusters,algs_res[[7]])
    # --------------------------------------------------------------------------
  }
  resMat = data.frame(ce_avg = apply(res_ce_dat,2,mean),
                      ce_sd  = apply(res_ce_dat,2,sd),
                      recovery_avg = apply(res_recovery_dat,2,mean),
                      recovery_sd  = apply(res_recovery_dat,2,sd))
  
  row.names(resMat) = c("cancer1","cancer2","cancer3","cancer4","Sum","Modularity_SA","TSCCA")
  
  return(resMat)
}

# ------------------------------------------------------------------------------
# get table result
source('Functions/Evaluation.R')
true_biclusters = list()
true_biclusters[[1]] = list(u=  1:100, v=1:10,  w=c(1,2))
true_biclusters[[2]] = list(u=101:200, v=11:20, w=c(1,3))
true_biclusters[[3]] = list(u=201:300, v=21:30, w=c(2,3))
# ------------------------------------------------------------------------------
res3D = get_3D_comparison_result(true_biclusters,algorithm_output_res)
colnames(res3D) = paste("3D", colnames(res3D)) 

res2D = get_2D_comparison_result(true_biclusters,algorithm_output_res)
colnames(res2D) = paste("2D", colnames(res2D)) 

res1D_gene = get_1D_comparison_result(true_biclusters,algorithm_output_res, dim=1)
colnames(res1D_gene) = paste("gene", colnames(res1D_gene)) 

res1D_miR  = get_1D_comparison_result(true_biclusters,algorithm_output_res, dim=2)
colnames(res1D_miR) = paste("miR", colnames(res1D_miR)) 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
get_table = function(res3D,res2D,res1D_gene,res1D_miR){
  raw_table = cbind(cbind(cbind(res1D_gene,res1D_miR),res2D),res3D)
  table = format(round(raw_table, 3), nsmall = 3)
  
  str_table = matrix(NA,7,8)
  rownames(str_table) = rownames(table)
  colnames(str_table) = c("1D-CE-gene", "1D-recovery-gene", 
                          "1D-CE-miR",  "1D-recovery-miR", 
                          "2D-CE",      "2D-recovery",
                          "3D-CE",      "3D-recovery")
  
  for(i in 1:nrow(str_table)){
    for(j in 1:ncol(str_table)){
      k = 2*(j-1) + 1
      str_table[i,j] = paste(table[i,k], "+", table[i,k+1],sep="")
    }
  }
  # str_table = str_table[c("TSCCA","Modularity_SA","Sum","cancer4","cancer3","cancer2","cancer1"), ]
  return(str_table)
}

str_table = get_table(res3D,res2D,res1D_gene,res1D_miR)
write.csv(str_table, file = "Res_Table.csv", row.names = T, quote = F)

td = Sys.time()
save.image("simulation_result.RData")
# ------------------------------------------------------------------------------
