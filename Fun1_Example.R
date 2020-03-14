############################################################################
# source('Functions/Fun0_setDefaultPath.R')
source('Functions/TSCCA.R') # TSCCA algorithm 
source('Functions/SCCA.R')
source('Functions/SA_modularity.R')
source('Functions/TSCCA_SA.R')
############################################################################
library("irlba")
library("ggplot2")
library("cowplot")
library("ComplexHeatmap")
require("reshape2")
library("circlize")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
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

get_jaccard = function(M1, M2){
  # M1 = true_bicluster[[1]]
  # M2 = est.biclusters[[2]]
  n1 = length(intersect(M1$u,M2$u))
  n2 = length(intersect(M1$v,M2$v))
  N = length(M1$u)*length(M1$v) + length(M2$u)*length(M2$v) - n1*n2
  jaccard_index = n1*n2/N 
}

get_consensus_score = function(true_bicluster,est_bicluster){
  # est_bicluster = get_SA_bicluster(out.SA)
  scoreMat = matrix(0,nrow=3,ncol=3)
  for(i in 1:3)
    for(j in 1:3){
      M1 = true_bicluster[[i]]
      M2 = est_bicluster[[j]]
      scoreMat[i,j] = get_jaccard(M1,M2) 
    }
  score = 0
  # 每一列取最�?
  score = score + max(scoreMat[,1]); flg1 = which.max(scoreMat[,1]); scoreMat[ ,1] = 0;
  score = score + max(scoreMat[,2]); flg2 = which.max(scoreMat[,2]); scoreMat[ ,2] = 0; 
  score = score + max(scoreMat[,3]);
  score = score/3
}
############################################################################
############################################################################
n_repeat_times = 50
res.dat = matrix(0, nrow = n_repeat_times, ncol = 7)
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
  # Runing TSCCA
  ku = 100; # set parameters of TSCCA 
  kv = 10;  # set parameters of TSCCA
  kw = 3;   # set parameters of TSCCA
  tscca.out = TSCCA(Tensor, ku, kv, kw, J=3)
  
  # Runing SCCA
  corMatSum = Tensor[,,1] + Tensor[,,2] + Tensor[,,3] + Tensor[,,4]
  out.scca = SCCA(corMatSum, ku, kv, J=3)
  
  # Runing SCCA for each cancer data
  out.cancer1 = SCCA(Tensor[,,1], ku, kv, J=3)
  out.cancer2 = SCCA(Tensor[,,2], ku, kv, J=3)
  out.cancer3 = SCCA(Tensor[,,3], ku, kv, J=3)
  out.cancer4 = SCCA(Tensor[,,4], ku, kv, J=3)
  
  # Runing SA for maximizing modulaity 
  # out.SA.rank1 = SA.modularity(Tensor, ku, kv, kw, niter=4000)
  # plot(out.SA.rank1$obj)
  Res = TSCCA_SA(Tensor, ku, kv, kw, J=3, SA_niter=8000)
  out.SA = Res$SA_res
  
  score1 = get_consensus_score(true_bicluster,get_tscca_bicluster(out.cancer1))
  score2 = get_consensus_score(true_bicluster,get_tscca_bicluster(out.cancer2))
  score3 = get_consensus_score(true_bicluster,get_tscca_bicluster(out.cancer3))
  score4 = get_consensus_score(true_bicluster,get_tscca_bicluster(out.cancer4))
  score5 = get_consensus_score(true_bicluster,get_tscca_bicluster(out.scca))
  score6 = get_consensus_score(true_bicluster,get_SA_bicluster(out.SA))
  score7 = get_consensus_score(true_bicluster,get_tscca_bicluster(tscca.out))
  
  scores = c(score1,score2,score3,score4,score5,score6,score7)
  print(scores)
  res.dat[kk,] = scores
}
################################################################################
################################################################################
df = data.frame(method = factor(x = c("cancer1","cancer2","cancer3","cancer4","Sum","Modularity_SA","TSCCA"), 
                                levels=c("cancer1","cancer2","cancer3","cancer4","Sum","Modularity_SA","TSCCA")),
                avg=apply(res.dat,2,mean),
                sd =apply(res.dat,2,sd))

fig = ggplot(df, aes(x=method, y=avg)) + theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +  
  geom_bar(stat="identity", color="black", fill="white", size=1.1, width=0.6) +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2, size=1.1) + 
  theme(axis.text = element_text(colour = "black", size = rel(1.5))) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = "black"))+
  # X和Y轴的线设
  theme(axis.line  = element_line(size = rel(1.5), colour = "black")) +
  # 坐标轴上的小线设
  theme(axis.ticks = element_line(size = rel(1.5), colour = "black")) + 
  theme(axis.title.y = element_text(size = rel(1.5))) +
  labs(title=NULL, x =NULL, y = "Jaccard index") 
ggsave("Fig_result_figure.png", width = 7.5, height = 4.5, units = "in")

save.image("ResultData_simulationStudy.RData")
################################################################################
################################################################################