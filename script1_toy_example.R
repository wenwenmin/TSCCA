# ------------------------------------------------------------------------------
# Step 1 generate data
set.seed(1) 
p = 300; q = 30; cancerTypes = 4
Tensor = array(0, c(p,q,cancerTypes))
sd0 = 0.2
Tensor = array(rnorm(p*q*cancerTypes, mean = 0, sd = sd0), c(p,q,cancerTypes))
Tensor[  1:100,  1:10, 1] <- matrix(rnorm(100*10, mean =  0.5, sd = sd0), nrow=100) 
Tensor[101:200, 11:20, 1] <- matrix(rnorm(100*10, mean = -0.5, sd = sd0), nrow=100) 

Tensor[  1:100,  1:10, 2] <- matrix(rnorm(100*10, mean = -0.5, sd = sd0), nrow=100) 
Tensor[201:300, 21:30, 2] <- matrix(rnorm(100*10, mean =  0.5, sd = sd0), nrow=100) 

Tensor[101:200, 11:20, 3] <- matrix(rnorm(100*10, mean =  0.5, sd = sd0), nrow=100) 
Tensor[201:300, 21:30, 3] <- matrix(rnorm(100*10, mean = -0.5, sd = sd0), nrow=100) 

true_cluster = list()
true_cluster[[1]] = list(u=  1:100, v=1:10)
true_cluster[[2]] = list(u=101:200, v=11:20)
true_cluster[[3]] = list(u=201:300, v=21:30)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Runing TSCCA
# sparse level
ku = 100; # extract top 100 genes
kv = 10;  # extract top 10  miRNAs
kw = 2;   # extract top 3   cancers 
# TSCCA algorithm 
source('Functions/TSCCA.R') 
out = TSCCA(Tensor, ku, kv, kw, J=1)
names(out)
# ------------------------------------------------------------------------------