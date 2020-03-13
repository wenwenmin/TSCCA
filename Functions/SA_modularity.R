SA.modularity = function(Tensor, ku, kv, kw, niter=1000, seed0=1){
  dat.dim = dim(Tensor) # miRNAs*genes*cancers
  #初始化 miR.x  gene.x  cancer.x
  set.seed(seed0)
  x1 = sample(1:dat.dim[1], ku, replace=F)
  x2 = sample(1:dat.dim[2], kv, replace=F)
  x3 = sample(1:dat.dim[3], kw, replace=F)
  
  obj = c()
  for(k in 1:niter){
    T0 = niter/k # Temperature
    # get.SA.newSolution
    # 随机选择更新miR,gene,cancer三个维度中的一个
    set.seed(seed0*k) # 这里的随机性非常重要
    if(kw<dat.dim[3]){
      flg = sample(c(1,2,3),1)
    }else{flg = sample(c(1,2),1)}#不用选择cancer维度
    Threshold.p = runif(1,0,1)
    if(flg==1){
      i = sample(x1,1)
      xx = 1:dat.dim[1]
      j = sample(xx[-c(x1)],1) # 更新点
      deta.E = sum(abs(Tensor[j,x2,x3])) - sum(abs(Tensor[i,x2,x3]))
      if(exp(deta.E/T0)>=Threshold.p){
        x1 = c(x1,j)
        x1 = x1[-which(x1==i)]}
    }else if(flg==2){
      i = sample(x2,1)
      xx = 1:dat.dim[2]
      j = sample(xx[-c(x2)],1) 
      deta.E = sum(abs(Tensor[x1,j,x3])) - sum(abs(Tensor[x1,i,x3]))
      if(exp(deta.E/T0)>=Threshold.p){
        x2 = c(x2,j)
        x2 = x2[-which(x2==i)]}
    }else{
      i = sample(x3,1)
      xx = 1:dat.dim[3]
      j = sample(xx[-c(x3)],1) 
      deta.E = sum(abs(Tensor[x1,x2,j])) - sum(abs(Tensor[x1,x2,i]))
      if(exp(deta.E/T0)>=Threshold.p){
        x3 = c(x3,j)
        x3 = x3[-which(x3==i)]
      }
    }
    obj = c(obj,mean(abs(Tensor[x1,x2,x3])))
  }
  return(list(u=sort(x1), v=sort(x2), w=sort(x3),obj=obj))
}
# ############################################################################
# ############################################################################
# 1 generate data
# p = 120
# q = 150
# cancerTypes = 4
# sd0 = 0.2
# 
# set.seed(1);
# Tensor = array(rnorm(p*q*cancerTypes, mean = 0, sd = sd0), c(p,q,cancerTypes))
# Tensor[1:40,1:50,1] <- matrix(rnorm(40*50, mean = -0.5, sd = sd0), nrow=40)
# Tensor[1:40,1:50,2] <- matrix(rnorm(40*50, mean =  0.5,   sd = sd0), nrow=40)
# out3 = SA.modularity(Tensor, 40, 50, 3, niter=4000)
# plot(out3$obj)
# ############################################################################
# ############################################################################