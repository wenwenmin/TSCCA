# ############################################################################
# ############################################################################
SA.modularity = function(Tensor, ku, kv, kw, niter=1000, flag=1, seed0=1){
  dat.dim = dim(Tensor) # miRNAs*genes*cancers
  # -------------------------------------------------------------
  # Initialize miR.x  gene.x  cancer.x
  sparse_level = list(ku=ku,kv=kv,kw=kw)
  if(flag==1){
    #print("Initialization with random point")
    Init_point = Initialize_random(Tensor, sparse_level, seed0)
    x1 = Init_point$x1
    x2 = Init_point$x2
    x3 = Init_point$x3
  }else if(flag==2){
    #print("Initialization with svd")
    Init_point = Initialize_svd(Tensor, sparse_level)
    x1 = Init_point$x1
    x2 = Init_point$x2
    x3 = Init_point$x3
  }else{
    #print("Initialization with cp")
    Init_point = Initialize_cp(Tensor, sparse_level)
    x1 = Init_point$x1
    x2 = Init_point$x2
    x3 = Init_point$x3
  }
  # -------------------------------------------------------------
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
Initialize_random = function(Tensor, sparse_level, seed0){
  dat.dim = dim(Tensor) # miRNAs*genes*cancers
  set.seed(seed0)
  x1 = sample(1:dat.dim[1], sparse_level$ku, replace=F)
  x2 = sample(1:dat.dim[2], sparse_level$kv, replace=F)
  x3 = sample(1:dat.dim[3], sparse_level$kw, replace=F)
  return(list(x1=x1,x2=x2,x3=x3))
}
# ############################################################################
# ############################################################################
Initialize_svd = function(Tensor, sparse_level){
  # Initialize
  DIM = dim(Tensor)
  # The two way
  w0 = matrix(1/sqrt(DIM[3]), nrow=DIM[3], ncol = 1)
  #XTY = compute.C(Tensor, w0, DIM[3])
  XTY = w0[1]*Tensor[,,1]
  for(i in 2:DIM[3]){
    XTY = XTY + w0[i]*Tensor[,,i]
  }
  ## Use irlba efficient partial svd for large problems
  library(irlba)
  if(min(dim(XTY)) < 500) S <- svd(XTY[],nu=1,nv=1)
  else S <- irlba(XTY,nu=1,nv=1)
  u0 <- S$u
  v0 <- S$v
  
  d1_order_index = order(abs(u0), decreasing = T)
  d2_order_index = order(abs(v0), decreasing = T)
  d3_order_index = order(abs(w0), decreasing = T)
  
  return(list(x1=d1_order_index[1:sparse_level$ku],
              x2=d2_order_index[1:sparse_level$kv],
              x3=d3_order_index[1:sparse_level$kw]))
}
# ############################################################################
# ############################################################################
Initialize_cp = function(Tensor, sparse_level){
  library(rTensor)
  cp_decomp = cp(as.tensor(Tensor), num_components = 1)
  d1 = cp_decomp$U[[1]]
  d2 = cp_decomp$U[[2]]
  d3 = cp_decomp$U[[3]]
  
  d1_order_index = order(abs(d1), decreasing = T)
  d2_order_index = order(abs(d2), decreasing = T)
  d3_order_index = order(abs(d3), decreasing = T)
  
  return(list(x1=d1_order_index[1:sparse_level$ku],
              x2=d2_order_index[1:sparse_level$kv],
              x3=d3_order_index[1:sparse_level$kw]))
  
}
# ############################################################################
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