SA_variance = function(Tensor, ku, kv, kw, niter=1000, flag=1, seed0=1){
  # Use the 3D Merit Function of Eq. (1) 
  # Reference: Triclustering algorithms for three-dimensional data analysis: a comprehensive survey. ACM Computing Surveys.
  
  ptm = proc.time()
  dat.dim = dim(Tensor) # miRNAs*genes*cancers

  # -------------------------------------------------------------
  # Initialize miR.x  gene.x  cancer.x
  sparse_level = list(ku=ku,kv=kv,kw=kw)
  s0 = Initialize_random(Tensor, sparse_level, seed0)
  
  obj = c()
  for(k in 1:niter){
    T0 = niter/k^3 # Temperature
    set.seed(seed0*k) # 这里的随机性非常重要
    if(kw<dat.dim[3]){
      flg = sample(c(1,2,3),1)
    }else{flg = sample(c(1,2),1)} #不用选择cancer维度
    
    Threshold.p = runif(1,0,1)
    if(flg==1){
      i = sample(s0$x1,1)
      xx = 1:dat.dim[1]
      j = sample(xx[-c(s0$x1)],1) #更新点
      
      # find a neighbour point
      x1 = c(s0$x1,j)
      x1 = x1[-which(x1==i)]
      sn = s0; sn$x1 = x1
      
      deta.E = energyFun(Tensor, sn) - energyFun(Tensor, s0)
      if(exp(-deta.E/T0)>=Threshold.p) s0 = sn
    }else if(flg==2){
      i = sample(s0$x2,1)
      xx = 1:dat.dim[2]
      j = sample(xx[-c(s0$x2)],1) 
      
      # find a neighbour point
      x2 = c(s0$x2,j)
      x2 = x2[-which(x2==i)]
      sn = s0; sn$x2 = x2
      
      deta.E = energyFun(Tensor, sn) - energyFun(Tensor, s0)
      if(exp(-deta.E/T0)>=Threshold.p) s0 = sn
    }else{
      i = sample(s0$x3,1)
      xx = 1:dat.dim[3]
      j = sample(xx[-c(s0$x3)],1) 
      x3 = c(s0$x3,j)
      x3 = x3[-which(x3==i)]
      sn = s0; sn$x3 = x3
      
      deta.E = energyFun(Tensor, sn) - energyFun(Tensor, s0)
      if(exp(-deta.E/T0)>=Threshold.p) s0 = sn
    }
    obj = c(obj, energyFun(Tensor,s0))
  }
  tim = proc.time() - ptm; # time = tim[1]
  return(list(u=sort(s0$x1), v=sort(s0$x2), w=sort(s0$x3), obj=obj, opt_obj=obj[niter], time = tim[1]))
}

Initialize_random = function(Tensor, sparse_level, seed0){
  dat.dim = dim(Tensor) # miRNAs*genes*cancers
  set.seed(seed0)
  x1 = sample(1:dat.dim[1], sparse_level$ku, replace=F)
  x2 = sample(1:dat.dim[2], sparse_level$kv, replace=F)
  x3 = sample(1:dat.dim[3], sparse_level$kw, replace=F)
  return(list(x1=x1,x2=x2,x3=x3))
}
# -----------------------------------------
# Modularity
# energyFun = function(Tensor, s){
#   subTensor = abs(Tensor[s$x1,s$x2,s$x3])
#   obj = mean(subTensor)
# }
# -----------------------------------------
energyFun = function(Tensor, s){
  subTensor = abs(Tensor[s$x1,s$x2,s$x3])
  a_IJK = mean(subTensor)
  var = mean((subTensor - a_IJK)^2)
}