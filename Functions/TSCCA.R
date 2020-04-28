TSCCA = function(Tensor, ku, kv, kw, J=1, niter=1000, err=0.0001){
  # Tensor sparse CCA (TSCCA)
  # Input
  # Tensor: An array in R can have three dimensions
  # Tensor[,,i] = crossprod(X^i, Y^i)
  # X^i \in R^{n_i \times p} (n_i:samples, p:variables), i = 1, ..., M
  # Y^i \in R^{n_i \times q} (n_i:samples, q:variables), i = 1, ..., M
  # Extrate J TSCCA sparse factors
  
  DIM = dim(Tensor)
  U = matrix(0, nrow=DIM[1], ncol=J)
  V = matrix(0, nrow=DIM[2], ncol=J)
  W = matrix(0, nrow=DIM[3], ncol=J)
  d = rep(0, J)
  d_iter = list()
  
  for(j in 1:J){
    res = Initialize_TSCCA(Tensor,DIM)
    w = w0 = res$w0; u = u0 = res$u0; v = v0 = res$v0
    temp_d = c()
    
    # Iterative algorithm
    for (i in 1:niter){
      XTY = compute.C(Tensor, w, DIM[3])
      
      z.u = XTY%*%v
      u = TSCCA.project(z.u, ku)
      
      z.v = crossprod(XTY,u)
      v = TSCCA.project(z.v, kv) 
      
      z.w = compute.zw(Tensor, u, v, DIM[3])
      w = TSCCA.project(z.w, kw)
      
      temp_d = c(temp_d, t(w)%*%z.w) 
      
      # Algorithm termination condition
      if (norm(u - u0,'E') + norm(v - v0,'E') + norm(w - w0,'E') < err){break}
      else {
        u0=u; v0=v; w0=w}
    }
    d[j] = t(w)%*%z.w; d_iter[[j]] = temp_d
    U[,j]=u; V[,j]=v; W[,j]=w
    Tensor = update.tensor(Tensor,u,v,w,d[j])
  }
  return (list(U=U, V=V, W=W, d=d, d_iter=d_iter))
}

Initialize_TSCCA = function(Tensor, DIM, flag=1, seed0=1){
  if(flag==1){
    #print("Initialization with irlba method")
    init_res = Initialize_tscca_irlba(Tensor,DIM)
  }else if(flag==2){
    #print("Initialization with random point")
    init_res = Initialize_tscca_random(Tensor,DIM,seed0)
    
  }else{
    #print("Initialization with cp")
    init_res = Initialize_tscca_cp(Tensor)
  }
  return(init_res)
}

Initialize_tscca_irlba = function(Tensor,DIM){
  # The two way
  w0 = matrix(1/sqrt(DIM[3]), nrow=DIM[3], ncol = 1)
  XTY = compute.C(Tensor, w0, DIM[3])
  
  ## Use irlba efficient partial svd for large problems
  library(irlba)
  if(min(dim(XTY)) < 500) S <- svd(XTY[],nu=1,nv=1)
  else S <- irlba(XTY,nu=1,nv=1)
  u0 <- S$u
  v0 <- S$v
  
  return(list(w0=w0,u0=u0,v0=v0))
}

Initialize_tscca_random = function(Tensor,DIM,seed0=100){
  set.seed(seed0)
  u0 = matrix(rnorm(DIM[1]),ncol=1);u0 = u0/norm(u0,'E')
  v0 = matrix(rnorm(DIM[2]),ncol=1);v0 = v0/norm(v0,'E')
  w0 = matrix(rnorm(DIM[3]),ncol=1);w0 = w0/norm(w0,'E')
  return(list(w0=w0,u0=u0,v0=v0))
}

Initialize_tscca_cp = function(Tensor){
  library(rTensor)
  cp_decomp = cp(as.tensor(Tensor), num_components = 1)
  u0 = cp_decomp$U[[1]]
  v0 = cp_decomp$U[[2]]
  w0 = cp_decomp$U[[3]]
  return(list(w0=w0,u0=u0,v0=v0))
}

# Compute matrix C
compute.C = function(Tensor, w, dim){
  C = w[1]*Tensor[,,1]
  for(i in 2:dim){
    C = C + w[i]*Tensor[,,i]
  }
  return(C)
}
# Compute vector zw
compute.zw = function(Tensor, u, v, dim){
  zw = matrix(0, nrow=dim, ncol=1)
  for(i in 1:dim){
    zw[i] = t(u)%*%Tensor[,,i]%*%v
  }
  return(zw)
}
# Update the tensor
update.tensor = function(Tensor,u,v,w,obj){
  for(k in 1:length(w)){
    Tensor[,,k] = Tensor[,,k] - obj*w[k]*u%*%t(v)
  }
  return(Tensor)
}
# Sparse project function
TSCCA.project = function(z, k){
  if(sum(z^2)==0) return(z)  
  u = abs(z);
  u[-order(u,decreasing=T)[1:k]] = 0
  u = u/sqrt(sum(u^2))
  u = sign(z)*u
  return(u) 
}
