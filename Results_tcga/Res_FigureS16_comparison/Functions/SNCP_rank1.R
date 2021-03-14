SNCP_rank1 = function(Tensor, ku, kv, kw, Init_flag=1, niter=50, seed0=1, err=0.0001){
  # Sparse Non-negative CP Decomposition
  # Ref.1: Sparse higher-order principal components analysis. Artificial Intelligence and Statistics. 2012.
  # Ref.2: From k-means to higher-way co-clustering: Multilinear decomposition with sparse latent factors. TSP 2013.
  
  Tensor = abs(Tensor)
  ptm = proc.time()
  
   # Initialize_TSCCA
  DIM = dim(Tensor)
  set.seed(seed0)
  u0 = matrix(rnorm(DIM[1]),ncol=1);u0 = u0/norm(u0,'E')
  v0 = matrix(rnorm(DIM[2]),ncol=1);v0 = v0/norm(v0,'E')
  w0 = matrix(rnorm(DIM[3]),ncol=1);w0 = w0/norm(w0,'E')
  
  u0 = abs(u0)
  v0 = abs(v0)
  w0 = abs(w0)
  u = u0; v = v0;w = w0
  
  temp_d = c()
  # Iterative algorithm
  for (i in 1:niter){
    XTY = compute.C(Tensor, w, DIM[3])
    
    z.u = XTY%*%v
    u = SNCP_project(z.u, ku)
    
    z.v = crossprod(XTY,u)
    v = SNCP_project(z.v, kv) 
    
    z.w = compute.zw(Tensor, u, v, DIM[3])
    w = SNCP_project(z.w, kw)
    
    temp_d = c(temp_d, t(w)%*%z.w) 
    
    # Algorithm termination condition
    if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)&(norm(w - w0,'E')<= err)){break}
    else {
      u0=u; v0=v; w0=w}
  }
  tim = proc.time() - ptm
  modularity = mean(abs(Tensor[which(u!=0), which(v!=0), which(w!=0)]))
  
  return (list(u=u, v=v, w=w, d_iter=temp_d, d=temp_d[i], modularity=modularity, time=tim[1]))
}

Initialize_uvw_irlba = function(Tensor,DIM){
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

# rank one 
SNCP_project = function(z, k){
  z[z<0] = 0
  lambda = sort(z, decreasing=T)[k+1]
  z = z - lambda
  z[z<0] = 0

  z = z/sqrt(sum(z^2))
  return(z)
}

SCP_project = function(z, k){
  lambda = sort(abs(z), decreasing=T)[k+1]
  u = abs(z) - lambda
  u[u<0] = 0
  u = sign(z)*u
  u = u/sqrt(sum(u^2))
}