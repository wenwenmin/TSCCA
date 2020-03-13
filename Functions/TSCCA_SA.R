TSCCA_SA = function(Tensor, ku, kv, kw, J=1, SA_niter=8000, niter=1000, err=0.0001){
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
  
  SA_res = list()
  for(j in 1:J){
    res = Initialize_TensorCCA(Tensor,DIM)
    
    w = w0 = res$w0; u = u0 = res$u0; v = v0 = res$v0
    temp_d = c()
    
    # Iterative algorithm
    for (i in 1:niter){
      XTY = compute.C(Tensor, w, DIM[3])
      
      SA_res[[j]] = SA.modularity(Tensor, ku, kv, kw, SA_niter) # SA_niter
      
      z.u = XTY%*%v
      u = TSCCA.project(z.u, ku)
      
      z.v = crossprod(XTY,u)
      v = TSCCA.project(z.v, kv) 
      
      z.w = compute.zw(Tensor, u, v, DIM[3])
      w = TSCCA.project(z.w, kw)
      
      temp_d = c(temp_d, t(w)%*%z.w) 
      
      # Algorithm termination condition
      if ((norm(u - u0,'E')<= err)&(norm(v - v0,'E')<= err)&(norm(w - w0,'E')<= err)){break}
      else {
        u0=u; v0=v; w0=w}
    }
    d[j] = t(w)%*%z.w; d_iter[[j]] = temp_d
    U[,j]=u; V[,j]=v; W[,j]=w
    Tensor = update.tensor(Tensor,u,v,w,d[j])
  }
  return (list(SA_res=SA_res, U=U, V=V, W=W, d=d, d_iter=d_iter))
}

# Initialize TSCCA
Initialize_TensorCCA = function(Tensor,DIM){
  # -----------------------------------------------------
  # The first way
  # set.seed(100)
  # u0 = matrix(rnorm(DIM[1]),ncol=1);u0 = u0/norm(u0,'E')
  # set.seed(200)
  # v0 = matrix(rnorm(DIM[2]),ncol=1);v0 = v0/norm(v0,'E')
  # set.seed(300)
  # w0 = matrix(rnorm(DIM[3]),ncol=1);w0 = w0/norm(w0,'E')
  # return(list(w0=w0,u0=u0,v0=v0))
  # -----------------------------------------------------
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
