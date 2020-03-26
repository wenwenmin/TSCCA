SCCA = function(XTY, ku=10, kv=100, J=3,  niter=100, err=0.0001){
  # --------------------------------------------------------------------------
  # XTY: miRNA-gene correlation matrix with the samples from all cancer types
  # J: Number of rank
  # --------------------------------------------------------------------------
  # Output
  p = nrow(XTY)
  n = ncol(XTY)
  U = matrix(0,p,J); D = rep(0,J); V = matrix(0,n,J)
  
  tX = XTY
  d_iter = list()
  for(i in 1:J){
    out = rank1.SCCA(tX, ku, kv, niter, err)
    U[,i] = out$u; V[,i] = out$v; D[i] = out$d; 
    
    d_iter[[i]] = out$ds
    # update X
    tX = tX - c(out$d)*out$u%*%t(out$v) 
  }
  return (list(U=U, D=D, V=V, d_iter = d_iter))
}

rank1.SCCA = function(X, ku, kv, niter=1000, err=0.0001){
  # Initialize u and v
  Res = svd(X,1,1)
  u0 = Res$u 
  v0 = Res$v
  
  # save singlar values
  ds = c()
  # Iterative algorithm to solve u and v values
  for (i in 1:niter){
    # update u and v 
    
    u = SCCA_l0_project(X%*%v0,   ku) # for miRNAs
    v = SCCA_l0_project(t(X)%*%u, kv)  #for genes
    
    d = t(u)%*%X%*%v
    ds = c(ds,d)
    # Termination condition
    if ((norm(u - u0,'2')<= err)&(norm(v - v0,'2')<= err)){break}
    else {
      u0 = u;v0 = v}
  }
  return (list(u=u, v=v, d=d, ds=ds))
}
##############################################################################
##############################################################################
# L0 constrained CCA project function
SCCA_l0_project = function(z, k){  
  u = abs(z);
  if(sum(u^2)==0) return(u)
  u[-order(u,decreasing=T)[1:k]] = 0
  u = u/sqrt(sum(u^2))
  u = sign(z)*u
  return(u)
}
##############################################################################
##############################################################################