SCCA_rank1 = function(X, ku, kv, seed0=100, niter=1000, err=0.0001){
  ptm = proc.time()
  # Initialize u and v
  # Res = svd(X,1,1)
  # u0 = Res$u 
  # v0 = Res$v
  DIM = dim(X)
  set.seed(seed0)
  u0 = matrix(rnorm(DIM[1]),ncol=1);u0 = u0/norm(u0,'E')
  v0 = matrix(rnorm(DIM[2]),ncol=1);v0 = v0/norm(v0,'E')
  
  # save singlar values
  ds = c()
  # Iterative algorithm to solve u and v values
  for (i in 1:niter){
    # update u and v 
    
    u = SCCA_rank1_l0_project(X%*%v0,   ku) # for miRNAs
    v = SCCA_rank1_l0_project(t(X)%*%u, kv)  #for genes
    
    d = t(u)%*%X%*%v
    ds = c(ds,d)
    # Termination condition
    if ((norm(u - u0,'2')<= err)&(norm(v - v0,'2')<= err)){break}
    else {
      u0 = u;v0 = v}
  }
  tim = proc.time() - ptm; # time = tim[1]
  return (list(u=u, v=v, d=d, ds=ds, time = tim[1]))
}
##############################################################################
##############################################################################
# L0 constrained CCA project function
SCCA_rank1_l0_project = function(z, k){  
  u = abs(z);
  if(sum(u^2)==0) return(u)
  u[-order(u,decreasing=T)[1:k]] = 0
  u = u/sqrt(sum(u^2))
  u = sign(z)*u
  return(u)
}
##############################################################################
##############################################################################