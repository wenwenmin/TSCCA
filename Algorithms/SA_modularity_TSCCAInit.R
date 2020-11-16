SA_modularity_TSCCAInit = function(tscca.out, Tensor, ku, kv, kw, J=1, SA_niter=8000){
  # For comparative fairness, SA.modularity and TSCCA have the same input data
  # ---------------------------------------------------------------------------
  # Input
  # Tensor: An array in R can have three dimensions
  # Tensor[,,i] = crossprod(X^i, Y^i)
  # X^i \in R^{n_i \times p} (n_i:samples, p:variables), i = 1, ..., M
  # Y^i \in R^{n_i \times q} (n_i:samples, q:variables), i = 1, ..., M
  # Extrate J sparse factors
  # ---------------------------------------------------------------------------
  
  DIM = dim(Tensor)
  U = tscca.out$U
  V = tscca.out$V
  W = tscca.out$W
  D = tscca.out$d
  
  res = list();  times_vec = rep(0,J)
  
  for(j in 1:J){
    ptm = proc.time()
    
    # ---------------------------------------------------------------------------
    # For comparative fairness, SA.modularity and TSCCA have the same input data
    res[[j]] = SA.modularity(Tensor, ku, kv, kw, SA_niter) # SA_niter
    # ---------------------------------------------------------------------------
    tim = proc.time() - ptm
    times_vec[j] = tim[1]
    
    u = U[,j]
    v = V[,j]
    w = W[,j]
    d = D[j]
    
    # ---------------------------------------------------------------------------
    # Tensor = update.tensor(Tensor, u, v, w, d)
    for(k in 1:length(w)){
      Tensor[,,k] = Tensor[,,k] - d*w[k]*u%*%t(v)
    }
  }
  return(list(res=res, times=times_vec))
}