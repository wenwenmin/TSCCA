SCP_TSCCAInit = function(tscca.out, Tensor, ku, kv, kw, J=1, niter=20){
  # For comparative fairness, SCP and TSCCA have the same input data
  
  DIM = dim(Tensor)
  U = tscca.out$U
  V = tscca.out$V
  W = tscca.out$W
  D = tscca.out$d
  
  res = list();  times_vec = rep(0,J)
  for(j in 1:J){
    ptm = proc.time()
    # --------------------------------------------------------------------------
    # For comparative fairness, SCP and TSCCA have the same input data
    res[[j]] = SCP_rank1(Tensor, ku, kv, kw, Init_flag=1, niter)
    
    # --------------------------------------------------------------------------
    tim = proc.time() - ptm
    times_vec[j] = tim[1]
    
    u = U[,j]
    v = V[,j]
    w = W[,j]
    d = D[ j]
    # --------------------------------------------------------------------------
    # Tensor = update.tensor(Tensor, u, v, w, d)
    for(k in 1:length(w)){
      Tensor[,,k] = Tensor[,,k] - d*w[k]*u%*%t(v)
    }
    # --------------------------------------------------------------------------
  }
  return (list(res=res, times=times_vec))
}
