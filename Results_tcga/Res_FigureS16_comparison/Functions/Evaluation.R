# get_consensus_score
# get_clustering_error_score
# get_relevance_score
# get_recovery_score

library("clue")
# solve_LSAP form "clue" to solve the linear sum assignment problem using the Hungarian method.
get_consensus_score = function(true_biclusters,pre_biclusters){
  # -------------------------------------------------------------------------------------
  # Reference: Hochreiter, Sepp, et al. FABIA: factor analysis for bicluster acquisition. 
  #            Bioinformatics, 2010, 26.12: 1520-1527.
  # -------------------------------------------------------------------------------------
  number_true_biclusters = length(true_biclusters)
  number_pre_biclusters = length(pre_biclusters)
  ## Step1: Compute jaccard_index matrix
  scoreMat = matrix(0, nrow=number_pre_biclusters, ncol=number_true_biclusters)
  for(i in 1:number_pre_biclusters){
    for(j in 1:number_true_biclusters){
      ## Compute jaccard_index score 
      T1 = true_biclusters[[i]]
      P1 = pre_biclusters[[j]]
      # intersect row and column for each bicluster
      n1 = length(intersect(T1$u,P1$u))
      n2 = length(intersect(T1$v,P1$v))
      # Number of Union between T1 and P1 bicluster.
      N = length(T1$u)*length(T1$v) + length(P1$u)*length(P1$v) - n1*n2
      jaccard_index = n1*n2/N 
      scoreMat[i,j] = jaccard_index
    }
  }
  ## Step2: the best matching between sets is found using the Hungarian algorithm
  # 确保scoreMat列的数目大于行
  if(number_true_biclusters < number_pre_biclusters){
    scoreMat2 = t(scoreMat)
  }else{scoreMat2 = scoreMat}
  # Hungarian algorithm
  y = solve_LSAP(scoreMat2, maximum = TRUE)
  
  ## Step3: The final score is the sum of similarities divided by the size of the larger set
  s1 = sum(scoreMat2[cbind(seq_along(y), y)])
  s2 = max(number_true_biclusters,number_pre_biclusters)
  consensus_score = s1/s2 
}
######################################################################################
######################################################################################
get_clustering_error_score = function(true_biclusters,pre_biclusters){
  # -------------------------------------------------------------------------------------
  # Reference: PATRIKAINEN et al. Comparing subspace clusterings. 
  #            IEEE TKDE, 2006, 18.7: 902-916.
  # -------------------------------------------------------------------------------------
  number_true_biclusters = length(true_biclusters)
  number_pre_biclusters = length(pre_biclusters)
  
  ## Step1: Compute intersect matrix between true_biclusters and pre_biclusters
  scoreMat = matrix(0, nrow=number_pre_biclusters, ncol=number_true_biclusters)
  for(i in 1:number_pre_biclusters){
    for(j in 1:number_true_biclusters){
      ## Compute jaccard_index score 
      T1 = true_biclusters[[i]]
      P1 = pre_biclusters[[j]]
      # intersect row and column for each bicluster
      n1 = length(intersect(T1$u,P1$u))
      n2 = length(intersect(T1$v,P1$v))
      # Number of intersect between T1 and P1 bicluster.
      scoreMat[i,j] = n1*n2
    }
  }
  if(number_true_biclusters < number_pre_biclusters){
    scoreMat2 = t(scoreMat)
  }else{scoreMat2 = scoreMat}
  # Hungarian algorithm
  y = solve_LSAP(scoreMat2, maximum = TRUE)
  score1 = sum(scoreMat2[cbind(seq_along(y), y)])
  
  ## Step2: Compute number of union between all biclusters 
  row_max_index = col_max_index = 0
  for(bicluster in true_biclusters){
    row_max_index = max(row_max_index, max(bicluster$u))
    col_max_index = max(col_max_index, max(bicluster$v))
  }
  for(bicluster in pre_biclusters){
    row_max_index = max(row_max_index, max(bicluster$u))
    col_max_index = max(col_max_index, max(bicluster$v))
  }  
  U_index_mat = matrix(0, nrow=row_max_index, ncol=col_max_index)
  for(bicluster in true_biclusters){
    U_index_mat[bicluster$u,bicluster$v] = U_index_mat[bicluster$u,bicluster$v] + 1
  }
  for(bicluster in pre_biclusters){
    U_index_mat[bicluster$u,bicluster$v] = U_index_mat[bicluster$u,bicluster$v] + 1
  }
  # 只要有一个非零元素，说明该元素属于并集合
  U_index_mat[U_index_mat>0]=1
  score2 = sum(U_index_mat)
  ### Step3: compute clustering_error_score
  clustering_error_score = score1/score2
}
######################################################################################
######################################################################################
# average_bicluster_relevance
get_relevance_score = function(true_biclusters,pre_biclusters){
  # true_biclusters = true_bicluster
  # pre_biclusters = get_tscca_bicluster(out.cancer1)
  # -------------------------------------------------------------------------------------
  # Reference: PRELI, Amela, et al. A systematic comparison and evaluation of 
  # biclustering methods for gene expression data. Bioinformatics, 2006, 22.9: 1122-1129.
  # -------------------------------------------------------------------------------------
  number_true_biclusters = length(true_biclusters)
  number_pre_biclusters = length(pre_biclusters)
  ## Step1: Compute jaccard_index matrix
  scoreMat = matrix(0, nrow=number_pre_biclusters, ncol=number_true_biclusters)
  for(i in 1:number_pre_biclusters){
    for(j in 1:number_true_biclusters){
      ## Compute jaccard_index score 
      T1 = true_biclusters[[i]]
      P1 = pre_biclusters[[j]]
      # intersect row and column for each bicluster
      n1 = length(intersect(T1$u,P1$u))
      n2 = length(intersect(T1$v,P1$v))
      # Number of Union between T1 and P1 bicluster.
      N = length(T1$u)*length(T1$v) + length(P1$u)*length(P1$v) - n1*n2
      jaccard_index = n1*n2/N 
      scoreMat[i,j] = jaccard_index
    }
  }
  average_match_score = sum(apply(scoreMat,1,max))/number_pre_biclusters
}
######################################################################################
######################################################################################
# average module recovery
get_recovery_score = function(true_biclusters,pre_biclusters){
  # -------------------------------------------------------------------------------------
  # Reference: PRELI, Amela, et al. A systematic comparison and evaluation of 
  # biclustering methods for gene expression data. Bioinformatics, 2006, 22.9: 1122-1129.
  # -------------------------------------------------------------------------------------
  number_true_biclusters = length(true_biclusters)
  number_pre_biclusters = length(pre_biclusters)
  ## Step1: Compute jaccard_index matrix
  scoreMat = matrix(0, nrow=number_pre_biclusters, ncol=number_true_biclusters)
  for(i in 1:number_pre_biclusters){
    for(j in 1:number_true_biclusters){
      ## Compute jaccard_index score 
      T1 = true_biclusters[[i]]
      P1 = pre_biclusters[[j]]
      # intersect row and column for each bicluster
      n1 = length(intersect(T1$u,P1$u))
      n2 = length(intersect(T1$v,P1$v))
      # Number of Union between T1 and P1 bicluster.
      N = length(T1$u)*length(T1$v) + length(P1$u)*length(P1$v) - n1*n2
      jaccard_index = n1*n2/N 
      scoreMat[i,j] = jaccard_index
    }
  }
  average_match_score = sum(apply(scoreMat,2,max))/number_true_biclusters
}
######################################################################################
######################################################################################