##CG(our method)
#This Function can be adjusted according to the number of culsters(k)
#Set k=3 be an example
CG = function(x){
  result = list()
  x = mvrnorm(n, rep(0,p), cov)
  s = cov(x)
  #obtain clusters by kmediods
  v_diag = as.vector(diag(s))
  km_result = pam(v_diag,k=3,nstart = 30)
  r_1 = which(km_result$cluster == 1)
  r_2 = which(km_result$cluster == 2)
  r_3 = which(km_result$cluster == 3)
  vec1 = x[,r_1]
  vec2 = x[,r_2]
  vec3 = x[,r_3]
  s1 = cov(vec1)
  s2 = cov(vec2)
  s3 = cov(vec3)
  d = c(max(diag(s1)),max(diag(s2)),max(diag(s3)))
  index1 = which(d == max(d))
  index2 = which(d == median(d))
  index3 = which(d == min(d))
  r_1 = which(km_result$cluster == index1)
  r_2 = which(km_result$cluster == index2)
  r_3 = which(km_result$cluster == index3)
  vec1 = x[,r_1]
  vec2 = x[,r_2]
  vec3 = x[,r_3]
  s1 = cov(vec1)
  s2 = cov(vec2)
  s3 = cov(vec3)
  precisionmatrix1 = precisionmatrix[r_1,r_1]
  precisionmatrix2 = precisionmatrix[r_2,r_2]
  precisionmatrix3 = precisionmatrix[r_3,r_3]
  #estimate respectively
  estimate_part1 = gelnet ( S = s1 , lambda = lambda1, alpha = alpha1)$Theta 
  estimate_part2 = gelnet ( S = s2 , lambda = lambda2, alpha = alpha2)$Theta  
  estimate_part3 = gelnet ( S = s3 , lambda = lambda3, alpha = alpha3)$Theta
  #combine different estimations
  estimate_sum = matrix(0,p,p)
  estimate_sum[r_1,r_1] = estimate_part1
  estimate_sum[r_2,r_2] = estimate_part2
  estimate_sum[r_3,r_3] = estimate_part3
  result[[1]] = c('part1','part2','part3','final')
  result[[2]] = estimate_part1
  result[[3]] = estimate_part2
  result[[4]] = estimate_part3
  result[[5]] = estimate_sum
  return(result)
}


##CGL(a competive method)
#This Function can be adjusted according to the number of culsters(k)
#Set k=3 be an example
CGL = function(method,s,p,rho1,rho2,rho3){
  s_star = generate_s_star(3,p,T = 200,s)
  s_star_dist = as.dist(s_star)
  SLC_hclust = hclust(s_star_dist,method = method)
  out_SLC = cutree(SLC_hclust,k = 3)
  ps = as.data.frame(table(out_SLC))
  x1 = x[,which(out_SLC == 1)]
  x2 = x[,which(out_SLC == 2)]
  x3 = x[,which(out_SLC == 3)]
  if(ps$Freq[1] != 1){
    s1 = cov(x1)
    precisionmatrix1 = precisionmatrix[which(out_SLC == 1),which(out_SLC == 1)]
    precision_estimation1 = glasso(s1, rho = rho1)$wi
  }else{
    precision_estimation1 = 1
  }
  if(ps$Freq[2] != 1){
    s2 = cov(x2)
    precisionmatrix2 = precisionmatrix[which(out_SLC == 2),which(out_SLC == 2)]
    precision_estimation2 = glasso(s2, rho = rho2)$wi
  }else{
    precision_estimation2 = 1
  }
  if(ps$Freq[3] != 1){
    s3 = cov(x3)
    precisionmatrix3 = precisionmatrix[which(out_SLC == 3),which(out_SLC == 3)]
    precision_estimation3 = glasso(s3, rho = rho3)$wi
  }else{
    precision_estimation3 = 1
  }
  precision_estimation = matrix(0,p,p)
  precision_estimation[1:ps$Freq[1],1:ps$Freq[1]] = precision_estimation1
  precision_estimation[(ps$Freq[1]+1):(ps$Freq[1]+ps$Freq[2]),
                       (ps$Freq[1]+1):(ps$Freq[1]+ps$Freq[2])] = precision_estimation2
  precision_estimation[(ps$Freq[1]+ps$Freq[2]+1):dim(precision_estimation)[1],
                       (ps$Freq[1]+ps$Freq[2]+1):dim(precision_estimation)[1]] = precision_estimation3
  num_order = c(which(out_SLC == 1),which(out_SLC == 2),which(out_SLC == 3))
  precision_estimation_new = cbind(precision_estimation,num_order)
  precision_estimation_new = precision_estimation_new[order(num_order),]
  precision_estimation_new = precision_estimation_new[,-(p+1)]
  precision_estimation_new = rbind(precision_estimation,num_order)
  precision_estimation_new = precision_estimation_new[,order(num_order)]
  precision_estimation_new = precision_estimation_new[-(p+1),]
  return(precision_estimation_new)
}

##select K
select_K = function(T = 200,k = 3:15,method = 'average',p,s){
  mk_all_T = NULL
  for (i in 1:T){
    s_star = generate_s_star(i,p,T,s)
    s_star_dist = as.dist(s_star)
    SLC_hclust = hclust(s_star_dist,method = method)
    plot(SLC_hclust)
    mk_all = NULL
    for(m in k){
      out_SLC = cutree(SLC_hclust,k = m)
      count_ck = as.data.frame(table(out_SLC))
      count_ck$count_square = count_ck$Freq^2
      cluster_zuhe = as.data.frame(t(combn(m,2)))
      cluster_zuhe = rbind(cluster_zuhe,as.data.frame(matrix(c(1:m,1:m),m,2)))
      colnames(cluster_zuhe) = c('i','j')
      B_out = NULL
      for(q in 1:length(cluster_zuhe[,1])){
        num_order = NULL
        if(cluster_zuhe[q,1] == cluster_zuhe[q,2]){
          Bij = Bij_same(out_SLC,cluster_zuhe,q,s_star)
          B_out = append(B_out,Bij)
        }
        if(cluster_zuhe[q,1] != cluster_zuhe[q,2]){
          Bij = Bij_difference(out_SLC,cluster_zuhe,q,s_star,count_ck)
          B_out = append(B_out,Bij)
        }
      }
      cluster_zuhe$B_out = B_out
      B = matrix(0,p,p)
      for(i in 1:p){
        for(j in 1:p){
          B[i,j] = cluster_zuhe[intersect(which((cluster_zuhe[,1] == out_SLC[i]) == TRUE),
                                          which((cluster_zuhe[,2] == out_SLC[i]) == TRUE)),3]
        }
      }
      B_sstar = B-s_star
      diag(B_sstar) = 0
      mk_all = append(mk_all,sum((B_sstar)^2)/(p*(p-1)/T))
    }
    mk_all_T = rbind(mk_all_T,mk_all)
  }
  colnames(mk_all_T) = paste("k",3:15,sep = '_')
  return(mk_all_T)
}
