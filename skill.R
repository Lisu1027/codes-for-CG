##Some functions used in simulations
#Calculate evaluation indicators involving NZ, TPR,FPR and F1 score
NZ <- function(precisionmatrixhat,threshold) {
  betahat = as.vector(precisionmatrixhat)
  NZ = length(which(betahat > threshold))
  return(NZ)
}

TPR <- function(precisionmatrixhat,precisionmatrix,threshold){
  betahat = as.vector(precisionmatrixhat)
  beta = as.vector(precisionmatrix)
  beta_NZ = which(beta > 0)
  betahat_NZ = which(betahat > threshold)
  TP = length(intersect(beta_NZ , betahat_NZ))
  NZ = length(beta_NZ)
  TPR = TP/NZ
  return(TPR)
}

FPR <- function(precisionmatrixhat,precisionmatrix,threshold){
  betahat = as.vector(precisionmatrixhat)
  beta = as.vector(precisionmatrix)
  beta_Z = which(beta == 0)
  betahat_NZ = which(betahat > threshold)
  FP = length(intersect(beta_Z , betahat_NZ))
  Z = length(beta_Z)
  FPR = FP/Z
  return(FPR)
}

F1_score = function(precisionmatrixhat,precisionmatrix,threshold){
  betahat = as.vector(precisionmatrixhat)
  beta = as.vector(precisionmatrix)
  beta_Z = which(beta == 0)
  beta_NZ = which(beta != 0)
  betahat_NZ = which(betahat > threshold)
  betahat_Z = which(betahat <= threshold)
  FP = length(intersect(beta_Z , betahat_NZ))
  TP = length(intersect(beta_NZ , betahat_NZ))
  FN = length(intersect(beta_NZ , betahat_Z))
  precision = TP/(TP+FP)
  recall = TP/(TP+FN)
  result = (2*(precision*recall))/(precision+recall)
  return(result)
}

##Caculate results of given indicators for an estimation
result_out = function(estimate,precisionmatrix,threshold){
  err_2 = sum((estimate - precisionmatrix)*(estimate - precisionmatrix))
  err_2 = sqrt(err_2)
  NZ_1 = NZ(estimate,threshold)
  FPR_1 =FPR(estimate,precisionmatrix,threshold)
  TPR_1 =TPR(estimate,precisionmatrix,threshold)
  F_1 =F1_score(estimate,precisionmatrix,threshold)
  all = c(err_1,err_2,NZ_1,FPR_1,TPR_1,F_1)
  return(all)
}

##Simulate 100 replicates and report the average results with their standard errors 
result_all = function(lambda1,lambda2,lambda3,lambda_ridge,lambda_lasso,lambda_en,lambda_CGL1,lambda_CGL2,lambda_CGL3,alpha1,alpha2,alpha3,alpha_ridge,alpha_lasso,alpha_en,c1,c2,c3,p1,p2,p3,n,p,precisionmatrix,threshold1,threshold2){
  all = matrix(0,100,6)
  all_1 = matrix(0,100,6)
  all_2 = matrix(0,100,6)
  all_3 = matrix(0,100,6)
  all_ridge = matrix(0,100,6)
  all_lasso = matrix(0,100,6)
  all_en = matrix(0,100,6)
  all_CGL = matrix(0,100,6)
  result = matrix(0,6,16)
  precisionmatrix = precision_matrix(c1,c2,c3,p1,p2,p3)
  precisionmatrix1 = precisionmatrix[1:sum(p1),1:sum(p1)]
  precisionmatrix2 = precisionmatrix[(sum(p1)+1):(sum(p1)+sum(p2)),(sum(p1)+1):(sum(p1)+sum(p2))]
  precisionmatrix3 = precisionmatrix[(sum(p1)+sum(p2)+1):p,(sum(p1)+sum(p2)+1):p]
  cov = solve(precisionmatrix)
  for(i in 1:100){
    x = mvrnorm(n, rep(0,p), cov)
    s = cov(x)
    #CG
    result_CG = CG(x)
    estimate_part1 = as.matrix(result_CG[[2]])
    result_part1 = result_out(estimate_part1,precisionmatrix1,threshold1)
    all_1[i,] = result_part1
    
    estimate_part2 = estimate_part1 = as.matrix(result_CG[[3]])  
    result_part2 = result_out(estimate_part2,precisionmatrix2,threshold1)
    all_2[i,] = result_part2
    
    estimate_part3 = estimate_part1 = as.matrix(result_CG[[4]]) 
    result_part3 = result_out(estimate_part3,precisionmatrix3,threshold1)
    all_3[i,] = result_part3
    
    estimate_sum = as.matrix(result_CG[[5]])
    result_sum = result_out(estimate_sum,precisionmatrix,threshold1)
    all[i,] = result_sum
    
    #Graphic Ridge 
    estimate_ridge = gelnet (S=s,lambda=lambda_ridge,alpha=alpha_ridge)$Theta
    result_ridge = result_out(estimate_ridge,precisionmatrix,threshold2)
    all_ridge[i,] = result_ridge
    
    #Graphic Lasso
    estimate_lasso = gelnet (S=s,lambda=lambda_lasso,alpha=alpha_lasso)$Theta
    result_lasso = result_out(estimate_lasso,precisionmatrix,threshold2)
    all_lasso[i,] = result_lasso
    
    #Graphic elastic net,
    estimate_en = gelnet (S=s,lambda =lambda_en,alpha=alpha_en)$Theta
    result_en = result_out(estimate_en,precisionmatrix,threshold2)
    all_en[i,] = result_en
    
    #CGL
    estimate_CGL = CGL(method="single",s,p,lambda_CGL1,lambda_CGL2,lambda_CGL3)
    result_CGL = result_out(estimate_CGL,precisionmatrix,threshold2)
    all_CGL[i,] = result_CGL
  }
  result[,1] = apply(all_1,2,mean)
  result[,2] = apply(all_1,2,sd)
  result[,3] = apply(all_2,2,mean)
  result[,4] = apply(all_2,2,sd)
  result[,5] = apply(all_3,2,mean)
  result[,6] = apply(all_3,2,sd)
  result[,7] = apply(all,2,mean)
  result[,8] = apply(all,2,sd)
  result[,9] = apply(all_ridge,2,mean)
  result[,10] = apply(all_ridge,2,sd)
  result[,11] = apply(all_lasso,2,mean)
  result[,12] = apply(all_lasso,2,sd)
  result[,13] = apply(all_en,2,mean)
  result[,14] = apply(all_en,2,sd)
  result[,15] = apply(all_CGL,2,mean)
  result[,16] = apply(all_CGL,2,sd)
  colnames(result) = c('mean_part1','sd_part1','mean_part2','sd_part2','mean_part3','sd_part3','mean_all','sd_all',
                       'mean_ridge','sd_ridge','mean_lasso','sd_lasso','mean_en','sd_en','mean_CGL','sd_CGL') 
  rownames(result) = c('err1','err2','NZ','FPR','TPR','F1 score')
  return(result)
}

