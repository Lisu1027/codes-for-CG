##Cross-validation for CGL
gelnet_cv = function(lambda,x,precisionmatrix){
  p = dim(x)[2]
  n = dim(x)[1]
  s = cov(x)
  if(lambda[1]==0 & p>n){
    precision = solve(1/(n-1)*(t(x - mean(x))%*%(x - mean(x))))
    err = sum((precision - precisionmatrix)*(precision - precisionmatrix))
  }else if(lambda[1]==0 & p<=n){
    err = 10000
  }else{
    precision = gelnet ( S = s , lambda = lambda[1], alpha = lambda[2])$Theta 
    err = sum((precision - precisionmatrix)*(precision - precisionmatrix))
  }
  return(err)
}

cv_err <- function(x,lambda,alpha = alpha,precisionmatrix,fold = 5){
  p = dim(x)[2]
  n = dim(x)[1]
  num = sample(n)
  l = n / fold
  lambda0 = matrix(NA, length(lambda) * length(alpha), 2)
  lambda0[, 1] = rep(lambda,length(alpha))
  lambda0[, 2] = c(matrix(rep(alpha, length(lambda)), length(lambda), byrow = TRUE))
  err<-matrix(NA,fold,length(lambda) * length(alpha))
  for (i in 1:fold) {
    xtest = x[num[((i - 1) * l + 1):(i * l)],]
    xtrain = x[-num[((i - 1) * l + 1):(i * l)],]
    result = apply(lambda0, 1, gelnet_cv, xtrain,precisionmatrix)
    err[i,] = result
  }
  err_mean = apply(err, 2, mean)
  j = which.min(err_mean)
  lambda_final = lambda0[j,]
  return(lambda_final)
}

#output final lambdas of different clusters based on CG
err_result = function(alpha,lambda,x1,x2,x3,s1,s2,s3,precisionmatrix1,precisionmatrix2,precisionmatrix3,precisionmatrix,fold=5){
  result = list()
  result_lambda1 = cv_err(x = x1,lambda = lambda ,alpha = alpha,precisionmatrix = precisionmatrix1,fold = fold)
  result[[1]] = result_lambda1
  result_lambda2 = cv_err(x = x2,lambda = lambda ,alpha = alpha,precisionmatrix = precisionmatrix2,fold = fold)
  result[[2]] = result_lambda2
  result_lambda3 = cv_err(x = x3,lambda = lambda ,alpha = 1,precisionmatrix = precisionmatrix3,fold = fold)
  result[[3]] = result_lambda3
  return(result)
}

##Cross-validation for Graphic elastic net, Graphic lasso, Graphic ridge and CGL
#Cross-validation for CGL
CGL_cv = function(rho,x,precisionmatrix){
  p = dim(x)[2]
  n = dim(x)[1]
  s = cov(x)
  if(rho==0 & p>n){
    precision = solve(1/(n-1)*(t(x - mean(x))%*%(x - mean(x))))
    err = sum((precision - precisionmatrix)*(precision - precisionmatrix))
  }else if(rho==0 & p<=n){
    err = 10000
  }else{
    precision = glasso(s, rho = rho)$wi
    err = sum((precision - precisionmatrix)*(precision - precisionmatrix))
  }
  return(err)
}

#output final lambdas of CGL
CGL_cv_err <- function(x,rho,precisionmatrix,fold = 5){
  p = dim(x)[2]
  n = dim(x)[1]
  num = sample(n)
  l = n / fold
  err = NULL
  for (i in 1:fold){
    print(i)
    xtest = x[num[((i - 1) * l + 1):(i * l)],]
    xtrain = x[-num[((i - 1) * l + 1):(i * l)],]
    err_fen = NULL
    for (j in 1:length(rho)){
      result =  glasso_cv(rho[j],xtrain,precisionmatrix)
      err_fen = append(err_fen,result)
    }
    err = rbind(err,err_fen)
  }
  err_mean = apply(err, 2, mean)
  j = which.min(err_mean)
  rho_final = rho[j]
  return(rho_final)
}

#output final lambdas of Graphic elastic net, Graphic lasso, Graphic ridge and CGL
all_result = function(x,lambda,alpha,precisionmatrix){
  result = list()
  ERR_en = c()
  #Graphic elastic net
  for(t in 1:length(alpha)){
    fitGelnetIdCV_en = crossvalidation(Y = x, lambda = lambda, alpha = alpha[t],penalize.diagonal = FALSE)
    w_en_part = fitGelnetIdCV_en$Theta
    err_en_part = sum((w_en_part - precisionmatrix)*(w_en_part - precisionmatrix))
    ERR_en = append(ERR_en,err_en_part)
    if(err_en_part<=min(ERR_en)){
      alpha_en = alpha[t]
      lambda_en = fitGelnetIdCV_en$optimal
      w_en = w_en_part
      err_en = err_en_part 
    }
  }
  
  #Graphic lasso
  fitGelnetIdCV_glasso = crossvalidation( Y = x , lambda = lambda, alpha =1)
  lambda_glasso = fitGelnetIdCV_glasso$optimal
  
  #Graphic ridge
  fitGelnetIdCV_ridge = crossvalidation( Y = x , lambda = lambda, alpha =0)
  lambda_ridge = fitGelnetIdCV_ridge$optimal
  
  #CGL
  lambda_CGL = CGL_cv_err(x,lambda,precisionmatrix,fold = 5)
  
  result[[1]] = c('elastic net','ridge','lasso','CGL')                      
  result[[2]] = c(lambda_en,alpha_en)                                      
  result[[3]] = c(lambda_ridge,0)                                           
  result[[4]] = c(lambda_glasso,1)                                          
  result[[5]] = lambda_CGL
  return(result)
}

