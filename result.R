## The adquate procedure for simulations 
source("data.R")
source('method.R')
source('CV.R')
source('skill.R')
library('MASS')
library('GLassoElnetFast')
library('cluster')
library('glasso')

#generate simulation data
p1 = c(5,5,5,5,5,5,5)
p2 = c(5,5,5,5,5,5,5)
p3 = 330
n = 200
p = sum(p1)+sum(p2)+sum(p3)
precisionmatrix = precision_matrix(c1=0.9,c2=0.6,c3=0,p1,p2,p3)
precisionmatrix1 = precisionmatrix[1:sum(p1),1:sum(p1)]
precisionmatrix2 = precisionmatrix[(sum(p1)+1):(sum(p1)+sum(p2)),(sum(p1)+1):(sum(p1)+sum(p2))]
precisionmatrix3 = precisionmatrix[(sum(p1)+sum(p2)+1):p,(sum(p1)+sum(p2)+1):p]
cov = solve(precisionmatrix)
x = mvrnorm(n, rep(0,p), cov)
s = cov(x)

##Choose lambdas by Cross-validation
#Output final lambdas of CG
alpha = seq(0,1,0.1)
lambda = 0.9^seq(-20,60,1)
lambda_CG = err_result(alpha,lambda,x1,x2,x3,s1,s2,s3,precisionmatrix1,precisionmatrix2,precisionmatrix3,precisionmatrix,fold=5)

#Output final lambdas of Graphic Ridge,Graphic elastic net,GLasso and CGL
alpha = seq(0.1,0.9,0.1)
lambda = 0.9^seq(-20,60,2)
lambda_competive = all_result(x,lambda,alpha,precisionmatrix)

lambda1 = lambda_CG[[1]][1]
alpha1 = lambda_CG[[1]][2]
lambda2 = lambda_CG[[2]][1]
alpha2 = lambda_CG[[2]][2]
lambda3 = lambda_CG[[3]][1]
alpha3 = lambda_CG[[3]][2]
lambda_ridge = lambda_competive[[3]][1] 
alpha_ridge = lambda_competive[[3]][2]
lambda_lasso = lambda_competive[[4]][1]
alpha_lasso = lambda_competive[[4]][2]
lambda_en = lambda_competive[[2]][1]
alpha_en = lambda_competivel[[2]][2]
lambda_CGL1 = lambda_competive[[5]][1]
lambda_CGL2 = lambda_competive[[5]][2]
lambda_CGL3 = lambda_competive[[5]][3]

#simulate 100 replicates for the given setting and report the average results with their standard errors
result = result_all(lambda1,lambda2,lambda3,lambda_ridge,lambda_lasso,lambda_en,lambda_CGL1,lambda_CGL2,lambda_CGL3,
                           alpha1,alpha2,alpha3,alpha_ridge,alpha_lasso,alpha_en,alpha_our,c1=0.9,c2=0.6,c3=0,p1,p2,p3,n,p,precisionmatrix,threshold)

