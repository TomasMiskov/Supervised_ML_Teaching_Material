#Course: Supervised Machine Learning
#Author: Tomas Miskov
#Purpose: Implement kernel ridge regression on Airline dataset
#----------------------------------------------------------------
# PREPARE THE ENVIRONMENT
rm(list=ls())
source("\\Functions_Kernel_Ridge_Regression.R")

#--------------------
# IMPORT AIRLINE DATA
#--------------------
vY <- as.vector(Airline$output)
mX <- as.matrix(Airline[,-4])                                         #drop output from the predictors
mX <- dummy_cols(mX,select_columns=c('airline'),remove_first_dummy=T) #airline 1 is the reference
mX <- mX[,2:10]                                                       #drop the original airline variable
mX_orig <- as.matrix(mX)                                              #save the unscaled X just in case
mX <- scale(mX)
mX <- as.matrix(mX)


#            - 1 -
#--------------------------------
# ESTIMATE MODEL USING RBF KERNEL
#--------------------------------
RBFresults <- KernelRidge(vY, mX, 0.1, kernel = "RBF")
RBFkrr_results <- krr(vY, mX, lambda = 0.1, kernel.type = "RBF")

# Comparison
cat("Mean absolute difference between our solution and krr package: ", mean(abs(RBFresults$y_hats - RBFkrr_results$yhat)))

#----------------------------
# RBF KERNEL CROSS VALIDATION
#----------------------------
set.seed(42)
cv_results <- KernelRidgeCV(vY, mX, kfolds = 10, dLambda = 0.1, kernel = "inpoly")

lambdas <- 10^seq(-8, 8, length.out = 100)
RMSEs <- rep(NA, 100)
for(i in seq(100)){
  cv_krr <- KernelRidgeCV(vY, mX, kfolds = 15, dLambda = lambdas[i], kernel = "inpoly")
  RMSEs[i] <- cv_krr$RMSE
}
plot(log(lambdas), RMSEs)
lambdas[which.min(RMSEs)]
#            - 2 -
#-------------------------------
# ESTIMATE MODEL USING NP KERNEL
#-------------------------------
InpolyResults <- KernelRidge(vY, mX, 0.1, kernel = "inpoly")
Inpolykrr_results <- krr(vY, mX, lambda = 0.1, kernel.type = "nonhompolynom")

# Comparison
cat("Mean absolute difference between our solution and krr package: ", mean(abs(InpolyResults$y_hats - Inpolykrr_results$yhat)))



