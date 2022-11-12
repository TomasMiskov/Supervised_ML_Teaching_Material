#Course: Supervised Machine Learning
#Author: Tomas Miskov
#Purpose: Function library for kernel ridge regression
#----------------------------------------------------------------
# PRELOAD PACKAGES
#------------------
# note: `dsmle` package needs to be installed manually

if (!require("pacman")) install.packages("pacman")
pacman::p_load(fastDummies, rdetools, docstring, dsmle, matlib, pracma, Ecdat)

#-----------------------------
# SPECIALIZED KERNEL FUCNTIONS
#-----------------------------

RBFkernel_matrix <- function(mX, dGamma, scaled = TRUE, mXu = NULL){
  #' Calculate the Radial Basis Function kernel matrix
  #' 
  #' @param mX The matrix of explanatory variables
  #' @param dGamma Float number that determines the size of the kernel
  #' @param scaled Indication whether mX is Z-scores. If TRUE, sets Gamma = 1/n
  #' @param mXu Matrix of out-of-sample variables used for out-of-sample predictions
  #' 
  #' @return mK: RBF kernel matrix corresponding to mX (or mX*mXu)
  
  if(scaled){
    dGamma <- 1/ncol(mX)           #optimal Gamma value from slide 44 (*find reference for report)
  }
  if(!is.null(mXu)){
    mK <- exp(as.matrix(-dGamma*distmat(mXu, mX)^2))
  } else {
    mK <- exp(as.matrix(-dGamma*dist(mX)^2))
  }
  
  return(mK)
}

InPolyKernel_matrix <- function(mX, dD, mXu = NULL){
  #' Calculate the Inhomogeneous polynomial kernel matrix
  #' 
  #' @param mX The matrix of explanatory variables
  #' @param dD Float number that determines the degree of the polynomial
  #' @param mXu Matrix of out-of-sample variables used for out-of-sample predictions
  #' 
  #' @return mK: Inhomogeneous Polynomial kernel matrix corresponding to mX
  
  if(is.null(mXu)){
    mK <- (1 + mX%*%t(mX)) ^ dD
  } else{
    mK <- (1 + mXu%*%t(mX)) ^ dD
  }
  return(mK)
}

#-------------------------
# KERNEL RIDGE REGRESSION
#-------------------------

KernelRidge <- function(vY, mX, dLambda, kernel = c("RBF", "inpoly"), 
                        dD = 2, dGamma = 1, scale = TRUE, eigen = TRUE){
  #' Solve Kernel Ridge Regression in dual formulation
  #' 
  #' @param vY Vector of the response variable values
  #' @param mX The matrix of explanatory variables
  #' @param dLambda Float coefficient determining the "strength" of the penalty
  #' @param kernel Type of kernel used
  #' @param dD Float number that determines the degree of the polynomial
  #' @param dGamma Float number determining the width of the RBF kernel
  #' @param scale Scales and centers the matrix mX of explanatory variables
  #' @param eigen If TRUE, the solution is arrived to through eigendecomposition of the kernel matrix
  #' 
  #' @return results: list of elements necessary for further analysis/prediction

  
  iN <- nrow(mX)
  iK <- ncol(mX)
  
  if(scale){
    mXmean <- colMeans(mX)
    mXsd <- apply(mX, 2, sd)
    # mX <- mX - matrix(mXmean, nrow = iN, ncol = iK, byrow = TRUE)     #centering
    # mX <- mX/matrix(mXsd, nrow = iN, ncol = iK, byrow = TRUE)         #scaling
    mX <- scale(mX)
  }
  
  # Compute the kernel matrix
  if(kernel == "RBF"){
    if(scale){
      mK <- RBFkernel_matrix(mX, dGamma, scaled = TRUE)
    } else {
      mK <- RBFkernel_matrix(mX, dGamma, scaled = FALSE)
    }
  }
  if(kernel == "inpoly"){
    mK <- InPolyKernel_matrix(mX, dD)
  }
  if(!(kernel %in% c("RBF", "inpoly"))){
    stop("Wrong kernel specified")
  }
  
  # Dual formulation analytical solution
  if(!eigen){
    iN <- nrow(mX)
    mJ <- diag(iN) - 1/iN
    tilde_mK <- mJ %*% mK %*% mJ               #centering the kernel space
    inverse_tilde_mK <- ifelse(det(tilde_mK) < 1e-10, Ginv(tilde_mK), solve(tilde_mK))
    dW0 <- mean(vY)
    tilde_q <- solve((diag(iN) + dLambda * inverse_tilde_mK)) %*% (mJ %*% vY)
    y_hats <- dW0 + tilde_q
  } else {
  
    # Dual formulation eigendecomposition
    iN <- nrow(mX)
    mJ <- diag(iN) - 1/iN
    tilde_mK <- mJ %*% mK %*% mJ               #centering the kernel space
    inverse_tilde_mK <- NULL                   #inverse not necessary
    dW0 <- mean(vY)
    eigK <- eigen(tilde_mK, symmetric = TRUE)
    mK_eigvec <- eigK$vectors
    vK_eigval <- eigK$values
    vDiag_val <- vK_eigval / (vK_eigval + dLambda)
    
    mD <- diag(vDiag_val)
    
    tilde_q <- mK_eigvec %*% mD %*% t(mK_eigvec) %*% scale(vY,center=T,scale=F)
    y_hats <- dW0 + tilde_q
  }

  
  results <- list(mX = mX, mXmean = mXmean, mXsd = mXsd, vY = vY, mK = mK, tilde_mK = tilde_mK, 
                  inverse_tilde_mK = inverse_tilde_mK, eigK = eigK, dW0 = dW0, tilde_q = tilde_q,
                  y_hats = y_hats, dGamma = dGamma, dD = dD, kernel = kernel)
  return(results)
}

#-------------------------
# KERNEL RIDGE PREDICTION
#-------------------------
KernelRidgePredict <- function(res, mXu){
  #' Predict the out-of-sample values from KRR model
  #' 
  #' @param res List of results from the Kernel Ridge estimation function
  #' @param mXu The matrix of out-of-sample explanatory variables
  #' 
  #' @return y_hatsu: vector of out-of-sample predictions
  
  iN <- nrow(res$mK)
  
  # Centering and scaling mXu such that our intercept is the same for both training and testing
  mXu <- mXu - matrix(res$mXmean, nrow = nrow(mXu), ncol = ncol(mXu), byrow = TRUE) 
  mXu <- mXu/matrix(res$mXsd, nrow = nrow(mXu), ncol = ncol(mXu), byrow = TRUE)
  
  if(res$kernel == "RBF"){
    mKu <- RBFkernel_matrix(mX, res$dGamma, mXu = mXu)
  } else if(res$kernel == "inpoly"){
    mKu <- InPolyKernel_matrix(mX, res$dD, mXu = mXu)
  }
  
  #Eigendecompose K
  vK_eigval <- eigen(res$tilde_mK, symmetric = TRUE)$values
  mD_sqrd_inv <- diag(ifelse(res$eigK$values <= 1e-10, 0, 1/res$eigK$values))
  mU <- res$eigK$vectors
  
  y_hatsu <- res$dW0 + ((mKu %*% mU) %*% (mD_sqrd_inv %*% t(mU))) %*% res$y_hats
  return(y_hatsu)
}

#------------------------------
# KERNEL RIDGE CROSS VALIDATION
#------------------------------
KernelRidgeCV <- function(vY, mX, kfolds = 10, dLambda, 
                          kernel = c("RBF", "inpoly"), dD = 2, dGamma = 1){
  #' K-fold cross-validation for a particular value of Lambda parameter
  #' 
  #' @param vY Vector of the response variable values
  #' @param mX The matrix of explanatory variables
  #' @param kfolds Number of folds
  #' @param dLambda Float coefficient determining the "strength" of the penalty
  #' @param kernel Type of kernel used
  #' @param dD Float number that determines the degree of the polynomial
  #' @param dGamma Float number determining the width of the RBF kernel
  #' 
  #' @return results: List of values necessary for further analysis

  iN <- nrow(mX)
  vFoldInd <- rep(sample(seq(kfolds)), ceiling(iN/kfolds))[1:iN] #vector of fold indices
  
  vMSEs <- rep(NA, length(kfolds))
  for(i in seq(kfolds)){
    mX.train <- mX[!(vFoldInd == i),]
    vY.train <- vY[!(vFoldInd == i)]
    mX.test <- mX[vFoldInd == i,]
    vY.test <- vY[vFoldInd == i]
    
    train_results <- KernelRidge(vY = vY, mX = mX, dLambda = dLambda,
                                 kernel = kernel, dD = dD, dGamma = dGamma)
    y_hatsu <- KernelRidgePredict(train_results, mX.test)
    
    # cat("The difference: ", mean(abs(y_hatsu - vY.test)))
    
    dMSE <- mean(t((vY.test - y_hatsu)) %*% (vY.test - y_hatsu))
    vMSEs[i] <- dMSE
  }
  dRootMSE <- ((1/length(vY.test)) * sum(vMSEs))^(0.5) #RMSE
  results <- list(RMSE = dRootMSE, MSEs = vMSEs, dLambda = dLambda)
  return(results)
}
