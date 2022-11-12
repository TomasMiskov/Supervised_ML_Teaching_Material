## ---------------------------
## Script Name: K_fold_cross_val
## Author: Tomas Miskov
## Date Created: 2022-11-03
## Purpose: Implement k-fold cross-validation from scratch
## ---------------------------

#--------
# SET UP |
#--------
rm(list=ls())                                         # clean the environment
if (!require("pacman")) install.packages("pacman")    # install pacman
pacman::p_load(ggplot2, tidyverse)                    # pre-load packages
options(scipen = 6, digits = 4)                       # clean numerical notation

#------
# DATA |
#------
set.seed(1001001)
X <- matrix(rnorm(900), ncol = 9)
b <- c(2,3,1,0,2,3,1,0,3)
y <- 1 + X %*% b + rnorm(100)

#------------------
# HELPER FUNCTIONS |
#------------------
mse <- function(y, y_hat) {
  n <- length(y)
  res <- y - y_hat
  return((t(res) %*% res)/n)
}

rmse <- function(mses) {
  k <- length(mses)
  return((1/k * sum(mses))^(1/2))
}

train_test_split <- function(X, y, prop) {
  n <- length(y)
  prop <- floor(prop * n)
  shuffled_indices <- sample(n, n)
  train_indices <- shuffled_indices[1:prop]
  test_indices <- shuffled_indices[(prop + 1):n]
  
  X_train <- X[train_indices, ]
  y_train <- y[train_indices]
  X_test <- X[test_indices, ]
  y_test <- y[test_indices]
  
  return(list(X_train = X_train, y_train = y_train, 
              X_test = X_test, y_test = y_test))
}

#---------
# K-FOLDS |
#---------
k_folds <- function(n, k) {
  shuffled_indices <- sample(n, n)
  fold_len <- floor(n/k)
  fold_start_indices <- seq(1, (n - fold_len), fold_len)
  
  folds_indices <- lapply(fold_start_indices, 
                          function(start_index) {
                            end_index <- (start_index + fold_len - 1)
                            indices <- shuffled_indices[start_index:end_index]
                            return(indices)
                          })
  
  remaining_obs <- n %% k
  if(remaining_obs > 0) {
    remaining_obs_indices <- shuffled_indices[(n - remaining_obs):n]
    folds_indices[[k]] <- c(folds_indices[[k]], remaining_obs_indices)
  }
  return(folds_indices)
}

#-----------
# CROSS-VAL |
#-----------
cross_val <- function(X, y, k) {
  folds_indices <- k_folds(length(y), k)
  
  mses <- lapply(folds_indices, function(fold_indices) {
    X_test <- X[fold_indices,]
    X_train <- X[-fold_indices,]
    y_test <- y[fold_indices]
    y_train <- y[-fold_indices]
    
    df_train <- data.frame(cbind(y_train, X_train))
    colnames(df_train) <- c('y', seq(9))
    df_test <- data.frame(X_test)
    colnames(df_test) <- seq(9)
    
    model <- lm(y ~ ., data = df_train)
    y_hat <- predict.lm(model, df_test)
    mse_fold <- mse(y_test, y_hat)
    
    return(mse_fold[[1]])
  })
  return(unlist(mses))
}

data <- train_test_split(X, y, 0.8)
mses <- cross_val(data$X_train, data$y_train, 10)
rmse(mses)
