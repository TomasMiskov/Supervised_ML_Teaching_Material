## ---------------------------
## Script Name: MM_From_Scratch.R
## Author: Tomas Miskov
## Date Created: 2022-10-16
## Purpose: Implement Minimization by Majorization from scratch
## ---------------------------

#--------
# SET UP |
#--------
rm(list = ls())                                       # clean the environment
if (!require("pacman")) install.packages("pacman")    # install pacman
pacman::p_load(ggplot2, tidyverse)                    # pre-load packages
# source("functions/packages.R")                      # load local libraries

options(scipen = 6, digits = 4)                       # clean numerical notation
# setwd("C:/Users/misko/")                            # set working directory

#----------------
# SYNTHETIC DATA |
#----------------
set.seed(1011011)
betas <- c(1, 2, 3)
X <- matrix(rnorm(200), ncol = 2)
X <- cbind(rep(1, 100), X)
y <- X %*% betas + rnorm(100)
res <- X %*% betas - y
ssr <- sum(res ^ 2)

#-----------------
# HELPER FUNCTION |
#-----------------

# Residual Sum of Squares
RSS <- function(X, y, b){
  e <- y - X %*% b
  rss <- t(e) %*% e
  return(rss)
}

#---------
# MM LOOP |
#---------
MM <- function(X, y, b0 = as.vector(runif(dim(X)[2])),
               epsilon = 1e-5){
  
  XtX <- t(X) %*% X
  Xty <- t(X) %*% y
  lambda <- max(eigen(XtX)$values)
  
  rssk <- RSS(X, y, b0)
  rsskm1 <- rssk * 0.9
  bkm1 <- b0
  
  k <- 1
  while(k == 1 | (rsskm1 - rssk)/rsskm1 > epsilon){
    k <- k + 1
    bk <- bkm1 - (1 / lambda) * XtX %*% bkm1 + (1 / lambda) * Xty
    rsskm1 <- rssk
    rssk <- RSS(X, y, bk)
    bkm1 <- bk
    
    # print(k)
    # print(all((rsskm1 - rssk) / rsskm1 > epsilon))
  }
  return(bk)
}

b <- MM(X, y)
summary(lm(y~X[, 2:3]))



