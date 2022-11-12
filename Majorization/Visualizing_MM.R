## ---------------------------
## Script Name: Visualizing_MM
## Author: Tomas Miskov
## Date Created: 2022-10-28
## Purpose: Visualize the majorization idea on linear regression
## ---------------------------

#--------
# SET UP |
#--------
rm(list=ls())                                         # clean the environment
if (!require("pacman")) install.packages("pacman")    # install pacman
pacman::p_load(ggplot2, tidyverse)                    # pre-load packages
# source("functions/packages.R")                      # load local libraries

options(scipen = 6, digits = 4)                       # clean numerical notation
# setwd("C:/Users/misko/")                            # set working directory

#------
# DATA |
#------
x <- seq(0, 99)
y <- x
lambda <- eigen(t(x)%*%x)$values[1] * 5

#-------------------------
# RSS as function of Beta |
#-------------------------
rss <- function(b) {
  return(t(b*x - y) %*% (b*x - y))
}

#---------------------
# Majorizing function |
#---------------------
mm <- function(b, b0) {
  return(lambda * b * b - 
           2 * lambda * b %*% (b0 - 1/lambda * t(x) %*% x * b0 + 
                               1/lambda * t(x) %*% y)
         + lambda * b0 * (1 - 1/lambda * t(x) %*% x) * b0 
         + t(y) %*% y)
}

test_betas <- seq(0, 2, 0.01)
test_rss <- mapply(rss, test_betas)
test_mm <- mapply(mm, test_betas, 0.75)

plot(test_betas, test_rss, type = 'l', lty = 1)
lines(test_betas, test_mm, type = 'l', lty = 1)
