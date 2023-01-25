## ---------------------------
## Script Name: decision_trees_tutorial.R
## Author: Tomas Miskov
## Date Created: 2022-11-24
## Purpose:Understanding decision trees and imbalanced data
## ---------------------------

#--------
# SET UP |
#--------
rm(list=ls())                                         # clean the environment
if (!require("pacman")) install.packages("pacman")    # install pacman
pacman::p_load(rpart, caret)                          # pre-load packages
options(scipen = 6, digits = 4)                       # clean numerical notation

#------
# DATA |
#------
set.seed(1101001)
x <- round(runif(100), 2)
y <- round(runif(100), 2)

# change the ratio "r" to create more balanced or less balanced dataset
# r = 0 should give a very balanced dataset, r > 0 increases the imbalance
# the closer the "r" is to 1, the more imbalanced the two classes will be

r <- 3/4
c <- ifelse(x + y >= 1 & runif(100) > r, "red3", "steelblue1")
data <- data.frame(x, y, c)
plot(x, y, col = c, pch = 19)
barplot(prop.table(table(data$c)),
        col = c("red3", "steelblue1"),
        ylim = c(0,1),
        names = c("red", "blue"))

data.train <- data[1:70, ]
data.test <- data[71:100, ]

plot(data.train$x, data.train$y, col = data.train$c, pch = 19)
plot(data.test$x, data.test$y, col = data.test$c, pch = 19)

#------
# GINI |
#------
giniTree <- rpart(c ~ ., data = data.train, parms = list(split = "gini"))
rpart.plot::rpart.plot(giniTree)

giniPred <- predict(giniTree, newdata = data.test, type = "class")
confusionMatrix(as.factor(data.test$c), giniPred)

#---------
# ENTROPY |
#---------
infoTree <- rpart(c ~ ., data = data.train, parms = list(split = "information"))
rpart.plot::rpart.plot(infoTree)

infoPred <- predict(infoTree, newdata = data.test, type = "class")
confusionMatrix(as.factor(data.test$c), infoPred)

#-------
# PRIOR |
#-------
priorTree <- rpart(c ~ ., data = data.train, parms = list(split = "information",
                                                          prior = c(0.5, 0.5)))
rpart.plot::rpart.plot(priorTree)
priorPred <- predict(priorTree, newdata = data.test, type = "class")
confusionMatrix(as.factor(data.test$c), priorPred)
