
# Stats packages

library(mlbench)
library(MASS)
library(MCMCpack)
library(caret)
library(rstan)
library(glmnet)

# Wrangling packages

library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(readr)

# Visualisation packages

library(ggplot2)
library(gridExtra)
library(scales)
library(GGally)

#### Set seed

set.seed(20190621)

#### Get data

data(BostonHousing) 

d <- data.table(BostonHousing) 

y <- d$medv

# Format chas as numeric

d[, chas := as.numeric(chas) - 1]

#### Standardise the y matrix

y <- log(y)

y.sd <- sd(y)
y.bar <- mean(y)

y <- (y - y.bar) / y.sd

#### Log transform some obvious candidates

d[, crim := log(crim)]
d[, indus := log(indus)]
d[, nox := log(nox)]
d[, dis := log(dis)]
d[, lstat := log(lstat)]
d[, ptratio := log(ptratio)]

#### Recode tax, rad and zn to binaries

d[, zone_over := ifelse(zn == 0, 0, 1)]
d[, tax_over := ifelse(tax > 500, 1, 0)]
d[, rad_over := ifelse(rad > 10, 1, 0)]
d[, c("tax", "rad", "zn") := NULL]


#### Prepare the design matrix

X <- as.matrix(d[, - "medv"])
X <- X[, order(colnames(X))]

#### Standardise the X matrix

X.means <- apply(X, 2, mean)
X.sd <- apply(X, 2, sd)

X <- sweep(X, 2, X.means)
X <- sweep(X, 2, X.sd, FUN = "/")

#### Setup the k-folds

### 20 fold partition

k.number <- 20

folds <- createFolds(y, k = k.number)

y.trains <- list()

y.tests <- list()

X.trains <- list()

X.tests <- list()

for (i in 1:k.number) {
  
  y.trains[[i]] <- y[- folds[[i]]]
  
  y.tests[[i]] <- y[folds[[i]]]
  
  X.trains[[i]] <- X[- folds[[i]], ]
  
  X.tests[[i]] <- X[folds[[i]], ]
  
}

## 5 fold partition

l.k.number <- 5

l.folds <- createFolds(y, k = l.k.number)

ly.trains <- list()

ly.tests <- list()

lX.trains <- list()

lX.tests <- list()

for (i in 1:l.k.number) {
  
  ly.trains[[i]] <- y[- l.folds[[i]]]
  
  ly.tests[[i]] <- y[l.folds[[i]]]
  
  lX.trains[[i]] <- X[- l.folds[[i]], ]
  
  lX.tests[[i]] <- X[l.folds[[i]], ]
  
}






