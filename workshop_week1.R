########################################################
## Workshop week 1
########################################################

# You can use the following command to clear your workspace

rm(list=ls())

##It's good practice to start all of your R scripts with this command in order to clear your workspace

########################################################
## Setting your working directory
########################################################

## To change the working directory, use "setwd"

setwd("~/Desktop") 

## You can also do this manually by clicking on the 
## R console, going to "File," and then clicking "Change dir"

## Now you can load files from your working directory

## To load an R data file:
# load("sample.RData")

##To load a data file saved as a text document:
# data <- read.table("sample.txt")

##To load a Stata file:
# data <- read.dta("sample.dta")

##To load a csv:
# data <- read.csv(file="sample.csv",head=TRUE,sep=",")

########################################################
## R packages
########################################################

## People write software for R all the time.
## These are called "packages" 

install.packages("rbounds")
## this must be done only once 

library(rbounds)
## this must be done every time you use the foreign library

## Other useful packages are
## foreign -- permits you to load data formatted for other software
## xtables -- helps you write up tables in LaTeX code
## arm -- applied regression and multi-level modeling
## ggplot2 -- powerful and versatile graphing functions
## more packages are at http://cran.r-project.org/web/packages/

## good practice when starting a new r script

rm(list=ls())
setwd("~/Desktop") 
library(rbounds)

########################################################
## Basic R functions
########################################################

# R has many preprogrammed functions that manipulate objects.
# To use a function, you type the function name followed by the
# arguments in parentheses

a <- c(1,3,6,5,9,22)

b <- c(4,5,6,5,2,1)

sum(a) ## sum of all elements of a

sum(b)

sum(a,b)

max(a) ## maximum element in a

min(a) ## minimum element in a

mean(a) ## mean of a

unique(b) ## unique values

length(a) ## number of elements in a,
## useful for when you need to calculate the sample size, n

sort(a) ## sorts a from lowest to highest

## you can store the output of a function, as a new object.

output <- length(a)
output

## These functions are useful for creating vectors

## creates a sequence of numbers
seq_num <- seq(from = 0, to = 5, by = .5) 
seq_num

## repeats the number "10" 27 times.
rep_num = rep(10, 27) 
rep_num

# matrices
matrix1 <- cbind(a,b)
matrix2 <- rbind(a,b)
matrix1
matrix2

## You can give your matrix column names and row names
matrix1
colnames(matrix1) <- c("col1","col2")
matrix1
rownames(matrix1) <- c("r1","r2","r3","r4","r5","r6")
matrix1

## We can extract particular elements of a matrix just like
## we could from a vector, though this time we have to specify
## two dimensions, the row and the column

matrix1[1,1] # first row and column

matrix1[,1] # first column
matrix1[,2] # second column

matrix1[1,] # first row
matrix1[4,] # fourth row

########################################################
## Let's create a sample data set
########################################################

# Set seed for simulation
set.seed(123)

# First, to create a sample dataset

## let's create a few independent variables with normal distributions (mean=0, sd=1)

n <- 1000 # number of observations 
x1 <- rnorm(n, 0,1) # independent variables 1
x2 <- rnorm(n, 0,1) + 2 # independent variables 2
x3 <- rnorm(n, 0,1) - 1 # independent variables 3
y <-  rnorm(n, 0,1) + 15 # dependent variables 

n
x1

## here we are binding the variables columns into a dataset

data <- data.frame(x1,x2,x3,y)
View(data)
names(data)
head(data)

# check variables using $

data$x1
data$y

########################################################
## Getting a feel for the data in R
########################################################

## Graphically

hist(data$x1)

plot(data$x1,data$y)

## Numerically

sum(data$x1)
mean(data$x1)
sd(data$x1)
cor(data$x1,data$y)
summary(data$x1)
summary(data)

# OLS

ols <- lm(data$y ~ data$x1 + data$x2 + data$x3)
summary(ols)

plot(data$x1, data$y, pch=16)
abline(lm(data$y ~ data$x1), lwd=3, col="red")

plot(data$x2, data$y, pch=16)
abline(lm(data$y ~ data$x2), lwd=3, col="red")

plot(data$x3, data$y, pch=16)
abline(lm(data$y ~ data$x3), lwd=3, col="red")


