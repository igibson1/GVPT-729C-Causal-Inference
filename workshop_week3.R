########################
# Workshop week 3
########################

rm(list=ls())    
library(ri) # need to go to CRAN 
library(foreign)    

##########################
# Randomization inference
##########################

titiunik <- read.dta("titiunik_data.dta")
Z <-  titiunik$term2year       # treatment is 2 year rather than 4 year term
Y <- titiunik$bills_introduced
block <- titiunik$texas0_arkansas1   # randomization occurs within each state

probs <- genprobexact(Z,blockvar=block)   # blocking is assumed when generating probability of treatment
                                          # different probabilities to be treated, you need to take this step so you do not have biased #s
table(probs)

ate <- estate(Y,Z,prob=probs)      # estimate the ATE
ate

perms <- genperms(Z,maxiter=1000,blockvar=block)   # set the number of simulated random assignments

Ys <- genouts(Y,Z,ate=0)    # create potential outcomes under the sharp null of no effect for any unit

distout <- gendist(Ys,perms,prob=probs)  # generate the sampling distribution  based on the schedule of potential outcomes implied by the null hypothesis

ate                             # estimated ATE
mean(abs(distout))        
mean(abs(distout) >= abs(ate))  # two-tailed comparison used to calculate p-value


dispdist(distout,ate)       # display p-values, 95% confidence interval, standard error under the null, and graph the sampling distribution under the null

# illustration that naive estimation of the ATE is misleading

summary(lm(Y ~ Z))

# note what happens when one includes both weights and block dummy
# notice that the p-value here matches RI quite closely
summary(lm(Y ~ Z + block))

###########################
# Kalla and Broockman 2015
###########################

#Read in overall data
data <- read.dta("kalla-broockman-donor-access-2013-data.dta", convert.underscore = TRUE)
data <- data[order(data$block),] #put data in order of blocks to make RI easier
data <- data[c(1:135,138:191,136,137),] #put block 46 at end since it only has two in it
head(data)

#Results
table(data$treat.donor, data$staffrank)

#Table 1 point estimates
round(table(data$treat.donor, data$staffrank)/rowSums(table(data$treat.donor, data$staffrank)), digits = 3)*100

#Cumulative probs
cumulative.probs <- matrix(NA, nrow = 2, ncol = 6)
for(i in 6:1){
  for(j in 1:2){
    cumulative.probs[j,7-i] <- sum(table(data$treat.donor, data$staffrank)[j,i:6])
  }
}
cumulative.probs[1,] <- cumulative.probs[1,] / cumulative.probs[1,6]
cumulative.probs[2,] <- cumulative.probs[2,] / cumulative.probs[2,6]
round(cumulative.probs, digits = 3)*100

#Generate permutations for randomization inference

Z <-  data$treat.donor  
Y <-  data$staffrank
block <- data$block  

# Outcome == 5

Y_1 = Y
Y_1[Y==5] = 1
Y_1[Y!=5] = 0
table(Y_1,Z)

probs <- genprobexact(Z,blockvar=block)   # blocking is assumed when generating probability of treatment
table(probs)

ate <- estate(Y_1,Z,prob=probs)      # estimate the ATE
ate

perms <- genperms(Z,maxiter=1000,blockvar=block)   # set the number of simulated random assignments

Ys <- genouts(Y_1,Z,ate=0)    # create potential outcomes under the sharp null of no effect for any unit

distout <- gendist(Ys,perms,prob=probs)  # generate the sampling distribution  based on the schedule of potential outcomes implied by the null hypothesis

ate                             # estimated ATE
mean(abs(distout) >= abs(ate))  # two-tailed comparison used to calculate p-value

dispdist(distout,ate) 

# Outcome == 5 or 4

Y_2 = Y
Y_2[Y==5|Y ==4] = 1
Y_2[Y<4] = 0
table(Y_2)

probs <- genprobexact(Z,blockvar=block)   # blocking is assumed when generating probability of treatment
table(probs)

ate <- estate(Y_2,Z,prob=probs)      # estimate the ATE
ate

perms <- genperms(Z,maxiter=1000,blockvar=block)   # set the number of simulated random assignments

Ys <- genouts(Y_2,Z,ate=0)    # create potential outcomes under the sharp null of no effect for any unit

distout <- gendist(Ys,perms,prob=probs)  # generate the sampling distribution  based on the schedule of potential outcomes implied by the null hypothesis

ate                             # estimated ATE
mean(abs(distout) >= abs(ate))  # two-tailed comparison used to calculate p-value

dispdist(distout,ate) 

# Outcome == 5 or 4 or 3

Y_3 = Y
Y_3[Y==5|Y ==4 |Y ==3] = 1
Y_3[Y<3] = 0
table(Y_3)

probs <- genprobexact(Z,blockvar=block)   # blocking is assumed when generating probability of treatment
table(probs)

ate <- estate(Y_3,Z,prob=probs)      # estimate the ATE
ate

perms <- genperms(Z,maxiter=100,blockvar=block)   # set the number of simulated random assignments

Ys <- genouts(Y_3,Z,ate=0)    # create potential outcomes under the sharp null of no effect for any unit

distout <- gendist(Ys,perms,prob=probs)  # generate the sampling distribution  based on the schedule of potential outcomes implied by the null hypothesis

ate                             # estimated ATE
mean(abs(distout) >= abs(ate))  # two-tailed comparison used to calculate p-value

dispdist(distout,ate) 

# Outcome == 5 or 4 or 3 or 2

Y_4 = Y
Y_4[Y==5 | Y ==4 | Y ==3 | Y ==2] = 1
Y_4[Y<2] = 0
table(Y_4)

probs <- genprobexact(Z,blockvar=block)   # blocking is assumed when generating probability of treatment
table(probs)

ate <- estate(Y_4,Z,prob=probs)      # estimate the ATE
ate

perms <- genperms(Z,maxiter=1000,blockvar=block)   # set the number of simulated random assignments

Ys <- genouts(Y_4,Z,ate=0)    # create potential outcomes under the sharp null of no effect for any unit

distout <- gendist(Ys,perms,prob=probs)  # generate the sampling distribution  based on the schedule of potential outcomes implied by the null hypothesis

ate                             # estimated ATE
mean(abs(distout) >= abs(ate))  # two-tailed comparison used to calculate p-value

dispdist(distout,ate) 

# Outcome == 5 or 4 or 3 or 2 or 1

Y_5 = Y
Y_5[Y==5 | Y ==4 | Y ==3 | Y ==2 | Y == 1] = 1
Y_5[Y<1] = 0
table(Y_5)

probs <- genprobexact(Z,blockvar=block)   # blocking is assumed when generating probability of treatment
table(probs)

ate <- estate(Y_5,Z,prob=probs)      # estimate the ATE
ate

perms <- genperms(Z,maxiter=1000,blockvar=block)   # set the number of simulated random assignments

Ys <- genouts(Y_5,Z,ate=0)    # create potential outcomes under the sharp null of no effect for any unit

distout <- gendist(Ys,perms,prob=probs)  # generate the sampling distribution  based on the schedule of potential outcomes implied by the null hypothesis

ate                             # estimated ATE
mean(abs(distout) >= abs(ate))  # two-tailed comparison used to calculate p-value

dispdist(distout,ate) 

###########################
# Exercise 
###########################

canvasser = c(0,1,1,0,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1)
older_50 = c(1,1,0,0,0,1,1,0,0,1,1,0,0,1,1,1,0,0,1,1)
turnout = c(0,0,1,0,0,1,1,0,0,0,1,1,1,0,0,0,1,1,0,0)

data_exercise <- data.frame(canvasser, older_50, turnout)

# Use an OLS to estimate the ATE and p-values 

summary(lm(turnout ~ canvasser))
#it does not account for block 

summary(lm(turnout ~ canvasser + older_50))

# Use randomization inference to estimate the ATE and p-values

# Define treatment and outcome variables
Z <- canvasser
Y <- turnout

# Generate random permutations (1000 iterations)
probs <-genprobexact(canvasser, blockvar = older_50)
table(probs)

ate <- estate(turnout, canvasser, prob = probs)
ate

#permutation matrix
perms <- genperms(canvasser, maxiter = 41580, blockvar = probs)

Ys <- genouts(turnout, canvasser, ate=0)

#distribution 
distout <- gendist(Ys, perms, prob = probs)

dispdist(distout, ate)
