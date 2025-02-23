########################
# Workshop week 4
########################

rm(list=ls())   

###########################
# Standard errors
###########################

## Compute Bell-McCaffrey Standard Errors
library(sandwich)
library(Matrix)

message1 <- paste0(
  'Bell-McCaffrey SE undefined. This happens, e.g., when a dummy regressor is 1 ',
  'for one cluster and 0 otherwise.'
)

MatSqrtInverse <- function(A) {
  ##  Compute the inverse square root of a matrix
  if (rankMatrix(A) < NROW(A)) stop(message1)
  ei <- eigen(A, symmetric = TRUE)
  d2 <- 1/sqrt(ei$values)
  ## diag(d2) is d2 x d2 identity if d2 is scalar, instead we want 1x1 matrix
  ei$vectors %*% (if (length(d2)==1) d2 else diag(d2)) %*% t(ei$vectors)
}

BMlmSE <- function(model, clustervar=NULL, ell=NULL, IK=TRUE) {
  X <- model.matrix(model)
  sum.model <- summary.lm(model)
  n <- sum(sum.model$df[1:2])
  K <- model$rank
  XXinv <- sum.model$cov.unscaled # XX^{-1}
  u <- residuals(model)
  
  df <- function(GG) {                # Compute DoF given G'*Omega*G
    sum(diag(GG))^2 / sum(GG * GG)
  }
  
  if(is.null(clustervar)) {           # no clustering
    Vhat <- vcovHC(model, type="HC2")
    Vhat.Stata <- Vhat*NA
    
    M <- diag(n)-X %*% XXinv %*% t(X)       # annihilator matrix
    GOG <- function(ell) {           # G'*Omega*G
      Xtilde <- drop(X %*% XXinv %*% ell / sqrt(diag(M)))
      crossprod(M * Xtilde)
    }
  } else {
    if(!is.factor(clustervar)) stop("'clustervar' must be a factor")
    
    ## Stata
    S <- length(levels(clustervar)) # number clusters
    uj <- apply(u*X, 2, function(x) tapply(x, clustervar, sum))
    Vhat.Stata <- S/(S-1) * (n-1)/(n-K) * sandwich(model, meat = crossprod(uj)/n)
    
    ## LZ2
    tXs <- function(s) {
      Xs <- X[clustervar==s, , drop=FALSE]
      MatSqrtInverse(diag(NROW(Xs))-Xs%*% XXinv %*% t(Xs)) %*% Xs
    }
    tX <- lapply(levels(clustervar), tXs) # list of matrices
    
    tu <- split(u, clustervar)
    tutX <- sapply(seq_along(tu),function(i) crossprod(tu[[i]],tX[[i]]))
    Vhat <- sandwich(model, meat = tcrossprod(tutX)/n)
    
    ## DOF adjustment
    tHs <- function(s) {
      Xs <- X[clustervar==s, , drop=FALSE]
      index <- which(clustervar==s)
      ss <- outer(rep(0,n),index)     # n x ns matrix of 0
      ss[cbind(index,1:length(index))] <- 1
      ss-X %*% XXinv %*% t(Xs)
    }
    tH <- lapply(levels(clustervar), tHs) # list of matrices
    
    Moulton <- function() {
      ## Moulton estimates
      ns <- tapply(u, clustervar,length)
      ssr <- sum(u^2)
      rho <- max((sum(sapply(seq_along(tu), function(i)
        sum(tu[[i]] %o% tu[[i]])))-ssr) / (sum(ns^2)-n), 0)
      c(sig.eps=max(ssr/n - rho, 0), rho=rho)
    }
    
    GOG <- function(ell) {
      G <- sapply(seq_along(tX),
                  function(i)  tH[[i]] %*% tX[[i]] %*% XXinv %*% ell)
      GG <- crossprod(G)
      
      if (IK==TRUE) {            # IK method
        Gsums <- apply(G, 2, function(x) tapply(x, clustervar, sum)) # Z'*G
        GG <- Moulton()[1]*GG + Moulton()[2]*crossprod(Gsums)
      }
      GG
    }
  }
  
  if (!is.null(ell)) {
    se <- drop(sqrt(crossprod(ell,Vhat) %*% ell))
    dof <- df(GOG(ell))
    se.Stata <- drop(sqrt(crossprod(ell,Vhat.Stata) %*% ell))
  } else {
    se <- sqrt(diag(Vhat))
    dof <- sapply(seq(K), function(k) df(GOG(diag(K)[,k])))
    se.Stata <- sqrt(diag(Vhat.Stata))
  }
  names(dof) <- names(se)
  
  return(list(vcov=Vhat, dof=dof, adj.se=se*qt(0.975,df=dof)/qnorm(0.975),
              se=se,
              se.Stata=se.Stata))
}

# Generate fake data (small sample)

set.seed(1234567)
treatment <- c(rep(1, 10), rep(0, 5))
covariate <- rnorm(15)
outcome <- treatment + rnorm(15)

# Use lm() to fit an OLS regression

ols.fit <- lm(outcome ~ treatment + covariate, singular.ok = FALSE)

# Use BMlmSE() to compute the BM SE and degrees of freedom

## The argument ell must have the same length as the vector of regression coefficients.
## In this example, ell = c(0, 1, 0) indicates that we only want to compute
## the Bell-McCaffrey SE and degrees of freedom for the linear combination
##   0 * intercept + 1 * (coef on treatment) + 0 * (coef on covariate),
## i.e., the coefficient on treatment.

bm <- BMlmSE(model = ols.fit, ell = c(0, 1, 0))

# Construct a 95% confidence interval for the average treatment effect

point.estimate <- coef(ols.fit)['treatment']

critical.value <- qt(0.975, df = bm$dof)  # Critical value for 95% CI
margin.of.error <- critical.value * bm$se

ci = c(point.estimate - margin.of.error, point.estimate + margin.of.error)
names(ci) = c('lower','upper')

point.estimate  # Estimated average treatment effect
bm$se           # HC2 robust SE
bm$dof          # Bell-McCaffrey degrees of freedom
ci              # Bell-McCaffrey confidence interval

# Standard regression
summary(ols.fit)
confint(ols.fit, "treatment", level=.95)

# Generate fake data (large sample)

set.seed(1234567)
treatment <- c(rep(1, 1000), rep(0, 500))
covariate <- rnorm(1500)
outcome <- treatment + rnorm(1500)

# Use lm() to fit an OLS regression

ols.fit <- lm(outcome ~ treatment + covariate, singular.ok = FALSE)

# Use BMlmSE() to compute the BM SE and degrees of freedom

## The argument ell must have the same length as the vector of regression coefficients.
## In this example, ell = c(0, 1, 0) indicates that we only want to compute
## the Bell-McCaffrey SE and degrees of freedom for the linear combination
##   0 * intercept + 1 * (coef on treatment) + 0 * (coef on covariate),
## i.e., the coefficient on treatment.

bm <- BMlmSE(model = ols.fit, ell = c(0, 1, 0))

# Construct a 95% confidence interval for the average treatment effect

point.estimate <- coef(ols.fit)['treatment']

critical.value <- qt(0.975, df = bm$dof)  # Critical value for 95% CI
margin.of.error <- critical.value * bm$se

ci = c(point.estimate - margin.of.error, point.estimate + margin.of.error)
names(ci) = c('lower','upper')

point.estimate  # Estimated average treatment effect
bm$se           # HC2 robust SE
bm$dof          # Bell-McCaffrey degrees of freedom
ci              # Bell-McCaffrey confidence interval

# Standard regression
summary(ols.fit)
confint(ols.fit, "treatment", level=.95)

################################
# Using covariates in analysis
################################

suppressMessages({
  library(randomizr)
  library(sandwich)
  library(lmtest)
})

N <- 200

# Make some covariates
X1 <- rnorm(N)
X2 <- rbinom(N, size = 1, prob = 0.5)

# Make some potential outcomes
Y0 <- .6*X1 + 3*X2 + rnorm(N)
Y1 <- Y0 + .4

# Conduct a random assignment and reveal outcomes
Z <- complete_ra(N, m= 100)
Y_obs <- Y1*Z + Y0*(1-Z)

# Mean-center the covariates
X1_c <- X1 - mean(X1)
X2_c <- X2 - mean(X2)

# Conduct Estimation
fit_adj <- lm(Y_obs ~ Z + Z*(X1_c + X2_c), singular.ok = FALSE)
# regress outcome and it has interaction
# the only thing that has a causal effect is Z, other variables no need to interpret 

# Robust Standard Errors
coeftest(fit_adj, vcov = vcovHC(fit_adj, type = "HC2"))

# Compare to unadjusted model - this is the ideal approach 
fit_unadj <- lm(Y_obs ~ Z, singular.ok = FALSE)
coeftest(fit_unadj, vcov = vcovHC(fit_unadj, type = "HC2"))

# Compare to adjusted model but not centered
fit_adj2 <- lm(Y_obs ~ Z + X1 + X2, singular.ok = FALSE)
coeftest(fit_adj2, vcov = vcovHC(fit_adj2, type = "HC2"))

#########################
# Missing Covariates
#########################

# Some covariate values are missing:
X1_obs <- X1
X2_obs <- X2

X1_obs[sample(1:N, size = 10)] <- NA
X2_obs[sample(1:N, size = 50)] <- NA

# Less than 10% of X1_obs is missing, so:
X1_obs[is.na(X1_obs)] <- mean(X1_obs, na.rm = TRUE)

# More than 10% of X2_obs is missing, so:
X2_missing <- is.na(X2_obs)
X2_obs[X2_missing] <- 0

# Mean-center the covariates
X1_obs_c <- X1_obs - mean(X1_obs)
X2_obs_c <- X2_obs - mean(X2_obs)
X2_missing_c <- X2_missing - mean(X2_missing)

# Conduct Estimation
fit_adj <- lm(Y_obs ~ Z + Z*(X1_obs_c + X2_obs_c + X2_missing_c), singular.ok = FALSE)

# Robust Standard Errors
coeftest(fit_adj, vcov = vcovHC(fit_adj, type = "HC2"))

########################
# Covariate Balance
########################

suppressMessages({
  library(randomizr)
  library(sandwich)
})

# Generate Covariates

set.seed(1234567)

N <- 1000

gender <- sample(c("M", "F"), N, replace=TRUE)
age <- sample(18:65, N, replace = TRUE)
lincome <- rnorm(N, 10, 3)
party <- sample(c("D", "R", "I"), N, prob=c(.45, .35,.2), replace=TRUE)
education <- sample(10:20, N, replace=TRUE)

# Conduct Random Assignment
Z <- complete_ra(N, 500)

# Regress treatment on covariates
fit <- lm(Z ~ gender + age + lincome + party + education, singular.ok = FALSE)
#check covariate balance, we are regressing treatment on the covariate 
#regress Z on each of the covariate

# Obtain observed heteroskedasticity-robust Wald statistic
# See Wooldridge (2010), p. 62
# Null hypothesis is that the slope coefficients are all zero, i.e.
#  R beta = 0
#  where beta is the 7 x 1 vector of coefficients, including the intercept
#  and R is the 6 x 7 matrix with all elements zero except
#   R[1,2] = R[2,3] = R[3,4] = R[4,5] = R[5,6] = R[6,7] = 1

Rbeta.hat <- coef(fit)[-1]
RVR <- vcovHC(fit, type <- 'HC0')[-1,-1]
W_obs <- as.numeric(Rbeta.hat %*% solve(RVR, Rbeta.hat))  # Wooldridge, equation (4.13)

# Compare to permutation distribution of W
#randomize the treatment again in different ways 
#construct multiple wall test to compare to our actual test and this will give us a P-value, we want a large p-value 
#compare the observed to the simulated one 
# we get a wall test, and compare it to a simulated wall test for permutation
sims <- 10000
W_sims <- numeric(sims)

for(i in 1:sims){
  Z_sim <- complete_ra(N, 500)
  fit_sim <- lm(Z_sim ~ gender + age + lincome + party + education, singular.ok = FALSE)
  
  Rbeta.hat <- coef(fit_sim)[-1]
  RVR <- vcovHC(fit_sim, type <- 'HC0')[-1,-1]
  W_sims[i] <- as.numeric(Rbeta.hat %*% solve(RVR, Rbeta.hat))
}

# Obtain p-value
p <- mean(W_sims >= W_obs)
p

hist(W_sims)
abline(v = W_obs, lwd=3, col="red")

########################
# Multiple comparisons
########################

set.seed(343)

# Generate 50 test statistics
# Half are drawn from a normal with mean 0
# The other half are drawn from a normal with mean 3
x <- rnorm(50, mean = c(rep(0, 25), rep(3, 25)))

# Obtain 50 p-values
p <- round(2*pnorm(sort(-abs(x))), 3)

# Choose alpha level
alpha <- 0.05

# Without any corrections
sig <- p < alpha
table(sig)

# Conduct three corrections
# and compare to target alpha
bonferroni_sig <- p.adjust(p, "bonferroni") < alpha
table(bonferroni_sig)
holm_sig <- p.adjust(p, "holm") < alpha
table(holm_sig)
BH_sig <- p.adjust(p, "BH") <alpha
table(BH_sig)

