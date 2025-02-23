#################################
# Workshop #2
#################################

rm(list=ls())

# This experiment attempts to reduce prejudice.
# Control group: no intervention.
# Treatment group: workshop on discrimination.
# Outcome: discrimination scale that goes from 1 to 10. 1: low discrimination, 10: high discrimination. 
# Education: primary, secondary, higher education. 
# Age (in years).

treatment = c(0,1,0,1,1,1,0,0,0,1,0,1,0,0,1,0,1,1,1,0,0,0,1,0,1,0,0,1,0,1,1,0,1,1,1,1,0,0,0,1)
discrimination = c(10,5,6,4,3,6,8,6,3,2,3,10,9,1,1,1,5,4,6,7,6,5,4,7,8,9,3,7,6,5,3,4,10,1,2,3,2,5,6,7) 
education = c(1,1,3,2,1,1,1,2,3,2,1,2,3,2,3,1,2,1,1,1,2,3,2,1,2,3,2,2,3,3,1,1,2,3,3,1,2,2,3,1)
age = c(40,33,65,63,22,12,44,18,19,22,45,32,22,65,55,23,44,55,66,76,90,43,21,45,67,87,54,32,19,45,43,18,68,76,59,30,90,44,56,34)

# Q.1) Check if the treatment and control groups are balanced in terms of observed covariates (education and age). Use a regression, explain what is the null hypothesis and interpret the p-values

# Data frame
data <- data.frame(treatment, education, age, discrimination)

# Check balance for education
balance_education <- lm(education ~ treatment, data = data)
summary(balance_education)
# There is high p-values for treatment, p = .442, the covariates are balanced 

# Check balance for age
balance_age <- lm(age ~ treatment, data = data)
summary(balance_age)
# There is high p-values for treatment, p = .546, the covariates are balanced

#Interpretation 
#The p-values for treatment age and education are greater than 0.05 we FAIL to reject the null Hypothesis H0, meaning my covariates are balanced 

# Q.2) What is the treatment effect? Use a difference-in-means estimator. 

# Compute mean discrimination for treatment and control groups
mean_treatment <- mean(data$discrimination[data$treatment == 1])
mean_control <- mean(data$discrimination[data$treatment == 0])

# Difference-in-means estimator
treatment_effect <- mean_treatment - mean_control
treatment_effect #diff in means

#the Treatment effect is -0.55, the average reduction in discrimination due to treatment (workshop on discrimination) 

# Q.3) Compute the standard error using the formula provided in class. 

# Sample sizes for control and for the treatment
n_treatment <- sum(data$treatment == 1)
n_control <- sum(data$treatment == 0)

# Variances
var_treatment <- var(data$discrimination[data$treatment == 1])
var_control <- var(data$discrimination[data$treatment == 0])

# Standard error formula
standard_error <- sqrt(var_treatment / n_treatment + var_control / n_control)
standard_error

#the standard error is 0.833114, The standard error provides an estimate of the variability in the treatment effect across samples and our varability is 
#0.833114
#Since 0.833114 is relatively moderate, there is some level of uncertainty in the estimate
#but it is small enough for reliable inference when combined with p-values and confidence intervals.
#where is this uncertainty come from? If the treatments were randomly assigned, where could the uncertainty coming from? I can allocate the treatment in different ways 


# Q.4) Compute the confidence intervals using the formula provided in class. Interpret them. 

#lower CI 

-.55-(1.96*0.833114)
#-2.182903
#upper CI
-.55+(1.96*0.833114)
#1.082903

#this CI ranges from [1.082903, -2.182903]
#this CI interval includes 0 which does not provide us with statistically significance 

# Q.5)  Run a regression to estimate the treatment effect. Interpret the intercept, estimate, p-value, and confidence intervals.
# Regression to estimate treatment effect
reg_treatment <- lm(discrimination ~ treatment, data = data)
summary(reg_treatment)

# Compute confidence intervals
confint(reg_treatment)

#Intercept is when everything = 0, which means this is the outcome for the control group Yi(0)
#The average outcome for the control group = 5.35 
#ATE is -0.55, p = 0.513, 
#The estimated effect of the workshop on discrimination is -0.55 units, meaning that, 
#on average, participants in the treatment group show a 0.55-unit reduction in discrimination compared to the control group.
#we fail to reject the null hypothesis 
#with p-value we assume the null is TRUE, 
#However, the p-value of 0.513 suggests that this effect is not statistically significant, as it is greater than the conventional threshold of 0.05.
