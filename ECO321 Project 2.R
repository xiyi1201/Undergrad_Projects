# Project 2: Relationship between the Income and Education
# Date: 2021-05-14
# Group: Yihan Wang, Xiyi Lin, Yechen Li

rm(list=ls())
options(stringsAsFactors=FALSE)

setwd("/Users/Cecilia/ECO321")
data <- read.csv("Project_Filter_Data.csv", header=TRUE)
attach(data)

library(AER)
library(lmtest)
library(sandwich)
library(ggplot2)
library(car)
library(MASS)      

# create the dummy variable for sex
# if sex==1 then male=1, otherwise, male=0 (sex==2)
data$male <- ifelse(data$sex == 1, 1, 0)
# create the dummy variables for regions, excluding the fourth region 'West'
Northeast <- ifelse(data$region==1,1,0)
North <- ifelse(data$region==2,1,0)
South <- ifelse(data$region==3,1,0)

### Regression model 1: linear regression
# hh_income = β0 + β1*highest_degree + u
regression1 <- lm(hh_income ~ highest_degree, data = data)
summary(regression1)
# using heteroskedasticity robust standard errors
sum.fit.het1 <- coeftest(regression1, vcov=vcovHC(regression1, "HC1"))
sum.fit.het1

###Regression model 2: linear with more add-ons
# hh_income = β0 + β1*highest_degree + β2*male + β3*age + β4*Northeast + β5*North + β6*South + u
regression2 <- lm(hh_income ~ highest_degree + male + age + Northeast + North + South, data = data)
summary(regression2)
sum.fit.het2 <- coeftest(regression2, vcov=vcovHC(regression2, "HC1"))
sum.fit.het2

### Regression model 3: non-linear regression: polynomial
# hh_income = β0 + β1*highest_degree + β2*highest_degree^2 + β3*highest_degree^3 + u
regression3 <- lm(hh_income ~ highest_degree + I(highest_degree^2) + I(highest_degree^3), data = data)
summary(regression3)
ttest3.cubic <- coeftest(regression3, vcov=vcovHC(regression3, "HC1"))
ttest3.cubic

# since |t-value| = 1.6015 < 1.96, we cannot reject H0 that squared highest_degree has no effect on TS.

# Hypothesis Testing: F-test
# H0: population coefficients on highest_degree^2 and highest_degree^3 = 0. (linearity)
# vs. H1: at least one of these coefficients in nonzero. (non-linearity)
waldtest(regression3, c("I(highest_degree^2)", "I(highest_degree^3)"), vcov = vcovHC(regression3, "HC1"))
# the p-value = 0.237 is larger than the significance level 1%, thus we cannot reject H0.

### Regression model 4: non-linear regression: log-transformation
# log(hh_income) = β0 + β1*highest_degree + u
regression4 <- lm(I(log(hh_income)) ~ highest_degree, data = data)
subset(data, hh_income==0)
# we cannot do the log-transformation in our model for all three ways. (x/y/both)
# because the data contains the value of 0, undefined, which cannot be log().

### Regression model 5: interactions between independent variables (b/w binary and continuous)
# hh_income = β0 + β1*highest_degree + β2*male + β3*(highest_degree x male) + u
regression5 <- lm(hh_income ~ highest_degree + male + I(highest_degree*male), data = data)
summary(regression5)
ttest5.int <- coeftest(regression5,vcov = vcovHC(regression5, "HC1"))
ttest5.int
# compare with linear regression model
anova(regression1, regression5)

### Regression model 6: interaction between one binary and non-linear transformation
# hh_income = β0+β1*highest_degree+β2*highest_degree^2+β3*highest_degree^3+β4*male
# +β5*(highest_degree x male)+β6*(highest_degree^2 x male)+β7*(highest_degree^3 x male)+u
regression6 <- lm(hh_income ~ highest_degree + I(highest_degree^2) + I(highest_degree^3) + male + I(highest_degree*male) + I(highest_degree^2*male) + I(highest_degree^3*male), data = data)
summary(regression6)
ttest6.int.nl <- coeftest(regression6,vcov = vcovHC(regression6, "HC1"))
ttest6.int.nl
# F-test: the interaction coefficients to be 0
anova(regression5, regression6)
# since the F = 1.2916 with p-value = 0.2708 > 0.05, then we conclude that they are not significantly different.

### TSLS: instrumental variables
cor(regression6$residuals,data)
cor(regression6$residuals,data$male*data$highest_degree)
