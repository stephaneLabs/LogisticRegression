#' Logit regression of diabete data
#'
#' Copyright Stéphane Lassalvy 2023
#' 
#' This R code in under GPL-3 Licence
#' 
#' Disclaimer of Warranty :
#' THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
#' EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
#' PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
#' INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#' FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.
#' SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
#' REPAIR OR CORRECTION.
#' 
#' Diabete data are available on Kaggle under licence CC0:
#' License : CC0: Public Domain. Smith, J.W., Everhart, J.E., Dickson, W.C., Knowler, W.C., & Johannes, R.S. (1988).
#' Using the ADAP learning algorithm to forecast the onset of diabetes mellitus. In Proceedings of the Symposium on
#' Computer Applications and Medical Care (pp. 261--265). IEEE Computer Society Press, and it is published t to reuse
#' in the google research dataset.
#' 
#' 
library(tidyverse)
library(FactoMineR)
library(missMDA)
library(InformationValue)
library(limma)
library(ISLR)
library(caret)

# Read the dataset
df <- read.csv("~/Documents/Info_Stats/R/Kaggle_Diabete/diabetes.csv")
df <- as_tibble(df)
df
dim(df)

# Source the R addons That I have coded for the logistic regression under R
source("~/Documents/Info_Stats/R/Kaggle_Diabete/logisticRegressionAddons_20220216.R")


# Quick format of the data
df <- df %>% mutate(Pregnancies      = Pregnancies %>% as.numeric(),
                    Glucose          = Glucose %>% as.numeric(),
                    SkinThickness    = SkinThickness %>% as.numeric(),
                    Insulin          = Insulin %>% as.numeric(),
                    BMI              = BMI %>% as.numeric(),
                    DiabetesPedigreeFunction = DiabetesPedigreeFunction %>% as.numeric(),
                    Age              = Age %>% as.numeric(),
                    Outcome = Outcome %>% as.factor())
dim(df)

# Strange paterns wich make think that 0 are NA for some variables
pairs(df[,1:8])

# Descriptive statistics on the continuous predictors
# Pregnancies
quantitativeVariableDescription(df = df, variableName = "Pregnancies")

# Glucose, seems a little O inflated and 0 is an outlier, 0 glucose in blood seems odd to me, put 0 to NA
quantitativeVariableDescription(df = df, variableName = "Glucose")
df$Glucose[df$Glucose==0] <- NA

quantitativeVariableDescription(df = df, variableName = "Glucose")


# SkinThickness, very 0 inflated, 0 for skin thickness seems odd too, put 0 to NA
quantitativeVariableDescription(df = df, variableName = "SkinThickness")
df$SkinThickness[df$SkinThickness==0] <- NA

quantitativeVariableDescription(df = df, variableName = "SkinThickness")

# Insulin, a little 0 inflated, 0 for Insulin seems a little odd to me, put 0 to NA
quantitativeVariableDescription(df = df, variableName = "Insulin")
df$Insulin[df$Insulin == 0] <- NA

quantitativeVariableDescription(df = df, variableName = "Insulin")

# BMI, 0 is very strange a BMI of 0 would imply a weight of 0 for the individual, , put 0 to NA
quantitativeVariableDescription(df = df, variableName = "BMI")
df$BMI[df$BMI==0] <- NA

quantitativeVariableDescription(df = df, variableName = "BMI")

# DiabetesPedigreeFunction, seems OK
quantitativeVariableDescription(df = df, variableName = "DiabetesPedigreeFunction")

# Age, seems OK
quantitativeVariableDescription(df = df, variableName = "Age")

# PCA of data, NA are replaced by the mean of the variable
# PCA replace missing values by the mean of the variable, some variable are highly correlated (colinearity problems ?)
df_PCA <- PCA(X = df, scale.unit = TRUE, ncp = 8, quali.sup = 9)

# 6 components explains 89% of the data
barplot(df_PCA$eig[,2])
df_PCA$eig

# Discrimination of the 2 diagnostics are not very good in the 1st plane of the PCA, maybe a linear discriminant analysis could do better
plot(df_PCA,axes = c(1,2), choix = "var")
plot(df_PCA,axes = c(1,2), choix = "ind", habillage = 9)
# 5 axes explain 82% of the inertia

# summary(df$Outcome)
# 
# estim_ncpPCA(df, ncp.max = 9, scale = FALSE, quali.sup = 9)
# Impute missing values with a better method following a PCA model
df_Imputation <- imputePCA(df, quali.sup = 9, ncp = 5, maxiter = 5000, scale = FALSE)

df_Complete <- df_Imputation$completeObs

# description between the log(odds) and the quantitative predictors (should be linear)
logOddsVsQuantitativePredictor(df = df_Complete, binaryResponse = "Outcome", method = "cluster", breaks = 10, quantitativePredictor = "Pregnancies")
logOddsVsQuantitativePredictor(df = df_Complete, binaryResponse = "Outcome", method = "cluster", breaks = 10, quantitativePredictor = "BMI")
logOddsVsQuantitativePredictor(df = df_Complete, binaryResponse = "Outcome", method = "cluster", breaks = 10, quantitativePredictor = "Glucose")
logOddsVsQuantitativePredictor(df = df_Complete, binaryResponse = "Outcome", method = "cluster", breaks = 10, quantitativePredictor = "SkinThickness")
logOddsVsQuantitativePredictor(df = df_Complete, binaryResponse = "Outcome", method = "cluster", breaks = 10, quantitativePredictor = "Insulin")

# For insulin the relationship is not so good, we choose do discretize this variable
df_Complete <- df_Complete %>% mutate(Insulin_grouped = as.factor(discretize(Insulin, method = "frequency", breaks = 3)))
attr(df_Complete$Insulin_grouped, "discretized:breaks") <- NULL
attr(df_Complete$Insulin_grouped, "discretized:method") <- NULL
binaryResponseVSCategoricalPredictor(df = df_Complete, binaryResponse = "Outcome", categoricalPredictor = "Insulin_grouped")
# seems to have a slightly influence on the outcome

# DiabetesPedigreeFunction discretized too
logOddsVsQuantitativePredictor(df = df_Complete, binaryResponse = "Outcome", digits = 1, quantitativePredictor = "DiabetesPedigreeFunction")
df_Complete <- df_Complete %>% mutate(DiabetesPedigreeFunction_grouped = as.factor(discretize(DiabetesPedigreeFunction, method = "frequency", breaks = 3)))
attr(df_Complete$DiabetesPedigreeFunction_grouped, "discretized:breaks") <- NULL
attr(df_Complete$DiabetesPedigreeFunction_grouped, "discretized:method") <- NULL
binaryResponseVSCategoricalPredictor(df = df_Complete, binaryResponse = "Outcome", categoricalPredictor = "DiabetesPedigreeFunction_grouped")

# Age discretized too
logOddsVsQuantitativePredictor(df = df_Complete, binaryResponse = "Outcome", digits = 1, quantitativePredictor = "Age")
df_Complete <- df_Complete %>% mutate(Age_grouped = as.factor(discretize(Age, method = "frequency", breaks = 3)))
attr(df_Complete$Age_grouped, "discretized:breaks") <- NULL
attr(df_Complete$Age_grouped, "discretized:method") <- NULL
binaryResponseVSCategoricalPredictor(df = df_Complete, binaryResponse = "Outcome", categoricalPredictor = "Age_grouped")

df_Complete <- droplevels(df_Complete)

#make this example reproducible
set.seed(1)

#Use 70% of dataset as training set and remaining 30% as testing set
sample_train   <- sample(c(TRUE, FALSE), nrow(df_Complete), replace = TRUE, prob = c(0.7,0.3))
df_train       <- df_Complete[sample_train, ]
df_test        <- df_Complete[!sample_train, ]

# Null model
m0 <- glm(Outcome ~ 1, family=binomial(link = logit),
          data = df_train)
AIC(m0)

# Initial model
m1 <- glm(Outcome ~ Pregnancies + Glucose + SkinThickness + Insulin_grouped + BMI + DiabetesPedigreeFunction_grouped + Age_grouped,
          family = binomial(link = logit), data = df_train)

plotVIF(m1)
summary(m1)
anova(m1, test = "Chisq")
mcFaddenR2(m1, m0)
adjMcFaddenR2(m1, m0)
sasR2(m1, m0)
adjSasR2(m1, m0)
devR2(m1, m0)
Cstat(m1)
HoslemTest(m1)
AIC(m1)
residualPlots(m1, binaryresponse = "Outcome")

# discard observation n°268
# df <- df[-268,]
# m1 <- update(m1, data = df)
# plotVIF(m1)
# summary(m1)
# anova(m1, test = "Chisq")
# mcFaddenR2(m1, m0)
# adjMcFaddenR2(m1, m0)
# sasR2(m1, m0)
# adjSasR2(m1, m0)
# devR2(m1, m0)
# Cstat(m1) 
# HoslemTest(m1)
# AIC(m1)
# residualPlots(m1)

# Backward elimination
m <- m1
m <- update(m, . ~ . - DiabetesPedigreeFunction_grouped)
summary(m)
anova(m, test = "Chisq")
mcFaddenR2(m, m0)
adjMcFaddenR2(m, m0)
sasR2(m, m0)
adjSasR2(m, m0)
devR2(m, m0)
Cstat(m)
HoslemTest(m)
AIC(m)
# the chi2 test says to discard  DiabetesPedigreeFunction_grouped but AIC disagree, Cstats is a little less than in m1, let's keep the model as is

# Interpretation des coefficients
summary(m1)
exp(coef(m1))
exp(confint(m1))

m_train <- m1

observations_train  <- df_train$Outcome
predictedProb_train <- predict(m_train, type = "response")

# Roc curves for the model fitted on the training set, then on the predictions based on the testing set
rocCurve(m_train, testDf = df_test)
