library(tidyverse)
df <- read.csv("~/Documents/Info_Stats/R/Kaggle_Diabete/diabetes.csv")
df <- as_tibble(df)
df
source("~/Documents/Info_Stats/R/Kaggle_Diabete/logisticRegressionAddons_20220216.R")


df$DiabeteDiagostic <- "yes"
df$DiabeteDiagostic[df$Outcome == 0] <- "no"
dim(df)
df <- df %>% mutate(Pregnancies = Pregnancies %>% as.numeric(),
                    tPregnancies = log(1 + Pregnancies),
                    Glucose = Glucose %>% as.numeric(),
                    SkinThickness = SkinThickness %>% as.numeric(),
                    tSkinThickness  = log(1 + SkinThickness),
                    Insulin = Insulin %>% as.numeric(),
                    tInsulin = log(1 + Insulin),
                    BMI = BMI %>% as.numeric(),
                    DiabetesPedigreeFunction = DiabetesPedigreeFunction %>% as.numeric(),
                    Age = Age %>% as.numeric(),
                    DiabeteDiagostic = DiabeteDiagostic %>% as.factor())
dim(df)
quantitativeVariableDescription(df = df, variableName = "Pregnancies")
quantitativeVariableDescription(df = df, variableName = "tPregnancies")
quantitativeVariableDescription(df = df, variableName = "Glucose")
quantitativeVariableDescription(df = df, variableName = "SkinThickness")
quantitativeVariableDescription(df = df, variableName = "tSkinThickness")
quantitativeVariableDescription(df = df, variableName = "Insulin")
quantitativeVariableDescription(df = df, variableName = "tInsulin")
quantitativeVariableDescription(df = df, variableName = "BMI")

df <- df[df$BMI != 0,]
quantile(df$BMI, probs = c(0.33, 0.66))
BMI_grouped <- ifelse(df$BMI <= 28.81, "]-Inf;28.81]", ifelse(df$BMI > 28.81 & df$BMI <= 34.60, "]28.81;34.60]", ifelse(df$BMI > 34.6,"]34.60;+Inf]", NA)))
df$BMI_grouped <- as.factor(BMI_grouped)

quantitativeVariableDescription(df = df, variableName = "DiabetesPedigreeFunction")

quantitativeVariableDescription(df = df, variableName = "Age")

par(mfrow = c(3,3))
logOddsVsQuantitativePredictor(df = df, binaryResponse = "DiabeteDiagostic", method = "cluster", breaks = 10, quantitativePredictor = "Pregnancies")
logOddsVsQuantitativePredictor(df = df, binaryResponse = "DiabeteDiagostic", method = "cluster", breaks = 10, quantitativePredictor = "Glucose")
df <- df[df$Glucose != 0,]
logOddsVsQuantitativePredictor(df = df, binaryResponse = "DiabeteDiagostic", method = "cluster", breaks = 10, quantitativePredictor = "SkinThickness")
logOddsVsQuantitativePredictor(df = df, binaryResponse = "DiabeteDiagostic", method = "cluster", breaks = 10, quantitativePredictor = "Insulin")
logOddsVsQuantitativePredictor(df = df, binaryResponse = "DiabeteDiagostic", method = "round", digits = 2, quantitativePredictor = "Insulin")
summary(df$Insulin)
Insulin_grouped <- ifelse(df$Insulin <= 0, "0.00",
                          ifelse(df$Insulin > 0 & df$Insulin <= 37.00, "]0.00;37.00]",
                                 ifelse(df$Insulin > 37.00 & df$Insulin <= 130,"]37.00;130]",
                                        ifelse(df$Insulin > 130,"]130;+Inf]", NA))))
df$Insulin_grouped <- as.factor(Insulin_grouped)
df$Insulin_grouped <- relevel(df$Insulin_grouped, ref = "0.00")
binaryResponseVSCategoricalPredictor(df = df, binaryResponse = "DiabeteDiagostic", categoricalPredictor = "Insulin_grouped")

logOddsVsQuantitativePredictor(df = df, binaryResponse = "DiabeteDiagostic", method = "cluster", breaks = 10, quantitativePredictor = "BMI")
df <- df[df$BMI != 0,]
quantile(df$BMI, probs = c(0.33, 0.66))
BMI_grouped <- ifelse(df$BMI <= 28.81, "]-Inf;28.81]", ifelse(df$BMI > 28.81 & df$BMI <= 34.60, "]28.81;34.60]", ifelse(df$BMI > 34.6,"]34.60;+Inf]", NA)))
df$BMI_grouped <- as.factor(BMI_grouped)
df$BMI_grouped <- relevel(df$BMI_grouped, ref = "]-Inf;28.81]")
binaryResponseVSCategoricalPredictor(df = df, binaryResponse = "DiabeteDiagostic", categoricalPredictor = "BMI_grouped")

logOddsVsQuantitativePredictor(df = df, binaryResponse = "DiabeteDiagostic", digits = 1, quantitativePredictor = "DiabetesPedigreeFunction")

logOddsVsQuantitativePredictor(df = df, binaryResponse = "DiabeteDiagostic", digits = 1, quantitativePredictor = "Age")
quantile(df$Age, probs = c(0.33, 0.66))
Age_grouped <- ifelse(df$Age <= 25.00, "]-Inf;25.00]", ifelse(df$Age > 25.00 & df$Age <= 36.00, "]25.00;36.00]", ifelse(df$Age > 36.00,"]36.00;+Inf]", NA)))
df$Age_grouped <- as.factor(Age_grouped)
df$Age_grouped <- relevel(df$Age_grouped, ref = "]-Inf;25.00]")
binaryResponseVSCategoricalPredictor(df = df, binaryResponse = "DiabeteDiagostic", categoricalPredictor = "Age_grouped")

# Null model
m0 <- glm(Outcome ~ 1, family=binomial(link = logit), data = df)
AIC(m0)
# m1 <- glm(Outcome ~ Pregnancies + Glucose + SkinThickness + Insulin + BMI + DiabetesPedigreeFunction + Age,  data = df)
m1 <- glm(Outcome ~ Pregnancies + Glucose + SkinThickness + Insulin_grouped + BMI_grouped + DiabetesPedigreeFunction + Age_grouped,  data = df)
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
residualPlots(m1)

# Backward elimination
m <- m1
m <- update(m, . ~ . - SkinThickness)
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


m <- update(m, . ~ . - Insulin_grouped)
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


summary(m)
exp(coef(m))