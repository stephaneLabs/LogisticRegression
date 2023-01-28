#' Tools for Logit regression
#'
#' Copyright Stéphane Lassalvy 2023
#' 
#' Licence GPL-3
#' Disclaimer of Warranty.
#'
#' THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
#' EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
#' PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
#' INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
#' FOR A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE 
#' PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
#' NECESSARY SERVICING, REPAIR OR CORRECTION.
#' 
#' 
#' 
library(ggplot2)
library(ggpubr)
library(DescTools)
library(ResourceSelection)
library(moments)
library(arules)
library(car)

#' Function : quantitativeVariableDescription
#' Histogram of a quantitative variable, with skewness and kurtosis information
#' Inputs :
#' @param df data frame including the variable of interest
#' @param variableName name of the variable of interest
#' 
#' Outputs :
#' Descriptive graphics for the variable of interest
#' 
quantitativeVariableDescription <- function(df, variableName){
  
  is_outlier <- function(x) {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
  }
  
  with(df, {
    variable <- eval(parse(text = paste("df$", variableName, sep="")));
    cat("Lenght of the variable :\n")
    cat(length(variable))
    cat("\n")
    cat("Variable summary :\n")
    print(summary(variable))
    
    df_temp <- data.frame(variable=variable, outlier = is_outlier(variable))
    df_temp <- subset(df_temp, !is.na(variable))
    
    figureHistogram <- ggplot(df, aes(x = variable)) + 
                              geom_histogram(aes(y = after_stat(density)), colour = 1, fill = "lightblue") +
                              geom_density(color="darkgreen") +
                              xlab(paste("Values, Skewness = ", round(skewness(variable),2), "Kurstosis = ", round(kurtosis(variable), 2))) +
                              ylab("Prob.") +
                              ggtitle(paste("Histogram of ", toupper(variableName), sep = "")) +
                              geom_vline(xintercept = mean(variable), color = "blue", linewidth = 1) +
                              geom_vline(xintercept = median(variable), color = "darkgreen", linewidth = 1)  +
                              stat_function(fun = dnorm, args = list(mean = mean(variable), sd = sd(variable)), col = "blue", linewidth = 1)

    figureBoxplot <- ggplot(df_temp, aes(x=eval(variableName), y=variable)) + 
                            geom_boxplot(outlier.size=2, outlier.colour="black") + 
                            xlab(toupper(variableName)) + ylab(paste("Values of ", toupper(variableName), sep="")) +
                            geom_text(data = subset(df_temp, outlier == TRUE), aes(label = variable), na.rm = TRUE, hjust = -0.3) + 
                            ggtitle(paste("Boxplot and outliers of ", toupper(variableName), sep=""))
    
    figure <- ggarrange(figureHistogram, figureBoxplot, nrow = 1, ncol = 2)
    figure
  })
}


#' Function : logOddsVsQuantitativePredictor
#' Relationship between the log odds and the quantitative predictor
#' Inputs :
#' @param df data frame including the variable of interest
#' @param binaryResponse name of the binary response"
#' @param quantitativePredictor name of the quantitative predictor
#' @param method method to discretize the quantitative predictor ("round" rounds the quantitative variable to "digits" digits,
#'        "interval" (equal interval width), "frequency" (equal frequency), "cluster" (k-means clustering) 
#'         and "fixed" (categories specifies interval boundaries).
#'         The function then use the rounded variable ("round" option) or compute the mean for each level of the discretization.
#' @param digits number of digits for the quantitative predictor
#' @param breaks option for the discretize function from package arules
#' @param variableName name of the variable of interest
#' 
#' Outputs :
#' Descriptive graphics for the variable of interest
#' 
logOddsVsQuantitativePredictor <- function(df, binaryResponse = "dm", quantitativePredictor, method = "round", digits = 2, breaks = 3){
  # Get the variables of interest
  binaryResp     <- as.factor(eval(parse(text=paste("df$", binaryResponse, sep = ""))))
  xQuantitative  <- as.numeric(eval(parse(text=paste("df$", quantitativePredictor, sep = ""))))
  
  # Discretize the quantitative predictor to observe the relationship between means of predictor values and log(odds) of the response
  if(method=="round"){
    quantPredictor  <- as.factor(round(xQuantitative, digits = digits))
  } else {
    xDiscretization <- as.factor(as.character(discretize(xQuantitative, method = method, breaks = breaks)))
    ddf <- data.frame(xQuantitative = xQuantitative, xDiscretization  = xDiscretization)
    quantPredictor <- aggregate(xQuantitative ~ xDiscretization, data = ddf, FUN = "mean", na.action = "na.omit")
    quantPredictor <- rename(quantPredictor, Xmeans = xQuantitative)
    ddf            <- full_join(ddf, quantPredictor, by = "xDiscretization")
    quantPredictor <- as.numeric(ddf$Xmeans)
  }
  if(any(binaryResp == "")){binaryResp[binaryResp == ""] <- NA}
  if(any(quantPredictor == "")){quantPredictor[quantPredictor == ""] <- NA}

  # Compute the odds for each value of predictor
  frequencyTable <- prop.table(table(quantPredictor , binaryResp), margin = 1)
  odds <- frequencyTable[, "yes"] / frequencyTable[, "no"]

  # Make a table with log(odds) and predictor
  oddsTable <- data.frame(predictor = as.numeric(row.names(frequencyTable)),
                          yes = frequencyTable[, "yes"],
                          no = frequencyTable[, "no"],
                          odds = odds,
                          logodds = log(odds))

  
  oddsTable <- oddsTable[oddsTable$logodds != -Inf & oddsTable$logodds != Inf, ]
  
  # regression and correlations between log(odds) and predictor
  coefOddsModel <- coef(lm(logodds ~ predictor, data = oddsTable))
  correlationPredictorVsLogOdds <- with(oddsTable, round(cor(predictor, odds, use = "complete.obs"), 3))

  # print relationship between log(odds) and predictor
  if(method=="round"){
    with(oddsTable, {figure <- ggplot(data = oddsTable, aes(x = predictor, y = logodds)) +
      geom_point() +
      geom_abline(intercept = coefOddsModel[1], slope = coefOddsModel[2], linetype = 2, color = "red") +
      xlab(paste("Predictor : ", toupper(quantitativePredictor), ", digits = ", digits)) +
      ylab("Log(Odds)") +
      ggtitle(paste("Predictor :", toupper(quantitativePredictor)," vs Log(Odds),\n corr = ", correlationPredictorVsLogOdds))
    figure
    })
  } else {
    with(oddsTable, {figure <- ggplot(data = oddsTable, aes(x = predictor, y = logodds)) +
                               geom_point() +
                               geom_abline(intercept = coefOddsModel[1], slope = coefOddsModel[2], linetype = 2, color = "red") +
                               xlab(paste("Predictor : ", toupper(quantitativePredictor), ", means of ", nrow(oddsTable), " groups")) +
                               ylab("Log(Odds)") +
                               ggtitle(paste("Predictor :", toupper(quantitativePredictor)," vs Log(Odds),\n corr = ", correlationPredictorVsLogOdds))
                     figure
    })
  }
}


#' Function : binaryResponseVSCategoricalPredictor
#' Percentage table and conditionnal percentage tables crossing the binary response with a categorical variable
#' Inputs :
#' @param binaryResponse : name of the binary response variable
#' @param categoricalPredictor : name of the categorical predictor
#' @param digits : number of digits in the text output
#' 
#' Output :
#' Only text outputs
#' 
binaryResponseVSCategoricalPredictor <- function(df, binaryResponse = "dm", categoricalPredictor, digits = 2){
  # Cross table of the two categorical variables
  binaryResp <- as.character(eval(parse(text=paste("df$", binaryResponse, sep = ""))))
  categoPred <- as.character(eval(parse(text=paste("df$", categoricalPredictor, sep = ""))))
  
  binaryResp[binaryResp == ""] <- NA
  categoPred[categoPred == ""] <- NA
  
  binaryResp <- as.factor(binaryResp)
  categoPred <- as.factor(categoPred)
  
  # Discard lines with missing data
  df <- na.omit(data.frame(binaryResp, categoPred))

  responseVsPredictorTable <- with(df, table(binaryResp, categoPred))
  
  # Proportion of observation in each cell of the table
  print(round(addmargins(100*prop.table(responseVsPredictorTable)), digits = digits))
  
  # Chi-square test of independance based on the observed frequencies
  print(chisq.test(responseVsPredictorTable))
  
  # Conditional proportions 
  print(round(addmargins(100*prop.table(responseVsPredictorTable, margin = 1), margin = 2), digits = digits))
  print(round(addmargins(100*prop.table(responseVsPredictorTable, margin = 2), margin = 1), digits = digits))
}

#' plotVIF
#' Plot Variable Inflation Factor for variables in a linear of logit model to assess multicolinearity among the predictors
#' Inputs :
#' @param model the considered model
#' 
#' Outputs :
#' Diagram of the VIF value withe 3 thresholds as a rule of thumb : 
#' 2.5 warning, 5 some concerns, 10 serious problem
#' 
plotVIF <- function(model){
  vifDf <- as.data.frame(vif(model))
  vifDf <- add_rownames(vifDf, var = "rowName")
  print(vifDf)
  barplotFigure <- ggplot(data = vifDf, aes(x = rowName, y = GVIF)) +
                   geom_col(color ="blue", fill = "lightblue") +
                   xlab("Variable/Factor") +
                   ylab("GVIF value") +
                   ggtitle("VIF Values for the considered model") +
                   geom_hline(yintercept = c(2.5, 5, 10), color = c("darkgreen", "blue", "red"))
  barplotFigure
}


#' Function mcFaddenR2
#' Inputs :
#' @param model : current model M
#' @param nullmondel : null model M0
#' 
#' Output :
#' Mc Fadden R2
#'  
mcFaddenR2 <- function(model, nullModel){
  R2 <- as.numeric(1 - logLik(model) / logLik(nullModel))
  return(R2)
}


#' Function adjMcFaddenR2
#' Inputs :
#' @param model : current model M
#' @param nullmondel : null model M0
#' 
#' Output :
#' R2 of Mc Fadden value 
#' 
adjMcFaddenR2 <- function(model, nullModel){
  R2 <-  as.numeric(1 - (logLik(model) - length(coef(model))) / (logLik(nullModel) - 1))
  return(R2)
}


#' Function sasR2
#' Inputs :
#' @param model : current model M
#' @param nullmondel : null model M0
#' 
#' Output :
#' SAS R2
#' 
sasR2 <- function(model, nullModel){n  <- nrow(na.omit(model$data))
                                    R2 <- as.numeric(1 - exp(-2 * (logLik(model) - logLik(nullModel)) / n))
                                    return(R2)
}

#' Function adjSasR2
#' Inputs :
#' @param model : current model M
#' @param nullmondel : null model M0
#' 
#' Output :
#' Ajusted SAS R2
#' 
adjSasR2 <- function(model, nullModel){
                                       n  <- nrow(na.omit(model$data))
                                       R2 <- as.numeric(1 - exp(-2 * (logLik(model) - logLik(nullModel)) / n))
                                       adjR2 <- as.numeric(R2 / (1 - exp(2 * logLik(nullModel) / n)))
                                       return(adjR2)
}

#' Function devR2
#' Inputs :
#' @param model : current model M
#' @param nullmondel : null model M0
#' 
#' Output :
#' dev R2
#' 
devR2 <- function(model, nullModel){
  n <- nrow(na.omit(model$data))
  nullDeviance  <- model$null.deviance
  modelDeviance <- model$deviance
  deltaDeviance <- nullDeviance - modelDeviance
  # it's strange that 2*(logL(M) - logL(M0)) does not give the same delta-deviance
  R2 <- deltaDeviance / nullDeviance
  return(R2)
}


#' Function HoslemTest
#' Performs the Hoslem Lemeshow test and the related charts
#' Inputs :
#' @param model   : the current model M
#' @param ngroups : number of groups to perform de test
#' 
#' Outputs :
#' Print the results of the test and displays the related graphics 
#' 
HoslemTest <- function(model, ngroups = 10){
  # Compute and display Hoslem Ledeshow test results
  hoslemTest <- hoslem.test(x = model$y, y = fitted(model), g = ngroups)
  print(hoslemTest)

  df <- as.tibble(cbind(hoslemTest$observed, hoslemTest$expected))
  
  # Display the related graphics
  # Observed cases vs Predicted cases for each of the 10 groups
  observed1VsPredicted1 <- ggplot(data = df, aes(x = y1, y = yhat1))  +
                           geom_point() +
                           geom_abline(intercept = 0, slope = 1, color="red") +
                           geom_smooth(method =lm, formula = y ~ x, color = "blue") +
                           xlab("Observed cases") + 
                           ylab("Predicted cases") +
                           ggtitle("Obs. cases Vs Pred. cases")

  # print(observed1VsPredicted1)
  observed0VsPredicted0 <- ggplot(data = df, aes(x = y0, y = yhat0))  +
                           geom_point() +
                           geom_abline(intercept = 0, slope = 1, color="red") +
                           geom_smooth(method =lm, formula = y ~ x, color = "blue") +
                           xlab("Observed cases") + 
                           ylab("Predicted cases") +
                           ggtitle("Obs. non cases Vs Pred. non cases")
  
  df <- df %>% mutate(total = y1 + y0, 
                      totalexpected = yhat1 + yhat0,
                      y1pct = y1 / total * 100,
                      yhat1pct= yhat1 / totalexpected * 100)
  
  # Observed cases vs Predicted cases for each of the 10 groups in %
  observed1VsPredicted1pct <- ggplot(data = df, aes(x = y1pct, y = yhat1pct))  +
                              geom_point() +
                              geom_abline(intercept = 0, slope = 1, color="red") +
                              geom_smooth(method =lm, formula = y ~ x, color = "blue") +
                              xlab("Observed cases") + 
                              ylab("Predicted cases") +
                              ggtitle("Obs. cases Vs Pred. cases in %")
  
  outputGraph <- ggarrange(observed1VsPredicted1, observed0VsPredicted0, observed1VsPredicted1, nrow = 1, ncol = 3)
  print(outputGraph)
  }

#' Residuals plots
#' 
residualPlots <- function(model, binaryresponse = "DiabeteDiagostic", quantitativePredictor = "none", label.size = 3){
 df <- model$data
 df$response <- as.factor(eval(parse(text = paste("df$", binaryresponse))))
 
 df$cooksDistances    <- cooks.distance(model)
 df$cooksIndex        <- seq(along.with = as.numeric(df$cooksDistances))
 
 df$leverages        <- influence(model)$hat
 df$indexLeverages   <- seq(along.with = as.numeric(df$leverages))
 
 df$pearsonResiduals  <- residuals(model, type = "pearson")
 df$pearsonIndex      <- seq(along.with = as.numeric(df$pearsonResiduals))
 
 df$devianceResiduals <- residuals(model, type = "deviance")
 df$devianceIndex     <- seq(along.with = as.numeric(df$devianceResiduals))

 df$devianceResiduals <- residuals(model, type = "deviance")
 df$devianceIndex     <- seq(along.with = as.numeric(df$devianceResiduals))

 df <- df %>% mutate(stdPearsonResiduals  = pearsonResiduals  / sqrt(1 - leverages),
                     stdDevianceResiduals = devianceResiduals / sqrt(1 - leverages))
 
 if(quantitativePredictor == "none"){
   # Cook's distances figure
   cooksTreshold <- 4 / nrow(df)
   CooksFigure  <- ggplot(data = df, aes(ymax = cooksDistances, ymin = 0, x=cooksIndex, color = response)) +
                          geom_linerange() +
                          geom_text(data = subset(df, cooksDistances > cooksTreshold),
                                    aes(y = cooksDistances, x=cooksIndex, label = cooksIndex), 
                                    hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                          geom_hline(yintercept = cooksTreshold, color = "red") +
                          xlab("Observation index") +
                          ylab("Cook's distances") +
                          ggtitle("Cook's distances") +
                          scale_colour_discrete(name = "Response")
   
   # Pearson residuals figure
   pearsonResFigure  <- ggplot(data = df, aes(y = pearsonResiduals, x=pearsonIndex, color = response)) +
                        geom_point() +
                        geom_text(data = subset(df, cooksDistances > cooksTreshold), 
                                  aes(y = pearsonResiduals, x=pearsonIndex, label = pearsonIndex),
                                  hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                        geom_hline(yintercept = 0, color = "black") +
                        xlab("Observation index") +
                        ylab("Pearson residuals") +
                        ggtitle("Pearson residuals") +
                        scale_colour_discrete(name = "Response") 
   
   # Pearson standardized residuals figure
   stdPearsonResFigure  <- ggplot(data = df, aes(y = stdPearsonResiduals, x=pearsonIndex, color = response)) +
                           geom_point() +
                           geom_text(data = subset(df, cooksDistances > cooksTreshold), 
                                     aes(y = stdPearsonResiduals, x=pearsonIndex, label = pearsonIndex), 
                                     hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                           geom_hline(yintercept = 0, color = "black") +
                           xlab("Observation index") +
                           ylab("Standardized Pearson residuals") +
                           ggtitle("Standardized Pearson residuals") +
                           scale_colour_discrete(name = "Response")
   
   # Leverages figures
   nParameters <- sum(df$leverages)
   leveragesTreshold <- 3 * nParameters / nrow(df)
   leveragesFigure   <- ggplot(data = df, aes(ymax = leverages, ymin = 0, x=indexLeverages, color = response)) +
                        geom_linerange() +
                        geom_hline(yintercept = leveragesTreshold, color = "red") +
                        geom_text(data = subset(df, leverages > leveragesTreshold), 
                                  aes(y = leverages, x=indexLeverages, label = indexLeverages), 
                                  hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                        xlab("Observation index") +
                        ylab("Leverages") +
                        ggtitle("Leverages")  +
                        scale_colour_discrete(name = "Response")
   
     
   # Deviance residuals figure
   devianceResFigure <- ggplot(data = df, aes(y = devianceResiduals, x=devianceIndex, color = response)) +
                        geom_point() +
                        geom_text(data = subset(df, cooksDistances > cooksTreshold), 
                                  aes(y = devianceResiduals, x=devianceIndex, label = devianceIndex), 
                                  hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                        geom_hline(yintercept = 0, color = "black") +
                        xlab("Observation index") +
                        ylab("Deviance residuals") +
                        ggtitle("Deviance residuals") +
                        scale_colour_discrete(name = "Response")       
   
   # Deviance standardized residuals figure
   stdDevianceResFigure <- ggplot(data = df, aes(y = stdDevianceResiduals, x=devianceIndex, color = response)) +
                           geom_point() +
                           geom_text(data = subset(df, cooksDistances > cooksTreshold), 
                                     aes(y = stdDevianceResiduals, x=devianceIndex, label = devianceIndex), 
                                     hjust = -0.1, vjust = 0.1, size = label.size, check_overlap = TRUE, color = "black") +
                           geom_hline(yintercept = 0, color = "black") +
                           xlab("Observation index") +
                           ylab("Standardized Deviance residuals") +
                           ggtitle("Standardized Deviance residuals")  +
                           scale_colour_discrete(name = "Response")
   
   figure <- ggarrange(CooksFigure, pearsonResFigure, stdPearsonResFigure, leveragesFigure, devianceResFigure, stdDevianceResFigure, nrow = 2, ncol = 3)
   figure
 }
}
