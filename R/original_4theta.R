#This code can be used to reproduce the forecasts submitted to the M4 competition for the 4Theta method

#Authors: E. Spiliotis and V. Assimakopoulos (2017) / Forecasting & Strategy Unit - NTUA

#Method Description: Generalizing the Theta model for automatic forecasting
#Method Type: Statistical - Decomposition

library(forecast) #requires version 8.2


FourTheta<- function(input, fh){
  #Used to automatically select the best Theta model

  #Scale
  base <- mean(input) ; input <- input/base

  molist <- c("M","A") ; trlist <- c("Lrl","Exp")

  #Check seasonality & Create list of models
  ppy <- frequency(input) ; ST <- F
  if (ppy>1){ ST <- SeasonalityTest(input, ppy) }
  if (ST==T){

    selist <- c("M","A")
    listnames <- c()
    for (i in 1:length(selist)){
      for (ii in 1:length(molist)){
        for (iii in 1:length(trlist)){
          listnames <- c(listnames,paste(selist[i], molist[ii], trlist[iii]))
        }
      }
    }

  }else{

    listnames <- c()
    for (ii in 1:length(molist)){
      for (iii in 1:length(trlist)){
        listnames <- c(listnames, paste("N", molist[ii], trlist[iii]))
      }
    }

  }

  modellist <- NULL
  for (i in 1:length(listnames)){
    modellist[length(modellist)+1] <- list(c(substr(listnames,1,1)[i], substr(listnames,3,3)[i],
                                             substr(listnames,5,7)[i]))
  }

  #Start validation
  errorsin <- c() ; models <- NULL

  #With this function determine opt theta per case
  optfun <- function(x, input, fh, curve, model, seasonality){
    mean(abs(Theta.fit(input=input, fh, theta=x, curve, model, seasonality , plot=FALSE)$fitted-input))
  }

  for (j in 1:length(listnames)){
    optTheta <- optimize(optfun, c(1:3),
                         input=input, fh=fh, curve=modellist[[j]][3], model=modellist[[j]][2],
                         seasonality=modellist[[j]][1])$minimum

    fortheta <- Theta.fit(input=input, fh=fh, theta=optTheta, curve=modellist[[j]][3], model=modellist[[j]][2],
                          seasonality=modellist[[j]][1], plot=F)
    models[length(models)+1] <- list(fortheta)
    errorsin <- c(errorsin, mean(abs(input-fortheta$fitted)))
  }

  #Select model and export
  selected.model <- models[[which.min(errorsin)]]
  description <- selected.model$model
  output <- list(fitted=selected.model$fitted*base,mean=selected.model$mean*base,
                 description=description)
  #Returns the fitted and forecasted values, as well as the model used (Type of seasonality, Type of Model, Type of Trend, Theta coef.)

  return(output)

}
