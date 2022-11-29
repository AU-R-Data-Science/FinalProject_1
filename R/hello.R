# Logistic Regression!
#
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


leastSquares <- function(y, X) {
  solve(t(X) %*% X) %*% t(X) %*% y
}

logistic_obj <- function(beta, y, X) {
  p <- 1/(1+exp(-X%*%beta))
  sum(-y * log(p) - (1 - y)*log(1-p))
}

beta_hat <- function(y, X) {

  beta_est <- optim(leastSquares(y, X), logistic_obj, y = y, X = X)$par
  output <- list("beta_hat" = beta_est, "response" = y, "predictors" = X)

  class(output) = "logistic_optim"

  return(output)

}


bootStrapLogisticRegression <- function(B = 20, alpha, y, X) {
  boot_mean <- list()
  for (x in 1:ncol(X)) {
    boot_mean[[x]] <- rep(NA, B)

  }

  for (i in 1:B){
    toSample <- (cbind(y, X))
    toSample <- toSample[sample(1:nrow(toSample), nrow(toSample), replace = TRUE), ]
    yBoot <- toSample[,1]
    xBoot <- toSample[,-c(1)]

    betaMat <- beta_hat(yBoot, xBoot)$beta_hat

    # Step 3
    for (o in 1:ncol(X)) {
      boot_mean[[o]][i] <- betaMat[o]
    }

  }

  for (k in 1:ncol(X)) {
    cat(colnames(X)[k], "\n")
    print(quantile(boot_mean[[k]], c(alpha/2, 1 - alpha/2)))
    cat("\n")
  }


}


plotRegressionCurves <- function(y, X) {
  par(mfrow=c(ncol(X), ncol(X)))
  for (x in 1:ncol(X)) {
    (X[,x])
    Predicted_data <- data.frame(varToPlot=seq(
      min(X[,x]), max(X[,x]),len=500))
    mat <- matrix(data = NA, nrow = 500, ncol = ncol(X))
    for (q in 1:ncol(X)) {

      if (q == x) {
        mat[, q] <- Predicted_data$varToPlot
      }
      else {
        mat[, q] <- mean(X[,q])
      }

    }
    (mat)
    Predicted_dataToGraph <- 1/(1+exp(-mat%*%output))

    finalDataFrame <- data.frame(Predicted_dataToGraph, Predicted_data$varToPlot)

    colnames(finalDataFrame) <- c('y','x')

    plot(y ~ X[,x], xlab=colnames(X)[x],ylab="y")
    lines(y ~ x, finalDataFrame, lwd=2, col="green")



    # mat <- matrix(data = c(rep(1, 500), rep(mean, 500), Predicted_data$varToPlot), nrow = 500, ncol = 3)


  }

}



#confusion matrix function

confusionMatrix <- function(y, X) {
  #X <- datatoTest$X
  #y <- datatoTest$y

  output <- beta_hat(y, X)$beta_hat
  (output)
  predBeforeRound <- 1/(1+exp(-X%*%output))
  (predBeforeRound)
  pred[predBeforeRound <= .5] <- 0
  pred[predBeforeRound > .5] <- 1
  (pred)
  TP <- 0
  TN <- 0
  FN <- 0
  FP <- 0

  for (i in 1:nrow(y)) {
    if (y[i] == 1 && pred[i] == 1) {
      TP <- TP + 1
    } else if (y[i] == 0 && pred[i] == 0) {
      TN <- TN + 1
    } else if (y[i] == 1 && pred[i] == 0) {
      FN <- FN + 1
    }
    else if (y[i] == 0 && pred[i] == 1) {
      FP <- FP + 1
    }

  }

  Prediction <- factor(c(0, 0, 1, 1))
  Target <- factor(c(0, 1, 0, 1))
  Y      <- c(TN, FP, FN, TP)
  df <- data.frame(Prediction, Target, Y)
  p <- ggplot(data =  df, mapping = aes(x = Target, y = Prediction)) +
    geom_tile(aes(fill = Y), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", Y)), vjust = 1) +
    scale_fill_gradient(low = "red", high = "green") +
    theme_bw() + theme(legend.position = "none")
  print(p)

  cat("Prevalence: ", getPrevalence(TP, TN, FN, FP), "\n")
  cat("Accuracy: ", getAccuracy(TP, TN, FN, FP), "\n")
  cat("Sensitivity: ", getSensitivity(TP, TN, FN, FP), "\n")
  cat("Specificity: ", getSpecificity(TP, TN, FN, FP), "\n")
  cat("False Discovery Rate: ", getFalseDiscoveryRate(TP, TN, FN, FP), "\n")
  cat("Diagnostic Odds Ratio: ", getDiagnosticOddsRatio(TP, TN, FN, FP), "\n")
}

getPrevalence <- function(TP, TN, FN, FP) {
  N <- FP + TN
  P <- TP + FN
  Prev <- P / P + N

  return(PT)
}

getAccuracy <- function(TP, TN, FN, FP) {
  ACC <- (TP + TN)/(TP + TN + FP + FN)
  return(ACC)
}

getSensitivity <- function(TP, TN, FN, FP) {
  TPR <- TP / (TP + FN)
  return(TPR)
}

getSpecificity <- function(TP, TN, FN, FP) {
  TNR <- TN / (TN + FP)
  return(TNR)
}

getFalseDiscoveryRate <- function(TP, TN, FN, FP) {

  FDR <- FP/(FP + TP)
  return(FDR)
}

getDiagnosticOddsRatio <- function(TP, TN, FN, FP) {

  DOR <- (TP/FN) / (FP/TN)
  return(DOR)
}



#' Gives ability to plot certain metrics from logistic regression Model
#' @description This function plots the metric given as a paremeter evaluted over a grid of cut-off values
#' for prediction going from .1 to .9 with steps of .1
#' @param y A \code{double} value of the 1 column matrix containing the response of interest.
#' @param X An \eqn{n \times p} \code{double} value of the matrix containing the values of the predictors.
#' @param metric An \code{string} value of the metric the user would like to graph, this can be ""Prevalence", "Accuracy",
#' "Sensitivity", "Specificity", "False Discovery Rate", "Diagnostic Odds Ratio"
#' @author
#' @export
plotMetric <- function(y, X, metric) {
  output <- beta_hat(y, X)$beta_hat
  
  whichMetric <- tolower(metric)
  
  if (whichMetric == "prevalence") {
    cutOffArr <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9)
    metricVals <- c()
    for (d in 1:length(cutOffArr)) {
      predBeforeRound <- 1/(1+exp(-X%*%output))
      pred[predBeforeRound <= cutOffArr[d]] <- 0
      pred[predBeforeRound > cutOffArr[d]] <- 1
      
      TP <- 0
      TN <- 0
      FN <- 0
      FP <- 0
      
      for (i in 1:nrow(y)) {
        if (y[i] == 1 && pred[i] == 1) {
          TP <- TP + 1
        } else if (y[i] == 0 && pred[i] == 0) {
          TN <- TN + 1
        } else if (y[i] == 1 && pred[i] == 0) {
          FN <- FN + 1
        }
        else if (y[i] == 0 && pred[i] == 1) {
          FP <- FP + 1
        }
        
      }
      metricVals <- c(metricVals, getPrevalence(TP, TN, FN, FP))
    }
    
    df <- data.frame(CutOff=cutOffArr,
                     Prevalence=metricVals)
    
    p<-ggplot(data=df, aes(x=CutOff, y=Prevalence)) +
      geom_bar(stat="identity")
    p
    
    
  } else if (whichMetric == "accuracy") {
    cutOffArr <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9)
    metricVals <- c()
    for (d in 1:length(cutOffArr)) {
      predBeforeRound <- 1/(1+exp(-X%*%output))
      pred[predBeforeRound <= cutOffArr[d]] <- 0
      pred[predBeforeRound > cutOffArr[d]] <- 1
      
      TP <- 0
      TN <- 0
      FN <- 0
      FP <- 0
      
      for (i in 1:nrow(y)) {
        if (y[i] == 1 && pred[i] == 1) {
          TP <- TP + 1
        } else if (y[i] == 0 && pred[i] == 0) {
          TN <- TN + 1
        } else if (y[i] == 1 && pred[i] == 0) {
          FN <- FN + 1
        }
        else if (y[i] == 0 && pred[i] == 1) {
          FP <- FP + 1
        }
        
      }
      metricVals <- c(metricVals, getAccuracy(TP, TN, FN, FP))
    }
    
    df <- data.frame(CutOff=cutOffArr,
                     Accuracy=metricVals)
    
    p<-ggplot(data=df, aes(x=CutOff, y=Accuracy)) +
      geom_bar(stat="identity")
    p
    
    
  } else if  (whichMetric == "sensitivity") {
    cutOffArr <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9)
    metricVals <- c()
    for (d in 1:length(cutOffArr)) {
      predBeforeRound <- 1/(1+exp(-X%*%output))
      pred[predBeforeRound <= cutOffArr[d]] <- 0
      pred[predBeforeRound > cutOffArr[d]] <- 1
      
      TP <- 0
      TN <- 0
      FN <- 0
      FP <- 0
      
      for (i in 1:nrow(y)) {
        if (y[i] == 1 && pred[i] == 1) {
          TP <- TP + 1
        } else if (y[i] == 0 && pred[i] == 0) {
          TN <- TN + 1
        } else if (y[i] == 1 && pred[i] == 0) {
          FN <- FN + 1
        }
        else if (y[i] == 0 && pred[i] == 1) {
          FP <- FP + 1
        }
        
      }
      metricVals <- c(metricVals, getSensitivity(TP, TN, FN, FP))
    }
    
    df <- data.frame(CutOff=cutOffArr,
                     Sensitivity=metricVals)
    
    p<-ggplot(data=df, aes(x=CutOff, y=Sensitivity)) +
      geom_bar(stat="identity")
    p
    
  } else if (whichMetric == "specificity") {
    cutOffArr <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9)
    metricVals <- c()
    for (d in 1:length(cutOffArr)) {
      predBeforeRound <- 1/(1+exp(-X%*%output))
      pred[predBeforeRound <= cutOffArr[d]] <- 0
      pred[predBeforeRound > cutOffArr[d]] <- 1
      
      TP <- 0
      TN <- 0
      FN <- 0
      FP <- 0
      
      for (i in 1:nrow(y)) {
        if (y[i] == 1 && pred[i] == 1) {
          TP <- TP + 1
        } else if (y[i] == 0 && pred[i] == 0) {
          TN <- TN + 1
        } else if (y[i] == 1 && pred[i] == 0) {
          FN <- FN + 1
        }
        else if (y[i] == 0 && pred[i] == 1) {
          FP <- FP + 1
        }
        
      }
      metricVals <- c(metricVals, getSpecificity(TP, TN, FN, FP))
    }
    
    df <- data.frame(CutOff=cutOffArr,
                     Specificity=metricVals)
    
    p<-ggplot(data=df, aes(x=CutOff, y=Specificity)) +
      geom_bar(stat="identity")
    p
    
    
  } else if (whichMetric == "false discovery rate") {
    cutOffArr <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9)
    metricVals <- c()
    for (d in 1:length(cutOffArr)) {
      predBeforeRound <- 1/(1+exp(-X%*%output))
      pred[predBeforeRound <= cutOffArr[d]] <- 0
      pred[predBeforeRound > cutOffArr[d]] <- 1
      
      TP <- 0
      TN <- 0
      FN <- 0
      FP <- 0
      
      for (i in 1:nrow(y)) {
        if (y[i] == 1 && pred[i] == 1) {
          TP <- TP + 1
        } else if (y[i] == 0 && pred[i] == 0) {
          TN <- TN + 1
        } else if (y[i] == 1 && pred[i] == 0) {
          FN <- FN + 1
        }
        else if (y[i] == 0 && pred[i] == 1) {
          FP <- FP + 1
        }
        
      }
      metricVals <- c(metricVals, getFalseDiscoveryRate(TP, TN, FN, FP))
    }
    
    df <- data.frame(CutOff=cutOffArr,
                     FalseDiscoveryRate=metricVals)
    
    p<-ggplot(data=df, aes(x=CutOff, y=FalseDiscoveryRate)) +
      geom_bar(stat="identity")
    p
    
    
  } else if (whichMetric == "Diagnostic Odds Ratio") {
    cutOffArr <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9)
    metricVals <- c()
    for (d in 1:length(cutOffArr)) {
      predBeforeRound <- 1/(1+exp(-X%*%output))
      pred[predBeforeRound <= cutOffArr[d]] <- 0
      pred[predBeforeRound > cutOffArr[d]] <- 1
      
      TP <- 0
      TN <- 0
      FN <- 0
      FP <- 0
      
      for (i in 1:nrow(y)) {
        if (y[i] == 1 && pred[i] == 1) {
          TP <- TP + 1
        } else if (y[i] == 0 && pred[i] == 0) {
          TN <- TN + 1
        } else if (y[i] == 1 && pred[i] == 0) {
          FN <- FN + 1
        }
        else if (y[i] == 0 && pred[i] == 1) {
          FP <- FP + 1
        }
        
      }
      metricVals <- c(metricVals, getDiagnosticOddsRatio(TP, TN, FN, FP))
    }
    
    df <- data.frame(CutOff=cutOffArr,
                     DiagnosticOddsRatio=metricVals)
    
    p<-ggplot(data=df, aes(x=CutOff, y=DiagnosticOddsRatio)) +
      geom_bar(stat="identity")
    p
    
    
  } else {
    print("Invalid metric.. Must be one of the following..")
  }
  
  
}



