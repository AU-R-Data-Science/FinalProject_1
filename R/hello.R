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
