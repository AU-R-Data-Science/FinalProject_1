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
