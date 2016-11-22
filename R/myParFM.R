#'
#' @export
getMedianCurves.ParFM <- function(data, validation){
  myModel <- parfm::parfm(survival::Surv(data, status) ~ 1, cluster="Sample", 
                          frailty="gamma", data=data, method="Nelder-Mead")
  
  hessian <- myModel$hessian
  myModel <- myModel$resmodel
  
  frailty_coefficients <- parfm::predict.parfm(myModel)
  estimates <- data.frame(myModel)$ESTIMATE
  SEs <- data.frame(myModel)$SE
  
  coefs <- unlist(mapply(rep, as.vector(frailty_coefficients), as.vector(table(validation$Sample))))
  ICDF <- weibull_gamma_fraily_ICDF(estimates[2], estimates[3], coefs, validation$data)
  matrix_medianCurves <- parfm_ci(estimates, myModel, hessian, validation)
  matrix_medianCurves[2,] <- ICDF
  
  return(matrix_medianCurves)
}

#'
#' @export
parfm_ci <- function(estimates, myModel, hessian, validation, J=500){
  fisher_IM <- solve(hessian)
  mySample <- exp(MASS::mvrnorm(J, log(estimates), Sigma=fisher_IM))
  
  k <- as.vector(table(validation$Sample))
  frailty_coefficients_mat <- (mapply(parfm::predict.parfm, mySample[,1], MoreArgs = list(object=myModel)))
  frailty_coefficients_mat <- t(apply(frailty_coefficients_mat, 2, myRep, y=k))
  frailty_coefficients_mat <- split(frailty_coefficients_mat, 1:J)
  
  X <- lapply(unique(validation$Sample), function(ss) t(matrix(subset(validation, Sample == ss)$data)))
  X <- do.call(cbind, X)
  myGrid <- matrix(rep(X, J), nrow=J, byrow = T)
  myGrid <- split(myGrid, 1:J)
  
  rho_mat <- split(mySample[, 2], 1:J)
  lambda_mat <- split(mySample[, 3], 1:J)
  
  curves <- t(mapply(weibull_gamma_fraily_ICDF, rho_mat, lambda_mat, frailty_coefficients_mat, myGrid))
  
  quantiles_curves <- apply(curves, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T)
  
  return(quantiles_curves)
  }

#'
#' @export
weibull_gamma_fraily_ICDF <- function(rho, lambda, frailty, t){
  hazard <- lambda*rho*t^(rho-1)
  toRet <- exp(-frailty*hazard*t)
  toRet
}

