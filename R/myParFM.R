#'
#' @export
parfm_survival_curves_estimate <- function(data, validation){
  myModel <- parfm::parfm(survival::Surv(data, status) ~ 1, cluster="Sample", 
                          frailty="gamma", data=data)
  
  hessian <- myModel$hessian
  myModel <- myModel$resmodel
  
  frailty_coefficients <- parfm::predict.parfm(myModel)
  estimates <- data.frame(myModel)$ESTIMATE
  SEs <- data.frame(myModel)$SE
  
  coefs <- unlist(mapply(rep, as.vector(frailty_coefficients), as.vector(table(validation$Sample))))
  ICDF <- weibull_gamma_fraily_ICDF(estimates[2], estimates[3], coefs, validation$data)
  ci <- parfm_ci(estimates, myModel, hessian, validation)
  
  toRet <- cbind(c(ci$upper_ICDF), c(ICDF), c(ci$lower_ICDF))
  toRet
}

#'
#' @export
parfm_ci <- function(estimates, myModel, hessian, validation){
  fisher_IM <- solve(hessian)
  mySample <- exp(MASS::mvrnorm(50, log(estimates), Sigma=fisher_IM))
  
  k <- as.vector(table(validation$Sample))
  frailty_coefficients_mat <- (mapply(parfm::predict.parfm, mySample[,1], MoreArgs = list(object=myModel)))
  frailty_coefficients_mat <- t(apply(frailty_coefficients_mat, 2, myRep, y=k))
  
  curves <- mapply(weibull_gamma_fraily_ICDF, mySample[, 2], mySample[, 3], frailty_coefficients_mat, 
                   MoreArgs = list(t=validation$data))
                                     
  quantiles_curves <- apply(curves, 1, quantile, c(0.025, 0.975))
  return(list(lower_ICDF=quantiles_curves[1,], upper_ICDF=quantiles_curves[2,]))
}

#'
#' @export
weibull_gamma_fraily_ICDF <- function(rho, lambda, frailty, t){
  hazard <- lambda*rho*t^(rho-1)
  toRet <- exp(-frailty*hazard*t)
  toRet
}

#'
#' @export
validation_parFM <- function(ICDF, status){
  RMSE <- sqrt(mean((1-ICDF[,2]-status)^2))
  MWCI <- mean(ICDF[,1]-ICDF[,3])
  toRet <- c(RMSE, MWCI)
  print(toRet)
  toRet
}

#'
#' @export
myRep <- function(x, y){
  unlist(mapply(rep, x, y))
}