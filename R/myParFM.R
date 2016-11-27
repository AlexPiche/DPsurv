#'
#' @export
fitParfm <- function(data, App=F){
  myModel <- parfm::parfm(survival::Surv(data, status) ~ 1, cluster="Sample", 
                          frailty="gamma", data=data@presentation, method="Nelder-Mead")
 
  params <- getSampleParameters.GFM(myModel, data@validation)
  
  matrix_medianCurves <- parfm_ci(params, data@validation)
  
  log_pred <- attributes(myModel$resmodel)$loglik
  
  score <- validate.Score(matrix_medianCurves, data, App)
  score <- c(score[1], log_pred, score[2])
  coverage <- validate.Coverage(matrix_medianCurves, data@validation)
  print(paste("GFM", score))
  return(list(score=score, coverage=coverage))
}


#'
#' @export
getSampleParameters.GFM <- function(myModel, validation, J=500){
  hessian <- myModel$hessian
  myModel <- myModel$resmodel
  
  frailty_coefficients <- parfm::predict.parfm(myModel)
  estimates <- data.frame(myModel)$ESTIMATE
  SEs <- data.frame(myModel)$SE
  
  coefs <- unlist(mapply(rep, as.vector(frailty_coefficients), as.vector(table(validation$Sample))))
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
  
  X <- lapply(unique(validation$Sample), function(ss) t(matrix(subset(validation, Sample == ss)$data)))
  X <- do.call(cbind, X)
  myStatus <- matrix(rep(X, J), nrow=J, byrow = T)
  myStatus <- split(myStatus, 1:J)
  
  
  rho_mat <- split(mySample[, 2], 1:J)
  lambda_mat <- split(mySample[, 3], 1:J)
  return(list(coefs=coefs, estimates=estimates, rho_mat=rho_mat, lambda_mat=lambda_mat, 
              frailty_coefficients_mat=frailty_coefficients_mat, myGrid=myGrid,
              myStatus=myStatus))
}

#'
#' @export
parfm_ci <- function(params, validation){
  curves <- t(mapply(weibull_gamma_fraily_ICDF, params$rho_mat, params$lambda_mat, 
                     params$frailty_coefficients_mat, params$myGrid))
  
  quantiles_curves <- apply(curves, 2, quantile, c(0.025, 0.5, 0.975), na.rm=T)
  quantiles_curves[2,] <- weibull_gamma_fraily_ICDF(params$estimates[2], params$estimates[3], params$coefs, validation$data)
  return(quantiles_curves)
  }

#'
#' @export
weibull_gamma_fraily_ICDF <- function(rho, lambda, frailty, t){
  hazard <- lambda*rho*t^(rho-1)
  toRet <- exp(-frailty*hazard*t)
  toRet
}


#   - p      : the parameters vector, in the form                              #
#              c( frailty distribution parameter(s),                           #
#                 baseline hazard parameter(s),                                #
#                 regression parameter(s) )                                    #


logLik <- function(){
  h_t <- lambda*rho*t^(rho-1)
  H_t <- lambda*t^(rho)
  loglik <- sum(as.numeric(loghaz[[1]]) + logSurv)
}

#'
#' @export
validate.logPredictiveAccuracy.GFM <- function(params){
  browser()
  Mloglikelihood(c(params$frailty_coefficients_mat[1], params$rho_mat[1], params$lambda_mat[1]),
                 list(time=params$myGrid[1], event=params$myStatus[1]),"weibull", "gamma")
  toRet <- mapply(Mloglikelihood, 
                  c(params$frailty_coefficients_mat, params$rho_mat, params$lambda_mat), 
                  list(time=params$myGrid, event=params$myStatus), 
                  MoreArgs = list("weibull", "gamma"))
  return(toRet)
}

#'
#' @export
weibull_gamma_fraily_PDF <- function(rho, lambda, frailty, t){
  hazard <- lambda*rho*t^(rho-1)
  toRet <- exp(-frailty*hazard*t)
  toRet
}

