#'
#' @export
fitParfm <- function(data, App=F){
  myModel <- parfm::parfm(survival::Surv(data, status) ~ 1, cluster="Sample", 
                          frailty="gamma", data=data@presentation, method="Nelder-Mead")
 
  params <- getSampleParameters.GFM(myModel, data@validation)
  matrix_medianCurves <- parfm_ci(params, data@validation)
  log_pred <- attributes(myModel$resmodel)$loglik/length(data@presentation$data)
  
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
  cov_mat <- solve(hessian)
  mySample <- exp(MASS::mvrnorm(J, log(estimates), Sigma=cov_mat))
  mySample_mat <- split(mySample, 1:nrow(mySample))
  
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
              myStatus=myStatus, mySample_mat=mySample_mat))
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

#'
#' @export
validate.logPredictiveAccuracy.GFM <- function(myModel, params){
  Mlog_lik <- mapply(Mloglikelihood, params$mySample_mat, MoreArgs = list(obs=myModel$obsdata, dist="weibull", frailty = "gamma"))
  if(length(is.infinite(c(Mlog_lik)))>0) Mlog_lik[is.infinite(Mlog_lik)] <- NA
  log_pred <- mean(-Mlog_lik, na.rm = T)
  plot(density(-Mlog_lik))
  return(log_pred)
}

#'
#' @export
weibull_gamma_fraily_PDF <- function(rho, lambda, frailty, t){
  hazard <- lambda*rho*t^(rho-1)
  toRet <- exp(-frailty*hazard*t)
  toRet
}


# functions from the parfm package
weibull <- function(pars, t, what){
  if (what == "H")
    return(pars[2] * t^(pars[1]))
  else if (what == "lh")
    return(log(pars[1]) + log(pars[2]) + ((pars[1] - 1) * log(t)))
}

fr.gamma <- function(k,
         s, 
         theta, 
         what="logLT"){
  if (what=="logLT") {
    res <- ifelse(k == 0, 
                  - 1 / theta  * log(1 + theta * s),
                  - (k + 1 / theta) * log(1 + theta * s) +
                    sum(log(1 + (seq(from=0, to=k-1, by=1) * theta))))
    return(res)
  }
  else if (what == "tau")
    return(theta / (theta + 2))
}


