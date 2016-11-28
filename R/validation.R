#'
#' @export
getMedianCurves <- function(DP, DataStorage, J=0){
  medianCurves <- validate.ICDF(DP=DP, validation=DataStorage@validation, quantiles=c(0.025,0.5,0.975), J=J)

  matrix_medianCurves <- matrix(medianCurves, nrow=3)
  matrix_medianCurves <- matrix_medianCurves[, matrix_medianCurves[2,]>0]
  
  return(matrix_medianCurves)
}

#'
#'@export
validate.Score <- function(matrix_medianCurves, data, App=F){
  mwci <- mean(matrix_medianCurves[3,]-matrix_medianCurves[1,])
  
  if(!App){
    probability <- 1 - matrix_medianCurves[2,]
    rmse <- RMSE(probability, data@validation$status)
  }else{
    rmse <- 0
  }
  score <- c(rmse, mwci)
  return(score)
}


#'
#'@export
validate.DP <- function(DP, data, App=F){
  J <- as.numeric(unique(data@validation$Sample))
  matrix_medianCurves <- getMedianCurves(DP, data, J = J)    
  log_pred <- validate.logPredictiveAccuracy(DP, data, J=J)
  
  score <- validate.Score(matrix_medianCurves, data, App)
  score <- c(score[1], log_pred, score[2])
  print(paste(class(DP), score))
  
  coverage <- validate.Coverage(matrix_medianCurves, data@validation)
  return(list(score=score, coverage=coverage))
}

#'
#'@export
validate.Coverage <- function(matrix_medianCurves, validation){
  validation$coverage <- mapply(isWithin, matrix_medianCurves[1,], matrix_medianCurves[3,], 1-validation$status)
  X <- lapply(unique(validation$Sample), function(ss) t(matrix(subset(validation, Sample == ss)$coverage)))
  coverage_matrix <- do.call(plyr::rbind.fill.matrix, X)
  coverage_array <- aperm(array(coverage_matrix, c(10,3,10)), c(1,3,2))
  coverage <- apply(coverage_array, c(2,3), mean)
  return(t(coverage))
}

#'
#' @export
validate.ICDF <- function(DP, validation, quantiles=0.5, i=0, J=0){
  params <- getParameters.ChainStorage(DP, validation, i, J)
  curves <- mapply(evaluate.ICDF, params$theta_mat, params$phi_mat,
                   params$weights_mat, params$myGrid_mat)
  curves[is.na(curves)] <- -100000000
  
  curvesReshape <- matrix(c(curves), ncol=dim(curves)[1], byrow=T)
  medianCurves <- apply(array(t(curvesReshape), c(dim(curves)[1], length(J), params$N/length(J))), c(1,2), quantile, quantiles, na.rm = T)
  return(medianCurves)
}

#'
#'@export
validate.logPredictiveAccuracy <- function(DP, DataStorage, i=0, J=0){
  params <- getParameters.ChainStorage(DP, DataStorage@presentation, i, J)
  log_pred <- mapply(evaluate.logPredictiveAccurracy, params$theta_mat, params$phi_mat,
                   params$weights_mat, params$myGrid_mat, params$myStatus_mat)
  
  log_pred_mat <- matrix(c(log_pred), ncol = dim(log_pred)[1]*length(J), byrow = T)
  if(length(is.infinite(c(log_pred_mat))) > 0) log_pred_mat[is.infinite(log_pred_mat)] <- NA
  log_pred_mean <- apply(log_pred_mat, 2, mean, na.rm=T)
  return(sum(log_pred_mean))
}

#'
#'@export
evaluate.logPredictiveAccurracy <- function(theta, phi, weights, x, status){
  
  # vectors of different length
  #status <- status[!is.na(status)]
  #x <- x[!is.na(x)]
  
  m <- length(theta)
  n <- length(x)

  log_lik <- rep(status, each=m)*dlnorm(rep(x, each=m), rep(theta, n), rep(phi, n), log = T) + 
    (1-rep(status, each=m))*plnorm(rep(x, each=m), rep(theta, n), rep(phi, n), log = T, lower.tail = F)
  
  log_lik_reshape <- matrix(log_lik, ncol=length(x))
  weighted_log_lik_reshape <- log(weights) + log_lik_reshape
  log_lik_sum <- log(colSums(exp(weighted_log_lik_reshape)))
  return(log_lik_sum)
}

#'
#' @export
RMSE <- function(probability, status){
  toRet <- (probability-status)^2
  return(sqrt(mean(toRet, na.rm=T)))
}

#'
#' @export
meanCIWidth <- function(upper, lower){
  toRet <- upper - lower
  return(mean(toRet), na.rm=T)
}

#'
#' @export
LogScore <- function(probability, status){
  toRet <- ifelse((1-status), log(1-probability), log(probability))
  return(mean(toRet, na.rm = T))
}