#' @useDynLib DPsurv
#' 

# utils functions call by the other files

#'
#' @export
stickBreaking <- function(beta){
  beta[length(beta)] <- 1
  pi <- sapply(1:length(beta), function(i) {beta[i]*prod(1-beta[0:(i-1)])})
  return(pi)
}

#'
#' @export
rNIG <- function(n, mu_0, n_0, v_0, vs2_0){

  phi <- rgamma(n=n, shape=0.5*v_0, rate=0.5*vs2_0)
  #phi <- invgamma::rinvgamma(n=n, shape=0.5*v_0, rate=0.5*vs2_0)
  theta <- rnorm(n=n, mean=mu_0, sd=(phi*n_0)^(-1/2))
  return(cbind(theta, phi))
}


#'
#' @export
remainingSum <- function(m){
  M <- length(m)
  sums <- sapply(1:M, function(i) {sum(m[(i+1):M])})
  sums[M] <- 0
  sums
}

#'
#' @export
evaluate.ICDF <- function(theta, phi, weights, grid, ...){
  L = length(theta)
  toRet <- rep(weights, length(grid)) * plnorm(rep(grid, each=L),
                                            rep(theta, length(grid)),
                                            rep(sqrt(phi), length(grid)),
                                            lower.tail = F)
  toRet <- matrix(toRet, ncol = length(grid), byrow = F)
  toRet <- colSums(toRet)
  #toRet[toRet > 1] <- 1
  toRet
}

#'
#' @export
evaluate.PDF <- function(theta, phi, weights, L, grid, ...){
  toRet <- rep(weights, length(grid)) * dlnorm(rep(grid, each=L),
                                               rep(theta, length(grid)),
                                               rep(sqrt(phi), length(grid)))
  
  toRet <- matrix(toRet, ncol = length(grid), byrow = F)
  toRet <- colSums(toRet)
  toRet
}


#'
#' @export
stabilize <- function(mat){
  maxPerCol = apply(mat, 2, max)
  mat <- sweep(mat, 2, maxPerCol, "-")
  toRet <- exp(mat)
  toRet
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

#'
#' @export
DP_sample <- function(n, size = n, replace = FALSE, prob = NULL, max_lik = F){
  if(!is.na(prob)[1]){
    prob = prob/sum(prob)
    ifelse(max_lik, which.max(prob), sample.int(n, size, replace, prob))
  }else{
    return(NA)
  }
}

#'
#' @export
validate.App <- function(DP, DataStorage){
  MCIW <- rep(NA, length(unique(DataStorage@presentation$zeta)))
  for(i in as.numeric(unique(DataStorage@presentation$zeta))){
    grid <- subset(DataStorage@presentation, zeta==i)$data
    medianCurves <- getICDF.ChainStorage(DP=DP, myGrid=grid, zeta=i, quantiles=c(0.05,0.95))
    MCIW[i] <- mean(abs(medianCurves[2,,1]-medianCurves[1,,1]))
  }
  mean(MCIW)
}

#'
#' @export
validate.DP <- function(DP, DataStorage){
  medianCurves <- getICDF.ChainStorage(DP=DP, myGrid=DataStorage@validation$data, zeta=1:DP@J, quantiles=c(0.025,0.5,0.975))
  #mean_diff <- mean(apply(medianCurves[c(1,3),,],2, function(vec){return(vec[2]-vec[1])}))
  indices <- c(rep(1:10, length(DataStorage@mask)/3), 
               rep(11:20, length(DataStorage@mask)/3), 
               rep(21:30, length(DataStorage@mask)/3))
  score <- validate(curves=medianCurves, 
                    indices = indices,
                    status=rep(DataStorage@validation$status, length(DataStorage@mask)/3),
                    zeta=rep(1:length(DataStorage@mask), each=10))
  print(paste(class(DP)[1], score))
  return(score)
}

#'
#' @export
validate <- function(zeta, indices, curves, status){
  probability <- rep(NA, length(zeta))
  width <- rep(NA, length(zeta))
  for(i in 1:length(zeta)){
    if(!is.na(zeta[i])){
      probability[i] <- 1-curves[2, indices[i], as.numeric(as.character(zeta[i]))]
      width[i] <- curves[3, indices[i], as.numeric(as.character(zeta[i]))][[1]] - curves[1, indices[i], as.numeric(as.character(zeta[i]))][[1]]
    }
  }
  rmse <- RMSE(probability, status)#, LogScore(probability, status))
  return(c(rmse, mean(width, na.rm = T)))
}

#'
#'@export
mappingMu <- function(xi, zeta, mu,...){
  toRet <- rep(NA, length(xi))
  for(i in 1:length(xi)){
    toRet[i] <- mu[xi[i],zeta[i]]
  }
  #toRet <- matrix(toRet, ncol=1)
  toRet <- as.numeric(toRet)
  return(toRet)
}
