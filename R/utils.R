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
  phi <- MCMCpack::rinvgamma(n=n, shape=0.5*v_0, scale=0.5*vs2_0)
  theta <- rnorm(n=n, mean=mu_0, sd=sqrt(phi/n_0))
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
  toRet[toRet > 1] <- 1
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
  maxpercol = apply(mat, 2, max)
  maxPerColMat = matrix(rep(maxpercol, each=nrow(mat)), ncol=ncol(mat))
  mat = mat - maxPerColMat
  toRet <- exp(mat)
  toRet
}

#'
#' @export
BrierScore <- function(probability, status){
  toRet <- (probability-status)^2
  return(mean(sqrt(toRet), na.rm=T))
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
validate <- function(zeta, curves, status){
  probability <- rep(NA, length(zeta))
  for(i in 1:length(zeta)){
    if(!is.na(zeta[i])){
      probability[i] <- 1-curves[i, as.numeric(as.character(zeta[i]))]
    }
  }
  toRet <- BrierScore(probability, status)#, LogScore(probability, status))
  return(toRet)
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
