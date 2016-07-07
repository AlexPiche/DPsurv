# utils functions call by the other files

stickBreaking <- function(beta){
  beta[length(beta)] <- 1
  pi <- sapply(1:length(beta), function(i) {beta[i]*prod(1-beta[0:(i-1)])})
  return(pi)
}

rNIG <- function(n, mu_0, n_0, v_0, vs2_0){
  phi <- MCMCpack::rinvgamma(n=n, shape=0.5*v_0, scale=0.5*vs2_0)
  theta <- rnorm(n=n, mean=mu_0, sd=sqrt(phi/n_0))
  return(cbind(theta, phi))
}

updateChain <- function(){
  #take matrix or vector and append it to an array
  return(-1)
}

mySums <- function(m){
  M <- length(m)
  sums <- sapply(1:M, function(i) {sum(m[(i+1):M])})
  sums[M] <- 0
  sums
}

evaluate.ICDF <- function(theta, phi, weights, grid, ...){
  L = length(theta)
  toRet <- rep(weights, length(grid)) * plnorm(rep(grid, each=L),
                                            rep(theta, length(grid)),
                                            rep(sqrt(phi), length(grid)),
                                            lower.tail = F)
  toRet <- matrix(toRet, ncol = length(grid), byrow = F)
  toRet <- colSums(toRet)
  toRet
}

evaluate.PDF <- function(theta, phi, weights, L, grid, ...){
  toRet <- rep(weights, length(grid)) * dlnorm(rep(grid, each=L),
                                               rep(theta, length(grid)),
                                               rep(sqrt(phi), length(grid)))
  
  toRet <- matrix(toRet, ncol = length(grid), byrow = F)
  toRet <- colSums(toRet)
  toRet
}


stabilize <- function(mat){
  maxpercol = apply(mat, 2, max)
  maxPerColMat = matrix(rep(maxpercol, each=nrow(mat)), ncol=ncol(mat))
  mat = mat - maxPerColMat
  toRet <- exp(mat)
  toRet
}

BrierScore <- function(probability, status){
  toRet <- mean((probability-status)^2)
  return(toRet)
}
