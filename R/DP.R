#'
#' @useDynLib DPsurv
#' @export
setClass("DP", representation(weights = 'matrix', phi = 'matrix', theta = 'matrix', details='list', clustering = 'logical',
                              prior = 'array', L = 'numeric', Chains='list', ChainStorage='ChainStorage'))

#'
#' @export
init.DP <- function(DP, prior, L, thinning, burnin, max_iter, clustering=F, ...){
  DP@L <- L
  DP@clustering <- clustering
  DP@prior <- array(rep(c(prior$mu, prior$n, prior$v, prior$vs2), each=(L)), c(L,1,4))
  DP <- update.DP(DP)
  DP@Chains <- list(theta=DP@theta, phi=DP@phi, weights=DP@weights)
  DP@details <- list(iteration=0, thinning=thinning, burnin=burnin, max_iter=max_iter)
  DP@ChainStorage <- init.ChainStorage(DP@Chains, max_iter-burnin, thinning)
  return(DP)
}

#'
#' @export
update.DP <- function(DP, ...){
  atoms <- rNIG(DP@L, c(DP@prior[,,1]), c(DP@prior[,,2]), c(DP@prior[,,3]), c(DP@prior[,,4]))
  DP@theta <- matrix(atoms[,1], nrow=DP@L)
  DP@phi <- matrix(atoms[,2], nrow=DP@L)
  sums <- mySums(round(c(DP@prior[,,2])))
  beta_0 <- rbeta(DP@L, shape1 = 1 + c(DP@prior[,,2]), shape2 = 1 + sums)
  DP@weights <- as.matrix(stickBreaking(beta_0), ncol=1)
  return(DP)
}

#'
#' @export
selectXi.DP <- function(DP, DataStorage, max_lik=F){
    xx <- eStep(DP@theta, DP@phi, DP@weights, DataStorage)
    xx <- stabilize(xx)
    zz <- array(xx, c(DP@L, max(DataStorage@mask), length(DataStorage@mask)))
    rr <- if(DP@clustering) apply(zz, 3, rowSums, na.rm=T) else xx
    xi <- apply(rr, 2, DP_sample, n=DP@L, size=1, replace=F, max_lik=max_lik)
    if(length(xi) != length(DataStorage@computation)) {
      xi <- rep(xi, each=max(DataStorage@mask))
      xi[is.na(DataStorage@computation)] <- NA
    }
    return(xi)
}

#'
#' @export
MCMC.DP <- function(DP, DataStorage, iter, ...){
  i <- 0
  while(DP@details[['iteration']] < DP@details[['max_iter']] & i < iter){

    DP@details[['iteration']] <- DP@details[['iteration']] + 1
    i <- i + 1
    xi <- selectXi.DP(DP, DataStorage)
    DataStorage@presentation$xi <- xi[!is.na(xi)] #rep(xi, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    DataStorage <- gibbsStep(DP=DP, DataStorage=DataStorage, xi=xi, zeta=rep(1, length(DataStorage@computation))) 
    DP <- mStep(DP, DataStorage, xi=xi,
                        zeta=rep(1, length(DataStorage@computation)))
    
    DP <- update.DP(DP)
    DP@Chains <- list(theta=DP@theta, phi=DP@phi, weights=DP@weights)
    if(DP@details[["iteration"]]>DP@details[["burnin"]] & (DP@details[["iteration"]] %% DP@details[["thinning"]])==0){
      DP@ChainStorage <- saveChain.ChainStorage(DP@Chains, (DP@details[["iteration"]]-DP@details[["burnin"]])/DP@details[["thinning"]], DP@ChainStorage)
    }
  }
  return(DP)
}

#'
#' @export
validate.DP <- function(DP, DataStorage){
  DP <- posterior.DP(DP, 0.5)
  if(DP@clustering){
    xi <- selectXi.DP(DP, DataStorage, max_lik = T)
    DataStorage@presentation$xi <- xi[!is.na(xi)] #rep(xi, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    DataStorage <- update_validation_set(DataStorage)
    score <- validate(data=DataStorage@validation$data, status=DataStorage@validation$status, zeta=DataStorage@validation$xi,
                      theta=t(DP@theta), phi=t(DP@phi), weights=diag(DP@L))
  }else{
    score <- validate(data=DataStorage@validation$data, status=DataStorage@validation$status, zeta=rep(1, length(DataStorage@validation$status)),
                      theta=DP@theta, phi=DP@phi, weights=DP@weights)
  }
  return(score)
}
