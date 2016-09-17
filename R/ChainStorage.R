#'
#' @export
setClass("ChainStorage", representation(chains='list'))

#'
#' @export
init.ChainStorage <- function(Chains, iterations, thinning, ...){
  ChainStorage <- new('ChainStorage')
  chains <- list()
  for(chain in names(Chains)){
    chains[[chain]] = array(dim=c(dim(Chains[[chain]]), iterations/thinning))
  }
  ChainStorage@chains <- chains
  return(ChainStorage)
}

#'
#' @export
saveChain.ChainStorage <- function(zeta, xi, Chains, iteration, ChainStorage, ...){
  for(chain in names(Chains)){
    if(dim(Chains[[chain]])[2] >= length(zeta)){
      if(chain == "weights") xi <- 1:dim(Chains[[chain]])[1]
      if(dim(Chains[[chain]])[1]==1) xi <- 1

      ChainStorage@chains[[chain]][xi, zeta, iteration] <- Chains[[chain]][xi, zeta]
    }else{
      #HDP case
      ChainStorage@chains[[chain]][,,iteration] <- Chains[[chain]]
    }
  }
  return(ChainStorage)
}


#'
#' @export
getQuantile.ChainStorage <- function(ChainStorage, quantiles, ...){
  lapply(ChainStorage@chains, apply, c(1,2), quantile, quantiles, na.rm=T)
}


#'
#' @export
getICDF.ChainStorage <- function(DP, myGrid, zeta, quantiles=0.5){
  theta_mat <- matrix(DP@ChainStorage@chains[["theta"]][,zeta,], nrow=DP@L)
  ind_na <- apply(theta_mat, 2, function(x) all(is.na(x)))
  theta_mat <- theta_mat[, !ind_na]
  theta_mat <- split(theta_mat, col(theta_mat))
  
  phi_mat <- matrix(DP@ChainStorage@chains[["phi"]][,zeta,], nrow=DP@L)
  phi_mat <- phi_mat[, !ind_na]
  phi_mat <- split(phi_mat, col(phi_mat))
  
  weights_mat <- matrix(DP@ChainStorage@chains[["weights"]][,zeta,], nrow=DP@L)
  weights_mat <- weights_mat[, !ind_na]
  weights_mat <- split(weights_mat, col(weights_mat))
  curves <- mapply(evaluate.ICDF, theta_mat, phi_mat, weights_mat, MoreArgs = list(grid=myGrid))
  curvesReshape <- array(curves, c(dim(curves)[1], length(zeta), sum(!ind_na)/length(zeta)))
  meanCurves <- apply(curvesReshape, c(1,2), quantile, quantiles, na.rm = T)
  return(meanCurves)
}

#'
#' @export
posterior.DP <- function(DP, quantiles){
  csMedian <- getQuantile.ChainStorage(DP@ChainStorage, quantiles)
  DP@theta <- csMedian[["theta"]]
  DP@theta[is.na(DP@theta)] <- 0
  DP@phi <- csMedian[["phi"]]
  DP@phi[is.na(DP@phi)] <- 100
  DP@weights <- apply(DP@ChainStorage@chains[["weights"]], c(1,2), median, na.rm=T)#csMedian[["weights"]]
  DP@weights[is.na(DP@weights)] <- rbeta(length(is.na(DP@weights)),1,1)
  DP@weights <- apply(DP@weights, 2, normalizeVec)
  #DP@weights <- csMedian[["weights"]]
  #DP@RE@computation <- csMedian[["RE"]]
  if(length(csMedian[["pi"]]) > 1){
    DP@pi <- apply(DP@ChainStorage@chains[["pi"]], c(1,2), median, na.rm=T)
    DP@pi[is.na(DP@pi)] <- 0
    DP@pi <- normalizeVec(DP@pi)
    #DP@pi <- normalizeVec(normalizeVec(csMedian[["pi"]]))
  }
  return(DP)
}

#'
#' @export
normalizeVec <- function(vec){
  toRet <- vec/sum(vec)
  return(toRet)
}

#'
#' @export
getCIICDF.ChainStorage <- function(DP){
  return(-1)
}