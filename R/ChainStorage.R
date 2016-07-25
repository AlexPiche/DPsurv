#'
#' @export
setClass("ChainStorage", representation(chains='list'))

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
saveChain.ChainStorage <- function(Chains, iteration, ChainStorage, ...){
  for(chain in names(Chains)){
    ChainStorage@chains[[chain]][,,iteration] <- Chains[[chain]]
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
getICDF.ChainStorage <- function(DP, max_iter, myGrid){
  theta_mat <- matrix(DP@ChainStorage@chains[["theta"]][,,1:max_iter], nrow=DP@L)
  theta_mat <- split(theta_mat, col(theta_mat))
  phi_mat <- matrix(DP@ChainStorage@chains[["phi"]][,,1:max_iter], nrow=DP@L)
  phi_mat <- split(phi_mat, col(phi_mat))
  weights_mat <- matrix(DP@ChainStorage@chains[["weights"]][,,1:max_iter], nrow=DP@L)
  weights_mat <- split(weights_mat, col(weights_mat))
  toRet <- mapply(evaluate.ICDF, theta_mat, phi_mat, weights_mat, MoreArgs = list(grid=myGrid))
  toRet <- array(toRet, c(dim(toRet)[1], dim(DP@theta)[2], max_iter))
  toRet
}

#'
#' @export
posterior.DP <- function(DP, quantiles){
  csMedian <- getQuantile.ChainStorage(DP@ChainStorage, quantiles)
  DP@theta <- csMedian[["theta"]]
  DP@phi <- csMedian[["phi"]]
  DP@weights <- csMedian[["weights"]]
  if(length(csMedian[["pi"]]) > 1) DP@pi <- csMedian[["pi"]]
  return(DP)
}

#'
#' @export
getCIICDF.ChainStorage <- function(DP){
  return(-1)
}