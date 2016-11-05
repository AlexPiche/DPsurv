#'
#' @export
setClass("ChainStorage", representation(chains='list'))

#'
#' @export
init.ChainStorage <- function(Chains, iterations, thinning, ...){
  ChainStorage <- methods::new('ChainStorage')
  chains <- list()
  for(chain in names(Chains)){
    chains[[chain]] = array(dim=c(dim(Chains[[chain]]), iterations/thinning))
  }
  ChainStorage@chains <- chains
  return(ChainStorage)
}

#'
#' @export
saveChain.ChainStorage <- function(zeta, Chains, iteration, ChainStorage, ...){
  for(chain in names(Chains)){
    for(i in 1:length(zeta)){
      if(chain=="prob"){
        ChainStorage@chains[[chain]][, i, iteration] <- Chains[[chain]][, i]
      } else if(dim(Chains[[chain]])[2] > 1){
        ChainStorage@chains[[chain]][, i, iteration] <- Chains[[chain]][, zeta[i]]
      } else {
        # HDP
        ChainStorage@chains[[chain]][, i, iteration] <- Chains[[chain]][, 1]
      }
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
getICDF.ChainStorage <- function(DP, validation, quantiles=0.5, i=0){
  J <- unique(validation$Sample)
  if(i==0) i <- 1:dim(DP@ChainStorage@chains[["theta"]])[3]
  N <- dim(DP@ChainStorage@chains[["theta"]])[3]*length(J)
  
  theta_mat <- matrix(DP@ChainStorage@chains[["theta"]][,J,i], nrow=DP@L)
  theta_mat <- split(t(theta_mat), 1:N)
    
  phi_mat <- matrix(DP@ChainStorage@chains[["phi"]][,J,i], nrow=DP@L)
  phi_mat <- split(t(phi_mat), 1:N)
  
  if(!is.null(DP@ChainStorage@chains[["weights"]])){
    # not the DP
    weights_mat <- matrix(DP@ChainStorage@chains[["weights"]][,J,i], nrow=DP@L)
  }else{
    weights_mat <- rep(1, N)
  }
  weights_mat <- split(t(weights_mat), 1:N)
  
  X <- lapply(unique(validation$Sample), function(ss) t(matrix(subset(validation, Sample == ss)$data)))
  myGrid <- do.call(plyr::rbind.fill.matrix, X)
  myGrid <- apply(myGrid, 2, rep, each=length(i)) 
  myGrid_mat <- split(myGrid, 1:N)
  
  curves <- mapply(evaluate.ICDF, theta_mat, phi_mat, weights_mat, myGrid_mat)
  curvesReshape <- array(curves, c(dim(curves)[1], length(J), N/length(J)))
  medianCurves <- apply(curvesReshape, c(1,2), quantile, quantiles, na.rm = T)
  return(medianCurves)
}

#'
#' @export
posterior.DP <- function(DP, quantiles){
  csMedian <- getQuantile.ChainStorage(DP@ChainStorage, quantiles)
  DP@theta <- csMedian[["theta"]]
  DP@phi <- csMedian[["phi"]]
  DP@weights <- apply(DP@ChainStorage@chains[["weights"]], c(1,2), median, na.rm=T)#csMedian[["weights"]]
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