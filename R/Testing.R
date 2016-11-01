#'
#' @export
Testing <- function(seed, data=NA, iterations=55000,burnin=5000,thinning=50,L=5,K=3, valid_prop = 0.0,
                    n=10,J=10,case="1",factor=1, Prior = c(0, 0.01, 2*1/3, 2*1/3, rep(c(5, 0.1), 2)), plotting=F){
  options(gsubfn.engine = "R")
  case <- toString(case)
  weights <- switch(case,
                    "0"=diag(3),
                    "1"=matrix(c(0.6, 0.4, 0, 0, 0.4, 0.6, 0.25, 0, 0.75), ncol=3),
                    "2"=matrix(c(rep(0,6), rep(1,3)), ncol=3, byrow = T))
  set.seed(seed)
  if(is.na(data)){
    data <- sim.data(n=n, J=J, weights=weights, factor=factor, validation_prop = valid_prop)
  }
  J <- length(unique(data@presentation$Sample))
  
  Prior[1]<- log(median(data@presentation$data))
  

  ICDF <- parfm_survival_curves_estimate(data@presentation, data@validation)
  score_parFM <- validation_parFM(ICDF, data@validation$status)
  
  save(data, file = "data.RData")
  G <- init.DP(prior=Prior, K=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations, DataStorage =  data)
  G1 <- MCMC.DP(G, data, iterations)
  save(G1, file = "G1.RData")
  
  if(plotting) {
    plotICDF.DP(G1, data)
    plotHeatmap(G1)
  }
  score_DP <- validate.DP(G1, data)
  
  G <- init.NDP(prior=Prior, K=K, L=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations)
  G2 <- MCMC.NDP(G, data, iterations)
  save(G2, file = "G2.RData")
  if(plotting) {
    plotICDF.DP(G2, data)
    plotHeatmap(G2)
  }
  score_NDP <- validate.DP(G2, data)
  
  G <- init.HDP(prior=Prior, L=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations)
  G3 <- MCMC.HDP(G, data, iterations)
  save(G3, file = "G3.RData")
  
  if(plotting) {
    plotICDF.DP(G3, data)
    plotHeatmap(G3)
  }
  score_HDP <- validate.DP(G3, data)
  #toRet <- c(score_DP, score_NDP, score_HDP)
  #print(toRet)
}

