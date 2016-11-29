#'
#' @export
Testing <- function(seed, data=NA, iterations=1000,burnin=500,thinning=50,L=55,K=35,valid_prop = 0.0, frailty=T, mySample=0,
                    n=20,J=10,case="1",factor=1, Prior = c(0, 0.1, 2*1/3, 2*1/3, rep(c(5, 0.1), 2)), plotting=F){
  options(gsubfn.engine = "R")
  case <- toString(case)
  weights <- switch(case,
                    "0"=matrix(c(rep(1,3), rep(0,6)), ncol=3, byrow = T),
                    "1"=matrix(c(0.6, 0.4, 0, 0, 0.4, 0.6, 0.25, 0, 0.75), ncol=3),
                    "2"=matrix(c(rep(0,6), rep(1,3)), ncol=3, byrow = T))
  set.seed(seed)
  if(is.na(data)){
    data <- sim.data(n=n, J=J, weights=weights, factor=factor, validation_prop = valid_prop, frailty)
  }
  J <- length(unique(data@presentation$Sample))
  
  Prior[1]<- log(median(data@presentation$data))
  
  parfm.model <- fitParfm(data)
  
  save(data, file = "data.RData")
  G <- init.DP(prior=Prior, K=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations, DataStorage =  data)
  G1 <- MCMC.DP(G, data, iterations)
  score <- validate.DP(G1, data)
  save(G1, file = "G1.RData")
  
  if(plotting) {
    plotICDF.DP(G1, data, mySample)
    plotHeatmap(G1)
  }
  
  G <- init.HDP(prior=Prior, L=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations)
  G3 <- MCMC.HDP(G, data, iterations)
  score <- validate.DP(G3, data)
  save(G3, file = "G3.RData")
  
  if(plotting) {
    plotICDF.DP(G3, data, mySample)
    plotHeatmap(G3)
  }
  
  G <- init.NDP(prior=Prior, K=K, L=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations)
  G2 <- MCMC.NDP(G, data, iterations)
  score <- validate.DP(G2, data)
  save(G2, file = "G2.RData")
  if(plotting) {
    plotICDF.DP(G2, data, mySample)
    plotHeatmap(G2)
  }
  
}

