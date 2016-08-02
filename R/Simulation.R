#' Illustration of crayon colors
#'
#' Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @param seed
#' 
#' @return None
#'
#' @examples
#' Simulation(seed=1)
#'
#' @export
Simulation <- function(seed,replications=35,iterations=55000,burnin=5000,thinning=20,L=35,K=20,
                       n=100,J=10, case=1, factor=1){
  options(gsubfn.engine = "R")
  weights <- switch(case,
                    "0"=matrix(c(0,1,0,0,1,0,0,1,0), ncol=3),
                    "1"=matrix(c(0.6,0.4,0, 0.25,.75,0, 0.25,0,0.75), ncol=3))
  data <- sim.data(n=n, J=J, weights=weights)
  filename = paste0("results","n", n, "J", J, "case", case, "factor", factor, "seed", seed, ".csv")
  write.table(t(c("Brier_DP", "Log_DP", "Brier_NDP", "Log_NDP", "Brier_HDP", "Log_HDP")), file = filename, sep = ",", col.names = F, row.names = F)
  set.seed(seed)
  scores <- parallel::mclapply(1:replications, function(i){
  #for(i in 1:replications){
    data <- sim.data(n=n, J=J, weights=weights, factor=factor)
    print(paste("DP", i, sep=" "))
    G <- new("DP")
    G <- init.DP(G, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, thinning, burnin, iterations, clustering=T)
    G <- MCMC.DP(G, data, iterations)
    score_DP <- validate.DP(G, data)
    
    print(paste("NDP", i, sep=" "))
    G <- new("NDP")
    G <- init.NDP(G, prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=K, L=L, thinning, burnin, iterations)
    G <- MCMC.NDP(G, data, iterations)
    score_NDP <- validate.NDP(G, data)
    
    print(paste("HDP", i, sep=" "))
    G <- new("HDP")
    G <- init.HDP(G, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, 
                   J=length(levels(data@presentation$Sample)), thinning, burnin, iterations)
    G <- MCMC.HDP(G, data, iterations)
    score_HDP <- validate.HDP(G, data)
    toRet <- c(score_DP, score_NDP, score_HDP)
    write.table(t(toRet), file = filename, sep = ",", col.names = F, row.names=F, append=TRUE)
    return(toRet)
  },
  mc.cores = 4
  )
}
