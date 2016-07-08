
init.Simulation <- function(replications=2, iterations=10, burnin=0, thinning=1, L=55, K=35){
  options(gsubfn.engine = "R")
  scores <- parallel::mclapply(1:replications, function(i){
    set.seed(i)
    data <- sim.data()
    G1 <- new("DP")
    G1 <- init.DP(G1, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, thinning, burnin, iterations)
    G1 <- MCMC.DP(G1, data, iterations)
    score_DP <- validate.DP(G1, data)
    
    G2 <- new("NDP")
    G2 <- init.NDP(G2, prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=K, L=L, thinning, burnin, iterations)
    G2 <- MCMC.NDP(G2, data, iterations)
    score_NDP <- validate.NDP(G2, data)
    
    G3 <- new("HDP")
    G3 <- init.HDP(G3, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, 
                   J=length(levels(data@presentation$Sample)), thinning, burnin, iterations)
    G3 <- MCMC.HDP(G3, data, iterations)
    score_HDP <- validate.HDP(G3, data)
    toRet <- c(score_DP, score_NDP, score_HDP)
    print(toRet)
    return(toRet)
  },
  mc.cores = 4
  )
  return(scores)
}
