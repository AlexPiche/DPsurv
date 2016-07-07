
init.Simulation <- function(n, iterations, burnin, thinning, L, K){
  parallel::mclapply(1:25, function(i){
    burnin <- 0
    thinning <- 1
    iterations <- 10
    L <- 55
    K <- 35
    set.seed(i)
    data <- sim.data()
    G1 <- new("DP")
    G1 <- init.DP(G1, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, thinning, burnin, iterations)
    G1 <- MCMC.DP(G1, data, iterations)
    
    G2 <- new("NDP")
    G2 <- init.NDP(G2, prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=K, L=L, thinning, burnin, iterations)
    G2 <- MCMC.NDP(G2, data, iterations)
    
    G3 <- new("HDP")
    G3 <- init.HDP(G3, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, 
                   J=length(levels(data@presentation$Sample)), thinning, burnin, iterations)
    G3 <- MCMC.HDP(G3, data, iterations)
  }
  )
}