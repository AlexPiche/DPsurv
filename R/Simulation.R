
init.Simulation <- function(filename, replications=35, iterations=55000, burnin=5000, thinning=20, L=35, K=20, n=500, J=20){
  options(gsubfn.engine = "R")
  write.table(t(c("score_DP", "score_NDP", "score_HDP")), file = paste0(filename, ".csv"), sep = ",", col.names = NA)
  scores <- parallel::mclapply(1:replications, function(i){
    #set.seed(i)
    data <- sim.data(n=n, J=J)
    print(paste("DP", i, sep=" "))
    G <- new("DP")
    G <- init.DP(G, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, thinning, burnin, iterations, clustering=F)
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
    write.table(t(toRet), file = paste0(filename, ".csv"), sep = ",", col.names = FALSE, append=TRUE)
    return(toRet)
  },
  mc.cores = 4
  )
  write.csv(scores, file="results_simulation_DP.csv")
  return(scores)
}
