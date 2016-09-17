#' Simulation
#'
#' 35 simulations
#'
#' @param seed
#' 
#' @return None
#'
#' @examples
#' Simulation(seed=1)
#'
#' @export
Simulation <- function(seed,replications=35,iterations=55000,burnin=5000,thinning=50,L=35,K=20,
                       n=100,J=10,case="1",factor=1){
  options(gsubfn.engine = "R")
  case <- toString(case)
  weights <- switch(case,
                    "0"=diag(3),
                    "1"=matrix(c(0.6,0.4,0, 0.25,.75,0, 0.25,0,0.75), ncol=3),
                    "2"=matrix(c(rep(1,3), rep(0,6)), ncol=3))
  filename = paste0("results","n", n, "J", J, "case", case, "factor", factor, "seed", seed, ".csv")
  write.table(t(c("Brier_DP", "Log_DP", "Brier_NDP", "Log_NDP", "Brier_HDP", "Log_HDP")), file = filename, sep = ",", col.names = F, row.names = F)
  set.seed(seed)

  scores <- parallel::mclapply(1:replications, function(i){
    data <- sim.data(n=n, J=J, weights=weights, factor=factor)
    print(paste("DP", i, sep=" "))
    G <- init.DP(prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, thinning=thinning, burnin=burnin, max_iter=iterations, DataStorage =  data)
    G1 <- MCMC.DP(G, data, iterations)
    score_DP <- validate.DP(G1, data)
    
    print(paste("NDP", i, sep=" "))
    G <- init.NDP(prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=K, L=L, thinning=thinning, burnin=burnin, max_iter=iterations)
    G2 <- MCMC.NDP(G, data, iterations)
    plotICDF.NDP(G2, data)
    
    score_NDP <- validate.NDP(G2, data)
    
    print(paste("HDP", i, sep=" "))
    G <- init.HDP(prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, 
                   J=length(levels(data@presentation$Sample)), thinning=thinning, burnin=burnin, max_iter=iterations)
    G3 <- MCMC.HDP(G, data, iterations)
    plotICDF.HDP(G3, data)
    score_HDP <- validate.HDP(G3, data)
    toRet <- c(score_DP, score_NDP, score_HDP)
    write.table(t(toRet), file = filename, sep = ",", col.names = F, row.names=F, append=TRUE)
  },mc.cores = 4 )
}

#' Simulation
#'
#' 35 simulations
#'
#' @param seed
#' 
#' @return None
#'
#' @examples
#' Simulation(seed=1)
#'
#' @export
Applications <- function(seed,iterations=55000,burnin=5000,thinning=20,L=35,K=20){
  options(gsubfn.engine = "R")
  filename = paste0("applications", "seed", seed, ".csv")
  write.table(t(c("Dataset", "Brier_DP", "Log_DP", "Brier_NDP", "Log_NDP", "Brier_HDP", "Log_HDP")), file = filename, sep = ",", col.names = F, row.names = F)

  data("performArt")
  myData <- c("performArt")
  
  set.seed(seed)
  
  scores <- parallel::mclapply(myData, function(dataset){
    data <- init.DataStorage.simple(get(dataset), 0.1)
    print(paste("DP", dataset, sep=" "))
    G <- new("DP")
    G <- init.DP(G, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, thinning=thinning, burnin=burnin, max_iter=iterations, DataStorage =  data)
    G <- MCMC.DP(G, data, iterations)
    score_DP <- validate.DP(G, data)
    
    print(paste("NDP", dataset, sep=" "))
    G <- new("NDP")
    G <- init.NDP(G, prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=K, L=L, thinning, burnin, iterations)
    G <- MCMC.NDP(G, data, iterations)
    score_NDP <- validate.NDP(G, data)
    
    print(paste("HDP", dataset, sep=" "))
    G <- new("HDP")
    G <- init.HDP(G, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, 
                   J=length(levels(data@presentation$Sample)), thinning, burnin, iterations)
    G <- MCMC.HDP(G, data, iterations)
    score_HDP <- validate.HDP(G, data)
    toRet <- c(dataset, score_DP, score_NDP, score_HDP)
    write.table(t(toRet), file = filename, sep = ",", col.names = F, row.names=F, append=TRUE)
    #return(toRet)
  },
  mc.cores = 4
  )
}