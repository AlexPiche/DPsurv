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
Simulation <- function(seed,replications=32,iterations=55000,burnin=5000,thinning=50,L=35,K=20,
                       n=100,J=10,case="1",factor=1, val_prop=0, nb_cores=4){
  options(gsubfn.engine = "R")
  case <- toString(case)
  weights <- switch(case,
                    "0"=diag(3),
                    "1"=matrix(c(0.6,0.4,0, 0.25,.75,0, 0.25,0,0.75), ncol=3),
                    "2"=matrix(c(rep(1,3), rep(0,6)), ncol=3))
  filename = paste0("results","n", n, "J", J, "case", case, "seed", seed, ".csv")
  write.table(t(c("Brier_DP", "CI_DP", "Brier_NDP", "CI_NDP", "Brier_HDP", "CI_HDP")), file = filename, sep = ",", col.names = F, row.names = F)
  set.seed(seed)

  scores <- parallel::mclapply(1:replications, function(i){
    data <- sim.data(n=n, J=J, weights=weights, factor=factor, validation_prop = val_prop)
    print(paste("DP", i, sep=" "))
    G <- init.DP(prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=L, J=3*J, thinning=thinning, burnin=burnin, max_iter=iterations, DataStorage =  data)
    G <- MCMC.DP(G, data, iterations)
    score_DP <- validate.DP(G, data)
    
    print(paste("NDP", i, sep=" "))
    G <- init.NDP(prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=K, L=L, J=3*J, thinning=thinning, burnin=burnin, max_iter=iterations)
    G <- MCMC.NDP(G, data, iterations)
    
    score_NDP <- validate.DP(G, data)
    
    print(paste("HDP", i, sep=" "))
    G <- init.HDP(prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, 
                   J=length(levels(data@presentation$Sample)), thinning=thinning, burnin=burnin, max_iter=iterations)
    G <- MCMC.HDP(G, data, iterations)
    score_HDP <- validate.DP(G, data)
    toRet <- c(score_DP, score_NDP, score_HDP)
    write.table(t(toRet), file = filename, sep = ",", col.names = F, row.names=F, append=TRUE)
  },mc.cores = nb_cores )
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
Applications <- function(seed,iterations=55000,burnin=5000,thinning=50,L=35,K=20){
  options(gsubfn.engine = "R")
  filename = paste0("applications", "seed", seed, ".csv")
  write.table(t(c("Dataset", "Log_DP", "Log_NDP", "Log_HDP")), file = filename, sep = ",", col.names = F, row.names = F)

  data("performArt")
  myData <- c("performArt")
  
  set.seed(seed)
  
  for(dataset in myData){
    data <- init.DataStorage.simple(get(dataset), 0, weights=0, application=T)
    
    print(paste("DP", dataset, sep=" "))
    G <- new("DP")
    J <- length(unique(data@presentation$zeta))
    G <- init.DP(G, prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations, DataStorage =  data)
    G1 <- MCMC.DP(G, data, iterations)
    save(G1, "G1pA.Rdata")
    #score_DP <- validate.DP(G1, data)[2]
    G1 <- NA
    
    print(paste("NDP", dataset, sep=" "))
    G <- new("NDP")
    G <- init.NDP(G, prior=list(mu=0, n=0.1, v=3, vs2=1*3), K=K, L=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations)
    G2 <- MCMC.NDP(G, data, iterations)
    save(G2, "G2pA.Rdata")
    #score_NDP <- validate.DP(G2, data)[2]
    G2 <- NA
    
    print(paste("HDP", dataset, sep=" "))
    G <- new("HDP")
    G <- init.HDP(G, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=L, 
                   J=length(levels(data@presentation$Sample)), thinning=thinning, burnin=burnin, max_iter=iterations)
    G3 <- MCMC.HDP(G, data, iterations)
    save(G3, "G3pA.Rdata")
    #score_HDP <- validate.DP(G3, data)[2]
    G3 <- NA
    
    #toRet <- c(dataset, score_DP, score_NDP, score_HDP, T)
    #write.table(t(toRet), file = filename, sep = ",", col.names = F, row.names=F, append=TRUE)
    #return(toRet)
  }
}