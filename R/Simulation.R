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
                       n=100,J=10,case="1",factor=1, Prior = c(0, 0.01, 2*1/3, 2*1/3, rep(c(5, 0.1), 2)), val_prop=0, nb_cores=4){
  options(gsubfn.engine = "R")
  
  case <- toString(case)
  
  weights <- switch(case,
                    "0"=diag(3),
                    "1"=matrix(c(0.6, 0.4, 0, 0, 0.4, 0.6, 0.25, 0, 0.75), ncol=3),
                    "2"=matrix(c(rep(0,6), rep(1,3)), ncol=3, byrow = T))
  
  filename = paste0("results","n", n, "J", J, "case", case, "seed", seed, ".csv")
  
  write.table(t(c("RMSE_DP", "CI_DP", "RMSE_NDP", "CI_NDP", "RMSE_HDP", "CI_HDP",
                  "RMSE_ParFM", "CI_ParFM")), file = filename, sep = ",", col.names = F,
              row.names = F)
  
  set.seed(seed)
  
  scores <- parallel::mclapply(1:replications, function(i){
    data <- sim.data(n=n, J=J, weights=weights, factor=factor, validation_prop = val_prop)
    Prior[1] <- log(median(data@presentation$data))
    
    #ICDF <- parfm_survival_curves_estimate(data@presentation, data@validation$data, 10)
    score_parFM <- c(0,0)#validation_parFM(ICDF, data@validation$status)
    
    print(paste("DP", i, sep=" "))
    G <- init.DP(prior=Prior, K=L, J=3*J, thinning=thinning, burnin=burnin, max_iter=iterations, DataStorage =  data)
    G <- MCMC.DP(G, data, iterations)
    score_DP <- validate.DP(G, data)
    
    print(paste("NDP", i, sep=" "))
    G <- init.NDP(prior=Prior, K=K, L=L, J=3*J, thinning=thinning, burnin=burnin, max_iter=iterations)
    G <- MCMC.NDP(G, data, iterations)
    
    score_NDP <- validate.DP(G, data)
    
    print(paste("HDP", i, sep=" "))
    G <- init.HDP(prior=Prior, L=L, J=3*J, thinning=thinning, burnin=burnin, max_iter=iterations)
    G <- MCMC.HDP(G, data, iterations)
    score_HDP <- validate.DP(G, data)
    toRet <- c(score_DP, score_NDP, score_HDP, score_parFM)
    write.table(t(toRet), file = filename, sep = ",", col.names = F, row.names=F, append=TRUE)
  }, mc.cores = nb_cores )
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
Application <- function(seed,iterations=55000,burnin=5000,thinning=50,L=55,K=35){
  options(gsubfn.engine = "R")
  filename = paste0("applications", "seed", seed, ".csv")
  write.table(t(c("Dataset", "Log_DP", "Log_NDP", "Log_HDP")), file = filename, sep = ",", col.names = F, row.names = F)

  data("performArt")
  myData <- c("performArt")
  
  set.seed(seed)
  
  #for(dataset in myData){
  data <- init.DataStorage.simple(performArt, 0, weights=0, application=T)
  data@validation$Sample <- data@validation$zeta
  Prior = c(median(log(data@presentation$data)), 0.01, 2*1/3, 2*1/3, 5, 0.1, 5, 0.1)
  parallel::mclapply(2, function(i){
  #for(i in 2:3){
    print(i)
    J <- length(unique(data@presentation$zeta))
    if(i == 1){
      print("DP")
      G <- init.DP(prior=Prior, K=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations, DataStorage =  data)
      G1 <- MCMC.DP(G, data, iterations)
      score_DP <- validate.DP(G1, data)
    }else if (i == 2){
      print("NDP")
      G <- init.NDP(prior=Prior, K=K, L=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations)
      G2 <- MCMC.NDP(G, data, iterations)
      save(G2, file="G2pa.RData")
      score_DP <- validate.DP(G2, data)
    } else if(i == 3) {
      print("HDP")
      G <- init.HDP(prior=Prior, L=L, 
                    J=J, thinning=thinning, burnin=burnin, max_iter=iterations)
      G3 <- MCMC.HDP(G, data, iterations)
      save(G3, file="G3pa.RData")
      score_DP <- validate.DP(G3, data)
    } else {
      data@presentation$Sample <- as.numeric(as.factor(data@presentation$Sample))
      ICDF <- parfm_survival_curves_estimate(data@presentation, data@validation)
      score_parFM <- validation_parFM(ICDF, data@validation$status)
    }
    #toRet <- c(dataset, score_DP, score_NDP, score_HDP, T)
    #write.table(t(toRet), file = filename, sep = ",", col.names = F, row.names=F, append=TRUE)
    #return(toRet)
  }, mc.cores = 4
  )
}