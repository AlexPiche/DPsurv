#' Simulation
#'
#' Simulate heterogeneous datasets of censored observations to be model by three different
#' Bayesian nonparametric methods. It writes the results in a csv.
#'
#' @param seed
#' 
#' @return None
#'
#' @examples
#' \dontrun{
#' Simulation(seed=1)
#' }
#'
#' @export
Simulation <- function(seed,replications=1,iterations=55000,burnin=5000,thinning=50,L=35,K=20,
                       n=10,J=10,case="1",factor=1, Prior = c(0, 0.01, 2*1/3, 2*1/3, rep(c(5, 0.1), 2)), val_prop=0, nb_cores=4){
  options(gsubfn.engine = "R")
  
  case <- toString(case)
  
  weights <- switch(case,
                    "0"=diag(3),
                    "1"=matrix(c(0.6, 0.4, 0, 0, 0.4, 0.6, 0.25, 0, 0.75), ncol=3),
                    "2"=matrix(c(rep(0,6), rep(1,3)), ncol=3, byrow = T))
  
  filename = paste0("n", n, "J", J, "case", case, "seed", seed, ".csv")
  
  colnames_cov <- c("model", "n", "curve", "0.05","0.15","0.25","0.35","0.45","0.55","0.65","0.75","0.85", "0.95" )
  colnames_res <- c("RMSE_DP", "CI_DP", "RMSE_NDP", "CI_NDP", "RMSE_HDP", "CI_HDP", "RMSE_ParFM", "CI_ParFM")
  write.table(t(colnames_res), file = paste0("res",filename), sep = ",", col.names = F, row.names = F)
  write.table(t(colnames_cov), file = paste0("cov",filename), sep = ",", col.names = F, row.names = F)
  
  
  set.seed(seed)
  
  scores <- parallel::mclapply(1:replications, function(i){
    data <- sim.data(n=n, J=J, weights=weights, factor=factor, validation_prop = val_prop)
    Prior[1] <- log(median(data@presentation$data))
    
    # This line is commented out, since it depends on a modified version of parfm
    # to get confidence interval on the survival curves
    #matrix_medianCurves <- getMedianCurves.ParFM(data@presentation, data@validation)
    score_parFM <- c(0,0) #validate.Score(matrix_medianCurves, data@validation)
    cov_parFM <- matrix(rep(0,3*10), nrow=3) #validate.Coverage(matrix_medianCurves, data@validation)
    print(paste("GFM", score_parFM))
    
    print(paste("DP", i, sep=" "))
    G <- init.DP(prior=Prior, K=L, J=3*J, thinning=thinning, burnin=burnin, max_iter=iterations, DataStorage =  data)
    G <- MCMC.DP(G, data, iterations)
    matrix_medianCurves <- getMedianCurves(G, data)    
    score_DP <- validate.Score(matrix_medianCurves, data@validation)
    cov_DP <- validate.Coverage(matrix_medianCurves, data@validation)
    
    print(paste("NDP", i, sep=" "))
    G <- init.NDP(prior=Prior, K=K, L=L, J=3*J, thinning=thinning, burnin=burnin, max_iter=iterations)
    G <- MCMC.NDP(G, data, iterations)
    matrix_medianCurves <- getMedianCurves(G, data)    
    score_NDP <- validate.Score(matrix_medianCurves, data@validation)
    cov_NDP <- validate.Coverage(matrix_medianCurves, data@validation)
    
    print(paste("HDP", i, sep=" "))
    G <- init.HDP(prior=Prior, L=L, J=3*J, thinning=thinning, burnin=burnin, max_iter=iterations)
    G <- MCMC.HDP(G, data, iterations)
    matrix_medianCurves <- getMedianCurves(G, data)    
    score_HDP <- validate.Score(matrix_medianCurves, data@validation)
    cov_HDP <- validate.Coverage(matrix_medianCurves, data@validation)
    
    toRet <- c(score_DP, score_NDP, score_HDP, score_parFM)
    write.table(t(toRet), file = paste0("res",filename), sep = ",", col.names = F, row.names=F, append=TRUE)
    
    toRet <- rbind(cov_DP, cov_NDP, cov_HDP, cov_parFM)
    
    write.table(cbind(rep(c("DP","NDP", "HDP", "GFM"),each=3), rep(n, 12), rep(c(1,2,3), 4), toRet), file = paste0("cov",filename), sep = ",", col.names = F, row.names=F, append=TRUE)
  }, mc.cores = nb_cores )
}

#' Application
#'
#' Produce the experiment on the performance art dataset.
#'
#' @param seed
#' 
#' @return None
#'
#' @examples
#' \dontrun{
#' Application(seed=1)
#' }
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
  parallel::mclapply(c(4), function(i){
    print(i)
    J <- length(unique(data@presentation$zeta))
    if(i == 1){
      print("DP")
      G <- init.DP(prior=Prior, K=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations, DataStorage =  data)
      G <- MCMC.DP(G, data, iterations)
      matrix_medianCurves <- getMedianCurves(G, data)    
      score_DP <- validate.Score(matrix_medianCurves, data@validation)
      print(paste("DP", score_DP))
    }else if (i == 2){
      print("NDP")
      G <- init.NDP(prior=Prior, K=K, L=L, J=J, thinning=thinning, burnin=burnin, max_iter=iterations)
      G <- MCMC.NDP(G, data, iterations)
      matrix_medianCurves <- getMedianCurves(G, data)    
      score_NDP <- validate.Score(matrix_medianCurves, data@validation)
      print(paste("NDP", score_NDP))
    } else if(i == 3) {
      print("HDP")
      G <- init.HDP(prior=Prior, L=L, 
                    J=J, thinning=thinning, burnin=burnin, max_iter=iterations)
      G <- MCMC.HDP(G, data, iterations)
      matrix_medianCurves <- getMedianCurves(G, data)    
      score_HDP <- validate.Score(matrix_medianCurves, data@validation)
      print(paste("HDP", score_HDP))
      
    } else {
      data@presentation$Sample <- as.numeric(as.factor(data@presentation$Sample))
      matrix_medianCurves <- getMedianCurves.ParFM(data@presentation, data@validation)
      score_parFM <- validate.Score(matrix_medianCurves, data@validation)
      print(paste("GFM", score_parFM))
      #ICDF <- parfm_survival_curves_estimate(data@presentation, data@validation)
      #score_parFM <- validation_parFM(ICDF, data@validation$status)
    }
    #toRet <- c(dataset, score_DP, score_NDP, score_HDP, T)
    #write.table(t(toRet), file = filename, sep = ",", col.names = F, row.names=F, append=TRUE)
    #return(toRet)
  }, mc.cores = 4
  )
}