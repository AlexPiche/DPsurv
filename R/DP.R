#'Dirichlet Process
#'DP applied to censored data
#'@examples
#'\dontrun{
#'weights <- matrix(c(1,0,0,0,1,0,0,0,1), ncol=3)
#'
#'data <- sim.data(n=100, J=10, weights)
#'
#'G1 <- init.DP(DataStorage=data, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=35, thinning=50,
#'           burnin = 5000, max_iter = 55000 )
#'           
#'G1 <- MCMC.DP(G1, data, 55000)
#'
#'validate.DP(G1, data)
#'  
#'plot.ICDF(G1@theta, G1@phi, G1@weights, G1@L, grid=0:500,
#'            distribution=data@presentation, xlim=500)
#'}
#'
#' @export
setClass("DP", representation(weights = 'matrix', phi = 'matrix', theta = 'matrix', details='list', RE='RE', conc_param='numeric',
                              prior = 'array', L = 'numeric', Chains='list', ChainStorage='ChainStorage'))

#'
#' @export
init.DP <- function(DP, DataStorage, prior, L, thinning, burnin, max_iter, ...){
  DP <- new("DP")
  DP@L <- L
  DP@prior <- array(rep(c(prior$mu, prior$n, prior$v, prior$vs2), each=(L)), c(L,1,4))
  DP@conc_param <- 1
  #a_conc <- 5
  #b_conc <- 0.1
  #DP@conc_param <- rgamma(1, a_conc, b_conc)
  DP <- update.DP(DP)
  DP@RE <- init.RE(DataStorage)
  DP@Chains <- list(theta=DP@theta, phi=DP@phi, weights=DP@weights, RE=DP@RE@computation)
  DP@details <- list(iteration=0, thinning=thinning, burnin=burnin, max_iter=max_iter)
  DP@ChainStorage <- init.ChainStorage(DP@Chains, max_iter-burnin, thinning)
  return(DP)
}

#'
#' @export
update.DP <- function(DP, ...){
  atoms <- rNIG(DP@L, c(DP@prior[,,1]), c(DP@prior[,,2]), c(DP@prior[,,3]), c(DP@prior[,,4]))
  DP@theta <- matrix(atoms[,1], nrow=DP@L)
  DP@phi <- matrix(atoms[,2], nrow=DP@L)
  sums <- remainingSum(round(c(DP@prior[,,2])))
  beta_0 <- rbeta(DP@L, shape1 = 1 + c(DP@prior[,,2]), shape2 = DP@conc_param + sums)
  DP@weights <- as.matrix(stickBreaking(beta_0), ncol=1)
  #a_conc <- 5
  #b_conc <- 0.1
  #DP@conc_param <- rgamma(1, a_conc + DP@L - 1, b_conc - sum(log(1-DP@weights[1:(DP@L-1)])))
  return(DP)
}

#'
#' @export
selectXi.DP <- function(DP, DataStorage, max_lik=F){
    log_prob <- eStep(data=DataStorage@computation-c(t(DP@RE@computation)), censoring=DataStorage@censoring, theta=DP@theta, phi=DP@phi, w=DP@weights)
    prob <- stabilize(log_prob)
    reshape_prob <- array(prob, c(DP@L, max(DataStorage@mask), length(DataStorage@mask)))
    sum_prob <- apply(reshape_prob, 3, rowSums, na.rm=T)
    xi <- apply(sum_prob, 2, DP_sample, n=DP@L, size=1, replace=F, max_lik=max_lik)
    if(length(xi) != length(DataStorage@computation)) {
      xi <- rep(xi, each=max(DataStorage@mask))
      xi[is.na(DataStorage@computation)] <- NA
    }
    return(xi)
}

#'
#' @export
MCMC.DP <- function(DP, DataStorage, iter, ...){
  i <- 0
  pb <- txtProgressBar(style = 3)
  while(DP@details[['iteration']] < DP@details[['max_iter']] & i < iter){

    DP@details[['iteration']] <- DP@details[['iteration']] + 1
    i <- i + 1
    xi <- selectXi.DP(DP, DataStorage)
    
    DataStorage@presentation$xi <- xi[!is.na(xi)] #rep(xi, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    DataStorage <- gibbsStep(DP=DP, DataStorage=DataStorage, xi=xi, zeta=rep(1, length(DataStorage@computation))) 
    
    DP@prior <- mStep(DP@prior, DataStorage@simulation-c(t(DP@RE@computation)), xi=xi, zeta=rep(1, length(DataStorage@computation)))
    DP <- update.DP(DP)
    if(F){
      map <- mappingMu(c(xi),rep(1,length(xi)), DP@theta)
      DP@RE <- MCMC.RE(DP@RE, DataStorage@simulation-map)
    }
    if(DP@details[["iteration"]]>DP@details[["burnin"]] & (DP@details[["iteration"]] %% DP@details[["thinning"]])==0){
      setTxtProgressBar(pb, i/iter)
      DP@Chains <- list(theta=DP@theta, phi=DP@phi, weights=DP@weights)#, RE=DP@RE@computation)
      DP@ChainStorage <- saveChain.ChainStorage(1, unique(xi[!is.na(xi)]), DP@Chains, (DP@details[["iteration"]]-DP@details[["burnin"]])/DP@details[["thinning"]], DP@ChainStorage)
    }
  }
  close(pb)
  return(DP)
}



#'
#' @export
diagonalize.DP <- function(DP){
  DP <- posterior.DP(DP, 0.5)
  theta_array <- array(0, dim=c(DP@L, DP@L, dim(DP@ChainStorage@chains[["theta"]])[3]))
  phi_array <- array(0, dim=c(DP@L, DP@L, dim(DP@ChainStorage@chains[["theta"]])[3]))
  weights_array <- array(0, dim=c(DP@L, DP@L, dim(DP@ChainStorage@chains[["theta"]])[3]))
  for(i in 1:dim(DP@ChainStorage@chains[["theta"]])[3]){
    theta_array[,,i] <- diag(DP@ChainStorage@chains[["theta"]][,,i])
    phi_array[,,i] <- diag(DP@ChainStorage@chains[["phi"]][,,i])
    weights_array[,,i] <- diag(DP@L)
  }
  DP@ChainStorage@chains[["theta"]] <- theta_array
  DP@ChainStorage@chains[["phi"]] <- phi_array
  DP@ChainStorage@chains[["weights"]] <- weights_array
  return(DP)
}

#'
#' @export
validate.DP <- function(DP, DataStorage){
  DP <- diagonalize.DP(DP)
  xi <- selectXi.DP(DP, DataStorage, max_lik = T)
  xi <- xi[!is.na(xi)]
  medianCurves <- getICDF.ChainStorage(DP, DataStorage@validation$data, zeta=unique(xi), quantiles=c(0.05,0.5,0.95))
  medianCurvesArranged <- matrix(NA, dim(medianCurves)[2], dim(DP@ChainStorage@chains[["theta"]])[2])
  ciCurvesArranged <- array(NA, c(2, dim(medianCurves)[2], dim(DP@ChainStorage@chains[["theta"]])[2]))
  for(i in 1:length(unique(xi))) medianCurvesArranged[,xi[i]] <- medianCurves[2,,i]
  for(i in 1:length(unique(xi))) ciCurvesArranged[1:2,1:27,xi[i]] <- medianCurves[c(1,3),,i]
  DataStorage@presentation$xi <- xi #rep(xi, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
  mapping <- unique(DataStorage@presentation[, c("xi", "Sample")])
  DataStorage@validation$xi <- plyr::mapvalues(DataStorage@validation$Sample, mapping$Sample, mapping$xi, warn_missing=F)
  score <- validate(curves=medianCurvesArranged, status=DataStorage@validation$status, zeta=DataStorage@validation$xi)
  mean_diff <- mean(apply(ciCurvesArranged[,,unique(xi)],2, function(vec){return(vec[2]-vec[1])}))
  score <- c(score, mean_diff)
  return(score)
}

#'
#' @export
plotICDF.DP <- function(DP, DataStorage){
  DP <- diagonalize.DP(DP)
  
  xi <- selectXi.DP(DP, DataStorage, max_lik = T)
  xi <- xi[!is.na(xi)]
  DataStorage@presentation$zeta <- xi
  for(zeta in unique(DataStorage@presentation$zeta)){
    p <- plot.ICDF(DP, zeta, DataStorage@presentation) + 
      ggplot2::ggtitle(paste("DP", zeta))
    print(p)
  }
}
