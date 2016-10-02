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
setClass("DP", representation(weights = 'matrix', phi = 'matrix', theta = 'matrix', details='list', conc_param='numeric', J='numeric', K='numeric',
                              prior='numeric', posterior = 'array', L = 'numeric', Chains='list', ChainStorage='ChainStorage'))

#'
#' @export
init.DP <- function(DP, DataStorage, prior, K, J, thinning, burnin, max_iter, ...){
  DP <- methods::new("DP")
  DP@L <- 1
  DP@K <- K
  DP@J <- J
  DP@prior <- prior
  DP@posterior <- array(rep(c(prior[1], prior[2], prior[3], prior[4]), each=(K)), c(K,1,4))
  DP@conc_param <- rgamma(1, prior[5], prior[7])
  DP <- update.DP(DP)
  #DP@RE <- init.RE(DataStorage)
  myProb <- matrix(NA, nrow=K, ncol=J)
  myParams <- matrix(NA, nrow=1, ncol=J)
  DP@Chains <- list(theta=myParams, phi=myParams, prob=myProb)
  DP@details <- list(iteration=0, thinning=thinning, burnin=burnin, max_iter=max_iter)
  DP@ChainStorage <- init.ChainStorage(DP@Chains, max_iter-burnin, thinning)
  return(DP)
}

#'
#' @export
update.DP <- function(DP, ...){
  atoms <- rNIG(DP@K, mu_0=c(DP@posterior[,,1]), n_0=c(DP@posterior[,,2]),
                v_0=c(DP@posterior[,,3]), vs2_0=c(DP@posterior[,,4]))
  DP@theta <- matrix(atoms[,1], nrow=DP@K)
  DP@phi <- matrix(atoms[,2], nrow=DP@K)
  sums <- remainingSum(round(c(DP@posterior[,,2])))
  beta_0 <- rbeta(DP@K, shape1 = 1 + c(DP@posterior[,,2]), shape2 = DP@conc_param + sums)
  DP@weights <- as.matrix(stickBreaking(beta_0), ncol=1)
  DP@conc_param <- rgamma(1, DP@prior[5] + DP@L - 1, DP@prior[7] - sum(log(1-DP@weights[1:(DP@L-1)])))
  return(DP)
}

#'
#' @export
selectXi.DP <- function(DP, DataStorage, max_lik=F){
  log_prob <- eStep(data=DataStorage@computation, censoring=DataStorage@censoring, theta=DP@theta, phi=DP@phi, w=DP@weights)
  prob <- stabilize(log_prob)
  reshape_prob <- array(prob, c(DP@K, max(DataStorage@mask), length(DataStorage@mask)))
  sum_prob <- apply(reshape_prob, 3, rowSums, na.rm=T)
  zeta <- apply(sum_prob, 2, DP_sample, n=DP@K, size=1, replace=F, max_lik=max_lik)
  if(length(zeta) != length(DataStorage@computation)) {
    xi <- rep(zeta, each=max(DataStorage@mask))
    xi[is.na(DataStorage@computation)] <- NA
  }
  return(list(zeta=zeta, xi=xi, prob=sum_prob))
}

#'
#' @export
MCMC.DP <- function(DP, DataStorage, iter, ...){
  i <- 0
  pb <- txtProgressBar(style = 3)
  while(DP@details[['iteration']] < DP@details[['max_iter']] & i < iter){
    
    DP@details[['iteration']] <- DP@details[['iteration']] + 1
    i <- i + 1
    xiZeta <- selectXi.DP(DP, DataStorage)
    xi <- xiZeta[["xi"]]
    zeta <- xiZeta[["zeta"]]
    prob <- xiZeta[["prob"]]
    DataStorage@presentation$xi <- xi[!is.na(xi)] #rep(xi, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    DataStorage <- gibbsStep(DP=DP, DataStorage=DataStorage, xi=xi, zeta=rep(1, length(DataStorage@computation))) 
    
    DP@posterior <- mStep(DP@prior, DP@posterior, DataStorage@simulation, xi=xi, zeta=rep(1, length(DataStorage@computation)))
    if(F){
      map <- mappingMu(c(xi),rep(1,length(xi)), DP@theta)
      #DP@RE <- MCMC.RE(DP@RE, DataStorage@simulation-map)
    }
    if(DP@details[["iteration"]]>DP@details[["burnin"]] & (DP@details[["iteration"]] %% DP@details[["thinning"]])==0){
      setTxtProgressBar(pb, i/iter)
      DP@Chains <- list(theta=t(DP@theta), phi=t(DP@phi), prob=prob)
      DP@ChainStorage <- saveChain.ChainStorage(zeta=zeta, Chains=DP@Chains, iteration=(DP@details[["iteration"]]-DP@details[["burnin"]])/DP@details[["thinning"]], ChainStorage=DP@ChainStorage)
    }
    DP <- update.DP(DP)
  }
  close(pb)
  return(DP)
}

