#'
#' @export
setClass("DP", representation(weights = 'matrix', phi = 'matrix', theta = 'matrix', details='list',
                              prior = 'array', L = 'numeric', Chains='list', ChainStorage='ChainStorage'))

init.DP <- function(DP, prior, L, thinning, burnin, max_iter, ...){
  DP@L <- L
  DP@prior <- array(rep(c(prior$mu, prior$n, prior$v, prior$vs2), each=(L)), c(L,1,4))
  DP <- update.DP(DP)
  DP@Chains <- list(theta=DP@theta, phi=DP@phi, weights=DP@weights)
  DP@details <- list(iteration=0, thinning=thinning, burnin=burnin, max_iter=max_iter)
  DP@ChainStorage <- init.ChainStorage(DP@Chains, max_iter-burnin, thinning)
  return(DP)
}

update.DP <- function(DP, ...){
  atoms <- rNIG(DP@L, c(DP@prior[,,1]), c(DP@prior[,,2]), c(DP@prior[,,3]), c(DP@prior[,,4]))
  DP@theta <- matrix(atoms[,1], nrow=DP@L)
  DP@phi <- matrix(atoms[,2], nrow=DP@L)
  sums <- mySums(round(c(DP@prior[,,2])))
  beta_0 <- rbeta(DP@L, shape1 = 1 + c(DP@prior[,,2]), shape2 = 1 + sums)
  DP@weights <- as.matrix(stickBreaking(beta_0), ncol=1)
  return(DP)
}

MCMC.DP <- function(DP, DataStorage, iter, ...){
  i <- 0
  while(DP@details[['iteration']] < DP@details[['max_iter']] & i < iter){

    DP@details[['iteration']] <- DP@details[['iteration']] + 1
    i <- i + 1
    xx <- DPsurv::eStep(DP@theta, DP@phi, DP@weights, DataStorage)
    xx <- stabilize(xx)
    zz <- array(xx, c(DP@L, max(DataStorage@mask), length(DataStorage@mask)))
    rr <- apply(zz, 3, rowSums, na.rm=T)
    xi <- apply(rr, 2, sample.int, n=DP@L, size=1, replace=F)
    DataStorage@presentation$xi <- rep(xi, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    DataStorage <- DPsurv::gibbsStep(DP=DP, DataStorage=DataStorage, 
                                     RealData=DataStorage@presentation$data, xi=rep(xi, each=max(DataStorage@mask)), 
                                     zeta=rep(1, length(DataStorage@computation)), censoring=DataStorage@presentation$status) 
    DP <- DPsurv::mStep(DP, DataStorage, xi=rep(xi, each=max(DataStorage@mask)),
                        zeta=rep(1, length(DataStorage@computation)))
    
    DataStorage@computation <- DataStorage@presentation$data
    DP <- update.DP(DP)
    DP@Chains <- list(theta=DP@theta, phi=DP@phi, weights=DP@weights)
    if(DP@details[["iteration"]]>DP@details[["burnin"]]){
      DP@ChainStorage <- saveChain.ChainStorage(DP@Chains, (DP@details[["iteration"]]-DP@details[["burnin"]])/DP@details[["thinning"]], DP@ChainStorage)
    }
    validate.DP(DP, DataStorage)
  }
  return(DP)
}

validate.DP <- function(DP, DataStorage){
  DataStorage <- update_validation_set(DataStorage)
  score <- validate(data=DataStorage@validation$data, status=DataStorage@validation$status, zeta=DataStorage@validation$xi,
                    theta=t(DP@theta), phi=t(DP@phi), weights=diag(DP@L))
}
