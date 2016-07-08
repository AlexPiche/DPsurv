#'
#' @export
setClass("HDP", representation(phi = 'matrix', theta='matrix', weights='matrix', Nmat='matrix', details='list',
                              prior = 'array', L = 'numeric', J = 'numeric', Chains='list', ChainStorage='ChainStorage'))

init.HDP <- function(HDP, prior, L, J, thinning, burnin, max_iter, ...){
  HDP@L <- L
  HDP@J <- J
  HDP@prior <- array(rep(c(prior$mu, prior$n, prior$v, prior$vs2), each=(L*J)), c(L,J,4))
  HDP@Nmat <- matrix(0, nrow=L, ncol=J)
  HDP@details <- list(iteration=0, thinning=thinning, burnin=burnin, max_iter=max_iter)
  HDP <- update.HDP(HDP)
  HDP@Chains <- list(theta=HDP@theta, phi=HDP@phi, weights=HDP@weights)
  HDP@ChainStorage <- init.ChainStorage(HDP@Chains, max_iter-burnin, thinning)
  return(HDP)
}

update.HDP <- function(HDP, ...){
  atoms <- rNIG(HDP@L, c(HDP@prior[,,1]), c(HDP@prior[,,2]), c(HDP@prior[,,3]), c(HDP@prior[,,4]))
  HDP@theta <- matrix(atoms[,1], nrow=HDP@L)
  HDP@phi <- matrix(atoms[,2], nrow=HDP@L)
  # 1st level
  sums <- mySums(round(c(HDP@prior[,,2])))
  beta_0 <- rbeta(HDP@L, shape1 = 1 + c(HDP@prior[,,2]), shape2 = 1 + sums)
  beta_0 <- stickBreaking(beta_0)
  
  # 2nd level
  zz <- matrix(apply(HDP@Nmat, 2, mySums), nrow=1)
  kk <- matrix(HDP@Nmat, nrow=1)
  alpha <- 1
  sums <- sapply(1:HDP@L, function(i) {alpha*(1-sum(beta_0[1:i]))})
  beta_j <- matrix(rbeta(HDP@J*HDP@L, rep(beta_0, HDP@J)+kk, 
                         zz+rep(sums, HDP@J)), ncol = HDP@J, byrow = F)
  
  HDP@weights <- apply(beta_j, 2, stickBreaking)
  return(HDP)
}

computeXi.HDP <- function(HDP, DataStorage, ...){
  xx <- DPsurv::eStep(theta=HDP@theta, phi=HDP@phi, w=rep(1, length(HDP@phi)), DataStorage=DataStorage)
  xx <- stabilize(xx)
  zz <- array(xx, c(HDP@L, max(DataStorage@mask), length(DataStorage@mask)))
  
  for(i in 1:dim(zz)[3]){
    zz[,,i] <- exp(log(zz[,,i]) + log(HDP@weights[,i]))
  }
  
  rr <- matrix(zz, nrow=HDP@L)
  xi <- apply(rr, 2, mySample, n=HDP@L, size=1, replace=F)
  return(xi)
}

MCMC.HDP <- function(HDP, DataStorage, iter, ...){
  i <- 0
  while(HDP@details[['iteration']] < HDP@details[['max_iter']] & i < iter){

    HDP@details[['iteration']] <- HDP@details[['iteration']] + 1
    i <- i + 1
    xi <- computeXi.HDP(HDP, DataStorage)
    DataStorage@presentation$zeta <- rep(1:HDP@J, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    Nmat <- matrix(factor(xi, levels = 1:HDP@L), ncol=HDP@J)
    Nmat <- apply(Nmat, 2, as.numeric)
    HDP@Nmat <- apply(Nmat, 2, test, L = HDP@L)
    DataStorage <- DPsurv::gibbsStep(DP=HDP, DataStorage=DataStorage, RealData=DataStorage@presentation$data, xi=xi, 
                                                 zeta=rep(1, length(xi)), censoring=DataStorage@presentation$status) 
    HDP <- DPsurv::mStep(HDP, DataStorage, xi, rep(1, length(xi))) 
    DataStorage@computation <- DataStorage@presentation$data
    HDP <- update.HDP(HDP)
    if(HDP@details[["iteration"]]>HDP@details[["burnin"]] & (HDP@details[["iteration"]] %% HDP@details[["thinning"]])==0){
      HDP@ChainStorage <- saveChain.ChainStorage(HDP@Chains, (HDP@details[["iteration"]]-HDP@details[["burnin"]])/HDP@details[["thinning"]], HDP@ChainStorage)
    }
  }
  return(HDP)
}

validate.HDP <- function(HDP, DataStorage){
  HDP <- posterior.DP(HDP, 0.5)
  DataStorage <- update_validation_set(DataStorage)
  score <- validate(data=DataStorage@validation$data, status=DataStorage@validation$status, zeta=DataStorage@validation$xi,
                    theta=matrix(rep(c(HDP@theta), HDP@J), ncol=HDP@J), phi=matrix(rep(c(HDP@phi), HDP@J), ncol=HDP@J), weights=HDP@weights)
  return(score)
}

test <- function(vec, L){
  vec <- factor(vec, 1:L)
  table(vec, useNA = 'no')
}

