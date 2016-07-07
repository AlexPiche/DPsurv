#'
#' @export
setClass("NDP", representation(DPs = 'list', K = 'numeric', phi = 'matrix', theta='matrix', weights='matrix', details='list',
                              prior = 'array', L = 'numeric', pi='matrix', Chains='list', ChainStorage='ChainStorage'))

init.NDP <- function(NDP, prior, K, L, thinning, burnin, max_iter, ...){
  NDP@K <- K
  NDP@L <- L
  NDP@prior <- array(rep(c(prior$mu, prior$n, prior$v, prior$vs2), each=(K*L)), c(L,K,4))
  NDP <- update.NDP(NDP)
  NDP@Chains <- list(theta=NDP@theta, phi=NDP@phi, weights=NDP@weights, pi=NDP@pi)
  NDP@details <- list(iteration=0, thinning=thinning, burnin=burnin, max_iter=max_iter)
  NDP@ChainStorage <- init.ChainStorage(NDP@Chains, max_iter-burnin, thinning)
  return(NDP)
}

update.NDP <- function(NDP, ...){
  atoms <- rNIG(NDP@L*NDP@K, c(NDP@prior[,,1]), c(NDP@prior[,,2]), c(NDP@prior[,,3]), c(NDP@prior[,,4]))
  NDP@theta <- matrix(atoms[,1], nrow=NDP@L)
  NDP@phi <- matrix(atoms[,2], nrow=NDP@L)
  sums <- apply(round(NDP@prior[,,2]), 2, mySums)
  beta_0 <- matrix(rbeta(NDP@L*NDP@K, shape1 = 1 + c(round(NDP@prior[,,2])), shape2 = 1 + c(sums)),
                   nrow=NDP@L)
  NDP@weights <- apply(beta_0, 2, stickBreaking)
  sums <- mySums(colSums(round(NDP@prior[,,2])))
  beta_0 <- rbeta(sums, 1 + colSums(round(NDP@prior[,,2])), 1 + sums)
  NDP@pi <- matrix(stickBreaking(beta_0), nrow=1)
  
  return(NDP)
}

selectZetaXi.NDP <- function(NDP, DataStorage, ...){
    zz <- DPsurv::eStep(theta=c(NDP@theta), phi=c(NDP@phi), w=c(NDP@weights), DataStorage=DataStorage)
    mega_cube <- array(zz, c(NDP@L, NDP@K, max(DataStorage@mask)*length(DataStorage@mask)))
    mega_cube <- aperm(mega_cube, c(1,3,2))
    jj <- array(mega_cube, c(NDP@L, max(DataStorage@mask), length(DataStorage@mask), NDP@K))
    zz <- apply(exp(jj), c(3,4), sum, na.rm=T)
    for(i in 1:dim(zz)[1]){
      zz[i,] <- exp(log(zz[i,]) + log(NDP@pi))
    }
    zeta <- apply(zz, 1, sample.int, n=NDP@K, size=1, replace=F)
    kk <- sapply(1:length(zeta), function(i) return(jj[,, i, zeta[i]]))
    ii <- matrix(kk, nrow=NDP@L)
    ii <- stabilize(ii)
    xi <- apply(ii, 2, mySample, n=NDP@L, size=1, replace=F)
    return(list(zeta=zeta, xi=xi))
}

MCMC.NDP <- function(NDP, DataStorage, iter, ...){
  j <- 0
  while(NDP@details[['iteration']] < NDP@details[['max_iter']] & j < iter){
    
    NDP@details[['iteration']] <- NDP@details[['iteration']] + 1
    j <- j + 1
    ZetaXi <- selectZetaXi.NDP(NDP, DataStorage)
    zeta <- ZetaXi[["zeta"]]
    xi <- ZetaXi[["xi"]]
    DataStorage@presentation$zeta <- rep(zeta, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    DataStorage <- DPsurv::gibbsStep(DP=NDP, DataStorage=DataStorage, RealData=DataStorage@presentation$data,
                                     xi=xi, zeta=rep(zeta, each = max(DataStorage@mask)),
                                     censoring=DataStorage@presentation$status)
    NDP <- DPsurv::mStep(NDP, DataStorage, xi=xi, zeta=rep(zeta, each = max(DataStorage@mask)))
    DataStorage@computation <- DataStorage@presentation$data
    NDP <- update.NDP(NDP)
    NDP@Chains <- list(theta=NDP@theta, phi=NDP@phi, weights=NDP@weights, pi=NDP@pi)
    if(NDP@details[["iteration"]]>NDP@details[["burnin"]] & (NDP@details[["iteration"]] %% NDP@details[["thinning"]])==0){
      NDP@ChainStorage <- saveChain.ChainStorage(NDP@Chains, (NDP@details[["iteration"]]-NDP@details[["burnin"]])/NDP@details[["thinning"]], NDP@ChainStorage)
    }
  }
  return(NDP)
}


validate.NDP <- function(NDP, DataStorage){
  NDP <- posterior.DP(NDP, 0.5)
  ZetaXi <- selectZetaXi.NDP(NDP, DataStorage)
  zeta <- ZetaXi[["zeta"]]
  DataStorage@presentation$zeta <- rep(zeta, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
  DataStorage <- update_validation_set(DataStorage)
  score <- validate(data=DataStorage@validation$data, status=DataStorage@validation$status, zeta=DataStorage@validation$xi,
                    theta=NDP@theta, phi=NDP@phi, weights=NDP@weights)
  return(score)
}

mySample <- function(n, size = n, replace = FALSE, prob = NULL){
  if(!is.na(prob)[1]){
    prob = prob/sum(prob)
    sample.int(n, size, replace, prob)
  }else{
    return(NA)
  }
}