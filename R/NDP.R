#'Nested Dirichlet Process
#'NDP applied to censored data
#'@examples
#'\dontrun{
#'weights <- matrix(c(0,1,0,0,1,0,0,1,0), ncol=3)
#'data <- sim.data(weights)
#'G2 <- new("NDP")
#'G2 <- init.NDP(G2, prior=list(mu=0, n=0.1, v=3, vs2=1*3),K=5, L=35, thinning=2,
#'               burnin = 0, max_iter = 5000 )
#'G2 <- MCMC.NDP(G2, data, 500)
#plot.ICDF(G2@theta[,which.max(G2@pi)], G2@phi[,which.max(G2@pi)], G2@weights[,which.max(G2@pi)],
#          G2@L, grid=0:500, distribution=data@presentation, xlim=500)
#'validate.NDP(G2, data)
#'}
#'
#' @export
setClass("NDP", representation(DPs = 'list', K = 'numeric', phi = 'matrix', theta='matrix', weights='matrix', details='list', conc_param = 'numeric',
                              prior = 'array', L = 'numeric', pi='matrix', Chains='list', ChainStorage='ChainStorage'))

#'
#' @export
init.NDP <- function(NDP, prior, K, L, thinning, burnin, max_iter, ...){
  NDP@K <- K
  NDP@L <- L
  NDP@prior <- array(rep(c(prior$mu, prior$n, prior$v, prior$vs2), each=(K*L)), c(L,K,4))
  NDP@conc_param <- c(1,1)#rgamma(2,1,1)
  NDP <- update.NDP(NDP)
  NDP@Chains <- list(theta=NDP@theta, phi=NDP@phi, weights=NDP@weights, pi=NDP@pi)
  NDP@details <- list(iteration=0, thinning=thinning, burnin=burnin, max_iter=max_iter)
  NDP@ChainStorage <- init.ChainStorage(NDP@Chains, max_iter-burnin, thinning)
  return(NDP)
}

#'
#' @export
update.NDP <- function(NDP, ...){
  print(NDP@conc_param)
  atoms <- rNIG(NDP@L*NDP@K, c(NDP@prior[,,1]), c(NDP@prior[,,2]), c(NDP@prior[,,3]), c(NDP@prior[,,4]))
  NDP@theta <- matrix(atoms[,1], nrow=NDP@L)
  NDP@phi <- matrix(atoms[,2], nrow=NDP@L)
  sums <- apply(round(NDP@prior[,,2]), 2, mySums)
  u_k <- matrix(rbeta(NDP@L*NDP@K, shape1 = 1 + c(round(NDP@prior[,,2])), shape2 = NDP@conc_param[1] + c(sums)),
                   nrow=NDP@L)
  NDP@weights <- apply(u_k, 2, stickBreaking)
  sums <- mySums(colSums(round(NDP@prior[,,2])))
  v_lk <- rbeta(sums, 1 + colSums(round(NDP@prior[,,2])), NDP@conc_param[2] + sums)
  NDP@pi <- matrix(stickBreaking(v_lk), nrow=1)
  
  #a_conc <- 1
  #b_conc <- 1
  #NDP@conc_param[1] <- rgamma(1, a_conc + NDP@K - 1, b_conc - sum(log(1-u_k[1:NDP@K-1])))
  #NDP@conc_param[2] <- rgamma(1, a_conc + NDP@K*(NDP@L-1) - 1, b_conc - sum(log(1-v_lk[-seq(NDP@L, NDP@K*NDP@L, NDP@L)])))
  
  return(NDP)
}

#'
#' @export
selectZetaXi.NDP <- function(NDP, DataStorage, max_lik=F, ...){
    zz <- DPsurv::eStep(theta=c(NDP@theta), phi=c(NDP@phi), w=c(NDP@weights), DataStorage=DataStorage)
    mega_cube <- array(zz, c(NDP@L, NDP@K, max(DataStorage@mask)*length(DataStorage@mask)))
    mega_cube <- aperm(mega_cube, c(1,3,2))
    jj <- array(mega_cube, c(NDP@L, max(DataStorage@mask), length(DataStorage@mask), NDP@K))
    zz <- apply(exp(jj), c(3,4), sum, na.rm=T)
    for(i in 1:dim(zz)[1]){
      zz[i,] <- exp(log(zz[i,]) + log(NDP@pi))
    }
    zeta <- apply(zz, 1, DP_sample, n=NDP@K, size=1, replace=F, max_lik=max_lik)
    #zeta <- apply(zz, 1, sample.int, n=NDP@K, size=1, replace=F)
    kk <- sapply(1:length(zeta), function(i) return(jj[,, i, zeta[i]]))
    ii <- matrix(kk, nrow=NDP@L)
    ii <- stabilize(ii)
    xi <- apply(ii, 2, DP_sample, n=NDP@L, size=1, replace=F, max_lik=max_lik)
    return(list(zeta=zeta, xi=xi))
}

#'
#' @export
MCMC.NDP <- function(NDP, DataStorage, iter, ...){
  j <- 0
  while(NDP@details[['iteration']] < NDP@details[['max_iter']] & j < iter){
    
    NDP@details[['iteration']] <- NDP@details[['iteration']] + 1
    j <- j + 1
    ZetaXi <- selectZetaXi.NDP(NDP, DataStorage)
    zeta <- ZetaXi[["zeta"]]
    xi <- ZetaXi[["xi"]]
    DataStorage@presentation$zeta <- rep(zeta, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    DataStorage <- DPsurv::gibbsStep(DP=NDP, DataStorage=DataStorage, 
                                     xi=xi, zeta=rep(zeta, each = max(DataStorage@mask)))
    NDP <- DPsurv::mStep(NDP, DataStorage, xi=xi, zeta=rep(zeta, each = max(DataStorage@mask)))
    NDP <- update.NDP(NDP)
    NDP@Chains <- list(theta=NDP@theta, phi=NDP@phi, weights=NDP@weights, pi=NDP@pi)
    if(NDP@details[["iteration"]]>NDP@details[["burnin"]] & (NDP@details[["iteration"]] %% NDP@details[["thinning"]])==0){
      NDP@ChainStorage <- saveChain.ChainStorage(NDP@Chains, (NDP@details[["iteration"]]-NDP@details[["burnin"]])/NDP@details[["thinning"]], NDP@ChainStorage)
    }
  }
  return(NDP)
}


#'
#' @export
validate.NDP <- function(NDP, DataStorage){
  NDP <- posterior.DP(NDP, 0.5)
  ZetaXi <- selectZetaXi.NDP(NDP, DataStorage, max_lik = T)
  zeta <- ZetaXi[["zeta"]]
  DataStorage@presentation$zeta <- rep(zeta, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
  DataStorage@validation$zeta <- as.numeric(as.character(plyr::mapvalues(DataStorage@validation$Sample, DataStorage@presentation$Sample, DataStorage@presentation$zeta,  warn_missing=F)))
  score <- validate(data=DataStorage@validation$data, status=DataStorage@validation$status, zeta=DataStorage@validation$zeta,
                    theta=NDP@theta, phi=NDP@phi, weights=NDP@weights)
  return(score)
}

