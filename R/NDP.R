#'Nested Dirichlet Process
#'NDP applied to censored data
#'@examples
#'\dontrun{
#'weights <- matrix(c(1,0,0,0,1,0,0,0,1), ncol=3)
#'
#'data <- sim.data(n=100, J=10, weights)
#'
#'G2 <- init.NDP(prior=list(mu=0, n=0.1, v=3, vs2=1*3),K=5, L=35, thinning=50,
#'               burnin = 0, max_iter = 5000 )
#'               
#'G2 <- MCMC.NDP(G2, data, 500)
#'
#'validate.NDP(G2, data)
#'
#'plot.ICDF(G2@theta[,which.max(G2@pi)], G2@phi[,which.max(G2@pi)], G2@weights[,which.max(G2@pi)],
#'             G2@L, grid=0:500, distribution=data@presentation, xlim=500)
#'}
#'
#' @export
setClass("NDP", representation(DPs = 'list', K = 'numeric', phi = 'matrix', theta='matrix', weights='matrix', details='list', conc_param = 'numeric', J = 'numeric',
                               m = 'numeric', prior = 'numeric', posterior = 'array', L = 'numeric', pi='matrix', Chains='list', ChainStorage='ChainStorage'))

#'
#' @export
init.NDP <- function(prior, J, K, L, thinning, burnin, max_iter, ...){
  NDP <- methods::new("NDP")
  NDP@K <- K
  NDP@L <- L
  NDP@J <- J
  NDP@prior <- prior
  NDP@posterior <- array(rep(c(prior[1], prior[2], prior[3], prior[4]), each=(K*L)), c(L,K,4))
  NDP@conc_param <- rgamma(2, prior[c(5,7)], prior[c(6,8)])
  NDP@m <- rep(0, NDP@K)
  NDP <- update.NDP(NDP)
  myProb <- matrix(NA, nrow=K, ncol=J)
  myParams <- matrix(NA, nrow=L, ncol=J)
  NDP@Chains <- list(theta=myParams, phi=myParams, weights=myParams, prob=myProb)
  NDP@details <- list(iteration=0, thinning=thinning, burnin=burnin, max_iter=max_iter)
  NDP@ChainStorage <- init.ChainStorage(NDP@Chains, max_iter-burnin, thinning)
  return(NDP)
}

#'
#' @export
update.NDP <- function(NDP, ...){
  atoms <- rNIG(NDP@L*NDP@K, c(NDP@posterior[,,1]), c(NDP@posterior[,,2]), c(NDP@posterior[,,3]), c(NDP@posterior[,,4]))
  NDP@theta <- matrix(atoms[,1], nrow=NDP@L)
  NDP@phi <- matrix(atoms[,2], nrow=NDP@L)
  sums <- remainingSum(NDP@m)
  u_k <- rbeta(sums, 1 + NDP@m, NDP@conc_param[1] + sums)
  NDP@pi <- matrix(stickBreaking(u_k), nrow=1)
  sums <- apply(round(NDP@posterior[,,2]), 2, remainingSum)
  v_lk <- matrix(rbeta(NDP@L*NDP@K, shape1 = 1 + c(round(NDP@posterior[,,2])), shape2 = NDP@conc_param[2] + c(sums)),nrow=NDP@L)
  NDP@weights <- apply(v_lk, 2, stickBreaking)
  
  NDP@conc_param[1] <- rgamma(1, NDP@prior[5] + NDP@K - 1, NDP@prior[6] - sum(log(1-u_k[1:NDP@K-1])))
  NDP@conc_param[2] <- rgamma(1, NDP@prior[7] + NDP@K*(NDP@L-1), NDP@prior[8] - sum(log(1-v_lk[-seq(NDP@L, NDP@K*NDP@L, NDP@L)])))
  return(NDP)
}

#'
#' @export
selectZetaXi.NDP <- function(NDP, DataStorage, max_lik=F, ...){
  log_prob <- DPsurv::eStep(data=DataStorage@computation, censoring=DataStorage@censoring, theta=c(NDP@theta), phi=c(NDP@phi), w=c(NDP@weights))
  log_prob <- array(log_prob, c(NDP@L, NDP@K, max(DataStorage@mask)*length(DataStorage@mask)))
  log_prob <- aperm(log_prob, c(1,3,2))
  reshape_log_prob <- array(log_prob, c(NDP@L, max(DataStorage@mask), length(DataStorage@mask), NDP@K))
  DP_prob <- apply(exp(reshape_log_prob), c(3,4), sum, na.rm=T)
  weighted_DP_prob <- exp(log(DP_prob)+log(matrix(rep(NDP@pi, dim(DP_prob)[1]), ncol=NDP@K, byrow = T)))
  zeta <- apply(weighted_DP_prob, 1, DP_sample, n=NDP@K, size=1, replace=F, max_lik=max_lik)
  atom_log_prob <- sapply(1:length(zeta), function(i) return(reshape_log_prob[,, i, zeta[i]]))
  atom_log_prob <- matrix(atom_log_prob, nrow=NDP@L)
  atom_prob <- stabilize(atom_log_prob)
  xi <- apply(atom_prob, 2, DP_sample, n=NDP@L, size=1, replace=F, max_lik=max_lik)
  return(list(zeta=zeta, xi=xi, prob=weighted_DP_prob))
}

#'
#' @export
selectZetaXi.NDP_c <- compiler::cmpfun(selectZetaXi.NDP)


#'
#' @export
MCMC.NDP <- function(NDP, DataStorage, iter, ...){
  j <- 0
  pb <- utils::txtProgressBar(style = 3)
  while(NDP@details[['iteration']] < NDP@details[['max_iter']] & j < iter){
    
    NDP@details[['iteration']] <- NDP@details[['iteration']] + 1
    j <- j + 1
    ZetaXi <- selectZetaXi.NDP(NDP, DataStorage)
    zeta <- ZetaXi[["zeta"]]
    xi <- ZetaXi[["xi"]]
    NDP@m <- c(table(factor(zeta, levels=1:NDP@K)))
    DataStorage@presentation$zeta <- rep(zeta, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    DataStorage <- DPsurv::gibbsStep(DP=NDP, DataStorage=DataStorage, 
                                     xi=xi, zeta=rep(zeta, each = max(DataStorage@mask)))
    NDP@posterior <- DPsurv::mStep(NDP@prior, NDP@posterior, DataStorage@simulation, xi=xi, zeta=rep(zeta, each = max(DataStorage@mask)))
    if(NDP@details[["iteration"]]>NDP@details[["burnin"]] & (NDP@details[["iteration"]] %% NDP@details[["thinning"]])==0){
      utils::setTxtProgressBar(pb, j/iter)
      NDP@Chains <- list(theta=NDP@theta, phi=NDP@phi, weights=NDP@weights, prob=t(ZetaXi[["prob"]]))
      NDP@ChainStorage <- saveChain.ChainStorage(zeta=zeta, Chains=NDP@Chains, iteration=(NDP@details[["iteration"]]-NDP@details[["burnin"]])/NDP@details[["thinning"]], ChainStorage=NDP@ChainStorage)
    }
    NDP <- update.NDP(NDP)
  }
  close(pb)
  return(NDP)
}

#'
#' @export
MCMC.NDP_c <- compiler::cmpfun(MCMC.NDP)
