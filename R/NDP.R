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
setClass("NDP", representation(DPs = 'list', K = 'numeric', phi = 'matrix', theta='matrix', weights='matrix', details='list', conc_param = 'numeric',
                               prior = 'array', L = 'numeric', pi='matrix', Chains='list', ChainStorage='ChainStorage'))

#'
#' @export
init.NDP <- function(prior, K, L, thinning, burnin, max_iter, ...){
  NDP <- new("NDP")
  NDP@K <- K
  NDP@L <- L
  NDP@prior <- array(rep(c(prior$mu, prior$n, prior$v, prior$vs2), each=(K*L)), c(L,K,4))
  NDP@conc_param <- c(1,1)
  #NDP@conc_param[1] <- rgamma(1,5,0.1)
  #NDP@conc_param[2] <- rgamma(1,0.1,0.1)
  NDP <- update.NDP(NDP)
  NDP@Chains <- list(theta=NDP@theta, phi=NDP@phi, weights=NDP@weights, pi=NDP@pi)
  NDP@details <- list(iteration=0, thinning=thinning, burnin=burnin, max_iter=max_iter)
  NDP@ChainStorage <- init.ChainStorage(NDP@Chains, max_iter-burnin, thinning)
  return(NDP)
}

#'
#' @export
update.NDP <- function(NDP, ...){
  atoms <- rNIG(NDP@L*NDP@K, c(NDP@prior[,,1]), c(NDP@prior[,,2]), c(NDP@prior[,,3]), c(NDP@prior[,,4]))
  NDP@theta <- matrix(atoms[,1], nrow=NDP@L)
  NDP@phi <- matrix(atoms[,2], nrow=NDP@L)
  sums <- remainingSum(colSums(round(NDP@prior[,,2])))
  u_k <- rbeta(sums, 1 + colSums(round(NDP@prior[,,2])), NDP@conc_param[1] + sums)
  NDP@pi <- matrix(stickBreaking(u_k), nrow=1)
  #NDP@pi <- matrix(c(1, rep(0, length(NDP@pi)-1)))
  sums <- apply(round(NDP@prior[,,2]), 2, remainingSum)
  v_lk <- matrix(rbeta(NDP@L*NDP@K, shape1 = 1 + c(round(NDP@prior[,,2])), shape2 = NDP@conc_param[2] + c(sums)),nrow=NDP@L)
  NDP@weights <- apply(v_lk, 2, stickBreaking)
  a_conc1 <- 5
  b_conc1 <- 0.1
  a_conc2 <- 0.1
  b_conc2 <- 0.1
  #NDP@conc_param[1] <- rgamma(1, a_conc1 + NDP@K - 1, b_conc1 - sum(log(1-u_k[1:NDP@K-1])))
  #NDP@conc_param[2] <- rgamma(1, a_conc2 + NDP@K*(NDP@L-1) - 1, b_conc2 - sum(log(1-v_lk[-seq(NDP@L, NDP@K*NDP@L, NDP@L)])))
  
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
  return(list(zeta=zeta, xi=xi))
}

#'
#' @export
selectZetaXi.NDP_c <- compiler::cmpfun(selectZetaXi.NDP)


#'
#' @export
MCMC.NDP <- function(NDP, DataStorage, iter, ...){
  j <- 0
  pb <- txtProgressBar(style = 3)
  while(NDP@details[['iteration']] < NDP@details[['max_iter']] & j < iter){
    
    NDP@details[['iteration']] <- NDP@details[['iteration']] + 1
    j <- j + 1
    ZetaXi <- selectZetaXi.NDP(NDP, DataStorage)
    zeta <- ZetaXi[["zeta"]]
    xi <- ZetaXi[["xi"]]
    DataStorage@presentation$zeta <- rep(zeta, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    DataStorage <- DPsurv::gibbsStep(DP=NDP, DataStorage=DataStorage, 
                                     xi=xi, zeta=rep(zeta, each = max(DataStorage@mask)))
    NDP@prior <- DPsurv::mStep(NDP@prior, DataStorage@simulation, xi=xi, zeta=rep(zeta, each = max(DataStorage@mask)))
    NDP <- update.NDP(NDP)
    if(NDP@details[["iteration"]]>NDP@details[["burnin"]] & (NDP@details[["iteration"]] %% NDP@details[["thinning"]])==0){
      setTxtProgressBar(pb, j/iter)
      NDP@Chains <- list(theta=NDP@theta, phi=NDP@phi, weights=NDP@weights, pi=NDP@pi)
      NDP@ChainStorage <- saveChain.ChainStorage(unique(zeta), 1:NDP@L, NDP@Chains, (NDP@details[["iteration"]]-NDP@details[["burnin"]])/NDP@details[["thinning"]], NDP@ChainStorage)
    }
  }
  close(pb)
  return(NDP)
}

#'
#' @export
MCMC.NDP_c <- compiler::cmpfun(MCMC.NDP)

#'
#' @export
posteriorZeta.NDP <- function(NDP, DataStorage){
  NDP <- posterior.DP(NDP, 0.5)
  ZetaXi <- selectZetaXi.NDP(NDP, DataStorage, max_lik = T)
  zeta <- ZetaXi[["zeta"]]
  DataStorage@presentation$zeta <- rep(zeta, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
  DataStorage@validation$zeta <- as.numeric(as.character(plyr::mapvalues(DataStorage@validation$Sample, DataStorage@presentation$Sample, DataStorage@presentation$zeta,  warn_missing=F)))
  return(DataStorage)
}

#'
#' @export
plotICDF.NDP <- function(NDP, DataStorage){
  DataStorage <- posteriorZeta.NDP(NDP, DataStorage)
  for(zeta in unique(DataStorage@presentation$zeta)){
    p <- plot.ICDF(NDP, zeta, DataStorage@presentation) + 
      ggplot2::ggtitle(paste("NDP", zeta))
    print(p)
  }
}

#'
#' @export
validate.NDP <- function(NDP, DataStorage){
  DataStorage <- posteriorZeta.NDP(NDP, DataStorage)
  zeta <- unique(DataStorage@presentation$zeta)
  medianCurves <- getICDF.ChainStorage(NDP, DataStorage@validation$data, zeta=zeta)
  medianCurvesArranged <- matrix(NA, dim(medianCurves)[1], dim(NDP@ChainStorage@chains[["theta"]])[2])
  for(i in 1:length(zeta)) medianCurvesArranged[,zeta[i]] <- medianCurves[,i]
  score <- validate(curves=medianCurvesArranged, status=DataStorage@validation$status, zeta=DataStorage@validation$zeta)
  return(score)
}

