#'Hierarchical Dirichlet Process
#'HDP applied to censored data
#'@examples
#'\dontrun{
#'weights <- matrix(c(1,0,0,0,1,0,0,0,1), ncol=3)
#'data <- sim.data(n=100, J=10, weights)
#'G3 <- new("HDP")
#'G3 <- init.HDP(G3, prior=list(mu=0, n=0.1, v=3, vs2=1*3), L=15, 
#'                    J=length(unique(data@presentation$Sample)), thinning=20,
#'                    burnin = 5000, max_iter = 55000)
#'G3 <- MCMC.HDP(G3, data, 55000)
#'plot.ICDF(G3@theta, G3@phi, G3@weights[,1], G3@L, grid=0:500,
#'            distribution=data@presentation, xlim=500)
#'validate.HDP(G3, data)
#'}
#'
#' @export
setClass("HDP", representation(phi = 'matrix', theta='matrix', weights='matrix', Nmat='matrix', details='list', conc_param = 'numeric',
                              prior = 'array', L = 'numeric', J = 'numeric', Chains='list', pi='matrix', ChainStorage='ChainStorage'))

#'
#' @export
init.HDP <- function(HDP, prior, L, J, thinning, burnin, max_iter, ...){
  HDP@L <- L
  HDP@J <- J
  HDP@prior <- array(rep(c(prior$mu, prior$n, prior$v, prior$vs2), each=(L*J)), c(L,J,4))
  HDP@Nmat <- matrix(0, nrow=L, ncol=J)
  HDP@details <- list(iteration=0, thinning=thinning, burnin=burnin, max_iter=max_iter)
  HDP@conc_param <- c(1,1)#rgamma(2,1,1)
  HDP <- update.HDP(HDP)
  HDP@Chains <- list(theta=HDP@theta, phi=HDP@phi, weights=HDP@weights, pi=HDP@pi)
  HDP@ChainStorage <- init.ChainStorage(HDP@Chains, max_iter-burnin, thinning)
  return(HDP)
}

#'
#' @export
update.HDP <- function(HDP, ...){
  #print(HDP@conc_param)
  atoms <- rNIG(HDP@L, c(HDP@prior[,,1]), c(HDP@prior[,,2]), c(HDP@prior[,,3]), c(HDP@prior[,,4]))
  HDP@theta <- matrix(atoms[,1], nrow=HDP@L)
  HDP@phi <- matrix(atoms[,2], nrow=HDP@L)
  # 1st level
  sums <- mySums(round(c(HDP@prior[,1,2])))
  u_k <- rbeta(HDP@L, shape1 = 1 + c(HDP@prior[,1,2]), shape2 = HDP@conc_param[1] + sums)
  HDP@pi <- as.matrix(stickBreaking(u_k))
  
  # 2nd level
  sumsNmat <- matrix(apply(HDP@Nmat, 2, mySums), nrow=dim(HDP@Nmat)[1])
  alpha <- 1
  sums <- sapply(1:HDP@L, function(i) {alpha*(1-sum(HDP@pi[1:i]))})
  sums[which(sums < 0)] <- 0
  
  v_lk <- rbeta(HDP@J*HDP@L, HDP@conc_param[2]*rep(HDP@pi, HDP@J)+c(HDP@Nmat), c(sumsNmat)+HDP@conc_param[2]*rep(sums, HDP@J))
  v_lk <- matrix(v_lk, ncol = HDP@J, byrow = F)
  
  HDP@weights <- apply(v_lk, 2, stickBreaking)
  #a_conc <- 1
  #b_conc <- 1
  #HDP@conc_param[1] <- rgamma(1, a_conc + HDP@L - 1, b_conc - sum(log(1-u_k[1:HDP@L-1])))
  #HDP@conc_param[2] <- rgamma(1, a_conc + HDP@J*(HDP@L-1), b_conc - sum(log(1-v_lk[-c(seq(HDP@L, HDP@J*HDP@L, HDP@L), which(v_lk==1))])))
  return(HDP)
}

#'
#' @export
computeXi.HDP <- function(HDP, DataStorage, max_lik=F, ...){
  xx <- DPsurv::eStep(theta=HDP@theta, phi=HDP@phi, w=rep(1, length(HDP@phi)), DataStorage=DataStorage)
  xx <- stabilize(xx)
  zz <- array(xx, c(HDP@L, max(DataStorage@mask), length(DataStorage@mask)))
  
  for(i in 1:dim(zz)[3]){
    zz[,,i] <- exp(log(zz[,,i]) + log(HDP@weights[,i]))
  }
  
  rr <- matrix(zz, nrow=HDP@L)
  xi <- apply(rr, 2, DP_sample, n=HDP@L, size=1, replace=F, max_lik=max_lik)
  return(xi)
}

#'
#' @export
MCMC.HDP <- function(HDP, DataStorage, iter, ...){
  i <- 0
  while(HDP@details[['iteration']] < HDP@details[['max_iter']] & i < iter){

    HDP@details[['iteration']] <- HDP@details[['iteration']] + 1
    i <- i + 1
    xi <- computeXi.HDP(HDP, DataStorage)
    DataStorage@presentation$zeta <- rep(1:HDP@J, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
    Nmat <- matrix(factor(DataStorage@presentation$xi, levels = 1:HDP@L), ncol=HDP@J)
    Nmat <- apply(Nmat, 2, as.numeric)
    HDP@Nmat <- apply(Nmat, 2, test, L = HDP@L)
    DataStorage <- DPsurv::gibbsStep(DP=HDP, DataStorage=DataStorage, xi=xi, zeta=rep(1, length(xi)))
    HDP <- DPsurv::mStep(HDP, DataStorage, xi, rep(1, length(xi))) 
    HDP <- update.HDP(HDP)
    HDP@Chains <- list(theta=HDP@theta, phi=HDP@phi, weights=HDP@weights, pi=HDP@pi)
    if(HDP@details[["iteration"]] > HDP@details[["burnin"]] & (HDP@details[["iteration"]] %% HDP@details[["thinning"]])==0){
      HDP@ChainStorage <- saveChain.ChainStorage(HDP@Chains, (HDP@details[["iteration"]]-HDP@details[["burnin"]])/HDP@details[["thinning"]], HDP@ChainStorage)
    }
  }
  HDP <- posterior.DP(HDP, 0.5)
  return(HDP)
}

#'
#' @export
validate.HDP <- function(HDP, DataStorage){
  xi <- computeXi.HDP(HDP, DataStorage, max_lik=T)
  DataStorage@presentation$xi <- xi[!is.na(xi)]
  DataStorage@presentation$zeta <- rep(1:HDP@J, as.vector(table(DataStorage@presentation$Sample, useNA = "no")))
  mapping <- unique(DataStorage@presentation[, c("zeta", "Sample")])
  DataStorage@validation$zeta <- plyr::mapvalues(DataStorage@validation$Sample, mapping$Sample, mapping$zeta, warn_missing=F)
  score <- validate(data=DataStorage@validation$data, status=DataStorage@validation$status, zeta=as.numeric(DataStorage@validation$zeta),#zeta=DataStorage@validation$xi,
                    theta=matrix(rep(c(HDP@theta), HDP@J), ncol=HDP@J), phi=matrix(rep(c(HDP@phi), HDP@J), ncol=HDP@J), weights=HDP@weights)
  return(score)
}

#'
#' @export
test <- function(vec, L){
  vec <- factor(vec, 1:L)
  table(vec, useNA = 'no')
}

