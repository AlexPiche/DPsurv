#' @export
setClass("RE", representation(RE='matrix', computation = 'matrix', n = 'numeric', individual='matrix',
                                         prior = 'array'))

#we need one per individual, then map it to a matrix

#' @export
init.RE <- function(DataStorage, ...){
  RE <- new("RE")
  DataStorage@presentation$individual <- 1:dim(DataStorage@presentation)[1]
  X <- lapply(unique(DataStorage@presentation$Sample), function(ss) t(matrix(subset(DataStorage@presentation, Sample == ss)$individual)))
  RE@individual <- do.call(plyr::rbind.fill.matrix, X)
  RE@individual <- matrix(as.numeric(RE@individual), ncol=ncol(RE@individual))
  RE@n <- length(unique(DataStorage@presentation$individual))
  prior <- list(mu=0, n=0.1, v=3, vs2=1*3)
  RE@prior <- array(rep(c(prior$mu, prior$n, prior$v, prior$vs2), each=(RE@n)), c(RE@n,1,4))
  RE <- update.RE(RE)
  return(RE)
}

#' @export
update.RE <- function(RE, ...){
  atoms <- rNIG(RE@n, c(RE@prior[,,1]), c(RE@prior[,,2]), c(RE@prior[,,3]), c(RE@prior[,,4]))
  mapping <- data.frame(ind=1:length(atoms[,1]), atoms=atoms[,1])
  RE@RE <- matrix(atoms[,1])
  if(T) mapping$atoms <- rep(0, length(atoms[,1]))
  temp <- plyr::mapvalues(c(RE@individual), mapping$ind, mapping$atoms, warn_missing=F)#
  RE@computation <- matrix(as.numeric(temp), ncol=ncol(RE@individual))
  return(RE)
}

#' @export
MCMC.RE <- function(RE, simulation){
  RE@prior <- mStep(RE@prior, data=simulation, xi=c(t(RE@individual)), zeta=rep(1, length(simulation)))
  RE <- update.RE(RE)
  return(RE)
}
