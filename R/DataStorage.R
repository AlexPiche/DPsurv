#'
#' @export
setClass("DataStorage", representation(presentation='data.frame', computation='numeric', simulation='numeric', validation='data.frame',
                                       grid='numeric', mask='numeric', censoring='numeric', episode='matrix'))#, xi='numeric'))



#'
#' @export
train_test_split <- function(dataset, fraction=0.1, DataStorage){
  idx <- sample.int(dim(dataset)[1], round(fraction*dim(dataset)[1]))
  validation <- dataset[sort(idx),]
  DataStorage@validation <- validation
  dataset <- dataset[-sort(idx),]
  DataStorage@presentation <- dataset
  return(DataStorage)
}

#'
#' @export
init.DataStorage.simple <- function(dataset, fraction_test, ...){
  dataset <- plyr::arrange(dataset, Sample)
  DataStorage <- new("DataStorage")
  DataStorage <- train_test_split(dataset, fraction_test, DataStorage)
  X <- lapply(unique(DataStorage@presentation$Sample), function(ss) t(matrix(subset(DataStorage@presentation, Sample == ss)$status)))
  censoring <- do.call(plyr::rbind.fill.matrix, X)
  DataStorage@censoring <- c(t(censoring))
  X <- lapply(unique(DataStorage@presentation$Sample), function(ss) t(matrix(subset(DataStorage@presentation, Sample == ss)$data)))
  computation <- do.call(plyr::rbind.fill.matrix, X)
  # log the data for computation purposes
  DataStorage@computation <- log(c(t(computation)))
  DataStorage@simulation <- DataStorage@computation
  max_grid <- ceiling(round(1.5*max(dataset$data)))
  DataStorage@grid <- seq(0, max_grid, round(max_grid/250))
  DataStorage@mask <-rowSums(!is.na(computation))
  #DataStorage@xi <- rep(0, 2)
  return(DataStorage)
}

init.DataStorage.recurrent <- function(dataset, ...){
  #X <- lapply(unique(dataset$Sample), function(ss) t(matrix(subset(dataset, Sample == ss)$real.episode)))
  #DataStorage@episode <- do.call(plyr::rbind.fill.matrix, X)
  return(-1)
}

findRt <- function(x){
  xx <- survsim::simple.surv.sim(1000, anc.ev = 3, foltime = 10000, anc.cens = 3, beta0.cens = x, beta0.ev = 4.5, dist.ev = "weibull")[,c('stop', 'status')]
  return(sum(xx$status)/1000-0.85)
}

#'
#' @export
simSurvMix <- function(N, prob, factor){
  # 15% censoring
    toRet <- data.frame(data=rep(NA,N), status=rep(NA,N))
    for (n in 1:N){
      i <- runif(1)
      if(i < prob[1]){
        # -log(1/uniroot(function(x) pweibull(x,a1,b1)-0.85, interval=c(0,exp(150)))$root)
        toRet[n,] <- survsim::simple.surv.sim(1, anc.ev = 1, foltime = 10000, anc.cens = 1, beta0.cens = 5.522906, beta0.ev = 3.75, dist.ev = "weibull")[,c('stop', 'status')]
      }else if(i < sum(prob[1:2])){
        # log(uniroot(function(x) plnorm(x,3.887,1)-0.85, interval=c(0,exp(15)))$root)
        toRet[n,] <- survsim::simple.surv.sim(1, anc.ev = 1, foltime = 10000, anc.cens = 1, beta0.cens = 5.336191, beta0.ev = 3.887, dist.ev = "lnorm", dist.cens = "lnorm")[,c('stop', 'status')]
      }else if(i < sum(prob[1:3])){
        #-log(1/(uniroot(function(x) pweibull(x,a3,b3)-0.85, interval=c(0,exp(150)))$root)^3)
        toRet[n,] <- survsim::simple.surv.sim(1, anc.ev = 3, foltime = 10000, anc.cens = 3, beta0.cens = 5.091579, beta0.ev = 4.5, dist.ev = "weibull")[,c('stop', 'status')]
      }else{
        toRet[n,] <- survsim::simple.surv.sim(1, anc.ev = 1, foltime = 10000, anc.cens = 1, beta0.cens = factor*8, beta0.ev = 8, dist.ev = "weibull")[,c('stop', 'status')]
      }
    }
    toRet
}

#'
#' @export
simRecSurvMix <- function(){
#  sim.data <- rec.ev.sim(n=5000, foltime=1825, dist.ev=c('weibull','weibull'),
#                         anc.ev=c(1, 1),beta0.ev=c(6.678, 4.430),,
#                       anc.cens=c(1, 1), beta0.cens=c(6.712, 4.399))
  return(-1)
}

#'
#' @export
getGrid <- function(weights){
  a1 <- 1
  b1 <- (1/exp(-1 * 3.75))^(1/1)
  a2 <- 1
  b2 <- factor*3.887
  a3 <- 3
  b3 <- (1/exp(-3 * (4.5)))^(1/3)
  myFct <- function(x, weights, pt){
    weights[1]*pweibull(x,a1,b1)+weights[2]*plnorm(x,b2,a2)+weights[3]*pweibull(x,a3,b3) - pt
  }
  rootT1 <- uniroot(function(x) myFct(x, weights=c(weights)[1:3], pt=1), interval = c(-1,5))$root
}

#'
#' @export
sim.data <- function(weights, n, J, validation_prop=0.1, factor=1){
  N <-  J*n
  
  T1 <- simSurvMix(N, c(weights)[1:3], factor=factor)
  T2 <- simSurvMix(N, c(weights)[4:6], factor=factor)
  T3 <- simSurvMix(N, c(weights)[7:9], factor=factor)
  

  
  T1$TrueDistribution <- "T1"
  T2$TrueDistribution <- "T2"
  T3$TrueDistribution <- "T3"
  
  mixture <- rbind(T1,T2,T3)
  
  mixture$TrueDistribution <- as.factor(mixture$TrueDistribution)
  mixture$xi <- rep(0, 3*N)
  mixture$zeta <- rep(0, 3*N)
  mixture$Sample <- paste("S", rep(1:(3*J), each=n),sep='')
  mixture$Sample <- as.factor(mixture$Sample)
  mixture$data <- mixture$data
  
  #names(mixture)[3] <- "data"
  data <- init.DataStorage.simple(mixture, validation_prop)
  data
}

#'
#' @export
update_validation_set <- function(DataStorage, ...){
  mapping <- unique(DataStorage@presentation[, c("xi", "Sample")])
  DataStorage@validation$xi <- plyr::mapvalues(DataStorage@validation$Sample, mapping$Sample, mapping$xi, warn_missing=F)
  #DataStorage@validation$xi[which(DataStorage@validation$xi == DataStorage@validation$Sample)] <- NA
  #DataStorage@validation$xi <- as.numeric(DataStorage@validation$xi)
  return(DataStorage)
}
