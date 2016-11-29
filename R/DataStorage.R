#'
#' @export
setClass("DataStorage", representation(presentation='data.frame', computation='numeric', simulation='numeric', validation='data.frame',
                                       grid='numeric', mask='numeric', censoring='numeric', episode='matrix'))#, xi='numeric'))


#'
#' @export
getGrid <- function(weights, myPt){
  a1 <- 1
  b1 <- (1/exp(-1 * 3.75))^(1/1)
  a2 <- 1
  b2 <- 3.887
  a3 <- 3
  b3 <- (1/exp(-3 * (4.5)))^(1/3)
  myFct <- function(x, weights, pt){
    weights[1]*pweibull(x,a1,b1)+weights[2]*plnorm(x,b2,a2)+weights[3]*pweibull(x,a3,b3) - pt
  }
  sln <- uniroot(function(x) myFct(x, weights=weights, pt=myPt), interval = c(0,500))$root
  return(sln)
}

getICDF.DataStorage <- function(x, weights, params=c(3.75, 3.887, 4.5)){
  a1 <- 1
  b1 <- (1/exp(-1 * params[1]))^(1/1)
  a2 <- 1
  b2 <- params[2]
  a3 <- 3
  b3 <- (1/exp(-3 * (params[3])))^(1/3)
  toRet <- weights[1]*pweibull(x,a1,b1)+weights[2]*plnorm(x,b2,a2)+weights[3]*pweibull(x,a3,b3)
  return(1-toRet)
}

#'
#' @export
train_test_split <- function(dataset, DataStorage, J, fraction=0, weights=rep(0,9)){
  if(fraction > 0){
    idx <- sample.int(dim(dataset)[1], round(fraction*dim(dataset)[1]))
    validation <- dataset[sort(idx),]
    DataStorage@validation <- validation
    dataset <- dataset[-sort(idx),]
  }else{
    #J <- length(unique(DataStorage@presentation$Sample))/3
    myQuantiles <- seq(0.05, 0.95, 0.1)
    S1 <- rep(mapply(getGrid, myQuantiles, MoreArgs = list(weights=c(weights)[1:3])), J)
    S2 <- rep(mapply(getGrid, myQuantiles, MoreArgs = list(weights=c(weights)[4:6])), J)
    S3 <- rep(mapply(getGrid, myQuantiles, MoreArgs = list(weights=c(weights)[7:9])), J)
    DataStorage@validation <- data.frame(data=c(S1,S2,S3), 
                                         status=rep(myQuantiles, 3*J), 
                                         Sample=rep(1:(3*J), each=length(myQuantiles)))
  }
  DataStorage@presentation <- dataset
  return(DataStorage)
}

#'
#' @export
init.DataStorage.simple <- function(dataset, fraction_test, weights, J, application=F, ...){
  dataset <- plyr::arrange(dataset, Sample)
  DataStorage <- methods::new("DataStorage")
  if(!application) {
    DataStorage <- train_test_split(dataset=dataset, DataStorage=DataStorage, J=J, fraction=fraction_test, weights=weights)
  }else{
    DataStorage@presentation <- dataset
    DataStorage@validation <- dataset
  }
  X <- lapply(unique(DataStorage@presentation$Sample), function(ss) t(matrix(subset(DataStorage@presentation, Sample == ss)$status)))
  censoring <- do.call(plyr::rbind.fill.matrix, X)
  DataStorage@censoring <- as.numeric(c(t(censoring)))
  X <- lapply(unique(DataStorage@presentation$Sample), function(ss) t(matrix(subset(DataStorage@presentation, Sample == ss)$data)))
  computation <- do.call(plyr::rbind.fill.matrix, X)
  # log the data for computation purposes
  DataStorage@computation <- log(c(t(computation)))
  DataStorage@simulation <- DataStorage@computation
  max_grid <- ceiling(round(2*max(dataset$data)))
  DataStorage@grid <- seq(0, max_grid, ceiling(max_grid/250))
  DataStorage@mask <-rowSums(!is.na(computation))
  #DataStorage@xi <- rep(0, 2)
  return(DataStorage)
}



findRt <- function(x){
  xx <- survsim::simple.surv.sim(1000, anc.ev = 3, foltime = 10000, anc.cens = 3, beta0.cens = x, beta0.ev = 4.5, dist.ev = "weibull")[,c('stop', 'status')]
  return(sum(xx$status)/1000-0.5)
}

#'
#' @export
simSurvMix <- function(N, prob, frailty=0, params=c(3.75, 3.887, 4.5)){
  factor <- 1
  foltime <- 10000
  toRet <- data.frame(data=rep(NA,N), status=rep(NA,N))
  for (n in 1:N){
    i <- runif(1)
    if(i < prob[1]){
      # -log(1/uniroot(function(x) pweibull(x,,b1)-0.85, interval=c(0,exp(150)))$root)
      toRet[n,] <- survsim::simple.surv.sim(1, anc.ev = 1, foltime = foltime, anc.cens = 1, beta0.cens = params[1]+log(frailty), beta0.ev = params[1]+log(frailty), dist.ev = "weibull")[,c('stop', 'status')]
    }else if(i < sum(prob[1:2])){
      # log(uniroot(function(x) plnorm(x,3.887,1)-0.85, interval=c(0,exp(15)))$root)
      toRet[n,] <- survsim::simple.surv.sim(1, anc.ev = 1, foltime = foltime, anc.cens = 1, beta0.cens = params[2]+log(frailty), beta0.ev = params[2]+log(frailty), dist.ev = "lnorm", dist.cens = "lnorm")[,c('stop', 'status')]
    }else if(i < sum(prob[1:3])){
      #-log(1/(uniroot(function(x) pweibull(x,a3,b3)-0.85, interval=c(0,exp(150)))$root)^3) 
      toRet[n,] <- survsim::simple.surv.sim(1, anc.ev = 3, foltime = foltime, anc.cens = 3, beta0.cens = params[3]+log(frailty), beta0.ev = params[3]+log(frailty), dist.ev = "weibull")[,c('stop', 'status')]
    }else{
      toRet[n,] <- survsim::simple.surv.sim(1, anc.ev = 1, foltime = 10000, anc.cens = 1, beta0.cens = factor*8, beta0.ev = 8, dist.ev = "weibull")[,c('stop', 'status')]
    }
  }
  toRet
}




#'
#' @export
sim.data <- function(weights, n, J, validation_prop=0, frailty=F){
  N <-  J*n
  mixture<- data.frame(data=numeric(), status=numeric())
  idx <- apply(matrix(1:9, ncol=3), 1, rep, each=J)
  for(i in 1:(3*J)){
    fr <- ifelse(frailty, rgamma(1,1,1), 1)
    temp <- simSurvMix(n, c(weights)[idx[i,]], frailty=1)
    temp$data <- fr * temp$data
    mixture <- rbind(mixture, temp)
  }
  
  mixture$xi <- rep(0, 3*N)
  mixture$zeta <- rep(1:(3*J), each=n)
  mixture$Sample <- rep(1:(3*J), each=n)
  
  data <- init.DataStorage.simple(dataset=mixture, fraction_test = validation_prop, weights=weights, J=J)
  return(data)
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
