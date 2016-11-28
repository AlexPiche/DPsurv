#'
#' @export
plot.PDF <- function(theta, phi, weights, L, grid, emp_data, xlim, ...){
  #emp_data <- ifelse(rep(survival, length(emp_data)), exp(emp_data), emp_data)
  
  #grid <- seq(2, xlim, 1)
  
  myPDF <- evaluate.PDF(theta, phi, weights, L, grid)
  
  mydensity <- data.frame(x = grid, y = myPDF)
  p <- ggplot2::ggplot() + ggplot2::geom_line(data=mydensity, ggplot2::aes(x=x,y=y), colour='blue')+ ggplot2::theme_linedraw()
  if(length(emp_data)>1){
    p <- p + ggplot2::geom_density(data=data.frame(x=(emp_data)), ggplot2::aes(x=x), linetype='dashed') + ggplot2::theme()
  }
 # p <- p + ggplot2::ggtitle(paste("DP", DPnb, "with probability:", round(theObject@toPlot_pi[[DPnb]], 3), sep=" ")) +
#    ggplot2::xlim(0, xlim)
  return(p) 
}

#'
#' @export
plot.ICDF <- function(DP, mySample, data, ...){
  data_pres <- subset(data@presentation, Sample == mySample)
  data_pres$status <- as.numeric(as.character(data_pres$status))
  temp <- data_pres[,c("status", "data", "Sample")]
  
  xlim <- 1.5*max(temp$data)
  grid <- subset(data@validation, Sample == mySample)$data
  
  x_cens <- subset(temp, status == 0)$data
  curves <- data.frame(grid, t(getMedianCurves(DP, data, J=mySample)))
  names(curves) <- c("grid", "0.025", "0.5", "0.975")
  curves <- reshape2::melt(curves, "grid")
  KM <- data.frame(cbind(grid, "True", 1-subset(data@validation, Sample == mySample)$status))
  names(KM) <- names(curves)
  KM[, c(1,3)] <- lapply(KM[, c(1,3)], as.character)
  KM[, c(1,3)] <- lapply(KM[, c(1,3)], as.numeric)
  curvestoPlot <- rbind(KM, curves)
  curves$variable <- factor(curves$variable)
  
  p <- ggplot2::ggplot(curvestoPlot, ggplot2::aes(x=grid, y=value, group=variable, colour=variable)) + ggplot2::geom_line() +
    ggplot2::xlim(0, xlim)+ ggplot2::theme_linedraw() +
    ggplot2::xlab("Time") + ggplot2::ylab("Survival") + ggplot2::theme(legend.title = ggplot2::element_blank())  
  p <- p + ggplot2::scale_color_manual(values=c("black","lightskyblue2", 'blue', "lightskyblue2")) + ggplot2::scale_linetype_manual(values = c("dashed", rep("solid", 3)))
  p <- p + ggplot2::ggtitle(paste(class(DP)[1], mySample)) + ggplot2::geom_vline(xintercept = x_cens, linetype=2)
  print(p)
  return(p)
}




#'
#' @export
grid_arrange_shared_legend <- function(plots, ...) {
  g <- ggplot2::ggplotGrob(plots[[1]] + ggplot2::theme(legend.position="right", legend.key.height = ggplot2::unit(1.5, "cm"),
                                     legend.key.width = ggplot2::unit(0.5, "cm"),legend.title = ggplot2::element_blank()))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  gridExtra::grid.arrange(
    do.call(gridExtra::arrangeGrob, lapply(plots, function(x)
      x + ggplot2::theme(legend.position="none"))),
    legend,
    ncol = 2,
    widths=grid::unit.c(grid::unit(5.5,"in"), grid::unit(1,"in")))
  #heights = grid::unit.c(unit(10,"npc"), unit(1,"npc")))#unit(1, "npc") - lheight, lheight))
}

#'
#' @export
plotICDF.DP <- function(DP, DataStorage, samples=0){
  if(samples[1]==0) samples <- 1:DP@J
  for(sample in samples){
    plot.ICDF(DP=DP, mySample=sample, data=DataStorage) 
  }
}

#'
#' @export
plotHeatmap <- function(DP){
  myLevels <- 1:DP@J
  if(class(DP)[1] == "HDP"){
    Z <- DP@ChainStorage@chains[["weights"]]
  }else{
    Z <- DP@ChainStorage@chains[["prob"]]
  }
  zz <- matrix(aperm(Z, c(1,3,2)), ncol=ncol(Z))
  
  crm1 <- reshape2::melt(cor(zz))#t(matrix(zz, ncol=ncol(zz)))))
  crm1$Var1 <- factor(crm1$Var1, levels = myLevels, ordered=T)
  crm1$Var2 <- factor(crm1$Var2, levels = myLevels, ordered=T)
  
  
  p <- ggplot2::ggplot(crm1, ggplot2::aes(Var1, Var2)) + ggplot2::geom_tile(ggplot2::aes(fill=value)) +
    ggplot2::scale_fill_gradient(low="white", high="steelblue") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6),
                   axis.text.y = ggplot2::element_text(size = 6),
                   axis.title.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.key.height = ggplot2::unit(2.5, "cm"),
                   legend.box.just = "right")
  print(p)
}