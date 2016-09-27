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
plot.ICDF <- function(DP, myZeta, data, ...){
  data <- subset(data, zeta == myZeta)
  if(sum(data$status) == 0) return(-1) # cannot be graph
  temp <- data[,c("status", "data", "Sample")]
  
  xlim <- 1.25*max(temp$data)
  grid <- seq(0, xlim, length.out = min(xlim, 1500))
  
  if(dim(temp[temp$status==1,])[1] > 1){
    s <- with(temp, survival::Surv(data, status))
    LN1 <- survival::survfit(s ~ 1, data=temp)
    xx <- summary(LN1)
    KM <- data.frame(Var1="KM", Var2=xx$time, value=xx$surv)
    KM <- rbind(data.frame(Var1="KM", Var2=0, value=1), KM)
  }else{
    KM <- data.frame(Var1=character(), Var2=numeric(), value=numeric())
  }
  
  curves <- getICDF.ChainStorage(DP, grid, myZeta, quantiles=c(0.05, 0.5, 0.95))
  curves <- reshape2::melt(curves)
  curves$Var2 <- rep(grid, each=3)
  KM$Var3 <- 1
  curves <- rbind(KM, curves)
  
  p <- ggplot2::ggplot(curves, ggplot2::aes(x=Var2, y=value, group=Var1, colour=Var1)) + ggplot2::geom_line() +
    ggplot2::xlim(0, xlim)+ ggplot2::theme_linedraw() +
    ggplot2::xlab("Time") + ggplot2::ylab("Survival") + ggplot2::theme(legend.title = ggplot2::element_blank())  
  p <- p + ggplot2::scale_color_manual(values=c("black","lightskyblue2", 'blue', "lightskyblue2")) + ggplot2::scale_linetype_manual(values = c("dashed", rep("solid", 3)))
  p <- p + ggplot2::ggtitle(paste(class(DP)[1], myZeta))
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
plotICDF.DP <- function(DP, DataStorage, zetas=1:DP@J){
  for(zeta in zetas){
    plot.ICDF(DP=DP, myZeta=zeta, data=DataStorage@presentation) 
  }
}

#'
#' @export
plotHeatmap <- function(DP, DataStorage){
  
  myLevels <- unique(DataStorage@presentation$zeta)
  zz <- getICDF.ChainStorage(DP, grid, myLevels, quantiles=c(0.05, 0.5, 0.95))
  browser()
  #zz <- theObject@HM[,which(!apply(theObject@HM,2,FUN = function(x){all(x == 0)}))]
  
  crm1 <- reshape2::melt(cor(t(zz)))
  
  mapping <- data.frame(from=1:length(myLevels),to=unique(theObject@T["Sample"]))#, stringsAsFactors=FALSE)
  names(mapping) <- c("from", "to")
  mapping$to <- lapply(mapping$to, as.character)
  crm1$Var1 <- plyr::mapvalues(crm1$Var1, mapping$from, mapping$to)
  crm1$Var2 <- plyr::mapvalues(crm1$Var2, mapping$from, mapping$to)
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
}