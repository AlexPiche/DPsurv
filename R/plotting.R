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
plot.ICDF <- function(theta, phi, weights, L, grid, xlim, distribution, ...){
  
  temp <- distribution[,c("status", "data", "Sample")]
  
  if(dim(temp[temp$status==1,])[1] > 1){
    #browser()
    s <- with(temp, survival::Surv(exp(data), status))
    LN1 <- survival::survfit(s ~ 1, data=temp)
    xx <- summary(LN1)
    KM <- data.frame(Var1="KM", Var2=xx$time, value=xx$surv)
    KM <- rbind(data.frame(Var1="KM", Var2=0, value=1), KM)
  }else{
    KM <- data.frame(Var1=character(), Var2=numeric(), value=numeric())
  }
  
  #if(posterior){
  #  weight <- (theObject@iteration - theObject@burnin)/theObject@thinning
  #  toPlot <- reshape2::melt(apply(theObject@MCMC_ICDF[,DPnb,1:weight], 1, quantile, c(0.025,0.5,0.975), na.rm=T))
  #  toPlot$Var2 <- rep(theObject@grid,each=3)
  #}else{
  toPlot <- evaluate.ICDF(theta, phi, weights, grid)
  toPlot <- data.frame(toPlot, grid)
  names(toPlot) <- c("value", "Var2")
  toPlot$Var1 <- "Estimation"
  #}
  if(dim(temp[temp$status==1,])[1] > 3) toPlot <- rbind(KM, toPlot)
  p <- ggplot2::ggplot(toPlot, ggplot2::aes(x=Var2, y=value, group=Var1, colour=Var1)) + ggplot2::geom_line() +
    ggplot2::xlim(0, xlim)+ ggplot2::theme_linedraw() +
    ggplot2::xlab("Time") + ggplot2::ylab("Survival")+ ggplot2::theme(legend.position="none")# + ggplot2::theme(legend.title = ggplot2::element_blank())
  #if(posterior){
  #  if(dim(temp[temp$status==1,])[1] > 3){
  #    p <- p + ggplot2::scale_color_manual(values=c("black","lightskyblue2", 'blue', "lightskyblue2")) + ggplot2::scale_linetype_manual(values = c("dashed", rep("solid", 3)))
  #  }else{
  #    p <- p + ggplot2::scale_color_manual(values=c("lightskyblue2", 'blue', "lightskyblue2"))
  #  }
  #}
  return(p)
}