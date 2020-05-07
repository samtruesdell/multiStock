### mType plot matrix
pm <- function(x){
  
  xc <- 1:ncol(x)
  yc <- rev(1:nrow(x))
  plot(0, type='n', xlim=c(0, ncol(x)+1), ylim=c(0, nrow(x)+1), las=1,
       xlab='', ylab='', yaxt='n')
  axis(2, at=pretty(yc), labels=rev(pretty(yc)), las=1)
  mtext(side=1:2, line=3, cex=1.5,
        c('Stock', 'Scenario'))
  col <- c('gray', 'cornflowerblue', 'firebrick1')
  for(i in 1:length(xc)){
    rect(xleft=xc[i]-0.5, ybottom=yc-0.5, xright=xc[i]+0.5, ytop=yc+0.5,
         col=col[x[,i]+1])
  }
  
}