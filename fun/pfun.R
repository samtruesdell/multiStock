# plotting function
pfun <- function(x, v=NULL, egX, nm, rgy, rgx, txtxp=0.1, txtyp=0.9,
                 xaxt='s', yaxs, leg=FALSE, legtxt=NULL, legpos=NULL,
                 labs=NULL){
  # if plotting against the scalar, turn it into a matrix
  # of identical rows
  if(!is.matrix(egX)){
    egX <- matrix(data=egX, nrow=dim(x)[1], ncol=dim(x)[2],
                  byrow=TRUE)
  }
  nms <- matrix(c('esc', 'harv', 'harvR', 'over', 'ext', 'pFailSub', 
                  'cvCatch', 'HRsurplus', 'Escapement', 'Harvest', 'Harvest Rate', 
                  'P(Overfishing)', 'P(Extirpation)',
                  'P(Fail Meet Sub. Needs)', 'CV Catch',
                  'Harv. Rate on Surplus'), ncol=2)
  if(!any(nms[,1]==nm))
    stop('check name')
  nmind <- which(nm == nms[,1])
  
  # establish colors
  nscen <- dim(x)[1]
  cr <- colorRampPalette(c('steelblue', 'firebrick1'))
  col <- cr(nscen)
  
  # set up reference line intervals
  if(rgy[2] < 10){
    rh <- seq(0, 10, 0.2)
  }else{
    rh <- seq(0, round(rgy[2], -3)+1000, 10000)
  }
  rv <- seq(0, round(rgx[2], -2)+10000, 10000)
  
  plot(NA, ylim=rgy, xlim=rgx, xlab='', ylab='',
       xaxt=xaxt, yaxt='n', type='n')
  abline(h=rh, v=rv, lty=3, col='gray80')
  for(j in 1:nscen){
    lines(x[j,,nmind] ~ egX[j,], col=col[j])
    if(!is.null(v)){
      mpv <- x[j,,nmind] + v[j,,nmind]
      mmv <- x[j,,nmind] - v[j,,nmind]
      lines(mpv ~ egX, lty=2, col=col[j])
      lines(mmv ~ egX, lty=2, col=col[j])
    }
  }
  axis(yaxs, las=1)
  # txty <- txtyp * (rgy[2]-rgy[1]) + rgy[1]
  # txtx <- txtxp * (rgx[2]-rgx[1]) + rgx[1]
  # text(txtx, txty, labels=nms[nmind], cex=2, pos=4)
  title(main=nms[nmind,2], cex=1.5)
  if(leg)
    legend(legpos, col=col, lty=1, 
           legend=paste(legtxt, labs), bty='n')
}
