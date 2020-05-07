# plot of ratio (x) vs harvest etc (y) for the different
# numbers of stocks
ratMetPlot <- function(nm, sfun){
  
  nms <- matrix(c('esc', 'harv', 'harvR', 'over', 'ext', 'pFailSub', 
                  'cvCatch', 'HRsurplus', #end abbrevs
                  'Escapement', 'Harvest', 
                  'Harvest Rate', 
                  'P(Overfishing)', 'P(Extirpation)',
                  'P(Fail Meet Sub. Needs)', 'CV Catch',
                  'Harv. Rate on Surplus'), ncol=2)
  if(!any(nms[,1]==nm))
    stop('check name')
  nmind <- which(nm == nms[,1])
  
  # What is the ratio of weir:aerial
  srat <- apply(metaType[[1]], 1, function(x)
    sum(x>1) / sum(x>=1))
  # list the number of stocks monitored
  nstock <- apply(metaType[[1]], 1, function(x)
    sum(x!=0))
  
  # find get the max or min or whatever you want
  # out of the aggregated data
  z <- apply(mres[,,nmind], 1, sfun, na.rm=TRUE)
  
  xdat <- split(srat, nstock)
  zdat <- split(z, nstock)
  
  cf <- colorRampPalette(c('firebrick1', 'slateblue'))
  cols <- cf(length(xdat))
  
  plot(NA, type='n', xlim=range(srat), ylim=range(z),
       xlab='', ylab='')
  for(i in 1:length(xdat)){
    lines(zdat[[i]] ~ xdat[[i]], type='o', col=cols[i],
          pch=16)
  }
  mtext(side=1:2, line=3, cex=1.5,
        c('% weired stocks', nms[nmind,2]))
  
  # make the legend
  xlmin <- min(srat)
  xlmax <- max(srat)
  nl <- 60
  
  xl <- seq(xlmin, xlmax, length.out=nl)[1:(nl-1)]
  xr <- seq(xlmin, xlmax, length.out=nl)[2:nl]
  lc <- cf(nl)
  
  yb <- max(z)+(max(z)-min(z))*0.08
  yt <- yb+(max(z)-min(z))*0.05
  
  rect(xleft=xl, ybottom=yb, xright=xr, ytop=yt,
       xpd=TRUE, col=lc, border=NA)
  rect(xleft=xl[1], ybottom=yb, xright=tail(xr,1), ytop=yt,
       xpd=TRUE, lwd=1)
  
  nlab <- 4
  labx <- seq(nstock[1], tail(nstock, 1), length.out=nlab)
  newr <- (labx - min(labx, na.rm=TRUE)) / 
    max(labx - min(labx, na.rm=TRUE), na.rm=TRUE) * 
    (max(srat, na.rm=TRUE) - min(srat, na.rm=TRUE)) + 
    min(srat, na.rm=TRUE)
  text(x=newr, y=yt, labels=labx, xpd=TRUE, pos=3, offset=0.2)
}
