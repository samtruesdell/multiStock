# plot function identifying min/max/etc of metric relative to
# the ratio of weir surveys and the numb of stocks sampled
nrPlot <- function(nm, legrd=-3, sfun,
                   xlmin=1, xlmax=13){
  
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
  rgx <- c(1, 13)
  rgy <- c(0, 1)
  
  # What is the ratio of weir:aerial
  srat <- apply(metaType[[1]], 1, function(x)
    sum(x>1) / sum(x>=1))
  # list the number of stocks monitored
  nstock <- apply(metaType[[1]], 1, function(x)
    sum(x!=0))
  
  # find get the max or min or whatever you want
  # out of the aggregated data
  # z <- apply(mres[,,nmind], 1, sfun, na.rm=TRUE)
  
  # extract value at maximum harvest
  z <- sapply(1:dim(mres)[1], function(x) mres[x,which.max(mres[x,,2]),nmind])
  z <- z-max(z)
  
  # standardize z on a scale to get the colors
  range01 <- function(x){
    (x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
  }
  cr <- colorRamp(c('slateblue', 'firebrick1'))
  stdz <- range01(z)
  
  if(any(is.nan(stdz) | any(is.na(stdz)))){
    stdz <- rep(1, length(stdz))
  }
  if(any(is.na(cr(stdz)))) browser()
  col <- rgb(cr(stdz)/255)
  
  # get distances between the ratio points to set up the
  # boxes on the plot
  # getmin <- function(x){
  #   yc <- combn(x, 2)
  #   yr0 <- abs(yc[1,]-yc[2,])
  #   yr <- min(yr0[yr0>0], na.rm=TRUE)
  #   yr
  # }
  u <- unique(nstock[nstock>1])
  # mind <- min(sapply(1:length(u), function(x) 
  #   getmin(srat[nstock==u[x]])), na.rm=TRUE)
  # 
  # # set up the plot
  # plot(srat ~ nstock, type='n', las=1, xlab='', ylab='',
  #      xlim=c(xlmin, xlmax), ylim=c(0, 1))
  # rect(xleft=nstock-0.5, xright=nstock+0.5,
  #      ybottom=srat-mind/2, ytop=srat+mind/2,
  #      col=col)
  # mtext(side=1:2, line=2.5, cex=1.5,
  #       c('Number of stocks', 'Proportion weired stocks'))
  # box()
  # browser()
  
  library(akima)
  vv <- interp(srat, nstock, stdz, duplicate='strip')
  filled.contour(vv$x, vv$y, vv$z,
                 plot.title = title(main=nm, 
                                    xlab = 'Percent weir sampling',
                                    ylab = 'Number of stocks monitored'))
  
  
  # # make the legend
  # xlmin <- 1
  # xlmax <- 13
  # nl <- 60
  # 
  # xl <- seq(xlmin, xlmax, length.out=nl)[1:(nl-1)]
  # xr <- seq(xlmin, xlmax, length.out=nl)[2:nl]
  # lc <- rgb(cr(seq(0, 1, length.out=nl))/255)
  # 
  # rect(xleft=xl, ybottom=1.1, xright=xr, ytop=1.3,
  #      xpd=TRUE, col=lc, border=NA)
  # rect(xleft=xl[1], ybottom=1.1, xright=xr[nl-1], ytop=1.3,
  #      xpd=TRUE, lwd=1)
  # 
  # nlab <- 4
  # labx <- seq(xlmin, xlmax, length.out=nlab)
  # newr <- (labx - min(labx, na.rm=TRUE)) / 
  #   max(labx - min(labx, na.rm=TRUE), na.rm=TRUE) * (max(z, na.rm=TRUE) - min(z, na.rm=TRUE)) + min(z, na.rm=TRUE)
  # ltext <- round(newr, legrd)
  # text(x=labx, y=1.35, labels=ltext, xpd=TRUE)
}
