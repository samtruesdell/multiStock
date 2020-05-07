

unix <- ifelse(Sys.info()[1] == 'Windows', FALSE, TRUE)

if(unix){
  wdir <- system('pwd', intern=TRUE)
}else
  wdir <- dirname(parent.frame(2)$ofile)
	



if(unix){
	source('ms_fun.R')
  load('out.Rdata')
}

mType <- metaType[[1]]

setwd(file.path(wdir, 'fig'))

yrg_eschar <- range(mres[,,1], mres[,,2], na.rm=TRUE)
yrg_cvCatch <- range(mres[,,7], na.rm=TRUE)
pmVplot <- NULL
legtxt <- paste0('s', 1:nrow(metaType[[1]]))
rgx <- c(0, max(mres[,,9], na.rm=TRUE))

# plot the PMs

# for(i in 2:length(metaNS[[1]])){
#   
#   zeros <- sapply(1:nrow(metaType[[1]]), function(x){
#     metaType[[1]][x,]!=0
#   })
#   ind <- which(apply(t(zeros), 1, sum) == metaNS[[1]][i])
#   
#   jpeg(paste0('nstock', metaNS[[1]][i], '.jpg'))
# 		par(mfrow=c(3,2), mar=c(0,0.5,1.5,0.5), oma=c(5,5,0,5))
# 		
# 		pfun(x=mres[ind,,,drop=FALSE], v=pmVplot, egX=mres[ind,,9], nm='esc', 
# 				 rgy=yrg_eschar, rgx=rgx, xaxt='n', yaxs=2, txtxp=0, txtyp=0.4)
# 		pfun(x=mres[ind,,,drop=FALSE], v=pmVplot, egX=mres[ind,,9], 
# 		     nm='harv', rgy=yrg_eschar, rgx=rgx, 
# 				 xaxt='n', yaxs=4, txtxp=0, txtyp=0.9)
# 		pfun(x=mres[ind,,,drop=FALSE], v=pmVplot, egX=mres[ind,,9], 
# 		     nm='harvR', rgy=c(0,1), rgx=rgx, 
# 				 xaxt='n', yaxs=2, txtxp=0, txtyp=0.2)
# 		pfun(x=mres[ind,,,drop=FALSE], v=pmVplot, egX=mres[ind,,9], 
# 		     nm='cvCatch', rgy=yrg_cvCatch, rgx=rgx, 
# 				 xaxt='n', yaxs=4, txtxp=0, txtyp=0.2)
# 		pfun(x=mres[ind,,,drop=FALSE], v=pmVplot, egX=mres[ind,,9], 
# 		     nm='ext', rgy=c(0,1), rgx=rgx, 
# 				 xaxt='s', yaxs=2, txtxp=0, txtyp=0.2)
# 		pfun(x=mres[ind,,,drop=FALSE], v=pmVplot, egX=mres[ind,,9], 
# 		     nm='pFailSub', rgy=c(0,1), rgx=rgx, 
# 				 xaxt='s', yaxs=4, leg=TRUE, legtxt=legtxt, legpos='topright',
# 				 txtxp=0, txtyp=0.1)
# 		mtext('Escapement Goal', side=1, line=3, cex=1.25, outer=TRUE)
# 	dev.off()
#   
# }


# number of stocks monitored as vector
ns <- apply(mType>0, 1, sum)

### Plot Smsy estimates
# version where colors represent # of samples > 7
jpeg('Smsy_nsamp.jpg')
  par(oma=c(0,0,0,0), mar=c(5,5,1,1), mfrow=c(1,1))
  
  # determine the proportion of weired stocks
  stockpa <- mType
  stockpa[stockpa>0] <- 1
  stockpa <- apply(stockpa, 1, sum)
  stockpab <- sapply(stockpa, function(x) ifelse(x<max(stockpa)/2,
                                                 'firebrick1',
                                                 'cornflowerblue'))
  
  
  boxplot(SR_stats[,'SMSY'] ~ SR_stats[,'scen'],
  xlab='Scenario', ylab='Smsy estimate', col=stockpab, outline=FALSE)
  
  # Old, incorrect method of calculating a multi-stock Smsy
  # abline(h=sum(log(alpha) / beta * (0.5 - 0.07 * log(alpha))), lty=3)
  
  # Correct method for estimating multistock Smsy
  U_range <- seq(0, 1, 0.01)
  Seq = sapply(U_range, function(x) eq_ricker(alpha, beta, x)$S)
  Ceq = sapply(U_range, function(x) eq_ricker(alpha, beta, x)$C)
  Smsy_multi <- Seq[which.max(Ceq)]
  abline(h = Smsy_multi)
  
  # divisions for the number of stocks
  nstockdiv <- cumsum(rle(ns)$lengths)
  abline(v=nstockdiv-0.5, lty=3)
dev.off()

# version where colors represent % weird samples
jpeg('Smsy_weir.jpg')
par(oma=c(0,0,0,0), mar=c(5,5,1,1), mfrow=c(1,1))

# determine the proportion of weired stocks
srat <- apply(metaType[[1]], 1, function(x)
  sum(x>1) / sum(x>=1))
# assign binary value to > 50% weired stocks
sratb <- sapply(srat, function(x) ifelse(x > 0.7, 'cornflowerblue', 
                                         'firebrick1'))


boxplot(SR_stats[,'SMSY'] ~ SR_stats[,'scen'],
        xlab='Scenario', ylab='Smsy estimate', col=sratb,ylim=c(0, 10000))
abline(h = Smsy_multi)
# divisions for the number of stocks
  abline(v=nstockdiv-0.5, lty=3)
dev.off()



### Plot the SR functions
jpeg('SREstimatePlots.jpg', width=480*1.5)
  escTemp <- seq(1, 70000, length.out=50)
  getval <- sapply(1:nrow(SR_stats), function(x)
    escTemp * exp(SR_stats[x,'a'] * (1-escTemp/SR_stats[x,'b'])))
  cr <- colorRampPalette(c('firebrick1', 'cornflowerblue'))
  cols <- cr(ncol(getval))
  matplot(x=escTemp, getval, type='l', col=cols, lty=1,
          xlab = 'Escapement', ylab = 'Recruits',
          main = 'blue = more stocks / more weirs')
  sapply(1:nrow(SR_stats), function(x)
    rug(SR_stats[x,'SMSY'], col=cols[x]))
dev.off()




### Plot the sampling distribution by stock / scenario
if(unix) mType <- mType_lst[[1]]
jpeg('samples.jpg')
  pm(mType)
dev.off()



toplot <- c('esc', 'harv', 'pFailSub', 'ext', 'cvCatch', 
            'harvR', 'over', 'HRsurplus')
tpfun <- list(max, max, min, min, min, min, min, min)
rdarg <- c(-3, -3, rep(2, 6))

for(i in 1:length(toplot)){
  jpeg(paste0('RNplt_', toplot[i], '.jpg'),
       width=480*3, height=480*2, res=72*3)
    par(mar=c(4,4,5,1), oma=rep(0,4))
    nrPlot(nm=toplot[i], legrd=rdarg[i], sfun=tpfun[[i]])
  dev.off()
  
  # jpeg(paste0('RATplt_', toplot[i], '.jpg'),
  #      width=480*3, height=480*2, res=72*3)
  #   par(mar=c(4,4,5,1), oma=rep(0,4))
  #   ratMetPlot(nm=toplot[i], sfun=tpfun[[i]])
  # dev.off()
}








# Plot example performance comparison for weirs vs non-weirs
jpeg('harvestByScalar.jpg',
     width=480*3, height=480*2, res=72*3)
  par(mar=c(4,4,1,1))
  nr <- nrow(mType)
  yl <- range(mres[1,,2], mres[nr,,2])
  plot(mres[1,,2] ~ egscalar, type='l', lwd=3, ylim=yl, col='cornflowerblue',
       xlab='', ylab='')
  lines(mres[nr,,2] ~ egscalar, lwd=3, col='firebrick1')
  mtext(side=1:2, line=c(2.5,2.5), cex=1.25,
       c('Escapement goal scalar', 'Harvest'))
  legend('bottomright', col=c('cornflowerblue', 'firebrick1'),
         lty=1, lwd=3, bty='n', cex=1.25,
         legend=c('n subpop = 2', 'n subpop = 13'))
dev.off()





# number of stocks monitored as vector
ns <- apply(mType>0, 1, sum)

#3 / 13 stocks
harv3 <- apply(mres[which(ns == 3),,2], 2, mean)
harv13 <- apply(mres[which(ns == 13),,2], 2, mean)

# Plot example performance comparison for avg Nstocks
jpeg('nstockByScalar.jpg',
     width=480*3, height=480*2, res=72*3)
par(mar=c(4,4,1,1))
yl <- range(harv3, harv13)
plot(harv3 ~ egscalar, type='l', lwd=3, ylim=yl, col='cornflowerblue',
     xlab='', ylab='')
lines(harv13 ~ egscalar, lwd=3, col='firebrick1')
mtext(side=1:2, line=c(2.5,2.5), cex=1.25,
      c('Escapement goal scalar', 'Harvest'))
legend('bottomright', col=c('cornflowerblue', 'firebrick1'),
       lty=1, lwd=3, bty='n', cex=1.25,
       legend=c('n subpop = 3', 'n subpop = 13'))
dev.off()




# stocks monitored by weir > 0.75
ws75 <- which(srat >= 0.75)

# stocks monitored by weir < 0.25
ws25 <- which(srat <= 0.25)

# 75/25 stocks
harv75 <- apply(mres[ws75,,2], 2, mean)
harv25 <- apply(mres[ws25,,2], 2, mean)

# Plot example performance comparison for avg Nstocks
jpeg('weirRatioByScalar.jpg',
     width=480*3, height=480*2, res=72*3)
par(mar=c(4,4,1,1))
yl <- range(harv75, harv25)
plot(harv25 ~ egscalar, type='l', lwd=3, ylim=yl, col='cornflowerblue',
     xlab='', ylab='')
lines(harv75 ~ egscalar, lwd=3, col='firebrick1')
mtext(side=1:2, line=c(2.5,2.5), cex=1.25,
      c('Escapement goal scalar', 'Harvest'))
legend('bottomright', col=c('cornflowerblue', 'firebrick1'),
       lty=1, lwd=3, bty='n', cex=1.25,
       legend=c('%weirs < 0.25', '%weirs > 0.75'))
dev.off()

