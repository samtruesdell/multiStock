#------------------------------------------------------------------------------#
# Multi-stock salmon simulation model with subsistence harvest control rule and 
# ability to track management outcomes. 
# Last updated October 14, 2016
#------------------------------------------------------------------------------#

################################################################################
## STEP 1. create functions
################################################################################

#------------------------------------------------------------------------------#
# Subsistence harvest control rule function
#------------------------------------------------------------------------------#
#sub <- subsistence requirement
#egfloor <- escapement goal
#run <- run size (i.e., pre-harvest recruitment) 
#com <- maximum commercial harvest

require(mvtnorm)
require(Matrix)

sub_hcr = function(sub, com, egfloor, run,for.error,OU, dg=FALSE){
	run.est <- run * (1 + rnorm(1, 0, for.error))
	  if(is.na(run.est)==TRUE){run.est <- run}
	if(run.est - egfloor <= 0){hr = 0}
	if(run.est > egfloor){
		if((run.est - egfloor) <= (sub)){hr = (run - egfloor)/run}
		if((run.est - egfloor) > (sub)){
			if((run.est - egfloor) > (sub + com)){hr = (sub + com)/run}
			if((run.est - egfloor) <= (sub + com)){hr = (sub + (run - egfloor - sub))/run}
			}
		}
	hr <- hr + rnorm(1,0,OU)
	if(hr < 0){hr=0}
	if(hr >1 ){hr=1}
	if(dg){
		if(hr > 0.9){
			hr <- 0.9
		}
	}
	hr
}

#------------------------------------------------------------------------------#
# Status function
# estimates equilibrium spwn abundance and catch plus whether over fished 
# or extinct
#------------------------------------------------------------------------------#
#U <- harvest rate
#a <- productivity
#b <- density dependence 
SC.eq <- function(U,a,b){
  a <- log(a)
  S.eq <- max(0,(a-(-log(1-U)))/b)
  C.eq <- max(0,((a-(-log(1-U)))/b)*exp(a-b*((a-(-log(1-U)))/b))-((a-(-log(1-U)))/b))
  OF <- ifelse(U>0.5*a-0.07*a^2,1,0)
  EX <- ifelse(S.eq==0,1,0)
  return(c(S.eq,C.eq,OF,EX))
  }

#------------------------------------------------------------------------------#
# Multi-stock simulation function
#------------------------------------------------------------------------------#
#ny <- the number of years
#Ro <- the sub-stock recruiment at time zero
#rho <- the expected correlation among stocks
#phi <- the expected correlation through time
#Preturn <- the expected proportion of spawners from age 4-7
#episd <- random recruitment variation
#alpha <- sub-stock productivity (not in log space)
#bet <- sub-stock density depedence 
#U <- finite annual exploitation rate
#sub <- subsistence requirement
#com <- maximum commercial harvest
#egfloor <- escapement goal
#pm.yr <- year of simulation that pms start to be calculated over
#for.erro <- forcast error (CV)
#OU <- outcome uncertainty (CV)
#TV <- whether or not there is a change in productivity over cource of simulations
process = function(ny,Ro,rho,phi,Preturn,episd,U,alpha,bet,sub,
                   com,egfloor,pm.yr,for.error,OU,TV,
                   sdlog_h, # lognormal sd for run size observations
                   sdlog_e, # lognormal sd for harvest observations
                   ess_paa_e, # effective sample size for run observations
                   ess_paa_h=NA,
                   nstockM,
                   mTypeV,  # type of stock monitoring.
                   cvAW,
                   # rhoAA,
                   # rhoAW,
                   # rhoWW,
                   cvH,
                   ESS,
                   abias,
                   # sirat_sd,
									 dg=FALSE # indicates if model in data-generating phase. If so don't
									          # let harvest rates get too high.
                   ){
  
  # Edit the escapement goal.  If a single number is given
  # use that throughout the time series.  If two numbers are
  # given create a linear scale (high to low) over the number
  # of years the simulation occurs.
  if(length(egfloor) > 2){
    stop('length of egfloor should be 1 or 2')
  }else if(length(egfloor) == 1){
    egfloorV <- rep(egfloor, 50+ny)
  }else{
    if(egfloor[1] > egfloor[2]){
      stop('egfloor[1] should be < egfloor[2]')
    }
    egfloorV <- c(rep(mean(egfloor), 50),        # vary over the ouput
                  sample(seq(from = egfloor[2],  # period only
                             to = egfloor[1], length.out = ny)))
  }## GOT RID OF SAMPLE FOR TEST ... high bias from 2:1; 
  

  # if(tstockM == 'lgs'){
  #   mstock <- order(alpha / (bet * exp(1)), 
  #                   decreasing=TRUE)[1:nstockM]
  # }else if(tstockM == 'sms'){
  #   mstock <- order(alpha / (bet * exp(1)), 
  #                   decreasing=FALSE)[1:nstockM]
  # }else if(tstockM == 'rnd'){
  #   mstock <- sample(1:length(alpha), size=nstockM)
  # }else stop('tstockM should be lgs, sms or rnd')
  
  
  ny = ny+50
  ns = length(Ro) #number of sub-stocks
  for.error = for.error
  OU = OU
  
  
  # establish CVs and ESS associated with each of the stocks. OK to have
  # 0 cv for non-monitored stocks because these will be excluded from
  # the summed indices.
  cvTab <- matrix(c(0:2, c(0, cvAW)), nrow=3, 
                  dimnames=list(1:3, c('mType', 'CV')))
  stockCV <- cvTab[match(mTypeV, cvTab[,1]),2]
  stockESS <- rep(0, ns)
  stockESS[mTypeV==2] <- ESS
  
  # what are the monitored stocks (both aerial and weired). Assume
  # this applies to the harvest as well as aerial/weir.
  mstock <- which(stockCV > 0)
  
  #Sets the Ricker curve parameter values
  #bet = log(alpha)/Ro
  #Creates correlation among stocks
  R <- matrix(rho,ns,ns) #correlation matrix
  for(i in 1:ns) R[i,i] <- 1
  V <- diag(rep(episd,ns)) # diagonal matrix of sqrt(variances)
  Sig <- V %*% R %*% V    #cov matrix
  epi = rmvnorm(ny, sigma=Sig)
  #Builds time series of Spawners (S), abundance of returning spawners pre-harvest
  # (N), and the component of the residual that is correlated throught time (v)
  R = t(matrix(Ro,ns,ny)) 
  S = R * (1-0)
  v = R; v[,]=0
  R[1:7,]=t(replicate(7,Ro,simplify=T))*exp(epi[1:7,]) 
  N = array(0,dim=c(ny,4,ns))
  Ntot = R; Ntot[,]=0
  H = Ntot
  S = Ntot
  predR = Ntot
  N[5:7,1,]=R[5:7-(4),]*Preturn[1]
  N[6:7,2,]=R[6:7-(5),]*Preturn[2]
  N[7,3,]=R[7-(6),]*Preturn[3]
  
  
  
  for(i in (7+1):ny){   
  	N[i,1,]=R[i-(4),]*Preturn[1]
  	N[i,2,]=R[i-(5),]*Preturn[2]
  	N[i,3,]=R[i-(6),]*Preturn[3]
  	N[i,4,]=R[i-(7),]*Preturn[4]
                 
  	Ntot[i,]=colSums(N[i,,])
  	# scale up to full basin via abundance
  	# rat <- sum(Ntot[i,]) / sum(Ntot[i,mstock])
  	# here the ratio is constant so it is never adjusted.
  	
  	# made changes here according to Jones' suggestion from
  	# 9/6 e-mail. sirat_sd is no longer a parameter but the
  	# CV for the ratio is calculated instead based on the
  	# number of stocks that are sampled.
  	# 
  	# ## Check that rat_p is the whole-stock Smsy based on
  	# Ben's material that he sent.
  	rat_p <- sum(log(alpha[mstock]) / beta[mstock]) / 
  	  sum(log(alpha) / beta)
  	         
  	rat_se <- sqrt(rat_p * (1 - rat_p)/length(alpha))/rat_p
  	sirat_error <- rlnorm(1, meanlog=log(rat_p)-rat_se^2/2, 
  	                      sdlog=rat_se)

  	egfloorV[i] <- egfloorV[i] * sirat_error
  	# get the harvest rate
  	HR= sub_hcr(sub,com,egfloorV[i],sum(Ntot[i,]),for.error,OU, dg=dg)
  	H[i,] = HR*Ntot[i,] # Harvest number
  	S[i,]= Ntot[i,]-H[i,] # escapement
  	
  	if(TV == "TRUE"){
  	  if(i>50){
  	    alpha <- alpha[]*c(0.95,0.95,0.95,0.95,0.95,0.95,
  	                       0.95,0.95,0.95,0.95,0.95,0.95,0.95)
  	  }
  	 }
            
  	R[i,] = alpha[]*S[i,]*exp(-bet[]*S[i,]+phi*v[i-1,]+epi[i,])
  	predR[i,] = alpha[]*S[i,]*exp(-bet[]*S[i,])
  	v[i,] = log(R[i,])-log(predR[i,])
  }

  #output
  S[S[,]=='NaN']=0
  Ntot[Ntot[,]=='NaN']=0
  #performance measures
  pms <- matrix(NA,1,8) # escapement, harvest, harvest rate, 
                        # overfished, extinct, prop years failed 
                        # to meet subsistence harvest
  over<- matrix(NA,length(alpha))
  ext<- matrix(NA,length(alpha))
  harvest_rate = (H[51:ny,]/Ntot[51:ny,])[,1]
  
  # changed "100" to ny-50 in all pm cases (100 used in BC code)
  for(j in 1:length(alpha)){
  	over[j] <- SC.eq(mean(harvest_rate[pm.yr:(ny-50)]),alpha[j],beta[j])[3]
  	ext[j] <- SC.eq(mean(harvest_rate[pm.yr:(ny-50)]),alpha[j],beta[j])[4]
 	}
  pms[,1] <- sum(S[pm.yr:(ny-50),])/((ny-50) - pm.yr +1)
  pms[,2] <- sum(H[pm.yr:(ny-50),])/((ny-50) - pm.yr +1)
  pms[,3] <- mean(harvest_rate[pm.yr:(ny-50)])
  pms[,4] <- sum(over)/length(alpha)
  pms[,5] <- sum(ext)/length(alpha)
  pms[,6] <- sum(rowSums(H[pm.yr:(ny-50),]) < 
                   (sub*0.95))/((ny-50) - pm.yr +1)
  
  HR_surp0 <- (H / (Ntot - S))
  HR_surp0[is.nan(HR_surp0)] <- 0
  HR_surplus <- HR_surp0[pm.yr:(ny-50),1]

  tempcv <- sd(tail(rowSums(H[pm.yr:(ny-50),]), 20)) / 
            mean(tail(rowSums(H[pm.yr:(ny-50),]), 20))
  pms[,7] <- ifelse(any(is.nan(tempcv)), 0, tempcv)
  pms[,8] <- mean(HR_surplus) # Harvest rate on surplus
  
  # pms[,1] <- sum(S[pm.yr:(ny),])/((ny) - pm.yr +1)
  # pms[,2] <- sum(H[pm.yr:(ny),])/((ny) - pm.yr +1)
  # pms[,3] <- mean(harvest_rate[pm.yr:(ny)])
  # pms[,4] <- sum(over)/length(alpha)
  # pms[,5] <- sum(ext)/length(alpha)
  # pms[,6] <- sum(rowSums(H[pm.yr:(ny),]) < 
  #                  (sub*0.95))/((ny) - pm.yr +1)
  
  colnames(pms) <- c('esc', 'harv', 'harvR', 'over', 
                     'ext', 'pFailSub', 'cvCatch', 'HRsurplus')
  
  # output the proportions-at-age N and harvest basin-wide
  naa <- apply(N, c(1,2),sum)
  paaB <- prop.table(naa, margin=1)[51:ny,]
  HB <- apply(H, 1, sum)[51:ny]    # basin harvest
  NB <- apply(Ntot, 1, sum)[51:ny] # basin total N
  SB <- apply(S, 1, sum)[51:ny]    # basin escapement

  # output harvest / run / paa data with error
  # apply errors to each of the harvest series
  HB_emat <- mapply(FUN=function(x,y) oe_t(H[51:ny,x], sdlog=y), 
                    x=1:ncol(H), y=cvH)
  HB_e <- apply(as.matrix(HB_emat[,mstock]), 1, sum)
  
  
  # determine which stocks need to have added bias
  asind <- which(stockCV == cvAW[1])
  biasv <- rep(1, ns)
  biasv[asind] <- 1 + abias
  
  
  # apply errors to each of the escapement estimates
  SB_emat <- mapply(FUN=function(x,y) oe_t(S[51:ny,x] * biasv[x], 
                                           sdlog=y),
                    x=1:ncol(S), y=stockCV)
  
  # cvs <- mTypeV
  # cvs[mTypeV == 1] <- cvAW[1]
  # cvs[mTypeV == 2] <- cvAW[2]
  # 
  # # build correlation matrix
  # cmat <- matrix(NA, ns, ns)
  # if(sum(mTypeV==1) > 1){
  #   icaa <- t(combn(which(mTypeV==1), 2))
  # }else{
  #   icaa <- NULL
  # }
  # if(sum(mTypeV==2) > 1){
  #   icww <- t(combn(which(mTypeV==2), 2))
  # }else{
  #   icww <- NULL
  # }
  # icaw <- as.matrix(expand.grid(which(mTypeV==1), which(mTypeV==2)))
  # if(dim(icaw)[1] != 0){
  #   icaw <- t(sapply(1:nrow(icaw), function(x) sort(icaw[x,])))
  # }
  # 
  # cmat[icaa] <- rhoAA
  # cmat[icww] <- rhoWW
  # cmat[icaw] <- rhoAW
  # diag(cmat) <- 1
  # cmat[is.na(cmat)] <- 0
  # cmat[lower.tri(cmat)] <- t(cmat)[lower.tri(cmat)]
  # 
  # 
  # # get the normal correlated errors
  # # check whether matrix is pd
  # SB_emat <- try(t(sapply(51:ny, function(x){
  #                covmat <- get_covmat(mu=S[x,], cvs=cvs, cormat=cmat)
  #                evs <- eigen(covmat)$values
  #                if(all(evs>0)) browser()
  #                if(is.complex(evs) || any(evs < 0)){
  #                  covmat <- as.matrix(nearPD(covmat)$mat)
  #                  # print(x)
  #                  # browser()
  #                }
  #                 rmvnorm(n = 1, mean=S[x,], sigma=covmat)
  #            })))
  # if(class(SB_emat)=='try-error') stop('mvnorm error')
  # # covmat <- get_covmat(mu=c(0,0), cvs=cvs, cormat=cmat)
  # # mve <- rmvnorm(n=ny-50, mean=0, sigma=cmat)
  # 

  SB_e <- apply(as.matrix(SB_emat[,mstock]), 1, sum)
  
  
  
  # apply errors for the age comps. For now don't separate between
  # escapement age comps and harvest age comps ... assume that
  # any corrections have been made and use just N.

  # get a list of paa according to ESS in each year for each stock
  paaB_earr <- mapply(FUN=function(x,y) oe_p(omat=N[51:ny,,x], ess=y),
                      x=1:(dim(N)[3]), y=stockESS, SIMPLIFY=FALSE)
  
  paaB_weights <- lapply(1:ns, function(x) SB_emat[,x] * paaB_earr[[x]])
  paaB_weightedSum <- Reduce('+', paaB_weights[stockESS>0])
  paaB_e <- prop.table(paaB_weightedSum, margin=1)

  
  
  list(S=S[51:ny,],N=Ntot[51:ny,], H=H[51:ny,], PMs=pms,
       paaB=paaB, HB=HB, NB=NB,
       paaB_e=paaB_e, HB_e=HB_e, SB_e=SB_e, mstock=mstock,
       HR_surplus=HR_surplus)
}



### Functions to add error to the output

# observation error for total numbers and catch
oe_t <- function(ovec, sdlog){
  out <- sapply(1:length(ovec),
  function(i) rlnorm(1, meanlog = log(ovec[i]) - sdlog^2/2,
                     sdlog = sdlog))
  # out <- sapply(1:length(ovec),
  #               function(i) rtnorm(1, mean = ovec[i],
  #                                  sd = ovec[i]*sdlog, #not really log
  #                                  lower=-1e-10))
                           
  return(out)
}



# observation error for paa
oe_p <- function(omat, ess){
  out <- sapply(1:nrow(omat),
                function(i) rmultinom(n=1, size=ess, 
                                      prob=omat[i,]))
  return(prop.table(t(out), margin=1))
}


# function to get the covariance matrix
get_covmat <- function(mu, cvs, cormat){
  ns <- length(mu)
  sds <- mu * cvs
  psdmat <- matrix(NA, ns, ns)
  psdmat[lower.tri(psdmat)] <- combn(sds, 2, FUN=prod)
  psdmat[upper.tri(psdmat)] <- t(psdmat)[upper.tri(psdmat)]
  diag(psdmat) <- sds^2
  # build covariance matrix
  covmat <- psdmat * cormat
  return(covmat)
}


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




# Plotting function for the linear version of the
# Ricker model
plRick <- function(x){
  sr <- x$SRDAT
  a <- x$a
  b <- x$b
  smsy <- x$SMSY
  rmsy <- smsy * exp(a * (1 - smsy/b))
  newx <- seq(min(sr[,1]), max(sr[,1]), length.out=100)
  newy <- newx * exp(a * (1 - newx / b))
  plot(sr[,2] ~ sr[,1], ylim=range(newy))
  lines(newy ~ newx)
  segments(x0=smsy, y0=-0.01*max(sr[,2]), y1=rmsy, lty=2)
  segments(x0=-0.01*max(sr[,1]), x1=smsy, y0=rmsy, lty=2)
}
# plRick(x=SR_stats[[5]])



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



# Staton function to help in the calculation of Smsy reference point
# for multiple populations.
# function to get BRPs
eq_ricker = function(alpha, beta, U) {
  # fished equilibrium states for each substock
  Req_u = log(alpha * (1 - U)) / (beta * (1 - U))
  Req_u = ifelse(Req_u < 0, 0, Req_u)
  Seq_u = Req_u * (1 - U)
  Seq_u = ifelse(Seq_u < 0, 0, Seq_u)
  Ceq_u = Req_u * U
  Ceq_u = ifelse(Ceq_u < 0, 0, Ceq_u)
  
  # sum across substocks
  Seq_tot_u = sum(Seq_u)
  Ceq_tot_u = sum(Ceq_u)
  
  # unfished equilibrium recruitment
  Req = (log(alpha) / beta)
  
  # approximation to Smsy
  Smsy_sum = sum(Req * (0.5 - 0.07 * log(alpha)))
  
  # output
  list(S = Seq_tot_u, C = Ceq_tot_u, Smsy_sum = Smsy_sum)
}

