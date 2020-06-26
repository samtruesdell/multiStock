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
    # MLJ: remove these calculations - they shouldn't affect
    #      the data generation phase
#    if (sirat_sd == 'yes') {
#      rat_p <- sum(log(alpha[mstock]) / beta[mstock]) / 
#        sum(log(alpha) / beta)
#      rat_se <- sqrt(rat_p * (1 - rat_p)/length(alpha))/rat_p
#      sirat_error <- rlnorm(1, meanlog=log(rat_p)-rat_se^2/2, 
#                            sdlog=rat_se)
#      }
#    else {sirat_error <-  1  } # turn off uncert due to scaling... 

    
#    egfloorV[i] <- egfloorV[i] * sirat_error
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