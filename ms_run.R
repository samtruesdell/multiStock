


unix <- ifelse(Sys.info()[1] == 'Windows', FALSE, TRUE)

if(unix){
  wdir <- system('pwd', intern=TRUE)
}else{
  wdir <- dirname(parent.frame(2)$ofile)
}

setwd(wdir)

# load in seed to use
rseed <- scan('ms_seed.txt')

# load the functions
funLst <- list.files('fun', full.names = TRUE)
sapply(funLst, source)


# number of stocks to be monitored
ns <- 1:13

# sample bias -- any bias for sampling high/low productivity streams?
# options are "high" "low or "none"
sbias <- 'none'

# CVs for aerial surveys [1] and weirs [2] and ESS for weirs
# and correlations btw aerial surveys and weirs
cvAW <- c(1, 0.1)
# rhoAA <- 0.1
# rhoAW <- 0.75
# rhoWW <- 0.01
ESS <- 50
abias <- 0 # percent bias for aerial surveys

# CV for the harvest of every stock
cvH <- 0.2

# escapement goal scalars
egscalar <- seq(from=0.01, to=2, length.out=5)

# number of repetitions to average the results over for each combination
# of mType and EG
nrepProcess <- 3

# number of simulation years (greater than 50)
ny <- 150

# number of years to use to fit SR model
sryrs <- 30

# uncertainty in the scalar of Smsy
# sirat_sd <- 0

# Do you want to include sd on the plots?
plotvar <- FALSE

# Initial value for commercial removals (affects the characteristics
# of the stock-recruitment estimates)
comInit <- 1e6

# alpha = c(5.5,1.9,6.4,2.6,5.8,5.5,2.4,2.1,8.0,3.7,1.7,4.3,2.3)
# beta = c(0.0002986,0.0001986, 0.0003486, 0.001986,
#          0.000168,0.00122,0.000131, 0.00286, 0.0001986, 0.0001986,
#          0.002095,0.000784,0.000911)

# method to randomize alpha and beta pairs based on Brendan's original
# values.
abounds <- matrix(c(5, 8,
                    3, 5,
                    1.5,3),
                  nrow = 3, byrow = TRUE)

bbounds <- matrix(c(0.001, 0.003,
                    0.00055, 0.001,
                    0.0001, 0.00055),
                  nrow = 3, byrow = TRUE)

# choose 13 alpha and betas
aval <- unlist(mapply(function(x,n) runif(n, abounds[x,1], abounds[x,2]),
                      1:nrow(abounds), c(4,5,4)))

bval <- unlist(mapply(function(x,n) runif(n, bbounds[x,1], bbounds[x,2]),
                      1:nrow(bbounds), c(4,5,4)))
# aval <- rep(aval[1], 13);bval <- rep(bval[1], 13)

# randomize the pairing of alphas and betas.
if(sbias == 'low' | sbias =='high'){
  set.seed(rseed)
}
alpha <- sample(aval)
beta <- sample(bval)
set.seed(NULL)
# Calculate the approximate Smsy for use in varying the initial conditions
  U_range <- seq(0, 1, 0.01)
  Seq = sapply(U_range, function(x) eq_ricker(alpha, beta, x)$S)
  Ceq = sapply(U_range, function(x) eq_ricker(alpha, beta, x)$C)
Smsy_approx <- Seq[which.max(Ceq)]


Ro = 2*log(alpha)/beta
rho = 0.4 # among-stock correlation
phi = 0.65 #0.85 # annual correlation
Preturn = c(.2, .39, .38, .03)
episd = 0.6
U = 0.1 #doesn't matteR???
sub = 0 #42500
com = 1e6#35000
# escapement for initial data generating only (it is updated
# after assessment model is fit in each run)
egfloor =  c(0.25*Smsy_approx, 3*Smsy_approx) #30000 #5000
pm.yr = 80
for.error = 0
OU = 0
TV = "FALSE"



require(mvtnorm)
require(Matrix)


# Determine which stocks will be monitored in this run
mType <- stock2mon(nstock = max(ns))


if(pm.yr > (ny-50)){
  stop('pm.yr > ny-50')
}

# baseline random seed for this run
roundrand <- round(runif(1, 1, 1e6))

# object to save the CVs to review.
cvSummary <- numeric(0)

# array to hold the results
# performance metric array dimensions: (1) scenario sampling distribution
# type; (2) escapement goal scalar value; (3) performance metric.
# Performance metrics are: esc, harv, harvR, over, 
#                          ext ,pFailSub, cvCatch, HRsurplus
pmU <- array(data=NA, dim=c(nrow(mType), length(egscalar), 8+1))
pmV <- array(data=NA, dim=c(nrow(mType), length(egscalar), 8+1))
SR_stats_lst <- list()
#get rid of loop over monitoring scenarios  (MJ)
#for(u in 1:nrow(mType)){
u <- 5  #pick an arbitrary scenario (ML)

  
  ## Generate some data. This will be used to fit a stock
  ## recruitment function to. There will be some error
  ## in the output for harvest.
  
  
  inputs <- list(
    alpha = alpha,
    bet = beta,
    Ro = log(alpha)/beta,
    ny = ny,
    rho = rho,
    phi = phi,
    Preturn = Preturn,
    episd = episd,
    U = U,
    sub = sub,
    com = com,
    egfloor =  egfloor,#5000,
    pm.yr = pm.yr,
    for.error = for.error,
    OU = OU,
    TV= TV,
    mTypeV = mType[u,],
    cvAW = cvAW,
    # rhoAA = rhoAA,
    # rhoAW = rhoAW,
    # rhoWW = rhoWW,
    cvH = cvH,
    ESS = ESS,
    abias = abias,
    # sirat_sd=sirat_sd,
		dg = TRUE
  )
  
  pmEG <- matrix(NA, length(egscalar), 8+1)
  pmsd <- matrix(NA, length(egscalar), 8+1)
  
  SR_stats_i <- list()
  for(d in seq_along(egscalar)){
    # set.seed(d + roundrand)   ## can't remember rational for this
    pmsave <- matrix(NA, nrow=nrepProcess, ncol=8+1) #+1 for EG
    
    i <- 1
    while(i <= nrepProcess){
      inputs$com <- comInit
      srdat <- try(do.call(process, inputs))
      if(class(srdat) == 'try-error'){
        break
      }
      
      #### fit an assessment model ####
      
      # spawning stock: just SB_e
      # related recruits have to match ages/years
      
      # Get the (estimated) recruits in every year
      R <- c()
      for(y in 1:(ny-7)){
        R[y] <- diag(srdat$paaB_e[(y+4):(y+7),]) %*% # estimated PAA
          (srdat$SB_e[(y+4):(y+7)] +           # sum of escapement and
             srdat$HB_e[(y+4):(y+7)])          # catch is total recruits
      }
      # print(length(R))
      # print(length(srdat$SB_e))
      # stop()
     
      
      
      # don't need to match S&R indices because just using
      # the end of both series
      
      R <- tail(R, sryrs)
      S <- tail(srdat$SB_e[1:(ny-7)], sryrs)
   

      
      # Calculate the parameters of the Ricker model
      lnRS <- log(R+1e-5) - log(S+1e-5)
      SRlm <- try(lm(lnRS ~ S))
      if(class(SRlm)=='try-error') browser()
      lRpar <- coef(SRlm)
      abase <- lRpar[1]
      bbase <- -lRpar[1] / lRpar[2]
      # get unbiased estimates (H&W p. 269)
      sdR <- summary(SRlm)$sigma
      aprime <- abase + sdR^2/2
      bprime <- aprime / abase * bbase
      Rpar <- c(aprime, bprime)
      
      # rp[[k]] <- Rpar; k <- k+1
      # Calculate Smsy
      Smsy <- Rpar[2] * (0.5 - 0.07 * Rpar[1])

      if(Smsy > 0){
        
        inputs2 <- inputs
        inputs2$egfloor <- Smsy * egscalar[d]
        inputs2$com <- com
				inputs2$dg <- FALSE # no longer data generating phase
    

        srdat2 <- try(do.call(process, inputs2))
        if(class(srdat2) == 'try-error') break
        
        # add the actual escapement to the performance measures.
        pmsave[i,] <- c(srdat2$PMs, inputs2$egfloor)
        if(d == 1){ # only bother with this in round 1 ... EG scalar
                    # does not play a role here.
          mstock <- which(mType[u,] > 0)
          
          # Determine the ratio of productivity for the sampled and
          # unsampled stocks
          
          # CV = K / sqrt(a)
          # Fix the maximum CV at 1.0 and define K from there. Then calculate
          # the CV for each of the combinations of sampled stocks.
          maxCV <- 0.5
          prodSamp <- sum(log(alpha[mstock]) / beta[mstock])
          K <- sqrt(min(log(alpha) / beta)) * maxCV
          cv <- K / sqrt(prodSamp)
          cvSummary <- c(cvSummary, cv)
          require(msm)
          prodSampEst <- rtnorm(1, mean=prodSamp, sd=cv*prodSamp, lower=0)
          
          rat <- sum(log(alpha) / beta) / prodSampEst

          SR_stats_i[[i]] <- c(u, i, # ~!!!!!unc[u] need name??
                               Rpar[1], Rpar[2],
                               Smsy*rat, use.names=FALSE)
        }
        i <- i+1
      }
      if(i > 4*nrepProcess){
        stop('i > 4*nrepProcess')
      }
      if(i %*% 5 == 0) print(paste(i, '/', nrepProcess))
    }
    
    # save results for that scalar
    pmEG[d,] <- apply(pmsave, 2, mean)
    pmsd[d,] <- apply(pmsave, 2, sd)
    
  }
  SR_stats_lst[[u]] <- do.call(rbind, SR_stats_i)
  pmU[u,,] <- pmEG
  pmV[u,,] <- pmsd
  
# }  end of former monitorting scenario loop (MJ)

SR_stats <- do.call(rbind, SR_stats_lst)
colnames(SR_stats) <- c('scen', 'i', 'a', 'b', 'SMSY')

#define y range for plotting harvest and escapement
if(plotvar){
  yrg_eschar <- range(pmU[,,1], pmU[,,2], 
                      pmU[,,1] + pmV[,,1], pmU[,,1] - pmV[,,1],
                      pmU[,,2] + pmV[,,2], pmU[,,2] - pmV[,,2])
  yrg_cvCatch <- range(pmU[,,7] + pmV[,,7],
                       pmU[,,7] - pmV[,,7])
  pmVplot <- pmV
}else{
  yrg_eschar <- range(pmU[,,1], pmU[,,2])
  yrg_cvCatch <- range(pmU[,,7])
  pmVplot <- NULL
}


legtxt <- paste0('s', 1:nrow(mType))
rgx <- c(0, max(pmU[,,9]))



mean(log(alpha)/beta*(0.5-0.07*log(alpha)))


rnn <- round(runif(1, min=1000000, max=9000000))
save.image(paste0(wdir, '/out', 'rnn', '.Rdata'))

dir.create('./fig', showWarnings = FALSE)

if(!unix){
  mres <- pmU
  metaType <- list(mType)
  metaNS <- list(ns)
#  source('ms_plot.r')   
#  source('../get_utility.r')
}




cat(
    'mType', mType,
    'egscalar', egscalar,
    'nrepProcess', nrepProcess,
    'ny', ny,
    'sryrs', sryrs,
    # 'sdlog_e', sdlog_e,
    # 'sdlog_h', sdlog_h,
    # 'ess_paa_e', ess_paa_e,
    # 'ess_paa_h', ess_paa_h,
    'alpha', alpha,
    'beta', beta,
    'Ro', Ro,
    'rho', rho,
    'phi', phi,
    'Preturn', Preturn,
    'episd', episd,
    'U', U,
    'sub', sub,
    'com', com,
    'egfloor', egfloor,
    'pm.yr', pm.yr,
    'for.error', for.error,
    'OU', OU,
    'TV', TV,
    sep='\n',
    file = paste0(wdir, '/param.txt'),
    append = FALSE
)
