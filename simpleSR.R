
#### Simulation parameters ####

EGInit <- 3000        # Escapement goal avg during burn in period
EGInitSD <- 1000      # Escapement goal sd durning burn in period
EGPar <- c(100, 4000) # Min & max EG during data collection period
EG_typ <- 'random'    # Trend for EG during data collection period
                      # (increase/decrease/random)

ny <- 200      # Total number of years (ny-nySR-7 will be burn in)
nySR <- 30     # Number of years for data collection
Preturn = c(0, 0, 0, .2, .39, .38, .03)  # P(return) by age class
alpha <- 20    # stock alpha
beta <- 0.001  # stock beta
sdR <- 0.2     # sd for realized recruitment

sdH <- 0.25    # sd for harvest
ESS_H <- 50    # ESS for harvest
sdS <- 0.25    # sd for escapement
ESS_S <- 50    # ESS for escapement


#### Packages ####
msm <- require(msm)
if(!msm){
  install.packages('msm')
  require(msm)
}


#### Set up containers ####

N <- matrix(0, nrow=ny, ncol=7)  # N-at-age
S <- numeric(ny)                 # escapement
R <- numeric(ny)                 # recruitment
H <- numeric(ny)                 # harvest

H_obs <- numeric(ny)                    # observed harvest total
HPAA_obs <- matrix(0, nrow=ny, ncol=7)  # observed PAA harvest
HAA_obs <- matrix(0, nrow=ny, ncol=7)   # observed harvest-at-age

S_obs <- numeric(ny)                    # observed escapement total
SPAA_obs <- matrix(0, nrow=ny, ncol=7)  # observed PAA escapement
SAA_obs <- matrix(0, nrow=ny, ncol=7)   # observed escapement-at-age

R_obs <- numeric(ny)                    # observed recruitment total


#### Initial values ####

# Initial values for recruitment to initialize model
R[1:7] <- 4000


#### Set up escapement goal ramp ####
EGdata <- seq(EGPar[1], EGPar[2], length.out=nySR)
if(EG_typ == 'increase'){
  EGdata <- sort(EGdata)
}else if(EG_typ == 'decrease'){
  EGdata <- sort(EGdata, decreasing=TRUE)
}else if(EG_typ == 'random'){
  EGdata <- sample(EGdata)
}else{
  stop('EG_typ not recognized')
}

# escapement goal vector for entire time-period
EG <- c(rtnorm(ny-nySR-7, EGInit, EGInitSD, lower=0),
        EGdata,
        rtnorm(7, EGInit, EGInitSD, lower=0))

# Loop to generate data for assessment
for(i in 8:ny){   

  # Fill the N-at-age matrix
  N[i,4] <- R[i-(4)] * Preturn[4]
  N[i,5] <- R[i-(5)] * Preturn[5]
  N[i,6] <- R[i-(6)] * Preturn[6]
  N[i,7] <- R[i-(7)] * Preturn[7]
  
  # Total run
  Run <- sum(N[i,])
  
  # Annual harvest and escapement
  H[i] <- ifelse(Run < EG[i], 0, Run - EG[i])
  S[i] <- Run - H[i]
  
  # This year's recruitment given escapement
  R[i] = alpha * S[i] * exp(-beta * S[i]) * rlnorm(1, 0, sdR)
  
  # Observations of harvest and escapement
  H_obs[i] <- rtnorm(1, mean=H[i], sd=sdH, lower=0)
  HPAA_obs[i,] <- c(rmultinom(n = 1, size = ESS_H, prob = N[i,])) / ESS_H
  HAA_obs[i,] <- H_obs[i] * HPAA_obs[i,]
  
  S_obs[i] <- rtnorm(1, mean=S[i], sd=sdS, lower=0)
  SPAA_obs[i,] <- c(rmultinom(n = 1, size = ESS_S, prob = N[i,])) / ESS_S
  SAA_obs[i,] <- S_obs[i] * SPAA_obs[i,]
  
}


# Get the (estimated) recruits in every year given harvest
# and escapement
for(i in 1:(ny-7)){
  
  R_obs[i] <- sum( diag( SAA_obs[(i+1):(i+7),] ) ) + 
              sum( diag( HAA_obs[(i+1):(i+7),] ) )
  
}



# Estimation
Estart <- ny - nySR - 7 + 1   # year 1 for model
Eend <- ny-7                  # year 30 for model
R_Edat <- R_obs[Estart:Eend]  # model recruitment data
S_Edat <- S_obs[Estart:Eend]  # model escapement data

# Calculate the parameters of the Ricker model (linear approach)
lnRS <- log(R_Edat) - log(S_Edat)
SRlm <- lm(lnRS ~ S_Edat)
lRpar <- coef(SRlm)
abase <- lRpar[1]
bbase <- -lRpar[1] / lRpar[2]

# get unbiased estimates (H&W p. 269)
lm_sdR <- summary(SRlm)$sigma
aprime <- abase + sdR^2/2
bprime <- aprime / abase * bbase

# Calculate Smsy
Smsy_est <- bprime * (0.5 - 0.07 * aprime)


#### Plotting ####

# Escapement history
plot(EG[Estart:Eend], type = 'o', pch=3, main='Escapement Goal')

# True dynamics
S_true <- seq(0, 5000, length.out = 1000)
R_true <- alpha * S_true * exp(-beta * S_true)
Smsy_true <- log(alpha)/beta * (0.5 - 0.07 * log(alpha))

# x and y limits for plots
xl <- range(S, S_true, S_obs)
yl <- range(R, R_true, R_obs)

# Plot of the true stock / recruit dynamics
plot(R_true ~ S_true, type='l', xlim=xl, ylim=yl, main='True dynamics')
abline(a=0, b=1, lty=3, col='cornflowerblue')
segments(x0=Smsy_true, 
         y0=Smsy_true, 
         y1=R_true[which.min(abs(Smsy_true-S_true))],
         lty=3)


# Actual recruitment
plot(R ~ S, pch=3, col='gray40', xlim=xl, ylim=yl, 
     main = 'Actual recruitment')
lines(R_true ~ S_true)


# Observations
plot(R ~ S, pch=3, col='gray80', xlim=xl, ylim=yl,
     main = 'Observed recruitment (rug is EG)')
lines(R_true ~ S_true)
points(R_Edat ~ S_Edat, pch=3, col='firebrick1')
legend('topright',
       col=c('gray40', 'firebrick1'),
       lty = c(1, NA),
       pch = c(NA, 3),
       bty = 'n',
       cex = 0.75,
       legend = c('True model', 'Observations'))
rug(EG[Estart:Eend], col='cornflowerblue')


# New model
R_mod <- S_true * exp(aprime * (1 - S_true / bprime))
plot(R_Edat ~ S_Edat, pch=3, col='firebrick1', xlim=xl, ylim=yl,
     main = 'Model comparison')
lines(R_true ~ S_true, type='l', col='gray40')
lines(R_mod ~ S_true, col='firebrick1')
abline(a=0, b=1, lty=3, col='cornflowerblue')
segments(x0=Smsy_true, 
         y0=Smsy_true, 
         y1=R_true[which.min(abs(Smsy_true-S_true))],
         lty=3, col='gray40')
segments(x0=Smsy_est, 
         y0=Smsy_est, 
         y1=R_true[which.min(abs(Smsy_est-S_true))],
         lty=3, col='firebrick1')
legend('topright',
       col=c('gray40', 'firebrick1', 'firebrick1', 'gray40', 'firebrick1'),
       lty = c(1,1,NA,3,3),
       pch = c(NA, NA, 3, NA, NA),
       bty = 'n',
       cex = 0.75,
       legend = c('True model', 'Estimate', 'Observations',
                  'True Smsy', 'Est Smsy'))
rug(EG[Estart:Eend], col='cornflowerblue')



