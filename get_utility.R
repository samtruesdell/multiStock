



# What is the ratio of weir:aerial
srat <- apply(metaType[[1]], 1, function(x)
  sum(x>1) / sum(x>=1))
# list the number of stocks monitored
nstock <- apply(metaType[[1]], 1, function(x)
  sum(x!=0))

# standardize the values for harvest and p(overfishing)
stdf <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}


# extract value at maximum harvest
UHMH <- sapply(1:dim(mres)[1], function(x) mres[x,which.max(mres[x,,2]),2])
UOMH <- sapply(1:dim(mres)[1], function(x) mres[x,which.max(mres[x,,2]),4])
# U2 <- U2-max(U2)

# Utility:harvest, utility:p(overfishing) and overall utility
UH <- stdf(UHMH)
UO <- stdf(UOMH)
U <- UH + (-UO)



# # standardize z on a scale to get the colors
# range01 <- function(x){
#   (x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))
# }
# cr <- colorRamp(c('slateblue', 'firebrick1'))
# stdz <- range01(U)
# 
# if(any(is.nan(stdz) | any(is.na(stdz)))){
#   stdz <- rep(1, length(stdz))
# }
# if(any(is.na(cr(stdz)))) browser()
# col <- rgb(cr(stdz)/255)


# ustock <- unique(nstock[nstock>1])


library(akima)
vv <- interp(srat, nstock, U, duplicate='strip')
filled.contour(vv$x, vv$y, vv$z, 
               plot.title = title(main='Utility', 
                                  xlab = 'Percent weir sampling',
                                  ylab = 'Number of stocks monitored'))







