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