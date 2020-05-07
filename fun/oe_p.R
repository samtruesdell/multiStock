# observation error for paa
oe_p <- function(omat, ess){
  out <- sapply(1:nrow(omat),
                function(i) rmultinom(n=1, size=ess, 
                                      prob=omat[i,]))
  return(prop.table(t(out), margin=1))
}
