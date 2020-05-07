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