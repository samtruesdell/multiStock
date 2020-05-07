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