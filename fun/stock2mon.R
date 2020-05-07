
stock2mon <- function(nstock){
  wstock <- t(sapply(1:nstock, function(x){
    n <- 1:nstock
    if(sbias == 'low'){
      set.seed(rseed)
      s <- n[order(alpha)][1:x]
      set.seed(NULL)
    }else if(sbias == 'high'){
      set.seed(rseed)
      s <- rev(n[order(alpha)])[1:x]
      set.seed(NULL)
    }else if(sbias == 'none'){
      s <- sample(n, size=x)
    }else{
      stop('error in sample bias assignment')
    }
    g <- rep(0, nstock)
    g[s] <- 1
    return(g)
  }))
  
  out <- list()
  k <- 0
  for(i in 1:length(ns)){
    v <- wstock[i,]
    k <- k+1
    out[[k]] <- v
    if(length(which(v>0)) > 1){
      whichOne <- sample(which(v > 0))
    }else{
      whichOne <- which(v > 0)
    }
    for(j in whichOne){
      v2 <- v
      v2[j] <- 2
      k <- k+1
      out[[k]] <- v2
      v <- v2
    }
  }
  
  mType_all <- do.call(rbind, out)
  
  
  # remove any case where there are zero weir surveys
  wNoWeir <- sapply(1:nrow(mType_all),
                    function(x) any(mType_all[x,] == 2))
  return(mType_all[wNoWeir,])
  
}
