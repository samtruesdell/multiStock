


unix <- ifelse(Sys.info()[1] == 'Windows', FALSE, TRUE)

if(unix){
  wdir <- system('pwd', intern=TRUE)
}else
  wdir <- dirname(parent.frame(2)$ofile)




setwd(file.path(wdir, 'runs'))

dirf <- dir()

# get out the results
res <- list()
SR_stats_lst2 <- list()
mType_lst <- list()
metaType <- list()
metaNS <- list()
for(i in seq_along(dirf)){
  e1 <- new.env() 
  load(dirf[i], e1)
  res[[i]] <- e1$pmU
  metaType[[i]] <- e1$mType
  metaNS[[i]] <- e1$ns
  SR_stats_lst2[[i]] <- e1$SR_stats
  cn <- colnames(e1$SR_stats)
  alpha <- e1$alpha
  beta <- e1$beta
  mType_lst[[i]] <- e1$mType
  egscalar <- e1$egscalar # make this available for plotting later
  rm('e1')
}

SR_stats <- do.call(rbind, SR_stats_lst2)
colnames(SR_stats) <- cn

# each result is an "average" of 1.0 so just need to average
# over all the elements in the list to get an average over
# nrep

mres <- Reduce("+", res) / length(res)

resmx <- apply(mres, c(1,3), max) #SOMETIMES NEEDS TO BE MIN!

save.image(paste0(wdir, '/out.Rdata'))

