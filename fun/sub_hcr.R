#------------------------------------------------------------------------------#
# Subsistence harvest control rule function
#------------------------------------------------------------------------------#
#sub <- subsistence requirement
#egfloor <- escapement goal
#run <- run size (i.e., pre-harvest recruitment) 
#com <- maximum commercial harvest



sub_hcr = function(sub, com, egfloor, run,for.error,OU, dg=FALSE){
  run.est <- run * (1 + rnorm(1, 0, for.error))
  if(is.na(run.est)==TRUE){run.est <- run}
  if(run.est - egfloor <= 0){hr = 0}
  if(run.est > egfloor){
    if((run.est - egfloor) <= (sub)){hr = (run - egfloor)/run}
    if((run.est - egfloor) > (sub)){
      if((run.est - egfloor) > (sub + com)){hr = (sub + com)/run}
      if((run.est - egfloor) <= (sub + com)){hr = (sub + (run - egfloor - sub))/run}
    }
  }
  hr <- hr + rnorm(1,0,OU)
  if(hr < 0){hr=0}
  if(hr >1 ){hr=1}
  if(dg){
    if(hr > 0.9){
      hr <- 0.9
    }
  }
  hr
}
