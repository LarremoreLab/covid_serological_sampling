# Inputs:
## samps = number of MCMC samples desired
## pos = number of positive tests pop
## n = number of neg tests in pop
## tp = true positive tests in the lab
## tn = true negative tests in the lab
## fp = false positive tests in the lab
## fn = false negative tests in the lab


# Output:
## (samps x 3) matrix of posterior samples: [r,se,sp]
sample_posterior_r_mcmc_testun <- function(samps,pos,n,tp,tn,fp,fn){
  
  ## Initial values
  sp <- (tn+1)/(tn+fp+2)
  se <- (tp+1)/(tp+fn+2)
  r <- (pos+1)/(n+2)
  
  ## Posterior samples
  r_post <- rep(NA,samps)
  se_post <- rep(NA,samps) 
  sp_post <- rep(NA,samps) 
  
  # MCMC tuning parameter (larger values decrease variability in proposals)
  delta_r <- 100*(1+floor(n/3000))
  delta_sp <- 100*(1+floor((tn+fp)/3000))
  delta_se <- 100*(1+floor((tp+fn)/3000))
  
  # push proposal for r slightly towards .5 to avoid boundary issues
  r_cent <- 1 + (1/n)*(-1)^(r>.5) # CHANGED 
  
  if(pos/n < 1-sp){
    delta_sp <- 100*(1+floor((n+tn+fp)/3000))
  }

  thin <- 50  
  burn_in <- 2*thin
  ac_r <- ac_se <- ac_sp <- 0
  for(s in 1:(samps*thin+burn_in))
  {
    
    #MH step to update r
    #propose r_prop | r ~ B((r*r_cent)*delta_r,(1-(r*r_cent))*delta_r)
      r_prop <- rbeta(1,(r*r_cent)*delta_r,(1-(r*r_cent))*delta_r) ##CHANGED
      ar_r <- dbinom(pos,n,r_prop*se+(1-r_prop)*(1-sp),log=TRUE)-
        dbinom(pos,n,r*se+(1-r)*(1-sp),log=TRUE)+
        dbeta(r,(r_prop*r_cent)*delta_r,(1-(r_prop*r_cent))*delta_r,log=TRUE)- ## CHANGED
        dbeta(r_prop,(r*r_cent)*delta_r,(1-(r*r_cent))*delta_r,log=TRUE) ## CHANGED
      if(log(runif(1))<ar_r){r <- r_prop;ac_r <- ac_r+1}
      
      
      #MH step to update se
      #propose se_prop | se ~ B(se*delta_se,(1-se)*delta_se)
      se_prop <- rbeta(1,se*delta_se,(1-se)*delta_se)
      ar_se <- dbinom(pos,n,r*se_prop+(1-r)*(1-sp),log=TRUE)-
        dbinom(pos,n,r*se+(1-r)*(1-sp),log=TRUE)+
        dbinom(tp,(tp+fn),se_prop,log=TRUE)-
        dbinom(tp,(tp+fn),se,log=TRUE)+
        dbeta(se,se_prop*delta_se,(1-se_prop)*delta_se,log=TRUE)-
        dbeta(se_prop,se*delta_se,(1-se)*delta_se,log=TRUE)
      if(log(runif(1))<ar_se){se <- se_prop;ac_se <- ac_se+1}
      
      
      ## MH step to update sp
      sp_prop <- rbeta(1,sp*delta_sp,(1-sp)*delta_sp)
      ar_sp <- dbinom(pos,n,r*se+(1-r)*(1-sp_prop),log=TRUE)-
        dbinom(pos,n,r*se+(1-r)*(1-sp),log=TRUE)+
        dbinom(tn,(fp+tn),sp_prop,log=TRUE)-
        dbinom(tn,(fp+tn),sp,log=TRUE)+
        dbeta(sp,sp_prop*delta_sp,(1-sp_prop)*delta_sp,log=TRUE)-
        dbeta(sp_prop,sp*delta_sp,(1-sp)*delta_sp,log=TRUE)
      if(log(runif(1))<ar_sp){sp <- sp_prop;ac_sp <- ac_sp+1}
    
    if(s%%thin==0 && s>burn_in) # problematic if burn_in is not multiple of thin
    {
      r_post[(s-burn_in)/thin] <- r
      se_post[(s-burn_in)/thin] <- se
      sp_post[(s-burn_in)/thin] <- sp
    }
  }
  #print(paste("Acceptance rates: ",round(c(ac_r,ac_se,ac_sp)/s,2)))
  param_samps <- cbind(r_post,se_post,sp_post)
  colnames(param_samps) <- c("r","se","sp")
  return(param_samps)
}



# tp <- 279
# fn <- 21
# fp <- 25
# tn <- 975
# 
# 
# n <- 1000
# pos <- 150
# samps <- 1000
# 
# res <- sample_posterior_r_mcmc_testun(samps,pos,n,tp,tn,fp,fn)
# 
# par(mfrow=c(3,1))
# plot(res[,1],type="l")
# plot(res[,2],type="l")
# plot(res[,3],type="l")
# 
# library(coda)
# apply(res,2,effectiveSize)
# 
# par(mfrow=c(3,1))
# plot(density(res[,1]),main="r")
# plot(density(res[,2]),main="se")
# plot(density(res[,3]),main="sp")
# 
# 



# tp <- 25+78+75
# fn <- 12+7
# fp <- 0+2
# tn <- 30+369
# 
# 
# neg <- 3280
# pos <- 50
# n <- neg+pos
# samps <- 1000
# 
# res <- sample_posterior_r_mcmc_testun(samps,pos,n,tp,tn,fp,fn)

# library(coda)
# apply(res,2,effectiveSize)
# 
# par(mfrow=c(3,1))
# plot(density(res[,1]),main="r")
# plot(density(res[,2]),main="se")
# plot(density(res[,3]),main="sp")