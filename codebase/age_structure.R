sample_posterior_r_mcmc_hyperR <- function(samps,posi,ni,se,sp,gam0){
  ## This function samples from the posterior distribution
  ## across age bins. It is the same as the Python code, but
  ## for whatever reason, it runs faster. Code by Bailey Fosdick.
  fn <- 1-se
  fp <- 1-sp
  
  nu <- 1 
  
  ## Initial values
  ri <- (posi+1)/(ni+2)
  r <- mean(ri) 
  gam <- gam0  
  
  ## Posterior samples
  ri_post <- matrix(NA,nrow=samps,ncol=length(ri))
  r_post <- rep(NA,samps)
  gam_post <- rep(NA,samps) 
  
  # MCMC tuning parameter (larger values decrease variability in proposals)
  delta_r <- 200
  delta_ri <- 100 
  delta_gam <- 10 
  thin <- 50  
  burn_in <- 2*thin
  
  for(s in 1:(samps*thin+burn_in))
  {
    
    ## MH step for gamma
    ## propose gam_prop | gamma ~ gamma(delta_gam,scale=gamma/delta_gam)
    gam_prop <- rgamma(1,delta_gam,scale=gam/delta_gam)
    ar_gam <- sum(dbeta(ri,r*gam_prop,(1-r)*gam_prop,log=TRUE))-
      sum(dbeta(ri,r*gam,(1-r)*gam,log=TRUE)) +
      dgamma(gam_prop,nu,scale=gam0/nu,log=TRUE) -
      dgamma(gam,nu,scale=gam0/nu,log=TRUE) +
      dgamma(gam,delta_gam,scale=gam_prop/delta_gam,log=TRUE) -
      dgamma(gam_prop,delta_gam,scale=gam/delta_gam,log=TRUE)
    if(log(runif(1))<ar_gam){gam <- gam_prop}
    
    
    # MH step to update r
    # propose r_prop | r ~ B(r*delta_r,(1-r)*delta_r)
    r_prop <- rbeta(1,r*delta_r,(1-r)*delta_r)
    ar_r <- sum(dbeta(ri,r_prop*gam,(1-r_prop)*gam,log=TRUE))-
      sum(dbeta(ri,r*gam,(1-r)*gam,log=TRUE)) +
      dbeta(r,r_prop*delta_r,(1-r_prop)*delta_r,log=TRUE)-
      dbeta(r_prop,r*delta_r,(1-r)*delta_r,log=TRUE)
    if(log(runif(1))<ar_r){r <- r_prop}
    
    #MH step to update each ri
    #propose ri_prop | ri ~ B(ri*delta_ri,(1-ri)*delta_ri)
    for(k in 1:length(ri))
    {
      ri_prop <- rbeta(1,ri[k]*delta_ri,(1-ri[k])*delta_ri)
      if(ri_prop==0 || ri_prop==1) 
      {
        # if error, sample 100 new values 
        ri_propMANY <- rbeta(100,ri[k]*delta_ri,(1-ri[k])*delta_ri)
        ri_prop <- ri_propMANY[which(!(ri_propMANY%in%c(0,1)))[1]] #grab first value not 0 or 1
      }
      ar_ri <- dbinom(posi[k],ni[k],ri_prop*(1-fn)+(1-ri_prop)*fp,log=TRUE)-
        dbinom(posi[k],ni[k],ri[k]*(1-fn)+(1-ri[k])*fp,log=TRUE)+
        dbeta(ri_prop,r*gam,(1-r)*gam,log=TRUE)-
        dbeta(ri[k],r*gam,(1-r)*gam,log=TRUE)+
        dbeta(ri[k],ri_prop*delta_ri,(1-ri_prop)*delta_ri,log=TRUE)-
        dbeta(ri_prop,ri[k]*delta_ri,(1-ri[k])*delta_ri,log=TRUE)
      if(log(runif(1))<ar_ri){ri[k] <- ri_prop}
    }
    
    if(s%%thin==0 && s>burn_in) # problematic if burn_in is not multiple of thin
    {
      gam_post[(s-burn_in)/thin] <- gam 
      r_post[(s-burn_in)/thin] <- r
      ri_post[(s-burn_in)/thin,] <- ri
    }
  }
  param_samps <- cbind(gam_post,r_post,ri_post)
  colnames(param_samps) <- c("gam","r",paste("r",c(1:length(ri)),sep=""))
  return(param_samps)
}
