#setwd('~/Google drive (jlewnard@berkeley.edu)/graham natural immunity')



library(truncdist)
library(parallel)
library(coxme)
library(survival)
library(truncnorm)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

good_odds_plot <- function(theta,family,params,negmult,numstrap,numcases,obstime){
  #num_cases <- 100
  #numstrap <- 100
  means <- c()
  vars = c()
  OBS_TIME <- obstime
  
  vals = seq(0,3,0.005)
  dbactrian = function(y,x,sigma,m){
    out = (1/(2*sigma*sqrt(2*pi*(1-m^2))))*(exp(-((y-x+m*sigma)^2)/(2*(1-m^2)*sigma^2)) + exp(-((y-x-m*sigma)^2)/(2*(1-m^2)*sigma^2))) 
    out[y<0]= 0
    return(out)
  }
  
  obs_vec <- c()
  
  for (j in 1:numstrap){
    
    # if (j%%5 == 0){
    #   sprintf("Sample no. %d",j)
    # }
    #
    num_sum <- 0
    den_sum <- 0
    
    Esum <- 0
    Fsum <- 0
    Gsum <- 0
    Hsum <- 0
    
    for (i in 1:numcases){
      
      if (family=="Truncated Normal"){
        lambda_n <- 1#1*negmult#rtruncnorm(1,a=0,b=Inf,mean=params[1],sd=params[2])
        lambda_tb <- rtruncnorm(1,a=0,b=Inf,mean=params[1],sd=params[2])
      } else if (family=="Exponential"){
        lambda_n <- 1#1*negmult#rexp(1, rate=params[1])
        lambda_tb <- rexp(1, rate=params[1])
      } else if (family=='Uniform'){
        lambda_n = 1#runif(1,params[1],params[2])
        lambda_tb = runif(1,params[1],params[2])
      } else if (family=='Truncated Cauchy'){
        lambda_n = 1#rtrunc(1,'cauchy',a=0,b=Inf,location=params[1],scale=params[2])
        lambda_tb = rtrunc(1,'cauchy',a=0,b=Inf,location=params[1],scale=params[2])
      } else if (family=='Gamma'){
        lambda_n = 1#rgamma(1,shape=params[1],scale=params[2])
        lambda_tb = rgamma(1,shape=params[1],scale=params[2])
      } else if (family=='Bactrian'){
        lambda_n = 1#sample(vals,1,prob=dbactrian(vals,x=params[1],sigma=params[2],m=params[3]))
        lambda_tb = sample(vals,1,prob=dbactrian(vals,x=params[1],sigma=params[2],m=params[3]))
      } else if (family=="Mixture of Gammas"){
        randnum <- runif(1)
        randnum2 <- runif(1)
        
        if (randnum <= 0.5) {
          lambda_n <- 1#gamma(1,shape=params[1],scale=params[2])
        } else {
          lambda_n <- 1#rgamma(1,shape=params[3],scale=params[4])
        }
        
        
        if (randnum2 <= 0.5) {
          lambda_tb <- rgamma(1,shape=params[1],scale=params[2])
        } else {
          lambda_tb <- rgamma(1,shape=params[3],scale=params[4])
        }
      } else if (family=='Mixture of Uniforms'){
        randnum = runif(1)
        randnum2 = runif(1)
        if (randnum<=0.5){
          lambda_n = 1#runif(1,params[1],params[2])
        } else{
          lambda_n = 1#runif(1,params[1],params[2])
        }
        
        if (randnum2<=0.5){
          lambda_tb = runif(1,params[1],params[2])
        } else{
          lambda_tb = runif(1,params[3],params[4])
        }
      } else{
        stop("Use a valid family!")
      }
      
      A <- rexp(1,lambda_tb)
      B <- rexp(1,lambda_n)
      C <- rexp(1,theta*lambda_tb)
      D <- rexp(1,lambda_n)
      
      # if(A+C<2|B+D<2){
      #   next
      # }
      
      if (A<B){
        Gi <- 1
        Hi <- 0
      } else {
        Gi <- 0
        Hi <- 1
      }
      
      if (C<D){
        Ei <- 1
        Fi <- 0
      } else {
        Ei <- 0
        Fi <- 1
      }
      
      # if ((A+C)<OBS_TIME){
      num_sum <- num_sum + Ei*Hi
      den_sum <- den_sum + Fi*Gi          
      Esum <- Esum + Ei
      Fsum <- Fsum + Fi
      Gsum <- Gsum + Gi
      Hsum <- Hsum + Hi
      # }
      
    }
    
    obs_vec <- c(obs_vec,(num_sum/den_sum))
    #bad_vec <- c(bad_vec,(Esum*Hsum)/(Fsum*Gsum))
    
  }
  means <- c(means,mean(obs_vec,na.rm=T))
  vars = c(vars,var(obs_vec))
  
  return(means)
}



alt_analysis <- function(theta,family,params,negmult,numstrap,numcases,obstime){
  
  #num_cases <- 100
  #numstrap <- 100
  
  meansBase = meansFixed = meansRand = meansFrail <- c()
  
  OBS_TIME <- obstime
  
  obs_vecBase = obs_vecFixed = obs_vecRand = obs_vecFrail <- c()
  
  vals = seq(0,3,0.005)
  dbactrian = function(y,x,sigma,m){
    out = (1/(2*sigma*sqrt(2*pi*(1-m^2))))*(exp(-((y-x+m*sigma)^2)/(2*(1-m^2)*sigma^2)) + exp(-((y-x-m*sigma)^2)/(2*(1-m^2)*sigma^2))) 
    out[y<0]= 0
    return(out)
  }
  
  frails = c()
  
  for (j in 1:numstrap){
    
    master.list <- c()
    
    for (i in 1:numcases){
      if (family=="Truncated Normal"){
        lambda_n <- 1#rtruncnorm(1,a=0,b=Inf,mean=params[1],sd=params[2])
        lambda_tb <- rtruncnorm(1,a=0,b=Inf,mean=params[1],sd=params[2])
      } else if (family=="Exponential"){
        lambda_n <- 1#rexp(1, rate=params[1])
        lambda_tb <- rexp(1, rate=params[1])
      } else if (family=='Uniform'){
        lambda_n = 1#runif(1,params[1],params[2])
        lambda_tb = runif(1,params[1],params[2])
      } else if (family=='Truncated Cauchy'){
        lambda_n = 1#rtrunc(1,'cauchy',a=0,b=Inf,location=params[1],scale=params[2])
        lambda_tb = rtrunc(1,'cauchy',a=0,b=Inf,location=params[1],scale=params[2])
      } else if (family=='Gamma'){
        lambda_n = 1#rgamma(1,shape=params[1],scale=params[2])
        lambda_tb = rgamma(1,shape=params[1],scale=params[2])
      } else if (family=='Bactrian'){
        lambda_n = 1#sample(vals,1,prob=dbactrian(vals,x=params[1],sigma=params[2],m=params[3]))
        lambda_tb = sample(vals,1,prob=dbactrian(vals,x=params[1],sigma=params[2],m=params[3]))
      } else if (family=="Mixture of Gammas"){
        randnum <- runif(1)
        randnum2 <- runif(1)
        
        if (randnum <= 0.5) {
          lambda_n <- 1#rgamma(1,shape=params[1],scale=params[2])
        } else {
          lambda_n <- 1#rgamma(1,shape=params[3],scale=params[4])
        }
        
        if (randnum2 <= 0.5) {
          lambda_tb <- rgamma(1,shape=params[1],scale=params[2])
        } else {
          lambda_tb <- rgamma(1,shape=params[3],scale=params[4])
        }
      } else if (family=='Mixture of Uniforms'){
        randnum = runif(1)
        randnum2 = runif(1)
        if (randnum<=0.5){
          lambda_n = 1#runif(1,params[1],params[2])
        } else{
          lambda_n = 1#runif(1,params[1],params[2])
        }
        
        if (randnum2<=0.5){
          lambda_tb = runif(1,params[1],params[2])
        } else{
          lambda_tb = runif(1,params[3],params[4])
        }
        
      } else{
        stop("Use a valid family!")
      }
      
      num_TB <- 0
      total_time <- 0
      j <- 1
      data.list <- c()
      
      while ((total_time<OBS_TIME)&num_TB<2){
        if (num_TB == 0){
          #print("First Infect")
          time <- rexp(1,lambda_tb)
          if (time + total_time>=OBS_TIME){
            event <- 0
            time <- OBS_TIME - total_time
            total_time <- OBS_TIME
          } else {
            event <- 1
            total_time <- total_time + time
          }
        } else {
          #print("Secondary Infect")
          time <- rexp(1,theta*lambda_tb)
          if (time + total_time>=OBS_TIME){
            event <- 0
            time <- OBS_TIME - total_time
            total_time <- OBS_TIME
          } else {
            event <- 1
            total_time <- total_time + time
          }
        }
        data.vec <- c("ID"=i,"dt"=time,"X"=num_TB>0,"Y"=event)
        data.list[[j]] <- data.vec
        j <- j + 1
        if (event == 1){
          num_TB <- num_TB + 1
        }
      }
      master.list <- c(master.list,data.list)
      
      
    }
    
    resultdfr <- as.data.frame(do.call("rbind", master.list)) 
    resultdfr = resultdfr[resultdfr$X<=1,]
    
    
    mBase = coxph(Surv(dt,Y)~X,data=resultdfr,control=coxph.control(iter.max=5,eps=1e-6))
    mFixed = coxph(Surv(dt,Y)~X+as.factor(ID),data=resultdfr,control=coxph.control(iter.max=5,eps=1e-6))
    mRand = coxme(Surv(dt,Y)~X+(1|ID),data=resultdfr,control=coxme.control(iter.max=5,eps=1e-6))
    mFrail = coxph(Surv(dt,Y)~X+frailty(ID,distribution='gamma'),data=resultdfr,control=coxph.control(iter.max=5,eps=1e-6))
    #frails = c(frails,mFrail$frail)

    estBase <- exp(coef(mBase)['X'])
    estFixed <- exp(coef(mFixed)['X'])
    estRand <- exp(coef(mRand)['X'])
    estFrail <- exp(coef(mFrail)['X'])
    
    obs_vecBase <- c(obs_vecBase,estBase)
    obs_vecFixed <- c(obs_vecFixed,estFixed)
    obs_vecRand <- c(obs_vecRand,estRand)
    obs_vecFrail <- c(obs_vecFrail,estFrail)
    
    #  if (length(mFrail$frail)>0){
    #    frails[k,j,1:length(mFrail$frail)] = exp(mFrail$frail)
    #  }
    
  }
  
  meansBase <- c(meansBase,mean(obs_vecBase,na.rm=T))
  meansFixed <- c(meansFixed,mean(obs_vecFixed,na.rm=T))
  meansRand <- c(meansRand,mean(obs_vecRand,na.rm=T))
  meansFrail <- c(meansFrail,mean(obs_vecFrail,na.rm=T))
  
  return(cbind(meansBase,meansFixed,meansRand,meansFrail))
#  return(frails)
}

get_effects <- function(theta,family,params,negmult,numstrap,numcases,obstime){
  
  #num_cases <- 100
  #numstrap <- 100
  
  meansBase = meansFixed = meansRand = meansFrail <- c()
  
  OBS_TIME <- obstime
  
  obs_vecBase = obs_vecFixed = obs_vecRand = obs_vecFrail <- c()
  
  vals = seq(0,3,0.005)
  dbactrian = function(y,x,sigma,m){
    out = (1/(2*sigma*sqrt(2*pi*(1-m^2))))*(exp(-((y-x+m*sigma)^2)/(2*(1-m^2)*sigma^2)) + exp(-((y-x-m*sigma)^2)/(2*(1-m^2)*sigma^2))) 
    out[y<0]= 0
    return(out)
  }
  
  frails = c()
  fixed_effs = c()
  rand_effs = c()
  
  for (j in 1:numstrap){
    
    master.list <- c()
    
    for (i in 1:numcases){
      if (family=="Truncated Normal"){
        lambda_n <- 1#rtruncnorm(1,a=0,b=Inf,mean=params[1],sd=params[2])
        lambda_tb <- rtruncnorm(1,a=0,b=Inf,mean=params[1],sd=params[2])
      } else if (family=="Exponential"){
        lambda_n <- 1#rexp(1, rate=params[1])
        lambda_tb <- rexp(1, rate=params[1])
      } else if (family=='Uniform'){
        lambda_n = 1#runif(1,params[1],params[2])
        lambda_tb = runif(1,params[1],params[2])
      } else if (family=='Truncated Cauchy'){
        lambda_n = 1#rtrunc(1,'cauchy',a=0,b=Inf,location=params[1],scale=params[2])
        lambda_tb = rtrunc(1,'cauchy',a=0,b=Inf,location=params[1],scale=params[2])
      } else if (family=='Gamma'){
        lambda_n = 1#rgamma(1,shape=params[1],scale=params[2])
        lambda_tb = rgamma(1,shape=params[1],scale=params[2])
      } else if (family=='Bactrian'){
        lambda_n = 1#sample(vals,1,prob=dbactrian(vals,x=params[1],sigma=params[2],m=params[3]))
        lambda_tb = sample(vals,1,prob=dbactrian(vals,x=params[1],sigma=params[2],m=params[3]))
      } else if (family=="Mixture of Gammas"){
        randnum <- runif(1)
        randnum2 <- runif(1)
        
        if (randnum <= 0.5) {
          lambda_n <- 1#rgamma(1,shape=params[1],scale=params[2])
        } else {
          lambda_n <- 1#rgamma(1,shape=params[3],scale=params[4])
        }
        
        if (randnum2 <= 0.5) {
          lambda_tb <- rgamma(1,shape=params[1],scale=params[2])
        } else {
          lambda_tb <- rgamma(1,shape=params[3],scale=params[4])
        }
      } else if (family=='Mixture of Uniforms'){
        randnum = runif(1)
        randnum2 = runif(1)
        if (randnum<=0.5){
          lambda_n = 1#runif(1,params[1],params[2])
        } else{
          lambda_n = 1#runif(1,params[1],params[2])
        }
        
        if (randnum2<=0.5){
          lambda_tb = runif(1,params[1],params[2])
        } else{
          lambda_tb = runif(1,params[3],params[4])
        }
        
      } else{
        stop("Use a valid family!")
      }
      
      num_TB <- 0
      total_time <- 0
      j <- 1
      data.list <- c()
      
      while ((total_time<OBS_TIME)&num_TB<2){
        if (num_TB == 0){
          #print("First Infect")
          time <- rexp(1,lambda_tb)
          if (time + total_time>=OBS_TIME){
            event <- 0
            time <- OBS_TIME - total_time
            total_time <- OBS_TIME
          } else {
            event <- 1
            total_time <- total_time + time
          }
        } else {
          #print("Secondary Infect")
          time <- rexp(1,theta*lambda_tb)
          if (time + total_time>=OBS_TIME){
            event <- 0
            time <- OBS_TIME - total_time
            total_time <- OBS_TIME
          } else {
            event <- 1
            total_time <- total_time + time
          }
        }
        data.vec <- c("ID"=i,"dt"=time,"X"=num_TB>0,"Y"=event)
        data.list[[j]] <- data.vec
        j <- j + 1
        if (event == 1){
          num_TB <- num_TB + 1
        }
      }
      master.list <- c(master.list,data.list)
      
      
    }
    
    resultdfr <- as.data.frame(do.call("rbind", master.list)) 
    resultdfr = resultdfr[resultdfr$X<=1,]
    
    
    mBase = coxph(Surv(dt,Y)~X,data=resultdfr,control=coxph.control(iter.max=5,eps=1e-6))
    mFixed = coxph(Surv(dt,Y)~X+as.factor(ID),data=resultdfr,control=coxph.control(iter.max=5,eps=1e-6))
    mRand = coxme(Surv(dt,Y)~X+(1|ID),data=resultdfr,control=coxme.control(iter.max=5,eps=1e-6))
    mFrail = coxph(Surv(dt,Y)~X+frailty(ID,distribution='gamma'),data=resultdfr,control=coxph.control(iter.max=5,eps=1e-6))
    
    frails = c(frails,mFrail$frail)
    fixed_effs = c(fixed_effs,coef(mFixed)[2:length(coef(mFixed))])
    rand_effs = c(rand_effs,mRand$frail$ID)
  
    
    
    # estBase <- exp(coef(mBase)['X'])
    # estFixed <- exp(coef(mFixed)['X'])
    # estRand <- exp(coef(mRand)['X'])
    # estFrail <- exp(coef(mFrail)['X'])
    # 
    # obs_vecBase <- c(obs_vecBase,estBase)
    # obs_vecFixed <- c(obs_vecFixed,estFixed)
    # obs_vecRand <- c(obs_vecRand,estRand)
    # obs_vecFrail <- c(obs_vecFrail,estFrail)
    
    #  if (length(mFrail$frail)>0){
    #    frails[k,j,1:length(mFrail$frail)] = exp(mFrail$frail)
    #  }
    
  }
  
  # meansBase <- c(meansBase,mean(obs_vecBase,na.rm=T))
  # meansFixed <- c(meansFixed,mean(obs_vecFixed,na.rm=T))
  # meansRand <- c(meansRand,mean(obs_vecRand,na.rm=T))
  # meansFrail <- c(meansFrail,mean(obs_vecFrail,na.rm=T))
  
  return(cbind(fixed_effs,rand_effs,frails))
  #  return(frails)
}



set.seed(1)
thetas = seq(0.01,1,0.01)
k = c(1/8,1/4,1,4,8)
badGamma = goodGamma = list()
for (i in 1:5){
  a = as.numeric(Sys.time())
  goodGamma[[i]] = unlist(mclapply(thetas, function(thetas) good_odds_plot(thetas,
                                                                           family='Gamma',params=c(1/k[i],1*k[i]),negmult=1,numstrap=250,numcases=200,obstime=5),
                                   mc.cores=20))
  print(c(i,'good',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  badGamma[[i]] = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                               family='Gamma',params=c(1/k[i],1*k[i]),negmult=1,numstrap=250,numcases=200,obstime=5),
                                         mc.cores=20)),nrow=100,byrow=T)
  print(c(i,'bad',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  save(goodGamma,file='goodGamma.Rdata'); save(badGamma,file='badGamma.Rdata')
  print(c(i,'saved',as.numeric(Sys.time())-a))
}

nums
xs = seq(0.01,1,0.01)
plot(y=badGamma[[5]][,4],x=xs,type='l',xlim=c(0,1),ylim=c(0,2)); lines(y=xs,x=xs,col='red')
plot(y=goodGamma[[5]],x=xs,type='l',xlim=c(0,1),ylim=c(0,2)); lines(y=xs,x=xs,col='red')
# 
 i = 5
 set.seed(1)
 temp1 = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                      family='Gamma',params=c(1/k[i],1*k[i]),negmult=1,numstrap=50,numcases=200,obstime=5),
                                mc.cores=20)),nrow=100,byrow=T)
 set.seed(2)
 temp2 = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                      family='Gamma',params=c(1/k[i],1*k[i]),negmult=1,numstrap=50,numcases=200,obstime=5),
                                mc.cores=20)),nrow=100,byrow=T)
 set.seed(3)
 temp3 = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                      family='Gamma',params=c(1/k[i],1*k[i]),negmult=1,numstrap=50,numcases=200,obstime=5),
                                mc.cores=20)),nrow=100,byrow=T)
 set.seed(8)
 temp4 = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                      family='Gamma',params=c(1/k[i],1*k[i]),negmult=1,numstrap=50,numcases=200,obstime=5),
                                mc.cores=20)),nrow=100,byrow=T)
 set.seed(9)
 temp5 = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                      family='Gamma',params=c(1/k[i],1*k[i]),negmult=1,numstrap=50,numcases=200,obstime=5),
                                mc.cores=20)),nrow=100,byrow=T)
# 
# badGamma[[5]] = array(NA,dim=dim(badGamma[[4]]))
# for (i in 1:dim(badGamma[[5]])[1]) for (j in 1:dim(badGamma[[5]])[2]){
#   badGamma[[5]][i,j] = mean(c(temp1[i,j],temp2[i,j],temp3[i,j],temp4[i,j],temp5[i,j]))
# }
# save(badGamma,file='badGamma.Rdata')


set.seed(1)
thetas = seq(0.01,1,0.01)
sds = c(0.25,0.5,1,2,4)
badTnorm = goodTnorm = list()
for (i in 1:length(sds)){
  a = as.numeric(Sys.time())
  goodTnorm[[i]] = unlist(mclapply(thetas, function(thetas) good_odds_plot(thetas,
                                                                           family='Truncated Normal',params=c(1,sds[i]),negmult=1,numstrap=250,numcases=200,obstime=5),
                                   mc.cores=20))
  print(c(i,'good',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  badTnorm[[i]] = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                               family='Truncated Normal',params=c(1,sds[i]),negmult=1,numstrap=250,numcases=200,obstime=5),
                                         mc.cores=20)),nrow=100,byrow=T)
  print(c(i,'bad',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  save(goodTnorm,file='goodTnorm.Rdata'); save(badTnorm,file='badTnorm.Rdata')
  print(c(i,'saved',as.numeric(Sys.time())-a))
}
xs = seq(0.01,1,0.01)
plot(y=badTnorm[[4]][,4],x=xs,type='l',xlim=c(0,1),ylim=c(0,1.5)); lines(y=xs,x=xs,col='red')
plot(y=goodTnorm[[5]],x=xs,type='l',xlim=c(0,1),ylim=c(0,1.5)); lines(y=xs,x=xs,col='red')


set.seed(1)
thetas = seq(0.01,1,0.01)
vals = c(0.125,0.25,0.5,0.75,1)
badUnif = goodUnif = list()
for (i in 1:length(vals)){
  a = as.numeric(Sys.time())
  goodUnif[[i]] = unlist(mclapply(thetas, function(thetas) good_odds_plot(thetas,
                                                                          family='Uniform',params=1+c(-1,1)*vals[i],negmult=1,numstrap=250,numcases=200,obstime=5),
                                  mc.cores=20))
  print(c(i,'good',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  badUnif[[i]] = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                              family='Uniform',params=1+c(-1,1)*vals[i],negmult=1,numstrap=250,numcases=200,obstime=5),
                                        mc.cores=20)),nrow=100,byrow=T)
  print(c(i,'bad',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  save(goodUnif,file='goodUnif.Rdata'); save(badUnif,file='badUnif.Rdata')
  print(c(i,'saved',as.numeric(Sys.time())-a))
}
xs = seq(0.01,1,0.01)
plot(y=badUnif[[5]][,4],x=xs,type='l',xlim=c(0,1),ylim=c(0,1.5)); lines(y=xs,x=xs,col='red')
plot(y=goodUnif[[5]],x=xs,type='l',xlim=c(0,1),ylim=c(0,1.5)); lines(y=xs,x=xs,col='red')



#a1 = 2; b1 = 4; a2 = 300/2; b2 = 0.01
mean(c(rgamma(1e4,shape=a1,scale=b1),rgamma(1e4,shape=a2,scale=b2)))
hist(c(rgamma(1e4,shape=a1,scale=b1),rgamma(1e4,shape=a2,scale=b2)),col='red')
sd(c(rgamma(1e4,shape=a1,scale=b1),rgamma(1e4,shape=a2,scale=b2)))


set.seed(1)
thetas = seq(0.01,1,0.01)

params = rbind(c(1/8,1,(15/8)/5,5),
               c(1/8,1,15/8,1),
               c(1/8,1,(15/8)/0.5,0.5),
               c(1/8,1,(15/8)/0.05,0.05),
               c(1/8,1,(15/8)/0.005,0.005))


badMixed = goodMixed = list()
for (i in 1:length(vals)){
  a = as.numeric(Sys.time())
  goodMixed[[i]] = unlist(mclapply(thetas, function(thetas) good_odds_plot(thetas,
                                                                           family='Mixture of Gammas',params=params[i,],negmult=1,numstrap=250,numcases=200,obstime=5),
                                   mc.cores=20))
  print(c(i,'good',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  badMixed[[i]] = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                               family='Mixture of Gammas',params=params[i,],negmult=1,numstrap=250,numcases=200,obstime=5),
                                         mc.cores=20)),nrow=100,byrow=T)
  print(c(i,'bad',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  save(goodMixed,file='goodMixed.Rdata'); save(badMixed,file='badMixed.Rdata')
  print(c(i,'saved',as.numeric(Sys.time())-a))
}




set.seed(1)
thetas = seq(0.01,1,0.01)
scales = c(0.25,0.5,1,2,4)
badTcauchy = goodTcauchy = list()
for (i in 1:length(vals)){
  a = as.numeric(Sys.time())
  goodTcauchy[[i]] = unlist(mclapply(thetas, function(thetas) good_odds_plot(thetas,
                                                                             family='Truncated Cauchy',params=c(1,scales[i]),negmult=1,numstrap=250,numcases=200,obstime=5),
                                     mc.cores=20))
  print(c(i,'good',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  badTcauchy[[i]] = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                                 family='Truncated Cauchy',params=c(1,scales[i]),negmult=1,numstrap=250,numcases=200,obstime=5),
                                           mc.cores=20)),nrow=100,byrow=T)
  print(c(i,'bad',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  save(goodTcauchy,file='goodTcauchy.Rdata'); save(badTcauchy,file='badTcauchy.Rdata')
  print(c(i,'saved',as.numeric(Sys.time())-a))
}
plot(y=badTcauchy[[1]][,3],x=xs,type='l',xlim=c(0,1),ylim=c(0,1.5)); lines(y=xs,x=xs,col='red')
plot(y=goodTcauchy[[1]],x=xs,type='l',xlim=c(0,1),ylim=c(0,1.5)); lines(y=xs,x=xs,col='red')



set.seed(1)
thetas = seq(0.01,1,0.01)

params = rbind(c(1/2,1,(3/2)/5,5),
               c(1/2,1,3/2,1),
               c(1/2,1,(3/2)/0.5,0.5),
               c(1/2,1,(3/2)/0.05,0.05),
               c(1/2,1,(3/2)/0.005,0.005))


badMixed2 = goodMixed2 = list()
for (i in 1:length(vals)){
  a = as.numeric(Sys.time())
  goodMixed2[[i]] = unlist(mclapply(thetas, function(thetas) good_odds_plot(thetas,
                                                                            family='Mixture of Gammas',params=params[i,],negmult=1,numstrap=250,numcases=200,obstime=5),
                                    mc.cores=20))
  print(c(i,'good',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  badMixed2[[i]] = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
                                                                                family='Mixture of Gammas',params=params[i,],negmult=1,numstrap=250,numcases=200,obstime=5),
                                          mc.cores=20)),nrow=100,byrow=T)
  print(c(i,'bad',as.numeric(Sys.time())-a))
  
  a = as.numeric(Sys.time())
  save(goodMixed2,file='goodMixed2.Rdata'); save(badMixed2,file='badMixed2.Rdata')
  print(c(i,'saved',as.numeric(Sys.time())-a))
}

j = 1; k = 4; plot(y=badMixed2[[j]][,k],x=xs,type='l',xlim=c(0,1),ylim=c(0,3)); lines(y=xs,x=xs,col='red'); lines(y=badMixed[[j]][,k],x=xs,col='blue')
plot(y=goodMixed2[[5]],x=xs,type='l',xlim=c(0,1),ylim=c(0,1.5)); lines(y=xs,x=xs,col='red')

params = rbind(c(1/2,1,(3/2)/5,5),
               c(1/2,1,3/2,1),
               c(1/2,1,(3/2)/0.5,0.5),
               c(1/2,1,(3/2)/0.05,0.05),
               c(1/2,1,(3/2)/0.005,0.005))

# set.seed(1)
# i = 1
# temp1 = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
#                                                                      family='Mixture of Gammas',params=params[i,],negmult=1,numstrap=100,numcases=200,obstime=5),
#                                mc.cores=20)),nrow=100,byrow=T)
# set.seed(2)
# temp2 = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
#                                                                      family='Mixture of Gammas',params=params[i,],negmult=1,numstrap=100,numcases=200,obstime=5),
#                                mc.cores=20)),nrow=100,byrow=T)
# set.seed(3)
# temp3 = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
#                                                                      family='Mixture of Gammas',params=params[i,],negmult=1,numstrap=100,numcases=500,obstime=5),
#                                mc.cores=20)),nrow=100,byrow=T)
# set.seed(4)
# temp4 = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
#                                                                      family='Mixture of Gammas',params=params[i,],negmult=1,numstrap=100,numcases=500,obstime=5),
#                                mc.cores=20)),nrow=100,byrow=T)
# set.seed(5)
# temp5 = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
#                                                                      family='Mixture of Gammas',params=params[i,],negmult=1,numstrap=100,numcases=500,obstime=5),
#                                mc.cores=20)),nrow=100,byrow=T)
# 
# badMixed2[[1]] = array(NA,dim=dim(badMixed2[[4]]))
# for (i in 1:dim(badMixed2[[1]])[1]) for (j in 1:dim(badMixed2[[1]])[2]){
#   badMixed2[[1]][i,j] = mean(c(temp1[i,j],temp2[i,j],temp3[i,j],temp4[i,j],temp5[i,j]))
# }
# save(badMixed2,file='badMixed2.Rdata')

vals = seq(0,10,0.01)
dbactrian = function(y,x,sigma,m){
  # if (y<0){
  #   out = 0
  # } else{
  out = (1/(2*sigma*sqrt(2*pi*(1-m^2))))*(exp(-((y-x+m*sigma)^2)/(2*(1-m^2)*sigma^2)) + exp(-((y-x-m*sigma)^2)/(2*(1-m^2)*sigma^2))) 
  # }
  out[y<0]= 0
  return(out)
}


### gamma [[5]]
### gamma mixed2 [[1]]

#set.seed(1)
#thetas = seq(0.01,1,0.01)
#ms = c(0.86,0.89,0.92,0.95,0.98)
#badBact = goodBact = list()
#for (i in 1:length(ms)){
#  a = as.numeric(Sys.time())
#  goodBact[[i]] = unlist(mclapply(thetas, function(thetas) good_odds_plot(thetas,
#                                                                          family='Bactrian',params=c(1,0.35,ms[i]),negmult=1,numstrap=250,numcases=500,obstime=5),
#                                  mc.cores=20))
#  print(c(i,'good',as.numeric(Sys.time())-a))
#  
#  a = as.numeric(Sys.time())
#  badBact[[i]] = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
#                                                                              family='Bactrian',params=c(1,0.35,ms[i]),negmult=1,numstrap=250,numcases=500,obstime=5),
#                                        mc.cores=20)),nrow=100,byrow=T)
#  print(c(i,'bad',as.numeric(Sys.time())-a))
#  
#  a = as.numeric(Sys.time())
#  save(goodBact,file='goodBact.Rdata'); save(badBact,file='badBact.Rdata')
#  print(c(i,'saved',as.numeric(Sys.time())-a))
#}


#k1 = 0.25; k2 = 0.5; hist(c(runif(1e4,0.25-k1,0.25+k1),runif(1e4,1.75-k2,1.75+k2)),breaks=1e2)
#1/8
#k1 = c(0.01,0.02,0.05,0.0625,0.125)
#k2 = rep(15/8,5)
#
#set.seed(1)
#thetas = seq(0.01,1,0.01)
#
#badMixUnif = goodMixUnif = list()
#for (i in 1:length(ms)){
#  a = as.numeric(Sys.time())
#  goodMixUnif[[i]] = unlist(mclapply(thetas, function(thetas) good_odds_plot(thetas,
#                                                                             family='Mixture of Uniforms',params=c(1/8-k1[i],1/8+k1[i],15/8-k2[i],15/8+k2[i]),negmult=1,numstrap=250,numcases=500,obstime=5),
#                                     mc.cores=20))
#  print(c(i,'good',as.numeric(Sys.time())-a))
#  
#  a = as.numeric(Sys.time())
#  badMixUnif[[i]] = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
#                                                                                 family='Mixture of Uniforms',params=c(1/8-k1[i],1/8+k1[i],15/8-k2[i],15/8+k2[i]),negmult=1,numstrap=250,numcases=500,obstime=5),
#                                           mc.cores=20)),nrow=100,byrow=T)
#  print(c(i,'bad',as.numeric(Sys.time())-a))
#  
#  a = as.numeric(Sys.time())
#  #save(goodMixUnif,file='goodMixUnif.Rdata'); save(badMixUnif,file='badMixUnif.Rdata')
#  print(c(i,'saved',as.numeric(Sys.time())-a))
#}
#
#plot(y=goodMixUnif[[1]],x=xs,type='l'); lines(y=xs,x=xs,col='red')
#
#k = 5; plot(y=badMixUnif[[k]][,4],x=xs,type='l',ylim=c(0,2)); lines(y=xs,x=xs,col='red'); lines(y=goodMixUnif[[k]],x=xs,col='blue')


#k1 = c(0.01,0.02,0.05,0.0625,0.125)
#k2 = rep(3/2,5)
#
#set.seed(1)
#thetas = seq(0.01,1,0.01)
#
#badMixUnif2 = goodMixUnif2 = list()
#for (i in 1:length(ms)){
#  a = as.numeric(Sys.time())
#  goodMixUnif2[[i]] = unlist(mclapply(thetas, function(thetas) good_odds_plot(thetas,
#                                                                              family='Mixture of Uniforms',params=c(1/2-k1[i],1/2+k1[i],3/2-k2[i],3/2+k2[i]),negmult=1,numstrap=250,numcases=500,obstime=5),
#                                      mc.cores=20))
#  print(c(i,'good',as.numeric(Sys.time())-a))
#  
#  a = as.numeric(Sys.time())
#  badMixUnif2[[i]] = matrix(unlist(mclapply(thetas, function(thetas) alt_analysis(thetas,
#                                                                                  family='Mixture of Uniforms',params=c(1/2-k1[i],1/2+k1[i],3/2-k2[i],3/2+k2[i]),negmult=1,numstrap=250,numcases=500,obstime=5),
#                                            mc.cores=20)),nrow=100,byrow=T)
#  print(c(i,'bad',as.numeric(Sys.time())-a))
#  
#  a = as.numeric(Sys.time())
#  save(goodMixUnif2,file='goodMixUnif2.Rdata'); save(badMixUnif2,file='badMixUnif2.Rdata')
#  print(c(i,'saved',as.numeric(Sys.time())-a))
#}
#
#
#k = 5; plot(y=badMixUnif[[k]][,4],x=xs,type='l',ylim=c(0,2)); lines(y=xs,x=xs,col='red'); lines(y=goodMixUnif[[k]],x=xs,col='blue')




### dist #### normal bias ### corrected (3 cols)
# FE, RE, frailty

### truncnorm
### cauchy
### Unif
### Unif mix 1
### Unif mix 2
### Gamma mix 1
### Gamma mix 2

load('goodTnorm.Rdata'); load('goodTcauchy.Rdata'); load('goodUnif.Rdata'); load('goodGamma.Rdata'); load('goodMixed.Rdata'); load('goodMixed2.Rdata')
load('badTnorm.Rdata'); load('badTcauchy.Rdata'); load('badUnif.Rdata'); load('badGamma.Rdata'); load('badMixed.Rdata'); load('badMixed2.Rdata')

goodObj = list(goodTnorm,goodTcauchy,goodUnif,goodGamma,goodMixed,goodMixed2)
badObj = list(badTnorm,badTcauchy,badUnif,badGamma,badMixed,badMixed2)

sds = c(0.25,0.5,1,2,4) ### sds (truncnorm)
vals = c(0.125,0.25,0.5,0.75,1) ## unif 
ks = c(1/8,1/4,1,4,8) ## gamma
paramsG1 = rbind(c(1/8,1,(15/8)/5,5),
                 c(1/8,1,15/8,1),
                 c(1/8,1,(15/8)/0.5,0.5),
                 c(1/8,1,(15/8)/0.05,0.05),
                 c(1/8,1,(15/8)/0.005,0.005)) ### gamma Mix
paramsG2 = rbind(c(1/2,1,(3/2)/5,5),
                 c(1/2,1,3/2,1),
                 c(1/2,1,(3/2)/0.5,0.5),
                 c(1/2,1,(3/2)/0.05,0.05),
                 c(1/2,1,(3/2)/0.005,0.005)) ### gamma Mix 2
k1u2 = c(0.01,0.02,0.05,0.0625,0.125); k2u2 = rep(3/2,5) ### unif mix 2
k1u1 = c(0.01,0.02,0.05,0.0625,0.125); k2u1 = rep(15/8,5) #### unif mix 
scales = c(0.25,0.5,1,2,4) ### cauchy

xs = seq(0.01,10,0.01)
out = list()
out[[1]] = out[[2]] = out[[3]] = out[[4]] = out[[5]] = out[[6]] = array(NA,dim=c(length(xs),5))
for (i in 1:5){
  out[[1]][,i] = dtruncnorm(xs,a=0,b=Inf,mean=1,sd=sds[i])
  out[[2]][,i] = dtrunc(xs,spec='cauchy',a=0,b=Inf,location=1,scale=scales[i])
  out[[3]][,i] = dunif(xs,1-vals[i],1+vals[i])
  out[[4]][,i] = dgamma(xs,shape=1/ks[i],scale=1*ks[i])
  out[[5]][,i] = dgamma(xs,shape=paramsG1[i,1],scale=paramsG1[i,2]) + dgamma(xs,shape=paramsG1[i,3],scale=paramsG1[i,4])
  out[[6]][,i] = dgamma(xs,shape=paramsG2[i,1],scale=paramsG2[i,2]) + dgamma(xs,shape=paramsG2[i,3],scale=paramsG2[i,4])
}




pdf("good plots.pdf",width=3,height=6.5)
layout(matrix(1:18,nrow=6,ncol=3,byrow=T))
thetas = seq(0.01,1,0.01)
xs = seq(0.01,10,0.01)

par(mar=c(2,1.85,1,1))
par(tck=-0.02)
par(mgp=c(3,0.3,0))
by = c(rep(0.5,2),1,0.5,1,1)

ylims1 = c(1.55,1.55,4,1.55,5,5)
ylims2 = c(rep(1.25,3),3,3,2.5)

count = 0

for (k in c(1:6)){
  
  
  count = count+1  
  plot(y=out[[k]][,1],x=xs,ylim=c(0,ylims1[k]),type='n',axes=F,ann=F,xlim=c(0,4));
  abline(v=1,col='dark grey',lwd=0.5)
  for (i in 1:5){
    lines(y=out[[k]][,i],col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,x=xs)
  }
  box(bty='l',lwd=0.5)
  axis(side=1,at=0:4,labels=NA,lwd=0,lwd.ticks=0.5)
  axis(side=2,at=seq(0,5,by[k]),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
  text(x=0:4,y=-0.1*ylims1[k],xpd=T,srt=45,0:4,adj=1,cex=0.65)
  mtext(side=2,line=1.125,'Density',cex=0.5)
  mtext(side=1,line=1,expression(lambda[italic(i)]^italic(P)),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[count],')',sep=''),cex=0.5,adj=0)
  
  count = count+1  
  plot(y=badObj[[k]][[1]][,1],x=thetas,type='n',axes=F,ann=F,ylim=c(0,ylims2[k])); abline(h=1,col='dark grey',lwd=0.5)
  for (i in 1:5){
    lines(y=badObj[[k]][[i]][,1],x=thetas,col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,xpd=T)
  }
  lines(y=c(0,1),x=c(0,1),col='dark grey',lwd=0.5)
  box(bty='l',lwd=0.5)
  axis(side=1,at=seq(0,1,0.25),labels=NA,lwd=0,lwd.ticks=0.5)
  labs = seq(0,3,0.25); if (k>3){labs[seq(2,length(labs),2)] = NA}
  axis(side=2,at=seq(0,3,0.25),labels=labs,las=1,cex.axis=0.65,lwd=0,lwd.ticks=0.5)
  text(x=seq(0,1,0.25),y=-0.1*ylims2[k],xpd=T,srt=45,seq(0,1,0.25),adj=1,cex=0.65)
  mtext(side=2,line=1.375,expression(hat(theta)[Naive]),cex=0.5)
  mtext(side=1,line=1,expression(theta),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[count],')',sep=''),cex=0.5,adj=0)
  
  
  count = count+1  
  plot(y=goodObj[[k]][[1]],x=thetas,type='n',axes=F,ann=F,ylim=c(0,1.25)); abline(h=1,col='dark grey',lwd=0.5)
  for (i in 1:5){
    lines(y=goodObj[[k]][[i]],x=thetas,col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,xpd=T)
  }
  lines(y=c(0,1),x=c(0,1),col='dark grey',lwd=0.5)
  box(bty='l',lwd=0.5)
  axis(side=1,at=seq(0,1,0.25),labels=NA,lwd=0,lwd.ticks=0.5)
  labs = seq(0,3,0.25);
  axis(side=2,at=seq(0,3,0.25),labels=labs,las=1,cex.axis=0.65,lwd=0,lwd.ticks=0.5)
  text(x=seq(0,1,0.25),y=-0.1*1.25,xpd=T,srt=45,seq(0,1,0.25),adj=1,cex=0.65)
  mtext(side=2,line=1.375,expression(hat(theta)[MH]),cex=0.5)
  mtext(side=1,line=1,expression(theta),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[count],')',sep=''),cex=0.5,adj=0)
}
dev.off()




pdf("bad plots.pdf",width=4,height=6.5)
layout(matrix(1:24,nrow=6,ncol=4,byrow=T))
thetas = seq(0.01,1,0.01)
xs = seq(0.01,10,0.01)

par(mar=c(2,1.85,1,1))
par(tck=-0.02)
par(mgp=c(3,0.3,0))
by = c(rep(0.5,2),1,0.5,1,1)

ylims1 = c(1.55,1.55,4,1.55,5,5)
ylims2 = c(1.25,1.25,1.25,2,2,1.25)

count = 0

for (k in c(1:6)){
  
  count = count+1  
  plot(y=out[[k]][,1],x=xs,ylim=c(0,ylims1[k]),type='n',axes=F,ann=F,xlim=c(0,4));
  abline(v=1,col='dark grey',lwd=0.5)
  for (i in 1:5){
    lines(y=out[[k]][,i],col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,x=xs)
  }
  box(bty='l',lwd=0.5)
  axis(side=1,at=0:4,labels=NA,lwd=0,lwd.ticks=0.5)
  axis(side=2,at=seq(0,5,by[k]),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
  text(x=0:4,y=-0.1*ylims1[k],xpd=T,srt=45,0:4,adj=1,cex=0.65)
  mtext(side=2,line=1.125,'Density',cex=0.5)
  mtext(side=1,line=1,expression(lambda[italic(i)]^italic(P)),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[count],')',sep=''),cex=0.5,adj=0)
  
  for (j in c(2,3,1)){
    count = count+1  
    plot(y=badObj[[k]][[1]][,j+1],x=thetas,type='n',axes=F,ann=F,ylim=c(0,ylims2[k])); abline(h=1,col='dark grey',lwd=0.5)
    for (i in 1:5){
      lines(y=badObj[[k]][[i]][,j+1],x=thetas,col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,xpd=T)
    }
    lines(y=c(0,1),x=c(0,1),col='dark grey',lwd=0.5)
    box(bty='l',lwd=0.5)
    axis(side=1,at=seq(0,1,0.25),labels=NA,lwd=0,lwd.ticks=0.5)
    labs = seq(0,3,0.25); if (k%in%c(4,5)){labs[seq(2,length(labs),2)] = NA}
    axis(side=2,at=seq(0,3,0.25),labels=labs,las=1,cex.axis=0.65,lwd=0,lwd.ticks=0.5)
    text(x=seq(0,1,0.25),y=-0.1*ylims2[k],xpd=T,srt=45,seq(0,1,0.25),adj=1,cex=0.65)
    
    if (j==1){
      mtext(side=2,line=1.375,expression(hat(theta)[FE]),cex=0.5)      
    } else{
      if (j==2){
        mtext(side=2,line=1.375,expression(hat(theta)[RE]),cex=0.5)
      } else{
    mtext(side=2,line=1.375,expression(hat(theta)[Frailty]),cex=0.5)
      }
    }

    
    mtext(side=1,line=1,expression(theta),cex=0.5)
    mtext(side=3,line=0,paste(LETTERS[count],')',sep=''),cex=0.5,adj=0)    
  }

}
dev.off()


pdf("all est plots.pdf",width=6,height=7)
layout(matrix(1:36,nrow=6,ncol=6,byrow=T))
thetas = seq(0.01,1,0.01)
xs = seq(0.01,10,0.01)

par(mar=c(2,1.85,2,0.5))
par(tck=-0.02)
par(mgp=c(3,0.3,0))
by = c(rep(0.5,2),1,0.5,1,1)
ylimsA = c(1.55,1.55,4,1.55,5,5)

ylims1 = c(1,1,1,4,3,2); bys1 = c(0.25,0.25,0.25,1,0.5,0.5)
ylims2 = c(1,1,1,2,2.5,1.5); bys2 =c(0.25,0.25,0.25,0.5,0.5,0.3)
ylims3 = c(1,1,1,1.2,2,1.2); bys3 =c(0.25,0.25,0.25,0.3,0.5,0.3)

count = 0

ords = list(1:5,1:5,1:5,1:5,5:1,5:1)

distnames = c('A) Truncated Normal distribution',
              'B) Truncated Cauchy distribution',
              'C) Uniform distribution',
              'D) Gamma distribution',
              'E) Gamma mixture distribtuion (i)',
              'F) Gamma mixture distribution (ii)')
badObj[[5]][[1]][1,2] = NA
for (k in c(1:6)){
  
  
  count = count+1  
  plot(y=out[[k]][,1],x=xs,ylim=c(0,ylimsA[k]),type='n',axes=F,ann=F,xlim=c(0,4));
  abline(v=1,col='dark grey',lwd=0.5)
  for (i in 1:5){
    lines(y=out[[k]][,ords[[k]][i]],col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,x=xs)
  }
  box(bty='l',lwd=0.5)
  axis(side=1,at=0:4,labels=NA,lwd=0,lwd.ticks=0.5)
  axis(side=2,at=seq(0,5,by[k]),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
  text(x=0:4,y=-0.1*ylimsA[k],xpd=T,srt=45,0:4,adj=1,cex=0.65)
  mtext(side=2,line=1.125,'Density',cex=0.5)
  mtext(side=1,line=1,expression(lambda[italic(i)]^italic(P)),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[k],'.1)',sep=''),cex=0.5,adj=0)
  mtext(side=3,line=1,distnames[k],cex=0.5,font=2,adj=0,xpd=T)
  
  count = count+1  
  yUB = ylims1[k]
  plot(y=badObj[[k]][[1]][,1],x=thetas,type='n',axes=F,ann=F,ylim=c(0,yUB)); abline(h=1,col='dark grey',lwd=0.5)
  for (i in 1:5){
    lines(y=badObj[[k]][[ords[[k]][i]]][,1],x=thetas,col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,xpd=T)
  }
  lines(y=c(0,1),x=c(0,1),col='dark grey',lwd=0.5)
  box(bty='l',lwd=0.5)
  axis(side=1,at=seq(0,1,0.25),labels=NA,lwd=0,lwd.ticks=0.5)
  labs = seq(0,ylims1[k],bys1[k]); #if (k>3){labs[seq(2,length(labs),2)] = NA}
  axis(side=2,at=labs,labels=labs,las=1,cex.axis=0.65,lwd=0,lwd.ticks=0.5)
  text(x=seq(0,1,0.25),y=-0.1*yUB,xpd=T,srt=45,seq(0,1,0.25),adj=1,cex=0.65)
  mtext(side=2,line=1.25,expression(hat(theta)[Naive]),cex=0.5)
  mtext(side=1,line=1,expression(theta),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[k],'.2)',sep=''),cex=0.5,adj=0)
  
  valspar = c('.5)','.3)','.4)')
  
  for (j in c(2,3,1)){
    count = count+1  
    if (j==2){
      yUB = ylims2[k]
      bys = bys2[k]
    } else{
      if (j==3){
        yUB = ylims3[k]
        bys = bys3[k]
      } else {
        yUB = 1
        bys = 0.25 
      }
    }
    plot(y=badObj[[k]][[1]][,j+1],x=thetas,type='n',axes=F,ann=F,ylim=c(0,yUB)); abline(h=1,col='dark grey',lwd=0.5)
    for (i in 1:5){
      lines(y=badObj[[k]][[ords[[k]][i]]][,j+1],x=thetas,col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,xpd=T)
    }
    lines(y=c(0,1),x=c(0,1),col='dark grey',lwd=0.5)
    box(bty='l',lwd=0.5)
    axis(side=1,at=seq(0,1,0.25),labels=NA,lwd=0,lwd.ticks=0.5)
    labs = seq(0,yUB,bys); #if (k%in%c(4:6)){labs[seq(2,length(labs),2)] = NA}
    axis(side=2,at=labs,labels=labs,las=1,cex.axis=0.65,lwd=0,lwd.ticks=0.5)
    text(x=seq(0,1,0.25),y=-0.1*yUB,xpd=T,srt=45,seq(0,1,0.25),adj=1,cex=0.65)
    
    if (j==1){
      mtext(side=2,line=1.25,expression(hat(theta)[FE]),cex=0.5)      
    } else{
      if (j==2){
        mtext(side=2,line=1.25,expression(hat(theta)[RE]),cex=0.5)
      } else{
        mtext(side=2,line=1.25,expression(hat(theta)[Frailty]),cex=0.5)
      }
    }
    mtext(side=1,line=1,expression(theta),cex=0.5)
    mtext(side=3,line=0,paste(LETTERS[k],valspar[j],sep=''),cex=0.5,adj=0)    
  }
  
  count = count+1  
  plot(y=goodObj[[k]][[1]],x=thetas,type='n',axes=F,ann=F,ylim=c(0,1)); abline(h=1,col='dark grey',lwd=0.5)
  for (i in 1:5){
    lines(y=goodObj[[k]][[ords[[k]][i]]],x=thetas,col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,xpd=T)
  }
  lines(y=c(0,1),x=c(0,1),col='dark grey',lwd=0.5)
  box(bty='l',lwd=0.5)
  axis(side=1,at=seq(0,1,0.25),labels=NA,lwd=0,lwd.ticks=0.5)
  labs = seq(0,3,0.25);
  axis(side=2,at=seq(0,3,0.25),labels=labs,las=1,cex.axis=0.65,lwd=0,lwd.ticks=0.5)
  text(x=seq(0,1,0.25),y=-0.1*1,xpd=T,srt=45,seq(0,1,0.25),adj=1,cex=0.65)
  mtext(side=2,line=1.25,expression(hat(theta)[MH]),cex=0.5)
  mtext(side=1,line=1,expression(theta),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[k],'.6)',sep=''),cex=0.5,adj=0)
}
dev.off()




sds = c(0.25,0.5,1,2,4) ### sds (truncnorm)
vals = c(0.125,0.25,0.5,0.75,1) ## unif 
ks = c(1/6,1/3,1,3,6) ## gamma
paramsG1 = rbind(c(1/8,1,(15/8)/5,5),
                 c(1/8,1,15/8,1),
                 c(1/8,1,(15/8)/0.5,0.5),
                 c(1/8,1,(15/8)/0.05,0.05),
                 c(1/8,1,(15/8)/0.005,0.005)) ### gamma Mix
paramsG2 = rbind(c(1/2,1,(3/2)/5,5),
                 c(1/2,1,3/2,1),
                 c(1/2,1,(3/2)/0.5,0.5),
                 c(1/2,1,(3/2)/0.05,0.05),
                 c(1/2,1,(3/2)/0.005,0.005)) ### gamma Mix 2
k1u2 = c(0.01,0.02,0.05,0.0625,0.125); k2u2 = rep(3/2,5) ### unif mix 2
k1u1 = c(0.01,0.02,0.05,0.0625,0.125); k2u1 = rep(15/8,5) #### unif mix 
scales = c(0.25,0.5,1,2,4) ### cauchy


set.seed(1)
thetas = c(0.2,0.5,0.8)
frails = list()
frails[[1]] = frails[[2]] = frails[[3]] = frails[[4]] = frails[[5]] = frails[[6]] = list()
for (i in 1:5){
  frails[[1]][[i]] = frails[[2]][[i]] = frails[[3]][[i]] = frails[[4]][[i]] = frails[[5]][[i]] = frails[[6]][[i]] = list()
  for (j in 1:3){
    frails[[1]][[i]][[j]] = alt_analysis(theta=thetas[j],family='Truncated Normal',params=c(1,sds[i]),negmult=1,numstrap=250,numcases=500,obstime=5)
    print(c(i,1))
    frails[[2]][[i]][[j]] = alt_analysis(theta=thetas[j],family='Truncated Cauchy',params=c(1,scales[i]),negmult=1,numstrap=250,numcases=500,obstime=5)
    print(c(i,2))
    frails[[3]][[i]][[j]] = alt_analysis(theta=thetas[j],family='Uniform',params=c(1-vals[i],1+vals[i]),negmult=1,numstrap=250,numcases=500,obstime=5)
    print(c(i,3))
    frails[[4]][[i]][[j]] = alt_analysis(theta=thetas[j],family='Gamma',params=c(1/ks[i],1*ks[i]),negmult=1,numstrap=250,numcases=500,obstime=5)
    print(c(i,4))
    frails[[5]][[i]][[j]] = alt_analysis(theta=thetas[j],family='Mixture of Gammas',params=paramsG1[i,],negmult=1,numstrap=250,numcases=500,obstime=5)
    print(c(i,5))
    frails[[6]][[i]][[j]] = alt_analysis(theta=thetas[j],family='Mixture of Gammas',params=paramsG2[i,],negmult=1,numstrap=250,numcases=500,obstime=5)
  }
  print(c(i,6))
}
hist(exp(frails[[6]][[1]]))

frailDens = list()
for (k in 1:6){
  frailDens[[k]] = list()
  for (i in 1:5){
    frailDens[[k]][[i]] = list()
    for (j in 1:3){
      temp = density((frails[[k]][[i]][[j]]),bw=0.5)
      frailDens[[k]][[i]][[j]] = cbind(exp(temp$x),temp$y)
    }
  }
}
plot(frailDens[[4]][[1]][[1]],type='l')

###get effects
sds = c(0.25,0.5,1,2,4) ### sds (truncnorm)
vals = c(0.125,0.25,0.5,0.75,1) ## unif 
ks = c(1/6,1/3,1,3,6) ## gamma
paramsG1 = rbind(c(1/8,1,(15/8)/5,5),
                 c(1/8,1,15/8,1),
                 c(1/8,1,(15/8)/0.5,0.5),
                 c(1/8,1,(15/8)/0.05,0.05),
                 c(1/8,1,(15/8)/0.005,0.005)) ### gamma Mix
paramsG2 = rbind(c(1/2,1,(3/2)/5,5),
                 c(1/2,1,3/2,1),
                 c(1/2,1,(3/2)/0.5,0.5),
                 c(1/2,1,(3/2)/0.05,0.05),
                 c(1/2,1,(3/2)/0.005,0.005)) ### gamma Mix 2
k1u2 = c(0.01,0.02,0.05,0.0625,0.125); k2u2 = rep(3/2,5) ### unif mix 2
k1u1 = c(0.01,0.02,0.05,0.0625,0.125); k2u1 = rep(15/8,5) #### unif mix 
scales = c(0.25,0.5,1,2,4) ### cauchy


set.seed(1)
thetas = c(0.2,0.5,0.8)
effects = list()
effects[[1]] = effects[[2]] = effects[[3]] = effects[[4]] = effects[[5]] = effects[[6]] = list()
for (i in 1:5){
  effects[[1]][[i]] = effects[[2]][[i]] = effects[[3]][[i]] = effects[[4]][[i]] = effects[[5]][[i]] = effects[[6]][[i]] = list()
  for (j in 1:3){
    effects[[1]][[i]][[j]] = get_effects(theta=thetas[j],family='Truncated Normal',params=c(1,sds[i]),negmult=1,numstrap=10,numcases=1000,obstime=5)
    print(c(i,1))
    effects[[2]][[i]][[j]] = get_effects(theta=thetas[j],family='Truncated Cauchy',params=c(1,scales[i]),negmult=1,numstrap=10,numcases=1000,obstime=5)
    print(c(i,2))
    effects[[3]][[i]][[j]] = get_effects(theta=thetas[j],family='Uniform',params=c(1-vals[i],1+vals[i]),negmult=1,numstrap=10,numcases=1000,obstime=5)
    print(c(i,3))
    effects[[4]][[i]][[j]] = get_effects(theta=thetas[j],family='Gamma',params=c(1/ks[i],1*ks[i]),negmult=1,numstrap=10,numcases=1000,obstime=5)
    print(c(i,4))
    effects[[5]][[i]][[j]] = get_effects(theta=thetas[j],family='Mixture of Gammas',params=paramsG1[i,],negmult=1,numstrap=10,numcases=1000,obstime=5)
    print(c(i,5))
    effects[[6]][[i]][[j]] = get_effects(theta=thetas[j],family='Mixture of Gammas',params=paramsG2[i,],negmult=1,numstrap=10,numcases=1000,obstime=5)
  }
  print(c(i,6))
}

save(effects,file='effects.Rdata')

sds = c(0.25,0.5,1,2,4) ### sds (truncnorm)
vals = c(0.125,0.25,0.5,0.75,1) ## unif 
ks = c(1/6,1/3,1,3,6) ## gamma
paramsG1 = rbind(c(1/8,1,(15/8)/5,5),
                 c(1/8,1,15/8,1),
                 c(1/8,1,(15/8)/0.5,0.5),
                 c(1/8,1,(15/8)/0.05,0.05),
                 c(1/8,1,(15/8)/0.005,0.005)) ### gamma Mix
paramsG2 = rbind(c(1/2,1,(3/2)/5,5),
                 c(1/2,1,3/2,1),
                 c(1/2,1,(3/2)/0.5,0.5),
                 c(1/2,1,(3/2)/0.05,0.05),
                 c(1/2,1,(3/2)/0.005,0.005)) ### gamma Mix 2
k1u2 = c(0.01,0.02,0.05,0.0625,0.125); k2u2 = rep(3/2,5) ### unif mix 2
k1u1 = c(0.01,0.02,0.05,0.0625,0.125); k2u1 = rep(15/8,5) #### unif mix 
scales = c(0.25,0.5,1,2,4) ### cauchy


set.seed(1)
thetas = c(0.2,0.5,0.8)
effects2 = list()
effects2[[1]] = effects2[[2]] = effects2[[3]] = effects2[[4]] = effects2[[5]] = effects2[[6]] = list()
for (i in 1:5){
  effects2[[1]][[i]] = effects2[[2]][[i]] = effects2[[3]][[i]] = effects2[[4]][[i]] = effects2[[5]][[i]] = effects2[[6]][[i]] = list()
  for (j in 1:3){
    effects2[[1]][[i]][[j]] = get_effects(theta=thetas[j],family='Truncated Normal',params=c(1,sds[i]),negmult=1,numstrap=10,numcases=1000,obstime=5)
    print(c(i,1))
    effects2[[2]][[i]][[j]] = get_effects(theta=thetas[j],family='Truncated Cauchy',params=c(1,scales[i]),negmult=1,numstrap=10,numcases=1000,obstime=5)
    print(c(i,2))
    effects2[[3]][[i]][[j]] = get_effects(theta=thetas[j],family='Uniform',params=c(1-vals[i],1+vals[i]),negmult=1,numstrap=10,numcases=1000,obstime=5)
    print(c(i,3))
    effects2[[4]][[i]][[j]] = get_effects(theta=thetas[j],family='Gamma',params=c(1/ks[i],1*ks[i]),negmult=1,numstrap=10,numcases=1000,obstime=5)
    print(c(i,4))
    effects2[[5]][[i]][[j]] = get_effects(theta=thetas[j],family='Mixture of Gammas',params=paramsG1[i,],negmult=1,numstrap=10,numcases=1000,obstime=5)
    print(c(i,5))
    effects2[[6]][[i]][[j]] = get_effects(theta=thetas[j],family='Mixture of Gammas',params=paramsG2[i,],negmult=1,numstrap=10,numcases=1000,obstime=5)
  }
  print(c(i,6))
}

save(effects2,file='effects2.Rdata')

fixedDens = list()
for (k in 1:6){
  fixedDens[[k]] = list()
  for (i in 1:5){
    fixedDens[[k]][[i]] = list()
    for (j in 1:3){
      temp = density((effects[[k]][[i]][[j]][[1]]),bw=0.5)
      fixedDens[[k]][[i]][[j]] = cbind(exp(temp$x),temp$y)
    }
  }
}

randDens = list()
for (k in 1:6){
  randDens[[k]] = list()
  for (i in 1:5){
    randDens[[k]][[i]] = list()
    for (j in 1:3){
      temp = density((effects[[k]][[i]][[j]][[2]]),bw=0.5)
      randDens[[k]][[i]][[j]] = cbind(exp(temp$x),temp$y)
    }
  }
}

frailDens = list()
for (k in 1:6){
  frailDens[[k]] = list()
  for (i in 1:5){
    frailDens[[k]][[i]] = list()
    for (j in 1:3){
      temp = density((effects[[k]][[i]][[j]][[3]]),bw=0.5)
      frailDens[[k]][[i]][[j]] = cbind(exp(temp$x),temp$y)
    }
  }
}

save(fixedDens,file='fixedDens.Rdata')
save(randDens,file='randDens.Rdata')
save(frailDens,file='frailDens2.Rdata')


pdf("fixed plots.pdf",width=4,height=7)
layout(matrix(1:24,nrow=6,ncol=4,byrow=F))
thetas = seq(0.01,1,0.01)
xs = seq(0.01,10,0.01)

par(mar=c(2,1.85,2,0.5))
par(tck=-0.02)
par(mgp=c(3,0.3,0))
by = c(rep(0.5,2),1,0.5,1,1)
ylimsA = c(1.55,1.55,4,1.55,5,5)

ylims1 = c(1,1,1,4,3,2); bys1 = c(0.25,0.25,0.25,1,0.5,0.5)
ylims2 = c(1,1,1,2,2.5,1.5); bys2 =c(0.25,0.25,0.25,0.5,0.5,0.3)
ylims3 = c(1,1,1,1.2,2,1.2); bys3 =c(0.25,0.25,0.25,0.3,0.5,0.3)
ylims4 = c(1,1,1,1,0.75,0.75); bys4 = c(0.5,0.5,0.5,0.5,0.25,0.25)

xbys1=c(1,1,1,1,1,5,1,5,5,1,5,5,1,1,1,5,5,1)

xlims1 = c(6,15,10,10,10,50,10,15,30,5,15,80,10,10,15,20,15,5)

count = 0

ords = list(1:5,1:5,1:5,1:5,5:1,5:1)

distnames = c('A) Truncated Normal distribution',
              'B) Truncated Cauchy distribution',
              'C) Uniform distribution',
              'D) Gamma distribution',
              'E) Gamma mixture distribtuion (i)',
              'F) Gamma mixture distribution (ii)')
icons = c('.2)','.3)','.4)')
for(k in 1:6){
  plot(y=out[[k]][,1],x=xs,ylim=c(0,ylimsA[k]),type='n',axes=F,ann=F,xlim=c(0,4));
  abline(v=1,col='dark grey',lwd=0.5)
  for (i in 1:5){
    lines(y=out[[k]][,ords[[k]][i]],col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,x=xs)
  }
  box(bty='l',lwd=0.5)
  axis(side=1,at=0:4,lwd=0,lwd.ticks=0.5,cex.axis=0.65)
  axis(side=2,at=seq(0,5,by[k]),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
  #text(x=0:4,y=-0.1*ylimsA[k],xpd=T,srt=45,0:4,adj=1,cex=0.65)
  mtext(side=2,line=1.125,'Density',cex=0.5)
  mtext(side=1,line=1,expression(lambda[italic(P)][italic(i)]),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[k],'.1)',sep=''),cex=0.5,adj=0)
  mtext(side=3,line=1,distnames[k],cex=0.5,font=2,adj=0,xpd=T)
  
  count = count+1  
  yUB = ylims1[k]
}
for (j in 1:3){
  for (k in 1:6){
    plot(fixedDens[[k]][[1]][[j]],col="red",xlab="frailty",type='l',lwd=0.5,axes=F,ann=F,xlim=c(0,xlims1[k+6*(j-1)]),ylim=c(0,1))
    for (i in 2:5){
      lines(fixedDens[[k]][[i]][[j]],type='l',col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5)
    }
    box(bty='l',lwd=0.5)
    axis(side=1,at=seq(0,100,xbys1[k+6*(j-1)]),lwd=0,lwd.ticks=0.5,cex.axis=0.65)
    axis(side=2,at=seq(0,1,bys4[k]),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
    #text(x=0:xs1[k],y=-0.1*ylimsA[k],xpd=T,srt=45,seq(0,xs1[k],xby1[k]),adj=1,cex=0.65)
    #mtext(side=2,line=1.125,'Density',cex=0.5)
    mtext(side=1,line=1,expression(lambda[italic(P)][italic(i)]),cex=0.5)
    mtext(side=3,line=0,paste(LETTERS[k],icons[j],sep=''),cex=0.5,adj=0)
    
    count = count+1  
    yUB = ylims1[k]
  }
}

dev.off()



pdf("rand plots.pdf",width=4,height=7)
layout(matrix(1:24,nrow=6,ncol=4,byrow=F))
thetas = seq(0.01,1,0.01)
xs = seq(0.01,10,0.01)

par(mar=c(2,1.85,2,0.5))
par(tck=-0.02)
par(mgp=c(3,0.3,0))
by = c(rep(0.5,2),1,0.5,1,1)
ylimsA = c(1.55,1.55,4,1.55,5,5)

ylims1 = c(1,1,1,4,3,2); bys1 = c(0.25,0.25,0.25,1,0.5,0.5)
ylims2 = c(1,1,1,2,2.5,1.5); bys2 =c(0.25,0.25,0.25,0.5,0.5,0.3)
ylims3 = c(1,1,1,1.2,2,1.2); bys3 =c(0.25,0.25,0.25,0.3,0.5,0.3)
ylims4 = c(1,1,1,1,0.75,0.75); bys4 = c(0.5,0.5,0.5,0.5,0.5,0.5)

xbys1=c(1,1,1,1,1,1,1,5,5,1,1,5,1,1,1,5,5,5)

xlims1 = c(6,10,10,5,10,5,3,15,30,5,5,15,6,20,6,40,30,30)

count = 0

ords = list(1:5,1:5,1:5,1:5,5:1,5:1)

distnames = c('A) Truncated Normal distribution',
              'B) Truncated Cauchy distribution',
              'C) Uniform distribution',
              'D) Gamma distribution',
              'E) Gamma mixture distribtuion (i)',
              'F) Gamma mixture distribution (ii)')
icons = c('.2)','.3)','.4)')
for(k in 1:6){
  plot(y=out[[k]][,1],x=xs,ylim=c(0,ylimsA[k]),type='n',axes=F,ann=F,xlim=c(0,4));
  abline(v=1,col='dark grey',lwd=0.5)
  for (i in 1:5){
    lines(y=out[[k]][,ords[[k]][i]],col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,x=xs)
  }
  box(bty='l',lwd=0.5)
  axis(side=1,at=0:4,lwd=0,lwd.ticks=0.5,cex.axis=0.65)
  axis(side=2,at=seq(0,5,by[k]),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
  #text(x=0:4,y=-0.1*ylimsA[k],xpd=T,srt=45,0:4,adj=1,cex=0.65)
  mtext(side=2,line=1.125,'Density',cex=0.5)
  mtext(side=1,line=1,expression(lambda[italic(P)][italic(i)]),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[k],'.1)',sep=''),cex=0.5,adj=0)
  mtext(side=3,line=1,distnames[k],cex=0.5,font=2,adj=0,xpd=T)
  
  count = count+1  
  yUB = ylims1[k]
}
for (j in 1:3){
  for (k in 1:6){
    plot(y=randDens[[k]][[1]][[j]][,2],x=randDens[[k]][[1]][[j]][,1],col="red",xlab="frailty",type='l',lwd=0.5,axes=F,ann=F,xlim=c(0,xlims1[k+6*(j-1)]),ylim=c(0,1))
    for (i in 2:5){
      lines(randDens[[k]][[i]][[j]],type='l',col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5)
    }
    box(bty='l',lwd=0.5)
    axis(side=1,at=seq(0,100,xbys1[k+6*(j-1)]),lwd=0,lwd.ticks=0.5,cex.axis=0.65)
    axis(side=2,at=seq(0,1,bys4[k]),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
    #text(x=0:xs1[k],y=-0.1*ylimsA[k],xpd=T,srt=45,seq(0,xs1[k],xby1[k]),adj=1,cex=0.65)
    #mtext(side=2,line=1.125,'Density',cex=0.5)
    mtext(side=1,line=1,expression(lambda[italic(P)][italic(i)]),cex=0.5)
    mtext(side=3,line=0,paste(LETTERS[k],icons[j],sep=''),cex=0.5,adj=0)
    
    count = count+1  
    yUB = ylims1[k]
  }
}

dev.off()


pdf("frail plots.pdf",width=4,height=7)
layout(matrix(1:24,nrow=6,ncol=4,byrow=F))
thetas = seq(0.01,1,0.01)
xs = seq(0.01,10,0.01)

par(mar=c(2,1.85,2,0.5))
par(tck=-0.02)
par(mgp=c(3,0.3,0))
by = c(rep(0.5,2),1,0.5,1,1)
ylimsA = c(1.55,1.55,4,1.55,5,5)

ylims1 = c(1,1,1,4,3,2); bys1 = c(0.25,0.25,0.25,1,0.5,0.5)
ylims2 = c(1,1,1,2,2.5,1.5); bys2 =c(0.25,0.25,0.25,0.5,0.5,0.3)
ylims3 = c(1,1,1,1.2,2,1.2); bys3 =c(0.25,0.25,0.25,0.3,0.5,0.3)
ylims4 = c(1,1,1,1,0.75,0.75); bys4 = c(0.5,0.5,0.5,0.5,0.5,0.5)

xbys1=c(1,1,1,1,1,5,1,1,1,1,1,5,1,1,1,5,5,1)

xlims1 = c(6,10,10,5,10,20,3,10,10,5,5,50,6,6,6,40,20,5)

count = 0

ords = list(1:5,1:5,1:5,1:5,5:1,5:1)

distnames = c('A) Truncated Normal distribution',
              'B) Truncated Cauchy distribution',
              'C) Uniform distribution',
              'D) Gamma distribution',
              'E) Gamma mixture distribtuion (i)',
              'F) Gamma mixture distribution (ii)')
icons = c('.2)','.3)','.4)')
for(k in 1:6){
  plot(y=out[[k]][,1],x=xs,ylim=c(0,ylimsA[k]),type='n',axes=F,ann=F,xlim=c(0,4));
  abline(v=1,col='dark grey',lwd=0.5)
  for (i in 1:5){
    lines(y=out[[k]][,ords[[k]][i]],col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5,x=xs)
  }
  box(bty='l',lwd=0.5)
  axis(side=1,at=0:4,lwd=0,lwd.ticks=0.5,cex.axis=0.65)
  axis(side=2,at=seq(0,5,by[k]),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
  #text(x=0:4,y=-0.1*ylimsA[k],xpd=T,srt=45,0:4,adj=1,cex=0.65)
  mtext(side=2,line=1.125,'Density',cex=0.5)
  mtext(side=1,line=1,expression(lambda[italic(P)][italic(i)]),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[k],'.1)',sep=''),cex=0.5,adj=0)
  mtext(side=3,line=1,distnames[k],cex=0.5,font=2,adj=0,xpd=T)
  
  count = count+1  
  yUB = ylims1[k]
}
for (j in 1:3){
  for (k in 1:6){
    plot(y=frailDens[[k]][[1]][[j]][,2],x=randDens[[k]][[1]][[j]][,1],col="red",xlab="frailty",type='l',lwd=0.5,axes=F,ann=F,xlim=c(0,xlims1[k+6*(j-1)]),ylim=c(0,1))
    for (i in 2:5){
      lines(frailDens[[k]][[i]][[j]],type='l',col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5)
    }
    box(bty='l',lwd=0.5)
    axis(side=1,at=seq(0,100,xbys1[k+6*(j-1)]),lwd=0,lwd.ticks=0.5,cex.axis=0.65)
    axis(side=2,at=seq(0,1,bys4[k]),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
    #text(x=0:xs1[k],y=-0.1*ylimsA[k],xpd=T,srt=45,seq(0,xs1[k],xby1[k]),adj=1,cex=0.65)
    #mtext(side=2,line=1.125,'Density',cex=0.5)
    mtext(side=1,line=1,expression(lambda[italic(P)][italic(i)]),cex=0.5)
    mtext(side=3,line=0,paste(LETTERS[k],icons[j],sep=''),cex=0.5,adj=0)
    
    count = count+1  
    yUB = ylims1[k]
  }
}

dev.off()


######### power analysis 


simtimes <- function(theta,family,params,numcases){

  res = resRaw = array(NA,dim=c(numcases,4))
  
  vals = seq(0,3,0.005)
  dbactrian = function(y,x,sigma,m){
    out = (1/(2*sigma*sqrt(2*pi*(1-m^2))))*(exp(-((y-x+m*sigma)^2)/(2*(1-m^2)*sigma^2)) + exp(-((y-x-m*sigma)^2)/(2*(1-m^2)*sigma^2))) 
    out[y<0]= 0
    return(out)
  }
  
  
    
    for (i in 1:numcases){
      
      if (family=="Truncated Normal"){
        lambda_n <- 1#1*negmult#rtruncnorm(1,a=0,b=Inf,mean=params[1],sd=params[2])
        lambda_tb <- rtruncnorm(1,a=0,b=Inf,mean=params[1],sd=params[2])
      } else if (family=="Exponential"){
        lambda_n <- 1#1*negmult#rexp(1, rate=params[1])
        lambda_tb <- rexp(1, rate=params[1])
      } else if (family=='Uniform'){
        lambda_n = 1#runif(1,params[1],params[2])
        lambda_tb = runif(1,params[1],params[2])
      } else if (family=='Truncated Cauchy'){
        lambda_n = 1#rtrunc(1,'cauchy',a=0,b=Inf,location=params[1],scale=params[2])
        lambda_tb = rtrunc(1,'cauchy',a=0,b=Inf,location=params[1],scale=params[2])
      } else if (family=='Gamma'){
        lambda_n = 1#rgamma(1,shape=params[1],scale=params[2])
        lambda_tb = rgamma(1,shape=params[1],scale=params[2])
      } else if (family=='Bactrian'){
        lambda_n = 1#sample(vals,1,prob=dbactrian(vals,x=params[1],sigma=params[2],m=params[3]))
        lambda_tb = sample(vals,1,prob=dbactrian(vals,x=params[1],sigma=params[2],m=params[3]))
      } else if (family=="Mixture of Gammas"){
        randnum <- runif(1)
        randnum2 <- runif(1)
        
        if (randnum <= 0.5) {
          lambda_n <- 1#gamma(1,shape=params[1],scale=params[2])
        } else {
          lambda_n <- 1#rgamma(1,shape=params[3],scale=params[4])
        }
        
        
        if (randnum2 <= 0.5) {
          lambda_tb <- rgamma(1,shape=params[1],scale=params[2])
        } else {
          lambda_tb <- rgamma(1,shape=params[3],scale=params[4])
        }
      } else if (family=='Mixture of Uniforms'){
        randnum = runif(1)
        randnum2 = runif(1)
        if (randnum<=0.5){
          lambda_n = 1#runif(1,params[1],params[2])
        } else{
          lambda_n = 1#runif(1,params[1],params[2])
        }
        
        if (randnum2<=0.5){
          lambda_tb = runif(1,params[1],params[2])
        } else{
          lambda_tb = runif(1,params[3],params[4])
        }
      } else{
        stop("Use a valid family!")
      }
      
      A <- rexp(1,lambda_tb)
      B <- rexp(1,lambda_n)
      C <- rexp(1,theta*lambda_tb)
      D <- rexp(1,lambda_n)
      
      resRaw[i,] = c(A,B,C,D)
      res[i,] = c(A,B,A+C,A+D)
    }
  return(list(resRaw,res))
}



simtimesCor <- function(theta,family,params,numcases){
  
  res = resRaw = array(NA,dim=c(numcases,4))
  
  vals = seq(0,3,0.005)
  dbactrian = function(y,x,sigma,m){
    out = (1/(2*sigma*sqrt(2*pi*(1-m^2))))*(exp(-((y-x+m*sigma)^2)/(2*(1-m^2)*sigma^2)) + exp(-((y-x-m*sigma)^2)/(2*(1-m^2)*sigma^2))) 
    out[y<0]= 0
    return(out)
  }
  
  
  
  for (i in 1:numcases){
    
    if (family=="Truncated Normal"){
      #lambda_n <- 1#1*negmult#rtruncnorm(1,a=0,b=Inf,mean=params[1],sd=params[2])
      lambda_tb = lambda_n <- rtruncnorm(1,a=0,b=Inf,mean=params[1],sd=params[2])
    } else if (family=="Exponential"){
      #lambda_n <- 1#1*negmult#rexp(1, rate=params[1])
      lambda_tb = lambda_n <- rexp(1, rate=params[1])
    } else if (family=='Uniform'){
      #lambda_n = 1#runif(1,params[1],params[2])
      lambda_tb = lambda_n = runif(1,params[1],params[2])
    } else if (family=='Truncated Cauchy'){
     # lambda_n = 1#rtrunc(1,'cauchy',a=0,b=Inf,location=params[1],scale=params[2])
      lambda_tb = lambda_n = rtrunc(1,'cauchy',a=0,b=Inf,location=params[1],scale=params[2])
    } else if (family=='Gamma'){
      #lambda_n = 1#rgamma(1,shape=params[1],scale=params[2])
      lambda_tb = lambda_n = rgamma(1,shape=params[1],scale=params[2])
    } else if (family=='Bactrian'){
    #  lambda_n = 1#sample(vals,1,prob=dbactrian(vals,x=params[1],sigma=params[2],m=params[3]))
      lambda_tb = lambda_n = sample(vals,1,prob=dbactrian(vals,x=params[1],sigma=params[2],m=params[3]))
    } else if (family=="Mixture of Gammas"){
      #randnum <- runif(1)
      randnum2 <- runif(1)
      
      #if (randnum <= 0.5) {
      #  lambda_n <- 1#gamma(1,shape=params[1],scale=params[2])
      #} else {
      #  lambda_n <- 1#rgamma(1,shape=params[3],scale=params[4])
      #}
      
      
      if (randnum2 <= 0.5) {
        lambda_tb = lambda_n <- rgamma(1,shape=params[1],scale=params[2])
      } else {
        lambda_tb = lambda_n <- rgamma(1,shape=params[3],scale=params[4])
      }
    } else if (family=='Mixture of Uniforms'){
     # randnum = runif(1)
      randnum2 = runif(1)
     # if (randnum<=0.5){
     #   lambda_n = 1#runif(1,params[1],params[2])
     # } else{
     #   lambda_n = 1#runif(1,params[1],params[2])
     # }
      
      if (randnum2<=0.5){
        lambda_tb = lambda_n = runif(1,params[1],params[2])
      } else{
        lambda_tb = lambda_n = runif(1,params[3],params[4])
      }
    } else{
      stop("Use a valid family!")
    }
    
    A <- rexp(1,lambda_tb)
    B <- rexp(1,lambda_n)
    C <- rexp(1,theta*lambda_tb)
    D <- rexp(1,lambda_n)
    
    resRaw[i,] = c(A,B,C,D)
    res[i,] = c(A,B,A+C,A+D)
  }
  return(list(resRaw,res))
}



set.seed(1)
thetas = c(seq(0.1,0.9,0.1))

sds = c(0.25,0.5,1,2,4) ### sds (truncnorm)
vals = c(0.125,0.25,0.5,0.75,1) ## unif 
ks = c(1/8,1/4,1,4,8) ## gamma
paramsG1 = rbind(c(1/8,1,(15/8)/5,5),
                 c(1/8,1,15/8,1),
                 c(1/8,1,(15/8)/0.5,0.5),
                 c(1/8,1,(15/8)/0.05,0.05),
                 c(1/8,1,(15/8)/0.005,0.005)) ### gamma Mix
paramsG2 = rbind(c(1/2,1,(3/2)/5,5),
                 c(1/2,1,3/2,1),
                 c(1/2,1,(3/2)/0.5,0.5),
                 c(1/2,1,(3/2)/0.05,0.05),
                 c(1/2,1,(3/2)/0.005,0.005)) ### gamma Mix 2
scales = c(0.25,0.5,1,2,4) ### cauchy

timesTnorm = timesUnif = timesTcauchy = timesGamma = timesMixGamma = timesMixGamma2 = list()
for (i in 1:5){
  timesTnorm[[i]] = (mclapply(thetas, function(thetas) simtimes(thetas,
                                                           family='Truncated Normal',
                                                           params=c(1,sds[i]),
                                                           numcases=1e5),
                         mc.cores=20))  
  timesTcauchy[[i]] = (mclapply(thetas, function(thetas) simtimes(thetas,
                                                                family='Truncated Cauchy',
                                                                params=c(1,scales[i]),
                                                                numcases=1e5),
                              mc.cores=20))  
  timesUnif[[i]] = (mclapply(thetas, function(thetas) simtimes(thetas,
                                                                family='Uniform',
                                                                params=1+c(-1,1)*vals[i],
                                                                numcases=1e5),
                              mc.cores=20))  
  timesGamma[[i]] = (mclapply(thetas, function(thetas) simtimes(thetas,
                                                               family='Gamma',
                                                               params=c(1/ks[i],1*ks[i]),
                                                               numcases=1e5),
                             mc.cores=20))  
  timesMixGamma[[i]] = (mclapply(thetas, function(thetas) simtimes(thetas,
                                                                   family='Mixture of Gammas',
                                                                   params=paramsG1[i,],
                                                                   numcases=1e5),
                                 mc.cores=20))  
  timesMixGamma2[[i]] = (mclapply(thetas, function(thetas) simtimes(thetas,
                                                                   family='Mixture of Gammas',
                                                                   params=paramsG2[i,],
                                                                   numcases=1e5),
                                 mc.cores=20))  
  print(i)
}

save(timesTnorm,file='timesTnorm.Rdata'); save(timesTcauchy,file='timesTcauchy.Rdata')
save(timesUnif,file='timesUnif.Rdata'); save(timesGamma,file='timesGamma.Rdata')
save(timesMixGamma,file='timesMixGamma.Rdata'); save(timesMixGamma2,file='timesMixGamma2.Rdata')


set.seed(1)
thetas = c(seq(0.1,0.9,0.1))



sds = c(0.25,0.5,1,2,4) ### sds (truncnorm)
vals = c(0.125,0.25,0.5,0.75,1) ## unif 
ks = c(1/8,1/4,1,4,8) ## gamma
paramsG1 = rbind(c(1/8,1,(15/8)/5,5),
                 c(1/8,1,15/8,1),
                 c(1/8,1,(15/8)/0.5,0.5),
                 c(1/8,1,(15/8)/0.05,0.05),
                 c(1/8,1,(15/8)/0.005,0.005)) ### gamma Mix
paramsG2 = rbind(c(1/2,1,(3/2)/5,5),
                 c(1/2,1,3/2,1),
                 c(1/2,1,(3/2)/0.5,0.5),
                 c(1/2,1,(3/2)/0.05,0.05),
                 c(1/2,1,(3/2)/0.005,0.005)) ### gamma Mix 2
scales = c(0.25,0.5,1,2,4) ### cauchy

timesTnormCor = timesUnifCor = timesTcauchyCor = timesGammaCor = timesMixGammaCor = timesMixGamma2Cor = list()
for (i in 1:5){
  timesTnormCor[[i]] = (mclapply(thetas, function(thetas) simtimesCor(thetas,
                                                                family='Truncated Normal',
                                                                params=c(1,sds[i]),
                                                                numcases=1e5),
                              mc.cores=20))  
  timesTcauchyCor[[i]] = (mclapply(thetas, function(thetas) simtimesCor(thetas,
                                                                  family='Truncated Cauchy',
                                                                  params=c(1,scales[i]),
                                                                  numcases=1e5),
                                mc.cores=20))  
  timesUnifCor[[i]] = (mclapply(thetas, function(thetas) simtimesCor(thetas,
                                                               family='Uniform',
                                                               params=1+c(-1,1)*vals[i],
                                                               numcases=1e5),
                             mc.cores=20))  
  timesGammaCor[[i]] = (mclapply(thetas, function(thetas) simtimesCor(thetas,
                                                                family='Gamma',
                                                                params=c(1/ks[i],1*ks[i]),
                                                                numcases=1e5),
                              mc.cores=20))  
  timesMixGammaCor[[i]] = (mclapply(thetas, function(thetas) simtimesCor(thetas,
                                                                   family='Mixture of Gammas',
                                                                   params=paramsG1[i,],
                                                                   numcases=1e5),
                                 mc.cores=20))  
  timesMixGamma2Cor[[i]] = (mclapply(thetas, function(thetas) simtimesCor(thetas,
                                                                    family='Mixture of Gammas',
                                                                    params=paramsG2[i,],
                                                                    numcases=1e5),
                                  mc.cores=20))  
  print(i)
}
save(timesTnormCor,file='timesTnormCor.Rdata'); save(timesTcauchyCor,file='timesTcauchyCor.Rdata')
save(timesUnifCor,file='timesUnifCor.Rdata'); save(timesGammaCor,file='timesGammaCor.Rdata')
save(timesMixGammaCor,file='timesMixGammaCor.Rdata'); save(timesMixGamma2Cor,file='timesMixGamma2Cor.Rdata')






power.fn = function(nums,timesInput,tmax){
    out = c()
    totTime = timesInput[,1]+timesInput[,3]
    sels = totTime<tmax
    timesInput = timesInput[sels,]
    for (i in 1:2e3){
      samp = sample(1:dim(timesInput)[1],nums,replace=F)
      or = rep(NA,200)
      xi = rep(NA,200)
      for (j in 1:200){
        boots = sample(samp,length(samp),replace=T)
        Ei = timesInput[boots,3]<timesInput[boots,4]
        Fi = timesInput[boots,4]<timesInput[boots,3]
        Gi = timesInput[boots,1]<timesInput[boots,2]
        Hi = timesInput[boots,2]<timesInput[boots,1]
        or[j] = sum(Ei*Hi)/sum(Fi*Gi)
      }
      ub = quantile(or,0.975,na.rm=T)
      out[i] = ub<1
    }
  return(out)
}

#analytic
ANApower.fn = function(nums,timesInput,tmax){
  out = c()
  totTime = timesInput[,1]+timesInput[,3]
  sels = totTime<tmax
  timesInput = timesInput[sels,]
  for (i in 1:2e3){
    samp = sample(1:dim(timesInput)[1],nums,replace=F)
    xi = rep(NA,200)
    for (j in 1:200){
      boots = sample(samp,length(samp),replace=T)
      Ei = timesInput[boots,3]<timesInput[boots,4]
      Fi = timesInput[boots,4]<timesInput[boots,3]
      Gi = timesInput[boots,1]<timesInput[boots,2]
      Hi = timesInput[boots,2]<timesInput[boots,1]
      xi[j] = sum(0.5*(Ei-Gi))^2/(0.25*sum((Ei+Gi)*(2-Ei+Gi)))
    }
    ub = quantile(xi,0.975,na.rm=T)
    out[i] = ub<5.02
  }
  return(out)
}



set.seed(1)
nums = c(25,50,75,100,150,200,250,300,350,400,450,500,600,700,800,1000)
normOut = gammaOut = mixGamma1out = array(NA,dim=c(5,9,16))
for (i in c(1,5)) for (j in 1:9){
  a = as.numeric(Sys.time())
  tempNorm = mclapply(nums, function(nums) power.fn(nums,
                                                timesInput=timesTnorm[[i]][[j]][[1]],10),
                  mc.cores=20)
  
  print(c(as.numeric(Sys.time())-a,'Norm'))
  
  a = as.numeric(Sys.time())
  tempGamma = mclapply(nums, function(nums) power.fn(nums,
                                                    timesInput=timesGamma[[i]][[j]][[1]],10),
                      mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Gamma'))
  
  a = as.numeric(Sys.time())
  tempMix = mclapply(nums, function(nums) power.fn(nums,
                                                    timesInput=timesMixGamma[[i]][[j]][[1]],10),
                      mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Mix'))
  
  for (k in 1:length(nums)){
    normOut[i,j,k] = mean(tempNorm[[k]])
    gammaOut[i,j,k] = mean(tempGamma[[k]])
    mixGamma1out[i,j,k] = mean(tempMix[[k]])
  }
  print(c(i,j))
}


#save(normOut,file='normOut.Rdata')
#save(gammaOut,file='gammaOut.Rdata')
#save(mixGamma1out,file='mixGamma1out.Rdata')#


set.seed(1)
nums = c(25,50,75,100,150,200,250,300,350,400,450,500,600,700,800,1000)
normOutCor = gammaOutCor = mixGamma1outCor = array(NA,dim=c(5,9,16))
for (i in c(1,5)) for (j in 1:9){
  a = as.numeric(Sys.time())
  tempNorm = mclapply(nums, function(nums) power.fn(nums,
                                                    timesInput=timesTnormCor[[i]][[j]][[1]]),
                      mc.cores=20)
  
  print(c(as.numeric(Sys.time())-a,'Norm'))
  
  a = as.numeric(Sys.time())
  tempGamma = mclapply(nums, function(nums) power.fn(nums,
                                                     timesInput=timesGammaCor[[i]][[j]][[1]]),
                       mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Gamma'))
  
  a = as.numeric(Sys.time())
  tempMix = mclapply(nums, function(nums) power.fn(nums,
                                                   timesInput=timesMixGammaCor[[i]][[j]][[1]]),
                     mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Mix'))
  
  for (k in 1:length(nums)){
    normOutCor[i,j,k] = mean(tempNorm[[k]])
    gammaOutCor[i,j,k] = mean(tempGamma[[k]])
    mixGamma1outCor[i,j,k] = mean(tempMix[[k]])
  }
  print(c(i,j))
}

#analytic

load('timesTnorm.Rdata')
load('timesGamma.Rdata')
load('timesMixGamma.Rdata')
load('timesTnormCor.Rdata')
load('timesGammaCor.Rdata')
load('timesMixGammaCor.Rdata')

set.seed(1)
nums = c(25,50,75,100,150,200,250,300,350,400,450,500,600,700,800,1000)
ANAnormOut = ANAgammaOut = ANAmixGamma1out = array(NA,dim=c(5,9,16))
for (i in c(1,5)) for (j in 1:9){
  a = as.numeric(Sys.time())
  tempNorm = mclapply(nums, function(nums) ANApower.fn(nums,
                                                    timesInput=timesTnorm[[i]][[j]][[1]],tmax=100),
                      mc.cores=20)
  
  print(c(as.numeric(Sys.time())-a,'Norm'))
  
  a = as.numeric(Sys.time())
  tempGamma = mclapply(nums, function(nums) ANApower.fn(nums,
                                                     timesInput=timesGamma[[i]][[j]][[1]],tmax=100),
                       mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Gamma'))
  
  a = as.numeric(Sys.time())
  tempMix = mclapply(nums, function(nums) ANApower.fn(nums,
                                                   timesInput=timesMixGamma[[i]][[j]][[1]],tmax=100),
                     mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Mix'))
  
  for (k in 1:length(nums)){
    ANAnormOut[i,j,k] = mean(tempNorm[[k]])
    ANAgammaOut[i,j,k] = mean(tempGamma[[k]])
    ANAmixGamma1out[i,j,k] = mean(tempMix[[k]])
  }
  print(c(i,j))
}

set.seed(1)
nums = c(25,50,75,100,150,200,250,300,350,400,450,500,600,700,800,1000)
ANAnormOutCor = ANAgammaOutCor = ANAmixGamma1outCor = array(NA,dim=c(5,9,16))
for (i in c(1,5)) for (j in 1:9){
  a = as.numeric(Sys.time())
  tempNorm = mclapply(nums, function(nums) ANApower.fn(nums,
                                                    timesInput=timesTnormCor[[i]][[j]][[1]],tmax=100),
                      mc.cores=20)
  
  print(c(as.numeric(Sys.time())-a,'Norm'))
  
  a = as.numeric(Sys.time())
  tempGamma = mclapply(nums, function(nums) ANApower.fn(nums,
                                                     timesInput=timesGammaCor[[i]][[j]][[1]],tmax=100),
                       mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Gamma'))
  
  a = as.numeric(Sys.time())
  tempMix = mclapply(nums, function(nums) ANApower.fn(nums,
                                                   timesInput=timesMixGammaCor[[i]][[j]][[1]],tmax=100),
                     mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Mix'))
  
  for (k in 1:length(nums)){
    ANAnormOutCor[i,j,k] = mean(tempNorm[[k]])
    ANAgammaOutCor[i,j,k] = mean(tempGamma[[k]])
    ANAmixGamma1outCor[i,j,k] = mean(tempMix[[k]])
  }
  print(c(i,j))
}
save(ANAnormOut,file='ANAnormOut2side.Rdata')
save(ANAgammaOut,file='ANAgammaOut2side.Rdata')
save(ANAmixGamma1out,file='ANAmixGamma1out2side.Rdata')#
save(ANAnormOutCor,file='ANAnormOutCor2side.Rdata')
save(ANAgammaOutCor,file='ANAgammaOutCor2side.Rdata')
save(ANAmixGamma1outCor,file='ANAmixGamma1outCor2side.Rdata')#

set.seed(1)
nums = c(25,50,75,100,150,200,250,300,350,400,450,500,600,700,800,1000)
normOut2y = gammaOut2y = mixGamma1out2y = array(NA,dim=c(5,9,16))
for (i in c(1,5)) for (j in 1:9){
  a = as.numeric(Sys.time())
  tempNorm = mclapply(nums, function(nums) power.fn(nums,
                                                    timesInput=timesTnorm[[i]][[j]][[1]],tmax=2),
                      mc.cores=20)
  
  print(c(as.numeric(Sys.time())-a,'Norm'))
  
  a = as.numeric(Sys.time())
  tempGamma = mclapply(nums, function(nums) power.fn(nums,
                                                     timesInput=timesGamma[[i]][[j]][[1]],tmax=2),
                       mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Gamma'))
  
  a = as.numeric(Sys.time())
  tempMix = mclapply(nums, function(nums) power.fn(nums,
                                                   timesInput=timesMixGamma[[i]][[j]][[1]],tmax=2),
                     mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Mix'))
  
  for (k in 1:length(nums)){
    normOut2y[i,j,k] = mean(tempNorm[[k]])
    gammaOut2y[i,j,k] = mean(tempGamma[[k]])
    mixGamma1out2y[i,j,k] = mean(tempMix[[k]])
  }
  print(c(i,j))
}
save(normOut2y,file='normOut2y.Rdata')
save(gammaOut2y,file='gammaOut2y.Rdata')
save(mixGamma1out2y,file='mixGamma1out2y.Rdata')#

set.seed(1)
nums = c(25,50,75,100,150,200,250,300,350,400,450,500,600,700,800,1000)
normOut5y = gammaOut5y = mixGamma1out5y = array(NA,dim=c(5,9,16))
for (i in c(1,5)) for (j in 1:9){
  a = as.numeric(Sys.time())
  tempNorm = mclapply(nums, function(nums) power.fn(nums,
                                                    timesInput=timesTnorm[[i]][[j]][[1]],tmax=5),
                      mc.cores=20)
  
  print(c(as.numeric(Sys.time())-a,'Norm'))
  
  a = as.numeric(Sys.time())
  tempGamma = mclapply(nums, function(nums) power.fn(nums,
                                                     timesInput=timesGamma[[i]][[j]][[1]],tmax=5),
                       mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Gamma'))
  
  a = as.numeric(Sys.time())
  tempMix = mclapply(nums, function(nums) power.fn(nums,
                                                   timesInput=timesMixGamma[[i]][[j]][[1]],tmax=5),
                     mc.cores=20)
  print(c(as.numeric(Sys.time())-a,'Mix'))
  
  for (k in 1:length(nums)){
    normOut5y[i,j,k] = mean(tempNorm[[k]])
    gammaOut5y[i,j,k] = mean(tempGamma[[k]])
    mixGamma1out5y[i,j,k] = mean(tempMix[[k]])
  }
  print(c(i,j))
}
save(normOut5y,file='normOut5y.Rdata')
save(gammaOut5y,file='gammaOut5y.Rdata')
save(mixGamma1out5y,file='mixGamma1out5y.Rdata')#



j = 9; plot(normOut[5,j,],type='l',ylim=c(0,1)); lines(normOutCor[5,j,],col='red')

plotObj = list(normOut,gammaOut,mixGamma1out)
plotCorObj = list(normOutCor,gammaOutCor,mixGamma1outCor)
plot5yObj = list(normOut5y,gammaOut5y,mixGamma1out5y)
plot2yObj = list(normOut2y,gammaOut2y,mixGamma1out2y)


pdf('power plots.pdf',height=3,width=6.5)
layout(matrix(1:12,nrow=2,ncol=6,byrow=T),widths=rep(c(1,0.9),3))
vs = rbind(c(1,5),c(1,5),c(5,1))
par(mgp=c(3,0.3,0)); par(tck=-0.02)
count = 0
for (i in 1:3){
  count = count+1
  for (j in 1:2){
    if (j==1){
      par(mar=c(3,3.5,1,0.5))
    } else{
      par(mar=c(3,1.5,1,1.5))
    }
    plot(y=plotObj[[i]][vs[i,j],1,],x=(nums),type='n',axes=F,ann=F,ylim=c(0,1))
    #polygon(x=c(-100,10000,10000,-100),y=c(0.8,0.8,1,1),col='darkgoldenrod2',lty=0)
    abline(h=seq(0,1,0.2),col='grey',lwd=0.25)
    for (k in 1:9){
      lines(y=plotObj[[i]][vs[i,j],k,],x=(nums),lwd=0.5,col=rgb((k-1)/8,0,1-(k-1)/8))
    }
    box(bty='l',lwd=0.5)
    axis(1,at=seq(0,1e3,200),labels=NA,lwd=0,lwd.ticks=0.5,tck=-0.04)
    axis(1,at=seq(100,900,200),labels=NA,lwd=0,lwd.ticks=0.5,tck=-0.02)
    axis(2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=0.65,lwd=0,lwd.ticks=0.5,las=1)
    text(x=seq(0,1e3,200),y=-0.1,seq(0,1e3,200),xpd=T,srt=45,adj=1,cex=0.65)
    mtext(side=1,'Sample size',cex=0.5,line=1.5)
    mtext(side=2,'Power (%)',cex=0.5,line=1.125)
    mtext(side=3,paste(LETTERS[i],'.',j,')',sep=''),cex=0.5,adj=0,line=0)
    #text(x=1e3,y=1.05,adj=1,expression(italic(t)[max]==infinity),xpd=T,cex=0.65)
    text(x=1e3,y=1.05,adj=1,ifelse(j==1,'Min. variance','Max. variance'),cex=0.65,xpd=T)
    #if (j==1){
    #  text(x=0,y=1.15,adj=0,
    #        ifelse(i==1,'Truncated normal',ifelse(i==2,'Gamma','Gamma mixture')),
    #        font=2,cex=0.65,xpd=T) 
    #   text(x=0,y=1.05,adj=0,expression(italic(t)[max]==infinity),xpd=T,cex=0.65)
    #}
  }
}


for (i in 1:3){
  for (j in 1:2){
    if (j==1){
      par(mar=c(3,3.5,1,0.5))
    } else{
      par(mar=c(3,1.5,1,1.5))
    }
    plot(y=plotCorObj[[i]][vs[i,j],j,],x=(nums),type='n',axes=F,ann=F,ylim=c(0,1))
    #polygon(x=c(-100,10000,10000,-100),y=c(0.8,0.8,1,1),col='darkgoldenrod2',lty=0)
    abline(h=seq(0,1,0.2),col='grey',lwd=0.25)
    for (k in 1:9){
      lines(y=plotCorObj[[i]][vs[i,j],k,],x=(nums),lwd=0.5,col=rgb((k-1)/8,0,1-(k-1)/8))
    }  
    box(bty='l',lwd=0.5)
    axis(1,at=seq(0,1e3,200),labels=NA,lwd=0,lwd.ticks=0.5,tck=-0.04)
    axis(1,at=seq(100,900,200),labels=NA,lwd=0,lwd.ticks=0.5,tck=-0.02)
    axis(2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=0.65,lwd=0,lwd.ticks=0.5,las=1)
    text(x=seq(0,1e3,200),y=-0.1,seq(0,1e3,200),xpd=T,srt=45,adj=1,cex=0.65)
    mtext(side=1,'Sample size',cex=0.5,line=1.5)
    mtext(side=2,'Power (%)',cex=0.5,line=1.125)
    mtext(side=3,paste(LETTERS[i],'.',j+2,')',sep=''),cex=0.5,adj=0)
    #text(x=1e3,y=1.05,adj=1,expression(italic(t)[max]==5),xpd=T,cex=0.65)
    text(x=1e3,y=1.05,adj=1,ifelse(j==1,'Min. variance','Max. variance'),cex=0.65,xpd=T)
  }
}


#for (i in 1:3){
#  for (j in 1:2){
#    if (j==1){
#      par(mar=c(2.5,3,2,0.5))
#    } else{
#      par(mar=c(2.5,2,2,0.5))
#    }
#    plot(y=plot2yObj[[i]][vs[i,j],j,],x=(nums),type='n',axes=F,ann=F,ylim=c(0,1))
#    #polygon(x=c(-100,10000,10000,-100),y=c(0.8,0.8,1,1),col='darkgoldenrod2',lty=0)
#    abline(h=seq(0,1,0.2),col='grey',lwd=0.25)
#    for (k in 1:9){
#      lines(y=plot2yObj[[i]][vs[i,j],k,],x=(nums),lwd=0.5,col=rgb((k-1)/8,0,1-(k-1)/8))
#    }  
#    box(bty='l',lwd=0.5)
#    axis(1,at=seq(0,1e3,200),labels=NA,lwd=0,lwd.ticks=0.5,tck=-0.04)
#    axis(1,at=seq(100,900,200),labels=NA,lwd=0,lwd.ticks=0.5,tck=-0.02)
#    axis(2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=0.65,lwd=0,lwd.ticks=0.5,las=1)
#    text(x=seq(0,1e3,200),y=-0.1,seq(0,1e3,200),xpd=T,srt=45,adj=1,cex=0.65)
#    mtext(side=1,'Sample size',cex=0.5,line=1.5)
#    mtext(side=2,'Power (%)',cex=0.5,line=1.125)
#    mtext(side=3,paste(LETTERS[i],'.',j+4,')',sep=''),cex=0.5,adj=0)
#    text(x=1e3,y=1.05,adj=1,expression(italic(t)[max]==2),xpd=T,cex=0.65)
#    text(x=1e3,y=1.15,adj=1,ifelse(j==1,'Min. variance','Max. variance'),cex=0.65,xpd=T)
#  }
#}

j = 9; plot(ANAnormOut[5,j,],type='l',ylim=c(0,1)); lines(ANAnormOutCor[5,j,],col='red')

ANAplotObj = list(ANAnormOut,ANAgammaOut,ANAmixGamma1out)
ANAplotCorObj = list(ANAnormOutCor,ANAgammaOutCor,ANAmixGamma1outCor)

pdf('ANA power plots2.pdf',height=3,width=6.5)
layout(matrix(1:12,nrow=2,ncol=6,byrow=T),widths=rep(c(1,0.9),3))
vs = rbind(c(1,5),c(1,5),c(5,1))
par(mgp=c(3,0.3,0)); par(tck=-0.02)
count = 0
for (i in 1:3){
  count = count+1
  for (j in 1:2){
    if (j==1){
      par(mar=c(3,3.5,1,0.5))
    } else{
      par(mar=c(3,1.5,1,1.5))
    }
    plot(y=1-ANAplotObj[[i]][vs[i,j],1,],x=(nums),type='n',axes=F,ann=F,ylim=c(0,1))
    #polygon(x=c(-100,10000,10000,-100),y=c(0.8,0.8,1,1),col='darkgoldenrod2',lty=0)
    abline(h=seq(0,1,0.2),col='grey',lwd=0.25)
    for (k in 1:9){
      lines(y=1-ANAplotObj[[i]][vs[i,j],k,],x=(nums),lwd=0.5,col=rgb((k-1)/8,0,1-(k-1)/8))
    }
    box(bty='l',lwd=0.5)
    axis(1,at=seq(0,1e3,200),labels=NA,lwd=0,lwd.ticks=0.5,tck=-0.04)
    axis(1,at=seq(100,900,200),labels=NA,lwd=0,lwd.ticks=0.5,tck=-0.02)
    axis(2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=0.65,lwd=0,lwd.ticks=0.5,las=1)
    text(x=seq(0,1e3,200),y=-0.1,seq(0,1e3,200),xpd=T,srt=45,adj=1,cex=0.65)
    mtext(side=1,'Sample size',cex=0.5,line=1.5)
    mtext(side=2,'Power (%)',cex=0.5,line=1.125)
    mtext(side=3,paste(LETTERS[i],'.',j,')',sep=''),cex=0.5,adj=0,line=0)
    #text(x=1e3,y=1.05,adj=1,expression(italic(t)[max]==infinity),xpd=T,cex=0.65)
    text(x=1e3,y=1.05,adj=1,ifelse(j==1,'Min. variance','Max. variance'),cex=0.65,xpd=T)
    #if (j==1){
    #  text(x=0,y=1.15,adj=0,
    #        ifelse(i==1,'Truncated normal',ifelse(i==2,'Gamma','Gamma mixture')),
    #        font=2,cex=0.65,xpd=T) 
    #   text(x=0,y=1.05,adj=0,expression(italic(t)[max]==infinity),xpd=T,cex=0.65)
    #}
  }
}


for (i in 1:3){
  for (j in 1:2){
    if (j==1){
      par(mar=c(3,3.5,1,0.5))
    } else{
      par(mar=c(3,1.5,1,1.5))
    }
    plot(y=1-ANAplotCorObj[[i]][vs[i,j],j,],x=(nums),type='n',axes=F,ann=F,ylim=c(0,1))
    #polygon(x=c(-100,10000,10000,-100),y=c(0.8,0.8,1,1),col='darkgoldenrod2',lty=0)
    abline(h=seq(0,1,0.2),col='grey',lwd=0.25)
    for (k in 1:9){
      lines(y=1-ANAplotCorObj[[i]][vs[i,j],k,],x=(nums),lwd=0.5,col=rgb((k-1)/8,0,1-(k-1)/8))
    }  
    box(bty='l',lwd=0.5)
    axis(1,at=seq(0,1e3,200),labels=NA,lwd=0,lwd.ticks=0.5,tck=-0.04)
    axis(1,at=seq(100,900,200),labels=NA,lwd=0,lwd.ticks=0.5,tck=-0.02)
    axis(2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=0.65,lwd=0,lwd.ticks=0.5,las=1)
    text(x=seq(0,1e3,200),y=-0.1,seq(0,1e3,200),xpd=T,srt=45,adj=1,cex=0.65)
    mtext(side=1,'Sample size',cex=0.5,line=1.5)
    mtext(side=2,'Power (%)',cex=0.5,line=1.125)
    mtext(side=3,paste(LETTERS[i],'.',j+2,')',sep=''),cex=0.5,adj=0)
    #text(x=1e3,y=1.05,adj=1,expression(italic(t)[max]==5),xpd=T,cex=0.65)
    text(x=1e3,y=1.05,adj=1,ifelse(j==1,'Min. variance','Max. variance'),cex=0.65,xpd=T)
  }
}

dev.off()
