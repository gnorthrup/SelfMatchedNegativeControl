setwd("~/Desktop/berkeley/Lewnard Lab/Non-Parametric Estimation NAI/Google Drive/new frails")
load("effects.Rdata")

library(truncdist)
library(parallel)
library(coxme)
library(survival)
library(truncnorm)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

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

fixedDens = list()
fixedInf = list()
for (k in 1:6){
  fixedDens[[k]] = list()
  fixedInf[[k]] = list()
  for (i in 1:5){
    fixedDens[[k]][[i]] = list()
    fixedInf[[k]][[i]] = list()
    for (j in 1:3){
      temp = density(effects[[k]][[i]][[j]][which(exp(effects[[k]][[i]][[j]][,1])<1000),1],bw=0.5)
      fixedDens[[k]][[i]][[j]] = cbind(exp(temp$x),temp$y)
      fixedInf[[k]][[i]][[j]] = length(which(exp(effects[[k]][[i]][[j]][,1])>1000))/10000
    }
  }
}

randDens = list()
for (k in 1:6){
  randDens[[k]] = list()
  for (i in 1:5){
    randDens[[k]][[i]] = list()
    for (j in 1:3){
      temp = density((effects[[k]][[i]][[j]][,2]),bw=0.5)
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
      temp = density((effects[[k]][[i]][[j]][,3]),bw=0.5)
      frailDens[[k]][[i]][[j]] = cbind(exp(temp$x),temp$y)
    }
  }
}

pdf("fixed plots test.pdf",width=5,height=7)
layout(matrix(1:30,nrow=6,ncol=5,byrow=F))
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
    #plot(fixedDens[[k]][[1]][[j]],col="red",xlab="frailty",type='l',lwd=0.5,axes=F,ann=F,xlim=c(0,10),ylim=c(0,1))
    plot(y=fixedDens[[k]][[1]][[j]][,2],x=fixedDens[[k]][[1]][[j]][,1],col="red",xlab="frailty",type='l',lwd=0.5,axes=F,ann=F,xlim=c(0,10),ylim=c(0,1))
    abline(v=1,col='dark grey',lwd=0.5)
    for (i in 2:5){
      lines(fixedDens[[k]][[i]][[j]],type='l',col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5)
    }
    box(bty='l',lwd=0.5)
    axis(side=1,at=seq(0,10,2),lwd=0,lwd.ticks=0.5,cex.axis=0.65)
    axis(side=2,at=seq(0,1,0.5),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
    #text(x=0:xs1[k],y=-0.1*ylimsA[k],xpd=T,srt=45,seq(0,xs1[k],xby1[k]),adj=1,cex=0.65)
    #mtext(side=2,line=1.125,'Density',cex=0.5)
    mtext(side=1,line=1,expression(lambda[italic(P)][italic(i)]),cex=0.5)
    mtext(side=3,line=0,paste(LETTERS[k],icons[j],sep=''),cex=0.5,adj=0)
    
    count = count+1  
    yUB = ylims1[k]
  }
}

for (k in 1:6){
  plot(y=fixedInf[[k]][[1]],x=c(0.2,0.5,0.8),type="p",pch=4,col="red",ylim=c(0,0.2),xlim=c(0,1),lwd=0.5,axes=F,ann=F)
  for (i in 2:5){
    points(y=fixedInf[[k]][[i]],x=c(0.2,0.5,0.8),type="p",pch=4,col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5)
  }
  box(bty='l',lwd=0.5)
  axis(side=1,at=c(0.2,0.5,0.8),lwd=0,lwd.ticks=0.5,cex.axis=0.65)
  axis(side=2,at=seq(0,0.2,0.1),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
  mtext(side=2,line=1.125,'Fraction',cex=0.5)
  mtext(side=1,line=1,expression(theta),cex=0.5)
  mtext(side=3,line=0,paste(LETTERS[k],'.5)',sep=''),cex=0.5,adj=0)
}

dev.off()



pdf("rand plots test.pdf",width=4,height=7)
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
    plot(y=randDens[[k]][[1]][[j]][,2],x=randDens[[k]][[1]][[j]][,1],col="red",xlab="frailty",type='l',lwd=0.5,axes=F,ann=F,xlim=c(0,10),ylim=c(0,1))
    abline(v=1,col='dark grey',lwd=0.5)
    for (i in 2:5){
      lines(randDens[[k]][[i]][[j]],type='l',col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5)
    }
    box(bty='l',lwd=0.5)
    axis(side=1,at=seq(0,10,2),lwd=0,lwd.ticks=0.5,cex.axis=0.65)
    axis(side=2,at=seq(0,1,0.5),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
    #text(x=0:xs1[k],y=-0.1*ylimsA[k],xpd=T,srt=45,seq(0,xs1[k],xby1[k]),adj=1,cex=0.65)
    #mtext(side=2,line=1.125,'Density',cex=0.5)
    mtext(side=1,line=1,expression(lambda[italic(P)][italic(i)]),cex=0.5)
    mtext(side=3,line=0,paste(LETTERS[k],icons[j],sep=''),cex=0.5,adj=0)
    
    count = count+1  
    yUB = ylims1[k]
  }
}

dev.off()


pdf("frail plots test.pdf",width=4,height=7)
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
    plot(y=frailDens[[k]][[1]][[j]][,2],x=frailDens[[k]][[1]][[j]][,1],col="red",xlab="frailty",type='l',lwd=0.5,axes=F,ann=F,xlim=c(0,10),ylim=c(0,1))
    abline(v=1,col='dark grey',lwd=0.5)
    for (i in 2:5){
      lines(frailDens[[k]][[i]][[j]],type='l',col=rgb(1-(i-1)/4,0,(i-1)/4,1),lwd=0.5)
    }
    box(bty='l',lwd=0.5)
    axis(side=1,at=seq(0,10,2),lwd=0,lwd.ticks=0.5,cex.axis=0.65)
    axis(side=2,at=seq(0,1,0.5),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
    #text(x=0:xs1[k],y=-0.1*ylimsA[k],xpd=T,srt=45,seq(0,xs1[k],xby1[k]),adj=1,cex=0.65)
    #mtext(side=2,line=1.125,'Density',cex=0.5)
    mtext(side=1,line=1,expression(lambda[italic(P)][italic(i)]),cex=0.5)
    mtext(side=3,line=0,paste(LETTERS[k],icons[j],sep=''),cex=0.5,adj=0)
    
    count = count+1  
    yUB = ylims1[k]
  }
}

dev.off()

