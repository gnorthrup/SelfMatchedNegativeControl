setwd("~/Desktop/berkeley/Lewnard Lab/Non-Parametric Estimation NAI/Rotavirus") 

library(DescToolsAddIns)
load('baseProtect.Rdata')
load('protect1prior.Rdata')
load('protect2prior.Rdata')

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

xs <- c(0.6,1.2,1.8)
points <- c(baseProtect[2],protect1prior[2],protect2prior[2])
lowers <- c(baseProtect[1],protect1prior[1],protect2prior[1])
uppers <- c(baseProtect[3],protect1prior[3],protect2prior[3])

V_points <- c(0.23,0.17,0.08)
V_lowers <- c(0.12,0.08,0.01)
V_uppers <- c(0.40,0.36,0.56)

G_points <- c(0.57,0.29,0.19)
G_lowers <- c(0.44,0.20,0.12)
G_uppers <- c(0.76,0.41,0.31)

L_points <- c(0.67,0.50,0.38)
L_lowers <- c(0.59,0.43,0.30)
L_uppers <- c(0.77,0.58,0.45)
# plot(x, points,
#      ylim=c(0, 2.5),
#      pch=19, xlab="Number of prior infections", ylab="OR of protection from next inf.",
#      main="RV Protection"
# )
# # hack: we draw arrows but with very special "arrowheads"
# arrows(x, lowers, x, uppers, length=0.05, angle=90, code=3)
# axis(side=1,at=0:4,labels=NA,lwd=0,lwd.ticks=0.5)
pdf("RotaRatiosColor.pdf")
par(mar=c(5, 4.1, 4.1, 2.1))
plot(y=points,x=xs-0.1,ylim=c(0,1),type='p',pch=19,axes=F,ann=F,xlim=c(0,2),col=cbbPalette[1])
points(0,1)
abline(h=1,lty='dotted',lwd=0.5)
box(bty='l',lwd=0.5)
axis(side=1,at=c(0,0.65,1.25,1.85),labels=NA,lwd=0,lwd.ticks=0.5)
text(x=c(0,0.65,1.25,1.85),y=-0.1,xpd=T,srt=45,c(0,1,2,3),adj=1,cex=0.65)
axis(side=2,at=seq(0,2.5,0.5),las=1,cex.axis=0.65,lwd.ticks=0.5,lwd=0)
#text(x=0:4,y=-0.1*1.5,xpd=T,srt=45,0:4,adj=1,cex=0.65)
mtext(side=2,line=2,expression(hat(theta)[MH]),cex=1)
mtext(side=1,line=2,'Number of Prior Infections',cex=0.8)
mtext(side=3,adj=0,line=0.5,'Protective effect of previous infections',cex=1)
arrows(xs-0.1, lowers,xs-0.1,uppers, length=0.05, angle=90, code=3,lwd=2,col=cbbPalette[1])

points(y=L_points,x=xs,pch=17,col=cbbPalette[2])
arrows(xs,L_lowers,xs,L_uppers,length=0.05,angle=90,code=3,col=cbbPalette[2],lwd=2)
points(x=-0.1,y=-0.2,xpd=T,pch=17,col=cbbPalette[2])
text(x=0.18,y=-0.2,xpd=T,"Lewnard et al. 2019",cex=0.65)

points(y=G_points,x=xs+0.1,pch=15,col=cbbPalette[3])
arrows(xs+0.1,G_lowers,xs+0.1,G_uppers,length=0.05,angle=90,code=3,col=cbbPalette[3],lwd=2)
points(x=0.6,y=-0.2,xpd=T,pch=15,col=cbbPalette[3])
text(x=0.9,y=-0.2,xpd=T,"Gladstone et al. 2011",cex=0.65)

points(y=V_points,x=xs+0.2,pch=18,cex=1.5,col=cbbPalette[4])
arrows(xs+0.2,V_lowers,xs+0.2,V_uppers,length=0.05,angle=90,code=3,col=cbbPalette[4],lwd=2)
points(x=1.3,y=-0.2,xpd=T,pch=15,col=cbbPalette[4])
text(x=1.6,y=-0.2,xpd=T,"Velazquez et al. 1996",cex=0.65)

dev.off()

