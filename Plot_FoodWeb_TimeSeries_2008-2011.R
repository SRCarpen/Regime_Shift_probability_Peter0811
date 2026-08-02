# Plot food web time series from Squeal 1
# 2026-02-22

rm(list = ls())
graphics.off()

library('colorspace')

# load data from Combine_Raft+Zoop+Fish_08-11_squeal1_v1_2026-02-10.R
print('saving all covariates + t-stats for Squeal 1 logistic regression',quote=F)
#save(all4reg,file='CombineCovarSplineSqueal1.Rdata')
load(file='CombineCovarSplineSqueal1.Rdata')

print(all4reg[1,])

dat0 = all4reg

print(c('data dim ',dim(dat0)),quote=F)

# extract some variates and transform food web
Tindx = dat0$Tindx
Tprob = dat0$Tprob 
Tstat = dat0$Tstat
zmix = dat0$zmix
medtemp = dat0$medtemp
dosat50 = dat0$dosat50
Chl = sqrt(dat0$chl)
ZB = log10(dat0$ZB)
Dbar = dat0$Dbar
pelratio = log10(dat0$pelratio + 1)
litmin = log10(dat0$litmin)
totmin = log10(dat0$totminnow)
lmb = log10(dat0$lmb)

# original units, not transformed
# lmb, totmin+litmin,ZB,Dbar,Chl
#
windows(height=12,width=8)
par(mfrow=c(6,1),mar=c(3,4.2,0.3,1)+0.1,cex.lab=1.8,cex.axis=1.8)
plot(Tindx,dat0$lmb,type='l',lwd=2,col='black',
     ylab='Bass Pop.',xlab='')
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')
plot(Tindx,dat0$litmin,type='l',lwd=2,col='darkgreen',
     ylab='Minnow CPE',xlab='')
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')
#points(Tindx,dat0$litmin,type='l',lwd=2,col='lightseagreen')
plot(Tindx,dat0$ZB,type='l',lwd=2,col='lightseagreen',
     ylab='Zoopl. Biom.',xlab='')
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')
plot(Tindx,dat0$Dbar,type='l',lwd=2,col='mediumblue',
     ylab='Daph. length',xlab='')
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')
LG = adjust_transparency('limegreen',alpha=0.5)
plot(Tindx,dat0$chl,type='l',lwd=1,col='black',
     ylab='Chl ug/L',xlab='')
points(Tindx,dat0$chl,type='l',lwd=3,col=LG)
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')
par(mar=c(4,4.2,0.2,1)+0.1)
plot(Tindx,dat0$Tprob,type='l',lwd=1,col='black',
     ylab='p(Low Chl)',xlab='Year')
points(Tindx,dat0$Tprob,type='l',lwd=3,col='deepskyblue')
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')

# omit probabilities (last panel)
windows(height=12,width=8)
par(mfrow=c(5,1),mar=c(3,4.2,0.3,1)+0.1,cex.lab=1.8,cex.axis=1.8)
plot(Tindx,dat0$lmb,type='l',lwd=2,col='black',
     ylab='Bass Pop.',xlab='')
text(x=2008.2,y=130,'A',cex=1.8)
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')
plot(Tindx,dat0$litmin,type='l',lwd=2,col='darkgreen',
     ylab='Minnow CPE',xlab='')
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')
text(x=2011.8,y=19,'B',cex=1.8)
#points(Tindx,dat0$litmin,type='l',lwd=2,col='lightseagreen')
plot(Tindx,dat0$ZB,type='l',lwd=2,col='lightseagreen',
     ylab='Zoopl. Biom.',xlab='')
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')
text(x=2011.8,y=0.21,'C',cex=1.8)
plot(Tindx,dat0$Dbar,type='l',lwd=2,col='mediumblue',
     ylab='Daph. length',xlab='')
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')
text(x=2008.2,y=1.4,'D',cex=1.8)
LG = adjust_transparency('limegreen',alpha=0.5)
plot(Tindx,dat0$chl,type='l',lwd=1,col='black',
     ylab='Chl mg/m^3',xlab='')
points(Tindx,dat0$chl,type='l',lwd=3,col=LG)
grid()
abline(v=c(2008:2012),lty=3,lwd=3,col='darkblue')
text(x=2008.2,y=12,'E',cex=1.8)
