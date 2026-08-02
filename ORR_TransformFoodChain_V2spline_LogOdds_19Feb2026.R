# logistic regression with more complete fish data & daphnia length
# SRC 16 Feb 2026

rm(list = ls())
graphics.off()

library('stats')
library('car')
library('olsrr')

# load data from Combine_Raft+Zoop+Fish_08-11_squeal1_v1_2026-02-10.R
print('saving all covariates + t-stats for Squeal 1 logistic regression',quote=F)
#save(all4reg,file='CombineCovarSplineSqueal1.Rdata')
load(file='CombineCovarSplineSqueal1.Rdata')

print(all4reg[1,])

dat0 = all4reg

# extract some variates
Tindx = dat0$Tindx
Tprob = dat0$Tprob 
Tstat = dat0$Tstat
zmix = dat0$zmix
medtemp = dat0$medtemp
dosat50 = dat0$dosat50
Chl = dat0$chl
ZB = dat0$ZB
Dbar = dat0$Dbar
pelratio = dat0$pelratio
litmin = dat0$litmin
totmin = dat0$totminnow
lmb = dat0$lmb
lmb.c = dat0$lmb - mean(dat0$lmb,na.rm=T)
#zmix.c = dat0$zmix - mean(dat0$zmix,na.rm=T)

# log odds ratio for the low attractor, below the unstable point
#lodrat = log(p.below/(1-p.below))  # we get infinity if pbelow = 1
p.below = Tprob
p.above = 1-Tprob
lodrat = log( (p.below + 0.001)/(p.above + 0.001))  # add a small constant to avoid 0
dat0$lodrat = lodrat  # add lodrat to the dataframe

print('summary of log odds ratio',quote=F)
print(summary(lodrat))

windows()
par(mfrow=c(1,1),mar=c(4.5,4.5,2,1)+0.1,cex.lab=1.6,cex.axis=1.6)
plot(lmb,lodrat,type='p',pch=20,cex=0.8,col='blue',
     xlab='LMB',
     ylab='Log(odds ratio) of Low-Chlorophyll State',
     main='Peter Lake 2008-2011')
grid(col='gray')

# Check steps of the food chain to lodrat =======================================
print('',quote=F)
print('steps of the food chain to lodrat =======================================',quote=F)
print('',quote=F)
fwmod0 = lm(Chl ~ Dbar*ZB)
print(summary(fwmod0))
print('',quote=F)
fwmod1a = lm(lodrat ~ Dbar)
print(summary(fwmod1a))
print('',quote=F)
fwmod1b = lm(lodrat ~ ZB)
print(summary(fwmod1b))
print('',quote=F)
fwmod1c = lm(lodrat ~ ZB*Dbar)
print(summary(fwmod1c))
print('',quote=F)
fwmod2a = lm(Dbar ~ pelratio*litmin)
print(summary(fwmod2a))
print('',quote=F)
fwmod2b = lm(ZB ~ pelratio*litmin)
print(summary(fwmod2b))
print('',quote=F)
fwmod2c = lm(Dbar ~ totmin)
print(summary(fwmod2a))
print('',quote=F)
fwmod2d = lm(ZB ~ totmin)
print(summary(fwmod2b))
print('',quote=F)
fwmod3a = lm(pelratio ~ lmb)
print(summary(fwmod3a))
print('',quote=F)
fwmod3b = lm(litmin ~ lmb)
print(summary(fwmod3b))
print('',quote=F)
fwmod3b = lm(totmin ~ lmb)
print(summary(fwmod3b))
print('',quote=F)
fwmod4a = lm(Dbar ~ lmb)
print(summary(fwmod4a))
print('',quote=F)
fwmod4b = lm(ZB ~ lmb)
print(summary(fwmod4b))
print('',quote=F)
fwmod5a = lm(lodrat ~ lmb, x=T, y=T)
print(summary(fwmod5a))
print('',quote=F)
loglmb = log10(lmb)
fwmod5b = lm(lodrat ~ loglmb, x=T)
print(summary(fwmod5b))
print('End of food chain regressions ============================================',quote=F)

# ===============================================================================

save(fwmod5a,file='LMB_model_final.Rdata')

# convert lodrat to probability of low attractor
plow = 1/(1 + exp(lodrat))

# MAKE A PLOT WITH lmb, model fwmod5a
# calculate predicted probs on fish gradient using lm5
X5 = unname(fwmod5a$X)
b5 = unname(fwmod5a$coefficients)
# lmb gradient and log10(lmbgradient)
rng.lmb = range(lmb,na.rm=T)
lmbgrad = seq(rng.lmb[1],rng.lmb[2],length.out=100)
#loglmbgrad = log10(lmbgrad)
x0 = rep(1,100)
xgrad = matrix(c(x0,lmbgrad),nr=100,nc=2,byrow=F)
# convert log(odds ratio) for low chl to prob(low chl)
lods = xgrad%*%b5  # prediction is log(odds ratio)
elods = exp(lods)  # odds ratio
phat = elods/(1 + elods)  # probability of low chl state

windows()
par(mfrow=c(1,1),mar=c(4.5,4.5,2,1)+0.1,cex.lab=1.6,cex.axis=1.6)
plot(lmbgrad,phat,type='l',lwd=2,col='blue',
     xlab='Piscivorous Largemouth Bass',
     ylab='Probability of Low-Chlorophyll State',
     main='model 5a, Peter Lake 2008-2011')
grid(col='gray')

windows()
par(mfrow=c(1,1),mar=c(4.5,4.5,2,1)+0.1,cex.lab=1.6,cex.axis=1.6)
plot(lmbgrad,phat,type='l',lwd=2,col='blue',
     xlab='Piscivorous Largemouth Bass',
     ylab='Probability of Low-Chlorophyll State',
     main='Peter Lake 2008-2011')
grid(col='gray')

# calculate probability of lower attractor from lmb data at each time step
plow.lmb = splinefun(x=lmbgrad,y=phat,method='fmm')
plow.hat = plow.lmb(lmb)

years = c(2008:2012)  # integer years
windows()
par(mfrow=c(1,1),mar=c(4.5,4.5,2,1)+0.1,cex.lab=1.6,cex.axis=1.6)
plot(Tindx,plow.hat,type='l',lwd=2,col='blue',xlab='Year',ylab='P(low Chl state')
grid()
abline(v=years,lty=3,lwd=3,col='darkred')

windows(height=9,width=6)
par(mfrow=c(3,1),mar=c(4.5,4.5,2,1)+0.1,cex.lab=1.6,cex.axis=1.6)
plot(lmbgrad,phat,type='l',lwd=2,col='blue',
     xlab='Piscivorous LMB',
     ylab='P(low Chl state')
grid()
#
plot(Tindx,lmb,type='l',lwd=2,col='blue',xlab='Year',ylab='LMB population')
grid()
abline(v=years,lty=3,lwd=3,col='darkred')
#
plot(Tindx,plow.hat,type='l',lwd=2,col='blue',xlab='Year',ylab='P(low Chl state')
grid()
abline(v=years,lty=3,lwd=3,col='darkred')

windows(height=5,width=9)
par(mfrow=c(1,2),mar=c(4.5,4.3,2,1)+0.1,cex.lab=1.8,cex.axis=1.8)
plot(lmbgrad,phat,type='l',lwd=2,col='blue',
     xlab='Piscivorous LMB',
     ylab='P(low Chl state)')
grid()
text(x=38,y=0.85,'A',cex=1.8)
#
plot(Tindx,plow.hat,type='l',lwd=2,col='blue',xlab='Year',
     ylab='P(low Chl state)')
grid()
abline(v=years,lty=3,lwd=3,col='darkred')
text(x=2008.5,y=0.85,'B',cex=1.8)
