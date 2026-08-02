# Final LMB model 2026-07-19

rm(list = ls())
graphics.off()

# Source of the bootstrapped error of the unstable equilibrium
# Peter_bootstrap_v0residuals/After-Boot_Chl-R08-11_using_DOsat_afterboot_2026-07-09.R
sd.xeq2 = 0.27929915
# source of the splined Chl and its prediction error:
# Read+Spline_Chl_MARSS-DLM_result_15March2026.R
# t probabilities of lower tail where vdxhat is varpred, the prediction variance
# Tstat = Tdist/sqrt(vdxhat)  
# Tprob = 1 - pt(Tstat,df=rdf.D1,lower.tail=T,log.p=F) # if Tstat < 0 Tprob > 0.5
rdf.D1 = 426  # spline's residual d.f. needed for Student t distribution
#
# source of the merged spline result and food web data:
# Combine_Raft+Zoop+Fish_08-11_squeal1_v1_2026-02-10.R
# load data from Combine_Raft+Zoop+Fish_08-11_squeal1_v1_2026-02-10.R
#save(all4reg,file='CombineCovarSplineSqueal1.Rdata')
load(file='CombineCovarSplineSqueal1.Rdata')

print(all4reg[1,])
dat0 = all4reg

# extract some variates
Tindx = dat0$Tindx  # year index
Tprob = dat0$Tprob  # student-t probabilities for nominal case
Tdist = dat0$Threshdist
Tscale = sqrt(dat0$varpred)
# T statistic is Tdist/Tscale
#zmix = dat0$zmix
#medtemp = dat0$medtemp
#dosat50 = dat0$dosat50
#Chl = dat0$chl
#ZB = dat0$ZB
#Dbar = dat0$Dbar
#pelratio = dat0$pelratio
#litmin = dat0$litmin
#totmin = dat0$totminnow
lmb = dat0$lmb

# source of the bootstrapped standard deviation pf EPF equilibria 
#   of the unstable equilibrium of chl in final units = 0.27929915
# After-Boot_Chl-R08-11_using_DOsat_afterboot_2026-07-09.R

# source of the model:
# ORR_TransformFoodChain_V2spline_LogOdds_19Feb2026.R
# save(fwmod5a,file='LMB_model_final.Rdata')
load(file='LMB_model_final.Rdata')

print(summary(fwmod5a))

# generate predix and predix s.e. for fwmod5a
p5a = predict(fwmod5a,se.fit=T,interval='none')
yhat = p5a$fit
#yse = p5a$se.fit  # regression s.e. of prediction
# yse = sd.xeq2  # s.d. of xeq2
yse = sqrt( (p5a$se.fit)^2 + sd.xeq2^2)  # sqrt sum of variances
predmat = matrix(0,nr=length(yhat),nc=3)
predmat[,1] = yhat - yse
predmat[,2] = yhat
predmat[,3] = yhat + yse
epred = exp(predmat)
ppred = epred/(1 + epred)
  
# extract predix of fwmod5a and convert them to prob(low chl state)
yhat = fwmod5a$fitted.values
eyhat = exp(yhat)
plowhat = eyhat/(1 + eyhat)  # prob(low chl) predicted by model 5a

# MAKE A PLOT WITH lmb, model fwmod5a
# calculate predicted probs on fish gradient using lm5
#X5 = unname(fwmod5a$X)
b5 = unname(fwmod5a$coefficients)
# lmb gradient and log10(lmbgradient)
rng.lmb = range(lmb,na.rm=T)
#lmbgrad = seq(rng.lmb[1],rng.lmb[2],length.out=100)
lmbgrad = seq(25,175,length.out=100)  # beyond the range of the data
#loglmbgrad = log10(lmbgrad)
x0 = rep(1,100)
xgrad = matrix(c(x0,lmbgrad),nr=100,nc=2,byrow=F)
# convert log(odds ratio) for low chl to prob(low chl)
lods = xgrad%*%b5  # prediction is log(odds ratio)
elods = exp(lods)  # odds ratio
phat0 = elods/(1 + elods)  # probability of low chl state
print('phat0 is p(low chl) for predix at mean threshold',quote=F)
print('phat1 will be p(low chl) for predix 1 sigma below threshold',quote=F)
# Tstat = Tdist/sqrt(vdxhat)  
# Tprob = 1 - pt(Tstat,df=rdf.D1,lower.tail=T,log.p=F) # if Tstat < 0 Tprob > 0.5
Tstat1 = (Tdist - sd.xeq2)/Tscale
Tprob1 = 1 - pt(Tstat1,df=rdf.D1,lower.tail=T,log.p=F) # if Tstat < 0 Tprob > 0.5
p.below = Tprob1
p.above = 1-Tprob1
lodrat1 = log( (p.below + 0.001)/(p.above + 0.001))  # add a small constant to avoid 0
mod5a1 = lm(lodrat1 ~ lmb, x=T)
print(summary(mod5a1))
b5a1 = mod5a1$coefficients
lods = xgrad%*%b5a1  # prediction is log(odds ratio)
elods = exp(lods)  # odds ratio
phat1 = elods/(1 + elods)  # probability of low chl state
print('',quote=F)
print('phat2 will be p(low chl) for predix 1 sigma above threshold',quote=F)
Tstat2 = (Tdist + sd.xeq2)/Tscale
Tprob2 = 1 - pt(Tstat2,df=rdf.D1,lower.tail=T,log.p=F) # if Tstat < 0 Tprob > 0.5
p.below = Tprob2
p.above = 1-Tprob2
lodrat2 = log( (p.below + 0.001)/(p.above + 0.001))  # add a small constant to avoid 0
mod5a2 = lm(lodrat2 ~ lmb, x=T)
print(summary(mod5a2))
b5a2 = mod5a2$coefficients
lods = xgrad%*%b5a2  # prediction is log(odds ratio)
elods = exp(lods)  # odds ratio
phat2 = elods/(1 + elods)  # probability of low chl state

phatrange = range(c(phat0,phat1,phat2))
windows()
par(mfrow=c(1,1),mar=c(4.5,4.5,2,1)+0.1,cex.lab=1.6,cex.axis=1.6)
plot(lmbgrad,phat0,ylim=phatrange,type='l',lwd=2,col='darkblue',
     xlab='Piscivorous Largemouth Bass',
     ylab='Probability of Low-Chlorophyll State',
     main='Peter Lake 2008-2011')
points(lmbgrad,phat1,type='l',lwd=2,col='deepskyblue')
points(lmbgrad,phat2,type='l',lwd=2,col='deepskyblue')
grid(col='gray')

years = c(2008:2012)  # integer years

windows(height=5,width=9)
par(mfrow=c(1,2),mar=c(4.5,4.3,2,1)+0.1,cex.lab=1.7,cex.axis=1.6)
plot(lmbgrad,phat0,ylim=phatrange,type='l',lwd=2,col='darkblue',
     xlab='Piscivorous Largemouth Bass',
     ylab='P(low Chl state)',
     #main='Peter Lake 2008-2011'
     )
points(lmbgrad,phat1,type='l',lwd=2,col='deepskyblue')
points(lmbgrad,phat2,type='l',lwd=2,col='deepskyblue')
grid(col='gray')
text(x=38,y=0.95,'A',cex=1.7)
#
plim = range(ppred,na.rm=T)
plot(Tindx,ppred[,2],ylim=c(0,1),type='l',lwd=2,col='darkblue',xlab='Year',
     ylab='P(low Chl state)')
points(Tindx,ppred[,1],type='l',lwd=2,col='deepskyblue')
points(Tindx,ppred[,3],type='l',lwd=2,col='deepskyblue')
grid(col='gray')
abline(v=years,lty=3,lwd=3,col='darkred')
text(x=2008.5,y=0.95,'B',cex=1.7)