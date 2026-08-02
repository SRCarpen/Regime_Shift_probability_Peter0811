# Read and spline MARSS-DLM results
# 27 Jan 2026

rm(list = ls()) 
graphics.off()

library(stats)
library(Boom)
#library(MARSS)
library(tictoc)
library(cubature)
library(numDeriv)
library(splines)
library(npreg)
library(tictoc)

source('Construct_Bspline_Multiyear_v1.R')
source('EPFunction+EQ_ssML_nKnot.R')  # use ss() from npreg

# FUNCTIONS ====================================================================

# function to regress-out nonstationarity
regstat = function(T,Z) {
  lm0 = lm(Z ~ T)
  err = lm0$residuals
  return(err)
}

# END FUNCTIONS ================================================================

# Save results for spline fit of Langevin
# Include a list of input data: xoriginal is untransformed, Xvar0 is sqrt(xoriginal),
#  mean and sd are mu.Xvar0 and sd.Xvar0, XZ is Z-score of Xvar0
# Xtrans.lst = list(xoriginal,Xvar0,mu.Xvar0,sd.Xvar0,XZ)
#save(dat0c,Xtrans.lst,                                    # input data & transform
#     mat1,xvar0,xvar1,Tindx,TT,DT,                         # inputs to MARSS
#     b0,b1,b0.se,b1.se,stdlevel,yhat,sm.y,sm.yse,   # outputs of MARSS
#     file=Fname)
# MARSS RESULTS ARE Z SCORES OF SQRT(CHL)
load(file="DLM-MARSS_result_Chl_PeterSqueal.Rdata")

# build data frame with MARSS-DLM result using mat1 as skeleton
mat1.original = mat1  # save the original mat1
# MARSS-DLM output vectors are same length as mat1, BUT we need to make dx to fit mat1
# Options:  Add a fake data point to output vectors OR delete one row of mat1.
# The MARSS-DLM output starts with a transient steep decline that could be deleted.
# Below I add one fake data point at the start of the output vector, then delete it later
#
# DLM outputs to consider:  b0, stdlevel, yhat, sm.y
Xdlm = c(0,sm.y)  # select DLM output and add a placeholder at the start
title = c('Spline of MARSS-DLM')
Fname = c('Spline_MARSSdlm_predix-tprobs_squeal1.Rdata')
ndlm = length(Xdlm)
Xdlm0 = Xdlm[1:(ndlm-1)]
Xdlm1 = Xdlm[2:ndlm]
# Plug the dlm result into mat1
mat1$X0 = Xdlm0
mat1$X1 = Xdlm1

# for stdlevel, the initial points larger than about 5 are outliers
# Remove the first 8 points:  NOT NECESSARY FOR sm.y
#mat2 = tail(mat1,-8)  # remove the first 8 points
# insert mat2 for mat1 for splining
#mat1 = mat2

# unpack mat1
year = mat1$year
doy = mat1$doy
x0 = mat1$X0
x1 = mat1$X1
dx = (x1-x0)/DT   # first moment
nx = length(x0)

uy = unique(year)

# make a doy variate with year and scaled doy within the year
doyrng = range(doy,na.rm=T)
Tindx = year + ((doy - doyrng[1] + 1)/(doyrng[2] - doyrng[1] + 1))

# plot X0 and diff(X0)
windows(height=8,width=5)
par(mfrow=c(2,1),mar=c(4, 4.2, 1, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tindx,x0,type='l',lwd=2,col='forestgreen',main=title)
abline(v=uy,lty=3,lwd=2)
plot(Tindx,dx,type='l',lwd=1,col='slateblue')
grid()
abline(h=0,lty=3,lwd=2,col='darkred')
abline(v=uy,lty=3,lwd=2)

# Number of knots for mu (deterministic core) and sigma (sqrt(conditional variance))
nk.mu = 5
nk.sig = 5
nk.ss = 7  # knots for EPF calculation via ss()

# Order of poynomial spline (usually 3 for cubic spline; 2 may dampen fluctuations)
npoly.mu = 3
npoly.sig = 3

# call Bspline.Langevin
# single year version:
#Bspline.Langevin = function(Xvar,nx,nk.mu,nk.sig,npoly.mu,npoly.sig)
# multiyear version
Bout = Bspline.Langevin(x0,x1,dx,nx,DT,nk.mu,nk.sig,npoly.mu,npoly.sig)

# output list 
#outlist=list(DF1,LSD1,LSsigma,xeq.D1,lmD1err,lmsigerr,
# basisD1,D1y,D1mat,D1D1,iD1D1,rdf.D1,
#basisig,sigmat,isig2)
DF1 = Bout[[1]]
LSD1 = Bout[[2]]
LSsigma = Bout[[3]]
xeq.D1 = Bout[[4]]
lmD1err = Bout[[5]]
lmsigerr = Bout[[6]]
basisD1 = Bout[[7]]  # recyclable spline basis
D1y = Bout[[8]]
D1mat = Bout[[9]]  # recyclable spline design matrix
D1D1 = Bout[[10]]  # squared spline design matrix
iD1D1 = Bout[[11]]  # recyclable inverse squared spline design matrix
rdf.D1 = Bout[[12]]  # residual d.f. for D1
basisig = Bout[[13]]
sigmat = Bout[[13]]
isig2 = Bout[[15]]
print(c('equilibria of D1: ',round(Bout[[4]],4)),quote=F)

windows(width=8,height=10)
par(mfrow=c(4,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(DF1$x0,LSD1,type='l',lwd=2,col='blue',
     xlab='x0',ylab='D1 by LS')
abline(h=0,lty=3,lwd=2)
plot(DF1$x0,lmD1err,type='p',pch=19,col='magenta',
     xlab='X0',ylab='D1 residual')
grid()
plot(DF1$x0,LSsigma,type='l',lwd=2,col='red',
     xlab='x0',ylab='sigma by LS')
plot(DF1$x0,lmsigerr,type='p',pch=19,col='black',
     xlab='X0',ylab='sigma residual')
grid()

# attempt to integrate EPF with EPFEQ
#EPFEQ = function(X,D1,sigma,nknots) 
#outlist = list(X.ep,EPF,dEPdx,xeq,D12)
EPout = EPFEQ(DF1$x0,LSD1,LSsigma,nk.ss) 
EPX = EPout[[1]]
EPF = EPout[[2]]
dEPdx = EPout[[3]]
xeq.epf = EPout[[4]]
NEPF = length(EPF)

windows(width=12,height=4)
par(mfrow=c(1,3),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.8,cex.lab=1.8)
plot(DF1$x0,LSD1,type='l',lwd=2,col='darkblue',
     xlim=c(-1.5,1),
     #ylim=c(-0.1,0.15),
     xlab='Z score sqrt(Chl)',ylab='Drift')
abline(h=0,lty=3,lwd=2,col='black')
grid()
plot(EPX[1:NEPF],EPF,type='l',lwd=2,col='darkblue',
     xlim=c(-1.5,1),ylim=c(0,4),
     xlab='Z score sqrt(Chl)',
     ylab='Effective Potential')
abline(v=xeq.epf,lty=3,lwd=2,col='darkred')
grid()
plot(EPX[1:length(dEPdx)],dEPdx,type='l',lwd=2,col='darkblue',
     xlim=c(-1.5,1),ylim=c(-3,3),
     xlab='Z score sqrt(Chl)',ylab='-Slope Effective Potential')
abline(v=xeq.epf,lty=3,lwd=2,col='darkred')
grid()
abline(h=0,lty=3,lwd=2,col='darkred')

print('equilibria from original data',quote=F)
print(c('deterministic, D1:  ',round(xeq.D1,4)),quote=F)
print(c('stochastic, EPF:  ',round(xeq.epf,4)),quote=F)

# MAKE SUMMARY SLIDEs ==========================================================
windows(height=8,width=8)
par(mfrow=c(2,1),mar=c(4, 4.2, 1, 2) + 0.1,cex.axis=1.7,cex.lab=1.7)
plot(EPX[1:NEPF],EPF,type='l',lwd=2,col='darkblue',
     xlim=c(-1.5,1),ylim=c(0,4),
     xlab='Z score sqrt(Chl)',
     ylab='Effective Potential')
abline(v=xeq.epf,lty=3,lwd=2,col='darkred')
grid()
plot(EPX[1:length(dEPdx)],dEPdx,type='l',lwd=2,col='darkblue',
     xlim=c(-1.5,1),ylim=c(-3,3),
     xlab='Z score sqrt(Chl)',ylab='-Slope Effective Potential')
abline(v=xeq.epf,lty=3,lwd=2,col='darkred')
grid()
abline(h=0,lty=3,lwd=2,col='darkred')

#
xoriginal = Xtrans.lst[[1]]
#plotmat = matrix(c(1,1,2,2))
#layout(plotmat)
windows(height=8,width=8)
par(mfrow=c(2,1),mar=c(2, 4.2, 2, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tindx,xoriginal[2:433],pch=20,cex=1,col='forestgreen',xlab='',
     ylab='Chlorophyll mg/m^3',main='Daily Chlorophyll During Summer')
grid()
abline(v=c(2008,2009,2010,2011,2012),lty=3,lwd=3,col='black')
par(mar=c(4, 4.2, 2, 2))
plot(Tindx,x0,type='l',lwd=2,col='blue',xlab='Year',
     ylab='Z score of sqrt(Chlorophyll)')
grid()
abline(v=c(2008,2009,2010,2011,2012),lty=3,lwd=3,col='black')
abline(h=0,lty=3,lwd=3,col='black')


# END OF SUMMARY SLIDES ==========================================================

# BACK-TRANSFORM X-AXES OF SPLINE AND EPF ***************************************
# Xtrans.lst = list(xoriginal,Xvar0,mu.Xvar0,sd.Xvar0,XZ)
mu.Xvar0 = Xtrans.lst[[3]]
sd.Xvar0 = Xtrans.lst[[4]]
# Z = (X - mu)/sigma implies Z*sigma + mu = X
x.spline = DF1$x0*sd.Xvar0 + mu.Xvar0  # units are sqrt(chl)
# from Carpenter & Pace 2018 appendix the back-transform of sqrt(x) is
#  sqrt(x)^2 + sd(sqrt(x))^2 in the case where sqrt(x) is predicted by regression.
# However our x-axes are observed not predicted.  We can just square x.spline.
x.chl.spl = x.spline^2
# the EPF x-axis is also an observation. 
x.EPX = EPX*sd.Xvar0 + mu.Xvar0  # reversal of z-scoring
x.chl.ep = x.EPX^2
# back-transform the equilibria
# reverse the z-transform
xeq.D1.z = xeq.D1*sd.Xvar0 + mu.Xvar0
xeq.EP.z = xeq.epf*sd.Xvar0 + mu.Xvar0
# reverse the sqrt
xeq.D1.chl = xeq.D1.z^2
xeq.EP.chl = xeq.EP.z^2
print('equilibria in Chl mg/m^3',quote=F)
print(c('deterministic, D1:  ',round(xeq.D1.chl,4)),quote=F)
print(c('stochastic, EPF:  ',round(xeq.EP.chl,4)),quote=F)

# plot D1, EPF, dEPF/dx versus Chl
# row format
windows(width=12,height=4)
par(mfrow=c(1,3),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.8,cex.lab=1.8)
plot(x.chl.spl,LSD1,type='l',lwd=2,col='darkblue',xlim=c(1.5,7),ylim=c(-0.1,0.2),
     xlab='Chlorophyll a, mg/m^3',ylab='Drift')
abline(h=0,lty=3,lwd=3,col='darkred')
grid()
text(x=2,y=0.18,'A',cex=1.8)
plot(x.chl.ep[1:NEPF],EPF,type='l',lwd=2,col='blue',xlim=c(1.5,7),ylim=c(0,6),
     xlab='Chlorophyll a, mg/m^3',
     ylab='Effective Potential')
grid()
abline(v=xeq.EP.chl,lty=3,lwd=2,col='darkred')
text(x=2,y=5.6,'B',cex=1.8)
plot(x.chl.ep[1:length(dEPdx)],dEPdx,type='l',lwd=2,col='purple',
     xlim=c(1.5,7),ylim=c(-6,10),
     xlab='Chlorophyll a, mg/m^3',ylab='- Slope Effective Potential',
    )
grid()
abline(h=0,lty=3,lwd=3,col='darkred')
abline(v=xeq.EP.chl,lty=3,lwd=2,col='darkred')
text(x=2,y=9,'C',cex=1.8)

# column format
windows(width=5,height=12)
par(mfrow=c(3,1),mar=c(1, 4.3, 2, 2) + 0.1,cex.axis=1.8,cex.lab=1.8)
plot(x.chl.spl,LSD1,type='l',lwd=2,col='darkblue',xlim=c(1.5,7),ylim=c(-0.1,0.15),
     #xlab='Chlorophyll a, mg/m^3',
     ylab='Drift')
abline(h=0,lty=3,lwd=3,col='darkred')
grid()
plot(x.chl.ep[1:NEPF],EPF,type='l',lwd=2,col='blue',xlim=c(1.5,7),ylim=c(0,6),
     #xlab='Chlorophyll a, mg/m^3',
     ylab='Effective Potential')
grid()
abline(v=xeq.EP.chl,lty=3,lwd=2,col='darkred')
par(mar=c(4, 4.3, 2, 2) + 0.1)
plot(x.chl.ep[1:length(dEPdx)],dEPdx,type='l',lwd=2,col='purple',
     xlim=c(1.5,7),ylim=c(-6,10),
     xlab='Chlorophyll a, mg/m^3',ylab='- Slope Effective Potential',
)
grid()
abline(h=0,lty=3,lwd=3,col='darkred')
abline(v=xeq.EP.chl,lty=3,lwd=2,col='darkred')

# STATISTICS CONTINUE IN THE TRANSFORMED UNITS
# calculate log(odds) ratio for low attractor
# We want to predict x1 = x0 + D1*DT which has error x1 - x0 - D1*DT
# residuals = DF1$x1 - DF1$x0 - (D1mat%*%wD1)*DT 
thresh = xeq.epf[2]  # threshold separating attractors
Tdist = mat1$X0 - thresh  # distance from threshold; positive above, negative below
modsig2D1 = sum(lmD1err^2)/length(lmD1err)    # prediction variance for dx
# use eq 1.4.9 of Draper & Smith
mux0 = mean(mat1$X0,na.rm=T)
devx0 = mat1$X0 - mux0
devx02 = devx0^2
ss.devx0 = sum(devx02)
# variances of dx predictions are variances of X1 predictions
vdxhat = modsig2D1 + ((devx02*modsig2D1)/ss.devx0) 
# variances of smoothed yhat from MARSS-DLM
varobs = sm.yse^2
print('',quote=F)
print('compare variances of spline predix and dlm-MARSS smoothed predix',quote=F)
print('spline',quote=F)
print(summary(vdxhat))
print('dlm-MARSS',quote=F)
print(summary(varobs))
print('',quote=F)
#
# t probabilities of lower tail
Tstat = Tdist/sqrt(vdxhat)  
Tprob = 1 - pt(Tstat,df=rdf.D1,lower.tail=T,log.p=F) # if Tstat < 0 Tprob > 0.5

# scatter plot of T result
windows()
par(mfrow=c(1,1),mar=c(4, 4.3, 0.5, 2) + 0.1,cex.axis=1.8,cex.lab=1.8)
plot(Tstat,Tprob,type='p',pch=20,cex=0.8,col='darkred',
     xlab='t statistic',ylab='probability below threshold')
grid()
abline(v=0,lty=3,lwd=2)
abline(h=0.5,lty=3,lwd=2)

# time plot of T result
windows(height=8,width=8)
par(mfrow=c(3,1),mar=c(1.3, 4.3, 0.7, 2) + 0.1,cex.axis=1.8,cex.lab=1.8)
plot(Tindx,Tdist,type='l',lwd=2,col='blue',#xlab='Year',
     ylab='Distance')
grid()
abline(v=uy,lty=3,lwd=2)
abline(h=0,lty=3,lwd=3,col='darkred')
text(x=2008.3,y=3.1,'A',cex=1.8)
plot(Tindx,Tstat,type='l',lwd=2,col='blue',#xlab='Year',
     ylab='Student t statistic')
grid()
abline(v=uy,lty=3,lwd=2)
abline(h=0,lty=3,lwd=3,col='darkred')
text(x=2008.3,y=10,'B',cex=1.8)
par(mar=c(4, 4.3, 0.7, 2) + 0.1)
plot(Tindx,Tprob,type='l',lwd=2,col='blue',xlab='Year',
     ylab='p(Low Chlorophyll)')
grid()
abline(v=uy,lty=3,lwd=2)
text(x=2008.3,y=0.92,'C',cex=1.8)

# add Tdist, Tdev, vdxhat, tprob to mat1
Stut = as.data.frame(cbind(Tdist,Tstat,vdxhat,Tprob))
colnames(Stut) = c('Threshdist','Tstat','varpred','Tprob')

# merge mat1 with Stut
matall = cbind(mat1,Stut)
# merge matall with dat0c, priority to year and doy in matall
dat4reg = merge(x=matall,y=dat0c,by=c('year','doy'))

# save final dat4reg plus input dataframes
save(dat0c,mat1,Stut,dat4reg,Tindx,file=Fname)
print(Fname,quote=F)
