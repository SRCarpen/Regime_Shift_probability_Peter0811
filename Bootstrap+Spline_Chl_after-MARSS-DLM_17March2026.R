# Nominal and bootstrap spline test
# SRC 2025-12-07

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

# Load data
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

# SET UP FOR BOOTSTRAP =========================================================

tstart = Sys.time() # start clock --------------------------------------------

#Set file name and title for output
Fname = c('Boot_10_Chl_Spline_Peter0811.Rdata')
title = c("Daily Chl")

# Bootstrap attempt 1 ==========================================================
# Remember that D2 is calculated from the residuals of D1
nboot = 10 # set number of samples
rboot = length(LSD1)
D1boot = matrix(0,nr=rboot,nc=nboot)  # save bootstraps of D1
sigboot = matrix(0,nr=rboot,nc=nboot)  # save bootstraps of sigma
booteqD1 = as.vector(rep(0,nboot))  # save bootstrapped number of eq of D1
EPXboot = matrix(0,nr=(NEPF+1),nc=nboot)
EPFboot = matrix(0,nr=NEPF,nc=nboot) # save bootstrapped EPF
booteqepf = as.vector(rep(0,nboot))  # save bootstrapped number of eq of EPF

# rationale:
# (1) The nominal D1 estimate LSD1 = D1mat%*%wD1
# (2) The spline for D1 predicts x1 = x0 + D1*DT which has error x1 - x0 - D1*DT
# (3) Therefore we calculate x1 pseudodata as
#       x1.pseudo =  DF1$x0 + (LSD1 + randomized(lmD1err))*DT
#           where the lmD1err values are randomized without replacement (all are used)
#       dx.pseudo = x1.pseudo - DF1$x0 
# (4) and fit 1 column of D1boot using DF1$x0 and dx.pseudo
# Then repeat for nboot cycles
# Because x0 stays the same the spline basis basisD1 
#   and the design matrix D1mat are unchanged and can be recycled
# 
nD1err = length(lmD1err)
for(i in 1:nboot) {      # start bootstrap cycle
  ranD1err = sample(lmD1err,size=nD1err,replace=F)
  ranx1 = DF1$x0 + (LSD1 + ranD1err)*DT
  randx = ranx1 - DF1$x0
  ranD1y = randx  
  ranwD1 = as.vector(iD1D1%*%t(D1mat)%*%ranD1y) # random weights for D1 spline
  D1boot[,i] = D1mat%*%ranwD1
  # find D1 equilibria
  sdrift = sign(D1boot[,i])
  dsdrift = c(0,-diff(sdrift))
  raneq = DF1$x0[which(!dsdrift == 0)]
  raneq.4 = round(raneq,4)
  #print(raneq)
  neq = length(raneq)
  booteqD1[i] = neq
  # find sigma equilibria
  eP = ranx1 - DF1$x0 - D1boot[,i]*DT  # which should match ranD1err!
  DF1$ep2 = 0.5*(eP^2)  # D2
  DF1$sep = sqrt(2*DF1$ep2)
  sigy = DF1$sep   # dependent variate for design matrix of sigma
  wsig = as.vector(isig2%*%t(sigmat)%*%sigy) # weights for sigma spline
  sigboot[,i] = sigmat%*%wsig 
  # attempt to integrate EPF with EPFEQ
  #EPFEQ = function(X,D1,sigma,nknots) 
  #outlist = list(X.ep,EPF,dEPdx,xeq,D12)
  epf.out = EPFEQ(DF1$x0,D1boot[,i],sigboot[,i],nk.ss)
  EPXboot[,i] = epf.out[[1]]  # this might be D1$x0 
  EPFboot[,i] = epf.out[[2]]
  eqepf = epf.out[[4]]
  eqepf.4 = round(eqepf,4)
  neq = length(eqepf)
  booteqepf[i] = neq
  print('',quote=F)
  print(c('complete cycles: ',i,' ----------------------------------------'),quote=F)
  print(c('D1 eq ',raneq.4),quote=F)
  print(c('EPF eq ',eqepf.4),quote=F)
}

# -----------------------------------------------------------------------------
tstop = Sys.time()
print('----------------------------------------------------------',quote=F)
runtime = difftime(tstop,tstart,units='mins')
print(c('runtime for bootstrapping, minutes ',runtime),quote=F)
print('----------------------------------------------------------',quote=F)
print('',quote=F)

print('frequencies of numbers of deterministic D1 equilibria',quote=F)
freq.table = table(booteqD1)
print('frequency table for number of equilibria',quote=F)
print('first row is number of equilibria ',quote=F)
print(c('second row is number of bootstrap cycles of total ',nboot),quote=F)
print(freq.table)
print(c('median number of equilibria = ',median(booteqD1)),quote=F)

print('frequencies of numbers of stochastic EPF equilibria',quote=F)
freq.table = table(booteqepf)
print('frequency table for number of equilibria',quote=F)
print('first row is number of equilibria ',quote=F)
print(c('second row is number of bootstrap cycles of total ',nboot),quote=F)
print(freq.table)
print(c('median number of equilibria = ',median(booteqepf)),quote=F)

D1range = range(c(LSD1,D1boot),na.rm=T)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5) 
plot(DF1$x0,LSD1,type='l',ylim=D1range,lwd=3,col='darkblue',xlab='X0',ylab='D1')
for(i in 1:nboot) {
  points(DF1$x0,D1boot[,i],type='l',lwd=1,col='skyblue2')
}
abline(h=0,lty=3,lwd=3,col='darkred')
grid()

sigrange = range(LSsigma,sigboot)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5) 
plot(DF1$x0,LSsigma,type='l',ylim=sigrange,lwd=3,col='darkred',xlab='X0',ylab='sigma')
for(i in 1:nboot) {
  points(DF1$x0,sigboot[,i],type='l',lwd=1,col='hotpink')
}
grid()

EPFrange = range(EPF,EPFboot,rm=T)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(EPX[1:NEPF],EPF,ylim=EPFrange,log='y',type='l',lwd=3,col='darkgreen',xlab='X',ylab='EPF')
for(i in 1:nboot)  {
  points(EPXboot[1:NEPF,i],EPFboot[,i],type='l',lwd=1,col='limegreen')
}
grid()
# nboot is number of bootstrap cycles, DF1 is data frame of data to bootstrap,
#  DF1$x0 is the X vector for plots of Langevin bootstraps,
#  D1boot is D1 bootstraps, booteqD1 is number of equilibria of D1 bootstraps, 
#  sigboot is sigma bootstraps, EPXboot is X axes of bootstrapped EPFboot,
#  booteqepf is number of equilibria of EPF bootstraps.
#  Nominal estimates from data:  LSD1, LSsigma, EPX, EPF, xeq.D1, xeq.epf
save(nboot,runtime,DF1,D1boot,booteqD1,sigboot,EPXboot,EPFboot,booteqepf,NEPF,
     LSD1,LSsigma,xeq.D1,EPX,EPF,xeq.epf,file=Fname)
print(c('results saved in ',Fname),quote=F)
