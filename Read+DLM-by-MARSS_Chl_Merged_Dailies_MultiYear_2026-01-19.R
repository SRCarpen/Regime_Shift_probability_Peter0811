# Read merged daily data Squeal 1 
# from Merge_daily_FoodWeb+sonde_Peter_squeal1.R
# MARSS the DLM for daily chl, and save results for spline fit to Langevin

rm(list = ls()) 
graphics.off()

library(stats)
library(Boom)
library(MARSS)
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

# read data
#save(datyd,FW,R.all.daily,file='Merged_sonde+foodweb_R_squeal1.Rdata')
load(file='Merged_sonde+foodweb_R_squeal1.Rdata')

# select data for analysis
dat0all = R.all.daily   # includes NA
dat0c = na.omit(dat0all)  #

#Set file name and title for output
Fname = c('DLM-MARSS_result_Chl_PeterSqueal.Rdata')
title = c("Daily Chl")

# Select data
xoriginal = dat0c$chl
#Xvar0 = dat0c$chl
#Xvar0 = log10(dat0c$chl)  # log transform if needed
Xvar0 = sqrt(dat0c$chl)  # sqrt transform if needed
#Xvar0 = 1/sqrt(dat0c$chl)  # reciprocal sqrt transform if needed
#Xvar0 = 1/dat0c$chl # reciprocal transform if needed
#
doy0 = dat0c$doy
year0 = dat0c$year
N0 = length(Xvar0)
# define DT
DT = 1

# Normal Probability plot -------------------------------------------
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
qqnorm(Xvar0,pch=20,col='red',xlab='Quantiles, sd',ylab='sqrt(Chl)',
       main='normal probability plot')
qqline(Xvar0,lwd=2,col='blue')
grid()
# --------------------------------------------------------------------

# Z-score for Xvar0
mu.Xvar0 = mean(Xvar0,na.rm=T)
sd.Xvar0 = sd(Xvar0,na.rm=T)
XZ = (Xvar0 - mu.Xvar0)/sd.Xvar0

# detrend before analysis?
#T0 = seq(1,N0,by=DT)
#resids = regstat(T0,Xvar0)
#Xvar0 = resids
# COMMENT OUT THIS BLOCK TO SKIP DETRENDING

# save original and transformed data for export
Xtrans.lst = list(xoriginal,Xvar0,mu.Xvar0,sd.Xvar0,XZ)

# CHOOSE Xvar0 OR XZ FOR DLM BY MARSS
# WHATEVER IS SELECTED WILL BE NAMED Xvar0 in further calculations
dat1 = as.data.frame(cbind(year0,doy0,XZ))
cnames = c('year0','doy0','Xvar0')
colnames(dat1)=cnames

# build dataset with x0, x1, dx
uy = unique(dat1$year0) 
nuy = length(uy)
NX = N0 - nuy  # allow for time step
mat0 = matrix(0,nrow=NX,ncol=4) # matrix for year0,doy0,x0,x1)
mat1 = as.data.frame(mat0)
cnames = c('year','doy','X0','X1')
colnames(mat1) = cnames
print(c('dim(mat1) = ',dim(mat1)),quote=F)
print(mat1[1,])

ncount =  0  # start counter
#
for(i in 1:nuy) {   
  dat1a = subset(dat1,subset=(dat1$year0 == uy[i]))
  dimd1a = dim(dat1a)
  nx1a = dimd1a[1]
  print(c('year = ',uy[i],', dim(dat1a) = ',dimd1a,', ncount = ',ncount),quote=F)
  mat1[(ncount+1):(ncount+nx1a-1),1:3] = dat1a[1:(nx1a - 1),1:3]
  mat1[(ncount+1):(ncount+nx1a-1),4] = dat1a[2:nx1a,3]
  ncount = ncount + nx1a - 1   # update ncount accounting for shift of x0 to x1
}

print(c('dim(mat1) = ',dim(mat1)),quote=F)
print(mat1[1,])

# unpack mat1
year = mat1$year
doy = mat1$doy
xvar0 = mat1$X0
xvar1 = mat1$X1
dx = (xvar1-xvar0)/DT   # first moment
nx = length(xvar0)

# make a doy variate with year and scaled doy within the year
doyrng = range(doy,na.rm=T)
Tindx = year + ((doy - doyrng[1] + 1)/(doyrng[2] - doyrng[1] + 1))

# plot x0 and dx
windows(height=8,width=5)
par(mfrow=c(2,1),mar=c(4, 4.2, 1, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(Tindx,xvar0,type='l',lwd=2,col='forestgreen',main=title)
abline(v=uy,lty=3,lwd=2)
plot(Tindx,dx,type='l',lwd=1,col='slateblue')
grid()
abline(h=0,lty=3,lwd=2,col='darkred')
abline(v=uy,lty=3,lwd=2)

# START CLOCK FOR MARSS
print('',quote=F)
print('Starting DLM by MARSS ================================================= ',quote=F)
tstart = Sys.time() # start clock --------------------------------------------

# set up for DLM by MARSS  
# IN THIS MODEL THE PREDICTOR IS SMOOTHED PREDICTION
# IN THE BAYES METHOD THE PREDICTOR IS THE PREVIOUS TIME STEP
DOY = Tindx
#yfull = matrix(XZ,nr=1,nc=nx) # y vector matching length DOY
#Nall = length(yfull)
#y0 = yfull[1:(Nall-1)]
#y1 = yfull[2:Nall]
ydat = matrix(t(xvar1),nr=1)  # response
xdat = matrix(t(xvar0),nr=1)  # predictor
TT = length(xdat)
x0 = matrix(rep(1,TT),nr=1)  # intercept dummy variate
m = 2 # number of regressors, xdat & x0

# for process eqn
B.mod = diag(m) # 2x2; Identity
U.mod = matrix(0,nrow=m,ncol=1) # 2x1; both elements = 0
Q.mod = matrix(list(0),m,m) # 2x2; all 0 for now
diag(Q.mod) = c("q1","q2") # 2x2; diag = (q1,q2)

# for observation eqn
Z.mod = array(NA, c(1,m,TT)) # NxMxT; empty for now
Z.mod[1,1,] = x0 #rep(1,TT) # Nx1; 1's for intercept dummy x0
Z.mod[1,2,] = xdat # Nx1; predictor variable
A.mod = matrix(0) # 1x1; scalar = 0
R.mod = matrix("r") # 1x1; scalar = r

# only need starting values for regr parameters
inits.list = list(x0=matrix(c(0, 0), nrow=m))
# list of model matrices & vectors
mod.list = list(B=B.mod, U=U.mod, Q=Q.mod, Z=Z.mod, A=A.mod, R=R.mod)

# fit univariate DLM
dlm1 = MARSS(ydat, inits=inits.list, model=mod.list)
# parameter series is in dlm1$states
# parameter s.e. are in dlm1$states.se

print('',quote=F)
print('Parameter estimates & std errors',quote=F)
parstats = MARSSparamCIs(dlm1)
print(parstats)


# extract intercept and predictor coef
b0 = dlm1$states[1,]
b1 = dlm1$states[2,]
b0.se = dlm1$states.se[1,]
b1.se = dlm1$states.se[2,]

# -----------------------------------------------------------------------------
tstop = Sys.time()
print('',quote=F)
print('----------------------------------------------------------',quote=F)
runtime = difftime(tstop,tstart,units='mins')
print(c('runtime for DLM, minutes ',runtime),quote=F)

windows(height=10,width=8)
par(mfrow=c(4,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(DOY[1:TT],b0,type='l',lwd=2,col='blue')
plot(DOY[1:TT],b0.se,type='l',lwd=2,col='blue')
plot(DOY[1:TT],b1,type='l',lwd=2,col='blue')
plot(DOY[1:TT],b1.se,type='l',lwd=2,col='blue')

# standardized level
stdlevel = b0/b0.se

windows(height=10,width=8)
par(mfrow=c(3,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(DOY[1:TT],b0,type='l',lwd=2,col='blue')
plot(DOY[1:TT],b0.se,type='l',lwd=2,col='blue')
plot(DOY[1:TT],stdlevel,type='l',lwd=2,col='blue')

windows()
par(mfrow=c(1,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(DOY[1:TT],ydat,type='l',lwd=2,col='seagreen',xlab='DOY',
     ylab='Chlorophyll, Z-scored',
     main='Chl Z-score green & b0 blue')
points(DOY[1:TT],b0,type='l',lwd=2,col='blue')
grid()

# model 1 step predix
yhat = dlm1$ytT
yhat.se = dlm1$ytT.se  # these are zeroes

# smoothed one-step predix
smdlm1 = tsSmooth(dlm1,type='ytt1')
sm.y = smdlm1$.estimate
sm.yse = smdlm1$.se  
# we can also get smoothed states, contemporaneous (xtt) or one-step-ahead (xtt1)
# smoothed results conditioned on all data are types ytT with states xtT

windows(height=10,width=8)
par(mfrow=c(3,1),mar=c(4,4,2,2)+0.1,cex.lab=1.5,cex.axis=1.5)
plot(DOY[1:TT],yhat,type='l',lwd=2,col='blue',ylab ='one-step pred. Chl',
     main='Chl Z score green, yhat blue')
points(DOY[1:TT],ydat,type='p',pch=19,col='forestgreen')
grid()
plot(DOY[1:TT],sm.y,type='l',lwd=2,col='blue',ylab='smoothed pred. Chl')
points(DOY[1:TT],ydat,type='p',pch=19,lwd=2,col='forestgreen')
grid()
plot(DOY[1:TT],sm.yse,type='l',lwd=2,col='blue',ylab='sm. pred. Chl s.e.')
grid()

windows(height=6,width=10)
par(mfrow=c(1,1),mar=c(4,4.2,2,2)+0.1,cex.lab=1.8,cex.axis=1.8)
plot(DOY[1:TT],sm.y,type='l',lwd=2,col='blue',ylab='Z score sqrt(Chlorophyll)',
     xlab='Year')
points(DOY[1:TT],ydat,type='p',pch=19,lwd=2,col='forestgreen')
grid()
legend('topright',legend=c('Smoothed Prediction','Observation'),pch=c(NA,19),
       lwd=c(2,NA),col=c('blue','forestgreen'),cex=1.5,bty='n')

# Density plots
dens.ydat = density(ydat,bw='SJ',window="epanechnikov",n=512,na.rm='T')
dens.lev = density(b0,bw='SJ',window="epanechnikov",n=512,na.rm='T')
dens.slev = density(stdlevel,bw='SJ',window="epanechnikov",n=512,na.rm='T')

windows(width=12,height=8)
par(mfrow=c(2,3),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.8)
plot(dens.ydat$x,dens.ydat$y,type='l',lwd=2,col='forestgreen',xlab='Z score of Observation',
     ylab='density')
plot(dens.lev$x,dens.lev$y,type='l',lwd=2,col='forestgreen',xlab='Level',
     ylab='density')
plot(dens.slev$x,dens.slev$y,type='l',lwd=2,col='forestgreen',xlab='Standardized Level',
     ylab='density')

# Density plots for predictions
dens.yhat = density(yhat,bw='SJ',window="epanechnikov",n=512,na.rm='T')
dens.smy = density(sm.y,bw='SJ',window="epanechnikov",n=512,na.rm='T')
# standardized smoothed y
stdysm = sm.y/sm.yse
dens.sysm = density(stdysm,bw='SJ',window="epanechnikov",n=512,na.rm='T')

#windows(width=12,height=4)
#par(mfrow=c(1,3),mar=c(4, 4.2, 3, 2) + 0.1,cex.axis=1.6,cex.lab=1.6)
plot(dens.yhat$x,dens.yhat$y,type='l',lwd=2,col='forestgreen',xlab='One-step Yhat',
     ylab='density')
plot(dens.smy$x,dens.smy$y,type='l',lwd=2,col='forestgreen',xlab='Smoothed Yhat',
     ylab='density')
plot(dens.sysm$x,dens.sysm$y,type='l',lwd=2,col='forestgreen',
     xlab='Standardized Smoothed Yhat',
     ylab='density')

# construct D1 and sigma from b0; use s.e. of b0 for parametric bootstrap
# Make data set with first and second moments
# 
# CHOOSE VARIATE TO SPLINE: all of the following have Z-scored Chl units;
#  XZ is Z-scored original data; b0 is DLM level; sm.y is smoothed y from DLM
# stdlevel from spline is also possible
#  Variates with estimated s.e. are b0 and sm.y
#x0 = sm.y[1:(TT-1)]
#x1 = sm.y[2:TT]

# Save results for spline fit of Langevin
# Include a list of input data: xoriginal is untransformed, Xvar0 is sqrt(xoriginal),
#  mean and sd are mu.Xvar0 and sd.Xvar0, XZ is Z-score of Xvar0
# Xtrans.lst = list(xoriginal,Xvar0,mu.Xvar0,sd.Xvar0,XZ)
save(dat0c,Xtrans.lst,                                    # input data & transform
     mat1,xvar0,xvar1,Tindx,TT,DT,                         # inputs to MARSS
     b0,b1,b0.se,b1.se,stdlevel,yhat,sm.y,sm.yse,   # outputs of MARSS
     file=Fname)
print(Fname)
