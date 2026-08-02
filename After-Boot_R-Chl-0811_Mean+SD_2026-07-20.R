# analyses of bootstraps
# SRC 2026-07-20

rm(list = ls())
graphics.off()

library(stats)
library(cubature)
library(numDeriv)
library(splines)
library(npreg)

# Function for equilibria of EPF ---------------------------------------------------
epfeq=function(epx,epf,nknots)  {
  # convert epf to a function
  xEPqt = quantile(epx,probs=seq(0,1,length.out=nknots),na.rm=T) # quantiles of X used for knots
  EPspline = ss(x=epx,y=epf,method='ML',m=2,knots=xEPqt)
  EPfun = function(x)  {
    yhat = predict(EPspline,x)$y
    return(yhat)
  }
  # take first derivative and find the roots
  dEPdx = grad(EPfun,epx,'Richardson')
  # Potential is a negative integral; reverse sign of derivative for plots
  # find roots of first derivative of effective potential
  sdrift = sign(dEPdx)
  dsdrift = c(0,-diff(sdrift))
  xeq = epx[which(!dsdrift == 0)]
  return(xeq)
}
# --------------------------------------------------------------------------------

# load bootstrap result
# nboot is number of bootstrap cycles, DF1 is data frame of data to bootstrap,
#  DF1$x0 is the X vector for plots of Langevin bootstraps,
#  D1boot is D1 bootstraps, booteqD1 is number of equilibria of D1 bootstraps, 
#  sigboot is sigma bootstraps, EPXboot is X axes of bootstrapped EPFboot,
#  booteqepf is number of equilibria of EPF bootstraps.
#  Nominal estimates from data:  LSD1, LSsigma, EPX, EPF, xeq.D1, xeq.epf
#save(nboot,DF1,D1boot,booteqD1,sigboot,EPXboot,EPFboot,booteqepf,NEPF,
#     LSD1,LSsigma,xeq.D1,EPX,EPF,xeq.epf,file=Fname)
load(file='Boot_1000_Chl_Spline_Peter0811.Rdata')
# file for results:
Fname = c('Density-mean-sd_3eq_Peter0811.Rdata')
# Analysis & plots from the bootstrap program

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

# Extract cases of 3 equilibria
# density(X.dlm,bw='SJ',window="epanechnikov",n=512,na.rm='T')
# epfeq=function(epx,epf,nknots)  # use the inline function - it's faster
# return(xeq) where xeq is equilibria of EPF
print('',quote=F)
print('Extract cases of 3 equilibria',quote=F)
tstart = Sys.time() # start clock --------------------------------------------
print(freq.table[4])
N3 = unname(freq.table[4])
print(c('Bootstraps with 3 equilibria found by the bootstrap program ',N3),quote=F)
# if there were 3 EPF equilibria then find them
EQ3mat = matrix(0,nr=nboot,nc=3)  # matrix to hold cases of 3 EPF equilibria
Ncase = 0  # counter for cases of 3 equilibria
i3eq = c(0)  # save row numbers with 3 equilibia
for(i in 1:nboot) {
  if(booteqepf[i] != 3) {next}  # if there are not 3 equilibria then skip
  else{
    # epfeq=function(epx,epf,nknots)
    xeqboot = epfeq(EPXboot[(1:NEPF),i],EPFboot[(1:NEPF),i],7) 
    #EPX = EPout[[1]]
    #EPF = EPout[[2]]
    #xeq.epf = EPout[[4]]
    #NEPF = length(EPF)
  }
  if(length(xeqboot) != 3){next}
  EQ3mat[i,] = xeqboot
  Ncase = Ncase + 1
  i3eq = c(i3eq,i)
  # print(c('row ',i,' done, total ',Ncase),quote=F) 
}

# -----------------------------------------------------------------------------
tstop = Sys.time()
print('----------------------------------------------------------',quote=F)
runtime = difftime(tstop,tstart,units='mins')
print(c('runtime for EPF equilibria, minutes ',runtime),quote=F)
print(c('confirmed cases of 3 stochastic (EPF) equilibria ',Ncase),quote=F)
print('----------------------------------------------------------',quote=F)
print('',quote=F)

# Select cases for density plots, means & sd, then make the plots
EQ3case = EQ3mat[i3eq,]  # save only the cases with 3 verified equilibria

# Do means of all bootstrap EPF have 3 equilibria?
mean.all = apply(EPFboot,1,mean)
x.all = apply(EPXboot,1,mean)
# find equilibria of mean 
xeq.mean = epfeq(x.all[1:199],mean.all,7)
print(c('equilibria of mean EPF ',xeq.mean),quote=F)
print('',quote=F)
#
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.8,cex.lab=1.8)
plot(x.all[1:199],mean.all,type='l',lwd=2,col='darkblue',xlab='Chlorophyll',
           ylab='EPF of mean, all boot cycles',main='dashed lines are eq of data')
grid()
abline(v=xeq.epf,lty=3,lwd=3,col='black')

# Do means of bootstrap EPF with 3 eq also have 3 eq?
# select EPFs for plots +/- s.d.
EPX3boot = EPXboot[,i3eq]
EPF3boot = EPFboot[,i3eq]
EPF3mean = apply(EPF3boot,1,mean)
EPF3sd = apply(EPF3boot,1,sd)
EPF3plus = EPF3mean + EPF3sd
EPF3minus = EPF3mean - EPF3sd
EPFxmean = apply(EPX3boot,1,mean)
# find equilibria of EPF3mean on EPFxmean
# epfeq=function(epx,epf,nknots)
xeq3 = epfeq(EPFxmean[2:200],EPF3mean,7)
print('Equilibria of bootstrap mean EPFs with 3 eq',quote=F)
print(xeq3)
# Mean and sd for equibria over all cases with 3 eq
meanxeq3 = colMeans(EQ3case)
sdxeq3 = apply(EQ3case,2,sd)
print('Mean equilibria for cases with 3 equilibria',quote=F)
print(meanxeq3)
print('SD of these means',quote=F)
print(sdxeq3)

# can we correct bias of bootstrap mean (when EPF has 3 eq) versus original data?
print('',quote=F)
print('=============================================================',quote=F)
print('Bootstrap bias estimate for equilibria at the mean',quote=F)
print('Deviations of xeq3 from xeq of the original data',quote=F)
print(round( (xeq3 - xeq.epf),3))  # note E&T use the order as xeq50 - xeq.EPF
print('Relative bootstrap bias for mean equilibria',quote=F)
RBB = (xeq3 - xeq.epf)/xeq.epf
print(round(RBB,3))
print('',quote=F)
print('bias-corrected equilibria ',quote=F)
BB = (xeq3 - xeq.epf)  # bootstrap bias by E&T text between eqs 10.40 & 10.41 p. 138
# On p 124 they define bias as estimator - population value 
# consistent with xeq.epf - xeq3:  
#BB = xeq.epf - xeq3
# However this correction moves the equilibria farther from peak & valleys of mean EPF
print(c('Bias based on mean EPF ',BB),quote=F)
print('Bias based on mean of all bootstrap EPF',quote=F)
print(meanxeq3 - xeq.epf)
print('compare to s.d. of mean eq and of EPFs in the plots',quote=F)
print(sdxeq3)
print('Note that SD of equilibria is larger than the bias',quote=F)
print('',quote=F)
#xeq.epf.adj = 2*xeq.epf - xeq3 # identical to xeq.epf - BB eq. 10.40
xeq.epf.adj = c(xeq.epf[1]-BB[1],xeq.epf[2],xeq.epf[3]+BB[3])
print(c('mean adjusted',round(xeq.epf.adj,3)),quote=F)
#print('the correction shifts the equilibria left of the peaks and valleys',quote=F)
print('',quote=F)
# plot result for cycles with 3 equilibria
# plot mean bootstrap EPF plus and minus sd
yrange = range(c(EPF3plus,EPF3minus),na.rm=T)
windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.8,cex.lab=1.8)
plot(EPFxmean[2:200],EPF3mean,ylim=yrange,type='l',lwd=3,col='darkblue',xlab='X',
     ylab='Effective Potential Mean +/- S.D.',
     main='Dash Lines:  Bias-corrected xeq')
points(EPFxmean[2:200],EPF3plus,type='l',lwd=2,col='deepskyblue')
points(EPFxmean[2:200],EPF3minus,type='l',lwd=2,col='deepskyblue')
grid()
abline(v=xeq.epf.adj,lty=3,lwd=3,col='black')
#abline(v=xeq.epf,lty=3,lwd=3,col='black')
#abline(v=xeq3,lty=3,lwd=3,col='red')
#abline(v=meanxeq3,lty=3,lwd=3,col='orchid')

# Plot densities for EQ3case, the cases of 3 equilibria

eqden1 = density(EQ3case[,1],bw='SJ',window="epanechnikov",n=512,na.rm='T')
eqden2 = density(EQ3case[,2],bw='SJ',window="epanechnikov",n=512,na.rm='T')
eqden3 = density(EQ3case[,3],bw='SJ',window="epanechnikov",n=512,na.rm='T')

windows(width=12,height=4)
par(mfrow=c(1,3),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.8,cex.lab=1.8)
plot(eqden1$x,eqden1$y,type='l',lwd=2,col='blue',xlab='Low Stable Eq',ylab='Density')
#abline(v=meanxeq3[1],lty=3,lwd=3,col='black')
plot(eqden2$x,eqden2$y,type='l',lwd=2,col='red',xlab='Unstable Eq',ylab='Density')
#abline(v=meanxeq3[2],lty=3,lwd=3,col='black')
plot(eqden3$x,eqden3$y,type='l',lwd=2,col='blue',xlab='High Stable Eq',ylab='Density')
#abline(v=meanxeq3[3],lty=3,lwd=3,col='black')

print('',quote=F)
print('Equilibria by 3 methods and s.d. of bootstraps',quote=F)
print('Data:',quote=F)
print(round(xeq.epf,4))
print('Mean EPF of boot cycles with 3 eq',quote=F)
print(round(xeq3,4))
print('Mean equilibria of boot cycles with 3 eq',quote=F)
print(round(meanxeq3,4))
print('SD of equilibria of boot cycles with 3 eq',quote=F)
print(round(sdxeq3,4))

# Data matrix for bar plot with errors
rnames = c('Stable Low','Unstable','Stable High')
cnames = c('Data','Mean EPF','Mean Eq')
meanbar = matrix(c(xeq.epf,xeq3,meanxeq3),nr=3,nc=3,byrow=F,
                 dimnames = list(rnames,cnames))
sdbar = matrix(rep(NA,9),nr=3,nc=3,byrow=F)
sdbar[,3] = sdxeq3

yrange = range(c(meanbar[,3]-sdbar[,3]-0.2,meanbar[3,]+sdbar[3,]+0.2),na.rm=T)

windows()
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.8,cex.lab=2)
bp_axes <- barplot(meanbar, 
                     beside = TRUE, 
                     col = c("deepskyblue", "red",'forestgreen'), 
                     ylim = yrange,
                     legend.text = TRUE,
                     args.legend = list(x = "topleft",cex=1.5),
                   xlab='',ylab='Equilibrium Chlorophyll',
                   main='Equilibria of Z(sqrt(Chlorophyll))')
xpos = bp_axes[,3]
yval = meanbar[,3]
err = sdbar[,3]
arrows(x0 = xpos, y0 = yval - err, 
       x1 = xpos, y1 = yval + err, 
       code = 3, angle = 90, length = 0.05, lwd = 2)
