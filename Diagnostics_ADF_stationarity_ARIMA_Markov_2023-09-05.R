# ADF stationarity test and ARIMA for Markov lag
# SRC 2023-09-05

rm(list = ls())
graphics.off()

library(stats)
library(tseries)

library(forecast)

library(parallel)
options(mc.cores = parallel::detectCores())

# load the data
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

#
# DLM outputs to consider:  b0, stdlevel, yhat, sm.y

# convert [Z-score of sqrt(chl)] to sqrt(Chl))
mu = Xtrans.lst[[3]]
sig = Xtrans.lst[[4]]
#x0 = sig*mat1$X0 + mu
#x1 = sig*mat1$X1 + mu

x0 = sm.y

# find optimal AR lags using auto.arima from forecast library
arfit = auto.arima(x0) # use defaults
#ncores=detectCores()
#arfit = auto.arima(Xvar0,stepwise=F,parallel=T,num.cores=ncores)  # use parallel processing
lagopt = arimaorder(arfit)
print('optimal order using autoarima() and arimaorder()',quote=F)
print(lagopt)
aropt = unname(lagopt[1])  # save optimal AR order
# if optimal lag is 0 then data are uncorrelated, use original data
aropt = ifelse(aropt==0,1,aropt)  
print(c('optimal order aropt = ',aropt),quote=F)

# test stationarity
ADF.result = adf.test(x0)
pvalue = ADF.result$p.value

print('ADF result',quote=F)
print('p value for HO: series NOT stationary',quote=F)
print(pvalue)
