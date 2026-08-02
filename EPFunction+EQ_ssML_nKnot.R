# Function returns effective potential and equilibria given X, D1, sigma
# SRC 2023-03-10

#require(stats)
#require(moments)
#require(bvpSolve)
require(cubature)
require(numDeriv)
require(npreg)

EPFEQ = function(X,D1,sigma,nknots)  {
  # convert xrate and sigma to functions of x
  #
  #D1mod = ss(D1x,D1y,method='ML',m=2,all.knots=T) # m=2 is cubic spline
  #nknots=7
  xqt = quantile(X,probs=seq(0,1,length.out=nknots),na.rm=T) # quantiles of X used for knots
  D1spline = ss(x=X,y=D1,method='ML',m=2,knots=xqt)
  D1fun = function(x) {
    yhat=predict(D1spline,x)$y
    return(yhat)
  }
  #
  # Make a D2 spline
  D2 = 0.5*sigma^2
  D2spline = ss(x=X,y=D2,method='ML',m=2,knots=xqt)
  D2fun = function(x)  {
    yhat = predict(D2spline,x)$y
    return(yhat)
  }
  #
  # Spline for D1/D2 ratio
  D12 = D1fun(X)/D2fun(X)
  D12spline = ss(x=X,y=D12,method='ML',m=2,knots=xqt)
  D12fun = function(x)  {
    yhat = predict(D12spline,x)$y
    return(yhat)
  }
  #
  # Calculate effective potential function
  X.ep = seq(min(X)+0.01,max(X)-0.01,length.out=200)
  EPF = rep(0,199)
  for(i in 1:199) {
    x0 = X.ep[1]
    x1 = X.ep[i+1]
    xhalf = (x0 + x1)/2 # midpoint of discrete interval
    integral = hcubature(f=D12fun,lowerLimit=x0,upperLimit=x1)$integral
    loghalf = ifelse(D2fun(xhalf)>=0,log(D2fun(xhalf)),0)  
    EPF[i] = -1*integral + loghalf # function from Science paper & Tabar eq 4.17
    #EPF[i] = -2*integral + 2*log(sigfun(xhalf)) # Babak 2023-03-07 email
    #EPF[i] = -1*integral + 2*log(sigfun(xhalf)) # Babak 2023-03-07 without the 2 
  }
  # See Proto_EffectivePotentialFunction_2023-03-09.R in /SRC/Cascade/Bandi_Mesh_Tau_March2023
  # The Science version (identical to Tabar eq 4.17) and 
  #   Babak's version without the 2 give the same result.
  # The standardization below has no effect on the roots 
  # Axis shift for alternate effective potential
  minEPF = min(EPF,na.rm=T)
  EPF = EPF - minEPF + 0.5
  
  Xep.mid = (X.ep[1:199] + X.ep[2:200])/2
  dmid = Xep.mid[100] - Xep.mid[99]  # x step at Xep.mid
  EPdat0 = as.data.frame(cbind(Xep.mid,EPF))
  EPdat = na.omit(EPdat0)
  
  # convert EPF to a numerical function
  #xEPqt = quantile(X.ep[2:100],probs=seq(0,1,length.out=nknots),na.rm=T) # quantiles of X used for knots
  #EPspline = ss(x=X.ep[2:100],y=EPF,method='ML',m=2,knots=xEPqt)
  #xEPqt = quantile(EPdat[,1],probs=seq(0,1,length.out=nknots),na.rm=T) # quantiles of X used for knots
  # EPspline = ss(x=EPdat[,1],y=EPdat[,2],method='ML',m=2,knots=xEPqt) # original knots
  EPspline = ss(x=EPdat[,1],y=EPdat[,2],method='ML',m=2,all.knots = T) # all points are knots
  EPfun = function(x)  {
    yhat = predict(EPspline,x)$y
    return(yhat)
  }
  
  # take first derivative and find the roots
  dEPdx = -1*grad(EPfun,X.ep,'Richardson')
  # Potential is a negative integral; reverse sign of derivative for plots
  # find roots of first derivative of effective potential
  sdrift = sign(dEPdx)
  dsdrift = c(0,-diff(sdrift))
  xeq = Xep.mid[which(!dsdrift == 0)] - dmid
  #
  outlist = list(X.ep,EPF,dEPdx,xeq,D12)
  return(outlist)
  #
} # end EPFEQ