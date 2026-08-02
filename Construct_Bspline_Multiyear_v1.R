# function to compute B spline and return D1 and sigma functions
# SRC 2025-11-30
require(splines)
#require(npreg)

# time series Xvar has been thinned, if necessary, to Markov time step
# X has been trimmed or transformed if necessary
# nx is length of Xvar
# nk is number of knots for mu or sigma function
# npoly is polynomial order for mu or sigma function
Bspline.Langevin = function(x0,x1,dx,nx,DT,nk.mu,nk.sig,npoly.mu,npoly.sig) {
  # Make variates for first and second moments
  #x0 = Xvar0[1:(nx-1)]
  #x1 = Xvar0[2:nx]
  dx = (x1-x0)/DT   # first moment
  dx2 = 0.5*(x1-x0)^2/DT  # second moment 
  dsd = sqrt(2*dx2)  # sigma?
  # Make data frame and sort by x0
  DF0 = as.data.frame(cbind(x0,x1,dx,dx2,dsd))
  DF1 = DF0[order(DF0$x0),]  # sort by x0; in earlier versions named DF0s 
  # SPLINE THE MU FUNCTION =============================================
  # Quantiles for knots of mu function
  xqt.mu = quantile(DF1$x0,probs=seq(0,1,length.out=nk.mu),na.rm=T)
  # Internal knots based on quantiles
  kmid.mu = xqt.mu[2:(nk.mu-1)]
  #
  # Defaults for bs() are order = 3, boundary knots are limits of data
  basisD1 = bs(DF1$x0,knots=kmid.mu,degree=npoly.mu,intercept=F,Boundary.knots = range(DF1$x0))
  lmD1 = lm(dx ~ basisD1 -1,data=DF1,x=T,y=T)  # least-squares fit of the B-spline
  D1mat = unname(lmD1$x) # remove row and col names from design matrix
  D1y = unname(lmD1$y) # dx conforming to spline design
  # Solve spline for weights with normal equations
  D1D1 = t(D1mat)%*%D1mat
  iD1D1 = solve(D1D1)
  # as.vector strips the row names from the result
  wD1 = as.vector(iD1D1%*%t(D1mat)%*%D1y) # weights for D1 spline
  print('Solution complete for D1 spline',quote=F)
  print('summary of linear spline fit of D1 function',quote=F)
  print(summary(lmD1))
  rdf.D1 = lmD1$df.residual  # residual d.f.
  print(c('residual df = ',rdf.D1),quote=F)
  print('dimension of D1mat',quote=F)
  print(dim(D1mat))
  #
  # Find errors of dx model
  lmD1err = lmD1$residuals  # from lm()
  # We want to predict x1 = x0 + D1*DT which has error x1 - x0 - D1*DT
  eP = DF1$x1 - DF1$x0 - (D1mat%*%wD1)*DT  # which should match lmD1err!
  # Plot shows that errors are identical by lm() & normal equations
  # Comment out the plot unless we need confirmation
  #windows()  
  #plot(lmD1err,eP,type='p',pch='+')
  
  # CONSTRUCT AND SPLINE THE SIGMA FUNCTION =======================================
  # Sigma vector from eP is sqrt(2*D2) 
  DF1$ep2 = 0.5*(eP^2)  # D2
  DF1$sep = sqrt(2*DF1$ep2)
  # Build spline for sqrt(error^2), i.e. sqrt(D2 corrected for mean) for effective potential
  # Quantiles
  xqt.sig = quantile(DF1$x0,probs=seq(0,1,length.out=nk.sig),na.rm=T)
  print('quantiles',quote=F)
  print(xqt.sig,quote=F)
  
  # Internal knots based on quantiles
  kmid.sig = xqt.sig[2:(nk.sig-1)]
  
  # Fit sigma 
  # Defaults for bs() are order = 3, boundary knots are limits of data
  basisig = bs(DF1$x0,knots=kmid.sig,degree=npoly.sig,intercept=F,Boundary.knots = range(DF1$x0))
  lmsig = lm(sep ~ basisig -1,data=DF1,x=T,y=T)
  print('summary of linear spline fit of sigma function',quote=F)
  print(summary(lmsig))
  lmsigerr = lmsig$residuals
  sigmat = unname(lmsig$x) # remove row and col names from sigma design matrix
  print('dimensions of design matrix for sigma',quote=F)
  print(dim(sigmat))
  sigy = unname(lmsig$y) # sigma conforming to spline design
  # Solve spline of sigma for weights using normal equations
  sig2 = t(sigmat)%*%sigmat
  isig2 = solve(sig2)
  # as.vector strips the row names from the result
  wsig = as.vector(isig2%*%t(sigmat)%*%sigy) # weights for sigma spline
  print('Solution complete for sigma spline',quote=F)
  
  # B-spline least-squares estimates of D1 and sigma
  LSD1 = D1mat%*%wD1
  LSsigma = sigmat%*%wsig 
  
  # Equilibria from spline fit of D1
  sdrift = sign(LSD1)
  dsdrift = c(0,-diff(sdrift))
  xeq = DF1$x0[which(!dsdrift == 0)]
  ixeq = which(!dsdrift == 0)  # indices of the equilibria
  
  print('',quote=F)
  print('equilibria of D1 on x0 axis',quote=F)
  print(xeq,quote=F)
  print(ixeq,quote=F)  
  
  xeq.D1 = xeq
  outlist=list(DF1,LSD1,LSsigma,xeq.D1,lmD1err,lmsigerr,
               basisD1,D1y,D1mat,D1D1,iD1D1,rdf.D1,
               basisig,sigmat,isig2)
  print('B spline result list contains DF1,LSD1,LSsigma,xeq.D1, & errors',quote=F)
  print('DF1 data frame includes x0, the x axis of the splines',quote=F)
  print('LSD1 is D1 by B-spline, LSsigma is sigma by B-spline',quote=F)
  print('xeq.D1 is equilibria of D1',quote=F)
  print('lmD1err and lmsigerr are residuals of D1 and sigma fits',quote=F)
  print('basisD1, D1mat, and iD1D1 are constant matrices for D1 spline',quote=F)
  print('D1y is dx vector,D1mat is design matrix X, D1D1 is t(X)X, iD1D1 is inverse',quote=F)
  print('rdf.D1 is residual degrees of freedom D1 spline',quote=F)
  print('basisig, sigmat, and isig2 are constant matrices for sigma spline',quote=F)
  return(outlist)
}