
.First.lib <- function(lib, pkg) {
  library.dynam("timeslab", pkg, lib)
}


acf <- function(x, m = 0) 
#------------------------------------------------------------------
# 
#   Function to find acf of lags 1 through m via 2 FFT's
#
#------------------------------------------------------------------
{
   if(m == 0) return(dot(x - mean(x), x - mean(x)) / length(x))

   n   <- length(x)
   Q   <- 2^(trunc(log(n + m, base = 2) + 0.05) + 1)
   if(Q<n+m) Q<-2*Q

   z   <- c(x - mean(x), rep(0, Q - n))

   fz  <- Re( fft( Mod(fft(z))^2 / n, inverse = TRUE)) / Q

   list(var = fz[1], corr = fz[2:(m + 1)]/fz[1])
}
acf1 <- function(x, m = 0) 
#-------------------------------------------------------------------
# 
#   Function to find acf for lags 1 through m via convolution
#
#-------------------------------------------------------------------
{
   var <- tsvar(x)
   if(m == 0) return(var)

   n   <- length(x)
   x   <- c( x - mean(x), rep(0,m) )
   y   <- x
   rho <- rep(0,m)

   for(i in 1:m)  rho[i] <- dot( x, (y <- crlag(y))) / (n * var)

   list(var = var, corr = rho)
}
araic <- function(x, maxp)
#---------------------------------------------------------------------
#
#   Function to find AIC for AR(p) processes of orders 0 through maxp
#          from data x
#
#   returns: $p    : best order
#            $aic0 : aic for order 0
#            $aic  : aic for orders 1, ..., maxp
#            $part : pacf for lags  1, ..., maxp
#
#---------------------------------------------------------------------
{
   lvar <- log(tsvar(x))
   n    <- length(x)
   part <- pacf(x,maxp)
   aic  <- n * (lvar + cumsum(log(1 - part^2))) + 2 * c(1:maxp)
   aic0 <- n * lvar

   list(p = order(c(aic0, aic))[1] - 1, aic = aic, aic0 = aic0, part = part)
}
arcorr <- function(alpha,rvar=1,m=0)
#---------------------------------------------------------------------
#
#   Function (Fortran) to find AR variance and autocorrelations.   
#
#----------------------------------------------------------------------
{
   return(armacorr(alpha,c(),rvar,m))
}
ardt <- function(alpha,rvar,n,seed=0)
#-----------------------------------------------------------------------
#
#   Function to simulate AR data.
#
#-----------------------------------------------------------------------
{
   x <- c()
   p <- length(alpha)
   z <- arcorr(alpha,rvar,p)
   if(z$ier!=0) return(list(x=x,ier=z$ier))
   x <- diffeq(alpha,p,n,c(corrdt(z$corr,z$var,p,seed),
                           sqrt(rvar)*wn(seed,n-p)))
   return(list(x=x,ier=0))
}
arfilt <- function(alpha,rvar,x)
#--------------------------------------------------------------------
#
#   Function to apply an AR filter to a matrix.
#
#   returns: $w :  new matrix
#            $ier: 0 means process stationary, j (> 0) means jth partial
#                  outside (-1,1)
#
#---------------------------------------------------------------------
{
   p <- length(alpha)
   if(is.matrix(x)) { nr <- nrow(x);   nc <- ncol(x)}
   if(is.vector(x)) { nr <- length(x); nc <- 1} 

   z <- .Fortran("arfilt", 
                  as.double(alpha), 
                  as.integer(p), 
                  as.double(rvar),
                  as.double(x), 
                  as.integer(nr), 
                  as.integer(nc), 
              w = as.double(matrix(0,nr,nc)), 
                  as.double(rep(0,p)),
                  as.double(matrix(0,p,p)), 
            ier = as.integer(ier <- 0))

   if(is.matrix(x)) z$w <- matrix(z$w,nr,nc)

   list(w = z$w, ier = z$ier)
}
arma <- function(alpha,beta,x,iopt,p,q,rvar,n,m,seed=0)
#--------------------------------------------------------------------------
#
#   Function to illustrate patterns in ARMA parametrizations.
#
#   iopt: 0  means use input alpha, beta, and x
#         1  means use input alpha and beta and simulate x
#         2  means simulate alpha, beta, and x
#
#--------------------------------------------------------------------------
{
   if(seed!=0) set.seed(seed)

   par(mfrow=c(4,1))
   
   if(iopt==2) {
      alpha <- simcfs(p)
      beta  <- simcfs(q) }

   if(iopt>=1) x <- armadt(alpha,beta,rvar,n,seed)$x

   options(digits=3)

   lab <- paste("ARMA(",p,",",q)
   if(p!=0) lab <- paste(lab,",(",paste(format(alpha),collapse=","),")")
   if(q!=0) lab <- paste(lab,",(",paste(format(beta),collapse=","),")")
   lab <- paste(lab,",",format(rvar),")")

   plot(x,ylab="value",type="l",main=lab)

   z      <-   armacorr(alpha,beta,rvar,m)
   rho    <-   z$corr
   r0     <-   z$var
   z      <-   acf(x,m)
   rhohat <-   z$corr
   r0hat  <-   z$var

   plot(rho,xlab="lag",ylab="value",ylim=c(-1,1),
        main="ACF (True is solid curve)",type="l")
   points(c(1:m),rhohat)

   theta  <-   armapart(alpha,beta,rvar,m)$theta
   thetahat <- pacf(x,m)

   plot(theta,xlab="lag",ylab="value",ylim=c(-1,1),
        main="PACF (True is solid curve)",type="l")
   points(c(1:m),thetahat)

   f    <- stdf(armasp(alpha,beta,rvar,n),r0,exp(-6),exp(6))
   fhat <- stdf(perdgm(x),r0hat,exp(-6),exp(6))

   matplot(freqs(n),cbind(log(f),log(fhat)),xlab="frequency",ylab="value",
           main="Log of Spectral Density",type="l",ylim=c(-6,6))
 
   return(list(alpha=alpha,beta=beta,x=x))
}





armacorr <- function(alpha,beta,rvar=1,m)
#---------------------------------------------------------------------
#
#   Function to find ARMA correlations.   
#
#----------------------------------------------------------------------
{

   p    <-  length(alpha)
   q    <-  length(beta)
   if(p==0) alpha <- c(0)
   if(q==0) beta  <- c(0)

   ma   <-  m+1
   mx   <-  max(p,q)
   mxp1 <-  mx+1 

   z <- .Fortran("wilson", 
                 as.double(-alpha), 
                 as.integer(p),
                 as.double(-beta),
                 as.integer(q),
             acf=as.double(rep(0,ma)),
                 as.integer(ma),
                 as.double(rep(0,mxp1)),
                 as.integer(mxp1),
                 as.double(rep(0,mxp1)),
                 as.integer(mx),
             ier=as.integer(ier <- 0))

   list(var=rvar*z$acf[1],corr=z$acf[2:ma]/z$acf[1],ier=z$ier)
}

armadt <- function(alpha,beta,rvar,n,seed=0)
#-----------------------------------------------------------------------
#
#   Function to simulate ARMA data.
#
#-----------------------------------------------------------------------
{
   p <- length(alpha)
   q <- length(beta)
   if(p==0) alpha <- c(0)
   if(q==0) beta <- c(0)
   z <- armacorr(alpha,beta,rvar,p)
   if(z$ier!=0) return(list(x=c(),ier=z$ier))
   x <- diffeq(alpha,p,n,c(corrdt(z$corr,z$var,p,seed),
                           madt(beta,rvar,n-p,seed)))
   return(list(x=x,ier=0))
}
armapart <- function(alpha,beta,rvar,m)
#-----------------------------------------------------------------------
#
#   Function to find ARMA partial autocorrelations.
#
#------------------------------------------------------------------------
{
   z <- armacorr(alpha,beta,rvar,m)
   if(z$ier!=0) return(list(theta=c(),ier=z$ier))
   z <- corrar(z$corr,z$var,m)
   return(list(theta=arpart(z$alpha)$theta,ier=0))
}
armapred <- function(x,alpha,beta,rvar,t1,t2,h1,h2)
#--------------------------------------------------------------------
#
#   Function to find ARMA predictors.
#
#   Find predictors and prediction error stndard deviation
#   of x(t+h) given x(t),...,x(1) for t=t1,...,t2 and h=h1,...,h2
#
#   Returns a list:
#      $xp, $se : predictors and standard deviations (those for t=t1
#                 come first, then those for t=t1+1, and so on)
#
#---------------------------------------------------------------------
{
   alph <- c(0)
   bet  <- c(0)
   p <- length(alpha)
   q <- length(beta)
   if(p+q == 0) stop("One of p or q must be positive in armapred()")

   if(p > 0) alph <- alpha
   if(q > 0) bet  <- beta

   n <- length(x)
   if(n < 1 || t1 > t2 || h1 > h2 || h1 < 1 || t1 < 1 || t2 > n || t1 > n)
      stop("Illegal input in armapred()")

   np  <- (t2-t1+1)*(h2-h1+1)
   nr  <- max(p,q+1)
   nr2 <- nr*nr

   z <- .Fortran("mxpd",
                 as.double(x),
                 as.double(alph),
                 as.double(bet),
                 as.integer(p),
                 as.integer(q),
                 as.double(rvar),
                 as.integer(t1),
                 as.integer(t2),
                 as.integer(h1),
                 as.integer(h2),
                 as.integer(nr),
                 as.double(st  <- rep(0,nr2)),
                 as.double(w   <- rep(0,nr2)),
                 as.double(a   <- rep(0,nr2)),
                 as.double(sh  <- rep(0,nr2)),
                 as.double(st1 <- rep(0,nr2)),
                 as.double(st2 <- rep(0,nr2)),
              xp=as.double(xp  <- rep(0,np)),
              se=as.double(se  <- rep(0,np)))

   return(list(xp=z$xp,se=z$se))

}
armasp <- function(alpha,beta,rvar=1,Q=256)
#------------------------------------------------------------------------
#
#   Function to find ARMA spectral density.
#
#   Find ARMA(p,q,alpha,beta,rvar) spectral density at
#   Q frequencies in [0,1)  and return first [Q/2] + 1
#
#------------------------------------------------------------------------
{
   f1 <- f2  <- rep(1, Q)

   p <- length(alpha)
   q <- length(beta)

   if(q != 0) f1 <- Re(Mod(fft(c(1, beta,  rep(0, Q - q - 1))))^2)

   if(p != 0) f2 <- Re(Mod(fft(c(1, alpha, rep(0, Q - p - 1))))^2)

   f <- rvar * f1 / f2
   f <- f[1:((Q / 2) + 1)]
}
arpart <- function(alpha)
#-------------------------------------------------------------------
#
#   Function (Fortran) to find AR partial autocorrelations.
#
#-------------------------------------------------------------------
{
   p <- length(alpha)
   if(p==0) return()

   z <- .Fortran("arpart", 
           alpha=as.double(alpha), 
                 as.integer(p),
             ier=as.integer(ier <- 0))

   list(theta=z$alpha,ier=z$ier)
}
arsp <- function(alpha,rvar=1,Q=256)
#------------------------------------------------------------------------
#
#   Function to find AR spectral density.
#
#   Find AR(p, alpha, rvar) spectral density at Q frequencies
#   in [0,1) and return first [Q/2] + 1
#
#-----------------------------------------------------------------------
{
   f <- rep(1, Q)
   p <- length(alpha)

   if( p != 0) f <- Re(Mod(fft(c(1, alpha, rep(0, Q - p - 1))))^2)

   f <- rvar / f[1:((Q/2)+1)]
}
 arsppeak <- function(alpha,rvar,n,start=0)
#----------------------------------------------------------------------
#
#   Function to determine a peak frequency of an AR process and its
#   standard error.
#
#   Input: 
#
#      alpha, rvar: AR parameters
#      n          : sample size
#      start      : initial guess at peak frequency (default value of
#                   0 will get starting value as largest relative max
#                   of f evaluated at 256 frequencies).
#
#   Returns:
#
#      $peakf : estimate of peak frequency 
#      $stderr: estimated standard error
#      $ier   : 0 means successful completion
#               1 means no peaks
#               2 means 0 second derivative encountered
#               3 means convergence to 0 or .5
#               4 means nonconvergence
#
#---------------------------------------------------------------------
{
   p    <- length(alpha)
   if(p < 1 || n < 1 || start < 0 || start >= .5 || rvar <= 0) 
      stop("illegal input to arsppeak()")
   f    <- arsp(alpha,rvar,256)

   z    <- .Fortran("arsppk",
                 as.double(alpha),
                 as.integer(p),
                 as.double(rvar),
                 as.double(start),
                 as.integer(n),
                 as.double(f),
                 as.double(rep(0,p*p)),
                 as.double(rep(0,p)),
                 as.double(rep(0,p*p)),
                 as.double(rep(0,p)),
                 as.double(rep(0,p*p)),
             ier=as.integer(ier <- 0),
           peakf=as.double(peakf <- 0),
              se=as.double(se <- 0))
 
   return(list(peakf=z$peakf, se=z$se, ier=z$ier))
}
clip <- function(x,low=-1.e20,up=1.e20)
#------------------------------------------------------------------------
#
#   Function to clip an array.
#
#   Replace any element of x that is less than low by low and any element
#   greater than up by up.
#
#------------------------------------------------------------------------
{
   ind <- c(1:length(x))
   x   <- replace(x,ind[x<low],low)
   x   <- replace(x,ind[x>up],up)
}
coeffcsd <- function(alpha,beta,n)
#----------------------------------------------------------------
#
#   Function to find asymptotic standard error of ARMA MLE's.
#
#   Returns:
#      $ier : error indicator (0 no error, 1 error)
#      $sd  : standard errors 
#
#----------------------------------------------------------------
{
   p <- length(alpha)
   q <- length(beta)
   if(p+q == 0) stop("One of p or q must be positive in coeffcsd()")

   if(p == 0 && q > 0 ) return(list(sd=sqrt(diag(schur(beta))/n),ier=0))
   if(p > 0  && q == 0) return(list(sd=sqrt(diag(schur(alpha))/n),ier=0))
   mdim <- 2*max(p,q)+1

   z <- .Fortran("mxcsd",
                 as.double(alpha),
                 as.double(beta),
                 as.integer(p),
                 as.integer(q),
                 as.integer(n),
                 as.integer(mdim),
                 as.double(matrix(rep(0,mdim*mdim),mdim,mdim)),
                 as.double(rep(0,mdim)),
                 as.double(matrix(rep(0,mdim*mdim),mdim,mdim)),
              sd=as.double(rep(0,p+q)),
             ier=as.integer(ier <- 0))

   list(sd=z$sd,ier=z$ier)
}
corrar <- function(rho,r0,p)
#---------------------------------------------------------------------
#
#   Function to find AR parameters from correlations.   
#
#----------------------------------------------------------------------
{
   z <- .Fortran("corrar", 
                 as.double(rho), 
                 as.double(r0),
                 as.integer(p),
           alpha=as.double(rep(0,p)),
            rvar=as.double(rvar <- 0))

   list(alpha=z$alpha,rvar=z$rvar)
}
corrarma <- function(rho,r0,p,q,maxit=100,del=1.e-5)
#---------------------------------------------------------------------
#
#   
#
#----------------------------------------------------------------------
{
   r    <- rho*r0
   ndim <- p+q
   if(length(r)<ndim) stop("length(rho) must be at least p+q")

   z <- .Fortran("cvmx1", 
                 as.double(r), 
                 as.double(r0),
                 as.integer(p),
                 as.integer(q),
                 as.integer(ndim),
                 as.integer(maxit),
                 as.double(del),
                 as.integer(ip <- rep(0,ndim)),
                 as.double(al <- matrix(rep(0,ndim*ndim),ndim,ndim)),
           alpha=as.double(alpha <- rep(0,ndim)),
            beta=as.double(beta <- rep(0,ndim)),
            rvar=as.double(rvar <- 0),
             ier=as.integer(ier <- 0))

   alpha <- c()
   beta  <- c()
   if(p>0) alpha <- z$alpha[1:p]
   if(q>0) beta  <- z$beta[1:q]

   list(alpha=alpha, beta=beta, rvar=z$rvar, ier=z$ier)
}


corrdt <- function(rho,r0,n,seed=0)
#----------------------------------------------------------------------
#
#   Function to simulate data having specified correlations.
#
#   Generate realization of length n from Gaussian process having
#   variance r0 and correlations rho.
#
#----------------------------------------------------------------------
{
   if(n==0) return(c())
   return(t(chol(toepl(rho*r0,r0,n))) %*% wn(seed,n))
}
corrma <- function(rho,r0,q,maxit=100,del=1.e-5)
#---------------------------------------------------------------------
#
#   
#
#----------------------------------------------------------------------
{
   r    <- rho*r0
   ndim <- q
   if(length(r)<ndim) stop("length(rho) must be at least q")

   z <- corrarma(rho,r0,0,q,maxit,del)   

   return(list(beta=z$beta,rvar=z$rvar,ier=z$ier))

}


crlag <- function(x) 
#------------------------------------------------------------------
# 
#   Function to apply circular shift operator to an array.
#
#------------------------------------------------------------------
{
   lx <- c(x[length(x)], x[1:(length(x) - 1)])
}
#delay <- function(seconds)
#----------------------------------------------------------------------
#
#   Function to delay execution of S for seconds seconds.
#
#---------------------------------------------------------------------
#{
#   unix(paste("sleep ",seconds))
#   invisible()
#}
diffeq <- function(alpha,p,n,e)
#---------------------------------------------------------------------
#
#   Function to perform a difference equation.
#
#----------------------------------------------------------------------
{
   if(p==0) alpha <- c(0)
   z <- .Fortran("diffeq", as.double(alpha), as.integer(p), as.integer(n),
               as.double(e), x = as.double(rep(0,n)))
z$x
}
divpoly <- function(num,den,n)
#---------------------------------------------------------------------
#
#   Function to divide polynomials.
#
#   Find first n coefficients of 
#
#      h(x) = num(x)/den(x)
#
#   where num(x) and den(x) are polynomials having zeroth coefficient
#   equal to 1.
#
#----------------------------------------------------------------------
{
   p <- length(den)
   q <- length(num)
   if(p+q == 0) return(rep(0,n))

   alpha <- c(0)
   beta  <- c(0)
   if(p > 0) alpha <- den
   if(q > 0) beta  <- num

   z <- .Fortran("mxma",
                 as.double(alpha),
                 as.double(beta),
                 as.integer(p),
                 as.integer(q),
                 as.integer(n),
           ratio=as.double(ratio <- rep(0,n)))

   return(z$ratio)
}
dot <- function(x,y) 
#-------------------------------------------------------------------
#
#   Function to find inner product of vectors.
#
#------------------------------------------------------------------
{
   dot <- sum(x * y)
}
dtar <- function(x, p)
#---------------------------------------------------------------------
#
#   dtar: Fit AR(p) (if p > 0) or AR(m) (if p < 0, m = AIC order) to x 
#
#   returns: $p      : order
#            $rvar   : estimate of white noise variance
#            $alpha  : coefficients
#
#----------------------------------------------------------------------
{

   if(p > 0)
   {
      part  <- pacf(x,p)
      alpha <- partar(part)
      rvar  <- exp(log(tsvar(x)) + sum(log(1 - part^2)))
      return(list(p = p, alpha = alpha, rvar = rvar))
   }

   else if(p < 0)
   {
      z  <- araic(x, -p)
      m  <- z$p

      if(m == 0) return(list(p = 0, alpha = c(0), rvar = tsvar(x)))

      if(m != 0) return(list(p = m, alpha = partar(z$part[1:m]),
           rvar = exp(log(tsvar(x)) + sum(log(1 - z$part[1:m]^2)))))
   }

}
dtarma <- function(x,alpha,beta,maxit=200,eps=.0001)
#-----------------------------------------------------------------
#
#   Function to find ARMA MLE's.
#
#   Arguments:
#      x  : data
#      alpha: starting values for AR parameters ( 0's usually work,
#             c() if doing MA)
#      beta : starting values for MA parameters ( 0's usually work,
#             c() if doing AR)
#
#   Return:
#      $ier : termination information ( <= 0 means convergence, 1 means 
#             nonconvergence, > 1 means error)
#
#      $alpha, $beta, $rvar : parameter values at termination
#      $m2ll : -2log(likelihhod) at terminating values
#      $var  : measure of degree of convergence at termination
#
#--------------------------------------------------------------------
#
{
   p <- length(alpha) 
   q <- length(beta)
   if((p+q)==0) stop("One of p or q must be nonzero in DTARMA")

   if(p > 0) {
      if(ap <- arpart(alpha)$ier != 0) 
         stop("Starting alpha not stationary") }

   if(q > 0) {
      if(bp <- arpart(beta)$ier != 0)  
         stop("Starting beta not invertible") }

   if(p > 0) ap <- -log((1-ap)/(1+ap))
   if(q > 0) bp <- -log((1-bp)/(1+bp))

   if(p == 0) ap <- c(0)
   if(q == 0) bp <- c(0)

   n <- length(x)

   if(n  == 0) stop("Data undefined in dtarma()")

   z <- .Fortran("dtarma",
                 as.double(c(ap,bp)),
                 as.double(eps),
                 as.integer(maxit),
                 as.double(x),
                 as.integer(n),
                 as.double(e <- rep(0,n)),
                 as.integer(p),
                 as.integer(q),
           alpha=as.double(alpha <- rep(0,p+q)),
            beta=as.double(beta  <- rep(0,p+q)),
                 as.integer(p+q),
                 as.double(pmat <- matrix(rep(0,(p+q)*(p+q+1)),p+q,p+q+1)),
              rv=as.double(rv <- 1),
            m2ll=as.double(m2ll <- 0),
             ier=as.integer(ier <- 0),
             var=as.double(var <- 0))

   alp <- c()
   bet <- c()
   if(p > 0) alp <- z$alpha[1:p]
   if(q > 0) bet <- z$beta[1:q]

   list(alpha=alp,beta=bet,ier=z$ier,rvar=z$rv,var=z$var,m2ll=z$m2ll)
}



filt <- function(beta,beta0,x)
#---------------------------------------------------------------------
#
#   Function filter an array.   
#
#----------------------------------------------------------------------
{
   m <- length(beta)
   if(m==0) return(beta0*x)
   n <- length(x)

   z <- .Fortran("filt", 
                 as.double(beta),
                 as.double(beta0),
                 as.integer(m),
                 as.integer(n),
                 as.double(x),
               y=as.double(rep(0,n-m)))

   return(z$y)
}
freqs <- function(n)
#-------------------------------------------------------------------
#
#   Function to form an array of frequencies in [0,.5]
#
#   Form the array of length (n / 2) + 1 :
#
#    f[i] = (i - 1) / n,   i = 1, ..., (n / 2 ) + 1
#
#-------------------------------------------------------------------
{
   f <- c(0:(n/2))/n
}
infqnt <- function(x)
#-----------------------------------------------------------------------
#
#   Function to plot informative quantile plot of x.
#
#-----------------------------------------------------------------------
{
   n <- length(x)
   xs <- sort(x)
   infq <- clip((xs - median(xs))/(2*(xs[3*n/4] - xs[n/4])),-1,1)
   vs <- ((1:n) - .5)/n
   plot(vs,infq,xlab="u",ylab="Value",ylim=c(-1,1),
         main="Informative Quantile Plot",type="l",lwd=3)
   lines(c(0,1),c(-.5,.5))
   lines(c(0,0,1,1,0),c(-1,1,1,-1,-1))
}
madt <- function(beta,rvar,n,seed=0)
#----------------------------------------------------------------------
#
#   Function to simulate MA data.
#
#   Simulate MA(q,beta,rvar) data of length n
#
#----------------------------------------------------------------------
{
   x <- filt(beta,1,sqrt(rvar) * wn(seed,n+length(beta)))
}
masmooth <- function(x,k)
#-------------------------------------------------------------------
#
#   Function to apply moving average smoother of length k.
#
#------------------------------------------------------------------
{
   z <- x
   if(k == 0) return(z)

   n <- length(x)
   w <- c(x[(k+1):2], x, x[(n-1):(n-k)])

   for(i in 1:k)  z <- z + w[(k+1-i):(k-i+n)] + w[(k+1+i):(k+i+n)]

   z <- z / (2 * k + 1)
}
macorr <- function(beta,rvar=1,m)
#---------------------------------------------------------------------
#
#   Function to find MA correlations.   
#
#----------------------------------------------------------------------
{
   return(armacorr(c(),beta,rvar,m))
}
masp <- function(beta,rvar=1,Q=256)
#------------------------------------------------------------------------
#
#   Function to find MA spectral density.
#
#   Find MA(q,beta,rvar) spectral density at Q frequencies
#   in [0,1) and return first [Q/2] + 1
#
#-----------------------------------------------------------------------
{
   f <- rep(1, Q)
   q <- length(beta)

   if( q != 0) f <- Re(Mod(fft(c(1, beta, rep(0, Q - q - 1))))^2)

   f <- rvar * f[1:((Q/2)+1)]
}
 movave <- function(x,k)
#---------------------------------------------------------------------
#
#   Function to apply moving average smoother (Fortran)
#
#
#----------------------------------------------------------------------
{
   k2 <- 2*k+1
   n  <- length(x)
   z  <- .Fortran("movave", 
                  as.double(x),  
                  as.integer(n),  
                  as.integer(k2),
              y = as.double(rep(0,n)))
z$y
}
movbox <- function(x,k)
#---------------------------------------------------------------------
#
#   Function to get information for moving box plot.
#
#----------------------------------------------------------------------
{
   n  <- length(x)
   k2 <- 2 * k + 1
   z <- .Fortran("movbox", 
                 as.double(x), 
                 as.integer(n),
                 as.integer(k2), 
                 as.integer(rep(0,k2)),
                 as.integer(rep(0,k2)),
          summ = as.double(matrix(rep(0,5*n),5,n)),
         nouts = as.integer(0),
        indout = as.integer(rep(0,n)),
        outval = as.double(rep(0,n)))

   if(z$nouts == 0) { z$indout <- NULL; z$outval <- NULL}

   if(z$nouts != 0) { z$indout <- z$indout[1:z$nouts]; 
                      z$outval <- z$outval[1:z$nouts]}

   list(summ=t(matrix(z$summ,5,n)), inds=z$indout, outs=z$outval)

}

movord <- function(x,nord,k)
#---------------------------------------------------------------------
#
#   Function to get moving order statistics.
#
#----------------------------------------------------------------------
{
   k2 <- 2 * k + 1
   n  <- length(x)
   z  <- .Fortran("movord", 
                  as.double(x),  
                  as.integer(n),  
                  as.integer(k2),
                  as.integer(nord),
                  as.integer(rep(0,k2)), 
                  as.integer(rep(0,k2)),
              y = as.double(rep(0,n)))
z$y
}
multpoly <- function(alpha,beta)
#--------------------------------------------------------------------
#
#   Function to multiply polynomials.
#
#   Return the first n coefficients of the product of the polynomials
#   having coefficients alpha and beta (these polynomials have zeroth
#   coefficient equal to 1).
#
#--------------------------------------------------------------------
{
   p <- length(alpha)
   q <- length(beta)
   if(p == 0 || q == 0) stop("p and q must be positive in multpoly()")

   z <- .Fortran("mtpoly",
                 as.double(alpha),
                 as.double(beta),
                 as.integer(p),
                 as.integer(q),
           gamma=as.double(gamma <- rep(0,p+q)))

   return(z$gamma)
}
odot <- function(x,y)
#--------------------------------------------------------------
#
#   Function to find the outer product of two vectors.
#
#--------------------------------------------------------------
{
   n <- length(x)
   m <- length(y)
   if(n != m) stop("Vectors not of same length in odot()")

   return(matrix(x,n,1) %*% matrix(y,1,n))
}
pacf <- function(x,m)
#------------------------------------------------------------------
#
#   Function to find sample partial autocorrelations.
#
#------------------------------------------------------------------
{
   n <- length(x)

   z <- .Fortran("pacf", 
                  as.double(x),  
                  as.integer(n),
                  as.integer(m),
                  as.double(e <- rep(0,n+m)),
                  as.double(f <- rep(0,n+m)),
                  as.double(temp <- rep(0,n+m)),
            theta=as.double(theta <- rep(0,m)))

   return(z$theta)
}
partar <- function(theta)
#-------------------------------------------------------------------
#
#   Function to find AR coefficients from partial autocorrelations.
#
#-------------------------------------------------------------------
{
   p <- length(theta)
   if(p==0) return()

   z <- .Fortran("partar", 
           theta=as.double(theta), 
                 as.integer(p))

   return(z$theta)
}
perdgm <- function(x)
#------------------------------------------------------------------------
#
#   Function to find periodogram of x.
#
#   result: periodogram at frequencies j / n, j = 0, 1, ..., [n / 2],
#           n = length(x)
#
#-----------------------------------------------------------------------
{
   ft <- Mod(fft(x - mean(x)))^2 / (n <- length(x))
   ft <- ft[1:((n / 2) + 1)]
}
plotsp <- function(f,n,div,main="Log Std Spectra")
#----------------------------------------------------------------------
#
#   Function to plot the log of the standardized version of f.
#   The array f is divided by div.
#
#-----------------------------------------------------------------------
{
   plot(freqs(n),log(stdf(f,div,exp(-6),exp(6))),xlab="frequency",
        ylab="value",ylim=c(-6,6),main,type="b")
}
tlpoly <- function(coeffs,x)
#-----------------------------------------------------------------------
#
#   Function to evaluate a polynomial having coefficients coeffs 
#   (arranged from degree 0, 1, etc) at the values in the array x
#
#-----------------------------------------------------------------------
{
   deg <- length(coeffs)-1
   n   <- length(x)

   if(deg==0) return(rep(coeffs[1],length(x)))

   z  <- .Fortran("poly",
                  as.double(coeffs),
                  as.integer(deg),
                  as.double(x),
                  as.integer(n),
               fx=as.double(rep(0,n)))

   return(z$fx)
}
polyrt <- function(coeffs,m=100,eps=1.e-6)
#---------------------------------------------------------------------
#
#   Function to find roots of polynomial 1 + sum(j=1,p) a(j)z^j
#           given its coefficients
#
#----------------------------------------------------------------------
{
   p   <- length(coeffs)
   ier <- 0
   z <- .Fortran("polyrt", 
                 as.double(coeffs), 
                 as.integer(p),
                 as.integer(m),
                 as.double(eps),
                 as.double(matrix(0,2,p+1)),
                 as.double(matrix(0,2,p+1)),
           roots=as.double(matrix(0,2,p)),
             ier=as.integer(ier))
           roots <- matrix(z$roots,2,p)
   list(roots=complex(real=roots[1,1:p],imag=roots[2,1:p]),ier=z$ier)
}
rtpoly <- function(roots)
#---------------------------------------------------------------------
#
#   Function to find coefficients of 1 + sum(j=1,p) a(j)z^j from its roots.
#
#----------------------------------------------------------------------
{
   p  <- length(roots)
   z <- .Fortran("rtpoly", 
                 as.double(rbind(Re(roots),Im(roots))),
                 as.integer(p),
                 as.double(matrix(0,2,p)),
           coeff=as.double(rep(0,p)))
   z$coeff
}

rw <- function(n,seed=0)
#----------------------------------------------------------------------
#
#   Function to simulate data from a Gaussian random walk.
#
#----------------------------------------------------------------------
{
   if(seed!=0) set.seed(seed)
   x <- cumsum(rnorm(n))
}
schur <- function(alpha)
#------------------------------------------------------------------
#
#   Function to form the schur matrix for AR coefficients.
#
#-------------------------------------------------------------------
{
   p <- length(alpha)
   if(p == 0) stop("Length of alpha must be positive in schur()")

   z <- .Fortran("schur",
                 as.double(alpha),
                 as.integer(p),
                 as.integer(p),
              SA=as.double(A <- matrix(rep(0,p*p),p,p)))

   return(matrix(z$SA,p,p))
}
seasest <- function(y,ords,coeffs,lags,back,maxit=50,eps=.000001)
#------------------------------------------------------------------
#
#
#
#------------------------------------------------------------------
{
   n    <- length(y)
   nn   <- sum(ords)
   ndim <- nn+1

   z <- .Fortran("marq",
                 as.double(y),
                 as.integer(n),
                 as.integer(ndim),
                 as.integer(maxit),
                 as.double(eps),
                 as.integer(ords),
          coeffs=as.double(coeffs),
                 as.integer(back),
                 as.integer(n+back),
                 as.integer(nn),
               e=as.double(rep(0,n+back)),
                 as.double(rep(0,n+back)),
                 as.double(rep(0,back+1)),
                 as.double(rep(0,(n+back)*nn)),
                 as.double(rep(0,ndim*ndim)),
                 as.double(rep(0,ndim*ndim)),
                 as.integer(lags),
             ier=as.integer(ier <- 0),
              rv=as.double(rv <- 0),
              se=as.double(se <- rep(0,ndim)))

   return(list(coeffs=z$coeffs,e=z$e[(back+1):(back+n)],ier=z$ier,
               rv=z$rv,se=z$se[1:nn]))
}





seaspred <- function(x,ords,coeffs,lags,rvar,tf,tl,hl,conf)
#-----------------------------------------------------------------
#
#
#-----------------------------------------------------------------
{

   n   <- length(x)
   nn  <- sum(ords[1:5])
   nd  <- ords[6]
   ndd <- ords[7]
   ns  <- ords[8]
   if(nd < 0 || ndd < 0 || ns < 0) stop("Error in d, D, or S")
   if(length(coeffs) < nn) stop("Error in coeffs")
   npl <- ords[2]
   nql <- ords[4]
   maxp <- ords[1]+nd+ndd*ns
   maxq <- ords[3]
   if(npl > 0) maxp <- maxp+lags[ords[2]]
   if(nql > 0) maxq <- maxq+lags[ords[2]+ords[4]]
   conf <- qnorm((1+conf)/2)

   ncf  <- maxp+maxq+1
   ncf1 <- (tl-tf+1)*hl

   z <- .Fortran("sspr",
                 as.double(x),
                 as.integer(n),
                 as.integer(ords),
                 as.double(coeffs),
                 as.integer(lags),
                 as.double(conf),
                 as.integer(tf),
                 as.integer(tl),
                 as.integer(hl),
               e=as.double(rep(0,n)),
              xx=as.double(rep(0,n+hl)),
              xp=as.double(rep(0,ncf1)),
             xpl=as.double(rep(0,ncf1)),
             xpu=as.double(rep(0,ncf1)),
             ier=as.integer(ier <- 0),
                 as.integer(npds <- 0),
                 as.double(rvar))


   return(xp=z$xp,xpl=z$xpl,xpu=z$xpu,ier=z$ier,check=z$e[1:hl],xx=z$xx[1:hl])
}






simcfs <- function(p,seed=0)
#-----------------------------------------------------------------------
#
#   Function to simulate coefficients of a polynomial 
#   having all of its zeros greater than one in modulus.
#
#-----------------------------------------------------------------------
{
   if(p==0) return()
   if(seed!=0) set.seed(seed)
   alpha <- partar(2 * (runif(p) - .5))
}
stdf <- function(f,fac,a,b)
#--------------------------------------------------------------------
#
#   Function to divide f by fac and then clip between a and b.
#
#--------------------------------------------------------------------
{
   f   <- f / fac
   f1  <- c(1:length(f))
   f   <- replace(f, f1[f > b], b)
   f   <- replace(f, f1[f < a], a)
}
swp <- function(a, k1, k2)
#---------------------------------------------------------------------
#
#   Function to sweep the matrix a on diagonals k1, ... , k2.
#
#----------------------------------------------------------------------
{
   n <- ncol(a)

   if( k1 < 1 || k2 > n || k1 > k2) stop("Illegal k1 or k2 in sweep()")

   z <- .Fortran("swpk12", 
               a1=as.matrix(a), 
                  as.integer(n),
                  as.integer(n), 
                  as.integer(k1),
                  as.integer(k2), 
              ier=as.integer(ier <- 0))

   list(A = z$a1, ier = z$ier)
}
toepl <- function(R,R0,M)
#-------------------------------------------------------------------
#
#   Function to form TOEPL(R0,...,R(M-1)).
#
#-------------------------------------------------------------------
{   
   G <- matrix(0,M,M)
   r <- c(R0, R[1:(M-1)])
   G <- matrix(r[abs(col(G) - row(G)) + 1], M, M)
}
tsvar <- function(x)
#------------------------------------------------------------------------
#
#   Function to find sample variance of a time series.
#
# x<-rnorm(100)
# tsvar(x)  
#------------------------------------------------------------------------
{
   dot(x - mean(x), x - mean(x)) / length(x)
}
windowf <- function(rho, R0, Q, ioptw, M, n, alpha=.05)
#---------------------------------------------------------------------
#
#   Function to find  nonparametric spectral estimator 
#   at Q frequencies in [0,1]
#   using window number ioptw, scale parameter M. 
#   The sample size n is also input. 
#
#   The spectral estimator and the constant c used in confidence
#   intervals are returned in a list.
#
#---------------------------------------------------------------------
{

   if(M >= n || Q < n || ioptw < 1 || ioptw > 8 
      || (ioptw <= 5 && length(rho) < M)
      || (ioptw > 5 && length(rho) < n - 1)
      || alpha <= 0. || alpha >= 1.)
   stop("Illegal Input to window()")

   z          <- rep(0, Q)
   if(ioptw <= 5) z[1:(M+1)] <- c(R0, R0 * rho[1:M])
   if(ioptw > 5)  z[1:n] <- c(R0, R0 * rho[1:(n-1)])

   z1 <- rep(0,Q)

   if(ioptw == 1) z1[1:(M+1)] <- rep(1, M+1)
      
   else if(ioptw == 2)  z1[1:(M+1)] <- 1 - (c(0:M) / M)

   else if(ioptw == 3)  z1[1:(M+1)] <- .54+.46*cos(pi*c(0:M)/M)

   else if(ioptw == 4) {
      u  <- c(1:M)/M
      u1 <- u[u<=.5]
      u2 <- u[u>.5]
      z1[1:(M+1)]  <- c(1,1-6*u1^2+6*u1^3,2*(1-u2)^3)
   }

   else if(ioptw == 5) {
      u   <- c(0:M) / M
      piu <- pi * u
      z1[1:(M+1)]  <- (1 - u) * cos(piu) + sin(piu) / pi
   }
 
   else if(ioptw == 6) {
      piu <- pi * (c(1:(n-1)) / M)
      z1[1:n] <- c(1,sin(piu)/piu)
   }

   else if(ioptw == 7) {
      piu <- pi * (c(1:(n-1)) / M)
      z1[1:n] <- c(1,3*((sin(piu)/piu) - cos(piu)) / (piu*piu))
   }

   else if(ioptw == 8) {
      u <- c(0:(n-1))/M
      z1[1:n] <- 1 /(1+u^4)
   }

   z <- c(z1[1] * z[1], 2 * z1[2:Q] * z[2:Q])

   z <- Re(fft(z))
   z <- z[1:((Q / 2) + 1)]

   ilam2 <- c(2,2./3.,.795,.539,.586,1,1.2,1.66)
   fac <- qnorm(1-(alpha/2)) * sqrt(M*ilam2[ioptw]/n)

   return(list(f=z,c=fac))

}





wn <- function(seed,n,dist=1)
#------------------------------------------------------------------------
#
#   Function to simulate white noise.
#
#   Generate a white noise time series of length n having distribution
#   specified by dist (see page 5 of the TIMESLAB text). seed is
#   a nonnegative integer seed for the random number generator (if
#   seed=0, then the output seed from last call is used).
#
# x<-wn(0,1000,2)
# plot(x,type="l")
#
#------------------------------------------------------------------------
{
   if(seed!=0) set.seed(seed)

   if(dist==1) return(rnorm(n))

   u <- runif(n)

   if(dist==2) return(u)

   if(dist==3) return(-log(1-u))

   if(dist==4) return(log(u/(1-u)))

   if(dist==5) return(tan(pi*(u-.5)))

   if(dist==6) return(log(-log(1-u)))

   if(dist==7) return(exp(rnorm(n)))

   if(dist==8) { xx <- (1:n)[u>.5]; u[xx] <- 1-u[xx]; u <- log(2*u); 
                 u[xx] <- -u[xx]; return(u);}

}
wntest <- function(x,m,alpha=.05)
#----------------------------------------------------------------------
#
#   Function to form the two plots for testing white noise.
#
# x<-wn(0,1000,2)
# wntest(x,20) 
#----------------------------------------------------------------------
{
   par(mfrow=c(2,1))
   n <- length(x)
   rho <- acf(x,m)$corr

   plot(rho,xlab="lag",ylab="value",main="Correlogram",ylim=c(-1,1))

   lim <- qnorm((1+(1-alpha)^(1./m))/2)/sqrt(n)

   lines(c(0,m),c(lim,lim))
   lines(c(0,m),c(-lim,-lim))
   result <- 0
   if(max(abs(rho)) > lim) result <-1

   cper <- perdgm(x)
   cper <- cumsum(cper)/sum(cper)
   plot(freqs(n),cper,xlab="frequency",ylab="value",
        main="Cumulative Periodogram")

   sq1 <- 1.36/sqrt(length(cper))
   sq2 <- 1-sq1
   sq3 <- .5*sq2
   sq4 <- .5*sq1

   lines(c(0,.5),c(0,1))
   lines(c(0,sq3),c(sq1,1))
   lines(c(sq4,.5),c(0,sq2))

   return(result)
}
