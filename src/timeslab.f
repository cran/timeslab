c**********************************************************************
      subroutine arfilt(alpha,p,rvar,x,n,m,w,sigma,a,ier)
c
c   Subroutine to multiply the (n x m) matrix  x by the inverse square
c   root of the (n x n) covariance matrix of an AR(p,alpha,rvar) process.
c
c***********************************************************************
c
      integer p, n, m
      double precision alpha(p),x(n,m),w(n,m),sigma(p),a(p,p)
      double precision al, al1, sr0, c, rvar
c
c   Find coefficients and residual variances for orders 1, ... , p:
c
      sigma(p) = rvar
      do 10 i = 1, p
 10   a(p,i) = alpha(i)
c
      if(p .gt. 1) then
         do 20 i = p - 1, 1, -1
         ip1 = i + 1
         ier = ip1
         al = a(i+1,i+1)
         if(abs(al) .ge. 1.) return
         al1 = 1.d0 - al * al
         sigma(i) = sigma(ip1) / al1
            do 15 j = 1, i
 15         a(i,j) = (a(ip1,j) - al * a(ip1, ip1-j)) / al1
c
 20      continue
c
      endif
c
      ier = 1 
      if(abs(a(1,1)) .ge. 1.) return
      sr0 = sqrt(sigma(1) / (1. - a(1,1) * a(1,1)))
      do 25 i = 1, p
 25   sigma(i) = sqrt(sigma(i))
      ier = 0
c
c   Multiply:
c
      do 100 k = 1, m
      w(1,k) = x(1,k) / sr0
         do 120 i = 2, n
         l = min0(i - 1, p)
c
            c = x(i,k)
            do 110 j = 1, l
 110        c = c + a(l,j)*x(i-j,k)
c
 120     w(i,k) = c / sigma(l)
 100  continue
c
c
      return
      end
c
c********************************************************************
c
	subroutine mxpd(x,alpha,beta,np,nq,rvar,ntf,ntl,nhf,nhl,
     1	nr,st,w,a,sh,st1,st2,xp,std)
c
c********************************************************************
c
        integer nr, np, nq, ntf, ntl, nhf, nhl
	double precision x(1),alpha(1),beta(1)
        double precision rvar
     	double precision xt(100),st(nr,nr),w(nr,nr),a(nr,nr)
        double precision xh(100)
      	double precision sh(nr,nr),xt1(100),st1(nr,nr),xp(1)
        double precision std(1),xt2(100)
      	double precision st2(nr,nr)
        double precision c1,c2,c,st11,z
c
c   initial conditions:
c
        iptstd=1
	mxpq=max0(np,nq)
	mxpqp1=mxpq+1
	if(np.gt.0) call vneg(alpha,np)
    	if(nq.gt.0) call vneg(beta,nq)
c
c   Form Initial State Vector (xh) and Covariance Matrix (sh):
c
	call wilson(alpha,np,beta,nq,xt,mxpqp1,xt1,mxpqp1,
     1	xt2,mxpq,ier)
	if(ier.eq.1) go to 99
	if(np.gt.0) call vneg(alpha,np)
    	if(nq.gt.0) call vneg(beta,nq)
c
c
	do 50 i=1,nr
	xh(i)=0.
	do 50 j=1,i
    	sh(i,j)=xt(i-j+1)*rvar
  50	sh(j,i)=sh(i,j)
	if(nr.gt.1) then
		call mxma(alpha,beta,np,nq,nr,xt1)
		do 51 j=2,nr
		do 51 k=j,nr
		c=sh(j,k)
			do 52 l=0,j-2
			c1=1.
			if(l.gt.0) c1=xt1(l)
			c2=1.D0
			l1=l+k-j
			if(l1.gt.0) c2=xt1(l1)
  52			c=c-c1*c2*rvar
		sh(j,k)=c
  51		sh(k,j)=c
	endif
c
c   form a and w:
c
	xt1(1)=1.
	if(nr.gt.1) then
		do 55 i=2,nr
		c=0.0
		if(i-1.le.nq) c=beta(i-1)
		if(np.gt.0) then
			do 54 j=1,min0(np,i-1)
  54			c=c-alpha(j)*xt1(i-j)
		endif
  55		xt1(i)=c
	endif
	do 56 i=1,nr
	do 56 j=1,nr
  56	w(i,j)=rvar*xt1(i)*xt1(j)
	do 57 i=1,nr
	do 57 j=1,nr
  57	a(i,j)=0.0
	if(np.gt.0) then
		do 58 i=1,np
  58		a(nr,nr-i+1)=-alpha(i)
	endif
	if(nr.gt.1) then
		do 59 i=1,nr-1
  59		a(i,i+1)=1.
	endif
c
c   iterate:
c
	npr=0
	nr2=nr*nr
	nr24=4*nr2
	nr4=4*nr
c
c
	do 100 it=0,ntl
c
c  x tilda and sigma tilda:
c
	call mnvar(xh,a,sh,w,nr,nr,xt,st)
c
c   Forecasts:
c
	if(it.ge.ntf.and.it.le.ntl) then
	call movxyr(xt1,xt,nr)
	call movxyr(st1,st,nr2)
 	do 60 ih=1,nhl
		if(ih.ge.nhf.and.ih.le.nhl) then
			npr=npr+1
			xp(npr)=xt1(1)
			if(iptstd.eq.1) std(npr)=sqrt(st1(1,1))
		endif
		if(ih.eq.nhl) go to 60
c
c
		call mnvar(xt1,a,st1,w,nr,nr,xt2,st2)
		call movxyr(xt1,xt2,nr)
		call movxyr(st1,st2,nr2)
  60	continue
	endif
	if(it.eq.ntl) go to 100
c
c   Update x hat and sigma hat:
c
	st11=st(1,1)
	do 65 i=1,nr
  65	xt1(i)=st(1,i)/st11
c
c   xhat:
c
	z=x(it+1)-xt(1)
	do 70 i=1,nr
  70	xh(i)=xt(i)+xt1(i)*z
c
c  sigma hat:
c
        do 68 i=1,nr
           do 68 j=1,nr
 68        sh(i,j)=st(i,j)-xt1(i)*st(1,j)
c	do 66 i=1,nr
c	do 66 j=1,nr
c  66	st1(i,j)=0.
c	do 67 i=1,nr
c    	st1(i,i)=1.
c  67	st1(i,1)=st1(i,1)-xt1(i)
c	do 68 i=1,nr
c	do 68 j=1,nr
c	c=0.
c		do 69 k=1,nr
c  69		c=c+st1(i,k)*st(k,j)
c  68	sh(i,j)=c
c
c
 100	continue
c
c
c
  99	continue
	return
	end
c
c	
c*******************************************************************
c
        subroutine mnvar(x,a,sig,w,ndim,n,y,b)
c
c   Subroutine to calculate
c
c   y=a*x     and     b=a*sig*a' + w
c
c   for nx1 x and nxn a and sig  in Kalman Filter for ARMAPRED
c
c*******************************************************************
c
        integer n, ndim
	double precision x(n),a(ndim,ndim),sig(ndim,ndim)
        double precision b(ndim,ndim),y(n)
	double precision w(ndim,ndim)
        double precision c,c1
c
c
	do 20 i=1,n
		if(i.lt.n) then
			c=x(i+1)
			go to 20
		endif
c
		c=0.0
		do 10 j=1,n
  10		c=c+a(i,j)*x(j)
  20	y(i)=c
c
c
	if(n.eq.1) go to 40
	do 30 i=1,n-1
	do 30 j=1,i
	c=w(i,j)+sig(i+1,j+1)
	b(i,j)=c
  30	b(j,i)=c
c
c   i=n:
c
  40	continue
	do 70 j=1,n
		c=w(n,j)
		do 60 k=1,n
	 	if(j.lt.n) then
			c1=sig(k,j+1)
			go to 60
	 	endif
			c1=0.
			do 50 l=1,n
  50			c1=c1+sig(k,l)*a(j,l)
  60		c=c+a(n,k)*c1
  	b(n,j)=c
  70	b(j,n)=c
c
c
	return
	end
c
c	
c********************************************************************
c
        subroutine vneg(x,n)
c
c   Subroutine to negate a vector
c
c*********************************************************************
c
        integer n
        double precision x(n)
c
	do 10 i=1,n
  10	x(i)=-x(i)
	return
	end
c
c
c********************************************************************
c
	subroutine mxma(alpha,beta,np,nq,n,d)
c
c   Subroutine to find the first n coefficients of the MA(infinity)
c   representation of an ARMA(alpha,beta,np,nq) process.
c
c********************************************************************
c
        integer n
	double precision alpha(1),beta(1),d(n)
        double precision c,c1
c
	do 10 i=1,n
  10	d(i)=0.0
	if(nq.gt.0) then
		do 20 i=1,nq
  20		d(i)=beta(i)
	endif
c
c
	if(np.gt.0) then
		do 40 i=1,n
		c=d(i)
		if(np.eq.0) go to 40
			do 30 j=1,min0(i,np)
			c1=1.d0
			if(i-j.gt.0) c1=d(i-j)
  30			c=c-alpha(j)*c1
  40		d(i)=c
	endif
c
c
	return
	end
c
c***********************************************************************
c
        subroutine movxyr(x,y,n)
c
c************************************************************************
c
        integer n
        double precision x(1),y(1)
c
        do 10 i=1,n
 10     x(i)=y(i)
c
        return
        end
c
c************************************************************************
c
        subroutine mtpoly(alpha,beta,p,q,gamma)
c
c************************************************************************
c
        integer p,q
        double precision alpha(1),beta(1),gamma(1)
c
c
        do 10 i=1,p+q
 10     gamma(i)=0.0
        do 20 i=1,p
 20     gamma(i)=gamma(i)+alpha(i)
        do 30 i=1,q
 30     gamma(i)=gamma(i)+beta(i)
        do 40 i=1,p
        do 40 j=1,q
 40     gamma(i+j)=gamma(i+j)+alpha(i)*beta(j)
c
c
        return
        end
c
c******************************************************************
c
        subroutine arpart(alpha,p,ier)
c
c   Subroutine to find partials from AR coefficients.
c
c******************************************************************
c
        integer p,ier
        double precision alpha(p)
        double precision temp,c,c2
c
        do 10 j=1,p
 10     alpha(j)=-alpha(j)
        ier=p
        if(abs(alpha(p)).ge.1.) return
        ier=0
        if(p.eq.1) return
c
c
        do 30 j=p-1,1,-1
           ier=j
c
           c=alpha(j+1)
           c2=1.d0-c*c
c
           do 20 k=1,(j+1)/2
              temp=(alpha(k)+c*alpha(j+1-k))/c2
              alpha(j+1-k)=(alpha(j+1-k)+c*alpha(k))/c2
 20           alpha(k)=temp
c
           if(abs(alpha(j)).ge.1.) return
 30     continue
c
        ier=0
c
c
        return
        end
c
c*************************************************************************
c
        subroutine arsppk(alpha,p,rvar,start,n,f,covri,ri,a,b,c,
     1                    ier,omega,se)
c
c*************************************************************************
c
        integer p,n,ier,oi
        double precision alpha(p),rvar,start,f
        double precision covri(p,p),ri(p),a(p,p),b(p),c(p,p)
        double precision omega,se,ri0,omega1
        double precision fm,twopi,fp,fpp,var,del
c
c
c   if start.eq.0 find largest relative max in f:
c
        if(start.eq.0.0) then
           fm=0.d0
           do 20 i=2,128
              if(f(i-1).lt.f(i).and.f(i).gt.f(i+1)) then
                 if(f(i).gt.fm) then
                    start=dble(i-1)/256.
                    fm=f(i)
                 endif
              endif
 20        continue
        endif
c
c   error return if no relative max's
c
	if(start.eq.0.0) then
           ier=1
           go to 99
        endif
c
c   Find Peak Frequency:
c
	call macv(alpha,p,1.d0,ri,ri0)
	twopi=8.*atan(1.0)
	it=1
	omega=start
c
c   start iterations:
c
 35        fp=0.d0
           fpp=0.d0
           do 40 i=1,p
              oi=i
              fp=fp+oi*ri(i)*sin(twopi*oi*omega)
 40        fpp=fpp+oi*oi*ri(i)*cos(twopi*oi*omega)
c
c   check for 0 second derivative:
c
           if(abs(fpp).lt.1.e-20) then
              ier=2
              go to 99
           endif
c
c   check for convergence:
c
           del=fp/(fpp*twopi)
           omega1=omega-del
           if(abs(del/omega).lt.1.e-5) then
c
c   yes:
c
c
c   check for convergence to 0 or .5:
c
	      if(omega1.lt.1.e-5.or.omega1.gt..49999) then
                   ier=3
                   go to 99
              endif
	      start=omega1
	      go to 60
           endif
c
c   no:
c
	   it=it+1
	   omega=omega1
c
c   check for max # of iterations:
c
	   if(it.gt.100) then
	      ier=4
	      go to 99
	   endif
	   go to 35
c
c   Asymptotic Standard Error:
c
  60	continue

        call freqcv(start,p,alpha,n,p,covri,ri,a,c,b,var)
        omega=start
	se=sqrt(var)
        ier=0
c
  99	continue
	return
	end
c
c*****************************************************************
c
      subroutine freqcv(freq,np,alpha,nobs,ndim,covri,
     1ri,a,c,b,var)
c
c   subroutine to find the asymptotic variance of a frequency
c   estimator freq for an ar(np) process having coeffs alpha.
c
c   input :
c            freq,np,alpha(1),alpha(np)
c            nobs: number of observations
c            ndim : row dimension of various arrays in calling prog
c
c   output :
c            var : asymptotic variance
c            covri,ri
c
c********************************************************************
c
      integer np, ndim
      double precision freq
      double precision alpha(np),covri(ndim,ndim),ri(np)
      double precision a(ndim,ndim), c(ndim,ndim),b(np)
      double precision cc,var,c1,c2,twopi
c
c
      call schur(alpha,ndim,np,a)
      twopi=8.*atan(1.0)
      do 10 iv=1,np
      do 10 ij=1,np
      cc=0.d0
      ivpij=iv+ij
      ijmiv=ij-iv
      if(ivpij.le.np) cc=cc+alpha(ivpij)
      if(ijmiv.eq.0) cc=cc+1.d0
      if(ijmiv.gt.0) cc=cc+alpha(ijmiv)
  10  c(iv,ij)=cc
      do 40 j=1,np
      do 40 k=1,j
         cc=0.0
         do 30 l=1,np
         c1=c(j,l)
         c2=0.0
            do 20 m=1,np
  20        c2=c2+a(l,m)*c(k,m)
  30     cc=cc+c1*c2
      covri(j,k)=cc
  40  covri(k,j)=cc
      do 60 iv=1,np
      cc=alpha(iv)
      if(iv.eq.np) go to 60
         npmiv=np-iv
         do 50 j=1,npmiv
  50     cc=cc+alpha(j)*alpha(j+iv)
  60  ri(iv)=cc
c
c   find vector b :
c
      do 90 iv=1,np
  90  b(iv)=dble(iv)*sin(iv*freq*twopi)
c
c   find var :
c
c
c   find -h''/4pi
c
	cc=0.d0
	do 95 i=1,np
  95	cc=cc+dble(i*i)*ri(i)*cos(twopi*dble(i)*freq)
	cc=twopi*cc
      var=0.0
      do 100 j=1,np
      do 100 k=1,np
 100  var=var+b(j)*covri(j,k)*b(k)
      var=var/(cc*cc)
	var=var/dble(nobs)

c
c
      return
      end
c
c**********************************************************
c
      subroutine mxcsd(alph,beta,np,nq,n,mdim,wk1,r,
     1 ainf,sd,ier)
c
c   subroutine to find the standard deviation of the estimated
c   coefficients (based on n observations) of an arma(np,nq,alph,beta)
c   process.
c
c   input :
c           np,nq,alph,beta
c           mdim : dimension of various arrays in
c           calling program (see dimensions below)
c           mdim.ge.(2*max(np,nq))+1
c
c   output :
c           ainf
c           ier : 0 is normal return, 1 means thereis
c           a singular matrix in the procedure
c           sd(1),...,sd(np+nq)
c
c   subroutines called : mxcvin,schur,swpk12,decomp,solv
c
c**********************************************************
c
      integer np, nq, mdim, n, ier
      double precision alph(np),beta(nq),wk1(mdim,mdim)
      double precision ainf(mdim,mdim),r(1),sd(1),on
c
c   iab :
c
	if(np*nq.eq.0) go to 99
	ier=0
	on=n
      m=max0(np,nq)
      mm=2*m+1
      nppnq=np+nq
      call mxcvin(alph,beta,np,nq,mm,mdim,ainf,r,ier)
	if(ier.eq.1) go to 99
      do 1 i=1,nppnq
      do 1 j=1,nppnq
   1  ainf(i,j)=0.
      do 2 i=1,np
      do 2 j=1,nq
      ainf(i,np+j)=-r(i-j+m+1)
   2  continue
c
c   iaa :
c
      call schur(alph,mdim,np,wk1)
	call swpk12(wk1,mdim,np,1,np,ier)
	if(ier.eq.1) go to 99
      do 3 i=1,np
      do 3 j=i,np
   3  ainf(i,j)=wk1(i,j)
c
c   ibb :
c
      call schur(beta,mdim,nq,wk1)
	call swpk12(wk1,mdim,nq,1,nq,ier)
	if(ier.eq.1) go to 99
      do 4 i=1,nq
      ii=np+i
      do 4 j=i,nq
   4  ainf(ii,np+j)=wk1(i,j)
c
      do 5 i=1,nppnq
      do 5 j=i,nppnq
   5  ainf(j,i)=ainf(i,j)
	call swpk12(ainf,mdim,nppnq,1,nppnq,ier)
	on=n
	do 10 i=1,nppnq
  10	sd(i)=sqrt(ainf(i,i)/on)
  99	continue
c
      return
      end
c

c**********************************************************
c
      subroutine mxcvin(alph,beta,np,nq,mm,ndim,a,r,ier)
c
c   subroutine to calculate the cross covariances for
c   lags (-max(np,nq),...,max(np,nq)), (stored in
c   r(1),...,r(mm),mm=(2*max(np,nq))+1) of the two
c   dimensional process (x1(t),x2(t)) where x1(.) and
c   x2(.) are defined by :
c
c      sum(j=0,np) alph(j)*x1(t-j)=e(t)
c      sum(j=0,nq)beta(j)*x2(t-j)=e(t)
c
c   input :
c           np,nq,mm=2*max(np,nq)+1
c           alph,beta
c           ndim : dimension of a in calling program
c           (ndim.ge.mm)
c
c   output :
c           r(1),...,r(mm) (r(j) is r12(j-max(np,nq)-1)
c           ier : 0 is normal return, 1 means a is singular
c
c   subroutines called : decomp,solv
c
c**********************************************************
c
      integer np, nq, ndim, mm, ip(100)
      double precision alph(np),beta(nq),a(ndim,ndim),r(mm)
c      data nout/0/
c
c
c
      m=max0(np,nq)
      mp1=m+1
      mp2=m+2
      do 1 i=1,mm
      do 1 j=1,mm
   1  a(i,j)=0.
c
c
      if(np.lt.nq) go to 50
c
c   for np.ge.nq :
c
      a(1,1)=1.
      do 5 i=1,np
   5  a(1,i+1)=alph(i)
      do 6 j=2,mp1
      jm1=j-1
      a(j,1)=a(jm1,mm)
         do 7 k=2,mm
   7     a(j,k)=a(jm1,k-1)
   6  continue
      m1=m-nq+1
      do 8 j=1,nq
   8  a(mp2,m1+j)=beta(nq-j+1)
      a(mp2,m1+nq+1)=1.
      mp3=m+3
      do 9 j=mp3,mm
      jm1=j-1
      a(j,1)=a(jm1,mm)
         do 10 k=2,mm
  10     a(j,k)=a(jm1,k-1)
   9  continue
      go to 100
c
c   for np.lt.nq :
c
  50  a(1,1)=1.
      do 51 i=1,np
  51  a(1,i+1)=alph(i)
      if(m.eq.1) go to 70
      do 52 j=2,m
      jm1=j-1
      a(j,1)=a(jm1,mm)
         do 53 k=2,mm
  53     a(j,k)=a(jm1,k-1)
  52  continue
  70  do 71 j=1,nq
  71  a(mp1,j)=beta(nq-j+1)
      a(mp1,nq+1)=1.
      do 72 j=mp2,mm
      jm1=j-1
      a(j,1)=a(jm1,mm)
         do 73 k=2,mm
  73     a(j,k)=a(jm1,k-1)
  72  continue
c
c   solve system of equations :
c
 100  do 101 i=1,mm
 101  r(i)=0.
      r(mp1)=1.
c
c
      call decomp(mm,ndim,a,ip)
      if(ip(mm).eq.0) go to 125
      ier=0
      call solv(mm,ndim,a,r,ip)
c
c   r contains r12(-m),...,r12(m)
c
      return
c
 125  ier=1
      return
      end


c************************************************************************
c
      subroutine corrar(rho,r0,p,alpha,rvar)
c
c************************************************************************
c
      integer p
      double precision rho(p),alpha(p)
      double precision rvar, r0, last, temp
c
c
      alpha(1)=-rho(1)
      rvar=1-rho(1)*rho(1)
      if(p.eq.1) go to 99
c
      do 30 j=2,p
c
         last=-rho(j)
         do 10 k=1,j-1
 10      last=last-alpha(k)*rho(j-k)
         last=last/rvar
c
         alpha(j)=last
         do 20 k=1,j/2
            temp=alpha(k)
            alpha(k)=alpha(k)+last*alpha(j-k)
            if(k.eq.(j-k)) go to 20
            alpha(j-k)=alpha(j-k)+last*temp
 20      continue
c
         rvar=rvar*(1-last*last)
c
 30   continue
c
 99   rvar=r0*rvar
c
c
      return
      end




c*******************************************************************
c
        subroutine cvmx1(r,r0,np,nq,ndim,maxit,del,ip,al,
     1  alpha,beta,rvar,ier)
c
c        subroutine cvmx1(r,r0,np,nq,ndim,maxit,del,ip,al,ra,ry,g,
c     1  alpha,beta,rvar,ier,t,f,fn,fp,d,dn,gam)
c   Subroutine to find arma paramters from covariances.
c
c   ndim.ge.np+nq
c
c******************************************************************
c
        integer ndim,nq,np,maxit,ip(ndim),ier
        double precision del,rvar
        double precision r(ndim),r0,al(ndim,ndim),ra(100),ra0
        double precision ry(100),alpha(ndim),ry0
        double precision beta(ndim),g(100)
        double precision t(100),f(100),fn(100),d(100)
        double precision dn(100),gam(100),fp(100)
        double precision c,c1

c        dimension r(ndim),al(ndim,ndim),ra(ndim),ry(ndim),alpha(ndim),
c     1  beta(ndim),ip(ndim),g(ndim),t(ndim),f(ndim),fn(ndim),d(ndim),
c     1	dn(ndim),gam(ndim),fp(ndim)
c
c   find alpha :
c
        call hiyw(r,r0,np,nq,ndim,ip,al,alpha,ier)
        if(ier.eq.1) then
		ier=3
		go to 99
	endif
c
c   find ra :
c
        call macv(alpha,np,1.d0,ra,ra0)
c
c   find ry :
c
        ry0=ra0*r0
        do 12 i=1,np
  12    ry0=ry0+2.*ra(i)*r(i)
        do 13 i=1,nq
                c=ra0*r(i)
                do 14 j=1,np
                if(i-j.eq.0) go to 16
                n1=iabs(i-j)
                c1=r(n1)
                go to 14
  16            c1=r0
  14            c=c+ra(j)*(c1+r(i+j))
  13    ry(i)=c
c
c   find corresponding betas :
c
	call cvmawl(ry,ry0,nq,maxit,del,t,f,fn,fp,g,d,dn,gam,
     1  beta,rvar,ier)
  99	continue
        return
        end
c

c******************************************************************
c
        subroutine hiyw(r,r0,np,nq,ndim,ip,al,alpha,ier)
c
c   Subroutine to solve high order Yule-Walker equations. ndim.ge.np+nq
c
c**********************************************************************
c
        integer ndim,np,nq,ier,ip(ndim)
        double precision r(ndim),r0
        double precision al(ndim,ndim),alpha(np)
c
        ier=0
        if(np.eq.1) go to 10
c
        do 1 i=1,np
        alpha(i)=-r(nq+i)
                do 1 j=i,np
                if(nq+i-j.eq.0) go to 3
                n1=iabs(nq+i-j)
                al(i,j)=r(n1)
                go to 1
   3            al(i,j)=r0
   1    al(j,i)=r(nq+j-i)
c
        call decomp(np,ndim,al,ip)
        if(ip(np).eq.0) go to 99
        call solv(np,ndim,al,alpha,ip)
        return
  10    alpha(1)=-r(nq+1)/r(nq)
        return
 99     ier=1
        return
        end
c

c******************************************************************
c
	subroutine cvmawl(r,r0,nq,maxit,del,t,f,fn,fp,g,d,dn,gam,
     1	beta,rvar,ier)
c
c   Subroutine to do Wilson's algorithm for MA processes.
c
c******************************************************************
c
        integer nq,maxit,ier
	double precision r(nq),r0, t(nq),f(nq),fn(nq),fp(nq)
        double precision fpk, fpk2, c,d0,e1,g0,t0,dn0, eps
        double precision g(nq),d(nq),dn(nq)
        double precision gam(nq),beta(nq),del,rvar
c
c   starting values:
c
	t0=sqrt(r0)
	do 10 j=1,nq
  10	t(j)=r(j)/t0
c
c   start iterations (t plays the role of tau, g plays role of r):
c
	do 150 it=1,maxit
c
c
		do 20 j=1,nq
  20		f(j)=t(j)/t0
		g0=t0*t0
		do 22 i=1,nq
  22		g0=g0+t(i)*t(i)
		do 24 i=1,nq
			c=t0*t(i)
			if(i.eq.nq) go to 26
			do 25 j=1,nq-i
  25			c=c+t(j)*t(j+i)
  26			g(i)=c
  24		continue
		g0=(g0+r0)/t0
		do 30 j=1,nq
  30		g(j)=(g(j)+r(j))/t0
c
c   work backward:
c
		do 70 k=nq,1,-1
		gam(k)=g(k)
		fp(k)=f(k)
		fpk=fp(k)
		fpk2=1.D0-fpk*fpk
		if(fpk2.le.0.) go to 299
		if(k.eq.1) go to 70
			do 40 j=1,k-1
  40			fn(j)=(f(j)-fpk*f(k-j))/fpk2
			do 50 j=1,k-1
  50			g(j)=g(j)-gam(k)*fn(k-j)
			do 60 j=1,k-1
  60			f(j)=fn(j)
  70		continue
c
c   work forward:
c
		d0=g0/2.
		do 110 k=0,nq-1
		d(k+1)=gam(k+1)
		fpk=fp(k+1)
		fpk2=1.D0-fpk*fpk
		dn0=(d0-fpk*d(k+1))/fpk2
		dn(k+1)=(d(k+1)-fpk*d0)/fpk2
		if(k.eq.0) go to 90
			do 80 j=1,k
  80			dn(j)=(d(j)-fpk*d(k+1-j))/fpk2
  90		d0=dn0
			do 100 j=1,k+1
 100			d(j)=dn(j)
 110		continue
c
c   check convergence:
c
		eps=abs(t0-d0)/abs(d0)
		do 120 j=1,nq
		e1=abs(t(j)-d(j))/abs(d(j))
		if(e1.gt.eps) eps=e1
 120		continue
		t0=d0
		do 130 j=1,nq
 130		t(j)=d(j)
c
c   yes:
c
		if(eps.lt.del) then
			ier=0
			rvar=t0*t0
			do 140 j=1,nq
 140			beta(j)=t(j)/t0
			go to 99
		endif
c
c   no:
c
 150	continue
c
c   nonconvergence in maxit iterations:
c
	ier=1
	go to 99
c
c   partial outside (-1,1):
c
 299    ier=2
c
c
  99	continue
	return
	end
c

c***********************************************************
c
      subroutine decomp(n,ndim,a,ip)
c
c   matrix triangularization by gaussian elimination.
c
c   input....
c     n = order of matrix
c     ndim = declared dimension of array a
c     a = matrix to be triangularized
c
c   output....
c     a(i,j), i.le.j = upper triangular factor, u.
c     a(i,j), i.gt.j = multipliers = lower triangular
c                      factor, i-l
c     ip(k),  k.lt.n = index of k-th pivot row
c     ip(n) = (-1)**(number of interchanges) or 0
c
c   use subroutine solve to obtain solution of linear system.
c   determ(a) = ip(n)*a(1,1)*a(2,2)*...*a(n,n)
c   if ip(n)=0, a is singular,  solve  will divide by zero
c   interchanges finished in u, only partly in l.
c
c   reference...algorithm 423 'linear equation solver'
c               by cleve b. moler, cacm, april 1972 p. 274
c
c***********************************************************
c
      integer ndim,ip(ndim)
      double precision a(ndim,ndim),t
c
      ip(n)=1
      do 6 k=1,n
         if(k.eq.n) go to 5
         kp1=k+1
         m=k
         do 1 i=kp1,n
            if(abs(a(i,k)).gt.abs(a(m,k))) m=i
 1       continue
         ip(k)=m
         if(m.ne.k) ip(n)=-ip(n)
         t=a(m,k)
         a(m,k)=a(k,k)
         a(k,k)=t
         if(t.eq.0.) go to 5
         do 2 i=kp1,n
 2          a(i,k)=-a(i,k)/t
         do 4 j=kp1,n
            t=a(m,j)
            a(m,j)=a(k,j)
            a(k,j)=t
            if(t.eq.0.) go to 4
            do 3 i=kp1,n
 3             a(i,j)=a(i,j)+a(i,k)*t
 4       continue
 5       if(a(k,k).eq.0.) ip(n)=0
 6    continue
      return
      end
c

c******************************************************************
c     
      subroutine solv(n,ndim,a,b,ip)
c     
c     solution of linear system  a*x=b
c     
c     do not use if  decomp  has set ip(n)=0
c     
c     input....
c     n = order of matrix a
c     ndim = declared dimension of array a
c     a = triangularized matrix obtained from subroutine |decomp|
c     b = right hand side vector
c     ip = pivot vector obtained from |decomp|
c     
c     output....
c     b = solution vector, x.
c     
c     reference...algorithm 423 'linear equation solver'
c     by cleve b. moler, cacm, april 1972 p. 274
c     
c******************************************************************
c     

      integer ndim,n,ip(ndim)
      double precision a(ndim,ndim),b(ndim),t

c     
      if(n.eq.1) go to 9
      nm1=n-1
      do 7 k=1,nm1
         kp1=k+1
         m=ip(k)
         t=b(m)
         b(m)=b(k)
         b(k)=t
         do 7 i=kp1,n
 7          b(i)=b(i)+a(i,k)*t
      do 8 kb=1,nm1
         km1=n-kb
         k=km1+1
         b(k)=b(k)/a(k,k)
         t=-b(k)
         do 8 i=1,km1
 8          b(i)=b(i)+a(i,k)*t
 9    b(1)=b(1)/a(1,1)
      return
      end
c     
c     
c*******************************************************
c     
      subroutine macv(beta,nq,sig,r,r0)
c     
c     subroutine to calculate the autocovariances r0,r(1),
c     ...,r(nq) for a moving average process of order nq
c     with parameters beta(1),...,beta(nq), and sig (res var)
c     
c     input :
c     nq,beta(1),...,beta(nq),sig
c     
c     output :
c     r0,r(1),...,r(nq)
c     
c     subroutines called : none
c     
c*******************************************************
c     
      integer nq
      double precision beta(nq),r(nq),r0,sig,c
c     
      c=1.D0
      do 1 i=1,nq
 1       c=c+beta(i)*beta(i)
      r0=c*sig
c     
      do 2 i=1,nq
         c=beta(i)
         if(i.eq.nq) go to 2
         nqmi=nq-i
         do 3 j=1,nqmi
 3          c=c+beta(j)*beta(j+i)
 2    r(i)=c*sig
c     
      return
      end
c
c**********************************************************************
c
      subroutine diffeq(alpha,np,n,e,x)
c
c**********************************************************************
c
      integer np,n
      double precision alpha(1),e(n),x(n),c
c
c
      if(np.eq.0) then
         do 5 i=1,n
 5          x(i)=e(i)
         return
      endif
c
      do 10 i=1,np
 10      x(i)=e(i)
      do 20 i=np+1,n
         c=e(i)
         do 30 j=1,np
 30         c=c-alpha(j)*x(i-j)
 20   x(i)=c
c
c
      return
      end
c
c*********************************************************************
c
      subroutine dtarma(start,eps,maxit,dat,ny,e,np,nq,
     1                  alpha,beta,nppnq,p,rv,am2ll,ier,var)
c
c
c   Subroutine to find exact ARMA MLE's as described in TIMESLAB.
c
c   Input:
c      start: starting values (those for alpha, then for beta)
c      eps  : convergence criterion
c      maxit: maximum number of function evaluations
c      dat  : data
c      ny   : number of data points
c      np,nq: order of ARMA
c      nppnq: np+nq
c
c   Output:
c      alpha, beta, rv: values of parameters when algorithm terminates
c      am2ll : -2log(likelihhod) at last values
c      ier   : <=0 means convergence, 1 means nonconvergence, >1 means 
c              error.
c      var   : variance of vertices upon termination
c
c   Auxilliary:
c      p(np+nq,np+nq+1) 
c      e(ny)
c
c   Used routines: armalk,partar,lkhood,wilson
c
c********************************************************************
c
        integer maxit,ny,np,nq,nppnq,ier
	double precision start(1),eps,min(20)
        double precision ynewlo,step(20)
      	double precision p(nppnq,1),pstar(20),p2star(20)
        double precision pbar(20),y(20)
        double precision dn,dnn,z,ylo,rcoeff,ystar,ecoeff
        double precision y2star,ccoeff,curmin,del,armalk
	double precision dat(ny),alpha(1),beta(1)
	double precision e(ny),vw(20),vk(20),vl(20)
	double precision sum,summ,cmin,am2ll,var,rv
c
	data rcoeff,ecoeff,ccoeff/1.,2.,.5/
c
c
        n=np+nq
        konvge=10
        do 310 i=1,n
 310    step(i)=1.0
	ier=0
	kcount=maxit
	maxit=0
	if(eps.le.0.) maxit=maxit-1
	if(n.gt.20) maxit=maxit-10
	if(konvge.le.0) maxit=maxit-100
	if(maxit.lt.0) return
c
	jcount=konvge
	dn=dble(n)
	nn=n+1
	dnn=dble(nn)
	del=1.
c
c   construction of initial simplex:
c
1001	do 1 i=1,n
   1	p(i,nn)=start(i)
	z=armalk(start,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	y(nn)=z
	sum=dble(z)
	summ=dble(z*z)
	do 2 j=1,n
	start(j)=start(j)+step(j)*del
		do 3 i=1,n
   3		p(i,j)=start(i)
	z=armalk(start,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	y(j)=z
	sum=sum+dble(z)
	summ=summ+dble(z)*dble(z)
   2	start(j)=start(j)-step(j)*del
c
c   simplex construction complete
c
c   find highest and lowest y values. ynewlo (=y(ihi)) indicates
c   the vertex of the simplex to be replaced.
c
1000	ylo=y(1)
	ynewlo=ylo
	ilo=1
	ihi=1
	do 5 i=2,nn
	if(y(i).ge.ylo) go to 4
	ylo=y(i)
	ilo=i
   4	if(y(i).le.ynewlo) go to 5
	ynewlo=y(i)
	ihi=i
   5	continue
	sum=sum-dble(ynewlo)
	summ=summ-dble(ynewlo)*dble(ynewlo)
c
c   calculate pbar. The centroid of the simplex vertices
c   excepting that with y value ynewlo.
c
	do 7 i=1,n
	z=0.
		do 6 j=1,nn
   6		z=z+p(i,j)
	z=z-p(i,ihi)
   7	pbar(i)=z/dn
c
c   reflection through the centroid
c
	do 8 i=1,n
   8	pstar(i)=(1.+rcoeff)*pbar(i)-rcoeff*p(i,ihi)
	ystar=armalk(pstar,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	maxit=maxit+1
	if(ystar.ge.ylo) go to 12
c
c   successful reflection, so extension
c
	do 9 i=1,n
   9	p2star(i)=ecoeff*pstar(i)+(1.-ecoeff)*pbar(i)
	y2star=armalk(p2star,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	maxit=maxit+1
c
c   retain extension or contraction
c
	if(y2star.ge.ylo) go to 19
  10	do 11 i=1,n
  11	p(i,ihi)=p2star(i)
	y(ihi)=y2star
	go to 900
c  no extension
  12	l=0
	do 13 i=1,nn
	if(y(i).gt.ystar) l=l+1
  13	continue
	if(l.gt.1) go to 19
	if(l.eq.0) go to 15
c
c  contraction on the reflection side of the centroid
c
	do 14 i=1,n
  14	p(i,ihi)=pstar(i)
	y(ihi)=ystar
c
c   contraction on the y(ihi) side of the centroid
c
  15	do 16 i=1,n
  16	p2star(i)=ccoeff*p(i,ihi)+(1.-ccoeff)*pbar(i)
	y2star=armalk(p2star,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	maxit=maxit+1
	if(y2star.le.y(ihi)) go to 10
c
c   contract whole simplex
c
	sum=0.d0
	summ=0.d0
	do 18 j=1,nn
		do 17 i=1,n
		p(i,j)=(p(i,j)+p(i,ilo))*0.5
  17		min(i)=p(i,j)
	y(j)=armalk(min,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	sum=sum+dble(y(j))
  18	summ=summ+dble(y(j))*dble(y(j))
	maxit=maxit+nn
	go to 901
c   retain reflection
  19	do 20 i=1,n
  20	p(i,ihi)=pstar(i)
	y(ihi)=ystar
 900	sum=sum+dble(y(ihi))
	summ=summ+dble(y(ihi))*dble(y(ihi))
 901	jcount=jcount-1
	if(jcount.ne.0) go to 1000
c
c   check to see if minimum reached
c
	if(maxit.gt.kcount) go to 22
	jcount=konvge
	cmin=(summ-(sum*sum)/dble(dnn))/dble(dn)
	curmin=cmin
c
c   curmin is the variance of the n+1 fn values at the vertices
c
	if(curmin.ge.eps) go to 1000
c
c   factorial tests to check that ynewlo is a local minimum
c
  22	do 23 i=1,n
  23	min(i)=p(i,ihi)
	ynewlo=y(ihi)
	if(maxit.gt.kcount) then
		ier=1
		return
	endif
	do 24 i=1,n
	del=step(i)*.001
	min(i)=min(i)+del
	z=armalk(min,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	if(z.lt.ynewlo) go to 25
	min(i)=min(i)-del-del
	z=armalk(min,np,nq,ny,dat,alpha,beta,e,vw,vl,vk,rv,ier)
	if(ier.gt.0) return
	if(z.lt.ynewlo) go to 25
  24	min(i)=min(i)+del
	am2ll=z
	var=curmin
	return
c
c   restart procedure
c
  25	do 26 i=1,n
  26	start(i)=min(i)
	del=.001
	maxit=maxit+1
	go to 1001
	end
c

c********************************************************************
c
	double precision function armalk(theta,np,nq,n,y,
     1  alpha,beta,e,vw,vl,vk,rv,ier)
c
c   Subroutine to calculate -2log likelihood for ARMA(np,nq) process
c
c*********************************************************************
c
        integer np,nq,n,ier
	double precision theta(1),y(n),alpha(1),beta(1)
        double precision e(n),vw(1),vl(1),vk(1)
        double precision sumsq, fact, tol,e1,al2pi,c,rv
c
c   find -alpha and -beta corresponding to input transformed partials:
c
	if(np.gt.0) then
           do 5 i=1,np
           e1=exp(-theta(i))
   5	   alpha(i)=(1.-e1)/(1.+e1)
           call partar(alpha,np)
	   do 10 i=1,np
  10	   alpha(i)=-alpha(i)
        endif
  	if(nq.gt.0) then
           do 20 i=1,nq
   	   e1=exp(-theta(np+i))
  20	   beta(i)=(1.-e1)/(1.+e1)
           call partar(beta,nq)
	   do 25 i=1,nq
  25	   beta(i)=-beta(i)
        endif
c
c   call evaluator:
c
	tol=.0001
	call lkhood(alpha,np,beta,nq,y,e,n,sumsq,fact,vw,vl,vk,tol,ier)
c
c   
	if(ier.gt.0) return
	if(np.gt.0) then
	   do 35 i=1,np
  35	   alpha(i)=-alpha(i)
	endif
	if(nq.gt.0) then
 	   do 40 i=1,nq
  40	   beta(i)=-beta(i)
	endif
	rv=sumsq/dble(n)
c	al2pi=alog(8.*atan(1.0))
c	c=float(n)*(1.+al2pi+alog(float(rv))+alog(float(fact)))
c
c
	al2pi=log(8.d0*atan(1.d0))
	c=dble(n)*(1.+al2pi+log(rv)+log(fact))
c
c
	armalk=c
	return
	end
c
c********************************************************************
c
	subroutine lkhood(alpha,np,beta,nq,y,e,n,sumsq,fact,vw,vl,
     1	vk,tol,ier)
c
c   Subroutine to find the 2 terms in the exact ARMA likelihood for
c   -alpha and -beta as in Melard (Applied Statistics, 1984, 104)
c
c   vw and vl are max(p,q+1)+1 vk is max(p,q+1)
c
c********************************************************************
c
        integer np,nq,ier,n, mr,mrp1
	double precision alpha(1),beta(1),y(n),e(n),vw(1)
        double precision vl(1),vk(1),sumsq,fact,tol
        double precision a,r,fn,detcar,detman
        double precision vl1,vw1,alf,flj,aor
	data epsil1/1.d-10/
c
c
	mr=max0(np,nq+1)
	mrp1=mr+1
	fact=0.D0
	detman=1.D0
	detcar=0.D0
	sumsq=0.D0
	mxpq=max0(np,nq)
	mxpqp1=mxpq+1
	mqp1=nq+1
	mpp1=np+1
c
c
	call wilson(alpha,np,beta,nq,vw,mxpqp1,vl,mxpqp1,vk,mxpq,ier)
	if(ier.gt.0) return
	vk(1)=vw(1)
	if(mr.eq.1) go to 15
	do 14 k=2,mr
		vk(k)=0.
		if(k.gt.np) go to 12
		do 11 j=k,np
  11			vk(k)=vk(k)+alpha(j)*vw(j+2-k)
  12		if(k.gt.mqp1) go to 14
		do 13 j=k,mqp1
  13			vk(k)=vk(k)-beta(j-1)*vl(j+1-k)
  14	continue
c
c
  15	r=vk(1)
	vl(mr)=0.
	do 16 j=1,mr
		vw(j)=0.
		if(j.ne.mr) vl(j)=vk(j+1)
		if(j.le.np) vl(j)=vl(j)+alpha(j)*r
  16    vk(j)=vl(j)
c
c
	last=mpp1-nq
	loop=np
	jfrom=mpp1
	vw(mpp1)=0.
	vl(mxpqp1)=0.
c
c
	if(n.le.0) go to 50
	do 29 i=1,n
c
c
		if(i.ne.last) go to 17
		loop=min0(np,nq)
		jfrom=loop+1
c
c
		if(nq.le.0) go to 30
  17		if(r.le.epsil1) go to 40
		if(abs(r-1.).lt.tol.and.i.gt.mxpq) go to 30
c
c
		detman=detman*r
  19		if(abs(detman).lt.1.) go to 20
		detman=detman*.0625
		detcar=detcar+4.
		go to 19
  20		if(abs(detman).ge..0625) go to 21
		detman=detman*16.
		detcar=detcar-4.
		go to 20
  21		vw1=vw(1)
		a=y(i)-vw1
		e(i)=a/sqrt(r)
		aor=a/r
		sumsq=sumsq+a*aor
		vl1=vl(1)
		alf=vl1/r
		r=r-alf*vl1
		if(loop.eq.0) go to 23
c
c
		do 22 j=1,loop
			flj=vl(j+1)+alpha(j)*vl1
			vw(j)=vw(j+1)+alpha(j)*vw1+aor*vk(j)
			vl(j)=flj-alf*vk(j)
  22 		vk(j)=vk(j)-alf*flj
  23		if(jfrom.gt.nq) go to 25
		do 24 j=jfrom,nq
			vw(j)=vw(j+1)+aor*vk(j)
			vl(j)=vl(j+1)-alf*vk(j)
  24		vk(j)=vk(j)-alf*vl(j+1)
  25		if(jfrom.gt.np) go to 27
		do 26 j=jfrom,np
  26		vw(j)=vw(j+1)+alpha(j)*y(i)
  27		continue
  29	continue
	go to 39
c
c
  30	nexti=i
	ier=-nexti
	do 31 i=nexti,n
  31	e(i)=y(i)
        if(np.eq.0) go to 34
        do 33 i=nexti,n
		do 32 j=1,np
  32		e(i)=e(i)-alpha(j)*y(i-j)
  33	continue
  34	if(nq.eq.0) go to 37
	do 36 i=nexti,n
		do 35 j=1,nq
  35		e(i)=e(i)+beta(j)*e(i-j)
  36	continue
c
c
  37	do 38 i=nexti,n
  38	sumsq=sumsq+e(i)*e(i)
  39	fn=n
	fact=detman**(1./fn)*2.**(detcar/fn)
	return
  40	ier=8
	return
  50	ier=9
	return
	end

c*************************************************************************
c
      subroutine filt(beta,beta0,m,n,x,y)
c
c   y(t) = beta0*x(t) + sum(j=1,m) beta(j)x(t-j),  t=m+1,...,n
c
c*************************************************************************
c
      integer t,tpm,m,n
      double precision beta(1),x(1),y(1),c,beta0
c
      do 20 t=1,n-m
         tpm=t+m
         c=beta0*x(tpm)
         do 10 j=1,m
 10      c=c+beta(j)*x(tpm-j)
 20   y(t)=c
c
      return
      end

c*********************************************************************
c
      subroutine movave(x,n,k2,xave)
c
c   Find moving averages of length k2 (except at edges where
c   as many as available are used).
c
c**********************************************************************
c
      integer n,k2,oi,ok
      double precision x(n),xave(n)
c
c   Get sums:
c
      k=(k2-1)/2
      xave(1)=x(1)
      xave(n)=x(n)
      ii=1
      jj=n
      do 10 i=2,k+1
         ii=ii+2
         jj=jj-2
         xave(n-i+1)=xave(n-i+2)+x(jj)+x(jj+1)
 10      xave(i)=xave(i-1)+x(ii)+x(ii-1)
c
      kp1=k+1
      do 20 i=k+2,n-k-1
 20   xave(i)=xave(i-1)+x(i+k)-x(i-kp1)
c
      do 30 i=1,k+1
         oi=dble(2*i-1)
         xave(i)=xave(i)/oi
 30      xave(n-i+1)=xave(n-i+1)/oi
c
      ok=2*k+1
      do 40 i=k+2,n-k-1
 40   xave(i)=xave(i)/ok
c
      return
      end
c
c***********************************************************************
c
      subroutine movbox(x,n,k2,ii,jj,xsumm,nouts,indout,outval)
c
c   Find moving five-pt summaries and extreme values (as in a box-plot)
c   for a time series x of length n. The moving window has k2 (odd) points
c   (except for ends where as many observations as possible are used).
c
c   xsumm is (5 x n): (m is number of points in the window for column i)
c
c       Row 1 is largest x in [u4, u4 + 1.5 (u4-l4)] or u4 if none
c       Row 2 is upper fourth (median of largest (m+1)/2)
c       Row 3 is median
c       Row 4 is lower fourth (median of smallest (m+1)/2)
c       Row 5 is smallest x in [l4, l4 - 1.5 (u4-l4)] or l4 if none
c
c   indout and outval are arrays (integer and real) of indices and 
c   corresponding values outside of [l4 - 1.5 (u4-l4), u4 + 1.5 (u4-l4)].
c
c   indout and outval can be no longer than n elements
c
c   The arrays ii and jj are integer work arrays of length k2.
c
c***********************************************************************
c
      integer n,k2,nouts,indout
      double precision x(n),ii(k2),jj(k2)
      double precision xsumm(5,n),outval(1)
c
      k=(k2-1)/2
      ii(1)=1
      jj(1)=1
      nl=1
      call summ5(ii,jj,nl,1,x,xsumm(1,1),nouts,indout,outval)
c
      do 10 m=2,k+1
      call medadd(ii,jj,nl,x)      
      nl=nl+1
      call medadd(ii,jj,nl,x)
      nl=nl+1
 10   call summ5(ii,jj,nl,m,x,xsumm(1,m),nouts,indout,outval)
c
      do 20 m=k+2,n-k
      call meddel(ii,jj,k2)
      call medadd(ii,jj,k2-1,x)
 20   call summ5(ii,jj,k2,m,x,xsumm(1,m),nouts,indout,outval)
c
      nl=k2
      do 30 m=n-k+1,n
      call meddel(ii,jj,nl)
      nl=nl-1
      call meddel(ii,jj,nl)
      nl=nl-1
 30   call summ5(ii,jj,nl,m,x,xsumm(1,m),nouts,indout,outval)
c
      return
      end

c*********************************************************************
c
      subroutine summ5(ii,jj,n,npos,x,xsumm,nouts,indout,outval)
c
c   Subroutine to 
c
c
c*********************************************************************
c
      integer n,npos,indout(1)
      double precision ii(n),jj(n),x(1),xsumm(5)
      double precision outval(1),xmid,flow,d,fup
c
c
      m=(n+1)/2
      call medii(ii,n,x,xsumm(3))
      call medii(ii(n-m+1),m,x,xsumm(2))
      call medii(ii,m,x,xsumm(4))
c
      d=1.d5*(xsumm(2)-xsumm(4))
      fup=xsumm(2)+d
      flow=xsumm(4)-d
      xsumm(1)=xsumm(2)
      xsumm(5)=xsumm(4)
c
      nn=0
      do 10 i=m,n
      if(x(ii(i)).gt.fup) go to 20
 10   if(x(ii(i)).gt.xsumm(2)) nn=i
 20   if(nn.ne.0) xsumm(1)=x(ii(nn))
c
      nn=0
      do 30 i=m,1,-1
      if(x(ii(i)).lt.flow) go to 40
 30   if(x(ii(i)).lt.xsumm(4)) nn=i
 40   if(nn.ne.0) xsumm(5)=x(ii(nn))
c
      xmid=x(npos)
      if(xmid.lt.flow.or.xmid.gt.fup) then
         nouts=nouts+1
         indout(nouts)=npos
         outval(nouts)=xmid
      endif
c
      return
      end
c
c*********************************************************************
c
      subroutine medii(ii,n,x,xmed)
c
c**********************************************************************
c
      integer n
      double precision ii(1),x(1),xmed
c
      xmed=(x(ii((n+1)/2))+x(ii((n+2)/2)))/2
c
      return
      end

c***********************************************************************
c
      subroutine movmed(x,n,k2,ii,jj,xmed)
c
c   Find moving medians (xmed(1),...xmed(n)) of x(1),...,x(n) where 
c   (except at ends) medians are of k2 x's. At end points, medians 
c   are of as many points as possible (for example, xmed(1) and xmed(n)
c   are medians of 1 point).
c
c   The arrays ii and jj are integer work arrays of length k2.
c
c***********************************************************************
c
      integer n,k2
      double precision x(n),ii(k2),jj(k2),xmed(n)
c
      k=(k2-1)/2
      ii(1)=1
      jj(1)=1
      xmed(1)=x(1)
c
      do 10 m=2,k+1
      call medadd(ii,jj,2*m-3,x)      
      call medadd(ii,jj,2*m-2,x)
 10   xmed(m)=x(ii(m))
c
      do 20 m=k+2,n-k
      call meddel(ii,jj,k2)
      call medadd(ii,jj,k2-1,x)
 20   xmed(m)=x(ii(k+1))
c
      nl=k2
      do 30 m=n-k+1,n
      call meddel(ii,jj,nl)
      nl=nl-1
      call meddel(ii,jj,nl)
      nl=nl-1
 30   xmed(m)=x(ii((nl+1)/2))
c
      return
      end
c
c***********************************************************************
c
      subroutine movord(x,n,k2,nord,ii,jj,xord)
c
c   Find moving nord order statistic (xord(1),...xord(n)) of 
c   x(1),...,x(n) where (except at ends) order statistics are of k2 x's. 
c   At end points, they are of as many points as possible (for example, 
c   xord(1) and xord(n) are for 1 point).
c
c   The arrays ii and jj are integer work arrays of length k2.
c
c***********************************************************************
c
      integer n,k2
      double precision x(n),ii(k2),jj(k2),xord(n)
c
      k=(k2-1)/2
      ii(1)=1
      jj(1)=1
      xord(1)=x(1)
c
      do 10 m=2,k+1
         call medadd(ii,jj,2*m-3,x)      
         call medadd(ii,jj,2*m-2,x)
         if(m.eq.k+1) xord(m)=x(ii(nord))
         if(m.lt.k+1) xord(m)=x(ii(m))
 10   continue
c
      do 20 m=k+2,n-k
      call meddel(ii,jj,k2)
      call medadd(ii,jj,k2-1,x)
 20   xord(m)=x(ii(nord))
c
      nl=k2
      do 30 m=n-k+1,n
      call meddel(ii,jj,nl)
      nl=nl-1
      call meddel(ii,jj,nl)
      nl=nl-1
 30   xord(m)=x(ii((nl+1)/2))
c
      return
      end
c
c***********************************************************************
c
      subroutine meddel(ii,jj,nl)
c   Delete a point in moving median
c
c**********************************************************************
c
      integer nl
      double precision ii(1),jj(1)
c
      nr=jj(1)
      noff=ii(nr)-1
c
      if(nr.lt.nl) then
         do 10 i=nr+1,nl
         i1=ii(i)-noff
         jj(i1)=jj(i1)-1
 10      ii(i-1)=ii(i)
      endif
c
      do 20 i=2,nl
 20   jj(i-1)=jj(i)
c
      return
      end
c
c***********************************************************************
c
      subroutine medadd(ii,jj,nl,x)
c
c   Add a point in moving median
c
c***********************************************************************
c
      integer bsrch,nl
      double precision xnew
      double precision x(1),ii(1),jj(1)
c
c
      xnew=x(ii(jj(nl))+1)
      noff=ii(jj(1))-1
c
      ns=bsrch(x,ii,nl,xnew)
c
      if(ns.le.nl) then
         do 30 i=ns,nl
         i1=ii(i)-noff
 30      jj(i1)=jj(i1)+1
         do 40 i=nl,ns,-1
 40      ii(i+1)=ii(i)
      endif
c
      jj(nl+1)=ns
      ii(ns)=nl+noff+1
c
      return
      end
c
c*************************************************************************
c
      integer function bsrch(x,ii,n,xnew)
c
c   Find location of xnew in x(ii(1)) <= x(ii(2)) <= ... <= x(ii(nl))
c
c*************************************************************************
c
      integer n
      double precision x(1),ii(1),xnew
c
      nl=0
      nr=n+1
c
 5    if(nr-nl.eq.1) go to 99
c
      mid=(nl+nr)/2
      if(xnew.gt.x(ii(mid))) then
         nl=mid
         go to 5
      else
         nr=mid
         go to 5
      endif
c
 99   bsrch=nr
      return
      end
c
c*************************************************************************
c
      subroutine pacf(x,n,m,e,f,temp,theta)
c
c   x of length n
c   e, f, and temp of length n+m
c   theta of length m
c
c*************************************************************************
c
      INTEGER n
      DOUBLE PRECISION x(1),e(1),f(1),temp(1),theta(1)
      double precision dtprod, apart,part,xbar
c
c
      call submns(x,n,xbar,e)
      npm=n+m
      do 10 i=n+1,npm
 10   e(i)=0.0
      call crlag(e,npm,f)
c
c
      do 30 i=1,m
         part=dtprod(f,e,npm)/dtprod(f,f,npm)
         apart=-part
         call movvec(e,npm,temp)
         call accvec(temp,apart,f,npm,e)
         call accvec(f,apart,temp,npm,f)
         call crlag(f,npm,f)
         theta(i)=part
 30   continue
c
c
      return
      end
c
c***********************************************************************
c
      subroutine accvec(x,a,y,n,z)
c
c   Find z=x+a*y
c
c************************************************************************
c
      integer n
      double precision x(1),y(1),z(1),a
c
c
      do 10 i=1,n
 10   z(i)=x(i)+a*y(i)
c
      return
      end
c
c*************************************************************************
c
      double precision function dtprod(x,y,n)
c
c*************************************************************************
c
      integer n
      double precision x(1),y(1)
c
      dtprod=0.0d0
      do 10 i=1,n
 10   dtprod=dtprod+dble(x(i))*dble(y(i))
c
      return
      end
c
c**************************************************************************
c
      subroutine movvec(x,n,y)
c
c***************************************************************************
c
      integer n
      double precision x(1),y(1)
c
      do 10 i=1,n
 10   y(i)=x(i)
c
      return
      end
c
c***************************************************************************
c
      subroutine crlag(x,n,y)
c
c
c***************************************************************************
c
      integer n
      double precision x(1),y(1),temp
c
c
      temp=x(n)
      do 10 i=n,2,-1
 10   y(i)=x(i-1)
      y(1)=temp
c
      return
      end
c
c***************************************************************************
c
      subroutine submns(x,n,xbar,y)
c
c***************************************************************************
c
      integer n
      double precision x(1),y(1),xbar
c
c
      xbar=0.0
      do 10 i=1,n
 10   xbar=xbar+x(i)
      xbar=xbar/dble(n)
c
      do 20 i=1,n
 20   y(i)=x(i)-xbar
c
      return
      end
c
c***********************************************************************
c
      subroutine partar(alpha,p)
c
c***********************************************************************
c
      integer p
      double precision alpha(p),temp,c
c
      if(p.gt.1) then
         do 20 j=2,p
            c=alpha(j)
            do 10 k=1,j/2
               temp=alpha(k)-c*alpha(j-k)
               alpha(j-k)=alpha(j-k)-c*alpha(k)
 10            alpha(k)=temp
 20      continue
      endif
c
      do 30 j=1,p
 30   alpha(j)=-alpha(j)
c
c
      return
      end
c
c**********************************************************************
c
      subroutine poly(coeffs,degree,x,n,fx)
c
c   Evaluate polynomial
c
c   f(x)=coeffs(1) + sum(j=1,degree) coeffs(j+1)*x^j
c
c   for x(1),...,x(n)
c
c***********************************************************************
c
      integer degree
      double precision coeffs(1),x(n),fx(n)
      double precision c,dc,xi
c
c
      dc=dble(coeffs(degree+1))
      do 20 i=1,n
         xi=dble(x(i))
         c=dc
         do 10 j=degree,1,-1
 10      c=c*xi+dble(coeffs(j))
 20   fx(i)=c
c
c
      return
      end
c
c********************************************************************
c
      subroutine rtpoly(r,np,wk,a)
c
c   Subroutine to find the coefficients a(1),...,a(np) of the polynomial
c
c      1 + Sum(j=1,np) a(j)z**j
c
c   given its roots r(1),...,r(np)
c
c   Input:  r (complex), np
c
c   Output: a (real)
c
c   Work:   wk (complex array of length np)
c
c**********************************************************************
c
      integer np
      double precision a(np)
      double precision r(2,np),wk(2,np),zi(2),cone(2)
      double precision ccreal
      data cone/1.0, 0.0/
c
      call ccdiv(cone,r(1,1),wk(1,1))
      call ccneg(wk(1,1),wk(1,1))
      if(np.eq.1) go to 30
c
      do 20 i=2,np
         call ccdiv(cone,r(1,i),zi)
         call ccneg(zi,zi)
         call ccmul(zi,wk(1,i-1),wk(1,i))
         if(i.eq.2) go to 15
         do 10 j=i-1,2,-1
 10      call ccacc(wk(1,j),zi,wk(1,j-1),wk(1,j))
 15      call ccacc(wk(1,1),zi,cone,wk(1,1))
 20   continue
c
 30   do 40 i=1,np
 40   a(i)=ccreal(wk(1,i))
c
      return
      end
c
c**********************************************************************
c
      subroutine ccadd(x,y,z)
c
c**********************************************************************
c
      double precision x(2),y(2),z(2)
c
      z(1)=x(1)+y(1)
      z(2)=x(2)+y(2)
      return
      end
c
c**********************************************************************
c
      subroutine ccsub(x,y,z)
c
c**********************************************************************
c
      double precision x(2),y(2),z(2)
c
      z(1)=x(1)-y(1)
      z(2)=x(2)-y(2)
      return
      end
c
c**********************************************************************
c
      subroutine ccneg(x,y)
c
c**********************************************************************
c
      double precision x(2),y(2)
c
      y(1)=-x(1)
      y(2)=-x(2)
      return
      end
c
c**********************************************************************
c
      subroutine ccmul(x,y,z)
c
c**********************************************************************
c
      double precision x(2),y(2),z(2),temp
c
      temp=x(1)*y(1)-x(2)*y(2)
      z(2)=x(1)*y(2)+x(2)*y(1)
      z(1)=temp
      return
      end
c
c**********************************************************************
c
      subroutine ccdiv(x,y,z)
c
c**********************************************************************
c
      double precision x(2),y(2),z(2),temp,y2
c
      y2=y(1)*y(1)+y(2)*y(2)
      temp=(x(1)*y(1)+x(2)*y(2))/y2
      z(2)=(x(2)*y(1)-x(1)*y(2))/y2
      z(1)=temp
      return
      end
c
c**********************************************************************
c
      subroutine ccacc(x,a,y,z)
c
c**********************************************************************
c
c
c   z = x + a * y
c
      double precision x(2),a(2),y(2),z(2),temp
c
      temp=x(1)+a(1)*y(1)-a(2)*y(2)
      z(2)=x(2)+a(2)*y(1)+a(1)*y(2)
      z(1)=temp
      return
      end
c
c**********************************************************************
c
      subroutine cccopy(x,y)
c
c**********************************************************************
c
      double precision x(2),y(2)
c
      x(1)=y(1)
      x(2)=y(2)
      return
      end
c
c**********************************************************************
c
      DOUBLE PRECISION function ccimag(x)
c
c**********************************************************************
c
      double precision x(2)
c
      ccimag=x(2)
      return
      end
c
c**********************************************************************
c
      DOUBLE PRECISION function ccreal(x)
c
c**********************************************************************
c
      double precision x(1)
c
      ccreal=x(1)
      return
      end
c
c**********************************************************************
c
      subroutine ccplx(a,b,z)
c
c**********************************************************************
c
      double precision a,b,z(2)
c
      z(1)=a
      z(2)=b
      return
      end
c
c**********************************************************************
c
      DOUBLE PRECISION function ccabs(x)
c
c**********************************************************************
c
      double precision x(2)
c
      ccabs=sqrt(x(1)*x(1)+x(2)*x(2))
      return
      end
c
c**********************************************************************
c
      subroutine crmul(a,x,y)
c
c**********************************************************************
c
      double precision a,x(2),y(2)
c
      y(1)=a*x(1)
      y(2)=a*x(2)
      return
      end
c
c**********************************************************************
c
      subroutine ccsqrt(x,y)
c
c**********************************************************************
c
      double precision x(2),y(2),xx,yy,r,w
c
      if(x(1).eq.0.0.and.x(2).eq.0.0) then
         y(1)=0.D0
         y(2)=0.D0
         return
      endif
      xx=abs(x(1))
      yy=abs(x(2))
      if(xx.ge.yy) then
         r=yy/xx
         w=sqrt(xx)*sqrt(0.5*(1.0+sqrt(1.0+r*r)))
      else
         r=xx/yy
         w=sqrt(yy)*sqrt(0.5*(r+sqrt(1.0+r*r)))
      endif
      if(x(1).ge.0.0) then
         y(1)=w
         y(2)=x(2)/(2.0*w)
      else
         y(2)=-w
         if(x(2).ge.0.0) y(2)=w
         y(1)=x(2)/(2.0*y(2))
      endif
      return
      end
c
c*****************************************************************
c
	subroutine polyrt(ra,np,maxit,eps,a,a1,roots,ier)
c
c   calculate roots of a(1)+a(2)z+a(3)z**2+...+a(np+1)z**np
c
c*****************************************************************
c
        integer np,maxit,ier,i
	double precision ra(1)
        double precision a(2,1),a1(2,1),roots(2,1),x(2),b(2),c(2)
        double precision cone(2),czero(2)
        double precision ccimag
        double precision ccreal
        data cone,czero/1.0, 0.0, 0.0, 0.0/
c
        call cccopy(a(1,1),cone)
        call cccopy(a1(1,1),cone)
c
	do 1 i=2,np+1
	a(1,i)=ra(i-1)
        a(2,i)=0.0
 1      call cccopy(a1(1,i),a(1,i))
	do 3 j=np,1,-1
           call cccopy(x,czero)
           call root1(a1,j,x,eps,0,maxit,ier)
           if(ier.eq.1) then
              ier=2
              go to 99
           endif
           if(abs(ccimag(x)).le.2.*eps**2*abs(ccreal(x))) then
              call ccplx(ccreal(x),0.d0,x)
           endif
           call cccopy(roots(1,j),x)
           call cccopy(b,a1(1,j+1))
              do 2 jj=j,1,-1
                 call cccopy(c,a1(1,jj))
                 call cccopy(a1(1,jj),b)
 2               call ccacc(c,b,x,b)
 3	continue
	do 4 j=1,np
    	call root1(a,np,roots(1,j),eps,1,maxit,ier)
   4	if(ier.eq.1) go to 99
	do 6 j=2,np
           call cccopy(x,roots(1,j))
           do 5 i=j-1,1,-1
              if(ccreal(roots(1,i)).le.ccreal(x)) go to 10
              call cccopy(roots(1,i+1),roots(1,i))
   5       continue
  	   i=0
  10	   call cccopy(roots(1,i+1),x)
   6	continue
  99	continue
	return
	end
c
c********************************************************************
c
      subroutine root1(a,np,x,eps,iopt,maxit,ier)
c
c   calculate one root.
c
c********************************************************************
c
      integer np,maxit,ier,iopt
      double precision a(2,1),x(2),dx(2),x1(2),b(2),d(2),f(2)
      double precision g(2),h(2),sq(2),gp(2),gm(2),g2(2),czero(2)
      double precision t1(2),t2(2)
      double precision ccabs,dxold,cdx
      data czero/0.0, 0.0/
c
      ier=0
      if(iopt.eq.1) then
         dxold=ccabs(x)
         npol=0
      endif
      do 2 it=1,maxit
         call cccopy(b,a(1,np+1))
         call cccopy(d,czero)
         call cccopy(f,czero)
         do 1 j=np,1,-1
            call ccacc(d,x,f,f)
            call ccacc(b,x,d,d)
 1          call ccacc(a(1,j),x,b,b)
         if(ccabs(b).le.1.e-15) then
            call cccopy(dx,czero)
         else if(ccabs(d).le.1.e-15.and.ccabs(f).le.1.e-15) then
            call ccdiv(b,a(1,np+1),t1)
            call ccplx(ccabs(t1)**(1./np),0.d0,dx)
         else
            call ccdiv(d,b,g)
            call ccmul(g,g,g2)
            call ccplx(-2.d0,0.d0,t1)
            call ccdiv(f,b,t2)
            call ccacc(g2,t1,t2,h)
            call crmul(dble(np),h,t1)
            call ccsub(t1,g2,t2)
            call crmul(dble(np-1),t2,t1)
            call ccsqrt(t1,sq)
            call ccadd(g,sq,gp)
            call ccsub(g,sq,gm)
            if(ccabs(gp).lt.ccabs(gm)) call cccopy(gp,gm)
            call ccplx(dble(np),0.d0,t1)
            call ccdiv(t1,gp,dx)
         endif
         call ccsub(x,dx,x1)
         if(x(1).eq.x1(1).and.x(2).eq.x1(2)) return
         call cccopy(x,x1)
         if(iopt.eq.1) then
            npol=npol+1
            cdx=ccabs(dx)
            if(npol.gt.9.and.cdx.ge.dxold) return
            dxold=cdx
         else
            if(ccabs(dx).le.eps*ccabs(x)) return
         endif
 2    continue
      ier=1
      return
      end
c

c*******************************************************
c
      subroutine schur(alph,ndim,np,a)
c
c   subroutine to form the schur matrix of the coefficients
c   alph(1),...,alph(np) (see pagano : when is an
c   autoregressive scheme stationary )
c
c   input :
c           ndim : dimension of matrix a in main program
c           np, alph(1),...,alph(np)
c
c   output :
c           a
c
c   subroutines called : none
c
c*******************************************************
c
      integer ndim,np
      double precision a(ndim,ndim),alph(np),c,d1,d2
c
      do 1 j=1,np
      do 1 k=j,np
      c=0.d0
         do 6 l=1,j
         if(l.eq.j) go to 2
         d1=alph(j-l)
         go to 3
   2     d1=1.d0
   3     if(l.eq.k) go to 4
         d2=alph(k-l)
         go to 6
   4     d2=1.d0
   6     c=c+d1*d2-alph(np+l-j)*alph(np+l-k)
      a(j,k)=c
   1  a(k,j)=a(j,k)
      return
      end
c
c******************************************************************
c
        subroutine marq(y,n,ndim,maxit,eps,nords,coeffs,nt,npnt,ntot,
     1  e,e1,xx,x,a,as,lags,ier,rvar,d)
c
c   subroutine to implement the algorithm on page 504 of box and jenkins
c   except in the language of timeslab.
c
c******************************************************************
c
        integer n,ndim,maxit,nt,npnt,ntot,ier,nords(10), lags
        double precision dtprod, eps, rvar, d(1)
        double precision y(2),coeffs(2),e(2),e1(2)
        double precision xx(2),x(npnt,2),a(ndim,ndim)
        double precision as(ndim,ndim)
        double precision g(100),h(100)
        double precision lp(100),lq(100),cc,gg,gg1,gg2
        integer nl(100)
        double precision alpha(100),beta(100),a1(100),beta1(100)
        double precision ss,ss1,ppi,f2,delta
        double precision seaslk
c
c   ppi, f2, and delta have the meaning in B+J
c
        ppi=0.d1
        f2=2.d0
        delta=0.d05
        npnt=n+nt
        ntp=ntot+1
c
c
c   start iterations:
c
c
        do 100 it=1,maxit
c
           ss=seaslk(y,n,nords,coeffs,lags,lp,lq,alpha,beta,
     1               e,a1,nl,xx,nt)

c
c   Get x matrix:
c
           do 10 i=1,ntot
              coeffs(i)=coeffs(i)+delta
              ss1=seaslk(y,n,nords,coeffs,lags,lp,lq,alpha,beta,
     1                   e1,a1,nl,xx,nt)
              coeffs(i)=coeffs(i)-delta
              do 20 j=1,npnt
  20          x(j,i)=(e(j)-e1(j))/delta
  10       continue
c
c   Get a, d, g:
c
           do 15 i=1,ntot
              do 30 j=1,i
              a(i,j)=dtprod(x(1,i),x(1,j),npnt)
  30          a(j,i)=a(i,j)
              g(i)=dtprod(x(1,i),e,npnt)
              d(i)=sqrt(a(i,i))
  15       continue
           if(maxit.eq.1) then
              ier=1
              go to 200
           endif
c
c
  25       do 40 i=1,ntot
              as(ntp,i)=g(i)/d(i)
              as(i,ntp)=as(ntp,i)
              do 40 j=1,ntot
              as(i,j)=a(i,j)/(d(i)*d(j))
  40       continue
           do 50 i=1,ntot
  50       as(i,i)=as(i,i)+ppi
           as(ntp,ntp)=1.
c
c   solve system:
c
           call swpk12(as,ntp,ntp,1,ntot,ier)
c
c   check for errors in solving system:
c
           if(ier.ne.0) then
              ier=2
              go to 99
           endif
c
c   update and get new sum of squares:
c
            do 60 i=1,ntot
               h(i)=as(i,ntp)/d(i)
  60        beta1(i)=coeffs(i)+h(i)
            ss1=seaslk(y,n,nords,beta1,lags,lp,lq,alpha,beta,
     1                 e1,a1,nl,xx,nt)
c
c   check for decrease in sum of squares:
c
            if(ss1.ge.ss) then
c
c   no:
c
               if(ppi.gt.200.) then
                  ier=3
                  go to 99
               endif
               ppi=ppi*f2
               go to 25
            endif
c
c   yes:
c
c   check for convergence:
c
            do 70 i=1,ntot
  70        if(abs(h(i)).gt.eps*abs(coeffs(i))) go to 75
c
c   yes:
c
            call movxyr(coeffs,beta1,ntot)
            ier=0
            go to 200
c
c   no:
c
  75        call movxyr(coeffs,beta1,ntot)
            if(ppi.lt.1.e-20) then
                ier=4
                go to 99
            endif
            ppi=ppi/f2
c
c   end iteration:
c
 100    continue
c
c   ran out of iterations:
c
        ier=1
 200    continue
  99    continue
c
c   get standard errors:
c
        rvar=ss/n
        call swpk12(a,ndim,ntot,1,ntot,ier)
        do 110 i=1,ntot
 110    d(i)=sqrt(rvar*a(i,i))
c
c   get constant:
c
        if(nords(5).eq.1) then
c
           cc=coeffs(ntot)
           gg1=1.
           if(nords(1).gt.0) then
              do 120 i=1,nords(1)
 120          gg1=gg1+coeffs(i)
           endif
           gg2=1.
           if(nords(2).gt.0) then
              do 130 i=1,nords(2)
 130          gg2=gg2+coeffs(nords(1)+i)
           endif
           gg=gg1*gg2
           cc=cc*gg
           coeffs(ntot)=cc
           d(ntot)=sqrt(rvar*d(ntot)*gg*gg)
        endif
c
        return
        end
c
c
c*********************************************************************
c
        double precision function seaslk(x,n,nords,coeffs,
     1  lags,lp,lq,alpha,beta,e,a,nl,xx,nt)
c
c
c   Function to evaluate a seasonal ARMA sum of squares.
c
c**********************************************************************
c
        integer n,nords,lags,nl,nt,ii,ij
        double precision x(1),coeffs(1),lp(1)
        double precision lq(1),alpha(1)
        double precision beta(1),e(1),a(1),xx(1)
        double precision c,c1,ss,amu
c
c
        call convt(nords,coeffs,lags,maxp,maxq,nlp,nlq,lp,lq,alpha,
     1  beta,a,nl,amu)
c
c
c
c   find eta:
c
        if(nt.eq.0) go to 715
        if(maxq.gt.0) then
           do 310 i=n-maxp+1,n-maxp+maxq
 310       e(i)=0.d0
        endif
        do 360 i=n-maxp,1,-1
           c=x(i)-amu
           if(nlp.gt.0) then
              do 320 j=1,nlp
 320          c=c+alpha(j)*(x(i+lp(j))-amu)
           endif
           if(nlq.gt.0) then
              do 340 k=1,nlq
 340          c=c-beta(k)*e(i+lq(k))
           endif
 360    e(i)=c
c
c   now get x(0),x(-1),... (x(i)=xx(nt+i), i=0,...,1-nt)
c
        do 400 i=0,-nt+1,-1
           c=0.0
           if(nlq.gt.0) then
                do 420 k=1,nlq
                ii=i+lq(k)
                if(ii.ge.1) c=c+beta(k)*e(ii)
 420            continue
           endif
           if(nlp.gt.0) then
                do 440 j=1,nlp
                ii=i+lp(j)
                if(ii.gt.0) c1=x(ii)-amu
                if(ii.le.0) c1=xx(nt+ii)-amu
 440            c=c-alpha(j)*c1
           endif
  400   xx(nt+i)=c+amu
c
c   now get e(-nt+1),...,e(n) and put them in indices 1,...,n+nt of e:
c
c
 715    ss=0.0
        ij=1-nt
        do 60 i=-nt+1,n
        if(i.gt.0) c=x(i)-amu
        if(i.le.0) c=xx(nt+i)-amu
           if(nlp.gt.0) then
                do 20 j=1,nlp
                ii=i-lp(j)
                c1=0.0
                if(ii.gt.0) c1=x(ii)-amu
                if(ii.le.0.and.ii.ge.ij) c1=xx(nt+ii)-amu
  20            c=c+alpha(j)*c1
           endif
           if(nlq.gt.0) then
              do 40 k=1,nlq
                 ii=i-lq(k)
                 if(ii.ge.ij) c=c-beta(k)*e(nt+ii)
  40          continue
           endif
        e(nt+i)=c
  60    ss=ss+c*c
c
c
        seaslk=ss
        return
        end
c
 
c*******************************************************************
c
       subroutine convt(nords,coeffs,lags,maxp,maxq,nlp,nlq,lp,lq,
     1  alpha,beta,a,nl,amu)
c
c   Subroutine to convert multiplicative stuff to subset arma.
c
c   Input:
c
c      nords  : p,P,q,Q,M
c      coeffs : full AR, subset AR, full MA, subset MA coefficients
c               and (if M=1) mean
c      lags   : subset AR and MA lags
c 
c   Output:
c 
c      maxp, maxq  : maximum AR and MA lags in subset form
c      nlp, nlq    : number of AR and MA lags in subset form
c      lp, lq      : arrays of length nlp, nlq having AR and MA lags
c      alpha, beta : arrays of length nlp, nlq having AR and MA coeffs
c      amu         : mean
c
c   Auxilliary:
c
c      a, nl : arrays of length max(maxp,maxq)
c      
c
c*******************************************************************
c
        integer nords(10),lags(1),nl(1),nlp,nlq,n1,n2,n3,n4
        integer maxp, maxq,jj
        double precision coeffs(1),lp(1),lq(1)
        double precision alpha(1),beta(1),a(1)
        double precision c1,amu
c
c
        np=nords(1)
        npl=nords(2)
        nq=nords(3)
        nql=nords(4)
        if(nords(5).eq.0) amu=0.d0
        if(nords(5).eq.1) amu=coeffs(np+npl+nq+nql+1)
        n1=0
        n2=n1+np
        n3=n2+npl
        n4=n3+nq
c
c
        maxp=np
        maxq=nq
        if(npl.gt.0) maxp=maxp+lags(npl)
        if(nql.gt.0) maxq=maxq+lags(npl+nql)
c
c
        nlp=0
        nlq=0
c
c
        if(maxp.eq.0) go to 50
c
c
        do 10 i=1,maxp
        nl(i)=0
  10    a(i)=0.0
c
c
        if(np.gt.0) then
                do 15 i=1,np
                nl(i)=1
  15            a(i)=coeffs(i)
        endif
c
c
        if(npl.eq.0) go to 30
c
        do 20 i=0,np
        c1=1.d0
        if(i.gt.0) c1=coeffs(i)
                do 20 j=1,npl
                jj=i+lags(j)
                nl(jj)=1
  20            a(jj)=a(jj)+c1*coeffs(n2+j)
c
  30    continue
c
        do 40 i=1,maxp
        if(nl(i).eq.0) go to 40
        nlp=nlp+1
        lp(nlp)=i
        alpha(nlp)=a(i)
  40    continue
c
c
  50    continue
c
c
        if(maxq.eq.0) go to 100
c
        do 60 i=1,maxq
        nl(i)=0
  60    a(i)=0.0
c
        if(nq.gt.0) then
                do 65 i=1,nq
                nl(i)=1
  65            a(i)=coeffs(n3+i)
        endif
c
        if(nql.eq.0) go to 80
c
        do 70 i=0,nq
        c1=1.d0
        if(i.gt.0) c1=coeffs(n3+i)
                do 70 j=1,nql
                jj=i+lags(npl+j)
                nl(jj)=1
  70            a(jj)=a(jj)+c1*coeffs(n4+j)
c
  80    continue
c
c
        do 90 i=1,maxq
        if(nl(i).eq.0) go to 90
        nlq=nlq+1
        lq(nlq)=i
        beta(nlq)=a(i)
  90    continue
c
c
 100    continue
        return
        end

c
c
c
c*******************************************************************
c
        subroutine sspr(x,n,nords,coeffs,lags,conf,ntf,ntl,nhl,
     1  e,xx,xp,xpl,xpu,ier,npds,rvar)
c
c
c*******************************************************************
c
        integer nords(10),lags,nl(100),ier,ntf,ntl,nhl,n
        integer nd, ndd, ns, maxq2, maxp2, ii
        double precision  rvar
        double precision fctr10
        double precision x(1),coeffs(1)
        double precision lp(100),lq(100)
        double precision alpha(100),beta(100)
        double precision a(100),a1(100),e(1),xx(1),xp(1)
        double precision xpl(1),xpu(1),xmin
        double precision c1,c2,alam,alam1,c,c3,amu,am
c
c
        ier=0
c
c   express model (without differencing) as subset ARMA:
c
        call convt(nords,coeffs,lags,maxp,maxq,nlp,nlq,lp,lq,
     1  alpha,beta,a,nl,amu)
c
c   adjust for differences:
c
        nd=nords(6)
        ndd=nords(7)
        ns=nords(8)
        maxq2=maxq
        maxp2=maxp+nd+ndd*ns
c
        if(maxp2.eq.0) go to 26
        do 10 i=1,maxp2
        a(i)=0.d0
  10    nl(i)=0
        do 20 i=0,nlp
           c1=1.d0
           ii=0
           if(i.gt.0) then
              ii=lp(i)
              c1=alpha(i)
           endif
           do 20 j=0,nd
              c2=((-1.d0)**j)*fctr10(nd,j)
              do 20 k=0,ndd
                 c3=((-1.d0)**k)*fctr10(ndd,k)
                 ijks=ii+j+k*ns
                 if((ijks).gt.0) then
                    nl(ijks)=1
                    a(ijks)=a(ijks)+c1*c2*c3
                 endif
 20     continue
        nlp=0
        do 25 i=1,maxp2
           if(nl(i).eq.1) then
              nlp=nlp+1
              lp(nlp)=i
              alpha(nlp)=a(i)
           endif
 25     continue
 26     continue
c
c   now nlp,nlq,alpha,beta,lp,lq are lags and coefficients of
c   full model
c
c   transform the series:
c
        ntot=nords(1)+nords(2)+nords(3)+nords(4)+nords(5)
        am=coeffs(ntot+1)
        alam=coeffs(ntot+2)
c
c   add the constant to make data positive:
c
        do 200 i=1,n
 200    x(i)=x(i)+am
c
c   power transform:
c
        iptlam=1
        if(alam.eq.1.0) go to 231
c 
c   check data + constant all positive:
c
         call vmin(x,n,xmin,imin)
         if(xmin.le.0.) then
            ier=1
            go to 99
         endif
c
c   do log transform:
c
         if(alam.eq.0.0) then
            iptlam=2
            do 210 i=1,n
 210        x(i)=log(x(i))
         endif
c
c   other power transforms:
c
         if(abs(alam).gt..0001.and.abs(alam).lt.1.) then
            iptlam=3
            do 220 i=1,n
 220        x(i)=x(i)**alam
         endif
c
 231     continue
c
c   forecasts: first one possible is x(maxp2+maxq2+1):
c
c
c   find one step ahead forecast errors for x(maxp2+1),...,x(n):
c
        do 250 i=1,n
 250    e(i)=0.0
        do 260 i=maxp2+1,n
           c=x(i)-amu
           if(nlp.gt.0) then
              do 270 j=1,nlp
                 if(i-lp(j).ge.1) c=c+alpha(j)*x(i-lp(j))
 270          continue
           endif
           if(nlq.gt.0) then
              do 280 j=1,nlq
                 if(i-lq(j).ge.1) c=c-beta(j)*e(i-lq(j))
 280          continue
           endif
 260    e(i)=c
c
c   find forecasts:
c
        npds=1
c
        do 350 nt=ntf,ntl
           call movxyr(xx,x,nt)
c
           do 360 nh=1,nhl
              ntph=nt+nh
              c=amu
                 if(maxq2.gt.0) then
                    do 370 j=1,nlq
                       jj=ntph-lq(j)
                       if((jj.lt.1).or.(jj.gt.nt)) go to 370
                       c=c+beta(j)*e(ntph-lq(j))
 370                continue
                 endif
                 if(maxp2.gt.0) then
                    do 380 j=1,nlp
 380                c=c-alpha(j)*xx(ntph-lp(j))
                 endif
 360             xx(ntph)=c
c
           call movxyr(xp(npds),xx(nt+1),nhl)
           npds=npds+nhl
 350    continue
        npds=npds-1
c
c   find psi weights (put them in xx):
c
        if(maxp2.gt.0) then
           do 410 i=1,maxp2
 410       a(i)=0.0
           do 420 i=1,nlp
 420       a(lp(i))=alpha(i)
        endif
        if(maxq2.gt.0) then
           do 430 i=1,maxq2
 430       a1(i)=0.0
           do 440 i=1,nlq
 440       a1(lq(i))=beta(i)
        endif
c
        do 450 i=1,nhl
 450    xx(i)=0.0
        if(maxq2.gt.0) then
           do 451 i=1,maxq2
 451       xx(i)=a1(i)
        endif
        if(maxp2.gt.0) then
           do 453 i=1,nhl
              c=xx(i)
              do 452 j=1,min0(i,maxp2)
                 c1=1.
                 if(i-j.gt.0) c1=xx(i-j)
 452          c=c-a(j)*c1
 453       xx(i)=c
        endif
c
c        do 450 i=1,nhl
c           c=0.
c           if(i.le.maxq2) c=a1(i)
c           if(maxp2.gt.0) then
c              do 460 j=1,min0(maxp2,i)
c                 c1=1.
c                 if(i-j.ne.0) c1=xx(i-j)
c 460          c=c-a(j)*c1
c           endif
c 450    xx(i)=c
c
        e(1)=rvar
        if(nhl.gt.1) then
           do 470 i=2,nhl
              c=1.
              do 480 j=1,i-1
 480          c=c+xx(j)*xx(j)
 470          e(i)=rvar*c
        endif
        do 471 i=1,nhl
 471    xx(i)=sqrt(e(i))
c
        do 490 i=1,nhl
 490    e(i)=conf*sqrt(e(i))
c
        do 500 nt=ntf,ntl
           ntt=(nt-ntf)*nhl
           do 500 nh=1,nhl
              nthh=ntt+nh
              xpl(nthh)=xp(nthh)-e(nh)
 500       xpu(nthh)=xp(nthh)+e(nh)
c
c   transform back:
c
        if(iptlam.eq.3) then
           call vmin(xpl,npds,xmin,imin)
           if(xmin.le.0.) then
              ier=1
              go to 99
           endif
           alam1=1./alam
           do 510 i=1,npds
              xp(i)=xp(i)**alam1
              xpl(i)=xpl(i)**alam1
 510          xpu(i)=xpu(i)**alam1
        endif
c
        if(iptlam.eq.2) then
           do 520 i=1,npds
              xp(i)=exp(xp(i))
              xpl(i)=exp(xpl(i))
 520          xpu(i)=exp(xpu(i))
        endif
c
 600    do 610 i=1,npds
           xp(i)=xp(i)-am
           xpl(i)=xpl(i)-am
 610       xpu(i)=xpu(i)-am
c
c
  99    continue
        return
        end
c


c*************************************************************************
c
        DOUBLE PRECISION function fctr10(n,k)
c
c   function to find the binomial coefficient { n choose k }.
c
c************************************************************************
c
        integer n,k
        double precision c,c1,c2
c
        if(k.eq.0.or.k.eq.n) then
                fctr10=1.d0
                return
        endif
        kk=min0(k,n-k)
        c=1.d0
        do 10 i=1,kk
        c1=dble(n-i+1)
        c2=dble(i)
  10    c=c*c1/c2
        fctr10=c
        return
        end
c
c*********************************************************
c
        subroutine vmin(x,n,xmin,ind)
c
c*********************************************************
c
        integer n,ind
        double precision x(1),xmin
c
        ind=1
        xmin=1
        do 10 i=1,n
           if(x(i).lt.xmin) then
              xmin=x(i)
              ind=i
           endif
 10     continue
        return
        end
c
c**********************************************************************
c
      subroutine swpk12(a,ndim,n,k1,k2,ier)
c
c   Subroutine to sweep the (n x n) matrix a on diagonals k1 through k2.
c
c   The output integer ier is 1 if an error occurs and 0 otherwise.
c
c*************************************************************************
c
      integer n,k1,k2,ier,ndim
      double precision a(ndim,1),diag
c
c
      ier = 1
c
      do 30 k = k1, k2
         diag = a(k,k)
         if(abs(diag).lt.1.e-20) return
c
         do 10 i = 1, n
         do 10 j = 1, n
 10      if(i.ne.k.and.j.ne.k) a(i,j) = a(i,j) - a(i,k)*a(k,j)/diag
c
         do 20 i = 1, n
         a(i,k) = -a(i,k) / diag
 20      a(k,i) =  a(k,i) / diag
c
         a(k,k) = 1 /diag
c
 30   continue
c
      ier = 0
c
c
      return
      end


c********************************************************************
c
	subroutine wilson(alpha,np,beta,nq,acf,ma,cvli,mxpqp1,alph,
     1	mxpq,ier)
c
c
c   subroutine to find mx covariances
c
c********************************************************************
c
        integer np,nq,mxpqp1,ier,mxpq,ma,j1
	double precision alpha(*),beta(*),acf(*)
        double precision cvli(*),alph(*)
        double precision div

c
	data epsil2/1.d-20/
c
c
	ier=0
	if(np.lt.0.or.nq.lt.0) ier=1
	if(mxpq.ne.max0(np,nq)) ier=2
	if(mxpqp1.ne.mxpq+1) ier=3
	if(ma.lt.mxpqp1) ier=4
	if(ier.gt.0) return
c
c
	acf(1)=1.d0
	cvli(1)=1.d0
	if(ma.eq.1) return
	do 1 i=2,ma
   1	acf(i)=0.d0
	if(mxpqp1.eq.1) return
	do 2 i=2,mxpqp1
   2	cvli(i)=0.d0
	do 9 k=1,mxpq
   9	alph(k)=0.d0
c
c
	if(nq.eq.0) go to 18
	do 13 k=1,nq
		cvli(k+1)=-beta(k)
		acf(k+1)=-beta(k)
		kc=nq-k
		if(kc.eq.0) go to 12
		do 11 j=1,kc
  11		acf(k+1)=acf(k+1)+beta(j)*beta(j+k)
  12		acf(1)=acf(1)+beta(k)*beta(k)
  13	continue
c
c
  18	if(np.eq.0) return
	do 19 k=1,np
	alph(k)=alpha(k)
  19	cvli(k)=alpha(k)
c
c
	do 29 k=1,mxpq
		kc=mxpq-k
		if(kc.ge.np) go to 24
		div=1.-alph(kc+1)*alph(kc+1)
		if(div.le.epsil2) go to 70
		if(kc.eq.0) go to 29
		do 23 j=1,kc
  23		alph(j)=(cvli(j)+alph(kc+1)*cvli(kc+1-j))/div
  24		if(kc.ge.nq) go to 26
		j1=max0(kc+1-np,1)
		do 25 j=j1,kc
  25		acf(j+1)=acf(j+1)+acf(kc+2)*alph(kc+1-j)
  26		if(kc.ge.np) go to 29
		do 27 j=1,kc
  27		cvli(j)=alph(j)
  29	continue
c
c
	acf(1)=.5*acf(1)
	do 33 k=1,mxpq
	if(k.gt.np) go to 33
	div=1.-alph(k)*alph(k)
		do 31 j=1,k+1
  31		cvli(j)=(acf(j)+alph(k)*acf(k+2-j))/div
		do 32 j=1,k+1
  32		acf(j)=cvli(j)
  33	continue
c
c
	do 43 i=1,ma
	miim1p=min0(i-1,np)
	if(miim1p.eq.0) go to 43
		do 42 j=1,miim1p
  42		acf(i)=acf(i)+alpha(j)*acf(i-j)
  43	continue
	acf(1)=acf(1)*2.
c
c
	cvli(1)=1.
	if(nq.le.0) go to 60
	do 53 k=1,nq
	cvli(k+1)=-beta(k)
	if(np.eq.0) go to 53
	mikp=min0(k,np)
		do 52 j=1,mikp
  52		cvli(k+1)=cvli(k+1)+alpha(j)*cvli(k+1-j)
  53	continue
c
  60	return
c
c
  70	ier=5
	return
	end

c
c********************************************************************
