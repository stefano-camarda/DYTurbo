      function alphasl(nq2)
c.....reference scale is factorization scale: muf2=muf**2
      implicit real*8(a-h,o-z)
      double complex xlambda,aa1,all,alphasl,qq,t,xlt,bstar,b,blog
      double complex log1xlambda
      double complex nq2,aa2
      include 'const.h'
      include 'sudakov_inc.f'
      include 'scales_inc.f'

c.....Here computes NLL expression for alphas
c     here nq2=b0^2/b^2 and the result is now  alpha(nq2)/alpha(Qres)  

c     To understand these formulas:
c     nq2 = a*b0^2/b^2 = a*qb2
c     Q2 = (m_ll/a)^2
c     alphas(qb2) = alphas(Q2) / (1 - beta0 * alphas(mur2) * log (Q2/qb2)) (why there is a mismatch between Q2 and mur2?)
c     with the changes:
c     IR cut off: b=b0p/nq -> bstar = b/sqrt(1+(b**2)/(blim**2))
c     modified sudakov: loq(Q2/qb2) -> loq(Q2/qb2 + 1)
c     The result is actually alphas(qb2)/alphas(Q2), where Q2 is the resummation scale

      
c     HERE CHANGE: order of alphas related to order of evolution
      if(iord.eq.1) then
      xlp=1
      elseif(iord.eq.0) then
      xlp=0
      endif

      b=sqrt((b0p**2/nq2))
      blim=rblim

      bstar=b

c.....choose bstar (b) for real axis (complex plane) integration
      if (flagrealcomplex.eq.0) bstar=b/sqrt(1+(b**2)/(blim**2))
      if (imod.eq.1) blog=log( (q*bstar/b0p)**2 + 1) !modified sudakov
      if (imod.eq.0) blog= log( (q*bstar/b0p)**2 )    !normal sudakov

      xlambda=beta0*aass*blog

c     I think it would be more correct to calculate alphas at the resummation scale, rather than aass which is alphas at the renormalisation scale
c     xlambda=beta0*dyalphas_lhapdf(q/a_param)/pi*blog
c     --> Not really, the running of alphas in this function reflects the definition of lambda in Eq. 25 of hep-ph/0508068.
      
c     print *,aass*pi,dyalphas_lhapdf(q),dyalphas_lhapdf(q/a_param),sqrt(nq2),xlambda
      
c     HERE now a dependence (without constant term)!
      log1xlambda=log(1-xlambda)
      aa1=log1xlambda+aass*xlp*
     .      (beta1/beta0*log1xlambda/(1-xlambda) 
     .       + beta0*xlambda/(1-xlambda)*rlogq2mur2
c     .       +   beta0*log(q2/muf2)
     .       -2d0*beta0*xlambda*rloga/(1-xlambda)   )
      alphasl=Exp(-aa1)
!      write(*,*) iord,alp,b,blim,flagrealcomplex,xlambda,as,a_param
!     .,mur,blog,nq2,aa1,alphasl
      
c.....Now compute the factors
c.....needed to resum the logs which multiply the N-dependent part
c.....of the C coefficients
c  the limit below implies xlambda<1/2 and then aa2<= 1      
c      blim=b0p*(1/q)*exp(1/(2*as*beta0)) 
c      blim=b0p*(1/q)*exp(1/(4*as*beta0))
C Set a limit to avoid very large values of b (= very small scales ~1/b)
!       blim=b0p*(1/q)*exp(1/(2*aass*beta0)) ! avoid Landau pole     
!       write(*,*) "blim",blim
      blim=0.5d0

      if (flagrealcomplex.eq.0) bstar=b/sqrt(1+(b**2)/(blim**2))
      if (imod.eq.1) blog=log( (q*bstar/b0p)**2 + 1) !modified sudakov
      if (imod.eq.0) blog= log( (q*bstar/b0p)**2 )    !normal sudakov

      xlambda=beta0*aass*blog

c     HERE now a dependence (without constant term)!
      log1xlambda=log(1-xlambda)
      aa1=log1xlambda+aass*xlp*
     .      (beta1/beta0*log1xlambda/(1-xlambda) 
     .       + beta0*xlambda/(1-xlambda)*rlogq2mur2
     .       -2d0*beta0*xlambda*rloga/(1-xlambda)   )
      aexp=Exp(-aa1)      

      aa2= xlambda/(1- xlambda)
      aexpB=Exp(aa2) 
      
      return
      end
      
c.....Sudakov form factor
      function S(b)
      implicit none
      complex *16 S,f0,f1,f2,b,bstar,blim,blog
      integer flag1
      COMMON/flag1/flag1
      real*8 g
      common/NP/g
      complex *16 y
      include 'sudakov_inc.f'
      include 'scales_inc.f'
      include 'const.h'
!      blim=(1/q)*exp(1/(2*aass*beta0))
c****************************
c mass dependence in blim
      blim=cblim

c     In reading these formulas, notice that L = q*bstar/b0p = (q/a_param)*bstar/b0 = Q * bstar/b0, according to Eq. (13) and (17) of hep-ph/0508068.
      
      bstar=b
c.....choose bstar (b) for real axis (complex plane) integration
c mass dependence in bstar
      if (flagrealcomplex.eq.0) bstar=b/sqrt(1+(b**2)/(blim**2))

c mass dependence in blog
      if (imod.eq.1) blog=log( (q*bstar/b0p)**2+1) !modified sudakov
      if (imod.eq.0) blog= log( (q*bstar/b0p)**2 )    !normal sudakov

c mass dependence in f0(y), f1(y), f2(y)
      y = beta0*aass*blog
      log1y=log(1-y)
      if (flag1.eq.0) then
         S=exp(blog*f0(y))
      elseif (flag1.eq.1) then
         S=exp(blog*f0(y)+f1(y))
      elseif (flag1.eq.2) then
         S=exp(blog*f0(y)+f1(y)+aass*f2(y))
      endif
!
      if(DBLE(S).gt.1d2) then
      write(*,*) "WARNING! LARGE SUDAKOV, S(b)=",S,"; for bstar=",bstar
      S=cmplx(0d0,0d0)
      endif
!
      S=S*exp(-g*b**2)

      return
c****************************
      end
      
      
c.....Soft-gluon-Resummation of LL
      function f0(y)
      implicit none
      complex *16 f0,y
      include 'const.h'
      include 'sudakov_inc.f'
      f0=(A1q/beta0)*(y+log1y)/(y)
      return
      end
      
c.....Soft-gluon-Resummation of NLL
c.....Now we have mu_r dependence!
      function f1(y)
      implicit none
      complex *16 f1,y
      include 'const.h'
      include 'sudakov_inc.f'
      f1=((A1q*beta1)/(beta0**3))*((1d0/2)*log1y*log1y +
     \     (y)/(1-y)+log1y/(1-y)) -
     \     (A2q/(beta0**2))*(log1y+(y)/(1-y)) + 
     \     (B1q/beta0)*log1y +
     \     (A1q/beta0)*(y/(1-y)+log1y)*rlogq2mur2
c    a dependence      
      f1=f1-2*rloga*A1q/beta0*y/(1-y)
      return
      end

c.....Soft-gluon-Resummation of NNLL
c.....Now we have mu_r dependence!
      function f2(y)
      implicit none
      complex *16 f2,y
      include 'const.h'
      include 'sudakov_inc.f'
      f2=((A2q*beta1)/(beta0**3))*((y/2)*((3*y-2d0)/(1-y)**2)-
     \     ((1-2*y)*log1y/(1-y)/(1-y))) - 
     \     (B2q/beta0)*((y)/(1-y))+
     \     (B1q*beta1/beta0**2)*((y)/(1-y)+log1y/(1-y))-
     \     (A3q/2/beta0**2)*(y)*(y)/(1-y)/(1-y) +
     \     A1q*((beta1**2/2/beta0**4)*(1-2*y)/(1-y)
     \     /(1-y)*log1y*log1y +
     \     log1y*((beta0*beta2-beta1**2)/(beta0**4)+
     \     beta1**2/beta0**4/(1-y)) +
c    \     beta1**4/beta0**4/(1-y)) +
     \     (y)/(2*beta0**4*(1-y)*(1-y))*
     \     (beta0*beta2*(2d0-3*y)+beta1**2*y)) -
     \     (A1q/2)*(y)*(y)/(1-y)/(1-y)*rlogq2mur2*rlogq2mur2 +
     \     rlogq2mur2*(B1q*y/(1-y)+A2q/beta0*y*y/(1-y)/(1-y)+
     \     A1q*beta1/beta0**2*(y/(1-y)+(1-2*y)/(1-y)/(1-y)*log1y)) 
     \     +2d0*C1qqn*((y)/(1-y))
c    a dependence (now without constant term)
      f2=f2+2*A1q*y*(y-2)/(1-y)**2*rloga**2-rloga
     \ *(2*B1q*y/(1-y)+2*y/beta0*A2q/(1-y)**2
     \ -2*A1q*beta1/beta0**2*y*log1y/(1-y)**2)
!
     \ + A1q*rloga*rlogq2mur2*y*2d0/(1-y)**2   
      return
      end