c*************************************************************
c.....Frequently occurring functions used in defining B,C,D,E
c*************************************************************
      subroutine utils_scales(q2)
      implicit none
      double precision q2
      include 'scales2_inc.f'
      include 'functions_inc.f'
      fm2=log(xmuf2/q2)
      fmu2=log(xmur2/q2)
      end

      subroutine utils_fu(uh,q2)
      implicit none
      double precision uh,q2
      include 'scales2_inc.f'
      include 'functions_inc.f'
      fu=log(-uh/q2)
      fa=fu !log(a/q2)
      end
      
      subroutine utils_dilog(sh,th,uh,q2)
      implicit none
      double precision sh,th,uh,q2
      double precision li2qs
      double precision Li2
      external li2
      include 'functions_inc.f'

      li2qs=Li2(q2/sh)
      f1t=Li2(q2/(q2-th))+1d0/2d0*(log(q2/(q2-th)))**2
      f2t=li2qs+1d0/2d0*fs**2+fs*log(-th/(sh-q2))
      f1u=Li2(q2/(q2-uh))+1d0/2d0*(log(q2/(q2-uh)))**2
      f2u=li2qs+1d0/2d0*fs**2+fs*log(-uh/(sh-q2))
      
      return
      end
      
      subroutine utils(sh,th,uh,q2,ss2)
      implicit none
      double precision sh,th,uh,q2,ss2
      include 'scales2_inc.f'
      include 'functions_inc.f'
      include 'fodyqt_inc.f'
c      INTEGER :: nthreads, myid
c      INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS

c     pass s2 as a parameter to set it to 0 when appropriate, and improve numerical precision
      s2=ss2                    !s2=sh+th+uh-q2
      
c     denominator factors
      dt=1d0/(s2-th)
      du=1d0/(s2-uh)
      ds=1d0/(sh+q2-s2)
      dst=1d0/(sh+th-s2)
      dsu=1d0/(sh+uh-s2)
      dtu=1d0/(th*uh-s2*q2)
      la=sqrt((uh+th)**2-4d0*s2*q2)

c     transcendental functions
      fs=log(sh/q2)
      ft=log(-th/q2)
!     fu=log(-uh/q2)
!     fa=log(a/q2)

      if (s2.gt.0d0) then
         fs2=log(s2/q2)
      else
         fs2 = 0d0              !when s2 is 0, fs2 is infinity!
      endif
         
      fst=log(sh*th**2/(q2*(s2-th)**2))
      fsu=log(sh*uh**2/(q2*(s2-uh)**2))
      fla=log((sh+q2-s2+la)/(sh+q2-s2-la))
      fstu=log(sh*q2/((s2-th)*(s2-uh)))
      ftu=log((th*uh-s2*q2)/((s2-th)*(s2-uh)))
      flat=log(sh*q2*(s2-th)**2/(s2*(2*q2-uh)-q2*th)**2)
      flau=log(sh*q2*(s2-uh)**2/(s2*(2*q2-th)-q2*uh)**2)
!      print *,OMP_GET_THREAD_NUM(),'mand',sh,uh,th,q2,s2
!      print *,OMP_GET_THREAD_NUM(),'dens',dt,du,ds,dst,dsu,dtu
!      print *,OMP_GET_THREAD_NUM(),'trans1',fs,ft,fu,fs2,fm2,fmu2,fa
!      print *,OMP_GET_THREAD_NUM(),'trans2',fst,fsu,fla,fstu,ftu,flat,flau
      return
      end
