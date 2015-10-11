      subroutine a06(xb,q2,PDFS,DPDFS,NPDF,NPAR)
c
c     This is a code for the NNLO parton distributions 
c     in the variable-flavor-number (VFN) schem with account 
c     of their experimental and theoretical uncertainties. 
c     The Q**2 range is 0.8d0 < Q**2 < 2d8, the x range is 1d-7 < x < 1d0
c     (for the values of PDFs and strong coupling constant at Q2 < 0.8 GeV^2 
c     (x < 1d-7) their values at Q^2 = 0.8 GeV^2 (x = 1d-7) are returned).
c
c  Output parameters:
c     The array PDFS contains fitted values of the strong coupling constant 
c     and the parton distributions at given x and Q:
c        PDFS(0) -- \alpha_s
c        PDFS(1) -- valence u-quarks 
c        PDFS(2) -- valence d-quarks
c        PDFS(3) -- gluons 
c        PDFS(4) -- sea u-quarks 
c        PDFS(5) -- s-quarks 
c        PDFS(6) -- sea d-quarks 
c        PDFS(7) -- c-quarks
c        PDFS(8) -- b-quarks
c        PDFS(9) -- t-quarks
c     NPDF is the number of PDFs returned (NPDF=9 for the VFN scheme).
c     Output array DPDFS(0:npdf,npar) contains derivatives of \alpha_s and
c     the PDFs on the fitted parameters with the number of the parameters 
c     returned in NPAR. With the derivatives of \alpha_s included one can take 
c     into account the correlations of the fitted PDFs with \alpha_s as well.
c     All derivatives are transformed to the orthonormal 
c     basis of eigenvectors of the parameters error matrix. For this reason 
c     the variation of the PDFs in the derivatives directions can be performed 
c     independently. For example the dispersion of the i-th PDF can be stored 
c     in DELPDF using the code 
c
c-----------------
c          DELPDF=0.
c          do k=1,npar
c            DELPDF=DELPDF+dpdfs(i,k)**2
c          end do
c-----------------
c     and its random value can be stored in RPDF using the code 
c-----------------
c          RPDF=pdfs(i)          
c          do k=1,npar
c            s=0.
c            do l=1,96
c              s=s+(2*rndm(xxx)-1)/sqrt(32.)
c            end do
c            RPDF=RPDF+s*dpdfs(i,k)
c          end do
c-----------------
c          
c         Reference: Phys. Rev. D74, 054033 (2006) [hep-ph/0606237]
c      
c         Comments: Sergey.Alekhin@ihep.ru                      
c                                                               
c     Initial version: Dec 2006      

      implicit none

      integer nxb,nq,np,nvar
      parameter(nxb=99,nq=20,np=9,nvar=23)

      integer k,i,n,m,kx,nxbb
      integer NPDF,NPAR,KORD,KSCHEM

      integer nexp(0:np) 

      real*8 f(nxb,nq+1,0:np),xx(nxb)
      real*8 fsp(nxb),bs(nxb),cs(nxb),ds(nxb)
      real*8 bsp(nxb,nq+1,0:np),csp(nxb,nq+1,0:np),dsp(nxb,nq+1,0:np)
      real*8 bspd(nvar,nxb,nq+1,0:np),cspd(nvar,nxb,nq+1,0:np)
     ,      ,dspd(nvar,nxb,nq+1,0:np)
      real*8 pdfs(0:np),dpdfs(0:np,nvar)
      real*8 df(nvar,0:np,nxb,nq+1)
      real*8 x,qsq,dels,delx,x1,delx1,xlog1,xd,b,aa,ss,f0,fp,fm
      real*8 xb,q2,df0,dfp,dfm

      character pdford*1
      dimension pdford(3)
      character pdfschem*3
      dimension pdfschem(0:1)

      real*8 xmin,xmax,qsqmin,qsqmax
      integer kords,kschems

c I/O channel to read the data
      integer nport
      character locdir*2
      data nport/1/
c put in your local address of the PDFs files in LOCDIR
      data locdir /'./'/
      data pdford/'1','2','3'/
      data pdfschem /'ffn','vfn'/
      data nexp / 0, 3, 4, 5, 5, 5, 5, 5, 5, 5 /
      data xmin,xmax,qsqmin,qsqmax/1d-7,1d0,0.8d0,2d8/
      data KORDS,KSCHEMS /-1,-1/
      data KORD,KSCHEM /3,1/

      save kords,kschems,f,df,dels,delx,x1,delx1,xlog1,nxbb,xx

      if (kschem.eq.0) then 
        npdf=6
      else 
        npdf=9
      end if
      npar=nvar
*
*...Reset arrays
*
      do i=0,npdf
        pdfs(i) = 0.d0
        do k=1,npar
          dpdfs(i,k)= 0.d0
        end do
      end do
*
      if((kords.eq.kord).and.
     +   (kschems.eq.kschem)) goto 10

      kords=kord
      kschems=kschem

      dels=(dlog(dlog(qsqmax/0.04d0))-
     +      dlog(dlog(qsqmin/0.04d0)))/dble(nq-1)

      nxbb=nxb/2
      x1=0.3d0
      xlog1=dlog(x1)
      delx=(dlog(x1)-dlog(xmin))/dble(nxbb-1)
      DELX1=(1.d0-x1)**2/dble(nxbb+1)

*...X GRID
      do kx=1,nxbb
        xx(kx)=dexp(dlog(xmin)+delx*dble(kx-1))
      end do
      do kx=nxbb+1,nxb-1
        xx(kx)=1.d0-dsqrt(dabs((1.d0-x1)**2-delx1*dble(kx-nxbb)))
      end do
      xx(nxb)=1d0

*...Read input tables
      print *,'***** Reading PDFs from tables *****'
      open(unit=nport,status='old'
     ,    ,file=locdir//'Pdfdata/a06.dpdfs_'//pdford(kord)//'_'
     /                              //pdfschem(kschem))
      do n=1,nxb-1
        do m=1,nq
          do i=0,npdf 
            read (nport,*) (df(k,i,n,m),k=1,npar)
          end do
        end do
      end do
      close(unit=nport)

      do k=1,npar
        do i=0,npdf
          do m=1,nq
            if (i.ne.0) then 
              df(k,i,nxb,m)=0d0
            else 
              df(k,i,nxb,m)=df(k,i,nxb-1,m)
            end if
            do n=1,nxb
              fsp(n)=df(k,i,n,m)
            end do
            call spline (nxb,xx,fsp,bs,cs,ds)
            do n=1,nxb
              bspd(k,n,m,i)=bs(n)
              cspd(k,n,m,i)=cs(n)
              dspd(k,n,m,i)=ds(n)
            end do
          end do
        end do
      end do

      open(unit=nport,status='old',err=199
     ,    ,file=locdir//'Pdfdata/a06.pdfs_'//pdford(kord)//'_'
     /     //pdfschem(kschem))
      do n=1,nxb-1
        do m=1,nq
          read(nport,*) (f(n,m,i),i=0,npdf)
        end do
      end do
      do i=0,npdf
        do m=1,nq
          if (i.ne.0) then 
            f(nxb,m,i)=0d0
          else 
            f(nxb,m,i)=f(nxb-1,m,i)
          end if
          do n=1,nxb-1
            f(n,m,i)=f(n,m,i)/(1d0-xx(n))**nexp(i)
          end do
          do n=1,nxb
            fsp(n)=f(n,m,i)
          end do
          call spline (nxb,xx,fsp,bs,cs,ds)
          do n=1,nxb
            bsp(n,m,i)=bs(n)
            csp(n,m,i)=cs(n)
            dsp(n,m,i)=ds(n)
          end do
        end do
      end do
      close(unit=nport)

  10  continue

      if((q2.lt.qsqmin).or.(q2.gt.qsqmax)) print 99,q2,qsqmin,qsqmax
      if((xb.lt.xmin).or.(xb.gt.xmax)) print 98,xb,xmin,xmax

  99  format('  A06 WARNING:  Q^2 VALUE IS OUT OF RANGE   ',3g12.3)
  98  format('  A06 WARNING:   X  VALUE IS OUT OF RANGE   ',3g12.3)

      x=max(xb,xmin)
      x=min(x,xmax)
      qsq=max(q2,qsqmin)
      qsq=min(qsq,qsqmax)

      if (x.gt.x1) then
        xd=(1d0-x1)**2-(1d0-x)**2
        n=int(xd/delx1)+nxbb
      else
        xd=dlog(x)-xlog1
        n=nxbb+int(xd/DELX)-1
      end if
      aa=x-xx(n)

      ss=dlog(dlog(qsq/0.04d0))-dlog(dlog(qsqmin/0.04d0))
      m=int(ss/dels)+1
      b=ss/dels-dble(m)+1.d0

      do i=0,npdf
        f0=f(n,m,i) + aa*bsp(n,m,i) + aa**2*csp(n,m,i) 
     +              + aa**3*dsp(n,m,i)
        fp=f(n,m+1,i) + aa*bsp(n,m+1,i) + aa**2*csp(n,m+1,i)
     +                + aa**3*dsp(n,m+1,i)
        if (m.ge.2) then 
          fm=f(n,m-1,i) + aa*bsp(n,m-1,i) + aa**2*csp(n,m-1,i)
     +                   +aa**3*dsp(n,m-1,i)
          pdfs(i)=fm*b*(b-1d0)/2d0 + f0*(1d0-b**2) + fp*b*(b+1d0)/2d0
        else 
          pdfs(i)= f0*(1d0-b) + fp*b
        end if
        pdfs(i) = pdfs(i)*(1d0-x)**nexp(i)
        do k=1,npar
          df0=df(k,i,n,m) + aa*bspd(k,n,m,i) + aa**2*cspd(k,n,m,i) 
     +                    + aa**3*dspd(k,n,m,i)
          dfp=df(k,i,n,m+1)+aa*bspd(k,n,m+1,i)+aa**2*cspd(k,n,m+1,i)
     +                      + aa**3*dspd(k,n,m+1,i)
          if (m.ge.2) then 
            dfm=df(k,i,n,m-1)+aa*bspd(k,n,m-1,i)+aa**2*cspd(k,n,m-1,i)
     +                        + aa**3*dspd(k,n,m-1,i)
            dpdfs(i,k)=dfm*b*(b-1d0)/2d0 
     +                + df0*(1d0-b**2) +dfp*b*(b+1d0)/2d0
          else 
            dpdfs(i,k) = df0*(1d0-b) + dfp*b
          end if
        end do
      end do

      return

 199  print *,'The PDF set is inavailable (FILE:'
     ,     ,'a06.pdfs_'//pdford(kord)//'_'
     /     //pdfschem(kschem),')'

      return
      end
