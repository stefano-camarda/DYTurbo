#include "dyres_interface.h"
#include "settings.h"
#include "interface.h"

void dyres::init()
{
    //Initialise some DYRES settings
    g_param_.g_param_ = opts.g_param;
    nnlo_.order_ = opts.order;            //order (0=LO, 1=NLO, 2=NNLO)
    opts_.fixedorder_  = opts.fixedorder; //fixed order/resummation switch
    qtsub_.xqtcut_= opts.xqtcut;          //Cut on qt/Q
    qtsub_.qtcut_= opts.qtcut;            //Cut on qt

    //move here the flaq.eq.0 initialisation part of resumm() in main2 instead of using this initialisation flag
    resummconst_();                      //Initialise beta and other QCD coefficients used in sudakov.f
    /*
!     ! ONE TIME INITIALIZATION
!print *, 'Warming up: first call to resum ... '
      write(*,'(A)',advance='no') 'First call to resum ... '
!write(*,'()',advance='no') 'First call to resum ... '
      call resummconst
      call reno2const           !this seems to be unused
      
C     Initialization for Gauss inversion
      if (approxpdf.eq.0) then
         call initmellingauss   ! higher order quadrature rule
      else
         call INITO             ! dyres original quadrature rule
      endif
      print *,'Done'
c     Initialization of redundant variables       
      flag1=order               ! set flag1 to order of calculation (carbon copy of order, 1=NLO+NLL, 2=NNLO+NNLL)
c     Choose pp or ppbar collider
      ih1=iih1                  !1
      ih2=iih2                  !1

c     Set factorization and renormalization scales (to work with dynamic scale need to move this outside init stage)
      mur=scale
      muf=facscale
C     Scales           
      mur2=mur**2
      muf2=muf**2

      q2mur=mur2
      q2muf=muf2

C     non-perturbative parameter
      g=g_param

c     precompute log of scales
c      loga=log(a_param)
c      rloga=log(a_param)

c     narrow width, no branching ratio, always 0
      fnwa=0
      brflag=0

c     choose real axis (complex plane) integration of bstar (b) (always 0)
c     flagrealcomplex=0

c     set flag5 = 3 for Z, 21 for W+, 22 for W- (if phot is on (-1) flag5 is 5)
      if(nproc.eq.3) then
         flag5=3-phot*2
      elseif(nproc.eq.1.or.nproc.eq.2) then
         flag5=nproc+20
      endif
      
c     narrow width approximation (never used)      
      if (fnwa.eq.1) then
         if (flag5.eq.3)  q=Mz
         if (flag5.eq.2.or.flag5.eq.21.or.flag5.eq.22) q=Mw
      endif

c     normal (imod=0) or modified (imod=1) sudakov      
c     imod=1

      nnf=nf                    !number of flavours (why duplicated?)
c      b0p=b0*a_param
      
c     check on remove branching ratio, not used
      if (brflag.eq.1.and.fnwa.eq.0) then
         write(*,*)"ERROR: Remove BR flag can be true, only if Narrow Width
     /Approximation flag is true!"
         stop
      endif

         
      if (mod.eq.0.and.approxpdf.eq.1) then
C     Call fitting routine
         CALL fiteador(x,muf2,etam)
C     CALL initialization subroutine for PDFS and Anomalous dimensions
         CALL INITOFIT
      else
C     initialization of exact PDF moments and anomalous dimensions
         call initmoments
      endif

C     define points for quadratures integration
      call intdeoini(lenaw, tiny, 1.0d-2, aw)

      return
      end
    */
}
