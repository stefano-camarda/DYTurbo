//interface to Pegasus-QCD
//The interface is a translation to C++ of the following Pegasus QCD files:
//initevol.f nparton.f initinp.f inplmom1.f
  
//The files are translated and modified on order to use LHAPDF as input for the PDF moments
//in N-space, and to use a different contour in the complex plane for the support points of the
//Mellin inversion
  
#include "pegasus.h"

#include "pdfevol.h"
#include "resint.h"
#include "resconst.h"
#include "besselint.h"
#include "mellinint.h"
#include "settings.h"
#include "LHAPDF/LHAPDF.h"

int pegasus::nff;
int pegasus::ivfns;

void pegasus::init()
{
  //No need to initialise if pegasus is not used
  if (opts.evolmode != 1 && opts.evolmode != 3 && opts.evolmode != 4)
    return;
  
  // From:
  // ..File: initevol.f    
  //                   
  //
  // ..The global initialization for the unpolarized parton evolution. 
  //
  // ..This routine initializes some constants, the array  NA  of complex
  //    support points (depending on  IFAST),  the corresponding weights
  //    WN  for the Gauss quadrature employed for the Mellin inversion, 
  //    and the simple harmonic sums S_i, i=1,...,6 on the support points.
  //
  // ..The corresponding arrays of N-space splitting functions, evolution
  //    matrix elements and flavour-threshold operator matrix elements are 
  //    then initialized (depending on  IMODEV, NPORD, IVFNS, NFF, FR2) 
  //    by the respective subroutines and stored in common blocks by them. 
  //
  // ..The default values of the initialization parameters speficied below
  //    can by overwritten either by reading the file  usrinit.dat  (for 
  //    IPAR = 1)  or by calling the subroutine  USRINIT  (for IPAR = 2).
  //
  // =====================================================================

  const int NFMIN = 3;
  const int NFMAX = 6;


  //Some constants and the two-dimensional Kronecker symbol
  double EMC = 0.57721566490153;

  rzeta_.zeta_[0] = EMC;
  rzeta_.zeta_[1] = 1.644934066848226;
  rzeta_.zeta_[2] = 1.202056903159594;
  rzeta_.zeta_[3] = 1.082323233711138;
  rzeta_.zeta_[4] = 1.036927755143370;
  rzeta_.zeta_[5] = 1.017343061984449;

  kron2d_.d_[0][0] =  fcx(complex <double> (1., 0.));
  kron2d_.d_[1][0] =  fcx(complex <double> (0., 0.));
  kron2d_.d_[0][1] =  fcx(complex <double> (0., 0.));
  kron2d_.d_[1][1] =  fcx(complex <double> (1., 0.));
  
  // QCD colour factors
  colour_.ca_ = 3.;
  colour_.cf_ = 4./3.;
  colour_.tr_ = 0.5;

  // Some default settings of the external initialization parameters
  // (standard-speed iterated VFNS evolution at NLO for mu_f/mu_r = 1)

  //reproduce DYRES evolution
  if (opts.evolmode == 1)
    {
      ivfns = 0; //FFN evolution
      nff = 5;   //5 light flavours
      order_.npord_ = opts.order - 1; //order of evolution is LO(NLO) for NLL(NNLL)
    }
  //reproduce in N-space the LHAPDF evolution in x-space
  else if (opts.evolmode == 3 || opts.evolmode == 4)
    {
      ivfns = 1; //VFN evolution (read from LHAPDF)
      nff = 4;   //in FFN evolution, number of flavours  (read from LHAPDF)
      order_.npord_ = LHAPDF::getOrderPDF(); //order of evolution (read from LHAPDF)
    }
  
  //mode of evolution in Pegasus-QCD
  //there are three available schemes for solving the evolution equations in N-space at NNLO
  // mode = 1 reproduces evolution in x-space
  evmod_.imodev_ = evolution_mode; 
  
  double FR2 = 1.;//ratio of muren2/mufac2
  frrat_.logfr_ = log(FR2);

  nnused_.nmax_ = mellinint::mdim;

  //Location of the Mellin inversion contour in the complex N plane 
  //Support points NA(K) for the gauss integration
  for (int i = 0; i < mellinint::mdim; i++)
    {
      moms_.na_[i] = fcx(mellinint::Np[i]);

      //The lowest simple harmonic sums on the support points
      fcomplex N1 = fcx(mellinint::Np[i] + 1.);
      int one = 1;
      int two = 2;
      int three = 3;
      int four = 4;
      int five = 5;
      hsums_.s_[0][i] = fcx(cx(psi_(N1)) + EMC);
      hsums_.s_[1][i] = fcx(rzeta_.zeta_[1] - cx(dpsi_(N1,one)));
      hsums_.s_[2][i] = fcx(rzeta_.zeta_[2] + 0.5 * cx(dpsi_(N1,two)));
      hsums_.s_[3][i] = fcx(rzeta_.zeta_[3] - 1./6. * cx(dpsi_(N1,three)));
      hsums_.s_[4][i] = fcx(rzeta_.zeta_[4] + 1./24. * cx(dpsi_(N1,four)));
      hsums_.s_[5][i] = fcx(rzeta_.zeta_[5] - 1./120. * cx(dpsi_(N1,five)));
    }

  //Some more (internal) evolution and initialization parameters 
  aspar_.naord_  = order_.npord_;
  aspar_.nastps_ = alphas_steps; //number of steps for the Runge-Kutta integration of alphas beyond LO

  
  if (evmod_.imodev_ == 1 || evmod_.imodev_ == 2)
    itord_.nuord_ = 15;
  else
    itord_.nuord_ = order_.npord_;

  if (ivfns == 0)
    {
      nfused_.nflow_  = nff;
      nfused_.nfhigh_ = nff;
    }
  else
    {
      nfused_.nflow_  = NFMIN;
      nfused_.nfhigh_ = NFMAX;
    }

  //The N-space splitting functions, evolution operators etc.
  betafct_();
  pns0mom_();
  psg0mom_();
  lsgmom_();

  if (order_.npord_ > 0)
    {
      pns1mom_();
      psg1mom_();
      usg1mom_();
    }
  if (itord_.nuord_ > 1) usg1hmom_();

  if (order_.npord_ > 1)
    {
      pns2mom_();
      psg2mom_();
      uns2mom_();
      usg2mom_();
    }
  if (itord_.nuord_ > 2) usg2hmom_();

  if (ivfns != 0)
    {
      ans2mom_();
      asg2mom_();
    }

  // ..File nparton.f
  //
  //
  // ..The subroutine  NPARTON  for one (in general complex) Mellin moment 
  //    N  of the parton densities evolved to the scale  M2 = mu_f^2  by 
  //    EVNFFN  (for IVFNS = 0)  or  EVNVFN.  
  //
  // ..The results are returned as  f_i(N,mu_f^2)  by the array  PDFN.
  //    See the header of evnffn.f for the notation for f_i depending
  //    on the input parameter  IPSTD.  For  IPOL = 0  we are dealing
  //    with the unpolarized, otherwise with the polarized case.
  //
  // ..Note that, unlike its x-space counterpart  XPARTON  the present
  //    routine includes calls of  INITPOL (.., EPAR)  and  INITINP (IPAR) 
  //    or  INITPINP (IPAR),  i.e., for the N-space evolution this routine 
  //    is the only one the user needs to call.  Input and initialization
  //    parameters are specified as usual depending on  EPAR  and  IPAR.
  // 
  // =====================================================================

  //Initialize the input parameters and initial moments
  //
  // ..File: initinp.f
  //
  //
  // ..The input initialization for the unpolarized parton evolution.
  //
  // ..This routine initializes the parton distributions at the initial
  //    scale  M20 (in GeV^2)  and, in the VFNS case, at the flavour
  //    thresholds,  Mi2, i = C, B, T, together with the corresponding
  //    values of the strong coupling constant. The corresponding routines
  //    are  INPLMOMj  and  EVNFTHR  where more information can be found.
  //    So far  INPLMOM1,2  are available and selected by  NFORM = 1, 2. 
  //
  // ..The flavour parameters  IVFNS, NFF  and the (fixed) scale log
  //    LOGFR = ln(mu_f^2/mu_r^2)  are provided by the common blocks
  //    VARFLV, NFFIX  and  FRRAT.  If called with  IPAR  unequal 1 or 2,
  //    the initial conditions for the evolution are those of the 2001
  //    Les Houches benchmark.  Other inputs (including the choice of
  //    NFORM  are specified by reading the file  usrinp.f  (for IPAR = 1)
  //    of by calling the routine  USRINP  (for IPAR = 2).  
  //

  //Some default settings of scales, couplings and input distributions

  //Initial scale  M20 = M_0^2 (in GeV^2)  and  ASI = alpha_s(M_0^2)

  //backward evolution from the factorisation scale
  if (opts.evolmode == 1)
    asinp_.m20_ = pow(opts.kmufac*opts.rmass,2);

  //input values from LHAPDF and forward evolution
  else if (opts.evolmode == 3 || opts.evolmode == 4)
    {
      LHAPDF::PDFInfo info(opts.LHAPDFset, opts.LHAPDFmember);
      double qmin = info.get_entry_as<double>("QMin", -1);
      asinp_.m20_ = pow(qmin,2);
    }

  double ASI = LHAPDF::alphasPDF(sqrt(asinp_.m20_));

  //The heavy quark masses squared, input values from LHAPDF
  asfthr_.m2c_ = pow(LHAPDF::getThreshold(4),2);
  asfthr_.m2b_ = pow(LHAPDF::getThreshold(5),2);
  asfthr_.m2t_ = pow(LHAPDF::getThreshold(6),2);
  
  //Stop some nonsense
  if (ivfns == 1 && asinp_.m20_ > asfthr_.m2c_)
    {
      cout << "Too high mu_0 for VFNS evolution. STOP" << endl;
      exit (-1);
    }

  if (ASI > 2. || ASI < 2.0e-2 )
    {
      cout << "alpha_s out of range. STOP" << endl;
      exit (-1);
    }

  if (ivfns == 1 && asfthr_.m2c_ > asfthr_.m2b_)
    {
      cout << "wrong charm-bottom mass hierarchy. stop" << endl;
      exit (-1);
    }
  if (ivfns == 1 && asfthr_.m2b_ > asfthr_.m2t_)
    {
      cout << "Wrong bottom-top mass hierarchy. STOP" << endl;
      exit (-1);
    }

  //For mu_r unequal mu_f  AS0  is different from the input parameter 
  //ASM = a_s(M_0^2). The VFNS evolution starts with n_f = 3 at M_0^2.

  double ASM = ASI / (4.* M_PI);
  double R20 = asinp_.m20_ * exp(-frrat_.logfr_);
  int NF = nff;
  if (ivfns != 0)  NF  = 3;
  asinp_.as0_ = as_(R20, asinp_.m20_, ASM, NF);

  //  cout<< " alphas at the starting scale " << asinp_.as0_ << endl;
  //  cout<< " alphas LHAPDF " << LHAPDF::alphasPDF(sqrt(asinp_.m20_))/ (4.* M_PI) << endl;
  
  //Input initialization (including threshold values for the VFNS case)
  // ..File: inplmom1.f 
  //
  //
  // ..The subroutine  INPLMOM1  (version 1) for the initial light-parton 
  //    distributions of a non-strange hadron in N-space.  The moments are 
  //    stored in the common-block  PAINP for an external NDIM-dimensional
  //    array  NA  of complex moments specified in the common-block  MOMS.
  //
  // ..The functional form of the distributions for this version reads
  //
  //        xf = Nf * x^Pf(2) * (1 - x)^Pf(3)  
  //             * ( 1 + x^Pf(4) * Pf(5) + x * Pf(6) )                 (1)
  //
  //    for  f = UV = u - ubar,  DV = d - dbar,  DL = dbar - ubar, 
  //    LS = 2 * (dbar + ubar),  SM = s - sbar,  SS = s + sbar,   GL = g. 
  //
  // ..The normalization factors  NUV and NDV  are fixed by the respective 
  //    quark numbers given by  PUV(1)  and  PDV(1)  (i.e.  PUV(1) = 2 and
  //    PDV(1) = 1  for the proton).
  //    NDL = PDL(1);  for  IMOMIN = 0  also  NLS = PLS(1), NSS = PSS(1).
  //    For non-zero  IMOMIN,  PLS(1) and PSS(1)  represent the respective
  //    momentum fractions.  NGL  is fixed by the momentum sum  PGL(1)  of
  //    all partons.  For non-zero  ISSIMP,  SS = PSS(1) * LS  is used for 
  //    the strange sea instead of (1), and  SM  is set to zero. Otherwise
  //    the coefficient of x^PSM(4) is such that the first moment is zero.
  //
  // ..The normalization factors  Nf  and the momentum fractions  Af  are 
  //    written to the common-block PANORM.  The flag  IINNEW  in  INPNEW
  //    is set to '1' at the end of this routine.
  //

  //Begin of the Mellin-N loop 
  fcomplex uval,dval,ubar,dbar,s,sbar,glu,charm,bot;
  for (int i = 0; i < mellinint::mdim; i++)
    {
      int hadron = 1; //opts.ih1;
      fcomplex N = fcx(mellinint::Np[i]); //compute positive branch only, the negative branch is obtained by complex conjugation
      double facscale = sqrt(asinp_.m20_);
      pdfmoments_(hadron,facscale,N,uval,dval,ubar,dbar,s,sbar,glu,charm,bot);
//      cout << "pegasus init " << facscale << "  "
//	   << i << "  " << mellinint::Np[i] << "  "
//	   << cx(glu) << "  "
//	   << cx(uval) << "  " << cx(dval) << "  "
//	   << cx(ubar) << "  " << cx(dbar) << "  "
//	   << cx(s) << "  "
//	   << cx(charm) << "  " << cx(bot) << endl;
      
      painp_.gli_[i] = glu;
      
      //Arrays of non-singlet and singlet quark combinations for N_f = 3 (and 4 and 5)
      //defined as in Eq. (2.16) of hep-ph/0408244
      //Pegasus evoltion: input PDFs at the starting scale, with 3 flavours
      //reproduce DYRES evolution: input PDFs at the factorisation scale, with 5 flavours
      if (opts.evolmode == 1)
	{
	  complex <double> qp[5];
	  qp[0] = cx(uval) + 2.*cx(ubar);
	  qp[1] = cx(dval) + 2.*cx(dbar);
	  qp[2] = cx(s) + cx(sbar);
	  qp[3] = 2.*cx(charm);
	  qp[4] = 2.*cx(bot);
	  
	  complex <double> qm[5];
	  qm[0] = cx(uval);
	  qm[1] = cx(dval);
	  qm[2] = cx(s) - cx(sbar);
	  qm[3] = 0.;
	  qm[4] = 0.;

	  painp_.vai_[i] = fcx(qm[0]+qm[1]+qm[2]+qm[3]+qm[4]);
	  painp_.m3i_[i] = fcx(qm[0]-qm[1]);
	  painp_.m8i_[i] = fcx(qm[0]+qm[1]-2.*qm[2]);
	  hfpainp_.m15i_[i] = fcx(qm[0]+qm[1]+qm[2]-3.*qm[3]);
	  hfpainp_.m24i_[i] = fcx(qm[0]+qm[1]+qm[2]+qm[3]-4.*qm[4]);
  
	  painp_.sgi_[i] = fcx(qp[0]+qp[1]+qp[2]+qp[3]+qp[4]);
	  painp_.p3i_[i] = fcx(qp[0]-qp[1]);
	  painp_.p8i_[i] = fcx(qp[0]+qp[1]-2.*qp[2]);
	  hfpainp_.p15i_[i] = fcx(qp[0]+qp[1]+qp[2]-3.*qp[3]);
	  hfpainp_.p24i_[i] = fcx(qp[0]+qp[1]+qp[2]+qp[3]-4.*qp[4]);
	}
      else if (opts.evolmode == 3 || opts.evolmode == 4)
	{
	  painp_.vai_[i] = fcx(cx(uval) + cx(dval) + (cx(s)-cx(sbar)));
	  painp_.m3i_[i] = fcx(cx(uval) - cx(dval));
	  painp_.m8i_[i] = fcx(cx(painp_.vai_[i]) - 3.*(cx(s)-cx(sbar)));

	  painp_.sgi_[i] = fcx(cx(uval) + cx(dval) + 2.* (cx(dbar) + cx(ubar)) + (cx(s)+cx(sbar)));
	  painp_.p3i_[i] = fcx(cx(painp_.m3i_[i]) - 2.* (cx(dbar)-cx(ubar)));
	  painp_.p8i_[i] = fcx(cx(uval) + cx(dval) + 2.* (cx(dbar) + cx(ubar)) - 2.* (cx(s)+cx(sbar)));
	}
    }
 
  //In VFN evolution calculate PDFs at the mc, mb, mt thresholds
  if (ivfns != 0)
    evnfthr_(asfthr_.m2c_, asfthr_.m2b_, asfthr_.m2t_);
  //  else
  //    dyevnfthr_(asfthr_.m2c_, asfthr_.m2b_, asfthr_.m2t_);
}

void pegasus::evolve()
{
  //The values of alphas: ASF for the evolution are obtained
  //using the Pegasus alphas function. In order to use the LHAPDF alphas function
  //the code which calculates the evolution at the thresholds should be changed

  //Values of a_s and N_f for the call of EVNFFN or EVNVFN below

  double M2 = pow(fabs(pdfevol::qbstar),2);   //qbstar = b0/bstar (without a_param) --> The factorisation scale of the PDFs should always be b0/b(star), and never scaled by a_param

  //check to reproduce bb->Z peak as in the old plots
  //M2 = pow(fabs(pdfevol::bstarscale),2); //bstarscale = a*b0/bstar
  
  //check PDFs at the starting scale
  //  double M2 = asinp_.m20_*10;
  //  pdfevol::bstarscale = sqrt(M2);

  //R2 is the final scale of the evolution, corrected by the ratio of muren2/mufac2
  double R2  = M2 * exp(-frrat_.logfr_);

  //Reproduce the L -> L~ redefinition by a redefinition of the b0/b scale (see Eq. 17 of hep-ph/0508068)
  //With the L -> L~ modification, alphas reaches asymptotically
  //alphas(mures/mufac/muren?) when bscale -> inf, so that the PDF evolution is frozen at muf as upper scale.
  //This modification is required to restore the fixed order cross section upon qt-integration

  double M2tilde = M2 * resint::mures2 / (M2 + resint::mures2);
  //  double M2tilde = M2 * resint::mufac2 / (M2 + resint::mufac2);
  //  double M2tilde = M2 * resint::muren2 / (M2 + resint::muren2);

  double R2tilde  = M2tilde * exp(-frrat_.logfr_);

  //use the modification L -> L~
  M2 = M2tilde;
  R2 = R2tilde;

  //Further modification of the final scale to account for the fact that alphas(muren)
  //instead of alphas(mures) is used as starting alphas(mu0) in alphasl(nq2) (variable aass)
  //This procedure should hold also for evolmode = 3 (VFN) as far as muf=mur? -> not really working, used only in evolmode = 1
  double M2prime = M2tilde * resint::muren2/resint::mures2; //scale for differences between muren and mures
  double R2prime = M2prime * exp(-frrat_.logfr_);

  /*
  //use the modification alphas(mures) -> alphas(muren)
  M2 = M2prime;
  R2 = R2prime;
  */

  //  R2 = R2 * resint::mures2 / resint::muren2;
    
  double ASI, ASF;
  int NF;

  //DYRES evolution: fixed number of flavours = 5 evolution from the factorisation scale downward
  if (opts.evolmode == 1)
    {
      NF  = nff;
      //set the factorisation scale as initial scale for the evolution 
      //ASI = asinp_.as0_;

      //set the resummation scale as initial scale for the evolution 
      ASI = resint::alpqres;

      //to set muf as initial scale, ASF should be multiplied by an additional alphas(muf)/alphas(mures) factor (i.e. change alpqres to alpqfac also in ASF)
      //-> consider to do it so that ASF is a valid final scale for the VFN evolution
      
      //set the factorisation scale as initial scale for the evolution 
      //ASI = resint::alpqfac;

      /*
      //calculate alphas at the final scale with pegasus, including the modification alpha(mures) -> alpha(muren)
      //      double R20 = asinp_.m20_ * R2/M2;
      //      double AS0 = asinp_.as0_;
      //Start the running of alphas from the renormalisation scale, to mimick alphasl(nq2)
      double R20 = resint::muren2 * R2/M2;
      double AS0 = resint::alpqren;
      //hence use the modification of the scale which accounts for alphas(mures) -> alphas(muren)
      M2 = M2prime;
      R2 = R2prime;
      //alphas from Pegasus (it is identical to the DYRES alphas at LO/LL running, but differs at low scales at NLO/NLL running)
      ASF = as_(R2, R20, AS0, NF)/resint::alpqren * resint::alpqres;
      */
  
      //calculate alphas at the final scale with DYRES (which has the modification L -> L~)
      //DYRES also has a different solution for the running of alphas beyond 1-loop
      fcomplex scale2 = fcx(pow(pdfevol::bscale,2));
      ASF = fabs(cx(alphasl_(scale2))) * resint::alpqres; //DYRES LL(NLL) running of alphas

      //---> Pegasus and DYRES running of alphas are identical only at LO/LL. At NLO/NLL DYRES evolution is different (it is not a Runge-Kutta solution),
      //and has dependences on the resummation and renormalisation scales
      //evolutionmode = 1 allows to check the impact of this difference.
      //To use the DYRES running of alphas in the VFN evolution set evolmode = 4

      //      cout << pdfevol::bscale << " dyres " << fabs(cx(alphasl_(scale2))) * resint::alpqres << " pegasus " << as_(R2, R20, AS0, NF)/resint::alpqren * resint::alpqres << endl;
      
      //      cout << pdfevol::bscale << "  " << resint::_m << "  " <<  as_(R2, R20, asinp_.as0_, NF) << "  " << fabs(cx(alphasl_(scale2))) * resint::alpqres << "  " << LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI) << endl;
      //      cout << pdfevol::qbstar << "  " << pdfevol::SALP << "  " << log(ASI/ASF) << endl;
      //      cout << "pegasus blog: " << pdfevol::bscale << "  " << resconst::beta0*resint::alpqres*4.*log(resint::mures2/M2tilde) << endl;
    }

  //Normal Pegasus evolution: VFNS evolution from the starting scale upward
  else if (opts.evolmode == 3)
    {
      //use the modification alphas(mures) -> alphas(muren)
      //M2 = M2prime;
      //R2 = R2prime;
      
      if (M2 > asfthr_.m2t_) //If M2 = M2tilde than the scale is frozen at muf, and never goes above the top
	{
	  NF = 6;
	  double R2T = asfthr_.m2t_ * R2/M2;
	  ASI = asfthr_.ast_;
	  //ASF = LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI);
	  ASF = as_(R2, R2T, asfthr_.ast_, NF);
	}
      else if (M2 > asfthr_.m2b_)
	{
	  NF = 5;
	  double R2B = asfthr_.m2b_ * R2/M2;
	  ASI = asfthr_.asb_;

	  //above the bottom mass, calculate alphas at the final scale with DYRES (which has the modification L -> L~)
	  //	  fcomplex fscale2 = fcx(pow(pdfevol::bscale,2));
	  //	  ASF = fabs(cx(alphasl_(fscale2))) * resint::alpqres; //DYRES LL(NLL) running of alphas.
	  
	  //ASF =  LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI);
	  ASF =  as_(R2, R2B, asfthr_.asb_, NF);
	}
      else if (M2 > asfthr_.m2c_)
	{
	  NF = 4;
	  double R2C = asfthr_.m2c_ * R2/M2;
	  ASI = asfthr_.asc_;
	  //ASF =  LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI);
	  ASF =  as_(R2, R2C, asfthr_.asc_, NF);
	}
      else
	{
	  NF = 3;
	  double R20 = asinp_.m20_ * R2/M2;
	  ASI = asinp_.as0_;
	  //ASF =  LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI);
	  ASF =  as_(R2, R20, asinp_.as0_, NF);
	}
      //      ASF = ASF/resint::alpqren * resint::alpqfac; //--> this is very crude
      
      //In the forward evolution need to rescale the final alphas ASF to account for the fact that the backward PDF evolution
      //is done from Qres to b0/b starting from PDFs at muf, while the residual evolution
      //from muf to Qres is factorised in the (H?) coefficients
      ASF = ASF/resint::alpqres * resint::alpqfac;
      
      //As a consequence of the ASF rescaling, recompute NF and ASI
      if (ASF < asfthr_.ast_)
	{
	  NF = 6;
	  ASI = asfthr_.ast_;
	}
      else if (ASF < asfthr_.asb_)
	{
	  NF = 5;
	  ASI = asfthr_.asb_;
	}
      else if (ASF < asfthr_.asc_)
	{
	  NF = 4;
	  ASI = asfthr_.asc_;
	}
      else
	{
	  NF = 3;
	  ASI = asinp_.as0_;
	}
    }
  
  //VFN forward evolution mode in which ASF is evaluated by DYRES
  else if (opts.evolmode == 4)
    {
      //calculate alphas at the final scale with DYRES (which has the modification L -> L~)
      //DYRES also has a different solution for the running of alphas beyond 1-loop
      fcomplex scale2 = fcx(pow(pdfevol::bscale,2));
      ASF = fabs(cx(alphasl_(scale2))) * resint::alpqres; //DYRES LL(NLL) running of alphas
      
      //Rescale the final alphas ASF for the fact that the PDF evolution is done from Qres to b0/b
      //starting from PDFs at muf, while the residual evolution from muf to Qres is factorised in the (H?) coefficients
      ASF = ASF/resint::alpqres * resint::alpqfac; 
      if (ASF < asfthr_.ast_) //If M2 = M2tilde than the scale is frozen at muf, and never goes above the top
	{
	  NF = 6;
	  ASI = asfthr_.ast_;
	}
      else if (ASF < asfthr_.asb_)
	{
	  NF = 5;
	  ASI = asfthr_.asb_;
	}
      else if (ASF < asfthr_.asc_)
	{
	  NF = 4;
	  ASI = asfthr_.asc_;
	}
      else
	{
	  NF = 3;
	  ASI = asinp_.as0_;
	}
    }
  
  //Calculation of the moments of the parton densities and output
  int nlow = 1;
  int nhigh = mellinint::mdim;
  int IPSTD = 1;
  fcomplex PDFN[13][144];
  if (ivfns == 0)
    dyevnffn_(PDFN, ASI, ASF, NF, nlow, nhigh, IPSTD); //modified ffn evolution with charm and bottom at the starting scale
  else
    evnvfn_(PDFN, ASI, ASF, NF, nlow, nhigh, IPSTD);


  for (int i = 0; i < mellinint::mdim; i++)
    {
      complex <double> fx[11];
      for (int p = -MAXNF; p <= MAXNF; p++)
	fx[p+MAXNF] = cx(PDFN[p+6][i]);
      pdfevol::storemoments(i, fx);
    }
  

  //Compare Pegasus QCD evolution with direct Mellin transform
  /*
  for (int i = 0; i < mellinint::mdim; i++)
    {
      pdfevol::calculate (i);
      cout << "direct " << fabs(pdfevol::bstarscale) << "  "
  	   << i << "  " << mellinint::Np[i] << "  "
  	   << cx(creno_.cfx1_[i][0+MAXNF]) << "  "
  	   << cx(creno_.cfx1_[i][1+MAXNF])-cx(creno_.cfx1_[i][-1+MAXNF]) << "  "
  	   << cx(creno_.cfx1_[i][2+MAXNF])-cx(creno_.cfx1_[i][-2+MAXNF]) << "  "
  	   << cx(creno_.cfx1_[i][-1+MAXNF]) << "  "
  	   << cx(creno_.cfx1_[i][-2+MAXNF]) << "  "
  	   << cx(creno_.cfx1_[i][3+MAXNF]) << "  "
	   << cx(creno_.cfx1_[i][4+MAXNF]) << "  "
	   << cx(creno_.cfx1_[i][5+MAXNF]) << endl;
  
      cout << "pegasus " << sqrt(M2) << "  "
  	   << i << "  " << mellinint::Np[i] << "  "
  	   << cx(PDFN[0+6][i]) << "  "
  	   << cx(PDFN[1+6][i])-cx(PDFN[-1+6][i]) << "  "
  	   << cx(PDFN[2+6][i])-cx(PDFN[-2+6][i]) << "  "
  	   << cx(PDFN[-1+6][i]) << "  "
  	   << cx(PDFN[-2+6][i]) << "  "
  	   << cx(PDFN[3+6][i]) << "  "
       	   << cx(PDFN[4+6][i]) << "  "
       	   << cx(PDFN[5+6][i]) << endl;
    }
  */
}
