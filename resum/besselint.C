#include "resconst.h"
#include "settings.h"
#include "besselint.h"
#include "pdfevol.h"
#include "pegasus.h"
#include "mesq.h"
#include "hcoefficients.h"
#include "hcoeff.h"
#include "expc.h"
#include "muf.h"
#include "phasespace.h"
#include "gaussrules.h"
#include "resint.h"
#include "isnan.h"
#include "pdf.h"
#include "sudakovff.h"
#include "npff.h"
#include "ccoeff.h"
#include "gint.h"

#include <LHAPDF/LHAPDF.h>

double C1 = 1;
double C3 = 1;

//fortran interface
void besselint_(double &b, double &qt, double &q2)
{
  resint::_qt = qt;
  resint::_m = sqrt(q2);
  besselint::bint(b);
};

complex <double> besselint::bint(complex <double> b)
{
  if (b == 0.)
    //if (opts.modlog && abs(b) < 1e-8)
    return 1.; //Should be correct, since J_0(0) = 1, and S(0) = 1 (Sudakov). May be different in the case that the L~ = L+1 matching is not used.
  
  double qt = resint::_qt; //better take this from phasespace::qt

  /***************************************************** Moved to pdfevol::scales (begin)
  //pdfevol::bscale is used for
  //alpq (C1)
  //pdfevol::alpr (C3)
  //pdfevol::XL (C3)
  //aexp (C1)
  //it also used in evolmode 1 and 4 as factorisation scale
  

  //The soft scale of the resummation is b0/b. The factor a = mll/mures is introduced because the b scale is used in alphasl
  //Introduce two scales:
  //mub = C1*b0/b is the lower integration limit of the Sudakov integral
  complex <double> mub = resint::a * resconst::b0/b * C1;
  //pdfevol::bscale = C3*b0/b is the "soft" factorisation scale at which the Wilson coefficient functions are evaluated (see arxiv:1309.1393 for details)
  pdfevol::bscale = resint::a * resconst::b0/b * C3;
  
  //convert to fortran complex number
  fcomplex fscale2_mub = fcx(pow(mub,2));
  fcomplex fscale2_mufb = fcx(pow(pdfevol::bscale,2));

  //freeze PDF evolution below a certain scale
  //  if (pdfevol::bscale < 5.)
  //    pdfevol::bscale = 5.;

  //pdfevol::bstarscale and pdfevol::bstartilde are not used (be careful, they do no have the C3 scale)

  //bstarscale = a*b0/bstar (final scale used for the PDF evolution)
  double bstar = real(b / sqrt(1.+(b*b)/(blimit_.rblim_*blimit_.rblim_)));
  //double bstar = real(b / sqrt(1.+(b*b)/(1.1229190*1.1229190)));  //// -->hard-coded!! --> should allow a different blim in the PDF evolution as a setting
  pdfevol::bstarscale = resconst::b0*resint::a/bstar;

  //qbstar is used in evolmode 3
  //qbstar = b0/bstar, it corresponds to pdfevol::bsstarcale but without a_param
  pdfevol::qbstar = resconst::b0/bstar*C3;

  //simulate pythia ISR factorisation scale, which is muf^2 = qt^2 (i.e. avoid the bstar prescription to freeze the factorisation scale)
  //pdfevol::qbstar = resconst::b0/b*C3;
  
  //bstartilde is bstarscale (qbstar) with the modification L -> L~, which freezes the scale at muf
  //bstartilde is used in evolmode 2, for the direct mellin transfrom at each scale
  //pdfevol::bstartilde = pdfevol::bstarscale * resint::mufac / sqrt((pow(pdfevol::bstarscale,2) + resint::mufac2));
  if (opts.modlog)
    pdfevol::bstartilde = pdfevol::qbstar * resint::mures / sqrt((pow(pdfevol::qbstar,2) + resint::mures2));
  else
    pdfevol::bstartilde = pdfevol::qbstar;


  //complex PDF scale to be used for the minimal prescription
  pdfevol::bcomplex = resconst::b0/b;
  
  //  cout << b << "  " << bstar << "  " << blimit_.rblim_ << endl;
  *****************************************************/ // Moved to pdfevol::scales (end)

  
  //The integration from b to qt space is done with the bstar prescription (real axis in the b space), and use the bessel function
  //The (complex) integration in the minimal prescription would require hankel functions
  //********************
  //qt and b dependence (bessel function) (is xj0 a jacobian? Probably yes, the Jacobian for the change of variable from cartesian to radial coordinate, which translates from Fourier to Hankel transform)
  complex <double> xj0;
  double qtb = qt*real(b);
  if (resint::_mode == 3 || resint::_mode == 4) //qt-integrated mode --> Needs to be implemented also for the other prescriptions
    {
      //bstar prescription
      if (opts.bprescription == 0)
	xj0 = 2.*fort_besj1_(qtb)/b;
      //Integrate up to blim with besche
      else if (opts.bprescription == 1)
	xj0 = 2./b;
      //Minimal prescription
      else if (opts.bprescription == 2 || opts.bprescription == 3)
	xj0 = 2./b;
      //Local bstar prescription
      else if (opts.bprescription == 4)
	xj0 = 2.*fort_besj1_(qtb)/b;
    }
  else
    {
      //bstar prescription
      if (opts.bprescription == 0)
	xj0 = 2.*fort_besj0_(qtb);
	//xj0 = 2.*fort_besj1_(qtb)/real(b); //--> this line is to check the qt-integral as a function of qt used in mode 3
      //Integrate up to blim with besche
      else if (opts.bprescription == 1)
	xj0 = 2.;
      //Minimal prescription
      else if (opts.bprescription == 2 || opts.bprescription == 3)
	xj0 = 2.;
      //Local bstar prescription
      else if (opts.bprescription == 4)
	xj0 = 2.*fort_besj0_(qtb);
    }
  //********************

  //The Sudakov is mass and b dependent

  //numerical Sudakov
  if (opts.numsud || opts.numexpc)
    gint::calc(b);
  
  //fortran
  //fcomplex fb = fcx(b);
  //complex <double> sudak=cx(s_(fb));

  //C++
  complex <double> sudak=sudakov::sff(b);
  if (sudak == 0.)
    return 0.;

  //non-perturbative form factor
  //sudak *= exp(-opts.g_param*pow(b,2));
  complex <double> ff = npff::S(b,resint::_m,resint::x1,resint::x2);
  if (!opts.flavour_kt)
    sudak *= ff;

  //********************
  //b, qt and mass dependence
  complex <double> factorfin = b*xj0*sudak;
  //cout << "where is the nan " << b << "  " << xj0 << "  " << sudak << endl;
  //********************

  //Do not need the Mellin transform for the LL case, because the HN coefficient is 1 --> Use PDFs in x space
  if (opts.order == 0 && opts.xspace && (opts.evolmode == 0 || opts.evolmode == 1)) //if (opts.order == 0 && opts.xspace && (opts.evolmode == 1))
    return factorfin;
  
  
  //**************************************
  //b-dependence
  // Set scales for evolution in pdfevol
  //alphasl gives the LL/NLL/NNLL evolution of alphas from Qres=Q/a_param to q2=b0^2/b^2
  pdfevol::alphasl(b);
  pdfevol::scales(b);

  /********************************************
  //pdfevol::alpqf is used as starting scale of the PDF evolution in pdfevol::evolution
  //according to Eq. 42 of arXiv:hep-ph/0508068 it should be the factorisation scale,
  //but in Eq. 98 in the resummation scale is used
  double alpqf = resint::alpqres;              //alphas at resummation scale   (as in dyres)
  //double alpqf = resint::alpqfac;            //alphas at factorisation scale (as in Eq. 42 of arXiv:hep-ph/0508068, but see also Eq. 98)
  ********************************************/
  
  /********************************************
  //alpq is used in hcoefficients::calcb, it is alpq = alphas(b0^2/b^2)
  //it is used only at NLL, at NNLL instead aexp and aexpb are used (aexp is the same as alphasl, but with a different blim)
  complex <double> alpq;
  alpq = resint::alpqres * cx(alphasl_(fscale2_mub)); //--> Attention! in DYRES it is alphas(qres), in DyQT it is alphas(mur)
  if (opts.evolmode == 2 || opts.evolmode == 3)
    //In order to have variable flavour evolution, use here a VFN definition of alpq
    //There is possibly an issue here when the renormalisation scale is (very) different from mll, since aass=alphas(muren) is used in xlambda = beta0*aass*blog
    {
      //      double M2 = pow(fabs(pdfevol::bstartilde),2);
      double M2 = pow(fabs(pdfevol::qbstar),2);                     //qbstar = b0/bstar (without a_param)
      double M2tilde;
      if (opts.modlog)
	M2tilde = M2 * resint::mures2 / (M2 + resint::mures2); //modified logs L -> L~
      else
	M2tilde = M2;
      double M2prime = M2tilde * resint::muren2/resint::mures2;     //scale for differences between muren and mures
      M2 = M2prime;
      double R2 = M2;
      double ASI;
      int NF;
      alpq = pegasus::alphas(M2, R2, ASI, NF);

      //alpq = LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI);


//      //tweak this to look as in dyres, where alphas(muren) is used:
//      double R20 = resint::_m;
//      alpq = as_(M2, R20, resint::alpqren, NF);
//      alpq = alpq/resint::alpqren*resint::alpqres;

      
      //Based on the definition of aexp, may be need to rescale alpq for as(qres)/as(muren)?
      //alpq = alpq/resint::alpqren*resint::alpqres;
      //cout << "alpq " << alpq << " M2 " << M2 << endl;
    }

  //  complex <double> alpq = LHAPDF::alphasPDF(fabs(pdfevol::bstartilde))/4./M_PI;        //alpq = alphas(b0^2/b^2)
  //  cout << pdfevol::bstarscale << "  " << resint::alpqres * cx(alphasl_(fscale2)) << "  " <<  LHAPDF::alphasPDF(fabs(pdfevol::bstartilde))/4./M_PI << endl;

  ********************************************

  //pdfevol::XL = pdfevol::alpqf / alpq; // = 1./cx(alphasl_(fscale2));

  double alpqf = resint::alpqres;              //alphas at resummation scale   (as in dyres)
  fcomplex fscale2_mufb = fcx(pow(resint::a * resconst::b0/b,2));


  pdfevol::XL = 1./cx(alphasl_(fscale2_mufb)); //XL = alphas(mures2)/alphas(b0^2/b^2)
  pdfevol::XL1 = 1.- pdfevol::XL;
  pdfevol::SALP = log(pdfevol::XL);
  //cout << pdfevol::bscale << "  " << pdfevol::XL << endl;
  // SELECT ORDER FOR EVOLUTION LO/NLO
  pdfevol::alpr = alpqf * cx(alphasl_(fscale2_mufb))*(double)(opts.order-1);
  //force LO evolution
  //pdfevol::alpr = alpqf * cx(alphasl_(fscale2_mufb))*(double)(0);
  //  cout << "old " << b << "  " << pdfevol::SALP << "  " << log(1./cx(alphasl_(fscale2_mufb))) << "  " << pdfevol::alpr << "  " << pdfevol::alpq <<  endl;

  **************************************/
  //Perform PDF evolution
  pdfevol::evolution();


  /*
  //original dyres evolution: Perform PDF evolution from muf to the scale b0/b ~ pt
  //the scales used in the evolution correspond to SALP and alpr
  if (opts.evolmode == 0)
    pdfevol::evolution();
  
  //PDF evolution with Pegasus QCD from muf to the scale b0/b ~ pt
  else if (opts.evolmode == 1)
    pegasus::evolve();

  //Calculate PDF moments by direct Mellin transformation at each value of bstarscale = b0p/bstar
  else if (opts.evolmode == 2)
    pdfevol::calculate();
  
  //PDF evolution with Pegasus QCD from the starting scale Q20, in VFN
  else if (opts.evolmode == 3)
    pegasus::evolve();

  //PDF evolution with Pegasus QCD from the starting scale Q20, in VFN, and using alphasl(nq2) to evaluate ASF at the final scale
  else if (opts.evolmode == 4)
    pegasus::evolve();
  */

  //Truncate moments from xmin to 1, with xmin = m/sqrt(s) e^-y0 (currently works only at LL where the HN coefficient is 1)
  //if (opts.mellin1d)
  //{
      //  pdfevol::truncate();
      //  ccoeff::truncate();
  //}
  //pdfevol::uppertruncate();
  
  //  for (int i = 0; i < mellinint::mdim; i++)
  //    cout << "C++ " << b << "  " << i << "  " << cx(creno_.cfx1_[i][5]) << endl;
  //**************************************

  /*  
  //aexp and aexpb are calculated in alphasl, they are used in hcoefficients::calcb only for the NNLL cross section
  alphasl_(fscale2_mub);
  complex <double> aexp = cx(exponent_.aexp_); //aexp is actually the ratio alphas(a*b0/b)/alphas(muren)

  //calculate aexp = alphas(a*b0/b)/alphas(muren) using LHAPDF
  //if (opts.evolmode == 2 || opts.evolmode == 3)
  //aexp = LHAPDF::alphasPDF(real(pdfevol::bstartilde))/LHAPDF::alphasPDF(resint::muren);
  
//  //In order to have variable flavour evolution, use here a VFN definition of aexp (aexpB instead, should be independent from NF)
//  //It seems that the pt-integrated cross section is invariant for variations of a_param if aexp is the ratio of alpq/as(muren) (and not alpq/as(Qres))
//  //Indeed, aexp is always multiplied by aass, which is alphas at the renormalisation scale
//  if (opts.evolmode == 2 || opts.evolmode == 3) //account for VFN evolution
//    aexp = alpq / resint::alpqren;
//    //aexp = alpq / resint::alpqres;
  
  //cout << pdfevol::bstartilde << "  " << cx(exponent_.aexp_) << "  " << alpq / resint::alpqres << "  " << alpq / resint::alpqren << endl;
  complex <double> aexpb = cx(exponent_.aexpb_);
  */

  //In case of truncation of hcoeff, need to recalculate the original coefficients before truncation
  //if (opts.mellin1d)
  //hcoeff::calc(resint::aass,resint::logmuf2q2,resint::logq2muf2,resint::logq2mur2,resint::loga);
  
  // Cache the positive and negative branch of coefficients which depend only on one I index
  expc::reset();
  if (opts.mellin1d)
    {
      expc::calc(b);
      if (opts.mufevol)
	{
	  muf::reset();
	  muf::calc(b);
	}
      //expc::reset();
      //hcoeff::calcb(resint::aass,resint::logmuf2q2,resint::loga,pdfevol::alpq,expc::aexp,expc::aexpB); // --> Need to access aass,logmuf2q2,loga,alpq,aexp,aexpb
    }
  else
    //--> Implement 2D expc here 
    {
      //complex <double> aexpb = cx(exponent_.aexpb_);
      //complex <double> aexp = cx(exponent_.aexp_); //aexp is actually the ratio alphas(a*b0/b)/alphas(muren)
      //hcoefficients::calcb(resint::aass,resint::logmuf2q2,resint::loga,pdfevol::alpq,expc::aexp,expc::aexpB); // --> Need to access aass,logmuf2q2,loga,alpq,aexp,aexpb

      expc::calc(b);
      //expc::reset();
      //hcoefficients::calcb(resint::aass,resint::logmuf2q2,resint::loga,pdfevol::alpq,expc::aexp,expc::aexpB); // --> Need to access aass,logmuf2q2,loga,alpq,aexp,aexpb
    }      
    
  //Inverting the HN coefficients from N to z space does not work, because it would miss the convolution with PDFs
  //double q2 = pow(phasespace::m,2);
  //if (opts.mellin1d)
  //hcoeff::invert(q2);

  //Truncate HN coefficients, to implement mellin1d also for rapdity intervals (currently not working, not sure if it is eventually possible. Need to figure out how to truncate the moments)
  //if (opts.mellin1d)
  //hcoeff::truncate();

  //Apply flavour dependent non-perturbative form factors
  if (opts.flavour_kt)
    pdfevol::flavour_kt();

      
  complex <double> fun = 0.;

  mellinint::reset();
    
  //Do not need the Mellin transform for the LL case, because the HN coefficient is 1 --> Use PDFs in x space
  if (opts.order == 0 && opts.xspace && (opts.evolmode == 2 || opts.evolmode == 3 || opts.evolmode == 4))
    //if (opts.order == 0 && opts.xspace && (opts.evolmode == 0 || opts.evolmode == 2 || opts.evolmode == 3 || opts.evolmode == 4))
    {
      double muf;
      
      //double muf = resconst::b0/b;
//      if (opts.evolmode == 2)
//	muf = fabs(pdfevol::qbstar);
//      if (opts.evolmode == 4)
//	muf = fabs(pdfevol::qbstar);
      //muf = fabs(pdfevol::bscale);

      if (opts.evolmode == 0 || opts.evolmode == 1)
	muf = resint::mufac;
      else
	muf = real(pdfevol::bstartilde);


      if (resint::_mode < 2 || resint::_mode == 4) //rapidity differential
	fun = mesq::loxs(resint::x1, resint::x2, muf);
      else                   //rapidity integrated
	fun = mesq::loxs(resint::tau, muf);      

      //cout << "besselint " << fun << endl;
    }
  else
    {
      //1d mellin
      if (opts.mellin1d)
	fun = mellinint::calc1d();
//	{
//	  //double q2 = resint::_m*resint::_m;
//	  //double bjx= q2/pow(opts.sroot,2);
//	  //double ax = log(bjx);
//	  for (int i = 0; i < mellinint::mdim; i++)
//	    {
//	      pdfevol::retrieve1d(i::positive);
//	      mellinint::pdf_mesq_expy(i,i,mesq::positive);
//	      double int1 = mellinint::integrand1d(i);
//      
//	      pdfevol::retrieve1d(i,mesq::negative); //--> both branches are negative
//	      mellinint::pdf_mesq_expy(i,i,mesq::negative);
//	      complex <double> int2 = mellinint::integrand1d(i,mesq::negative);
//	      //fun += 0.5*real(int1-int2);
//	      fun += 0.5*(int1-int2); //--> complex result for bprescription 2
//	      
//	      //fun += 0.5*(int1-int2);
//	      //cout << "C++ " << setprecision(16) << i1 << "  " << i2 << "  " << int1 << "  " << int2 << endl;
//	      
//	      //do not actually need the negative branch, because it is equal to the complex conjugate of the positive branch
//	      //so the sum of the two branches is equal to the real part of the positive branch: 1/2(a+conj(a)) = real(a)
//	      fun += int1;
//	      //cout << "mellin 1d " << setprecision(16) << i << "  " << int1 << endl;
//	      
//	    }
//	  //cout << "besselint 1d inversion " << fun << endl;
//          }
      else
	fun = mellinint::calc2d();
//	{
//	  //#pragma omp parallel for reduction(+:fun) num_threads(opts.mellincores) copyin(creno_,mesq::mesqij_expy,hcoefficients::Hqqb,hcoefficients::Hqg_1,hcoefficients::Hqg_2,hcoefficients::Hgg,hcoefficients::Hqq,hcoefficients::Hqq_1,hcoefficients::Hqq_2,hcoefficients::Hqqp_1,hcoefficients::Hqqp_2)
//	  for (int i1 = 0; i1 < mellinint::mdim; i1++)
//	    {
//	      pdfevol::retrieve_beam1(i1);
//	      for (int i2 = 0; i2 < mellinint::mdim; i2++)
//		{
//		  // here scale2 is fixed (b-dependent), and the function is called many times at I1 I2 points       
//		  // part of the coefficients calculation is hoisted in the previous i loop
//	      
//		  //  -->   merge positive and negative branch?
//	      
//		  //     In Rapidity integrated mode:
//		  //     sigmaij are fatorised from HCRN and numerical integration in y is performed in rapintegrals
//		  //     the full expression is HCRN(I1,I2)_ij * ccex(I1,I2) * sigma_ij
//		  //     HCRN(I1,I2)_ij is only b dependent
//		  //     ccex(I1,I2) is rapidity and mass dependent
//		  //     sigma_ij is costh and mass dependent, but becomes rapidity dependent after integration of the costh moments
//		  //     The integrals are solved analitically when no cuts on the leptons are applied
//	      
//		  //pdfevol::retrieve(i1,i2,mesq::positive);
//		  pdfevol::retrieve_beam2_pos(i2);
//		  mellinint::pdf_mesq_expy(i1,i2,mesq::positive);
//		  double int1 = real(mellinint::integrand2d(i1,i2,mesq::positive));
//	      
//		  //pdfevol::retrieve(i1,i2,mesq::negative);
//		  pdfevol::retrieve_beam2_neg();
//		  mellinint::pdf_mesq_expy(i1,i2,mesq::negative);
//		  double int2 = real(mellinint::integrand2d(i1,i2,mesq::negative));
//
//		  //fun += -real(0.5*(int1-int2));
//		  fun += -0.5*(int1-int2);
//		  //cout << "C++ " << setprecision(16) << i1 << "  " << i2 << "  " << int1 << "  " << int2 << endl;
//		}
//	    }
//	  //invres = fun*factorfin;
//	}
    }

  complex <double> invres;
  //if (fun == 0. || factorfin == 0.)
  //  invres = 0;
  //else
  invres = fun*factorfin;

  //complex <double> invres = fun*real(factorfin); //with bstar prescription factorfin is real
  //complex <double> invres = factorfin/resint::_m*(8./3.)*pow(opts.sroot,2)/2.; // --> Check unitarity of Sudakov integral
  //complex <double> invres = factorfin;

  
  //cout << "invres " << b << "  " << fun << "  " << factorfin << "  " << invres << endl;
  if (isnan_ofast(real(invres)) || isnan_ofast(imag(invres)))
    {
      cout << "Warning, invres = " << invres << ", qt = " << qt << ", b = "  << b << ", pdf*mesq = " << fun << ", S*bj0 = " << factorfin << endl;
      cout << "Warning, invres = " << invres << ", qt = " << qt << ", b = "  << b << ", pdf = " << pdfevol::fx1[6] << ", mesq = " << mesq::mesqij_expy[mesq::index(0,0,0,mesq::positive)] << endl;
      invres = 0;
    }
  //cout << "invres = " << invres << ", qt = " << qt << ", b = "  << b << ", pdf*mesq = " << fun << ", S " << sudakov::S << " bj0 = " << b*xj0 << endl;
  //cout << "invres = " << invres << ", qt = " << qt << ", b = "  << b << ", pdf = " << pdfevol::fx1[6] << ", mesq = " << mesq::mesqij_expy[mesq::index(0,0,0,mesq::positive)] << endl;
  //cout << setprecision(16) << "C++ " << b << "  " << invres << "  " << fun << "  " << factorfin << endl;
  //cout << setprecision(16) << " b " << b << " bstar " << bstar << " besselint " << invres << " J0 " << fort_besj0_(qtb) << " sud " << sudak << " fun " << fun << endl;
  return invres;
}
