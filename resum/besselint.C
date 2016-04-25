#include "resconst.h"
#include "settings.h"
#include "besselint.h"
#include "pdfevol.h"
#include "pegasus.h"
#include "mesq.h"
#include "hcoefficients.h"
#include "resint.h"

#include <LHAPDF/LHAPDF.h>

//fortran interface
void besselint_(double &b, double &qt, double &q2)
{
  resint::_qt = qt;
  resint::_m = sqrt(q2);
  besselint::bint(b);
};

double besselint::bint(double b)
{
  complex <double> bb = b;
  double qt = resint::_qt;

  //bscale = a*b0/b
  pdfevol::bscale = resconst::b0*resint::a/b;
  //freeze PDF evolution below a certain scale
  //  if (pdfevol::bscale < 5.)
  //    pdfevol::bscale = 5.;
  fcomplex fscale2 = fcx(pow(pdfevol::bscale,2));
  
  //bstarscale = a*b0/bstar (final scale used for the PDF evolution)
  double bstar = b / sqrt(1+(b*b)/(blimit_.rblim_*blimit_.rblim_));
  pdfevol::bstarscale = resconst::b0*resint::a/bstar;

  //bstarscale with the modification L -> L~, which freezes the scale at muf
  pdfevol::bstartilde = pdfevol::bstarscale * resint::mufac / sqrt((pow(pdfevol::bstarscale,2) + resint::mufac2));
  
  //qbstar = b0/bstar (without a_param)
  pdfevol::qbstar = resconst::b0/bstar;
  
  //  cout << b << "  " << bstar << "  " << blimit_.rblim_ << endl;
  
  //The integration from b to qt space is done with the bstar prescription (real axis in the b space), and use the bessel function
  //The (complex) integration in the minimal prescription would require hankel functions
  //********************
  //qt and b dependence (bessel function) (is xj0 a jacobian?)
  double qtb = qt*b;
  double xj0= 2.*fort_besj0_(qtb);
  //********************

  //The Sudakov is mass and b dependent
  fcomplex fbb = fcx(bb);
  complex <double> sudak=cx(s_(fbb));
  if (sudak == 0.)
    return 0.;

  //********************
  //b, qt and mass dependence
  complex <double> factorfin = bb*xj0*sudak;
  //********************

  //**************************************
  //b-dependence
  // Set scales for evolution in pdfevol
  //alphasl gives the LL/NLL evolution of alpha from Qres=Q/a_param to  q2=b0^2/b^2

  /********************************************/
  //pdfevol::alpqf is used as starting scale in pdfevol::evolution
  //according to Eq. 42 of arXiv:hep-ph/0508068 it should be the factorisation scale,
  //but in Eq. 98 in the resummation scale is used
  double alpqf = resint::alpqres;              //alphas at resummation scale   (as in dyres)
  //double alpqf = resint::alpqfac;            //alphas at factorisation scale (as in Eq. 42 of arXiv:hep-ph/0508068, but see also Eq. 98)
  /********************************************/
  
  /********************************************/
  //alpq is used in hcoefficients::calcb, it is alpq = alphas(b0^2/b^2)
  //it is used only at NLL, at NNLL instead aexp and aexpb are used (aexp is the same as alphasl, but with a different blim)
  complex <double> alpq;
  alpq = resint::alpqres * cx(alphasl_(fscale2));
  if (opts.evolmode == 2 || opts.evolmode == 3)
    //In order to have variable flavour evolution, use here a VFN definition of alpq
    //There is possibly an issue here when the renormalisation scale is (very) different from mll, since aass=alphas(muren) is used in xlambda = beta0*aass*blog
    {
      //      double M2 = pow(fabs(pdfevol::bstartilde),2);

      double M2 = pow(fabs(pdfevol::qbstar),2);   //qbstar = b0/bstar (without a_param)
      double M2tilde = M2 * resint::mures2 / (M2 + resint::mures2);
      double M2prime = M2tilde * resint::muren2/resint::mures2; //scale for differences between muren and mures
      M2 = M2prime;

      if (M2 > asfthr_.m2t_) //The scale is frozen at muf, and never goes above the top
	{
	  int NF = 6;
	  double R2T = asfthr_.m2t_;
	  alpq = as_(M2, R2T, asfthr_.ast_, NF);
	}
      else if (M2 > asfthr_.m2b_)
	{
	  int NF = 5;
	  double R2B = asfthr_.m2b_;
	  alpq = as_(M2, R2B, asfthr_.asb_, NF);

	  /*
	  //tweak this to look as in dyres, where alphas(muren) is used:
	  double R20 = resint::_m;
	  alpq = as_(M2, R20, resint::alpqren, NF);
	  alpq = alpq/resint::alpqren*resint::alpqres;
	  */
	}
      else if (M2 > asfthr_.m2c_)
	{
	  int NF = 4;
	  double R2C = asfthr_.m2c_;
	  alpq = as_(M2, R2C, asfthr_.asc_, NF);
	}
      else
	{
	  int NF = 3;
	  double R20 = asinp_.m20_;
	  alpq = as_(M2, R20, asinp_.as0_, NF);
	}
      //Based on the definition of aexp, may be need to rescale alpq for as(qres)/as(muren)?
      //alpq = alpq/resint::alpqren*resint::alpqres;
    }

  //  complex <double> alpq = LHAPDF::alphasPDF(fabs(pdfevol::bstartilde))/4./M_PI;        //alpq = alphas(b0^2/b^2)
  //  cout << pdfevol::bstarscale << "  " << resint::alpqres * cx(alphasl_(fscale2)) << "  " <<  LHAPDF::alphasPDF(fabs(pdfevol::bstartilde))/4./M_PI << endl;
  /********************************************/

  //pdfevol::XL = pdfevol::alpqf / alpq; // = 1./cx(alphasl_(fscale2));
  pdfevol::XL = 1./cx(alphasl_(fscale2)); //XL = alphas(mures2)/alphas(b0^2/b^2)
  pdfevol::XL1 = 1.- pdfevol::XL;
  pdfevol::SALP = log(pdfevol::XL);
  //cout << pdfevol::bscale << "  " << pdfevol::XL << endl;
  // SELECT ORDER FOR EVOLUTION LO/NLO
  pdfevol::alpr = alpqf * cx(alphasl_(fscale2))*(double)(opts.order-1);
  //force LO evolution
  //pdfevol::alpr = alpqf * cx(alphasl_(fscale2))*(double)(0);
  //cout << b << "  " << scale2 << "  " << pdfevol::SALP << "  " << log(1./cx(alphasl_(fscale2))) << "  " << pdfevol::alpr << "  " << alpq <<  endl;
  //**************************************

  //**************************************
  //Perform PDF evolution
  
  //original dyres evolution: Perform PDF evolution from muf to the scale b0/b ~ pt
  //the scales used in the evolution correspond to SALP and alpr
  if (opts.evolmode == 0)
    for (int i = 0; i < mellinint::mdim; i++)
      pdfevol::evolution (i);

  //PDF evolution with Pegasus QCD from muf to the scale b0/b ~ pt
  else if (opts.evolmode == 1)
    pegasus::evolve();

  //Calculate PDF moments by direct Mellin transformation at each value of bstarscale = b0p/bstar
  else if (opts.evolmode == 2)
    for (int i = 0; i < mellinint::mdim; i++)
      pdfevol::calculate (i);

  //PDF evolution with Pegasus QCD from the starting scale Q20, in VFN
  else if (opts.evolmode == 3)
    pegasus::evolve();

  //PDF evolution with Pegasus QCD from the starting scale Q20, in VFN, and using alphasl(nq2) to evaluate ASF at the final scale
  else if (opts.evolmode == 4)
    pegasus::evolve();

  //for (int i = 0; i < mellinint::mdim; i++)
  //  cout << cx(creno_.cfx1_[i][0]) << "  " << cx(creno_.cfx2p_[i][5]) << "  " <<  cx(creno_.cfx2m_[i][5]) << endl;
  //**************************************

  //aexp and aexpb are calculated in alphasl, they are used in hcoefficients::calcb only for the NNLL cross section
  complex <double> aexp = cx(exponent_.aexp_); //aexp is actually the ratio alphas(a*b0/b)/alphas(muren)
  
  //In order to have variable flavour evolution, use here a VFN definition of aexp (aexpB instead, should be independent from NF)
  //It seems that the pt-integrated cross section is invariant for variations of a_param if aexp is the ratio of alpq/as(muren) (and not alpq/as(Qres))
  //Indeed, aexp is always multiplied by aass, which is alphas at the renormalisation scale
  if (opts.evolmode == 2 || opts.evolmode == 3) //account for VFN evolution
    aexp = alpq / resint::alpqren;
    //aexp = alpq / resint::alpqres;
  
  //  cout << pdfevol::bstartilde << "  " << cx(exponent_.aexp_) << "  " << alpq / resint::alpqres << "  " << alpq / resint::alpqren << endl;
  complex <double> aexpb = cx(exponent_.aexpb_);
  
  // Cache the positive and negative branch of coefficients which depend only on one I index
  hcoefficients::calcb(resint::aass,resint::logmuf2q2,resint::loga,alpq,aexp,aexpb); // --> Need to access aass,logmuf2q2,loga,alpq,aexp,aexpb
  
  double fun = 0.;
  for (int i1 = 0; i1 < mellinint::mdim; i1++)
    for (int i2 = 0; i2 < mellinint::mdim; i2++)
      {
	// here scale2 is fixed (b-dependent), and the function is called many times at I1 I2 points       
	// part of the coefficients calculation is hoisted in the previous i loop

	//  -->   merge positive and negative branch?

	//     In Rapidity integrated mode:
	//     sigmaij are fatorised from HCRN and numerical integration in y is performed in rapintegrals
	//     the full expression is HCRN(I1,I2)_ij * ccex(I1,I2) * sigma_ij
	//     HCRN(I1,I2)_ij is only b dependent
	//     ccex(I1,I2) is rapidity and mass dependent
	//     sigma_ij is costh and mass dependent, but becomes rapidity dependent after integration of the costh moments
	//     The integrals are solved analitically when no cuts on the leptons are applied
	pdfevol::retrieve(i1,i2,mesq::positive);
	mellinint::pdf_mesq_expy(i1,i2,mesq::positive);
	complex <double> int1 = mellinint::integrand(i1,i2,mesq::positive);
	pdfevol::retrieve(i1,i2,mesq::negative);
	mellinint::pdf_mesq_expy(i1,i2,mesq::negative);
	complex <double> int2 = mellinint::integrand(i1,i2,mesq::negative);
	//complex <double> FZ=-0.5*(real(int1)-real(int2));
	fun = fun -real(0.5*(int1-int2));
	//	cout << "C++ " << setprecision(16) << i1 << "  " << i2 << "  " << int1 << "  " << int2 << endl;
      }

  double invres = fun*real(factorfin);
  if (invres != invres)
    //  if (isnan(invres))
    {
      cout << "Warning, invres = " << invres << ", qt = " << qt << ", b = "  << b << ", pdf*mesq = " << fun << ", S*bj0 = " << factorfin << endl;
      invres = 0;
    }
  //  cout << setprecision(16) << "C++ " << b << "  " << invres << "  " << fun << "  " << factorfin << endl;
  return invres;
}
