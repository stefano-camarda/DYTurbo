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

  complex <double> scale2 = pow(resconst::b0*opts.a_param/b,2);
  //freeze PDF evolution below a certain scale
  //  if (fabs(sqrt(scale2)) < 5.)
  //    scale2 = 5.*5.;
  pdfevol::bscale = sqrt(scale2);
  fcomplex fscale2 = fcx(scale2);

  double bstar= b / sqrt(1+(b*b)/(blimit_.rblim_*blimit_.rblim_));
  complex <double> bstarscale2 = pow(resconst::b0*opts.a_param/bstar,2);
  pdfevol::bstarscale = sqrt(bstarscale2);

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
  //alpq is used in hcoefficients::calcb, it is alphas(res scale) * alphas(b0^2/b^2)
  //it is used only at NLL, at NNLL instead aexp and aexpb are used
  complex <double> alpq = resint::alpqres * cx(alphasl_(fscale2));        //alphas at the resummation scale times alphas at 1/b
  //complex <double> alpq = resint::alpqfac * cx(alphasl_(fscale2));          //alphas at the factorisation scale times alphas at 1/b
  /********************************************/

  //pdfevol::XL = pdfevol::alpqf / alpq; // = 1./cx(alphasl_(fscale2));
  pdfevol::XL = 1./cx(alphasl_(fscale2));
  pdfevol::XL1 = 1.- pdfevol::XL;
  pdfevol::SALP = log(pdfevol::XL);

  // SELECT ORDER FOR EVOLUTION LO/NLO
  pdfevol::alpr = alpqf * cx(alphasl_(fscale2))*(double)(opts.order-1);
  //force LO evolution
  //pdfevol::alpr = alpqf * cx(alphasl_(fscale2))*(double)(0);
  //cout << b << "  " << scale2 << "  " << pdfevol::SALP << "  " << log(1./cx(alphasl_(fscale2))) << "  " << pdfevol::alpr << "  " << alpq <<  endl;
  //**************************************

  //**************************************

  //original dyres evolution: Perform PDF evolution from muf to the scale corresponding to the impact parameter b
  //the scales used in the evolution correspond to SALP and alpr
  if (opts.evolmode == 0)
    for (int i = 0; i < mellinint::mdim; i++)
      pdfevol::evolution (i);

  //Calculate PDF moments by direct Mellin transformation at each value of bstarscale = b0p/bstar
  if (opts.evolmode == 2)
    for (int i = 0; i < mellinint::mdim; i++)
      pdfevol::calculate (i);

  //PDF evolution with Pegasus QCD from the starting scale Q20 
  if (opts.evolmode == 1)
    pegasus::evolve();

  //for (int i = 0; i < mellinint::mdim; i++)
  //  cout << cx(creno_.cfx1_[i][0]) << "  " << cx(creno_.cfx2p_[i][5]) << "  " <<  cx(creno_.cfx2m_[i][5]) << endl;
  //**************************************

  //aexp and aexpb are calculated in alphasl
  complex <double> aexp = cx(exponent_.aexp_);
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
