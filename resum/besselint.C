#include "resconst.h"
#include "settings.h"
#include "besselint.h"
#include "pdfevol.h"
#include "mesq.h"
#include <LHAPDF/LHAPDF.h>

//fortran interface
void besselint_(double &b, double &qt, double &q2)
{
  besselint::bint(b, qt, q2);
};

double besselint::bint(double b, double qt, double q2)
{
  //      double complex loga,logmuf2q2,logq2muf2,logq2mur2
  //      common/clogs/loga,logmuf2q2,logq2muf2,logq2mur2
      
  //      double complex aexp,aexpB
  //      COMMON/exponent/aexp,aexpB

  //      double precision aass
  //      COMMON/aass/aass

  complex <double> scale2 = pow(resconst::b0*opts.a_param/b,2);
  complex <double> bb = b;
  //     USES BESSEL FUNCTION SINCE INTEGRATION IS DONE ALONG THE REAL AXIS
  //     ********************
  //     qt and b dependence (bessel function) (jacobian?)
  double qtb = qt*b;
  double xj0= 2.*fort_besj0_(qtb);
  //     ********************
  //     Sudakov is only mass and b dependent
  fcomplex fbb = fcx(bb);
  complex <double> sudak=cx(s_(fbb));
  if (sudak == 0.)
    return 0.;
  //     ********************

  //     ********************
  //     qt and mass dependence
  complex <double> factorfin = bb*xj0*sudak;
  //     ********************

  //**************************************
  //     b-dependence
  //...  alphasl gives the LL/NLL evolution of alpha from Qres=Q/a_param to
  //     q2=bo^2/b^2
//**************************************

  double q2s = q2/pow(opts.a_param,2);                   //resummation scale
  pdfevol::alpqf = LHAPDF::alphasPDF(sqrt(q2s))/4./M_PI; //alphas at resummation scale
  fcomplex fscale2 = fcx(scale2);
  pdfevol::alpq = pdfevol::alpqf * cx(alphasl_(fscale2));    //alphas at the resummation scale times alphas at 1/b

  pdfevol::XL = pdfevol::alpqf / pdfevol::alpq; // = 1./cx(alphasl_(fscale2));
  pdfevol::XL1 = 1.- pdfevol::XL;
  pdfevol::SALP = log(pdfevol::XL);

  // SELECT ORDER FOR EVOLUTION LO/NLO
  pdfevol::alpr = pdfevol::alpq *(double)(opts.order-1);

  for (int i = 0; i < mellinint::mdim; i++)
    {
      //perform PDF evolution from muf to the scale corresponding to the impact parameter b
      pdfevol::evolution (i, mesq::positive, 1);
      pdfevol::evolution (i, mesq::positive, 2);
      pdfevol::evolution (i, mesq::negative, 2);
    }
  // Cache the positive and negative branch of coefficients which depend only on one I index
  //  double aass,logmuf2q2,loga,aexp,aexpb;
  //  hcoefficients::calcb(aass,logmuf2q2,loga,real(pdfevol::alpq),aexp,aexpb); // --> Need to access aass,logmuf2q2,loga,alpq,aexp,aexpb
  
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
	pdfevol::evolve(i1,i2,mesq::positive);
	mellinint::pdf_mesq_expy(i1,i2,mesq::positive);
	complex <double> int1 = mellinint::integrand(i1,i2,mesq::positive);
	pdfevol::evolve(i1,i2,mesq::negative);
	mellinint::pdf_mesq_expy(i1,i2,mesq::negative);
	complex <double> int2 = mellinint::integrand(i1,i2,mesq::negative);
	//complex <double> FZ=-0.5*(real(int1)-real(int2));
	fun = fun -real(0.5*(int1-int2));
	//	cout << "C++ " << setprecision(16) << i1 << "  " << i2 << "  " << int1 << "  " << int2 << endl;
      }

  double invres = fun*real(factorfin);
  if (invres != invres)
    {
      cout << "Warning, invres = " << invres << " b = "  << b << " qt = " << qt << endl;
      cout << fun << "  " << factorfin << endl;
      invres = 0;
    }
  //  cout << setprecision(16) << "C++ " << b << "  " << fun << "  " << factorfin << endl;
  return invres;
}
