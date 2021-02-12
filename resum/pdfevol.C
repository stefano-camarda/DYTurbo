#include "pdfevol.h"
#include "interface.h"
#include "settings.h"
#include "mesq.h"
#include "anomalous.h"
#include "pegasus.h"
#include "resconst.h"
#include "resint.h"
#include "chebyshev.h"
#include "phasespace.h"
#include "mellinpdf.h"
//#include "cerespdf.h"
#include "scales.h"
#include "parton.h"
#include "npff.h"
#include "blim.h"
#include "pdf.h"
#include "evolnative.h"
#include "evolnum.h"
#include "gint.h"
#include "clock_real.h"

#include <LHAPDF/LHAPDF.h>

#include <iostream>
#include <complex>

//Main interface for PDF evolution, which should:
// - manage ancillary evolution modules: pegasus, native-evolution (mela in the future)
// - Call mellinpdf for the transform from x to N space
// - provide evolved PDFs to mellinint
// - pdfevol should store only evolved PDFs, not PDFs at the starting scale (whether muf or Q0)

//init
//update
//evolve

complex <double> *pdfevol::fx1;
complex <double> *pdfevol::fx2;

complex <double> pdfevol::fn1[2*MAXNF+1];
complex <double> pdfevol::fn2[2*MAXNF+1];

//scales (obsolete)
complex <double> pdfevol::bscale;
complex <double> pdfevol::bstarscale;
complex <double> pdfevol::bstartilde;
complex <double> pdfevol::qbstar;
complex <double> pdfevol::bcomplex;

//scales
complex <double> pdfevol::bstar;
complex <double> pdfevol::mub_a;
complex <double> pdfevol::mub;
complex <double> pdfevol::mubstar_a;
complex <double> pdfevol::mubstar;
complex <double> pdfevol::mubstartilde;

complex <double> pdfevol::logasl;
complex <double> pdfevol::asl;

using namespace parton;
using namespace resconst;
using namespace resint;

//fortran interface
void pdfevol_(int& i1, int& i2, int& sign)
{
  pdfevol::retrieve(i1-1, i2-1, sign-1);
};

//evolved moments
void pdfevol::allocate_fx()
{
  fx1 = new complex <double>[mellinint::mdim*(2*MAXNF+1)*2];
  fx2 = new complex <double>[mellinint::mdim*(2*MAXNF+1)*2];
}
void pdfevol::free_fx()
{
  delete[] fx1;
  delete[] fx2;
}

void pdfevol::init()
{
  //set up the x-spacing for the x-to-N transform
  mellinpdf::init();

  //Always init all the modules
  evolnative::init(); //--> Needed for the counterterm by ctmellin
  pegasus::init();    //--> Needed for the gamma1,2,3 (moments of the splitting functions)
  //evolnum::init();
  
  //switch (opts.evolmode)
  //  {
  //  case 0: evolnative::init(); break;
  //  case 1: pegasus::init();    break;
  //  case 2: evolnum::init();    break;
  //  case 3: pegasus::init();    break;
  //  case 4: pegasus::init();    break;
  //  }
}

void pdfevol::release()
{
  pegasus::release();
  evolnative::release();
  mellinpdf::release();
 
  //switch (opts.evolmode)
  //  {
  //  case 0: evolnative::release(); break;
  //  case 1: pegasus::release();    break;
  //  case 2: evolnum::release();     break;
  //  case 3: pegasus::release();    break;
  //  case 4: pegasus::release();    break;
  //  }
}

//compute the LL,NLL,NNLL,NNNLL evolution of alphas
void pdfevol::alphasl(complex <double> b)
{
  double as = resint::aass; // = pdf::alphas(muren)/4./M_PI;
  double Q = scales::res;
  double LQR = LR-LQ;//real(logq2mur2-2.*loga); //2*loga = logq2mures2;
  /*
  if (opts.mufevol)  //should be actually if (opts.mufevol || opts.evolmode == 2 || opts.evolmode == 3)
    {
      //      cout << endl;
      //      cout << as << "  " << Q << "  " << LQR << endl;
      Q = scales::fac;
      LQR = resint::LR - resint::LF;

      //Q = phasespace::m;
      //as = pdf::alphas(Q)/M_PI;
      //LQR = resint::LR;

      //Q = phasespace::m;
      //as = pdf::alphas(Q)/M_PI;
      //LQR = 0.;
      
      //Q = scales::fac;
      //as = pdf::alphas(Q)/M_PI;
      //LQR = 0.;
      //      cout << as << "  " << Q << "  " << LQR << endl;
    }
  */
  //artificially remove mur scale variations
  //Q = scales::fac;
  //as = pdf::alphas(Q)/M_PI;
  //LQR = 0.;

  double as2 = pow(as,2);

  double blim = blim::pdf;

  //Set b according to bstar or other prescriptions
  complex <double> bstar;
  if (opts.bprescription == 0 || opts.bprescription == 4 || opts.bstar_pdf)
    bstar = real(b)/sqrt(1.+pow(real(b)/blim,2));
  else
    bstar = b;

  complex <double> blog;
  if (!opts.modlog)
    blog = log(pow(Q*bstar/b0,2));   //normal sudakov
  else if (opts.p == 1)
    blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  else
    blog = 1./opts.p*log(pow(Q*bstar/b0,2*opts.p) + 1.); //modified sudakov with exponent p
  

  //  if (opts.modlog)
  //    blog = log(pow(Q*bstar/b0,2) + 1.); //modified sudakov
  //  else
  //    blog = log(pow(Q*bstar/b0,2));   //normal sudakov

  //lambda as defined in Eq. (25) of hep-ph/0508068.
  complex <double> xlambda = beta0*as*blog;
  complex <double> log1xlambda = log(1.-xlambda);
  //cout << " b " << b << " bstar " << bstar << " as " << as << " beta0 " << beta0 << " blog " << blog << " xlambda " << xlambda << endl;
  //  aexp = -Log[1 - lam]
  //    - aS/Pi*(beta1*Log[1 - lam] + beta0^2*lam*Log[Q2/muR2])/(beta0*(1 - lam))
  //    + aS^2/Pi^2*(-2*beta1^2*lam + 2*beta0*beta2*lam - 2*beta1^2*Log[1 - lam] + beta1^2*Log[1 - lam]^2 - 2*beta0^2*beta1*lam*Log[Q2/muR2] + 2*beta0^2*beta1*Log[1 - lam]*Log[Q2/muR2])/(2*beta0^2*(-1 + lam)^2)

  //LL (no evolution)
  logasl = 0.;
  
  //NLL (LO evolution)
  if (opts.order >= 1)
    logasl += log1xlambda;

  //NNLL (NLO evolution)
  if (opts.order >= 2)
    logasl += as* beta1/beta0*log1xlambda/(1.-xlambda);

  //NNNLL (NNLO evolution)
  if (opts.order >= 3)
    logasl += as2* ((pow(beta1/beta0,2)-beta2/beta0) *xlambda/pow(1.-xlambda,2)
		      + pow(beta1/beta0,2)             *log1xlambda/pow(1.-xlambda,2)
		      - pow(beta1/beta0,2)             *pow(log1xlambda,2)/(2.*pow(1.-xlambda,2)));

  //QCD coupling scale dependence
  if (opts.order >= 2)
    logasl += as*(LQR)
      *beta0*xlambda/(1.-xlambda);

  if (opts.order >= 3)
    logasl += as2*(+LQR*beta1                   *(xlambda-log1xlambda)/pow(1.-xlambda,2)
		   +LQR*beta1                   *xlambda/(1.-xlambda)                      //missing piece
		   +0.5*pow(LQR,2)*pow(beta0,2) *xlambda*(xlambda-2.)/pow(1.-xlambda,2));  //missing piece

  //LF scale dependence --> This is not formally correct, but improves muF variation bands, at least at NNLL
  //double LFQ = LQ-LF;
  //if (opts.order >=2)
  //  logasl += as*(LFQ)
  //    *beta0*xlambda/(1.-xlambda);
  //
  //if (opts.order >= 3)
  //  logasl += as2*(+LFQ*beta1                   *(xlambda-log1xlambda)/pow(1.-xlambda,2)
  //		   +LFQ*beta1                   *xlambda/(1.-xlambda)                      //missing piece
  //		   +0.5*pow(LFQ,2)*pow(beta0,2) *xlambda*(xlambda-2.)/pow(1.-xlambda,2));  //missing piece

  //cout << endl;
  //cout << "analytic logasl " << logasl        << endl;
  //cout << "numeric  logasl " << gint::logasl_pdf  << endl;

  if (opts.numsud || opts.order >= 4)
    logasl = gint::logasl_pdf;
  
  asl = exp(-logasl);
}

void pdfevol::scales(complex <double> b)
{
  //bb = b;
  
  mub_a = resint::a * resconst::b0/b;

  //complex PDF scale to be used for the minimal prescription
  mub = resconst::b0/b;
  
  double blim = blim::pdf; //1.1229190 -->hard-coded!! --> should allow a different blim in the PDF evolution as a setting
  bstar = real(b)/sqrt(1.+(pow(real(b)/blim,2)));

  //mubstar_a = a*b0/bstar is the final scale used for the PDF evolution (previously called bstarscale)
  mubstar_a  = resconst::b0*resint::a/bstar; //bstarscale = resconst::b0*resint::a/bstar;
  
  //mubstar is used in evolmode 3 (previously called qbstar), it corresponds to mubstar_a but without a_param
  mubstar = resconst::b0/bstar;

  //simulate pythia ISR factorisation scale, which is muf^2 = qt^2 (i.e. avoid the bstar prescription to freeze the factorisation scale)
  //if (opts.bprescription_pdf = ...) //--> minimal prescription for pdfs
  //mubstar = resconst::b0/b;
  
  //mubstartilde is mubstar (qbstar) with the modification L -> L~, which freezes the scale at mures (or mufac)
  //it is used in evolmode 2, for the direct mellin transfrom at each scale

  double Q = resint::mures;
  //if (opts.mufevol || opts.evolmode == 2)
  //Q = resint::mufac; //--> instead of usingg mufac in mubstartilde, apply an evolution operator U(Q,muF) to the PDFs also in evolmode = 2

  bool useC2 = false;
  if (useC2) //convert muF variations to C2 variations
    Q = resint::mures;

  if (!opts.modlog)
    mubstartilde = mubstar;
  else if (opts.p == 1)
    mubstartilde = mubstar * Q / sqrt((pow(mubstar,2) + pow(Q,2)));
  //mubstartilde = mubstar / sqrt(1.+pow(mubstar/Q,2));
  else
    mubstartilde = mubstar * Q / pow((pow(mubstar,2*opts.p) + pow(Q,2*opts.p)),0.5/opts.p);

  //if (opts.modlog)
  //  mubstartilde = mubstar * Q / sqrt((pow(mubstar,2) + pow(Q,2)));
  ////mubstartilde = mubstar / sqrt(1.+pow(mubstar/Q,2));
  //else
  //  mubstartilde = mubstar;

  if (useC2) //convert muF to a C2 variation
    mubstartilde *= resint::mufac/resint::mures;
    
  /***********************  C1 and C3  ******************
  //Resbos scales C1 and C3 see arxiv:1309.1393 for details (not really correctly implemented)
  double C1 = 1;
  double C3 = 1;

  //bscale is used for
  //alpq (C1)
  //alpr (C3)
  //XL (C3)
  //aexp in cexp (C1)
  //it also used in evolmode 1 and 4 as factorisation scale
  

  //The soft scale of the resummation is b0/b. The factor a = mll/mures is introduced because the b scale is used in alphasl
  //Introduce two scales:
  //mubC1 = C1*b0/b is the lower integration limit of the Sudakov integral
  mubC1 = mub * C1;

  //mubC3 = C3*b0/b is the "soft" factorisation scale at which the Wilson coefficient functions are evaluated (see arxiv:1309.1393 for details)
  mubC3 = mub * C3;
  
  //convert to fortran complex number
  fcomplex fscale2_mubC1 = fcx(pow(mubC1,2));
  fcomplex fscale2_mubC3 = fcx(pow(mubC3,2));

  //mubbstar and mubstartilde should have the C3 scale
  mubstarC3 = mubstar*C3;
  mubstartildeC3 = mubstartilde;
  ***********************  end C1 and C3  ******************/
}

//Main evolution selector
void pdfevol::evolution()
{
  switch (opts.evolmode)
    {
    //original dyres evolution: Perform PDF evolution from Q to b0/b
    case 0: evolnative::evolve(); break;
      
    //PDF evolution with Pegasus QCD from Q to b0/b
    case 1: pegasus::evolve();    break;
      
    //Calculate PDF moments by direct Mellin transformation at each value of mubstartilde ~ b0/bstar
    case 2: evolnum::calculate();          break;
      
    //PDF evolution with Pegasus QCD from the starting scale Q20, in VFN
    case 3: pegasus::evolve();    break;
      
    //PDF evolution with Pegasus QCD from the starting scale Q20, in VFN, and using alphasl(nq2) to evaluate ASF at the final scale
    case 4: pegasus::evolve();    break;
    }
}

void pdfevol::allocate()
{
  if (opts.evolmode == 0) //--> always allocate evolnative for the counterterm in Mellin space
    evolnative::allocate();

  if (opts.evolmode == 1 || opts.evolmode == 3)
    pegasus::allocate();
}

void pdfevol::free()
{
  if (opts.evolmode == 0) //--> always free evolnative for the counterterm in Mellin space
    evolnative::free();

  if (opts.evolmode == 1 || opts.evolmode == 3)
    pegasus::free();
}


//Update the Mellin transform from x to N at the factorisation scale --> probably no need for this function, evolutions modules should deal with this inside their evolution routines
//This function is needed by ctmellin!!!
void pdfevol::update()
{
  if (opts.evolmode == 0)
    evolnative::update();

  if (opts.evolmode == 1)
    pegasus::update();

  if (opts.evolmode == 3)
    pegasus::update();
  

  /*
//Deal with this inside pegasus
  if (opts.evolmode == 1)
    .....
  

  //Restrict the integration of moments to xmin = m/sqrt(s)*exp(-ymax) = (m/sqrt(s))^2
  //double xmin = pow(phasepace::m/opts.sroot,2);

  //No need to init mellinpdf as it was already initialised
  //mellinpdf::init();
  
  //scales::set(phasespace::m); //assume scales were already set

  double muf = -1;
  if (opts.evolmode == 3)
    muf = pdf::qmin;
  
  mellinpdf::allocate();
  mellinpdf::update(muf);
  
  for (int k = 0; k < mellinint::mdim; k++)
    {
      UVP[k] = mellinpdf::UV[k];
      DVP[k] = mellinpdf::DV[k];
      USP[k] = mellinpdf::US[k];
      DSP[k] = mellinpdf::DS[k];
      SSP[k] = mellinpdf::SM[k];
      GLP[k] = mellinpdf::GL[k];
      CHP[k] = mellinpdf::CP[k];
      BOP[k] = mellinpdf::BP[k];
      SVP[k] = mellinpdf::SP[k] - mellinpdf::SM[k];
      CVP[k] = mellinpdf::CP[k] - mellinpdf::CM[k];
      BVP[k] = mellinpdf::BP[k] - mellinpdf::BM[k];
      //cout << "moment " << k << "  " << mellinint::Np[k] << "  ";
      //cout << "uval  " << mellinpdf::UV[k] << endl;
      //cout << "dval  " << mellinpdf::DV[k] << endl;
      //cout << "usea  " << mellinpdf::US[k] << endl;
      //cout << "dsea  " << mellinpdf::DS[k] << endl;
      //cout << "ssea  " << mellinpdf::SM[k] << endl;
      //cout << "gluon " << mellinpdf::GL[k] << endl;
      //cout << "charm " << mellinpdf::CP[k] << endl;
      //cout << "bottom" << mellinpdf::BP[k] << endl;
    }

  mellinpdf::free();
  */
}

void pdfevol::storemoments(int i, complex <double> fx[11])
{
  //Save the evolved PDFs into the fx1 and fx2 arrays
  int negidx = mellinint::mdim*11;
  int nf = 2*MAXNF+1;

  copy(fx, fx+nf, fx1+i*nf);
  copy(fx, fx+nf, fx2+i*nf);

  if (opts.ih1 == -1)
    {
      fx1[i*nf+bb] = fx[b ];
      fx1[i*nf+cb] = fx[c ];
      fx1[i*nf+sb] = fx[s ];
      fx1[i*nf+db] = fx[d ];
      fx1[i*nf+ub] = fx[u ];
      fx1[i*nf+u ] = fx[ub];
      fx1[i*nf+d ] = fx[db];
      fx1[i*nf+s ] = fx[sb];
      fx1[i*nf+c ] = fx[cb];
      fx1[i*nf+b ] = fx[bb];
    }
  if (opts.ih2 == -1)
    {
      fx2[i*nf+bb] = fx[b ];
      fx2[i*nf+cb] = fx[c ];
      fx2[i*nf+sb] = fx[s ];
      fx2[i*nf+db] = fx[d ];
      fx2[i*nf+ub] = fx[u ];
      fx2[i*nf+u ] = fx[ub];
      fx2[i*nf+d ] = fx[db];
      fx2[i*nf+s ] = fx[sb];
      fx2[i*nf+c ] = fx[cb];
      fx2[i*nf+b ] = fx[bb];
    }

//  //beam 1 positive
//  fx1[i*nf+bb] = fx[bb];
//  fx1[i*nf+cb] = fx[cb];
//  fx1[i*nf+sb] = fx[sb];
//  fx1[i*nf+g ] = fx[g ];
//  fx1[i*nf+s ] = fx[s ];
//  fx1[i*nf+c ] = fx[c ];
//  fx1[i*nf+b ] = fx[b ];
//  if (opts.ih1 == 1)
//    {
//      fx1[i*nf+db] = fx[db];
//      fx1[i*nf+ub] = fx[ub];
//      fx1[i*nf+u ] = fx[u ];
//      fx1[i*nf+d ] = fx[d ];
//    }
//  else if (opts.ih1 == -1)
//    {
//      fx1[i*nf+db] = fx[d ];
//      fx1[i*nf+ub] = fx[u ];
//      fx1[i*nf+u ] = fx[ub];
//      fx1[i*nf+d ] = fx[db];
//    }
//
//  //beam 1 negative (never used)
//  fx1[negidx+i*nf+bb] = conj(fx[bb]);
//  fx1[negidx+i*nf+cb] = conj(fx[cb]);
//  fx1[negidx+i*nf+sb] = conj(fx[sb]);
//  fx1[negidx+i*nf+g ] = conj(fx[g ]);
//  fx1[negidx+i*nf+s ] = conj(fx[s ]);
//  fx1[negidx+i*nf+c ] = conj(fx[c ]);
//  fx1[negidx+i*nf+b ] = conj(fx[b ]);
//  if (opts.ih1 == 1)
//    {
//      fx1[negidx+i*nf+db] = conj(fx[db]);
//      fx1[negidx+i*nf+ub] = conj(fx[ub]);
//      fx1[negidx+i*nf+u ] = conj(fx[u ]);
//      fx1[negidx+i*nf+d ] = conj(fx[d ]);
//    }
//  else if (opts.ih1 == -1)
//    {
//      fx1[negidx+i*nf+db] = conj(fx[d ]);
//      fx1[negidx+i*nf+ub] = conj(fx[u ]);
//      fx1[negidx+i*nf+u ] = conj(fx[ub]);
//      fx1[negidx+i*nf+d ] = conj(fx[db]);
//    }
//
//  //beam 2 positive
//  fx2[i*nf+bb] = fx[bb];
//  fx2[i*nf+cb] = fx[cb];
//  fx2[i*nf+sb] = fx[sb];
//  fx2[i*nf+g ] = fx[g ];
//  fx2[i*nf+s ] = fx[s ];
//  fx2[i*nf+c ] = fx[c ];
//  fx2[i*nf+b ] = fx[b ];
//  if (opts.ih2 == 1)
//    {
//      fx2[i*nf+db] = fx[db];
//      fx2[i*nf+ub] = fx[ub];
//      fx2[i*nf+u ] = fx[u ];
//      fx2[i*nf+d ] = fx[d ];
//    }
//  else if (opts.ih2 == -1)
//    {
//      fx2[i*nf+db] = fx[d ];
//      fx2[i*nf+ub] = fx[u ];
//      fx2[i*nf+u ] = fx[ub];
//      fx2[i*nf+d ] = fx[db];
//    }

//  if (opts.mellin1d && opts.bprescription != 2) return;
  
  //beam 2 negative
  fx2[negidx+i*nf+g ] = conj(fx[g ]);
  if (opts.ih2 == 1)
    {
      fx2[negidx+i*nf+bb] = conj(fx[bb]);
      fx2[negidx+i*nf+cb] = conj(fx[cb]);
      fx2[negidx+i*nf+sb] = conj(fx[sb]);
      fx2[negidx+i*nf+db] = conj(fx[db]);
      fx2[negidx+i*nf+ub] = conj(fx[ub]);
      fx2[negidx+i*nf+u ] = conj(fx[u ]);
      fx2[negidx+i*nf+d ] = conj(fx[d ]);
      fx2[negidx+i*nf+s ] = conj(fx[s ]);
      fx2[negidx+i*nf+c ] = conj(fx[c ]);
      fx2[negidx+i*nf+b ] = conj(fx[b ]);
    }
  else if (opts.ih2 == -1)
    {
      fx2[negidx+i*nf+bb] = conj(fx[b ]);
      fx2[negidx+i*nf+cb] = conj(fx[c ]);
      fx2[negidx+i*nf+sb] = conj(fx[s ]);
      fx2[negidx+i*nf+db] = conj(fx[d ]);
      fx2[negidx+i*nf+ub] = conj(fx[u ]);
      fx2[negidx+i*nf+u ] = conj(fx[ub]);
      fx2[negidx+i*nf+d ] = conj(fx[db]);
      fx2[negidx+i*nf+s ] = conj(fx[sb]);
      fx2[negidx+i*nf+c ] = conj(fx[cb]);
      fx2[negidx+i*nf+b ] = conj(fx[bb]);
    }
}  

void pdfevol::storemoments_1(int i, complex <double> fx[11])
{
  //Save the evolved PDFs into the fx1 and fx2 arrays
  int negidx = mellinint::mdim*11;
  int nf = 2*MAXNF+1;

  copy(fx, fx+nf, fx1+i*nf);

  if (opts.ih1 == -1)
    {
      fx1[i*nf+bb] = fx[b ];
      fx1[i*nf+cb] = fx[c ];
      fx1[i*nf+sb] = fx[s ];
      fx1[i*nf+db] = fx[d ];
      fx1[i*nf+ub] = fx[u ];
      fx1[i*nf+u ] = fx[ub];
      fx1[i*nf+d ] = fx[db];
      fx1[i*nf+s ] = fx[sb];
      fx1[i*nf+c ] = fx[cb];
      fx1[i*nf+b ] = fx[bb];
    }

//  //beam 1 negative
//  fx1[negidx+i*nf+bb] = conj(fx[bb]);
//  fx1[negidx+i*nf+cb] = conj(fx[cb]);
//  fx1[negidx+i*nf+sb] = conj(fx[sb]);
//  fx1[negidx+i*nf+g ] = conj(fx[g ]);
//  fx1[negidx+i*nf+s ] = conj(fx[s ]);
//  fx1[negidx+i*nf+c ] = conj(fx[c ]);
//  fx1[negidx+i*nf+b ] = conj(fx[b ]);
//  if (opts.ih1 == 1)
//    {
//      fx1[negidx+i*nf+db] = conj(fx[db]);
//      fx1[negidx+i*nf+ub] = conj(fx[ub]);
//      fx1[negidx+i*nf+u ] = conj(fx[u ]);
//      fx1[negidx+i*nf+d ] = conj(fx[d ]);
//    }
//  else if (opts.ih1 == -1)
//    {
//      fx1[negidx+i*nf+db] = conj(fx[d ]);
//      fx1[negidx+i*nf+ub] = conj(fx[u ]);
//      fx1[negidx+i*nf+u ] = conj(fx[ub]);
//      fx1[negidx+i*nf+d ] = conj(fx[db]);
//    }
//
}  

void pdfevol::storemoments_2(int i, complex <double> fx[11])
{
  //Save the evolved PDFs into the fx1 and fx2 arrays
  int negidx = mellinint::mdim*11;
  int nf = 2*MAXNF+1;

  copy(fx, fx+nf, fx2+i*nf);

  if (opts.ih2 == -1)
    {
      fx2[i*nf+bb] = fx[b ];
      fx2[i*nf+cb] = fx[c ];
      fx2[i*nf+sb] = fx[s ];
      fx2[i*nf+db] = fx[d ];
      fx2[i*nf+ub] = fx[u ];
      fx2[i*nf+u ] = fx[ub];
      fx2[i*nf+d ] = fx[db];
      fx2[i*nf+s ] = fx[sb];
      fx2[i*nf+c ] = fx[cb];
      fx2[i*nf+b ] = fx[bb];
    }

  //beam 2 negative
  fx2[negidx+i*nf+g ] = conj(fx[g ]);
  if (opts.ih2 == 1)
    {
      fx2[negidx+i*nf+bb] = conj(fx[bb]);
      fx2[negidx+i*nf+cb] = conj(fx[cb]);
      fx2[negidx+i*nf+sb] = conj(fx[sb]);
      fx2[negidx+i*nf+db] = conj(fx[db]);
      fx2[negidx+i*nf+ub] = conj(fx[ub]);
      fx2[negidx+i*nf+u ] = conj(fx[u ]);
      fx2[negidx+i*nf+d ] = conj(fx[d ]);
      fx2[negidx+i*nf+s ] = conj(fx[s ]);
      fx2[negidx+i*nf+c ] = conj(fx[c ]);
      fx2[negidx+i*nf+b ] = conj(fx[b ]);
    }
  else if (opts.ih2 == -1)
    {
      fx2[negidx+i*nf+bb] = conj(fx[b ]);
      fx2[negidx+i*nf+cb] = conj(fx[c ]);
      fx2[negidx+i*nf+sb] = conj(fx[s ]);
      fx2[negidx+i*nf+db] = conj(fx[d ]);
      fx2[negidx+i*nf+ub] = conj(fx[u ]);
      fx2[negidx+i*nf+u ] = conj(fx[ub]);
      fx2[negidx+i*nf+d ] = conj(fx[db]);
      fx2[negidx+i*nf+s ] = conj(fx[sb]);
      fx2[negidx+i*nf+c ] = conj(fx[cb]);
      fx2[negidx+i*nf+b ] = conj(fx[bb]);
    }
}  

void pdfevol::storemoments_fortran(int i, complex <double> fx[11])
{
  //Save the evolved PDFs into the fortran common block
  //beam 1
  creno_.cfx1_[i][-5+MAXNF] = fcx(fx[-5+MAXNF]);
  creno_.cfx1_[i][-4+MAXNF] = fcx(fx[-4+MAXNF]);
  creno_.cfx1_[i][-3+MAXNF] = fcx(fx[-3+MAXNF]);
  creno_.cfx1_[i][ 0+MAXNF] = fcx(fx[ 0+MAXNF]);
  creno_.cfx1_[i][ 3+MAXNF] = fcx(fx[ 3+MAXNF]);
  creno_.cfx1_[i][ 4+MAXNF] = fcx(fx[ 4+MAXNF]);
  creno_.cfx1_[i][ 5+MAXNF] = fcx(fx[ 5+MAXNF]);
  if (opts.ih1 == 1)
    {
      creno_.cfx1_[i][-2+MAXNF] = fcx(fx[-2+MAXNF]);
      creno_.cfx1_[i][-1+MAXNF] = fcx(fx[-1+MAXNF]);
      creno_.cfx1_[i][ 1+MAXNF] = fcx(fx[ 1+MAXNF]);
      creno_.cfx1_[i][ 2+MAXNF] = fcx(fx[ 2+MAXNF]);
    }
  else if (opts.ih1 == -1)
    {
      creno_.cfx1_[i][-2+MAXNF] = fcx(fx[ 2+MAXNF]);
      creno_.cfx1_[i][-1+MAXNF] = fcx(fx[ 1+MAXNF]);
      creno_.cfx1_[i][ 1+MAXNF] = fcx(fx[-1+MAXNF]);
      creno_.cfx1_[i][ 2+MAXNF] = fcx(fx[-2+MAXNF]);
    }

  //beam 2 positive
  creno_.cfx2p_[i][-5+MAXNF] = fcx(fx[-5+MAXNF]);
  creno_.cfx2p_[i][-4+MAXNF] = fcx(fx[-4+MAXNF]);
  creno_.cfx2p_[i][-3+MAXNF] = fcx(fx[-3+MAXNF]);
  creno_.cfx2p_[i][ 0+MAXNF] = fcx(fx[ 0+MAXNF]);
  creno_.cfx2p_[i][ 3+MAXNF] = fcx(fx[ 3+MAXNF]);
  creno_.cfx2p_[i][ 4+MAXNF] = fcx(fx[ 4+MAXNF]);
  creno_.cfx2p_[i][ 5+MAXNF] = fcx(fx[ 5+MAXNF]);
  if (opts.ih2 == 1)
    {
      creno_.cfx2p_[i][-2+MAXNF] = fcx(fx[-2+MAXNF]);
      creno_.cfx2p_[i][-1+MAXNF] = fcx(fx[-1+MAXNF]);
      creno_.cfx2p_[i][ 1+MAXNF] = fcx(fx[ 1+MAXNF]);
      creno_.cfx2p_[i][ 2+MAXNF] = fcx(fx[ 2+MAXNF]);
    }
  else if (opts.ih2 == -1)
    {
      creno_.cfx2p_[i][-2+MAXNF] = fcx(fx[ 2+MAXNF]);
      creno_.cfx2p_[i][-1+MAXNF] = fcx(fx[ 1+MAXNF]);
      creno_.cfx2p_[i][ 1+MAXNF] = fcx(fx[-1+MAXNF]);
      creno_.cfx2p_[i][ 2+MAXNF] = fcx(fx[-2+MAXNF]);
    }

  //beam 2 negative
  creno_.cfx2m_[i][-5+MAXNF] = fcx(conj(fx[-5+MAXNF]));
  creno_.cfx2m_[i][-4+MAXNF] = fcx(conj(fx[-4+MAXNF]));
  creno_.cfx2m_[i][-3+MAXNF] = fcx(conj(fx[-3+MAXNF]));
  creno_.cfx2m_[i][ 0+MAXNF] = fcx(conj(fx[ 0+MAXNF]));
  creno_.cfx2m_[i][ 3+MAXNF] = fcx(conj(fx[ 3+MAXNF]));
  creno_.cfx2m_[i][ 4+MAXNF] = fcx(conj(fx[ 4+MAXNF]));
  creno_.cfx2m_[i][ 5+MAXNF] = fcx(conj(fx[ 5+MAXNF]));
  if (opts.ih2 == 1)
    {
      creno_.cfx2m_[i][-2+MAXNF] = fcx(conj(fx[-2+MAXNF]));
      creno_.cfx2m_[i][-1+MAXNF] = fcx(conj(fx[-1+MAXNF]));
      creno_.cfx2m_[i][ 1+MAXNF] = fcx(conj(fx[ 1+MAXNF]));
      creno_.cfx2m_[i][ 2+MAXNF] = fcx(conj(fx[ 2+MAXNF]));
    }
  else if (opts.ih2 == -1)
    {
      creno_.cfx2m_[i][-2+MAXNF] = fcx(conj(fx[ 2+MAXNF]));
      creno_.cfx2m_[i][-1+MAXNF] = fcx(conj(fx[ 1+MAXNF]));
      creno_.cfx2m_[i][ 1+MAXNF] = fcx(conj(fx[-1+MAXNF]));
      creno_.cfx2m_[i][ 2+MAXNF] = fcx(conj(fx[-2+MAXNF]));
    }
}

//Apply the flavour dependent form factors
void pdfevol::flavour_kt()
{
  int nf = 2*MAXNF+1;
  int negidx = mellinint::mdim*nf;
  complex <double> uv, dv, us, ds;
  for (int i = 0; i < mellinint::mdim; i++)
    {
      //beam 1 positive
      fx1[i*nf+bb] *= npff::boff;
      fx1[i*nf+cb] *= npff::chff;
      fx1[i*nf+sb] *= npff::ssff;
      fx1[i*nf+g ] *= npff::glff;
      fx1[i*nf+s ] *= npff::ssff;
      fx1[i*nf+c ] *= npff::chff;
      fx1[i*nf+b ] *= npff::boff;
      if (opts.ih1 == 1)
	{
	  uv = (fx1[i*nf+u] - fx1[i*nf+ub]);
	  dv = (fx1[i*nf+d] - fx1[i*nf+db]);
	  us = fx1[i*nf+ub];
	  ds = fx1[i*nf+db];
	  fx1[i*nf+u] = uv * npff::uvff + us * npff::usff;
	  fx1[i*nf+d] = dv * npff::dvff + ds * npff::dsff;
	  fx1[i*nf+ub] = us * npff::usff;
	  fx1[i*nf+db] = ds * npff::dsff;
	}
      else if (opts.ih1 == -1)
	{      
	  uv = (fx1[i*nf+ub] - fx1[i*nf+u]);
	  dv = (fx1[i*nf+db] - fx1[i*nf+d]);
	  us = fx1[i*nf+u];
	  ds = fx1[i*nf+d];
	  fx1[i*nf+ub] = uv * npff::uvff + us * npff::usff;
	  fx1[i*nf+db] = dv * npff::dvff + ds * npff::dsff;
	  fx1[i*nf+u] = us * npff::usff;
	  fx1[i*nf+d] = ds * npff::dsff;
	}

      //beam 2 positive
      fx2[i*nf+bb] *= npff::boff;
      fx2[i*nf+cb] *= npff::chff;
      fx2[i*nf+sb] *= npff::ssff;
      fx2[i*nf+g ] *= npff::glff;
      fx2[i*nf+s ] *= npff::ssff;
      fx2[i*nf+c ] *= npff::chff;
      fx2[i*nf+b ] *= npff::boff;
      if (opts.ih2 == 1)
	{
	  uv = (fx2[i*nf+u] - fx2[i*nf+ub]);
	  dv = (fx2[i*nf+d] - fx2[i*nf+db]);
	  us = fx2[i*nf+ub];
	  ds = fx2[i*nf+db];
	  fx2[i*nf+u] = uv * npff::uvff + us * npff::usff;
	  fx2[i*nf+d] = dv * npff::dvff + ds * npff::dsff;
	  fx2[i*nf+ub] = us * npff::usff;
	  fx2[i*nf+db] = ds * npff::dsff;
	}
      else if (opts.ih2 == -1)
	{      
	  uv = (fx2[i*nf+ub] - fx2[i*nf+u]);
	  dv = (fx2[i*nf+db] - fx2[i*nf+d]);
	  us = fx2[i*nf+u];
	  ds = fx2[i*nf+d];
	  fx2[i*nf+ub] = uv * npff::uvff + us * npff::usff;
	  fx2[i*nf+db] = dv * npff::dvff + ds * npff::dsff;
	  fx2[i*nf+u] = us * npff::usff;
	  fx2[i*nf+d] = ds * npff::dsff;
	}

      //beam 2 negative
      fx2[negidx+i*nf+bb] *= npff::boff;
      fx2[negidx+i*nf+cb] *= npff::chff;
      fx2[negidx+i*nf+sb] *= npff::ssff;
      fx2[negidx+i*nf+g ] *= npff::glff;
      fx2[negidx+i*nf+s ] *= npff::ssff;
      fx2[negidx+i*nf+c ] *= npff::chff;
      fx2[negidx+i*nf+b ] *= npff::boff;
      if (opts.ih2 == 1)
	{
	  uv = (fx2[negidx+i*nf+u] - fx2[negidx+i*nf+ub]);
	  dv = (fx2[negidx+i*nf+d] - fx2[negidx+i*nf+db]);
	  us = fx2[negidx+i*nf+ub];
	  ds = fx2[negidx+i*nf+db];
	  fx2[negidx+i*nf+u] = uv * npff::uvff + us * npff::usff;
	  fx2[negidx+i*nf+d] = dv * npff::dvff + ds * npff::dsff;
	  fx2[negidx+i*nf+ub] = us * npff::usff;
	  fx2[negidx+i*nf+db] = ds * npff::dsff;
	}
      else if (opts.ih2 == -1)
	{      
	  uv = (fx2[negidx+i*nf+ub] - fx2[negidx+i*nf+u]);
	  dv = (fx2[negidx+i*nf+db] - fx2[negidx+i*nf+d]);
	  us = fx2[negidx+i*nf+u];
	  ds = fx2[negidx+i*nf+d];
	  fx2[negidx+i*nf+ub] = uv * npff::uvff + us * npff::usff;
	  fx2[negidx+i*nf+db] = dv * npff::dvff + ds * npff::dsff;
	  fx2[negidx+i*nf+u] = us * npff::usff;
	  fx2[negidx+i*nf+d] = ds * npff::dsff;
	}
    }
}

void pdfevol::retrieve(int i1, int i2, int sign)
{
  int nf = 2*MAXNF+1;
  int negidx = mellinint::mdim*nf;

  fn1[bb] = fx1[i1*nf+bb];
  fn1[cb] = fx1[i1*nf+cb];
  fn1[sb] = fx1[i1*nf+sb];
  fn1[db] = fx1[i1*nf+db];
  fn1[ub] = fx1[i1*nf+ub];
  fn1[g ] = fx1[i1*nf+g ];
  fn1[u ] = fx1[i1*nf+u ];
  fn1[d ] = fx1[i1*nf+d ];
  fn1[s ] = fx1[i1*nf+s ];
  fn1[c ] = fx1[i1*nf+c ];
  fn1[b ] = fx1[i1*nf+b ];
  if (sign == mesq::positive)
    {
      fn2[bb] = fx2[i2*nf+bb];
      fn2[cb] = fx2[i2*nf+cb];
      fn2[sb] = fx2[i2*nf+sb];
      fn2[db] = fx2[i2*nf+db];
      fn2[ub] = fx2[i2*nf+ub];
      fn2[g ] = fx2[i2*nf+g ];
      fn2[u ] = fx2[i2*nf+u ];
      fn2[d ] = fx2[i2*nf+d ];
      fn2[s ] = fx2[i2*nf+s ];
      fn2[c ] = fx2[i2*nf+c ];
      fn2[b ] = fx2[i2*nf+b ];
    }
  else if (sign == mesq::negative)
    {
      fn2[bb] = fx2[negidx+i2*nf+bb];
      fn2[cb] = fx2[negidx+i2*nf+cb];
      fn2[sb] = fx2[negidx+i2*nf+sb];
      fn2[db] = fx2[negidx+i2*nf+db];
      fn2[ub] = fx2[negidx+i2*nf+ub];
      fn2[g ] = fx2[negidx+i2*nf+g ];
      fn2[u ] = fx2[negidx+i2*nf+u ];
      fn2[d ] = fx2[negidx+i2*nf+d ];
      fn2[s ] = fx2[negidx+i2*nf+s ];
      fn2[c ] = fx2[negidx+i2*nf+c ];
      fn2[b ] = fx2[negidx+i2*nf+b ];
    }
}  

void pdfevol::retrieve_fortran(int i1, int i2, int sign)
{
  //  cout << i1 << endl;
  //  cout << creno_.cfx1_[i1][5].real << "  " << creno_.cfx1_[i1][5].imag << endl;
  fn1[-5+MAXNF] = cx(creno_.cfx1_[i1][-5+MAXNF]);
  fn1[-4+MAXNF] = cx(creno_.cfx1_[i1][-4+MAXNF]);
  fn1[-3+MAXNF] = cx(creno_.cfx1_[i1][-3+MAXNF]);
  fn1[-2+MAXNF] = cx(creno_.cfx1_[i1][-2+MAXNF]);
  fn1[-1+MAXNF] = cx(creno_.cfx1_[i1][-1+MAXNF]);
  fn1[ 0+MAXNF] = cx(creno_.cfx1_[i1][ 0+MAXNF]);
  fn1[ 1+MAXNF] = cx(creno_.cfx1_[i1][ 1+MAXNF]);
  fn1[ 2+MAXNF] = cx(creno_.cfx1_[i1][ 2+MAXNF]);
  fn1[ 3+MAXNF] = cx(creno_.cfx1_[i1][ 3+MAXNF]);
  fn1[ 4+MAXNF] = cx(creno_.cfx1_[i1][ 4+MAXNF]);
  fn1[ 5+MAXNF] = cx(creno_.cfx1_[i1][ 5+MAXNF]);
  if (sign == mesq::positive)
    {
      fn2[-5+MAXNF] = cx(creno_.cfx2p_[i2][-5+MAXNF]);
      fn2[-4+MAXNF] = cx(creno_.cfx2p_[i2][-4+MAXNF]);
      fn2[-3+MAXNF] = cx(creno_.cfx2p_[i2][-3+MAXNF]);
      fn2[-2+MAXNF] = cx(creno_.cfx2p_[i2][-2+MAXNF]);
      fn2[-1+MAXNF] = cx(creno_.cfx2p_[i2][-1+MAXNF]);
      fn2[ 0+MAXNF] = cx(creno_.cfx2p_[i2][ 0+MAXNF]);
      fn2[ 1+MAXNF] = cx(creno_.cfx2p_[i2][ 1+MAXNF]);
      fn2[ 2+MAXNF] = cx(creno_.cfx2p_[i2][ 2+MAXNF]);
      fn2[ 3+MAXNF] = cx(creno_.cfx2p_[i2][ 3+MAXNF]);
      fn2[ 4+MAXNF] = cx(creno_.cfx2p_[i2][ 4+MAXNF]);
      fn2[ 5+MAXNF] = cx(creno_.cfx2p_[i2][ 5+MAXNF]);
    }
  else if (sign == mesq::negative)
    {
      fn2[-5+MAXNF] = cx(creno_.cfx2m_[i2][-5+MAXNF]);
      fn2[-4+MAXNF] = cx(creno_.cfx2m_[i2][-4+MAXNF]);
      fn2[-3+MAXNF] = cx(creno_.cfx2m_[i2][-3+MAXNF]);
      fn2[-2+MAXNF] = cx(creno_.cfx2m_[i2][-2+MAXNF]);
      fn2[-1+MAXNF] = cx(creno_.cfx2m_[i2][-1+MAXNF]);
      fn2[ 0+MAXNF] = cx(creno_.cfx2m_[i2][ 0+MAXNF]);
      fn2[ 1+MAXNF] = cx(creno_.cfx2m_[i2][ 1+MAXNF]);
      fn2[ 2+MAXNF] = cx(creno_.cfx2m_[i2][ 2+MAXNF]);
      fn2[ 3+MAXNF] = cx(creno_.cfx2m_[i2][ 3+MAXNF]);
      fn2[ 4+MAXNF] = cx(creno_.cfx2m_[i2][ 4+MAXNF]);
      fn2[ 5+MAXNF] = cx(creno_.cfx2m_[i2][ 5+MAXNF]);
    }
}

void pdfevol::retrieve_beam1(int i1)
{
  int nf = 2*MAXNF+1;
  //memcpy(fn1, &(fx1[i1*nf]), nf*sizeof(complex<double>));
  copy(fx1+i1*nf, fx1+(i1+1)*nf, fn1);
}

void pdfevol::retrieve_beam2_pos(int i2)
{
  int nf = 2*MAXNF+1;
  //memcpy(fn2, &(fx2[i2*nf]), nf*sizeof(complex<double>));
  copy(fx2+i2*nf, fx2+(i2+1)*nf, fn2);
}

void pdfevol::retrieve_beam2_neg()
{
  fn2[bb] = conj(fn2[bb]);
  fn2[cb] = conj(fn2[cb]);
  fn2[sb] = conj(fn2[sb]);
  fn2[db] = conj(fn2[db]);
  fn2[ub] = conj(fn2[ub]);
  fn2[g ] = conj(fn2[g ]);
  fn2[u ] = conj(fn2[u ]);
  fn2[d ] = conj(fn2[d ]);
  fn2[s ] = conj(fn2[s ]);
  fn2[c ] = conj(fn2[c ]);
  fn2[b ] = conj(fn2[b ]);
}


void pdfevol::retrieve1d_pos(int i)
{
  int nf = 2*MAXNF+1;

  //memcpy(fn1, &(fx1[i*nf]), nf*sizeof(complex<double>));
  //memcpy(fn2, &(fx2[i*nf]), nf*sizeof(complex<double>));

  copy(fx1+i*nf, fx1+(i+1)*nf, fn1);
  copy(fx2+i*nf, fx2+(i+1)*nf, fn2);
}

void pdfevol::retrieve1d_neg()
{
  //the negative branch is needed in mellin1d only for bprescription = 2
  //!!! This may be wrong when PDFs are evolved to complex scale (bprescription = 2) !!!
  fn1[bb] = conj(fn1[bb]);
  fn1[cb] = conj(fn1[cb]);
  fn1[sb] = conj(fn1[sb]);
  fn1[db] = conj(fn1[db]);
  fn1[ub] = conj(fn1[ub]);
  fn1[g ] = conj(fn1[g ]);
  fn1[u ] = conj(fn1[u ]);
  fn1[d ] = conj(fn1[d ]);
  fn1[s ] = conj(fn1[s ]);
  fn1[c ] = conj(fn1[c ]);
  fn1[b ] = conj(fn1[b ]);

  fn2[bb] = conj(fn2[bb]);
  fn2[cb] = conj(fn2[cb]);
  fn2[sb] = conj(fn2[sb]);
  fn2[db] = conj(fn2[db]);
  fn2[ub] = conj(fn2[ub]);
  fn2[g ] = conj(fn2[g ]);
  fn2[u ] = conj(fn2[u ]);
  fn2[d ] = conj(fn2[d ]);
  fn2[s ] = conj(fn2[s ]);
  fn2[c ] = conj(fn2[c ]);
  fn2[b ] = conj(fn2[b ]);

  //  int negidx = mellinint::mdim*nf;
  //  copy(fx1+negidx+i*nf, fx1+negidx+(i+1)*nf, fn1);
  //  copy(fx2+negidx+i*nf, fx2+negidx+(i+1)*nf, fn2);
}

void pdfevol::retrieve1d_fortran(int i, int sign)
{
  //cout << i << endl;
  //cout << creno_.cfx1_[i][5].real << "  " << creno_.cfx1_[i][5].imag << endl;
  if (sign == mesq::positive)
    {
      fn1[-5+MAXNF] = cx(creno_.cfx1_[i][-5+MAXNF]);
      fn1[-4+MAXNF] = cx(creno_.cfx1_[i][-4+MAXNF]);
      fn1[-3+MAXNF] = cx(creno_.cfx1_[i][-3+MAXNF]);
      fn1[-2+MAXNF] = cx(creno_.cfx1_[i][-2+MAXNF]);
      fn1[-1+MAXNF] = cx(creno_.cfx1_[i][-1+MAXNF]);
      fn1[ 0+MAXNF] = cx(creno_.cfx1_[i][ 0+MAXNF]);
      fn1[ 1+MAXNF] = cx(creno_.cfx1_[i][ 1+MAXNF]);
      fn1[ 2+MAXNF] = cx(creno_.cfx1_[i][ 2+MAXNF]);
      fn1[ 3+MAXNF] = cx(creno_.cfx1_[i][ 3+MAXNF]);
      fn1[ 4+MAXNF] = cx(creno_.cfx1_[i][ 4+MAXNF]);
      fn1[ 5+MAXNF] = cx(creno_.cfx1_[i][ 5+MAXNF]);
    }
  else if (sign == mesq::negative)
    {
      fn1[-5+MAXNF] = conj(cx(creno_.cfx1_[i][-5+MAXNF]));
      fn1[-4+MAXNF] = conj(cx(creno_.cfx1_[i][-4+MAXNF]));
      fn1[-3+MAXNF] = conj(cx(creno_.cfx1_[i][-3+MAXNF]));
      fn1[-2+MAXNF] = conj(cx(creno_.cfx1_[i][-2+MAXNF]));
      fn1[-1+MAXNF] = conj(cx(creno_.cfx1_[i][-1+MAXNF]));
      fn1[ 0+MAXNF] = conj(cx(creno_.cfx1_[i][ 0+MAXNF]));
      fn1[ 1+MAXNF] = conj(cx(creno_.cfx1_[i][ 1+MAXNF]));
      fn1[ 2+MAXNF] = conj(cx(creno_.cfx1_[i][ 2+MAXNF]));
      fn1[ 3+MAXNF] = conj(cx(creno_.cfx1_[i][ 3+MAXNF]));
      fn1[ 4+MAXNF] = conj(cx(creno_.cfx1_[i][ 4+MAXNF]));
      fn1[ 5+MAXNF] = conj(cx(creno_.cfx1_[i][ 5+MAXNF]));
    }
  if (sign == mesq::positive)
    {
      fn2[-5+MAXNF] = cx(creno_.cfx2p_[i][-5+MAXNF]);
      fn2[-4+MAXNF] = cx(creno_.cfx2p_[i][-4+MAXNF]);
      fn2[-3+MAXNF] = cx(creno_.cfx2p_[i][-3+MAXNF]);
      fn2[-2+MAXNF] = cx(creno_.cfx2p_[i][-2+MAXNF]);
      fn2[-1+MAXNF] = cx(creno_.cfx2p_[i][-1+MAXNF]);
      fn2[ 0+MAXNF] = cx(creno_.cfx2p_[i][ 0+MAXNF]);
      fn2[ 1+MAXNF] = cx(creno_.cfx2p_[i][ 1+MAXNF]);
      fn2[ 2+MAXNF] = cx(creno_.cfx2p_[i][ 2+MAXNF]);
      fn2[ 3+MAXNF] = cx(creno_.cfx2p_[i][ 3+MAXNF]);
      fn2[ 4+MAXNF] = cx(creno_.cfx2p_[i][ 4+MAXNF]);
      fn2[ 5+MAXNF] = cx(creno_.cfx2p_[i][ 5+MAXNF]);
    }
  else if (sign == mesq::negative)
    {
      fn2[-5+MAXNF] = cx(creno_.cfx2m_[i][-5+MAXNF]);
      fn2[-4+MAXNF] = cx(creno_.cfx2m_[i][-4+MAXNF]);
      fn2[-3+MAXNF] = cx(creno_.cfx2m_[i][-3+MAXNF]);
      fn2[-2+MAXNF] = cx(creno_.cfx2m_[i][-2+MAXNF]);
      fn2[-1+MAXNF] = cx(creno_.cfx2m_[i][-1+MAXNF]);
      fn2[ 0+MAXNF] = cx(creno_.cfx2m_[i][ 0+MAXNF]);
      fn2[ 1+MAXNF] = cx(creno_.cfx2m_[i][ 1+MAXNF]);
      fn2[ 2+MAXNF] = cx(creno_.cfx2m_[i][ 2+MAXNF]);
      fn2[ 3+MAXNF] = cx(creno_.cfx2m_[i][ 3+MAXNF]);
      fn2[ 4+MAXNF] = cx(creno_.cfx2m_[i][ 4+MAXNF]);
      fn2[ 5+MAXNF] = cx(creno_.cfx2m_[i][ 5+MAXNF]);
    }
}

//Retrieve PDFs at the starting scale (muf), only positive branch
// --> Check that this is doing the right thing, as it now depends on evolnative
void pdfevol::retrievemuf_1d(int i)
{
  // i is the index of the complex mellin moment in the z-space for the gaussian quadrature used for the mellin inversion

  //N flavour dependence
  int nf = resconst::NF;

  //XP[i] are moments of PDFs at the starting scale (factorisation scale)
  complex <double> fx[11];
  fx[0+MAXNF] = evolnative::GLP[i];
  fx[1+MAXNF] = evolnative::UVP[i] + evolnative::USP[i];
  fx[-1+MAXNF] = evolnative::USP[i];
  fx[2+MAXNF] = evolnative::DVP[i] + evolnative::DSP[i];
  fx[-2+MAXNF] = evolnative::DSP[i];
  fx[3+MAXNF] = evolnative::SVP[i] + evolnative::SSP[i];
  fx[-3+MAXNF] = evolnative::SSP[i];
  if (nf >= 4)
    {
      fx[4+MAXNF] = evolnative::CVP[i] + evolnative::CHP[i];
      fx[-4+MAXNF] = evolnative::CHP[i];
    }
  else
    {
      fx[4+MAXNF] = 0.;
      fx[-4+MAXNF] = 0.;
    }
  if (nf >= 5)
    {
      fx[5+MAXNF] = evolnative::BVP[i] + evolnative::BOP[i];
      fx[-5+MAXNF] = evolnative::BOP[i];
    }
  else
    {
      fx[5+MAXNF] = 0.;
      fx[-5+MAXNF] = 0.;
    }

  storemoments(i, fx);
  retrieve1d_pos(i);
  //  cout << i << "  " << GLP[i] << "  " << fx[0+MAXNF] << "  " << fn1[MAXNF] << "  " << fn2[MAXNF] << endl;
  return;
}

void pdfevol::retrievemuf_2d(int i)
{
  // i is the index of the complex mellin moment in the z-space for the gaussian quadrature used for the mellin inversion

  //N flavour dependence
  int nf = resconst::NF;

  //XP[i] are moments of PDFs at the starting scale (factorisation scale)
  complex <double> fx_1[11];
  fx_1[-5+MAXNF] = evolnative::BOP[i];
  fx_1[-4+MAXNF] = evolnative::CHP[i];
  fx_1[-3+MAXNF] = evolnative::SSP[i];
  fx_1[-2+MAXNF] = evolnative::DSP[i];
  fx_1[-1+MAXNF] = evolnative::USP[i];
  fx_1[ 0+MAXNF] = evolnative::GLP[i];
  fx_1[ 1+MAXNF] = evolnative::UVP[i] + evolnative::USP[i];
  fx_1[ 2+MAXNF] = evolnative::DVP[i] + evolnative::DSP[i];
  fx_1[ 3+MAXNF] = evolnative::SVP[i] + evolnative::SSP[i];
  fx_1[ 4+MAXNF] = evolnative::CVP[i] + evolnative::CHP[i];
  fx_1[ 5+MAXNF] = evolnative::BVP[i] + evolnative::BOP[i];

  complex <double> fx_2[11];
  fx_2[-5+MAXNF] = evolnative::BOP[i+mellinint::mdim];
  fx_2[-4+MAXNF] = evolnative::CHP[i+mellinint::mdim];
  fx_2[-3+MAXNF] = evolnative::SSP[i+mellinint::mdim];
  fx_2[-2+MAXNF] = evolnative::DSP[i+mellinint::mdim];
  fx_2[-1+MAXNF] = evolnative::USP[i+mellinint::mdim];
  fx_2[ 0+MAXNF] = evolnative::GLP[i+mellinint::mdim];
  fx_2[ 1+MAXNF] = evolnative::UVP[i+mellinint::mdim] + evolnative::USP[i+mellinint::mdim];
  fx_2[ 2+MAXNF] = evolnative::DVP[i+mellinint::mdim] + evolnative::DSP[i+mellinint::mdim];
  fx_2[ 3+MAXNF] = evolnative::SVP[i+mellinint::mdim] + evolnative::SSP[i+mellinint::mdim];
  fx_2[ 4+MAXNF] = evolnative::CVP[i+mellinint::mdim] + evolnative::CHP[i+mellinint::mdim];
  fx_2[ 5+MAXNF] = evolnative::BVP[i+mellinint::mdim] + evolnative::BOP[i+mellinint::mdim];
  
  if (nf < 4)
    fx_1[4+MAXNF] = fx_1[-4+MAXNF] = fx_2[4+MAXNF] = fx_2[-4+MAXNF] = 0.;
  if (nf < 5)
    fx_1[5+MAXNF] = fx_1[-5+MAXNF] = fx_2[5+MAXNF] = fx_2[-5+MAXNF] = 0.;

  storemoments_1(i, fx_1);
  storemoments_2(i, fx_2);
  //  cout << i << "  " << GLP[i] << "  " << fx[0+MAXNF] << "  " << fn1[MAXNF] << "  " << fn2[MAXNF] << endl;
  return;
}

void pdfevol::truncate()
{
  //Calculate truncated moments
  double x1 = phasespace::m/opts.sroot*exp(phasespace::ymin);
  double x2 = phasespace::m/opts.sroot*exp(-phasespace::ymax);

  //double x1 = 1e-8;//pow(phasespace::m/opts.sroot,2);
  //double x2 = 1e-8;//pow(phasespace::m/opts.sroot,2);

  double lx1 = log(x1);
  double lx2 = log(x2);
  
  //truncated moments (output)
  complex <double> ffx1[mellinint::mdim][2*MAXNF+1] = {0.};
  complex <double> ffx2[mellinint::mdim][2*MAXNF+1] = {0.};

  //cache x^(N) values
  complex <double> x1n[mellinint::mdim];
  complex <double> x2n[mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    {
      x1n[n] = pow(x1,mellinint::Np[n]);
      x2n[n] = pow(x2,mellinint::Np[n]);
    }

  //Normalisation times Jacobian
  complex <double> facp = mellinint::CCp/2./M_PI/complex <double>(0.,1);
  complex <double> facm = mellinint::CCm/2./M_PI/complex <double>(0.,1);
  
  //original moments times prefactor and weight
  int nf = 2*MAXNF+1;
  int negidx = mellinint::mdim*11;
  complex <double> fm1p[mellinint::mdim][2*MAXNF+1];
  complex <double> fm2p[mellinint::mdim][2*MAXNF+1];
  complex <double> fm1m[mellinint::mdim][2*MAXNF+1];
  complex <double> fm2m[mellinint::mdim][2*MAXNF+1];
  for (int m = 0; m < mellinint::mdim; m++)
    for (int f = 0; f < 2*MAXNF+1; f++)
      {
	//fm1p[m][f] = facp * cx(creno_.cfx1_[m][f] ) * mellinint::wn[m];
	//fm2p[m][f] = facp * cx(creno_.cfx2p_[m][f]) * mellinint::wn[m];
	//fm1m[m][f] = facm * conj(cx(creno_.cfx1_[m][f]) ) * mellinint::wn[m];
	//fm2m[m][f] = facm * conj(cx(creno_.cfx2p_[m][f])) * mellinint::wn[m];
	fm1p[m][f] = facp * fx1[m*nf+f] * mellinint::wn[m];
	fm2p[m][f] = facp * fx2[m*nf+f] * mellinint::wn[m];
	fm1m[m][f] = facm * conj(fx1[m*nf+f]) * mellinint::wn[m];
	fm2m[m][f] = facm * conj(fx2[m*nf+f]) * mellinint::wn[m];
      }

  //cache factor (1-x^(N-M))/(N-M) which limit is ln(x) when N-M -> 0
  complex <double> llx1p[mellinint::mdim][mellinint::mdim];
  complex <double> llx2p[mellinint::mdim][mellinint::mdim];
  complex <double> llx1m[mellinint::mdim][mellinint::mdim];
  complex <double> llx2m[mellinint::mdim][mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	llx1p[n][m] = (1.-x1n[n]/x1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llx2p[n][m] = (1.-x2n[n]/x2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llx1m[n][m] = (1.-x1n[n]/conj(x1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	llx2m[n][m] = (1.-x2n[n]/conj(x2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
      }

  //overwrite divergent diagonal part
  for (int n = 0; n < mellinint::mdim; n++)
    {
      llx1p[n][n] = -lx1;
      llx2p[n][n] = -lx2;
    }
  
  for (int n = 0; n < mellinint::mdim; n++)
    {
      //positive branch
      for (int m = 0; m < mellinint::mdim; m++)
	for (int f = 0; f < 2*MAXNF+1; f++)
	  {
	    /*
	    if (m == n)
	      {
		ffx1[n][f] += fm1p[m][f] * (-lx1);
		ffx2[n][f] += fm2p[m][f] * (-lx2);
	      }
	    else
	      {
		ffx1[n][f] += fm1p[m][f] * (1.-x1n[n]/x1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
		ffx2[n][f] += fm2p[m][f] * (1.-x2n[n]/x2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	      }
	    */
	    //	    ffx1[n][f] += fm1p[m][f]*llx1p[n][m]; 
	    //	    ffx2[n][f] += fm2p[m][f]*llx2p[n][m];
	    ffx1[n][f] += fm1p[m][f]*llx1p[n][m] - fm1m[m][f]*llx1m[n][m]; 
	    ffx2[n][f] += fm2p[m][f]*llx2p[n][m] - fm2m[m][f]*llx2m[n][m];
	  }

      /*
      //negative branch
      for (int m = 0; m < mellinint::mdim; m++)
	for (int f = 0; f < 2*MAXNF+1; f++)
	  {
	    //	    ffx1[n][f] -= fm1m[m][f]*(1.-x1n[n]/conj(x1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    //	    ffx2[n][f] -= fm2m[m][f]*(1.-x2n[n]/conj(x2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    ffx1[n][f] -= fm1m[m][f]*llx1m[n][m];
	    ffx2[n][f] -= fm2m[m][f]*llx2m[n][m];
	  }
      */

  
      //cout << "truncated " << n << ffx1[5] << endl;
    }

  //replace moments
  for (int n = 0; n < mellinint::mdim; n++)
    storemoments(n, ffx1[n]);
  //storemoments(n, ffx1[n], ffx2[n]);
}


void pdfevol::uppertruncate()
{
  //Calculate truncated moments
  double xmin1 = pow(phasespace::m/opts.sroot,2);
  double xmin2 = pow(phasespace::m/opts.sroot,2);

  double xmax1 = phasespace::m/opts.sroot*exp(phasespace::ymax);
  double xmax2 = phasespace::m/opts.sroot*exp(-phasespace::ymin);

  double lx1 = log(xmin1/xmax1);
  double lx2 = log(xmin2/xmax2);
  
  //truncated moments (output)
  complex <double> fx1[mellinint::mdim][2*MAXNF+1] = {0.};
  complex <double> fx2[mellinint::mdim][2*MAXNF+1] = {0.};

  //cache x^(N) values
  complex <double> xmin1n[mellinint::mdim];
  complex <double> xmin2n[mellinint::mdim];
  complex <double> xmax1n[mellinint::mdim];
  complex <double> xmax2n[mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    {
      xmin1n[n] = pow(xmin1,mellinint::Np[n]);
      xmin2n[n] = pow(xmin2,mellinint::Np[n]);
      xmax1n[n] = pow(xmax1,mellinint::Np[n]);
      xmax2n[n] = pow(xmax2,mellinint::Np[n]);
    }

  //Normalisation times Jacobian
  complex <double> facp = mellinint::CCp/2./M_PI/complex <double>(0.,1);
  complex <double> facm = mellinint::CCm/2./M_PI/complex <double>(0.,1);
  
  //original moments times prefactor and weight
  complex <double> fm1p[mellinint::mdim][2*MAXNF+1];
  complex <double> fm2p[mellinint::mdim][2*MAXNF+1];
  complex <double> fm1m[mellinint::mdim][2*MAXNF+1];
  complex <double> fm2m[mellinint::mdim][2*MAXNF+1];
  for (int m = 0; m < mellinint::mdim; m++)
    for (int f = 0; f < 2*MAXNF+1; f++)
      {
	fm1p[m][f] = facp * cx(creno_.cfx1_[m][f] ) * mellinint::wn[m];
	fm2p[m][f] = facp * cx(creno_.cfx2p_[m][f]) * mellinint::wn[m];
	fm1m[m][f] = facm * conj(cx(creno_.cfx1_[m][f]) ) * mellinint::wn[m];
	fm2m[m][f] = facm * conj(cx(creno_.cfx2p_[m][f])) * mellinint::wn[m];
      }

  //cache factor (xmax^(N-M)-xmin^(N-M))/(N-M) which limit is ln(xmin/xmax) when N-M -> 0
  complex <double> llx1p[mellinint::mdim][mellinint::mdim];
  complex <double> llx2p[mellinint::mdim][mellinint::mdim];
  complex <double> llx1m[mellinint::mdim][mellinint::mdim];
  complex <double> llx2m[mellinint::mdim][mellinint::mdim];
  for (int n = 0; n < mellinint::mdim; n++)
    for (int m = 0; m < mellinint::mdim; m++)
      {
	llx1p[n][m] = (xmax1n[n]/xmax1n[m]-xmin1n[n]/xmin1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llx2p[n][m] = (xmax2n[n]/xmax2n[m]-xmin2n[n]/xmin2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	llx1m[n][m] = (xmax1n[n]/conj(xmax1n[m])-xmin1n[n]/conj(xmin1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	llx2m[n][m] = (xmax2n[n]/conj(xmax2n[m])-xmin2n[n]/conj(xmin2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
      }

  //overwrite divergent diagonal part
  for (int n = 0; n < mellinint::mdim; n++)
    {
      llx1p[n][n] = -lx1;
      llx2p[n][n] = -lx2;
    }
  
  for (int n = 0; n < mellinint::mdim; n++)
    {
      //positive branch
      for (int m = 0; m < mellinint::mdim; m++)
	for (int f = 0; f < 2*MAXNF+1; f++)
	  {
	    /*
	    if (m == n)
	      {
		fx1[n][f] += fm1p[m][f] * (-lx1);
		fx2[n][f] += fm2p[m][f] * (-lx2);
	      }
	    else
	      {
		fx1[n][f] += fm1p[m][f] * (1.-x1n[n]/x1n[m])/(mellinint::Np[n]-mellinint::Np[m]);
		fx2[n][f] += fm2p[m][f] * (1.-x2n[n]/x2n[m])/(mellinint::Np[n]-mellinint::Np[m]);
	      }
	    */
	    //	    fx1[n][f] += fm1p[m][f]*llx1p[n][m]; 
	    //	    fx2[n][f] += fm2p[m][f]*llx2p[n][m];
	    fx1[n][f] += fm1p[m][f]*llx1p[n][m] - fm1m[m][f]*llx1m[n][m]; 
	    fx2[n][f] += fm2p[m][f]*llx2p[n][m] - fm2m[m][f]*llx2m[n][m];
	  }

      /*
      //negative branch
      for (int m = 0; m < mellinint::mdim; m++)
	for (int f = 0; f < 2*MAXNF+1; f++)
	  {
	    //	    fx1[n][f] -= fm1m[m][f]*(1.-x1n[n]/conj(x1n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    //	    fx2[n][f] -= fm2m[m][f]*(1.-x2n[n]/conj(x2n[m]))/(mellinint::Np[n]-mellinint::Nm[m]);
	    fx1[n][f] -= fm1m[m][f]*llx1m[n][m];
	    fx2[n][f] -= fm2m[m][f]*llx2m[n][m];
	  }
      */

  
      //cout << "truncated " << n << fx1[5] << endl;
    }

  //replace moments
  for (int n = 0; n < mellinint::mdim; n++)
    storemoments(n, fx1[n]);
    //storemoments(n, fx1[n], fx2[n]);
}

/*
void pdfevol::invert(double x, double fx[2*MAXNF+1])
{
  fx = {0.};

  complex <double> fm[2*MAXNF+1];
  complex <double> fm2[2*MAXNF+1];

  double lx = log(x);
  
  //positive branch
  for (int m = 0; m < mellinint::mdim; m++)
    {
      complex <double> cexp = mellinint::CCp/2./M_PI/complex <double>(0.,1) * exp(-mellinint::Np[m] * lx);
      for (int f = 0; f < 2*MAXNF+1; f++)
	{
	  fm[f] = cx(creno_.cfx1_[m][f]);
	  fx[f] += cexp * fm[f] * mellinint::wn[m];
	}

  //negative branch
  for (int m = 0; m < mellinint::mdim; m++)
    {
      complex <double> cexm = mellinint::CCm/2./M_PI/complex <double>(0.,1) * exp(-mellinint::Nm[m] * lx);
      for (int f = 0; f < 2*MAXNF+1; f++)
	{
	  fm[f] = cx(creno_.cfx2m_[m][f]);
	fx1[f] -= mellinint::CCm/2./M_PI/complex <double>(0.,1) * conj(fm1[f])*(1.-pow(x1,mellinint::Np[n]-mellinint::Nm[m]))/(mellinint::Np[n]-mellinint::Nm[m]) * mellinint::wn[m];
	fx2[f] -= mellinint::CCm/2./M_PI/complex <double>(0.,1) * conj(fm2[f])*(1.-pow(x2,mellinint::Np[n]-mellinint::Nm[m]))/(mellinint::Np[n]-mellinint::Nm[m]) * mellinint::wn[m];
      }
  
  cout << "truncated " << n << fx1[5] << endl;


}
*/
