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
#include "mellinpdf.h"
#include "evolnum.h"
#include "settings.h"
#include "pdf.h"
#include "alphas.h"
#include "constants.h"
#include "scales.h"
#include "psi.h"
#include "phasespace.h"
#include "mesq.h"

#include "LHAPDF/LHAPDF.h"

using namespace constants;
using namespace resconst;

int pegasus::nff;
int pegasus::ivfns;
int pegasus::dim;

complex <double> *pegasus::gli;
complex <double> *pegasus::vai;
complex <double> *pegasus::m3i;
complex <double> *pegasus::m8i;
complex <double> *pegasus::m15i;
complex <double> *pegasus::m24i;
complex <double> *pegasus::sgi;
complex <double> *pegasus::p3i;
complex <double> *pegasus::p8i;
complex <double> *pegasus::p15i;
complex <double> *pegasus::p24i;

void pegasus::init_const()
{
  //Some constants and the two-dimensional Kronecker symbol
  rzeta_.zeta_[0] = constants::euler; //0.57721566490153;
  rzeta_.zeta_[1] = constants::zeta2; //1.644934066848226;
  rzeta_.zeta_[2] = constants::zeta3; //1.202056903159594;
  rzeta_.zeta_[3] = constants::zeta4; //1.082323233711138;
  rzeta_.zeta_[4] = constants::zeta5; //1.036927755143370;
  rzeta_.zeta_[5] = constants::zeta6; //1.017343061984449;

  kron2d_.d_[0][0] =  fcx(complex <double> (1., 0.));
  kron2d_.d_[1][0] =  fcx(complex <double> (0., 0.));
  kron2d_.d_[0][1] =  fcx(complex <double> (0., 0.));
  kron2d_.d_[1][1] =  fcx(complex <double> (1., 0.));
  
  // QCD colour factors
  colour_.ca_ = constants::CA; //3.;
  colour_.cf_ = constants::CF; //4./3.;
  colour_.tr_ = constants::TF; //0.5;

  //Initialise coefficients of the QCD beta function
  //betafct_();
  for (int nf = NFMIN; nf <= NFMAX; nf++)
    {
      pgbeta_.pgbeta0_[nf-NFMIN] = 4.  * (33.-2.*nf)/12.;
      pgbeta_.pgbeta1_[nf-NFMIN] = 16. * (153.-19.*nf)/24.;
      pgbeta_.pgbeta2_[nf-NFMIN] = 64. * (2857./128.-5033.*nf/1152.+325.*(nf*nf)/3456.);
      pgbeta_.pgbeta3_[nf-NFMIN] = 256.* ((149753./6. + 3564.*zeta3 + nf*(-1078361./162.-6508./27.*zeta3)+(nf*nf)*(50065./162.+6472./81.*zeta3)+1093./729.*(nf*nf*nf))/256.);
    }
  
  //number of steps for the Runge-Kutta integration of alphas beyond LO  
  aspar_.nastps_ = alphas_steps;
}

//Allocate PDF moments at the starting scale
void pegasus::allocate()
{
  gli  = new complex <double> [2*dim];
  vai  = new complex <double> [2*dim];
  m3i  = new complex <double> [2*dim];
  m8i  = new complex <double> [2*dim];
  m15i = new complex <double> [2*dim];
  m24i = new complex <double> [2*dim];
  sgi  = new complex <double> [2*dim];
  p3i  = new complex <double> [2*dim];
  p8i  = new complex <double> [2*dim];
  p15i = new complex <double> [2*dim];
  p24i = new complex <double> [2*dim];
}

void pegasus::free()
{
  //if (opts.evolmode == 0 || opts.evolmode == 2)
  //return;
  delete[] gli;
  delete[] vai;
  delete[] m3i;
  delete[] m8i;
  delete[] m15i;
  delete[] m24i;
  delete[] sgi;
  delete[] p3i;
  delete[] p8i;
  delete[] p15i;
  delete[] p24i;
}

void pegasus::release()
{
  if (opts.melup <= 1)
    free();
}

void pegasus::init_alphas()
{
  //Set the initial scale for alphas M20 = M_0^2 (in GeV^2) and ASI = alpha_s(M_0^2)

  //In the backward evolution mode the Pegasus alphas is not used
  if (opts.evolmode == 1)
    return;

  //In the forward evolution mode set input values from LHAPDF

  //LHAPDF::PDFInfo info(opts.LHAPDFset, opts.LHAPDFmember);
  //double qmin = info.get_entry_as<double>("QMin", -1);
  //asinp_.m20_ = pow(qmin,2);
  asinp_.m20_ = pow(pdf::qmin,2);
      
  //The heavy quark masses squared, input values from LHAPDF (kmux can be used to modify the matching scales)
  asfthr_.m2c_ = pow(pdf::mc*opts.kmuc,2);
  asfthr_.m2b_ = pow(pdf::mb*opts.kmub,2);
  asfthr_.m2t_ = pow(pdf::mt*opts.kmut,2);

  //Stop some nonsense
  if (asinp_.m20_ > asfthr_.m2c_) //--> Allow evolution with starting scale above mcharm
    {
      //cout << "Too high mu_0 for VFNS evolution. STOP" << endl;
      //exit (-1);
    }

  if (asfthr_.m2c_ > asfthr_.m2b_)
    {
      cout << "wrong charm-bottom mass hierarchy. stop" << endl;
      exit (-1);
    }
  if (asfthr_.m2b_ > asfthr_.m2t_)
    {
      cout << "Wrong bottom-top mass hierarchy. STOP" << endl;
      exit (-1);
    }

  //double ASI = pdf::alphas(sqrt(asinp_.m20_)); // / (4.* M_PI);
  double ASI = pdf::alphas(sqrt(asinp_.m20_)) / (4.* M_PI);

  if ((ASI*4.* M_PI) > 2. || (ASI*4.* M_PI) < 2.0e-2 )
    {
      cout << "alpha_s out of range: " << ASI << " STOP" << endl;
      exit (-1);
    }

  //For mu_r unequal mu_f  AS0  is different from the input parameter 
  //ASM = a_s(M_0^2). The VFNS evolution starts with n_f = 3 at M_0^2.

  double ASM = ASI;// / (4.* M_PI);
  double R20 = asinp_.m20_ * exp(-frrat_.logfr_);
  int NF = nff;
  if (ivfns != 0)
    {
      if (asinp_.m20_ > asfthr_.m2b_)
	NF = 5;
      else if (asinp_.m20_ > (asfthr_.m2c_-1e-3))
	NF = 4;
      else
	NF = 3;
    }
  asinp_.as0_ = as_(R20, asinp_.m20_, ASM, NF);

  //  cout<< " alphas at the starting scale " << asinp_.as0_ << endl;
  //  cout<< " alphas LHAPDF " << LHAPDF::alphasPDF(sqrt(asinp_.m20_))/ (4.* M_PI) << endl;
}


//Calculate Mellin moments
void pegasus::calc_mellin()
{
  //nnused_.nmax_ = dim;

  //Location of the Mellin inversion contour in the complex N plane 
  //Support points NA(K) for the gauss integration
  complex <double> cxn[dim*2] = {0.};
  if (opts.mellin1d)
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int i = 0; i < mellinint::mdim; i++)
	{
	  int idx = index(i,sign);
	  if (sign == mesq::positive)
	    cxn[idx] = mellinint::Np[i];
	  else
	    cxn[idx] = mellinint::Nm[i];
	}
  else
    for (int sign = mesq::positive; sign <= mesq::negative; sign++)
      for (int beam = 1; beam <= 2; beam++)
	for (int i = 0; i < mellinint::mdim; i++)
	  {
	    int idx = index(i,beam,sign);
	    if (beam == 1)
	      if (sign == mesq::positive)
		cxn[idx] = mellinint::Np_1[i];
	      else
		cxn[idx] = mellinint::Nm_1[i];
	    else
	      if (sign == mesq::positive)
		cxn[idx] = mellinint::Np_2[i];
	      else
		cxn[idx] = mellinint::Nm_2[i];
	  }

  for (int i = 0; i < dim*2; i++)
    {
//complex <double> cxn;
//if (opts.mellin1d)
//	cxn = mellinint::Np[i];
//else
//	if (i < mellinint::mdim)
//	  cxn = mellinint::Np_1[i];
//	else
//	  cxn = mellinint::Np_2[i-mellinint::mdim];
	  
      moms_.na_[i] = fcx(cxn[i]);

      //The lowest simple harmonic sums on the support points
      fcomplex N1 = fcx(cxn[i] + 1.);
      int one = 1;
      int two = 2;
      int three = 3;
      int four = 4;
      int five = 5;
      hsums_.s_[0][i] = fcx(cx(psi_(N1)) + constants::euler);
      hsums_.s_[1][i] = fcx(rzeta_.zeta_[1] - cx(dpsi_(N1,one)));
      hsums_.s_[2][i] = fcx(rzeta_.zeta_[2] + 0.5 * cx(dpsi_(N1,two)));
      hsums_.s_[3][i] = fcx(rzeta_.zeta_[3] - 1./6. * cx(dpsi_(N1,three)));
      hsums_.s_[4][i] = fcx(rzeta_.zeta_[4] + 1./24. * cx(dpsi_(N1,four)));
      hsums_.s_[5][i] = fcx(rzeta_.zeta_[5] - 1./120. * cx(dpsi_(N1,five)));

      //--> Check this simplified code before using
      //complex <double> N1 = mellinint::Np[i] + 1.;
      //hsums_.s_[0][i] = fcx(constants::euler + cpsi0(N1));
      //hsums_.s_[1][i] = fcx(constants::zeta2 - cpsi(1,N1));
      //hsums_.s_[2][i] = fcx(constants::zeta3 + 0.5 * cpsi(2,N1));
      //hsums_.s_[3][i] = fcx(constants::zeta4 - 1./6. * cpsi(3,N1));
      //hsums_.s_[4][i] = fcx(constants::zeta5 + 1./24. * cpsi(4,N1));
      //hsums_.s_[5][i] = fcx(constants::zeta6 - 1./120. * cpsi(5,N1));
    }

  //The N-space splitting functions, evolution operators etc.
  pns0mom_();
  psg0mom_();
  lsgmom_();

  if (order_.npord_ > 0 || opts.order >= 2)
    {
      pns1mom_();
      psg1mom_();
      usg1mom_();
    }
  if (itord_.nuord_ > 1) usg1hmom_();

  if (order_.npord_ > 1 || opts.order >= 3)
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
//      cout << setprecision(16) << endl;
//      cout << "pega " << mellinint::Np[0]
//	   << "  " << cx(asg2_.a2sg_[0][0][0])
//	   << "  " << cx(asg2_.a2sg_[0][1][0])
//	   << "  " << cx(asg2_.a2sg_[1][0][0])
//	   << endl;
//      asg2mom_form_();
//      cout << "form " << mellinint::Np[0]
//	   << "  " << cx(asg2_.a2sg_[0][0][0])
//	   << "  " << cx(asg2_.a2sg_[0][1][0])
//	   << "  " << cx(asg2_.a2sg_[1][0][0])
//	   << endl;
    }
}

void pegasus::init_pdf()
{
  //evolmode = 3,4 -> Pegasus forward evolution: input PDFs at the starting scale, with 3 flavours
  //evolmode = 1 -> reproduce DYRES evolution: input PDFs at the factorisation scale, with 5 flavours
  
  // --> change to the new code for the mellin moments (which allows asymmetric charm and bottom)
  // --> Need to figure out possible memory issues with allocate and free
  //  mellinpdf::allocate();
  //  double facscale = sqrt(asinp_.m20_);
  //  mellinpdf::evalpdfs(facscale);
  //  if (opts.mellininv == 1 || opts.phi > 0.5)
  //    mellinpdf::laguerre_ipol();
  //  else
  //    mellinpdf::gauss_quad();
  //
  //for (int i = 0; i < mellinint::mdim; i++)
  // ....
  //
  //  mellinpdf::free();

  // --> also need to update PDFs at the starting scale in evolmode = 1 when muf = mll

  /// Do not separate anymore the fmufac = 0 case !!!
  //if (opts.evolmode == 1 && opts.fmufac == 0)
  //  {
  //    update();
  //    store();
  //  }
  
  
  //Pegasus forward evolution: input PDFs at the starting scale Q0, with 3 (or 4 or 5) flavours
  bool threval = false; //if true then perform explicit Mellin transform at each threshold
  if (opts.evolmode == 3 || opts.evolmode == 4)
    {
      //No need to init mellinpdf as this was called by pdfevol
      //mellinpdf::init();
      
      double facscale = sqrt(asinp_.m20_);

      mellinpdf::allocate();
      mellinpdf::evalpdfs(facscale);
      mellinpdf::transform();
      
      complex <double> qp[5];
      complex <double> qm[5];
      for (int sign = mesq::positive; sign <= mesq::negative; sign++)
	for (int n = 0; n < dim; n++)
	  {
	    int i = index(n,sign);
	    if (opts.mellin1d)
	      {
		//gluon
	      painp_.gli_[i] = fcx(mellinpdf::GL[n]);
      
	      //Arrays of non-singlet and singlet quark combinations for N_f = 3 (and 4 and 5)
	      //defined as in Eq. (2.16) of hep-ph/0408244
	      qp[0] = mellinpdf::UP[n] + mellinpdf::UB[n];
	      qp[1] = mellinpdf::DO[n] + mellinpdf::DB[n];
	      qp[2] = mellinpdf::ST[n] + mellinpdf::SB[n];
	      qp[3] = mellinpdf::CH[n] + mellinpdf::CB[n];
	      qp[4] = mellinpdf::BO[n] + mellinpdf::BB[n];
	      
	      qm[0] = mellinpdf::UP[n] - mellinpdf::UB[n];
	      qm[1] = mellinpdf::DO[n] - mellinpdf::DB[n];
	      qm[2] = mellinpdf::ST[n] - mellinpdf::SB[n];
	      qm[3] = mellinpdf::CH[n] - mellinpdf::CB[n];
	      qm[4] = mellinpdf::BO[n] - mellinpdf::BB[n];
	    }
	  else
	    if (n < mellinint::mdim)
	      {
		painp_.gli_[i] = fcx(mellinpdf::GL_1[n]);
		qp[0] = mellinpdf::UP_1[n] + mellinpdf::UB_1[n];
		qp[1] = mellinpdf::DO_1[n] + mellinpdf::DB_1[n];
		qp[2] = mellinpdf::ST_1[n] + mellinpdf::SB_1[n];
		qp[3] = mellinpdf::CH_1[n] + mellinpdf::CB_1[n];
		qp[4] = mellinpdf::BO_1[n] + mellinpdf::BB_1[n];
		qm[0] = mellinpdf::UP_1[n] - mellinpdf::UB_1[n];
		qm[1] = mellinpdf::DO_1[n] - mellinpdf::DB_1[n];
		qm[2] = mellinpdf::ST_1[n] - mellinpdf::SB_1[n];
		qm[3] = mellinpdf::CH_1[n] - mellinpdf::CB_1[n];
		qm[4] = mellinpdf::BO_1[n] - mellinpdf::BB_1[n];
	      }
	    else
	      {
		painp_.gli_[i] = fcx(mellinpdf::GL_2[n-mellinint::mdim]);
		qp[0] = mellinpdf::UP_2[n-mellinint::mdim] + mellinpdf::UB_2[n-mellinint::mdim];
		qp[1] = mellinpdf::DO_2[n-mellinint::mdim] + mellinpdf::DB_2[n-mellinint::mdim];
		qp[2] = mellinpdf::ST_2[n-mellinint::mdim] + mellinpdf::SB_2[n-mellinint::mdim];
		qp[3] = mellinpdf::CH_2[n-mellinint::mdim] + mellinpdf::CB_2[n-mellinint::mdim];
		qp[4] = mellinpdf::BO_2[n-mellinint::mdim] + mellinpdf::BB_2[n-mellinint::mdim];
		qm[0] = mellinpdf::UP_2[n-mellinint::mdim] - mellinpdf::UB_2[n-mellinint::mdim];
		qm[1] = mellinpdf::DO_2[n-mellinint::mdim] - mellinpdf::DB_2[n-mellinint::mdim];
		qm[2] = mellinpdf::ST_2[n-mellinint::mdim] - mellinpdf::SB_2[n-mellinint::mdim];
		qm[3] = mellinpdf::CH_2[n-mellinint::mdim] - mellinpdf::CB_2[n-mellinint::mdim];
		qm[4] = mellinpdf::BO_2[n-mellinint::mdim] - mellinpdf::BB_2[n-mellinint::mdim];
	      }

	  int nf = 3;
	  if (asinp_.m20_ > (asfthr_.m2c_-1e-3))
	    nf = 4;
	  if (asinp_.m20_ > asfthr_.m2b_)
	    nf = 5;

	  complex <double> sumqm = 0.;
	  for (int f = 0; f < nf; f++)
	    sumqm += qm[f];

	  complex <double> sumqp = 0.;
	  for (int f = 0; f < nf; f++)
	    sumqp += qp[f];
	  
	  painp_.vai_[i] = fcx(sumqm);
	  painp_.m3i_[i] = fcx(qm[0]-qm[1]);
	  painp_.m8i_[i] = fcx(qm[0]+qm[1]-2.*qm[2]);
	  if (nf >= 4)
	    hfpainp_.m15i_[i] = fcx(qm[0]+qm[1]+qm[2]-3.*qm[3]);
	  if (nf >= 5)
	    hfpainp_.m24i_[i] = fcx(qm[0]+qm[1]+qm[2]+qm[3]-4.*qm[4]);
  
	  painp_.sgi_[i] = fcx(sumqp);
	  painp_.p3i_[i] = fcx(qp[0]-qp[1]);
	  painp_.p8i_[i] = fcx(qp[0]+qp[1]-2.*qp[2]);
	  if (nf >= 4)
	    hfpainp_.p15i_[i] = fcx(qp[0]+qp[1]+qp[2]-3.*qp[3]);
	  if (nf >= 5)
	    hfpainp_.p24i_[i] = fcx(qp[0]+qp[1]+qp[2]+qp[3]-4.*qp[4]);
	}
      //store();

      //Perform explicit Mellin transform at each threshold --> not supported
      if (threval)
	{
	  //c threshold
	  facscale = sqrt(asfthr_.m2c_);
	  mellinpdf::evalpdfs(facscale);
	  mellinpdf::transform();
	  
	  for (int i = 0; i < dim; i++)
	    {
	      if (opts.mellin1d)
		{
		  //gluon
		  pacthr_.glc_[i] = fcx(mellinpdf::GL[i]);
      
		  //Arrays of non-singlet and singlet quark combinations for N_f = 3 (and 4 and 5)
		  //defined as in Eq. (2.16) of hep-ph/0408244
		  qp[0] = mellinpdf::UP[i] + mellinpdf::UB[i];
		  qp[1] = mellinpdf::DO[i] + mellinpdf::DB[i];
		  qp[2] = mellinpdf::ST[i] + mellinpdf::SB[i];
		  qp[3] = mellinpdf::CH[i] + mellinpdf::CB[i];
		  qp[4] = mellinpdf::BO[i] + mellinpdf::BB[i];
	      
		  qm[0] = mellinpdf::UP[i] - mellinpdf::UB[i];
		  qm[1] = mellinpdf::DO[i] - mellinpdf::DB[i];
		  qm[2] = mellinpdf::ST[i] - mellinpdf::SB[i];
		  qm[3] = mellinpdf::CH[i] - mellinpdf::CB[i];
		  qm[4] = mellinpdf::BO[i] - mellinpdf::BB[i];
		}

	      int nf = 4;
	      
	      complex <double> sumqm = 0.;
	      for (int f = 0; f < nf; f++)
		sumqm += qm[f];
	      
	      complex <double> sumqp = 0.;
	      for (int f = 0; f < nf; f++)
		sumqp += qp[f];
	  
	      pacthr_.vac_[i] = fcx(sumqm);
	      pacthr_.m3c_[i] = fcx(qm[0]-qm[1]);
	      pacthr_.m8c_[i] = fcx(qm[0]+qm[1]-2.*qm[2]);
	      pacthr_.m15c_[i] = fcx(qm[0]+qm[1]+qm[2]-3.*qm[3]);
		  
	      pacthr_.sgc_[i] = fcx(sumqp);
	      pacthr_.p3c_[i] = fcx(qp[0]-qp[1]);
	      pacthr_.p8c_[i] = fcx(qp[0]+qp[1]-2.*qp[2]);
	      pacthr_.p15c_[i] = fcx(qp[0]+qp[1]+qp[2]-3.*qp[3]);
	    }

	  //b threshold
	  facscale = sqrt(asfthr_.m2b_);
	  mellinpdf::evalpdfs(facscale);
	  mellinpdf::transform();

	  for (int i = 0; i < dim; i++)
	    {
	      if (opts.mellin1d)
		{
		  //gluon
		  pabthr_.glb_[i] = fcx(mellinpdf::GL[i]);
      
		  //Arrays of non-singlet and singlet quark combinations for N_f = 3 (and 4 and 5)
		  //defined as in Eq. (2.16) of hep-ph/0408244
		  qp[0] = mellinpdf::UP[i] + mellinpdf::UB[i];
		  qp[1] = mellinpdf::DO[i] + mellinpdf::DB[i];
		  qp[2] = mellinpdf::ST[i] + mellinpdf::SB[i];
		  qp[3] = mellinpdf::CH[i] + mellinpdf::CB[i];
		  qp[4] = mellinpdf::BO[i] + mellinpdf::BB[i];
	      
		  qm[0] = mellinpdf::UP[i] - mellinpdf::UB[i];
		  qm[1] = mellinpdf::DO[i] - mellinpdf::DB[i];
		  qm[2] = mellinpdf::ST[i] - mellinpdf::SB[i];
		  qm[3] = mellinpdf::CH[i] - mellinpdf::CB[i];
		  qm[4] = mellinpdf::BO[i] - mellinpdf::BB[i];
		}
	      int nf = 5;

	      complex <double> sumqm = 0.;
	      for (int f = 0; f < nf; f++)
		sumqm += qm[f];

	      complex <double> sumqp = 0.;
	      for (int f = 0; f < nf; f++)
		sumqp += qp[f];
	  
	      pabthr_.vab_[i] = fcx(sumqm);
	      pabthr_.m3b_[i] = fcx(qm[0]-qm[1]);
	      pabthr_.m8b_[i] = fcx(qm[0]+qm[1]-2.*qm[2]);
	      pabthr_.m15b_[i] = fcx(qm[0]+qm[1]+qm[2]-3.*qm[3]);
	      pabthr_.m24b_[i] = fcx(qm[0]+qm[1]+qm[2]+qm[3]-4.*qm[4]);
  
	      pabthr_.sgb_[i] = fcx(sumqp);
	      pabthr_.p3b_[i] = fcx(qp[0]-qp[1]);
	      pabthr_.p8b_[i] = fcx(qp[0]+qp[1]-2.*qp[2]);
	      pabthr_.p15b_[i] = fcx(qp[0]+qp[1]+qp[2]-3.*qp[3]);
	      pabthr_.p24b_[i] = fcx(qp[0]+qp[1]+qp[2]+qp[3]-4.*qp[4]);
	    }
	  //Evaluate also the alphas thresholds from LHAPDF
	  asfthr_.asc_ = pdf::alphas(sqrt(asfthr_.m2c_)) /(4.*M_PI);
	  asfthr_.asb_ = pdf::alphas(sqrt(asfthr_.m2b_)) /(4.*M_PI);
	  asfthr_.ast_ = pdf::alphas(sqrt(asfthr_.m2t_)) /(4.*M_PI);
	  evnasthr_(asfthr_.m2c_, asfthr_.m2b_, asfthr_.m2t_);
	}
      
      mellinpdf::free();
      //xmin = pow(bins.mbins.front()/opts.sroot,2); //Restrict the integration of moments to xmin = m/sqrt(s)*exp(-ymax) = (m/sqrt(s))^2
      //mellinpdf::init(xmin);

      /*      
      fcomplex uval,dval,ubar,dbar,s,sbar,glu,charm,bot;
      for (int i = 0; i < mellinint::mdim; i++)
	{
	  int hadron = 1; //opts.ih1;
	  fcomplex N = fcx(mellinint::Np[i]); //compute positive branch only, the negative branch is obtained by complex conjugation
	  double facscale = sqrt(asinp_.m20_);
	  double xmin = 1e-8;
	  pdfmoments_(hadron,facscale,N,uval,dval,ubar,dbar,s,sbar,glu,charm,bot,xmin);
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
	  painp_.vai_[i] = fcx(cx(uval) + cx(dval) + (cx(s)-cx(sbar)));
	  painp_.m3i_[i] = fcx(cx(uval) - cx(dval));
	  painp_.m8i_[i] = fcx(cx(painp_.vai_[i]) - 3.*(cx(s)-cx(sbar)));

	  painp_.sgi_[i] = fcx(cx(uval) + cx(dval) + 2.* (cx(dbar) + cx(ubar)) + (cx(s)+cx(sbar)));
	  painp_.p3i_[i] = fcx(cx(painp_.m3i_[i]) - 2.* (cx(dbar)-cx(ubar)));
	  painp_.p8i_[i] = fcx(cx(uval) + cx(dval) + 2.* (cx(dbar) + cx(ubar)) - 2.* (cx(s)+cx(sbar)));
	}
      */
      
    }

  //cout << "m0 " << asinp_.m20_  << "  " << asinp_.as0_  *4.*M_PI<< "  " << pdf::alphas(sqrt(asinp_.m20_)) << endl;
  //cout << "mc " << asfthr_.m2c_ << "  " << asfthr_.asc_ *4.*M_PI<< "  " << pdf::alphas(sqrt(asfthr_.m2c_)) << endl;
  //cout << "mb " << asfthr_.m2b_ << "  " << asfthr_.asb_ *4.*M_PI<< "  " << pdf::alphas(sqrt(asfthr_.m2b_)) << endl;
  //cout << "mt " << asfthr_.m2t_ << "  " << asfthr_.ast_ *4.*M_PI<< "  " << pdf::alphas(sqrt(asfthr_.m2t_)) << endl;

  //In VFN evolution calculate PDFs at the mc, mb, mt thresholds
  if (ivfns != 0 && !threval)
    evnfthr_(asfthr_.m2c_, asfthr_.m2b_, asfthr_.m2t_);
  //  else
  //    dyevnfthr_(asfthr_.m2c_, asfthr_.m2b_, asfthr_.m2t_);
  
}

void pegasus::init()
{
  //No need to initialise if pegasus is not used --> Now need to initialise gamma matrices in all cases
  //if (opts.evolmode != 1 && opts.evolmode != 3 && opts.evolmode != 4)
  //if (opts.evolmode == 0 || opts.evolmode == 2)
  //return;

  if (opts.mellin1d)
    dim = mellinint::mdim;
  else
    dim = 2*mellinint::mdim;

  nnused_.nmax_ = 2*dim;
  
  //if (opts.melup <= 1)
  //allocate();
  
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


  // Some default settings of the external initialization parameters
  // (standard-speed iterated VFNS evolution at NLO for mu_f/mu_r = 1)

  //mode of evolution in Pegasus-QCD
  //there are three available schemes for solving the evolution equations in N-space at NNLO, IMODEV = 1, 2, 3
  evmod_.imodev_ = evolution_mode;

  if (opts.evolmode == 1)
    {
      ivfns = 0; //FFN evolution
      nff = 5;   //5 light flavours
      order_.npord_ = opts.order_evol - 1; //order of evolution is LO, NLO, NNLO for NLL, NNLL, NNNLL
      //if (opts.order == 3) order_.npord_ = 1; //---> Use NLO evolution at N3LL
      //evmod_.imodev_ = 3; //reproduces DYRes evolution
      evmod_.imodev_ = 4; //reproduces DYRes evolution (even better)
    }
  //reproduce in N-space the LHAPDF evolution in x-space
  else if (opts.evolmode == 3 || opts.evolmode == 4)
    {
      ivfns = 1; //VFN evolution (read from LHAPDF)
      nff = 4;   //in FFN evolution, number of flavours  --> should read from LHAPDF
      order_.npord_ = pdf::order; //order of evolution read from LHAPDF
      evmod_.imodev_ = 1; //reproduces evolution in x-space --> actually for NNPDF options 3 or 4 may be better
    }
  
  double FR2 = 1.;//ratio of muren2/mufac2
  frrat_.logfr_ = log(FR2);

  //Some more (internal) evolution and initialization parameters 
  aspar_.naord_  = order_.npord_;
  
  //maximal power of alphas in the U-matrix solution, Eq. (2.23) of https://arxiv.org/pdf/hep-ph/0408244.pdf
  if (evmod_.imodev_ == 1 || evmod_.imodev_ == 2)
    itord_.nuord_ = 15;            //Iterative solution with 15 iterations
  else
    itord_.nuord_ = order_.npord_; //Truncated solution

  /*
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
  */
  nfused_.nflow_  = NFMIN;
  nfused_.nfhigh_ = NFMAX;

  if (opts.melup == 0)
    calc_mellin();

  //No need to go beyond this point if Pegasus evolution is not used
  if (opts.evolmode == 0 || opts.evolmode == 2)
    return;
  
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
  init_alphas();
  
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

  if (opts.melup <= 1)
    {
    //init_pdf();
      allocate();
      update();
    }
}

//update PDFs at the starting scale
void pegasus::update()
{
  //No need to update PDFs at the starting scale for forward evolution modes, where the starting scale is Q0
  //if (opts.evolmode != 1 && opts.evolmode != 3)
  //return;

  //if (opts.evolmode == 1 && opts.fmufac == 0) //This would not work when the Mellin inversion points are updated
  //return;

  double facscale;
  if (opts.evolmode == 1)
    facscale = scales::fac;
  else if (opts.evolmode == 3 || opts.evolmode == 4)
    facscale = sqrt(asinp_.m20_);
  
  //Assume the following was already called
  //pdfevol::allocate();
  //scales::set(opts.rmass);
  //pdfevol::update();
  mellinpdf::allocate();
  mellinpdf::evalpdfs(facscale, phasespace::m, phasespace::y);
  mellinpdf::update_mellin();
  mellinpdf::transform();

  complex <double> qp[5];
  complex <double> qm[5];
  complex <double> gl;
  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
    for (int n = 0; n < dim; n++)
      {
	int nn;
	if (opts.mellin1d)
	  nn = n;
	else
	  if (n < mellinint::mdim)
	    nn = n;
	  else
	    nn = n-mellinint::mdim;
	
	if (opts.mellin1d)
	  if (sign == mesq::positive)
	    {
	      gl    = mellinpdf::GL[nn];
	      qp[0] = mellinpdf::UP[nn] + mellinpdf::UB[nn];
	      qp[1] = mellinpdf::DO[nn] + mellinpdf::DB[nn];
	      qp[2] = mellinpdf::ST[nn] + mellinpdf::SB[nn];
	      qp[3] = mellinpdf::CH[nn] + mellinpdf::CB[nn];
	      qp[4] = mellinpdf::BO[nn] + mellinpdf::BB[nn];
	      qm[0] = mellinpdf::UP[nn] - mellinpdf::UB[nn];
	      qm[1] = mellinpdf::DO[nn] - mellinpdf::DB[nn];
	      qm[2] = mellinpdf::ST[nn] - mellinpdf::SB[nn];
	      qm[3] = mellinpdf::CH[nn] - mellinpdf::CB[nn];
	      qm[4] = mellinpdf::BO[nn] - mellinpdf::BB[nn];
	    }
	  else
	    {
	      gl    = conj(mellinpdf::GL[nn]);
	      qp[0] = conj(mellinpdf::UP[nn] + mellinpdf::UB[nn]);
	      qp[1] = conj(mellinpdf::DO[nn] + mellinpdf::DB[nn]);
	      qp[2] = conj(mellinpdf::ST[nn] + mellinpdf::SB[nn]);
	      qp[3] = conj(mellinpdf::CH[nn] + mellinpdf::CB[nn]);
	      qp[4] = conj(mellinpdf::BO[nn] + mellinpdf::BB[nn]);
	      qm[0] = conj(mellinpdf::UP[nn] - mellinpdf::UB[nn]);
	      qm[1] = conj(mellinpdf::DO[nn] - mellinpdf::DB[nn]);
	      qm[2] = conj(mellinpdf::ST[nn] - mellinpdf::SB[nn]);
	      qm[3] = conj(mellinpdf::CH[nn] - mellinpdf::CB[nn]);
	      qm[4] = conj(mellinpdf::BO[nn] - mellinpdf::BB[nn]);
	    }
	else
	  if (n < mellinint::mdim)
	    if (sign == mesq::positive)
	      {
		gl    = mellinpdf::GL_1[nn];
		qp[0] = mellinpdf::UP_1[nn] + mellinpdf::UB_1[nn];
		qp[1] = mellinpdf::DO_1[nn] + mellinpdf::DB_1[nn];
		qp[2] = mellinpdf::ST_1[nn] + mellinpdf::SB_1[nn];
		qp[3] = mellinpdf::CH_1[nn] + mellinpdf::CB_1[nn];
		qp[4] = mellinpdf::BO_1[nn] + mellinpdf::BB_1[nn];
		qm[0] = mellinpdf::UP_1[nn] - mellinpdf::UB_1[nn];
		qm[1] = mellinpdf::DO_1[nn] - mellinpdf::DB_1[nn];
		qm[2] = mellinpdf::ST_1[nn] - mellinpdf::SB_1[nn];
		qm[3] = mellinpdf::CH_1[nn] - mellinpdf::CB_1[nn];
		qm[4] = mellinpdf::BO_1[nn] - mellinpdf::BB_1[nn];
	      }
	    else
	      {
		gl    = conj(mellinpdf::GL_1[nn]);
		qp[0] = conj(mellinpdf::UP_1[nn] + mellinpdf::UB_1[nn]);
		qp[1] = conj(mellinpdf::DO_1[nn] + mellinpdf::DB_1[nn]);
		qp[2] = conj(mellinpdf::ST_1[nn] + mellinpdf::SB_1[nn]);
		qp[3] = conj(mellinpdf::CH_1[nn] + mellinpdf::CB_1[nn]);
		qp[4] = conj(mellinpdf::BO_1[nn] + mellinpdf::BB_1[nn]);
		qm[0] = conj(mellinpdf::UP_1[nn] - mellinpdf::UB_1[nn]);
		qm[1] = conj(mellinpdf::DO_1[nn] - mellinpdf::DB_1[nn]);
		qm[2] = conj(mellinpdf::ST_1[nn] - mellinpdf::SB_1[nn]);
		qm[3] = conj(mellinpdf::CH_1[nn] - mellinpdf::CB_1[nn]);
		qm[4] = conj(mellinpdf::BO_1[nn] - mellinpdf::BB_1[nn]);
	      }
	  else
	    if (sign == mesq::positive)
	      {
		gl    = mellinpdf::GL_2[nn];
		qp[0] = mellinpdf::UP_2[nn] + mellinpdf::UB_2[nn];
		qp[1] = mellinpdf::DO_2[nn] + mellinpdf::DB_2[nn];
		qp[2] = mellinpdf::ST_2[nn] + mellinpdf::SB_2[nn];
		qp[3] = mellinpdf::CH_2[nn] + mellinpdf::CB_2[nn];
		qp[4] = mellinpdf::BO_2[nn] + mellinpdf::BB_2[nn];
		qm[0] = mellinpdf::UP_2[nn] - mellinpdf::UB_2[nn];
		qm[1] = mellinpdf::DO_2[nn] - mellinpdf::DB_2[nn];
		qm[2] = mellinpdf::ST_2[nn] - mellinpdf::SB_2[nn];
		qm[3] = mellinpdf::CH_2[nn] - mellinpdf::CB_2[nn];
		qm[4] = mellinpdf::BO_2[nn] - mellinpdf::BB_2[nn];
	      }
	    else
	      {
		gl    = conj(mellinpdf::GL_2[nn]);
		qp[0] = conj(mellinpdf::UP_2[nn] + mellinpdf::UB_2[nn]);
		qp[1] = conj(mellinpdf::DO_2[nn] + mellinpdf::DB_2[nn]);
		qp[2] = conj(mellinpdf::ST_2[nn] + mellinpdf::SB_2[nn]);
		qp[3] = conj(mellinpdf::CH_2[nn] + mellinpdf::CB_2[nn]);
		qp[4] = conj(mellinpdf::BO_2[nn] + mellinpdf::BB_2[nn]);
		qm[0] = conj(mellinpdf::UP_2[nn] - mellinpdf::UB_2[nn]);
		qm[1] = conj(mellinpdf::DO_2[nn] - mellinpdf::DB_2[nn]);
		qm[2] = conj(mellinpdf::ST_2[nn] - mellinpdf::SB_2[nn]);
		qm[3] = conj(mellinpdf::CH_2[nn] - mellinpdf::CB_2[nn]);
		qm[4] = conj(mellinpdf::BO_2[nn] - mellinpdf::BB_2[nn]);
	      }
	
	//arrays of non-singlet and singlet quark combinations for N_f = 3 (and 4 and 5)
	//defined as in Eq. (2.16) of hep-ph/0408244
	int i = index(n,sign);
	painp_.gli_[i] = fcx(gl);	
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
  mellinpdf::free();
  if (opts.evolmode == 3 || opts.evolmode == 4)
    store();

  //In VFN evolution calculate PDFs at the mc, mb, mt thresholds
  if (ivfns != 0)// && !threval)
    evnfthr_(asfthr_.m2c_, asfthr_.m2b_, asfthr_.m2t_);
}

//store PDFs at the starting scale
void pegasus::store()
{
  for (int i = 0; i < 2*dim; i++)
    {
      gli[i]  = cx(painp_.gli_[i]);
      vai[i]  = cx(painp_.vai_[i]);
      m3i[i]  = cx(painp_.m3i_[i]);
      m8i[i]  = cx(painp_.m8i_[i]);
      m15i[i] = cx(hfpainp_.m15i_[i]);
      m24i[i] = cx(hfpainp_.m24i_[i]);
      sgi[i]  = cx(painp_.sgi_[i]);
      p3i[i]  = cx(painp_.p3i_[i]);
      p8i[i]  = cx(painp_.p8i_[i]);
      p15i[i] = cx(hfpainp_.p15i_[i]);
      p24i[i] = cx(hfpainp_.p24i_[i]);
    }
}

//retrieve PDFs at the starting scale
void pegasus::retrieve()
{
  for (int i = 0; i < 2*dim; i++)
    {
      painp_.gli_[i]    = fcx(gli[i]);
      painp_.vai_[i]    = fcx(vai[i]);
      painp_.m3i_[i]    = fcx(m3i[i]);
      painp_.m8i_[i]    = fcx(m8i[i]);
      hfpainp_.m15i_[i] = fcx(m15i[i]);
      hfpainp_.m24i_[i] = fcx(m24i[i]);
      painp_.sgi_[i]    = fcx(sgi[i]);
      painp_.p3i_[i]    = fcx(p3i[i]);
      painp_.p8i_[i]    = fcx(p8i[i]);
      hfpainp_.p15i_[i] = fcx(p15i[i]);
      hfpainp_.p24i_[i] = fcx(p24i[i]);
    }
}


void pegasus::evolve()
{
  if (opts.evolmode == 3 || opts.evolmode == 4)
    retrieve();

  //--> this part is not working, should avoid evolmode 1 at LL, could automatically switch to evolmode 0
  if (opts.evolmode == 1 && opts.order_evol == 0)
    {
      cout << "at LL should use evolmode = 0" << endl;
      //N flavour dependence
      int nf = resconst::NF;

      complex <double> fx[11];
      //At LL there is no PDF evolution, PDFs are evaluated at the factorisation scale
      for (int i = 0; i < mellinint::mdim; i++)
	if (opts.mellin1d)
	  pdfevol::retrievemuf_1d(i);
	else
	  pdfevol::retrievemuf_2d(i);

      return;
    }

  //The values of alphas: ASF for the evolution are obtained
  //using the Pegasus alphas function. In order to use the LHAPDF alphas function
  //the code which calculates the evolution at the thresholds should be changed

  //Values of a_s and N_f for the call of EVNFFN or EVNVFN below

  double M2 = real(pow(pdfevol::mubstar,2));   //qbstar = b0/bstar (without a_param) --> The factorisation scale of the PDFs should always be b0/b(star), and not scaled by a_param

  //check to reproduce bb->Z peak as in the old plots
  //M2 = pow(fabs(pdfevol::mubstar_a),2); //bstarscale = a*b0/bstar
  
  //check PDFs at the starting scale
  //  double M2 = asinp_.m20_*10;
  //  pdfevol::bstarscale = sqrt(M2);

  //R2 is the final scale of the evolution, corrected by the ratio of muren2/mufac2
  double R2  = M2 * exp(-frrat_.logfr_);

  //Reproduce the L -> L~ redefinition by a redefinition of the b0/b scale (see Eq. 17 of hep-ph/0508068)
  //With the L -> L~ modification, alphas reaches asymptotically
  //alphas(mures/mufac/muren?) when bscale -> inf, so that the PDF evolution is frozen at muf as upper scale.
  //This modification is required to restore the fixed order cross section upon qt-integration

  double M2tilde = real(pow(pdfevol::mubstartilde,2));
  //  double M2tilde = M2 * resint::mures2 / (M2 + resint::mures2);
  //  double M2tilde = M2 * resint::mufac2 / (M2 + resint::mufac2);
  //  double M2tilde = M2 * resint::muren2 / (M2 + resint::muren2);

  double R2tilde  = M2tilde * exp(-frrat_.logfr_);

  //use the modification L -> L~
  M2 = M2tilde;
  R2 = R2tilde;

  //Further modification of the final scale to account for the fact that alphas(muren)
  //instead of alphas(mures) is used as starting alphas(mu0) in alphasl(nq2) (variable aass)
  //This procedure should hold also for evolmode = 3 (VFN) as far as muf=mur? -> not really working, used only in evolmode = 1
  //double M2prime = M2tilde * resint::muren2/resint::mures2; //scale for differences between muren and mures
  //double R2prime = M2prime * exp(-frrat_.logfr_);

  /*
  //use the modification alphas(mures) -> alphas(muren)
  M2 = M2prime;
  R2 = R2prime;
  */

  //  R2 = R2 * resint::mures2 / resint::muren2;
    
  double asi;
  complex <double> asf;
  int NF;

  
  //DYRES evolution: fixed number of flavours = 5 evolution from the factorisation scale downward
  if (opts.evolmode == 1)
    {
      NF  = nff;
      //set the factorisation scale as initial scale for the evolution 
      //asi = asinp_.as0_;

      //set the resummation scale as initial scale for the evolution 
      asi = resint::alpqres;
      //if (opts.mufevol)
      //asi = resint::alpqfac;

      //to set muf as initial scale, ASF should be multiplied by an additional alphas(muf)/alphas(mures) factor (i.e. change alpqres to alpqfac also in ASF)
      //-> consider to do it so that ASF is a valid final scale for the VFN evolution
      
      //set the factorisation scale as initial scale for the evolution 
      //asi = resint::alpqfac;

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
      asf = as_(R2, R20, AS0, NF)/resint::alpqren * resint::alpqres;
      */
  
      //calculate alphas at the final scale with DYRES (which has the modification L -> L~)
      //DYRES also has a different solution for the running of alphas beyond 1-loop

      //Fortran
      //fcomplex scale2 = fcx(pow(pdfevol::bscale,2));
      //asf = fabs(cx(alphasl_(scale2))) * resint::alpqres; //DYRES LL(NLL) running of alphas

      //C++
      complex <double> b = resconst::b0/pdfevol::bcomplex;
      /*
      if (!opts.mufevol)
	asf = fabs(pdfevol::asl)*resint::alpqres;
      else
	asf = fabs(pdfevol::asl)*resint::alpqfac;
      */

      //asf = fabs(pdfevol::asl)*resint::alpqres;
      asf = pdfevol::asl*resint::alpqres;

      //---> Pegasus and DYRES running of alphas are identical only at LO/LL. At NLO/NLL DYRES evolution is different (it is not a Runge-Kutta solution),
      //and has dependences on the resummation and renormalisation scales
      //evolutionmode = 1 allows to check the impact of this difference.
      //To use the DYRES running of alphas in the VFN evolution set evolmode = 4

      //      cout << pdfevol::bscale << " dyres " << fabs(cx(alphasl_(scale2))) * resint::alpqres << " pegasus " << as_(R2, R20, AS0, NF)/resint::alpqren * resint::alpqres << endl;
      
      //      cout << pdfevol::bscale << "  " << resint::_m << "  " <<  as_(R2, R20, asinp_.as0_, NF) << "  " << fabs(cx(alphasl_(scale2))) * resint::alpqres << "  " << LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI) << endl;
      //      cout << pdfevol::qbstar << "  " << pdfevol::SALP << "  " << log(asi/asf) << endl;
      //      cout << "pegasus blog: " << pdfevol::bscale << "  " << resconst::beta0*resint::alpqres*4.*log(resint::mures2/M2tilde) << endl;
    }

  //Normal Pegasus evolution: VFNS evolution from the starting scale upward
  else if (opts.evolmode == 3)
    {
      //-->force opts.mufevol = true
      
      //use the modification alphas(mures) -> alphas(muren)
      //M2 = M2prime;
      //R2 = R2prime;

      //At Q0 lhapdf and pegasus alphas are equal by construction
      //cout << pdf::alphas(sqrt(asinp_.m20_))/(4.* M_PI) << "  " << alphas(asinp_.m20_, asinp_.m20_, asi, NF) << endl;

      //asf = alphas(pdfevol::mubstartilde, asi, NF);  // --> evolve alphas with Pegasus Runge-Kutta (allow complex scales)
      //asf = pdf::alphas(sqrt(M2))/(4.* M_PI);        // --> evolve alphas with the LHAPDF interpolation

      if  (M2 >= asinp_.m20_)
      	asf = pdf::alphas(sqrt(M2))/(4.* M_PI);// --> above Q0 evolve alphas with the LHAPDF interpolation
      else
	{
	  //asf = alphas(pdfevol::mubstartilde, asi, NF);  // --> below Q0 evolve alphas with Pegasus Runge-Kutta (allow complex scales)
	  // --> below Q0 use the iterative solution for alphas
	  double Q;
	  if (asinp_.m20_ > asfthr_.m2b_)
	    {
	      NF = 5;
	      asi = asfthr_.asb_;
	      Q = sqrt(asfthr_.m2b_);
	    }	      
	  else if (asinp_.m20_ > (asfthr_.m2c_-1e-3))
	    {
	      NF = 4;
	      asi = asfthr_.asc_;
	      Q = sqrt(asfthr_.m2c_);
	    }
	  else
	    {
	      NF = 3;
	      asi = asinp_.as0_;
	      Q = sqrt(asinp_.m20_);
	    }
	  //double QR = scales::ren/scales::res;
	  //double LQR = log(pow(QR,2));
	  //Q *= QR;
	  complex <double> blog = log(pow(Q/pdfevol::mubstartilde,2));
	  double as = asi;
	  double as2 = as*as;
	  complex <double> xlambda = alphas::bet0[NF-NFMIN]*as*blog;
	  complex <double> log1xlambda = log(1.-xlambda);
	  complex <double> logas = 0.;
	  //opts.order should be opts.order_evol?
	  if (opts.order >= 1)
	    logas += log1xlambda;
	  if (opts.order >= 2)
	    logas += as* alphas::bet1[NF-NFMIN]/alphas::bet0[NF-NFMIN]*log1xlambda/(1.-xlambda);
	  if (opts.order >= 3)
	    logas += as2* ((pow(alphas::bet1[NF-NFMIN]/alphas::bet0[NF-NFMIN],2)-alphas::bet2[NF-NFMIN]/alphas::bet0[NF-NFMIN]) *xlambda/pow(1.-xlambda,2)
			   + pow(alphas::bet1[NF-NFMIN]/alphas::bet0[NF-NFMIN],2)             *log1xlambda/pow(1.-xlambda,2)
			   - pow(alphas::bet1[NF-NFMIN]/alphas::bet0[NF-NFMIN],2)             *pow(log1xlambda,2)/(2.*pow(1.-xlambda,2)));
	  //if (opts.order >= 2)
	  //	logas += as*(LQR)
	  //	  *alphas::bet0[NF-NFMIN]*xlambda/(1.-xlambda);
	  //if (opts.order >= 3)
	  //	logas += as2*(+LQR*alphas::bet1[NF-NFMIN]                   *(xlambda-log1xlambda)/pow(1.-xlambda,2)
	  //		      +LQR*alphas::bet1[NF-NFMIN]                   *xlambda/(1.-xlambda)                      //missing piece
	  //		      +0.5*pow(LQR,2)*pow(alphas::bet0[NF-NFMIN],2) *xlambda*(xlambda-2.)/pow(1.-xlambda,2));  //missing piece
	  asf = exp(-logas)*asi;
	}
      
      //asf = asf/resint::alpqren * resint::alpqfac; //--> this is very crude
      
      //In the forward evolution need to rescale the final alphas ASF to account for the fact that the backward PDF evolution
      //is done from Qres to b0/b starting from PDFs at muf, while the residual evolution
      //from muf to Qres is factorised in the (H?) coefficients
      //asf = asf/resint::alpqres * resint::alpqfac; // --> not needed for mufevol = true

      //cout << "scale " << sqrt(M2) << endl;

      //when alphas is evolved with LHAPDF, recompute NF and asi based on the values of ASF
      if  (M2 >= asinp_.m20_)
	if (real(asf) < asfthr_.ast_)
	  {
	    NF = 6;
	    asi = asfthr_.ast_;
	  }
	else if (real(asf) < asfthr_.asb_)
	  {
	    NF = 5;
	    asi = asfthr_.asb_;
	  }
	else if (real(asf) < asfthr_.asc_ || asinp_.m20_ > (asfthr_.m2c_-1e-3))  //if the starting is above mcharm never switch to 3 flavour (intrinsic charm)
	  {
	    NF = 4;
	    asi = asfthr_.asc_;
	  }
	else if (asinp_.m20_ < (asfthr_.m2c_-1e-3))
	  {
	    NF = 3;
	    asi = asinp_.as0_;
	  }
    }
  
  //VFN forward evolution mode in which ASF is evaluated by DYRES
  else if (opts.evolmode == 4)
    {
      //calculate alphas at the final scale with DYRES (which has the modification L -> L~)
      //DYRES also has a different solution for the running of alphas beyond 1-loop
      fcomplex scale2 = fcx(pow(pdfevol::bscale,2));
      asf = fabs(cx(alphasl_(scale2))) * resint::alpqres; //DYRES LL(NLL) running of alphas
      
      //cout << pdfevol::bscale << " dyres " << fabs(cx(alphasl_(scale2))) << endl;
      //cout << pdfevol::bscale << " dyres " << fabs(cx(alphasl_(scale2))) * resint::alpqfac << " LHAPDF " << resint::alpqfac << endl;

      //Rescale the final alphas ASF for the fact that the PDF evolution is done from Qres to b0/b
      //starting from PDFs at muf, while the residual evolution from muf to Qres is factorised in the (H?) coefficients
      asf = asf/resint::alpqres * resint::alpqfac; 
      if (real(asf) < asfthr_.ast_) //If M2 = M2tilde than the scale is frozen at muf, and never goes above the top
	{
	  NF = 6;
	  asi = asfthr_.ast_;
	}
      else if (real(asf) < asfthr_.asb_)
	{
	  NF = 5;
	  asi = asfthr_.asb_;
	}
      else if (real(asf) < asfthr_.asc_)
	{
	  NF = 4;
	  asi = asfthr_.asc_;
	}
      else
	{
	  NF = 3;
	  asi = asinp_.as0_;
	}
    }
  
  //Calculation of the moments of the parton densities and output
  int nlow = 1;
  int nhigh = 2*dim;
  int IPSTD = 1; //IPSTD = 1 return PDFN as valence sea, IPSTD = = return PDFN as q, qb (as in LHAPDF but with u and d inverted)
  fcomplex PDFN[13][ndim];

  //Do not allow the PDFs to evolve below q0 --> freeze the PDF evolution at q0
  //cout << " alphas " << ASF << " qbstar " << pdfevol::qbstar << endl;
  //if (real(ASF) > real(ASI) || fabs(pdfevol::qbstar) < sqrt(asinp_.m20_))
  //ASF = ASI;

  fcomplex ASI = fcx(asi);
  fcomplex ASF = fcx(asf);

  //  double ASIR = asi;
  //  double ASFR = real(asf);

  if (ivfns == 0)
    dyevnffn_(PDFN, ASI, ASF, NF, nlow, nhigh, IPSTD); //modified ffn evolution with charm and bottom at the starting scale
  else
    evnvfn_(PDFN, ASI, ASF, NF, nlow, nhigh, IPSTD);

  //Evolve from Q to muF, to compensate the LQF*gamma terms in the H coefficients
  //if (false)
  if (opts.evolmode == 3 && (opts.kmufac != opts.kmures))
    {
      NF = 5;
      /*
      //Fixed flavour evolution from Q to muF
      ASI = resint::alpqres;
      ASF = resint::alpqfac;
      //order_.npord_ = pdf::order; //order of evolution read from LHAPDF
      //evmod_.imodev_ = 1; //reproduces evolution in x-space
      order_.npord_ = opts.order_evol-1;
      evmod_.imodev_ = 4;
      */
      //cout << "alpqres " << resint::alpqres << endl;
      //cout << "alpqfac " << resint::alpqfac << endl;
      //cout << " ratio " << resint::alpqfac/resint::alpqres << endl;
      
      
      //Fixed flavour evolution from Q to muF
      order_.npord_ = opts.order_evol-1;
      //evmod_.imodev_ = 1; //reproduces evolution in x-space
      double blog = log(pow(scales::res/scales::fac,2));
      double as = resint::aass;
      double LQR = resint::LR-resint::LQ;
      double as2 = as*as;
      double xlambda = beta0*as*blog;
      double log1xlambda = log(1.-xlambda);
      double logas = 0.;
      //opts.order should be opts.order_evol?
      if (opts.order >= 1)
	logas += log1xlambda;
      if (opts.order >= 2)
	logas += as* beta1/beta0*log1xlambda/(1.-xlambda);
      if (opts.order >= 3)
	logas += as2* ((pow(beta1/beta0,2)-beta2/beta0) *xlambda/pow(1.-xlambda,2)
		       + pow(beta1/beta0,2)             *log1xlambda/pow(1.-xlambda,2)
		       - pow(beta1/beta0,2)             *pow(log1xlambda,2)/(2.*pow(1.-xlambda,2)));
      if (opts.order >= 2)
	logas += as*(LQR)
	  *beta0*xlambda/(1.-xlambda);
      if (opts.order >= 3)
	logas += as2*(+LQR*beta1                   *(xlambda-log1xlambda)/pow(1.-xlambda,2)
		       +LQR*beta1                   *xlambda/(1.-xlambda)                      //missing piece
		       +0.5*pow(LQR,2)*pow(beta0,2) *xlambda*(xlambda-2.)/pow(1.-xlambda,2));  //missing piece
      ASI = fcx(resint::alpqres);
      ASF = fcx(exp(-logas)*resint::alpqres);
      
      //double FR2 = pow(scales::ren/scales::fac,2);//ratio of muren2/mufac2
      //frrat_.logfr_ = log(FR2);
      //double R2  = M2 * exp(-frrat_.logfr_);

      //Need to allocate() store() and retrieve() moments at the starting scale,
      //because pegasus::update() is called at most once per each phase space point,
      //while pegasus::evolve is called at each value of b inside the Bessel inversion      
      //allocate();
      //store();
      for (int i = 0; i < 2*dim; i++)
	{
	  //gluon
	  painp_.gli_[i] = PDFN[6][i];
      
	  complex <double> qp[5];
	  qp[0] = (cx(PDFN[6+1][i]) + cx(PDFN[6-1][i]));
	  qp[1] = (cx(PDFN[6+2][i]) + cx(PDFN[6-2][i]));
	  qp[2] = (cx(PDFN[6+3][i]) + cx(PDFN[6-3][i]));
	  qp[3] = (cx(PDFN[6+4][i]) + cx(PDFN[6-4][i]));
	  qp[4] = (cx(PDFN[6+5][i]) + cx(PDFN[6-5][i]));

	  complex <double> qm[5];
	  qm[0] = (cx(PDFN[6+1][i]) - cx(PDFN[6-1][i]));
	  qm[1] = (cx(PDFN[6+2][i]) - cx(PDFN[6-2][i]));
	  qm[2] = (cx(PDFN[6+3][i]) - cx(PDFN[6-3][i]));
	  qm[3] = (cx(PDFN[6+4][i]) - cx(PDFN[6-4][i]));
	  qm[4] = (cx(PDFN[6+5][i]) - cx(PDFN[6-5][i]));

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

      //cout << endl;
      //cout << endl;
      //for (int i = 0; i < mellinint::mdim; i++)
      //	for (int p = -MAXNF; p <= MAXNF; p++)
      //	  cout << i << "  " << p << "  " << cx(PDFN[p+6][i]) << endl;
      //      fcomplex PDFN_evQF[13][ndim];
      dyevnffn_(PDFN, ASI, ASF, NF, nlow, nhigh, IPSTD); //modified ffn evolution with charm and bottom at the starting scale
      //      cout << endl;
      //      for (int i = 0; i < mellinint::mdim; i++)
      //	for (int p = -MAXNF; p <= MAXNF; p++)
      //	  cout << i << "  " << p << "  " << cx(PDFN[p+6][i]) << endl;
      //retrieve();
      //free();
      order_.npord_ = pdf::order; //order of evolution from LHAPDF
      //evmod_.imodev_ = 1;         //reproduces evolution in x-space
      //frrat_.logfr_ = 0.;
    }
  

  complex <double> fx[11];
  for (int sign = mesq::positive; sign <= mesq::negative; sign++)
    for (int n = 0; n < dim; n++)
    {
      int i = index(n,sign);
      for (int p = -MAXNF; p <= MAXNF; p++)
	fx[p+MAXNF] = cx(PDFN[p+6][i]);
      if (opts.mellin1d)
	pdfevol::storemoments(n, sign, fx);
      else
	if (n < mellinint::mdim)
	  pdfevol::storemoments_1(n, sign, fx);
	else
	  pdfevol::storemoments_2(n-mellinint::mdim, sign, fx);
    }

  /*
  cout << endl;
  for (int i = 0; i < dim; i++)
    {
      //int i = 0;
      cout << "pegasus " << sqrt(M2) << "  "
  	   << i << "  " << cx(moms_.na_[i]) << "  "
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

  
  /*
  //Compare Pegasus QCD evolution with direct Mellin transform
  //  if ( real(pdfevol::mubstartilde)  > 2  && real(pdfevol::mubstartilde)  < 4)
    //if ( real(pdfevol::mubstartilde)  > 10)
    {
      evolnum::calculate();

  //for (int i = 0; i < mellinint::mdim; i++)
  //{
      int i = 0;
      //cout << "direct " << fabs(pdfevol::bstarscale) << "  "
      cout << "direct " << real(pdfevol::mubstartilde) << "  "
  	   << i << "  " << mellinint::Np[i] << "  "
  	   << pdfevol::fx1[i*11+(0+MAXNF )] << "  "
  	   << pdfevol::fx1[i*11+(1+MAXNF )]-pdfevol::fx1[i*11+(-1+MAXNF)] << "  "
	   << pdfevol::fx1[i*11+(2+MAXNF )]-pdfevol::fx1[i*11+(-2+MAXNF)] << "  "
  	   << pdfevol::fx1[i*11+(-1+MAXNF)] << "  "
  	   << pdfevol::fx1[i*11+(-2+MAXNF)] << "  "
  	   << pdfevol::fx1[i*11+(3+MAXNF )] << "  "
	   << pdfevol::fx1[i*11+(4+MAXNF )] << "  "
	   << pdfevol::fx1[i*11+(5+MAXNF )] << endl;
  
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

//interface to the VFN alphas of pegasus
//double pegasus::alphas(double M2, double R2, double &ASI, int &NF)
complex <double> pegasus::alphas(complex <double> m, double &ASI, int &NF)
{
  complex <double> ASF = 0.;
  //double QR2 = pow(scales::res/scales::ren,2);
  double QR2 = 1.;
  complex <double> r = m;//*scales::res/scales::ren;

  double M2 = pow(real(m),2);
  double R2 = pow(real(r),2);
  
  if (M2 > asfthr_.m2t_) //If M2 = M2tilde than the scale is frozen at muf, and never goes above the top
    {
      NF = 6;
      double R2T = asfthr_.m2t_ * QR2;//R2/M2;;
      ASI = asfthr_.ast_;
      ASF = as_(R2, R2T, asfthr_.ast_, NF);
    }
  else if (M2 > asfthr_.m2b_)
    {
      NF = 5;
      double R2B = asfthr_.m2b_ * QR2;//R2/M2;;
      ASI = asfthr_.asb_;
      ASF =  as_(R2, R2B, asfthr_.asb_, NF);
    }
  else if (M2 > asfthr_.m2c_ || asinp_.m20_ > (asfthr_.m2c_-1e-3))   //if the starting scale is above mcharm never switch to 3 flavour (intrinsic charm)
    {
      NF = 4;
      double R2C = asfthr_.m2c_ * QR2;//R2/M2;;
      ASI = asfthr_.asc_;
      //ASF =  as_(R2, R2C, asfthr_.asc_, NF);

      //alphas::rgkt(r, sqrt(R2C), asfthr_.asc_*4*M_PI, NF);
      alphas::iter(r, sqrt(R2C), asfthr_.asc_*4*M_PI, NF);
      if (aspar_.naord_ == 0)
	ASF = alphas::asLO/4./M_PI;
      else if (aspar_.naord_ == 1)
	ASF = alphas::asNLO/4./M_PI;
      else if (aspar_.naord_ == 2)
	ASF = alphas::asNNLO/4./M_PI;
    }
  else
    {
      NF = 3;
      double R20 = asinp_.m20_ * QR2;//R2/M2;;
      ASI = asinp_.as0_;
      
      //ASF =  as_(R2, R20, asinp_.as0_, NF);

      //alphas::rgkt(r, sqrt(R20), asinp_.as0_*4*M_PI, NF);
      alphas::iter(r, sqrt(R20), asinp_.as0_*4*M_PI, NF);
      if (aspar_.naord_ == 0)
	ASF = alphas::asLO/4./M_PI;
      else if (aspar_.naord_ == 1)
	ASF = alphas::asNLO/4./M_PI;
      else if (aspar_.naord_ == 2)
	ASF = alphas::asNNLO/4./M_PI;
    }

  ////QCD coupling scale dependence
  //double LQR = -log(QR2); //log(pow(scales::ren/scales::res,2));
  //if (aspar_.naord_ >= 1)  
  //  ASF += 1/4./M_PI*(-pow(ASI*4*M_PI,2)*alphas::bet0[NF-NFMIN]/M_PI*LQR);
  //if (aspar_.naord_ >= 2)  
  //  ASF += 1/4./M_PI*(-pow(ASI*4*M_PI,3)*LQR*(alphas::bet1[NF-NFMIN] - pow(alphas::bet0[NF-NFMIN],2)*LQR)/pi2);
  //if (aspar_.naord_ >= 3)
  //  ASF += 1/4./M_PI*(-pow(ASI*4*M_PI,4)*LQR*(alphas::bet2[NF-NFMIN] - 5./2.*alphas::bet0[NF-NFMIN]*alphas::bet1[NF-NFMIN]*LQR + pow(alphas::bet0[NF-NFMIN],3)*pow(LQR,2))/(pi2*M_PI));
  //if (aspar_.naord_ >= 4)
  //  ASF += 1/4./M_PI*(-pow(ASI*4*M_PI,5)*LQR*(alphas::bet3[NF-NFMIN] - 3./2.*pow(alphas::bet1[NF-NFMIN],2)*LQR - 3.*alphas::bet0[NF-NFMIN]*alphas::bet2[NF-NFMIN]*LQR + 13./3.*pow(alphas::bet0[NF-NFMIN],2)*alphas::bet1[NF-NFMIN]*pow(LQR,2) - pow(alphas::bet0[NF-NFMIN],4)*pow(LQR,3))/pi4);

  //cout << "m " << m << " as " << ASF << endl;
  
  //above the bottom mass, calculate alphas at the final scale with DYRES (which has the modification L -> L~)
  //	  fcomplex fscale2 = fcx(pow(pdfevol::bscale,2));
  //	  ASF = fabs(cx(alphasl_(fscale2))) * resint::alpqres; //DYRES LL(NLL) running of alphas.

  //ASF = LHAPDF::alphasPDF(sqrt(M2)) / (4.* M_PI);
  return ASF;
}

