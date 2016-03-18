#include "mellinint.h"
#include "mesq.h"
#include "settings.h"
#include "gaussrules.h"
#include "hcoefficients.h"
#include "pdfevol.h"
#include "interface.h"
#include "anomalous.h"

#include <complex>
#include <iostream>
#include <iomanip>

double *mellinint::wn;
complex <double> *mellinint::Np;
complex <double> *mellinint::Nm;
complex <double> mellinint::CCp;
complex <double> mellinint::CCm;

complex <double> mellinint::GGN;
complex <double> mellinint::QGN_1;
complex <double> mellinint::QGN_2;
complex <double> mellinint::QQBN_1;
complex <double> mellinint::QQBN_2;
complex <double> mellinint::QQBN_3;
complex <double> mellinint::QQBN_4;

int mellinint::mdim;

//fortran interfaces
void mellinint_pdf_mesq_expy_(int& i1, int& i2, int& sign)
{
  mellinint::pdf_mesq_expy(i1-1, i2-1, sign-1);
};
fcomplex mellinint_integrand_(int& i1, int& i2, int& sign)
{
  return fcx(mellinint::integrand(i1-1, i2-1, sign-1));
};

void mellinint::initgauss()
{
  //set up weights and nodes for gaussian quadrature
  //cpoint is the starting point on the real axis for the positive and negative part of the integration path
  //phi is the angle in the complex plane of the positive part of the integration path
       
  double cpoint = 1.;
  double phi = M_PI * 1./2.;
      
  double min = 0;
  double max = 27.; //upper limit for the mellin integration in the complex plane (z). Above 50 the integral becomes unstable (blim issue with Landau pole?)
      
  mdim = opts.mellinintervals*opts.mellinrule;

  //allocate memory
  wn = new double [mdim];
  Np = new complex <double> [mdim];
  Nm = new complex <double> [mdim];
  
  //initialise to 0
  for (int i=0; i < mdim; i++)
    {
      Np[i]=cpoint+1.;
      Nm[i]=cpoint+1.;
      wn[i]=0;
    }

  // positive branch      
  CCp = complex <double> (cos(phi), sin(phi));
  for (int i=0; i < opts.mellinintervals; i++)
    {
      double a = 0.+(1.-0.)*i/opts.mellinintervals;
      double b = 0.+(1.-0.)*(i+1)/opts.mellinintervals;
      double c = 0.5*(a+b);
      double m = 0.5*(b-a);
      for (int j=0; j < opts.mellinrule; j++)
	{
	  double x = c+m*gr::xxx[opts.mellinrule-1][j];
	  double t = min+(max-min)*x;
	  double jac = max-min;
	  Np[j+i*opts.mellinrule]=complex <double> (cpoint+cos(phi)*t+1.,sin(phi)*t);
	  wn[j+i*opts.mellinrule]=gr::www[opts.mellinrule-1][j]*m*jac;
	  //	  cout << setprecision(16) <<  t << " " << Np[j+i*opts.mellinrule] << "  " << wn[j+i*opts.mellinrule] << endl;
	}
    }      
  // negative branch
  CCm = complex <double> (cos(phi), -sin(phi));
  for (int i=0; i < opts.mellinintervals; i++)
    {
      double a = 0.+(1.-0.)*i/opts.mellinintervals;
      double b = 0.+(1.-0.)*(i+1)/opts.mellinintervals;
      double c = 0.5*(a+b);
      double m = 0.5*(b-a);
      for (int j=0; j < opts.mellinrule; j ++)
	{
	  double x = c+m*gr::xxx[opts.mellinrule-1][j];
	  double t = min+(max-min)*x;
	  Nm[j+i*opts.mellinrule]=complex <double> (cpoint+cos(phi)*t+1.,-sin(phi)*t);
	}
    }
}


//input : fn1 fn2, sigmaij
//output: ggn, qgn qqbn

// l'integrando della trasformata di Mellin e' composto di vari pezzi:
// A) le luminosita' partoniche nello spazio di mellin: fn1 * fn2 -> GGN, QGN, QQBN  funzione di (i,j,I1,I2)
// B) le ampiezze EW sigmaij, funzione di m, costh (e se integrate in costh, i momenti di costh sono funzione di y, pt, m in presenza di tagli sui leptoni) funzione di (i,j)
// C) la dipendenza esplicita' della rapidita' cex1 cex2 funzione di (I1, I2)
// D) i coefficienti di Wilson Hgg Hqqb, etc.. funzione di (I1,I2)

//L'integrazione in costh phi entra solamente in B
//L'integrazione in costh phi e rapidita' entra in B e C, per cui B e C devono essere calcolati insieme
//Il prodotto con D si puo' fare alla fine separatemente

//This function performs the product of PDF, born level amplitudes and expy piece, and sums over partonic channels (i,j)
//It is a function of z1 and z2 in Mellin space
void mellinint::pdf_mesq_expy(int i1, int i2, int sign)
{
  GGN=0;
  QGN_1=0;
  QGN_2=0;
  QQBN_1=0;
  QQBN_2=0;
  QQBN_3=0;
  QQBN_4=0;

  complex<double>* fn1 = pdfevol::fn1;
  complex<double>* fn2 = pdfevol::fn2;

  //DYRES convention
  // bb cb sb db ub g u d s c b
  // -5 -4 -3 -2 -1 0 1 2 3 4 5
    
  //mesqij[n] -> mesqij_expy[mesq::index(n,i1,i2,sign)]
  if (opts.nproc == 3)
    {
      complex <double> mesq_uub = mesq::mesqij_expy[mesq::index(0,i1,i2,sign)];
      complex <double> mesq_ddb = mesq::mesqij_expy[mesq::index(1,i1,i2,sign)];
      complex <double> mesq_ubu = mesq::mesqij_expy[mesq::index(2,i1,i2,sign)];
      complex <double> mesq_dbd = mesq::mesqij_expy[mesq::index(3,i1,i2,sign)];

      //NLL part
      QGN_1 = fn1[g]*((fn2[ub]+fn2[cb])*mesq_uub
		      + (fn2[db]+fn2[sb]+fn2[bb])*mesq_ubu
		      + (fn2[u]+fn2[c])*mesq_ddb
		      + (fn2[d]+fn2[s]+fn2[b])*mesq_dbd);

      QGN_2 = fn2[g]*((fn1[u]+fn1[c])*mesq_uub
		      + (fn1[d]+fn1[s]+fn1[b])*mesq_ubu
		      + (fn1[ub]+fn1[cb])*mesq_ddb
		      + (fn1[db]+fn1[sb]+fn1[bb])*mesq_dbd);

      QQBN_1 = (fn1[u]*fn2[ub]+fn1[c]*fn2[cb])*mesq_uub
	+ (fn1[d]*fn2[db]+fn1[s]*fn2[sb]+fn1[b]*fn2[bb])*mesq_ubu
	+ (fn1[ub]*fn2[u]+fn1[cb]*fn2[c])*mesq_ddb
	+ (fn1[db]*fn2[d]+fn1[sb]*fn2[s]+fn1[bb]*fn2[b])*mesq_dbd;

      //NNLL
      if(opts.order >= 2)
	{
	  GGN = fn1[g]*fn2[g]*
	    (2.*mesq_uub
	     + 3.*mesq_ubu
	     + 2.*mesq_ddb
	     + 3.*mesq_dbd);
            
	  QQBN_2 = (fn1[u]*fn2[u]+fn1[c]*fn2[c])*mesq_uub
	    + (fn1[d]*fn2[d]+fn1[s]*fn2[s]+fn1[b]*fn2[b])*mesq_ubu
	    + (fn1[ub]*fn2[ub]+fn1[cb]*fn2[cb])*mesq_ddb
	    + (fn1[db]*fn2[db]+fn1[sb]*fn2[sb]+fn1[bb]*fn2[bb])*mesq_dbd;

	  QQBN_3 = fn1[u]*(fn2[db]*mesq_ubu
			   +fn2[sb]*mesq_ubu
			   +fn2[cb]*mesq_uub
			   +fn2[bb]*mesq_ubu
			   +fn2[d]*mesq_dbd
			   +fn2[s]*mesq_dbd
			   +fn2[c]*mesq_ddb
			   +fn2[b]*mesq_dbd)
	    + fn1[d]*(fn2[ub]*mesq_uub
		      +fn2[sb]*mesq_ubu
		      +fn2[cb]*mesq_uub
		      +fn2[bb]*mesq_ubu
		      +fn2[u]*mesq_ddb
		      +fn2[s]*mesq_dbd
		      +fn2[c]*mesq_ddb
		      +fn2[b]*mesq_dbd)
	    + fn1[s]*(fn2[ub]*mesq_uub
		      +fn2[db]*mesq_ubu
		      +fn2[cb]*mesq_uub
		      +fn2[bb]*mesq_ubu
		      +fn2[u]*mesq_ddb
		      +fn2[d]*mesq_dbd
		      +fn2[c]*mesq_ddb
		      +fn2[b]*mesq_dbd)
	    + fn1[c]*(fn2[ub]*mesq_uub
		      +fn2[db]*mesq_ubu
		      +fn2[sb]*mesq_ubu
		      +fn2[bb]*mesq_ubu
		      +fn2[u]*mesq_ddb
		      +fn2[d]*mesq_dbd
		      +fn2[s]*mesq_dbd
		      +fn2[b]*mesq_dbd)
	    + fn1[b]*(fn2[ub]*mesq_uub
		      +fn2[db]*mesq_ubu
		      +fn2[sb]*mesq_ubu
		      +fn2[cb]*mesq_uub
		      +fn2[u]*mesq_ddb
		      +fn2[d]*mesq_dbd
		      +fn2[s]*mesq_dbd
		      +fn2[c]*mesq_ddb)
	    + fn1[ub]*(fn2[db]*mesq_ubu
		       +fn2[sb]*mesq_ubu
		       +fn2[cb]*mesq_uub
		       +fn2[bb]*mesq_ubu
		       +fn2[d]*mesq_dbd
		       +fn2[s]*mesq_dbd
		       +fn2[c]*mesq_ddb
		       +fn2[b]*mesq_dbd)
	    + fn1[db]*(fn2[ub]*mesq_uub
		       +fn2[sb]*mesq_ubu
		       +fn2[cb]*mesq_uub
		       +fn2[bb]*mesq_ubu
		       +fn2[u]*mesq_ddb
		       +fn2[s]*mesq_dbd
		       +fn2[c]*mesq_ddb
		       +fn2[b]*mesq_dbd)
	    + fn1[sb]*(fn2[ub]*mesq_uub
		       +fn2[db]*mesq_ubu
		       +fn2[cb]*mesq_uub
		       +fn2[bb]*mesq_ubu
		       +fn2[u]*mesq_ddb
		       +fn2[d]*mesq_dbd
		       +fn2[c]*mesq_ddb
		       +fn2[b]*mesq_dbd)
	    + fn1[cb]*(fn2[ub]*mesq_uub
		       +fn2[db]*mesq_ubu
		       +fn2[sb]*mesq_ubu
		       +fn2[bb]*mesq_ubu
		       +fn2[u]*mesq_ddb
		       +fn2[d]*mesq_dbd
		       +fn2[s]*mesq_dbd
		       +fn2[b]*mesq_dbd)
	    + fn1[bb]*(fn2[ub]*mesq_uub
		       +fn2[db]*mesq_ubu
		       +fn2[sb]*mesq_ubu
		       +fn2[cb]*mesq_uub
		       +fn2[u]*mesq_ddb
		       +fn2[d]*mesq_dbd
		       +fn2[s]*mesq_dbd
		       +fn2[c]*mesq_ddb);

	  QQBN_4 = fn2[u]*(fn1[d]*mesq_ubu
			   +fn1[s]*mesq_ubu
			   +fn1[c]*mesq_uub
			   +fn1[b]*mesq_ubu
			   +fn1[db]*mesq_dbd
			   +fn1[sb]*mesq_dbd
			   +fn1[cb]*mesq_ddb
			   +fn1[bb]*mesq_dbd)
	    + fn2[d]* (fn1[u]*mesq_uub
		       +fn1[s]*mesq_ubu
		       +fn1[c]*mesq_uub
		       +fn1[b]*mesq_ubu
		       +fn1[ub]*mesq_ddb
		       +fn1[sb]*mesq_dbd
		       +fn1[cb]*mesq_ddb
		       +fn1[bb]*mesq_dbd)
	    + fn2[s]* (fn1[u]*mesq_uub
		       +fn1[d]*mesq_ubu
		       +fn1[c]*mesq_uub
		       +fn1[b]*mesq_ubu
		       +fn1[ub]*mesq_ddb
		       +fn1[db]*mesq_dbd
		       +fn1[cb]*mesq_ddb
		       +fn1[bb]*mesq_dbd)
	    + fn2[c]* (fn1[u]*mesq_uub
		       +fn1[d]*mesq_ubu
		       +fn1[s]*mesq_ubu
		       +fn1[b]*mesq_ubu
		       +fn1[ub]*mesq_ddb
		       +fn1[db]*mesq_dbd
		       +fn1[sb]*mesq_dbd
		       +fn1[bb]*mesq_dbd)
	    + fn2[b]* (fn1[u]*mesq_uub
		       +fn1[d]*mesq_ubu
		       +fn1[s]*mesq_ubu
		       +fn1[c]*mesq_uub
		       +fn1[ub]*mesq_ddb
		       +fn1[db]*mesq_dbd
		       +fn1[sb]*mesq_dbd
		       +fn1[cb]*mesq_ddb)
	    + fn2[ub]*(fn1[d]*mesq_ubu
		       +fn1[s]*mesq_ubu
		       +fn1[c]*mesq_uub
		       +fn1[b]*mesq_ubu
		       +fn1[db]*mesq_dbd
		       +fn1[sb]*mesq_dbd
		       +fn1[cb]*mesq_ddb
		       +fn1[bb]*mesq_dbd)
	    + fn2[db]*(fn1[u]*mesq_uub
		       +fn1[s]*mesq_ubu
		       +fn1[c]*mesq_uub
		       +fn1[b]*mesq_ubu
		       +fn1[ub]*mesq_ddb
		       +fn1[sb]*mesq_dbd
		       +fn1[cb]*mesq_ddb
		       +fn1[bb]*mesq_dbd)
	    + fn2[sb]*(fn1[u]*mesq_uub
		       +fn1[d]*mesq_ubu
		       +fn1[c]*mesq_uub
		       +fn1[b]*mesq_ubu
		       +fn1[ub]*mesq_ddb
		       +fn1[db]*mesq_dbd
		       +fn1[cb]*mesq_ddb
		       +fn1[bb]*mesq_dbd)
	    + fn2[cb]*(fn1[u]*mesq_uub
		       +fn1[d]*mesq_ubu
		       +fn1[s]*mesq_ubu
		       +fn1[b]*mesq_ubu
		       +fn1[ub]*mesq_ddb
		       +fn1[db]*mesq_dbd
		       +fn1[sb]*mesq_dbd
		       +fn1[bb]*mesq_dbd)
	    + fn2[bb]*(fn1[u]*mesq_uub
		       +fn1[d]*mesq_ubu
		       +fn1[s]*mesq_ubu
		       +fn1[c]*mesq_uub
		       +fn1[ub]*mesq_ddb
		       +fn1[db]*mesq_dbd
		       +fn1[sb]*mesq_dbd
		       +fn1[cb]*mesq_ddb);
	}
    }
}

complex <double> mellinint::integrand(int i1, int i2, int sign)
{
  int i = hcoefficients::index(i1,i2,sign);
  //  cout << i1 << "  " << i2 << "  " << GGN << "  " << hcoefficients::Hgg[i] <<  endl;
  return GGN*hcoefficients::Hgg[i] + QGN_1*hcoefficients::Hqg_1[i] + QGN_2*hcoefficients::Hqg_2[i]
    + QQBN_1*hcoefficients::Hqqb[i] + QQBN_2*hcoefficients::Hqq[i] + QQBN_3*hcoefficients::Hqqp_1[i] + QQBN_4*hcoefficients::Hqqp_2[i];
}
