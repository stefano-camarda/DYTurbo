#include "hcoeff_check.h"
#include "resconst.h"
#include "anomalous.h"
#include "ccoeff.h"
#include "mesq.h"
#include "besselint.h"
#include "resint.h"
#include "settings.h"
#include "interface.h"
#include "phasespace.h"
#include "parton.h"
#include "pmom.h"
#include <ctime>

#include <iostream>
#include <map>

complex <double> *hcoeff_check::Hqqb;
complex <double> *hcoeff_check::Hqg;
complex <double> *hcoeff_check::Hqq;
complex <double> *hcoeff_check::Hqqp;
complex <double> *hcoeff_check::Hgg;
complex <double> *hcoeff_check::Hqqbp;
complex <double> *hcoeff_check::Hqbg;
complex <double> *hcoeff_check::Hqpg;
complex <double> *hcoeff_check::Hqbpg;

using namespace anomalous;
using namespace resconst;
using namespace ccoeff;
using namespace resint;
using namespace parton;


//Flavour-to-parton mapping class
template <class T> class fpmap
{
private:
  enum part {q, qb, g, qp, qbp};

  //map <part, map <part, complex <double>> > m;
  map <part, map <part, T> > m;

  //parton mapping logic
  //f has 11 flavours, a is the compressed flavour (q,qb,g,qp,qbp)
  //Assume the Z/gamma* case of ab -> qqb with q=u
  void pmap(parton::partid f, part &a)
  {
    if (f == parton::u)
      a = part::q;
    else if (f == parton::ub)
      a = part::qb;
    else if (f == parton::g)
      a = part::g;
    else if (isq(f) == 1)
      a = part::qp;
    else if (isq(f) == -1)
      a = part::qbp;
    return;
  }

  void pmap(parton::partid A, parton::partid B, part &a, part &b)
  {
    if (A == parton::g)
      {
	a = part::g;
	pmap(B,b);
	return;
      }
	
    if (B == parton::g)
      {
	b = part::g;
	pmap(A,a);
	return;
      }

    if (A == B)
      {
	a = part::q;
	b = part::q;
	return;
      }

    if (A == charge_conjn(B))
      {
	a = part::q;
	b = part::qb;
	return;
      }

    if (isq(A)*isq(B) == 1)
      {
	a = part::q;
	b = part::qp;
	return;
      }

    if (isq(A)*isq(B) == -1)
      {
	a = part::q;
	b = part::qbp;
	return;
      }
  }
  
public:
  //complex <double> operator()(parton::partid A, parton::partid B)
  T operator()(parton::partid A, parton::partid B)
  {
    part a, b;
    //pmap(A, a);
    //pmap(B, b);
    pmap(A,B,a,b);
    return m[a][b];
  };

  //complex <double> operator()(int A, parton::partid B)
  T operator()(int A, parton::partid B)
  {
    part a, b;
    //pmap(partid(A), a);
    //pmap(B, b);
    pmap(partid(A),B,a,b);
    return m[a][b]*sqrt(11.);
  };
  

  //complex <double> operator()(parton::partid A, int B)
  T operator()(parton::partid A, int B)
  {
    part a, b;
    //pmap(A, a);
    //pmap(partid(B), b);
    pmap(A,partid(B),a,b);
    return m[a][b]*sqrt(11.);
  };

  //complex <double> operator()(int A, int B)
  T operator()(int A, int B)
  {
    part a, b;
    //pmap(partid(A), a);
    //pmap(partid(B), b);
    pmap(partid(A),partid(B),a,b);
    return m[a][b]*11.;
  };
  

  //complex <double> *get(parton::partid A, parton::partid B)
  T *get(parton::partid A, parton::partid B)
  {
    part a, b;
    pmap(A, a);
    pmap(B, b);
    //pmap(A,B,a,b);
    return &(m[a][b]);
  };
};

void hcoeff_check::allocate()
{
  if (opts.order == 0)
    return;

  //allocate memory
  //LL
  Hqqb = new complex <double> [mellinint::mdim];

  //NLL
  Hqg = new complex <double> [mellinint::mdim];

  //NNLL
  Hqq   = new complex <double> [mellinint::mdim];
  Hqqp  = new complex <double> [mellinint::mdim];
  Hqqbp = new complex <double> [mellinint::mdim];
  Hgg   = new complex <double> [mellinint::mdim];

  //NNNLL
  Hqbg  = new complex <double> [mellinint::mdim];
  Hqpg  = new complex <double> [mellinint::mdim];
  Hqbpg = new complex <double> [mellinint::mdim];
}

void hcoeff_check::reset()
{
  if (opts.order == 0)
    return;

  fill(Hqqb, Hqqb+mellinint::mdim, 0.);
  fill(Hqg , Hqg +mellinint::mdim, 0.);
  fill(Hqq , Hqq +mellinint::mdim, 0.);
  fill(Hqqp, Hqqp+mellinint::mdim, 0.);
  fill(Hqqbp, Hqqbp+mellinint::mdim, 0.);
  fill(Hgg , Hgg +mellinint::mdim, 0.);
  fill(Hqbg , Hqbg +mellinint::mdim, 0.);
  fill(Hqpg , Hqpg +mellinint::mdim, 0.);
  fill(Hqbpg , Hqbpg +mellinint::mdim, 0.);
}

void hcoeff_check::free()
{
  if (opts.order == 0)
    return;

  delete[] Hqqb;
  delete[] Hqg;
  delete[] Hqq;
  delete[] Hqqp;  
  delete[] Hqqbp;
  delete[] Hgg;
  delete[] Hqbg;
  delete[] Hqpg;
  delete[] Hqbpg;
}

void hcoeff_check::calc()
{
  double LQ2 = pow(LQ,2);
  double LQ3 = pow(LQ,3);
  double LQ4 = pow(LQ,4);
  double LQ5 = pow(LQ,5);
  double LQ6 = pow(LQ,6);

  double LR2 = pow(LR,2);

  double LF2 = pow(LF,2);
  double LF3 = pow(LF,3);
  
  for (int i = 0; i < mellinint::mdim; i++)
    {
      int idx = anomalous::index(i,mesq::positive);
      parton::partid q   = u;
      parton::partid qb  = ub;
      parton::partid qp  = d;
      parton::partid qbp = db;
	  
      fpmap <double> del;
      *del.get(q,q)    = 1.;
      *del.get(qb,qb)  = 1.;
      *del.get(g,g)    = 1.;
      *del.get(qp,qp)  = 1.;
      *del.get(qbp,qbp)= 1.;
	  
      fpmap <complex <double>> C1;
      *C1.get(q,q)   = *C1.get(qb,qb) = *C1.get(qp,qp)  = *C1.get(qbp,qbp) = C1qq[idx];
      *C1.get(q,g)   = *C1.get(g,q)   = *C1.get(qb,g)   = *C1.get(g,qb)    = C1qg[idx];
      *C1.get(qp,g)  = *C1.get(g,qp)  = *C1.get(qbp,g)  = *C1.get(g,qbp)   = C1qg[idx];
      *C1.get(q,qb)  = *C1.get(qb,q)  = *C1.get(qp,qbp) = *C1.get(qbp,qp)  = 0.;
      *C1.get(q,qp)  = *C1.get(qp,q)  = *C1.get(qb,qbp) = *C1.get(qbp,qb)  = 0.;
      *C1.get(q,qbp) = *C1.get(qbp,q) = *C1.get(qb,qp)  = *C1.get(qp,qb)   = 0.;
      //*C1.get(g,g) = 2. * resconst::C1ggn;
	  
      fpmap <complex <double>> gamma1;
      *gamma1.get(q,q)   = *gamma1.get(qb,qb) = *gamma1.get(qp,qp) = *gamma1.get(qbp,qbp) = pmom::gamma1qq[idx];
      *gamma1.get(q,g)   = *gamma1.get(qb,g)  = *gamma1.get(qp,g)  = *gamma1.get(qbp,g)   = pmom::gamma1qg[idx];
      *gamma1.get(g,q)   = *gamma1.get(g,qb)  = *gamma1.get(g,qp)  = *gamma1.get(g,qbp)   = pmom::gamma1gq[idx];
      *gamma1.get(q,qb)  = *gamma1.get(qb,q)  = *gamma1.get(qp,qbp) = *gamma1.get(qbp,qp) = 0.;
      *gamma1.get(q,qp)  = *gamma1.get(qp,q)  = *gamma1.get(qb,qbp) = *gamma1.get(qbp,qb) = 0.;
      *gamma1.get(q,qbp) = *gamma1.get(qbp,q) = *gamma1.get(qb,qp)  = *gamma1.get(qp,qb)  = 0.;
      *gamma1.get(g,g)   = pmom::gamma1gg[idx];
      
      fpmap <complex <double>> H1st;

      vector <partid> pa, pb;
      pa.push_back(q); pb.push_back(qb);
      pa.push_back(q); pb.push_back(g);

      for (int p = 0; p < pa.size(); p++)
	{
	  parton::partid a = pa[p];
	  parton::partid b = pb[p];
	  *H1st.get(a,b) = (-(B1q*LQ) - (A1q*pow(LQ,2))/2.)*del(q,a)*del(qb,b) + del(qb,b)*(C1(q,a) + (LF - LQ)*gamma1(q,a)) + del(q,a)*(C1(qb,b) + (LF - LQ)*gamma1(qb,b));
	}

      double as = aass;
      Hqqb[i] += as*H1st(q,qb);
      Hqg[i]  += as*H1st(q,g);

      //cout << endl;
      //cout << "Hst_scale2.txt " << endl;
      //cout << "H1st(q,qb)  " << H1st(q,qb)  << endl;
      //cout << "H1st(q,g)   " << H1st(q,g)   << endl;
      
      if (opts.order == 1)
	continue;
	  
      fpmap <complex <double>> C2;
      *C2.get(q,q)   = *C2.get(qb,qb) = *C2.get(qp,qp)  = *C2.get(qbp,qbp) = C2qq[idx];
      *C2.get(q,g)   = *C2.get(g,q)   = *C2.get(qb,g)   = *C2.get(g,qb)    = C2qg[idx];
      *C2.get(qp,g)  = *C2.get(g,qp)  = *C2.get(qbp,g)  = *C2.get(g,qbp)   = C2qg[idx];
      *C2.get(q,qb)  = *C2.get(qb,q)  = *C2.get(qp,qbp) = *C2.get(qbp,qp)  = C2qqb[idx];
      *C2.get(q,qp)  = *C2.get(qp,q)  = *C2.get(qb,qbp) = *C2.get(qbp,qb)  = C2qqp[idx];
      *C2.get(q,qbp) = *C2.get(qbp,q) = *C2.get(qb,qp)  = *C2.get(qp,qb)   = C2qqbp[idx];

      fpmap <complex <double>> gamma2;
      //*gamma2.get(q,q)   = *gamma2.get(qb,qb) = *gamma2.get(qp,qp)  = *gamma2.get(qbp,qbp) = 1./4.*(gamma2qqV[idx]+gamma2qqS[idx]);
      //*gamma2.get(q,g)   = *gamma2.get(qb,g)  = *gamma2.get(qp,g)   = *gamma2.get(qbp,g)   = 1./4.*gamma2qg[idx];
      //*gamma2.get(g,q)   = *gamma2.get(g,qb)  = *gamma2.get(g,qp)   = *gamma2.get(g,qbp)   = 1./4.*gamma2gq[idx];
      //*gamma2.get(q,qb)  = *gamma2.get(qb,q)  = *gamma2.get(qp,qbp) = *gamma2.get(qbp,qp)  = 1./4.*(gamma2qqbV[idx]+gamma2qqbS[idx]);
      //*gamma2.get(q,qp)  = *gamma2.get(qp,q)  = *gamma2.get(qb,qbp) = *gamma2.get(qbp,qb)  = 1./4.*gamma2qqbS[idx];
      //*gamma2.get(q,qbp) = *gamma2.get(qbp,q) = *gamma2.get(qb,qp)  = *gamma2.get(qp,qb)   = 1./4.*gamma2qqbS[idx];
      //*gamma2.get(g,g)   = 1./4.*gamma2gg[idx];
      
      *gamma2.get(q,q)   = *gamma2.get(qb,qb) = *gamma2.get(qp,qp)  = *gamma2.get(qbp,qbp) = pmom::gamma2qq[idx];
      *gamma2.get(q,g)   = *gamma2.get(qb,g)  = *gamma2.get(qp,g)   = *gamma2.get(qbp,g)   = pmom::gamma2qg[idx];
      *gamma2.get(g,q)   = *gamma2.get(g,qb)  = *gamma2.get(g,qp)   = *gamma2.get(g,qbp)   = pmom::gamma2gq[idx];
      *gamma2.get(q,qb)  = *gamma2.get(qb,q)  = *gamma2.get(qp,qbp) = *gamma2.get(qbp,qp)  = pmom::gamma2qqb[idx];
      *gamma2.get(q,qp)  = *gamma2.get(qp,q)  = *gamma2.get(qb,qbp) = *gamma2.get(qbp,qb)  = pmom::gamma2qqp[idx];
      *gamma2.get(q,qbp) = *gamma2.get(qbp,q) = *gamma2.get(qb,qp)  = *gamma2.get(qp,qb)   = pmom::gamma2qqbp[idx];
      *gamma2.get(g,g)   = pmom::gamma2gg[idx];

      fpmap <complex <double>> H2st;

      pa.push_back(q); pb.push_back(q);
      pa.push_back(q); pb.push_back(qp);
      pa.push_back(q); pb.push_back(qbp);
      pa.push_back(g); pb.push_back(g);

      vector <parton::partid> allf;
      allf.push_back(parton::bb);
      allf.push_back(parton::cb);
      allf.push_back(parton::sb);
      allf.push_back(parton::db);
      allf.push_back(parton::ub);
      allf.push_back(parton::g);
      allf.push_back(parton::u);
      allf.push_back(parton::d);
      allf.push_back(parton::s);
      allf.push_back(parton::c);
      allf.push_back(parton::b);
	  
      for (int p = 0; p < pa.size(); p++)
	{
	  parton::partid a = pa[p];
	  parton::partid b = pb[p];
	  *H2st.get(a,b) = 0;

	  for (int pa1 = 0; pa1 < allf.size(); pa1++)
	    for (int pb1 = 0; pb1 < allf.size(); pb1++)
	      {
		int a1 = allf[pa1];
		int b1 = allf[pb1];

		*H2st.get(a,b) += C1(q,a)*C1(qb,b) + (-(B2q*LQ) + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8. + beta0*(-(B1q*pow(LQ,2))/2. - (A1q*pow(LQ,3))/3.))*del(q,a)*del(qb,b) + (LF - LQ)*C1(qb,b)*gamma1(q,a) + ((LF - LQ)*C1(q,a) + (pow(LF,2) - 2*LF*LQ + pow(LQ,2))*gamma1(q,a))*gamma1(qb,b) - beta0*LR*((-(B1q*LQ) - (A1q*pow(LQ,2))/2.)*del(q,a)*del(qb,b) + del(qb,b)*(C1(q,a) + (LF - LQ)*gamma1(q,a)) + del(q,a)*(C1(qb,b) + (LF - LQ)*gamma1(qb,b))) + del(qb,b)*((-(B1q*LQ) + beta0*LQ - (A1q*pow(LQ,2))/2.)*C1(q,a) + C2(q,a) + (B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2. + beta0*(pow(LF,2)/2. - pow(LQ,2)/2.) + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.))*gamma1(q,a) + (pow(LF,2)*gamma1(a1,a)*gamma1(q,a1))/2. + (pow(LQ,2)*gamma1(a1,a)*gamma1(q,a1))/2. + LQ*(-(C1(q,a1)*gamma1(a1,a)) - gamma2(q,a)) + LF*(C1(q,a1)*gamma1(a1,a) - LQ*gamma1(a1,a)*gamma1(q,a1) + gamma2(q,a))) + del(q,a)*((-(B1q*LQ) + beta0*LQ - (A1q*pow(LQ,2))/2.)*C1(qb,b) + C2(qb,b) + (B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2. + beta0*(pow(LF,2)/2. - pow(LQ,2)/2.) + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.))*gamma1(qb,b) + (pow(LF,2)*gamma1(b1,b)*gamma1(qb,b1))/2. + (pow(LQ,2)*gamma1(b1,b)*gamma1(qb,b1))/2. + LQ*(-(C1(qb,b1)*gamma1(b1,b)) - gamma2(qb,b)) + LF*(C1(qb,b1)*gamma1(b1,b) - LQ*gamma1(b1,b)*gamma1(qb,b1) + gamma2(qb,b)));
	      }

	  *H2st.get(a,b) /= pow(11.,2);
	      /*
	      for (int p1 = 0; p1 < allf.size(); p1++)
		{
		  int a1 = allf[p1];
		  int b1 = allf[p1];
		  *H2st.get(a,b) += C1(q,a)*C1(qb,b)
		    + (-(B2q*LQ) + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8. + beta0*(-(B1q*pow(LQ,2))/2. - (A1q*pow(LQ,3))/3.))*del(q,a)*del(qb,b)
		    + (LF-LQ)*C1(qb,b)*gamma1(q,a)
		    + (0.
		       + (LF-LQ)*C1(q,a)
		       + (pow(LF,2) - 2*LF*LQ + pow(LQ,2))*gamma1(q,a)
		       )*gamma1(qb,b)
		    - beta0*LR*((-(B1q*LQ) - (A1q*pow(LQ,2))/2.)*del(q,a)*del(qb,b) + del(qb,b)*(C1(q,a) + (LF - LQ)*gamma1(q,a)) + del(q,a)*(C1(qb,b) + (LF - LQ)*gamma1(qb,b)))
		    + del(qb,b)*(
				 +(0.
				   + beta0*LQ
				   -(B1q*LQ)
				   - (A1q*pow(LQ,2))/2.
				   )*C1(q,a)
				 + C2(q,a)
				 + (0.
				    + B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2.
				    + beta0*(pow(LF,2)/2. - pow(LQ,2)/2.)
				    + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.)
				    )*gamma1(q,a)
				 + (pow(LF,2)*gamma1(a1,a)*gamma1(q,a1))/2.
				 + (pow(LQ,2)*gamma1(a1,a)*gamma1(q,a1))/2.
				 + LQ*(0.
				       -(C1(q,a1)*gamma1(a1,a))
				       - gamma2(q,a)
				       )
				 + LF*(0.
				       + C1(q,a1)*gamma1(a1,a)
				       - LQ*gamma1(a1,a)*gamma1(q,a1)
				       + gamma2(q,a)
				       ))
		    + del(q,a)*(
				+(0.
				  + beta0*LQ
				  -(B1q*LQ)
				  - (A1q*pow(LQ,2))/2.
				  )*C1(qb,b)
				+ C2(qb,b)
				+ (0.
				   + B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2.
				   + beta0*(pow(LF,2)/2. - pow(LQ,2)/2.)
				   + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.)
				   )*gamma1(qb,b)
				+ (pow(LF,2)*gamma1(b1,b)*gamma1(qb,b1))/2.
				+ (pow(LQ,2)*gamma1(b1,b)*gamma1(qb,b1))/2.
				+ LQ*(0.
				      -(C1(qb,b1)*gamma1(b1,b))
				      - gamma2(qb,b)
				      )
				+ LF*(0.
				      + C1(qb,b1)*gamma1(b1,b)
				      - LQ*gamma1(b1,b)*gamma1(qb,b1)
				      + gamma2(qb,b)
				      ));
		}
	      *H2st.get(a,b) /= 11.;
	      */
	}

      double as2 = pow(aass,2);
      Hqqb[i]  += as2*H2st(q,qb);
      Hqg[i]   += as2*H2st(q,g);
      Hqq[i]   += as2*H2st(q,q);
      Hqqp[i]  += as2*H2st(q,qp);
      Hqqbp[i] += as2*H2st(q,qbp);
      Hgg[i]   += as2*H2st(g,g);

      //cout << endl;
      //cout << "Hst_scale2.txt " << endl;
      //cout << "H2st(q,qb)  " << H2st(q,qb)  << endl;
      //cout << "H2st(q,g)   " << H2st(q,g)   << endl;
      //cout << "H2st(q,q)   " << H2st(q,q)   << endl;
      //cout << "H2st(q,qp)  " << H2st(q,qp)  << endl;
      //cout << "H2st(q,qbp) " << H2st(q,qbp) << endl;
      //cout << "H2st(g,g)   " << H2st(g,g)   << endl;
      
      if (opts.order == 2)
	continue;

      fpmap <complex <double>> C3;
      *C3.get(q,q)   = *C3.get(qb,qb) = *C3.get(qp,qp)  = *C3.get(qbp,qbp) = C3qq[idx];
      *C3.get(q,g)   = *C3.get(g,q)   = *C3.get(qb,g)   = *C3.get(g,qb)    = C3qg[idx];
      *C3.get(qp,g)  = *C3.get(g,qp)  = *C3.get(qbp,g)  = *C3.get(g,qbp)   = C3qg[idx];
      *C3.get(q,qb)  = *C3.get(qb,q)  = *C3.get(qp,qbp) = *C3.get(qbp,qp)  = C3qqb[idx];
      *C3.get(q,qp)  = *C3.get(qp,q)  = *C3.get(qb,qbp) = *C3.get(qbp,qb)  = C3qqp[idx];
      *C3.get(q,qbp) = *C3.get(qbp,q) = *C3.get(qb,qp)  = *C3.get(qp,qb)   = C3qqbp[idx];

      fpmap <complex <double>> gamma3;
      *gamma3.get(q,q)   = *gamma3.get(qb,qb) = *gamma3.get(qp,qp)  = *gamma3.get(qbp,qbp) = pmom::gamma3qq[idx];
      *gamma3.get(q,g)   = *gamma3.get(qb,g)  = *gamma3.get(qp,g)   = *gamma3.get(qbp,g)   = pmom::gamma3qg[idx];
      *gamma3.get(g,q)   = *gamma3.get(g,qb)  = *gamma3.get(g,qp)   = *gamma3.get(g,qbp)   = pmom::gamma3gq[idx];
      *gamma3.get(q,qb)  = *gamma3.get(qb,q)  = *gamma3.get(qp,qbp) = *gamma3.get(qbp,qp)  = pmom::gamma3qqb[idx];
      *gamma3.get(q,qp)  = *gamma3.get(qp,q)  = *gamma3.get(qb,qbp) = *gamma3.get(qbp,qb)  = pmom::gamma3qqp[idx];
      *gamma3.get(q,qbp) = *gamma3.get(qbp,q) = *gamma3.get(qb,qp)  = *gamma3.get(qp,qb)   = pmom::gamma3qqbp[idx];
      *gamma3.get(g,g)   = pmom::gamma3gg[idx];

      fpmap <complex <double>> H3st;

      pa.push_back(qb);  pb.push_back(g);
      pa.push_back(qp);  pb.push_back(g);
      pa.push_back(qbp); pb.push_back(g);
	  
      clock_t begin_time, end_time;

      begin_time = clock();
      for (int p = 0; p < pa.size(); p++)
	{
	  parton::partid a = pa[p];
	  parton::partid b = pb[p];
	  *H3st.get(a,b) = 0;

	  for (int p1 = 0; p1 < allf.size(); p1++)
	    for (int p2 = 0; p2 < allf.size(); p2++)
	      for (int p3 = 0; p3 < allf.size(); p3++)
		for (int p4 = 0; p4 < allf.size(); p4++)
		  {
		    int a1 = allf[p1];
		    int a2 = allf[p2];
		    int a3 = allf[p3];
		    int a4 = allf[p4];
		    int aa1 = allf[p1];
		    int b1 = allf[p1];
		    int b2 = allf[p2];
		    int b3 = allf[p3];
		    int b4 = allf[p4];
		    int bb1 = allf[p1];

		    *H3st.get(a,b) +=
		      (-(B3q*LQ) + (-A3q/2. + B1q*B2q - (B1q*beta1)/2.)*pow(LQ,2) + ((A2q*B1q)/2. - pow(B1q,3)/6. + (A1q*B2q)/2. - (A1q*beta1)/3.)*pow(LQ,3) + A1q*(A2q/4. - pow(B1q,2)/4.)*pow(LQ,4) - (pow(A1q,2)*B1q*pow(LQ,5))/8. - (pow(A1q,3)*pow(LQ,6))/48. + pow(beta0,2)*(-(B1q*pow(LQ,3))/3. - (A1q*pow(LQ,4))/4.) + beta0*(-(B2q*pow(LQ,2)) + ((-2*A2q)/3. + pow(B1q,2)/2.)*pow(LQ,3) + (7.*A1q*B1q*pow(LQ,4))/12. + (pow(A1q,2)*pow(LQ,5))/6.))*del(q, a)*del(qb, b) + (-(beta1*LR) + pow(beta0,2)*pow(LR,2))*((-(B1q*LQ) - (A1q*pow(LQ,2))/2.)*del(q, a)*del(qb, b) + del(qb, b)*(C1(q, a) + (LF - LQ)*gamma1(q, a)) + del(q, a)*(C1(qb, b) + (LF - LQ)*gamma1(qb, b))) + C1(qb, b)*((-(B1q*LQ) + 2.*beta0*LQ - (A1q*pow(LQ,2))/2.)*C1(q, a) + C2(q, a) + (pow(LF,2)*gamma1(a3, a)*gamma1(q, a3))/2. + (pow(LQ,2)*gamma1(a3, a)*gamma1(q, a3))/2. + LQ*(-(C1(q, a1)*gamma1(a1, a)) - gamma2(q, a)) + LF*(C1(q, a2)*gamma1(a2, a) - LQ*gamma1(a2, a)*gamma1(q, a2) + gamma2(q, a))) + gamma1(qb, b)*((B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2. + beta0*(pow(LF,2)/2. + LF*LQ - (3.*pow(LQ,2))/2.) + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.))*C1(q, a) - LQ*C2(q, a) + (-(B1q*pow(LQ,3)) - (A1q*pow(LQ,4))/2. + pow(LF,2)*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.) + beta0*(pow(LF,3) - pow(LF,2)*LQ - LF*pow(LQ,2) + pow(LQ,3)) + LF*(2.*B1q*pow(LQ,2) + A1q*pow(LQ,3)))*gamma1(q, a) + (pow(LF,3)*gamma1(a3, a)*gamma1(q, a3))/2. - (pow(LQ,3)*gamma1(a3, a)*gamma1(q, a3))/2. + LF*(C2(q, a) + pow(LQ,2)*(gamma1(a2, a)*gamma1(q, a2) + (gamma1(a3, a)*gamma1(q, a3))/2.) + LQ*(-(C1(q, a1)*gamma1(a1, a)) - C1(q, a2)*gamma1(a2, a) - 2.*gamma2(q, a))) + pow(LQ,2)*(C1(q, a1)*gamma1(a1, a) + gamma2(q, a)) + pow(LF,2)*(C1(q, a2)*gamma1(a2, a) + LQ*(-(gamma1(a2, a)*gamma1(q, a2)) - (gamma1(a3, a)*gamma1(q, a3))/2.) + gamma2(q, a))) + C1(q, a)*(C2(qb, b) + (pow(LF,2)*gamma1(b3, b)*gamma1(qb, b3))/2. + (pow(LQ,2)*gamma1(b3, b)*gamma1(qb, b3))/2. + LQ*(-(C1(qb, b1)*gamma1(b1, b)) - gamma2(qb, b)) + LF*(C1(qb, b2)*gamma1(b2, b) - LQ*gamma1(b2, b)*gamma1(qb, b2) + gamma2(qb, b))) + gamma1(q, a)*((B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2. + beta0*(pow(LF,2)/2. + LF*LQ - (3.*pow(LQ,2))/2.) + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.))*C1(qb, b) - LQ*C2(qb, b) + (pow(LF,3)*gamma1(b3, b)*gamma1(qb, b3))/2. - (pow(LQ,3)*gamma1(b3, b)*gamma1(qb, b3))/2. + LF*(C2(qb, b) + pow(LQ,2)*(gamma1(b2, b)*gamma1(qb, b2) + (gamma1(b3, b)*gamma1(qb, b3))/2.) + LQ*(-(C1(qb, b1)*gamma1(b1, b)) - C1(qb, b2)*gamma1(b2, b) - 2.*gamma2(qb, b))) + pow(LQ,2)*(C1(qb, b1)*gamma1(b1, b) + gamma2(qb, b)) + pow(LF,2)*(C1(qb, b2)*gamma1(b2, b) + LQ*(-(gamma1(b2, b)*gamma1(qb, b2)) - (gamma1(b3, b)*gamma1(qb, b3))/2.) + gamma2(qb, b))) - 2.*beta0*LR*(C1(q, a)*C1(qb, b) + (-(B2q*LQ) + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8. + beta0*(-(B1q*pow(LQ,2))/2. - (A1q*pow(LQ,3))/3.))*del(q, a)*del(qb, b) + (LF - LQ)*C1(qb, b)*gamma1(q, a) + ((LF - LQ)*C1(q, a) + (pow(LF,2) - 2.*LF*LQ + pow(LQ,2))*gamma1(q, a))*gamma1(qb, b) + del(qb, b)*((-(B1q*LQ) + beta0*LQ - (A1q*pow(LQ,2))/2.)*C1(q, a) + C2(q, a) + (B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2. + beta0*(pow(LF,2)/2. - pow(LQ,2)/2.) + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.))*gamma1(q, a) + (pow(LF,2)*gamma1(a1, a)*gamma1(q, a1))/2. + (pow(LQ,2)*gamma1(a1, a)*gamma1(q, a1))/2. + LQ*(-(C1(q, a1)*gamma1(a1, a)) - gamma2(q, a)) + LF*(C1(q, a1)*gamma1(a1, a) - LQ*gamma1(a1, a)*gamma1(q, a1) + gamma2(q, a))) + del(q, a)*((-(B1q*LQ) + beta0*LQ - (A1q*pow(LQ,2))/2.)*C1(qb, b) + C2(qb, b) + (B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2. + beta0*(pow(LF,2)/2. - pow(LQ,2)/2.) + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.))*gamma1(qb, b) + (pow(LF,2)*gamma1(b1, b)*gamma1(qb, b1))/2. + (pow(LQ,2)*gamma1(b1, b)*gamma1(qb, b1))/2. + LQ*(-(C1(qb, b1)*gamma1(b1, b)) - gamma2(qb, b)) + LF*(C1(qb, b1)*gamma1(b1, b) - LQ*gamma1(b1, b)*gamma1(qb, b1) + gamma2(qb, b)))) + del(qb, b)*(((-B2q + beta1)*LQ + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + pow(beta0,2)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8. + beta0*((-3*B1q*pow(LQ,2))/2. - (5.*A1q*pow(LQ,3))/6.))*C1(q, a) + C3(q, a) + ((beta1*pow(LF,2))/2. + (B2q - beta1/2.)*pow(LQ,2) + (A2q/2. - pow(B1q,2)/2.)*pow(LQ,3) - (A1q*B1q*pow(LQ,4))/2. - (pow(A1q,2)*pow(LQ,5))/8. + pow(beta0,2)*(pow(LF,3)/3. - pow(LQ,3)/3.) + LF*(-(B2q*LQ) + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8.) + beta0*(B1q*pow(LQ,3) + (7.*A1q*pow(LQ,4))/12. + pow(LF,2)*(-(B1q*LQ)/2. - (A1q*pow(LQ,2))/4.) + LF*(-(B1q*pow(LQ,2))/2. - (A1q*pow(LQ,3))/3.)))*gamma1(q, a) - (A1q*pow(LQ,4)*gamma1(a3, a)*gamma1(q, a3))/4. + (pow(LF,3)*gamma1(a3, a4)*gamma1(a4, a)*gamma1(q, a3))/6. + pow(LF,2)*((C1(q, a2)*gamma1(a2, a3)*gamma1(a3, a))/2. - (A1q*pow(LQ,2)*gamma1(a3, a)*gamma1(q, a3))/4. + LQ*(-(gamma1(a2, a3)*gamma1(a3, a)*gamma1(q, a2))/2. - (B1q*gamma1(a3, a)*gamma1(q, a3))/2.) + gamma1(q, a3)*gamma2(a3, a)) + pow(LQ,3)*(-(B1q*gamma1(a3, a)*gamma1(q, a3))/2. - (gamma1(a3, a4)*gamma1(a4, a)*gamma1(q, a3))/6. + A1q*((C1(q, a1)*gamma1(a1, a))/2. + gamma2(q, a)/2.)) + pow(LQ,2)*(-(A1q*C2(q, a))/2. + (C1(q, a1)*gamma1(a1, a3)*gamma1(a3, a))/2. + gamma1(q, a3)*gamma2(a3, a) + B1q*(C1(q, a1)*gamma1(a1, a) + gamma2(q, a))) + beta0*(LQ*(C1(a1, a)*C1(q, a1) - C1(aa1, a)*C1(q, aa1) + 2.*C2(q, a)) + LF*(LQ*C1(q, a2)*gamma1(a2, a) - (pow(LQ,2)*gamma1(a2, a)*gamma1(q, a2))/2.) + (pow(LF,3)*gamma1(a3, a)*gamma1(q, a3))/2. + (pow(LQ,3)*gamma1(a3, a)*gamma1(q, a3))/2. + pow(LQ,2)*(-(C1(q, a1)*gamma1(a1, a))/2. - C1(q, a3)*gamma1(a3, a) - gamma2(q, a)) + pow(LF,2)*((C1(q, a2)*gamma1(a2, a))/2. - (LQ*gamma1(a2, a)*gamma1(q, a2))/2. + gamma2(q, a))) + LQ*(-(B1q*C2(q, a)) - C2(q, a1)*gamma1(a1, a) - C1(q, a1)*gamma2(a1, a) - gamma3(q, a)) + LF*(C2(q, a2)*gamma1(a2, a) + (A1q*pow(LQ,3)*gamma1(a2, a)*gamma1(q, a2))/2. + C1(q, a2)*gamma2(a2, a) + pow(LQ,2)*(B1q*gamma1(a2, a)*gamma1(q, a2) + (gamma1(a2, a)*gamma1(a3, a2)*gamma1(q, a3))/2. + A1q*(-(C1(q, a2)*gamma1(a2, a))/2. - gamma2(q, a)/2.)) + LQ*(-(C1(q, a1)*gamma1(a1, a2)*gamma1(a2, a)) - gamma1(q, a2)*gamma2(a2, a) + B1q*(-(C1(q, a2)*gamma1(a2, a)) - gamma2(q, a)) - gamma1(a2, a)*gamma2(q, a2)) + gamma3(q, a))) + del(q, a)*(((-B2q + beta1)*LQ + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + pow(beta0,2)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8. + beta0*((-3*B1q*pow(LQ,2))/2. - (5.*A1q*pow(LQ,3))/6.))*C1(qb, b) + C3(qb, b) + ((beta1*pow(LF,2))/2. + (B2q - beta1/2.)*pow(LQ,2) + (A2q/2. - pow(B1q,2)/2.)*pow(LQ,3) - (A1q*B1q*pow(LQ,4))/2. - (pow(A1q,2)*pow(LQ,5))/8. + pow(beta0,2)*(pow(LF,3)/3. - pow(LQ,3)/3.) + LF*(-(B2q*LQ) + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8.) + beta0*(B1q*pow(LQ,3) + (7.*A1q*pow(LQ,4))/12. + pow(LF,2)*(-(B1q*LQ)/2. - (A1q*pow(LQ,2))/4.) + LF*(-(B1q*pow(LQ,2))/2. - (A1q*pow(LQ,3))/3.)))*gamma1(qb, b) - (A1q*pow(LQ,4)*gamma1(b3, b)*gamma1(qb, b3))/4. + (pow(LF,3)*gamma1(b3, b4)*gamma1(b4, b)*gamma1(qb, b3))/6. + pow(LF,2)*((C1(qb, b2)*gamma1(b2, b3)*gamma1(b3, b))/2. - (A1q*pow(LQ,2)*gamma1(b3, b)*gamma1(qb, b3))/4. + LQ*(-(gamma1(b2, b3)*gamma1(b3, b)*gamma1(qb, b2))/2. - (B1q*gamma1(b3, b)*gamma1(qb, b3))/2.) + gamma1(qb, b3)*gamma2(b3, b)) + pow(LQ,3)*(-(B1q*gamma1(b3, b)*gamma1(qb, b3))/2. - (gamma1(b3, b4)*gamma1(b4, b)*gamma1(qb, b3))/6. + A1q*((C1(qb, b1)*gamma1(b1, b))/2. + gamma2(qb, b)/2.)) + pow(LQ,2)*(-(A1q*C2(qb, b))/2. + (C1(qb, b1)*gamma1(b1, b3)*gamma1(b3, b))/2. + gamma1(qb, b3)*gamma2(b3, b) + B1q*(C1(qb, b1)*gamma1(b1, b) + gamma2(qb, b))) + beta0*(LQ*(C1(b1, b)*C1(qb, b1) - C1(bb1, b)*C1(qb, bb1) + 2.*C2(qb, b)) + LF*(LQ*C1(qb, b2)*gamma1(b2, b) - (pow(LQ,2)*gamma1(b2, b)*gamma1(qb, b2))/2.) + (pow(LF,3)*gamma1(b3, b)*gamma1(qb, b3))/2. + (pow(LQ,3)*gamma1(b3, b)*gamma1(qb, b3))/2. + pow(LQ,2)*(-(C1(qb, b1)*gamma1(b1, b))/2. - C1(qb, b3)*gamma1(b3, b) - gamma2(qb, b)) + pow(LF,2)*((C1(qb, b2)*gamma1(b2, b))/2. - (LQ*gamma1(b2, b)*gamma1(qb, b2))/2. + gamma2(qb, b))) + LQ*(-(B1q*C2(qb, b)) - C2(qb, b1)*gamma1(b1, b) - C1(qb, b1)*gamma2(b1, b) - gamma3(qb, b)) + LF*(C2(qb, b2)*gamma1(b2, b) + (A1q*pow(LQ,3)*gamma1(b2, b)*gamma1(qb, b2))/2. + C1(qb, b2)*gamma2(b2, b) + pow(LQ,2)*(B1q*gamma1(b2, b)*gamma1(qb, b2) + (gamma1(b2, b)*gamma1(b3, b2)*gamma1(qb, b3))/2. + A1q*(-(C1(qb, b2)*gamma1(b2, b))/2. - gamma2(qb, b)/2.)) + LQ*(-(C1(qb, b1)*gamma1(b1, b2)*gamma1(b2, b)) - gamma1(qb, b2)*gamma2(b2, b) + B1q*(-(C1(qb, b2)*gamma1(b2, b)) - gamma2(qb, b)) - gamma1(b2, b)*gamma2(qb, b2)) + gamma3(qb, b)));		    
		  }
	  *H3st.get(a,b) /= pow(11.,4);
	}
      end_time = clock();
	  
      cout << endl;
      cout << "Hst_scale2.txt " << endl;
      cout << H3st(q,qb) << endl;
      cout << H3st(q,g) << endl;
      cout << H3st(q,q) << endl;
      cout << H3st(q,qp) << endl;
      cout << H3st(q,qbp) << endl;
      cout << H3st(g,g) << endl;
      cout << H3st(qb,g)  << endl;
      cout << H3st(qp,g)  << endl;
      cout << H3st(qbp,g) << endl;

      cout << "tot time: " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << "s" << endl;

      /******************** Fast version ********************/
      begin_time = clock();
      for (int p = 0; p < pa.size(); p++)
	{
	  parton::partid a = pa[p];
	  parton::partid b = pb[p];
	  *H3st.get(a,b) =
	    (-(B3q*LQ) + (-A3q/2. + B1q*B2q - (B1q*beta1)/2.)*pow(LQ,2) + ((A2q*B1q)/2. - pow(B1q,3)/6. + (A1q*B2q)/2. - (A1q*beta1)/3.)*pow(LQ,3)
	     + A1q*(A2q/4. - pow(B1q,2)/4.)*pow(LQ,4) - (pow(A1q,2)*B1q*pow(LQ,5))/8.
	     - (pow(A1q,3)*pow(LQ,6))/48.
	     + pow(beta0,2)*(-(B1q*pow(LQ,3))/3. - (A1q*pow(LQ,4))/4.)
	     + beta0*(-(B2q*pow(LQ,2)) + ((-2*A2q)/3. + pow(B1q,2)/2.)*pow(LQ,3) + (7.*A1q*B1q*pow(LQ,4))/12. + (pow(A1q,2)*pow(LQ,5))/6.))*del(q, a)*del(qb, b)
	    + (-(beta1*LR) + pow(beta0,2)*pow(LR,2))*((-(B1q*LQ) - (A1q*pow(LQ,2))/2.)*del(q, a)*del(qb, b)
						      + del(qb, b)*(C1(q, a) + (LF - LQ)*gamma1(q, a))
						      + del(q, a)*(C1(qb, b) + (LF - LQ)*gamma1(qb, b)))
	    + del(qb, b)*(((-B2q + beta1)*LQ + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + pow(beta0,2)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2.
			   + (pow(A1q,2)*pow(LQ,4))/8. + beta0*((-3*B1q*pow(LQ,2))/2. - (5.*A1q*pow(LQ,3))/6.))*C1(q, a)
			  + C3(q, a)
			  + ((beta1*LF2)/2. + (B2q - beta1/2.)*pow(LQ,2) + (A2q/2. - pow(B1q,2)/2.)*pow(LQ,3) - (A1q*B1q*pow(LQ,4))/2. - (pow(A1q,2)*pow(LQ,5))/8.
			     + pow(beta0,2)*(LF3/3. - pow(LQ,3)/3.) + LF*(-(B2q*LQ) + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8.)
			     + beta0*(B1q*pow(LQ,3) + (7.*A1q*pow(LQ,4))/12. + LF2*(-(B1q*LQ)/2. - (A1q*pow(LQ,2))/4.) + LF*(-(B1q*pow(LQ,2))/2. - (A1q*pow(LQ,3))/3.)))*gamma1(q, a))
	    + del(q, a)*(((-B2q + beta1)*LQ + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + pow(beta0,2)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8.
			  + beta0*((-3*B1q*pow(LQ,2))/2. - (5.*A1q*pow(LQ,3))/6.))*C1(qb, b) + C3(qb, b)
			 + ((beta1*LF2)/2. + (B2q - beta1/2.)*pow(LQ,2) + (A2q/2. - pow(B1q,2)/2.)*pow(LQ,3) - (A1q*B1q*pow(LQ,4))/2. - (pow(A1q,2)*pow(LQ,5))/8.
			    + pow(beta0,2)*(LF3/3. - pow(LQ,3)/3.)
			    + LF*(-(B2q*LQ) + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8.)
			    + beta0*(B1q*pow(LQ,3) + (7.*A1q*pow(LQ,4))/12. + LF2*(-(B1q*LQ)/2. - (A1q*pow(LQ,2))/4.) + LF*(-(B1q*pow(LQ,2))/2. - (A1q*pow(LQ,3))/3.)))*gamma1(qb, b))
	    + del(qb, b)*(LF*(gamma3(q, a)))
	    + del(q, a)*(LF*(+ gamma3(qb, b)))
	    + del(qb, b)*(+ pow(LQ,2)*(-(A1q*C2(q, a))/2.))
	    + del(q, a)*(+ pow(LQ,2)*(-(A1q*C2(qb, b))/2.));

	  *H3st.get(a,b) *= pow(11.,2);

	  for (int p1 = 0; p1 < allf.size(); p1++)
	    {
	      int a1 = allf[p1];
	      int a2 = allf[p1];
	      int a3 = allf[p1];
	      int aa1 = allf[p1];
	      int b1 = allf[p1];
	      int b2 = allf[p1];
	      int b3 = allf[p1];
	      int bb1 = allf[p1];
	      *H3st.get(a,b) +=
		+ C1(qb, b)*((-(B1q*LQ) + 2.*beta0*LQ - (A1q*pow(LQ,2))/2.)*C1(q, a)
			     + C2(q, a)
			     + (LF2*gamma1(a3, a)*gamma1(q, a3))/2.
			     + (pow(LQ,2)*gamma1(a3, a)*gamma1(q, a3))/2.
			     + LQ*(-(C1(q, a1)*gamma1(a1, a)) - gamma2(q, a))
			     + LF*(C1(q, a2)*gamma1(a2, a) - LQ*gamma1(a2, a)*gamma1(q, a2) + gamma2(q, a)))
		*11.;
	      *H3st.get(a,b) +=
		+ gamma1(qb, b)*((B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2.
				  + beta0*(LF2/2. + LF*LQ - (3.*pow(LQ,2))/2.)
				  + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.))*C1(q, a)
				 - LQ*C2(q, a)
				 + (-(B1q*pow(LQ,3)) - (A1q*pow(LQ,4))/2.
				    + LF2*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.)
				    + beta0*(LF3 - LF2*LQ - LF*pow(LQ,2) + pow(LQ,3))
				    + LF*(2.*B1q*pow(LQ,2) + A1q*pow(LQ,3)))*gamma1(q, a)
				 + (LF3*gamma1(a3, a)*gamma1(q, a3))/2.
				 - (pow(LQ,3)*gamma1(a3, a)*gamma1(q, a3))/2.
				 + LF*(C2(q, a) + pow(LQ,2)*(gamma1(a2, a)*gamma1(q, a2) + (gamma1(a3, a)*gamma1(q, a3))/2.)
				       + LQ*(-(C1(q, a1)*gamma1(a1, a)) - C1(q, a2)*gamma1(a2, a) - 2.*gamma2(q, a)))
				 + pow(LQ,2)*(C1(q, a1)*gamma1(a1, a) + gamma2(q, a))
				 + LF2*(C1(q, a2)*gamma1(a2, a) + LQ*(-(gamma1(a2, a)*gamma1(q, a2)) - (gamma1(a3, a)*gamma1(q, a3))/2.) + gamma2(q, a)))
		*11.;
	      *H3st.get(a,b) +=
		+ C1(q, a)*(C2(qb, b)
			    + (LF2*gamma1(b3, b)*gamma1(qb, b3))/2.
			    + (pow(LQ,2)*gamma1(b3, b)*gamma1(qb, b3))/2.
			    + LQ*(-(C1(qb, b1)*gamma1(b1, b)) - gamma2(qb, b))
			    + LF*(C1(qb, b2)*gamma1(b2, b) - LQ*gamma1(b2, b)*gamma1(qb, b2) + gamma2(qb, b)))
		*11.;
	      *H3st.get(a,b) +=
		+ gamma1(q, a)*((B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2.
				 + beta0*(LF2/2. + LF*LQ - (3.*pow(LQ,2))/2.)
				 + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.))*C1(qb, b)
				- LQ*C2(qb, b)
				+ (LF3*gamma1(b3, b)*gamma1(qb, b3))/2. - (pow(LQ,3)*gamma1(b3, b)*gamma1(qb, b3))/2.
				+ LF*(C2(qb, b) + pow(LQ,2)*(gamma1(b2, b)*gamma1(qb, b2) + (gamma1(b3, b)*gamma1(qb, b3))/2.)
				      + LQ*(-(C1(qb, b1)*gamma1(b1, b)) - C1(qb, b2)*gamma1(b2, b) - 2.*gamma2(qb, b)))
				+ pow(LQ,2)*(C1(qb, b1)*gamma1(b1, b) + gamma2(qb, b))
				+ LF2*(C1(qb, b2)*gamma1(b2, b) + LQ*(-(gamma1(b2, b)*gamma1(qb, b2)) - (gamma1(b3, b)*gamma1(qb, b3))/2.) + gamma2(qb, b)))
		*11.;
	      *H3st.get(a,b) +=
		- 2.*beta0*LR*(C1(q, a)*C1(qb, b)
			       + (-(B2q*LQ) + (-A2q/2. + pow(B1q,2)/2.)*pow(LQ,2) + (A1q*B1q*pow(LQ,3))/2. + (pow(A1q,2)*pow(LQ,4))/8.
				  + beta0*(-(B1q*pow(LQ,2))/2. - (A1q*pow(LQ,3))/3.))*del(q, a)*del(qb, b)
			       + (LF - LQ)*C1(qb, b)*gamma1(q, a)
			       + ((LF - LQ)*C1(q, a) + (LF2 - 2.*LF*LQ + pow(LQ,2))*gamma1(q, a))*gamma1(qb, b)
			       + del(qb, b)*(
					     (-(B1q*LQ) + beta0*LQ - (A1q*pow(LQ,2))/2.)*C1(q, a)
					     + C2(q, a)
					     + (B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2. + beta0*(LF2/2. - pow(LQ,2)/2.) + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.))*gamma1(q, a)
					     + (LF2*gamma1(a1, a)*gamma1(q, a1))/2.
					     + (pow(LQ,2)*gamma1(a1, a)*gamma1(q, a1))/2.
					     + LQ*(-(C1(q, a1)*gamma1(a1, a)) - gamma2(q, a))
					     + LF*(C1(q, a1)*gamma1(a1, a) - LQ*gamma1(a1, a)*gamma1(q, a1) + gamma2(q, a)))
			       + del(q, a)*(
					    (-(B1q*LQ) + beta0*LQ - (A1q*pow(LQ,2))/2.)*C1(qb, b)
					    + C2(qb, b)
					    + (B1q*pow(LQ,2) + (A1q*pow(LQ,3))/2. + beta0*(LF2/2. - pow(LQ,2)/2.) + LF*(-(B1q*LQ) - (A1q*pow(LQ,2))/2.))*gamma1(qb, b)
					    + (LF2*gamma1(b1, b)*gamma1(qb, b1))/2.
					    + (pow(LQ,2)*gamma1(b1, b)*gamma1(qb, b1))/2.
					    + LQ*(-(C1(qb, b1)*gamma1(b1, b)) - gamma2(qb, b))
					    + LF*(C1(qb, b1)*gamma1(b1, b) - LQ*gamma1(b1, b)*gamma1(qb, b1) + gamma2(qb, b))))
		*11.;
	      *H3st.get(a,b) +=
		+ del(qb, b)*(- (A1q*pow(LQ,4)*gamma1(a3, a)*gamma1(q, a3))/4.)
		*11.;
	      *H3st.get(a,b) +=
		+ del(qb, b)*(+ beta0*(LQ*(C1(a1, a)*C1(q, a1) - C1(aa1, a)*C1(q, aa1) + 2.*C2(q, a)) + LF*(LQ*C1(q, a2)*gamma1(a2, a) - (pow(LQ,2)*gamma1(a2, a)*gamma1(q, a2))/2.)
				       + (LF3*gamma1(a3, a)*gamma1(q, a3))/2. + (pow(LQ,3)*gamma1(a3, a)*gamma1(q, a3))/2.
				       + pow(LQ,2)*(-(C1(q, a1)*gamma1(a1, a))/2. - C1(q, a3)*gamma1(a3, a) - gamma2(q, a))
				       + LF2*((C1(q, a2)*gamma1(a2, a))/2. - (LQ*gamma1(a2, a)*gamma1(q, a2))/2. + gamma2(q, a)))
			      + LQ*(-(B1q*C2(q, a)) - C2(q, a1)*gamma1(a1, a) - C1(q, a1)*gamma2(a1, a) - gamma3(q, a)))
		*11.;
	      *H3st.get(a,b) +=
		+ del(q, a)*(- (A1q*pow(LQ,4)*gamma1(b3, b)*gamma1(qb, b3))/4.)	
		*11.;
	      *H3st.get(a,b) +=
		+ del(q, a)*(+ beta0*(LQ*(C1(b1, b)*C1(qb, b1) - C1(bb1, b)*C1(qb, bb1) + 2.*C2(qb, b))
				      + LF*(LQ*C1(qb, b2)*gamma1(b2, b) - (pow(LQ,2)*gamma1(b2, b)*gamma1(qb, b2))/2.)
				      + (LF3*gamma1(b3, b)*gamma1(qb, b3))/2.
				      + (pow(LQ,3)*gamma1(b3, b)*gamma1(qb, b3))/2.
				      + pow(LQ,2)*(-(C1(qb, b1)*gamma1(b1, b))/2. - C1(qb, b3)*gamma1(b3, b) - gamma2(qb, b))
				      + LF2*((C1(qb, b2)*gamma1(b2, b))/2. - (LQ*gamma1(b2, b)*gamma1(qb, b2))/2. + gamma2(qb, b)))
			     + LQ*(-(B1q*C2(qb, b)) - C2(qb, b1)*gamma1(b1, b) - C1(qb, b1)*gamma2(b1, b) - gamma3(qb, b)))
		*11.;
	      *H3st.get(a,b) +=
		+ del(qb, b)*(+ LF*(C2(q, a2)*gamma1(a2, a) + (A1q*pow(LQ,3)*gamma1(a2, a)*gamma1(q, a2))/2.)
			      + LF*(C1(q, a2)*gamma2(a2, a)))
		*11.;
	      *H3st.get(a,b) +=
		+ del(q, a)*(+ LF*(C2(qb, b2)*gamma1(b2, b) + (A1q*pow(LQ,3)*gamma1(b2, b)*gamma1(qb, b2))/2.)
			     + LF*(+ C1(qb, b2)*gamma2(b2, b)))
		*11.;
	      
	      *H3st.get(a,b) +=
		+ del(qb, b)*(+ pow(LQ,3)*(-(B1q*gamma1(a3, a)*gamma1(q, a3))/2.))
		*11.;
	      *H3st.get(a,b) +=
		+ del(qb, b)*(+ pow(LQ,3)*(+ A1q*((C1(q, a1)*gamma1(a1, a))/2. + gamma2(q, a)/2.)))
		*11.;
	      *H3st.get(a,b) +=
		  + del(q, a)*(+ pow(LQ,3)*(-(B1q*gamma1(b3, b)*gamma1(qb, b3))/2.))
		*11.;
	      *H3st.get(a,b) +=
		+ del(q, a)*(+ pow(LQ,3)*(+ A1q*((C1(qb, b1)*gamma1(b1, b))/2. + gamma2(qb, b)/2.)))
	      *11.;

	      *H3st.get(a,b) +=
		  + del(qb, b)*(+ LF2*(- (A1q*pow(LQ,2)*gamma1(a3, a)*gamma1(q, a3))/4.))
		*11.;
	      *H3st.get(a,b) +=
		  + del(qb, b)*(+ LF2*(+ LQ*(- (B1q*gamma1(a3, a)*gamma1(q, a3))/2.)))
		*11.;
	      *H3st.get(a,b) +=
		  + del(qb, b)*(+ LF2*(+ gamma1(q, a3)*gamma2(a3, a)))
		*11.;
	      *H3st.get(a,b) +=

		  + del(qb, b)*(+ pow(LQ,2)*(+ gamma1(q, a3)*gamma2(a3, a)))
		*11.;
	      *H3st.get(a,b) +=
		  + del(qb, b)*(+ pow(LQ,2)*(+ B1q*(C1(q, a1)*gamma1(a1, a) + gamma2(q, a))))
		*11.;
	      *H3st.get(a,b) +=

		  + del(q, a)*(+ LF2*(- (A1q*pow(LQ,2)*gamma1(b3, b)*gamma1(qb, b3))/4.))
		*11.;
	      *H3st.get(a,b) +=
		  + del(q, a)*(+ LF2*(+ LQ*(- (B1q*gamma1(b3, b)*gamma1(qb, b3))/2.)))
		*11.;
	      *H3st.get(a,b) +=
		  + del(q, a)*(+ LF2*(+ gamma1(qb, b3)*gamma2(b3, b)))

		*11.;
	      *H3st.get(a,b) +=
		  + del(q, a)*(+ pow(LQ,2)*(+ gamma1(qb, b3)*gamma2(b3, b)))
		*11.;
	      *H3st.get(a,b) +=
		  + del(q, a)*(+ pow(LQ,2)*(+ B1q*(C1(qb, b1)*gamma1(b1, b) + gamma2(qb, b))))
		*11.;

	      *H3st.get(a,b) +=
		  + del(qb, b)*(+ LF*(pow(LQ,2)*(B1q*gamma1(a2, a)*gamma1(q, a2))))
		*11.;
	      *H3st.get(a,b) +=
		  + del(qb, b)*(+ LF*(pow(LQ,2)*(+ A1q*(-(C1(q, a2)*gamma1(a2, a))/2. - gamma2(q, a)/2.))))
		*11.;
	      *H3st.get(a,b) +=
		  + del(qb, b)*(+ LF*(LQ*(- gamma1(q, a2)*gamma2(a2, a))))
		*11.;
	      *H3st.get(a,b) +=
		  + del(qb, b)*(+ LF*(LQ*(+ B1q*(-(C1(q, a2)*gamma1(a2, a)) - gamma2(q, a)))))
		*11.;
	      *H3st.get(a,b) +=
		  + del(qb, b)*(+ LF*(LQ*(- gamma1(a2, a)*gamma2(q, a2))))
		*11.;
	      *H3st.get(a,b) +=

		  + del(q, a)*(+ LF*(+ pow(LQ,2)*(B1q*gamma1(b2, b)*gamma1(qb, b2))))
		*11.;
	      *H3st.get(a,b) +=
		  + del(q, a)*(+ LF*(+ pow(LQ,2)*(+ A1q*(-(C1(qb, b2)*gamma1(b2, b))/2. - gamma2(qb, b)/2.))))
		*11.;
	      *H3st.get(a,b) +=
		  + del(q, a)*(+ LF*(+ LQ*(- gamma1(qb, b2)*gamma2(b2, b))))
		*11.;
	      *H3st.get(a,b) +=
		  + del(q, a)*(+ LF*(+ LQ*(+ B1q*(-(C1(qb, b2)*gamma1(b2, b)) - gamma2(qb, b)))))
		*11.;
	      *H3st.get(a,b) +=
		  + del(q, a)*(+ LF*(+ LQ*(- gamma1(b2, b)*gamma2(qb, b2))))
		*11.;
	    }
		
	  for (int p1 = 0; p1 < allf.size(); p1++)
	    for (int p2 = 0; p2 < allf.size(); p2++)
	      {
		int a3 = allf[p1];
		int a4 = allf[p2];
		int b3 = allf[p1];
		int b4 = allf[p2];
		*H3st.get(a,b) +=
		  + del(qb, b)*(+ (LF3*gamma1(a3, a4)*gamma1(a4, a)*gamma1(q, a3))/6.)
		  + del(qb, b)*(+ LQ3*(- (gamma1(a3, a4)*gamma1(a4, a)*gamma1(q, a3))/6.))
		  + del(q, a)*(+ (LF3*gamma1(b3, b4)*gamma1(b4, b)*gamma1(qb, b3))/6.)
		  + del(q, a)*(+ LQ3*(- (gamma1(b3, b4)*gamma1(b4, b)*gamma1(qb, b3))/6.))
		  ;
		  
		  
	      }

	  for (int p1 = 0; p1 < allf.size(); p1++)
	    for (int p2 = 0; p2 < allf.size(); p2++)
	      {
		int a1 = allf[p1];
		int a2 = allf[p1];
		int a3 = allf[p2];
		int b1 = allf[p1];
		int b2 = allf[p1];
		int b3 = allf[p2];
		*H3st.get(a,b) +=
		  + del(qb, b)*(+ LF2*((C1(q, a2)*gamma1(a2, a3)*gamma1(a3, a))/2.))
		  + del(qb, b)*(+ LF2*(+ LQ*(-(gamma1(a2, a3)*gamma1(a3, a)*gamma1(q, a2))/2.)))
		  + del(qb, b)*(+ LQ2*(+ (C1(q, a1)*gamma1(a1, a3)*gamma1(a3, a))/2.))
		  + del(q, a)*(+ LF2*((C1(qb, b2)*gamma1(b2, b3)*gamma1(b3, b))/2.))
		  + del(q, a)*(+ LF2*(+ LQ*(-(gamma1(b2, b3)*gamma1(b3, b)*gamma1(qb, b2))/2.)))
		  + del(q, a)*(+ LQ2*(+ (C1(qb, b1)*gamma1(b1, b3)*gamma1(b3, b))/2.))
		  ;

	      }
	      
	  for (int p1 = 0; p1 < allf.size(); p1++)
	    for (int p2 = 0; p2 < allf.size(); p2++)
	      {
		int a1 = allf[p1];
		int a2 = allf[p2];
		int a3 = allf[p1];
		int b1 = allf[p1];
		int b2 = allf[p2];
		int b3 = allf[p1];
		  
		*H3st.get(a,b) +=
		  + del(qb, b)*(+ LF*(LQ2*(+ (gamma1(a2, a)*gamma1(a3, a2)*gamma1(q, a3))/2.)))
		  + del(qb, b)*(+ LF*(LQ*(-(C1(q, a1)*gamma1(a1, a2)*gamma1(a2, a)))))
		  + del(q, a)*(+ LF*(+ LQ2*(+ (gamma1(b2, b)*gamma1(b3, b2)*gamma1(qb, b3))/2.)))
		  + del(q, a)*(+ LF*(+ LQ*(-(C1(qb, b1)*gamma1(b1, b2)*gamma1(b2, b)))))
		  ;
		      
	      }
	  *H3st.get(a,b) /= pow(11.,2);
	}
      end_time = clock();
	  
      cout << endl;
      cout << "Hst_scale2.txt fast version" << endl;
      cout << "H3st(q,qb)  " << H3st(q,qb)  << endl;
      cout << "H3st(q,g)   " << H3st(q,g)   << endl;
      cout << "H3st(q,q)   " << H3st(q,q)   << endl;
      cout << "H3st(q,qp)  " << H3st(q,qp)  << endl;
      cout << "H3st(q,qbp) " << H3st(q,qbp) << endl;
      cout << "H3st(g,g)   " << H3st(g,g)   << endl;
      cout << "H3st(qb,g)  " << H3st(qb,g)  << endl;
      cout << "H3st(qp,g)  " << H3st(qp,g)  << endl;
      cout << "H3st(qbp,g) " << H3st(qbp,g) << endl;

      complex <double> gamma1_test_2 = 0.;
      //      parton::partid a = q;
      //      parton::partid b = g;

//for (int p1 = 0; p1 < allf.size(); p1++)
//	{
//	  int b3 = allf[p1];
//	  gamma1_test_2 +=
//	    + gamma1(q, a)*LF3*gamma1(b3, b)*gamma1(qb, b3)/2.
//	    *11.;
//	  //cout << allf[p1] << "  " << gamma1(q, a)*LF3*gamma1(b3, b)*gamma1(qb, b3)/2./11. << endl;
//	}
//      
//      for (int p1 = 0; p1 < allf.size(); p1++)
//	for (int p2 = 0; p2 < allf.size(); p2++)
//	  {
//	    int b3 = allf[p1];
//	    int b4 = allf[p2];
//	    gamma1_test_2 +=
//	      + del(q, a)*LF3*gamma1(b3, b4)*gamma1(b4, b)*gamma1(qb, b3)/6.;
//	    ;
//	    cout << allf[p1] << "  " << allf[p2] << "  " << del(q, a)*LF3*gamma1(b3, b4)*gamma1(b4, b)*gamma1(qb, b3)/6./11./11. << endl;
//	  }
//

//      parton::partid a = q;
//      parton::partid b = q;
//      
//      for (int p1 = 0; p1 < allf.size(); p1++)
//	{
//	  int b2 = allf[p1];
//	  gamma1_test_2 +=
//	    + del(q, a)*(LF*(C2(qb, b2)*gamma1(b2, b)))*11.;
//	}
//
//      gamma1_test_2 +=
//	+ gamma1(q, a)*(LF*(C2(qb, b)))*11.*11.;
//
//	      
//      gamma1_test_2 /= pow(11.,2);
//      cout << "gamma1_test_2 " << gamma1_test_2 << endl;
      
      cout << "tot time: " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << "s" << endl;

      double as3 = pow(aass,3);
      Hqqb[i]  += as3*H3st(q,qb);
      Hqg[i]   += as3*H3st(q,g);
      Hqq[i]   += as3*H3st(q,q);
      Hqqp[i]  += as3*H3st(q,qp);
      Hqqbp[i] += as3*H3st(q,qbp);
      Hgg[i]   += as3*H3st(g,g);
      Hqbg[i]  += as3*H3st(qb,g);
      Hqpg[i]  += as3*H3st(qp,g);
      Hqbpg[i] += as3*H3st(qbp,g);

      //if (opts.order == 3)
      //continue;
    }
}
