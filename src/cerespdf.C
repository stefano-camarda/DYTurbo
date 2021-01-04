#include "cerespdf.h"

#include <math.h>
#include <stdio.h>
#include <string.h>  // For NULL
#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/GridPDF.h>


#include <ceres/ceres.h>

#include "pdf.h"
#include "settings.h"
#include "gaussrules.h"
#include "zeta.h"

using ceres::SizedCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;

using namespace std;

double *cerespdf::xp;
double *cerespdf::logx;
double *cerespdf::log1mx;
double *cerespdf::bber;

//int cerespdf::N;
//int cerespdf::M;
//int cerespdf::P;

int cerespdf::np;
int cerespdf::num_observations;

cerespdf::approxpdf cerespdf::U;
cerespdf::approxpdf cerespdf::UB;
cerespdf::approxpdf cerespdf::D;
cerespdf::approxpdf cerespdf::DB;
cerespdf::approxpdf cerespdf::S;
cerespdf::approxpdf cerespdf::SB;
cerespdf::approxpdf cerespdf::C;
cerespdf::approxpdf cerespdf::CB;
cerespdf::approxpdf cerespdf::B;
cerespdf::approxpdf cerespdf::BB;
cerespdf::approxpdf cerespdf::G;

Solver::Options cerespdf::options;

void cerespdf::bernstein(int n, double t, double *bern)
{
  int i;
  int j;

  if ( n == 0 )
  {
    bern[0] = 1.0;
  }
  else if ( 0 < n )
  {
    bern[0] = 1.0 - t;
    bern[1] = t;
 
    for ( i = 2; i <= n; i++ )
    {
      bern[i] = t * bern[i-1];
      for ( j = i - 1; 1 <= j; j-- )
      {
        bern[j] =         t   * bern[j-1] 
                + ( 1.0 - t ) * bern[j];
      }
      bern[0] = ( 1.0 - t ) * bern[0];
    }
  }
  return;
}

bool cerespdf::ResAnalytic::Evaluate(double const* const* parameters,
					     double* residuals,
					     double** jacobians) const

//bool cerespdf::ResAnalytic::Evaluate(const double* const* parameters,
//					     double* residuals,
//					     double** jacobians) const

{
  //This function will have problems if x == 1 or C == 0;

  double x = x_;
  double y = y_;
  int i = i_;
       
  //Approximating function: f(x) = C * x^(a-1) * (1-x)^(b-1) * *l(x) * p(x)
  double a = parameters[0][0];
  double b = parameters[0][1];
  double C = parameters[0][2];

  double fun = C*pow(x,a-1.)*pow(1-x,b-1.);

  //polynomial in bernstein basis
  double pol = 0.;
  for (int n = 0; n < N; n++)
    pol += parameters[2][n] * bber[i*(N+1)+n+1]; //skip first bernstein polynomial (N,0)

  //log polynomial in power basis l(x) = ( 1 + l0*x + l1*x^2 + ... )
  double lx = logx[i]; //log(x)
  double lp = 1.;

  for (int m = 0; m < M; m++)
    lp += parameters[1][m] * pow(lx,m+1);
    
  double ff = fun*pol*lp;
  //double ff = fun*pol;

  //weight function
  //double w = 1.; //absolute error
  double w = 1./y; //relative error
    
  residuals[0] = (y - ff)*w;
    
  if (!jacobians) return true;

  double* jacobian = jacobians[0];
  if (jacobian)
    {
      jacobian[0] = - lx*ff*w;
      jacobian[1] = - log1mx[i]*ff*w;
      jacobian[2] = - ff/C*w;
    }

  jacobian = jacobians[1];
  if (jacobian)
    for (int m = 0; m < M; m++)
      jacobian[m] =  - fun*pol*pow(lx,m+1)*w;

  jacobian = jacobians[2];
  if (jacobian)
    for (int n = 0; n < N; n++)
      jacobian[n] =  - fun*bber[i*(N+1)+n+1]*lp*w;

  //cout << "residual " << i << "  " << residuals[0] << " y " << y << " ff " << ff << " w " << w << endl;
  return true;
}

void cerespdf::init(double xmin)
{
  //PDF interpolation points
  const LHAPDF::GridPDF& gridpdf = * dynamic_cast<const LHAPDF::GridPDF*>(pdf::lhapdf);
  std::vector <double> xknots = gridpdf.xKnots();

  //Drop knots below xmin to fit PDFs only in the region [xmin,1]
  int k = 0;
  for (int i = 0; i < xknots.size()-1; i++)
    if (xknots[i] < xmin)
      k = i;
    else
      break;
  std::vector<double>::iterator it = xknots.begin();
  xknots.erase(it, it+k);
  
  //for (int i = 0; i < xknots.size(); i++)
  //std::cout << " x " << xknots[i] << std::endl;

  //Define Bjorken-x interpolation points and cache everything which depends on x
  np = 1;                       //Density of interpolation points with respect to the knots --> make this an option
  num_observations = (xknots.size()-1)*np; //Total number of residuals

  xp     = new double[num_observations];
  logx   = new double[num_observations];
  log1mx = new double[num_observations];
  bber   = new double[num_observations*(N+1)];

  for (int i = 0; i < xknots.size()-1; i++)
    {
      double xmn = xknots[i];
      double xmx = xknots[i+1];
      double ll = log(xmx/xmn);
      for (int k = 0; k < np; k++)
	{
	  double t = xmn*exp(ll*double(k)/double(np));
	  xp[i*np+k] = t;
	  logx[i*np+k] = log(t);
	  log1mx[i*np+k] = log1p(-t);

	  //bernstein basis
	  double ber[N+1];
	  bernstein(N,t,ber);
	  for (int n = 0; n <= N; n++)
	    bber[(i*np+k)*(N+1)+n] = ber[n];
	}
    }

  //double data[num_observations];
  double data[2*num_observations];

  B.allocate();
  C.allocate();
  S.allocate();
  D.allocate();
  U.allocate();
  G.allocate();
  UB.allocate();
  DB.allocate();
  SB.allocate();
  CB.allocate();
  BB.allocate();

  setoptions();

  //cout << psi(0,0.5) << "  " << digamma(0.5) << endl;
}

void cerespdf::update(double muf)
{
//  const LHAPDF::GridPDF& gridpdf = * dynamic_cast<const LHAPDF::GridPDF*>(pdf::lhapdf);
//  std::vector <double> q2knots = gridpdf.q2Knots();
//
//  for (int q = 0; q < q2knots.size(); q++)
//    {
//      updatepdfs(sqrt(q2knots[q]));
//      cout << "Q = " << sqrt(q2knots[q]) << endl;
//      
//      //cout << "b PDF ";  B.fit();  		
//      //cout << "c PDF ";  C.fit();			
//      //cout << "s PDF ";  S.fit();			
//      //cout << "d PDF ";  D.fit();			
//      //cout << "u PDF ";  U.fit();			
//      cout << "g PDF ";  G.fit();			
//      //cout << "ub PDF ";  UB.fit();
//      //cout << "db PDF ";  DB.fit();
//      //cout << "sb PDF ";  SB.fit();
//      //cout << "cb PDF ";  CB.fit();
//      //cout << "bb PDF ";  BB.fit();
//    }

  
  updatepdfs(muf);

  //cout << "b PDF ";  B.fit();  		
  //cout << "c PDF ";  C.fit();			
  //cout << "s PDF ";  S.fit();			
  //cout << "d PDF ";  D.fit();			
  //cout << "u PDF ";  U.fit();			
  //cout << "g PDF ";  G.fit();			
  //cout << "ub PDF "; UB.fit();
  //cout << "db PDF "; DB.fit();
  //cout << "sb PDF "; SB.fit();
  //cout << "cb PDF "; CB.fit();
  //cout << "bb PDF "; BB.fit();

  if (muf > pdf::mb) B.fit(); else B.resetpar();
  if (muf > pdf::mc) C.fit(); else C.resetpar();
  S.fit();			
  D.fit();			
  U.fit();			
  G.fit();			
  UB.fit();
  DB.fit();
  SB.fit();
  if (muf > pdf::mc) CB.fit(); else CB.resetpar();
  if (muf > pdf::mb) BB.fit(); else BB.resetpar();
  
  
  //now refit
  //cout << endl << endl << endl;
  //cout << "b PDF " ;  B.fit();
  //cout << "c PDF " ;  C.fit();
  //cout << "s PDF " ;  S.fit();
  //cout << "d PDF " ;  D.fit();
  //cout << "u PDF " ;  U.fit();
  //cout << "g PDF " ;  G.fit();
  //cout << "ub PDF ";  UB.fit();
  //cout << "db PDF ";  DB.fit();
  //cout << "sb PDF ";  SB.fit();
  //cout << "cb PDF ";  CB.fit();
  //cout << "bb PDF ";  BB.fit();
  
  return;

  //Check approximation
  approxpdf PDF(G);
  for (int i = 0; i < num_observations; i++)
    {
      double t = xp[i];
//double fun = PDF.p[2]*pow(t,PDF.p[0]-1.)*pow(1-t,PDF.p[1]-1.);
//
////polynomial
//double pol = 0.;
//
//double ber[N+1];
////bernstein basis
//bernstein(N,t,ber);
//for (int n = 0; n < N; n++)
//	pol += PDF.c[n] * ber[n+1]; //skip first bernstein polynomial (N,0)
//
////log polynomial
//double lp = 1.;
//
////power basis:  p(x) = ( c0 + c1*x + c2*x^2 + ... )
//for (int m = 1; m <= M; m++)
//	lp += PDF.l[m-1] * pow(log(t),m);
      
      double approx = PDF.approx(t);
      
      std::cout << std::setw(10) << "x "       << std::setw(12) << t
		<< std::setw(10) << " pdf "    << std::setw(12) << PDF.f[i]
		<< std::setw(10) << " approx " << std::setw(12) << approx
		<< std::setw(10) << " error "  << std::setw(12) << (PDF.f[i]-approx)/PDF.f[i]
		<< std::endl;
    }
}

void cerespdf::free()
{
  delete[] xp;
  delete[] logx;
  delete[] log1mx;
  delete[] bber;

  B.free();
  C.free();
  S.free();
  D.free();
  U.free();
  G.free();
  UB.free();
  DB.free();
  SB.free();
  CB.free();
  BB.free();
}

void cerespdf::updatepdfs(double muf)
{
  for (int i = 0; i < num_observations; i++)
    {
      double t = xp[i];
      double fx[13];
      int ih = 1;
      fdist_(ih, t, muf, fx);

      //cout << " updatepdfs " << i << " " << t << "  " << fx[parton::B] << endl;
      
      B.f[i]  = fx[parton::B];
      C.f[i]  = fx[parton::C];
      S.f[i]  = fx[parton::S];
      D.f[i]  = fx[parton::D];
      U.f[i]  = fx[parton::U];
      G.f[i]  = fx[parton::G];
      UB.f[i] = fx[parton::Ub];
      DB.f[i] = fx[parton::Db];
      SB.f[i] = fx[parton::Sb];
      CB.f[i] = fx[parton::Cb];
      BB.f[i] = fx[parton::Bb];
    }
}

void cerespdf::setoptions()
{
  //Ceres options
  options.max_num_iterations = 100000;

  options.minimizer_type = ceres::TRUST_REGION; //ceres::LINE_SEARCH;
  options.linear_solver_type = ceres::DENSE_QR; //ceres::SPARSE_NORMAL_CHOLESKY;
  options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT; //ceres::DOGLEG;
  //options.dogleg_type = ceres::SUBSPACE_DOGLEG;
  //options.use_nonmonotonic_steps = true;
  
  options.minimizer_progress_to_stdout = false;
  options.logging_type = ceres::SILENT;
  //options.max_num_consecutive_invalid_steps = 10;
  if (opts.cubacores > 0)
    options.num_threads = opts.cubacores;

  options.dense_linear_algebra_library_type = ceres::EIGEN; //ceres::LAPACK
}


void cerespdf::approxpdf::fit()
{
  Problem problem;

  //analytic differentiation
  for (int i = 0; i < num_observations; ++i) {
    CostFunction* cost_function = new ResAnalytic(xp[i], f[i], i);

    //DynamicCostFunction* cost_function = new ResAnalytic(xp[i], data[i], i);
    //cost_function->AddParameterBlock(P);
    //cost_function->AddParameterBlock(M);
    //cost_function->AddParameterBlock(N);
    //cost_function->SetNumResiduals(num_observations);

    problem.AddResidualBlock(cost_function, NULL, p,l,c);
  }
  
  Solve(options, &problem, &summary);
  //std::cout << summary.BriefReport() << "\n";
  //std::cout << summary.FullReport() << "\n";

  //cout << "final cost " << summary.final_cost
  //     << " iterations " << summary.num_successful_steps+summary.num_unsuccessful_steps
  //     << " time " << summary.total_time_in_seconds << "s"
  //     << endl;
  
  //print parameters
  //printf("al= %g; be= %g; C= %g;\n", p[0], p[1], p[2]);
  //for (int n = 0; n < min(N,12); n++)
  //  printf("c%d= %g\t", n, c[n]);
  //printf("\n");
  //for (int m = 0; m < M; m++)
  //  printf("l%d= %g\t", m, l[m]);
  //printf("\n");

}

//beta function for complex argument :
complex <double> lngam(complex <double> x)
{
  return (x-0.5)*log(x) - x + 0.91893853 + 1./(12.*x)*(1.-1./(30.*x*x)*(1.-1./(3.5*x*x)*(1.-4./(3.*x*x))));
}


complex <double> cbeta (complex <double> z1, complex <double> z2)
{
  complex <double> sub = complex <double> (0., 0.);
  complex <double> zz1 = z1;
  while (real(zz1) < 15.)
    {
      sub = sub + log ((zz1+z2) / zz1);
      zz1 = zz1 + 1.;
    }
  complex <double> zz2 = z2;
  while (real (zz2) < 15.)
    {
      sub = sub + log ((zz1+zz2) / zz2);
      zz2 = zz2 + 1.;
    }
  complex <double> lg1 = lngam (zz1);
  complex <double> lg2 = lngam (zz2);
  complex <double> lg12 = lngam (zz1 + zz2);
  return exp (lg1 + lg2 - lg12 + sub);
}

//int binom(int n, int k)
//{
//  int i,j;
//
//  if ((k == 0) || (k == n))
//    return 1;
//
//  j = min(k, n - k);
//  int binom = 1;
//  for (i = 0; i < j; i++)
//    binom = (binom * (n - i)) / (i + 1);
//
//  return binom;
//}      

int binom(int n, int k)
{
  int C[n+1][k+1];
  int i, j;
  for (i = 0; i <= n; i++)
    for (j = 0; j <= min(i, k); j++)
      if (j == 0 || j == i) C[i][j] = 1;
      else C[i][j] = C[i-1][j-1] + C[i-1][j];
  
  return C[n][k];
}

extern "C"
{
  void psi0_(fcomplex &zz, fcomplex &res);
  void psi1_(fcomplex &zz, fcomplex &res);
  void psi2_(fcomplex &zz, fcomplex &res);
  void psi3_(fcomplex &zz, fcomplex &res);

  //psi function and derivatives from Pegasus
  fcomplex psi_ (fcomplex &z);
  fcomplex dpsi_ (fcomplex &z, int &m);
}

complex <double> psi(int n, complex <double> z)
{
  //To generalise psi to arbitrary n can use the Hurwitz zeta function from https://github.com/andrewfowlie/thermal_funcs/blob/master/src/zeta.cpp (need also zeta.h and bernoulli.h)
  //and the relation: psi(n,z) = (-1)^(n+1) * (n+1)! * zeta(n+1,z) (see https://gcc.gnu.org/onlinedocs/libstdc++/libstdc++-html-USERS-4.4/a01203.html)
  
  fcomplex zz = fcx(z);
  complex <double> result;

  /************ Hurwitz zeta function *************/
  //fcomplex res;
  //switch (n)
  //  {
  //  case 0: psi0_(zz,res); result = cx(res); break;
  //  case 1: psi1_(zz,res); result = cx(res); break;
  //  case 2: psi2_(zz,res); result = cx(res); break;
  //  case 3: psi3_(zz,res); result = cx(res); break;
  //  default:
  //    {
  //	result = hurwitz_zeta(double(n+1), z);
  //	result *= exp(lngam(double(n+1)));
  //	if ((n+1) % 2 == 1)
  //	  result = -result;
  //    }
  //  }
  /***********************************************/
  //cout << "order " << n << " result " << result << " res " << cx(res) << "  " << exp(lngam(double(n+1))) << "  " << hurwitz_zeta(double(n+1), z) << endl;
  //cout << "order " << n << " result " << result << endl;

  /************ Pegasus *************/
  //can use Pegasus functions instead (check if it faster) --> will save the dependency on thermal_funcs, and on gls and gslblas
  if (n == 0)
    result = cx(psi_(zz));
  else
    result = cx(dpsi_(zz, n));
  /***********************************/

  //cout << "Pegasus " << n << " result " << result << endl;
  //cout << endl;
  
  return result;
}

//See Eq. (8) of http://www.iaeng.org/IJAM/issues_v44/issue_4/IJAM_44_4_06.pdf
complex <double> h(int p, int l, complex <double> x, complex <double> y)
{
  complex <double> res = 0.;
  if (l == 0)
    res = psi(p-1,x)-psi(p-1,y);
  if (l == 1)
    for (int j = 1; j < p; j++)
      res += double(binom(p-1,j))*(psi(p-1-j,x)-psi(p-1-j,y))*(psi(j-1,x)-psi(j-1,y));
  else
    for (int j = 1; j < p; j++)
      res += double(binom(p-1,j))*(psi(p-1-j,x)-psi(p-1-j,y))*h(j,l-1,x,y);
  return res;
}

complex <double> h(int p, complex <double> b, complex <double> c)
{
  complex <double> res = 0.;
  switch (p)
    {
    case 0: res = 1.; break;
    case 1: res = psi(0,b) - psi(0,b+c); break;
    case 2: res = pow(psi(0,b) - psi(0,b+c),2) + psi(1,b) - psi(1,b+c); break;
    case 3: res = pow(psi(0,b) - psi(0,b+c),3) + psi(2,b) - psi(2,b+c) + 3.*(psi(0,b) - psi(0,b+c))*(psi(1,b) - psi(1,b+c)); break;
    default:
      for (int l = 0; l < p; l++)
	res += h(p,l,b,b+c);
    }
  return res;
}

complex <double> cerespdf::approxpdf::mellin(complex <double> s)
{
  //Derivatives of the beta function, see:
  //http://www.iaeng.org/IJAM/issues_v44/issue_4/IJAM_44_4_06.pdf
  //https://arxiv.org/pdf/1902.11125.pdf
  complex <double> res = 0.;
  for (int m = 0; m <= M; m++)
    //for (int m = 0; m <= 0; m++)
    {
      complex <double> logp = 0.;
      for (int n = 1; n <= N; n++)
	{
	  //beta function
	  complex <double> bb = p[0]+s-1.+double(n);
	  complex <double> cc = p[1]+double(N-n);
	  complex <double> bet = cbeta(bb,cc);
	  bet *= double(binom(N,n));
	  bet *= h(m,bb,cc);
	  bet *= c[n-1];
	  logp += bet;
	}
      if (m >= 1)
	logp *= l[m-1];
      res += logp;
    }
  res *= p[2];
  return res;
}

double cerespdf::approxpdf::approx(double x)
{
  double fun = p[2]*pow(x,p[0]-1.)*pow(1-x,p[1]-1.);

  //polynomial
  double pol = 0.;
      
  double ber[N+1];
  //bernstein basis
  bernstein(N,x,ber);
  for (int n = 0; n < N; n++)
    pol += c[n] * ber[n+1]; //skip first bernstein polynomial (N,0)
      
  //log polynomial
  double lp = 1.;

  //power basis:  p(x) = ( c0 + c1*x + c2*x^2 + ... )
  for (int m = 1; m <= M; m++)
    lp += l[m-1] * pow(log(x),m);
      
  return fun*pol*lp;
}

complex <double> cerespdf::approxpdf::mellin_num(complex <double> s)
{
  complex <double> res = 0.;

  double xmin = 1e-8;
  double xmax = 1;
  int rule = 500;
  
  double cx = 0.5;
  double mx = 0.5;
  double ll = log(xmax/xmin);
  for (int i = 0; i < rule; i++)
    {
      double x = cx+mx*gr::xxx[rule-1][i];
      double t = xmin*pow(xmax/xmin,x);
      double jac = mx*t*ll;
      double f = approx(t);
      res += pow(t,s-1.)*f*jac*gr::www[rule-1][i];
    }
  return res;
}
