#include <iostream>
#include <LHAPDF/LHAPDF.h>
//#include "DEIntegrator.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <ctime>
#include <cuba.h>
#include <iomanip>

#include "integr.h"
#include "settings.h"
#include "interface.h"

using namespace std;

void yline();
void mline();
void costhline();
void ptline();
void ptavar();
void ptgvar();

void integr2d(double &res, double &err);
void integr3d(double &res, double &err);
void integr4d(double &res, double &err);

int main()
{
  cout << "hello world" << endl;
  clock_t begin_time, end_time;

  //  lhapdfs_.lhapdfs_ = true;
  //  LHAPDF::initPDFSet("CT10nnlo");
  //  double zmass = 91.1876;
  //  couple_.amz_ = LHAPDF::alphasPDF(zmass);
  //  nlooprun_.nlooprun_ = 2;
  //  scale_.scale_ = zmass;
  //  setup_();
  dyinit_();
  //  cout << a_param_.a_param_ << endl;
  //  std::cout << std::setprecision(15);

  //initialise settings
  opts.init();
  
  double costh, m, qt, y;
  double value, error;
  int mode = 0;
  
  costh = 0.3; m = 91; qt = 1; y = 0;
  begin_time = clock();
  value = resumm_(costh,m,qt,y,mode);
  end_time = clock();
  cout << setw(10) << "Result" << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  costh = 0.1; m = 91; qt = 5; y = 0.2;
  begin_time = clock();
  value = resumm_(costh,m,qt,y,mode);
  end_time = clock();
  cout << setw(10) << "Result" << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  costh = 0.5; m = 70; qt = 10; y = 1.5;
  begin_time = clock();
  value = resumm_(costh,m,qt,y,mode);
  end_time = clock();
  cout << setw(10) << "Result" << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  costh = -1.0; m = 110; qt = 20; y = -2.5;
  begin_time = clock();
  value = resumm_(costh,m,qt,y,mode);
  end_time = clock();
  cout << setw(10) << "Result" << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  costh = 0.1; m = 91; qt = 5; y = 1.5;
  begin_time = clock();
  value = resumm_(costh,m,qt,y,mode);
  end_time = clock();
  cout << setw(10) << "Result" << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  costh = 0.1; m = 91; qt = 5; y = 0.2;
  begin_time = clock();
  value = resumm_(costh,m,qt,y,mode);
  end_time = clock();
  cout << setw(10) << "Result" << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  costh = 0.59172160171850774; m = 80.887838188405482;
  qt = 4.0234498520010531; y = 1.0630039379722271;
  begin_time = clock();
  value = resumm_(costh,m,qt,y,mode);
  end_time = clock();
  cout << setw(10) << "Result" << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  costh = 0.76061949271840290; m = 91.557029972501027;
  qt = 4.6388123004847825; y = 0.73501935633317217;
  begin_time = clock();
  value = resumm_(costh,m,qt,y,mode);
  end_time = clock();
  cout << setw(10) << "Result" << setw(15) << value
       << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  
  //mline
  int nocuts = (int)true;
  double m1 = 0;
  double m2 = 1;
  int nm = 1;
  setbounds(opts.mlow, opts.mhigh, 40, 60, opts.ylow, opts.yhigh);

  ofstream mf("mline.C");
  mf << "{" << endl;
  mf << "TGraph *gm = new TGraph();" << endl;
  double hm=(m2-m1)/nm;
  for(int i=0;i<=nm;i++)
    {
      double m = i*hm+m1;
      double ymin = 0.;
      double ymax = 1.;
      //      rapintegrals_(ymin,ymax,m,nocuts);
      double x[3];
      double f[1];
      const int ncomp = 1;
      const int ndim = 3;
      x[0] = m;
      x[1] = 0.5;
      x[2] = 0.5;
      //      resintegrand2d(ndim, x, ncomp, f);
      resintegrand3d(ndim, x, ncomp, f);
      mf << "gm->SetPoint(gm->GetN(), " << i*hm+m1 << ", " << f[0] << ");" << endl;
    }
  mf << "gm->Draw();" << endl;
  mf << "}" << endl;
  
  //ptline();
  //yline();
  //mline();
  //ptavar();
  //ptgvar();

  //return 0;
  // Cuba integration
  bins.init();

  cout << endl;
  cout << "Start integration" << endl;
  begin_time = clock();
  for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
    {
      //Set integration boundaries
      setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
      clock_t b_time = clock();
      if (opts.int2d) integr2d(value, error);
      if (opts.int3d) integr3d(value, error);
      if (opts.int4d) integr4d(value, error);
      clock_t e_time = clock();
      value = value / (*(qit+1) - *qit);
      error = error / (*(qit+1) - *qit);
      cout << setw(3) << "bin" << setw(5) << *qit << setw(2) << "-" << setw(5) << *(qit+1)
	   << setw(10) << value << setw(5) << "+/-" << setw(10) << error
	   << setw(6) << "time" << setw(10) <<  float(e_time - b_time) / CLOCKS_PER_SEC << endl;
    }
  end_time = clock();
  cout << endl;
  cout << setw(10) << "time "  << setw(15) << float(end_time - begin_time) / CLOCKS_PER_SEC << endl;

  return 0;
  /*
  //DEquadrature
  FunctResm fm;
  int evaluations;
  double errorEstimate;
  std::cout << std::setprecision(15);

  std::cout << "Start integration" << std::endl;
    
  result = DEIntegrator<FunctResm>::Integrate(fm, a, b, 10, evaluations, err);
  std::cout << integral << ", " << errorEstimate << ", " << evaluations << "\n";
  */

  /*
  //GSL
  size_t limit = 1000;
  gsl_integration_workspace * gslw = gsl_integration_workspace_alloc (limit);
  gsl_function gslfm;
  gslfm.function = &resm;
  gsl_integration_qag(&gslfm, a, b, 0, 1E-5, limit, 1, gslw, &result, &err);
  //size_t neval;
  //gsl_integration_qng (&gslfm, a, b, 0, 0.1, &result, &err, &neval);
  */

  /*
  //trapezoidal rule
  int n = 30;
  double h,r=0;
  h=(b-a)/n;
  for(int i=1;i<n;i++)
    r=(2*resm(i*h+a, NULL))+r;
  result=(h/2)*(resm(a, NULL)+resm(b, NULL)+r);
  */


  
  /*
  begin_time = clock();
  //Integrate
  qt = 5;
  set(costh,m,qt,y);
  double a = 66;
  double b = 116;
  double result = 0;
  double  err = 0;
  double pta = 0;
  double ptb = 30;

  int n = 120;
  double h,r=0;
  h=(ptb-pta)/n;
  for(int i=1;i<=n;i++)
    {
      setqt(i*h+pta);
      //cos theta integral
      //GSL
      a = -1;
      b = +1;
      size_t limit = 1000;
      gsl_integration_workspace * gslw = gsl_integration_workspace_alloc (limit);
      gsl_function gslfth;
      gslfth.function = &resth;
      begin_time = clock();
      gsl_integration_qag(&gslfth, a, b, 0, 1E-2, limit, 1, gslw, &result, &err);
      //size_t neval;
      //gsl_integration_qng (&gslfth, a, b, 0, 1E-3, &result, &err, &neval);
      end_time = clock();
      cout << "pt " << i*h+pta << "  " << result << "  " << err << endl;
      cout << "time " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << endl;
    }
  */

  //rapidity integral

  /*
  //GSL
  a = 0;
  b = 5;
  size_t limit = 1000;
  gsl_integration_workspace * gslw = gsl_integration_workspace_alloc (limit);
  gsl_function gslfy;
  gslfy.function = &resy;
  gsl_integration_qag(&gslfy, a, b, 0, 1E-1, limit, 1, gslw, &result, &err);
  end_time = clock();
  cout << result << "  " << err << endl;
  cout << "time " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << endl;
  */
}

void integr2d(double &res, double &err)
{
  const int ndim = 2;     //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[1];
  double error[1];
  double prob[1];
  const int flags = 0+opts.cubaverbosity;
  const int mineval = 65+2*65*opts.niter;
  const int maxeval = 65+2*65*opts.niter;
  const int key = 13;
  int nregions;
  Cuhre(ndim, ncomp,
	(integrand_t) resintegrand2d, userdata, nvec,
	epsrel, epsabs,
	flags,
	mineval, maxeval,
	key, statefile, NULL,
	&nregions, &neval, &fail,
  	integral, error, prob);

  res = integral[0];
  err = error[0];
  return;
}

void integr3d(double &res, double &err)
{
  const int ndim = 3;     //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[1];
  double error[1];
  double prob[1];
  const int flags = 0+opts.cubaverbosity;
  const int mineval = 127+2*127*opts.niter;
  const int maxeval = 127+2*127*opts.niter;
  const int key = 13;
  int nregions;
  Cuhre(ndim, ncomp,
	(integrand_t) resintegrand3d, userdata, nvec,
	epsrel, epsabs,
	flags,
	mineval, maxeval,
	key, statefile, NULL,
	&nregions, &neval, &fail,
  	integral, error, prob);

  res = integral[0];
  err = error[0];
  return;
}

//original integration, can use this to sample phase space and fill histograms
void integr4d(double &res, double &err)
{
  const int ndim = 4;   //dimensions of the integral
  const int ncomp = 1;  //components of the integrand
  void *userdata;
  const int nvec = 1;
  const double epsrel = 0.;
  const double epsabs = 0.;
  const char *statefile = "";
  void *spin=NULL;
  int neval;
  int fail;
  double integral[1];
  double error[1];
  double prob[1];
  const int flags = 4+opts.cubaverbosity;
  const int seed = 1;
  const int mineval = opts.vegasncalls;
  const int maxeval = opts.vegasncalls;
  const int nstart = 1000;
  const int nincrease = 1000;
  const int nbatch = 10000;
  const int gridno = 1;
  Vegas(ndim, ncomp, (integrand_t)resintegrand4d, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, NULL,
	&neval, &fail,
	integral, error, prob);
  res = integral[0];
  err = error[0];
  return;
}

//y line plot
void yline()
{
  double costh = 0.;
  double m = 80.385;
  double qt = 5.;
  double y = 0.0;
  int mode = 0;

  double y1 = -5;
  double y2 = 5;
  int ny = 200;

  ofstream yf("yline.C");
  yf << "{" << endl;
  yf << "TGraph *gy = new TGraph();" << endl;

  double hy=(y2-y1)/ny;
  for(int i=0;i<=ny;i++)
    {
      double y = i*hy+y1;
      set(m, qt, y);//set global variables to m, qt, y
      genV4p(m, qt, y, 0.);//generate boson 4-momentum, with m, qt, y and phi=0
      yf << "gy->SetPoint(gy->GetN(), " << i*hy+y1 << ", " << resumm_(costh,m,qt,y,mode) << ");" << endl;
    }
  yf << "gy->Draw();" << endl;
  yf << "}" << endl;
}

//m line plot
void mline()
{
  double costh = 0.1;
  double m = 80.385;
  double qt = 5;
  double y = 1.0;
  int mode = 0;

  double m1 = 10;
  double m2 = 200;
  int nm = 100;

  ofstream mf("mline.C");
  mf << "{" << endl;
  mf << "TGraph *gm = new TGraph();" << endl;
  double hm=(m2-m1)/nm;
  for(int i=0;i<=nm;i++)
    {
      double m = i*hm+m1;
      set(m, qt, y);//set global variables to m, qt, y
      genV4p(m, qt, y, 0.);//generate boson 4-momentum, with m, qt, y and phi=0
      mf << "gm->SetPoint(gm->GetN(), " << i*hm+m1 << ", " << resumm_(costh,m,qt,y,mode) << ");" << endl;
    }
  mf << "gm->Draw();" << endl;
  mf << "}" << endl;
}

/*
//costh line plot
void costhline()
{
  double costh = 0.;
  double m = 91;
  double qt = 5;
  double y = 0.0;

  double costh_CS;

  double c1 = -1;
  double c2 = +1;
  int nc = 2000;

  ofstream cf("cline.C");
  cf << "{" << endl;
  cf << "TGraph *gc = new TGraph();" << endl;
  double hc=0;
  hc=(c2-c1)/nc;
  for(int i=0;i<=nc;i++)
    {
      double costh_CS = i*hc+c1;
      if (cuts(costh, m, qt, y, 0., 0., costh_CS))
	cf << "gc->SetPoint(gc->GetN(), " << i*hc+c1 << ", " << resumm_(costh_CS,m,qt,y) << ");" << endl;
      else
	cf << "gc->SetPoint(gc->GetN(), " << i*hc+c1 << ", " << 0 << ");" << endl;
    }
  cf << "gc->Draw();" << endl;
  cf << "}" << endl;
}  
*/
//pt line plot
void ptline()
{
  double costh = 0.;
  double m = 80.385;
  double qt = 5.;
  double y = 0.0;
  int mode = 0;

  double p1 = 0.1;
  double p2 = 100;
  int np = 1998;

  ofstream pf("ptline.C");
  pf << "{" << endl;
  pf << "TGraph *gp = new TGraph();" << endl;
  double hp=(p2-p1)/np;
  for(int i=0;i<=np;i++)
    {
      double qt = i*hp+p1;
      set(m, qt, y);//set global variables to m, qt, y
      genV4p(m, qt, y, 0.);//generate boson 4-momentum, with m, qt, y and phi=0
      pf << "gp->SetPoint(gp->GetN(), " << i*hp+p1 << ", " << resumm_(costh,m,qt,y,mode) << ");" << endl;
    }
  pf << "gp->Draw();" << endl;
  pf << "}" << endl;
}
