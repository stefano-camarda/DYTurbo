#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <ctime>

#include "settings.h"
#include "interface.h"
#include "finintegr.h"
#include "finitemapping.h"
#include "integr.h"

using namespace std;

void yline();
void mline();
void ptline();
void ptavar();
void ptgvar();
void lowintegr(double &res, double &err);
void realintegr(double &res, double &err);
void virtintegr(double &res, double &err);
void ctintegr(double &res, double &err);

int main( int argc , const char * argv[])
{

  //initialise settings
  //opts.init();
  string conf_file = "input/CT10nlo_settings.in";
  if (argc>1) {
      conf_file = argv[1];
  }
  opts.readfromfile(conf_file.c_str());
  opts.initDyresSettings();
  dyinit_();

  //born level variables (6 dimensions)
  double m, qt, y, costh;
  double phicm, phiZ;
  costh = 0.0;  m = 91.1876;  qt = 5.;  y = 0.5;  phicm = 0.0;  phiZ = 0.;

  //variables to be integrated (4 dimensions in total)
  //1 collinear PDF dimension, common to real, virtual and lowint
  double zcth;
  zcth = 0.5;   //polar angle of Z in the CM frame
  //3 dimensions of real radiation
  double mjj, phijj, costhjj;
  mjj = 10.;    //invariant mass of the dijet (set to 0 to have the correct virtual phase space mapping)
  phijj = 0.3;  //phi of the dijet in dijet rest frame
  costhjj = 0.1;//costh of the dijet in dijet rest frame
  //1 dimension of virtual
  double vz;
  vz = sqrt(0.95);
  //2 dimensions for the counterterm
  double alpha,beta;
  beta = 0.1;
  alpha = 0.1;

  
  //call function wrappers, which map the variables into the unity hypercube of the vegas integration
  cout << " check phase space mapping " << endl;
  dyreal(m, y, qt, phicm, phiZ, costh, zcth, mjj, costhjj, phijj);
  dyvirt(m, y, qt, phicm, phiZ, costh, zcth, vz);
  dyct(m, y, qt, phicm, phiZ, costh, alpha, beta);

  //work in progress for rewriting the counterterm with the same integration structure as the resummed part
  /*
  double cthmom0 = 0;
  double cthmom1 = 0;
  double cthmom2 = 0;
  cout << countterm_(costh,m,qt,y,alpha,beta,cthmom0,cthmom1,cthmom2) << endl;
  */

  //  return 0;
  //lines
  /*
  //mass line
  costh = 0.;
  m = 91;
  qt = 5;
  y = 0.0;
  zcth = 0.5;
  mjj = 0.1;
  costhjj = 0.1;
  phijj = 0.1;
  double mm1 = 20;
  double mm2 = 200;
  int nm = 2000;
  ofstream mf("mline.C");
  mf << "{" << endl;
  mf << "TGraph *gm = new TGraph();" << endl;
  double hm=(mm2-mm1)/nm;
  for(int i=0;i<=nm;i++)
    {
      double m = i*hm+mm1;
      mf << "gm->SetPoint(gm->GetN(), " << i*hm+mm1 << ", " << dyreal(m, y, qt, 0., 0., costh, zcth, mjj, costhjj, phijj) << ");" << endl;
    }
  mf << "gm->Draw();" << endl;
  mf << "}" << endl;

  //pt line
  costh = 0.;
  m = 91;
  qt = 5;
  y = 0.0;
  zcth = 0.;
  mjj = 0.1;
  costhjj = 0.1;
  phijj = 0.1;
  double p1 = 0.0;
  double p2 = 100;
  int np = 2000;
  ofstream pf("ptline.C");
  pf << "{" << endl;
  pf << "TGraph *gp = new TGraph();" << endl;
  double hp=(p2-p1)/np;
  for(int i=0;i<=np;i++)
    {
      double qt = i*hp+p1;
      pf << "gp->SetPoint(gp->GetN(), " << i*hp+p1 << ", " << dyreal(m, y, qt, 0., 0., costh, zcth, mjj, costhjj, phijj) << ");" << endl;
    }
  pf << "gp->Draw();" << endl;
  pf << "}" << endl;

  //mjj line
  costh = 0.;
  m = 91;
  qt = 5;
  y = 0.0;
  zcth = 0.;
  mjj = 0.1;
  costhjj = 0.;
  phijj = 0.;
  double mj1 = 0.0;
  double mj2 = 7000;
  int nmj = 2000;
  ofstream mjf("mjjline.C");
  mjf << "{" << endl;
  mjf << "TGraph *gmj = new TGraph();" << endl;
  double hmj=(mj2-mj1)/nmj;
  for(int i=0;i<=nmj;i++)
    {
      double mjj = i*hmj+mj1;
      mjf << "gmj->SetPoint(gmj->GetN(), " << i*hmj+mj1 << ", " << dyreal(m, y, qt, 0., 0., costh, zcth, mjj, costhjj, phijj) << ");" << endl;
    }
  mjf << "gmj->Draw();" << endl;
  mjf << "}" << endl;

  //cjj line
  costh = 0.;
  m = 91;
  qt = 5;
  y = 0.0;
  zcth = 0.;
  mjj = 0.1;
  costhjj = 0.;
  phijj = 0.;
  double cj1 = 0.0;
  double cj2 = 1.0;
  int ncj = 2000;
  ofstream cjf("costhjjline.C");
  cjf << "{" << endl;
  cjf << "TGraph *gcj = new TGraph();" << endl;
  double hcj=(cj2-cj1)/ncj;
  for(int i=0;i<=ncj;i++)
    {
      double costhjj = i*hcj+cj1;
      cjf << "gcj->SetPoint(gcj->GetN(), " << i*hcj+cj1 << ", " << dyreal(m, y, qt, 0., 0., costh, zcth, mjj, costhjj, phijj) << ");" << endl;
    }
  cjf << "gcj->Draw();" << endl;
  cjf << "}" << endl;
  
  //mass line
  costh = 0.;
  m = 91;
  qt = 5;
  y = 0.0;
  beta = 0.1; //2 dimensions for the counterterm
  alpha = 0.1;
  double mm1 = 20;
  double mm2 = 200;
  int nm = 2000;
  ofstream mf("mline.C");
  mf << "{" << endl;
  mf << "TGraph *gm = new TGraph();" << endl;
  double hm=(mm2-mm1)/nm;
  for(int i=0;i<=nm;i++)
    {
      double m = i*hm+mm1;
      mf << "gm->SetPoint(gm->GetN(), " << i*hm+mm1 << ", " << dyct(m, y, qt, 0., 0., costh, alpha, beta) << ");" << endl;
    }
  mf << "gm->Draw();" << endl;
  mf << "}" << endl;

  //pt line
  costh = 0.;
  m = 91;
  qt = 5;
  y = 0.0;
  beta = 0.5; //2 dimensions for the counterterm
  alpha = 0.5;
  double p1 = 0.0;
  double p2 = 100;
  int np = 1998;
  ofstream pf("ptline.C");
  pf << "{" << endl;
  pf << "TGraph *gp = new TGraph();" << endl;
  double hp=(p2-p1)/np;
  for(int i=0;i<=np;i++)
    {
      double qt = i*hp+p1;
      pf << "gp->SetPoint(gp->GetN(), " << i*hp+p1 << ", " << dyct(m, y, qt, 0., 0., costh, alpha, beta) << ");" << endl;
    }
  pf << "gp->Draw();" << endl;
  pf << "}" << endl;

  //y line
  costh = 0.;
  m = 91;
  qt = 5.0;
  y = 0.0;
  beta = 0.5; //2 dimensions for the counterterm
  alpha = 0.5;
  double y1 = 0.0;
  double y2 = 5.;
  int ny = 2000;
  ofstream yf("yline.C");
  yf << "{" << endl;
  yf << "TGraph *gy = new TGraph();" << endl;
  double hy=(y2-y1)/ny;
  for(int i=0;i<=ny;i++)
    {
      double y = i*hy+y1;
      yf << "gy->SetPoint(gy->GetN(), " << i*hy+y1 << ", " << dyct(m, y, qt, 0., 0., costh, alpha, beta) << ");" << endl;
    }
  yf << "gy->Draw();" << endl;
  yf << "}" << endl;

  //costh line
  costh = 0.;
  m = 91;
  qt = 5.0;
  y = 0.0;
  beta = 0.5; //2 dimensions for the counterterm
  alpha = 0.5;
  double c1 = -1.0;
  double c2 = 1.;
  int nc = 2000;
  ofstream cf("cline.C");
  cf << "{" << endl;
  cf << "TGraph *gc = new TGraph();" << endl;
  double hc=(c2-c1)/nc;
  for(int i=0;i<=nc;i++)
    {
      double costh = i*hc+c1;
      cf << "gc->SetPoint(gc->GetN(), " << i*hc+c1 << ", " << dyct(m, y, qt, 0., 0., costh, alpha, beta) << ");" << endl;
    }
  cf << "gc->Draw();" << endl;
  cf << "}" << endl;

  //alpha line
  costh = 0.;
  m = 91;
  qt = 5;
  y = 0.0;
  beta = 0.; //2 dimensions for the counterterm
  alpha = 0.;
  double a1 = 0.0;
  double a2 = 1.0;
  int na = 2000;
  ofstream af("aline.C");
  af << "{" << endl;
  af << "TGraph *ga = new TGraph();" << endl;
  double ha=(a2-a1)/na;
  for(int i=0;i<=na;i++)
    {
      double alpha = i*ha+a1;
      af << "ga->SetPoint(ga->GetN(), " << i*ha+a1 << ", " << dyct(m, y, qt, 0., 0., costh, alpha, beta) << ");" << endl;
    }
  af << "ga->Draw();" << endl;
  af << "}" << endl;

  //beta line
  costh = 0.;
  m = 91;
  qt = 5;
  y = 0.0;
  beta = 0.; //2 dimensions for the counterterm
  alpha = 0.;
  double b1 = 0.0;
  double b2 = 1.0;
  int nb = 2000;
  ofstream bf("bline.C");
  bf << "{" << endl;
  bf << "TGraph *gb = new TGraph();" << endl;
  double hb=(b2-b1)/nb;
  for(int i=0;i<=nb;i++)
    {
      double beta = i*hb+b1;
      bf << "gb->SetPoint(gb->GetN(), " << i*hb+b1 << ", " << dyct(m, y, qt, 0., 0., costh, alpha, beta) << ");" << endl;
    }
  bf << "gb->Draw();" << endl;
  bf << "}" << endl;
  */



  bins.readfromfile(conf_file.c_str());
  clock_t begin_time, end_time;
  double value, error;

  cout << endl;
  cout << "Start integration of counterterm" << endl;
  begin_time = clock();
  for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
    {
      //Set integration boundaries
      setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
      clock_t b_time = clock();
      ctintegr(value, error);
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

  if (opts.order == 1)
    {
      cout << "Start integration of Z+j LO" << endl;
      begin_time = clock();
      for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
	{
	  //Set integration boundaries
	  setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
	  clock_t b_time = clock();
	  lowintegr(value, error);
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
    }
  if (opts.order == 2)
    {
      cout << endl;
      cout << "Start integration of real" << endl;
      begin_time = clock();
      for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
	{
	  //Set integration boundaries
	  setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
	  clock_t b_time = clock();
	  realintegr(value, error);
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

      cout << endl;
      cout << "Start integration of virtual" << endl;
      begin_time = clock();
      for (vector<double>::iterator qit = bins.qtbins.begin(); qit != bins.qtbins.end()-1; qit++)
	{
	  //Set integration boundaries
	  setbounds(opts.mlow, opts.mhigh, *qit, *(qit+1), opts.ylow, opts.yhigh);
	  clock_t b_time = clock();
	  virtintegr(value, error);
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
    }

  return 0;
}

//Cuba integration of Z+j LO
void lowintegr(double &res, double &err)
{
  const int ndim = 7;   //dimensions of the integral
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
  const int flags = 8+4+opts.cubaverbosity;
  const int seed = 1;
  const int mineval = 10000000;
  const int maxeval = 10000000;
  const int nstart = 100000;
  const int nincrease = 100000;
  const int nbatch = 1000;
  const int gridno = 1;
  Vegas(ndim, ncomp, (integrand_t)lowintegrand, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res = integral[0];
  err = error[0];

  return;
}

//Cuba integration of the real part
void realintegr(double &res, double &err)
{
  const int ndim = 10;   //dimensions of the integral
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
  const int flags = 8+4+opts.cubaverbosity;
  const int seed = 1;
  const int mineval = 10000000;
  const int maxeval = 10000000;
  const int nstart = 100000;
  const int nincrease = 100000;
  const int nbatch = 1000;
  const int gridno = 1;
  Vegas(ndim, ncomp, (integrand_t)realintegrand, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res = integral[0];
  err = error[0];

  return;
}

//Cuba integration of the virtual part
void virtintegr(double &res, double &err)
{
  const int ndim = 8;   //dimensions of the integral
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
  const int flags = 8+4+opts.cubaverbosity;
  const int seed = 1;
  const int mineval = 10000000;
  const int maxeval = 10000000;
  const int nstart = 100000;
  const int nincrease = 100000;
  const int nbatch = 1000;
  const int gridno = 1;
  Vegas(ndim, ncomp, (integrand_t)virtintegrand, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res = integral[0];
  err = error[0];

  return;
}

//Cuba integration of the counterterm
void ctintegr(double &res, double &err)
{
  const int ndim = 8;   //dimensions of the integral
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
  const int flags = 8+4+opts.cubaverbosity;
  const int seed = 1;
  const int mineval = 1000000;
  const int maxeval = 1000000;
  const int nstart = 100000;
  const int nincrease = 100000;
  const int nbatch = 1000;
  const int gridno = 1;
  Vegas(ndim, ncomp, (integrand_t)ctintegrand, userdata, nvec,
	epsrel, epsabs,
	flags, seed,
	mineval, maxeval,
	nstart, nincrease, nbatch,
	gridno, statefile, spin,
	&neval, &fail,
	integral, error, prob);
  res = integral[0];
  err = error[0];

  return;
}
