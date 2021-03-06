#include "dyres_interface.h"
#include "mcfm_interface.h"
#include "propagator.h"
#include "mesq.h"
#include "vjint.h"
#include "vjloint.h"

//y line plot
void yline()
{
  double costh = 0.;
  double m = opts.rmass;
  double qt = (phasespace::qtmin + phasespace::qtmax)/2.; //0.01;//0.850673; //0.01;
  double y = 0.;
  int mode = 1;
  double f[opts.totpdf];

  double y1 = phasespace::ymin;
  double y2 = phasespace::ymax;
  int ny = 100;

  ofstream yf("yline.C");
  yf << std::setprecision(15);
  yf << "{" << endl;
  yf << "TGraph *gy = new TGraph();" << endl;
  yf << "TGraph *gy1 = new TGraph();" << endl;
  yf << "TGraph *gy2 = new TGraph();" << endl;

  double hy=(y2-y1)/ny;
  for(int i=0;i<=ny;i++)
    {
      double y = i*hy+y1;
      double ym = -y;
      phasespace::set_mqtyphi(m, qt, y);//set global variables to costh, m, qt, y
      phasespace::set_cth(costh);//set global variables to costh, m, qt, y
      phasespace::calcexpy();
      phasespace::calcmt();
      omegaintegr::genV4p();//generate boson 4-momentum, with m, qt, y and phi=0
      //yf << "gy->SetPoint(gy->GetN(), " << i*hy+y1 << ", " << resumm_(costh,m,qt,y,mode) << ");" << endl;
      //      yf << "gy->SetPoint(gy->GetN(), " << i*hy+y1 << ", " << resint::rint(costh,m,qt,y,mode) << ");" << endl;
      //      if (vjfo_(m,qt,y) != 0)
      //yf << "gy->SetPoint(gy->GetN(), " << i*hy+y1 << ", " << (-ctint_(costh,m,qt,y,mode,f)*2*qt)/vjfo_(m,qt,y) << ");" << endl;
      //      yf << "gy1->SetPoint(gy1->GetN(), " << i*hy+y1 << ", " << -ctint_(costh,m,qt,y,mode,f)*2*qt << ");" << endl;
      //      yf << "gy2->SetPoint(gy2->GetN(), " << i*hy+y1 << ", " << vjfo_(m,qt,y) << ");" << endl;
      // check of y asymmetry
      // yf << "gy1->SetPoint(gy1->GetN(), " << i*hy+y1 << ", " << -ctint_(costh,m,qt,y,mode,f)*2*qt+ctint_(costh,m,qt,ym,mode,f)*2*qt << ");" << endl;
      // yf << "gy2->SetPoint(gy2->GetN(), " << i*hy+y1 << ", " << vjfo_(m,qt,y)-vjfo_(m,qt,ym) << ");" << endl;

      //      double vj = vjint::vint(m,qt,y);
      
      dofill_.doFill_ = 1;
      //qtint::calc(m,qt,0,1);
      //omegaintegr::genV4p();//generate boson 4-momentum, with m, qt, y and phi=0
      mode = 2;
      qtint::calc(m,max(opts.qtcut,opts.xqtcut*phasespace::m),1e10,mode);
      ctint::calc(costh,m,qt,y,mode,f);

      //      cout << vj << "  " << f[0] << endl;

      //yf << "gy->SetPoint(gy->GetN(), " << i*hy+y1 << ", " << (-f[0]*2*qt)/vj << ");" << endl;
      yf << "gy1->SetPoint(gy1->GetN(), " << i*hy+y1 << ", " << f[0]*2*qt << ");" << endl;
      //yf << "gy->SetPoint(gy->GetN(), " << i*hy+y1 << ", " << vjint::vint(m,qt,y) << ");" << endl;
      
    }
  yf << "//gy->Draw();" << endl;
  yf << "gy1->Draw();" << endl;
  yf << "//gy2->Draw(\"same\");" << endl;
  yf << "}" << endl;
}

//m line plot
void mline()
{
  double costh = 0.;
  double m = opts.rmass;
  double qt = (phasespace::qtmin + phasespace::qtmax)/2.;
  double y =  (phasespace::ymin + phasespace::ymax)/2.;
  int mode = 4; //3;
  double f[opts.totpdf];

  double m1 = phasespace::mmin;
  double m2 =  phasespace::mmax;
  int nm = 150;

  ofstream mf("mline.C");
  mf << std::setprecision(15);
  mf << "{" << endl;
  mf << "TGraph *gm = new TGraph();" << endl;
  mf << "TGraph *gm1 = new TGraph();" << endl;
  mf << "TGraph *gm2 = new TGraph();" << endl;
  double hm=(m2-m1)/nm;
  for(int i=0;i<=nm;i++)
    {
      double m = i*hm+m1;
      phasespace::set_mqtyphi(m, qt, y);//set global variables to costh, m, qt, y
      phasespace::set_cth(costh);//set global variables to costh, m, qt, y
      phasespace::calcexpy();
      phasespace::calcmt();
      omegaintegr::genV4p();//generate boson 4-momentum, with m, qt, y and phi=0
      //      double vj = vjint::vint(m,qt,y);
      //      mf << "gm->SetPoint(gm->GetN(), " << i*hm+m1 << ", " << resumm_(costh,m,qt,y,mode) << ");" << endl;

//      rapint::cache(phasespace::ymin, phasespace::ymax);
//      rapint::allocate();
//      rapint::integrate(phasespace::ymin,phasespace::ymax,(phasespace::mmin+phasespace::mmax)/2.);
      resint::rint(costh,m,qt,y,mode,f);
      mf << "gm->SetPoint(gm->GetN(), " << i*hm+m1 << ", " << f[0] << ");" << endl;
      //mf << "gm->SetPoint(gm->GetN(), " << i*hm+m1 << ", " << -ctint_(costh,m,qt,y,mode,f)*2*qt/vjint::vint(m,qt,y) << ");" << endl;
      //mf << "gm1->SetPoint(gm1->GetN(), " << i*hm+m1 << ", " << -ctint_(costh,m,qt,y,mode,f)*2*qt << ");" << endl;
      //mf << "gm2->SetPoint(gm2->GetN(), " << i*hm+m1 << ", " << vjint::vint(m,qt,y) << ");" << endl;

      dofill_.doFill_ = 1;

      //mode = 2;
      //qtint::calc(m,max(opts.qtcut,opts.xqtcut*phasespace::m),1e10,mode);
      //omegaintegr::genV4p();//generate boson 4-momentum, with m, qt, y and phi=0
      //ctint::calc(costh,m,qt,y,mode,f);
      //mf << "gm->SetPoint(gm->GetN(), " << i*hm+m1 << ", " << f[0]*2*qt << ");" << endl;

      //      cout << vj << "  " << f[0]*2*qt << endl;
      
      //      mf << "gm->SetPoint(gm->GetN(), " << i*hm+m1 << ", " << (-f[0]*2*qt)/vj << ");" << endl;
    }
  mf << "gm->Draw();" << endl;
  mf << "//gm1->Draw();" << endl;
  mf << "//gm2->Draw(\"same\");" << endl;
  mf << "}" << endl;
}

void mlinebw()
{
  //mline after bret wigner unweighting
  int nocuts = (int)true;
  double m1 = 0;
  double m2 = 1;
  int nm = 100;
  phasespace::setbounds(phasespace::mmin, phasespace::mmax, 0, 2, phasespace::ymin, phasespace::ymax);

  ofstream mf("mlinebw.C");
  mf << "{" << endl;
  mf << "TGraph *gm = new TGraph();" << endl;
  double hm=(m2-m1)/nm;
  for(int i=0;i<=nm;i++)
    {
      double m = i*hm+m1;
      double ymin = 0.;
      double ymax = 1.;
      //      rapintegrals_(ymin,ymax,m,nocuts);
      double x[4];
      double f[1];
      const int ncomp = 1;
      const int ndim = 4; //3; //2;
      x[0] = 0.5;
      x[1] = 0.1;
      x[2] = m;
      x[3] = 0.5;
      //resintegrand2d(ndim, x, ncomp, f);
      //      resintegrand3d(ndim, x, ncomp, f);
      ctintegrand3d(ndim, x, ncomp, f);
      //void* userdata; int nvec; int core; double weight; int iter; resintegrand4d(ndim, x, ncomp, f, userdata, nvec, core, weight, iter);
      mf << "gm->SetPoint(gm->GetN(), " << i*hm+m1 << ", " << f[0] << ");" << endl;
    }
  mf << "gm->Draw();" << endl;
  mf << "}" << endl;
}

void xline()
{
  //xline after unweighting
  int nocuts = (int)true;
  double x1 = 0.;
  double x2 = 1.;
  int nx = 100;
  phasespace::setbounds(phasespace::mmin, phasespace::mmax, phasespace::qtmin, phasespace::qtmax, phasespace::ymin, phasespace::ymax);

  ofstream xf("xline.C");
  xf << "{" << endl;
  xf << "TGraph *gx = new TGraph();" << endl;
  double hx=(x2-x1)/nx;
  for(int i=0;i<=nx;i++)
    {
      double xx = i*hx+x1;
      double ymin = 0.;
      double ymax = 1.;
      //    rapintegrals_(ymin,ymax,m,nocuts);
      const int ncomp = 1;
      const int ndim = 2; //4; //3; //10;
      double x[ndim];
      double f[ncomp],g[ncomp];
      x[0] = 0.5;
      x[1] = xx;
      //x[2] = 0.5;
      //x[3] = 0.5;
      //x[4] = 0.5;
      //x[5] = 0.5;
      //x[6] = 0.5;
      //x[7] = 0.5;
      //x[8] = 0.5;
      //x[9] = 0.5;
      void* userdata = 0;
      const int nvec = 1;
      const int core = 0;
      double weight;
      const int iter = 0;

      //rapint::cache(phasespace::ymin, phasespace::ymax);
      //rapint::allocate();
      //rapint::integrate(phasespace::ymin,phasespace::ymax,(phasespace::mmin+phasespace::mmax)/2.);

      //resintegrand1d(ndim, x, ncomp, f);
      //resintegrand2d(ndim, x, ncomp, f);
      //vjlointegrand(ndim, x, ncomp, f);
      //resintegrand3d(ndim, x, ncomp, f);
      //ctintegrand3d(ndim, x, ncomp, f);
      ctintegrand2d(ndim, x, ncomp, f);
      //                  x[0] = 0.5;
      //                  x[1] = 0.5;
      //                  x[2] = xx;
      //vjintegrand(ndim, x, ncomp, g);
      //lointegrand2d(ndim, x, ncomp, f);
      //lointegrandMC(ndim, x, ncomp, f, userdata, nvec, core, weight, iter);
      //resintegrandMC(ndim, x, ncomp, f, userdata, nvec, core, weight, iter);
      //ctintegrandMC(ndim, x, ncomp, f, userdata, nvec, core, weight, iter);
      //vjlointegrandMC(ndim, x, ncomp, f, userdata, nvec, core, weight, iter);
      //void* userdata; int nvec; int core; double weight; int iter; resintegrand4d(ndim, x, ncomp, f, userdata, nvec, core, weight, iter);
      //realintegrand(ndim, x, ncomp, f, userdata, nvec, core, weight, iter);
      //xf << "gx->SetPoint(gx->GetN(), " << i*hx+x1 << ", " << f[0]/g[0] << ");" << endl;
      //xf << "gx->SetPoint(gx->GetN(), " << i*hx+x1 << ", " << f[0]+g[0] << ");" << endl;
      //xf << "gx->SetPoint(gx->GetN(), " << i*hx+x1 << ", " << g[0] << ");" << endl;
      xf << "gx->SetPoint(gx->GetN(), " << i*hx+x1 << ", " << f[0] << ");" << endl;

    }
  xf << "gx->Draw();" << endl;
  xf << "}" << endl;
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
  //double m = (phasespace::mmin + phasespace::mmax)/2.;
  double m = opts.rmass;
  double qt = 5.;
  double y = 0.0;
  int mode = 1;
  double f[2];
  /*
  phasespace::setbounds(phasespace::mmin, phasespace::mmax, 0, 100, phasespace::ymin, phasespace::ymax);
  cacheyrapint_(phasespace::ymin, phasespace::ymax);
  if (mode >= 2)
    if (opts.resumcpp)
      {
	rapint::cache(phasespace::ymin, phasespace::ymax);
	rapint::allocate();
	rapint::integrate(phasespace::ymin,phasespace::ymax,m);
      }
  */

  double p1 =  0.01;
  double p2 =  1;
  int np = 99;

  ofstream pf("ptline.C");
  pf << std::setprecision(15);
  pf << "{" << endl;
  pf << "TGraph *gp = new TGraph();" << endl;
  pf << "TGraph *gp1 = new TGraph();" << endl;
  pf << "TGraph *gp2 = new TGraph();" << endl;
  double hp=(p2-p1)/np;

  for(int i=0;i<=np;i++)
    {
      double qt = i*hp+p1;
      phasespace::set_mqtyphi(m, qt, y);//set global variables to costh, m, qt, y
      phasespace::set_cth(costh);//set global variables to costh, m, qt, y
      phasespace::calcexpy();
      phasespace::calcmt();
      omegaintegr::genV4p();//generate boson 4-momentum, with m, qt, y and phi=0
      //double vj = vjint::vint(m,qt,y);

      //pf << "gp->SetPoint(gp->GetN(), " << i*hp+p1 << ", " << resumm_(costh,m,qt,y,mode) << ");" << endl;
      //phasespace::setbounds(phasespace::mmin, phasespace::mmax, 0, qt, phasespace::ymin, phasespace::ymax);
      resint::rint(costh,m,qt,y,mode,f);
      pf << "gp->SetPoint(gp->GetN(), " << i*hp+p1 << ", " << f[0] << ");" << endl;
      //pf << "gp->SetPoint(gp->GetN(), " << i*hp+p1 << ", " << vjint::vint(m,qt,y) << ");" << endl;
      /*      pf << "gp->SetPoint(gp->GetN(), "   << i*hp+p1 << ", " << vjfo_(m,qt,y)+ctint_(costh,m,qt,y,mode,f)*2*qt << ");" << endl;
      pf << "gp1->SetPoint(gp1->GetN(), " << i*hp+p1 << ", " << -ctint_(costh,m,qt,y,mode,f)*2*qt << ");" << endl;
      pf << "gp2->SetPoint(gp2->GetN(), " << i*hp+p1 << ", " << vjfo_(m,qt,y) << ");" << endl;*/

      //dofill_.doFill_ = 1;
      //qtint::calc(m,qt,0,1);
      //omegaintegr::genV4p();//generate boson 4-momentum, with m, qt, y and phi=0
      //mode = 1;
      //ctint::calc(costh,m,qt,y,mode,f);

      //cout << qt << "  " << vj << "  " << f[0]*2*qt << "  " << (f[0]*2*qt)+vj << endl;
      //pf << "gp->SetPoint(gp->GetN(), " << i*hp+p1 << ", " << -((f[0]*2*qt)+vj)*qt << ");" << endl;
      //pf << "gp->SetPoint(gp->GetN(), " << i*hp+p1 << ", " << 1.+(f[0]*2*qt)/vj << ");" << endl;
    }
  pf << "gp->Draw();" << endl;
  pf << "//gp1->Draw();" << endl;
  pf << "//gp2->Draw(\"same\");" << endl;
  pf << "}" << endl;
}

void ptomline()
{
  double costh = 0.;
  //double m = (phasespace::mmin + phasespace::mmax)/2.;
  double m = opts.rmass;
  double qt = 5.;
  double y = 0.0;
  int mode = 2;
  double f[opts.totpdf];

  double p1 = 0.01;
  double p2 = 0.5;
  int np = 49;

  ofstream pf("ptomline.C");
  pf << std::setprecision(15);
  pf << "{" << endl;
  pf << "TGraph *gp = new TGraph();" << endl;
  pf << "TGraph *gp1 = new TGraph();" << endl;
  pf << "TGraph *gp2 = new TGraph();" << endl;
  double hp=(p2-p1)/np;

  for(int i=0;i<=np;i++)
    {
      double qt;

      //phasespace::setbounds(phasespace::mmin, phasespace::mmax, 0, qt, phasespace::ymin, phasespace::ymax);
      opts.nproc = 3;
      dyres::init();
      mcfm::init();
      iniflavreduce_();                      //need to call this after nproc_.nproc_ is set
      prop::init();                          //Initialise mass and width used in the propagator
      mesq::init();                          //EW couplings for born amplitudes
      vjint::init();
      vjloint::init();
      m = opts.rmass;
      qt = (i*hp+p1)*m;
      phasespace::set_mqtyphi(m, qt, y);//set global variables to costh, m, qt, y
      phasespace::set_cth(costh);//set global variables to costh, m, qt, y
      phasespace::calcexpy();
      phasespace::calcmt();
      omegaintegr::genV4p();//generate boson 4-momentum, with m, qt, y and phi=0
      phasespace::setbounds(phasespace::mmin, phasespace::mmax, 0, 100, phasespace::ymin, phasespace::ymax);
      cacheyrapint_(phasespace::ymin, phasespace::ymax);
      if (mode >= 2)
	if (opts.resumcpp)
	  {
	    rapint::cache(phasespace::ymin, phasespace::ymax);
	    rapint::allocate();
	    rapint::integrate(phasespace::ymin,phasespace::ymax,m);
	  }
      resint::rint(costh,m,qt,y,mode,f);
      double z = f[0];
      //double z = vjint::vint(m,qt,y);
      
      opts.nproc = 1;
      dyres::init();
      mcfm::init();
      iniflavreduce_();                      //need to call this after nproc_.nproc_ is set
      prop::init();                          //Initialise mass and width used in the propagator
      mesq::init();                          //EW couplings for born amplitudes
      vjint::init();
      vjloint::init();
      m = opts.rmass;
      qt = (i*hp+p1)*m;
      phasespace::set_mqtyphi(m, qt, y);//set global variables to costh, m, qt, y
      phasespace::set_cth(costh);//set global variables to costh, m, qt, y
      phasespace::calcexpy();
      phasespace::calcmt();
      omegaintegr::genV4p();//generate boson 4-momentum, with m, qt, y and phi=0
      phasespace::setbounds(phasespace::mmin, phasespace::mmax, 0, 100, phasespace::ymin, phasespace::ymax);
      cacheyrapint_(phasespace::ymin, phasespace::ymax);
      if (mode >= 2)
	if (opts.resumcpp)
	  {
	    rapint::cache(phasespace::ymin, phasespace::ymax);
	    rapint::allocate();
	    rapint::integrate(phasespace::ymin,phasespace::ymax,m);
	  }
      resint::rint(costh,m,qt,y,mode,f);
      double wp = f[0];
      //double wp = vjint::vint(m,qt,y);

      cout << i*hp+p1 << " " << wp << "  " << z << endl;
      pf << "gp->SetPoint(gp->GetN(), " << i*hp+p1 << ", " << wp/z << ");" << endl;

    }
  pf << "gp->Draw();" << endl;
  pf << "//gp1->Draw();" << endl;
  pf << "//gp2->Draw(\"same\");" << endl;
  pf << "}" << endl;
}

//lines for finite order part
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
  setcthmqty(costh,m,qt,y);
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
