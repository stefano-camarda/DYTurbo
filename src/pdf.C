#include "pdf.h"
#include "dyres_interface.h"
#include "interface.h"
#include "scales.h"
#include "settings.h"
#include "coupling.h"
#include "alphas.h"
#include "isnan.h"
#include "vjint.h"
#include "mcfm_interface.h"

#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/LHAGlue.h>
#include <math.h>

LHAPDF::PDF* pdf::lhapdf = 0;

int pdf::order;
double pdf::xmin;
double pdf::qmin;
double pdf::mc;
double pdf::mb;
double pdf::mt;
double pdf::g;

void (*pdf::xfxq)(const double &x, const double &Q, double *fPDF) = 0;
double (*pdf::alphas)(const double &Q) = 0;

void pdf::lhaxfxq(const double &x, const double &Q, double *r)
{
  vector <double> fPDF; fPDF.resize(13);
  lhapdf->xfxQ(x,Q,fPDF);
  copy(fPDF.begin(),fPDF.end(),r);
}

double pdf::lhaalphas(const double &Q)
{
  return lhapdf->alphasQ(Q);
}

double pdf::rgktalphas(double q)
{
  //Interface to the Runge-Kutta alphas solution, with order matched to opts.order
  double q0 = coupling::zmass;
  double as0 = alphas(q0);
  alphas::rgkt(q,q0,as0);
  
  switch (opts.order)
    {
    case 0: return real(alphas::asLO);     break;
    case 1: return real(alphas::asNLO);    break;
    case 2: return real(alphas::asNNLO);   break;
    case 3: return real(alphas::asNNNLO);  break;
    case 4: return real(alphas::asNNNNLO); break;
    default: return 0.;
    }
}

void pdf::init()
{
  if (!opts.externalpdf)
    {
      //printf(" ==Initialize PDF set from LHAPDF==\n\n");
      //printf("\n");
      LHAPDF::Info& cfg = LHAPDF::getConfig();

      cfg.set_entry("Verbosity"             , 0              );
      cfg.set_entry("Interpolator"          , "logcubic"     );
      cfg.set_entry("Extrapolator"          , "continuation" );
      cfg.set_entry("ForcePositive"         , 0              );
      cfg.set_entry("AlphaS_Type"           , "analytic"     );
      cfg.set_entry("MZ"                    , 91.1876        );
      cfg.set_entry("MUp"                   , 0.002          );
      cfg.set_entry("MDown"                 , 0.005          );
      cfg.set_entry("MStrange"              , 0.10           );
      cfg.set_entry("MCharm"                , 1.29           );
      cfg.set_entry("MBottom"               , 4.19           );
      cfg.set_entry("MTop"                  , 172.9          );
      cfg.set_entry("Pythia6LambdaV5Compat" , true           );

      //Old interface
      LHAPDF::initPDFSet(opts.LHAPDFset); //LHAPDF::initPDFSetByName(opts.LHAPDFset);
      LHAPDF::initPDF(opts.LHAPDFmember);

      //New interface
      lhapdf = LHAPDF::mkPDF(opts.LHAPDFset, opts.LHAPDFmember);

      xfxq = lhaxfxq;
      //if (opts.lhaalphas)
      alphas = lhaalphas;
      //else
      //alphas = rgktalphas
      if (opts.PDFerrors && LHAPDF::numberPDF() > 1)
	{
	  opts.totpdf = LHAPDF::numberPDF()+1;
	  pdferropts_.pdferr_ = true;
	  pdferropts_.totpdf_ = LHAPDF::numberPDF()+1;
	}
      else
	{
	  opts.totpdf = 1;
	  pdferropts_.pdferr_ = false;
	  pdferropts_.totpdf_ = 1;
	}

      order = LHAPDF::getOrderPDF(); //order of evolution

      LHAPDF::PDFInfo info(opts.LHAPDFset, opts.LHAPDFmember);
      qmin = info.get_entry_as<double>("QMin", -1);
      xmin = info.get_entry_as<double>("XMin", -1);

      mc = LHAPDF::getThreshold(4);
      mb = LHAPDF::getThreshold(5);
      mt = LHAPDF::getThreshold(6);
    }
  else
    {
      opts.totpdf = 1;
      pdferropts_.pdferr_ = false;
      pdferropts_.totpdf_ = 1;
    }

  // initialization of alphas
  setalphas();

  //read g from the PDF
  setg();
  
  //take the cmass and b mass from the PDF
  //      cmass=dsqrt(mcsq)
  //      bmass=dsqrt(mbsq)
}

//set value of alphas at the renormalization scale in the fortran common blocks
void pdf::setalphas()
{
  //Old LHAPDF interface
  //double asmz = LHAPDF::alphasPDF(coupling::zmass);

  //New LHAPDF interface
  //double asmz = lhapdf->alphasQ(coupling::zmass);

  //Allows external alphas
  //couple_.amz_ = alphas(dymasses_.zmass_);
  double asmz = alphas(coupling::zmass);
  
  //run alphas
  double as;
  if (opts_.approxpdf_ == 1)
    {
      int nloop = 3;
      as = dyalphas_mcfm_(scales::ren,couple_.amz_,nloop);
    }
  else
    //as = dyalphas_lhapdf_(scales::ren);
    //as = LHAPDF::alphasPDF(scales::ren); //Old LHAPDF interface
    //as = lhapdf->alphasQ(scales::ren);   //New LHAPDF interface
    if (opts.alphaslha)
      as = alphas(scales::ren);		  //Allows external alphas
    else
      as = rgktalphas(scales::ren);

  //Set alphas and strong coupling in the fortran common blocks

  //MCFM
  qcdcouple_.as_ = as;
  qcdcouple_.ason2pi_ = as/(2*M_PI);
  qcdcouple_.ason4pi_ = as/(4*M_PI);
  qcdcouple_.gsq_= 4*M_PI*as;

  //vjet analytical calculation
  asnew_.as_ = as/M_PI;
  asp_.asp_ = as;
}

//set the value of the g-parameter of the non perturbative form factor
void pdf::setg()
{
  if (!opts.externalpdf)
    {
      LHAPDF::PDFInfo info(opts.LHAPDFset, opts.LHAPDFmember);
      g = info.get_entry_as<double>("g", -1);
    }
  if (g >= 0)
    {
      cout << "g form factor: input from PDF member: " << g << endl;
      opts.g1 = g;
      g_param_.g_param_ = g;
      np_.g_ = g;
    }
}


void dysetpdf_(int& member)
{
  if (opts.externalpdf) return;
  
  if (member == 0)
    {
      if (opts.PDFerrors && opts.totpdf > 1)
	LHAPDF::initPDF(0);
      else
	LHAPDF::initPDF(opts.LHAPDFmember);
    }
  else
    LHAPDF::initPDF(member);
  
  pdf::setalphas();
  //  setg(); //set g separately when setting a different PDF in the resummed part
}


void setmellinpdf_(int &member){
    // if member is still same than dont do anything
    //if (lastMember==member) return;
    //if (v_mellinpdf.size()<member) v_mellinpdf.
    // test current flag
    //if (v_mellinpdf[lastMember].isInitialized)
    // weights
    // moms
    // set init flag
    //lastMember=member;
}

void fdist_(int& ih, double& x, double& xmu, double fx[2*MAXNF+1])
{
  //set to zero if x out of range
  //check for nans
  if (isnan_ofast(x))
    {
      cout << "Bjorken-x is nan in fdist_" << endl;
      for (int i = -MAXNF; i <= MAXNF; i++)
	fx[MAXNF+i]=0.;
      return;
    }
  if (x > 1. || x <= 0.)
    {
      for (int i = -MAXNF; i <= MAXNF; i++)
	fx[MAXNF+i]=0.;
      return;
    }

  //Old LHAPDF interface
  //double fPDF[13];
  //LHAPDF::xfx(x,xmu,fPDF);

  //New LHAPDF interface
  //vector <double> fPDF; fPDF.resize(13);
  //pdf::lhapdf->xfxQ(x,xmu,fPDF);

  //Allows external PDF
  double fPDF[13];
  pdf::xfxq(x,xmu,fPDF);

  //cout << "fdist " << setprecision(16) << fPDF[6] << endl;
  //vector<int> pids = pdf::lhapdf->flavors();
  //for (int i = 0; i < pids.size(); i++)
  //cout << pids[i] << endl;

  //cout << " fdist " << endl;
  //for (int i = -MAXNF; i <= MAXNF; i++)
  //cout << fPDF[6+i]/x << "  ";
  //cout << endl;
  
  //xmu *= xmu;
  //fPDF[6-5] = pdf::lhapdf->xfxQ(-5,x,xmu);
  //fPDF[6-4] = pdf::lhapdf->xfxQ(-4,x,xmu);
  //fPDF[6-3] = pdf::lhapdf->xfxQ(-3,x,xmu);
  //fPDF[6-2] = pdf::lhapdf->xfxQ(-2,x,xmu);
  //fPDF[6-1] = pdf::lhapdf->xfxQ(-1,x,xmu);
  //fPDF[6+0] = pdf::lhapdf->xfxQ(21,x,xmu);
  //fPDF[6+1] = pdf::lhapdf->xfxQ(1 ,x,xmu);
  //fPDF[6+2] = pdf::lhapdf->xfxQ(2 ,x,xmu);
  //fPDF[6+3] = pdf::lhapdf->xfxQ(3 ,x,xmu);
  //fPDF[6+4] = pdf::lhapdf->xfxQ(4 ,x,xmu);
  //fPDF[6+5] = pdf::lhapdf->xfxQ(5 ,x,xmu);

  //  for (int i = -MAXNF; i <= MAXNF; i++)
  //    cout << fPDF[6+i]/x << "  ";
  //  cout << endl;
  

  if (ih == 1) //proton
    for (int i = -MAXNF; i <= MAXNF; i++)
      fx[MAXNF+i]=fPDF[6+i]/x;
  else if (ih == -1) //antiproton
    for (int i = -MAXNF; i <= MAXNF; i++)
      fx[MAXNF+i]=fPDF[6-i]/x;

  //switch off flavours
  //fx[MAXNF-5]=0.; //bbar
  //fx[MAXNF-4]=0.; //cbar
  //fx[MAXNF-3]=0.; //sbar
  //fx[MAXNF-2]=0.;
  //fx[MAXNF-1]=0.;
  //fx[MAXNF+0]=0.; //gluon
  //fx[MAXNF+1]=0.;
  //fx[MAXNF+2]=0.;
  //fx[MAXNF+3]=0.; //s
  //fx[MAXNF+4]=0.; //c
  //fx[MAXNF+5]=0.; //b

  //impose positivity
  //fx[MAXNF-5]=max(0.,fx[MAXNF-5]);
  //fx[MAXNF-4]=max(0.,fx[MAXNF-4]);
  //fx[MAXNF-3]=max(0.,fx[MAXNF-3]);
  //fx[MAXNF-2]=max(0.,fx[MAXNF-2]);
  //fx[MAXNF-1]=max(0.,fx[MAXNF-1]);
  //fx[MAXNF+0]=max(0.,fx[MAXNF+0]);
  //fx[MAXNF+1]=max(0.,fx[MAXNF+1]);
  //fx[MAXNF+2]=max(0.,fx[MAXNF+2]);
  //fx[MAXNF+3]=max(0.,fx[MAXNF+3]);
  //fx[MAXNF+4]=max(0.,fx[MAXNF+4]);
  //fx[MAXNF+5]=max(0.,fx[MAXNF+5]);

  //make u = d and ubar = dbar
  //fx[MAXNF-1]=fx[MAXNF-2];
  //fx[MAXNF+1]=fx[MAXNF+2];

  /*
  if (x < 90./13000.*exp(-1.))
    {
      fx[MAXNF-5]=0.;
      fx[MAXNF-4]=0.;
      fx[MAXNF-3]=0.;
      fx[MAXNF-2]=0.;
      fx[MAXNF-1]=0.;
      fx[MAXNF+0]=0.;
      fx[MAXNF+1]=0.;
      fx[MAXNF+2]=0.;
      fx[MAXNF+3]=0.;
      fx[MAXNF+4]=0.;
      fx[MAXNF+5]=0.;
    }
  */

  //flat PDFs
  //fx[MAXNF-5]=1.; //bbar
  //fx[MAXNF-4]=1.; //cbar
  //fx[MAXNF-3]=1.; //sbar
  //fx[MAXNF-2]=1.;
  //fx[MAXNF-1]=1.;
  //fx[MAXNF+0]=1.; //gluon
  //fx[MAXNF+1]=1.;
  //fx[MAXNF+2]=1.;
  //fx[MAXNF+3]=1.; //s
  //fx[MAXNF+4]=1.; //c
  //fx[MAXNF+5]=1.; //b
}
