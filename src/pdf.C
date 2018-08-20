#include "pdf.h"
#include "dyres_interface.h"
#include "interface.h"
#include "settings.h"

#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/LHAGlue.h>
#include <math.h>

void pdf::init()
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

  LHAPDF::initPDFSet(opts.LHAPDFset);
  //LHAPDF::initPDFSetByName(opts.LHAPDFset);
  LHAPDF::initPDF(opts.LHAPDFmember);

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

  // initialization of alphas
  setalphas();

  // read g from the PDF
  setg();

  //take the cmass and b mass from the PDF
  //      cmass=dsqrt(mcsq)
  //      bmass=dsqrt(mbsq)
}

//set value of alphas
void pdf::setalphas()
{
  couple_.amz_=LHAPDF::alphasPDF(dymasses_.zmass_);
  double scale = fabs(scale_.scale_);

  if (opts_.approxpdf_ == 1)
    {
      int nloop = 3;
      qcdcouple_.as_=dyalphas_mcfm_(scale,couple_.amz_,nloop);
    }
  else
    qcdcouple_.as_=dyalphas_lhapdf_(scale);
  
  qcdcouple_.ason2pi_=qcdcouple_.as_/(2*M_PI);
  qcdcouple_.ason4pi_=qcdcouple_.as_/(4*M_PI);
  qcdcouple_.gsq_=4*M_PI*qcdcouple_.as_;
}

//set the value of the g-parameter of the non perturbative form factor
void pdf::setg()
{
  LHAPDF::PDFInfo info(opts.LHAPDFset, opts.LHAPDFmember);
  double gformfactor = info.get_entry_as<double>("g", -1);
  if (gformfactor >= 0)
    {
      cout << "g form factor: input from PDF member: " << gformfactor << endl;
      opts.g_param = gformfactor;
      g_param_.g_param_ = gformfactor;
      np_.g_ = opts.g_param;
    }
}


void dysetpdf_(int& member)
{
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
  double fPDF[13];

  //set to zero if x out of range
  if (x > 1. || x <= 0.)
    {
      for (int i = -MAXNF; i <= MAXNF; i++)
	fx[MAXNF+i]=0.;
      return;
    }
  
  LHAPDF::xfx(x,xmu,fPDF);
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

}
