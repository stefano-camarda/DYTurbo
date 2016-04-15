#include "pdf.h"
#include "interface.h"
#include "settings.h"

#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/LHAGlue.h>
#include <math.h>

void pdfini_()
{
  //printf(" ==Initialize PDF set from LHAPDF==\n\n");
  //printf("\n");
  LHAPDF::Info& cfg = LHAPDF::getConfig();
  cfg.set_entry("Verbosity", 0);
  LHAPDF::initPDFSetByName(opts.LHAPDFset);
  LHAPDF::initPDF(opts.LHAPDFmember);
  // initialization of alphas
  couple_.amz_=LHAPDF::alphasPDF(dymasses_.zmass_) ;
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

  //setalphas();
  //setg();
}

void setpdf_(int& member)
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
  
  setalphas();
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

//set value of alphas
void setalphas()
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

//set value of g
void setg()
{
  LHAPDF::PDFInfo info(opts.LHAPDFset, opts.LHAPDFmember);
  double gformfactor = info.get_entry_as<double>("g", -1);
  if (gformfactor >= 0)
  {
    cout << "g form factor: input from PDF member: " << gformfactor << endl;
    opts.g_param = gformfactor;
    g_param_.g_param_ = gformfactor;
  }
}

void fdist_(int& ih, double& x, double& xmu, double fx[11])
{
  double fPDF[13];

  //set to zero if x out of range
  if (x > 1.)
    for (int i = -5; i <=5; i++)
      fx[5+i]=0.;
  //!!!! should return here??? !!!!

  
  LHAPDF::xfx(x,xmu,fPDF);
  if (ih == 1) //proton
    for (int i = -5; i <=5; i++)
      fx[5+i]=fPDF[6+i]/x;
  else if (ih == -1) //antiproton
    for (int i = -5; i <=5; i++)
      fx[5+i]=fPDF[6-i]/x;

  //fx[5-5]=0.;
  //fx[5-4]=0.;
  //fx[5-3]=0.;
  //fx[5-2]=0.;
  //fx[5-1]=0.;
  //fx[5+0]=0.;
  //fx[5+1]=0.;
  //fx[5+2]=0.;
  //fx[5+3]=0.;
  //fx[5+4]=0.;
  //fx[5+5]=0.;
}
