#include "evolnum.h"

#include "mellinpdf.h"
#include "mellinint.h"
#include "pdfevol.h"
#include "resconst.h"
#include "parton.h"
#include "phasespace.h"
#include "settings.h"

using namespace parton;

//Mellin transform from x to N space at the scale b0/b
void evolnum::calculate()
{
  clock_t begin_time, end_time;
  begin_time = clock();  

  double facscale = fabs(pdfevol::mubstartilde);  //mubstartilde = mubstar * Q / sqrt((pow(mubstar,2) + pow(Q,2))), with Q = mures or mufac
  
  mellinpdf::allocate();
  //mellinpdf::evalpdfs(facscale);
  mellinpdf::evalpdfs(facscale,phasespace::m,phasespace::y);
  mellinpdf::update_mellin();
  mellinpdf::transform();
  
  end_time = clock();  
  //cout << "facscale is " << facscale << " x to N done in "  << float(end_time - begin_time)/CLOCKS_PER_SEC*1000. << "ms" << endl;

  complex <double> fx[11] = {0};
  if (opts.mellin1d)
    for (int k = 0; k < mellinint::mdim; k++)
      {
	fx[g]  = mellinpdf::GL[k];
	fx[u]  = mellinpdf::UP[k];
	fx[ub] = mellinpdf::UB[k];
	fx[d]  = mellinpdf::DO[k];
	fx[db] = mellinpdf::DB[k];
	fx[s]  = mellinpdf::ST[k];
	fx[sb] = mellinpdf::SB[k];
	if (resconst::NF >= 4)
	  {
	    fx[c]  = mellinpdf::CH[k];
	    fx[cb] = mellinpdf::CB[k];
	  }
	if (resconst::NF >= 5)
	  {
	    fx[b]  = mellinpdf::BO[k];
	    fx[bb] = mellinpdf::BB[k];
	  }
	//cout << k << "  " << fx[g] << endl;
	pdfevol::storemoments(k, fx);
      }
  else
    for (int k = 0; k < mellinint::mdim; k++)
      {
	fx[g]  = mellinpdf::GL_1[k];
	fx[u]  = mellinpdf::UP_1[k];
	fx[ub] = mellinpdf::UB_1[k];
	fx[d]  = mellinpdf::DO_1[k];
	fx[db] = mellinpdf::DB_1[k];
	fx[s]  = mellinpdf::ST_1[k];
	fx[sb] = mellinpdf::SB_1[k];
	if (resconst::NF >= 4)
	  {
	    fx[c]  = mellinpdf::CH_1[k];
	    fx[cb] = mellinpdf::CB_1[k];
	  }
	if (resconst::NF >= 5)
	  {
	    fx[b]  = mellinpdf::BO_1[k];
	    fx[bb] = mellinpdf::BB_1[k];
	  }
	//cout << k << "  " << fx[g] << endl;
	pdfevol::storemoments_1(k, fx);
	
	fx[g]  = mellinpdf::GL_2[k];
	fx[u]  = mellinpdf::UP_2[k];
	fx[ub] = mellinpdf::UB_2[k];
	fx[d]  = mellinpdf::DO_2[k];
	fx[db] = mellinpdf::DB_2[k];
	fx[s]  = mellinpdf::ST_2[k];
	fx[sb] = mellinpdf::SB_2[k];
	if (resconst::NF >= 4)
	  {
	    fx[c]  = mellinpdf::CH_2[k];
	    fx[cb] = mellinpdf::CB_2[k];
	  }
	if (resconst::NF >= 5)
	  {
	    fx[b]  = mellinpdf::BO_2[k];
	    fx[bb] = mellinpdf::BB_2[k];
	  }
	//cout << k << "  " << fx[g] << endl;
	pdfevol::storemoments_2(k, fx);
      }

  mellinpdf::free();
}


/*
//Old fortran code
void evolnum::calculate(int i)
{
  //N flavour dependence
  int nf = resconst::NF;

  int hadron = 1;
  //double facscale = fabs(bscale);
  //double facscale = fabs(bstarscale);
  //double facscale = fabs(opts.muf);
  double facscale = fabs(mubstartilde); //mubstartilde = mubstar * Q / sqrt((pow(mubstar,2) + pow(Q,2))), with Q = mures or mufac
  fcomplex XN = fcx(mellinint::Np[i]);
  double xmin = 1e-8;

  complex <double> fx[11];
  
  fcomplex uval,dval,usea,dsea,s,sbar,glu,charm,bot;
  pdfmoments_(hadron,facscale,XN,uval,dval,usea,dsea,s,sbar,glu,charm,bot,xmin);

  fx[0+MAXNF] = cx(glu);
  fx[1+MAXNF] = cx(uval) + cx(usea);
  fx[-1+MAXNF] = cx(usea);
  fx[2+MAXNF] = cx(dval) + cx(dsea);
  fx[-2+MAXNF] = cx(dsea);
  fx[3+MAXNF] = cx(s);
  fx[-3+MAXNF] = cx(sbar);
  if (nf >= 4)
    {
      fx[4+MAXNF] = cx(charm);
      fx[-4+MAXNF] = cx(charm);
    }
  else
    {
      fx[4+MAXNF] = 0.;
      fx[-4+MAXNF] = 0.;
    }
  
  if (nf >= 5)
    {
      fx[5+MAXNF] = cx(bot);
      fx[-5+MAXNF] = cx(bot);
    }
  else
    {
      fx[5+MAXNF] = 0.;
      fx[-5+MAXNF] = 0.;
    }

  storemoments(i, fx);
}
*/
