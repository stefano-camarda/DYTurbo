#include "luminosity.h"

#include "parton.h"
#include "vjint.h"
#include "coupling.h"
#include "settings.h"
#include "pdf.h"
#include <iostream>

double luminosity::fh1[2*MAXNF+1];
double luminosity::fh2[2*MAXNF+1];
double luminosity::sumckm[MAXNF];

double luminosity::mesq[MAXNF][MAXNF];
double luminosity::mesq_nointerf[MAXNF][MAXNF];

double luminosity::sum_mesq;
double luminosity::sum_mesq_nointerf;
double luminosity::flqg[MAXNF]; // -> rename sum_mesqi
double luminosity::flgq[MAXNF];

double luminosity::sw;
double luminosity::cw;

void luminosity_calc_()
{
  luminosity::calc();
}

void luminosity_pdf1_(double & x)
{
  luminosity::pdf1(x);
}
void luminosity_pdf2_(double & x)
{
  luminosity::pdf2(x);
}

void luminosity::init()
{
  //--> bug fix in dyqt, sumckm are attributed incorrectly (make a difference only when switching off the diagonal ckm elements)
  //sumckm[0]=pow(quarks_.ckm_[1][0],2)+pow(quarks_.ckm_[1][3],2);                           //sum dx -> apply to u
  //sumckm[1]=pow(quarks_.ckm_[1][0],2)+pow(quarks_.ckm_[2][0],2)+pow(quarks_.ckm_[4][0],2); //sum ux -> apply to d
  //sumckm[2]=pow(quarks_.ckm_[1][3],2)+pow(quarks_.ckm_[2][3],2)+pow(quarks_.ckm_[4][3],2); //sum cx -> apply to s
  //sumckm[3]=pow(quarks_.ckm_[2][0],2)+pow(quarks_.ckm_[2][3],2);                           //sum sx -> apply to c
  //sumckm[4]=pow(quarks_.ckm_[1][5],2)+pow(quarks_.ckm_[2][5],2)+pow(quarks_.ckm_[4][5],2); //sum tx -> apply to b

  sumckm[1]=pow(quarks_.ckm_[1][0],2)+pow(quarks_.ckm_[1][3],2);                           //sum dx -> apply to d
  sumckm[0]=pow(quarks_.ckm_[1][0],2)+pow(quarks_.ckm_[2][0],2)+pow(quarks_.ckm_[4][0],2); //sum ux -> apply to u
  sumckm[3]=pow(quarks_.ckm_[1][3],2)+pow(quarks_.ckm_[2][3],2)+pow(quarks_.ckm_[4][3],2); //sum cx -> apply to c
  sumckm[2]=pow(quarks_.ckm_[2][0],2)+pow(quarks_.ckm_[2][3],2);                           //sum sx -> apply to s
  sumckm[4]=pow(quarks_.ckm_[4][0],2)+pow(quarks_.ckm_[4][3],2);                           //sum bx -> apply to b
  //sumckm[5]=pow(quarks_.ckm_[1][5],2)+pow(quarks_.ckm_[2][5],2)+pow(quarks_.ckm_[4][5],2); //sum tx -> apply to t
  
  sw = sqrt(coupling::xw); // sin_w
  cw = sqrt(1.-coupling::xw); // cos_w
}

void luminosity::cache()
{
  //mesq[i][j] are the costh integrated born amplitudes
  if (opts.nproc == 1 || opts.nproc == 2)
    {
      for (int i = 0; i < MAXNF; i++)
	for (int j = 0; j < MAXNF; j++)
	  mesq[i][j] = 1./(2.*coupling::xw)*pow(quarks_.ckm_[j][i],2)/2.*sigs_.sigw_;
    }
  if (opts.nproc == 3)
    {
      if (opts.useGamma)
	{
	  for (int i = 0; i < MAXNF; i++)
	    for (int j = 0; j < MAXNF; j++)
	      {
		mesq[i][j] = 2.*pow(quarks_.eq_[j],2)*quarks_.delta_[j][i]*sigs_.siggamma_+
		  (pow(1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2)+
		   pow(-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2))*sigs_.sigz_ +
		  (((1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw)+
		    (-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw))*2.*quarks_.eq_[j])*sigs_.sigint_;
		
		mesq_nointerf[i][j] = 2.*pow(quarks_.eq_[j],2)*quarks_.delta_[j][i]*sigs_.siggamma_+
		  (pow(1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2)+
		   pow(-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2))*sigs_.sigz_;
	      }

	}
      else
	{
	  for (int i = 0; i < MAXNF; i++)
	    for (int j = 0; j < MAXNF; j++)
	      {
		mesq[i][j] = (pow(1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2)+
			      pow(-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2))*sigs_.sigz_;
		mesq_nointerf[i][j] = mesq[i][j];
	      }
	}
    }




  
  sum_mesq = 0.;
  sum_mesq_nointerf = 0.;
  for (int i = 0; i < MAXNF; i++)
    flqg[i] = 0;
  
  if (opts.nproc == 1)
    {
      for (int i = 0; i < MAXNF; i++)
	for (int j = 0; j < MAXNF; j++)
	  sum_mesq += mesq[i][j];//1./(2.*coupling::xw)*pow(quarks_.ckm_[j][i],2)/2.*sigs_.sigw_;

      //(assuming CKM unitarity: 1./(2.*coupling::xw)*sigs_.sigw_ = sum_mesq/2.;)
      //sum_mesq = 1./(2.*coupling::xw)*sigs_.sigw_;
    }


  if (opts.nproc == 2)
    {
      for (int i = 0; i < MAXNF; i++)
	for (int j = 0; j < MAXNF; j++)
	  sum_mesq += mesq[i][j];//1./(2.*coupling::xw)*pow(quarks_.ckm_[j][i],2)/2.*sigs_.sigw_; // -> ckm unitarity implies sum_ij ckm[i][j]^2 / 2 = 1

      //(assuming CKM unitarity: 1./(2.*coupling::xw)*sigs_.sigw_ = sum_mesq/2.;)
      //sum_mesq = 1./(2.*coupling::xw)*sigs_.sigw_;
    }

  if (opts.nproc == 3)
    {
      if (opts.useGamma)
	{
	  for (int i = 0; i < MAXNF; i++)
	    for (int j = 0; j < MAXNF; j++)
	      {
		sum_mesq += mesq[i][j];//2.*pow(quarks_.eq_[j],2)*quarks_.delta_[j][i]*sigs_.siggamma_+
		//(pow(1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2)+
		//pow(-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2))*sigs_.sigz_ +
		//(((1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw)+
		//(-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw))*2.*quarks_.eq_[j])*sigs_.sigint_; // -> equal to sum_ij mesq[i][j] integrated wrt costh?
	      
		sum_mesq_nointerf +=  mesq_nointerf[i][j]; //2.*pow(quarks_.eq_[j],2)*quarks_.delta_[j][i]*sigs_.siggamma_+
		//(pow(1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2)+
		//pow(-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2))*sigs_.sigz_;
	      }
	  for (int i = 0; i < MAXNF; i++)
	    for (int j = 0; j < MAXNF; j++)
	      flqg[i] +=  mesq[i][j];//2.*pow(quarks_.eq_[j],2)*quarks_.delta_[j][i]*sigs_.siggamma_+
	  //(pow(1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2)+
	  //pow(-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2))*sigs_.sigz_ +
	  //(((1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw)+
	  //(-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw))*2.*quarks_.eq_[j])*sigs_.sigint_; // -> equal to sum_i mesq[i][j] =  mesq[i][-i] integrated wrt costh? 

	}
      else
	{
	  for (int i = 0; i < MAXNF; i++)
	    for (int j = 0; j < MAXNF; j++)
	      {
		sum_mesq +=  mesq[i][j];//(pow(1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2)+
		//pow(-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2))*sigs_.sigz_;
		sum_mesq_nointerf += mesq[i][j];
	      }

	  for (int i = 0; i < MAXNF; i++)
	    for (int j = 0; j < MAXNF; j++)
	      flqg[i] +=  mesq[i][j];//(pow(1./(2.*sw*cw)*quarks_.tau3_[j][i]-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2)+
	  //pow(-quarks_.delta_[j][i]*quarks_.eq_[j]*sw/cw,2))*sigs_.sigz_;
	}
    }
}

void luminosity::pdf1(double x)
{
  //set to zero if x out of range
  if (x > 1.)
    {
      for (int i = 0; i < 2*MAXNF+1; i++)
	fh1[i]=0.;
      return;
    }

  //call fdist instead of evolvePDF, so that all PDF calls in x space
  //pass through the same interface, and allows easily to switch off PDFs
  fdist_(opts.ih1,x,scales2_.xmuf_,fh1);

  // Change u and d quarks (convert from LHAPDF to DYRES convention)
  double temp1=fh1[parton::D];
  fh1[parton::u]=fh1[parton::U];
  fh1[parton::d]=temp1;
  double temp2=fh1[parton::Db];
  fh1[parton::ub]=fh1[parton::Ub];
  fh1[parton::db]=temp2;
}      

void luminosity::pdf2(double x)
{
  //set to zero if x out of range
  if (x > 1.)
    {
      for (int i = 0; i < 2*MAXNF+1; i++)
	fh2[i] = 0.;
      return;
    }

  //call fdist instead of evolvePDF, so that all PDF calls in x space
  //pass through the same interface, and allows easily to switch off PDFs
  fdist_(opts.ih2,x,scales2_.xmuf_,fh2);

  // Change u and d quarks (convert from LHAPDF to DYRES convention)
  double temp1=fh2[parton::D];
  fh2[parton::u]=fh2[parton::U];
  fh2[parton::d]=temp1;
  double temp2=fh2[parton::Db];
  fh2[parton::ub]=fh2[parton::Ub];
  fh2[parton::db]=temp2;
}      

void luminosity::calc()
{
  int nf,ih1,ih2,ic,isetproton,nloop,ord,prodflag;
  //      include 'scales2_inc.f'
  //      include 'internal_inc.f'
  //      include 'luminosities_inc.f'
  //      include 'fodyqt_inc.f'
  //      common/quarks/eq(5),alq(5),arq(5),ckm(6,6),delta(5,5),tau3(5,5)

  //initialize
  luminosities_.xlumgg_ = 0.;
  luminosities_.xlumqg_ = 0.;
  luminosities_.xlumgq_ = 0.;
  luminosities_.xlumqgtr_ = 0.;
  luminosities_.xlumgqtr_ = 0.;
  luminosities_.xlumqqb_ = 0.;
  luminosities_.xlumqqbtr_ = 0.;
  luminosities_.xlumqqbdbb_ = 0.;
  luminosities_.xlumqqbdbc_ = 0.;
  luminosities_.xlumqqbdcc_ = 0.;
  luminosities_.xlumqqbddd_ = 0.;
  luminosities_.xlumqqbLL_ = 0.;
  luminosities_.xlumqqbLR_ = 0.;
  luminosities_.xlumqq_ = 0.;
  luminosities_.xlumqqead_ = 0.;
  luminosities_.xlumqqeaa_ = 0.; //not used
  luminosities_.xlumqqebb_ = 0.; //not used
  luminosities_.xlumqqLL_ = 0.;  //not used
  luminosities_.xlumqqLR_ = 0.;  //not used

  
  if (opts.nproc == 3)
    {
      for (int i = 0; i < MAXNF; i++)
	{
	  //qqb luminosity
	  luminosities_.xlumqqb_ +=
	    (fh1[i+MAXNF+1]*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]
	     +fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]*fh2[i+MAXNF+1])
	    *mesq[i][i];

	  //qg,gq luminosities
	  luminosities_.xlumqg_ += (fh1[i+MAXNF+1]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))])*fh2[parton::g] * flqg[i];
	  luminosities_.xlumgq_ += (fh2[i+MAXNF+1]+fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))])*fh1[parton::g] * flqg[i];
	}

      if (opts.order <= 1)
	return;

      for (int i = 0; i < MAXNF; i++)
	{
	  //triangular loop qg luminosity (only for Z production)
	  luminosities_.xlumqgtr_ += (fh1[i+MAXNF+1]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))])*fh2[parton::g] * quarks_.tau3_[i][i] * (-(1./(2.*sw*cw))*(1./(2.*sw*cw))*sigs_.sigz_);
	  luminosities_.xlumqgtr_ += (fh2[i+MAXNF+1]+fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))])*fh1[parton::g] * quarks_.tau3_[i][i] * (-(1./(2.*sw*cw))*(1./(2.*sw*cw))*sigs_.sigz_);

	  //triangular loop and Dab qqb luminosity (only for Z production)
	  luminosities_.xlumqqbtr_ += (fh1[i+MAXNF+1]*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]*fh2[i+MAXNF+1]) * (-(1./(2.*sw*cw))*(1./(2.*sw*cw))*sigs_.sigz_);
	}

      //gg luminosity
      luminosities_.xlumgg_ = fh1[parton::g]*fh2[parton::g]*sum_mesq;
  
      double sumpdf = 0.;
      for (int i = 0; i < MAXNF; i++)
	sumpdf += (fh1[i+MAXNF+1]*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]*fh2[i+MAXNF+1]);

      luminosities_.xlumqqbdbb_ = sumpdf*sum_mesq_nointerf;

      for (int i = 0; i < MAXNF; i++)
	luminosities_.xlumqqbdbc_ += (fh1[i+MAXNF+1]*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]*fh2[i+MAXNF+1]) * sum_mesq_nointerf;

      double sumpdf1 = 0.;
      double sumpdf2 = 0.;
      for (int i = 0; i < MAXNF; i++)
	{
	  sumpdf1 += (fh1[i+MAXNF+1]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]);
	  sumpdf2 += (fh2[i+MAXNF+1]+fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]);
	}

      for (int i = 0; i < MAXNF; i++)
	{
	  luminosities_.xlumqqbdcc_ += (fh1[i+MAXNF+1]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]) * sumpdf2 * mesq[i][i];
	  luminosities_.xlumqqbddd_ += (fh2[i+MAXNF+1]+fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]) * sumpdf1 * mesq[i][i];
	}

      //LL and LR qqb luminosity
      if (opts.useGamma)
	{
	  double sumpdf1g = 0.;
	  double sumpdf2g = 0.;
	  for (int i = 0; i < MAXNF; i++)
	    {
	      sumpdf1g += quarks_.eq_[i]*(fh1[i+MAXNF+1]-fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]);
	      sumpdf2g += quarks_.eq_[i]*(fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]-fh2[i+MAXNF+1]);
	    }
	  luminosities_.xlumqqbLL_ += 2.*sumpdf1g*sumpdf2g*sigs_.siggamma_;
	  luminosities_.xlumqqbLR_ += 2.*sumpdf1g*sumpdf2g*sigs_.siggamma_;
	}

      double sumpdf1az = 0.;
      double sumpdf2bz = 0.;
      double sumpdf1bz = 0.;
      double sumpdf2az = 0.;
      for (int i = 0; i < MAXNF; i++)
	{
	  sumpdf1az += (1./(2.*sw*cw)*quarks_.tau3_[i][i]-quarks_.eq_[i]*sw/cw)*fh1[i+MAXNF+1]+(quarks_.eq_[i]*sw/cw)*fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))];
	  sumpdf2bz += (1./(2.*sw*cw)*quarks_.tau3_[i][i]-quarks_.eq_[i]*sw/cw)*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]+(quarks_.eq_[i]*sw/cw)*fh2[i+MAXNF+1];
	  sumpdf1bz += (1./(2.*sw*cw)*quarks_.tau3_[i][i]-quarks_.eq_[i]*sw/cw)*fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]+(quarks_.eq_[i]*sw/cw)*fh1[i+MAXNF+1];
	  sumpdf2az += (1./(2.*sw*cw)*quarks_.tau3_[i][i]-quarks_.eq_[i]*sw/cw)*fh2[i+MAXNF+1]+(quarks_.eq_[i]*sw/cw)*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))];
	}
      luminosities_.xlumqqbLL_ += (sumpdf1az*sumpdf2bz+sumpdf1bz*sumpdf2az)*sigs_.sigz_;

      double sumpdf2bzm = 0.;
      double sumpdf2azm = 0.;
      for (int i = 0; i < MAXNF; i++)
	{
	  sumpdf2bzm += -(1./(2.*sw*cw)*quarks_.tau3_[i][i]-quarks_.eq_[i]*sw/cw)*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]+(-quarks_.eq_[i]*sw/cw)*fh2[i+MAXNF+1];
	  sumpdf2azm += -(1./(2.*sw*cw)*quarks_.tau3_[i][i]-quarks_.eq_[i]*sw/cw)*fh2[i+MAXNF+1]+(-quarks_.eq_[i]*sw/cw)*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))];
	}
      luminosities_.xlumqqbLR_ += (sumpdf1az*sumpdf2azm+sumpdf1bz*sumpdf2bzm)*sigs_.sigz_;

      //qq luminosity
      //for (int i = 0; i < MAXNF; i++)
      for (int i = 2; i < 3; i++)
	luminosities_.xlumqq_ += (fh1[i+MAXNF+1]*fh2[i+MAXNF+1]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]) *mesq_nointerf[i][i];
      
      //Ead and Ebc qq luminosity (here q = q')
      luminosities_.xlumqqead_ = luminosities_.xlumqq_;
    }
  if (opts.nproc == 1)
    {
      //qqb luminosity
      luminosities_.xlumqqb_ =
	(fh1[parton::u]  *fh2[parton::db]*pow(quarks_.ckm_[1][0],2) 
	 +fh1[parton::u] *fh2[parton::sb]*pow(quarks_.ckm_[2][0],2) 
	 +fh1[parton::u] *fh2[parton::bb]*pow(quarks_.ckm_[4][0],2) 
	 +fh1[parton::c] *fh2[parton::db]*pow(quarks_.ckm_[1][3],2) 
	 +fh1[parton::c] *fh2[parton::sb]*pow(quarks_.ckm_[2][3],2) 
	 +fh1[parton::c] *fh2[parton::bb]*pow(quarks_.ckm_[4][3],2) 
	 +fh1[parton::db]*fh2[parton::u] *pow(quarks_.ckm_[0][1],2)
	 +fh1[parton::db]*fh2[parton::c] *pow(quarks_.ckm_[3][1],2)
	 +fh1[parton::sb]*fh2[parton::u] *pow(quarks_.ckm_[0][2],2)
	 +fh1[parton::sb]*fh2[parton::c] *pow(quarks_.ckm_[3][2],2)
	 +fh1[parton::bb]*fh2[parton::u] *pow(quarks_.ckm_[0][4],2)
	 +fh1[parton::bb]*fh2[parton::c] *pow(quarks_.ckm_[3][4],2))
	*1./(2.*coupling::xw)*sigs_.sigw_;

      //the use of sumckm looks wrong! should exchange u <-> d, c <-> s, and t <-> b
      //qg,gq luminosities
      //(assuming CKM unitarity: 1./(2.*coupling::xw)*sigs_.sigw_ = sum_mesq/2.;)
      luminosities_.xlumqg_ = (fh1[parton::u]*sumckm[0]+fh1[parton::c]*sumckm[3]+fh1[parton::db]*sumckm[1]+fh1[parton::sb]*sumckm[2]+fh1[parton::bb]*sumckm[4])*fh2[parton::g]*1./(2.*coupling::xw)*sigs_.sigw_;
      luminosities_.xlumgq_ = (fh2[parton::u]*sumckm[0]+fh2[parton::c]*sumckm[3]+fh2[parton::db]*sumckm[1]+fh2[parton::sb]*sumckm[2]+fh2[parton::bb]*sumckm[4])*fh1[parton::g]*1./(2.*coupling::xw)*sigs_.sigw_;

      if (opts.order <= 1)
	return;
      
      //gg luminosity
      luminosities_.xlumgg_ = fh1[parton::g]*fh2[parton::g]*sum_mesq;
      
      double sumpdf = 0.;
      for (int i = 0; i < MAXNF; i++)
	sumpdf += (fh1[i+MAXNF+1]*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]*fh2[i+MAXNF+1]);

      luminosities_.xlumqqbdbb_ = sumpdf*sum_mesq;

      luminosities_.xlumqqbdbc_ =
      	(fh1[parton::db]*fh2[parton::d]*sumckm[1]
      	 +fh1[parton::sb]*fh2[parton::s]*sumckm[2]
      	 +fh1[parton::bb]*fh2[parton::b]*sumckm[4]
      	 +fh2[parton::ub]*fh1[parton::u]*sumckm[0]
      	 +fh2[parton::cb]*fh1[parton::c]*sumckm[3])
      	*1./(2.*coupling::xw)*sigs_.sigw_;

      double sumpdf2 = 0.;
      for (int i = 0; i < MAXNF; i++)
	sumpdf2 += (fh2[i+MAXNF+1]+fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]);
      
      luminosities_.xlumqqbdcc_ =
	(fh1[parton::u]*sumckm[0]
	 +fh1[parton::c]*sumckm[3]
	 +fh1[parton::db]*sumckm[1]
	 +fh1[parton::sb]*sumckm[2]
	 +fh1[parton::bb]*sumckm[4])
	*sumpdf2
	*1./(2.*coupling::xw)*sigs_.sigw_;

      double sumpdf1 = 0.;
      for (int i = 0; i < MAXNF; i++)
	sumpdf1 += (fh1[i+MAXNF+1]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]);
      
      luminosities_.xlumqqbddd_ =
	(fh2[parton::u]*sumckm[0]
	 +fh2[parton::c]*sumckm[3]
	 +fh2[parton::db]*sumckm[1]
	 +fh2[parton::sb]*sumckm[2]
	 +fh2[parton::bb]*sumckm[4])
	*sumpdf1
	*1./(2.*coupling::xw)*sigs_.sigw_;
      
      //qq luminosity
      luminosities_.xlumqq_ =
	(fh1[parton::u]  *fh2[parton::d]*pow(quarks_.ckm_[1][0],2) 
	 +fh1[parton::u] *fh2[parton::s]*pow(quarks_.ckm_[2][0],2) 
	 +fh1[parton::u] *fh2[parton::b]*pow(quarks_.ckm_[4][0],2) 
	 +fh1[parton::c] *fh2[parton::d]*pow(quarks_.ckm_[1][3],2) 
	 +fh1[parton::c] *fh2[parton::s]*pow(quarks_.ckm_[2][3],2) 
	 +fh1[parton::c] *fh2[parton::b]*pow(quarks_.ckm_[4][3],2) 
	 +fh1[parton::db]*fh2[parton::ub] *pow(quarks_.ckm_[0][1],2)
	 +fh1[parton::db]*fh2[parton::cb] *pow(quarks_.ckm_[3][1],2)
	 +fh1[parton::sb]*fh2[parton::ub] *pow(quarks_.ckm_[0][2],2)
	 +fh1[parton::sb]*fh2[parton::cb] *pow(quarks_.ckm_[3][2],2)
	 +fh1[parton::bb]*fh2[parton::ub] *pow(quarks_.ckm_[0][4],2)
	 +fh1[parton::bb]*fh2[parton::cb] *pow(quarks_.ckm_[3][4],2))
	*1./(2.*coupling::xw)*sigs_.sigw_;

      //Ead and Ebc qq luminosity (here q = q')
      luminosities_.xlumqqead_ =
      	(fh1[parton::u]*fh2[parton::u]*sumckm[0]
      	 +fh1[parton::c]*fh2[parton::c]*sumckm[3]
      	 +fh1[parton::db]*fh2[parton::db]*sumckm[1]
      	 +fh2[parton::sb]*fh1[parton::sb]*sumckm[2]
      	 +fh2[parton::bb]*fh1[parton::bb]*sumckm[4])
      	*1./(2.*coupling::xw)*sigs_.sigw_;
    }

  if (opts.nproc == 2)
    {
      //qqb luminosity
      luminosities_.xlumqqb_ =
        (fh1[parton::ub] *fh2[parton::d] *pow(quarks_.ckm_[1][0],2) 
	 +fh1[parton::ub]*fh2[parton::s] *pow(quarks_.ckm_[2][0],2) 
	 +fh1[parton::ub]*fh2[parton::b] *pow(quarks_.ckm_[4][0],2) 
	 +fh1[parton::cb]*fh2[parton::d] *pow(quarks_.ckm_[1][3],2) 
	 +fh1[parton::cb]*fh2[parton::s] *pow(quarks_.ckm_[2][3],2) 
	 +fh1[parton::cb]*fh2[parton::b] *pow(quarks_.ckm_[4][3],2) 
	 +fh1[parton::d] *fh2[parton::ub]*pow(quarks_.ckm_[0][1],2) 
	 +fh1[parton::d] *fh2[parton::cb]*pow(quarks_.ckm_[3][1],2) 
	 +fh1[parton::s] *fh2[parton::ub]*pow(quarks_.ckm_[0][2],2) 
	 +fh1[parton::s] *fh2[parton::cb]*pow(quarks_.ckm_[3][2],2) 
	 +fh1[parton::b] *fh2[parton::ub]*pow(quarks_.ckm_[0][4],2) 
	 +fh1[parton::b] *fh2[parton::cb]*pow(quarks_.ckm_[3][4],2))
	*1./(2.*coupling::xw)*sigs_.sigw_;

      //qg,gq luminosities
      luminosities_.xlumqg_ = (fh1[parton::ub]*sumckm[0]+fh1[parton::cb]*sumckm[3]+fh1[parton::d]*sumckm[1]+fh1[parton::s]*sumckm[2]+fh1[parton::b]*sumckm[4])*fh2[parton::g]*1./(2.*coupling::xw)*sigs_.sigw_;
      luminosities_.xlumgq_ = (fh2[parton::ub]*sumckm[0]+fh2[parton::cb]*sumckm[3]+fh2[parton::d]*sumckm[1]+fh2[parton::s]*sumckm[2]+fh2[parton::b]*sumckm[4])*fh1[parton::g]*1./(2.*coupling::xw)*sigs_.sigw_;

      if (opts.order <= 1)
	return;

      //gg luminosity
      luminosities_.xlumgg_ = fh1[parton::g]*fh2[parton::g]*sum_mesq;
      
      double sumpdf = 0.;
      for (int i = 0; i < MAXNF; i++)
	sumpdf += (fh1[i+MAXNF+1]*fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]*fh2[i+MAXNF+1]);

      luminosities_.xlumqqbdbb_ = sumpdf*sum_mesq;

      luminosities_.xlumqqbdbc_ =
      	(fh1[parton::d]*fh2[parton::db]*sumckm[1]
      	 +fh1[parton::s]*fh2[parton::sb]*sumckm[2]
      	 +fh1[parton::b]*fh2[parton::bb]*sumckm[4]
      	 +fh2[parton::u]*fh1[parton::ub]*sumckm[0]
      	 +fh2[parton::c]*fh1[parton::cb]*sumckm[3])
      	*1./(2.*coupling::xw)*sigs_.sigw_;

      double sumpdf2 = 0.;
      for (int i = 0; i < MAXNF; i++)
	sumpdf2 += (fh2[i+MAXNF+1]+fh2[parton::charge_conj(parton::pdgid(i+MAXNF+1))]);
      
      luminosities_.xlumqqbdcc_ =
	(fh1[parton::ub]*sumckm[0]
	 +fh1[parton::cb]*sumckm[3]
	 +fh1[parton::d]*sumckm[1]
	 +fh1[parton::s]*sumckm[2]
	 +fh1[parton::b]*sumckm[4])
	*sumpdf2
	*1./(2.*coupling::xw)*sigs_.sigw_;

      double sumpdf1 = 0.;
      for (int i = 0; i < MAXNF; i++)
	sumpdf1 += (fh1[i+MAXNF+1]+fh1[parton::charge_conj(parton::pdgid(i+MAXNF+1))]);

      luminosities_.xlumqqbddd_ =
	(fh2[parton::ub]*sumckm[0]
	 +fh2[parton::cb]*sumckm[3]
	 +fh2[parton::d]*sumckm[1]
	 +fh2[parton::s]*sumckm[2]
	 +fh2[parton::b]*sumckm[4])
	*sumpdf1
	*1./(2.*coupling::xw)*sigs_.sigw_;
      
      luminosities_.xlumqq_ =
        (fh1[parton::d] *fh2[parton::u] *pow(quarks_.ckm_[1][0],2) 
	 +fh1[parton::s]*fh2[parton::u] *pow(quarks_.ckm_[2][0],2) 
	 +fh1[parton::b]*fh2[parton::u] *pow(quarks_.ckm_[4][0],2) 
	 +fh1[parton::d]*fh2[parton::c] *pow(quarks_.ckm_[1][3],2) 
	 +fh1[parton::s]*fh2[parton::c] *pow(quarks_.ckm_[2][3],2) 
	 +fh1[parton::b]*fh2[parton::c] *pow(quarks_.ckm_[4][3],2) 
	 +fh1[parton::db] *fh2[parton::ub]*pow(quarks_.ckm_[0][1],2) 
	 +fh1[parton::db] *fh2[parton::cb]*pow(quarks_.ckm_[3][1],2) 
	 +fh1[parton::sb] *fh2[parton::ub]*pow(quarks_.ckm_[0][2],2) 
	 +fh1[parton::sb] *fh2[parton::cb]*pow(quarks_.ckm_[3][2],2) 
	 +fh1[parton::bb] *fh2[parton::ub]*pow(quarks_.ckm_[0][4],2) 
	 +fh1[parton::bb] *fh2[parton::cb]*pow(quarks_.ckm_[3][4],2))
	*1./(2.*coupling::xw)*sigs_.sigw_;

      //Ead and Ebc qq luminosity (here q = q')
      luminosities_.xlumqqead_ =
      	(fh1[parton::d]*fh2[parton::d]*sumckm[1]
      	 +fh1[parton::s]*fh2[parton::s]*sumckm[2]
      	 +fh1[parton::b]*fh2[parton::b]*sumckm[4]
      	 +fh2[parton::ub]*fh1[parton::ub]*sumckm[0]
      	 +fh2[parton::cb]*fh1[parton::cb]*sumckm[3])
      	*1./(2.*coupling::xw)*sigs_.sigw_;
    }
  
  //luminosities_.xlumqqb_ = 0.;
  //luminosities_.xlumqg_ = 0.;
  //luminosities_.xlumgq_ = 0.;
  //luminosities_.xlumgg_ = 0.;
  //luminosities_.xlumqgtr_ = 0.;
  //luminosities_.xlumgqtr_ = 0.;
  //luminosities_.xlumqqbtr_ = 0.;
  //luminosities_.xlumqqbdbb_ = 0.;
  //luminosities_.xlumqqbdbc_ = 0.;
  //luminosities_.xlumqqbdcc_ = 0.;
  //luminosities_.xlumqqbddd_ = 0.;
  //luminosities_.xlumqqbLL_ = 0.;
  //luminosities_.xlumqqbLR_ = 0.;
  //luminosities_.xlumqq_ = 0.;
  //luminosities_.xlumqqead_ = 0.;
  //luminosities_.xlumqqeaa_ = 0.;//not used
  //luminosities_.xlumqqebb_ = 0.;//not used
  //luminosities_.xlumqqLL_ = 0.;	//not used
  //luminosities_.xlumqqLR_ = 0.;	//not used
}
