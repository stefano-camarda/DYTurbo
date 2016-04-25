#include <iostream>
#include <LHAPDF/LHAPDF.h>
//#include <gsl/gsl_integration.h>
//#include <gsl/gsl_math.h>
#include <ctime>

#include "settings.h"
#include "interface.h"
#include "finintegr.h"
#include "ctintegr.h"
#include "finitemapping.h"


using namespace std;

//these function wrappers map the variables into the unity hypercube of the vegas integration
double dyreal(double m, double y, double qt, double phicm, double phiZ, double cos_th, double zcth, double mjj, double phijj, double costhjj)
{
  //set up constants
  double rmass = opts.rmass;
  double rwidth = opts.rwidth;
  double minmass = sqrt(taumin_.taumin_)*energy_.sroot_;
  double sqrts = energy_.sroot_;;
  double s = sqrts*sqrts;
  double mmin = opts.mlow;
  double mmax = opts.mhigh;
  double xqtcut = qtcut_.xqtcut_;

  //  std::cout << std::setprecision(15);
  double begin_time, end_time;
  double value;
  double rre[22];
  double wgt = 1;

  /*
  //variables to be integrated (4 dimensions in total)
  //1 collinear PDF dimension
  double zcth;
  zcth = 0.5;   //polar angle of Z in the CM frame
  //3 dimensions of real radiation
  double mjj, phijj, costhjj;
  mjj = 10.;    //invariant mass of the dijet (set to 0 to have the correct virtual phase space mapping)
  phijj = 0.3;  //phi of the dijet in dijet rest frame
  costhjj = 0.1;//costh of the dijet in dijet rest frame
  */
  
  //evaluate four momentum of Z in lab frame, given y, m, qt
  double zql, zqt, zp, ze;
  zql = 0.5*sqrt(m*m+qt*qt)*(exp(y)-exp(-y));
  zqt = qt;
  zp = sqrt(zqt*zqt + zql*zql);
  ze = sqrt(m*m+zqt*zqt+zql*zql);

  //Compute Z four momentum in CM, given y, m, qt and zcth
  double zqtcm, zpcm, zqlcm, zecm, zycm;
  zqtcm = qt;
  zpcm = zqtcm / sin(acos(zcth)); //absolute momentum of Z in CM
  zqlcm = zpcm * zcth;         //longitudinal momentum of Z in CM
  zecm = sqrt(m*m + zpcm*zpcm);//Energy of Z in CM
  zycm = 0.5 * log((zecm + zqlcm)/(zecm - zqlcm)); //rapidity of Z in CM

  //boost from CM to lab
  double yboost;
  yboost = y - zycm;

  //Compute recoil four momentum in CM (massless recoil with mjj = 0)
  double rqtcm, rqlcm, recm, rycm;
  rqtcm = zqtcm;
  rqlcm = zqlcm;
  recm = sqrt(rqtcm*rqtcm+rqlcm*rqlcm+mjj*mjj);
  rycm = 0.5*log((recm-rqlcm)/(recm+rqlcm));

  //now apply yboost to the recoil to have the longitudinal momentum of the recoil in the lab frame
  double rqt, ry, rql, re;
  rqt = rqtcm;
  ry = rycm+yboost;
  rql = 0.5*sqrt(rqtcm*rqtcm+mjj*mjj)*(exp(ry)-exp(-ry));
  re = sqrt(rql*rql+rqt*rqt+mjj*mjj);
  
  //compute tau = mtot^2/s
  double mtot, tau;
  mtot = sqrt(pow(ze+re,2)-pow(zql+rql,2));
  tau = mtot*mtot/s;

  /*
  //cross check that rapidity of the system is the rapidity of the boost
  double ytot;
  ytot = 0.5* log((ze+re + zql+rql)/(ze+re - (zql+rql)));
  cout << ytot << "  " << yboost << endl;
  */

  //convert invariant mass to 0-1 for Breit-Wigner weighting
  double m2, m1, s3max, almin, almax, tanal, al, x1;
  m2 = mjj;
  m1 = mtot;
  s3max=(m2-m1)*(m2-m1);
  almin=atan((0.-rmass*rmass)/rmass/rwidth);
  almax=atan((s3max-rmass*rmass)/rmass/rwidth);
  tanal = (m*m - rmass*rmass)/rmass/rwidth;
  al = atan(tanal);
  x1 = (al-almin)/(almax-almin);

  //************** REAL RADIATION ***************
  //phase space mapping
  rre[1] = x1;                                     //mll
  rre[3] = phicm;                                  //phi of Z in CM frame
  rre[6] = (cos_th + 1.)/2.;                       //cos_th of dilepton in Z rest frame
  rre[7] = phiZ;                                   //phi of dilepton in Z rest frame
  rre[8] = log(tau)/log(minmass*minmass/s);        //tau of the system
  rre[9] = 0.5 - yboost/log(tau);                  //y of the system
  //dimension to be integrated
  rre[2] = (zcth + 1.)/2.;                         //cos th of Z in CM frame
  //dimensions of real radiation to be integrated
  rre[0] = mjj*mjj / (mtot*mtot);                  //mjj
  rre[4] = phijj;                                  //phi jj
  rre[5] = costhjj;                                //costh jj
  begin_time = clock_real();
  double f[opts.totpdf];
  value = realint_(rre,wgt,f);
  end_time = clock_real();
  cout << "Real: " << value << "  " << "time " << end_time - begin_time << "s" << endl;
  //******************************************

  return value;
}

double dyvirt(double m, double y, double qt, double phicm, double phiZ, double cos_th, double zcth, double vz)
{
  //set up constants
  double rmass = opts.rmass;
  double rwidth = opts.rwidth;
  double minmass = sqrt(taumin_.taumin_)*energy_.sroot_;
  double sqrts = energy_.sroot_;;
  double s = sqrts*sqrts;
  double mmin = opts.mlow;
  double mmax = opts.mhigh;
  double xqtcut = qtcut_.xqtcut_;

  //  std::cout << std::setprecision(15);
  double begin_time, end_time;
  double value;
  double rvi[22];
  double wgt = 1;

  /*
  //variables to be integrated (4 dimensions in total)
  //1 collinear PDF dimension
  double zcth;
  zcth = 0.5;   //polar angle of Z in the CM frame
  //1 dimension of virtual
  double vz;
  vz = 0.;
  */
  //invariant mass of the recoil (set to 0 to have the correct virtual phase space mapping)
  double mjj = 0.;
  
  //evaluate four momentum of Z in lab frame, given y, m, qt
  double zql, zqt, zp, ze;
  zql = 0.5*sqrt(m*m+qt*qt)*(exp(y)-exp(-y));
  zqt = qt;
  zp = sqrt(zqt*zqt + zql*zql);
  ze = sqrt(m*m+zqt*zqt+zql*zql);

  //Compute Z four momentum in CM, given y, m, qt and zcth
  double zqtcm, zpcm, zqlcm, zecm, zycm;
  zqtcm = qt;
  zpcm = zqtcm / sin(acos(zcth)); //absolute momentum of Z in CM
  zqlcm = zpcm * zcth;         //longitudinal momentum of Z in CM
  zecm = sqrt(m*m + zpcm*zpcm);//Energy of Z in CM
  zycm = 0.5 * log((zecm + zqlcm)/(zecm - zqlcm)); //rapidity of Z in CM

  //boost from CM to lab
  double yboost;
  yboost = y - zycm;

  //Compute recoil four momentum in CM (massless recoil with mjj = 0)
  double rqtcm, rqlcm, recm, rycm;
  rqtcm = zqtcm;
  rqlcm = zqlcm;
  recm = sqrt(rqtcm*rqtcm+rqlcm*rqlcm+mjj*mjj);
  rycm = 0.5*log((recm-rqlcm)/(recm+rqlcm));

  //now apply yboost to the recoil to have the longitudinal momentum of the recoil in the lab frame
  double rqt, ry, rql, re;
  rqt = rqtcm;
  ry = rycm+yboost;
  rql = 0.5*sqrt(rqtcm*rqtcm+mjj*mjj)*(exp(ry)-exp(-ry));
  re = sqrt(rql*rql+rqt*rqt+mjj*mjj);
  
  //compute tau = mtot^2/s
  double mtot, tau;
  mtot = sqrt(pow(ze+re,2)-pow(zql+rql,2));
  tau = mtot*mtot/s;

  /*
  //cross check that rapidity of the system is the rapidity of the boost
  double ytot;
  ytot = 0.5* log((ze+re + zql+rql)/(ze+re - (zql+rql)));
  cout << ytot << "  " << yboost << endl;
  */

  //convert invariant mass to 0-1 for Breit-Wigner weighting
  double m2, m1, s3max, almin, almax, tanal, al, x1;
  m2 = mjj;
  m1 = mtot;
  s3max=(m2-m1)*(m2-m1);
  almin=atan((0.-rmass*rmass)/rmass/rwidth);
  almax=atan((s3max-rmass*rmass)/rmass/rwidth);
  tanal = (m*m - rmass*rmass)/rmass/rwidth;
  al = atan(tanal);
  x1 = (al-almin)/(almax-almin);

  //************** VIRTUAL ***************
  //phase space mapping
  rvi[0] = x1;                                       //mll
  rvi[2] = phicm;                                    //phi of Z in CM frame
  rvi[3] = (cos_th + 1.)/2.;                         //cos_th of dilepton in Z rest frame
  rvi[4] = phiZ;                                     //phi of dilepton in Z rest frame
  rvi[5] = log(tau)/log(minmass*minmass/s);          //tau of the system
  rvi[6] = 0.5 - yboost/log(tau);                    //y of the system
  //dimension to be integrated
  rvi[1] = (zcth + 1.)/2.;                           //cos th of Z in CM frame
  //virtual dimension to be integrated
  rvi[9] = vz;
  begin_time = clock_real();
  double f[opts.totpdf];
  value = virtint_(rvi,wgt,f);
  end_time = clock_real();
  cout << "Virt: " << value << "  " << "time " << end_time - begin_time << "s" << endl;
  //******************************************
}
double dylow(double m, double y, double qt, double phicm, double phiZ, double cos_th, double zcth)
{
  //set up constants
  double rmass = opts.rmass;
  double rwidth = opts.rwidth;
  double minmass = sqrt(taumin_.taumin_)*energy_.sroot_;
  double sqrts = energy_.sroot_;;
  double s = sqrts*sqrts;
  double mmin = opts.mlow;
  double mmax = opts.mhigh;
  double xqtcut = qtcut_.xqtcut_;

  //  std::cout << std::setprecision(15);
  double begin_time, end_time;
  double value;
  double rlo[22];
  double wgt = 1;

  /*
  //variables to be integrated (1 dimension in total)
  //1 collinear PDF dimension
  double zcth;
  zcth = 0.5;   //polar angle of Z in the CM frame
  */
  //invariant mass of the recoil (set to 0 to have the correct LO phase space mapping)
  double mjj = 0.;
  
  //evaluate four momentum of Z in lab frame, given y, m, qt
  double zql, zqt, zp, ze;
  zql = 0.5*sqrt(m*m+qt*qt)*(exp(y)-exp(-y));
  zqt = qt;
  zp = sqrt(zqt*zqt + zql*zql);
  ze = sqrt(m*m+zqt*zqt+zql*zql);

  //Compute Z four momentum in CM, given y, m, qt and zcth
  double zqtcm, zpcm, zqlcm, zecm, zycm;
  zqtcm = qt;
  zpcm = zqtcm / sin(acos(zcth)); //absolute momentum of Z in CM
  zqlcm = zpcm * zcth;         //longitudinal momentum of Z in CM
  zecm = sqrt(m*m + zpcm*zpcm);//Energy of Z in CM
  zycm = 0.5 * log((zecm + zqlcm)/(zecm - zqlcm)); //rapidity of Z in CM

  //boost from CM to lab
  double yboost;
  yboost = y - zycm;

  //Compute recoil four momentum in CM (massless recoil with mjj = 0)
  double rqtcm, rqlcm, recm, rycm;
  rqtcm = zqtcm;
  rqlcm = zqlcm;
  recm = sqrt(rqtcm*rqtcm+rqlcm*rqlcm+mjj*mjj);
  rycm = 0.5*log((recm-rqlcm)/(recm+rqlcm));

  //now apply yboost to the recoil to have the longitudinal momentum of the recoil in the lab frame
  double rqt, ry, rql, re;
  rqt = rqtcm;
  ry = rycm+yboost;
  rql = 0.5*sqrt(rqtcm*rqtcm+mjj*mjj)*(exp(ry)-exp(-ry));
  re = sqrt(rql*rql+rqt*rqt+mjj*mjj);
  
  //compute tau = mtot^2/s
  double mtot, tau;
  mtot = sqrt(pow(ze+re,2)-pow(zql+rql,2));
  tau = mtot*mtot/s;

  /*
  //cross check that rapidity of the system is the rapidity of the boost
  double ytot;
  ytot = 0.5* log((ze+re + zql+rql)/(ze+re - (zql+rql)));
  cout << ytot << "  " << yboost << endl;
  */

  //convert invariant mass to 0-1 for Breit-Wigner weighting
  double m2, m1, s3max, almin, almax, tanal, al, x1;
  m2 = mjj;
  m1 = mtot;
  s3max=(m2-m1)*(m2-m1);
  almin=atan((0.-rmass*rmass)/rmass/rwidth);
  almax=atan((s3max-rmass*rmass)/rmass/rwidth);
  tanal = (m*m - rmass*rmass)/rmass/rwidth;
  al = atan(tanal);
  x1 = (al-almin)/(almax-almin);

  //************** VIRTUAL ***************
  //phase space mapping
  rlo[0] = x1;                                       //mll
  rlo[2] = phicm;                                    //phi of Z in CM frame
  rlo[3] = (cos_th + 1.)/2.;                         //cos_th of dilepton in Z rest frame
  rlo[4] = phiZ;                                     //phi of dilepton in Z rest frame
  rlo[5] = log(tau)/log(minmass*minmass/s);          //tau of the system
  rlo[6] = 0.5 - yboost/log(tau);                    //y of the system
  //dimension to be integrated
  rlo[1] = (zcth + 1.)/2.;                           //cos th of Z in CM frame
  begin_time = clock_real();
  double f[opts.totpdf];
  value = lowint_(rlo,wgt,f);
  end_time = clock_real();
  cout << "Virt: " << value << "  " << "time " << end_time - begin_time << "s" << endl;
  //******************************************
}

double dyct(double m, double y, double qt, double phicm, double phiZ, double cos_th, double alpha, double beta)
{
  //set up constants
  double rmass = opts.rmass; //91.1876;
  double rwidth = opts.rwidth; //2.495;
  double minmass = sqrt(taumin_.taumin_)*energy_.sroot_; //40.;
  double sqrts = energy_.sroot_; //7000.;
  double s = sqrts*sqrts;
  double mmin = opts.mlow; //66.;
  double mmax = opts.mhigh; //116.;
  double xqtcut = qtcut_.xqtcut_; //0.008;

  //  std::cout << std::setprecision(15);
  double begin_time, end_time;
  double value;
  double rct[22];
  double wgt = 1;

  /*
  //variables to be integrated
  //2 dimensions for the counterterm
  double alpha,beta;
  beta = 0.1;
  alpha = 0.1;
  */
  
  //phase space for the counterterm
  double qtcut,xth,tau0,phict;
  qtcut = xqtcut*m;
  xth = 1/(log((qt*qt)/(qtcut*qtcut))+1);
  tau0 = m*m/s;
  phict = -phicm + 0.25;
  if (phict < 0.)
    phict = phict+1.;

  //************** COUNTERTERM ***************
  //phase space mapping
  rct[5] = (m*m - mmin*mmin)/(mmax*mmax-mmin*mmin);     //q2 (why not breit wigner?)
  rct[2] = xth;                                         //costh of the system (qt)
  rct[6] = 0.5 - y/log(tau0);                           //y of the Z boson
  rct[0] = phict;                                       //phi of the system
  rct[3] = (cos_th + 1.)/2.;                            //costh of dilepton in Z rest frame
  rct[4] = phiZ;                                        //phi of dilepton in Z rest frame
  //dummy dimension
  rct[1] = 0.;
  //dimensions to be integrated (alpha and beta of the counterterm)
  rct[7] = beta;
  rct[8] = alpha;
  begin_time = clock_real();
  value = countint_(rct,wgt);
  end_time = clock_real();
  cout << "Counterterm: " << value << "  " << "time " << end_time - begin_time << "s" << endl;
  //******************************************

  return value;
}

#ifdef HAVE_DYRES
void dyres(double costh,double m,double qt,double y)
{
  double rmass = opts.rmass;
  double rwidth = opts.rwidth;
  double sqrts = energy_.sroot_; //7000.;
  double s = sqrts*sqrts;

  double begin_time, end_time;
  double value;
  double rr[22];
  double wgt = 1;
  
  //convert invariant mass to 0-1 for Breit-Wigner weighting
  double mmin2 = pow(opts.mlow,2);
  double mmax2 = pow(opts.mhigh,2);
  double almin, almax, tanal, al, x1;
  almin=atan((mmin2-rmass*rmass)/rmass/rwidth);
  almax=atan((mmax2-rmass*rmass)/rmass/rwidth);
  tanal = (m*m - rmass*rmass)/rmass/rwidth;
  al = atan(tanal);
  x1 = (al-almin)/(almax-almin);

  //rapidity
  double ymax=0.5*log(s/m/m);
  double x2 = (y + ymax)/(2.*ymax);

  //qt
  double qtmin=0.1;
  double cosh2y=pow((exp(y)+exp(-y))*0.5,2);
  double qtmax = sqrt(pow(s+m*m,2)/(4.*s*cosh2y) - m*m);
  double x3 = (qt-qtmin)/qtmax;

  //costh
  double x4 = (costh+1.)/2.;
  
  //phase space mapping
  const int ncomp = 1;
  const int ndim = 6; //4;
  double f[ncomp];
  double x[ndim];

  //************** COUNTERTERM ***************
  //phase space mapping
  rr[0] = x1;                                           //q2
  rr[1] = x2;                                          //y of the Z boson
  rr[2] = x3;                                         //costh of the system (qt)
  rr[3] = x4;                                          //costh of dilepton in Z rest frame
  begin_time = clock_real();
  value = lowinthst_(rr,wgt);
  end_time = clock_real();
  cout << "Resummed part: " << value << "  " << "time " << end_time - begin_time << "s" << endl;
  //******************************************
  
  return;
}
#endif
