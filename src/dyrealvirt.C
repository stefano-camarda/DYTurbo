#include <iostream>
#include <LHAPDF/LHAPDF.h>
#include <dyfin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <ctime>

using namespace std;

double dyreal(double m, double y, double qt, double phicm, double phiZ, double cos_th, double zcth, double mjj, double phijj, double costhjj)
{
  //set up constants
  double rmass = 91.1876;
  double rwidth = 2.495;
  double minmass = 40.;
  double sqrts = 7000.;
  double s = sqrts*sqrts;
  double mmin = 66.;
  double mmax = 116.;
  double xqtcut = 0.008;

  //  std::cout << std::setprecision(15);
  clock_t begin_time, end_time;
  double value;
  double rre[22];
  double rvi[22];
  double rct[22];
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
  begin_time = clock();
  value = realint_(rre,wgt);
  end_time = clock();
  cout << "Real: " << value << "  " << "time " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << endl;
  //******************************************

  return value;
}

double dyvirt(double m, double y, double qt, double phicm, double phiZ, double cos_th, double zcth, double vz)
{
  //set up constants
  double rmass = 91.1876;
  double rwidth = 2.495;
  double minmass = 40.;
  double sqrts = 7000.;
  double s = sqrts*sqrts;
  double mmin = 66.;
  double mmax = 116.;
  double xqtcut = 0.008;

  //  std::cout << std::setprecision(15);
  clock_t begin_time, end_time;
  double value;
  double rre[22];
  double rvi[22];
  double rct[22];
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
  begin_time = clock();
  value = virtint_(rvi,wgt);
  end_time = clock();
  cout << "Virt: " << value << "  " << "time " << float( end_time - begin_time ) /  CLOCKS_PER_SEC << endl;
  //******************************************
}
