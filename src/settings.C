#include <math.h>

#include "settings.h"
#include "interface.h"

settings opts;
binning bins;

void settings::init()
{
  //read input settings from file

  
  //Z
  //resonance mass and width (used for breit wigner unweighting)
  rmass = 91.1876;
  rwidth = 2.495;
  //Set boundaries of integration
  ylow = 2;
  yhigh = 2.4;
  mlow = 66.;
  mhigh = 116.;

  /*
  //W
  //resonance mass and width (used for breit wigner unweighting)
  rmass = 80.385;
  rwidth = 2.091;
  //Set boundaries of integration
  ylow = -5.;
  yhigh = 5.;
  mlow = 10.;
  mhigh = 1000.;
  */
  
  //type of integration
  int2d = false;
  int3d = true;
  int4d = false;
  
  //Cuba settings
  cubaverbosity = 0;   //Cuba info messsages, from 0 to 3
  niter = 0;           //only for 2d and 3d cuhre integration
  vegasncalls = 10000; //only for 4d vegas integration
  
  //total or with lepton cuts
  makelepcuts = true;

  //integration types and settings for costh phi_lep phase space
  cubaint = true;    //integration with Cuba Suave
  trapezint = false; //trapezoidal rule for the phi_lep integration and semi-analytical for costh
  quadint  = false;  //quadrature rule for the phi_lep integration and semi-analytical for costh
  
  suavepoints = 1000000;  //number of points for suave integration, newpoints is set to suavepoints/10;
  nphitrape = 1000;       //number of steps for trapezoidal rule of phi_lep integration
  quadnphi = 20;          //number of segments for quadrature rule of phi_lep integration
  ncstart = 1000;         //starting sampling for the costh semi-analytical integration (common settings for the trapezoidal and quadrature rules)

  //qt-recoil prescriptions
  qtrec_naive = false;
  qtrec_cs = false;
  qtrec_kt0 = true;
  
  //debug settings
  timeprofile = false; //debug and time profile resummation integration
  verbose = false; //debug and time profile costh phi_lep integration

  //Set to 1 to use the dyres approximation of PDFs and integration contour in the complex plane for the Mellin inversion
  //Set to 0 to use exact PDFs and straight line contour in the complex plane for the Mellin inversion
  opts_.approxpdf_ = 0;
}

bool cuts(double p3[4], double p4[4])
{
  if (!opts.makelepcuts)
    return true;
  double pt3 = sqrt((float)pow(p3[0],2)+pow(p3[1],2));
  if (pt3 < 20)
    return false;
  double pt4 = sqrt((float)pow(p4[0],2)+pow(p4[1],2));
  if (pt4 < 20)
    return false;
  double y3 = 0.5 *log((p3[3] + p3[2]) / (p3[3] - p3[2]));
  if (fabs(y3) > 2.4)
    return false;
  double y4 = 0.5 *log((p4[3] + p4[2]) / (p4[3] - p4[2]));
  if (fabs(y4) > 2.4)
    return false;

  //  if (y3-y4 < 0.2)
  //    return false;

  //  if (66 < m < 116
  //      && pt3 > 20 && pt4 > 20
  //      && fabs (y3) < 2.4 && fabs(y4) < 2.4)
  //    cut = true;
  return true;
}

void binning::init()
{
  qtbins.push_back(0  ); 
  qtbins.push_back(2  ); 
  qtbins.push_back(4  ); 
  qtbins.push_back(6  ); 
  qtbins.push_back(8  ); 
  qtbins.push_back(10 ); 
  qtbins.push_back(12 ); 
  qtbins.push_back(14 ); 
  qtbins.push_back(16 ); 
  qtbins.push_back(18 ); 
  qtbins.push_back(22 ); 
  qtbins.push_back(26 ); 
  qtbins.push_back(30 ); 
  qtbins.push_back(34 ); 
  qtbins.push_back(38 ); 
  qtbins.push_back(42 ); 
  qtbins.push_back(46 ); 
  qtbins.push_back(50 ); 
  qtbins.push_back(54 ); 
  qtbins.push_back(60 ); 
  qtbins.push_back(70 ); 
  qtbins.push_back(80 ); 
  qtbins.push_back(100); 
  qtbins.push_back(150); 
  qtbins.push_back(200); 
  qtbins.push_back(300);
  qtbins.push_back(800);
  /*
  qtbins.push_back(0  );       qtbins.push_back(0.5  );  
  qtbins.push_back(1  );       qtbins.push_back(1.5  );  
  qtbins.push_back(2  );       qtbins.push_back(2.5  );  
  qtbins.push_back(3  );       qtbins.push_back(3.5  );  
  qtbins.push_back(4  );       qtbins.push_back(4.5  );  
  qtbins.push_back(5  );       qtbins.push_back(5.5  );  
  qtbins.push_back(6  );       qtbins.push_back(6.5  );  
  qtbins.push_back(7  );       qtbins.push_back(7.5  );  
  qtbins.push_back(8  );       qtbins.push_back(8.5  );  
  qtbins.push_back(9  );       qtbins.push_back(9.5  );  
  qtbins.push_back(10 );       qtbins.push_back(10.5 ); 
  qtbins.push_back(11 );       qtbins.push_back(11.5 ); 
  qtbins.push_back(12 );       qtbins.push_back(12.5 ); 
  qtbins.push_back(13 );       qtbins.push_back(13.5 ); 
  qtbins.push_back(14 );       qtbins.push_back(14.5 ); 
  qtbins.push_back(15 );       qtbins.push_back(15.5 ); 
  qtbins.push_back(16 );       qtbins.push_back(16.5 ); 
  qtbins.push_back(17 );       qtbins.push_back(17.5 ); 
  qtbins.push_back(18 );       qtbins.push_back(18.5 ); 
  qtbins.push_back(19 );       qtbins.push_back(19.5 ); 
  qtbins.push_back(20 );       qtbins.push_back(20.5 ); 
  qtbins.push_back(21 );       qtbins.push_back(21.5 ); 
  qtbins.push_back(22 );       qtbins.push_back(22.5 ); 
  qtbins.push_back(23 );       qtbins.push_back(23.5 ); 
  qtbins.push_back(24 );       qtbins.push_back(24.5 ); 
  qtbins.push_back(25 );       qtbins.push_back(25.5 ); 
  qtbins.push_back(26 );       qtbins.push_back(26.5 ); 
  qtbins.push_back(27 );       qtbins.push_back(27.5 ); 
  qtbins.push_back(28 );       qtbins.push_back(28.5 ); 
  qtbins.push_back(29 );       qtbins.push_back(29.5 ); 
  qtbins.push_back(30 );       qtbins.push_back(30.5 ); 
  qtbins.push_back(31 );       qtbins.push_back(31.5 ); 
  qtbins.push_back(32 );       qtbins.push_back(32.5 ); 
  qtbins.push_back(33 );       qtbins.push_back(33.5 ); 
  qtbins.push_back(34 );       qtbins.push_back(34.5 ); 
  qtbins.push_back(35 );       qtbins.push_back(35.5 ); 
  qtbins.push_back(36 );       qtbins.push_back(36.5 ); 
  qtbins.push_back(37 );       qtbins.push_back(37.5 ); 
  qtbins.push_back(38 );       qtbins.push_back(38.5 ); 
  qtbins.push_back(39 );       qtbins.push_back(39.5 ); 
  qtbins.push_back(40 );       qtbins.push_back(40.5 ); 
  qtbins.push_back(41 );       qtbins.push_back(41.5 ); 
  qtbins.push_back(42 );       qtbins.push_back(42.5 ); 
  qtbins.push_back(43 );       qtbins.push_back(43.5 ); 
  qtbins.push_back(44 );       qtbins.push_back(44.5 ); 
  qtbins.push_back(45 );       qtbins.push_back(45.5 ); 
  qtbins.push_back(46 );       qtbins.push_back(46.5 ); 
  qtbins.push_back(47 );       qtbins.push_back(47.5 ); 
  qtbins.push_back(48 );       qtbins.push_back(48.5 ); 
  qtbins.push_back(49 );       qtbins.push_back(49.5 ); 
  qtbins.push_back(50 );       qtbins.push_back(50.5 ); 
  */
}
