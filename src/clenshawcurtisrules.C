#include "clenshawcurtisrules.h"
//nodes and weights for the Clenshaw Curtis quadrature rules
//from https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_clenshaw_curtis/quadrature_rules_clenshaw_curtis.html
#include "sandia_rules.hpp"
#include "settings.h"
#include <string.h>
#include <iostream>
#include <math.h>
#include <quadmath.h>

using namespace std;
double cc::xxx[CCNMAX][CCNMAX];
double cc::www[CCNMAX][CCNMAX];
double cc::cosw[65][65];
double cc::sinw[65][65];

const double cc::xxx9[9]={-1.0000000000000000,
			  -0.9238795325112867,    
			  -0.7071067811865475,    
			  -0.3826834323650897,    
			  0.000000000000000,    
			  0.3826834323650898,    
			  0.7071067811865475,    
			  0.9238795325112867,    
			  1.0000000000000000};
const double cc::www9[9]={0.1587301587301588E-01,
			  0.1462186492160182,    
			  0.2793650793650794,    
			  0.3617178587204897,    
			  0.3936507936507936,    
			  0.3617178587204898,    
			  0.2793650793650794,    
			  0.1462186492160182,    
			  0.1587301587301588E-01};

const double cc::xxx33[33]={-1.0000000000000000,
			    -0.9951847266721968,    
			    -0.9807852804032304,    
			    -0.9569403357322088,    
			    -0.9238795325112867,    
			    -0.8819212643483549,    
			    -0.8314696123025453,    
			    -0.7730104533627370,    
			    -0.7071067811865475,    
			    -0.6343932841636454,    
			    -0.5555702330196020,    
			    -0.4713967368259977,    
			    -0.3826834323650897,    
			    -0.2902846772544622,    
			    -0.1950903220161282,    
			    -0.9801714032956065E-01,
			    0.000000000000000,
			    0.9801714032956077E-01,
			    0.1950903220161283,    
			    0.2902846772544623,    
			    0.3826834323650898,    
			    0.4713967368259978,    
			    0.5555702330196023,    
			    0.6343932841636455,    
			    0.7071067811865475,    
			    0.7730104533627370,    
			    0.8314696123025452,    
			    0.8819212643483550,    
			    0.9238795325112867,    
			    0.9569403357322088,    
			    0.9807852804032304,    
			    0.9951847266721969,    
			    1.0000000000000000};
const double cc::www33[33]={0.9775171065493659E-03,
			    0.9393197962955013E-02,
			    0.1923424513268114E-01,
			    0.2845791667723369E-01,
			    0.3759434191404722E-01,
			    0.4626276283775175E-01,
			    0.5455501630398032E-01,
			    0.6227210954529399E-01,
			    0.6942757563043547E-01,
			    0.7588380044138848E-01,
			    0.8163481765493850E-01,
			    0.8657753844182743E-01,
			    0.9070611286772098E-01,
			    0.9394324443876872E-01,
			    0.9629232594548820E-01,
			    0.9769818820805558E-01,
			    0.9817857778176831E-01,
			    0.9769818820805558E-01,
			    0.9629232594548819E-01,
			    0.9394324443876871E-01,
			    0.9070611286772098E-01,
			    0.8657753844182743E-01,
			    0.8163481765493850E-01,
			    0.7588380044138850E-01,
			    0.6942757563043547E-01,
			    0.6227210954529399E-01,
			    0.5455501630398031E-01,
			    0.4626276283775177E-01,
			    0.3759434191404721E-01,
			    0.2845791667723370E-01,
			    0.1923424513268119E-01,
			    0.9393197962955048E-02,
			    0.9775171065493659E-03};

const double cc::xxx65[65]={-1.0000000000000000,
			    -0.9987954562051724,    
			    -0.9951847266721968,    
			    -0.9891765099647810,    
			    -0.9807852804032304,    
			    -0.9700312531945440,    
			    -0.9569403357322088,    
			    -0.9415440651830207,    
			    -0.9238795325112867,    
			    -0.9039892931234433,    
			    -0.8819212643483549,    
			    -0.8577286100002720,    
			    -0.8314696123025453,    
			    -0.8032075314806448,    
			    -0.7730104533627370,    
			    -0.7409511253549589,    
			    -0.7071067811865475,    
			    -0.6715589548470184,    
			    -0.6343932841636454,    
			    -0.5956993044924334,    
			    -0.5555702330196020,    
			    -0.5141027441932217,    
			    -0.4713967368259977,    
			    -0.4275550934302819,    
			    -0.3826834323650897,    
			    -0.3368898533922199,    
			    -0.2902846772544622,    
			    -0.2429801799032639,    
			    -0.1950903220161282,    
			    -0.1467304744553616,    
			    -0.9801714032956065E-01,
			    -0.4906767432741801E-01,
			    0.000000000000000,    
			    0.4906767432741813E-01,
			    0.9801714032956077E-01,
			    0.1467304744553617,    
			    0.1950903220161283,    
			    0.2429801799032640,    
			    0.2902846772544623,    
			    0.3368898533922201,    
			    0.3826834323650898,    
			    0.4275550934302822,    
			    0.4713967368259978,    
			    0.5141027441932218,    
			    0.5555702330196023,    
			    0.5956993044924335,    
			    0.6343932841636455,    
			    0.6715589548470184,    
			    0.7071067811865475,    
			    0.7409511253549592,    
			    0.7730104533627370,    
			    0.8032075314806448,    
			    0.8314696123025452,    
			    0.8577286100002721,    
			    0.8819212643483550,    
			    0.9039892931234433,    
			    0.9238795325112867,    
			    0.9415440651830208,    
			    0.9569403357322088,    
			    0.9700312531945440,    
			    0.9807852804032304,    
			    0.9891765099647810,    
			    0.9951847266721969,    
			    0.9987954562051724,    
			    1.0000000000000000};
const double cc::www65[65]={0.2442002442002449E-03,
			    0.2351490675311702E-02,
			    0.4831465448790911E-02,
			    0.7192693161736115E-02,
			    0.9582338795283791E-02,
			    0.1192339471421277E-01,
			    0.1425206043235199E-01,
			    0.1653498765728959E-01,
			    0.1878652974179578E-01,
			    0.2098627442973744E-01,
			    0.2314069493435819E-01,
			    0.2523506498175476E-01,
			    0.2727225714146840E-01,
			    0.2924065319746835E-01,
			    0.3114129710406762E-01,
			    0.3296454656997634E-01,
			    0.3471049818092511E-01,
			    0.3637092028663919E-01,
			    0.3794545992128482E-01,
			    0.3942698871295609E-01,
			    0.4081501340035782E-01,
			    0.4210333111141810E-01,
			    0.4329151496169082E-01,
			    0.4437417923925731E-01,
			    0.4535110955166067E-01,
			    0.4621766751092559E-01,
			    0.4697395904661415E-01,
			    0.4761604458525018E-01,
			    0.4814443257251221E-01,
			    0.4855584485714105E-01,
			    0.4885125664306610E-01,
			    0.4902801843102554E-01,
			    0.4908762351494248E-01,
			    0.4902801843102556E-01,
			    0.4885125664306610E-01,
			    0.4855584485714105E-01,
			    0.4814443257251221E-01,
			    0.4761604458525019E-01,
			    0.4697395904661414E-01,
			    0.4621766751092559E-01,
			    0.4535110955166067E-01,
			    0.4437417923925733E-01,
			    0.4329151496169082E-01,
			    0.4210333111141811E-01,
			    0.4081501340035782E-01,
			    0.3942698871295609E-01,
			    0.3794545992128483E-01,
			    0.3637092028663919E-01,
			    0.3471049818092511E-01,
			    0.3296454656997635E-01,
			    0.3114129710406762E-01,
			    0.2924065319746836E-01,
			    0.2727225714146839E-01,
			    0.2523506498175477E-01,
			    0.2314069493435821E-01,
			    0.2098627442973743E-01,
			    0.1878652974179578E-01,
			    0.1653498765728961E-01,
			    0.1425206043235200E-01,
			    0.1192339471421278E-01,
			    0.9582338795283809E-02,
			    0.7192693161736120E-02,
			    0.4831465448790926E-02,
			    0.2351490675311698E-02,
			    0.2442002442002449E-03};

void cc::init()
{
  /*
  for (int i = 0; i < 65; i++)
    for (int j = 0; j < 65; j++)
      {
	xxx[i][j] = 0;
	www[i][j] = 0;
      }
  for (int i = 0; i < 65; i++)
    {
      if (i < 9) xxx[9-1][i] = xxx9[i];
      if (i < 33) xxx[33-1][i] = xxx33[i];
      if (i < 65) xxx[65-1][i] = xxx65[i];

      if (i < 9) www[9-1][i] = www9[i];
      if (i < 33) www[33-1][i] = www33[i];
      if (i < 65) www[65-1][i] = www65[i];
    }
  */

  //fstream iff;
  //iff = fstream(string(SHAREDIR)+"/ccr.bin", ios::in | ios::binary);
  fstream iff(string(SHAREDIR)+"/ccr.bin", ios::in | ios::binary);
  if (iff)
    {
      iff.read((char*)&xxx,CCNMAX*CCNMAX*sizeof(double));
      iff.read((char*)&www,CCNMAX*CCNMAX*sizeof(double));
      iff.close();
      return;
    }
 
  //generation of weights and nodes
  memset(xxx, 0, sizeof(xxx));
  memset(www, 0, sizeof(xxx));
  for (int n = 1; n <= CCNMAX; n++)
    {
      double *w = new double[n];
      double *x = new double[n];
      webbur::clenshaw_curtis_compute(n, x, w);
      for ( int i = 1; i <= n; i++ )
	{
	  xxx[n-1][i-1] = x[i-1];
	  www[n-1][i-1] = w[i-1];
	}
      delete[] w;
      delete[] x;      
    }

  //fstream off;
  //off = fstream(string(SHAREDIR)+"/ccr.bin", ios::out | ios::binary);
  fstream off(string(SHAREDIR)+"/ccr.bin", ios::out | ios::binary);
  off.write((char*)&xxx,CCNMAX*CCNMAX*sizeof(double));
  off.write((char*)&www,CCNMAX*CCNMAX*sizeof(double));
  off.close();
  
  /*
  //calculate nodes and weights and cross check result
  int n = opts.mellinrule;
  int N = n-1;
    {
      double *w = new double[n];
      double *x = new double[n];
      webbur::clenshaw_curtis_compute(n, x, w);
      for ( int j = 0; j <= N; j++ )
	{
	  cout << "cc   nodes and weights: " << j << "  " << x[j] << "  " << w[j] << endl;

	  //calculate
	  double x = cos(M_PI*double(N-j)/double(N));

	  double sumjk = 0;
	  for ( int k = 0; k <= N; k++ )
	    {
	      double jktk;
	      if (k == 1)
		jktk = 0.;
	      else if (k == 0 || k == N)
		jktk = (-cos(M_PI*k)-1.)/(pow(k,2) - 1)/2.*cos(k*M_PI*j/double(N));
	      else
		jktk = (-cos(M_PI*k)-1.)/(pow(k,2) - 1)*cos(k*M_PI*j/double(N));
	      sumjk += jktk;
	    }

	    
	  double w = 2./double(N) * sumjk;
	  if (j == 0 || j == N)
	    w/= 2.;
	  cout << "calc nodes and weights: " << j << "  " << x << "  " << w << endl;
	}
    }
    cout << endl;
    setw(1e-2);
  */
}

//tabulate jkcoswx = int_0^pi cos(k*t)*cos(w*cos(t))*sin(t) dt
double jkcoswx(int k, double W)
{
  __float128 res;
  __float128 w = W;
  if (k%2 == 1) return 0.;
  if (k == 0) res =  2*sinq(w)/w;
  if (k == 2) res =  2*(4*w*cosq(w) + (powq(w,2) - 4)*sinq(w))/powq(w,3);
  if (k == 4) res =  2*(16*(powq(w,3) - 12*w)*cosq(w) + (powq(w,4) - 80*powq(w,2) + 192)*sinq(w))/powq(w,5);
  if (k == 6) res =  2*(12*(3*powq(w,5) - 224*powq(w,3) + 1920*w)*cosq(w) + (powq(w,6) - 420*powq(w,4) + 10368*powq(w,2) - 23040)*sinq(w))/powq(w,7);
  if (k == 8) res =  2.Q*(64.Q*(powq(w,7) - 252.Q*powq(w,5) + 10560.Q*powq(w,3) - 80640.Q*w)*cosq(w) + (powq(w,8) - 1344.Q*powq(w,6) + 126720.Q*powq(w,4) - 2396160.Q*powq(w,2) + 5160960.Q)*sinq(w))/powq(w,9);
  if (k == 10) res =  2.Q*(20.Q*(5.Q*powq(w,9) - 3168.Q*powq(w,7) + 384384.Q*powq(w,5) - 12902400.Q*powq(w,3) + 92897280.Q*w)*cosq(w) + (powq(w,10) - 3300.Q*powq(w,8) + 823680.Q*powq(w,6) - 52416000.Q*powq(w,4) + 877363200.Q*powq(w,2) - 1857945600.Q)*sinq(w))/powq(w,11);
  if (k == 12) res =  2.Q*(48.Q*(3.Q*powq(w,11) - 4004.Q*powq(w,9) + 1098240.Q*powq(w,7) - 98703360.Q*powq(w,5) + 2941747200.Q*powq(w,3) - 20437401600.Q*w)*cosq(w) + (powq(w,12) - 6864.Q*powq(w,10) + 3706560.Q*powq(w,8) - 570286080.Q*powq(w,6) + 30005821440.Q*powq(w,4) - 468202291200.Q*powq(w,2) + 980995276800.Q)*sinq(w))/powq(w,13);
  if (k == 14) res =  2.Q*(28.Q*(7.Q*powq(w,13) - 17472.Q*powq(w,11) + 9335040.Q*powq(w,9) - 1786060800.Q*powq(w,7) + 135908720640.Q*powq(w,5) - 3760481894400.Q*powq(w,3) + 25505877196800.Q*w)*cosq(w) + (powq(w,14) - 12740.Q*powq(w,12) + 13069056.Q*powq(w,10) - 4063288320.Q*powq(w,8) + 490095083520.Q*powq(w,6) - 23032951603200.Q*powq(w,4) + 343348346880000.Q*powq(w,2) - 714164561510400.Q)*sinq(w))/powq(w,15);
  if (k == 16) res =  2*(256*(powq(w,15) - 4284*powq(w,13) + 4031040*powq(w,11) - 1432569600*powq(w,9) + 223278612480*powq(w,7) - 15276957696000*powq(w,5) + 401717565849600*powq(w,3) - 2678117105664000*w)*cosq(w) + (powq(w,16) - 21760*powq(w,14) + 38697984*powq(w,12) - 21670871040*powq(w,10) + 5060981882880*powq(w,8) - 526467465216000*powq(w,6) + 22955289477120000*powq(w,4) - 331372356540825600*powq(w,2) + 685597979049984000)*sinq(w))/powq(w,17);
  if (k == 18) res =  2*(108*(3*powq(w,17) - 20672.Q*powq(w,15) + 31834880.Q*powq(w,13) - 19170385920.Q*powq(w,11) + 5375225856000.Q*powq(w,9) - 733293969408000.Q*powq(w,7) + 46599237638553600.Q*powq(w,5) - 1180752075030528000.Q*powq(w,3) + 7770110429233152000.Q*w)*cosq(w) + (powq(w,18) - 34884.Q*powq(w,16) + 100465920.Q*powq(w,14) - 93455631360.Q*powq(w,12) + 37957364121600.Q*powq(w,10) - 7424601440256000.Q*powq(w,8) + 698988564578304000.Q*powq(w,6) - 28891527335903232000.Q*powq(w,4) + 407245199555690496000.Q*powq(w,2) - 839171926357180416000.Q)*sinq(w))/powq(w,19);
  if (k == 20) res =  2.Q*(80.Q*(5.Q*powq(w,19) - 52668.Q*powq(w,17) + 125520384.Q*powq(w,15) - 119814912000.Q*powq(w,13) + 55413692006400.Q*powq(w,11) - 13290953195520000.Q*powq(w,9) + 1650944419194470400.Q*powq(w,7) - 99360287113818931200.Q*powq(w,5) + 2447584785208442880000.Q*powq(w,3) - 15944266600786427904000.Q*w)*cosq(w) + (powq(w,20) - 53200.Q*powq(w,18) + 235350720.Q*powq(w,16) - 342328320000.Q*powq(w,14) + 224293515264000.Q*powq(w,12) - 74429337894912000.Q*powq(w,10) + 12898003274956800000.Q*powq(w,8) - 1129979735804215296000.Q*powq(w,6) + 44872387728821452800000.Q*powq(w,4) - 620987225504313507840000.Q*powq(w,2) + 1275541328062914232320000.Q)*sinq(w))/powq(w,21);
  double val = res;
  //cout << val << endl;
  return val;
}


//tabulate jksinwx = int_0^pi cos(k*t)*sin(w*cos(t))*sin(t) dt
double jksinwx(int k, double W)
{
  __float128 res;
  __float128 w = W;
  if (k%2 == 0) return 0.;
  if (k == 1) res =  -2*(w*cosq(w) - sinq(w))/powq(w,2);
  if (k == 3) res =  -2*((powq(w,3) - 24*w)*cosq(w) - 3*(3*powq(w,2) - 8)*sinq(w))/powq(w,4);
  if (k == 5) res =  -2*((powq(w,5) - 200*powq(w,3) + 1920*w)*cosq(w) - 5*(5*powq(w,4) - 168*powq(w,2) + 384)*sinq(w))/powq(w,6);
  if (k == 7) res =  -2*((powq(w,7) - 784*powq(w,5) + 40320*powq(w,3) - 322560*w)*cosq(w) - 7*(7*powq(w,6) - 1008*powq(w,4) + 21120*powq(w,2) - 46080)*sinq(w))/powq(w,8);
  if (k == 9) res =  -2.Q*((powq(w,9) - 2160.Q*powq(w,7) + 342144.Q*powq(w,5) - 12579840.Q*powq(w,3) + 92897280.Q*w)*cosq(w) - 27.Q*(3.Q*powq(w,8) - 1232.Q*powq(w,6) + 91520.Q*powq(w,4) - 1612800.Q*powq(w,2) + 3440640.Q)*sinq(w))/powq(w,10);
  if (k == 11) res =  -2.Q*((powq(w,11) - 4840.Q*powq(w,9) + 1812096.Q*powq(w,7) - 184504320.Q*powq(w,5) + 5790597120.Q*powq(w,3) - 40874803200.Q*w)*cosq(w) - 11.Q*(11.Q*powq(w,10) - 10296.Q*powq(w,8) + 1921920.Q*powq(w,6) - 109670400.Q*powq(w,4) + 1765048320.Q*powq(w,2) - 3715891200.Q)*sinq(w))/powq(w,12);
  if (k == 13) res =  -2.Q*((powq(w,13) - 9464.Q*powq(w,11) + 7138560.Q*powq(w,9) - 1588654080.Q*powq(w,7) + 130025226240.Q*powq(w,5) - 3719607091200.Q*powq(w,3) + 25505877196800.Q*w)*cosq(w) - 13.Q*(13.Q*powq(w,12) - 24024.Q*powq(w,10) + 9335040.Q*powq(w,8) - 1250242560.Q*powq(w,6) + 61776691200.Q*powq(w,4) - 940120473600.Q*powq(w,2) + 1961990553600.Q)*sinq(w))/powq(w,14);
  if (k == 15) res =  -2*((powq(w,15) - 16800*powq(w,13) + 22913280*powq(w,11) - 9674496000*powq(w,9) + 1650320179200*powq(w,7) - 118455179673600*powq(w,5) + 3188234649600000*powq(w,3) - 21424936845312000*w)*cosq(w) - 15*(15*powq(w,14) - 49504*powq(w,12) + 35473152*powq(w,10) - 9376819200*powq(w,8) + 1041966858240*powq(w,6) - 47006023680000*powq(w,4) + 688658684313600*powq(w,2) - 1428329123020800)*sinq(w))/powq(w,16);
  //if (k == 17) res =  -2*((powq(w,17) - 27744*powq(w,15) + 63256320*powq(w,13) - 46050600960*powq(w,11) + 14339448668160*powq(w,9) - 2077666246656000*powq(w,7) + 136583972388864000*powq(w,5) - 3520831288246272000*powq(w,3) + 23310331287699456000ULL*w)*cosq(w) - 17*(17*powq(w,16) - 93024*powq(w,14) + 112869120*powq(w,12) - 52718561280*powq(w,10) + 11163930624000*powq(w,8) - 1099940954112000*powq(w,6) + 46599237638553600*powq(w,4) - 664173042204672000*powq(w,2) + 1371195958099968000)*sinq(w))/powq(w,18);
  //  if (k == 19) res =  -2*((powq(w,19) - 43320*powq(w,17) + 155536128*powq(w,15) - 182118666240*powq(w,13) + 94701706444800*powq(w,11) - 24382024482816000*powq(w,9) + 3162091125473280000*powq(w,7) - 195178318002546278400*powq(w,5) + 4871859239129186304000*powq(w,3) - 31888533201572855808000*w)*cosq(w) - 19*(19*powq(w,18) - 162792*powq(w,16) + 313800960*powq(w,14) - 239629824000*powq(w,12) + 87078658867200*powq(w,10) - 15949143834624000*powq(w,8) + 1444576366795161600*powq(w,6) - 58447227714011136000*powq(w,4) + 815861595069480960000*powq(w,2) - 1678343852714360832000)*sinq(w))/powq(w,20);

  double val = res;
  return val;
}


void cc::setw(double omega)
{
  //calculate nodes and weights with w(x) = cos(x) and w(x) = sin(x)
  int n = opts.mellinrule;
  int N = n-1;
    {
      //double *w = new double[n];
      //double *x = new double[n];
      //webbur::clenshaw_curtis_compute(n, x, w);
      for ( int j = 0; j <= N; j++ )
	{
	  //	  cout << "cc   nodes and weights: " << j << "  " << x[j] << "  " << w[j] << endl;

	  double sumjkcosw = 0;
	  double sumjksinw = 0;
	  for ( int k = 0; k <= N; k++ )
	    {
	      double jktk;
	      if (k == 0 || k == N)
		jktk = jkcoswx(k,omega)/2.*cos(k*M_PI*j/double(N));
	      else
		jktk = jkcoswx(k,omega)*cos(k*M_PI*j/double(N));
	      sumjkcosw += jktk;

	      if (k == 0 || k == N)
		jktk = jksinwx(k,omega)/2.*cos(k*M_PI*j/double(N));
	      else
		jktk = jksinwx(k,omega)*cos(k*M_PI*j/double(N));
	      sumjksinw += jktk;
	    }

	    
	  cosw[N][j] = 2./double(N) * sumjkcosw;
	  if (j == 0 || j == N)
	    cosw[N][j]/= 2.;

	  sinw[N][j] = 2./double(N) * sumjksinw;
	  if (j == 0 || j == N)
	    sinw[N][j]/= 2.;
	  
	  //	  cout << "calc nodes and weights: " << j << "  " << x[j] << "  " << cosw[N][j] << "  " << sinw[N][j] << endl;
	}
      //delete[] w;
      //delete[] x;      
    }

}
