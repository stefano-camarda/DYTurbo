#ifndef intde2_c_h
#define intde2_c_h

extern "C"
{
 void intdeini(int lenaw, double tiny, double eps, double *aw);
 void intde(double (*f)(double), double a, double b, double *aw, 
	   double *i, double *err);
 void intdeiini(int lenaw, double tiny, double eps, double *aw);
 void intdei(double (*f)(double), double a, double *aw, double *i, 
	    double *err);
 void intdeoini(int lenaw, double tiny, double eps, double *aw);
 void intdeo(double (*f)(double), double a, double omega, double *aw, 
	    double *i, double *err);
}
#endif
