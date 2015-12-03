/* test of intde1.c */

#include <math.h>
#include <stdio.h>

#include "intde1.c"

int nfunc;

main()
{
    double f1(double), f2(double), f3(double), f4(double), 
        f5(double), f6(double);
    void intde(double (*f)(double), double a, double b, double eps, 
        double *i, double *err);
    void intdei(double (*f)(double), double a, double eps, 
        double *i, double *err);
    void intdeo(double (*f)(double), double a, double omega, double eps, 
        double *i, double *err);
    extern int nfunc;
    double i, err;
    
    nfunc = 0;
    intde(f1, 0.0, 1.0, 1.0e-15, &i, &err);
    printf("I_1=int_0^1 1/sqrt(x) dx\n");
    printf(" I_1= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
    nfunc = 0;
    intde(f2, 0.0, 2.0, 1.0e-15, &i, &err);
    printf("I_2=int_0^2 sqrt(4-x*x) dx\n");
    printf(" I_2= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
    
    nfunc = 0;
    intdei(f3, 0.0, 1.0e-15, &i, &err);
    printf("I_3=int_0^infty 1/(1+x*x) dx\n");
    printf(" I_3= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
    nfunc = 0;
    intdei(f4, 0.0, 1.0e-15, &i, &err);
    printf("I_4=int_0^infty exp(-x)/sqrt(x) dx\n");
    printf(" I_4= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
    
    nfunc = 0;
    intdeo(f5, 0.0, 1.0, 1.0e-15, &i, &err);
    printf("I_5=int_0^infty sin(x)/x dx\n");
    printf(" I_5= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
    nfunc = 0;
    intdeo(f6, 0.0, 1.0, 1.0e-15, &i, &err);
    printf("I_6=int_0^infty cos(x)/sqrt(x) dx\n");
    printf(" I_6= %lg\t, err= %lg\t, N= %d\n", i, err, nfunc);
}


double f1(double x)
{
    extern int nfunc;
    
    nfunc++;
    return 1 / sqrt(x);
}


double f2(double x)
{
    extern int nfunc;
    
    nfunc++;
    return sqrt(4 - x * x);
}


double f3(double x)
{
    extern int nfunc;
    
    nfunc++;
    return 1 / (1 + x * x);
}


double f4(double x)
{
    extern int nfunc;
    
    nfunc++;
    return exp(-x) / sqrt(x);
}


double f5(double x)
{
    extern int nfunc;
    
    nfunc++;
    return sin(x) / x;
}


double f6(double x)
{
    extern int nfunc;
    
    nfunc++;
    return cos(x) / sqrt(x);
}

