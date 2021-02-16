# ifndef smolpack_h
# define smolpack_h

# define maxdim 40

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

extern double int_smolyak(int, int, double(*ff)(int, double x[]), int);
extern double cc_int_smolyak(int, int, double(*ff)(int, double x[]), int);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */
  
#endif
