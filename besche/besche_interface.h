#ifndef besche_interface_h
#define besche_interface_h
extern "C"
{
  void besche_(double(*F)(double&),double &C, double ALFA[], int &NUM, int &NU, int &N, double RESULT[], int INFO[]);

}
#endif
