#include "pdfevol.h"
#include "interface.h"
#include <iostream>

complex <double> *pdfevol::fn1;
complex <double> *pdfevol::fn2;

//fortran interface
void pdfevol_(int& i1, int& i2, int& sign)
{
  pdfevol::evolve(i1-1, i2-1, sign-1);
};

void pdfevol::init()
{
  fn1 = new complex <double>[11];
  fn2 = new complex <double>[11];
}

void pdfevol::evolve(int i1, int i2, int sign)
{
  //  cout << i1 << endl;
  //  cout << creno_.cfx1_[i1][5].real << "  " << creno_.cfx1_[i1][5].imag << endl;
  fn1[0] = cx(creno_.cfx1_[i1][0]);
  fn1[1] = cx(creno_.cfx1_[i1][1]);
  fn1[2] = cx(creno_.cfx1_[i1][2]);
  fn1[3] = cx(creno_.cfx1_[i1][3]);
  fn1[4] = cx(creno_.cfx1_[i1][4]);
  fn1[5] = cx(creno_.cfx1_[i1][5]);
  fn1[6] = cx(creno_.cfx1_[i1][6]);
  fn1[7] = cx(creno_.cfx1_[i1][7]);
  fn1[8] = cx(creno_.cfx1_[i1][8]);
  fn1[9] = cx(creno_.cfx1_[i1][9]);
  fn1[10] = cx(creno_.cfx1_[i1][10]);
  if (sign == mesq::positive)
    {
      fn2[0] = cx(creno_.cfx2p_[i2][0]);
      fn2[1] = cx(creno_.cfx2p_[i2][1]);
      fn2[2] = cx(creno_.cfx2p_[i2][2]);
      fn2[3] = cx(creno_.cfx2p_[i2][3]);
      fn2[4] = cx(creno_.cfx2p_[i2][4]);
      fn2[5] = cx(creno_.cfx2p_[i2][5]);
      fn2[6] = cx(creno_.cfx2p_[i2][6]);
      fn2[7] = cx(creno_.cfx2p_[i2][7]);
      fn2[8] = cx(creno_.cfx2p_[i2][8]);
      fn2[9] = cx(creno_.cfx2p_[i2][9]);
      fn2[10] = cx(creno_.cfx2p_[i2][10]);
    }
  else if (sign == mesq::negative)
    {
      fn2[0] = cx(creno_.cfx2m_[i2][0]);
      fn2[1] = cx(creno_.cfx2m_[i2][1]);
      fn2[2] = cx(creno_.cfx2m_[i2][2]);
      fn2[3] = cx(creno_.cfx2m_[i2][3]);
      fn2[4] = cx(creno_.cfx2m_[i2][4]);
      fn2[5] = cx(creno_.cfx2m_[i2][5]);
      fn2[6] = cx(creno_.cfx2m_[i2][6]);
      fn2[7] = cx(creno_.cfx2m_[i2][7]);
      fn2[8] = cx(creno_.cfx2m_[i2][8]);
      fn2[9] = cx(creno_.cfx2m_[i2][9]);
      fn2[10] = cx(creno_.cfx2m_[i2][10]);
    }
}
