#include "interface.h"
//#include "LHAPDF/LHAPDF.h"
#include "settings.h"
#include "init.h"

//rewrite initialisation functions
void dyturboinit()
{
  // Initialize efficiency variables (get rid of these)
  efficiency_.njetzero_ = 0;
  efficiency_.ncutzero_ = 0;
  efficiency_.ntotzero_ = 0;
  efficiency_.ntotshot_ = 0;
  
  // Set-up incoming beams and PS integration cut-offs
  rtsmin_.rtsmin_ = min (rtsmin_.rtsmin_, sqrt(limits_.wsqmin_ + cutoff_.cutoff_));

  if (zerowidth_.zerowidth_)
    {
      rtsmin_.rtsmin_ = dymasses_.wmass_;
      if (nproc_.nproc_ == 3) rtsmin_.rtsmin_ = dymasses_.zmass_;
    }

  taumin_.taumin_=pow((rtsmin_.rtsmin_/energy_.sroot_),2);
  taumin_.logtaumin_=log(taumin_.taumin_);
  xmin_.xmin_=taumin_.taumin_;

  pext_.p1ext_[3]=-0.5*energy_.sroot_;
  pext_.p1ext_[0]=0.;
  pext_.p1ext_[1]=0.;
  pext_.p1ext_[2]=-0.5*energy_.sroot_;

  pext_.p2ext_[3]=-0.5*energy_.sroot_;
  pext_.p2ext_[0]=0.;
  pext_.p2ext_[1]=0.;
  pext_.p2ext_[2]=+0.5*energy_.sroot_;
}
