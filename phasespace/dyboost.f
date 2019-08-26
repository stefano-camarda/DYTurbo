      subroutine dyboost(gamma,beta,pin,pout)
c     boost 4-momentum pin into 4-momentum pout
      implicit none
      double precision pin(4),pout(4)
      double precision gamma,beta(3),bdotp
      double precision one
      parameter(one=1d0)
      bdotp=pin(1)*beta(1)+pin(2)*beta(2)+pin(3)*beta(3);
      pout(4)=gamma*(pin(4)-bdotp)
      pout(1)=pin(1)+gamma*beta(1)*(gamma/(gamma+one)*bdotp-pin(4))
      pout(2)=pin(2)+gamma*beta(2)*(gamma/(gamma+one)*bdotp-pin(4))
      pout(3)=pin(3)+gamma*beta(3)*(gamma/(gamma+one)*bdotp-pin(4))
      return
      end      
      
