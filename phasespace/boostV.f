      subroutine boostV(m,p,gamma,beta)
      implicit none
      double precision m,p(4),gamma,beta(3)
c     Calculate the boost 4-vector from the rest frame
c     of a particle with mass m to the
c     frame in which the particle has 4-momentum p
      gamma=p(4)/m
      beta(1)=-p(1)/p(4);
      beta(2)=-p(2)/p(4);
      beta(3)=-p(3)/p(4);
      return
      end
