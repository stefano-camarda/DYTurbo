      subroutine genp(costh,phi,m,p)
c     generate 4-momentum from costh, phi, m
c     Important: phi must be in [-pi,pi]
      implicit none
      double precision costh, phi, m
      double precision p(4)
      double precision sintheta,cosphi,sinphi
      sintheta = sqrt(max(0d0,1d0-costh**2))
      cosphi = cos(phi);
      sinphi = sqrt(max(0d0,1d0-cosphi**2))
      if (phi.lt.0) sinphi = -sinphi
      p(4)=m/2d0                !E
      p(1)=p(4)*sintheta*sinphi !px
      p(2)=p(4)*sintheta*cosphi !py
      p(3)=p(4)*costh           !pz
      return
      end      
      
