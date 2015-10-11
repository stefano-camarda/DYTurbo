      subroutine branch(brwen,brzee,brtau,brtop)
      implicit none
C     Returns the lowest order branching ratios for 
C     1) W   --> e nu
C     2) Z   --> e e
C     3) tau --> e nu nubar
C     4) t   --> b W
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      
      double precision facz,facw,factau,factop,pwidth_e
c      double precision pwidth_u,pwidth_d,pwidth_n,width
      double precision brwen,brzee,brtau,brtop,xsq

      facz=esq/4d0*zmass/(6d0*pi)
      facw=gwsq/8d0*wmass/(6d0*pi)
      factau=gwsq**2/32d0/wmass**4*mtau**5/192d0/pi**3
      xsq=(wmass/mt)**2
      factop=(gw/wmass)**2*mt**3/(64d0*pi)*(1d0-xsq)**2*(1d0+2d0*xsq)

      pwidth_e=facz*(le**2+re**2)
c      pwidth_d=3*facz*(L(1)**2+R(1)**2)
c      pwidth_u=3*facz*(L(2)**2+R(2)**2)
c      pwidth_n=facz*(ln**2)
c calculated zwidth=3*pwidth_d+2*pwidth_u+3*pwidth_e+3*pwidth_n
      brzee=pwidth_e/zwidth
      brwen=facw/wwidth
      brtau=factau/tauwidth
c      brtop=factop/twidth

c      write(*,*)brzee,brwen,brtau,brtop
      return
      end
