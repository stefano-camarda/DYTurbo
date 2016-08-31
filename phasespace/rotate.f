      subroutine rotate(vin, c, s, ax, vout)
c     Rotate the 3-vector "vin" around the axis "ax"
c     by and angle phi with "c" = cos(phi) and "s" = sin(phi)
c     "vout" is the result
      implicit none
      double precision vin(3)
      double precision c,s
      double precision ax(3)
      double precision vout(3)
      vout(1)=(c+ax(1)*ax(1)*(1-c))      *vin(1)
     .     +  (ax(1)*ax(2)*(1-c)-ax(3)*s)*vin(2)
     .     +  (ax(1)*ax(3)*(1-c)+ax(2)*s)*vin(3)
      vout(2)=(ax(2)*ax(1)*(1-c)+ax(3)*s)*vin(1)
     .     +  (c+ax(2)*ax(2)*(1-c))      *vin(2)
     .     +  (ax(2)*ax(3)*(1-c)-ax(1)*s)*vin(3)
      vout(3)=(ax(3)*ax(1)*(1-c)-ax(2)*s)*vin(1)
     .     +  (ax(3)*ax(2)*(1-c)+ax(1)*s)*vin(2)
     .     +  (c+ax(3)*ax(3)*(1-c))      *vin(3)
      return
      end
