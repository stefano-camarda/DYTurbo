      function myLI3(x)
      implicit none
      double precision myLI3,xlog,x,PI,Z3
      PI=3.14159265358979312D0
      Z3=1.20205690315959429D0
      
      
      if (x.lt.-20) then ! asymptotic expansion, suggested by Mathematica
      xlog=dlog(-x)
      myLI3= -1d0/6*(xlog**3+Pi**2*xlog) +1/x +1d0/8/x**2 +1d0/27/x**3
      elseif (x.lt.0.5d0) then
      xlog=dlog(1d0-x)
      myLI3= -xlog -(3*xlog**2)/8. -(17*xlog**3)/216. -(5*xlog**4)/576
     .   -(7*xlog**5)/54000. +(7*xlog**6)/86400. +19*xlog**7/5556600
     .   -xlog**8/752640 -11*xlog**9/127008000 +11*xlog**10/435456000
      elseif (x.lt.1d0) then
      xlog=dlog(x)
      myLI3=Z3 +(Pi**2*xlog)/6 +(3d0/4-dlog(-xlog)/2)*xlog**2 
     .   -xlog**3/12-xlog**4/288 +xlog**6/86400 -xlog**8/10160640  
      elseif (x.eq.1d0) then
      myLI3=Z3
      else
      write(6,*)'wrong argument of Li3!!' 
      endif
      return
      end   
