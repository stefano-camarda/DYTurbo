      function pt(p,k)
      implicit none
      double precision pt
      double precision p(4)
      double precision k(4)

      pt = sqrt( (p(1)+k(1))**2+(p(2)+k(2))**2 )
      return
      end

      function mjj(p,k)
      implicit none
      double precision mjj
      double precision p(4)
      double precision k(4)

      mjj = sqrt( (p(4)+k(4))**2
     .     -(p(1)+k(1))**2
     .     -(p(2)+k(2))**2
     .     -(p(3)+k(3))**2 )
      return
      end
