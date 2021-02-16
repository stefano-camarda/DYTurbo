      function fort_besj0(x)
      double precision x,fort_besj0
      fort_besj0 = dbesj0(x)
      return
      end

      function fort_besy0(x)
      double precision x,fort_besy0
      fort_besy0 = dbesy0(x)
      return
      end

      function fort_besj1(x)
      double precision x,fort_besj1
      fort_besj1 = dbesj1(x)
      return
      end
      
      function fort_besy1(x)
      double precision x,fort_besy1
      fort_besy1 = dbesy1(x)
      return
      end

      function fort_besjn(n,x)
      double precision x,fort_besjn
      integer n
      fort_besjn = dbesjn(n,x)
      return
      end

      function fort_besyn(n,x)
      double precision x,fort_besyn
      integer n
      fort_besyn = dbesyn(n,x)
      return
      end
      
