
      if(k1.ne.0.and.k2.eq.0)then
         if(k1.eq.-5)then
            include 'sm50000.f'
         elseif(k1.eq.-4)then
            include 'sm40000.f'
         elseif(k1.eq.-3)then
            include 'sm30000.f'
         elseif(k1.eq.-2)then
            include 'sm20000.f'
         elseif(k1.eq.-1)then
            include 'sm10000.f'
         elseif(k1.eq. 1)then
            include 'sp10000.f'
         elseif(k1.eq. 2)then
            include 'sp20000.f'
         elseif(k1.eq. 3)then
            include 'sp30000.f'
         elseif(k1.eq. 4)then
            include 'sp40000.f'
         elseif(k1.eq. 5)then
            include 'sp50000.f'
         endif
      elseif(k1.ne.0.and.k2.ne.0.and.k3.eq.0)then
         if(k1.eq.-4.and.k2.eq.-1)then
            include 'sm4m1000.f'
         elseif(k1.eq.-4.and.k2.eq. 1)then
            include 'sm4p1000.f'
         elseif(k1.eq.-3.and.k2.eq.-2)then
            include 'sm3m2000.f'
         elseif(k1.eq.-3.and.k2.eq.-1)then
            include 'sm3m1000.f'
         elseif(k1.eq.-3.and.k2.eq. 1)then
            include 'sm3p1000.f'
         elseif(k1.eq.-3.and.k2.eq. 2)then
            include 'sm3p2000.f'
         elseif(k1.eq.-2.and.k2.eq.-3)then
            include 'sm2m3000.f'
         elseif(k1.eq.-2.and.k2.eq.-2)then
            include 'sm2m2000.f'
         elseif(k1.eq.-2.and.k2.eq.-1)then
            include 'sm2m1000.f'
         elseif(k1.eq.-2.and.k2.eq. 1)then
            include 'sm2p1000.f'
         elseif(k1.eq.-2.and.k2.eq. 2)then
            include 'sm2p2000.f'
         elseif(k1.eq.-2.and.k2.eq. 3)then
            include 'sm2p3000.f'
         elseif(k1.eq.-1.and.k2.eq.-4)then
            include 'sm1m4000.f'
         elseif(k1.eq.-1.and.k2.eq.-3)then
            include 'sm1m3000.f'
         elseif(k1.eq.-1.and.k2.eq.-2)then
            include 'sm1m2000.f'
         elseif(k1.eq.-1.and.k2.eq.-1)then
            include 'sm1m1000.f'
         elseif(k1.eq.-1.and.k2.eq. 1)then
            include 'sm1p1000.f'
         elseif(k1.eq.-1.and.k2.eq. 2)then
            include 'sm1p2000.f'
         elseif(k1.eq.-1.and.k2.eq. 3)then
            include 'sm1p3000.f'
         elseif(k1.eq.-1.and.k2.eq. 4)then
            include 'sm1p4000.f'
         elseif(k1.eq. 1.and.k2.eq.-4)then
            include 'sp1m4000.f'
         elseif(k1.eq. 1.and.k2.eq.-3)then
            include 'sp1m3000.f'
         elseif(k1.eq. 1.and.k2.eq.-2)then
            include 'sp1m2000.f'
         elseif(k1.eq. 1.and.k2.eq.-1)then
            include 'sp1m1000.f'
         elseif(k1.eq. 1.and.k2.eq. 1)then
            include 'sp1p1000.f'
         elseif(k1.eq. 1.and.k2.eq. 2)then
            include 'sp1p2000.f'
         elseif(k1.eq. 1.and.k2.eq. 3)then
            include 'sp1p3000.f'
         elseif(k1.eq. 1.and.k2.eq. 4)then
            include 'sp1p4000.f'
         elseif(k1.eq. 2.and.k2.eq.-3)then
            include 'sp2m3000.f'
         elseif(k1.eq. 2.and.k2.eq.-2)then
            include 'sp2m2000.f'
         elseif(k1.eq. 2.and.k2.eq.-1)then
            include 'sp2m1000.f'
         elseif(k1.eq. 2.and.k2.eq. 1)then
            include 'sp2p1000.f'
         elseif(k1.eq. 2.and.k2.eq. 2)then
            include 'sp2p2000.f'
         elseif(k1.eq. 2.and.k2.eq. 3)then
            include 'sp2p3000.f'
         elseif(k1.eq. 3.and.k2.eq.-2)then
            include 'sp3m2000.f'
         elseif(k1.eq. 3.and.k2.eq.-1)then
            include 'sp3m1000.f'
         elseif(k1.eq. 3.and.k2.eq. 1)then
            include 'sp3p1000.f'
         elseif(k1.eq. 3.and.k2.eq. 2)then
            include 'sp3p2000.f'
         elseif(k1.eq. 4.and.k2.eq.-1)then
            include 'sp4m1000.f'
         elseif(k1.eq. 4.and.k2.eq. 1)then
            include 'sp4p1000.f'
         endif
      elseif(k1.ne.0.and.k2.ne.0.and.k3.ne.0.and.k4.eq.0)then
         if(k1.eq.-3.and.k2.eq.-1.and.k3.eq.-1)then
            include 'sm3m1m100.f'
         elseif(k1.eq.-3.and.k2.eq.-1.and.k3.eq. 1)then
            include 'sm3m1p100.f'
         elseif(k1.eq.-3.and.k2.eq. 1.and.k3.eq.-1)then
            include 'sm3p1m100.f'
         elseif(k1.eq.-3.and.k2.eq. 1.and.k3.eq. 1)then
            include 'sm3p1p100.f'
         elseif(k1.eq.-2.and.k2.eq.-2.and.k3.eq.-1)then
            include 'sm2m2m100.f'
         elseif(k1.eq.-2.and.k2.eq.-2.and.k3.eq. 1)then
            include 'sm2m2p100.f'
         elseif(k1.eq.-2.and.k2.eq.-1.and.k3.eq.-2)then
            include 'sm2m1m200.f'
         elseif(k1.eq.-2.and.k2.eq.-1.and.k3.eq.-1)then
            include 'sm2m1m100.f'
         elseif(k1.eq.-2.and.k2.eq.-1.and.k3.eq. 1)then
            include 'sm2m1p100.f'
         elseif(k1.eq.-2.and.k2.eq.-1.and.k3.eq. 2)then
            include 'sm2m1p200.f'
         elseif(k1.eq.-2.and.k2.eq. 1.and.k3.eq.-2)then
            include 'sm2p1m200.f'
         elseif(k1.eq.-2.and.k2.eq. 1.and.k3.eq.-1)then
            include 'sm2p1m100.f'
         elseif(k1.eq.-2.and.k2.eq. 1.and.k3.eq. 1)then
            include 'sm2p1p100.f'
         elseif(k1.eq.-2.and.k2.eq. 1.and.k3.eq. 2)then
            include 'sm2p1p200.f'
         elseif(k1.eq.-2.and.k2.eq. 2.and.k3.eq.-1)then
            include 'sm2p2m100.f'
         elseif(k1.eq.-2.and.k2.eq. 2.and.k3.eq. 1)then
            include 'sm2p2p100.f'
         elseif(k1.eq.-1.and.k2.eq.-3.and.k3.eq.-1)then
            include 'sm1m3m100.f'
         elseif(k1.eq.-1.and.k2.eq.-3.and.k3.eq. 1)then
            include 'sm1m3p100.f'
         elseif(k1.eq.-1.and.k2.eq.-2.and.k3.eq.-2)then
            include 'sm1m2m200.f'
         elseif(k1.eq.-1.and.k2.eq.-2.and.k3.eq.-1)then
            include 'sm1m2m100.f'
         elseif(k1.eq.-1.and.k2.eq.-2.and.k3.eq. 1)then
            include 'sm1m2p100.f'
         elseif(k1.eq.-1.and.k2.eq.-2.and.k3.eq. 2)then
            include 'sm1m2p200.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-3)then
            include 'sm1m1m300.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-2)then
            include 'sm1m1m200.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-1)then
            include 'sm1m1m100.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 1)then
            include 'sm1m1p100.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 2)then
            include 'sm1m1p200.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 3)then
            include 'sm1m1p300.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-3)then
            include 'sm1p1m300.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-2)then
            include 'sm1p1m200.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-1)then
            include 'sm1p1m100.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 1)then
            include 'sm1p1p100.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 2)then
            include 'sm1p1p200.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 3)then
            include 'sm1p1p300.f'
         elseif(k1.eq.-1.and.k2.eq. 2.and.k3.eq.-2)then
            include 'sm1p2m200.f'
         elseif(k1.eq.-1.and.k2.eq. 2.and.k3.eq.-1)then
            include 'sm1p2m100.f'
         elseif(k1.eq.-1.and.k2.eq. 2.and.k3.eq. 1)then
            include 'sm1p2p100.f'
         elseif(k1.eq.-1.and.k2.eq. 2.and.k3.eq. 2)then
            include 'sm1p2p200.f'
         elseif(k1.eq.-1.and.k2.eq. 3.and.k3.eq.-1)then
            include 'sm1p3m100.f'
         elseif(k1.eq.-1.and.k2.eq. 3.and.k3.eq. 1)then
            include 'sm1p3p100.f'
         elseif(k1.eq. 1.and.k2.eq.-3.and.k3.eq.-1)then
            include 'sp1m3m100.f'
         elseif(k1.eq. 1.and.k2.eq.-3.and.k3.eq. 1)then
            include 'sp1m3p100.f'
         elseif(k1.eq. 1.and.k2.eq.-2.and.k3.eq.-2)then
            include 'sp1m2m200.f'
         elseif(k1.eq. 1.and.k2.eq.-2.and.k3.eq.-1)then
            include 'sp1m2m100.f'
         elseif(k1.eq. 1.and.k2.eq.-2.and.k3.eq. 1)then
            include 'sp1m2p100.f'
         elseif(k1.eq. 1.and.k2.eq.-2.and.k3.eq. 2)then
            include 'sp1m2p200.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-3)then
            include 'sp1m1m300.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-2)then
            include 'sp1m1m200.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-1)then
            include 'sp1m1m100.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 1)then
            include 'sp1m1p100.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 2)then
            include 'sp1m1p200.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 3)then
            include 'sp1m1p300.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-3)then
            include 'sp1p1m300.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-2)then
            include 'sp1p1m200.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-1)then
            include 'sp1p1m100.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 1)then
            include 'sp1p1p100.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 2)then
            include 'sp1p1p200.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 3)then
            include 'sp1p1p300.f'
         elseif(k1.eq. 1.and.k2.eq. 2.and.k3.eq.-2)then
            include 'sp1p2m200.f'
         elseif(k1.eq. 1.and.k2.eq. 2.and.k3.eq.-1)then
            include 'sp1p2m100.f'
         elseif(k1.eq. 1.and.k2.eq. 2.and.k3.eq. 1)then
            include 'sp1p2p100.f'
         elseif(k1.eq. 1.and.k2.eq. 2.and.k3.eq. 2)then
            include 'sp1p2p200.f'
         elseif(k1.eq. 1.and.k2.eq. 3.and.k3.eq.-1)then
            include 'sp1p3m100.f'
         elseif(k1.eq. 1.and.k2.eq. 3.and.k3.eq. 1)then
            include 'sp1p3p100.f'
         elseif(k1.eq. 2.and.k2.eq.-2.and.k3.eq.-1)then
            include 'sp2m2m100.f'
         elseif(k1.eq. 2.and.k2.eq.-2.and.k3.eq. 1)then
            include 'sp2m2p100.f'
         elseif(k1.eq. 2.and.k2.eq.-1.and.k3.eq.-2)then
            include 'sp2m1m200.f'
         elseif(k1.eq. 2.and.k2.eq.-1.and.k3.eq.-1)then
            include 'sp2m1m100.f'
         elseif(k1.eq. 2.and.k2.eq.-1.and.k3.eq. 1)then
            include 'sp2m1p100.f'
         elseif(k1.eq. 2.and.k2.eq.-1.and.k3.eq. 2)then
            include 'sp2m1p200.f'
         elseif(k1.eq. 2.and.k2.eq. 1.and.k3.eq.-2)then
            include 'sp2p1m200.f'
         elseif(k1.eq. 2.and.k2.eq. 1.and.k3.eq.-1)then
            include 'sp2p1m100.f'
         elseif(k1.eq. 2.and.k2.eq. 1.and.k3.eq. 1)then
            include 'sp2p1p100.f'
         elseif(k1.eq. 2.and.k2.eq. 1.and.k3.eq. 2)then
            include 'sp2p1p200.f'
         elseif(k1.eq. 2.and.k2.eq. 2.and.k3.eq.-1)then
            include 'sp2p2m100.f'
         elseif(k1.eq. 2.and.k2.eq. 2.and.k3.eq. 1)then
            include 'sp2p2p100.f'
         elseif(k1.eq. 3.and.k2.eq.-1.and.k3.eq.-1)then
            include 'sp3m1m100.f'
         elseif(k1.eq. 3.and.k2.eq.-1.and.k3.eq. 1)then
            include 'sp3m1p100.f'
         elseif(k1.eq. 3.and.k2.eq. 1.and.k3.eq.-1)then
            include 'sp3p1m100.f'
         elseif(k1.eq. 3.and.k2.eq. 1.and.k3.eq. 1)then
            include 'sp3p1p100.f'
         endif
      elseif(k1.ne.0.and.k2.ne.0.and.k3.ne.0
     .        .and.k4.ne.0.and.k5.eq.0) then
         if(k1.eq.-2.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sm2m1m1m10.f'
         elseif(k1.eq.-2.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sm2m1m1p10.f'
         elseif(k1.eq.-2.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sm2m1p1m10.f'
         elseif(k1.eq.-2.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sm2m1p1p10.f'
         elseif(k1.eq.-2.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sm2p1m1m10.f'
         elseif(k1.eq.-2.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sm2p1m1p10.f'
         elseif(k1.eq.-2.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sm2p1p1m10.f'
         elseif(k1.eq.-2.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sm2p1p1p10.f'
         elseif(k1.eq.-1.and.k2.eq.-2.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sm1m2m1m10.f'
         elseif(k1.eq.-1.and.k2.eq.-2.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sm1m2m1p10.f'
         elseif(k1.eq.-1.and.k2.eq.-2.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sm1m2p1m10.f'
         elseif(k1.eq.-1.and.k2.eq.-2.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sm1m2p1p10.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-2.and.k4.eq.-1)then
            include 'sm1m1m2m10.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-2.and.k4.eq. 1)then
            include 'sm1m1m2p10.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq.-2)then
            include 'sm1m1m1m20.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sm1m1m1m10.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sm1m1m1p10.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq. 2)then
            include 'sm1m1m1p20.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq.-2)then
            include 'sm1m1p1m20.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sm1m1p1m10.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sm1m1p1p10.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq. 2)then
            include 'sm1m1p1p20.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 2.and.k4.eq.-1)then
            include 'sm1m1p2m10.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 2.and.k4.eq. 1)then
            include 'sm1m1p2p10.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-2.and.k4.eq.-1)then
            include 'sm1p1m2m10.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-2.and.k4.eq. 1)then
            include 'sm1p1m2p10.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq.-2)then
            include 'sm1p1m1m20.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sm1p1m1m10.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sm1p1m1p10.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq. 2)then
            include 'sm1p1m1p20.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq.-2)then
            include 'sm1p1p1m20.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sm1p1p1m10.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sm1p1p1p10.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq. 2)then
            include 'sm1p1p1p20.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 2.and.k4.eq.-1)then
            include 'sm1p1p2m10.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 2.and.k4.eq. 1)then
            include 'sm1p1p2p10.f'
         elseif(k1.eq.-1.and.k2.eq. 2.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sm1p2m1m10.f'
         elseif(k1.eq.-1.and.k2.eq. 2.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sm1p2m1p10.f'
         elseif(k1.eq.-1.and.k2.eq. 2.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sm1p2p1m10.f'
         elseif(k1.eq.-1.and.k2.eq. 2.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sm1p2p1p10.f'
         elseif(k1.eq. 1.and.k2.eq.-2.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sp1m2m1m10.f'
         elseif(k1.eq. 1.and.k2.eq.-2.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sp1m2m1p10.f'
         elseif(k1.eq. 1.and.k2.eq.-2.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sp1m2p1m10.f'
         elseif(k1.eq. 1.and.k2.eq.-2.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sp1m2p1p10.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-2.and.k4.eq.-1)then
            include 'sp1m1m2m10.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-2.and.k4.eq. 1)then
            include 'sp1m1m2p10.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq.-2)then
            include 'sp1m1m1m20.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sp1m1m1m10.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sp1m1m1p10.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq. 2)then
            include 'sp1m1m1p20.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq.-2)then
            include 'sp1m1p1m20.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sp1m1p1m10.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sp1m1p1p10.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq. 2)then
            include 'sp1m1p1p20.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 2.and.k4.eq.-1)then
            include 'sp1m1p2m10.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 2.and.k4.eq. 1)then
            include 'sp1m1p2p10.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-2.and.k4.eq.-1)then
            include 'sp1p1m2m10.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-2.and.k4.eq. 1)then
            include 'sp1p1m2p10.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq.-2)then
            include 'sp1p1m1m20.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sp1p1m1m10.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sp1p1m1p10.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq. 2)then
            include 'sp1p1m1p20.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq.-2)then
            include 'sp1p1p1m20.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sp1p1p1m10.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sp1p1p1p10.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq. 2)then
            include 'sp1p1p1p20.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 2.and.k4.eq.-1)then
            include 'sp1p1p2m10.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 2.and.k4.eq. 1)then
            include 'sp1p1p2p10.f'
         elseif(k1.eq. 1.and.k2.eq. 2.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sp1p2m1m10.f'
         elseif(k1.eq. 1.and.k2.eq. 2.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sp1p2m1p10.f'
         elseif(k1.eq. 1.and.k2.eq. 2.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sp1p2p1m10.f'
         elseif(k1.eq. 1.and.k2.eq. 2.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sp1p2p1p10.f'
         elseif(k1.eq. 2.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sp2m1m1m10.f'
         elseif(k1.eq. 2.and.k2.eq.-1.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sp2m1m1p10.f'
         elseif(k1.eq. 2.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sp2m1p1m10.f'
         elseif(k1.eq. 2.and.k2.eq.-1.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sp2m1p1p10.f'
         elseif(k1.eq. 2.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq.-1)then
            include 'sp2p1m1m10.f'
         elseif(k1.eq. 2.and.k2.eq. 1.and.k3.eq.-1.and.k4.eq. 1)then
            include 'sp2p1m1p10.f'
         elseif(k1.eq. 2.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq.-1)then
            include 'sp2p1p1m10.f'
         elseif(k1.eq. 2.and.k2.eq. 1.and.k3.eq. 1.and.k4.eq. 1)then
            include 'sp2p1p1p10.f'
         endif
      else
         if(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-1
     .        .and.k4.eq.-1.and.k5.eq.-1)then
            include 'sm1m1m1m1m1.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-1
     .           .and.k4.eq.-1.and.k5.eq. 1)then
            include 'sm1m1m1m1p1.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-1
     .       .and.k4.eq. 1.and.k5.eq.-1)then
            include 'sm1m1m1p1m1.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq.-1
     .       .and.k4.eq. 1.and.k5.eq. 1)then
            
            include 'sm1m1m1p1p1.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 1
     .       .and.k4.eq.-1.and.k5.eq.-1)then
            
            include 'sm1m1p1m1m1.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 1
     .       .and.k4.eq.-1.and.k5.eq. 1)then
            
            include 'sm1m1p1m1p1.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 1
     .       .and.k4.eq. 1.and.k5.eq.-1)then
            
            include 'sm1m1p1p1m1.f'
         elseif(k1.eq.-1.and.k2.eq.-1.and.k3.eq. 1
     .       .and.k4.eq. 1.and.k5.eq. 1)then
            
            include 'sm1m1p1p1p1.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-1
     .       .and.k4.eq.-1.and.k5.eq.-1)then
            
            include 'sm1p1m1m1m1.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-1
     .       .and.k4.eq.-1.and.k5.eq. 1)then
            
            include 'sm1p1m1m1p1.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-1
     .       .and.k4.eq. 1.and.k5.eq.-1)then
            
            include 'sm1p1m1p1m1.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq.-1
     .       .and.k4.eq. 1.and.k5.eq. 1)then
            
            include 'sm1p1m1p1p1.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 1
     .       .and.k4.eq.-1.and.k5.eq.-1)then
            
            include 'sm1p1p1m1m1.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 1
     .       .and.k4.eq.-1.and.k5.eq. 1)then
            
            include 'sm1p1p1m1p1.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 1
     .       .and.k4.eq. 1.and.k5.eq.-1)then
            
            include 'sm1p1p1p1m1.f'
         elseif(k1.eq.-1.and.k2.eq. 1.and.k3.eq. 1
     .       .and.k4.eq. 1.and.k5.eq. 1)then
            
            include 'sm1p1p1p1p1.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-1
     .       .and.k4.eq.-1.and.k5.eq.-1)then
            
            include 'sp1m1m1m1m1.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-1
     .       .and.k4.eq.-1.and.k5.eq. 1)then
            
            include 'sp1m1m1m1p1.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-1
     .       .and.k4.eq. 1.and.k5.eq.-1)then
            
            include 'sp1m1m1p1m1.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq.-1
     .       .and.k4.eq. 1.and.k5.eq. 1)then
            
            include 'sp1m1m1p1p1.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 1
     .       .and.k4.eq.-1.and.k5.eq.-1)then
            
            include 'sp1m1p1m1m1.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 1
     .       .and.k4.eq.-1.and.k5.eq. 1)then
            
            include 'sp1m1p1m1p1.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 1
     .       .and.k4.eq. 1.and.k5.eq.-1)then
            
            include 'sp1m1p1p1m1.f'
         elseif(k1.eq. 1.and.k2.eq.-1.and.k3.eq. 1
     .       .and.k4.eq. 1.and.k5.eq. 1)then
            
            include 'sp1m1p1p1p1.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-1
     .       .and.k4.eq.-1.and.k5.eq.-1)then
            
            include 'sp1p1m1m1m1.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-1
     .       .and.k4.eq.-1.and.k5.eq. 1)then
            
            include 'sp1p1m1m1p1.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-1
     .       .and.k4.eq. 1.and.k5.eq.-1)then
            
            include 'sp1p1m1p1m1.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq.-1
     .       .and.k4.eq. 1.and.k5.eq. 1)then
            
            include 'sp1p1m1p1p1.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 1
     .       .and.k4.eq.-1.and.k5.eq.-1)then
            
            include 'sp1p1p1m1m1.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 1
     .       .and.k4.eq.-1.and.k5.eq. 1)then
            
            include 'sp1p1p1m1p1.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 1
     .       .and.k4.eq. 1.and.k5.eq.-1)then
            
            include 'sp1p1p1p1m1.f'
         elseif(k1.eq. 1.and.k2.eq. 1.and.k3.eq. 1
     .       .and.k4.eq. 1.and.k5.eq. 1)then
            
            include 'sp1p1p1p1p1.f'
         endif
      endif
