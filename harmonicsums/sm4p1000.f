        t1 = 5.d-1*gammaEplnn
        t2 = -1.d0*gammaEplnn
        t3 = 2.5d0*gammaEplnn
        t4 = -1.4d1*gammaEplnn
        t5 = 1.275d2*gammaEplnn
        t6 = -1.705d3*gammaEplnn
        t7 = 3.14405d4*gammaEplnn
        t8 = -1.52908d6*gammaEplnn
        t9 = 1.299718d7*gammaEplnn
        t10 = -6.9566601d7*gammaEplnn
        t11 = 2.904630795d8*gammaEplnn
        t12 = t10
        t13 = invnp1*(-6.817133989815981388092041d8 + t11)
        t14 = t12 + t13
        t15 = t9
        t16 = invnp1*(1.468218141769341826438904d8 + t14)
        t17 = t15 + t16
        t18 = t8
        t19 = invnp1*(-2.412561496454934030771255d7 + t17)
        t20 = t18 + t19
        t21 = invnp1*(2.333859023498723283410072d6 + t20)
        t22 = 5.8950948046398043516092d4 + t21
        t23 = t7 + invnp1*t22
        t24 = -4.326428061868686927482486d4 + t23
        t25 = invnp1*t24
        t26 = -2.770628787878787989029661d3 + t25
        t27 = t6 + invnp1*t26
        t28 = 2.09324424603174611547729d3 + t27
        t29 = invnp1*t28
        t30 = 1.753145833333333314385527d2 + t29
        t31 = t5 + invnp1*t30
        t32 = -1.343085317460317469340225d2 + t31
        t33 = invnp1*t32
        t34 = -1.575198412698412653298874d1 + t33
        t35 = t4 + invnp1*t34
        t36 = 1.178333333333333321490954d1 + t35
        t37 = invnp1*t36
        t38 = 2.191666666666666873908298d0 + t37
        t39 = t3 + invnp1*t38
        t40 = -1.416666666666666740681535d0 + t39
        t41 = -6.666666666666666296592325d-1 + invnp1*t40
        t42 = t2 + invnp1*t41
        t43 = t1
        t44 = invnp1*(5.d-1 + t42)
        t45 = n**(-4)
        t46 = t43 + t44
        t47 = t45
        t48 = eta*t46*t47
        HSexpand = -9.231833733969403432695344d-1 + t48