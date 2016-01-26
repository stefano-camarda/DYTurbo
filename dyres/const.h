         integer nf
	   real* 8 Euler,Z2,Z3
	   real *8 beta0,beta1,beta2,A1g,A2g,A3g,B1g,B2g,
     #	         A1q,A2q,A3q,B1q,B2q,CA,Cf,TR,b0,C1ggn,C1qqn
	   double precision pi
c	   parameter(pi=3.14159265358979d0)
	   double precision pisq329
c	   parameter(pisq329=2d0*pi**2d0/3d0-16d0/3d0)
c	   parameter(pisq329=2d0*3.14159265358979d0**2/3d0-16/3d0)
	   complex *16 II
	   common/nf/nf
	   common/Euler/Euler
	   common/Zeta/Z2,Z3
	   common/const/beta0,beta1,beta2,A1g,A2g,A3g,B1g,B2g,
     #                A1q,A2q,A3q,B1q,B2q,CA,Cf,TR,b0,C1ggn,C1qqn
	   common/data/pi,II,pisq329
	   double precision eequ,eeqd
	   parameter(eequ=2d0/3d0, eeqd=-1d0/3d0)
	   double precision gf,gZ,gW,Mz,Mw,zw,ww,sw2,cw2,aem,aweak,
     #                gLZu,gLZd,gRZu,gRZd,fLZ,fRZ,gLW,fLW
	   common/vectboson/gf,gZ,gW,Mz,Mw,zw,ww,sw2,cw2,aem,aweak,
     #                gLZu,gLZd,gRZu,gRZd,fLZ,fRZ,gLW,fLW
