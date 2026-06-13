      function FLpppp(x1,y,lx1,ly)       
      implicit none
      complex*16 FLpppp
      double precision x1,y,lx1,ly
      double precision Zeta3,pi
      double precision ldli2,ldli3,ldli4

      external ldli2,ldli3,ldli4
      Zeta3=1.202056903d0
      pi=dacos(-1.d0)

      FLpppp=-(x1*x1+y*y)*(4.0d0*ldli4(-x1)+(ly-3.0d0*lx1)*ldli3(-x1)
     . +lx1*lx1*(ldli2(-x1)-pi*pi/12.0d0)+1.0d0/48.0d0*((lx1+ly)**4.0d0)
     .        -109.0d0/720.0d0*(pi**4.0d0))
     . -2.0d0*x1*y*(ldli3(-x1/y)-(lx1-ly)*ldli2(-x1/y)-Zeta3
     .        +0.5d0*ly*((lx1-ly)*(lx1-ly)+pi*pi))
     .  -0.5d0*x1*x1*(ldli3(-x1/y)-(lx1-ly)*ldli2(-x1/y)-Zeta3
     . -0.5d0*((lx1-ly)**3.0d0)-ly*((lx1-ly)*(lx1-ly)+pi*pi))
     . +1.0d0/8.0d0*(y*y+12.0d0*x1*y-27.0d0*x1*x1-8.0d0/y+9.0d0/y/y)*
     . lx1*lx1-1.0d0/8.0d0*(38.0d0*x1*y-13.0d0)*lx1*ly +pi*pi/48.0d0*
     . (114.0d0*x1*y-43.0d0)-9.0d0/4.0d0*(1.0d0/y+2.0d0*x1)*lx1+0.25d0
     . +dcmplx(0,1)*pi*(-(x1*x1+y*y)*(ldli3(-x1/y)-(lx1-ly)*ldli2(-x1/y)
     .  -Zeta3 +0.5d0*(ly-0.75d0)*((lx1-ly)*(lx1-ly)+pi*pi))
     . +0.25d0*(14.0d0*(x1-y)-8.0d0/y+9.0d0/y/y)*lx1-9.0d0/4.0d0*
     . (1.0d0/y-1.0d0))-(y*y+x1*x1)*(4.0d0*ldli4(-y)+(lx1-3.0d0*ly)*
     . ldli3(-y)+ly*ly*(ldli2(-y)-pi*pi/12.0d0)+1.0d0/48.0d0*
     .   ((ly+lx1)**4.0d0)        -109.0d0/720.0d0*(pi**4.0d0))
     .    -2.0d0*y*x1*(ldli3(-y/x1)-(ly-lx1)*ldli2(-y/x1)-Zeta3
     .           +0.5d0*lx1*((ly-lx1)*(ly-lx1)+pi*pi))
     .    -0.5d0*y*y*(ldli3(-y/x1)-(ly-lx1)*ldli2(-y/x1)-Zeta3
     .  -0.5d0*((ly-lx1)**3.0d0)-lx1*((ly-lx1)*(ly-lx1)+pi*pi))
     . +1.0d0/8.0d0*(x1*x1+12.0d0*y*x1-27.0d0*y*y-8.0d0/x1+9.0d0/x1/x1)
     . *ly*ly-1.0d0/8.0d0*(38.0d0*y*x1-13.0d0)*ly*lx1 +pi*pi/48.0d0*
     . (114.0d0*y*x1-43.0d0)-9.0d0/4.0d0*(1.0d0/x1+2.0d0*y)*ly+0.25d0
     . +dcmplx(0,1)*pi*(-(y*y+x1*x1)*(ldli3(-y/x1)-(ly-lx1)*ldli2(-y/x1)
     .   -Zeta3+0.5d0*(lx1-0.75d0)*((ly-lx1)*(ly-lx1)+pi*pi))
     . +0.25d0*(14.0d0*(y-x1)-8.0d0/x1+9.0d0/x1/x1)*ly-9.0d0/4.0d0*
     . (1.0d0/x1-1.0d0))

      return
      end

      function FSLmmpp(x1,y,lx1,ly)
      implicit none
      double precision x1,y,lx1,ly
      complex*16 FSLmmpp
      FSLmmpp=dcmplx(-1.5d0,0.)
      return
      end

      function FLmmpp(x1,y,lx1,ly)
      implicit none
      double precision x1,y,lx1,ly
      complex*16 FLmmpp
      FLmmpp=dcmplx(0.5d0,0.)
      return
      end

      function FSLpppp(x1,y,lx1,ly)
      implicit none
      double precision x1,y,lx1,ly
      complex*16 FSLpppp
      double precision Zeta3,pi
      double precision ldli2,ldli3,ldli4

      external ldli2,ldli3,ldli4
      Zeta3=1.202056903d0
      pi=dacos(-1.d0)

      FSLpppp=-2.0d0*(x1*x1+y*y)*ldli4(-x1)-(x1-y)*ldli4(-x1/y)
     . +2.0d0*x1*
     . x1*lx1*(ldli3(-x1)+ldli3(-y))
     . +pi*pi/6.0d0*(x1-y)*ldli2(-x1)-x1*x1*(1.0d0/6.0d0*(lx1**4.0d0)
     . -2.0d0/3.0d0*lx1*lx1*lx1*ly
     . +pi*pi/6.0d0*lx1*ly-4.0d0/45.0d0*(pi**4.0d0))
     . -x1*(2.0d0*ldli3(-x1)-2.0d0*lx1*ldli2(-x1)-ldli3(-x1/y)+(lx1-ly)*
     . ldli2(-x1/y)-3.0d0*Zeta3
     . -2.0d0/3.0d0*lx1*(lx1*lx1+pi*pi)+5.0d0/12.0d0*(lx1-ly)*
     . ((lx1-ly)*(lx1-ly)+pi*pi)+lx1*lx1*(lx1-ly))
     . -0.25d0*(2.0d0-3.0d0*x1*x1+(x1**4.0d0)/y/y)*lx1*lx1+0.25d0*
     . (2.0d0*x1*y+3.0d0)*lx1*ly
     . -1.0d0/24.0d0*(6.0d0*x1*y+7.0d0)*pi*pi+(0.5d0/y+x1)*lx1-0.25d0
     . +dcmplx(0,1)*pi*(-2.0d0*x1*x1*(ldli3(-y/x1)-Zeta3)+2.0d0*x1*
     . (ldli2(-y/x1)-0.25d0*(lx1-ly)*(lx1-ly)
     . -pi*pi/12.0d0)-0.5d0*(2.0d0*x1*x1-1.0d0)/y/y*lx1+0.5d0/y-0.5d0)
     . -2.0d0*(y*y+x1*x1)*ldli4(-y)-(y-x1)*ldli4(-y/x1)+2.0d0*y*y*ly*
     . (ldli3(-y)+ldli3(-x1))
     . +pi*pi/6.0d0*(y-x1)*ldli2(-y)-y*y*(1.0d0/6.0d0*(ly**4.0d0)-
     . 2.0d0/3.0d0*ly*ly*ly*lx1
     . +pi*pi/6.0d0*ly*lx1-4.0d0/45.0d0*(pi**4.0d0))
     . -y*(2.0d0*ldli3(-y)-2.0d0*ly*ldli2(-y)-ldli3(-y/x1)+(ly-lx1)*
     . ldli2(-y/x1)-3.0d0*Zeta3
     . -2.0d0/3.0d0*ly*(ly*ly+pi*pi)+5.0d0/12.0d0*(ly-lx1)*((ly-lx1)*
     . (ly-lx1)+pi*pi)+ly*ly*(ly-lx1))
     . -0.25d0*(2.0d0-3.0d0*y*y+(y**4.0d0)/x1/x1)*ly*ly+0.25d0*
     . (2.0d0*y*x1+3.0d0)*ly*lx1
     . -1.0d0/24.0d0*(6.0d0*y*x1+7.0d0)*pi*pi+(0.5d0/x1+y)*ly-0.25d0
     . +dcmplx(0,1)*pi*(-2.0d0*y*y*(ldli3(-x1/y)-Zeta3)+2.0d0*y*
     . (ldli2(-x1/y)-0.25d0*(ly-lx1)*(ly-lx1)
     . -pi*pi/12.0d0)-0.5d0*(2.0d0*y*y-1.0d0)/x1/x1*ly+0.5d0/x1-0.5d0)

      return
      end

      function H1gYYint (yy,s) !yy es coseno de tita
      implicit none
      double precision yy,s,t,u,mur,poles,Nf,H1gYYint,factors
      double precision Sum_eq2
      double precision H1gYYint_0back1sig, H1gYYint_1back0sig,
     .  M1back0sig,M0back1sig
      double precision ldli2,ldli3,ldli4,pi,Zeta3,average,identical
      double precision x1,y,lx1,ly
      double precision Minterf

      double precision s_common,mur_common
      complex*16 Msig,M1sig

      external ldli2,ldli3,ldli4

      double precision et,gh,formfactor,etau,Mh,mtau
      double precision Nc,v,color
      complex*16 MsigAA,MsigBB
      double precision xt,xb,xc,xw,xtau
      double precision eb,ec

      double precision Mu,Md,Ms,Mc,Mb,Mtop,Mz,mw

      complex*16  Martinf,MartinF1,MartinF12
      complex*16 FLmmpp,FLpmpp,FLmmmp,FLpppp,FLpmmp,FLpmpm
      complex*16 FSLmmpp,FSLpmpp,FSLmmmp,FSLpppp,FSLpmmp,FSLpmpm
      complex*16 Mqback_pppp,Mqback_mmpp,Mqback_mmmm,Mqback_ppmm
      complex*16  AggH,AyyH,AyyHt
      complex*16 Mtilde1,M1back,M0back
      complex*16 M1back_pppp,M1back_mmpp
      complex*16 M0back1sig_mmpp,M0back1sig_pppp,M0back_mmpp,
     .  M0back_pppp,M1back0sig_mmpp,M1back0sig_pppp
	  double precision H1gYYintDY


      


	  mur=sqrt(s)                 !no puede depender de la escala tiene que ser 

c       write(*,*) 's = ',s
      ! puesto a 1 este ratio
      Zeta3=1.202056903d0
      pi=dacos(-1d0)
      u = -1d0/2d0*s*(1d0-yy)
      t = -1d0/2d0*s*(1d0+yy)
      x1=t/s
      y=u/s

      lx1=dlog(-x1)
      ly=dlog(-y)

      GH=4.165d-3
      et=2d0/3d0

      MH=125.d0
      mtau=1.776d0
      Nc = 3.d0
      v = 174.103d0
   
      MW=80.385d0
      MZ=91.1876d0
            
      mtop=173.2d0
      mb=4.66d0
      mc=1.275d0

      eb = -1.d0/3.d0
      ec =  2.d0/3.d0
      
      xt=4d0*mtop**2d0/s
      xb=4d0*mb**2d0/s
      xc=4d0*mc**2d0/s
      xtau=4d0*mtau**2d0/s
      xw=4d0*mw**2d0/s

      etau = 1.d0
      
      FormFactor = 3.0d0*pi**2.0d0/2.0d0 +11.0d0/2.0d0 !H1g Hard Scheme
c ecuacion 3 de 1208.1533 (sin alphas, esta en otro lugar del codigo)
        AggH = -(s)/(8d0*dsqrt(2d0)*Pi*v)*(MartinF12(xt)+MartinF12(xb)
     .   +MartinF12(xc))
c Without the b and c contributions to compare with Lean
c        AggH = -(s)/(8d0*dsqrt(2d0)*Pi*v)*(MartinF12(xt))
c ecuacion 4 de 1208.1533 (sin alpha, esta en otro lugar del codigo)
        AyyH = -(s)/(4d0*dsqrt(2d0)*Pi*v)*(MartinF1(xw)+Nc*et**2d0*
     .    MartinF12(xt)+Nc*eb**2d0*MartinF12(xb)+Nc*ec**2d0*
     .    MartinF12(xc)+etau**2d0*MartinF12(xtau))
c Without the b and c contributions to compare with Lean
c        AyyH = -(s)/(4d0*dsqrt(2d0)*Pi*v)*(MartinF1(xw)+Nc*et**2d0*
c     .    MartinF12(xt)+etau**2d0*MartinF12(xtau))
c I also write the top contribution separately because it is the only one we
c include O(alpha_s) corrections for
	    AyyHt = -(s)/(4d0*dsqrt(2d0)*Pi*v)*(Nc*et**2d0*MartinF12(xt))
c Primer término de la ecuacion 2 de 1208.1533 (sin las helicidades ni color)
c           write(*,*)'(AggHH) ',(AggH)
c           write(*,*)'(AyyH) ',(AyyH)
c           write(*,*)' MartinF12(xt) ',MartinF12(xt)
c           write(*,*)'(s-Mh**2d0+dcmpl ',(s-Mh**2d0+dcmplx(0.,1.)*Mh*Gh)
            Msig = -(AggH*AyyH)/(s-Mh**2d0+dcmplx(0.,1.)*Mh*Gh)
          MsigAA = -(AggH*AyyH)*(s-Mh**2) /
     .       ((s-MH**2)**2 + (Mh*Gh)**2)
          MsigBB = -(AggH*AyyH)*dcmplx(0,-1)*Mh*Gh /
     .      ((s-MH**2)**2 + (Mh*Gh)**2)
	    FormFactor = 3.0d0*pi**2.0d0/2.0d0 +11.0d0/2.0d0
c	M1sig = FormFactor*Msig -2d0 *( -(AggH*AyyHt)/(s-MH**2d0+dcmplx(0,1)*Mh*Gh) )
        M1sig = FormFactor*Msig
     . -2d0 *( -(AggH*AyyHt)/(s-MH**2d0+dcmplx(0,1)*Mh*Gh) )
c Multiplico por el conjugado para tener la señal:
c 	HARD SCHEME
      FLmmpp=0.5d0
c The minus sign in mppp and ppmp comes from spinor phases.
c      FLpmpp=-1.0d0/8.0d0*((2.0d0+4.0d0*x1/y/y-5.0d0*x1*x1/y/y)*
c     .  (lx1*lx1+2.0d0*dcmplx(0,1)*pi*lx1)
c     .  -(1.0d0-x1*y)*((lx1-ly)*(lx1-ly)+pi*pi)
c     .  +2.0d0*(9.0d0/y-10.0d0*x1)*(lx1+dcmplx(0,1)*pi)
c     .  +(2.0d0+4.0d0*y/x1/x1-5.0d0*y*y/x1/x1)*
c     .  (ly*ly+2.0d0*dcmplx(0,1)*pi*ly)
c     .  -(1.0d0-y*x1)*((lx1-ly)*(lx1-ly)+pi*pi)
c     .  +2.0d0*(9.0d0/x1-10.0d0*y)*(ly+dcmplx(0,1)*pi))
c      FLmmmp=-1.0d0/8.0d0*(  (2.0d0+6.0d0*x1/y/y-3.0d0*x1*x1/y/y)*
c     .  (lx1*lx1+2.0d0*dcmplx(0,1)*pi*lx1)
c     .             -(x1-y)*(x1-y)*((lx1-ly)*(lx1-ly)+pi*pi)
c     .              +2.0d0*(9.0d0/y-8.0d0*x1)*(lx1+dcmplx(0,1)*pi)
c     .              +(2.0d0+6.0d0*y/x1/x1-3.0d0*y*y/x1/x1)*
c     .  (ly*ly+2.0d0*dcmplx(0,1)*pi*ly)
c     .              -(x1-y)*(x1-y)*((lx1-ly)*(lx1-ly)+pi*pi)
c     .              +2.0d0*(9.0d0/x1-8.0d0*y)*(ly+dcmplx(0,1)*pi))
      FLpppp=-(x1*x1+y*y)*(4.0d0*ldli4(-x1)+(ly-3.0d0*lx1)*ldli3(-x1)
     . +lx1*lx1*(ldli2(-x1)-pi*pi/12.0d0)+1.0d0/48.0d0*((lx1+ly)**4.0d0)
     .        -109.0d0/720.0d0*(pi**4.0d0))
     . -2.0d0*x1*y*(ldli3(-x1/y)-(lx1-ly)*ldli2(-x1/y)-Zeta3
     .        +0.5d0*ly*((lx1-ly)*(lx1-ly)+pi*pi))
     .  -0.5d0*x1*x1*(ldli3(-x1/y)-(lx1-ly)*ldli2(-x1/y)-Zeta3
     . -0.5d0*((lx1-ly)**3.0d0)-ly*((lx1-ly)*(lx1-ly)+pi*pi))
     . +1.0d0/8.0d0*(y*y+12.0d0*x1*y-27.0d0*x1*x1-8.0d0/y+9.0d0/y/y)*
     . lx1*lx1-1.0d0/8.0d0*(38.0d0*x1*y-13.0d0)*lx1*ly +pi*pi/48.0d0*
     . (114.0d0*x1*y-43.0d0)-9.0d0/4.0d0*(1.0d0/y+2.0d0*x1)*lx1+0.25d0
     . +dcmplx(0,1)*pi*(-(x1*x1+y*y)*(ldli3(-x1/y)-(lx1-ly)*ldli2(-x1/y)
     .  -Zeta3 +0.5d0*(ly-0.75d0)*((lx1-ly)*(lx1-ly)+pi*pi))
     . +0.25d0*(14.0d0*(x1-y)-8.0d0/y+9.0d0/y/y)*lx1-9.0d0/4.0d0*
     . (1.0d0/y-1.0d0))-(y*y+x1*x1)*(4.0d0*ldli4(-y)+(lx1-3.0d0*ly)*
     . ldli3(-y)+ly*ly*(ldli2(-y)-pi*pi/12.0d0)+1.0d0/48.0d0*
     .   ((ly+lx1)**4.0d0)        -109.0d0/720.0d0*(pi**4.0d0))
     .    -2.0d0*y*x1*(ldli3(-y/x1)-(ly-lx1)*ldli2(-y/x1)-Zeta3
     .           +0.5d0*lx1*((ly-lx1)*(ly-lx1)+pi*pi))
     .    -0.5d0*y*y*(ldli3(-y/x1)-(ly-lx1)*ldli2(-y/x1)-Zeta3
     .  -0.5d0*((ly-lx1)**3.0d0)-lx1*((ly-lx1)*(ly-lx1)+pi*pi))
     . +1.0d0/8.0d0*(x1*x1+12.0d0*y*x1-27.0d0*y*y-8.0d0/x1+9.0d0/x1/x1)
     . *ly*ly-1.0d0/8.0d0*(38.0d0*y*x1-13.0d0)*ly*lx1 +pi*pi/48.0d0*
     . (114.0d0*y*x1-43.0d0)-9.0d0/4.0d0*(1.0d0/x1+2.0d0*y)*ly+0.25d0
     . +dcmplx(0,1)*pi*(-(y*y+x1*x1)*(ldli3(-y/x1)-(ly-lx1)*ldli2(-y/x1)
     .   -Zeta3+0.5d0*(lx1-0.75d0)*((ly-lx1)*(ly-lx1)+pi*pi))
     . +0.25d0*(14.0d0*(y-x1)-8.0d0/x1+9.0d0/x1/x1)*ly-9.0d0/4.0d0*
     . (1.0d0/x1-1.0d0))
c      FLpmmp=2.0d0*(x1*x1+1.0d0)/y/y*(ldli4(-x1/y)-ldli4(-y)+0.5d0*
c     .  (lx1-2.0d0*ly)*(ldli3(-x1)-Zeta3)
c     .    -1.0d0/48.0d0*(lx1**4.0d0)-1.0d0/6.0d0*lx1*ly*ly*ly+1.0d0/
c     .  24.0d0*(ly**4.0d0)
c     .    +pi*pi/24.0d0*(7.0d0*lx1*lx1+2.0d0*ly*ly)+7.0d0/360.0d0*
c     .  (pi**4.0d0))
c     .    +4.0d0*x1*(x1-3.0d0)/y/y*(ldli4(-x1)+ldli4(-x1/y)-ldli4(-y)-ly*
c     .  (ldli3(-x1)-Zeta3)
c     .  +pi*pi/6.0d0*(ldli2(-x1)+0.5d0*ly*ly)-1.0d0/6.0d0*lx1*ly*ly*ly+
c     .  1.0d0/24.0d0*(ly**4.0d0)
c     .    -7.0d0/360.0d0*(pi**4.0d0))
c     .  +2.0d0/3.0d0*(8.0d0-x1+30.0d0*x1/y)*(ldli3(-x1)-ldli3(-y)-(lx1+ly)
c     .  *ldli2(-x1)-0.5d0*lx1*ly*ly)
c     .    -1.0d0/6.0d0*(47.0d0+154.0d0*x1/y-4.0d0*x1*x1/y/y)*(ldli3(-x1)-
c     .  lx1*ldli2(-x1)-Zeta3)
c     .    +1.0d0/12.0d0*(3.0d0-2.0d0/y-12.0d0*x1/y/y)*lx1*
c     .  (lx1*lx1+3.0d0*pi*pi)-1.0d0/3.0d0*y*lx1*ly*ly
c     .    +pi*pi/3.0d0*(y-x1+2.0d0*x1/y)*ly
c     .    +5.0d0/9.0d0*pi*pi*x1*(1.0d0-2.0d0/y-10.0d0*x1/y/y)*lx1+
c     .  2.0d0*(1.0d0+2.0d0/y)*Zeta3
c     .    -1.0d0/24.0d0*(15.0d0-14.0d0*x1/y-48.0d0*x1/y/y)*lx1*lx1
c     .  +1.0d0/24.0d0*(y*y-24.0d0*y+44.0d0-8.0d0*x1*x1*x1/y)*((lx1-ly)
c     .  *(lx1-ly)+pi*pi)
c     .    +4.0d0/9.0d0*pi*pi*x1/y
c     .  +1.0d0/24.0d0*(8.0d0*x1/y+60.0d0-24.0d0*y/x1+27.0d0*y*y/x1/x1)
c     .  *ly*ly
c     .  +1.0d0/12.0d0*(2.0d0*x1*x1-54.0d0*x1-27.0d0*y*y)*(lx1/y+ly/x1)
c     .  +dcmplx(0,1)*pi*((18.0d0*x1/y/y+4.0d0*x1/y-1.0d0)*(ldli3(-x1)
c     .  -Zeta3)
c     .   -(x1*x1+1.0d0)/6.0d0/y/y*lx1*(2.0d0*lx1*lx1-pi*pi)
c     .  -1.0d0/6.0d0*(13.0d0-8.0d0*x1+78.0d0*x1/y+4.0d0/y/y)*ldli2(-x1)
c     .  -1.0d0/3.0d0*x1*(3.0d0+2.0d0*x1/y-10.0d0*x1/y/y)*lx1*lx1
c     .  -1.0d0/3.0d0*y*ly*(ly+2.0d0*lx1)-pi*pi/36.0d0*(7.0d0+4.0d0*x1+
c     .  18.0d0*x1/y-4.0d0/y/y)
c     .  -1.0d0/12.0d0*(15.0d0-14.0d0*x1/y-48.0d0*x1/y/y)*lx1
c     .  +1.0d0/12.0d0*(8.0d0*x1/y+60.0d0-24.0d0*y/x1+27.0d0*y*y/x1/x1)
c     .  *ly
c     .    -1.0d0/12.0d0*(2.0d0*x1/y-54.0d0/y-27.0d0*y/x1))
c      FLpmpm=2.0d0*(y*y+1.0d0)/x1/x1*(ldli4(-y/x1)-ldli4(-x1)+0.5d0*
c     .  (ly-2.0d0*lx1)*(ldli3(-y)-Zeta3)
c     .  -1.0d0/48.0d0*(ly**4.0d0)-1.0d0/6.0d0*ly*lx1*lx1*lx1+1.0d0/
c     .  24.0d0*(lx1**4.0d0)
c     .  +pi*pi/24.0d0*(7.0d0*ly*ly+2.0d0*lx1*lx1)+7.0d0/360.0d0*
c     .  (pi**4.0d0))
c     .  +4.0d0*y*(y-3.0d0)/x1/x1*(ldli4(-y)+ldli4(-y/x1)-ldli4(-x1)-lx1*
c     .  (ldli3(-y)-Zeta3)
c     .  +pi*pi/6.0d0*(ldli2(-y)+0.5d0*lx1*lx1)-1.0d0/6.0d0*ly*lx1*lx1*
c     .  lx1+1.0d0/24.0d0*(lx1**4.0d0)
c     . -7.0d0/360.0d0*(pi**4.0d0))
c     . +2.0d0/3.0d0*(8.0d0-y+30.0d0*y/x1)*(ldli3(-y)-ldli3(-x1)-(ly+lx1)
c     .  *ldli2(-y)-0.5d0*ly*lx1*lx1)
c     . -1.0d0/6.0d0*(47.0d0+154.0d0*y/x1-4.0d0*y*y/x1/x1)*(ldli3(-y)-
c     .  ly*ldli2(-y)-Zeta3)
c     . +1.0d0/12.0d0*(3.0d0-2.0d0/x1-12.0d0*y/x1/x1)*ly*(ly*ly+3.0d0*
c     .  pi*pi)-1.0d0/3.0d0*x1*ly*lx1*lx1
c     . +pi*pi/3.0d0*(x1-y+2.0d0*y/x1)*lx1
c     . +5.0d0/9.0d0*pi*pi*y*(1.0d0-2.0d0/x1-10.0d0*y/x1/x1)*ly+2.0d0*
c     .  (1.0d0+2.0d0/x1)*Zeta3
c     . -1.0d0/24.0d0*(15.0d0-14.0d0*y/x1-48.0d0*y/x1/x1)*ly*ly
c     . +1.0d0/24.0d0*(x1*x1-24.0d0*x1+44.0d0-8.0d0*y*y*y/x1)*((ly-lx1)
c     .  *(ly-lx1)+pi*pi)
c     . +4.0d0/9.0d0*pi*pi*y/x1
c     . +1.0d0/24.0d0*(8.0d0*y/x1+60.0d0-24.0d0*x1/y+27.0d0*x1*x1/y/y)
c     .  *lx1*lx1
c     . +1.0d0/12.0d0*(2*y*y-54.0d0*y-27.0d0*x1*x1)*(ly/x1+lx1/y)
c     . +dcmplx(0,1)*pi*((18.0d0*y/x1/x1+4.0d0*y/x1-1.0d0)*(ldli3(-y)
c     .  -Zeta3)
c     . -(y*y+1.0d0)/6.0d0/x1/x1*ly*(2.0d0*ly*ly-pi*pi)
c     . -1.0d0/6.0d0*(13.0d0-8.0d0*y+78.0d0*y/x1+4.0d0/x1/x1)*ldli2(-y)
c     . -1.0d0/3.0d0*y*(3.0d0+2.0d0*y/x1-10.0d0*y/x1/x1)*ly*ly
c     . -1.0d0/3.0d0*x1*lx1*(lx1+2.0d0*ly)-pi*pi/36.0d0*(7.0d0+4.0d0*y+
c     .  18.0d0*y/x1-4.0d0/x1/x1)
c     . -1.0d0/12.0d0*(15.0d0-14.0d0*y/x1-48.0d0*y/x1/x1)*ly
c     . +1.0d0/12.0d0*(8.0d0*y/x1+60.0d0-24.0d0*x1/y+27*x1*x1/y/y)*lx1
c     . -1.0d0/12.0d0*(2.0d0*y/x1-54.0d0/x1-27.0d0*x1/y))
      FSLmmpp=-1.5d0
c The minus sign in mppp and ppmp comes from spinor phases.
c      FSLpmpp=-1.0d0/8.0d0*( (x1*x1+1.0d0)/y/y*(lx1*lx1+2.0d0*
c     .  dcmplx(0,1)*pi*lx1)
c     .  +0.5d0*(x1*x1+y*y)*((lx1-ly)*(lx1-ly)+pi*pi)
c     .  -4.0d0*(1.0d0/y-x1)*(lx1+dcmplx(0,1)*pi)
c     .  +(y*y+1.0d0)/x1/x1*(ly*ly+2.0d0*dcmplx(0,1)*pi*ly)
c     .  +0.5d0*(x1*x1+y*y)*((lx1-ly)*(lx1-ly)+pi*pi)
c     .  -4.0d0*(1.0d0/x1-y)*(ly+dcmplx(0,1)*pi))
c      FSLmmmp=FSLpmpp
      FSLpppp=-2.0d0*(x1*x1+y*y)*ldli4(-x1)-(x1-y)*ldli4(-x1/y)
     . +2.0d0*x1*
     . x1*lx1*(ldli3(-x1)+ldli3(-y))
     . +pi*pi/6.0d0*(x1-y)*ldli2(-x1)-x1*x1*(1.0d0/6.0d0*(lx1**4.0d0)
     . -2.0d0/3.0d0*lx1*lx1*lx1*ly
     . +pi*pi/6.0d0*lx1*ly-4.0d0/45.0d0*(pi**4.0d0))
     . -x1*(2.0d0*ldli3(-x1)-2.0d0*lx1*ldli2(-x1)-ldli3(-x1/y)+(lx1-ly)*
     . ldli2(-x1/y)-3.0d0*Zeta3
     . -2.0d0/3.0d0*lx1*(lx1*lx1+pi*pi)+5.0d0/12.0d0*(lx1-ly)*
     . ((lx1-ly)*(lx1-ly)+pi*pi)+lx1*lx1*(lx1-ly))
     . -0.25d0*(2.0d0-3.0d0*x1*x1+(x1**4.0d0)/y/y)*lx1*lx1+0.25d0*
     . (2.0d0*x1*y+3.0d0)*lx1*ly
     . -1.0d0/24.0d0*(6.0d0*x1*y+7.0d0)*pi*pi+(0.5d0/y+x1)*lx1-0.25d0
     . +dcmplx(0,1)*pi*(-2.0d0*x1*x1*(ldli3(-y/x1)-Zeta3)+2.0d0*x1*
     . (ldli2(-y/x1)-0.25d0*(lx1-ly)*(lx1-ly)
     . -pi*pi/12.0d0)-0.5d0*(2.0d0*x1*x1-1.0d0)/y/y*lx1+0.5d0/y-0.5d0)
     . -2.0d0*(y*y+x1*x1)*ldli4(-y)-(y-x1)*ldli4(-y/x1)+2.0d0*y*y*ly*
     . (ldli3(-y)+ldli3(-x1))
     . +pi*pi/6.0d0*(y-x1)*ldli2(-y)-y*y*(1.0d0/6.0d0*(ly**4.0d0)-
     . 2.0d0/3.0d0*ly*ly*ly*lx1
     . +pi*pi/6.0d0*ly*lx1-4.0d0/45.0d0*(pi**4.0d0))
     . -y*(2.0d0*ldli3(-y)-2.0d0*ly*ldli2(-y)-ldli3(-y/x1)+(ly-lx1)*
     . ldli2(-y/x1)-3.0d0*Zeta3
     . -2.0d0/3.0d0*ly*(ly*ly+pi*pi)+5.0d0/12.0d0*(ly-lx1)*((ly-lx1)*
     . (ly-lx1)+pi*pi)+ly*ly*(ly-lx1))
     . -0.25d0*(2.0d0-3.0d0*y*y+(y**4.0d0)/x1/x1)*ly*ly+0.25d0*
     . (2.0d0*y*x1+3.0d0)*ly*lx1
     . -1.0d0/24.0d0*(6.0d0*y*x1+7.0d0)*pi*pi+(0.5d0/x1+y)*ly-0.25d0
     . +dcmplx(0,1)*pi*(-2.0d0*y*y*(ldli3(-x1/y)-Zeta3)+2.0d0*y*
     . (ldli2(-x1/y)-0.25d0*(ly-lx1)*(ly-lx1)
     . -pi*pi/12.0d0)-0.5d0*(2.0d0*y*y-1.0d0)/x1/x1*ly+0.5d0/x1-0.5d0)
c      FSLpmmp=-2.0d0*(x1*x1+1.0d0)/y/y*(ldli4(-x1/y)-ldli4(-y)-ly*
c     . (ldli3(-x1)-Zeta3)
c     . +0.5d0*lx1*(ldli3(-x1)-Zeta3)+1.0d0/24.0d0*(lx1**4.0d0)-
c     . 1.0d0/6.0d0*lx1*ly*ly*ly
c     . +1.0d0/24.0d0*(ly**4.0d0)+pi*pi/12.0d0*ly*ly+7.0d0/360.0d0*
c     . (pi**4.0d0))
c     . -2.0d0*(x1-1.0d0)/y*(ldli4(-x1)-0.5d0*lx1*(ldli3(-x1)-Zeta3)
c     . +pi*pi/6.0d0*(ldli2(-x1)-0.5d0*lx1*lx1)
c     . -1.0d0/48.0d0*(lx1**4.0d0)-7.0d0/180.0d0*(pi**4.0d0))
c     . +(2.0d0*x1/y-1.0d0)*(ldli3(-x1)-lx1*ldli2(-x1)+Zeta3-1.0d0/6.0d0*
c     . lx1*lx1*lx1
c     . -pi*pi/3.0d0*(lx1+ly))+2.0d0*(2.0d0*x1/y+1.0d0)*(ldli3(-y)+ly*
c     . ldli2(-x1)-Zeta3
c     . +0.25d0*lx1*(2.0d0*ly*ly+pi*pi)-1.0d0/8.0d0*lx1*lx1*lx1)
c     . -0.25d0*(2.0d0*x1*x1-y*y)*
c     . ((lx1-ly)*(lx1-ly)+pi*pi)-0.25d0*(3.0d0+2.0d0*x1/y/y)*lx1*lx1
c     . -(2.0d0-y*y)/4.0d0/x1/x1*ly*ly+pi*pi/6.0d0+0.5d0*(2.0d0*x1+y*y)
c     . *(lx1/y+ly/x1)-0.5d0
c     . +dcmplx(0,1)*pi*(2.0d0/y/y*(ldli3(-x1)-Zeta3)
c     . -(x1*x1+1.0d0)/6.0d0/y/y*lx1*lx1*lx1
c     . -(2.0d0/y-1.0d0)*ldli2(-x1)-3.0d0/4.0d0*(x1-1.0d0)/y*lx1*lx1
c     . -(x1/y/y+1.5d0)*lx1-(2.0d0-y*y)/2.0d0/x1/x1*ly
c     . -1.0d0/y-y/2.0d0/x1)
c      FSLpmpm=-2.0d0*(y*y+1.0d0)/x1/x1*(ldli4(-y/x1)-ldli4(-x1)-
c     . lx1*(ldli3(-y)-Zeta3)
c     . +0.5d0*ly*(ldli3(-y)-Zeta3)+1.0d0/24.0d0*(ly**4.0d0)
c     . -1.0d0/6.0d0*ly*lx1*lx1*lx1
c     . +1.0d0/24.0d0*(lx1**4.0d0)+pi*pi/12.0d0*lx1*lx1
c     . +7.0d0/360.0d0*(pi**4.0d0))
c     . -2.0d0*(y-1.0d0)/x1*(ldli4(-y)-0.5d0*ly*(ldli3(-y)-Zeta3)
c     . +pi*pi/6.0d0*(ldli2(-y)-0.5d0*ly*ly)
c     . -1.0d0/48.0d0*(ly**4.0d0)-7.0d0/180.0d0*(pi**4.0d0))
c     . +(2.0d0*y/x1-1.0d0)*(ldli3(-y)-ly*ldli2(-y)+Zeta3
c     . -1.0d0/6.0d0*ly*ly*ly
c     . -pi*pi/3.0d0*(ly+lx1))+2.0d0*(2.0d0*y/x1+1.0d0)*
c     . (ldli3(-x1)+lx1*ldli2(-y)-Zeta3
c     . +0.25d0*ly*(2.0d0*lx1*lx1+pi*pi)-1.0d0/8.0d0*ly*ly*ly)
c     . -0.25d0*(2.0d0*y*y-x1*x1)*
c     . ((ly-lx1)*(ly-lx1)+pi*pi)-0.25d0*(3.0d0+2.0d0*y/x1/x1)*ly*ly
c     . -(2.0d0-x1*x1)/4.0d0/y/y*lx1*lx1+pi*pi/6.0d0+0.5d0*
c     . (2.0d0*y+x1*x1)*(ly/x1+lx1/y)-0.5d0
c     . +dcmplx(0,1)*pi*(2.0d0/x1/x1*(ldli3(-y)-Zeta3)
c     . -(y*y+1.0d0)/6.0d0/x1/x1*ly*ly*ly
c     . -(2.0d0/x1-1.0d0)*ldli2(-y)-3.0d0/4.0d0*(y-1.0d0)/x1*ly*ly
c     . -(y/x1/x1+1.5d0)*ly-(2.0d0-x1*x1)/2.0d0/y/y*lx1
c     . -1.0d0/x1-x1/2.0d0/y)
	Poles = 1.0d0/6.0d0*(9.0d0*pi**2.0d0+(33.0d0-2.0d0*Nf)*
     .   dlog(mur**2.d0/s))
c         write(*,*)'Poles : ',Poles
c         write(*,*)'dlog(mur**2.d0/s) = ',dlog(mur**2.d0/s)
c Repito aca la definicion de Mqback_pppp,mmmm,mmpp y ppmm
	Mqback_pppp = -1d0 + yy*dlog((1d0+yy)/(1d0-yy))  
     . -( 1d0 + yy**2d0 )/4d0*((dlog((1d0+yy)/(1d0-yy)))**2d0+Pi**2d0)
	Mqback_mmpp = 1d0
	Mqback_mmmm = Mqback_pppp
	Mqback_ppmm = Mqback_mmpp
	M1back_pppp = Poles*Mqback_pppp+3.d0*FLpppp-1.0d0/3.0d0*FSLpppp
	M1back_mmpp = Poles*Mqback_mmpp+3.d0*FLmmpp-1.0d0/3.0d0*FSLmmpp
	M0back_pppp = Mqback_pppp
	M0back_mmpp = Mqback_mmpp
c        write(*,*)' Mqback_pppp = ', Mqback_pppp
c        write(*,*)' M1back_pppp = ', M1back_pppp 
c        write(*,*)' M1back_mmpp = ', M1back_mmpp
c Minterf es el output de la función BOXyy 
c(incluye todos los factores de norm, incluido el 4pi^2)
c agrego todos los factores de sumar y promediar, y el 4pi^2 q va a mano 
c (alpha y alphas estan en otro lado).
	average = 1d0/(2d0*2d0*8d0*8d0)
	identical = 1d0/2d0
	color = 8d0
	Sum_eq2 = dsqrt(121d0/81d0)
	factors = (4d0*Sum_eq2)**2d0
c Msig_l1l2l3l4 =  Msig si l1=l2 y l3 =l4 o 0 si no.
	M1back0sig_pppp = 2d0*dble(dconjg(M1back_pppp)*Msig)
	M1back0sig_mmpp = 2d0*dble(dconjg(M1back_mmpp)*Msig)
	M1back0sig = 2.d0 * (M1back0sig_pppp + M1back0sig_mmpp) ! times 2 
c because pppp = mmmm and mmpp = ppmm
	M1back0sig = M1back0sig*average*identical*color*dsqrt(factors)*
     .   4d0*Pi**2d0
c        write(*,*)'M1back0sig_pppp =',M1back0sig_pppp
c        write(*,*)'M1back0sig_mmpp =',M1back0sig_mmpp
c        write(*,*)'M1back0sig = ',M1back0sig
c        write(*,*)'M1back0sig = ',M1back0sig
	M0back1sig_pppp = 2d0*dble(dconjg(M1sig)*M0back_pppp)
	M0back1sig_mmpp = 2d0*dble(dconjg(M1sig)*M0back_mmpp)
c        write(*,*)'M0back_pppp = ',M0back_pppp
c        write(*,*)'dconjg(M1sig) = ',dconjg(M1sig)
c         write(*,*)'M0back1sig_pppp = ',M0back1sig_pppp
c          write(*,*)' M0back1sig_mmpp = ', M0back1sig_mmpp
	M0back1sig = 2.d0 * (M0back1sig_pppp + M0back1sig_mmpp) ! times 2 
c because pppp = mmmm and mmpp = ppmm
	M0back1sig = M0back1sig*average*identical*color*dsqrt(factors)*
     .   4d0*Pi**2d0
c Repito la definición de Minterf acá:
	Minterf = 2d0*dble(dconjg(Msig)*(Mqback_pppp+Mqback_mmmm
     .   +Mqback_ppmm+Mqback_mmpp))
	Minterf = Minterf*average*identical*color*dsqrt(factors)*
     .   4d0*Pi**2d0
	H1gYYint_1back0sig = 1.0d0/2.0d0*M1back0sig/Minterf
	H1gYYint_0back1sig = 1.0d0/2.0d0*M0back1sig/Minterf
	H1gYYintDY = 4d0/3d0*(Pi**2/2d0 - 4d0)
	H1gYYint = H1gYYint_1back0sig + H1gYYint_0back1sig 
     .  -0.d0*H1gYYintDY ! El H1DY no tiene nada que ver con el canal gg
        ! Hay un bug en resolve
        H1gYYint =2.d0* H1gYYint ! to adjust alpha_s/2/Pi expansion in DYTurbo
c        write(*,*) 'H1gYYint = ', H1gYYint        
c        
c$$$        write(*,*) ' F AggH = ', AggH
c$$$        write(*,*) ' F AyyH = ', AyyH
c$$$        write(*,*) ' F  AyyHt = ',  AyyHt
c$$$        write(*,*) ' F  Msig = ',  Msig
c$$$        write(*,*) ' F MsigAA = ', MsigAA
c$$$        write(*,*) ' F MsigBB = ', MsigBB
c$$$        write(*,*) ' F M1sig = ', M1sig
c        
c$$$        write(*,*)' Mqback_pppp = ',   Mqback_pppp 
c$$$        write(*,*)' Mqback_mmpp = ',  Mqback_mmpp  
c$$$        write(*,*)' Mqback_mmmm = ',  Mqback_mmmm 
c$$$        write(*,*)'  Mqback_ppmm = ',   Mqback_ppmm 
c$$$        write(*,*)' M1back_pppp = ',  M1back_pppp 
c$$$        write(*,*)' M1back_mmpp = ',  M1back_mmpp 
c$$$        write(*,*)' M0back_pppp = ',  M0back_pppp 
c$$$        write(*,*)'  M0back_mmpp = ',   M0back_mmpp 
c               
c$$$        write(*,*)' M1back0sig_pppp  ', M1back0sig_pppp  
c$$$        write(*,*)' M1back0sig_mmpp  ',  M1back0sig_mmpp  
c$$$        write(*,*)' M1back0sig  ', M1back0sig  
c$$$        write(*,*)'M1back0sig  ', M1back0sig  
c
c$$$        write(*,*)'  M1back0sig = ', M1back0sig
c$$$        write(*,*)' M0back1sig_pppp ', M0back1sig_pppp  
c$$$        write(*,*)'  M0back1sig_mmpp ', M0back1sig_mmpp 
c$$$        write(*,*)'   M0back1sig ', M0back1sig 
c$$$        write(*,*)' Minterf ', Minterf  
c$$$        write(*,*)' H1gYYint_1back0sig ', H1gYYint_1back0sig 
c$$$        write(*,*)'  H1gYYint_0back1sig ', H1gYYint_0back1sig 
c$$$        write(*,*)'   H1gYYint ', H1gYYint  
      return
      end


c------------------------------------------------------------------c
        complex*16 function Martinf(x)
	implicit none
        double precision x,Pi

        Pi = 3.14159265358979323846d0
        
        if (x.ge.1d0) then
            Martinf = (dasin(dsqrt(1d0/x)))**2d0
        else
            Martinf = -1d0/4d0*(dlog((1d0+dsqrt(1d0-x))/
     .      (1d0-dsqrt(1d0-x)))- dcmplx(0,1)*Pi )**2d0
	endif
      	return
        end
c------------------------------------------------------------------c
        complex*16 function MartinF1(x)
	implicit none
        double precision x
	complex*16 Martinf
        external Martinf
	
        MartinF1 = 2d0 + 3d0*x*(1d0+(2d0-x)*Martinf(x))
	return        
	end
c------------------------------------------------------------------c
        complex*16 function MartinF12(x)
	implicit none
        double precision x
	complex*16 Martinf
        external Martinf
 
        MartinF12 = -2d0*x*(1d0+(1d0-x)*Martinf(x))
	return
        end
c------------------------------------------------------------------c

