c Begin of CT contribution Massimiliano
c
      double precision function ctnnlo(q2,xx10,xx20,x1p
     &  ,x2p,xmur,xmuf,xmio,iset,a2,order,LL1,LL2,LL3,LL4,ah1)
     
      implicit double precision (a-z)
      
c      double precision wght,sig0,sfp1(7),sfp2(7),sfp12(7)
c      double precision p(4,6),pw(4),p5(4),p6(4),pout(4,6)
c      double precision sf0(7)
      double precision Z3
C
      integer order,iset
      parameter(Z3=1.20205690315959429d0)
     
      double precision fx10(-5:5),fx20(-5:5),fx1p(-5:5),fx2p(-5:5)
      double precision eq2(5),msqc(-5:5,-5:5)
     
      logical noglue,ggonly,gqonly
c      integer n1,n2
c      integer jproc,i,j,ni,nj,flag,finwidth,ict,
      integer i,j,k,nf,l
c      integer isetproton1,isetproton2
      integer nl,ih1,ih2
      parameter (pi=3.14159265358979312D0)
      
      integer iset2g,ih12g,ih22g
      common/parpdfs2gamma/iset2g,ih12g,ih22g
      iset = iset2g
      
c      parameter (n1=3)
        ih1=ih12g
        ih2=ih22g
        
        vvev_b=0d0
        vvev_1=0d0
        vvev_2=0d0
        sigct1=0d0
        sig_1=0d0
        sig_2=0d0
        sigborn=0d0
        nf=5
        nl=nf
        
c Eliminate particle i=3
        xiilo=0.d0
        yilo=1d0
        phiilo=0d0

c Set g^2 and scales  
c       zg2=zgmu2()      
c Call pdfs
       x1=xx10
       x2=xx20
       
         z1=x1/x1p
         z2=x2/x2p
         
         
     
c
c Q_i^2 for each flavor  
      eq2(1)=4.d0/9.d0
      eq2(4)=4.d0/9.d0
      eq2(2)=1.d0/9.d0
      eq2(3)=1.d0/9.d0
      eq2(5)=1.d0/9.d0
      
      xmur2=xmur**2.d0
      xmuf2=xmuf**2.d0 
      
C SCALE LOGS
      LR=dlog(q2/xmur2)
      LF=dlog(q2/xmuf2) 


        x10=xx10
        x20=xx20
      
      call mlmpdf(iset,ih1,xmuf2,xx10,fx10,5)
      call mlmpdf(iset,ih2,xmuf2,xx20,fx20,5)         

      
      call mlmpdf(iset,ih1,xmuf2,x1p,fx1p,5)
      call mlmpdf(iset,ih2,xmuf2,x2p,fx2p,5)         

      
      noglue=.false.
      ggonly=.false.
      gqonly=.false.

             

CC Switch off gluon !!
      if(noglue) then
        fx10(0)=0d0
        fx20(0)=0d0
        fx1p(0)=0d0
        fx2p(0)=0d0
      endif

CC Gluon only !
      if(ggonly) then
       do j=1,5
       fx10(j)=0d0
       fx10(-j)=0d0
       fx1p(j)=0d0
       fx1p(-j)=0d0
       fx20(j)=0d0
       fx20(-j)=0d0
       fx2p(j)=0d0
       fx2p(-j)=0d0
       enddo
      endif

      flgq=1
CC QG only
       if(gqonly)flgq=0

c         write(6,*)flgq
      
      
      
      xxi1=LL1
      xxi2=LL2
      xxi3=LL3
      xxi4=LL4

        tdelta=0d0
        tH1st=0d0
        tH1stF=0d0
        tgaga=0d0
        tcga=0d0
        tgamma2=0d0

        diffc10=0d0
        diffc1f=0d0
        diffc20=0d0
        diffc2f=0d0

        diffg10=0d0
        diffg1f=0d0
        diffg20=0d0
        diffg2f=0d0

        sig1=0d0
        sig2=0d0

        sig11=0d0
        sig12=0d0
        sig21=0d0      
        sig22=0d0
        sig23=0d0
        sig24=0d0

          do i=-nl,nl 
         do j=-nl,nl    
         msqc(i,j)= 0d0
         enddo
         enddo 
c          write(*,*) a2,LL1,LL2,LL3,LL4

     
c EN MI CODIGO LAS CARGAS NO ESTAN EN LA AMPLITUD ASIQUE HAY QUE PONERLAS
      do i=1,nl    
       msqc(i,-i)= eq2(i)**2 * a2 
       msqc(-i,i)= eq2(i)**2 * a2 
      enddo


CC    Resummation coefficients

      beta0=(33-2*nf)/12d0
      beta1=(153d0-19*nf)/24d0
      Kappa=67/6d0-(pi**2)/2d0-5d0/9d0*nf

      A1q=4d0/3
      A2q=0.5d0*A1q*Kappa
      B1q=-2d0
      B2q=4d0/9*(pi**2-3d0/4-12*Z3)+(11d0/9*pi**2-193d0/12+6*Z3)
     & +nf/6d0*(17d0/3-4d0/9*pi**2)

      
C     Delta term in c1qq coefficient
c      C1qqdelta=(pi**2-8)/3d0   ! Drell-Yan
       C1qqdelta= ah1/2d0/2d0  ! REVISAR  2gamma !!!!!!!!!!!!!!!!!
c      C1qqdelta= ah1/2d0  ! REVISAR  2gamma !!!!!!!!!!!!!!!!!
c       write(*,*) ah1

C     Delta term in P2qq splitting function (as/pi normalization)
      Delta2qq=16d0/9*(3d0/8-pi**2/2+6*Z3)
     &   +4*(17d0/24+11d0*pi**2/18-3*Z3)-2d0/3*nf*(1d0/6+2*pi**2/9d0)
      Delta2qq=Delta2qq/4d0

CC    Coefficients of D0 and D1 in P*P (as/pi normalization)
      D0qqqq=8d0/3
      D1qqqq=32d0/9

CC    Coefficients of delta(1-z) in P*P
      Deltaqqqq=4d0/9*(9d0/4-2*pi**2/3d0)

C     H2qq contribution: coefficient of delta(1-z)
      H2qqdelta=-2561d0/144+127d0*nf/72+3*pi**2/2-19d0*nf*Pi**2/81+
     &          49d0*Pi**4/324 +58d0*Z3/9 + 8d0*nf*Z3/27

C     H2qq contribution: coefficient of D0(z)
      H2qqD0=-404d0/27+(56d0*nf)/81+14*Z3

     

      do j=-nf,nf
      do k=-nf,nf

      if(msqc(j,k).eq.0d0) goto 75


C     Simplest term without convolutions
  
      tdelta=tdelta+fx10(j)*fx20(k)*msqc(j,k)*flgq
      
C     Start H1st: to be used later

C     H1st delta term            
      tH1st=tH1st+2*C1qqdelta*fx10(j)*fx20(k)*msqc(j,k)*flgq

C     H1st: non delta terms, first leg
      tH1st=tH1st+(fx1p(j)*Cqq(z1)*flgq+fx1p(0)*Cqg(z1))
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)

C     H1st: non delta terms, second leg
      tH1st=tH1st+(fx2p(k)*Cqq(z2)*flgq+fx2p(0)*Cqg(z2))         
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
      

C     H1st: muf dependence (LF factor to be added at the end)


c     gammaqq and gammaqg: first leg      
      diff=-dlog(xx10)
     &  *((fx1p(j)-fx10(j)*z1)*Pqq(z1)*flgq+fx1p(0)*Pqg(z1))
      tH1stF=tH1stF+diff*fx20(k)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx10)*fx10(j)*fx20(k)*msqc(j,k)*flgq

c     gammaqq and gammaqg: second leg   
      diff=-dlog(xx20)
     &  *((fx2p(k)-fx20(k)*z2)*Pqq(z2)*flgq+fx2p(0)*Pqg(z2))
      tH1stF=tH1stF+diff*fx10(j)*msqc(j,k)
      tH1stF=tH1stF-Pqqint(xx20)*fx10(j)*fx20(k)*msqc(j,k)*flgq      
      

CC    End of H1st

      if(order.eq.1) goto 75

CC    Now (gamma+gamma)*(gamma+gamma) term: to be used later

C     First part: one gamma for each leg: FLGQ here is non trivial ! DONE
      diffg1f=-dlog(xx10)*(fx1p(j)-fx10(j)*z1)*Pqq(z1)
     &  - Pqqint(xx10)*fx10(j)
      diffg10=-dlog(xx10)*fx1p(0)*Pqg(z1)
      diffg2f=-dlog(xx20)*(fx2p(k)-fx20(k)*z2)*Pqq(z2)
     &  - Pqqint(xx20)*fx20(k)
      diffg20=-dlog(xx20)*fx2p(0)*Pqg(z2)

      tgaga=tgaga+2*
     #   (flgq*diffg10*diffg20+flgq*diffg1f*diffg2f
     #   +diffg10*diffg2f+diffg1f*diffg20)*msqc(j,k)


CC     Second part: gamma*gamma terms
c     Pij * Pjk = D1ijjk (log(1-z)/(1-z))_+ + D0ijjk/(1-z)_+ 
c              + Pijjk(z) + Deltaijjk delta(1-z)
C     First leg
      diff1=-dlog(xx10)*(flgq*(fx1p(j)-fx10(j)*z1)
     &    *(D0qqqq/(1-z1)+D1qqqq*dlog(1-z1)/(1-z1))
     &    +fx1p(j)*Pqqqq(z1)*flgq+fx1p(0)*(Pqqqg(z1)+Pqggg(z1)))
     &    +(Deltaqqqq-D0qqqq*D0int(xx10)-D1qqqq*D1int(xx10))
     &    *fx10(j)*flgq

C    Second leg      
      diff2=-dlog(xx20)*(flgq*(fx2p(k)-fx20(k)*z2)
     &    *(D0qqqq/(1-z2)+D1qqqq*dlog(1-z2)/(1-z2))
     &    +fx2p(k)*Pqqqq(z2)*flgq+fx2p(0)*(Pqqqg(z2)+Pqggg(z2)))
     &    +(Deltaqqqq-D0qqqq*D0int(xx20)-D1qqqq*D1int(xx20))
     &    *fx20(k)*flgq


C     Include Pqggq

      do l=1,nf
      diff1=diff1-dlog(xx10)*(fx1p(l)+fx1p(-l))*Pqggq(z1)*flgq
      diff2=diff2-dlog(xx20)*(fx2p(l)+fx2p(-l))*Pqggq(z2)*flgq
      enddo

      tgaga=tgaga+diff1*fx20(k)*msqc(j,k)
      tgaga=tgaga+diff2*fx10(j)*msqc(j,k)



C    End of (gamma+gamma)*(gamma+gamma) term: FLGQ non trivial here ! DONE

C    Start  (C+C)*(gamma+gamma) term

c    gamma first leg, C second leg
      diffc2f=-dlog(xx20)*fx2p(k)*Cqq(z2)+C1qqdelta*fx20(k)
      diffc20=-dlog(xx20)*fx2p(0)*Cqg(z2)

      tcga=tcga+msqc(j,k)*
     # (flgq*diffg10*diffc20+flgq*diffg1f*diffc2f
     #          +diffg10*diffc2f+diffg1f*diffc20)


c    C first leg, gamma second leg
      diffc1f=-dlog(xx10)*fx1p(j)*Cqq(z1)+C1qqdelta*fx10(j)
      diffc10=-dlog(xx10)*fx1p(0)*Cqg(z1)

      tcga=tcga+msqc(j,k)*
     # (flgq*diffc10*diffg20+flgq*diffc1f*diffg2f
     #          +diffc10*diffg2f+diffc1f*diffg20)
    

c    C*gamma: first leg (ignore delta term in Cqq: taken into account with tH1stF)
      tcga=tcga
     &     +(fx1p(j)*CqqPqq(z1)*flgq+fx1p(0)*(CqqPqg(z1)+CqgPgg(z1)))
     &     *(-dlog(xx10))*fx20(k)*msqc(j,k) 

c    C*gamma: second leg (ignore delta term in Cqq: taken into account with tH1stF)
      tcga=tcga
     &     +(fx2p(k)*CqqPqq(z2)*flgq+fx2p(0)*(CqqPqg(z2)+CqgPgg(z2)))
     &     *(-dlog(xx20))*fx10(j)*msqc(j,k) 

c    Add Cqg*Pgq contribution

      do l=1,nf
      tcga=tcga+(fx1p(l)+fx1p(-l))*CqgPgq(z1)
     &           *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq 
      tcga=tcga+(fx2p(l)+fx2p(-l))*CqgPgq(z2)
     &           *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq 
      enddo

CC  Start 2-loop AP

C   Gluon + pure singlet
      do l=-nf,nf
      if(l.eq.0) then
      tgamma2=tgamma2+fx1p(0)*P2qg(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)
      tgamma2=tgamma2+fx2p(0)*P2qg(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)
      else
      tgamma2=tgamma2+fx1p(l)*P2qqS(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2=tgamma2+fx2p(l)*P2qqS(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      endif
      enddo


C   P2qq non-singlet: regular part
      tgamma2=tgamma2+fx1p(j)*P2qqV(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      tgamma2=tgamma2+fx2p(k)*P2qqV(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq


C   P2qq non-singlet: 1/(1-z)_+
      diff=-dlog(xx10)
     &  *(fx1p(j)-fx10(j)*z1)/(1-z1)
     &  - D0int(xx10)*fx10(j)      
  
      tgamma2=tgamma2+2d0/3*Kappa*diff*fx20(k)*msqc(j,k)*flgq

      diff=-dlog(xx20)
     &  *(fx2p(k)-fx20(k)*z2)/(1-z2)
     &  - D0int(xx20)*fx20(k)      
  
      tgamma2=tgamma2+2d0/3*Kappa*diff*fx10(j)*msqc(j,k)*flgq

      

C   P2qqb non singlet
      tgamma2=tgamma2+fx1p(-j)*P2qqbV(z1)
     & *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

      tgamma2=tgamma2+fx2p(-k)*P2qqbV(z2)
     & *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq

 75   continue


      enddo
      enddo


c!!!!!!!!!!!!!!!CANAL GG
         
c xfactor= (zg2/4/pi**2)**2 *a2*eq2(i)**2

c Leg 1 and 2  gg->qqb
        
      
             xfc=0d0
             do i=1,nl
             xfc=xfc+(-dlog(x10))*(-dlog(x20))*fx1p(0)*fx2p(0)*
     &    a2*eq2(i)**2*( pqg(z1)*pqg(z2)*xxi2 +
     &        (-Cqg(z1)*pqg(z2)-Cqg(z2)*pqg(z1))*xxi1) *2  ! extra 2 for qqbar + qbarq
            enddo
            sigctgg=xfc

C!!!!!!!!!!!!!!!!!!!!!!!








CC   First order      
      sig12=-0.5d0*A1q*tdelta
      sig11=-B1q*tdelta-tH1stF


CC   Second order
      sig24=(A1q)**2/8*tdelta
      sig23=-beta0*A1q/3*tdelta-0.5d0*A1q*sig11
      sig22=0.5d0*(beta0*A1q*LR-A2q)*tdelta
     &     -0.5d0*A1q*(tH1st+LF*tH1stF)
     &     -0.5d0*(B1q-beta0)*sig11
     &     +0.5d0*B1q*tH1stF
     &     +0.5d0*tgaga
      sig21=-beta0*LR*sig11-B1q*(tH1st+LF*tH1stF)
     &     -LF*tgaga-B2q*tdelta+beta0*tH1st-tcga-tgamma2
c     Include missing delta term from C*gamma (no factor 2 here !)
      sig21=sig21-C1qqdelta*tH1stF
C     Include missing term from contact term in 2 loop AP
      sig21=sig21-2*Delta2qq*tdelta

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CC Include as/pi factors and sum O(as) and O(as^2) contributions


      sig1=sig12*xxi2+sig11*xxi1
      sig2=sig24*xxi4+sig23*xxi3+sig22*xxi2+sig21*xxi1


      aspi= (zgmu2(xmur2,iset)/4.d0/pi**2) 
      
      if (order.eq.1) then
         vvev_1=(sig1 *aspi )
      write(*,*)'El orden 1 debe ser implementado solo desde C++!!!!'
      write(*,*)'Poner el flag en el archivo CT.cpp en 2'
      
      else
         vvev_1= (sig2+ sigctgg) *aspi**2
      endif
      
      
       
      

c       EL MENOS VIENE DEL MENOS DEL XJAC DE DAN
c       CHEQUEAR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ctnnlo=-vvev_1  
c       
       return
       end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


CC qq splitting function (with asopi normalization)

      function Pqq(z)
      implicit none
      real *8 Pqq,z
      Pqq=2d0/3*(1+z**2)/(1-z)
      return
      end

CC qg splitting function (with asopi normalization)

      function Pqg(z)
      implicit none
      real *8 Pqg,z
      Pqg=0.25d0*(1-2*z*(1-z))
      return
      end

CC Non delta term in Cqq coefficient (with asopi normalization)

      function Cqq(z)
      implicit none
      real *8 Cqq,z
      Cqq=2d0/3*(1-z)
      return
      end


CC Cqg coefficient (with asopi normalization)

      function Cqg(z)
      implicit none
      real *8 Cqg,z
      Cqg=0.5d0*z*(1-z)
      return
      end


CC Integral of Pqq=1/2 CF (1+x^2)/(1-x) from 0 to z

      function Pqqint(z)
      implicit none
      real *8 Pqqint,z
      Pqqint=-2d0/3*(z+z**2/2+2*dlog(1-z))
      return
      end

CC Integral of 1/(1-x) from 0 to z

      function D0int(z)
      implicit none
      real *8 D0int,z
      D0int=-dlog(1-z)
      return
      end

CC Integral of log(1-x)/(1-x) from 0 to z

      function D1int(z)
      implicit none
      real *8 D1int,z
      D1int=-0.5d0*dlog(1-z)**2
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                P*P convolutions
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Regular part of Pqq*Pqq (checked !)

      function Pqqqq(z)
      implicit none
      real *8 Pqqqq,z
      Pqqqq=4d0/9*(-4*dlog(z)/(1-z)-2*(1-z)
     &  +(1+z)*(3*dlog(z)-4*dlog(1-z)-3))
      return
      end


CC Pqq*Pqg (checked !)

      function Pqqqg(z) 
      implicit none
      real *8 Pqqqg,z
      Pqqqg=1d0/3*((z**2+(1-z)**2)*dlog((1-z)/z)
     &  -(z-0.5d0)*dlog(z)+z-0.25d0)
      return
      end

CC Pqg*Pgq (checked !)

      function Pqggq(z)
      implicit none
      real *8 Pqggq,z
      Pqggq=1d0/3*(2d0/3/z+(1+z)*dlog(z)-2d0/3*z**2-0.5d0*(z-1))
      return
      end


CC Full Pqg*Pgg (checked !)

      function Pqggg(z)
      implicit none
      real *8 Pqggg,z,beta0,Pqg
      integer nf
      external Pqg
      nf=5
      beta0=(33-2*nf)/12d0
      Pqggg=1.5d0*(1/3d0/z+(z**2-z+0.5d0)*dlog(1-z)
     &     +(2*z+0.5d0)*dlog(z)+0.25d0+2*z-31d0/12*z**2)

      Pqggg=Pqggg+beta0*Pqg(z)
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C                C*P convolutions
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Cqq*Pqq (without delta term in Cqq) (checked !)

      function CqqPqq(z)
      implicit none
      real *8 CqqPqq,z
      CqqPqq=2d0/9*(1-z)*(4*dlog(1-z)-2*dlog(z)-1)
      return
      end

CC Cqq*Pqg (without delta term in Cqq) (checked !)

      function CqqPqg(z)
      implicit none
      real *8 CqqPqg,z
      CqqPqg=(-2+z+z**2-(1+2*z)*dlog(z))/6d0
      return
      end

CC Cqg*Pgq (checked !)

      function CqgPgq(z) 
      implicit none
      real *8 CqgPgq,z
      CqgPgq=(1d0/3/z-1+2*z**2/3-z*dlog(z))/3d0
      return
      end

CC Cqg*Pgg (checked !)

      function CqgPgg(z)
      implicit none
      real *8 CqgPgg,z,beta0
      integer nf
      nf=5
      beta0=(33-2*nf)/12d0
      CqgPgg=3d0/4*(2*z*(1-z)*dlog(1-z)-4*z*dlog(z)
     &      +1d0/3/z-1-5*z+17d0*z**2/3)+beta0/2*z*(1-z)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C           Two loop AP:  pqq of ESW is my 3/2 Pqq
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C     Pqq NS: Eq. (4.107) ESW (no 1/(1-x)_+ and delta term)

      function P2qqV(x)
      implicit none
      real *8 x,P2qqV,Pqq,pi
      integer nf
      external Pqq

      pi=3.14159265358979d0
      nf=5

      P2qqV=16d0/9*(-(2*dlog(x)*dlog(1-x)+1.5d0*dlog(x))*3d0/2*Pqq(x)
     &     -(1.5d0+3.5d0*x)*dlog(x)-0.5d0*(1+x)*dlog(x)**2-5*(1-x))
     &     +4*((0.5d0*dlog(x)**2+11d0/6*dlog(x))*3d0/2*Pqq(x)
     &     -(67d0/18-pi**2/6)*(1+x)
     &     +(1+x)*dlog(x)+20d0/3*(1-x))
     &     +2d0/3d0*nf*(-dlog(x)*Pqq(x)+10d0/9*(1+x)-4d0/3*(1-x))

c     Change to as/pi normalization

      P2qqV=P2qqV/4  

      return
      end


C    Pqqb NS: Eq. (4.108) ESW

      function P2qqbV(x)
      implicit none
      real *8 x,P2qqbV,Pqq,S2
      external Pqq,S2

      P2qqbV=-2d0/9*(3d0*Pqq(-x)*S2(x)+2*(1+x)*dlog(x)+4*(1-x))
      
c     Change to as/pi normalization

      P2qqbV=P2qqbV/4 

      return
      end



C    Pqg Singlet: Eq. (4.110) ESW (ESW Pqg is 4 times my Pqg)

      function P2qg(x)
      implicit none
      real *8 x,P2qg,Pqg,pi,S2,logx,logomxsx
      external Pqg,S2

      pi=3.14159265358979d0
      logx=dlog(x)
      logomxsx=dlog((1-x)/x)

      P2qg=2d0/3*(4-9*x-(1-4*x)*logx-(1-2*x)*logx**2+4*dlog(1-x)
     &    +(2*logomxsx**2-4*logomxsx-2d0/3*pi**2+10d0)*4*Pqg(x))
     &    +1.5d0*(182d0/9+14d0/9*x+40d0/9/x+(136d0/3*x-38d0/3)*logx
     &    -4*dlog(1-x)-(2+8*x)*logx**2+8*Pqg(-x)*S2(x)
     &    +(-logx**2+44d0/3*logx-2*dlog(1-x)**2+4*dlog(1-x)+pi**2/3
     &    -218d0/9)*4*Pqg(x))

c     Change to as/pi normalization

      P2qg=P2qg/4d0
  
c     Divide by 2 to eliminate 2nf factor

      P2qg=P2qg/2d0

      return
      end

C     Pqq Pure Singlet appearing in ESW Eq. (4.95)
C     PSqq=PSqqb
C     Obtained through Eq.(4.101)
C     PSqq=1/2/nf (P2qq-P2qqbV-P2qqV) (contains only CF TR=2/3)

      function P2qqS(x)
      implicit none
      real *8 P2qqS,x

      P2qqS=2d0/3*(20 - 18*x + 54*x**2 - 56*x**3
     &    +3*x*(3 + 15*x + 8*x**2)*dlog(x) 
     &    - 9*x*(1 + x)*dlog(x)**2)/(9*x)
      
      P2qqS=P2qqS/4

      return
      end


C    S2: Eq. (4.114) ESW

      function S2(x)
      implicit none
      real *8 x,pi,S2,myli2
      external myli2      
      pi=3.14159265358979d0

      S2=-2*myli2(-x)+0.5d0*dlog(x)**2-2*dlog(x)*dlog(1+x)-pi**2/6
      return
      end

