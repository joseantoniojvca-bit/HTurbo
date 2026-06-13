c Begin of H2 sin el H2f pero con el H1f
c
      double precision function h2(q2,xx10,xx20,z1
     &  ,z2,xmur,xmuf,iset,a2,order,ah1on2pi)
     
      implicit real*8 (a-z)
      external cqq,cqg
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
      integer i,j,k,nf,l
c      integer ict
      integer nl,ih1,ih2,iaa
      parameter (pi=3.14159265358979312D0)
      integer iset2g,ih12g,ih22g
      common/parpdfs2gamma/iset2g,ih12g,ih22g
      iset = iset2g
      ih1=ih12g
      ih2=ih22g
      
c      parameter (n1=3)
      
       
      
      
      
        
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
       
       
        x1p=x1/z1
        x2p=x2/z2
         
         
c         write(*,*) a2,ah1,order,iset
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
c       write(*,*) xmur,xmuf,q2
      
      call mlmpdf(iset,ih1,xmuf2,xx10,fx10,5)
      call mlmpdf(iset,ih2,xmuf2,xx20,fx20,5)         

      
      call mlmpdf(iset,ih1,xmuf2,x1p,fx1p,5)
      call mlmpdf(iset,ih2,xmuf2,x2p,fx2p,5)         

      
      noglue=.false.
      ggonly=.false.
      gqonly=.false.

c                     pepa=0d0
c              do i=-nl,nl 
c               pepa= eq2(i)
c              enddo
c              write(*,*) pepa


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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CC Start calculation

      tdelta=0d0
        tH1st=0d0
        tH1stF=0d0
        tgaga=0d0
        tcga=0d0
        tgamma2=0d0
        tH2st=0d0
        diffc10=0d0
        diffc1f=0d0
        diffc20=0d0
        diffc2f=0d0
        diff=0d0 
        diffg10=0d0
        diffg1f=0d0
        diffg20=0d0
        diffg2f=0d0
        diff2=0d0
        diff1=0d0
        diffc2f=0d0        
        diffc0f=0d0


        sig1=0d0
        sig2=0d0

        sig11=0d0
        sig12=0d0
        sig21=0d0      
        sig22=0d0
        sig23=0d0
        sig24=0d0
        xmsq=0d0 
        

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
      

      
C     Delta term in c1qq coefficient
c      C1qqdelta=(pi**2-8)/3d0   ! Drell-Yan
         ! REVISAR  2gamma !!!!!!!!!!!!!!!!!
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
      H2qqdelta=0.d0
C     H2qq contribution: coefficient of D0(z)
      H2qqD0=-404d0/27+(56d0*nf)/81+14*Z3

     
c         write(*,*) order
      do j=-nf,nf
      do k=-nf,nf

      
      if(msqc(j,k).eq.0d0) goto 75

          


C     Simplest term without convolutions
  
      tdelta=tdelta+fx10(j)*fx20(k)*msqc(j,k)*flgq

c      if(order.eq.0) goto 75

C     Start H1st: to be used later

C     H1st delta term
c      C1qqdelta= ah1on2pi/2d0/2d0

      C1qqdelta = 0.d0 
      ah1=ah1on2pi/2.0d0
c       ah1 = 0.d0
      tH1st=tH1st+ah1*fx10(j)*fx20(k)*msqc(j,k)*flgq
c      tH1st=tH1st+2*C1qqdelta*fx10(j)*fx20(k)*msqc(j,k)*flgq
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

      if(order.eq.1) goto 75

CC    End of H1st

CC    Start H2 contribution

CC    H2st gg contribution

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)*flgq

CC    H2st qqbar contribution from C1*C1 (without delta term)

C     regular*regular

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)*flgq

C     regular-delta

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx20(k)*C1qqdelta*msqc(j,k)*flgq       

      tH2st=tH2st+fx2p(k)*Cqq(z2)*(-dlog(xx20))*
     &            fx10(j)*C1qqdelta*msqc(j,k)*flgq       


CC    H2st qg contribution from C1*C1

C     regular*regular

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx2p(k)*Cqq(z2)*(-dlog(xx20))*msqc(j,k)

      tH2st=tH2st+fx1p(j)*Cqq(z1)*(-dlog(xx10))*
     &            fx2p(0)*Cqg(z2)*(-dlog(xx20))*msqc(j,k)


C     regular-delta

      tH2st=tH2st+fx1p(0)*Cqg(z1)*(-dlog(xx10))*
     &            fx20(k)*C1qqdelta*msqc(j,k)       

      tH2st=tH2st+fx2p(0)*Cqg(z2)*(-dlog(xx20))*
     &            fx10(j)*C1qqdelta*msqc(j,k)    

CC    H2st qqbar channel: D0(z), first leg

      diff=-dlog(xx10)*(fx1p(j)-fx10(j)*z1)*H2qqD0/(1-z1)

      tH2st=tH2st+0.5d0*diff*fx20(k)*msqc(j,k)*flgq
      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx10)
     &            *fx10(j)*fx20(k)*msqc(j,k)*flgq

CC    H2st, qqbar channel: D0(z), second leg
      
      diff=-dlog(xx20)*(fx2p(k)-fx20(k)*z2)*H2qqD0/(1-z2)

      tH2st=tH2st+0.5d0*diff*fx10(j)*msqc(j,k)*flgq
      tH2st=tH2st-0.5d0*H2qqD0*D0int(xx20)
     &            *fx10(j)*fx20(k)*msqc(j,k)*flgq

CC    C2qq, regular part, first leg

      tH2st=tH2st+fx1p(j)*C2qqreg(z1)
     &                   *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

CC    C2qq, regular part, second leg

      tH2st=tH2st+fx2p(k)*C2qqreg(z2)
     &                   *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq

CC    C2qg, first leg

      tH2st=tH2st+fx1p(0)*C2qg(z1)*(-dlog(xx10))*fx20(k)*msqc(j,k)

CC    C2qg, second leg

      tH2st=tH2st+fx2p(0)*C2qg(z2)*(-dlog(xx20))*fx10(j)*msqc(j,k)

CC    Cqqbar contribution: first leg

      tH2st=tH2st+fx1p(-j)*C2qqb(z1)
     &                    *(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq

CC    Cqqbar contribution: second leg
  
      tH2st=tH2st+fx2p(-k)*C2qqb(z2)
     &                    *(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq  

      do a=1,nf

CC    Cqqp contribution: first leg

      if(a.ne.abs(j)) then      
       tH2st=tH2st+(fx1p(a)+fx1p(-a))*
     &       C2qqp(z1)*(-dlog(xx10))*fx20(k)*msqc(j,k)*flgq
      endif

CC    Cqqp contribution: second leg

      if(a.ne.abs(k)) then      
       tH2st=tH2st+(fx2p(a)+fx2p(-a))*
     &       C2qqp(z2)*(-dlog(xx20))*fx10(j)*msqc(j,k)*flgq
      endif

      enddo

CCCC Terms needed for NNLO scale dependence  CCCCCC


CC    (gamma+gamma)*(gamma+gamma) term

C     First part: one gamma for each leg


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



C    End of (gamma+gamma)*(gamma+gamma) term

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


CCCCCCCCCCCC   End of NNLO scale dependence CCCCCCCCCCCCCCCCC


 75   continue

      enddo
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


       asopi= (zgmu2(xmur2,iset)/4.d0/pi**2) 
      

 72   
cxmsq=tdelta
c                write(*,*) tgamma2,flgq
      xmsq=0d0
      if(order.eq.1)then
       xmsq=xmsq+asopi*(tH1st+LF*tH1stF)
                write(*,*) tgamma2,flgq
      elseif(order.eq.2)then
       xmsq=xmsq+asopi**2*(tH2st)
c       xmsq=xmsq+asopi**2*(tdelta*H2qqdelta)

CC     add scale dependence at NNLO
c               write(*,*) xmsq
       xmsq=xmsq+asopi**2*(0.5d0*beta0*LF**2*tH1stF
     &                   +tgamma2*LF
     &                   -beta0*LR*(tH1st+LF*tH1stF)
     &                   +LF*tcga+0.5d0*LF**2*tgaga)

c               write(*,*) xmsq
C     Include missing delta term from C*gamma (no factor 2 here !)

c      xmsq=xmsq+asopi**2*(LF*C1qqdelta*tH1stF)

C     Include missing term from contact term in 2 loop AP

      xmsq=xmsq+asopi**2*(2*Delta2qq*tdelta)*LF

C     Include missing terms from Hard Scheme
      xmsq=xmsq+asopi**2*(ah1*tH1st)
      xmsq=xmsq+asopi**2*(ah1*LF*tH1stF)

      endif     


      
c              write(*,*) asopi**2
        h2=xmsq  
c        write(*,*) C1qqdelta
      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC





      double precision function C2qqreg(z)
      implicit none
      real *8 Pi,Z2,Z3,myli2,myli3,z,CA,CF
      integer nf

      external myli2,myli3

      Pi=3.14159265358979d0
      z2=Pi**2/6

      Z3=1.20205690316d0

      CF=4d0/3
      CA=3d0
      nf=5

      C2qqreg=   
     & (CF*(-344+24*Pi**2+974*z-1600*CA*z+2052*CF*z+148*nf*z-60*Pi**2*z+
     & 108*CF*Pi**2*z-1188*z**2+1584*CA*z**2-4104*CF*z**2-72*nf*z**2+
     & 72*Pi**2*z**2-216*CF*Pi**2*z**2+830*z**3+16*CA*z**3+2052*CF*z**3-
     & 76*nf*z**3-60*Pi**2*z**3+108*CF*Pi**2*z**3-
     & 272*z**4+24*Pi**2*z**4+
     & 324*CA*z*z2-1728*CF*z*z2 - 648*CA*z**2*z2 + 3456*CF*z**2*z2+
     & 324*CA*z**3*z2-1728*CF*z**3*z2 + 1188*CA*z*Z3 + 864*CF*z*Z3-
     & 324*CA*z**3*Z3-864*CF*z**3*Z3 - 108*CA*z**2*dlog(1-z) +
     & 108*CF*z**2*dlog(1-z)+108*CA*z**3*dlog(1-z)-
     & 108*CF*z**3*dlog(1-z)-
     & 216*CA*z*z2*dlog(1-z) + 216*CF*z*z2*dlog(1-z) -
     & 216*CA*z**3*z2*dlog(1-z)+216*CF*z**3*z2*dlog(1-z)-252*z*dlog(z)+
     & 348*CA*z*dlog(z)-540*CF*z*dlog(z)-
     & 60*nf*z*dlog(z)+612*z**2*dlog(z)-
     & 432*CA*z**2*dlog(z)+1404*CF*z**2*dlog(z)-744*z**3*dlog(z)+
     & 996*CA*z**3*dlog(z)-1728*CF*z**3*dlog(z)-60*nf*z**3*dlog(z)+
     & 384*z**4*dlog(z)-144*dlog(1-z)*dlog(z)+360*z*dlog(1-z)*dlog(z)-
     & 216*CA*z*dlog(1-z)*dlog(z) + 648*CF*z*dlog(1-z)*dlog(z) -
     & 432*z**2*dlog(1-z)*dlog(z) + 432*CA*z**2*dlog(1-z)*dlog(z) -
     & 1296*CF*z**2*dlog(1-z)*dlog(z) + 360*z**3*dlog(1-z)*dlog(z) -
     & 216*CA*z**3*dlog(1-z)*dlog(z) + 648*CF*z**3*dlog(1-z)*dlog(z) -
     & 144*z**4*dlog(1-z)*dlog(z)+216*CA*z*dlog(1-z)**2*dlog(z) -
     & 324*CF*z*dlog(1-z)**2*dlog(z)+216*CA*z**3*dlog(1-z)**2*dlog(z) -
     & 324*CF*z**3*dlog(1-z)**2*dlog(z)+27*z*dlog(z)**2+
     & 99*CA*z*dlog(z)**2-
     & 162*CF*z*dlog(z)**2-18*nf*z*dlog(z)**2+108*CA*z**2*dlog(z)**2 -
     & 108*CF*z**2*dlog(z)**2+45*z**3*dlog(z)**2-9*CA*z**3*dlog(z)**2+
     & 108*CF*z**3*dlog(z)**2-18*nf*z**3*dlog(z)**2-72*z**4*dlog(z)**2-
     & 108*CF*z*dlog(1-z)*dlog(z)**2-108*CF*z**3*dlog(1-z)*dlog(z)**2-
     & 18*z*dlog(z)**3 + 18*CA*z*dlog(z)**3-
     & 18*CF*z*dlog(z)**3+18*z**3*dlog(z)**3 +
     & 18*CA*z**3*dlog(z)**3+18*CF*z**3*dlog(z)**3 -
     & 72*((-1 + z)**2*(2 + (-1+3*CA-6*CF)*z + 2*z**2) -
     & 3*(CA-CF)*z*(1 + z**2)*dlog(1-z)-3*(CA-3*CF)*z*(1+z**2)*dlog(z))*
     & myli2(z) + 216*CA*z*myli3(1-z)-216*CF*z*myli3(1-z) +
     & 216*CA*z**3*myli3(1-z)-216*CF*z**3*myli3(1-z) -
     & 432*CA*z*myli3(z) + 1080*CF*z*myli3(z) -
     & 432*CA*z**3*myli3(z)+1080*CF*z**3*myli3(z)-1944*CF*z*Z3-
     & 216*CF*z**3*Z3))/(864*(-1 + z)*z)
     & -1/4.d0*CF**2*(Pi**2/2.d0-4.d0)*(1-z)
c ultima linea a remover si no es Hard scheme     




        return
        end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function C2qqp(z)
      implicit none
      real *8 Pi,myli2,z,CF
c      integer nf

      external myli2

      Pi=3.14159265358979d0


      CF=4d0/3
     

      C2qqp=(CF*(2*(-1+z)*(-172+143*z-136*z**2+6*Pi**2*(2-z+2*z**2))-
     &   12*(z*(-21 + 30*z - 32*z**2)+
     &  6*(-2+3*z-3*z**2+2*z**3)*dlog(1-z))*
     &  dlog(z)-9*z*(3+3*z+8*z**2)*dlog(z)**2+18*z*(1 + z)*dlog(z)**3-
     &  72*(-2 + 3*z - 3*z**2 + 2*z**3)*myli2(z)))/(864d0*z)

 
      return
      end   

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision function C2qqb(z)
      implicit none
      real *8 Pi,Z3,myli2,myli3,z,CF,CA,C2qqp
    

      external myli2,myli3,C2qqp

      Pi=3.14159265358979d0
      Z3=1.20205690316d0

      CF=4d0/3   
      CA=3d0


       C2qqb=C2qqp(z)+
     &  (CF*(-CA+2*CF)*(45-3*Pi**2-2*Pi**2*z-45*z**2+Pi**2*z**2+9*
     &  dlog(z)+42*z*dlog(z)+33*z**2*dlog(z)+12*dlog(1-z)*dlog(z)-
     &  12*z**2*dlog(1-z)*dlog(z)-dlog(z)**3-z**2*
     &  dlog(z)**3+2*Pi**2*dlog(1+z) +
     &  2*Pi**2*z**2*dlog(1+z)-12*dlog(z)*
     &  dlog(1+z)-24*z*dlog(z)*dlog(1+z)-
     &  12*z**2*dlog(z)*dlog(1+z)+6*dlog(z)**2*dlog(1+z)+
     &  6*z**2*dlog(z)**2*dlog(1+z)-4*dlog(1+z)**3-4*z**2*dlog(1+z)**3-
     &  12*((1+z)**2+(1+z**2)*dlog(z))*myli2(-z)-
     &  12*(-1+z**2+dlog(z)+z**2*dlog(z))*myli2(z)+36*myli3(-z)+
     &  36*z**2*myli3(-z)+24*myli3(z)+24*z**2*myli3(z)+
     &  24*myli3(1d0/(1+z))+24*z**2*myli3(1d0/(1+z))-
     &  18*Z3-18*z**2*Z3))/(48*(1+z))


      return
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     NEW version August 2013
C

      double precision function C2qg(z)
      implicit none
      real *8 z,CF,CA,Pi,Z3,myli2,myli3
      external myli2,myli3


      Pi=3.14159265358979d0
      Z3=1.20205690316d0
      CF=4d0/3
      CA=3d0


      C2qg=(688*CA-1260*CA*z-702*CF*z+1548*CA*z**2+2322*CF*z**2+ 
     -    72*CA*Pi**2*z**2+144*CF*Pi**2*z**2-1192*CA*z**3-2160*CF*z**3- 
     -    144*CF*Pi**2*z**3+324*CA*z**2*Log(1-z)-324*CF*z**2*Log(1-z)- 
     -    432*CA*z**3*Log(1-z)+432*CF*z**3*Log(1-z)- 
     -    216*CA*z**2*Log(1-z)**2+216*CF*z**2*Log(1-z)**2+ 
     -    216*CA*z**3*Log(1-z)**2-216*CF*z**3*Log(1-z)**2 + 
     -    36*CA*z*Log(1 - z)**3 - 36*CF*z*Log(1 - z)**3 - 
     -    72*CA*z**2*Log(1-z)**3 + 72*CF*z**2*Log(1-z)**3 + 
     -    72*CA*z**3*Log(1-z)**3-72*CF*z**3*Log(1-z)**3 + 
     -    504*CA*z*Log(z)+432*CF*z*Log(z)-720*CA*z**2*Log(z) + 
     -    810*CF*z**2*Log(z)+144*CA*Pi**2*z**2*Log(z)
     -    +1632*CA*z**3*Log(z)- 
     -    432*CF*z**3*Log(z) - 432*CF*z**2*Log(1-z)*Log(z) + 
     -    432*CF*z**3*Log(1-z)*Log(z)+108*CF*z*Log(1-z)**2*Log(z)- 
     -    216*CF*z**2*Log(1-z)**2*Log(z)+
     -    216*CF*z**3*Log(1-z)**2*Log(z)- 
     -    54*CA*z*Log(z)**2+27*CF*z*Log(z)**2+216*CA*z**2*Log(z)**2+ 
     -    324*CF*z**2*Log(z)**2 - 792*CA*z**3*Log(z)**2 - 
     -    216*CF*z**3*Log(z)**2 + 108*CF*z*Log(1 - z)*Log(z)**2- 
     -    864*CA*z**2*Log(1-z)*Log(z)**2-
     -    216*CF*z**2*Log(1-z)*Log(z)**2+ 
     -    216*CF*z**3*Log(1-z)*Log(z)**2+36*CA*z*Log(z)**3 - 
     -    18*CF*z*Log(z)**3+72*CA*z**2*Log(z)**3+36*CF*z**2*Log(z)**3- 
     -    72*CF*z**3*Log(z)**3 + 36*CA*Pi**2*z*Log(1 + z) + 
     -    72*CA*Pi**2*z**2*Log(1 + z) + 72*CA*Pi**2*z**3*Log(1+z)+ 
     -    432*CA*z**2*Log(z)*Log(1 + z) + 432*CA*z**3*Log(z)*Log(1+z)+ 
     -    108*CA*z*Log(z)**2*Log(1+z)+216*CA*z**2*Log(z)**2*Log(1+z)+ 
     -    216*CA*z**3*Log(z)**2*Log(1+z)-72*CA*z*Log(1+z)**3 - 
     -    144*CA*z**2*Log(1 + z)**3 - 144*CA*z**3*Log(1 + z)**3 - 
     -    72*(3*(CA - CF)*z*(1 - 2*z + 2*z**2)*Log(1 - z) + 
     -       2*CA*(2 - 3*z + 12*z**2 - 11*z**3 + 6*z**2*Log(z)))*
     -     myLi2(1-z) - 216*CA*z*
     -     (-2*z*(1 + z) + (1 + 2*z + 2*z**2)*Log(z))*myli2(-z)+ 
     -    216*CF*z*Log(z)*myli2(z)-
     -    1728*CA*z**2*Log(z)*myli2(z)- 
     -    432*CF*z**2*Log(z)*myli2(z)+432*CF*z**3*Log(z)*myli2(z)+ 
     -    216*CA*z*myli3(1-z)-216*CF*z*myli3(1-z)- 
     -    432*CA*z**2*myli3(1-z)+432*CF*z**2*myli3(1-z) + 
     -    432*CA*z**3*myli3(1-z) - 432*CF*z**3*myli3(1-z) + 
     -    648*CA*z*myli3(-z) + 1296*CA*z**2*myli3(-z)+ 
     -    1296*CA*z**3*myli3(-z) - 216*CF*z*myli3(z) + 
     -    1728*CA*z**2*myli3(z) + 432*CF*z**2*myli3(z)- 
     -    432*CF*z**3*myli3(z) + 432*CA*z*myli3(1/(1+z))+ 
     -    864*CA*z**2*myli3(1/(1+z)) +
     -    864*CA*z**3*myli3(1/(1+z)) - 
     -    648*CA*z*Z3 + 1728*CF*z*Z3-3456*CF*z**2*Z3 - 
     -    1296*CA*z**3*Z3 + 3456*CF*z**3*Z3)/(1728d0*z)
     -   -1/4.d0*CF*(Pi**2/2.d0-4.d0)*z*(1-z)
c ultima linea a remover si no es en el Hard Scheme

      return
      end



    



