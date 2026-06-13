!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine as3vj_qqb_z1jet_v(p,msq)
      implicit none
************************************************************************
*     Authors: R.K. Ellis and John Campbell                            *
*     May, 2001.                                                       *
*     Matrix element for Z + jet production                            *
*     in order alpha_s^2                                               *
*     averaged over initial colours and spins                          *
*     q(-p1)+qbar(-p2)-->Z^+(l(p3)+a(p4))+g(p5)                        *
************************************************************************
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'as3vjet_masses.f'
      include 'as3vjet_qcdcouple.f'
      include 'as3vjet_ewcharge.f'
      include 'as3vjet_zcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'nflav.f'
      include 'cplx.h'
      include 'zcouple_cms.f'
      include 'blha.f'
      include 'as3vjet_constants.f'
      include 'toploops.f'
      integer:: j,k,nu,om
      real(dp):: msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     & p(as3vj_mxpart,4),fac,sz,virt5,subuv
      complex(dp):: zsin2winv
      complex(dp):: qqbZgLL(2),qqbZgRR(2),qqbZgLR(2),qqbZgRL(2)
      complex(dp):: gqZqLL(2),gqZqRR(2),gqZqLR(2),gqZqRL(2)
      complex(dp):: qgZqLL(2),qgZqRR(2),qgZqLR(2),qgZqRL(2)
      complex(dp):: qbqZgLL(2),qbqZgRR(2),qbqZgLR(2),qbqZgRL(2)
      complex(dp):: gqbZqbLL(2),gqbZqbRR(2),gqbZqbLR(2),gqbZqbRL(2)
      complex(dp):: qbgZqbLL(2),qbgZqbRR(2),qbgZqbLR(2),qbgZqbRL(2)
      complex(dp):: prop,virt5ax,cLL(nf),cLR(nf),cRL(nf),cRR(nf)
!      real(dp):: phi,muk,rho,ssig,csig,theta,mtsq,musq,
!     & p1true(4),p2true(4),p3true(4),p4true(4),p5true(4)
      logical,parameter:: includeanom=.true.
      integer,parameter::
     & iqqbgLL(5)=(/1,2,3,4,5/),
     & iqqbgLR(5)=(/1,2,4,3,5/),
     & iqgqLL(5)=(/1,5,3,4,2/),
     & iqgqLR(5)=(/1,5,4,3,2/),
     & igqqLL(5)=(/2,5,3,4,1/),
     & igqqLR(5)=(/2,5,4,3,1/)
      common/virt5ax/virt5ax
!$omp threadprivate(/virt5ax/)

      scheme='dred'

c--set msq=0 to initialize
      msq(:,:)=zero

c--calculate spinor and dot-products (using BDK type notation)
      call mcfm_spinoru(5,p,za,zb)
!      call spinoru6(p,za,zb)

c BEGIN checking code
c        include 'kinpoint5.f'
c        mu=1._dp
c        musq=mu**2
c        mt=0.4255266775_dp
c!        mtsq=mt**2
c        write(6,*) 'mt=',mt
c        do nu=1,4
c        om=nu-1
c        if (nu==1) om=4
c        p(1,om)=p1true(nu)
c        p(2,om)=p2true(nu)
c        p(3,om)=p4true(nu)
c        p(4,om)=p5true(nu)
c        p(5,om)=p3true(nu)
c        enddo
c        do nu=1,5
c        write(6,'(a4,i2,4f18.12)') 'p_',nu,
c     &   p(nu,4),p(nu,1),p(nu,2),p(nu,3)
c        enddo
c        call spinorz(5,p,za,zb)
c END checking code

c--- calculate lowest order
      call as3vj_qqb_z1jet(p,msq0)

c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      subuv=ason2pi*xn
     & *(epinv*(11._dp-2._dp*real(nflav,dp)/xn)-1._dp)/6._dp   

c--   calculate propagator
      sz=s(3,4)
      prop=sz/cplx2((sz-zmass**2),zmass*zwidth)

      fac=8._dp*cf*xnsq*abs(zesq)**2*gsq

c     first letter L/R is the index of the fermion line
c     second letter L/R is the index of the Z/gamma
c     already summed over gluons
      qqbZgLL(1)=aveqq*fac*cplx1(virt5(iqqbgLL,za,zb,.true.))
      qqbZgLL(2)=-aveqq*fac*virt5ax
      qqbZgLR(1)=aveqq*fac*cplx1(virt5(iqqbgLR,za,zb,.true.))
      qqbZgLR(2)=-aveqq*fac*virt5ax
      qqbZgRL(1)=aveqq*fac*cplx1(virt5(iqqbgLR,zb,za,.true.))
      qqbZgRL(2)=aveqq*fac*virt5ax
      qqbZgRR(1)=aveqq*fac*cplx1(virt5(iqqbgLL,zb,za,.true.))
      qqbZgRR(2)=aveqq*fac*virt5ax

      if (onlyaxial) then
        qqbZgLL(1)=0._dp
        qqbZgLR(1)=0._dp
        qqbZgRL(1)=0._dp
        qqbZgRR(1)=0._dp
        subuv=0._dp
      endif

      qbqZgLL(1)=qqbZgRL(1)
      qbqZgLR(1)=qqbZgRR(1)
      qbqZgRL(1)=qqbZgLL(1)
      qbqZgRR(1)=qqbZgLR(1)
      qbqZgLL(2)=-qqbZgRL(2)
      qbqZgLR(2)=-qqbZgRR(2)
      qbqZgRL(2)=-qqbZgLL(2)
      qbqZgRR(2)=-qqbZgLR(2)

      if (useblha==0) then
      gqZqLL(1)=aveqg*fac*cplx1(virt5(igqqLL,za,zb,.true.))
      gqZqLL(2)=-aveqg*fac*virt5ax
      gqZqLR(1)=aveqg*fac*cplx1(virt5(igqqLR,za,zb,.true.))
      gqZqLR(2)=-aveqg*fac*virt5ax
      gqZqRL(1)=aveqg*fac*cplx1(virt5(igqqLR,zb,za,.true.))
      gqZqRL(2)=aveqg*fac*virt5ax
      gqZqRR(1)=aveqg*fac*cplx1(virt5(igqqLL,zb,za,.true.))
      gqZqRR(2)=aveqg*fac*virt5ax
      endif

      if (onlyaxial) then
        gqZqLL(1)=0._dp
        gqZqLR(1)=0._dp
        gqZqRL(1)=0._dp
        gqZqRR(1)=0._dp
        subuv=0._dp
      endif

      gqbZqbRL(1)=gqZqLL(1)
      gqbZqbRR(1)=gqZqLR(1)
      gqbZqbLL(1)=gqZqRL(1)
      gqbZqbLR(1)=gqZqRR(1)
      gqbZqbRL(2)=-gqZqLL(2)
      gqbZqbRR(2)=-gqZqLR(2)
      gqbZqbLL(2)=-gqZqRL(2)
      gqbZqbLR(2)=-gqZqRR(2)


      if (useblha==0) then
      qgZqLL(1)=aveqg*fac*cplx1(virt5(iqgqLL,za,zb,.true.))
      qgZqLL(2)=-aveqg*fac*virt5ax
      qgZqLR(1)=aveqg*fac*cplx1(virt5(iqgqLR,za,zb,.true.))
      qgZqLR(2)=-aveqg*fac*virt5ax
      qgZqRL(1)=aveqg*fac*cplx1(virt5(iqgqLR,zb,za,.true.))
      qgZqRL(2)=aveqg*fac*virt5ax
      qgZqRR(1)=aveqg*fac*cplx1(virt5(iqgqLL,zb,za,.true.))
      qgZqRR(2)=aveqg*fac*virt5ax
      endif

      if (onlyaxial) then
        qgZqLL(1)=0._dp
        qgZqLR(1)=0._dp
        qgZqRL(1)=0._dp
        qgZqRR(1)=0._dp
        subuv=0._dp
      endif

      qbgZqbRL(1)=qgZqLL(1)
      qbgZqbRR(1)=qgZqLR(1)
      qbgZqbLL(1)=qgZqRL(1)
      qbgZqbLR(1)=qgZqRR(1)
      qbgZqbRL(2)=-qgZqLL(2)
      qbgZqbRR(2)=-qgZqLR(2)
      qbgZqbLL(2)=-qgZqRL(2)
      qbgZqbLR(2)=-qgZqRR(2)

      zsin2winv=zL(2)-zR(2)
      do j=1,nf
      cLL(j)=Q(j)*q1+zL(j)*zl1*prop
      cRR(j)=Q(j)*q1+zR(j)*zr1*prop
      cLR(j)=Q(j)*q1+zL(j)*zr1*prop
      cRL(j)=Q(j)*q1+zR(j)*zl1*prop
      enddo

      do j=-nflav,nflav
      do k=-nflav,nflav
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j == 0) .and. (k == 0)) then
         msq(j,k)=0._dp
      elseif ((j > 0) .and. (k < 0)) then
         msq(j,k)=+real(cLL(j)*conjg(cLL(j))*qqbZgLL(1),dp)
     &            +real(cRR(j)*conjg(cRR(j))*qqbZgRR(1),dp)
     &            +real(cLR(j)*conjg(cLR(j))*qqbZgLR(1),dp)
     &            +real(cRL(j)*conjg(cRL(j))*qqbZgRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(zsin2winv*zl1*prop*(qqbZgLL(2)*conjg(cLL(j))+qqbZgRL(2)*conjg(cRL(j))),dp)
     &            +real(zsin2winv*zr1*prop*(qqbZgRR(2)*conjg(cRR(j))+qqbZgLR(2)*conjg(cLR(j))),dp)
         endif
      elseif ((j < 0) .and. (k > 0)) then
         msq(j,k)=+real(cLL(k)*conjg(cLL(k))*qbqZgLL(1),dp)
     &            +real(cRR(k)*conjg(cRR(k))*qbqZgRR(1),dp)
     &            +real(cLR(k)*conjg(cLR(k))*qbqZgLR(1),dp)
     &            +real(cRL(k)*conjg(cRL(k))*qbqZgRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(zsin2winv*zl1*prop*(qbqZgLL(2)*conjg(cLL(k))+qbqZgRL(2)*conjg(cRL(k))),dp)
     &            +real(zsin2winv*zr1*prop*(qbqZgRR(2)*conjg(cRR(k))+qbqZgLR(2)*conjg(cLR(k))),dp)
         endif
      elseif ((j > 0) .and. (k == 0)) then
         msq(j,k)=+real(cLL(j)*conjg(cLL(j))*qgZqLL(1),dp)
     &            +real(cRR(j)*conjg(cRR(j))*qgZqRR(1),dp)
     &            +real(cLR(j)*conjg(cLR(j))*qgZqLR(1),dp)
     &            +real(cRL(j)*conjg(cRL(j))*qgZqRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(zsin2winv*zl1*prop*(qgZqLL(2)*conjg(cLL(j))+qgZqRL(2)*conjg(cRL(j))),dp)
     &            +real(zsin2winv*zr1*prop*(qgZqRR(2)*conjg(cRR(j))+qgZqLR(2)*conjg(cLR(j))),dp)
         endif
      elseif ((j < 0) .and. (k == 0)) then
         msq(j,k)=+real(cLL(-j)*conjg(cLL(-j))*qbgZqbLL(1),dp)
     &            +real(cRR(-j)*conjg(cRR(-j))*qbgZqbRR(1),dp)
     &            +real(cLR(-j)*conjg(cLR(-j))*qbgZqbLR(1),dp)
     &            +real(cRL(-j)*conjg(cRL(-j))*qbgZqbRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(zsin2winv*zl1*prop*(qbgZqbLL(2)*conjg(cLL(-j))+qbgZqbRL(2)*conjg(cRL(-j))),dp)
     &            +real(zsin2winv*zr1*prop*(qbgZqbRR(2)*conjg(cRR(-j))+qbgZqbLR(2)*conjg(cLR(-j))),dp)
         endif
      elseif ((j == 0) .and. (k > 0)) then
         msq(j,k)=+real(cLL(k)*conjg(cLL(k))*gqZqLL(1),dp)
     &            +real(cRR(k)*conjg(cRR(k))*gqZqRR(1),dp)
     &            +real(cLR(k)*conjg(cLR(k))*gqZqLR(1),dp)
     &            +real(cRL(k)*conjg(cRL(k))*gqZqRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(zsin2winv*zl1*prop*(gqZqLL(2)*conjg(cLL(k))+gqZqRL(2)*conjg(cRL(k))),dp)
     &            +real(zsin2winv*zr1*prop*(gqZqRR(2)*conjg(cRR(k))+gqZqLR(2)*conjg(cLR(k))),dp)
         endif
      elseif ((j == 0) .and. (k < 0)) then
         msq(j,k)=+real(cLL(-k)*conjg(cLL(-k))*gqbZqbLL(1),dp)
     &            +real(cRR(-k)*conjg(cRR(-k))*gqbZqbRR(1),dp)
     &            +real(cLR(-k)*conjg(cLR(-k))*gqbZqbLR(1),dp)
     &            +real(cRL(-k)*conjg(cRL(-k))*gqbZqbRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(zsin2winv*zl1*prop*(gqbZqbLL(2)*conjg(cLL(-k))+gqbZqbRL(2)*conjg(cRL(-k))),dp)
     &            +real(zsin2winv*zr1*prop*(gqbZqbRR(2)*conjg(cRR(-k))+gqbZqbLR(2)*conjg(cLR(-k))),dp)
         endif
      endif

   19 continue
      enddo
      enddo

      return
      end
