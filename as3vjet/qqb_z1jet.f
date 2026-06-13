!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine as3vj_qqb_z1jet(p,msq)
      implicit none
      include 'types.f'
c-----Authors: John Campbell, Keith Ellis
c-----June 2000 and December 2001
c----Matrix element for Z production
c----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))+g(p5)
c---
      include 'as3vjet_constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'as3vjet_masses.f'
      include 'as3vjet_qcdcouple.f'
      include 'as3vjet_zcouple.f'
      include 'as3vjet_ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'nflav.f'
      include 'zcouple_cms.f'
      integer:: j,k,hq,hl
      real(dp):: msq(-nf:nf,-nf:nf),p(as3vj_mxpart,4),fac
      complex(dp):: prop
      real(dp):: AqqbZg2(2,2),AqbqZg2(2,2),AqgZq2(2,2),
     &               AqbgZqb2(2,2),AgqbZqb2(2,2),AgqZq2(2,2)
      integer,parameter:: swap(2)=(/2,1/)
      include 'cplx.h'

      msq(:,:)=zip

      call mcfm_spinoru(5,p,za,zb)

c-----Protect from photon pole by cutting off at some value about 10 GeV
c      if (s(3,4) < 4._dp*mbsq) return

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)
      fac=four*V*abs(zesq)**2*gsq

c      qqbZg= +aveqq*s(3,4)**2*fac*z1jet(1,2,3,4,5)
c      gqbZqb=-aveqg*s(3,4)**2*fac*z1jet(5,2,3,4,1)
c      qgZq=  -aveqg*s(3,4)**2*fac*z1jet(1,5,3,4,2)
c      qbqZg= +aveqq*s(3,4)**2*fac*z1jet(2,1,3,4,5)
c      qbgZqb=-aveqg*s(3,4)**2*fac*z1jet(5,1,3,4,2)
c      gqZq=  -aveqg*s(3,4)**2*fac*z1jet(2,5,3,4,1)

      call as3vj_zgamps2(1,2,3,4,5,za,zb,AqqbZg2)
      call as3vj_zgamps2(1,5,3,4,2,za,zb,AqgZq2)
      call as3vj_zgamps2(2,5,3,4,1,za,zb,AgqZq2)
      do hq=1,2
      do hl=1,2
      AqbqZg2(hq,hl)=AqqbZg2(hq,swap(hl))
      AqbgZqb2(hq,hl)=AqgZq2(hq,swap(hl))
      AgqbZqb2(hq,hl)=AgqZq2(hq,swap(hl))
      enddo
      enddo

c      call as3vj_zgamps2(2,1,3,4,5,za,zb,AqbqZg2)
c      call as3vj_zgamps2(5,1,3,4,2,za,zb,AqbgZqb2)
c      call as3vj_zgamps2(5,2,3,4,1,za,zb,AgqbZqb2)

      do j=-nflav,nflav
      do k=-nflav,nflav

      if( (j  /=  0) .and. (k  /=  0) .and. (j  /=  -k)) goto 20

      if     ((j == 0) .and. (k == 0)) then
          msq(j,k)=zip
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=abs(Q(j)*q1+zL(j)*zl1*prop)**2*AqqbZg2(1,1)
     &            +abs(Q(j)*q1+zL(j)*zr1*prop)**2*AqqbZg2(1,2)
     &            +abs(Q(j)*q1+zR(j)*zl1*prop)**2*AqqbZg2(2,1)
     &            +abs(Q(j)*q1+zR(j)*zr1*prop)**2*AqqbZg2(2,2)
          msq(j,k)=msq(j,k)*aveqq*fac/s(3,4)**2
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=abs(Q(k)*q1+zL(k)*zl1*prop)**2*AqbqZg2(1,1)
     &            +abs(Q(k)*q1+zL(k)*zr1*prop)**2*AqbqZg2(1,2)
     &            +abs(Q(k)*q1+zR(k)*zl1*prop)**2*AqbqZg2(2,1)
     &            +abs(Q(k)*q1+zR(k)*zr1*prop)**2*AqbqZg2(2,2)
          msq(j,k)=msq(j,k)*aveqq*fac/s(3,4)**2
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=abs(Q(j)*q1+zL(j)*zl1*prop)**2*AqgZq2(1,1)
     &            +abs(Q(j)*q1+zL(j)*zr1*prop)**2*AqgZq2(1,2)
     &            +abs(Q(j)*q1+zR(j)*zl1*prop)**2*AqgZq2(2,1)
     &            +abs(Q(j)*q1+zR(j)*zr1*prop)**2*AqgZq2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=abs(Q(-j)*q1+zL(-j)*zl1*prop)**2*AqbgZqb2(1,1)
     &            +abs(Q(-j)*q1+zL(-j)*zr1*prop)**2*AqbgZqb2(1,2)
     &            +abs(Q(-j)*q1+zR(-j)*zl1*prop)**2*AqbgZqb2(2,1)
     &            +abs(Q(-j)*q1+zR(-j)*zr1*prop)**2*AqbgZqb2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=abs(Q(k)*q1+zL(k)*zl1*prop)**2*AgqZq2(1,1)
     &            +abs(Q(k)*q1+zL(k)*zr1*prop)**2*AgqZq2(1,2)
     &            +abs(Q(k)*q1+zR(k)*zl1*prop)**2*AgqZq2(2,1)
     &            +abs(Q(k)*q1+zR(k)*zr1*prop)**2*AgqZq2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=abs(Q(-k)*q1+zL(-k)*zl1*prop)**2*AgqbZqb2(1,1)
     &            +abs(Q(-k)*q1+zL(-k)*zr1*prop)**2*AgqbZqb2(1,2)
     &            +abs(Q(-k)*q1+zR(-k)*zl1*prop)**2*AgqbZqb2(2,1)
     &            +abs(Q(-k)*q1+zR(-k)*zr1*prop)**2*AgqbZqb2(2,2)
          msq(j,k)=msq(j,k)*aveqg*fac/s(3,4)**2
      endif

 20   continue
      enddo
      enddo
      return
      end


      subroutine as3vj_zgamps2(j1,j2,j3,j4,j5,za,zb,amps2)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: amps(2,2,2)
      real(dp):: amps2(2,2)
      integer:: hq,hl,j1,j2,j3,j4,j5
c-- amplitude helicities are amps(quark,lepton,gluon)

      amps(1,1,1)=za(j2,j3)/za(j1,j5)/za(j2,j5)
     &             *(za(j2,j1)*zb(j4,j1)+za(j2,j5)*zb(j4,j5))

      amps(1,1,2)=zb(j4,j1)/zb(j1,j5)/zb(j2,j5)
     &             *(za(j2,j3)*zb(j2,j1)+za(j3,j5)*zb(j1,j5))

      amps(1,2,1)=za(j2,j4)/za(j1,j5)/za(j2,j5)
     &             *(za(j2,j1)*zb(j3,j1)+za(j2,j5)*zb(j3,j5))

      amps(1,2,2)=zb(j3,j1)/zb(j1,j5)/zb(j2,j5)
     &             *(za(j2,j4)*zb(j2,j1)+za(j4,j5)*zb(j1,j5))

c obtained from above by  (za <-> zb)
      amps(2,2,2)=zb(j2,j3)/zb(j1,j5)/zb(j2,j5)
     &             *(zb(j2,j1)*za(j4,j1)+zb(j2,j5)*za(j4,j5))

      amps(2,2,1)=za(j4,j1)/za(j1,j5)/za(j2,j5)
     &             *(zb(j2,j3)*za(j2,j1)+zb(j3,j5)*za(j1,j5))

      amps(2,1,2)=zb(j2,j4)/zb(j1,j5)/zb(j2,j5)
     &             *(zb(j2,j1)*za(j3,j1)+zb(j2,j5)*za(j3,j5))

      amps(2,1,1)=za(j3,j1)/za(j1,j5)/za(j2,j5)
     &             *(zb(j2,j4)*za(j2,j1)+zb(j4,j5)*za(j1,j5))

c      write(6,*) j1,j2,j3,j4,j5
c      do hl=1,2
c      do hg=1,2
c        write(6,*) amps(2,hl,hg)/(-conjg(amps(1,3-hl,3-hg)))
c      enddo
c      enddo
c      pause

      do hq=1,2
      do hl=1,2
        amps2(hq,hl)=abs(amps(hq,hl,1))**2+abs(amps(hq,hl,2))**2
      enddo
      enddo


      return
      end

