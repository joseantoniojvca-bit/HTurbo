!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_z2jet_g(p,msq)
!          use MCFMStorage, only: selectpdfs
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     March, 2001.                                                     *
c     Tested by JMC, 6/7/01 and found to agree with Madgraph           *
c***********************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  Z^0 + f(p5)+f(p6)+g(p7)
c                           |
c                            --> e^-(p3)+e^+(p4)

      include 'as3vjet_constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'as3vjet_masses.f'
      include 'as3vjet_qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com_7.f'
      include 'as3vjet_zcouple.f'
      include 'as3vjet_ewcharge.f'
      include 'flags.f'
      include 'lc.f'
      include 'first.f'
      include 'mpicommon.f'
      include 'zcouple_cms.f'
      integer:: j,k,nquark
      real(dp):: P(as3vj_mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: mmsq_gg(2,2),mmsq_qqb(2,2),mmsq_qbq(2,2),
     & mmsq_qg(2,2),mmsq_gq(2,2),mmsq_gqb(2,2),mmsq_qbg(2,2)
      real(dp):: fac
      real(dp):: msqi_qq(2),msqi_qbqb(2),
     &                 msqi_qqb(2),msqi_qbq(2),
     &                 msqi_qqbs(2),msqi_qbqs(2),
     &                 msqi_qg(2),msqi_qbg(2),
     &                 msqi_gqb(2),msqi_gq(2),ggtemp
      real(dp):: msqn_qq(2,2),msqn_qbqb(2,2),
     &                 msqn_qqb(2,2),msqn_qbq(2,2),
     &                 msqn_qqbs(2,2),msqn_qbqs(2,2),
     &                 msqn_qg(2,2),msqn_qbg(2,2),
     &                 msqn_gqb(2,2),msqn_gq(2,2)
      complex(dp):: prop
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      include 'cplx.h'
      logical :: filled(1:8)
      first=.false.
      if (first) then
        first=.false.
!$omp master
        if (rank == 0) then
        if ((Gflag) .or. (QandGflag)) then
          write(*,*) 'Using QQGG+G (REAL) matrix elements'
c          write(*,*) '[LC is     N   ]'
c          write(*,*) '[SLC is   1/N  ]'
c          write(*,*) '[SSLC is 1/N**3]'
        endif
        if ((Qflag) .or. (QandGflag)) then
          write(*,*) 'Using QQBQQB+G (REAL) matrix elements'
c          write(*,*) '[LC is   1 ]'
c          write(*,*) '[SLC is 1/N]'
        endif
        if     (colourchoice == 1) then
          write(*,*) 'Leading colour only in REAL'
        elseif (colourchoice == 2) then
          write(*,*) 'Sub-leading colour only in REAL'
        elseif (colourchoice == 3) then
          write(*,*) 'Sub-sub-leading colour only in REAL'
        elseif (colourchoice == 0) then
          write(*,*) 'Total of all colour structures in REAL'
        else
          write(*,*) 'Bad colourchoice'
          stop
        endif
        endif
!$omp end master
      endif

      msq(:,:)=0._dp

c---call spinor routine and load common block twopij
      call mcfm_spinoru(7,p,za,zb)

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

      if (Gflag) then
************************************************************************
*     Calculate contributions from the QQGGG matrix elements            *
************************************************************************
      filled = .false.
c***********************************************************************
      do j=-nf,nf
      do k=-nf,nf
c     Calculate contributions from the QQGGG matrix elements            *
c***********************************************************************
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) cycle

      if     ((j == 0) .and. (k == 0)) then
c-- matrix elements for gg -> qbq
           if (.not. filled(1)) then
             call xzqqggg(5,1,2,7,6,3,4,mmsq_gg)
             filled(1) = .true.
           endif

         ggtemp=0._dp
         do nquark=1,nf

           ggtemp=ggtemp
     &            +abs(Q(nquark)*q1+zL(nquark)*zl1*prop)**2*mmsq_gg(1,1)
     &            +abs(Q(nquark)*q1+zR(nquark)*zr1*prop)**2*mmsq_gg(2,2)
     &            +abs(Q(nquark)*q1+zL(nquark)*zr1*prop)**2*mmsq_gg(1,2)
     &            +abs(Q(nquark)*q1+zR(nquark)*zl1*prop)**2*mmsq_gg(2,1)
          enddo
          msq(j,k)=ggtemp
      elseif ((j > 0) .and. (k < 0)) then

          if (.not. filled(2)) then
            call xzqqggg(1,5,6,7,2,3,4,mmsq_qqb)
            filled(2) = .true.
          endif

          msq(j,k)=+abs(Q(j)*q1+zL(j)*zl1*prop)**2*mmsq_qqb(1,1)
     &             +abs(Q(j)*q1+zR(j)*zr1*prop)**2*mmsq_qqb(2,2)
     &             +abs(Q(j)*q1+zL(j)*zr1*prop)**2*mmsq_qqb(1,2)
     &             +abs(Q(j)*q1+zR(j)*zl1*prop)**2*mmsq_qqb(2,1)

c---Statistical factor
          msq(j,k)=aveqq/avegg*msq(j,k)/6._dp
      elseif ((j < 0) .and. (k > 0)) then
          if (.not. filled(3)) then
            call xzqqggg(2,5,6,7,1,3,4,mmsq_qbq)
            filled(3) = .true.
          endif
          msq(j,k)=+abs(Q(k)*q1+zL(k)*zl1*prop)**2*mmsq_qbq(1,1)
     &             +abs(Q(k)*q1+zR(k)*zr1*prop)**2*mmsq_qbq(2,2)
     &             +abs(Q(k)*q1+zL(k)*zr1*prop)**2*mmsq_qbq(1,2)
     &             +abs(Q(k)*q1+zR(k)*zl1*prop)**2*mmsq_qbq(2,1)
c---Statistical factor
          msq(j,k)=aveqq/avegg*msq(j,k)/6._dp
      elseif ((j > 0) .and. (k == 0)) then
          if (.not. filled(4)) then
            call xzqqggg(1,2,6,7,5,3,4,mmsq_qg)
            filled(4) = .true.
          endif
          msq(j,k)=+abs(Q(j)*q1+zL(j)*zl1*prop)**2*mmsq_qg(1,1)
     &             +abs(Q(j)*q1+zR(j)*zr1*prop)**2*mmsq_qg(2,2)
     &             +abs(Q(j)*q1+zL(j)*zr1*prop)**2*mmsq_qg(1,2)
     &             +abs(Q(j)*q1+zR(j)*zl1*prop)**2*mmsq_qg(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
      elseif ((j < 0) .and. (k == 0)) then
          if (.not. filled(5)) then
            call xzqqggg(5,2,6,7,1,3,4,mmsq_qbg)
            filled(5) = .true.
          endif
          msq(j,k)=+abs(Q(-j)*q1+zL(-j)*zl1*prop)**2*mmsq_qbg(1,1)
     &             +abs(Q(-j)*q1+zR(-j)*zr1*prop)**2*mmsq_qbg(2,2)
     &             +abs(Q(-j)*q1+zL(-j)*zr1*prop)**2*mmsq_qbg(1,2)
     &             +abs(Q(-j)*q1+zR(-j)*zl1*prop)**2*mmsq_qbg(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
      elseif ((j == 0) .and. (k > 0)) then
          if (.not. filled(6)) then
            call xzqqggg(2,1,6,7,5,3,4,mmsq_gq)
            filled(6) = .true.
          endif
          msq(j,k)=+abs(Q(k)*q1+zL(k)*zl1*prop)**2*mmsq_gq(1,1)
     &             +abs(Q(k)*q1+zR(k)*zr1*prop)**2*mmsq_gq(2,2)
     &             +abs(Q(k)*q1+zL(k)*zr1*prop)**2*mmsq_gq(1,2)
     &             +abs(Q(k)*q1+zR(k)*zl1*prop)**2*mmsq_gq(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
      elseif ((j == 0) .and. (k < 0)) then
          if (.not. filled(7)) then
            call xzqqggg(5,1,6,7,2,3,4,mmsq_gqb)
            filled(7) = .true.
          endif
          msq(j,k)=+abs(Q(-k)*q1+zL(-k)*zl1*prop)**2*mmsq_gqb(1,1)
     &             +abs(Q(-k)*q1+zR(-k)*zr1*prop)**2*mmsq_gqb(2,2)
     &             +abs(Q(-k)*q1+zL(-k)*zr1*prop)**2*mmsq_gqb(1,2)
     &             +abs(Q(-k)*q1+zR(-k)*zl1*prop)**2*mmsq_gqb(2,1)
c---Statistical factor
          msq(j,k)=half*aveqg/avegg*msq(j,k)
      endif

      enddo
      enddo
      endif

      if (Qflag) then
c--- note the factor of 4._dp*xw**2 relative to wbb
      fac=4._dp*gsq**3*abs(zesq)**2*8._dp
c--- extra factor of 2**3=8 to compensate for Ta normalization

      filled = .false.

      do j=-nf,nf
      do k=-nf,nf

       
      if ((j > 0) .and. (k < 0)) then
!-qqb
           if(.not. filled(1)) then
               !--- note: the following array ends up being overall 1<->2 symmetric
              call msq_ZqqQQg(5,1,2,6,7,4,3,msqn_qqbs,msqi_qqbs)
              call msq_ZqqQQg(2,1,5,6,7,4,3,msqn_qqb,msqi_qqb)
              filled(1) = .true.
           endif

           if (k==-j) then
            msq(j,k)=msq(j,k)+fac*aveqq*(msqi_qqbs(jj(j))
     &      +real(1+jj(j),dp)*msqn_qqb(jj(j),1)
     &      +real(3-jj(j),dp)*msqn_qqb(jj(j),2))
           else
            msq(j,k)=msq(j,k)+fac*aveqq*msqn_qqbs(jj(j),-kk(k))
           endif
      elseif ((j < 0) .and. (k > 0)) then
!-qbq
          if (.not. filled(2)) then
            !--- note: the following array ends up being overall 1<->2 symmetric
            call msq_ZqqQQg(1,6,5,2,7,4,3,msqn_qbqs,msqi_qbqs)
            call msq_ZqqQQg(1,2,5,6,7,4,3,msqn_qbq,msqi_qbq)
            filled(2) = .true.
          endif

          if (j ==-k) then
            msq(j,k)=msq(j,k)+fac*aveqq*(msqi_qbqs(kk(k))
     &      +real(1+kk(k),dp)*msqn_qbq(kk(k),1)
     &      +real(3-kk(k),dp)*msqn_qbq(kk(k),2))
          else
           msq(j,k)=msq(j,k)+fac*aveqq*msqn_qbqs(-jj(j),kk(k))
          endif
      elseif ((j > 0) .and. (k == 0)) then
!-qg
          if (.not. filled(3)) then
            call msq_ZqqQQg(7,1,5,6,2,4,3,msqn_qg,msqi_qg)
            filled(3) = .true.
          endif

          msq(j,k)=msq(j,k)
     &     +fac*aveqg*(half*msqi_qg(jj(j))
     &      +real(1+jj(j),dp)*msqn_qg(jj(j),1)
     &      +real(3-jj(j),dp)*msqn_qg(jj(j),2)
     & )
!-qbg
      elseif ((j < 0) .and. (k == 0)) then

          if (.not. filled(4)) then
            call msq_ZqqQQg(1,7,5,6,2,4,3,msqn_qbg,msqi_qbg)
            filled(4) = .true.
          endif

          msq(j,k)=msq(j,k)
     &     +fac*aveqg*(half*msqi_qbg(-jj(j))
     &      +real(1-jj(j),dp)*msqn_qbg(-jj(j),1)
     &      +real(3+jj(j),dp)*msqn_qbg(-jj(j),2)
     & )
      elseif ((j == 0) .and. (k > 0)) then
!-gq
          if (.not. filled(5)) then
            call msq_ZqqQQg(7,2,5,6,1,4,3,msqn_gq,msqi_gq)
            filled(5) = .true.
          endif

          msq(j,k)=msq(j,k)
     &      +fac*aveqg*(half*msqi_gq(kk(k))
     &      +real(1+kk(k),dp)*msqn_gq(kk(k),1)
     &      +real(3-kk(k),dp)*msqn_gq(kk(k),2)
     & )
      elseif ((j == 0) .and. (k < 0)) then
!-gqb
          if (.not. filled(6)) then
            call msq_ZqqQQg(2,7,5,6,1,4,3,msqn_gqb,msqi_gqb)
            filled(6) = .true.
          endif

          msq(j,k)=msq(j,k)
     &      +fac*aveqg*(half*msqi_gqb(-kk(k))
     &      +real(1-kk(k),dp)*msqn_gqb(-kk(k),1)
     &      +real(3+kk(k),dp)*msqn_gqb(-kk(k),2)
     & )
      elseif ((j > 0) .and. (k > 0)) then
!-qq
          if (.not. filled(7)) then
            call msq_ZqqQQg(5,1,6,2,7,4,3,msqn_qq,msqi_qq)
            filled(7) = .true.
          endif

          if (j==k) then
          msq(j,k)=msq(j,k)+half*fac*aveqq*msqi_qq(jj(j))
          else
          msq(j,k)=msq(j,k)+fac*aveqq*msqn_qq(jj(j),kk(k))
         endif
      elseif ((j < 0) .and. (k < 0)) then
!-qbqb
          if (.not. filled(8)) then
            call msq_ZqqQQg(1,5,2,6,7,4,3,msqn_qbqb,msqi_qbqb)
            filled(8) = .true.
          endif

          if (j==k) then
          msq(j,k)=msq(j,k)+half*fac*aveqq*msqi_qbqb(-jj(j))
          else
          msq(j,k)=msq(j,k)+fac*aveqq*msqn_qbqb(-jj(j),-kk(k))
          endif
      endif


      enddo
      enddo

      endif


      return
      end
