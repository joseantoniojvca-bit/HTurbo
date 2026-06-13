!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

module nnlo_z1jet
    use types
    use SCET
    use SCET_Jet
    implicit none

    public :: passed_taucut_1jet
    public :: lumxmsq_z1jet
    public :: genps_z1jet_r ! for single real emission (snlo above)
    public :: genps_z1jet_rr ! for double real emission (nnlo above)
    public :: qqb_z1jet_vbis

    private
      
    private :: Zampqqbgsq
    private :: Zampgggsq
    private :: ampgggppp, ampgggpmm
    
    logical, parameter :: comparedecay = .false.

    contains

    function genps_z1jet_r(r,p,wt)
        use types
        use Multichannel
        implicit none

        include 'mxpart.f'
        include 'kpart.f'
        include 'taucut.f' ! usescet
        include 'npart.f'
        include 'vegas_common.f' ! ndim
        include 'mxdim.f'
        include 'ipsgen.f'

        logical :: genps_z1jet_r
        real(dp), intent(in) :: r(mxdim)
        real(dp), intent(out) :: p(as3vj_mxpart,4)
        real(dp), intent(out) :: wt

        real(dp) :: wtdip, plo(as3vj_mxpart,4), pone(as3vj_mxpart,4), pswt
        real(dp), allocatable :: dipconfig(:,:)

        ! select between beam and jet regions here
        ! both have to be added eventually
        integer, parameter :: taups = 0

        genps_z1jet_r = .false.

#define PSCHOICE 0

#if (PSCHOICE == 0)
        ! default gen_realps
        npart=4
        call gen_njets(r,2,p,pswt,*999) 
        wt = pswt

#elif (PSCHOICE == 1)
        ! trying one multichan emission

        npart = 3
        call gen_njets(r,1,pone,pswt,*999)
        wt = pswt

        ! adjust emitters/spectators here
        allocate(dipconfig(6,3))
        dipconfig(1,:) = [5,6,1]
        dipconfig(2,:) = [5,6,2]
        dipconfig(3,:) = [1,6,5]
        dipconfig(4,:) = [2,6,5]
        dipconfig(5,:) = [1,6,2]
        dipconfig(6,:) = [2,6,1]


        npart = npart + 1
        if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                    pone,p,pswt, dipconfig_in=dipconfig)) then
            genps_z1jet_r = .false.
            return
        endif
        wt = wt*pswt

! this two-staged approach works also well
! not better than one-staged though, but about the same
#elif (PSCHOICE == 2)

        npart = 2
        call gen2(r,plo,pswt,*999)
        wt = pswt

        !!! first emission
        allocate(dipconfig(2,3))
        dipconfig(1,:) = [1,5,2]
        dipconfig(2,:) = [2,5,1]

        npart = npart + 1
        if (.not. multichan(r(ndim-5),r(ndim-4),r(ndim-3),r(ndim+1), &
                    plo,pone,pswt, dipconfig_in=dipconfig)) then
            genps_z1jet_r = .false.
            return
        endif
        wt = wt*pswt
        deallocate(dipconfig)

        !!! second emission
        allocate(dipconfig(6,3))
        dipconfig(1,:) = [5,6,1]
        dipconfig(2,:) = [5,6,2]
        dipconfig(3,:) = [1,6,5]
        dipconfig(4,:) = [2,6,5]
        dipconfig(5,:) = [1,6,2]
        dipconfig(6,:) = [2,6,1]

        npart = npart + 1
        if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                    pone,p,pswt, dipconfig_in=dipconfig)) then
            genps_z1jet_r = .false.
            return
        endif
        wt = wt*pswt

!! first the two-staged approach, then select by
!! clustering with beam or jet

!! at first sight this doesn't seem to be any better
!! but I could imagine that for a very high precision
!! the adaptability is improved, resulting in a overall net benefit
#elif (PSCHOICE == 3)

        ! it turns out that the beam clustering region
        ! works *much* better with plain gen_njets than
        ! any dipole emission
        if (taups == 1) then
            npart=4
            call gen_njets(r,2,p,pswt,*999) 
            wt = pswt
        else
       ! on the other hand the jet clustering region
       ! really needs the dipole emissions
            npart = 2
            call gen2(r,plo,pswt,*999)
            wt = pswt

            !!! first emission
            allocate(dipconfig(2,3))
            dipconfig(1,:) = [1,5,2]
            dipconfig(2,:) = [2,5,1]

            npart = npart + 1
            if (.not. multichan(r(ndim-5),r(ndim-4),r(ndim-3),r(ndim+1), &
                       plo,pone,pswt, dipconfig_in=dipconfig)) then
               genps_z1jet_r = .false.
               return
            endif
            wt = wt*pswt
            deallocate(dipconfig)

            !!! second emission
            allocate(dipconfig(6,3))
            dipconfig(1,:) = [5,6,1]
            dipconfig(2,:) = [5,6,2]
            dipconfig(3,:) = [1,6,5]
            dipconfig(4,:) = [2,6,5]
            dipconfig(5,:) = [1,6,2]
            dipconfig(6,:) = [2,6,1]

            npart = npart + 1
            if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                       pone,p,pswt, dipconfig_in=dipconfig)) then
               genps_z1jet_r = .false.
               return
            endif
            wt = wt*pswt
        endif

#define E(i) p(i,4)
#define pz(i) (p(i,3))
        tauphasesnlo: block
            real(dp) :: p56(4)
            real(dp) :: taubeam, taujet
            real(dp) :: wtbeam, wtjet

            p56(:) = p(5,:) + p(6,:)

            ! extra emission clustered with a beam
            taubeam = min(E(5)-abs(pz(5)), E(6)-abs(pz(6)))
            ! extra emission clustered with the jet
            taujet =   E(5)+E(6) - sqrt((p56(1)**2 + p56(2)**2 + p56(3)**2))

            wtbeam = 1._dp/taubeam
            wtjet = 1._dp/taujet

            if (taups == 1) then
                wt = wt * wtbeam/(wtbeam + wtjet)
            elseif (taups == 2) then
                wt = wt * wtjet/(wtbeam + wtjet)
            else
                continue
                !error stop "wrong taups in phase space generation"
            endif

        end block tauphasesnlo
#undef E
#undef pz
    
#else
        write (*,*) "INVALID PSCHOICE"
        stop
#endif

#undef PSCHOICE

        genps_z1jet_r = .true.

        return
 999    continue

        genps_z1jet_r = .false.
        return

    end function genps_z1jet_r


    function genps_z1jet_rr(r,p,wt)
        use types
        use Multichannel
        implicit none

        include 'mxpart.f'
        include 'kpart.f'
        include 'taucut.f' ! usescet
        include 'npart.f'
        include 'vegas_common.f' ! ndim
        include 'mxdim.f'
        include 'ipsgen.f'

        logical :: genps_z1jet_rr
        real(dp), intent(in) :: r(mxdim)
        real(dp), intent(out) :: p(as3vj_mxpart,4)
        real(dp), intent(out) :: wt

        real(dp) :: wtdip, plo(as3vj_mxpart,4), pone(as3vj_mxpart,4), pswt
        real(dp) :: ptwo(as3vj_mxpart,4)
        real(dp), allocatable :: dipconfig(:,:)

        integer :: taups

        taups = ipsgen

        genps_z1jet_rr = .false.

#define PSCHOICE 0

#if (PSCHOICE == 0)
        ! default gen_realps
        npart=5
        call gen_njets(r,3,p,pswt,*999) 
        wt = pswt

#elif (PSCHOICE == 1)
        
        !!! starting with zero emissions does not seem to work
        !npart = 2
        !call gen2(r,plo,pswt,*999)
        !wt = pswt

        !!!! first emission
        !allocate(dipconfig(2,3))
        !dipconfig(1,:) = [1,5,2]
        !dipconfig(2,:) = [2,5,1]

        !npart = npart + 1
        !if (.not. multichan(r(ndim-8),r(ndim-7),r(ndim-6),r(ndim+1), &
        !            plo,pone,pswt, dipconfig_in=dipconfig)) then
        !    genps_z1jet_rr = .false.
        !    return
        !endif
        !wt = wt*pswt
        !deallocate(dipconfig)
        !!! end start with zero emissions

        npart = 3
        call gen_njets(r,1,pone,pswt,*999)
        wt = pswt

        !!! second emission
        allocate(dipconfig(6,3))
        dipconfig(1,:) = [5,6,1]
        dipconfig(2,:) = [5,6,2]
        dipconfig(3,:) = [1,6,5]
        dipconfig(4,:) = [2,6,5]
        dipconfig(5,:) = [1,6,2]
        dipconfig(6,:) = [2,6,1]

        npart = npart + 1
        if (.not. multichan(r(ndim-5),r(ndim-4),r(ndim-3),r(ndim+2), &
                    pone,ptwo,pswt, dipconfig_in=dipconfig)) then
            genps_z1jet_rr = .false.
            return
        endif
        wt = wt*pswt
        deallocate(dipconfig)


        !!! third emission
        allocate(dipconfig(12,3))
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        dipconfig(3,:)= (/ 1,7,5 /)
        dipconfig(4,:)= (/ 2,7,6 /)
        dipconfig(5,:)= (/ 1,7,6 /)
        dipconfig(6,:)= (/ 2,7,5 /)
        dipconfig(7,:)= (/ 5,7,1 /)
        dipconfig(8,:)= (/ 6,7,2 /)
        dipconfig(9,:)= (/ 6,7,1 /)
        dipconfig(10,:)=(/ 5,7,2 /)
        dipconfig(11,:)=(/ 5,7,6 /)
        dipconfig(12,:)=(/ 6,7,5 /)

        npart = npart + 1
        if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+3), &
                    ptwo,p,pswt, dipconfig_in=dipconfig)) then
            genps_z1jet_rr = .false.
            return
        endif
        wt = wt*pswt

! let's try some crazy three-staged emission
#elif (PSCHOICE == 2)

        npart = 2
        call gen2(r,plo,pswt,*999)
        wt = pswt

        !!! first emission
        allocate(dipconfig(2,3))
        dipconfig(1,:) = [1,5,2]
        dipconfig(2,:) = [2,5,1]

        npart = npart + 1
        if (.not. multichan(r(ndim-8),r(ndim-7),r(ndim-6),r(ndim+1), &
                    plo,pone,pswt, dipconfig_in=dipconfig)) then
            genps_z1jet_rr = .false.
            return
        endif
        wt = wt*pswt
        deallocate(dipconfig)

        !!! second emission
        allocate(dipconfig(6,3))
        dipconfig(1,:) = [5,6,1]
        dipconfig(2,:) = [5,6,2]
        dipconfig(3,:) = [1,6,5]
        dipconfig(4,:) = [2,6,5]
        dipconfig(5,:) = [1,6,2]
        dipconfig(6,:) = [2,6,1]

        npart = npart + 1
        if (.not. multichan(r(ndim-5),r(ndim-4),r(ndim-3),r(ndim+2), &
                    pone,ptwo,pswt, dipconfig_in=dipconfig)) then
            genps_z1jet_rr = .false.
            return
        endif
        wt = wt*pswt
        deallocate(dipconfig)


        !!! third emission
        allocate(dipconfig(12,3))
        dipconfig(1,:)= (/ 1,7,2 /)
        dipconfig(2,:)= (/ 2,7,1 /)
        dipconfig(3,:)= (/ 1,7,5 /)
        dipconfig(4,:)= (/ 2,7,6 /)
        dipconfig(5,:)= (/ 1,7,6 /)
        dipconfig(6,:)= (/ 2,7,5 /)
        dipconfig(7,:)= (/ 5,7,1 /)
        dipconfig(8,:)= (/ 6,7,2 /)
        dipconfig(9,:)= (/ 6,7,1 /)
        dipconfig(10,:)=(/ 5,7,2 /)
        dipconfig(11,:)=(/ 5,7,6 /)
        dipconfig(12,:)=(/ 6,7,5 /)

        npart = npart + 1
        if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+3), &
                    ptwo,p,pswt, dipconfig_in=dipconfig)) then
            genps_z1jet_rr = .false.
            return
        endif
        wt = wt*pswt

#elif (PSCHOICE == 3)

!       npart = 2
!       call gen2(r,plo,pswt,*999)
!       wt = pswt

!       !!! first emission
!       allocate(dipconfig(2,3))
!       dipconfig(1,:) = [1,5,2]
!       dipconfig(2,:) = [2,5,1]

!       npart = npart + 1
!       if (.not. multichan(r(ndim-8),r(ndim-7),r(ndim-6),r(ndim+1), &
!                   plo,pone,pswt, dipconfig_in=dipconfig)) then
!           genps_z1jet_rr = .false.
!           return
!       endif
!       wt = wt*pswt
!       deallocate(dipconfig)

!       npart = 3
!       call gen_njets(r,1,pone,pswt,*999)
!       wt = pswt

!       !!! second emission
!       allocate(dipconfig(6,3))
!       dipconfig(1,:) = [5,6,1]
!       dipconfig(2,:) = [5,6,2]
!       dipconfig(3,:) = [1,6,5]
!       dipconfig(4,:) = [2,6,5]
!       dipconfig(5,:) = [1,6,2]
!       dipconfig(6,:) = [2,6,1]

!       npart = npart + 1
!       if (.not. multichan(r(ndim-5),r(ndim-4),r(ndim-3),r(ndim+2), &
!                   pone,ptwo,pswt, dipconfig_in=dipconfig)) then
!           genps_z1jet_rr = .false.
!           return
!       endif
!       wt = wt*pswt
!       deallocate(dipconfig)

        ! the following is the result of some trial and error
        ! using various levels of multichan radiation and choices of pT integral parametrization
        if (taups == 2) then
            npart = 4
            call gen_njets(r,2,ptwo,pswt,*999)
            wt = pswt

            !!! third emission
            allocate(dipconfig(12,3))
            dipconfig(1,:)= (/ 1,7,2 /)
            dipconfig(2,:)= (/ 2,7,1 /)
            dipconfig(3,:)= (/ 1,7,5 /)
            dipconfig(4,:)= (/ 2,7,6 /)
            dipconfig(5,:)= (/ 1,7,6 /)
            dipconfig(6,:)= (/ 2,7,5 /)
            dipconfig(7,:)= (/ 5,7,1 /)
            dipconfig(8,:)= (/ 6,7,2 /)
            dipconfig(9,:)= (/ 6,7,1 /)
            dipconfig(10,:)=(/ 5,7,2 /)
            dipconfig(11,:)=(/ 5,7,6 /)
            dipconfig(12,:)=(/ 6,7,5 /)

            npart = npart + 1
            if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+3), &
                        ptwo,p,pswt, dipconfig_in=dipconfig)) then
                genps_z1jet_rr = .false.
                return
            endif
            wt = wt*pswt
        else
            npart=5
            call gen_njets(r,3,p,pswt,*999) 
            wt = pswt
        endif

#define E(i) p(i,4)
#define pz(i) (p(i,3))
        tauphasennlo: block
            real(dp) :: p56(4), p57(4), p67(4), p567(4)
            real(dp) :: tauallbeam, taubeamjet, taualljet
            real(dp) :: wtallbeam, wtbeamjet, wtalljet

            p56(:) = p(5,:) + p(6,:)
            p57(:) = p(5,:) + p(7,:)
            p67(:) = p(6,:) + p(7,:)
            p567(:) = p(5,:) + p(6,:) + p(7,:)

            tauallbeam = min(E(5)+E(6)-abs(pz(5))-abs(pz(6)), &
                             E(5)+E(7)-abs(pz(5))-abs(pz(7)), &
                             E(6)+E(7)-abs(pz(6))-abs(pz(7)))
            taubeamjet = min(E(5)+E(6)+E(7) - sqrt(p56(1)**2+p56(2)**2+p56(3)**2) - abs(pz(7)), &
                             E(5)+E(6)+E(7) - sqrt(p57(1)**2+p57(2)**2+p57(3)**2) - abs(pz(6)), &
                             E(5)+E(6)+E(7) - sqrt(p67(1)**2+p67(2)**2+p67(3)**2) - abs(pz(5)))
            taualljet = E(5)+E(6)+E(7) - sqrt(p567(1)**2+p567(2)**2+p567(3)**2)

            wtallbeam = 1._dp/tauallbeam
            wtbeamjet = 1._dp/taubeamjet
            wtalljet = 1._dp/taualljet

            if (taups == 1) then
                wt = wt * wtallbeam/(wtalljet+wtallbeam+wtbeamjet)
            elseif (taups == 2) then
                wt = wt * wtbeamjet/(wtalljet+wtallbeam+wtbeamjet)
            elseif (taups == 3) then
                wt = wt * wtalljet/(wtalljet+wtallbeam+wtbeamjet)
            else
                continue
                !error stop "wrong taups in phase space generation"
            endif

        end block tauphasennlo
#undef E
#undef pz
    
#else
        write (*,*) "INVALID PSCHOICE"
        stop
#endif

        genps_z1jet_rr = .true.

        return
 999    continue

        genps_z1jet_rr = .false.
        return

    end function genps_z1jet_rr

      subroutine qqb_z1jet_vbis(p,msq,order)
          use types
      implicit none
      include 'as3vjet_constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'as3vjet_masses.f'
      include 'as3vjet_qcdcouple.f'
      include 'as3vjet_ewcouple.f'
      include 'as3vjet_zcouple.f'
      include 'as3vjet_ewcharge.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'nflav.f'

      real(dp), intent(in) :: p(as3vj_mxpart,4)
      real(dp), intent(out) :: msq(-nf:nf,-nf:nf)
      integer, intent(in) :: order

      integer:: j,k
      real(dp):: fac,facgg
      complex(dp):: prop
      real(dp):: qqbZg(2,2),qbqZg(2,2),qgZq(2,2), &
       qbgZqb(2,2),gqbZqb(2,2),gqZq(2,2),ggZg(2)
      real(dp):: &
       Hqqb1(2,2),Hqg1(2,2),Hgqb1(2,2),Hqbq1(2,2),Hqbg1(2,2),Hgq1(2,2), &
       Hqqb2(2,2),Hqg2(2,2),Hgqb2(2,2),Hqbq2(2,2),Hqbg2(2,2),Hgq2(2,2), &
       Hqqb2NFV(2,2),Hqg2NFV(2,2),Hgqb2NFV(2,2), &
       Hqbq2NFV(2,2),Hqbg2NFV(2,2),Hgq2NFV(2,2),Qsumg,QsumZ

      msq(:,:)=0._dp
      
      call mcfm_spinoru(5,p,za,zb)
!      print *,'za',za
!      print *,'zb',zb
      prop=s(3,4)/cmplx((s(3,4)-zmass**2),zmass*zwidth,dp)
      fac=four*V*esq**2*gsq
      facgg=(40._dp/3._dp)*esq**2*gsq

      if (comparedecay) then
!--- remove photon contribution for debugging with/without decay
        q1=zip
      endif

!--- compute matrix elements
!--- for Zggg, note that 4,3 are flipped to reflect different notation in amp
      call Zampgggsq(order,1,2,5,4,3,za,zb,ggZg)

      call Zampqqbgsq(order,2,1,5,3,4,zip,za,zb,qqbZg,Hqqb1,Hqqb2)
      call Zampqqbgsq(order,5,1,2,3,4,zip,za,zb,qgZq,Hqg1,Hqg2)
      call Zampqqbgsq(order,2,5,1,3,4,zip,za,zb,gqbZqb,Hgqb1,Hgqb2)

      call Zampqqbgsq(order,1,2,5,3,4,zip,za,zb,qbqZg,Hqbq1,Hqbq2)
      call Zampqqbgsq(order,1,5,2,3,4,zip,za,zb,qbgZqb,Hqbg1,Hqbg2)
      call Zampqqbgsq(order,5,2,1,3,4,zip,za,zb,gqZq,Hgq1,Hgq2)

!--- extract term proportional to NFV piece of 2-loop calculation
      if ((order == 2) .and. (comparedecay .eqv. .false.)) then
      call Zampqqbgsq(order,2,1,5,3,4,one,za,zb,qqbZg,Hqqb1,Hqqb2NFV)
      call Zampqqbgsq(order,5,1,2,3,4,one,za,zb,qgZq,Hqg1,Hqg2NFV)
      call Zampqqbgsq(order,2,5,1,3,4,one,za,zb,gqbZqb,Hgqb1,Hgqb2NFV)

      call Zampqqbgsq(order,1,2,5,3,4,one,za,zb,qbqZg,Hqbq1,Hqbq2NFV)
      call Zampqqbgsq(order,1,5,2,3,4,one,za,zb,qbgZqb,Hqbg1,Hqbg2NFV)
      call Zampqqbgsq(order,5,2,1,3,4,one,za,zb,gqZq,Hgq1,Hgq2NFV)
!      print *,'qgZq',qgZq
!      print *,'Hqg1',Hqg1
!      print *,'Hqg2NFV',Hqg2NFV
!      print *,'gqZq',gqZq
!      print *,'Hgq1',Hgq1
!      print *,'Hgq2NFV',Hgq2NFV
      
!      print *,'Hqqb2NFV',Hqqb2NFV(:,:)
!      print *,'Hqg2NFV ',Hqg2NFV(:,:) 
!      print *,'Hgqb2NFV',Hgqb2NFV(:,:)
!      print *,'Hqbq2NFV',Hqbq2NFV(:,:)
!      print *,'Hqbg2NFV',Hqbg2NFV(:,:)
!      print *,'Hgq2NFV ',Hgq2NFV(:,:) 
      
      Hqqb2NFV(:,:)=Hqqb2NFV(:,:) - Hqqb2(:,:)
      Hqg2NFV(:,:) =Hqg2NFV(:,:)  - Hqg2(:,:) 
      Hgqb2NFV(:,:)=Hgqb2NFV(:,:) - Hgqb2(:,:)
      Hqbq2NFV(:,:)=Hqbq2NFV(:,:) - Hqbq2(:,:)
      Hqbg2NFV(:,:)=Hqbg2NFV(:,:) - Hqbg2(:,:)
      Hgq2NFV(:,:) =Hgq2NFV(:,:)  - Hgq2(:,:) 
       
      else
       
      Hqqb2NFV(:,:)=zip
      Hqg2NFV(:,:)=zip
      Hgqb2NFV(:,:)=zip
      Hqbq2NFV(:,:)=zip
      Hqbg2NFV(:,:)=zip
      Hgq2NFV(:,:)=zip 
       
      endif
!--- apply overall factor
      qqbZg= aveqq*fac*qqbZg
      qgZq  =aveqg*fac*qgZq
      gqbZqb=aveqg*fac*gqbZqb
      qbqZg =aveqq*fac*qbqZg
      qbgZqb=aveqg*fac*qbgZqb
      gqZq  =aveqg*fac*gqZq

      ggZg(:)=avegg*facgg*(ason2pi)**2*ggZg(:)
      
!--- pick out higher-order coefficients
      if (order == 1) then
        qqbZg(:,:) =ason2pi*Hqqb1(:,:)*qqbZg(:,:)
        qgZq(:,:)  =ason2pi*Hqg1(:,:)*qgZq(:,:)
        gqbZqb(:,:)=ason2pi*Hgqb1(:,:)*gqbZqb(:,:)
        qbqZg(:,:) =ason2pi*Hqbq1(:,:)*qbqZg(:,:)
        qbgZqb(:,:)=ason2pi*Hqbg1(:,:)*qbgZqb(:,:)
        gqZq(:,:)  =ason2pi*Hgq1(:,:)*gqZq(:,:)
      endif
      if (order == 2) then
        Hqqb2NFV(:,:)=ason2pi**2*Hqqb2NFV(:,:)*qqbZg(:,:)
        Hqg2NFV(:,:) =ason2pi**2*Hqg2NFV(:,:)*qgZq(:,:)
        Hgqb2NFV(:,:)=ason2pi**2*Hgqb2NFV(:,:)*gqbZqb(:,:)
        Hqbq2NFV(:,:)=ason2pi**2*Hqbq2NFV(:,:)*qbqZg(:,:)
        Hqbg2NFV(:,:)=ason2pi**2*Hqbg2NFV(:,:)*qbgZqb(:,:)
        Hgq2NFV(:,:) =ason2pi**2*Hgq2NFV(:,:)*gqZq(:,:)
       
        qqbZg(:,:) =ason2pi**2*Hqqb2(:,:)*qqbZg(:,:)
        qgZq(:,:)  =ason2pi**2*Hqg2(:,:)*qgZq(:,:)
        gqbZqb(:,:)=ason2pi**2*Hgqb2(:,:)*gqbZqb(:,:)
        qbqZg(:,:) =ason2pi**2*Hqbq2(:,:)*qbqZg(:,:)
        qbgZqb(:,:)=ason2pi**2*Hqbg2(:,:)*qbgZqb(:,:)
        gqZq(:,:)  =ason2pi**2*Hgq2(:,:)*gqZq(:,:)

      endif

! couplings for gg->Zg according to Eqs.(2.19) and (2.20) of 1302.2630;
! note that this is also the coupling for the NFV pieces as well,
!  c.f. Eqs. (3.8) and (3.9) of 1112.1531
!--- five light quark flavors are included in loop
      Qsumg=2._dp*Qu+3._dp*Qd
      QsumZ=half*(2._dp*(L(2)+R(2))+3._dp*(L(1)+R(1)))

      do j=-nflav,nflav
      do k=-nflav,nflav

      if( (j .ne. 0) .and. (k .ne. 0) .and. (j .ne. -k)) goto 20

      if     ((j == 0) .and. (k == 0)) then
          msq(j,k)=abs(Qsumg*q1+QsumZ*l1*prop)**2*ggZg(1) &
                  +abs(Qsumg*q1+QsumZ*r1*prop)**2*ggZg(2)
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=abs(Q(j)*q1+L(j)*l1*prop)**2*qqbZg(1,1) &
                  +abs(Q(j)*q1+L(j)*r1*prop)**2*qqbZg(1,2) &
                  +abs(Q(j)*q1+R(j)*l1*prop)**2*qqbZg(2,1) &
                  +abs(Q(j)*q1+R(j)*r1*prop)**2*qqbZg(2,2)
          
          msq(j,k)=msq(j,k) &
           +real(conjg(Q(j)*q1+L(j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqqb2NFV(1,1) &
           +real(conjg(Q(j)*q1+L(j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqqb2NFV(1,2) &
           +real(conjg(Q(j)*q1+R(j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqqb2NFV(2,1) &
           +real(conjg(Q(j)*q1+R(j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqqb2NFV(2,2)
     
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=abs(Q(k)*q1+L(k)*l1*prop)**2*qbqZg(1,1) &
                  +abs(Q(k)*q1+L(k)*r1*prop)**2*qbqZg(1,2) &
                  +abs(Q(k)*q1+R(k)*l1*prop)**2*qbqZg(2,1) &
                  +abs(Q(k)*q1+R(k)*r1*prop)**2*qbqZg(2,2)
          
          msq(j,k)=msq(j,k) &
           +real(conjg(Q(k)*q1+L(k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqbq2NFV(1,1) &
           +real(conjg(Q(k)*q1+L(k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqbq2NFV(1,2) &
           +real(conjg(Q(k)*q1+R(k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqbq2NFV(2,1) &
           +real(conjg(Q(k)*q1+R(k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqbq2NFV(2,2)
     
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=abs(Q(j)*q1+L(j)*l1*prop)**2*qgZq(1,1) &
                  +abs(Q(j)*q1+L(j)*r1*prop)**2*qgZq(1,2) &
                  +abs(Q(j)*q1+R(j)*l1*prop)**2*qgZq(2,1) &
                  +abs(Q(j)*q1+R(j)*r1*prop)**2*qgZq(2,2)
          
          msq(j,k)=msq(j,k) &
           +real(conjg(Q(j)*q1+L(j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqg2NFV(1,1) &
           +real(conjg(Q(j)*q1+L(j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqg2NFV(1,2) &
           +real(conjg(Q(j)*q1+R(j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqg2NFV(2,1) &
           +real(conjg(Q(j)*q1+R(j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqg2NFV(2,2)
!          print *,j,k,qgZq,Hqg2NFV
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=abs(Q(-j)*q1+L(-j)*l1*prop)**2*qbgZqb(1,1) &
                  +abs(Q(-j)*q1+L(-j)*r1*prop)**2*qbgZqb(1,2) &
                  +abs(Q(-j)*q1+R(-j)*l1*prop)**2*qbgZqb(2,1) &
                  +abs(Q(-j)*q1+R(-j)*r1*prop)**2*qbgZqb(2,2)
          
          msq(j,k)=msq(j,k) &
           +real(conjg(Q(-j)*q1+L(-j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqbg2NFV(1,1) &
           +real(conjg(Q(-j)*q1+L(-j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqbg2NFV(1,2) &
           +real(conjg(Q(-j)*q1+R(-j)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hqbg2NFV(2,1) &
           +real(conjg(Q(-j)*q1+R(-j)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hqbg2NFV(2,2)
!          print *,j,k,qbgZqb,Hqbg2NFV
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=abs(Q(k)*q1+L(k)*l1*prop)**2*gqZq(1,1) &
                  +abs(Q(k)*q1+L(k)*r1*prop)**2*gqZq(1,2) &
                  +abs(Q(k)*q1+R(k)*l1*prop)**2*gqZq(2,1) &
                  +abs(Q(k)*q1+R(k)*r1*prop)**2*gqZq(2,2)
          
          msq(j,k)=msq(j,k) &
           +real(conjg(Q(k)*q1+L(k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hgq2NFV(1,1) &
           +real(conjg(Q(k)*q1+L(k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hgq2NFV(1,2) &
           +real(conjg(Q(k)*q1+R(k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hgq2NFV(2,1) &
           +real(conjg(Q(k)*q1+R(k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hgq2NFV(2,2)
!          print *,j,k,gqZq,Hgq2NFV
          msq(j,k)=msq(j,k)/s(3,4)**2
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=abs(Q(-k)*q1+L(-k)*l1*prop)**2*gqbZqb(1,1) &
                  +abs(Q(-k)*q1+L(-k)*r1*prop)**2*gqbZqb(1,2) &
                  +abs(Q(-k)*q1+R(-k)*l1*prop)**2*gqbZqb(2,1) &
                  +abs(Q(-k)*q1+R(-k)*r1*prop)**2*gqbZqb(2,2)
          
          msq(j,k)=msq(j,k) &
           +real(conjg(Q(-k)*q1+L(-k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hgqb2NFV(1,1) &
           +real(conjg(Q(-k)*q1+L(-k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hgqb2NFV(1,2) &
           +real(conjg(Q(-k)*q1+R(-k)*l1*prop)*(Qsumg*q1+QsumZ*l1*prop),dp)*Hgqb2NFV(2,1) &
           +real(conjg(Q(-k)*q1+R(-k)*r1*prop)*(Qsumg*q1+QsumZ*r1*prop),dp)*Hgqb2NFV(2,2)
!          print *,j,k,gqbZqb,Hgqb2NFV
          msq(j,k)=msq(j,k)/s(3,4)**2
      endif

20    continue
      enddo
      enddo
      return
      end
 

      subroutine lumxmsq_z1jet(p,xx,z1,z2,QB,order,xmsq,central)
          use types
!          use LHAPDF
          use SCET, only: getdynamictau
      implicit none
      include 'as3vjet_constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'as3vjet_masses.f'
      include 'as3vjet_ewcouple.f'
      include 'as3vjet_zcouple.f'
      include 'as3vjet_ewcharge.f'
      include 'as3vjet_qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'taucut.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'toploops.f'
      include 'beamtype.f'

      real(dp), intent(in) :: p(as3vj_mxpart,4), xx(2), z1, z2
      integer, intent(in) :: order
      real(dp), intent(out) :: xmsq
      logical, intent(in) :: central

      ! to omit gg->Zg 1L^2 contribution, c.f. NNLOJET
      logical, parameter :: omitggZg=.true.

      integer:: j,k,m,n

      real(dp):: fac,facgg,hard(2), &
       soft1qg(-1:1),soft2qg(-1:3),soft2qg_nab(-1:3), &
       soft1gq(-1:1),soft2gq(-1:3),soft2gq_nab(-1:3), &
       soft1qa(-1:1),soft2qa(-1:3),soft2qa_nab(-1:3), &
       beama0(-5:5),beamb0(-5:5), &
       beama1(-5:5,-1:1),beamb1(-5:5,-1:1), &
       beama2(-5:5,-1:3),beamb2(-5:5,-1:3), &
       jet1q(-1:1),jet2q(-1:3),jet1g(-1:1),jet2g(-1:3), &
       QB(2),bit, &
       assemblejet,y12,y15,y25,Iijm(5),msqsave(-nf:nf,-nf:nf)

      real(dp):: msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)

      real(dp) :: tauc, ptpure, puremass
      real(dp) :: twomass, pttwo, gettaucut

      tauc=gettaucut(twomass(3,4,p),pttwo(3,4,p))
!      if (dynamictau) then
!         tauc = getdynamictau(p, taucut)
!      else
!          tauc = taucut
!      endif

      call qqb_z1jet_vbis(p,msq0,0)
      if (order > 0) then
        call qqb_z1jet_vbis(p,msq1,1)
!        write(6,*) 'vbis msq1(2,-2)',msq1(2,-2)
!        msqsave(:,:)=msq1(:,:)
! The alternative code below can be used to account for the axial
! contribution at NLO, given by the 'N^a_V' term in Eq. (70) of 1309.3245;
! this term is otherwise not included (absent in qqb_z1jet_vbis)
!         epinv=0._dp
!         epinv2=0._dp
!         call as3vj_qqb_z1jet_v(p,msq1)
!         if (onlyaxial .eqv. .false.) then
! ! translation from DR (MCFM) to tH-V (MG)
!         msq1(:,:)=msq1(:,:)-(CF+XN/six)*ason2pi*msq0(:,:)
! ! extra finite terms from (1+e^2*pisq/12)
!         msq1(:,:)=msq1(:,:)+(2*CF+XN)*pisq/12d0*ason2pi*msq0(:,:)
!         endif
! !        write(6,*) 'v    msq1(2,-2)',msq1(2,-2)
! !        write(6,*) '     msq0(2,-2)',msq0(2,-2)
! !        write(6,*) 'diff/LO/ason2pi',(msq1(2,-2)-msqsave(2,-2))/pisq/msq0(2,-2)/ason2pi
! !        pause
      endif
      if (order > 1) then
        if (onlyaxial) then
          msq2=0._dp
        else
          call qqb_z1jet_vbis(p,msq2,2)
        endif
      endif
!      print *,'msq0',msq0
!      print *,'msq1',msq1
!      do j=-nf,nf
!         do k=-nf,nf
!            print *,'msq2',j,k,msq2(j,k)
!         enddo
!      enddo
      if (comparedecay) then
        call averageoverZ(p,msq0,0)
        call averageoverZ(p,msq1,1)
        call averageoverZ(p,msq2,2)
        write(6,*) 'msq0(2,-2)',msq0(2,-2)
        write(6,*) 'msq1(2,-2)',msq1(2,-2)
        write(6,*) 'msq2(2,-2)',msq2(2,-2)
        write(6,*) 'msq0(2,0)',msq0(2,0)
        write(6,*) 'msq1(2,0)',msq1(2,0)
        write(6,*) 'msq2(2,0)',msq2(2,0)
        write(6,*) 'msq0(0,2)',msq0(0,2)
        write(6,*) 'msq1(0,2)',msq1(0,2)
        write(6,*) 'msq2(0,2)',msq2(0,2)
        pause
      endif

      if (order >= 0) then
!        call fdist(ih1,xx(1),facscale,beama0,1)
!        call fdist(ih2,xx(2),facscale,beamb0,2)
        call pdf_fdist(1,xx(1),facscale,beama0)
        call pdf_fdist(2,xx(2),facscale,beamb0)
     endif
!     print *,'beam1',ih1,xx(1),facscale,beama0
!     print *,'beam2',ih2,xx(2),facscale,beamb0
      if (order >= 1) then
        call xbeam1bis(ih1,z1,xx(1),QB(1),beama1,1)
        call xbeam1bis(ih2,z2,xx(2),QB(2),beamb1,2)
      endif
      if (order >= 2) then
        call xbeam2bis(ih1,z1,xx(1),QB(1),beama2,1)
        call xbeam2bis(ih2,z2,xx(2),QB(2),beamb2,2)
      endif

      call jetq(order,two*p(5,4),jet1q,jet2q, disttau=.true.)
      call jetg(order,two*p(5,4),jet1g,jet2g)

      y12=s(1,2)/p(1,4)/p(2,4)/four
      y15=s(1,5)/p(1,4)/p(5,4)/four
      y25=s(2,5)/p(2,4)/p(5,4)/four
      
      call computeIijm(y12,y25,y15,Iijm)
      call soft_ab_qgq(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft1qg,soft2qg)
      call soft_ab_qgq(order,y12,y15,y25,Iijm,2,1,3,5,4,3,soft1gq,soft2gq)
      call soft_ab_qag(order,y12,y15,y25,Iijm,1,2,3,4,5,3,soft1qa,soft2qa)
      if (order > 1) then
        call soft_nab_qgq(order,y12,y25,y15,Iijm,1,2,3,4,5,3,soft2qg_nab)
        call soft_nab_qgq(order,y12,y15,y25,Iijm,2,1,3,5,4,3,soft2gq_nab)
        call soft_nab_qag(order,y12,y15,y25,Iijm,1,2,3,4,5,3,soft2qa_nab)
        soft2qg(:)=soft2qg(:)+soft2qg_nab(:)
        soft2gq(:)=soft2gq(:)+soft2gq_nab(:)
        soft2qa(:)=soft2qa(:)+soft2qa_nab(:)
      endif

      xmsq=zip
      do j=-nf,nf; do k=-nf,nf
          if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) cycle

          if (order > 0) then
            hard(1)=msq1(j,k)/msq0(j,k)/ason2pi
          else
            hard(1)=zip
          endif
          if (order > 1) then
            hard(2)=msq2(j,k)/msq0(j,k)/ason2pi**2
          else
            hard(2)=zip
          endif
          
          if ((j == 0) .and. (k == 0)) then
!--- g    g->Zg only has a hard(2) contribution
            if ((order > 1) .and. (omitggZg .eqv. .false.)) then
              bit=beama0(j)*beamb0(k)*msq2(j,k)
            else
              bit=zip
            endif
            goto 20
          elseif (j*k == 0) then ! initial state gluon -> jet is a quark
            if (j == 0) then
              bit=assemblejet(order,tauc, &
               beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
               beama2(j,:),beamb2(k,:),soft1gq,soft2gq,jet1q,jet2q,hard)
            elseif (k == 0) then
              bit=assemblejet(order,tauc, &
               beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
               beama2(j,:),beamb2(k,:),soft1qg,soft2qg,jet1q,jet2q,hard)
            endif
          else               ! initial state quarks -> jet is a gluon
              bit=assemblejet(order,tauc, &
               beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
               beama2(j,:),beamb2(k,:),soft1qa,soft2qa,jet1g,jet2g,hard)
          endif
!          print *,j,k,bit,beama0(j),beamb0(k)
          bit=bit*msq0(j,k)

 20       continue
!          print *,j,k,msq0(j,k),msq1(j,k),msq2(j,k),beama0(j),beamb0(k)
          xmsq=xmsq+bit

      enddo; enddo

      if (central .and. doMultitaucut) then
          scetreweight(:) = 0._dp
          if (xmsq /= 0._dp) then
              do m=1,size(tcutarray)
                 tauc=gettaucut(twomass(3,4,p),pttwo(3,4,p))
!                  if (dynamictau) then
!                      tauc = getdynamictau(p, tcutarray(m))
!                  else
!                      tauc = tcutarray(m)
!                  endif

                  scetreweight(m) = 0._dp !getxmsq_z(p,xx,order,soft1,soft2,hard, &
     !                beama0,beamb0,beama1,beamb1,beama2,beamb2,fac,qqb,qbq,prop)
     

!!!! BEGIN COPIED BLOCK FOR SCALEREWEIGHT 
      do j=-nf,nf; do k=-nf,nf
          if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) cycle

          if (order > 0) then
            hard(1)=msq1(j,k)/msq0(j,k)/ason2pi
          else
            hard(1)=zip
          endif
          if (order > 1) then
            hard(2)=msq2(j,k)/msq0(j,k)/ason2pi**2
          else
            hard(2)=zip
          endif
          
          if ((j == 0) .and. (k == 0)) then
!--- g    g->Zg only has a hard(2) contribution
            if ((order > 1) .and. (omitggZg .eqv. .false.)) then
              bit=beama0(j)*beamb0(k)*msq2(j,k)
            else
              bit=zip
            endif
            goto 30
          elseif (j*k == 0) then ! initial state gluon -> jet is a quark
            if (j == 0) then
              bit=assemblejet(order,tauc, &
               beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
               beama2(j,:),beamb2(k,:),soft1gq,soft2gq,jet1q,jet2q,hard)
            elseif (k == 0) then
              bit=assemblejet(order,tauc, &
               beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
               beama2(j,:),beamb2(k,:),soft1qg,soft2qg,jet1q,jet2q,hard)
            endif
          else               ! initial state quarks -> jet is a gluon
              bit=assemblejet(order,tauc, &
               beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
               beama2(j,:),beamb2(k,:),soft1qa,soft2qa,jet1g,jet2g,hard)
          endif
          
          bit=bit*msq0(j,k)

 30       continue

          scetreweight(m)=scetreweight(m)+bit

      enddo; enddo
!!!! END COPIED BLOCK FOR SCALEREWEIGHT 

              enddo
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif
      
      end subroutine lumxmsq_z1jet

      function passed_taucut_1jet(pparton,scetreweight_local,taucut_in)
          use SCET
          use ieee_arithmetic
      implicit none
      include 'as3vjet_constants.f'
      include 'mxpart.f'
      include 'npart.f'
      include 'nqcdjets.f'
      include 'taucut.f'
      include 'plabel.f'
!      include 'first.f'
      include 'kprocess.f'
      include 'kpart.f'

      
      logical :: passed_taucut_1jet
      real(dp), intent(inout), optional :: scetreweight_local(:)
      real(dp), intent(in) :: pparton(as3vj_mxpart,4)
      real(dp), intent(in), optional :: taucut_in

      integer i,j
      real(dp) :: tau,taua,taub,tauj,Q(4),tauc
      real(dp) :: p56(4), p57(4), p67(4), p567(4)

      real(dp) :: puremass, dot, ptpure, dotvec

      real(dp) :: taubeam, taujet
      real(dp) :: tauallbeam, taubeamjet, taualljet

      real(dp) :: Qsq

      logical :: bin
      common/bin/bin

      if (present(scetreweight_local)) then
          scetreweight_local(:) = 0._dp
      endif

      passed_taucut_1jet = .false.

      if (present(taucut_in)) then
          tauc = taucut_in
      else
          tauc = taucut
      endif

      !Qsq = -dotvec(pparton(2,:)+pparton(3,:)+pparton(4,:)+pparton(5,:), &
      !                  pparton(2,:)+pparton(3,:)+pparton(4,:)+pparton(5,:))

#define E(i) pparton(i,4)
#define pz(i) (pparton(i,3))
      p56(:) = pparton(5,:) + pparton(6,:)

      if (origKpart == ksnlo .or. (present(taucut_in) .and. kpart == kvirt)) then
          tau = min( &
              ! extra emission clustered with the beam
              min(E(5)-abs(pz(5)), E(6)-abs(pz(6))), &
              ! extra emission clustered with the jet
              E(5)+E(6) - sqrt((p56(1)**2 + p56(2)**2 + p56(3)**2)) )

      elseif (origKpart == knnlo .or. (present(taucut_in) .and. kpart == kreal)) then
          p57(:) = pparton(5,:) + pparton(7,:)
          p67(:) = pparton(6,:) + pparton(7,:)
          p567(:) = pparton(5,:) + pparton(6,:) + pparton(7,:)

          tauallbeam = min(E(5)+E(6)-abs(pz(5))-abs(pz(6)), &
                           E(5)+E(7)-abs(pz(5))-abs(pz(7)), &
                           E(6)+E(7)-abs(pz(6))-abs(pz(7)))
          taubeamjet = min(E(5)+E(6)+E(7) - sqrt(p56(1)**2+p56(2)**2+p56(3)**2) - abs(pz(7)), &
                           E(5)+E(6)+E(7) - sqrt(p57(1)**2+p57(2)**2+p57(3)**2) - abs(pz(6)), &
                           E(5)+E(6)+E(7) - sqrt(p67(1)**2+p67(2)**2+p67(3)**2) - abs(pz(5)))
          taualljet = E(5)+E(6)+E(7) - sqrt(p567(1)**2+p567(2)**2+p567(3)**2)

          tau = min(tauallbeam,taubeamjet,taualljet)
      else
          error stop "unknown kpart in maketaucut_singletop"
      endif

      if (ieee_is_nan(tau)) then
        write(6,*) 'maketaucut.f:  tau=',tau
        stop
      endif

      if (dynamictau) then
          tau = tau**2/(getdynamictau(pparton, tau))
      endif


      if (bin .and. doMultitaucut .and. present(scetreweight_local)) then
          ! no sampled value for taus will pass cuts, if smaller than the smallest taucut
          if (tau < smallestTaucut*(tauc/taucut)) then
              scetreweight_local(:) = 0._dp
              return
          endif

          ! otherwise compute "weights" for other taucuts
          do j=1,size(tcutarray)
              if (tau < tcutarray(j)*(tauc/taucut)) then
                  scetreweight_local(j) = 0._dp
              else
                  scetreweight_local(j) = 1._dp
              endif
          enddo

          ! and postpone this cut for later in the *int routines
          ! and in the plotting
          if (tau < tauc) then
              return
          endif
      else
         if (tau < tauc) return

      endif

      passed_taucut_1jet = .true.

      end function passed_taucut_1jet

      subroutine Zampqqbgsq(order,p1,p2,p3,p4,p5,NFV,za,zb,amps0,hard1,hard2)
          use types
          use nnlo_z1jet_hfun
          use constants
!--- returns msq0, hard1, hard2 which can be used to reconstruct the
!--- complete hard function for ab -> Z+c using:
!---  H = amps0 * (1 + [as/2/pi]*hard1 + [as/2/pi]^2*hard2)
!---
!--- note that amps0, hard1 and hard2 are arrays wjere
!---    1st index is helicity of quark line
!---    2nd index is helicity of lepton line
      implicit none
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'hpls.f'
      include 'scale.f'
      include 'as3vjet_ewcharge.f'
      include 'as3vjet_ewcouple.f'
      include 'as3vjet_qcdcouple.f'
      include 'as3vjet_zcouple.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> q^+(1)+qb^1(2)+g^+(3)+l^-(4)+lb^+(5)
!     according to Eq.22 of 1309.3245v3
!     1st index is helicity quark line
!     2nd index is helicity lepton line
!     helicity of gluon is summe.e-_dpover
!
!     Order of calculation is passed in via integer "order"
!     Region is passed in via integer "region"
!
      integer, intent(in) :: order,p1,p2,p3,p4,p5
      real(dp), intent(in) :: NFV
      !      complex(dp), intent(in) :: za(as3vj_mxpart,as3vj_mxpart), zb(as3vj_mxpart,as3vj_mxpart)
      complex(dp), intent(in) :: za(6,6), zb(6,6)
      real(dp), intent(out) :: amps0(2,2), hard1(2,2), hard2(2,2)

      real(dp), parameter :: three = 3._dp, four = 4._dp, eight = 8._dp
      real(dp), parameter :: one = 1._dp, two = 2._dp, six = 6._dp
      real(dp), parameter :: be0 = 11/three*CA-4/three*TR*NF
      real(dp), parameter :: be1= 34/three*CA**2-(20/three*CA+4*CF)*TR*NF

      integer:: qq1,qq2,qq3,hq,hl,k,iperm,region
      complex(dp):: al,be,ga,de
      complex(dp):: amps(0:2,2,2,2),ampqqbgll
      real(dp):: hn,s12,s13,s23,qsq,res0,res1,res2,msq0,msq1,msq2, &
       evolve1Lqqb,evolve2Lqqb,c1
      real(dp):: uu,vv,amps1(2,2),amps2(2,2),xlf
      complex(dp):: alp0,bet0,gam0,alp1,bet1,gam1,alp2,bet2,gam2
      real(dp):: alR(0:2),alI(0:2)
      real(dp):: beR(0:2),beI(0:2)
      real(dp):: gaR(0:2),gaI(0:2)
      real(dp):: LgamHqq0,LgamHqq1,Lrat

      complex(dp) :: T1(6)
      complex(dp) :: schemeconvC0, schemeconv2lM0

      logical, parameter:: checkBLLM=.false. ! use to check result in Eq. (70)

      logical, parameter:: peter=.false. ! use original PeTeR functions in C++
      
    real(dp) :: alpha2a_c0re,alpha2a_c0im,alpha2a_c1re,alpha2a_c1im,alpha2a_c2re,alpha2a_c2im
    real(dp) :: beta2a_c0re,beta2a_c0im,beta2a_c1re,beta2a_c1im,beta2a_c2re,beta2a_c2im
    real(dp) :: gamma2a_c0re,gamma2a_c0im,gamma2a_c1re,gamma2a_c1im,gamma2a_c2re,gamma2a_c2im
    real(dp) :: alpha3a_c0re,alpha3a_c0im,alpha3a_c1re,alpha3a_c1im,alpha3a_c2re,alpha3a_c2im
    real(dp) :: beta3a_c0re,beta3a_c0im,beta3a_c1re,beta3a_c1im,beta3a_c2re,beta3a_c2im
    real(dp) :: gamma3a_c0re,gamma3a_c0im,gamma3a_c1re,gamma3a_c1im,gamma3a_c2re,gamma3a_c2im
    real(dp) :: alpha4a_c0re,alpha4a_c0im,alpha4a_c1re,alpha4a_c1im,alpha4a_c2re,alpha4a_c2im
    real(dp) :: beta4a_c0re,beta4a_c0im,beta4a_c1re,beta4a_c1im,beta4a_c2re,beta4a_c2im
    real(dp) :: gamma4a_c0re,gamma4a_c0im,gamma4a_c1re,gamma4a_c1im,gamma4a_c2re,gamma4a_c2im

      real(dp), parameter:: &
       UGamHqq0= + 2*CA + 4*CF, &
       UGamHqq1= &
           - 40/9._dp*CA*TR*nf &
           + 134/9._dp*CA**2 &
           - 2/3._dp*CA**2*pisq &
           - 80/9._dp*CF*TR*nf &
           + 268/9._dp*CF*CA &
           - 4/3._dp*CF*CA*pisq

!--- statement function for hn of Eq.(29)
      hn(s12,s23,s13,qsq,alp0,bet0,gam0,alp1,bet1,gam1) &
       =one/(s13*s23*qsq)*( &
       +real(alp0*conjg(alp1),dp)*s12*(2*s12*qsq+s13*s23) &
       +real(bet0*conjg(bet1),dp)*s23*(2*s23*qsq+s12*s13) &
       +real(gam0*conjg(gam1),dp)*s13*(2*s13*qsq+s12*s23) &
       +real((alp0*conjg(bet1)+alp1*conjg(bet0)),dp)*s12*s23*(2*qsq-s13) &
       -real((alp0*conjg(gam1)+alp1*conjg(gam0)),dp)*s12*s13*s23 &
       -real((bet0*conjg(gam1)+bet1*conjg(gam0)),dp)*s13*s23*(2*qsq-s12))

!--- statement function for one-loop evolution
      evolve1Lqqb(musq)= &
          -UGamHqq0*log(abs(s(p1,p2)/musq))**2/two &
          -LgamHqq0*log(abs(s(p1,p2)/musq))

!--- statement function for two-loop evolution
      evolve2Lqqb(musq,c1)= &
       + log(abs(s(p1,p2)/musq))**4/eight * ( &
          + UGamHqq0**2 &
          ) &
       + log(abs(s(p1,p2)/musq))**3/six * ( &
          + 3*UGamHqq0*LgamHqq0 &
          + be0*UGamHqq0 &
          ) &
       + log(abs(s(p1,p2)/musq))**2/two * ( &
          + LgamHqq0**2 &
          - UGamHqq1 &
          - c1*UGamHqq0 &
          + be0*LgamHqq0 &
          ) &
       + log(abs(s(p1,p2)/musq)) * ( &
          - LgamHqq1 &
          - c1*LgamHqq0 &
          - be0*c1 &
          )

      xlf=real(nf,dp)

      if (checkBLLM) then
!---- ! Check Becher et al. result at a single PS point
        if ((p1 == 2) .and. (p2 == 1)) then
        s(p1,p2)=1._dp
        s(p1,p3)=-0.4_dp
        s(p4,p5)=0.1_dp**2
        s(p2,p3)=s(p4,p5)-s(p1,p2)-s(p1,p3)
        s(p2,p1)=s(p1,p2)
        s(p3,p1)=s(p1,p3)
        s(p3,p2)=s(p2,p3)
        endif
        if ((p1 == 5) .and. (p2 == 1)) then
        s(p2,p3)=1._dp
        s(p1,p3)=-0.4_dp
        s(p4,p5)=0.1_dp**2
        s(p1,p2)=s(p4,p5)-s(p1,p3)-s(p2,p3)
        s(p2,p1)=s(p1,p2)
        s(p3,p1)=s(p1,p3)
        s(p3,p2)=s(p2,p3)
        endif
        scale=0.6_dp
        musq=scale**2
      endif

!      write(6,*) 'region',region
!      write(6,*) 's(p1,p2)',s(p1,p2)
!      write(6,*) 's(p1,p3)',s(p1,p3)
!      write(6,*) 's(p2,p3)',s(p2,p3)

!--- logarithm appearing in evolution expressions
!--- note: take absolute value in order to work for all crossings
      Lrat=log(abs(s(p1,p2)**2/(s(p1,p3)*s(p2,p3))))
!--- constants for evolution components
      LgamHqq0 = &
           - 2*CA*Lrat &
           - 6*CF
      LgamHqq1 = &
           + 256/27._dp*CA*TR*nf &
           - 2/9._dp*CA*TR*nf*pisq &
           + 40/9._dp*CA*Lrat*TR*nf &
           - 692/27._dp*CA**2 &
           + 2*CA**2*zeta3 &
           + 11/18._dp*CA**2*pisq &
           - 134/9._dp*CA**2*Lrat &
           + 2/3._dp*CA**2*Lrat*pisq &
           + 368/27._dp*CF*TR*nf &
           + 4/3._dp*CF*TR*nf*pisq &
           - 961/27._dp*CF*CA &
           + 52*CF*CA*zeta3 &
           - 11/3._dp*CF*CA*pisq &
           - 3*CF**2 &
           - 48*CF**2*zeta3 &
           + 4*CF**2*pisq &
           + be1

      res0=0._dp
      res1=0._dp
      res2=0._dp
      do iperm=1,2
      if (iperm == 1) then
        qq1=p1
        qq2=p2
        qq3=p3
      else
        qq1=p2
        qq2=p1
        qq3=p3
      endif

!--- determine correct region from invariants
      if (s(qq1,qq2) > 0) then
        region=2
      elseif (s(qq1,qq3) > 0) then
        region=3
      elseif (s(qq2,qq3) > 0) then
        region=4
      else
        write(6,*) 'region not found:'
        write(6,*) 's12,s13,s23',s(qq1,qq2),s(qq1,qq3),s(qq2,qq3)
      endif
!      write(6,*) 'iperm,region',iperm,region

      if (region == 2) then
        uu=-s(qq1,qq3)/s(qq1,qq2)
        vv=+s(p4,p5)/s(qq1,qq2)
      elseif (region == 3) then
        uu=-s(qq2,qq3)/s(qq1,qq3)
        vv=+s(p4,p5)/s(qq1,qq3)
      elseif (region == 4) then
        uu=-s(qq1,qq3)/s(qq2,qq3)
        vv=+s(p4,p5)/s(qq2,qq3)
      endif
            
      if (order > 0) then
!--- fill arrays for 2DHPLs
        call tdhpl(uu,vv,4,G1,G2,G3,G4,H1,H2,H3,H4)
        call sethplh(H1,H2,H3,H4)
        call sethplg(G1,G2,G3,G4)
      else
!--- fill with dummy values
        G1=0._dp
        G2=0._dp
        G3=0._dp
        G4=0._dp
        H1=0._dp
        H2=0._dp
        H3=0._dp
        H4=0._dp
      endif

      if (.not.peter) then
!--- Coefficients for Region 2
      if (region == 2) then
        alR(0)=Alpha_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        alI(0)=Alpha_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beR(0)=Beta_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beI(0)=Beta_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaR(0)=Gamma_2a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaI(0)=Gamma_2a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)

        if (order > 0) then
          alR(1)=Alpha_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(1)=Alpha_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(1)=Beta_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(1)=Beta_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(1)=Gamma_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(1)=Gamma_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif
      
        if (order > 1) then
          alR(2)=Alpha_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(2)=Alpha_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(2)=Beta_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(2)=Beta_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(2)=Gamma_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(2)=Gamma_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif



        if (order > 0) then

#define CHECK_OMEGA 0
#if CHECK_OMEGA == 1
        call omega_V3j_2a(uu,vv,G1,G2,G3,G4,H1,H2,H3,H4, xlf, nfv, T1)

        schemeconvC0 = (18*CA*im*Pi - 36*CF*im*Pi + CA*Pi**2 + 2*CF*Pi**2 + 6*CA*im**2*Pi**2 - 12*CF*im**2*Pi**2 +  &
          6*(CA - 2*CF)*(3 + 2*im*Pi)*Log(s(p4,p5)/s(qq1,qq2)) + 6*(CA - 2*CF)*Log(s(p4,p5)/s(qq1,qq2))**2 +  &
          (-20*CA + 4*NF*TR)*Log(-(s(p4,p5)/s(qq1,qq3))) - 6*CA*Log(-(s(p4,p5)/s(qq1,qq3)))**2 -  &
          20*CA*Log(-(s(p4,p5)/s(qq2,qq3))) + 4*NF*TR*Log(-(s(p4,p5)/s(qq2,qq3))) - 6*CA*Log(-(s(p4,p5)/s(qq2,qq3)))**2)/24._dp

        schemeconv2lM0 = (-2*(CA*(-67 + 3*Pi**2) + 20*NF*TR)*(72*CA*im*Pi - 144*CF*im*Pi + CA*Pi**2 + 2*CF*Pi**2 + 24*CA*im**2*Pi**2 -  &
             48*CF*im**2*Pi**2 + 24*(CA - 2*CF)*(3 + 2*im*Pi)*Log(s(p4,p5)/s(qq1,qq2)) +  &
             24*(CA - 2*CF)*Log(s(p4,p5)/s(qq1,qq2))**2 + (-80*CA + 16*NF*TR)*Log(-(s(p4,p5)/s(qq1,qq3))) -  &
             24*CA*Log(-(s(p4,p5)/s(qq1,qq3)))**2 - 80*CA*Log(-(s(p4,p5)/s(qq2,qq3))) + 16*NF*TR*Log(-(s(p4,p5)/s(qq2,qq3))) -  &
             24*CA*Log(-(s(p4,p5)/s(qq2,qq3)))**2) +  &
          3*(18*CA*im*Pi - 36*CF*im*Pi + CA*Pi**2 + 2*CF*Pi**2 + 6*CA*im**2*Pi**2 - 12*CF*im**2*Pi**2 +  &
              6*(CA - 2*CF)*(3 + 2*im*Pi)*Log(s(p4,p5)/s(qq1,qq2)) + 6*(CA - 2*CF)*Log(s(p4,p5)/s(qq1,qq2))**2 +  &
              (-20*CA + 4*NF*TR)*Log(-(s(p4,p5)/s(qq1,qq3))) - 6*CA*Log(-(s(p4,p5)/s(qq1,qq3)))**2 -  &
              20*CA*Log(-(s(p4,p5)/s(qq2,qq3))) + 4*NF*TR*Log(-(s(p4,p5)/s(qq2,qq3))) - 6*CA*Log(-(s(p4,p5)/s(qq2,qq3)))**2)**2 &
           - 2*(11*CA - 4*NF*TR)*(11*CA*Pi**2 + 18*CF*Pi**2 - 108*CA*im**2*Pi**2 + 216*CF*im**2*Pi**2 - 6*CA*im*Pi**3 +  &
             12*CF*im*Pi**3 - 24*CA*im**3*Pi**3 + 48*CF*im**3*Pi**3 - 4*NF*Pi**2*TR + 36*CA*Zeta3 + 72*CF*Zeta3 -  &
             6*(CA - 2*CF)*Pi*(36*im + Pi + 12*im**2*Pi)*Log(s(p4,p5)/s(qq1,qq2)) -  &
             36*(CA - 2*CF)*(3 + 2*im*Pi)*Log(s(p4,p5)/s(qq1,qq2))**2 - 24*(CA - 2*CF)*Log(s(p4,p5)/s(qq1,qq2))**3 +  &
             6*CA*Pi**2*Log(-(s(p4,p5)/s(qq1,qq3))) + 24*(5*CA - NF*TR)*Log(-(s(p4,p5)/s(qq1,qq3)))**2 +  &
             24*CA*Log(-(s(p4,p5)/s(qq1,qq3)))**3 + 6*CA*Pi**2*Log(-(s(p4,p5)/s(qq2,qq3))) +  &
             120*CA*Log(-(s(p4,p5)/s(qq2,qq3)))**2 - 24*NF*TR*Log(-(s(p4,p5)/s(qq2,qq3)))**2 +  &
             24*CA*Log(-(s(p4,p5)/s(qq2,qq3)))**3))/3456._dp


        ! for the two-loop results C0 multiplies the 1-loop formfactor
        ! and there is another piece schemeconv2lM0 that multiplies 
        ! the tree formfactor, which is zero for gamma

        ! alpha2
        T1(4) = T1(4) + schemeconvC0*T1(1) + schemeconv2lM0

        ! beta2
        T1(5) = T1(5) + schemeconvC0*T1(2) + schemeconv2lM0

        ! gamma2
        T1(6) = T1(6) + schemeconvC0*T1(3)


        ! alpha1
        T1(1) = T1(1) + schemeconvC0
        ! beta1
        T1(2) = T1(2) + schemeconvC0
        ! gamma1
        ! T1(3) stays, since the first order of gamma = 0


        write (*,*) "s12", s(qq1,qq2), s(qq1,qq3), s(qq2,qq3)
        write (*,*) "schemeconv", 2._dp*schemeconvC0
        write (*,*) "schemeconv2lM0", 4._dp*schemeconv2lM0

        ! the factors of 2 at NLO and 4 at NNLO are (likely) due to some as/4/pi vs. as/2/pi mismatch

        write (*,*) "order", order
        write (*,*) "u,v", uu, vv, xlf, nfv
        write (*,*) "alpha 1", alR(1), alI(1)
        write (*,*) "alpha 1 difference", 2._dp*T1(1)-complex(alR(1),alI(1))
        write (*,*) "alpha 2", alR(2), alI(2)
        write (*,*) "alpha 2 difference", 4._dp*T1(4)-complex(alR(2),alI(2))
        write (*,*) ""

        write (*,*) "beta", beR(1), beI(1)
        write (*,*) "beta 1", T1(2)
        write (*,*) "beta 1 difference", 2._dp*T1(2)-complex(beR(1),beI(1))
        write (*,*) "beta 2", beR(2), beI(2)
        write (*,*) "beta 2 difference", 4._dp*T1(5)-complex(beR(2),beI(2))
        write (*,*) ""

        write (*,*) "gamma", gaR(0), gaI(0)
        write (*,*) "gamma 1", T1(3)
        write (*,*) "gamma 1 difference", 2._dp*T1(3)-complex(gaR(1),gaI(1))
        write (*,*) "gamma 2", gaR(2), gaI(2)
        write (*,*) "gamma 2 difference", 4._dp*T1(6)-complex(gaR(2),gaI(2))
        write (*,*) ""

        pause
#endif

        endif
      
!--- Coefficients for Region 3
      elseif (region == 3) then
        alR(0)=Alpha_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        alI(0)=Alpha_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beR(0)=Beta_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beI(0)=Beta_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaR(0)=Gamma_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaI(0)=Gamma_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)

        if (order > 0) then
          alR(1)=Alpha_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(1)=Alpha_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(1)=Beta_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(1)=Beta_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(1)=Gamma_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(1)=Gamma_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif

        if (order > 1) then
          alR(2)=Alpha_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(2)=Alpha_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(2)=Beta_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(2)=Beta_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(2)=Gamma_3a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(2)=Gamma_3a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif

!--- Coefficients for Region 4
      elseif (region == 4) then
        alR(0)=Alpha_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        alI(0)=Alpha_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beR(0)=Beta_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        beI(0)=Beta_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaR(0)=Gamma_4a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        gaI(0)=Gamma_4a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)

        if (order > 0) then
          alR(1)=Alpha_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(1)=Alpha_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(1)=Beta_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(1)=Beta_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(1)=Gamma_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(1)=Gamma_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif

        if (order > 1) then
          alR(2)=Alpha_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          alI(2)=Alpha_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beR(2)=Beta_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          beI(2)=Beta_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaR(2)=Gamma_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
          gaI(2)=Gamma_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
        endif

      else
        write(6,*) 'Zampqqbgsq: region should be 2,3 or 4: ',region
        stop
      endif

   else
        ! Coefficients from PeTeR C functions
        !--- Coefficients for Region 2
        if (region == 2) then
            alR(0)=alpha2a_c0re(uu,vv)
            alI(0)=alpha2a_c0im(uu,vv)
            beR(0)=beta2a_c0re(uu,vv) 
            beI(0)=beta2a_c0im(uu,vv) 
            gaR(0)=gamma2a_c0re(uu,vv)
            gaI(0)=gamma2a_c0im(uu,vv)
            alR(1)=alpha2a_c1re(uu,vv)
            alI(1)=alpha2a_c1im(uu,vv)
            beR(1)=beta2a_c1re(uu,vv) 
            beI(1)=beta2a_c1im(uu,vv) 
            gaR(1)=gamma2a_c1re(uu,vv)
            gaI(1)=gamma2a_c1im(uu,vv)
            alR(2)=alpha2a_c2re(uu,vv)
            alI(2)=alpha2a_c2im(uu,vv)
            beR(2)=beta2a_c2re(uu,vv) 
            beI(2)=beta2a_c2im(uu,vv) 
            gaR(2)=gamma2a_c2re(uu,vv)
            gaI(2)=gamma2a_c2im(uu,vv)
            
         !--- Coefficients for Region 3
         elseif (region == 3) then
            alR(0)=alpha3a_c0re(uu,vv)
            alI(0)=alpha3a_c0im(uu,vv)
            beR(0)=beta3a_c0re(uu,vv) 
            beI(0)=beta3a_c0im(uu,vv) 
            gaR(0)=gamma3a_c0re(uu,vv)
            gaI(0)=gamma3a_c0im(uu,vv)
            alR(1)=alpha3a_c1re(uu,vv)
            alI(1)=alpha3a_c1im(uu,vv)
            beR(1)=beta3a_c1re(uu,vv) 
            beI(1)=beta3a_c1im(uu,vv) 
            gaR(1)=gamma3a_c1re(uu,vv)
            gaI(1)=gamma3a_c1im(uu,vv)
            alR(2)=alpha3a_c2re(uu,vv)
            alI(2)=alpha3a_c2im(uu,vv)
            beR(2)=beta3a_c2re(uu,vv) 
            beI(2)=beta3a_c2im(uu,vv) 
            gaR(2)=gamma3a_c2re(uu,vv)
            gaI(2)=gamma3a_c2im(uu,vv)
            
         !--- Coefficients for Region 4
         elseif (region == 4) then
            alR(0)=alpha4a_c0re(uu,vv)
            alI(0)=alpha4a_c0im(uu,vv)
            beR(0)=beta4a_c0re(uu,vv) 
            beI(0)=beta4a_c0im(uu,vv) 
            gaR(0)=gamma4a_c0re(uu,vv)
            gaI(0)=gamma4a_c0im(uu,vv)
            alR(1)=alpha4a_c1re(uu,vv)
            alI(1)=alpha4a_c1im(uu,vv)
            beR(1)=beta4a_c1re(uu,vv) 
            beI(1)=beta4a_c1im(uu,vv) 
            gaR(1)=gamma4a_c1re(uu,vv)
            gaI(1)=gamma4a_c1im(uu,vv)
            alR(2)=alpha4a_c2re(uu,vv)
            alI(2)=alpha4a_c2im(uu,vv)
            beR(2)=beta4a_c2re(uu,vv) 
            beI(2)=beta4a_c2im(uu,vv) 
            gaR(2)=gamma4a_c2re(uu,vv)
            gaI(2)=gamma4a_c2im(uu,vv)
         endif
      endif
!      print *,'region',region,alR(2),alI(2),beR(2),beI(2),gaR(2),gaI(2),NFV
      do k=0,2
      al=cmplx(alR(k),alI(k),dp)
      be=cmplx(beR(k),beI(k),dp)
      ga=cmplx(gaR(k),gaI(k),dp)
      de=s(qq1,qq2)*(al-be-ga)/(2._dp*(s(p1,p2)+s(p2,p3)+s(p3,p1)))

!     1st index is order of amplitude
!     2nd index is helicity quark line
!     3rd index is helicity lepton line
!     4th index is helicity of outgoing gluon line
      if (iperm == 2) then
      amps(k,1,1,2)=-ampqqbgll(p2,p1,p3,p4,p5,al,be,ga,de,za,zb)
      amps(k,2,2,1)=-ampqqbgll(p2,p1,p3,p4,p5,al,be,ga,de,zb,za)
      amps(k,1,2,2)=-ampqqbgll(p2,p1,p3,p5,p4,al,be,ga,de,za,zb)
      amps(k,2,1,1)=-ampqqbgll(p2,p1,p3,p5,p4,al,be,ga,de,zb,za)
      endif
      
      if (iperm == 1) then
      amps(k,1,1,1)=+ampqqbgll(p1,p2,p3,p5,p4,al,be,ga,de,zb,za)
      amps(k,2,2,2)=+ampqqbgll(p1,p2,p3,p5,p4,al,be,ga,de,za,zb)
      amps(k,1,2,1)=+ampqqbgll(p1,p2,p3,p4,p5,al,be,ga,de,zb,za)
      amps(k,2,1,2)=+ampqqbgll(p1,p2,p3,p4,p5,al,be,ga,de,za,zb)
      endif
      
      enddo
 
!---- This bit of code for testing squared MEs against Becher et al.
      alp0=cmplx(alR(0),alI(0),dp)
      bet0=cmplx(beR(0),beI(0),dp)
      gam0=cmplx(gaR(0),gaI(0),dp)
      alp1=cmplx(alR(1),alI(1),dp)
      bet1=cmplx(beR(1),beI(1),dp)
      gam1=cmplx(gaR(1),gaI(1),dp)
      alp2=cmplx(alR(2),alI(2),dp)
      bet2=cmplx(beR(2),beI(2),dp)
      gam2=cmplx(gaR(2),gaI(2),dp)

      res0=res0 &
       +two*hn(s(qq1,qq2),s(qq2,qq3),s(qq1,qq3),s(p4,p5), &
               alp0,bet0,gam0,alp0,bet0,gam0)
     
      res1=res1 &
       +two*hn(s(qq1,qq2),s(qq2,qq3),s(qq1,qq3),s(p4,p5), &
              alp0,bet0,gam0,alp1,bet1,gam1)
     
      res2=res2 &
       +two*hn(s(qq1,qq2),s(qq2,qq3),s(qq1,qq3),s(p4,p5), &
              alp0,bet0,gam0,alp2,bet2,gam2) &
           +hn(s(qq1,qq2),s(qq2,qq3),s(qq1,qq3),s(p4,p5), &
              alp1,bet1,gam1,alp1,bet1,gam1)
     
      enddo ! end of iperm loop
      
!---- This bit of code for testing squared MEs against Becher et al.
      msq0=res0
      if (order > 0) then
        msq1=two*res1
!--- RG running to restore scale dependence from musq=s(p4,p5)
        msq1=msq1+(evolve1Lqqb(musq)-evolve1Lqqb(s(p4,p5)))*res0
      endif
      if (order > 1) then
        msq2=two*res2
!--- RG running to restore scale dependence from musq=s(p4,p5)
        c1=two*res1/res0-evolve1Lqqb(s(p4,p5))
        msq2=msq2 &
         +(evolve2Lqqb(musq,c1) &
          -evolve2Lqqb(s(p4,p5),c1))*res0
      endif
      
      if (checkBLLM .and. (order == 2)) then
        write(6,*) 'qqbg p1,p2,p3,p4,p5,NFV = ',p1,p2,p3,p4,p5,NFV
        write(6,*) 'qqbg LO          ',res0
        write(6,*) 'qqbg 1L/LO/fourpi',msq1/res0/fourpi
        write(6,*) 'qqbg 2L/LO/fourpi**2',msq2/res0/fourpi**2
        pause
      endif
      if (comparedecay .and. (order == 2)) then
        if ((p1<3) .and. (p2<3)) then
        write(6,*) 'LO not decayed', &
         msq0*esq*(L(2)**2+R(2)**2)/two*gsq*8._dp/4._dp/3._dp**2 !V*aveqq
        write(6,*) '1L not decayed', &
         msq1*esq*(L(2)**2+R(2)**2)/two*gsq*ason4pi*8._dp/4._dp/3._dp**2 !V*aveqq
        write(6,*) '2L not decayed', &
         msq2*esq*(L(2)**2+R(2)**2)/two*gsq*ason4pi**2*8._dp/4._dp/3._dp**2 !V*aveqq
        else
        write(6,*) 'LO not decayed', &
         -msq0*esq*(L(2)**2+R(2)**2)/two*gsq*8._dp/4._dp/3._dp/8._dp !V*aveqg
        write(6,*) '1L not decayed', &
         -msq1*esq*(L(2)**2+R(2)**2)/two*gsq*ason4pi*8._dp/4._dp/3._dp/8._dp !V*aveqg
        write(6,*) '2L not decayed', &
         -msq2*esq*(L(2)**2+R(2)**2)/two*gsq*ason4pi**2*8._dp/4._dp/3._dp/8._dp !V*aveqg
        endif
      endif
      
      do hq=1,2
      do hl=1,2
        amps0(hq,hl)=abs(amps(0,hq,hl,1))**2+abs(amps(0,hq,hl,2))**2
        if (order > 0) then
          amps1(hq,hl)=two*( &
           +real(amps(0,hq,hl,1)*conjg(amps(1,hq,hl,1)),dp) &
           +real(amps(0,hq,hl,2)*conjg(amps(1,hq,hl,2)),dp))
!--- RG running to restore scale dependence from musq=s(p4,p5)
          amps1(hq,hl)=amps1(hq,hl) &
           +(evolve1Lqqb(musq)-evolve1Lqqb(s(p4,p5)))*amps0(hq,hl)
        else
          amps1(hq,hl)=0._dp
        endif
        if (order > 1) then
          amps2(hq,hl)=two*( &
             +real(amps(0,hq,hl,1)*conjg(amps(2,hq,hl,1)),dp) &
             +real(amps(0,hq,hl,2)*conjg(amps(2,hq,hl,2)),dp)) &
           +real(amps(1,hq,hl,1)*conjg(amps(1,hq,hl,1)),dp) &
           +real(amps(1,hq,hl,2)*conjg(amps(1,hq,hl,2)),dp)
!--- RG running to restore scale dependence from musq=s(p4,p5)
          c1=amps1(hq,hl)/amps0(hq,hl)-evolve1Lqqb(musq)
          amps2(hq,hl)=amps2(hq,hl) &
           +(evolve2Lqqb(musq,c1) &
            -evolve2Lqqb(s(p4,p5),c1))*amps0(hq,hl)
        else
          amps2(hq,hl)=0._dp
        endif
      enddo
      enddo
      
!--- convert from coefficient of [as/4/pi]^n to coefficient of [as/2/pi]^n
      amps1(:,:)=amps1(:,:)/two
      amps2(:,:)=amps2(:,:)/four
      msq1=msq1/two
      msq2=msq2/four

!--- normalize coefficients to LO
      hard1(:,:)=amps1(:,:)/amps0(:,:)
      hard2(:,:)=amps2(:,:)/amps0(:,:)
      
      return
      end

      subroutine Zampgggsq(order,p1,p2,p3,p4,p5,za,zb,amps2)
          use types
          use nnlo_z1jet_hfun
          use constants
!--- returns msq0, hard1, hard2 which can be used to reconstruct the
!--- complete hard function for gg -> Z+g using:
!---  H = [as/2/pi]^2 * amps2
!---
!--- note that amps2 is an array where index is helicity of lepton line
      implicit none
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'hpls.f'
      include 'scale.f'
!     This routine calculates the amplitude for
!     0--> g(1)+g(2)+g(3)+l^-(4)+lb^+(5)
!     according to Eq.22 of 1309.3245v3

      integer, intent(in) :: order,p1,p2,p3,p4,p5
!      complex(dp), intent(in) :: za(as3vj_mxpart,as3vj_mxpart), zb(as3vj_mxpart,as3vj_mxpart)
      complex(dp), intent(in) :: za(6,6), zb(6,6)
      real(dp), intent(out) :: amps2(2)

      logical, parameter:: checkBLLM=.false. ! use to check result in Eq. (70)

      real(dp), parameter :: zip = 0._dp, one = 1._dp, &
          four = 4._dp, two = 2._dp, eight = 8._dp
      complex(dp), parameter :: czip = (0._dp,0._dp)
!
!     1st index is helicity quark line
!     2nd index is helicity lepton line
!     helicity of gluon is summed over
!
!     Order of calculation is passed in via integer "order"
!
      integer:: q1,q2,q3,hl,ha,hb,hc,j,k,iperm,region
      complex(dp):: al,be,ga,de,alC,beC,gaC,deC
      complex(dp):: amps(2,2,2,2)
      real(dp):: hnpmm,hnppp,s12,s13,s23,qsq,res0,res1,res2, &
       msq0,msq1,msq2,c1
      real(dp):: uu,vv
      complex(dp):: alpa0,alpa1,alpb0,alpb1,alpc0,alpc1
      complex(dp):: beta0,beta1,betb0,betb1,betc0,betc1
      real(dp):: al1R,al1I,al2R,al2I,al3R,al3I,be1R,be1I,be2R,be2I,be3R,be3I

!--- statement function for hn of Eqs.(42), (43)
      hnpmm(s12,s23,s13,qsq,alpa0,alpb0,alpc0,alpa1,alpb1,alpc1) &
       =one/(eight*qsq*s12*s13)*( &
       +real(alpa0*conjg(alpa1),dp)*qsq*(s23*qsq+2*s12*s13) &
       +real(alpb0*conjg(alpb1),dp)*s23*(s12+s23)**2 &
       +real(alpc0*conjg(alpc1),dp)*s23*(s13+s23)**2 &
       -real((alpa0*conjg(alpb1)+alpa1*conjg(alpb0)),dp) &
          *s23*(s12+s23)*qsq &
       +real((alpa0*conjg(alpc1)+alpa1*conjg(alpc0)),dp) &
          *s23*(s13+s23)*qsq &
       -real((alpb0*conjg(alpc1)+alpb1*conjg(alpc0)),dp) &
          *s23*(s23*qsq-s12*s13))
      hnppp(s12,s23,s13,qsq,beta0,betb0,betc0,beta1,betb1,betc1) &
       =one/(eight*qsq*s12)*( &
       +real(beta0*conjg(beta1),dp)*s13*(s12+s13)**2/s23 &
       +real(betb0*conjg(betb1),dp)*s23*(s12+s23)**2/s13 &
       +real(betc0*conjg(betc1),dp)*(2*s12*qsq+s13*s23) &
       -real((beta0*conjg(betb1)+beta1*conjg(betb0)),dp) &
          *(s12*qsq-s13*s23) &
       +real((beta0*conjg(betc1)+beta1*conjg(betc0)),dp) &
          *s13*(s12+s13) &
       +real((betb0*conjg(betc1)+betb1*conjg(betc0)),dp) &
          *s23*(s12+s23))

!--- contribution only enters at NNLO
      if (order < 2) then
        amps2(:)=zip
        return
      endif

      if (checkBLLM) then
!---- ! Check Becher et al. result at a single PS point
        s(p1,p2)=1._dp
        s(p1,p3)=-0.4_dp
        s(p4,p5)=0.1_dp**2
        s(p2,p3)=s(p4,p5)-s(p1,p2)-s(p1,p3)
        s(p2,p1)=s(p1,p2)
        s(p3,p1)=s(p1,p3)
        s(p3,p2)=s(p2,p3)
        scale=0.6_dp
        musq=scale**2
      endif

!      write(6,*) 'region',region
!      write(6,*) 's(p1,p2)',s(p1,p2)
!      write(6,*) 's(p1,p3)',s(p1,p3)
!      write(6,*) 's(p2,p3)',s(p2,p3)

      res2=zip
      amps(:,:,:,:)=czip
      do iperm=1,3
      if     (iperm == 1) then
        q1=p1
        q2=p2
        q3=p3
      elseif (iperm == 2) then
        q1=p2
        q2=p1
        q3=p3
      elseif (iperm == 3) then
        q1=p3
        q2=p2
        q3=p1
      endif

!--- determine correct region from invariants
      if (s(q1,q2) > 0) then
        region=2
      elseif (s(q1,q3) > 0) then
        region=3
      elseif (s(q2,q3) > 0) then
        region=4
      else
        write(6,*) 'region not found:'
        write(6,*) 's12,s13,s23',s(q1,q2),s(q1,q3),s(q2,q3)
      endif
!      write(6,*) 'iperm,region',iperm,region

      if (region == 2) then
        uu=-s(q1,q3)/s(q1,q2)
        vv=+s(p4,p5)/s(q1,q2)
      elseif (region == 3) then
        write(6,*) 'Region 3 not anticipated in Zampgggsq'
        stop
        uu=-s(q2,q3)/s(q1,q3)
        vv=+s(p4,p5)/s(q1,q3)
      elseif (region == 4) then
        uu=-s(q1,q3)/s(q2,q3)
        vv=+s(p4,p5)/s(q2,q3)
      endif
      
!--- fill arrays for 2DHPLs
      call tdhpl(uu,vv,4,G1,G2,G3,G4,H1,H2,H3,H4)

!--- Coefficients for Region 2: need alpha and beta
      if (region == 2) then
        al1R=GGAlpha1_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al1I=GGAlpha1_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al2R=GGAlpha2_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al2I=GGAlpha2_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al3R=GGAlpha3_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al3I=GGAlpha3_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)

        be1R=GGBeta1_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        be1I=GGBeta1_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        be2R=GGBeta2_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        be2I=GGBeta2_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        be3R=GGBeta3_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        be3I=GGBeta3_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
      
!--- Coefficients for Region 4: only need alpha
      elseif (region == 4) then
        al1R=GGAlpha1_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al1I=GGAlpha1_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al2R=GGAlpha2_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al2I=GGAlpha2_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al3R=GGAlpha3_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
        al3I=GGAlpha3_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4)
      
      else
        write(6,*) 'Zampgggsq: region should be 2 or 4: ',region
        stop

      endif

      alpa1=cmplx(al1R,al1I,dp)
      alpb1=cmplx(al2R,al2I,dp)
      alpc1=cmplx(al3R,al3I,dp)
      beta1=cmplx(be1R,be1I,dp)
      betb1=cmplx(be2R,be2I,dp)
      betc1=cmplx(be3R,be3I,dp)

! result without Z decay
      if (iperm == 1) then
      res2=res2 &
       +two*hnppp(s(q1,q2),s(q2,q3),s(q1,q3),s(p4,p5), &
              beta1,betb1,betc1,beta1,betb1,betc1)
      endif  
      res2=res2 &
       +two*hnpmm(s(q1,q2),s(q2,q3),s(q1,q3),s(p4,p5), &
              alpa1,alpb1,alpc1,alpa1,alpb1,alpc1)
      
!     1st index is helicity of lepton line
!     2nd,3rd,4th indices are helicities of gluons
      if (iperm == 1) then
      amps(2,2,2,2)=ampgggppp(q1,q2,q3,p4,p5,beta1,betb1,betc1,za,zb)
      amps(1,2,2,2)=ampgggppp(q1,q2,q3,p5,p4,beta1,betb1,betc1,za,zb)
      amps(1,1,1,1)=ampgggppp(q1,q2,q3,p4,p5,beta1,betb1,betc1,zb,za)
      amps(2,1,1,1)=ampgggppp(q1,q2,q3,p5,p4,beta1,betb1,betc1,zb,za)

      amps(2,2,1,1)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,za,zb)
      amps(1,2,1,1)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,za,zb)
      amps(1,1,2,2)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,zb,za)
      amps(2,1,2,2)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,zb,za)
      endif

      if (iperm == 2) then
      amps(2,1,2,1)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,za,zb)
      amps(1,1,2,1)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,za,zb)
      amps(1,2,1,2)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,zb,za)
      amps(2,2,1,2)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,zb,za)
      endif
      
      if (iperm == 3) then
      amps(2,1,1,2)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,za,zb)
      amps(1,1,1,2)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,za,zb)
      amps(1,2,2,1)=ampgggpmm(q1,q2,q3,p4,p5,alpa1,alpb1,alpc1,zb,za)
      amps(2,2,2,1)=ampgggpmm(q1,q2,q3,p5,p4,alpa1,alpb1,alpc1,zb,za)
      endif
 
      enddo ! end of iperm loop
      
!---- This bit of code for testing squared MEs against Becher et al.
      msq2=res2
      if (checkBLLM) then
        write(6,*) 'ggg 2L/fourpi**2',msq2/fourpi**2
        pause
      endif
!      write(6,*) '2L not decayed',
!     &   (avegg*msq2*(40._dp/3._dp)*esq*gsq*ason4pi**2)
!     &  *(half*(2._dp*(L(2)+R(2))+3._dp*(L(1)+R(1))))**2
      
      do hl=1,2
        amps2(hl)=zip
        do ha=1,2
        do hb=1,2
        do hc=1,2
        amps2(hl)=amps2(hl) &
         +real(amps(hl,ha,hb,hc)*conjg(amps(hl,ha,hb,hc)),dp)
!        write(6,*) hl,ha,hb,hc,amps(hl,ha,hb,hc)
        enddo
        enddo
        enddo
      enddo
      
!--- convert from coefficient of [as/4/pi]^n to coefficient of [as/2/pi]^n
      amps2(:)=amps2(:)/four

      return
      end

      function ampgggpmm(p1,p2,p3,p5,p6,alpa,alpb,alpc,za,zb)
      implicit none
      include 'types.f'
      include 'as3vjet_constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> g^+(1)+g^-(2)+g^-(3)+l^+(5)+lb^-(6)
!     according to Eq.(36) of 1309.3245v3
      integer:: p1,p2,p3,p5,p6
      complex(dp)::ampgggpmm
      complex(dp)::alpa,alpb,alpc
      ampgggpmm=za(p2,p3)/(za(p1,p2)*za(p1,p3)*zb(p2,p3))*( &
       +alpa*za(p2,p5)*za(p3,p5)*zb(p5,p6) &
       +alpb*za(p2,p3)*za(p2,p5)*zb(p2,p6) &
       +alpc*za(p2,p3)*za(p3,p5)*zb(p3,p6))/rt2
     
! Overall factor of (-im) present in Eq. (2.18) of 1302.2630
!      ampgggpmm=(-im)*ampgggpmm

      return
      end

      function ampgggppp(p1,p2,p3,p5,p6,beta,betb,betc,za,zb)
      implicit none
      include 'types.f'
      include 'as3vjet_constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> g^+(1)+g^+(2)+g^+(3)+l^+(5)+lb^-(6)
!     according to Eq.(37) of 1309.3245v3
      integer:: p1,p2,p3,p5,p6
      complex(dp)::ampgggppp
      complex(dp)::beta,betb,betc
      ampgggppp=( &
       +beta*zb(p1,p3)*za(p1,p5)*zb(p1,p6)/(za(p1,p2)*za(p2,p3)) &
       +betb*zb(p2,p3)*za(p2,p5)*zb(p2,p6)/(za(p1,p2)*za(p1,p3)) &
       +betc*zb(p2,p3)*za(p2,p5)*zb(p1,p6)/(za(p1,p2)*za(p2,p3)))/rt2

! Overall factor of (-im) present in Eq. (2.18) of 1302.2630
!      ampgggppp=(-im)*ampgggppp

      return
      end


end module
