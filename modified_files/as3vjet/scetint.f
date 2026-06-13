!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function scetint(r,wgt,f)
          use ieee_arithmetic
          use types
          use SCET
          use SafetyCuts, only : passed_smallnew
      implicit none
      real(dp):: scetint
      double precision f(*)
      double precision dot,dytpt,mjj
      include 'as3vjet_constants.f'
      include 'mxpart.f'
      include 'vegas_common.f'
      include 'mxdim.f'
      include 'sprods_com.f'
      include 'energy.f'
      include 'npart.f'
      include 'taucut.f'
      include 'as3vjet_qcdcouple.f'
      include 'x1x2.f'
      include 'xmin.f'
      include 'nproc.f'
      include 'options.f'
      real(dp):: p(as3vj_mxpart,4),pjet(as3vj_mxpart,4),r(mxdim),W,xmsq,
     & val,val2,pswt,wgt,z1,z2,flux,BrnRat
      integer:: j
      logical:: as3vjet_includedipole
      real(dp) :: xjac
      common/BrnRat/BrnRat

      real(dp) :: scet_xmsq
      
      logical as3vj_smalls
      external as3vj_smalls

!      integer:: ndim
      ndim = 9
      
      f(1) = 0._dp
      scetint=0._dp

      W=sqrts**2

      p(:,:)=0._dp
!      pjet(:,:)=0._dp

      call gen_lops(r,p,pswt,*999)
!      if (debug) then
!          if (all(.not. ieee_is_nan(p(1:npart,:))) .eqv. .false.) then
!             write(6,*) 'Discarding NaN or infinite phase space point'
!             goto 999
!          endif
!      endif
      !print *,'npart',npart
      call as3vj_dotem(npart+2,p,s)

      !print *,'p'
      !print *,p(1,:)
      !print *,p(2,:)
      !print *,p(3,:)
      !print *,p(4,:)
      !print *,p(5,:)
      !print *,'s'
      !print *,s(:,:)      
!      print *,'try smallnew'
c      if (ntau == 0) then
cc----reject event if any tau is too small -- important for precision
c         call smalltau(p,npart,*999)
c      else
c small safety cuts
!         if (.not. passed_smallnew(p,npart)) then
!             goto 999
!         endif

      if (.not. as3vj_smalls(s,p,npart)) then
         goto 999
      endif
c      endif
!      print *,'smallnew done'

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (as3vjet_includedipole(0,p) .eqv. .false.) then
        goto 999
      endif
!      print *,'as3vjet_includedipole done'

      z1=r(ndim-1)**2
      z2=r(ndim)**2

      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts

      if ( (xx(1) > one)  .or. (xx(2) > one)
     & .or.(xx(1) < xmin) .or. (xx(2) < xmin)) goto 999

      call scaleset_mcfm(sqrt(2*dot(p,3,4)),
     .     dytpt(p(4,:),p(3,:)),
     .     mjj(p(5,:),p(6,:)))
      xmsq = scet_xmsq(z1,z2,p,.true.)
      xjac = four*sqrt(z1*z2)
      flux = fbGeV2/(two*xx(1)*xx(2)*W)
      scetint = flux*xjac*pswt*xmsq/BrnRat
!      print *,'scetint',xmsq,xjac,flux,pswt,BrnRat
      call getptildejet(0,pjet)
      val=scetint*wgt
      val2=val**2

      if (ieee_is_nan(val) .or. (.not. ieee_is_finite(val))) then
!          if (debug) then
!              write(6,*) 'Discarded NaN, val=',val
!          endif
          goto 999
      endif

      includeTaucutgrid(0) = .true.
      
C     Fill only if it's last iteration
      if (doFill.ne.0) then
          call hists_fill(p(3,:),p(4,:),scetint*wgt)
      endif
      
      f(1) = scetint
!      print *,'phasespace'
!      do j=1,12
!         print *,'x[',j-1,']',r(j)
!      enddo
!      print *,'scetint',scetint
      return

 999  continue

      scetint = 0._dp
      end

