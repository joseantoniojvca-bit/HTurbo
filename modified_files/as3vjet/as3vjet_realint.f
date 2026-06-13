!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function as3vj_realint(vector,wgt)
          use ieee_arithmetic
          use SafetyCuts, only : passed_smallnew
          use nnlo_z1jet, only: genps_z1jet_rr
      implicit none
      include 'types.f'
      include 'as3vjet_constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'noglue.f'
      include 'nflav.f'
      include 'mxdim.f'
      include 'vegas_common.f'
      include 'ptilde.f'
      include 'xmin.f'
      include 'new_pspace.f'
      include 'npart.f'
      include 'scale.f'
      include 'facscale.f'
      include 'maxwt.f'
      include 'kprocess.f'
      include 'as3vjet_masses.f'
      include 'dipolescale.f'
      include 'flags.f'
      include 'ipsgen.f'
      include 'outputoptions.f'
      include 'breit.f'
      include 'runstring.f'
      include 'energy.f'
      include 'incldip.f'
      include 'nproc.f'
      include 'initialscales.f'
      include 'nqcdjets.f'
      include 'taucut.f'
      include 'as3vjet_qcdcouple.f'
      include 'couple.f'
      include 'nlooprun.f'
      include 'cutoff.f'
      include 'ewcorr.f'
      include 'kpart.f'
      include 'beamtype.f'
      include 'plabel.f'
      include 'options.f'
      include 'sprods_com.f'
      real(dp):: as3vj_realint
      integer:: j,k,nd,nvec
      real(dp):: vector(mxdim),W,ptmp
      real(dp):: fx1(-nf:nf),fx2(-nf:nf),
     & dipfx1(0:maxd,-nf:nf),dipfx2(0:maxd,-nf:nf)
      real(dp):: p(as3vj_mxpart,4),pjet(as3vj_mxpart,4),p1ext(4),p2ext(4)
      real(dp):: pswt
!     real(dp):: s(as3vj_mxpart,as3vj_mxpart),wgt,msq(-nf:nf,-nf:nf)
      real(dp):: wgt,msq(-nf:nf,-nf:nf)
      real(dp):: msqc(maxd,-nf:nf,-nf:nf),xmsq(0:maxd)
      real(dp):: flux,BrnRat
      real(dp):: xx1,xx2,q(as3vj_mxpart,4)
      real(dp):: xmsqjk
      logical:: as3vjet_includedipole,includereal
      real(dp):: QandGint
      external qqb_w2jet_g,qqb_w2jet_gs,qqb_z2jet_g,qqb_z2jet_gs
      common/Pext/p1ext,p2ext
      common/BrnRat/BrnRat

      real(dp) :: plo(as3vj_mxpart,4)
      real(dp) :: pswtdip

      real(dp):: wt34,wt345,wt347,wt3457,wtprop,s34,s345,s347,s3457,wtips(4)

      double precision dot
      external dot

      double precision dytpt
      external dytpt

      double precision mjj
      external mjj

      logical isnan_ofast
      external isnan_ofast

      logical as3vj_smalls
      external as3vj_smalls

c--- statement function
      wtprop(s34,wmass,wwidth)=(s34-wmass**2)**2+(wmass*wwidth)**2

      as3vj_realint = 0._dp
      ndim = 13
c      call init_is_functions;
      
      QandGflag = .false.
      !p(1:npart+2,:) = 0._dp
      !pjet(1:npart+2,:) = 0._dp
      p = 0._dp
      pjet = 0._dp
      pswt = 0._dp

      W=sqrts**2

      if ( genps_z1jet_rr(vector,p,pswt) .eqv. .false. ) then
         goto 999
      endif

      if (all(.not. ieee_is_nan(p(1:npart+2,4))) .eqv. .false.) then
          goto 999
      endif

      nvec=npart+2
      call as3vj_dotem(nvec,p,s)

c----calculate the x's for the incoming partons from generated momenta
      xx1=two*(p(1,4)*p2ext(4)-p(1,3)*p2ext(3))/W
      xx2=two*(p(2,4)*p1ext(4)-p(2,3)*p1ext(3))/W

c      if (debug) write(*,*) 'Reconstructed x1,x2 ',xx1,xx2

      if ((xx1 >  1._dp) .or. (xx2 >  1._dp)
     &.or.(xx1 < xmin) .or. (xx2 < xmin)) then
         goto 999
      endif

c--- (moved to includedipole) impose cuts on final state
c      call masscuts(p,*999)

c small safety cuts
!      if (.not. passed_smallnew(p,npart)) then
!         goto 999
!      endif

      if (.not. as3vj_smalls(s,p,npart)) then
         goto 999
      endif
      
c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead set to zero
      includereal=as3vjet_includedipole(0,p)
      incldip(0)=includereal

c     it is always a good idea to initialize all elements as zero
      msq(:,:) = 0
      msqc(:,:,:) = 0

c--- test to see whether we need Gflag and Qflag together
      if ( (Qflag) .and. (Gflag) ) then
        QandGflag=.true.
        QandGint=0._dp
c--- first pass: Gflag
        Gflag=.true.
        Qflag=.false.
      endif

c--- restart from here when calculating with Qflag and Gflag
c--- (W+2 jet and Z+2 jet processes only)
   44 continue

      call scaleset_mcfm(sqrt(2*dot(p,3,4)),
     .     dytpt(p(4,:),p(3,:)),
     .     mjj(p(5,:),p(6,:)))  ! should compute invariant mass of 5,6,7
      dipscale(0)=facscale
      
      as3vj_realint = 0._dp

c---  Calculate the required matrix elements
      if (kcase==kW_2jet) then
c        if (includereal) call singcheck(qqb_w2jet_g,qqb_w2jet_gs_new,p) ! Re-checked June 09
        if (includereal)  call qqb_w2jet_g(p,msq)
        call qqb_w2jet_gs_new(p,msqc)
      elseif (kcase==kZ_2jet) then
c        if (includereal) call singcheck(qqb_z2jet_g,qqb_z2jet_gs_new,p) ! Re-checked June 09
        if (includereal) call qqb_z2jet_g(p,msq)
        call qqb_z2jet_gs_new(p,msqc)
c        call qqb_z2jet_gs(p,msqc)
c        write(6,*) 'Gflag,Qflag',Gflag,Qflag
c        do j=-nf,nf
c        do k=-nf,nf
c        bit1=0._dp
c        bit2=0._dp
c        bit1=sum(msqc(1:ndmax,j,k))
c        bit2=sum(msqc_new(1:ndmax,j,k))
c        if ((abs(bit1)  >  1.e-100_dp) .or. (abs(bit2)  >  1.e-100_dp)) then
c          if (abs(bit1/bit2-1._dp)  >  1.e-12_dp) then
c            write(6,*) 'j,k,bit1,bit1/bit2',j,k,bit1,bit2,bit1/bit2-1._dp
c          endif
c        endif
c        enddo
c        enddo
c        pause
      endif

      xmsq = 0._dp

      flux=fbGeV2/(two*xx1*xx2*W)

      xmsq = 0._dp

c--- calculate PDF's
      do nd=ndmax,0,-1          ! so that fx1,fx2 correct for real kinematics
         call pdf_fdist(1,xx1,dipscale(nd),fx1)
         call pdf_fdist(2,xx2,dipscale(nd),fx2)
         dipfx1(nd,:)=fx1(:)
         dipfx2(nd,:)=fx2(:)
      enddo

      do j=-nflav,nflav
         do k=-nflav,nflav
            xmsqjk=fx1(j)*fx2(k)*msq(j,k)
!            print *,j,k,fx1(j),fx2(k),msq(j,k)
            xmsq(0)=xmsq(0)+xmsqjk
! add dipole contributions
            do nd=1,ndmax
               xmsqjk=dipfx1(nd,j)*dipfx2(nd,k)*(-msqc(nd,j,k))
!      print *,j,k,nd,dipfx1(nd,j),dipfx2(nd,k),(-msqc(nd,j,k))
               xmsq(nd)=xmsq(nd)+xmsqjk
            enddo
         enddo
      enddo

! summation and histogramming for all dipole contributions
      do nd=0,ndmax

         xmsq(nd)=xmsq(nd)*flux*pswt/BrnRat
         if (xmsq(nd) /= 0._dp) then
            q(:,:) = ptilde(nd,:,:)
! check cuts
            incldip(nd)=as3vjet_includedipole(nd,q)
            if (incldip(nd) .eqv. .false.) then
               call as3vj_dotem(nvec,p,s)
               xmsq(nd) = 0._dp
            endif
         endif
!      print *,'includedipole',nd,incldip(nd),
!     .         xmsq(nd),flux,pswt

         if (ieee_is_nan(xmsq(nd))) then
!     if (debug) write(6,*) 'discarding point with weight NaN: pswt=',pswt
            goto 999
         endif

         if (.not. ieee_is_finite(xmsq(nd))) then
!     if (debug) write(6,*) 'discarding point with infinite weight: pswt=',pswt
            goto 999
         endif

         if (isnan_ofast(xmsq(nd))) then
            print *,'NaN in as3vj_realint, dipole',nd,xmsq(nd)
            print *,'phasespace',vector
!     xmsq(nd) = 0
            goto 999
         endif
         
! histogramming
         as3vj_realint = as3vj_realint + xmsq(nd)

         call getptildejet(nd,pjet)
         call as3vj_dotem(nvec,pjet,s)
!         print *,'q3',q(3,:)
!         print *,'ptilde3',ptilde(nd,3,:)
!         print *,'q4',q(4,:)
!         print *,'ptilde4',ptilde(nd,4,:)
C     Fill only if it's last iteration
         if (doFill.ne.0) then
            call hists_real_dipole(q(3,:),q(4,:),xmsq(nd)*wgt,nd)
         endif
      enddo                     ! loop over dipole contributions

C---  Fill only if it's last iteration
      if (doFill.ne.0) then
C---  fill the dipole contribution to each bin separatelly
         call hists_real_event()
      endif
      
      if (abs(as3vj_realint*wgt) > wtmax) then
         wtmax=abs(as3vj_realint*wgt)
      endif

c--- handle special case of Qflag and Gflag
      if (QandGflag) then
        QandGint=QandGint+as3vj_realint
        if ((Gflag) .and. (.not.(Qflag))) then
c--- go back for second pass (Qflag), calling includedipole first to reset
c--- the values of "jets" and "jetlabel"
          includereal=as3vjet_includedipole(0,p)
          Qflag=.true.
          Gflag=.false.
          goto 44
        else
c--- return both to .true. and assign value to as3vj_realint (to return to VEGAS)
          Qflag=.true.
          Gflag=.true.
          as3vj_realint=QandGint
        endif
      endif

      return

      ! any errors where the whole point is set to zero
 999  continue

      as3vj_realint = 0._dp

      if (QandGflag) then
        Qflag=.true.
        Gflag=.true.
      endif

      return
      end
