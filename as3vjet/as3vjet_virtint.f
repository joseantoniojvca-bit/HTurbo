!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function as3vj_virtint(r,wgt)
        use ieee_arithmetic
        use SCET
        use SafetyCuts, only : passed_smallnew
      implicit none
      include 'types.f'
      include 'as3vjet_constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'noglue.f'
      include 'mxdim.f'
      include 'vegas_common.f'
      include 'sprods_com.f'
      include 'npart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'agq.f'
      include 'PR_new.f'
      include 'PR_cs_new.f'
      include 'PR_h2j.f'
      include 'PR_twojet.f'
      include 'PR_stop.f'
      include 'PR_mix.f'
      include 'msq_cs.f'
      include 'msq_struc.f'
      include 'msq_mix.f'
      include 'as3vjet_qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'lc.f'
      include 'kprocess.f'
      include 'maxwt.f'
      include 'limits.f'
      include 'nflav.f'
      include 'b0.f'
      include 'as3vjet_masses.f'
      include 'wts_bypart.f'
      include 'nores.f'
      include 'stopbmass.f'
      include 'as3vjet_ewcouple.f'
      include 'flags.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'xmin.f'
      include 'outputoptions.f'
      include 'runstring.f'
      include 'x1x2.f'
      include 'bypart.f'
      include 'energy.f'
      include 'first.f'
      include 'initialscales.f'
      include 'taucut.f'
      include 'ppmax.f'
      include 'kpart.f'
      include 'hbbparams.f'
      include 'ewcorr.f'
      include 'mpicommon.f'
      include 'couple.f'
      include 'nlooprun.f'
      include 'badpoint.f'
      include 'zcouple_cms.f'
      include 'cutoff.f'
      include 'nproc.f'
      include 'beamtype.f'
      include 'anomcoup.f'
      include 'TRbadpoint.f' 
      include 'toploops.f'
      include 'options.f'
      include 'plabel.f'
      real(dp):: as3vj_virtint
      real(dp), intent(in) :: wgt, r(mxdim)

      real(dp):: mqq(0:2,-nf:nf,-nf:nf),
     & msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: ppmsqx(0:2,ppmax)
      real(dp):: msqx_cs(0:2,-nf:nf,-nf:nf)
      real(dp):: AP(-1:1,-1:1,3),APqg_mass,AP_mix(-1:1,-1:1,0:3,3),APqed(-1:1,-1:1,3)

      integer:: j,k,m,n,cs,ics,csmax,nvec,is,iq,ia
      integer itrial
      real(dp):: xmsqvar(8)
      real(dp):: p(as3vj_mxpart,4),pjet(as3vj_mxpart,4),W,xmsq,pttwo,
     & val,val2,fx1(-nf:nf),fx2(-nf:nf),fx1z(-nf:nf),fx2z(-nf:nf),xmsqt,
     & fx1_H(-nf:nf),fx2_H(-nf:nf),fx1_L(-nf:nf),fx2_L(-nf:nf),
     & fx1z_H(-nf:nf),fx2z_H(-nf:nf),fx1z_L(-nf:nf),fx2z_L(-nf:nf)
      real(dp):: xmsq_noew,msq_noew(-nf:nf,-nf:nf),virtint_noew
      real(dp):: pswt,xjac,msqqcd(-nf:nf,-nf:nf),
     & msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),msqvdk(-nf:nf,-nf:nf),
     & msqvdkW(-nf:nf,-nf:nf),
     & msq_qq,msq_aa,msq_aq,msq_qa,msq_qg,msq_gq,epcorr,msq_up,msq_dn,
     & msq1(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),msqm1(-nf:nf,-nf:nf),xlog,
     & bit2(-nf:nf,-nf:nf),bit1(-nf:nf,-nf:nf),bit0(-nf:nf,-nf:nf)
      real(dp) :: msqv_light(-nf:nf,-nf:nf), msqv_heavy(-nf:nf,-nf:nf) ! for ktopanom
      real(dp):: z,x1onz,x2onz,flux,omz,ptthree,
     & BrnRat,xmsq_old,tmp,wwidth_save,zwidth_save
      real(dp):: scaleup,scaledn,facscaleup,facscaledn
      integer:: rvcolourchoice
      logical:: bin,as3vjet_includedipole,WWjetcheckpiDpjk
      real(dp):: QandGint

      common/bin/bin
      common/BrnRat/BrnRat
      common/rvcolourchoice/rvcolourchoice
c      common/ggZZunstable/ggZZunstable
c      data p/56*0._dp/
      integer, save:: nshot=1
      external gg_ZZ,qqb_w1jet_vbis
!$omp threadprivate(nshot)

      real(dp) fxa
      common/photonpdf/fxa
!$omp threadprivate(/photonpdf/)

      double precision dot
      external dot

      double precision dytpt
      external dytpt
      
      logical as3vj_smalls
      external as3vj_smalls

      ndim = 11
c      call init_is_functions;
      
      QandGflag=.false.
      first=.false.
      nshot=3
      if (first) then
         first=.false.
c         write(*,*) case
         nshot=1
      endif

      as3vj_virtint=0._dp
c--- ensure isolation code does not think this is fragmentation piece
      z_frag=0._dp

      p(:,:)=0._dp
      pjet(:,:)=0._dp

      W=sqrts**2

      call gen_lops(r,p,pswt,*999)

      if (all(.not. ieee_is_nan(p(1:npart+2,:))) .eqv. .false.) then
          goto 999
      endif

      nvec=npart+2
      call as3vj_dotem(nvec,p,s)

c small safety cuts
!      if (.not. passed_smallnew(p,npart)) then
!         goto 999
!      endif

      if (.not. as3vj_smalls(s,p,npart)) then
         goto 999
      endif

      as3vj_virtint=0._dp

c--- see whether this point will pass cuts - if it will not, do not
c--- bother calculating the matrix elements for it, instead bail out
      if (as3vjet_includedipole(0,p) .eqv. .false.) then
        goto 999
      endif
!      print *,'passed as3vjet_includedipole'
      
c--- test to see whether we need Gflag and Qflag together
      if ( ((kcase==kW_2jet) .or. (kcase==kZ_2jet))
     &.and. (Qflag) .and. (Gflag) ) then
        QandGflag=.true.
        QandGint=0._dp
c--- first pass: Gflag
        Gflag=.true.
        Qflag=.false.
      endif

!      print *,Qflag,Gflag
!      print*,'first pass: Gflag'
      
c--- restart from here when calculating with Qflag and Gflag
c--- (W+2 jet and Z+2 jet processes only)
   44 continue

   66 continue

 771  continue

      call scaleset_mcfm(sqrt(2*dot(p,3,4)),
     .     dytpt(p(4,:),p(3,:)),sqrt(2*dot(p,5,6)))
      
      xx(1)=-2._dp*p(1,4)/sqrts
      xx(2)=-2._dp*p(2,4)/sqrts
      if ((xx(1) >  1._dp) .or. (xx(2) >  1._dp)
     &.or.(xx(1) < xmin)   .or. (xx(2) < xmin)) then
         goto 999
      endif

      z=r(11)**2
! Exclude exceptional points
      if ((z == 0._dp) .or. (z == 1._dp)) goto 999
c      if (nshot == 1) z=0.95_dp
      xjac=two*sqrt(z)

      omz=1._dp-z

      flux=fbGeV2/(2._dp*xx(1)*xx(2)*W)
c--- for mlm study, divide by (Ecm)**2=W
c      if (runstring(1:3) == 'mlm') then
c      flux=flux/W
c      endif

c--- to test poles, we need colourchoice=0, but save real value
      if (nshot == 1) then
        rvcolourchoice=colourchoice
        colourchoice=0
      endif

   12 continue
c--- point to restart from when checking epsilon poles

c--- correction to epinv from AP subtraction when mu_FAC != mu_REN,
c--- corresponding to subtracting -1/epinv*Pab*log(musq_REN/musq_FAC)
      epcorr=epinv+2._dp*log(scale/facscale)

      AP(q,q,1)=+ason2pi*Cf*1.5_dp*epcorr
      AP(q,q,2)=+ason2pi*Cf*(-1._dp-z)*epcorr
      AP(q,q,3)=+ason2pi*Cf*2._dp/omz*epcorr
      AP(a,a,1)=+ason2pi*Cf*1.5_dp*epcorr
      AP(a,a,2)=+ason2pi*Cf*(-1._dp-z)*epcorr
      AP(a,a,3)=+ason2pi*Cf*2._dp/omz*epcorr

      AP(q,g,1)=0._dp
      AP(q,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
      AP(q,g,3)=0._dp
      AP(a,g,1)=0._dp
      AP(a,g,2)=ason2pi*Tr*(z**2+omz**2)*epcorr
      AP(a,g,3)=0._dp


      AP(g,q,1)=0._dp
      AP(g,q,2)=ason2pi*Cf*(1._dp+omz**2)/z*epcorr
      AP(g,q,3)=0._dp
      AP(g,a,1)=0._dp
      AP(g,a,2)=ason2pi*Cf*(1._dp+omz**2)/z*epcorr
      AP(g,a,3)=0._dp

      AP(g,g,1)=+ason2pi*b0*epcorr
      AP(g,g,2)=+ason2pi*xn*2._dp*(1._dp/z+z*omz-2._dp)*epcorr
      AP(g,g,3)=+ason2pi*xn*2._dp/omz*epcorr

      Q1 = zip
      Q2 = zip
      H1 = zip
      H2 = zip
      R1 = zip
      R2 = zip
      S1 = zip
      S2 = zip
      M1 = zip
      M2 = zip

      msqv = zip
      msq = zip
      msq_qq = zip
      msq_aq = zip
      msq_qa = zip

c---  Calculate the required matrix elements
!      print *,kcase,nproc
      if (kcase==kW_2jet) then
         if (nproc==22) call qqb_wp2jetx_new(p,msq,mqq,ppmsqx,msqx_cs)
         if (nproc==27) call qqb_wm2jetx_new(p,msq,mqq,ppmsqx,msqx_cs)
         call qqb_w2jet_v(p,msqv)
         call qqb_w2jet_z(p,z)
      elseif (kcase==kZ_2jet) then
c        call qqb_z2jetx(p,msq,mqq,msqx,msqx_cs)
         call qqb_z2jetx_new(p,msq,mqq,ppmsqx,msqx_cs)
         call qqb_z2jet_v(p,msqv)
         if (onlyaxial) then
            msq=0._dp
            mqq=0._dp
            ppmsqx=0._dp
            msqx_cs=0._dp
            msq_cs=0._dp
         else
!            call qqb_z2jetx_new(p,msq,mqq,ppmsqx,msqx_cs)
            call qqb_z2jet_z(p,z)
         endif
      endif
!      R1 = zip
!      R2 = zip
!      AP = zip
!      msqv = zip
!      msq = zip
!      ppmsqx = zip
!      msqx_cs = zip
!      msq_cs = zip
!      print *,'msq',msq      
!      print *,'msq_cs',msq_cs
!      print *,'mqq',mqq
!      print *,'ppmsqx',ppmsqx
!      print *,'msqx_cs',msqx_cs
!      do j=-nf,nf
!      do k=-nf,nf
!         print *,'mqq',j,k,mqq(0,j,k)
!         print *,'msq',j,k,mqq(0,j,k)
!         print *,'msqx_cs',j,k,msqx_cs(0,j,k)
!         print *,'msqv',j,k,msqv(j,k)
!      enddo
!      enddo
!      print *,'z',z
!      print *,'R1',R1
!      print *,'R2',R2
c--- explicitly remove factor of LO if we are only interested in coefficient
      if (coeffonly) then
        msqv(:,:)=msqv(:,:)-msq(:,:)
      endif

      currentNd = 0

c--- initialize a PDF set here, if calculating errors
  777 continue
      xmsq=0._dp
      fx1z(:)=0._dp
      fx2z(:)=0._dp

      if (z > xx(1)) x1onz=xx(1)/z
      if (z > xx(2)) x2onz=xx(2)/z
      xmsq_noew=0._dp

c--- calculate PDF's
      call pdf_fdist(1,xx(1),facscale,fx1)
      call pdf_fdist(2,xx(2),facscale,fx2)
      if (z > xx(1)) call pdf_fdist(1,x1onz,facscale,fx1z)
      if (z > xx(2)) call pdf_fdist(2,x2onz,facscale,fx2z)

      do j=-nflav,nflav
      do k=-nflav,nflav


      tmp=xmsq

c--- The variables R1 and R2 provide the Regular and Plus pieces associated
c--- with radiation from leg 1 (R1(a,b,c,cs,is)) and leg 2 (R2(a,b,c,cs,is))
c--- In each case the parton labelling is using the normal QM notation of
c--- putting everything backward
c---       emitted line after emission =    a
c---       emitter before emission     =    b
c---       spectator                   =    c
c--- There is no label for he or she who is emitted.
c--- Note that in general each piece will be composed of many different
c--- dipole contributions

c--- SUM BY COLOUR STRUCTURES: H+2jets only
      if  ( (kcase==kggfus2) .or. (kcase==kgagajj)
     &     .or. (kcase==kHWW2jt) .or. (kcase==kHZZ2jt)) then
       xmsq=xmsq+fx1(j)*fx2(k)*(
     & msqv(j,k)+msq(j,k))
c      write(6,*) j,k,'-> msqv = ',fx1(j)*fx2(k)*(
c     & msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))
c      tmp=xmsq

c--- quark-quark or antiquark-antiquark
      if (  ((j > 0).and.(k > 0))
     & .or. ((j < 0).and.(k < 0))) then
c      write(6,*) 'as3vj_virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      if (j == k) then
        m=+1
        n=+1
        csmax=6
      else
        m=abs(j)
        n=abs(k)
        csmax=6
      endif
      xmsqt=0._dp
      do cs=1,csmax
      xmsqt=xmsqt
     & +msq_struc(cs,m,n)*(AP(q,q,1)-AP(q,q,3)
     &                 +H1(q,q,q,cs,1)-H1(q,q,q,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +H2(q,q,q,cs,1)-H2(q,q,q,cs,3))*fx1(j)*fx2(k)
     & +(msq_struc(cs,m,n)*(AP(q,q,2)+AP(q,q,3)
     &                  +H1(q,q,q,cs,2)+H1(q,q,q,cs,3))
     & +msq_struc(cs,g,k)*(AP(g,q,2)+H1(g,q,q,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_struc(cs,m,n)*(AP(q,q,2)+AP(q,q,3)
     &                  +H2(q,q,q,cs,2)+H2(q,q,q,cs,3))
     & +msq_struc(cs,j,g)*(AP(g,q,2)+H2(g,q,q,cs,2)))*fx1(j)*fx2z(k)/z
      enddo
      xmsq=xmsq+xmsqt
c      write(6,*) 'as3vj_virtint: j,k,  ct=',j,k,xmsqt

c--- quark-antiquark or antiquark-quark
      elseif (  ((j > 0).and.(k < 0))
     &     .or. ((j < 0).and.(k > 0))) then
c      write(6,*) 'as3vj_virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      if (j == -k) then
        m=+1
        n=-1
        csmax=7
      else
        m=abs(j)
        n=-abs(k)
        csmax=6
      endif
      xmsqt=0._dp
      do cs=1,csmax
c      if ((cs > 3) .and. (cs < 7)) goto 67
c      do cs=4,6
      xmsqt=xmsqt
     & +msq_struc(cs,m,n)*(AP(q,q,1)-AP(q,q,3)
     &                 +H1(q,q,a,cs,1)-H1(q,q,a,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +H2(a,a,q,cs,1)-H2(a,a,q,cs,3))*fx1(j)*fx2(k)
     & +(msq_struc(cs,m,n)*(AP(q,q,2)+AP(q,q,3)
     &                  +H1(q,q,a,cs,2)+H1(q,q,a,cs,3))
     & + msq_struc(cs,g,k)*(AP(g,q,2)+H1(g,q,a,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_struc(cs,m,n)*(AP(a,a,2)+AP(a,a,3)
     &                  +H2(a,a,q,cs,2)+H2(a,a,q,cs,3))
     & + msq_struc(cs,j,g)*(AP(g,a,2)+H2(g,a,q,cs,2)))*fx1(j)*fx2z(k)/z
c   67 continue
      enddo
      xmsq=xmsq+xmsqt
c      write(6,*) 'as3vj_virtint: j,k,  ct=',j,k,xmsqt

c--- gluon-gluon
      elseif ((j == g) .and. (k == g)) then
c      write(6,*) 'as3vj_virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      xmsqt=0._dp
c--- loop up to 3 to cancel poles from gggg
c      do cs=1,3
      do cs=1,6
      msq_qg=msq_struc(cs,+5,g)+msq_struc(cs,+4,g)+msq_struc(cs,+3,g)
     &      +msq_struc(cs,+2,g)+msq_struc(cs,+1,g)
     &      +msq_struc(cs,-5,g)+msq_struc(cs,-4,g)+msq_struc(cs,-3,g)
     &      +msq_struc(cs,-2,g)+msq_struc(cs,-1,g)
      msq_gq=msq_struc(cs,g,+5)+msq_struc(cs,g,+4)+msq_struc(cs,g,+3)
     &      +msq_struc(cs,g,+2)+msq_struc(cs,g,+1)
     &      +msq_struc(cs,g,-5)+msq_struc(cs,g,-4)+msq_struc(cs,g,-3)
     &      +msq_struc(cs,g,-2)+msq_struc(cs,g,-1)
      xmsqt=xmsqt
     & +msq_struc(cs,g,g)*(AP(g,g,1)-AP(g,g,3)
     &                 +H1(g,g,g,cs,1)-H1(g,g,g,cs,3)
     &                 +AP(g,g,1)-AP(g,g,3)
     &                 +H2(g,g,g,cs,1)-H2(g,g,g,cs,3))*fx1(g)*fx2(g)
     & +(msq_struc(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +H1(g,g,g,cs,2)+H1(g,g,g,cs,3))
     &  +msq_qg*(AP(q,g,2)+H1(q,g,g,cs,2)))*fx1z(g)/z*fx2(g)
     & +(msq_struc(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +H2(g,g,g,cs,2)+H2(g,g,g,cs,3))
     &  +msq_gq*(AP(q,g,2)+H2(q,g,g,cs,2)))*fx1(g)*fx2z(g)/z
      enddo
      xmsq=xmsq+xmsqt
c      write(6,*) 'as3vj_virtint: j,k,  ct=',j,k,xmsqt

c--- quark-gluon and anti-quark gluon
      elseif ((j  /=  0) .and. (k == g)) then
c      write(6,*) 'as3vj_virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      m=+1
      n=0
      xmsqt=0._dp
      do cs=4,6
      xmsqt=xmsqt
     &+ msq_struc(cs,m,g)*(AP(q,q,1)-AP(q,q,3)
     &                 +H1(q,q,g,cs,1)-H1(q,q,g,cs,3)
     &                 +H2(g,g,q,cs,1)-H2(g,g,q,cs,3)
     &                 +AP(g,g,1)-AP(g,g,3))*fx1(j)*fx2(g)
     &+(msq_struc(cs,m,g)*(AP(q,q,2)+AP(q,q,3)
     &                    +H1(q,q,g,cs,2)+H1(q,q,g,cs,3))
     & +msq_struc(cs,g,g)*(AP(g,q,2)+H1(g,q,g,cs,2)))*fx1z(j)/z*fx2(g)
     &+(msq_struc(cs,m,g)*(AP(g,g,2)+AP(g,g,3)
     &                    +H2(g,g,q,cs,2)+H2(g,g,q,cs,3))
     & +msq_struc(cs,m,-m)*(AP(a,g,2)+H2(a,g,q,cs,2)))*fx1(j)*fx2z(g)/z
      enddo
      do cs=1,3
      msq_qa=msq_struc(cs,m,-1)+msq_struc(cs,m,-2)+msq_struc(cs,m,-3)
     &      +msq_struc(cs,m,-4)+msq_struc(cs,m,-5)
      msq_qq=msq_struc(cs,m,+1)+msq_struc(cs,m,+2)+msq_struc(cs,m,+3)
     &      +msq_struc(cs,m,+4)+msq_struc(cs,m,+5)
      xmsqt=xmsqt
     &+(msq_struc(cs,g,g)*(AP(g,q,2)+H1(g,q,g,cs,2)))*fx1z(j)/z*fx2(g)
     &+(+msq_qa*(AP(a,g,2)+H2(a,g,q,cs,2))
     &  +msq_qq*(AP(q,g,2)+H2(q,g,q,cs,2)))*fx1(j)*fx2z(g)/z
      enddo
      xmsqt=xmsqt
     & +msq_struc(iqr,m,-m)*(AP(a,g,2)+H2(a,g,q,iqr,2))*fx1(j)*fx2z(g)/z
      xmsq=xmsq+xmsqt
c      write(6,*) 'as3vj_virtint: j,k,  ct=',j,k,xmsqt
c      write(6,*) 'as3vj_virtint: j,k, SUM=',j,k,xmsqt+fx1(j)*fx2(k)*msqv(j,k)

c--- gluon-quark and gluon anti-quark
      elseif ((j == 0) .and. (k  /=  0)) then
c      write(6,*) 'as3vj_virtint: j,k,msqv=',j,k,fx1(j)*fx2(k)*msqv(j,k)
      m=0
      n=+1
      xmsqt=0._dp
      do cs=4,6
      xmsqt=xmsqt
     & +msq_struc(cs,g,n)*(AP(g,g,1)-AP(g,g,3)
     &                 +H1(g,g,q,cs,1)-H1(g,g,q,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +H2(q,q,g,cs,1)-H2(q,q,g,cs,3))*fx1(g)*fx2(k)
     & +(msq_struc(cs,g,n)*(AP(g,g,2)+AP(g,g,3)
     &                 +H1(g,g,q,cs,2)+H1(g,g,q,cs,3))
     &  +msq_struc(cs,-n,n)*(AP(a,g,2)+H1(a,g,q,cs,2)))*fx1z(g)/z*fx2(k)
     & +(msq_struc(cs,g,n)*(AP(q,q,2)+AP(q,q,3)
     &                 +H2(q,q,g,cs,2)+H2(q,q,g,cs,3))
     &  +msq_struc(cs,g,g)*(AP(g,q,2)+H2(g,q,g,cs,2)))*fx1(g)*fx2z(k)/z
      enddo
      do cs=1,3
      msq_aq=msq_struc(cs,-1,n)+msq_struc(cs,-2,n)+msq_struc(cs,-3,n)
     &      +msq_struc(cs,-4,n)+msq_struc(cs,-5,n)
      msq_qq=msq_struc(cs,+1,n)+msq_struc(cs,+2,n)+msq_struc(cs,+3,n)
     &      +msq_struc(cs,+4,n)+msq_struc(cs,+5,n)
      xmsqt=xmsqt
     & +(+msq_aq*(AP(a,g,2)+H1(a,g,q,cs,2))
     &   +msq_qq*(AP(q,g,2)+H1(q,g,q,cs,2)))*fx1z(g)/z*fx2(k)
     & +(+msq_struc(cs,g,g)*(AP(g,q,2)+H2(g,q,g,cs,2)))*fx1(g)*fx2z(k)/z
      enddo
      xmsqt=xmsqt
     & +msq_struc(iqr,n,-n)*(AP(a,g,2)+H1(a,g,q,iqr,2))*fx1z(g)/z*fx2(k)
      xmsq=xmsq+xmsqt
c      write(6,*) 'as3vj_virtint: j,k,  ct=',j,k,xmsqt

      endif

      elseif ((kcase==kW_2jet) .or. (kcase==kZ_2jet)) then
c--- SUM BY COLOUR STRUCTURES: W/Z + 2 jet only

      xmsq=xmsq+fx1(j)*fx2(k)*(
     & msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))
c      write(6,*) j,k,'-> msqv = ',fx1(j)*fx2(k)*(
c     & msqv(j,k)+msq_cs(0,j,k)+msq_cs(1,j,k)+msq_cs(2,j,k))
c      tmp=xmsq
!      print *,j,k,xmsq,msqv(j,k),msq_cs(0,j,k),msq_cs(1,j,k),msq_cs(2,j,k)
      if ((j > 0) .and. (k>0)) then
      do cs=0,2
      ics=cs
      xmsq=xmsq
     & +msq_cs(ics,j,k)*(AP(q,q,1)-AP(q,q,3)
     &                  +R1(q,q,q,cs,1)-R1(q,q,q,cs,3)
     &                  +AP(q,q,1)-AP(q,q,3)
     &                  +R2(q,q,q,cs,1)-R2(q,q,q,cs,3))*fx1(j)*fx2(k)
     & +(msq_cs(ics,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                   +R1(q,q,q,cs,2)+R1(q,q,q,cs,3))
     & + msq_cs(ics,g,k)*(AP(g,q,2)+R1(g,q,q,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(ics,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                   +R2(q,q,q,cs,2)+R2(q,q,q,cs,3))
     & + msq_cs(ics,j,g)*(AP(g,q,2)+R2(g,q,q,cs,2)))*fx1(j)*fx2z(k)/z
      enddo
      elseif ((j < 0) .and. (k<0)) then
      do cs=0,2
      ics=cs
      xmsq=xmsq
     & +msq_cs(ics,j,k)*(AP(a,a,1)-AP(a,a,3)
     &                 +R1(a,a,a,cs,1)-R1(a,a,a,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +R2(a,a,a,cs,1)-R2(a,a,a,cs,3))*fx1(j)*fx2(k)
     & +(msq_cs(ics,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R1(a,a,a,cs,2)+R1(a,a,a,cs,3))
     & + msq_cs(ics,g,k)*(AP(g,a,2)+R1(g,a,a,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(ics,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R2(a,a,a,cs,2)+R2(a,a,a,cs,3))
     & + msq_cs(ics,j,g)*(AP(g,a,2)+R2(g,a,a,cs,2)))*fx1(j)*fx2z(k)/z

      enddo
      elseif ((j > 0) .and. (k<0)) then
      iq=q
      ia=a
      do cs=0,2
      xmsq=xmsq
     & +msq_cs(cs,j,k)*(AP(q,q,1)-AP(q,q,3)
     &                 +R1(iq,iq,ia,cs,1)-R1(iq,iq,ia,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +R2(ia,ia,iq,cs,1)-R2(ia,ia,iq,cs,3)
     &                   )*fx1(j)*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                  +R1(iq,iq,ia,cs,2)+R1(iq,iq,ia,cs,3))
     & + msq_cs(cs,g,k)*(AP(g,q,2)+R1(g,q,a,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R2(ia,ia,iq,cs,2)+R2(ia,ia,iq,cs,3))
     & + msq_cs(cs,j,g)*(AP(g,a,2)+R2(g,a,q,cs,2)))*fx1(j)*fx2z(k)/z

      enddo
      elseif ((j < 0) .and. (k>0)) then
      iq=q
      ia=a
      do cs=0,2
      xmsq=xmsq
     & +msq_cs(cs,j,k)*(AP(a,a,1)-AP(a,a,3)
     &                 +R1(ia,ia,iq,cs,1)-R1(ia,ia,iq,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +R2(iq,iq,ia,cs,1)-R2(iq,iq,ia,cs,3)
     &                   )*fx1(j)*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R1(ia,ia,iq,cs,2)+R1(ia,ia,iq,cs,3))
     & + msq_cs(cs,g,k)*(AP(g,a,2)+R1(g,a,q,cs,2)))*fx1z(j)/z*fx2(k)
     & +(msq_cs(cs,j,k)*(AP(q,q,2)+AP(q,q,3)
     &                  +R2(iq,iq,ia,cs,2)+R2(iq,iq,ia,cs,3))
     & + msq_cs(cs,j,g)*(AP(g,q,2)+R2(g,q,a,cs,2)))*fx1(j)*fx2z(k)/z

      enddo
      elseif ((j == g) .and. (k == g)) then
      do cs=0,2
      msq_qg=msq_cs(cs,+5,g)+msq_cs(cs,+4,g)+msq_cs(cs,+3,g)
     &      +msq_cs(cs,+2,g)+msq_cs(cs,+1,g)
     &      +msq_cs(cs,-5,g)+msq_cs(cs,-4,g)+msq_cs(cs,-3,g)
     &      +msq_cs(cs,-2,g)+msq_cs(cs,-1,g)
      msq_gq=msq_cs(cs,g,+5)+msq_cs(cs,g,+4)+msq_cs(cs,g,+3)
     &      +msq_cs(cs,g,+2)+msq_cs(cs,g,+1)
     &      +msq_cs(cs,g,-5)+msq_cs(cs,g,-4)+msq_cs(cs,g,-3)
     &      +msq_cs(cs,g,-2)+msq_cs(cs,g,-1)
      xmsq=xmsq
     & +msq_cs(cs,g,g)*(AP(g,g,1)-AP(g,g,3)
     &                 +R1(g,g,g,cs,1)-R1(g,g,g,cs,3)
     &                 +AP(g,g,1)-AP(g,g,3)
     &                 +R2(g,g,g,cs,1)-R2(g,g,g,cs,3))*fx1(g)*fx2(g)
     & +(msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +R1(g,g,g,cs,3)+R1(g,g,g,cs,2))
     & + msq_qg*(AP(q,g,2)+R1(q,g,g,cs,2)))*fx1z(g)/z*fx2(g)
     & +(msq_cs(cs,g,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +R2(g,g,g,cs,3)+R2(g,g,g,cs,2))
     & + msq_gq*(AP(q,g,2)+R2(q,g,g,cs,2)))*fx1(g)*fx2z(g)/z
      enddo
      elseif ((j == g) .and. (k > 0)) then
c--- special case for W+bj - remove b-PDF contribution
      if ((kcase==kW_bjet) .and. (k  /=  5)) then
      xmsq=xmsq+(msq(-5,k)+msq(+5,k))*APqg_mass*fx1z(g)/z*fx2(k)
      else
      do cs=0,2
      msq_aq=msq_cs(cs,-1,k)+msq_cs(cs,-2,k)+msq_cs(cs,-3,k)
     &      +msq_cs(cs,-4,k)+msq_cs(cs,-5,k)
      msq_qq=msq_cs(cs,+1,k)+msq_cs(cs,+2,k)+msq_cs(cs,+3,k)
     &      +msq_cs(cs,+4,k)+msq_cs(cs,+5,k)
      xmsq=xmsq
     & +msq_cs(cs,g,k)*(AP(g,g,1)-AP(g,g,3)
     &                 +R1(g,g,q,cs,1)-R1(g,g,q,cs,3)
     &                 +AP(q,q,1)-AP(q,q,3)
     &                 +R2(q,q,g,cs,1)-R2(q,q,g,cs,3))*fx1(g)*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(g,g,2)+AP(g,g,3)
     &                  +R1(g,g,q,cs,2)+R1(g,g,q,cs,3))
     & + msq_aq*(AP(a,g,2)+R1(a,g,q,cs,2))
     & + msq_qq*(AP(q,g,2)+R1(q,g,q,cs,2)))*fx1z(g)/z*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(q,q,2)+AP(q,q,3)
     &                +R2(q,q,g,cs,2)+R2(q,q,g,cs,3))
     & + msq_cs(cs,g,g)*(AP(g,q,2)+R2(g,q,g,cs,2)))*fx1(g)*fx2z(k)/z

      enddo
      endif

      elseif ((j == g) .and. (k < 0)) then
c--- special case for W+bj - remove b-PDF contribution
      if ((kcase==kW_bjet) .and. (k  /=  -5)) then
      xmsq=xmsq+(msq(-5,k)+msq(+5,k))*APqg_mass*fx1z(g)/z*fx2(k)
      else
      do cs=0,2
      msq_qa=msq_cs(cs,+1,k)+msq_cs(cs,+2,k)+msq_cs(cs,+3,k)
     &      +msq_cs(cs,+4,k)+msq_cs(cs,+5,k)
      msq_aa=msq_cs(cs,-1,k)+msq_cs(cs,-2,k)+msq_cs(cs,-3,k)
     &      +msq_cs(cs,-4,k)+msq_cs(cs,-5,k)
      xmsq=xmsq
     & +msq_cs(cs,g,k)*(AP(g,g,1)-AP(g,g,3)
     &                 +R1(g,g,a,cs,1)-R1(g,g,a,cs,3)
     &                 +AP(a,a,1)-AP(a,a,3)
     &                 +R2(a,a,g,cs,1)-R2(a,a,g,cs,3))*fx1(g)*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(g,g,2)+AP(g,g,3)
     &                  +R1(g,g,a,cs,2)+R1(g,g,a,cs,3))
     & + msq_qa*(AP(q,g,2)+R1(q,g,a,cs,2))
     & + msq_aa*(AP(a,g,2)+R1(a,g,a,cs,2)))*fx1z(g)/z*fx2(k)
     & +(msq_cs(cs,g,k)*(AP(a,a,2)+AP(a,a,3)
     &                  +R2(a,a,g,cs,2)+R2(a,a,g,cs,3))
     & + msq_cs(cs,g,g)*(AP(g,a,2)+R2(g,a,g,cs,2)))*fx1(g)*fx2z(k)/z

      enddo
      endif

      elseif ((j > 0) .and. (k == g)) then
c--- special case for W+bj - remove b-PDF contribution
      if ((kcase==kW_bjet) .and. (j  /=  5)) then
      xmsq=xmsq+(msq(j,-5)+msq(j,+5))*APqg_mass*fx1(j)*fx2z(g)/z
      else
      do cs=0,2
      msq_qa=msq_cs(cs,j,-1)+msq_cs(cs,j,-2)+msq_cs(cs,j,-3)
     &      +msq_cs(cs,j,-4)+msq_cs(cs,j,-5)
      msq_qq=msq_cs(cs,j,+1)+msq_cs(cs,j,+2)+msq_cs(cs,j,+3)
     &      +msq_cs(cs,j,+4)+msq_cs(cs,j,+5)
       xmsq=xmsq
     &+ msq_cs(cs,j,g)*(AP(q,q,1)-AP(q,q,3)
     &                 +R1(q,q,g,cs,1)-R1(q,q,g,cs,3)
     &                 +R2(g,g,q,cs,1)-R2(g,g,q,cs,3)
     &                 +AP(g,g,1)-AP(g,g,3))*fx1(j)*fx2(g)
     &+(msq_cs(cs,j,g)*(AP(q,q,2)+AP(q,q,3)
     &                 +R1(q,q,g,cs,2)+R1(q,q,g,cs,3))
     &+ msq_cs(cs,g,g)*(AP(g,q,2)+R1(g,q,g,cs,2)))*fx1z(j)/z*fx2(g)
     &+(msq_cs(cs,j,g)*(AP(g,g,2)+AP(g,g,3)
     &                 +R2(g,g,q,cs,2)+R2(g,g,q,cs,3))
     &+ msq_qa*(AP(a,g,2)+R2(a,g,q,cs,2))
     &+ msq_qq*(AP(q,g,2)+R2(q,g,q,cs,2)))*fx1(j)*fx2z(g)/z

      enddo
      endif

      elseif ((j < 0) .and. (k == g)) then
c--- special case for W+bj - remove b-PDF contribution
      if ((kcase==kW_bjet) .and. (j  /=  -5)) then
      xmsq=xmsq+(msq(j,-5)+msq(j,+5))*APqg_mass*fx1(j)*fx2z(g)/z
      else
      do cs=0,2
      msq_aq=msq_cs(cs,j,+1)+msq_cs(cs,j,+2)+msq_cs(cs,j,+3)
     &      +msq_cs(cs,j,+4)+msq_cs(cs,j,+5)
      msq_aa=msq_cs(cs,j,-1)+msq_cs(cs,j,-2)+msq_cs(cs,j,-3)
     &      +msq_cs(cs,j,-4)+msq_cs(cs,j,-5)
       xmsq=xmsq
     & + msq_cs(cs,j,g)*(AP(a,a,1)-AP(a,a,3)
     &                  +R1(a,a,g,cs,1)-R1(a,a,g,cs,3)
     &                  +AP(g,g,1)-AP(g,g,3)
     &                  +R2(g,g,a,cs,1)-R2(g,g,a,cs,3))*fx1(j)*fx2(g)
     & +(msq_cs(cs,j,g)*(AP(a,a,2)+AP(a,a,3)
     &                  +R1(a,a,g,cs,2)+R1(a,a,g,cs,3))
     & + msq_cs(cs,g,g)*(AP(g,a,2)+R1(g,a,g,cs,2)))*fx1z(j)/z*fx2(g)
     & +(msq_cs(cs,j,g)*(AP(g,g,2)+AP(g,g,3)
     &                  +R2(g,g,a,cs,2)+R2(g,g,a,cs,3))
     & + msq_aq*(AP(q,g,2)+R2(q,g,a,cs,2))
     & + msq_aa*(AP(a,g,2)+R2(a,g,a,cs,2)))*fx1(j)*fx2z(g)/z

      enddo
      endif
      endif

c--- END 2 JET

      endif
c--- subtract off LO (we don't want it) for Wcs_ms case and for
c--- the comparison with C. Oleari's e+e- --> QQbg calculation
c      if ((kcase==kWcs_ms) .or. (runstring(1:5) == 'carlo')) then
c--- (MCFM_original)  if (kcase==kWcs_ms) then

c     if ((kcase==kWcs_ms) .or. (purevirt)) then
c       xmsq=xmsq-msq(j,k)*fx1(j)*fx2(k)
c     endif

      if (kewcorr /= knone) xmsq_noew=xmsq_noew+fx1(j)*fx2(k)*msq_noew(j,k)

      enddo
      enddo

!      print *,flux,xjac,pswt,xmsq,BrnRat
      as3vj_virtint=flux*xjac*pswt*xmsq/BrnRat

c--- code to check that epsilon poles cancel
      if (nshot == 1) then
        if ((xmsq == 0._dp) .or. ieee_is_nan(xmsq)) goto 999
        xmsq_old=xmsq
        nshot=nshot+1
        epinv=0._dp
        epinv2=0._dp
        goto 12
      elseif (nshot == 2) then
        nshot=nshot+1

        if ((abs(xmsq_old/xmsq-1._dp) > 1.e-7_dp) .or. ieee_is_nan(xmsq_old/xmsq)) then
!$omp master
          if (rank == 0) then
            write(6,*) 'epsilon fails to cancel'
            write(6,*) 'xmsq (epinv=large) = ',xmsq_old
            write(6,*) 'xmsq (epinv=zero ) = ',xmsq
            write(6,*) 'fractional difference',xmsq/xmsq_old-1d0
            call flush(6)
            !pause
          endif
!$omp end master
          colourchoice=rvcolourchoice
        else
!$omp master
          if (rank == 0) then
            write(6,*) 'Poles cancelled!'
c            write(6,*) 'xmsq (epinv=large) = ',xmsq_old
c            write(6,*) 'xmsq (epinv=zero ) = ',xmsq
c            write(6,*) 'fractional difference',xmsq/xmsq_old-1d0
            call flush(6)
          endif
!$omp end master
          colourchoice=rvcolourchoice
        endif
      endif

      call getptildejet(0,pjet)

      call dotem(nvec,pjet,s)

      if (ieee_is_nan(as3vj_virtint) .or. (.not. ieee_is_finite(as3vj_virtint))) then
         write(6,*) 'Found as3vj_virtint=',as3vj_virtint
         !write(6,*) 'random #s',r
         as3vj_virtint=zip
         goto 999
      endif

      val=as3vj_virtint*wgt
      val2=val**2

      if (abs(val) > wtmax) then
        wtmax=abs(val)
      endif

      ! QandGflag part is also affected by this, so have it before
c      if (includeTaucutgrid(0) .eqv. .false.) then
c          as3vj_virtint = 0._dp
c      endif

c--- handle special case of Qflag and Gflag
      if (QandGflag) then
        QandGint=QandGint+as3vj_virtint
        if ((Gflag) .and. (.not.(Qflag))) then
c--- go back for second pass (Qflag)
          Qflag=.true.
          Gflag=.false.
!          print*,'second pass: Qflag'
          goto 44
        else
c--- return both to .true. and assign value to as3vj_virtint (to return to VEGAS)
          Qflag=.true.
          Gflag=.true.
          as3vj_virtint=QandGint
        endif
      endif
      
C     Fill only if it's last iteration
      if (doFill.ne.0) then
         call hists_fill(p(3,:),p(4,:),as3vj_virtint*wgt)
!         call hists_fill(pjet(3,:),pjet(4,:),xmsq*wgt)
      endif
      
!      print *,'virtint',as3vj_virtint
      return

 999  continue
      as3vj_virtint=0._dp

c--- safety catch
      if (QandGflag) then
        Qflag=.true.
        Gflag=.true.
      endif

      return
      end

