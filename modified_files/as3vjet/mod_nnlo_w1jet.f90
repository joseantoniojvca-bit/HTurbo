!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

module nnlo_w1jet
    use types
    use SCET
    use SCET_Jet
    implicit none

    public :: lumxmsq_w1jet

    private

    private :: Wampqqbgsq
      
    contains

      subroutine lumxmsq_w1jet(p,xx,z1,z2,QB,order,xmsq,central)
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
      include 'beamtype.f'

      real(dp), intent(in) :: p(as3vj_mxpart,4), xx(2), z1, z2
      integer, intent(in) :: order
      real(dp), intent(out) :: xmsq
      logical, intent(in) :: central

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
       assemblejet,y12,y15,y25,Iijm(5)

      real(dp) :: qqbWg,qbqWg,qgWq,qbgWqb,gqbWqb,gqWq
      real(dp) :: Hqqb(2),Hqg(2),Hgqb(2),Hqbq(2),Hqbg(2),Hgq(2)
      real(dp) :: propsq

      real(dp) :: tauc, ptpure, puremass
      real(dp) :: twomass, pttwo, gettaucut

      tauc=gettaucut(twomass(3,4,p),pttwo(3,4,p))
!      if (dynamictau) then
!          tauc = getdynamictau(p, taucut)
!      else
!          tauc = taucut
!      endif

      call mcfm_spinoru(5,p,za,zb)
      propsq = (s(3,4)-wmass**2)**2+(wmass*wwidth)**2
      fac = gwsq**2*gsq*V/propsq
!      print *,propsq,fac
      call Wampqqbgsq(order,2,1,5,3,4,za,zb,qqbWg,Hqqb(1),Hqqb(2))
      call Wampqqbgsq(order,5,1,2,3,4,za,zb,qgWq,Hqg(1),Hqg(2))
      call Wampqqbgsq(order,2,5,1,3,4,za,zb,gqbWqb,Hgqb(1),Hgqb(2))

      call Wampqqbgsq(order,1,2,5,3,4,za,zb,qbqWg,Hqbq(1),Hqbq(2))
      call Wampqqbgsq(order,1,5,2,3,4,za,zb,qbgWqb,Hqbg(1),Hqbg(2))
      call Wampqqbgsq(order,5,2,1,3,4,za,zb,gqWq,Hgq(1),Hgq(2))

      qqbWg  = aveqq*fac*qqbWg
      qgWq   = aveqg*fac*qgWq
      gqbWqb = aveqg*fac*gqbWqb
      qbqWg  = aveqq*fac*qbqWg
      qbgWqb = aveqg*fac*qbgWqb
      gqWq   = aveqg*fac*gqWq


      if (order >= 0) then
!        call fdist(ih1,xx(1),facscale,beama0,1)
!        call fdist(ih2,xx(2),facscale,beamb0,2)
        call pdf_fdist(1,xx(1),facscale,beama0)
        call pdf_fdist(2,xx(2),facscale,beamb0)
      endif
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
!      print *,'Vsq',Vsq
      xmsq=zip
      do j=-nf,nf; do k=-nf,nf
          if ((j == 0) .and. (k == 0)) cycle
          if (j*k > 0) cycle
      
      if (j*k == 0) then ! initial state gluon -> jet is a quark
        if (j == 0) then
          if (k > 0) then
            hard(:)=Hgq(:)
          else
            hard(:)=Hgqb(:)
          endif
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1gq,soft2gq,jet1q,jet2q,hard)
        elseif (k == 0) then
          if (j > 0) then
            hard(:)=Hqg(:)
          else
            hard(:)=Hqbg(:)
          endif
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1qg,soft2qg,jet1q,jet2q,hard)
        endif
      else               ! initial state quarks -> jet is a gluon
          if (j > 0) then
            hard(:)=Hqqb(:)
          else
            hard(:)=Hqbq(:)
          endif
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1qa,soft2qa,jet1g,jet2g,hard)
      endif
!      print *
!      print *,tauc
!      print *,'beama0(j)    ',beama0(j)    
!      print *,'beamb0(k)    ',beamb0(k)
!      print *,'beama1(j,:)  ',beama1(j,:)
!      print *,'beamb1(k,:)  ',beamb1(k,:)
!      print *,'beama2(j,:)  ',beama2(j,:)
!      print *,'beamb2(k,:)  ',beamb2(k,:)
!      print *,'soft1gq      ',soft1gq
!      print *,'soft2gq      ',soft2gq
!      print *,'jet1q        ',jet1q
!      print *,'jet2q        ',jet2q
!      print *,'hard         ',hard          
!      print *,j,k,bit
      
      if     ((j > 0) .and. (k < 0)) then
        bit=bit*Vsq(j,k)*qqbWg
      elseif ((j < 0) .and. (k > 0)) then
        bit=bit*Vsq(j,k)*qbqWg
      elseif ((j > 0) .and. (k == 0)) then
        bit=bit* &
         (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWq
      elseif ((j < 0) .and. (k == 0)) then
        bit=bit* &
          (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqb
      elseif ((j == 0) .and. (k > 0)) then
        bit=bit* &
          (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWq
      elseif ((j == 0) .and. (k < 0)) then
        bit=bit* &
          (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqb
      else
        bit=zip
      endif

      xmsq=xmsq+bit

      enddo; enddo

!      print *,'doMultitaucut',doMultitaucut
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
          if ((j == 0) .and. (k == 0)) cycle
          if (j*k > 0) cycle
      
      if (j*k == 0) then ! initial state gluon -> jet is a quark
        if (j == 0) then
          if (k > 0) then
            hard(:)=Hgq(:)
          else
            hard(:)=Hgqb(:)
          endif
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1gq,soft2gq,jet1q,jet2q,hard)
        elseif (k == 0) then
          if (j > 0) then
            hard(:)=Hqg(:)
          else
            hard(:)=Hqbg(:)
          endif
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1qg,soft2qg,jet1q,jet2q,hard)
        endif
      else               ! initial state quarks -> jet is a gluon
          if (j > 0) then
            hard(:)=Hqqb(:)
          else
            hard(:)=Hqbq(:)
          endif
          bit=assemblejet(order,tauc, &
           beama0(j),beamb0(k),beama1(j,:),beamb1(k,:), &
           beama2(j,:),beamb2(k,:),soft1qa,soft2qa,jet1g,jet2g,hard)
      endif
      
      if     ((j > 0) .and. (k < 0)) then
        bit=bit*Vsq(j,k)*qqbWg
      elseif ((j < 0) .and. (k > 0)) then
        bit=bit*Vsq(j,k)*qbqWg
      elseif ((j > 0) .and. (k == 0)) then
        bit=bit* &
         (Vsq(j,-1)+Vsq(j,-2)+Vsq(j,-3)+Vsq(j,-4)+Vsq(j,-5))*qgWq
      elseif ((j < 0) .and. (k == 0)) then
        bit=bit* &
          (Vsq(j,+1)+Vsq(j,+2)+Vsq(j,+3)+Vsq(j,+4)+Vsq(j,+5))*qbgWqb
      elseif ((j == 0) .and. (k > 0)) then
        bit=bit* &
          (Vsq(-1,k)+Vsq(-2,k)+Vsq(-3,k)+Vsq(-4,k)+Vsq(-5,k))*gqWq
      elseif ((j == 0) .and. (k < 0)) then
        bit=bit* &
          (Vsq(+1,k)+Vsq(+2,k)+Vsq(+3,k)+Vsq(+4,k)+Vsq(+5,k))*gqbWqb
      else
        bit=zip
      endif

      scetreweight(m)=scetreweight(m)+bit

      enddo; enddo
!!!! END COPIED BLOCK FOR SCALEREWEIGHT 


              enddo
              scetreweight(:) = scetreweight(:) / xmsq
          endif
      endif
      
      end subroutine lumxmsq_w1jet

    subroutine Wampqqbgsq(order,p1,p2,p3,p4,p5,za,zb,msq0,hard1,hard2)
        use types
        use nnlo_z1jet_hfun
        use constants
!--- returns msq0, hard1, hard2 which can be used to reconstruct the
!--- complete hard function for ab -> W+c using:
!---  H = msq0 * (1 + [as/2/pi]*hard1 + [as/2/pi]^2*hard2)
    implicit none
    include 'nf.f'
    include 'mxpart.f'
    include 'zprods_decl.f'
    include 'sprods_com.f'
    include 'hpls.f'
    include 'scale.f'
    include 'as3vjet_ewcouple.f'
    include 'as3vjet_qcdcouple.f'
    include 'as3vjet_masses.f'
!     This routine calculates the amplitude for a right-handed quark
!     0--> q**+(1)+qb**1(2)+g**+(3)+l**-(4)+lb**+(5)
!     according to Eq.22 of 1309.3245v3
!     index is helicity of outgoing gluon line
    real(dp):: W1Lampqqbgsq,uu,vv,xlf,NFV,hard1,hard2
    integer:: order,p1,p2,p3,p4,p5,q1,q2,q3,iperm,region,k
    complex(dp):: al,be,ga,de
    complex(dp):: amp(0:2,2),ampqqbgll,C0,lnrat
    real(dp):: hn,s12,s13,s23,qsq,res0,res1,res2,msq0,msq1,msq2, &
    evolve1Lqqb,evolve2Lqqb,c1
    complex(dp):: alp0,bet0,gam0,alp1,bet1,gam1,alp2,bet2,gam2
    real(dp):: alR(0:2),alI(0:2)
    real(dp):: beR(0:2),beI(0:2)
    real(dp):: gaR(0:2),gaI(0:2)
    real(dp):: LgamHqq0,LgamHqq1,Lrat

    real(dp) :: alpha2a_c0re,alpha2a_c0im,alpha2a_c1re,alpha2a_c1im,alpha2a_c2re,alpha2a_c2im
    real(dp) :: beta2a_c0re,beta2a_c0im,beta2a_c1re,beta2a_c1im,beta2a_c2re,beta2a_c2im
    real(dp) :: gamma2a_c0re,gamma2a_c0im,gamma2a_c1re,gamma2a_c1im,gamma2a_c2re,gamma2a_c2im
    real(dp) :: alpha3a_c0re,alpha3a_c0im,alpha3a_c1re,alpha3a_c1im,alpha3a_c2re,alpha3a_c2im
    real(dp) :: beta3a_c0re,beta3a_c0im,beta3a_c1re,beta3a_c1im,beta3a_c2re,beta3a_c2im
    real(dp) :: gamma3a_c0re,gamma3a_c0im,gamma3a_c1re,gamma3a_c1im,gamma3a_c2re,gamma3a_c2im
    real(dp) :: alpha4a_c0re,alpha4a_c0im,alpha4a_c1re,alpha4a_c1im,alpha4a_c2re,alpha4a_c2im
    real(dp) :: beta4a_c0re,beta4a_c0im,beta4a_c1re,beta4a_c1im,beta4a_c2re,beta4a_c2im
    real(dp) :: gamma4a_c0re,gamma4a_c0im,gamma4a_c1re,gamma4a_c1im,gamma4a_c2re,gamma4a_c2im

    logical, parameter:: peter=.false. ! use original PeTeR functions in C++
    
      real(dp), parameter :: zip = 0._dp
      real(dp), parameter :: three = 3._dp, four = 4._dp, eight = 8._dp
      real(dp), parameter :: one = 1._dp, two = 2._dp, six = 6._dp
      real(dp), parameter :: be0 = 11/three*CA-4/three*TR*NF
      real(dp), parameter :: be1= 34/three*CA**2-(20/three*CA+4*CF)*TR*NF

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
    NFV=zip

    if (1 == 2) then
    !---- ! Check Becher et al. result at a single PS point
        if ((p1 == 2) .AND. (p2 == 1)) then
            s(p1,p2)=1._dp
            s(p1,p3)=-0.4_dp
            s(p4,p5)=0.1_dp**2
            s(p2,p3)=s(p4,p5)-s(p1,p2)-s(p1,p3)
            s(p2,p1)=s(p1,p2)
            s(p3,p1)=s(p1,p3)
            s(p3,p2)=s(p2,p3)
        endif
        if ((p1 == 5) .AND. (p2 == 1)) then
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

    res0=zip
    res1=zip
    res2=zip
    do iperm=1,2
        if (iperm == 1) then
            q1=p1
            q2=p2
            q3=p3
        else
            q1=p2
            q2=p1
            q3=p3
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
            hard1 = 0._dp
            hard2 = 0._dp
            return
        endif
    !      write(6,*) 'iperm,region',iperm,region

        if (region == 2) then
            uu=-s(q1,q3)/s(q1,q2)
            vv=+s(p4,p5)/s(q1,q2)
        elseif (region == 3) then
            uu=-s(q2,q3)/s(q1,q3)
            vv=+s(p4,p5)/s(q1,q3)
        elseif (region == 4) then
            uu=-s(q1,q3)/s(q2,q3)
            vv=+s(p4,p5)/s(q2,q3)
        endif
                    
        if (order > 0) then
        !--- fill arrays for 2DHPLs
           call tdhpl(uu,vv,4,G1,G2,G3,G4,H1,H2,H3,H4)
           call sethplh(H1,H2,H3,H4)
           call sethplg(G1,G2,G3,G4)
        else
        !--- fill with dummy values
            G1=zip
            G2=zip
            G3=zip
            G4=zip
            H1=zip
            H2=zip
            H3=zip
            H4=zip
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
            if (order >= 1) then
                alR(1)=Alpha_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                alI(1)=Alpha_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                beR(1)=Beta_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                beI(1)=Beta_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                gaR(1)=Gamma_2a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                gaI(1)=Gamma_2a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            endif
            if (order >= 2) then
                alR(2)=Alpha_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                alI(2)=Alpha_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                beR(2)=Beta_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                beI(2)=Beta_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                gaR(2)=Gamma_2a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                gaI(2)=Gamma_2a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            endif

        !--- Coefficients for Region 3
        elseif (region == 3) then
            alR(0)=Alpha_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            alI(0)=Alpha_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            beR(0)=Beta_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            beI(0)=Beta_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            gaR(0)=Gamma_3a0re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            gaI(0)=Gamma_3a0im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            if (order >= 1) then
                alR(1)=Alpha_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                alI(1)=Alpha_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                beR(1)=Beta_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                beI(1)=Beta_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                gaR(1)=Gamma_3a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                gaI(1)=Gamma_3a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            endif
            if (order >= 2) then
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
            if (order >= 1) then
                alR(1)=Alpha_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                alI(1)=Alpha_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                beR(1)=Beta_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                beI(1)=Beta_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                gaR(1)=Gamma_4a1re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                gaI(1)=Gamma_4a1im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            endif
            if (order >= 2) then
                alR(2)=Alpha_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                alI(2)=Alpha_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                beR(2)=Beta_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                beI(2)=Beta_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                gaR(2)=Gamma_4a2re(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
                gaI(2)=Gamma_4a2im(uu,vv,H1,H2,H3,H4,G1,G2,G3,G4,xlf,NFV)
            endif
                        
            !print *
            !print *,'Region 4'
            !print *,uu,vv,xlf,NFV
            !print *,'alpha4a_c0re', alR(0),alpha4a_c0re(uu,vv)
            !print *,'alpha4a_c0im', alI(0),alpha4a_c0im(uu,vv)
            !print *,'beta4a_c0re',  beR(0),beta4a_c0re(uu,vv)
            !print *,'beta4a_c0im',  beI(0),beta4a_c0im(uu,vv)
            !print *,'gamma4a_c0re', gaR(0),gamma4a_c0re(uu,vv)
            !print *,'gamma4a_c0im', gaI(0),gamma4a_c0im(uu,vv)
            !print *,'alpha4a_c1re', alR(1),alpha4a_c1re(uu,vv)
            !print *,'alpha4a_c1im', alI(1),alpha4a_c1im(uu,vv)
            !print *,'beta4a_c1re',  beR(1),beta4a_c1re(uu,vv)
            !print *,'beta4a_c1im',  beI(1),beta4a_c1im(uu,vv)
            !print *,'gamma4a_c1re', gaR(1),gamma4a_c1re(uu,vv)
            !print *,'gamma4a_c1im', gaI(1),gamma4a_c1im(uu,vv)
            !print *,'alpha4a_c2re', alR(2),alpha4a_c2re(uu,vv)
            !print *,'alpha4a_c2im', alI(2),alpha4a_c2im(uu,vv)
            !print *,'beta4a_c2re',  beR(2),beta4a_c2re(uu,vv)
            !print *,'beta4a_c2im',  beI(2),beta4a_c2im(uu,vv)
            !print *,'gamma4a_c2re', gaR(2),gamma4a_c2re(uu,vv)
            !print *,'gamma4a_c2im', gaI(2),gamma4a_c2im(uu,vv)
           
        else
            write(6,*) 'Wampqqbgsq: region should be 2,3 or 4: ',region
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
        do k=0,order
            al=cmplx(alR(k),alI(k),dp)
            be=cmplx(beR(k),beI(k),dp)
            ga=cmplx(gaR(k),gaI(k),dp)
            de=s(p1,p2)*(al-be-ga)/(2._dp*(s(p1,p2)+s(p2,p3)+s(p3,p1)))

            if (iperm == 1) then
            ! H gluon
                amp(k,1)=+ampqqbgll(p1,p2,p3,p5,p4,al,be,ga,de,zb,za)
            endif
            if (iperm == 2) then
            ! H gluon
                amp(k,2)=-ampqqbgll(p2,p1,p3,p4,p5,al,be,ga,de,za,zb)
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
        +two*hn(s(p1,p2),s(q2,q3),s(q1,q3),s(p4,p5), &
        alp0,bet0,gam0,alp0,bet0,gam0)
             
        res1=res1 &
        +two*hn(s(p1,p2),s(q2,q3),s(q1,q3),s(p4,p5), &
        alp0,bet0,gam0,alp1,bet1,gam1)
             
        res2=res2 &
        +two*hn(s(p1,p2),s(q2,q3),s(q1,q3),s(p4,p5), &
        alp0,bet0,gam0,alp2,bet2,gam2) &
        +hn(s(p1,p2),s(q2,q3),s(q1,q3),s(p4,p5), &
        alp1,bet1,gam1,alp1,bet1,gam1)
             
    enddo ! end of iperm loop
          
!---- This bit of code for testing squared MEs against Becher et al.
    msq0=res0
!        write(6,*) 'LO not decayed',msq0*gwsq*gsq*V*aveqq/two
    if (order > 0) then
        msq1=two*res1
    !--- RG running to restore scale dependence from musq=s(p4,p5)
        msq1=msq1 &
        +(evolve1Lqqb(musq)-evolve1Lqqb(s(p4,p5)))*res0
    !        write(6,*) '1L/LO/fourpi',msq1/msq0/fourpi
    !        write(6,*) '1L not decayed',msq1*gwsq*gsq*V*aveqq/two*ason4pi
    endif
    if (order > 1) then
        msq2=two*res2
    !--- RG running to restore scale dependence from musq=s(p4,p5)
        c1=two*res1/res0-evolve1Lqqb(s(p4,p5))
        msq2=msq2 &
        +(evolve2Lqqb(musq,c1) &
        -evolve2Lqqb(s(p4,p5),c1))*res0
    !        write(6,*) '2L/LO/fourpi**2',msq2/msq0/fourpi**2
    !        write(6,*) '2L not decayed',msq2*gwsq*gsq*V*aveqq/two*ason4pi**2
    endif
          
!      write(6,*) 'amp(0,1)',amp(0,1)
!      write(6,*) 'amp(0,2)',amp(0,2)
!      write(6,*) 'amp(1,1)',amp(1,1)
!      write(6,*) 'amp(1,2)',amp(1,2)

    msq0=( &
    +real(amp(0,1)*conjg(amp(0,1)),dp) &
    +real(amp(0,2)*conjg(amp(0,2)),dp))
    if (order > 0) then
        msq1=two*( &
        +real(amp(0,1)*conjg(amp(1,1)),dp) &
        +real(amp(0,2)*conjg(amp(1,2)),dp))
    !--- RG running to restore scale dependence from musq=s(p4,p5)
        msq1=msq1 &
        +(evolve1Lqqb(musq)-evolve1Lqqb(s(p4,p5)))*msq0
    else
        msq1=zip
    endif
    if (order > 1) then
        msq2=two*( &
        +real(amp(0,1)*conjg(amp(2,1)),dp) &
        +real(amp(0,2)*conjg(amp(2,2)),dp)) &
        +real(amp(1,1)*conjg(amp(1,1)),dp) &
        +real(amp(1,2)*conjg(amp(1,2)),dp)
    !--- RG running to restore scale dependence from musq=s(p4,p5)
        c1=msq1/msq0-evolve1Lqqb(musq)
        msq2=msq2 &
        +(evolve2Lqqb(musq,c1) &
        -evolve2Lqqb(s(p4,p5),c1))*msq0
    else
        msq2=zip
    endif

!--- convert from coefficient of [as/4/pi]^n to coefficient of [as/2/pi]^n
    msq1=msq1/two
    msq2=msq2/four

!--- normalize coefficients to LO
    hard1=msq1/msq0
    hard2=msq2/msq0
          
    return
    end subroutine Wampqqbgsq

end module
