!-----------------------------------------------------------------------
! Written by D.Tommasini: Mellin of H2gg and H2gq
!-----------------------------------------------------------------------

      COMPLEX*16 FUNCTION MellinH2gq(ZN)
      IMPLICIT none
      COMPLEX*16 ZN,ZN2,ZN3,ZN4,ZN5,ZN6,ZNp1,ZNm1,ZNp2,ZNp3
      COMPLEX*16 ZNb,plinomZN,plinomZNm, temp
      COMPLEX*16 PG0ZN,PG1ZN,PG1ZNmez,PG1ZNp1mez,PG2ZN,PG2ZNmezp1
      COMPLEX*16 PG2ZNp1mez,PG0ZNmezp1,PG0ZNp1mez
      COMPLEX*16 DACG1ZN,DACG1ZNp1,ACG4ZN,ACG4ZNp1,ACG4ZNp2,DACG1ZNp2
      COMPLEX*16 DACG1,ACG4
      Real*8 GE,GE2,GE3,ZET3,Pi,Pi2,DL
      Real*8 CA,CF,nf
      CA=3.0D0
      CF=4.0D0/3.0D0
      nf=5.0D0
      GE   = 0.57721566490153D0
      GE2  =  GE*GE
      GE3  =  GE2*GE
      ZET3 = 1.20205690315959428540D0
      Pi = 3.141592653589793238462643D0
      Pi2 = Pi*Pi
      DL = DLOG(2.0D0)

      ZN2=ZN*ZN
      ZN3=ZN2*ZN
      ZN4=ZN3*ZN
      ZN5=ZN4*ZN
      ZN6=ZN5*ZN
      ZNp1 = ZN+(1.0D0,0.0D0)
      ZNm1 = ZN+(-1.0D0,0.0D0)
      ZNp2 = ZN+(2.0D0,0.0D0)
      ZNp3 = ZN+(3.0D0,0.0D0)
      ZNb = (-1.0D0 + ZN2)
      plinomZN = 2.0D0 + ZN + ZN2
      plinomZNm=-2.0D0 + ZN + ZN2
!      write(*,*) ZN,ZN2,zn3,zn4,zn5,zn6,znp1,znm1,znp2,znp3,znb
!----------------------------------------------------------------------
      CALL PSI0(ZN,PG0ZN)                        !PolyGamma[0, ZN]
      temp=ZN/2.0D0+1.0D0
      CALL PSI0(temp,PG0ZNmezp1)                 !PolyGamma[0, 1 + ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI0(temp,PG0ZNp1mez)                 !PolyGamma[0, ZNp1/2]
      CALL PSI1(ZN,PG1ZN)                        !PolyGamma[1, ZN]
      temp=ZN/2.0D0
      CALL PSI1(temp,PG1ZNmez)                   !PolyGamma[1, ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI1(temp,PG1ZNp1mez)                 !PolyGamma[1, ZNp1/2]
      CALL PSI2(ZN,PG2ZN)                        !PolyGamma[2, ZN]
      temp=ZN/2.0D0+1.0D0
      CALL PSI2(temp,PG2ZNmezp1)                 !PolyGamma[2, 1 + ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI2(temp,PG2ZNp1mez)                 !PolyGamma[2, ZNp1/2]
!--------------------------------------------------------------------
! NB:The Blumlein's convention to calculate the Mellin transform is z^N (not z^(N-1) )
!--------------------------------------------------------------------
      DACG1ZN=DACG1(ZNm1)
      DACG1ZNp1=DACG1(ZN)
      DACG1ZNp2=DACG1(ZNp1)
      ACG4ZN=ACG4(ZNm1)
      ACG4ZNp1=ACG4(ZN)
      ACG4ZNp2=ACG4(ZNp1)

!      write (*,*) 'H2gqZN',ZN

      MellinH2gq=(CF*((4.0D0*ZN*nf*(72.0D0 + ZN*ZNp1*(-12.0D0 +
     # 18.0D0*GE2*ZN*ZNp1*plinomZN +
     # 4.0D0*ZN*(-13.0D0 + ZN*(144.0D0 + ZN*(74.0D0 + 13.0D0*ZN))) -
     # 12.0D0*GE*
     # (-6.0D0+ZN*(7.0D0 + ZN*(27.0D0 +2.0D0*ZN*(5.0D0+ZN))))+3.0D0*
     # ZN*ZNp1*plinomZN*
     #    Pi2) + 6.0D0*ZN*ZNp1*
     #  (2.0D0*(6.0D0 + ZN*(-7.0D0 + 3.0D0*GE*ZNp1*plinomZN -
     #       ZN*(27.0D0 + 2.0D0*ZN*(5.0D0 + ZN))))*PG0ZN +
     #   3.0D0*ZN*ZNp1*plinomZN*PG0ZN**2 -
     #   3.0D0*ZN*ZNp1*plinomZN*PG1ZN)))/(ZNm1*ZNp1**3) +
     #    (CA*(8.0D0*ZNm1**2*(5184.0D0 +
     #  ZN*(24048.0D0 + ZN*(36744.0D0 + ZN*(19532.0D0+ZN*(20342.0D0+
     # ZN*(57126.0D0 +
     # ZN*(48911.0D0 + ZN*(-18621.0D0 +ZN*(-54864.0D0+ZN*(-34946.0D0+
     # ZN*(-9377.0D0 +
     #                     ZN*(-807.0D0 + 40.0D0*ZN)))))))))))) -
     # 3.0D0*ZN*ZNp1*plinomZNm**2*(576.0D0*ZN3*ZNp1**3*plinomZNm*
     #    DACG1ZN + 576.0D0*ZNm1**2*ZN2*ZNp1**3*ZNp2*DACG1ZNp1 -
     #   576.0D0*PG0ZN + ZN*(288.0D0*ZN2*ZNp2*ZNb**2*DACG1ZNp2 +
     #     ZNp1*(576.0D0*ZN2*ZNp1**2*plinomZNm*ACG4ZN +
     #       576.0D0*ZN*ZNp2*ZNb**2*ACG4ZNp1 + 576.0D0*ZN2*ACG4ZNp2 -
     #    288.0D0*ZN3*ACG4ZNp2 - 864.0D0*ZN4*ACG4ZNp2 +288.0D0*ZN5*
     # ACG4ZNp2+
     #       288.0D0*ZN6*ACG4ZNp2 - 8.0D0*(240.0D0 + ZN*(-484.0D0 +
     # ZN*(-721.0D0 +
     #        ZN*(211.0D0 + ZN*(462.0D0 + ZN*(213.0D0 + 43.0D0*
     # ZN))))))*PG0ZN +
     #       12.0D0*ZN*ZNp1*ZNp2*(-28.0D0 + ZN*(29.0D0 + ZN*(18.0D0 +
     # 5.0D0*ZN)))*
     #        PG0ZN**2 - 24.0D0*ZNm1*ZN*ZNp1*ZNp2*plinomZN*
     #     PG0ZN**3 - 36.0D0*ZNp2*(-4.0D0 + ZN2*(15.0D0 +ZN*(3.0D0 +
     # ZN + ZN2)))*
     #        PG1ZNmez - 12.0D0*(-48.0D0 + ZN*ZNp1*(-248.0D0 +
     # ZN*(-30.0D0 + ZN*(17.0D0 + ZN*(52.0D0 + 17.0D0*ZN)))))*PG1ZN -
     #       72.0D0*ZNm1*ZN*ZNp1*ZNp2*plinomZN*PG0ZN*
     # PG1ZN + 36.0D0*ZNp2*(-4.0D0 + ZN2*(15.0D0 + ZN*(3.0D0 +ZN+
     # ZN2)))*
     #        PG1ZNp1mez + 9.0D0*ZNm1*ZN*ZNp1*ZNp2*
     #        plinomZN*PG2ZNmezp1 + 120.0D0*ZNm1*ZN*ZNp1*
     #        ZNp2*plinomZN*PG2ZN - 9.0D0*ZNm1*ZN*ZNp1*
     #        ZNp2*plinomZN*PG2ZNp1mez)))))/
     # (ZNp2**3*ZNb**4) +
     # (6.0D0*CA*ZN*(2.0D0*GE*(72.0D0 + ZN*(204.0D0 + ZN*(-346.0D0 +
     # 6.0D0*GE2*ZNp1**2*
     #        (-2.0D0 + ZN + ZN3) - 3.0D0*GE*ZNp1**2*
     #        (-28.0D0 + ZN*(29.0D0 + ZN*(18.0D0 + 5.0D0*ZN))) +
     #    2.0D0*ZN*(-516.0D0 + ZN*(3.0D0 + ZN*(335.0D0+ZN*(170.0D0+
     # 43.0D0*ZN))))))) +
     # ZN*ZNp1**2*(-24.0D0 + ZN*(76.0D0 + 6.0D0*GE*(-2.0D0 + ZN + ZN3) +
     #     ZN*(85.0D0 + ZN*(-174.0D0 + 61.0D0*ZN))))*Pi2 + 6.0D0*ZN2*
     # ZNp1**2*
     #  ((-2.0D0 + ZN + ZN3)*Pi2*PG0ZNmezp1 +
     #   (2.0D0*GE*(28.0D0 + 3.0D0*GE*(-2.0D0 + ZN + ZN3) -
     # ZN*(29.0D0 + ZN*(18.0D0 + 5.0D0*ZN))) + (-2.0D0+ ZN+ZN3)*Pi2)*
     # PG0ZN +
     #   6.0D0*GE*(-2.0D0 + ZN + ZN3)*PG0ZN**2 -
     #   (-2.0D0 + ZN + ZN3)*(Pi2*PG0ZNp1mez - 6.0D0*GE*
     #      PG1ZN - 46.0D0*ZET3))))/(ZNm1**2*ZNp1**3) -
     # (18.0D0*CF*(24.0D0 - 2.0D0*ZN*(-2.0D0*GE3*ZN2*ZNp1**3*plinomZN -
     # 6.0D0*GE*ZNp1*(4.0D0 + ZN2*ZNp1*(1.0D0 + ZN*(36.0D0 +5.0D0*ZN*
     # (4.0D0+ZN)))) +
     # 3.0D0*GE2*ZN*ZNp1**2*(-4.0D0 + ZN*(4.0D0 + ZN*(17.0D0 + ZN*
     # (6.0D0+ZN)))) +
     # 3.0D0*(8.0D0 + ZN*(3.0D0 + ZN*(-99.0D0 + ZN*(-107.0D0 + ZN*
     # (63.0D0 + ZN*(113.0D0 +
     # (26.0D0 - 3.0D0*ZN)*
     #                ZN))))))) + ZN2*ZNp1**2*
     # (4.0D0 + ZN*(-4.0D0 + 2.0D0*GE*ZNp1*plinomZN - ZN*(17.0D0+ ZN*
     # (6.0D0+ ZN))))*
     #  Pi2 + 2.0D0*ZN*ZNp1*(3.0D0*ZN*ZNp1*
     #    (4.0D0 + ZN*(-4.0D0 + 2.0D0*GE*ZNp1*plinomZN -
     # ZN*(17.0D0 + ZN*(6.0D0 + ZN))))*PG0ZN**2 + 2.0D0*ZN2*ZNp1**2*
     #    plinomZN*PG0ZN**3 + PG0ZN*
     #    (24.0D0 + ZN*ZNp1*(6.0D0*GE2*ZN*ZNp1*plinomZN +
     #       6.0D0*ZN*(1.0D0 + ZN*(36.0D0 + 5.0D0*ZN*(4.0D0 + ZN))) -
     # 6.0D0*GE*
     #  (-4.0D0 + ZN*(4.0D0 + ZN*(17.0D0 + ZN*(6.0D0 + ZN)))) + ZN*
     # ZNp1*plinomZN*
     #        Pi2) - 6.0D0*ZN2*ZNp1**2*plinomZN*PG1ZN) +
     # ZN*ZNp1*(-3.0D0*(4.0D0 + ZN*(-4.0D0 + 2.0D0*GE*ZNp1*plinomZN -
     #         ZN*(17.0D0 + ZN*(6.0D0 + ZN))))*PG1ZN + 2.0D0*ZN*ZNp1*
     #      plinomZN*(PG2ZN + 2.0D0*ZET3)))))/
     # (ZNm1*ZNp1**4)))/(1728.0D0*ZN4)
!      write (*,*) 'H2gqZN',MellinH2gq,ZN
      RETURN
      END

      COMPLEX*16 FUNCTION MellinH2gg(ZN)
      IMPLICIT none
      COMPLEX*16 ZN,ZN2,ZN3,ZN4,ZN5,ZN6,ZNp1,ZNm1,ZNp2,ZNp3,ZNm2
      COMPLEX*16 ZNb,plinomZN,plinomZNm, temp
      COMPLEX*16 PG0ZNm1,PG0ZNmez,PG1ZNm1,PG0ZNp1,PG0ZNp2mez
      COMPLEX*16 PG1ZNmezp1,PG1ZNp1,PG0ZNp2,PG1ZNp2,PG0ZNp3,PG1ZNp3
      COMPLEX*16 PG2ZNm1,PG2ZNmez,PG2ZNp3,PG3ZNmezp1,PG3ZNm1,PG3ZN
      COMPLEX*16 PG3ZNp1mez,PG0ZNm1mez
      COMPLEX*16 PG0ZN,PG1ZN,PG1ZNmez,PG1ZNp1mez,PG2ZN,PG2ZNmezp1
      COMPLEX*16 PG2ZNp1mez,PG0ZNmezp1,PG0ZNp1mez
      COMPLEX*16 DDACG1ZNm1,DACG2ZNm1,DACG2ZN,DACG2ZNp1,DACG2ZNp2
      COMPLEX*16 DACG2ZNp3,DACG4ZNm1,ACG1ZNp1,ACG1ZNp2,ACG1ZNp3
      COMPLEX*16 ACG13ZNm1,ACG13ZN,ACG13ZNp1,ACG13ZNp2,ACG13ZNp3
      COMPLEX*16 ACG20ZNm1,ACG22ZNm1,ACG6ZNm1,ACG7ZNm1,ACG9ZNm1
      COMPLEX*16 DDACG1,DACG2,DACG4,ACG1,ACG13,ACG20,ACG22,ACG5ZNm1
      COMPLEX*16 ACG1ZNm1
      COMPLEX*16 ACG5,ACG6,ACG7,ACG9
      Real*8 GE,GE2,GE3,ZET3,Pi,Pi2,DL
      Real*8 CA,CF,nf,LT
      common/LT/LT
      CA=3.0D0
      CF=4.0D0/3.0D0
      nf=5.0D0
      GE   = 0.57721566490153
      GE2  =  GE*GE
      GE3  =  GE2*GE
      ZET3 = 1.20205690315959428540
      Pi = 3.141592653589793238462643
      Pi2 = Pi*Pi
      DL = DLOG(2.0D0)
!      write (*,*) 'H2ggZN',ZN
      ZN2=ZN*ZN
      ZN3=ZN2*ZN
      ZN4=ZN3*ZN
      ZN5=ZN4*ZN
      ZN6=ZN5*ZN
      ZNp1 = ZN+(1.0D0,0.0D0)
      ZNm1 = ZN+(-1.0D0,0.0D0)
      ZNp2 = ZN+(2.0D0,0.0D0)
      ZNm2 = ZN+(-2.0D0,0.0D0)
      ZNp3 = ZN+(3.0D0,0.0D0)
      ZNb = (-1.0D0 + ZN2)
      plinomZN = 2.0D0 + ZN + ZN2
      plinomZNm=-2.0D0 + ZN + ZN2
!      write(*,*) zn,zn2,zn3,zn4,zn5,zn6,znp1,znm1,znp2,znp3,znb
!------------------------------------------------------------------
      CALL PSI0(ZN,PG0ZN)                        !PolyGamma[0, ZN]
      temp=ZN/2.0D0+1.0D0
      CALL PSI0(temp,PG0ZNmezp1)                 !PolyGamma[0, 1 + ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI0(temp,PG0ZNp1mez)                 !PolyGamma[0, ZNp1/2]
      CALL PSI1(ZN,PG1ZN)                        !PolyGamma[1, ZN]
      temp=ZN/2.0D0
      CALL PSI1(temp,PG1ZNmez)                   !PolyGamma[1, ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI1(temp,PG1ZNp1mez)                 !PolyGamma[1, ZNp1/2]
      CALL PSI2(ZN,PG2ZN)                        !PolyGamma[2, ZN]
      temp=ZN/2.0D0+1.0D0
      CALL PSI2(temp,PG2ZNmezp1)                 !PolyGamma[2, 1 + ZN/2]
      temp=ZNp1/2.0D0
      CALL PSI2(temp,PG2ZNp1mez)                 !PolyGamma[2, ZNp1/2]
      temp=ZNm1
      CALL PSI0(temp,PG0ZNm1)                    !PolyGamma[0, -1 + ZN]
      temp=ZN/2.0D0
      CALL PSI0(temp,PG0ZNmez)                   !PolyGamma[0, ZN/2]
      temp=ZNm1
      CALL PSI1(temp,PG1ZNm1)                    !PolyGamma[1, -1 + ZN]
      temp=ZNp1
      CALL PSI0(temp,PG0ZNp1)                    !PolyGamma[0, 1 + ZN]
      temp=ZNp2/2.0D0
      CALL PSI0(temp,PG0ZNp2mez)                 !PolyGamma[0, ZNp2/2]
      temp=ZN/2.0D0+1.0D0
      CALL PSI1(temp,PG1ZNmezp1)                 !PolyGamma[1, 1 + ZN/2]
      CALL PSI1(ZNp1,PG1ZNp1)                    !PolyGamma[1, 1 + ZN]
      CALL PSI0(ZNp2,PG0ZNp2)                    !PolyGamma[0, 2 + ZN]
      CALL PSI1(ZNp2,PG1ZNp2)                    !PolyGamma[1, 2 + ZN]
      CALL PSI0(ZNp3,PG0ZNp3)                    !PolyGamma[0, 3 + ZN]
      CALL PSI1(ZNp3,PG1ZNp3)                    !PolyGamma[1, 3 + ZN]
      CALL PSI2(ZNm1,PG2ZNm1)                    !PolyGamma[2, -1 + ZN]
      temp=ZN/2.0D0
      CALL PSI2(temp,PG2ZNmez)                   !PolyGamma[2, ZN/2]
      CALL PSI2(ZNp3,PG2ZNp3)                    !PolyGamma[2, 3 + ZN]
      temp=ZN/2.0D0+1.0D0
      CALL PSI3(temp,PG3ZNmezp1)                 !PolyGamma[3, 1 + ZN/2]
      CALL PSI3(ZNm1,PG3ZNm1)                    !PolyGamma[3, -1 + ZN]
      CALL PSI3(ZN,PG3ZN)                        !PolyGamma[3, ZN]
      temp=ZNp1/2.0D0
      CALL PSI3(temp,PG3ZNp1mez)                 !PolyGamma[3, ZNp1/2]
      temp=ZNm1/2.0D0
      CALL PSI0(temp,PG0ZNm1mez)                 !PolyGamma[0, ZNm1/2]
!--------------------------------------------------------------------
! NB:The Blumlein's convention to calculate the Mellin transform is z^N (not z^(N-1) )
!--------------------------------------------------------------------
      DDACG1ZNm1=DDACG1(ZNm2)!DDACG1[-1 + ZN]
      DACG2ZNm1=DACG2(ZNm2)!DACG2[-1 + ZN]
      DACG2ZN=DACG2(ZNm1)!DACG2[ZN]
      DACG2ZNp1=DACG2(ZN)!DACG2[1 + ZN]
      DACG2ZNp2=DACG2(ZNp1)!DACG2[2 + ZN]
      DACG2ZNp3=DACG2(ZNp2)!DACG2[3 + ZN]
      DACG4ZNm1=DACG4(ZNm2)!DACG4[-1 + ZN]
      DACG4ZNm1=DACG4(ZNm2)!DACG4[-1 + ZN]
      ACG1ZNp1=ACG1(ZN)!ACG1[1 + ZN]
      ACG1ZNm1=ACG1(ZNm2)!ACG1[-1 + ZN]
      ACG1ZNp2=ACG1(ZNp1)!ACG1[2 + ZN]
      ACG1ZNp3=ACG1(ZNp2)!ACG1[3 + ZN]
      ACG5ZNm1=ACG5(ZNm2)!ACG5[-1 + ZN]
      ACG6ZNm1=ACG6(ZNm2)!ACG6[-1 + ZN]
      ACG7ZNm1=ACG7(ZNm2)!ACG7[-1 + ZN]
      ACG9ZNm1=ACG9(ZNm2)!ACG9[-1 + ZN]
      ACG13ZNm1=ACG13(ZNm2)!ACG13[-1 + ZN]
      ACG13ZN=ACG13(ZNm1)!ACG13[ZN]
      ACG13ZNp1=ACG13(ZN)!ACG13[1 + ZN]
      ACG13ZNp2=ACG13(ZNp1)!ACG13[2 + ZN]
      ACG13ZNp3=ACG13(ZNp2)!ACG13[3 + ZN]
      ACG20ZNm1=ACG20(ZNm2)!ACG20[-1 + ZN]
      ACG22ZNm1=ACG22(ZNm2)!ACG22[-1 + ZN]


!      write(*,*) ZN,DDACG1ZNm1,DACG2ZNm1,DACG2ZN


      MellinH2gg=(9.0D0*CF**2)/4.0D0-(CA**2*(8.0D0+ ZN*(24.0D0 +
     # ZN*(31.0D0 + 2.0D0*ZN*
     # (9.0D0 + 2.0D0*ZN))))*DL**2)/
     # (ZN2*ZNp1**2*ZNp2**2) + CA**2*(-DDACG1ZNm1/2.0D0 + DACG2ZNm1 +
     # 2.0D0*DACG2ZN +3.0D0*DACG2ZNp1 +2.0D0*DACG2ZNp2 + DACG2ZNp3 +
     # DACG4ZNm1+
     # (4.0D0*ACG1ZNp1)/ZN2 + (2.0D0*ACG1ZNp2)/ZNp1**2 +
     # (2.0D0*ACG1ZNp3)/ZNp2**2 +
     # 2.0D0*ACG13ZNm1 + 4.0D0*ACG13ZN + 6.0D0*ACG13ZNp1 +
     # 4.0D0*ACG13ZNp2 +
     # 2.0D0*ACG13ZNp3 + 3.0D0*ACG20ZNm1 - 2.0D0*ACG22ZNm1 +
     # ACG5ZNm1 -
     # 2.0D0*ACG6ZNm1 - 3.0D0*ACG7ZNm1 + 2.0D0*ACG9ZNm1) +
     # (CA**2*(32.0D0*ZN*ZNb*(-54.0D0 +
     # ZN2*(162.0D0+ZN*(-47.0D0 + ZN*(-648.0D0 + ZN*(465.0D0 + ZN*
     # (-108.0D0+101.0D0*ZN*
     # (-3.0D0 + ZN2)))))))*PG0ZNm1 +
     # (186624.0D0+ZN*(312192.0D0 +ZN*(-463296.0D0 +ZN*(-948160.0D0 +
     # ZN*(707600.0D0 + 756.0D0*LT*ZNm1**4*ZNp1**4*ZNp2**4 +
     # ZN*(1232712.0D0 + ZN*(-1566444.0D0 + ZN*(-2736488.0D0 + ZN*
     # (204381.0D0 +
     # ZN*(2124608.0D0 + ZN*(390184.0D0 + ZN*(-1291248.0D0 + ZN*
     # (-807822.0D0 +
     # ZN*(86264.0D0 + ZN*(236748.0D0 + ZN*(83096.0D0 + 9561.0D0*
     # ZN))))))))))))))) + 432.0D0*ZNm1**4*ZN*ZNp1*ZNp2*
     # (16.0D0+ZN*(72.0D0+ZN*(132.0D0+ZN*(119.0D0 + ZN*(57.0D0 +
     # ZN*(15.0D0+2.0D0*ZN))))))*
     # PG0ZNmez)/ZNp2**4))/(864.0D0*ZN4*ZNb**4) +
     # (CA**2*Pi2*(-288.0D0- 816.0D0*ZN + 356.0D0*ZN2+3312.0D0*ZN3 +
     # 3986.0D0*ZN4 +
     # 1884.0D0*ZN5 +
     # 314.0D0*ZN6+52.0D0*ZN2*Pi2+156.0D0*ZN3*Pi2+169.0D0*ZN4*Pi2 +
     # 78.0D0*ZN5*Pi2+
     # 13.0D0*ZN6*Pi2 + 192.0D0*ZN*DL+720.0D0*ZN2*DL+936.0D0*ZN3*DL +
     # 504.0D0*ZN4*DL+96.0D0*ZN5*DL+24.0D0*ZN2*(2.0D0+3.0D0*ZN+
     # ZN2)**2*ACG1ZNm1-
     # 12.0D0*ZN*(8.0D0 + 22.0D0*ZN + 23.0D0*ZN2 + 11.0D0*ZN3 +
     # 2.0D0*ZN4)*PG0ZNmez +
     # 96.0D0*ZN*PG0ZNp1mez + 264.0D0*ZN2*PG0ZNp1mez +
     # 276.0D0*ZN3*PG0ZNp1mez + 132.0D0*ZN4*PG0ZNp1mez +
     # 24.0D0*ZN5*PG0ZNp1mez - 48.0D0*ZN2*PG1ZNm1 -
     # 144.0D0*ZN3*PG1ZNm1 - 156.0D0*ZN4*PG1ZNm1 -
     # 72.0D0*ZN5*PG1ZNm1 - 12.0D0*ZN6*PG1ZNm1))/
     # (144.0D0*ZN2*ZNp1**2*ZNp2**2) +
     # CA**2*(-1.0D0/(12.0D0*ZNp1**2)+(-2.0D0/ZN3 + ZNp1**(-3)-1.0D0/
     # (12.0D0*ZNp1) -
     # ZNp2**(-3))*PG0ZNp1 +
     # ((-2.0D0/ZN3 + ZNp1**(-3) - ZNp2**(-3))*PG0ZNp2mez)/2.0D0 -
     # (PG0ZNm1**2*PG1ZNm1)/2.0D0 + PG1ZNm1**2 -
     # ((8.0D0 + ZN*(24.0D0 + ZN*(23.0D0+2.0D0*ZN*(5.0D0 + ZN))))*
     # PG1ZNmez)/
     # (4.0D0*ZN2*ZNp1**2*ZNp2**2) + PG0ZN**2*PG1ZN) +
     # (CA**2*((-2.0D0/ZN2 + ZNp1**(-2) - ZNp2**(-2))*PG1ZNmezp1 +
     # 2.0D0*(-4.0D0*PG1ZN**2 + ((8.0D0+ZN*(24.0D0+ZN*(23.0D0 +
     # 2.0D0*ZN*(5.0D0 + ZN))))*
     # PG1ZNp1mez)/(ZN2*ZNp1**2*ZNp2**2) +
     # (2.0D0*(-8.0D0 + ZN2*(47.0D0 + ZN*(56.0D0+ZN*(30.0D0+7.0D0*
     # ZN))))*PG1ZNp1)/
     # (ZNm1*ZN2*ZNp1**2*ZNp2**2) - 3.0D0*PG0ZNp1**2*
     # PG1ZNp1 + 6.0D0*PG1ZNp1**2 +
     # 2.0D0*PG0ZNp2**2*PG1ZNp2 - 4.0D0*PG1ZNp2**2 -
     # PG0ZNp3**2*PG1ZNp3 + 2.0D0*PG1ZNp3**2 +
     # PG2ZNm1/ZNb)))/4.0D0 +
     # (CA**2*(-96.0D0*(4.0D0 + ZN2*(7.0D0 + 4.0D0*ZN))*PG2ZNm1 +
     # ZNm1*(12.0D0*(4.0D0 + ZN*(5.0D0 + 2.0D0*ZN))*PG2ZNmez -
     # 12.0D0*(4.0D0 + ZN*(5.0D0 + 2.0D0*ZN))*PG2ZNp1mez +
     # ZN*ZNp1*ZNp2*(192.0D0*PG0ZNp3*PG2ZNp3 +
     # PG3ZNmezp1 - 80.0D0*PG3ZNm1 +
     # 16.0D0*PG3ZN - PG3ZNp1mez))))/
     # (192.0D0*ZN*ZNp2*ZNb) +
     # (CA**2*(504.0D0+2.0D0*ZN*(299.0D0+ZN*(55.0D0 - ZN*(173.0D0 +
     # 55.0D0*ZN))) +
     # 9.0D0*ZNm1*ZN*ZNp1*ZNp2*(PG0ZNm1mez -
     # PG0ZNmez - 14.0D0*PG0ZN))*ZET3)/
     # (36.0D0*ZNm1*ZN*ZNp1*ZNp2) +
     # (CA*(-45.0D0 -
     # (2.0D0*(144.0D0 + ZN*(48.0D0 + ZN*(-1964.0D0+ZN*(-6382.0D0+
     # 4.0D0*GE*ZNm1*ZNp1**2*
     # ZNp2*(47.0D0 + 56.0D0*ZN) + ZN*(-8669.0D0 + ZN*(-3142.0D0 +
     # ZN*(3748.0D0+
     # ZN*(3740.0D0 + 861.0D0*ZN))))))))*nf)/(ZN3*ZNp1**3*
     # plinomZNm) - 448.0D0*nf*PG0ZN + (72.0D0*nf*PG0ZNp2)/
     # ZNp1-12.0D0*(CF*(435.0D0+99.0D0*LT+54.0D0*Pi2) + 2.0D0*nf*
     # (5.0D0*Pi2 +16.0D0
     # *ZET3))))/
     # 864.0D0 + (CA**2*GE*(-2592.0D0 + ZN*ZNp1*
     # (-1296.0D0 + 54.0D0*GE*ZNm1**3*ZNp2*
     # (8.0D0 + ZN*(24.0D0 + ZN*(23.0D0 + 2.0D0*ZN*(5.0D0 + ZN)))) +
     # ZN*(11016.0D0 + ZN*(-9640.0D0+ ZN*(-14936.0D0+ ZN*(10346.0D0+
     # ZN*(6269.0D0+
     # ZN*(-9818.0D0+ ZN*(-6304.0D0+ ZN*(4.0D0+ ZN)*(395.0D0+
     # 404.0D0*ZN))))))))) -
     # 54.0D0*ZN*ZNp2**3*ZNb*(2.0D0*(-1.0D0+ ZN2*(5.0D0 + 2.0D0*ZNm1*
     # ZN))*
     # PG0ZNm1 + ZNm1*ZN*ZNp1*
     # ((-2.0D0 + ZN*(2.0D0 - 4.0D0*ZN + GE*ZNb))*PG1ZNm1 +
     # ZN*ZNb*(2.0D0*PG0ZNp3*PG1ZNp3 -
     # 2.0D0*PG2ZNm1 + 7.0D0*ZET3)))))/
     # (108.0D0*ZN3*ZNp2**3*ZNb**3) +
     # (CF*(-2.0D0 +
     # (nf*(24.0D0+ ZN*ZNp1*(48.0D0+ ZN*(-18.0D0+ ZN*(-72.0D0+ ZN*
     # (-62.0D0+211.0D0*ZN +
     # 142.0D0*ZN2-152.0D0*ZN3-164.0D0*ZN4-41.0D0*ZN5+12.0D0*LT*
     # ZNp1**3*plinomZNm +
     # 24.0D0*ZNp1**3*plinomZNm*ZET3))))))/
     # (ZN4*ZNp1**4*plinomZNm)))/24.0D0
!      write (*,*) 'H2ggZN',MellinH2gg,ZN

      RETURN
      END

      complex*16 function MellinH2qq(ZN)
      implicit none
      complex*16 ZN
      real *8 CF

      CF=4d0/3
      MellinH2qq=-(CF**2*(-(1+ZN)**(-2)+(4d0*(-3+2*ZN))/(-1+ZN)**2
     &           - (4 + 8*ZN)/ZN**2))/4d0
!      write (*,*) 'H2qqZN',MellinH2qq,ZN
      return
      end




      COMPLEX*16 FUNCTION ACG1(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,T
      REAL*8 AK2,DL,ZERO,ONE
      INTEGER L,K
C
      DIMENSION AK2(10)
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/

!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACCON1/ DL
!      COMMON/VAL   / ZERO,ONE
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
C
      T=DCMPLX(ZERO,ZERO)
      DO 1 L=2,11
      K=L-1
      T=T+AK2(K)/(ZN+DBLE(K+1))
1     CONTINUE
C
      ACG1=(DL*DL- ZN*T)/2.0D0
C
      RETURN
      END
      COMPLEX*16   FUNCTION DACG1(ZN)

C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT none
      COMPLEX*16 ZN,T1,T2
      REAL*8 AK2,DL,ZERO
      INTEGER L,K

      DIMENSION AK2(10)
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/
      DL = DLOG(2.0D0)
      ZERO=0.0D0

      T1=DCMPLX(ZERO,ZERO)
      T2=DCMPLX(ZERO,ZERO)
      DO 1 L=2,11
      K=L-1
      T1=T1+AK2(K)/(ZN+DBLE(L))
      T2=T2+AK2(K)/(ZN+DBLE(L))/(ZN+DBLE(L))
1     CONTINUE

      DACG1=(ZN*T2 - T1)/(2.0D0)

      RETURN
      END





C --------------------------------------------------------------------------------------------------------------------------------

      COMPLEX*16   FUNCTION DDACG1(ZN)

C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT none
      COMPLEX*16 ZN,T1,T2,ZNPKP1
      REAL*8 AK2,DL,ZERO
      INTEGER L,K

      DIMENSION AK2(10)
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/
      DL = DLOG(2.0D0)
      ZERO=0.0D0

      T1=DCMPLX(ZERO,ZERO)
      T2=DCMPLX(ZERO,ZERO)
      DO 1 L=2,11
      K=L-1
      ZNPKP1=ZN+DBLE(L)
      T1=T1+AK2(K)/ZNPKP1/ZNPKP1
      T2=T2+AK2(K)/ZNPKP1/ZNPKP1/ZNPKP1
1     CONTINUE

      DDACG1=(T1-ZN*T2)
c      write(*,*) ZN,T1,T2,DACG1
      RETURN
      END




      COMPLEX*16 FUNCTION DACG2(ZN)
C     ----------------------------
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X)*LOG(1+X)**2/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,T1,T2
      REAL*8 AK3,DL,ZERO,ONE
      INTEGER K,L
C
      DIMENSION AK3(11)
      DATA AK3/9.99999989322696D-1,
     &        -1.49999722020708D+0,
     &         1.74988008499745D+0,
     &        -1.87296689068405D+0,
     &         1.91539974617231D+0,
     &        -1.85963744001295D+0,
     &         1.62987195424434D+0,
     &        -1.17982353224299D+0,
     &         6.28710122994999D-1,
     &        -2.11307487211713D-1,
     &         3.28953352932140D-2/
!      COMMON/ACLOG3/ AK3(11)
!      COMMON/ACCON1/ DL
!      COMMON/VAL   / ZERO,ONE
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0

C
      T1=DCMPLX(ZERO,ZERO)
      T2=DCMPLX(ZERO,ZERO)
      DO 1 L=3,13
      K=L-2
      T1=T1+AK3(K)/(ZN+DBLE(L))
      T2=T2+AK3(K)/(ZN+DBLE(L))/(ZN+DBLE(L))
1     CONTINUE
C
      DACG2=(ZN*T2 - T1)/3.0D0
C
      RETURN
      END


      COMPLEX*16 FUNCTION ACG4(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI2(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZL,ZNL,ZNL1,T,V1
      REAL*8 AK1,DL,ZERO,ONE,ZET2,PI
      INTEGER L


C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/VAL   / ZERO,ONE
      DIMENSION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/

      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0


C
      T=DCMPLX(-ZET2/2.0D0*DL,ZERO)
C
      DO 1 L=1,9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL =ZN+ZL
      ZNL1=ZNL+DCMPLX(ONE,ZERO)
      CALL BET(ZNL1,V1)
C
      T=T+AK1(L)*(ZN/ZNL*ZET2/2.0D0+ZL/ZNL**2*(DL-V1))
1     CONTINUE
C
      ACG4=T
C
      RETURN
      END

      COMPLEX*16 FUNCTION DACG4(ZN)
C     ----------------------------
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: Log(X)*LI2(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZL,ZNL,ZNL1,T,V1,V2
      REAL*8 AK1,DL,ZERO,ONE,ZET2,PI
      INTEGER L

C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/

      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0

      T=DCMPLX(ZERO,ZERO)
C
      DO 1 L=1,9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL =ZN+ZL
      ZNL1=ZNL+DCMPLX(ONE,ZERO)
      CALL BET(ZNL1,V1)
      CALL BET1(ZNL1,V2)
C
      T=T+AK1(L)*L*(ZET2/2.0D0-V2-2*(DL-V1)/ZNL)/ZNL/ZNL
1     CONTINUE
C
      DACG4=T
C
      RETURN
      END


      COMPLEX*16 FUNCTION ACG5(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X) LI2(X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZL,ZNL,ZL1,T,Z,PS,S1P,S1
      REAL*8 AK1,DL,ZERO,ONE,ZET2,PI,GE
      INTEGER L

C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/VAL   / ZERO,ONE
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON4/ GE
C
      DIMENSION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/

      GE   = 0.57721566490153D+0
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0

      T=DCMPLX(ZERO,ZERO)
      Z=ZN
      DO 1 L=1,9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL=Z+ZL
      ZL1=ZNL+ONE
      CALL PSI0(ZL1,PS)
      CALL PSI1(ZL1,S1P)
      S1=PS+GE
      T=T-AK1(L)*ZL/ZNL**2*(ZET2+S1P-2.0D0*S1/ZNL)
1     CONTINUE
C
      ACG5=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG6(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI3(X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZL,ZNL,ZNL1,T,V1,S1
      REAL*8 AK1,DL,ZERO,ONE,ZET3,GE,ZET2,PI
      INTEGER L

C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/ACCON4/ GE
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/

      ZET3 = 1.20205690315959428540D0
      GE   = 0.57721566490153D+0
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0


      T=DCMPLX(DL*ZET3,ZERO)
      DO 1 L=1,9
      ZL=DCMPLX(DBLE(L),ZERO)
      ZNL=ZN+ZL
      ZNL1=ZNL+ONE
      CALL PSI0(ZNL1,V1)
      S1=V1+GE
C
      T=T-AK1(L)*(ZN/ZNL*ZET3+ZL/ZNL**2*(ZET2-S1/ZNL))
1     CONTINUE
C
      ACG6=T
C
      RETURN
      END
      COMPLEX*16 FUNCTION ACG7(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LI3(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZL,ZNL,ZNL1,T,V1
      REAL*8 AK1,DL,ZERO,ONE,ZET2,ZET3,PI
      INTEGER L

C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION AK1(9)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/

      ZET3 = 1.20205690315959428540D0
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0


      T=DCMPLX(-3.0D0*ZET3/4.0D0*DL,ZERO)
C
      DO 1 L=1,9
      ZL =DCMPLX(DBLE(L),ZERO)
      ZNL =ZN+ZL
      ZNL1=ZNL+DCMPLX(ONE,ZERO)
      CALL BET(ZNL1,V1)
C
      T=T+AK1(L)*(ZN/ZNL*3.0D0*ZET3/4.0D0+ZL/ZNL**2/2.0D0*ZET2
     & -ZL/ZNL**3*(DL-V1))
1     CONTINUE
C
      ACG7=T
C
      RETURN
      END


      COMPLEX*16 FUNCTION ACG9(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: S12(-X)/(1+X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,T,T1,T2,ZNKL,ZK,ZNK
      REAL*8 AK1,AK2,AK3,DL,ZERO,ONE,ZET3,PI
      INTEGER K,L,L1
C
!      COMMON/ACLOG1/ AK1(9)
!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACLOG3/ AK3(11)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON3/ ZET3
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION AK1(9),AK2(10),AK3(11)
      DATA AK1/0.999999974532240D+0,
     &        -0.499995525890027D+0,
     &         0.333203435554182D+0,
     &        -0.248529457735332D+0,
     &         0.191451164493502D+0,
     &        -0.137466222203386D+0,
     &         0.792107405737825D-1,
     &        -0.301109652783781D-1,
     &         0.538406198111749D-2/
c
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/
C
      DATA AK3/9.99999989322696D-1,
     &        -1.49999722020708D+0,
     &         1.74988008499745D+0,
     &        -1.87296689068405D+0,
     &         1.91539974617231D+0,
     &        -1.85963744001295D+0,
     &         1.62987195424434D+0,
     &        -1.17982353224299D+0,
     &         6.28710122994999D-1,
     &        -2.11307487211713D-1,
     &         3.28953352932140D-2/

      ZET3 = 1.20205690315959428540D0
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0


      T=DCMPLX(ZET3*DL/8.0D0,ZERO)
      DO 1 K=1,9
      T1=DCMPLX(ZERO,ZERO)
      DO 2 L=2,11
      L1=L-1
      ZNKL=ZN+DCMPLX(DBLE(K+L),ZERO)
      T1=T1+AK2(L1)/ZNKL
2     CONTINUE
      ZK=DCMPLX(DBLE(K),ZERO)
      ZNK=ZN+ZK
      T=T-AK1(K)*ZN/ZNK*(ZET3/8.0D0-T1/2.0D0)
1     CONTINUE
      T2=DCMPLX(ZERO,ZERO)
      DO 3 K=3,13
      L=K-2
      ZK=DCMPLX(DBLE(K),ZERO)
      ZNK=ZN+ZK
      T2=T2+AK3(L)/ZNK
3     CONTINUE
      T=T-T2/2.0D0
C
      ACG9=T
C
      RETURN
      END


      COMPLEX*16 FUNCTION ACG13(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(1+X)/(1+X)*LI2(-X)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,T0,T1,T2,V1,ZNK1
      REAL*8 AK2,AK3,DL,ZERO,ONE,ZET2,PI
      INTEGER L,K

C
!      COMMON/ACLOG2/ AK2(10)
!      COMMON/ACLOG3/ AK3(11)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION AK2(10),AK3(11)
      DATA AK2/0.999999980543793D+0,
     &        -0.999995797779624D+0,
     &         0.916516447393493D+0,
     &        -0.831229921350708D+0,
     &         0.745873737923571D+0,
     &        -0.634523908078600D+0,
     &         0.467104011423750D+0,
     &        -0.261348046799178D+0,
     &         0.936814286867420D-1,
     &        -0.156249375012462D-1/
C
      DATA AK3/9.99999989322696D-1,
     &        -1.49999722020708D+0,
     &         1.74988008499745D+0,
     &        -1.87296689068405D+0,
     &         1.91539974617231D+0,
     &        -1.85963744001295D+0,
     &         1.62987195424434D+0,
     &        -1.17982353224299D+0,
     &         6.28710122994999D-1,
     &        -2.11307487211713D-1,
     &         3.28953352932140D-2/


      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0

      T0=DCMPLX(-1.0D0/4.0D0*ZET2*DL**2,ZERO)
C
      T1=DCMPLX(ZERO,ZERO)
      DO 1 L=3,13
      K=L-2
      T1=T1+AK3(K)/(ZN+DBLE(L))
1     CONTINUE
C
      T2=DCMPLX(ZERO,ZERO)
      DO 2 L=2,11
      K=L-1
      ZNK1=ZN+DBLE(L+1)
      CALL BET(ZNK1,V1)
      T2=T2+AK2(K)*ZN/(ZN+DBLE(L))*(ZET2/2.0D0-(DL-V1)/(ZN+DBLE(L)))
2     CONTINUE
C
      ACG13=T0+(T1+T2)/2.0D0
C
      RETURN
      END

      COMPLEX*16 FUNCTION ACG20(ZZ)
C     -----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: (LI3(X)-ZETA3)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,S1,S2,ZN1,ZZ,ZNK,ZNK1,R,T,T1,T2,PS,PS1
      REAL*8 DL,ZERO,ONE,ZETA2,ZET2,ZETA3,ZET3,GE,P24,P34,P22,CK2,CK4,PI
      INTEGER L,K

C
!      COMMON/ACLOG7/ CK2(13)
!      COMMON/ACLOG9/ CK4(13)
!      COMMON/POLY2 / P22(4)
!      COMMON/POLY4 / P24(5),P34(5)
!      COMMON/ACCON1/ DL
!      COMMON/ACCON2/ ZET2
!      COMMON/ACCON3/ ZET3
!      COMMON/ACCON4/ GE
!      COMMON/VAL   / ZERO,ONE
C
      DIMENSION CK2(13),CK4(13),P22(4),P24(5),P34(5)
      DL = DLOG(2.0D0)
      ZERO=0.0D0
      ONE =1.0D0
      ZET3 = 1.20205690315959428540D0
      GE   = 0.57721566490153D+0
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0
      ZETA2=ZET2
      ZETA3=ZET3

      DATA CK2/ 1.9891961855478257D-8,!0.480322239287449D+0,
     &         -6.050451121009104D-6,!-0.168480825099580D+1,
     &          2.1253046159349207,!0.209270571620726D+1,
     &         -1.0523034804772378D0,!-0.101728150275998D+1,
     &          0.160179976133047D+0,
     &         -0.351982236677917D+0,
     &          0.141033316846244D+1,
     &         -0.353343997843391D+1,
     &          0.593934696819832D+1,
     &         -0.660019784804042D+1,
     &          0.466330349413063D+1,
     &         -0.189825467489058D+1,
     &          0.339772909487512D+0/
C
      DATA CK4/-1.844613928270178D-8,!0.192962504274437D+0,
     &          2.2150086978693064D0,!0.000005641557253D+0,
     &         -0.9133677154535804D0,!-0.196891075399448D+1,
     &          3.4783104357500143D0,!0.392919138747074D+1,
     &         -2.823955592989266D0,!-0.290306105685546D+1,
     &          0.992890266001707D+0,
     &         -0.130026190226546D+1,
     &          0.341870577921103D+1,
     &         -0.576763902370864D+1,
     &          0.645554138192407D+1,
     &         -0.459405622046138D+1,
     &          0.188510809558304D+1,
     &         -0.340476080290674D+0/



C----------------------------------------------------------------------

!      P12(1)=ZETA3-11.0D0/6.0*ZETA2+4.0D0/3.0D0!=-0.48032221939548725
!      P12(2)=3.0D0*ZETA2-13.0D0/4.0D0!=1.684802200544679
!      P12(3)=-3.0D0/2.0D0*ZETA2+5.0D0/2.0D0!=0.0325988997276605
!      P12(4)=1.0D0/3.0D0*ZETA2-7.0D0/12.0D0!=-0.03502197771725779
C-------------------------------
      P22(1)=-1.0D0
      P22(2)=5.0D0/2.0D0
      P22(3)=-2.0D0
      P22(4)=1.0D0/2.0D0

C-------------------------------
!      P14(1)=257.D0/144.0D0-205.0D0/72.0D0*ZET2+ZET2**2
!      P14(2)=-167.0D0/36.0D0+25.0D0/6.0D0*ZET2
!      P14(3)=101.0D0/24.0D0-23.0D0/12.0D0*ZET2
!      P14(4)=-59.0D0/36.0D0+13.0D0/18.0D0*ZET2
!      P14(5)=41.0D0/144.0D0-ZET2/8.0D0
C-------------------------------
      P24(1)=-167.0D0/36.0D0+25.0D0/6.0D0*ZET2
      P24(2)=235.0D0/18.0D0-8.0D0*ZET2
      P24(3)=-40.0D0/3.0D0+6.0D0*ZET2
      P24(4)=109.0D0/18.0D0-8.0D0/3.0D0*ZET2
      P24(5)=-41.0D0/36.0D0+ZET2/2.0D0
C-------------------------------
      P34(1)=35.0D0/12.0D0
      P34(2)=-26.0D0/3.0D0
      P34(3)=19.0D0/2.0D0
      P34(4)=-14.0D0/3.0D0
      P34(5)=11.0D0/12.0D0
C----------------------------------------------------------------------
C
C >>> ACCOUNT FOR POLYNOM PARTS
C
!      CK2(1)=CK2(1)+P12(1)
!      CK2(2)=CK2(2)+P12(2)
!      CK2(3)=CK2(3)+P12(3)
!      CK2(4)=CK2(4)+P12(4)
C
!      write(*,*) CK2(1),CK2(2),CK2(3),CK2(4),CK2(5)
!      CK4(1)=CK4(1)+P14(1)
!      CK4(2)=CK4(2)+P14(2)
!      CK4(3)=CK4(3)+P14(3)
!      CK4(4)=CK4(4)+P14(4)
!      CK4(5)=CK4(5)+P14(5)
!      write(*,*) CK4(1),CK4(2),CK4(3),CK4(4),CK4(5),CK4(6)
C
      T=DCMPLX(ZET2**2/2.0D0,ZERO)
C
      ZN = ZZ
      ZN1=ZN+ONE
      CALL PSI0(ZN1,PS)
      S1=PS+GE
      T=T-ZET3*S1
C
      DO 1 K=1,13
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
      R=ZN/ZNK
C
      T=T+CK2(K)*R*S1
1     CONTINUE
C
      DO 3 K=1,4
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      S1=PS+GE
      S2=-PS1+ZET2
      R=ZN/ZNK
      T1=S1**2+S2
C
      T=T-P22(K)*R*T1
3     CONTINUE
C
      DO 4 K=1,13
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      S1=PS+GE
C
      T=T-CK4(K)/ZNK*ZN/2.0D0
4     CONTINUE
C
      DO 5 K=1,5
      L=K-1
      ZNK=ZN+DBLE(L)
      ZNK1=ZNK+ONE
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      S1=PS+GE
      S2=-PS1+ZET2
      T1=-S1/ZNK
      T2=(S1**2+S2)/ZNK
C
      T=T-(P24(K)*T1+P34(K)*T2)*ZN/2.0D0
5     CONTINUE
C
C
      ACG20=T

C
      RETURN
      END


      COMPLEX*16 FUNCTION ACG22(ZN)
C     ----------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  MELLIN TRANSFORM FOR: LOG(X) LI2(X)/(X-1)
C---  FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT COMPLEX*16 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZN,ZNK1,PS,PS1,PS2,T,TA
      REAL*8 CK1,P21,GE,ZERO,ONE,ZET2,PI
      INTEGER L
C
!      COMMON/ACLOG6/ CK1(12)
!      COMMON/POLY1 / P21(4)
!      COMMON/ACCON4/ GE
!      COMMON/VAL   / ZERO,ONE
C

      DIMENSION CK1(12),P21(4)
      PI = 3.141592653589793238462643D0
      ZET2 = PI**2/6.0D0
      DATA CK1/ 2.201218318731435D-8,!-0.283822933724932D+0,
     &          2.833327652357064D0,!0.999994319023731D+0,
     &         -1.8330909624101532D0,!-0.124975762907682D+1,
     &          0.7181879191200942,!0.607076808008983D+0,
     &         -0.280403220046588D-1,
     &         -0.181869786537805D+0,
     &          0.532318519269331D+0,
     &         -0.107281686995035D+1,
     &          0.138194913357518D+1,
     &         -0.111100841298484D+1,
     &          0.506649587198046D+0,
     &         -0.100672390783659D+0/
C----------------------------------------------------------------------
!      P11(1)=-49.0D0/36.0+ZET2
!      P11(2)=11.0D0/6.0D0
!      P11(3)=-7.0D0/12.0D0
!      P11(4)=1.0D0/9.0D0
C---------------------------------------------------------
      P21(1)=11.0D0/6.0
      P21(2)=-3.0D0
      P21(3)=3.0D0/2.0D0
      P21(4)=-1.0D0/3.0D0
C
!      CK1(1)=CK1(1)+P11(1)
!      CK1(2)=CK1(2)+P11(2)
!      CK1(3)=CK1(3)+P11(3)
!      CK1(4)=CK1(4)+P11(4)
!      write(*,*) CK1(1),CK1(2),CK1(3),CK1(4),CK1(5)

      ZERO=0.0D0
      ONE =1.0D0
      GE   = 0.57721566490153D+0

      T=DCMPLX(ZERO,ZERO)
C
      DO 1 L=1,12
      ZNK1=ZN+DBLE(L)
      CALL PSI1(ZNK1,PS1)
!      write(*,*)T,CK1(L)
      T=T+CK1(L)*PS1
1     CONTINUE
C
      DO 2 L=1,4
      ZNK1=ZN+DBLE(L)
      CALL PSI0(ZNK1,PS)
      CALL PSI1(ZNK1,PS1)
      CALL PSI2(ZNK1,PS2)
      TA=(PS+GE)*PS1-PS2/2.0D0
      T=T-P21(L)*TA
2     CONTINUE
C
      ACG22=T
C
      RETURN
      END
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------

      SUBROUTINE GAMMAL(ZZ,RES)
C     -------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LOG(GAMMA(Z)) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
!      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT none
      COMPLEX*16 ZZ,Z,RES,T,T1,T2,X,X2,ONE
      REAL*8 PI,R
C
      Z=ZZ
      PI = 3.141592653589793238462643D0
C
      ONE=DCMPLX(1.0D0,0.0D0)
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T-LOG(Z)
      Z=Z+ONE
      GOTO 2
1     CONTINUE
C
      T1=Z*(LOG(Z)-1.0D0)+LOG(2.0D0*PI/Z)/2.0D0
C
      X=ONE/Z
      X2=X*X
      T2 = (1.D0/12.D0+(-1.D0/360.D0+(1.D0/1260.D0+(-1.D0/1680.D0+(1.D0/
     #1188.D0+(-691.D0/360360.D0+(1.D0/156.D0-3617.D0/122400.D0*X2)*X2
     #)*X2)*X2)*X2)*X2)*X2)*X
C
      RES=T1+T2+T
C
      RETURN
      END
      SUBROUTINE GAMMA(ZZ,RES)
C     ------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  GAMMA(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT none
      COMPLEX*16 ZZ,T,RES
C
      CALL GAMMAL(ZZ,T)
C
      RES=EXP(T)
C
      RETURN
      END
      SUBROUTINE BETA(AA,BB,RES)
C     --------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C---  BETA(A,B) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT none
      COMPLEX*16 T1,T2,T3,T,RES,AA,BB
C
      CALL GAMMAL(AA,T1)
      CALL GAMMAL(BB,T2)
      CALL GAMMAL(AA+BB,T3)
      T=T1+T2-T3
C
      RES=EXP(T)
C
      RETURN
      END
      SUBROUTINE PSI0(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C
      IMPLICIT none
!      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ,Y2   ! Dava Y2 come reale, ma e'complessa? Il risultato cambia, ma molto poco
      REAL*8 R
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T-ONE/Z
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
!      write(*,*) Y2
      T0 = (-1.D0/2.D0+(-1.D0/12.D0+(1.D0/120.D0+(-1.D0/252.D0+(1.D0/240
     #.D0+(-1.D0/132.D0+(691.D0/32760.D0+(-ONE/12.0D0+ONE*3617.0D0
     #             /8160.D0*Y2
     #  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y)*Y-LOG(Y) !formula(142)
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE PSI1(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI'(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T+ONE/Z**2
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
      T0 = (1.D0+(1.D0/2.D0+(1.D0/6.D0+(-1.D0/30.D0+
     &(1.D0/42.D0+(-1.D0/30.D0+(5.D0/66.D0-691.D0/2730.D0*Y2)
     &*Y2)*Y2)*Y2)*Y2)*Y)*Y)*Y
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE PSI2(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI''(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
      TWO=ONE*2.0D0
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T-TWO/Z**3
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
      T0 =(-1.D0+(-1.D0+(-1.D0/2.D0+(1.D0/6.D0+(-1.D0/6.D0+(3.D0/
     &10.D0+(-5.D0/6.D0+691.D0/210.D0*Y2)*Y2)*Y2)*Y2)*Y2)*Y)*Y)*Y2
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE PSI3(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  PSI'''(Z) FOR COMPLEX ARGUMENT
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,T,ONE,Y,T0,RES,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
      SIX=ONE*6.0D0
C
      Z=ZZ
      T=DCMPLX(0.0D0,0.0D0)
2     R=SQRT(DREAL(Z)**2+DIMAG(Z)**2)
      IF(R.GT.10.0D0) GOTO 1
      T=T+SIX/Z**4
      Z=Z+ONE
      GOTO 2
1     Y=ONE/Z
      Y2=Y*Y
C
      T0 = (2.D0+(3.D0+(2.D0+(-1.D0+(4.D0/3.D0+(-3.D0+(10.D0+(-691.D0/15
     #.D0+(280.D0-10851.D0/5.D0*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2  )*Y2
     #)*Y)*Y)*Y2*Y
C
      RES=T+T0
C
      RETURN
      END
      SUBROUTINE BET(ZZ,RES)
C     ----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  BET(Z) FOR COMPLEX ARGUMENT  $\beta(z)$
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,ONE,RES,V1,V2,Z1,Z2,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      Z1=(Z+ONE)/2.0D0
      Z2=Z/2.0D0
      CALL PSI0(Z1,V1)
      CALL PSI0(Z2,V2)
C
      RES=(V1-V2)/2.0D0
C
      RETURN
      END
      SUBROUTINE BET1(ZZ,RES)
C     -----------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  BET'(Z) FOR COMPLEX ARGUMENT  $\beta'(z)$
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8 (A-H,O-Z)
      COMPLEX*16 Z,ONE,RES,V1,V2,Z1,Z2,ZZ
C
      ONE=DCMPLX(1.0D0,0.0D0)
C
      Z=ZZ
      Z1=(Z+ONE)/2.0D0
      Z2=Z/2.0D0
      CALL PSI1(Z1,V1)
      CALL PSI1(Z2,V2)
C
      RES=(V1-V2)/4.0D0
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION FLI4(X)
C     -----------------------------------
C************************************************************************
C
C     J. BLUEMLEIN:  01.10.1999 (1.00)
C
C************************************************************************
C
C---  LI4(X) FOR -1. LE . X . LE .+1
C
C************************************************************************
C      IMPLICIT none
      IMPLICIT REAL*8(A-H,O-Z)
C
      A=1D0
      F=0D0
      AN=0D0
      TCH=1D-16
1     AN=AN+1D0
      A=A*X
      B=A/AN**4
      F=F+B
      IF(ABS(B)-TCH)2,2,1
2     FLI4 = F
      END


