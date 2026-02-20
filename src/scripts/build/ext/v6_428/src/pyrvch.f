 
C*********************************************************************
 
C...PYRVCH
C...Calculates R-violating chargino decay widths.
C...P. Z. Skands
 
      SUBROUTINE PYRVCH(KFIN,XLAM,IDLAM,LKNT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
      COMMON/PYMSRV/RVLAM(3,3,3), RVLAMP(3,3,3), RVLAMB(3,3,3)
C...Local variables.
      DOUBLE PRECISION XLAM(0:400)
      INTEGER IDLAM(400,3), PYCOMP
C...Information from main routine to PYRVGW
      COMMON/PYRVNV/AB(2,16,2),RMS(0:3),RES(6,2),INTRES(6,3),IDR,IDR2
     &     ,DCMASS,KFR(3)
C...Auxiliary variables needed for BV (RV Gauge STOre)
      COMMON/RVGSTO/XRESI,XRESJ,XRESK,XRESIJ,XRESIK,XRESJK,RVLIJK,RVLKIJ
     &     ,RVLJKI,RVLJIK
C...Running quark masses
      DOUBLE PRECISION RMQ(6)
C...Decay product masses on/off
      LOGICAL DCMASS
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/,/PYSSMT/,/PYMSRV/,/PYRVNV/,
     &     /RVGSTO/
 
 
C...IF R-VIOLATION ON.
      IF ((IMSS(51).GE.1).OR.(IMSS(52).GE.1).OR.(IMSS(53).GE.1)) THEN
        KFSM=KFIN-KSUSY1
        IF(KFSM.EQ.24.OR.KFSM.EQ.37) THEN
C...WHICH CHARGINO ?
          NCHI = 1
          IF (KFSM.EQ.37) NCHI = 2
 
C...Useful parameters for calculating the A and B constants.
C...SIGN OF MASS (Opposite convention as HERWIG)
          ISM  = 1
          IF (SMW(NCHI).LT.0D0) ISM = -1
          WMASS   = PMAS(PYCOMP(24),1)
          COSB    = 1/(SQRT(1+RMSS(5)**2))
          SINB    = RMSS(5)/SQRT(1+RMSS(5)**2)
          GW2     = 4*PARU(103)*PARU(1)/PARU(102)
          C1U     = UMIX(NCHI,2)/(SQRT(2D0)*COSB*WMASS)
          C1V     = VMIX(NCHI,2)/(SQRT(2D0)*SINB*WMASS)
          C2      = UMIX(NCHI,1)
          C3      = VMIX(NCHI,1)
C...Running masses at Q^2=MCHI^2.
          SQMCHI  = PMAS(PYCOMP(KFSM),1)**2
          DO 100 I=1,6
            RMQ(I)=PYMRUN(I,SQMCHI)
  100     CONTINUE
 
C... AB(x,y,z) coefficients:
C       x=1-2  : A or B coefficient  (1:A ; 2:B)
C       y=1-16 : Sparticle's SM code (1-6:d,u,s,c,b,t ;
C                                    11-16:e,nu_e,mu,...)
C       z=1-2  : Mass eigenstate number
          DO 110 I = 11,15,2
C...Intermediate sleptons
            AB(1,I,1)   = 0D0
            AB(1,I,2)   = 0D0
            AB(2,I,1)   = -PMAS(PYCOMP(I),1)*C1U*SFMIX(I,2) +
     &           SFMIX(I,1)*C2
            AB(2,I,2)   = -PMAS(PYCOMP(I),1)*C1U*SFMIX(I,4) +
     &           SFMIX(I,3)*C2
C...Intermediate sneutrinos
            AB(1,I+1,1) = -PMAS(PYCOMP(I),1)*C1U
            AB(1,I+1,2) = 0D0
            AB(2,I+1,1) = ISM*C3
            AB(2,I+1,2) = 0D0
C...Intermediate sdown
            J=I-10
            AB(1,J,1)   = -RMQ(J+1)*C1V*SFMIX(J,1)
            AB(1,J,2)   = -RMQ(J+1)*C1V*SFMIX(J,3)
            AB(2,J,1)   = -ISM*(RMQ(J)*C1U*SFMIX(J,2) - SFMIX(J,1)*C2)
            AB(2,J,2)   = -ISM*(RMQ(J)*C1U*SFMIX(J,4) - SFMIX(J,3)*C2)
C...Intermediate sup
            J=J+1
            AB(1,J,1)   = -RMQ(J-1)*C1U*SFMIX(J,1)
            AB(1,J,2)   = -RMQ(J-1)*C1U*SFMIX(J,3)
            AB(2,J,1)   = -ISM*(RMQ(J)*C1V*SFMIX(J,2) - SFMIX(J,1)*C3)
            AB(2,J,2)   = -ISM*(RMQ(J)*C1V*SFMIX(J,4) - SFMIX(J,3)*C3)
  110     CONTINUE
 
C...LLE TYPE R-VIOLATION
          IF (IMSS(51).GE.1) THEN
C...LOOP OVER DECAY MODES
            DO 140 ISC=0,26
 
C...CHI+ -> NUBAR_I + LEPTON+_J + NU_K.
              IF(MOD(ISC/9,3).NE.MOD(ISC/3,3)) THEN
                LKNT = LKNT+1
                IDLAM(LKNT,1) = -12 -2*MOD(ISC/9,3)
                IDLAM(LKNT,2) = -11 -2*MOD(ISC/3,3)
                IDLAM(LKNT,3) =  12 +2*MOD(ISC,3)
                XLAM(LKNT)    =  0D0
C...Set coupling, and decay product masses on/off
                RVLAMC        = GW2 * 5D-1 *
     &               RVLAM(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)+1)
     &               **2
                DCMASS=.FALSE.
                IF (IDLAM(LKNT,2).EQ.-15) DCMASS = .TRUE.
C...Resonance KF codes (1=I,2=J,3=K).
                KFR(1) = 0
                KFR(2) = 0
                KFR(3) = -IDLAM(LKNT,3)+1
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XLAM(LKNT))
                XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...KINEMATICS CHECK
                IF (XLAM(LKNT).EQ.0D0) THEN
                  LKNT=LKNT-1
                ENDIF
 
C * CHI+ -> NU_I + NU_J + LEPTON+_K. (NOTE: SYMM. IN I AND J)
  120           IF (MOD(ISC/9,3).LT.MOD(ISC/3,3)) THEN
                  LKNT = LKNT+1
                  IDLAM(LKNT,1) = 12 +2*MOD(ISC/9,3)
                  IDLAM(LKNT,2) = 12 +2*MOD(ISC/3,3)
                  IDLAM(LKNT,3) =-11 -2*MOD(ISC,3)
                  XLAM(LKNT)    = 0D0
C...Set coupling, and decay product masses on/off
                  RVLAMC = GW2 * 5D-1 *
     &              RVLAM(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)+1)**2
C...I,J SYMMETRY => FACTOR 2
                  RVLAMC=2*RVLAMC
                  DCMASS=.FALSE.
                  IF (IDLAM(LKNT,3).EQ.-15) DCMASS = .TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
                  KFR(1)=IDLAM(LKNT,1)-1
                  KFR(2)=IDLAM(LKNT,2)-1
                  KFR(3)=0
C...Calculate width.
                  CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &                 IDLAM(LKNT,3),XLAM(LKNT))
                 XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...KINEMATICS CHECK
                  IF (XLAM(LKNT).EQ.0D0) THEN
                    LKNT=LKNT-1
                  ENDIF

C * CHI+ -> LEPTON+_I + LEPTON+_J + LEPTON-_K (NOTE: SYMM. IN I AND J)
C * 19/04 2010: Bug corrected. Moved channel inside the I < J IF statement 
C *             from above, thanks to N.-E. Bomark.
                  LKNT = LKNT+1
                  IDLAM(LKNT,1) =-11 -2*MOD(ISC/9,3)
                  IDLAM(LKNT,2) =-11 -2*MOD(ISC/3,3)
                  IDLAM(LKNT,3) = 11 +2*MOD(ISC,3)
                  XLAM(LKNT)    = 0D0
C...Set coupling, and decay product masses on/off
                  RVLAMC = GW2 * 5D-1 *
     &              RVLAM(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)+1)**2
C...I,J SYMMETRY => FACTOR 2
                  RVLAMC=2*RVLAMC
                  DCMASS=.FALSE.
                  IF (IDLAM(LKNT,1).EQ.-15.OR.IDLAM(LKNT,2).EQ.-15
     &                 .OR.IDLAM(LKNT,3).EQ.15) DCMASS = .TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
                  KFR(1) =-IDLAM(LKNT,1)+1
                  KFR(2) =-IDLAM(LKNT,2)+1
                  KFR(3) = 0
C...Calculate width.
                  CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &                 IDLAM(LKNT,3),XLAM(LKNT))
                  XLAM(LKNT)=XLAM(LKNT)*RVLAMC
     &                 /((2*PARU(1)*RMS(0))**3*32)
C...KINEMATICS CHECK
                  IF (XLAM(LKNT).EQ.0D0) THEN
                    LKNT=LKNT-1
                  ENDIF
                ENDIF
              ENDIF
 140        CONTINUE
          ENDIF
 
C...LQD TYPE R-VIOLATION
          IF (IMSS(52).GE.1) THEN
C...LOOP OVER DECAY MODES
            DO 180 ISC=0,26
 
C...CHI+ -> NUBAR_I + DBAR_J + U_K
              LKNT = LKNT+1
              IDLAM(LKNT,1) =-12 -2*MOD(ISC/9,3)
              IDLAM(LKNT,2) = -1 -2*MOD(ISC/3,3)
              IDLAM(LKNT,3) =  2 +2*MOD(ISC,3)
              XLAM(LKNT)    =  0D0
C...Set coupling, and decay product masses on/off
              RVLAMC = 3. * GW2 * 5D-1 *
     &           RVLAMP(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)+1)**2
              DCMASS=.FALSE.
              IF (IDLAM(LKNT,2).EQ.-5.OR.IDLAM(LKNT,3).EQ.6)
     &             DCMASS = .TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)=0
              KFR(2)=0
              KFR(3)=-IDLAM(LKNT,3)+1
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XLAM(LKNT))
              XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...KINEMATICS CHECK
              IF (XLAM(LKNT).EQ.0D0) THEN
                LKNT=LKNT-1
              ENDIF
 
C * CHI+ -> LEPTON+_I + UBAR_J + U_K.
  150         LKNT = LKNT+1
              IDLAM(LKNT,1) =-11 -2*MOD(ISC/9,3)
              IDLAM(LKNT,2) = -2 -2*MOD(ISC/3,3)
              IDLAM(LKNT,3) =  2 +2*MOD(ISC,3)
              XLAM(LKNT)    =  0D0
C...Set coupling, and decay product masses on/off
              RVLAMC = 3. * GW2 * 5D-1 *
     &             RVLAMP(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)+1)**2
              DCMASS=.FALSE.
              IF (IDLAM(LKNT,1).EQ.-11.OR.IDLAM(LKNT,2).EQ.-6
     &             .OR.IDLAM(LKNT,3).EQ.6) DCMASS = .TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)=0
              KFR(2)=0
              KFR(3)=-IDLAM(LKNT,3)+1
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XLAM(LKNT))
              XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...KINEMATICS CHECK
              IF (XLAM(LKNT).EQ.0D0) THEN
                LKNT=LKNT-1
              ENDIF
 
C * CHI+ -> LEPTON+_I + DBAR_J + D_K.
  160         LKNT = LKNT+1
              IDLAM(LKNT,1) =-11 -2*MOD(ISC/9,3)
              IDLAM(LKNT,2) = -1 -2*MOD(ISC/3,3)
              IDLAM(LKNT,3) =  1 +2*MOD(ISC,3)
              XLAM(LKNT)    =  0D0
C...Set coupling, and decay product masses on/off
              RVLAMC = 3. * GW2 * 5D-1 *
     &             RVLAMP(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)+1)**2
              DCMASS = .FALSE.
              IF (IDLAM(LKNT,1).EQ.-15.OR.IDLAM(LKNT,2).EQ.-5
     &             .OR.IDLAM(LKNT,3).EQ.5) DCMASS = .TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)=-IDLAM(LKNT,1)+1
              KFR(2)=-IDLAM(LKNT,2)+1
              KFR(3)=0
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XLAM(LKNT))
              XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...KINEMATICS CHECK
              IF (XLAM(LKNT).EQ.0D0) THEN
                LKNT=LKNT-1
              ENDIF
 
C * CHI+ -> NU_I + U_J + DBAR_K.
  170         LKNT = LKNT+1
              IDLAM(LKNT,1) = 12 +2*MOD(ISC/9,3)
              IDLAM(LKNT,2) =  2 +2*MOD(ISC/3,3)
              IDLAM(LKNT,3) = -1 -2*MOD(ISC,3)
              XLAM(LKNT)    =  0D0
C...Set coupling, and decay product masses on/off
              DCMASS = .FALSE.
              RVLAMC = 3. * GW2 * 5D-1 *
     &             RVLAMP(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)+1)**2
              IF (IDLAM(LKNT,2).EQ.6.OR.IDLAM(LKNT,3).EQ.-5)
     &             DCMASS = .TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)=IDLAM(LKNT,1)-1
              KFR(2)=IDLAM(LKNT,2)-1
              KFR(3)=0
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XLAM(LKNT))
              XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...KINEMATICS CHECK
              IF (XLAM(LKNT).EQ.0D0) THEN
                LKNT=LKNT-1
              ENDIF
 
  180       CONTINUE
          ENDIF
 
C...UDD TYPE R-VIOLATION
C...These decays need special treatment since more than one BV coupling
C...contributes (with interference). Consider e.g. (symbolically)
C      |M|^2 = |l''_{ijk}|^2*(PYRVI1(RES_I) + PYRVI2(RES_I))
C             +|l''_{jik}|^2*(PYRVI1(RES_J) + PYRVI2(RES_J))
C             +l''_{ijk}*l''_{jik}*PYRVI3(PYRVI4(RES_I,RES_J))
C...The problem is that a single call to PYRVGW would evaluate all
C...these terms and sum them, but without the different couplings. The
C...way out is to call PYRVGW three times, once for the first line, once
C...for the second line, and then once for all the lines (it is
C...impossible to get just the last line out) without multiplying by
C...couplings. The last line is then obtained as the result of the third
C...call minus the results of the two first calls. Each term is then
C...multiplied by its respective coupling before the whole thing is
C...summed up in XLAM.
C...Note that with three interfering resonances, this procedure becomes
C...more complicated, as can be seen in the CHI+ -> 3*DBAR mode.
 
          IF (IMSS(53).GE.1) THEN
C...LOOP OVER DECAY MODES
            DO 190 ISC=1,25
 
C...CHI+ -> U_I + U_J + D_K
C...Decay mode I<->J symmetric.
              IF (MOD(ISC/9,3).LE.MOD(ISC/3,3).AND.ISC.NE.13) THEN
                LKNT = LKNT+1
                IDLAM(LKNT,1) =  2 +2*MOD(ISC/9,3)
                IDLAM(LKNT,2) =  2 +2*MOD(ISC/3,3)
                IDLAM(LKNT,3) =  1 +2*MOD(ISC,3)
                XLAM(LKNT)    =  0D0
C...Set coupling, and decay product masses on/off
                RVLAMC= 6. * GW2 * 5D-1
                RVLJIK= RVLAMB(MOD(ISC/3,3)+1,MOD(ISC/9,3)+1,MOD(ISC,3)
     &               +1)
                RVLIJK= RVLAMB(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)
     &               +1)
                IF (MOD(ISC/9,3).EQ.MOD(ISC/3,3)) RVLAMC = 5D-1
     &               * RVLAMC
                DCMASS=.FALSE.
                IF (IDLAM(LKNT,1).EQ.6.OR.IDLAM(LKNT,2).EQ.6
     &               .OR.IDLAM(LKNT,3).EQ.5) DCMASS =.TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1) = -IDLAM(LKNT,1)+1
                KFR(2) = 0
                KFR(3) = 0
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XRESI)
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1) = 0
                KFR(2) = -IDLAM(LKNT,2)+1
                KFR(3) = 0
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XRESJ)
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1) = -IDLAM(LKNT,1)+1
                KFR(2) = -IDLAM(LKNT,2)+1
                KFR(3) = 0
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XRESIJ)
                IF (ABS(XRESI+XRESJ-XRESIJ).GT.1D-4*XRESIJ) THEN
                  XRESIJ = XRESIJ-XRESI-XRESJ
                ELSE
                  XRESIJ = 0D0
                ENDIF
C...CALCULATE TOTAL WIDTH
                XLAM(LKNT) = RVLJIK**2 * XRESI + RVLIJK**2 * XRESJ
     &               + RVLJIK*RVLIJK * XRESIJ
                XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...KINEMATICS CHECK
                IF (XLAM(LKNT).EQ.0D0) THEN
                  LKNT=LKNT-1
                ENDIF
              ENDIF
C...CHI+ -> DBAR_I + DBAR_J + DBAR_K
C...Symmetry I<->J<->K.
              IF ((MOD(ISC/9,3).LE.MOD(ISC/3,3)).AND.(MOD(ISC/3,3).LE
     &             .MOD(ISC,3)).AND.ISC.NE.13) THEN
                LKNT = LKNT+1
                IDLAM(LKNT,1) = -1 -2*MOD(ISC/9,3)
                IDLAM(LKNT,2) = -1 -2*MOD(ISC/3,3)
                IDLAM(LKNT,3) = -1 -2*MOD(ISC,3)
                XLAM(LKNT)    =  0D0
C...Set coupling, and decay product masses on/off
                RVLAMC = 6. * GW2 * 5D-1
                RVLIJK = RVLAMB(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)
     &               +1)
                RVLKIJ = RVLAMB(MOD(ISC,3)+1,MOD(ISC/9,3)+1,MOD(ISC/3,3)
     &               +1)
                RVLJKI = RVLAMB(MOD(ISC/3,3)+1,MOD(ISC,3)+1,MOD(ISC/9,3)
     &               +1)
                DCMASS = .FALSE.
                IF (IDLAM(LKNT,1).EQ.-5.OR.IDLAM(LKNT,2).EQ.-5
     &               .OR.IDLAM(LKNT,3).EQ.-5) DCMASS = .TRUE.
C...Collect symmetry factors
                IF (MOD(ISC/9,3).EQ.MOD(ISC/3,3).OR.MOD(ISC/3,3).EQ
     &               .MOD(ISC,3).OR.MOD(ISC/9,3).EQ.MOD(ISC,3))
     &               RVLAMC = 5D-1 * RVLAMC
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1) = IDLAM(LKNT,1)-1
                KFR(2) = 0
                KFR(3) = 0
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XRESI)
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1) = 0
                KFR(2) = IDLAM(LKNT,2)-1
                KFR(3) = 0
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XRESJ)
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1) = 0
                KFR(2) = 0
                KFR(3) = IDLAM(LKNT,3)-1
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XRESK)
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1) = IDLAM(LKNT,1)-1
                KFR(2) = IDLAM(LKNT,2)-1
                KFR(3) = 0
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XRESIJ)
                IF (ABS(XRESI+XRESJ-XRESIJ).GT.1D-4*(XRESI+XRESJ)) THEN
                  XRESIJ = XRESI+XRESJ-XRESIJ
                ELSE
                  XRESIJ = 0D0
                ENDIF
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1) = 0
                KFR(2) = IDLAM(LKNT,2)-1
                KFR(3) = IDLAM(LKNT,3)-1
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XRESJK)
                IF (ABS(XRESJ+XRESK-XRESJK).GT.1D-4*(XRESJ+XRESK)) THEN
                  XRESJK = XRESJ+XRESK-XRESJK
                ELSE
                  XRESJK = 0D0
                ENDIF
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1) = IDLAM(LKNT,1)-1
                KFR(2) = 0
                KFR(3) = IDLAM(LKNT,3)-1
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XRESIK)
                IF (ABS(XRESI+XRESK-XRESIK).GT.1D-4*(XRESI+XRESK)) THEN
                  XRESIK = XRESI+XRESK-XRESIK
                ELSE
                  XRESIK = 0D0
                ENDIF
C...CALCULATE TOTAL WIDTH
                XLAM(LKNT) =
     &                 RVLIJK**2 * XRESI
     &               + RVLJKI**2 * XRESJ
     &               + RVLKIJ**2 * XRESK
     &               + RVLIJK*RVLJKI * XRESIJ
     &               + RVLIJK*RVLKIJ * XRESIK
     &               + RVLJKI*RVLKIJ * XRESJK
                XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2.*PARU(1)*RMS(0))**3*32)
C...KINEMATICS CHECK
                IF (XLAM(LKNT).EQ.0D0) THEN
                  LKNT=LKNT-1
                ENDIF
              ENDIF
  190       CONTINUE
          ENDIF
        ENDIF
      ENDIF
 
      RETURN
      END
