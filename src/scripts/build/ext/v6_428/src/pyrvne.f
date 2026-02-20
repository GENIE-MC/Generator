 
C*********************************************************************
 
C...PYRVNE
C...Calculates R-violating neutralino decay widths (pure 1->3 parts).
C...P. Z. Skands
 
      SUBROUTINE PYRVNE(KFIN,XLAM,IDLAM,LKNT)
 
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
      COMMON/PYRVNV/AB(2,16,2),RMS(0:3),RES(6,2),INTRES(6,3),IDR,IDR2
     &     ,DCMASS,KFR(3)
      DOUBLE PRECISION XLAM(0:400)
      DOUBLE PRECISION ZPMIX(4,4), NMIX(4,4), RMQ(6)
      INTEGER IDLAM(400,3), PYCOMP
      LOGICAL DCMASS
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/,/PYSSMT/,/PYMSRV/,/PYRVNV/
 
C...R-VIOLATING DECAYS
      IF ((IMSS(51).GE.1).OR.(IMSS(52).GE.1).OR.(IMSS(53).GE.1)) THEN
        KFSM=KFIN-KSUSY1
        IF(KFSM.EQ.22.OR.KFSM.EQ.23.OR.KFSM.EQ.25.OR.KFSM.EQ.35) THEN
C...WHICH NEUTRALINO ?
          NCHI=1
          IF (KFSM.EQ.23) NCHI=2
          IF (KFSM.EQ.25) NCHI=3
          IF (KFSM.EQ.35) NCHI=4
C...SIGN OF MASS (Opposite convention as HERWIG)
          ISM = 1
          IF (SMZ(NCHI).LT.0D0) ISM = -ISM
 
C...Useful parameters for the calculation of the A and B constants.
          WMASS = PMAS(PYCOMP(24),1)
          ECHG = 2*SQRT(PARU(103)*PARU(1))
          COSB=1/(SQRT(1+RMSS(5)**2))
          SINB=RMSS(5)/SQRT(1+RMSS(5)**2)
          COSW=SQRT(1-PARU(102))
          SINW=SQRT(PARU(102))
          GW=2D0*SQRT(PARU(103)*PARU(1))/SINW
C...Run quark masses to neutralino mass squared (for Higgs-type
C...couplings)
          SQMCHI=PMAS(PYCOMP(KFIN),1)**2
          DO 100 I=1,6
            RMQ(I)=PYMRUN(I,SQMCHI)
  100     CONTINUE
C...EXPRESS NEUTRALINO MIXING IN (photino,Zino,~H_u,~H_d) BASIS
            DO 110 NCHJ=1,4
              ZPMIX(NCHJ,1)= ZMIX(NCHJ,1)*COSW+ZMIX(NCHJ,2)*SINW
              ZPMIX(NCHJ,2)=-ZMIX(NCHJ,1)*SINW+ZMIX(NCHJ,2)*COSW
              ZPMIX(NCHJ,3)= ZMIX(NCHJ,3)
              ZPMIX(NCHJ,4)= ZMIX(NCHJ,4)
  110       CONTINUE
            C1=GW*ZPMIX(NCHI,3)/(2D0*COSB*WMASS)
            C1U=GW*ZPMIX(NCHI,4)/(2D0*SINB*WMASS)
            C2=ECHG*ZPMIX(NCHI,1)
            C3=GW*ZPMIX(NCHI,2)/COSW
            EU=2D0/3D0
            ED=-1D0/3D0
C... AB(x,y,z):
C       x=1-2  : Select A or B constant     (1:A ; 2:B)
C       y=1-16 : Sparticle's SM code (1-6:d,u,s,c,b,t ;
C                                    11-16:e,nu_e,mu,...)
C       z=1-2  : Mass eigenstate number
C...CALCULATE COUPLINGS
          DO 120 I = 11,15,2
            CMS=PMAS(PYCOMP(I),1)
C...Intermediate sleptons
            AB(1,I,1)=ISM*(CMS*C1*SFMIX(I,1) + SFMIX(I,2)
     &           *(C2-C3*SINW**2))
            AB(1,I,2)=ISM*(CMS*C1*SFMIX(I,3) + SFMIX(I,4)
     &           *(C2-C3*SINW**2))
            AB(2,I,1)= CMS*C1*SFMIX(I,2) - SFMIX(I,1)*(C2+C3*(5D-1-SINW
     &           **2))
            AB(2,I,2)=CMS*C1*SFMIX(I,4) - SFMIX(I,3)*(C2+C3*(5D-1-SINW
     &           **2))
C...Inermediate sneutrinos
            AB(1,I+1,1)=0D0
            AB(2,I+1,1)=5D-1*C3
            AB(1,I+1,2)=0D0
            AB(2,I+1,2)=0D0
C...Inermediate sdown
            J=I-10
            CMS=RMQ(J)
            AB(1,J,1)=ISM*(CMS*C1*SFMIX(J,1) - SFMIX(J,2)
     &           *ED*(C2-C3*SINW**2))
            AB(1,J,2)=ISM*(CMS*C1*SFMIX(J,3) - SFMIX(J,4)
     &           *ED*(C2-C3*SINW**2))
            AB(2,J,1)=CMS*C1*SFMIX(J,2) + SFMIX(J,1)
     &           *(ED*C2-C3*(1D0/2D0+ED*SINW**2))
            AB(2,J,2)=CMS*C1*SFMIX(J,4) + SFMIX(J,3)
     &           *(ED*C2-C3*(1D0/2D0+ED*SINW**2))
C...Inermediate sup
            J=J+1
            CMS=RMQ(J)
            AB(1,J,1)=ISM*(CMS*C1U*SFMIX(J,1) - SFMIX(J,2)
     &           *EU*(C2-C3*SINW**2))
            AB(1,J,2)=ISM*(CMS*C1U*SFMIX(J,3) - SFMIX(J,4)
     &           *EU*(C2-C3*SINW**2))
            AB(2,J,1)=CMS*C1U*SFMIX(J,2) + SFMIX(J,1)
     &           *(EU*C2+C3*(1D0/2D0-EU*SINW**2))
            AB(2,J,2)=CMS*C1U*SFMIX(J,4) + SFMIX(J,3)
     &           *(EU*C2+C3*(1D0/2D0-EU*SINW**2))
  120     CONTINUE
 
          IF (IMSS(51).GE.1) THEN
C...LAMBDA COUPLINGS (LLE TYPE R-VIOLATION)
C * CHI0_I -> NUBAR_I + LEPTON+_J + lEPTON-_K.
C...STEP IN I,J,K USING SINGLE COUNTER
            DO 130 ISC=0,26
C...LAMBDA COUPLING ASYM IN I,J
              IF(MOD(ISC/9,3).NE.MOD(ISC/3,3)) THEN
                LKNT = LKNT+1
                IDLAM(LKNT,1) =-12 -2*MOD(ISC/9,3)
                IDLAM(LKNT,2) =-11 -2*MOD(ISC/3,3)
                IDLAM(LKNT,3) = 11 +2*MOD(ISC,3)
                XLAM(LKNT)    = 0D0
C...Set coupling, and decay product masses on/off
                RVLAMC        = RVLAM(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1
     &               ,MOD(ISC,3)+1)**2
                DCMASS=.FALSE.
                IF (IDLAM(LKNT,2).EQ.-15.OR.IDLAM(LKNT,3).EQ.15)
     &               DCMASS = .TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1)=-IDLAM(LKNT,1)
                KFR(2)=-IDLAM(LKNT,2)
                KFR(3)=-IDLAM(LKNT,3)
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XLAM(LKNT))
                XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...Charge conjugate mode.
                LKNT=LKNT+1
                IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
                IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
                IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
                XLAM(LKNT)=XLAM(LKNT-1)
C...KINEMATICS CHECK
                IF (XLAM(LKNT).EQ.0D0) THEN
                  LKNT=LKNT-2
                ENDIF
              ENDIF
  130       CONTINUE
          ENDIF
 
          IF (IMSS(52).GE.1) THEN
C...LAMBDA' COUPLINGS. (LQD TYPE R-VIOLATION)
C * CHI0 -> NUBAR_I + DBAR_J + D_K
            DO 140 ISC=0,26
              LKNT = LKNT+1
              IDLAM(LKNT,1) =-12 -2*MOD(ISC/9,3)
              IDLAM(LKNT,2) = -1 -2*MOD(ISC/3,3)
              IDLAM(LKNT,3) =  1 +2*MOD(ISC,3)
              XLAM(LKNT)    =  0D0
C...Set coupling, and decay product masses on/off
              RVLAMC        = 3 * RVLAMP(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1
     &             ,MOD(ISC,3)+1)**2
              DCMASS=.FALSE.
              IF (IDLAM(LKNT,2).EQ.-5.OR.IDLAM(LKNT,3).EQ.5)
     &             DCMASS = .TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)=-IDLAM(LKNT,1)
              KFR(2)=-IDLAM(LKNT,2)
              KFR(3)=-IDLAM(LKNT,3)
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XLAM(LKNT))
              XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...Charge conjugate mode.
              LKNT=LKNT+1
              IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
              IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
              IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
              XLAM(LKNT)=XLAM(LKNT-1)
C...KINEMATICS CHECK
              IF (XLAM(LKNT).EQ.0D0) THEN
                LKNT=LKNT-2
              ENDIF
 
C * CHI0 -> LEPTON_I+ + UBAR_J + D_K
              LKNT = LKNT+1
              IDLAM(LKNT,1) =-11 -2*MOD(ISC/9,3)
              IDLAM(LKNT,2) = -2 -2*MOD(ISC/3,3)
              IDLAM(LKNT,3) =  1 +2*MOD(ISC,3)
              XLAM(LKNT)    =  0D0
C...Set coupling, and decay product masses on/off
              RVLAMC        = 3 * RVLAMP(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1
     &             ,MOD(ISC,3)+1)**2
              DCMASS=.FALSE.
              IF (IDLAM(LKNT,1).EQ.-15.OR.IDLAM(LKNT,2).EQ.-6
     &             .OR.IDLAM(LKNT,3).EQ.5) DCMASS=.TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)=-IDLAM(LKNT,1)
              KFR(2)=-IDLAM(LKNT,2)
              KFR(3)=-IDLAM(LKNT,3)
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XLAM(LKNT))
              XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...Charge conjugate mode.
              LKNT=LKNT+1
              IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
              IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
              IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
              XLAM(LKNT)=XLAM(LKNT-1)
C...KINEMATICS CHECK
              IF (XLAM(LKNT).EQ.0D0) THEN
                LKNT=LKNT-2
              ENDIF
  140       CONTINUE
          ENDIF
 
          IF (IMSS(53).GE.1) THEN
C...LAMBDA'' COUPLINGS. (UDD TYPE R-VIOLATION)
C * CHI0 -> UBAR_I + DBAR_J + DBAR_K
            DO 150 ISC=0,26
C...Symmetry J<->K. Also, LAMB antisymmetric in J and K, so no J=K.
              IF (MOD(ISC/3,3).LT.MOD(ISC,3)) THEN
                LKNT = LKNT+1
                IDLAM(LKNT,1) = -2 -2*MOD(ISC/9,3)
                IDLAM(LKNT,2) = -1 -2*MOD(ISC/3,3)
                IDLAM(LKNT,3) = -1 -2*MOD(ISC,3)
                XLAM(LKNT)    =  0D0
C...Set coupling, and decay product masses on/off
                RVLAMC        = 6. * RVLAMB(MOD(ISC/9,3)+1,MOD(ISC/3,3)
     &               +1,MOD(ISC,3)+1)**2
                DCMASS=.FALSE.
                IF (IDLAM(LKNT,1).EQ.-6.OR.IDLAM(LKNT,2).EQ.-5
     &               .OR.IDLAM(LKNT,3).EQ.-5) DCMASS=.TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
                KFR(1) = IDLAM(LKNT,1)
                KFR(2) = IDLAM(LKNT,2)
                KFR(3) = IDLAM(LKNT,3)
C...Calculate width.
                CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),
     &               IDLAM(LKNT,3),XLAM(LKNT))
                XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...Charge conjugate mode.
                LKNT=LKNT+1
                IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
                IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
                IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
                XLAM(LKNT)=XLAM(LKNT-1)
C...KINEMATICS CHECK
                IF (XLAM(LKNT).EQ.0D0) THEN
                  LKNT=LKNT-2
                ENDIF
              ENDIF
  150       CONTINUE
          ENDIF
        ENDIF
      ENDIF
 
      RETURN
      END
