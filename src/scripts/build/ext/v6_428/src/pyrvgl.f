 
C*********************************************************************
 
C...PYRVGL
C...Calculates R-violating gluino decay widths.
C...See BV part of PYRVCH for comments about the way the BV decay width
C...is calculated. Same comments apply here.
C...P. Z. Skands
 
      SUBROUTINE PYRVGL(KFIN,XLAM,IDLAM,LKNT)
 
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
 
C...IF LQD OR UDD TYPE R-VIOLATION ON.
      IF (IMSS(52).GE.1.OR.IMSS(53).GE.1) THEN
        KFSM=KFIN-KSUSY1
 
C... AB(x,y,z):
C       x=1-2  : Select A or B coupling     (1:A ; 2:B)
C       y=1-16 : Sparticle's SM code (1-6:d,u,s,c,b,t ;
C                                    11-16:e,nu_e,mu,... not used here)
C       z=1-2  : Mass eigenstate number
        DO 100 I = 1,6
C...A Couplings
          AB(1,I,1) = SFMIX(I,2)
          AB(1,I,2) = SFMIX(I,4)
C...B Couplings
          AB(2,I,1) = -SFMIX(I,1)
          AB(2,I,2) = -SFMIX(I,3)
  100   CONTINUE
        GSTR2 = 4D0*PARU(1) * PYALPS(PMAS(PYCOMP(KFIN),1)**2)
C...LQD DECAYS.
        IF (IMSS(52).GE.1) THEN
C...STEP IN I,J,K USING SINGLE COUNTER
          DO 120 ISC=0,26
C * GLUINO -> NUBAR_I + DBAR_J + D_K.
            LKNT          = LKNT+1
            IDLAM(LKNT,1) =-12 -2*MOD(ISC/9,3)
            IDLAM(LKNT,2) = -1 -2*MOD(ISC/3,3)
            IDLAM(LKNT,3) =  1 +2*MOD(ISC,3)
            XLAM(LKNT)=0D0
C...Set coupling, and decay product masses on/off
            RVLAMC=RVLAMP(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)+1)**2
     &           * 5D-1 * GSTR2
            DCMASS        = .FALSE.
            IF (IDLAM(LKNT,2).EQ.-5.OR.IDLAM(LKNT,3).EQ.5) DCMASS=.TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
            KFR(1)        = 0
            KFR(2)        = -IDLAM(LKNT,2)
            KFR(3)        = -IDLAM(LKNT,3)
C...Calculate width.
            CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &           ,XLAM(LKNT))
C...Normalize
            XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...Charge conjugate mode.
  110       LKNT          = LKNT+1
            IDLAM(LKNT,1) =-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2) =-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3) =-IDLAM(LKNT-1,3)
            XLAM(LKNT)    = XLAM(LKNT-1)
C...KINEMATICS CHECK
            IF (XLAM(LKNT).EQ.0D0) THEN
              LKNT=LKNT-2
            ENDIF
 
C * GLUINO -> LEPTON+_I + UBAR_J + D_K
            LKNT = LKNT+1
            IDLAM(LKNT,1) =-11 -2*MOD(ISC/9,3)
            IDLAM(LKNT,2) = -2 -2*MOD(ISC/3,3)
            IDLAM(LKNT,3) =  1 +2*MOD(ISC,3)
            XLAM(LKNT)=0D0
C...Set coupling, and decay product masses on/off
            RVLAMC = RVLAMP(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)+1)
     &           **2* 5D-1 * GSTR2
            DCMASS        = .FALSE.
            IF (IDLAM(LKNT,1).EQ.-15.OR.IDLAM(LKNT,2).EQ.-6
     &           .OR.IDLAM(LKNT,3).EQ.5) DCMASS = .TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
            KFR(1)        = 0
            KFR(2)        = -IDLAM(LKNT,2)
            KFR(3)        = -IDLAM(LKNT,3)
C...Calculate width.
            CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &           ,XLAM(LKNT))
            XLAM(LKNT)=XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...Charge conjugate mode.
            LKNT=LKNT+1
            IDLAM(LKNT,1) = -IDLAM(LKNT-1,1)
            IDLAM(LKNT,2) = -IDLAM(LKNT-1,2)
            IDLAM(LKNT,3) = -IDLAM(LKNT-1,3)
            XLAM(LKNT)    =  XLAM(LKNT-1)
C...KINEMATICS CHECK
            IF (XLAM(LKNT).EQ.0D0) THEN
              LKNT=LKNT-2
            ENDIF
 
  120     CONTINUE
        ENDIF
 
C...UDD DECAYS.
        IF (IMSS(53).GE.1) THEN
C...STEP IN I,J,K USING SINGLE COUNTER
          DO 130 ISC=0,26
C * GLUINO -> UBAR_I + DBAR_J + DBAR_K.
            IF (MOD(ISC/3,3).LT.MOD(ISC,3)) THEN
              LKNT          = LKNT+1
              IDLAM(LKNT,1) = -2 -2*MOD(ISC/9,3)
              IDLAM(LKNT,2) = -1 -2*MOD(ISC/3,3)
              IDLAM(LKNT,3) = -1 -2*MOD(ISC,3)
              XLAM(LKNT)=0D0
C...Set coupling, and decay product masses on/off. A factor of 2 for
C...(N_C-1) has been used to cancel a factor 0.5.
              RVLAMC=RVLAMB(MOD(ISC/9,3)+1,MOD(ISC/3,3)+1,MOD(ISC,3)+1)
     &             **2 * GSTR2
              DCMASS        = .FALSE.
              IF (IDLAM(LKNT,1).EQ.-6.OR.IDLAM(LKNT,2).EQ.-5
     &             .OR.IDLAM(LKNT,3).EQ.-5) DCMASS=.TRUE.
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)        = IDLAM(LKNT,1)
              KFR(2)        = 0
              KFR(3)        = 0
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XRESI)
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)        = 0
              KFR(2)        = IDLAM(LKNT,2)
              KFR(3)        = 0
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XRESJ)
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)        = 0
              KFR(2)        = 0
              KFR(3)        = IDLAM(LKNT,3)
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XRESK)
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)        = IDLAM(LKNT,1)
              KFR(2)        = IDLAM(LKNT,2)
              KFR(3)        = 0
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XRESIJ)
C...Calculate interference function. (Factor -1/2 to make up for factor
C...-2 in PYRVGW.
              IF (ABS(XRESI+XRESJ-XRESIJ).GT.1D-4*XRESIJ) THEN
                XRESIJ = 5D-1 * (XRESI+XRESJ-XRESIJ)
              ELSE
                XRESIJ = 0D0
              ENDIF
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)        = 0
              KFR(2)        = IDLAM(LKNT,2)
              KFR(3)        = IDLAM(LKNT,3)
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XRESJK)
              IF (ABS(XRESJ+XRESK-XRESJK).GT.1D-4*XRESJK) THEN
                XRESJK = 5D-1 * (XRESJ+XRESK-XRESJK)
              ELSE
                XRESJK = 0D0
              ENDIF
C...Resonance KF codes (1=I,2=J,3=K)
              KFR(1)        = IDLAM(LKNT,1)
              KFR(2)        = 0
              KFR(3)        = IDLAM(LKNT,3)
C...Calculate width.
              CALL PYRVGW(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),IDLAM(LKNT,3)
     &             ,XRESIK)
              IF (ABS(XRESI+XRESK-XRESIK).GT.1D-4*XRESIK) THEN
                XRESIK = 5D-1 * (XRESI+XRESK-XRESIK)
              ELSE
                XRESIK = 0D0
              ENDIF
C...Calculate total width (factor 1/2 from 1/(N_C-1))
              XLAM(LKNT) = XRESI + XRESJ + XRESK
     &             + 5D-1 * (XRESIJ + XRESIK + XRESJK)
C...Normalize
              XLAM(LKNT) = XLAM(LKNT)*RVLAMC/((2*PARU(1)*RMS(0))**3*32)
C...Charge conjugate mode.
              LKNT          = LKNT+1
              IDLAM(LKNT,1) =-IDLAM(LKNT-1,1)
              IDLAM(LKNT,2) =-IDLAM(LKNT-1,2)
              IDLAM(LKNT,3) =-IDLAM(LKNT-1,3)
              XLAM(LKNT)    = XLAM(LKNT-1)
C...KINEMATICS CHECK
              IF (XLAM(LKNT).EQ.0D0) THEN
                LKNT=LKNT-2
              ENDIF
            ENDIF
  130     CONTINUE
        ENDIF
      ENDIF
      RETURN
      END
