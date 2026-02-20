 
C*********************************************************************
 
C...PYRVGW
C...Generalized Matrix Element for R-Violating 3-body widths.
C...P. Z. Skands
      SUBROUTINE PYRVGW(KFIN,ID1,ID2,ID3,XLAM)
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
      PARAMETER (EPS=1D-4)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYRVNV/AB(2,16,2),RMS(0:3),RES(6,2),INTRES(6,3),IDR,IDR2
     &     ,DCMASS,KFR(3)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     & SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
      DOUBLE PRECISION XLIM(3,3)
      INTEGER KC(0:3), PYCOMP
      LOGICAL DCMASS, DCHECK(6)
      SAVE /PYDAT2/,/PYRVNV/,/PYSSMT/
 
      XLAM   = 0D0
 
      KC(0)  = PYCOMP(KFIN)
      KC(1)  = PYCOMP(ID1)
      KC(2)  = PYCOMP(ID2)
      KC(3)  = PYCOMP(ID3)
      RMS(0) = PMAS(KC(0),1)
      RMS(1) = PYMRUN(ID1,PMAS(KC(1),1)**2)
      RMS(2) = PYMRUN(ID2,PMAS(KC(2),1)**2)
      RMS(3) = PYMRUN(ID3,PMAS(KC(3),1)**2)
C...INITIALIZE OUTER INTEGRATION LIMITS AND KINEMATICS CHECK
      XLIM(1,1)=(RMS(1)+RMS(2))**2
      XLIM(1,2)=(RMS(0)-RMS(3))**2
      XLIM(1,3)=XLIM(1,2)-XLIM(1,1)
      XLIM(2,1)=(RMS(2)+RMS(3))**2
      XLIM(2,2)=(RMS(0)-RMS(1))**2
      XLIM(2,3)=XLIM(2,2)-XLIM(2,1)
      XLIM(3,1)=(RMS(1)+RMS(3))**2
      XLIM(3,2)=(RMS(0)-RMS(2))**2
      XLIM(3,3)=XLIM(3,2)-XLIM(3,1)
C...Check Phase Space
      IF (XLIM(1,3).LT.0D0.OR.XLIM(2,3).LT.0D0.OR.XLIM(3,3).LT.0D0) THEN
        RETURN
      ENDIF
 
C...INITIALIZE RESONANCE INFORMATION
      DO 110 JRES = 1,3
        DO 100 IMASS = 1,2
          IRES = 2*(JRES-1)+IMASS
          INTRES(IRES,1) = 0
          DCHECK(IRES)   =.FALSE.
C...NO RIGHT-HANDED NEUTRINOS
          IF (((IMASS.EQ.2).AND.((IABS(KFR(JRES)).EQ.12).OR
     &         .(IABS(KFR(JRES)).EQ.14).OR.(IABS(KFR(JRES)).EQ.16))).OR
     &         .KFR(JRES).EQ.0) GOTO 100
          RES(IRES,1) = PMAS(PYCOMP(IMASS*KSUSY1+IABS(KFR(JRES))),1)
          RES(IRES,2) = PMAS(PYCOMP(IMASS*KSUSY1+IABS(KFR(JRES))),2)
          INTRES(IRES,1) = IABS(KFR(JRES))
          INTRES(IRES,2) = IMASS
          IF (KFR(JRES).LT.0) INTRES(IRES,3) = 1
          IF (KFR(JRES).GT.0) INTRES(IRES,3) = 0
  100   CONTINUE
  110 CONTINUE
 
C...SUM OVER DIAGRAMS AND INTEGRATE OVER PHASE SPACE
 
C...RESONANCE CONTRIBUTIONS
C...(Only sum contributions where the resonance is off shell).
C...Store whether diagram on/off in DCHECK.
C...LOOP OVER MASS STATES
      DO 120 J=1,2
        IDR=J
        IF(INTRES(IDR,1).NE.0) THEN

        TMIX = SFMIX(INTRES(IDR,1),2*J+INTRES(IDR,3)-1)**2
        IF ((RMS(0).LT.(RMS(1)+RES(IDR,1)).OR.(RES(IDR,1).LT.(RMS(2)
     &       +RMS(3)))).AND.TMIX.GT.EPS.AND.INTRES(IDR,1).NE.0) THEN
          DCHECK(IDR) =.TRUE.
          XLAM = XLAM + TMIX * PYRVI1(2,3,1)
        ENDIF
        ENDIF
 
        IDR=J+2
        IF(INTRES(IDR,1).NE.0) THEN
        TMIX = SFMIX(INTRES(IDR,1),2*J+INTRES(IDR,3)-1)**2
        IF ((RMS(0).LT.(RMS(2)+RES(IDR,1)).OR.(RES(IDR,1).LT.(RMS(1)
     &       +RMS(3)))).AND.TMIX.GT.EPS.AND.INTRES(IDR,1).NE.0) THEN
          DCHECK(IDR) =.TRUE.
          XLAM = XLAM + TMIX * PYRVI1(1,3,2)
        ENDIF
        ENDIF
 
        IDR=J+4
        IF(INTRES(IDR,1).NE.0) THEN
        TMIX = SFMIX(INTRES(IDR,1),2*J+INTRES(IDR,3)-1)**2
        IF ((RMS(0).LT.(RMS(3)+RES(IDR,1)).OR.(RES(IDR,1).LT.(RMS(1)
     &       +RMS(2)))).AND.TMIX.GT.EPS.AND.INTRES(IDR,1).NE.0) THEN
          DCHECK(IDR) =.TRUE.
          XLAM = XLAM + TMIX * PYRVI1(1,2,3)
        ENDIF
        ENDIF
  120 CONTINUE
C... L-R INTERFERENCES
C... (Only add contributions where both contributing diagrams
C... are non-resonant).
      IDR=1
      IF (DCHECK(1).AND.DCHECK(2)) THEN
C...Bug corrected 11/12 2001. Skands.
        XLAM  = XLAM + 2D0 * PYRVI2(2,3,1)
     &     * SFMIX(INTRES(1,1),2+INTRES(1,3)-1)
     &     * SFMIX(INTRES(2,1),4+INTRES(2,3)-1)
      ENDIF
 
      IDR=3
      IF (DCHECK(3).AND.DCHECK(4)) THEN
        XLAM  = XLAM + 2D0 * PYRVI2(1,3,2)
     &     * SFMIX(INTRES(3,1),2+INTRES(3,3)-1)
     &     * SFMIX(INTRES(4,1),4+INTRES(4,3)-1)
      ENDIF
 
      IDR=5
      IF (DCHECK(5).AND.DCHECK(6)) THEN
        XLAM  = XLAM + 2D0 * PYRVI2(1,2,3)
     &     * SFMIX(INTRES(5,1),2+INTRES(5,3)-1)
     &     * SFMIX(INTRES(6,1),4+INTRES(6,3)-1)
      ENDIF
C... TRUE INTERFERENCES
C... (Only add contributions where both contributing diagrams
C... are non-resonant).
      PREF=-2D0
      IF ((KFIN-KSUSY1).EQ.24.OR.(KFIN-KSUSY1).EQ.37) PREF=2D0
      DO 140 IKR1 = 1,2
        DO 130 IKR2 = 1,2
          IDR  = IKR1+2
          IDR2 = IKR2
          IF (DCHECK(IDR).AND.DCHECK(IDR2)) THEN
            XLAM = XLAM + PREF*PYRVI3(1,3,2) *
     &           SFMIX(INTRES(IDR,1),2*IKR1+INTRES(IDR,3)-1)
     &           *SFMIX(INTRES(IDR2,1),2*IKR2+INTRES(IDR2,3)-1)
          ENDIF
 
          IDR  = IKR1+4
          IDR2 = IKR2
          IF (DCHECK(IDR).AND.DCHECK(IDR2)) THEN
            XLAM = XLAM + PREF*PYRVI3(1,2,3) *
     &           SFMIX(INTRES(IDR,1),2*IKR1+INTRES(IDR,3)-1)
     &           *SFMIX(INTRES(IDR2,1),2*IKR2+INTRES(IDR2,3)-1)
          ENDIF
 
          IDR  = IKR1+4
          IDR2 = IKR2+2
          IF (DCHECK(IDR).AND.DCHECK(IDR2)) THEN
            XLAM = XLAM + PREF*PYRVI3(2,1,3) *
     &           SFMIX(INTRES(IDR,1),2*IKR1+INTRES(IDR,3)-1)
     &           *SFMIX(INTRES(IDR2,1),2*IKR2+INTRES(IDR2,3)-1)
          ENDIF
  130   CONTINUE
  140 CONTINUE
 
      RETURN
      END
