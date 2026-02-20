 
C*********************************************************************
 
C...PYXTEE
C...Calculates total cross-section, including initial state
C...radiation effects.
 
      SUBROUTINE PYXTEE(KFL,ECM,XTOT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 
C...Status, (optimized) Q^2 scale, alpha_strong.
      PARJ(151)=ECM
      MSTJ(119)=10*MSTJ(102)+KFL
      IF(MSTJ(111).EQ.0) THEN
        Q2R=ECM**2
      ELSEIF(MSTU(111).EQ.0) THEN
        PARJ(168)=MIN(1D0,MAX(PARJ(128),EXP(-12D0*PARU(1)/
     &  ((33D0-2D0*MSTU(112))*PARU(111)))))
        Q2R=PARJ(168)*ECM**2
      ELSE
        PARJ(168)=MIN(1D0,MAX(PARJ(128),PARU(112)/ECM,
     &  (2D0*PARU(112)/ECM)**2))
        Q2R=PARJ(168)*ECM**2
      ENDIF
      ALSPI=PYALPS(Q2R)/PARU(1)
 
C...QCD corrections factor in R.
      IF(MSTJ(101).EQ.0.OR.MSTJ(109).EQ.1) THEN
        RQCD=1D0
      ELSEIF(IABS(MSTJ(101)).EQ.1.AND.MSTJ(109).EQ.0) THEN
        RQCD=1D0+ALSPI
      ELSEIF(MSTJ(109).EQ.0) THEN
        RQCD=1D0+ALSPI+(1.986D0-0.115D0*MSTU(118))*ALSPI**2
        IF(MSTJ(111).EQ.1) RQCD=MAX(1D0,RQCD+(33D0-2D0*MSTU(112))/12D0*
     &  LOG(PARJ(168))*ALSPI**2)
      ELSEIF(IABS(MSTJ(101)).EQ.1) THEN
        RQCD=1D0+(3D0/4D0)*ALSPI
      ELSE
        RQCD=1D0+(3D0/4D0)*ALSPI-(3D0/32D0+0.519D0*MSTU(118))*ALSPI**2
      ENDIF
 
C...Calculate Z0 width if default value not acceptable.
      IF(MSTJ(102).GE.3) THEN
        RVA=3D0*(3D0+(4D0*PARU(102)-1D0)**2)+6D0*RQCD*(2D0+
     &  (1D0-8D0*PARU(102)/3D0)**2+(4D0*PARU(102)/3D0-1D0)**2)
        DO 100 KFLC=5,6
          VQ=1D0
          IF(MOD(MSTJ(103),2).EQ.1) VQ=SQRT(MAX(0D0,1D0-
     &    (2D0*PYMASS(KFLC)/ ECM)**2))
          IF(KFLC.EQ.5) VF=4D0*PARU(102)/3D0-1D0
          IF(KFLC.EQ.6) VF=1D0-8D0*PARU(102)/3D0
          RVA=RVA+3D0*RQCD*(0.5D0*VQ*(3D0-VQ**2)*VF**2+VQ**3)
  100   CONTINUE
        PARJ(124)=PARU(101)*PARJ(123)*RVA/(48D0*PARU(102)*
     &  (1D0-PARU(102)))
      ENDIF
 
C...Calculate propagator and related constants for QFD case.
      POLL=1D0-PARJ(131)*PARJ(132)
      IF(MSTJ(102).GE.2) THEN
        SFF=1D0/(16D0*PARU(102)*(1D0-PARU(102)))
        SFW=ECM**4/((ECM**2-PARJ(123)**2)**2+(PARJ(123)*PARJ(124))**2)
        SFI=SFW*(1D0-(PARJ(123)/ECM)**2)
        VE=4D0*PARU(102)-1D0
        SF1I=SFF*(VE*POLL+PARJ(132)-PARJ(131))
        SF1W=SFF**2*((VE**2+1D0)*POLL+2D0*VE*(PARJ(132)-PARJ(131)))
        HF1I=SFI*SF1I
        HF1W=SFW*SF1W
      ENDIF
 
C...Loop over different flavours: charge, velocity.
      RTOT=0D0
      RQQ=0D0
      RQV=0D0
      RVA=0D0
      DO 110 KFLC=1,MAX(MSTJ(104),KFL)
        IF(KFL.GT.0.AND.KFLC.NE.KFL) GOTO 110
        MSTJ(93)=1
        PMQ=PYMASS(KFLC)
        IF(ECM.LT.2D0*PMQ+PARJ(127)) GOTO 110
        QF=KCHG(KFLC,1)/3D0
        VQ=1D0
        IF(MOD(MSTJ(103),2).EQ.1) VQ=SQRT(1D0-(2D0*PMQ/ECM)**2)
 
C...Calculate R and sum of charges for QED or QFD case.
        RQQ=RQQ+3D0*QF**2*POLL
        IF(MSTJ(102).LE.1) THEN
          RTOT=RTOT+3D0*0.5D0*VQ*(3D0-VQ**2)*QF**2*POLL
        ELSE
          VF=SIGN(1D0,QF)-4D0*QF*PARU(102)
          RQV=RQV-6D0*QF*VF*SF1I
          RVA=RVA+3D0*(VF**2+1D0)*SF1W
          RTOT=RTOT+3D0*(0.5D0*VQ*(3D0-VQ**2)*(QF**2*POLL-
     &    2D0*QF*VF*HF1I+VF**2*HF1W)+VQ**3*HF1W)
        ENDIF
  110 CONTINUE
      RSUM=RQQ
      IF(MSTJ(102).GE.2) RSUM=RQQ+SFI*RQV+SFW*RVA
 
C...Calculate cross-section, including QCD corrections.
      PARJ(141)=RQQ
      PARJ(142)=RTOT
      PARJ(143)=RTOT*RQCD
      PARJ(144)=PARJ(143)
      PARJ(145)=PARJ(141)*86.8D0/ECM**2
      PARJ(146)=PARJ(142)*86.8D0/ECM**2
      PARJ(147)=PARJ(143)*86.8D0/ECM**2
      PARJ(148)=PARJ(147)
      PARJ(157)=RSUM*RQCD
      PARJ(158)=0D0
      PARJ(159)=0D0
      XTOT=PARJ(147)
      IF(MSTJ(107).LE.0) RETURN
 
C...Virtual cross-section.
      XKL=PARJ(135)
      XKU=MIN(PARJ(136),1D0-(2D0*PARJ(127)/ECM)**2)
      ALE=2D0*LOG(ECM/PYMASS(11))-1D0
      SIGV=ALE/3D0+2D0*LOG(ECM**2/(PYMASS(13)*PYMASS(15)))/3D0-4D0/3D0+
     &1.526D0*LOG(ECM**2/0.932D0)
 
C...Soft and hard radiative cross-section in QED case.
      IF(MSTJ(102).LE.1) THEN
        SIGV=1.5D0*ALE-0.5D0+PARU(1)**2/3D0+2D0*SIGV
        SIGS=ALE*(2D0*LOG(XKL)-LOG(1D0-XKL)-XKL)
        SIGH=ALE*(2D0*LOG(XKU/XKL)-LOG((1D0-XKU)/(1D0-XKL))-(XKU-XKL))
 
C...Soft and hard radiative cross-section in QFD case.
      ELSE
        SZM=1D0-(PARJ(123)/ECM)**2
        SZW=PARJ(123)*PARJ(124)/ECM**2
        PARJ(161)=-RQQ/RSUM
        PARJ(162)=-(RQQ+RQV+RVA)/RSUM
        PARJ(163)=(RQV*(1D0-0.5D0*SZM-SFI)+RVA*(1.5D0-SZM-SFW))/RSUM
        PARJ(164)=(RQV*SZW**2*(1D0-2D0*SFW)+RVA*(2D0*SFI+SZW**2-
     &  4D0+3D0*SZM-SZM**2))/(SZW*RSUM)
        SIGV=1.5D0*ALE-0.5D0+PARU(1)**2/3D0+((2D0*RQQ+SFI*RQV)/
     &  RSUM)*SIGV+(SZW*SFW*RQV/RSUM)*PARU(1)*20D0/9D0
        SIGS=ALE*(2D0*LOG(XKL)+PARJ(161)*LOG(1D0-XKL)+PARJ(162)*XKL+
     &  PARJ(163)*LOG(((XKL-SZM)**2+SZW**2)/(SZM**2+SZW**2))+
     &  PARJ(164)*(ATAN((XKL-SZM)/SZW)-ATAN(-SZM/SZW)))
        SIGH=ALE*(2D0*LOG(XKU/XKL)+PARJ(161)*LOG((1D0-XKU)/
     &  (1D0-XKL))+PARJ(162)*(XKU-XKL)+PARJ(163)*
     &  LOG(((XKU-SZM)**2+SZW**2)/((XKL-SZM)**2+SZW**2))+
     &  PARJ(164)*(ATAN((XKU-SZM)/SZW)-ATAN((XKL-SZM)/SZW)))
      ENDIF
 
C...Total cross-section and fraction of hard photon events.
      PARJ(160)=SIGH/(PARU(1)/PARU(101)+SIGV+SIGS+SIGH)
      PARJ(157)=RSUM*(1D0+(PARU(101)/PARU(1))*(SIGV+SIGS+SIGH))*RQCD
      PARJ(144)=PARJ(157)
      PARJ(148)=PARJ(144)*86.8D0/ECM**2
      XTOT=PARJ(148)
 
      RETURN
      END
