 
C*********************************************************************
 
C...PYXKFL
C...Selects flavour for produced qqbar pair.
 
      SUBROUTINE PYXKFL(KFL,ECM,ECMC,KFLC)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 
C...Calculate maximum weight in QED or QFD case.
      IF(MSTJ(102).LE.1) THEN
        RFMAX=4D0/9D0
      ELSE
        POLL=1D0-PARJ(131)*PARJ(132)
        SFF=1D0/(16D0*PARU(102)*(1D0-PARU(102)))
        SFW=ECMC**4/((ECMC**2-PARJ(123)**2)**2+(PARJ(123)*PARJ(124))**2)
        SFI=SFW*(1D0-(PARJ(123)/ECMC)**2)
        VE=4D0*PARU(102)-1D0
        HF1I=SFI*SFF*(VE*POLL+PARJ(132)-PARJ(131))
        HF1W=SFW*SFF**2*((VE**2+1D0)*POLL+2D0*VE*(PARJ(132)-PARJ(131)))
        RFMAX=MAX(4D0/9D0*POLL-4D0/3D0*(1D0-8D0*PARU(102)/3D0)*HF1I+
     &  ((1D0-8D0*PARU(102)/3D0)**2+1D0)*HF1W,1D0/9D0*POLL+2D0/3D0*
     &  (-1D0+4D0*PARU(102)/3D0)*HF1I+((-1D0+4D0*PARU(102)/3D0)**2+
     &  1D0)*HF1W)
      ENDIF
 
C...Choose flavour. Gives charge and velocity.
      NTRY=0
  100 NTRY=NTRY+1
      IF(NTRY.GT.100) THEN
        CALL PYERRM(14,'(PYXKFL:) caught in an infinite loop')
        KFLC=0
        RETURN
      ENDIF
      KFLC=KFL
      IF(KFL.LE.0) KFLC=1+INT(MSTJ(104)*PYR(0))
      MSTJ(93)=1
      PMQ=PYMASS(KFLC)
      IF(ECM.LT.2D0*PMQ+PARJ(127)) GOTO 100
      QF=KCHG(KFLC,1)/3D0
      VQ=1D0
      IF(MOD(MSTJ(103),2).EQ.1) VQ=SQRT(MAX(0D0,1D0-(2D0*PMQ/ECMC)**2))
 
C...Calculate weight in QED or QFD case.
      IF(MSTJ(102).LE.1) THEN
        RF=QF**2
        RFV=0.5D0*VQ*(3D0-VQ**2)*QF**2
      ELSE
        VF=SIGN(1D0,QF)-4D0*QF*PARU(102)
        RF=QF**2*POLL-2D0*QF*VF*HF1I+(VF**2+1D0)*HF1W
        RFV=0.5D0*VQ*(3D0-VQ**2)*(QF**2*POLL-2D0*QF*VF*HF1I+VF**2*HF1W)+
     &  VQ**3*HF1W
        IF(RFV.GT.0D0) PARJ(171)=MIN(1D0,VQ**3*HF1W/RFV)
      ENDIF
 
C...Weighting or new event (radiative photon). Cross-section update.
      IF(KFL.LE.0.AND.RF.LT.PYR(0)*RFMAX) GOTO 100
      PARJ(158)=PARJ(158)+1D0
      IF(ECMC.LT.2D0*PMQ+PARJ(127).OR.RFV.LT.PYR(0)*RF) KFLC=0
      IF(MSTJ(107).LE.0.AND.KFLC.EQ.0) GOTO 100
      IF(KFLC.NE.0) PARJ(159)=PARJ(159)+1D0
      PARJ(144)=PARJ(157)*PARJ(159)/PARJ(158)
      PARJ(148)=PARJ(144)*86.8D0/ECM**2
 
      RETURN
      END
