 
C*********************************************************************
 
C...PYXDIF
C...Gives the angular orientation of events.
 
      SUBROUTINE PYXDIF(NC,NJET,KFL,ECM,CHI,THE,PHI)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
 
C...Charge. Factors depending on polarization for QED case.
      QF=KCHG(KFL,1)/3D0
      POLL=1D0-PARJ(131)*PARJ(132)
      POLD=PARJ(132)-PARJ(131)
      IF(MSTJ(102).LE.1.OR.MSTJ(109).EQ.1) THEN
        HF1=POLL
        HF2=0D0
        HF3=PARJ(133)**2
        HF4=0D0
 
C...Factors depending on flavour, energy and polarization for QFD case.
      ELSE
        SFF=1D0/(16D0*PARU(102)*(1D0-PARU(102)))
        SFW=ECM**4/((ECM**2-PARJ(123)**2)**2+(PARJ(123)*PARJ(124))**2)
        SFI=SFW*(1D0-(PARJ(123)/ECM)**2)
        AE=-1D0
        VE=4D0*PARU(102)-1D0
        AF=SIGN(1D0,QF)
        VF=AF-4D0*QF*PARU(102)
        HF1=QF**2*POLL-2D0*QF*VF*SFI*SFF*(VE*POLL-AE*POLD)+
     &  (VF**2+AF**2)*SFW*SFF**2*((VE**2+AE**2)*POLL-2D0*VE*AE*POLD)
        HF2=-2D0*QF*AF*SFI*SFF*(AE*POLL-VE*POLD)+2D0*VF*AF*SFW*SFF**2*
     &  (2D0*VE*AE*POLL-(VE**2+AE**2)*POLD)
        HF3=PARJ(133)**2*(QF**2-2D0*QF*VF*SFI*SFF*VE+(VF**2+AF**2)*
     &  SFW*SFF**2*(VE**2-AE**2))
        HF4=-PARJ(133)**2*2D0*QF*VF*SFW*(PARJ(123)*PARJ(124)/ECM**2)*
     &  SFF*AE
      ENDIF
 
C...Mass factor. Differential cross-sections for two-jet events.
      SQ2=SQRT(2D0)
      QME=0D0
      IF(MSTJ(103).GE.4.AND.IABS(MSTJ(101)).LE.1.AND.MSTJ(102).LE.1.AND.
     &MSTJ(109).NE.1) QME=(2D0*PYMASS(KFL)/ECM)**2
      IF(NJET.EQ.2) THEN
        SIGU=4D0*SQRT(1D0-QME)
        SIGL=2D0*QME*SQRT(1D0-QME)
        SIGT=0D0
        SIGI=0D0
        SIGA=0D0
        SIGP=4D0
 
C...Kinematical variables. Reduce four-jet event to three-jet one.
      ELSE
        IF(NJET.EQ.3) THEN
          X1=2D0*P(NC+1,4)/ECM
          X2=2D0*P(NC+3,4)/ECM
        ELSE
          ECMR=P(NC+1,4)+P(NC+4,4)+SQRT((P(NC+2,1)+P(NC+3,1))**2+
     &    (P(NC+2,2)+P(NC+3,2))**2+(P(NC+2,3)+P(NC+3,3))**2)
          X1=2D0*P(NC+1,4)/ECMR
          X2=2D0*P(NC+4,4)/ECMR
        ENDIF
 
C...Differential cross-sections for three-jet (or reduced four-jet).
        XQ=(1D0-X1)/(1D0-X2)
        CT12=(X1*X2-2D0*X1-2D0*X2+2D0+QME)/SQRT((X1**2-QME)*(X2**2-QME))
        ST12=SQRT(1D0-CT12**2)
        IF(MSTJ(109).NE.1) THEN
          SIGU=2D0*X1**2+X2**2*(1D0+CT12**2)-QME*(3D0+CT12**2-X1-X2)-
     &    QME*X1/XQ+0.5D0*QME*((X2**2-QME)*ST12**2-2D0*X2)*XQ
          SIGL=(X2*ST12)**2-QME*(3D0-CT12**2-2.5D0*(X1+X2)+X1*X2+QME)+
     &    0.5D0*QME*(X1**2-X1-QME)/XQ+0.5D0*QME*((X2**2-QME)*CT12**2-
     &    X2)*XQ
          SIGT=0.5D0*(X2**2-QME-0.5D0*QME*(X2**2-QME)/XQ)*ST12**2
          SIGI=((1D0-0.5D0*QME*XQ)*(X2**2-QME)*ST12*CT12+
     &    QME*(1D0-X1-X2+0.5D0*X1*X2+0.5D0*QME)*ST12/CT12)/SQ2
          SIGA=X2**2*ST12/SQ2
          SIGP=2D0*(X1**2-X2**2*CT12)
 
C...Differential cross-sect for scalar gluons (no mass effects).
        ELSE
          X3=2D0-X1-X2
          XT=X2*ST12
          CT13=SQRT(MAX(0D0,1D0-(XT/X3)**2))
          SIGU=(1D0-PARJ(171))*(X3**2-0.5D0*XT**2)+
     &    PARJ(171)*(X3**2-0.5D0*XT**2-4D0*(1D0-X1)*(1D0-X2)**2/X1)
          SIGL=(1D0-PARJ(171))*0.5D0*XT**2+
     &    PARJ(171)*0.5D0*(1D0-X1)**2*XT**2
          SIGT=(1D0-PARJ(171))*0.25D0*XT**2+
     &    PARJ(171)*0.25D0*XT**2*(1D0-2D0*X1)
          SIGI=-(0.5D0/SQ2)*((1D0-PARJ(171))*XT*X3*CT13+
     &    PARJ(171)*XT*((1D0-2D0*X1)*X3*CT13-X1*(X1-X2)))
          SIGA=(0.25D0/SQ2)*XT*(2D0*(1D0-X1)-X1*X3)
          SIGP=X3**2-2D0*(1D0-X1)*(1D0-X2)/X1
        ENDIF
      ENDIF
 
C...Upper bounds for differential cross-section.
      HF1A=ABS(HF1)
      HF2A=ABS(HF2)
      HF3A=ABS(HF3)
      HF4A=ABS(HF4)
      SIGMAX=(2D0*HF1A+HF3A+HF4A)*ABS(SIGU)+2D0*(HF1A+HF3A+HF4A)*
     &ABS(SIGL)+2D0*(HF1A+2D0*HF3A+2D0*HF4A)*ABS(SIGT)+2D0*SQ2*
     &(HF1A+2D0*HF3A+2D0*HF4A)*ABS(SIGI)+4D0*SQ2*HF2A*ABS(SIGA)+
     &2D0*HF2A*ABS(SIGP)
 
C...Generate angular orientation according to differential cross-sect.
  100 CHI=PARU(2)*PYR(0)
      CTHE=2D0*PYR(0)-1D0
      PHI=PARU(2)*PYR(0)
      CCHI=COS(CHI)
      SCHI=SIN(CHI)
      C2CHI=COS(2D0*CHI)
      S2CHI=SIN(2D0*CHI)
      THE=ACOS(CTHE)
      STHE=SIN(THE)
      C2PHI=COS(2D0*(PHI-PARJ(134)))
      S2PHI=SIN(2D0*(PHI-PARJ(134)))
      SIG=((1D0+CTHE**2)*HF1+STHE**2*(C2PHI*HF3-S2PHI*HF4))*SIGU+
     &2D0*(STHE**2*HF1-STHE**2*(C2PHI*HF3-S2PHI*HF4))*SIGL+
     &2D0*(STHE**2*C2CHI*HF1+((1D0+CTHE**2)*C2CHI*C2PHI-2D0*CTHE*S2CHI*
     &S2PHI)*HF3-((1D0+CTHE**2)*C2CHI*S2PHI+2D0*CTHE*S2CHI*C2PHI)*HF4)*
     &SIGT-2D0*SQ2*(2D0*STHE*CTHE*CCHI*HF1-2D0*STHE*(CTHE*CCHI*C2PHI-
     &SCHI*S2PHI)*HF3+2D0*STHE*(CTHE*CCHI*S2PHI+SCHI*C2PHI)*HF4)*SIGI+
     &4D0*SQ2*STHE*CCHI*HF2*SIGA+2D0*CTHE*HF2*SIGP
      IF(SIG.LT.SIGMAX*PYR(0)) GOTO 100
 
      RETURN
      END
