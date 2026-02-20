 
C*********************************************************************
 
C...PYGLUI
C...Calculates gluino decay modes.
 
      SUBROUTINE PYGLUI(KFIN,XLAM,IDLAM,IKNT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
CC     &SFMIX(16,4),
C      COMMON/PYINTS/XXM(20)
      COMPLEX*16 CXC
      COMMON/PYINTC/XXC(10),CXC(8)
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/,/PYSSMT/,/PYINTC/
 
C...Local variables
      COMPLEX*16 ZMIXC(4,4),VMIXC(2,2),UMIXC(2,2),OLPP,ORPP,GLIJ,GRIJ
      DOUBLE PRECISION XMI,XMJ,XMF,AXMJ,AXMI
      DOUBLE PRECISION XMI2,XMI3,XMA2,XMB2,XMFP
      DOUBLE PRECISION PYLAMF,XL
      DOUBLE PRECISION TANW,XW,AEM,C1,AS,S12MAX,S12MIN
      DOUBLE PRECISION CA,CB,AL,AR,BL,BR
      DOUBLE PRECISION XLAM(0:400)
      INTEGER IDLAM(400,3)
      INTEGER LKNT,IX,ILR,I,IKNT,IFL
      DOUBLE PRECISION SR2
      DOUBLE PRECISION GAM
      DOUBLE PRECISION PYALEM,PI,PYALPS,EI,T3I
      EXTERNAL PYGAUS,PYXXZ6
      DOUBLE PRECISION PYGAUS,PYXXZ6
      DOUBLE PRECISION PREC
      INTEGER KFNCHI(4),KFCCHI(2)
      DATA PI/3.141592654D0/
      DATA SR2/1.4142136D0/
      DATA PREC/1D-2/
      DATA KFNCHI/1000022,1000023,1000025,1000035/
      DATA KFCCHI/1000024,1000037/
 
C...COUNT THE NUMBER OF DECAY MODES
      LKNT=0
      IF(KFIN.NE.KSUSY1+21) RETURN
      KCIN=PYCOMP(KFIN)
 
      XW=PARU(102)
      TANW = SQRT(XW/(1D0-XW))
 
      XMI=PMAS(KCIN,1)
      AXMI=ABS(XMI)
      XMI2=XMI**2
      AEM=PYALEM(XMI2)
      AS =PYALPS(XMI2)
      C1=AEM/XW
      XMI3=AXMI**3
 
      XMI=SIGN(XMI,RMSS(3))
 
C...2-BODY DECAYS OF GLUINO -> GRAVITINO GLUON
 
      IF(IMSS(11).EQ.1) THEN
        XMP=RMSS(29)
        IDG=39+KSUSY1
        XMGR=PMAS(PYCOMP(IDG),1)
        XFAC=(XMI2/(XMP*XMGR))**2*AXMI/48D0/PI
        IF(AXMI.GT.XMGR) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=21
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC
        ENDIF
      ENDIF
 
C...2-BODY DECAYS OF GLUINO -> QUARK SQUARK
 
      DO 110 IFL=1,6
        DO 100 ILR=1,2
          XMJ=PMAS(PYCOMP(ILR*KSUSY1+IFL),1)
          AXMJ=ABS(XMJ)
          XMF=PMAS(IFL,1)
          IF(AXMI.GE.AXMJ+XMF) THEN
C...Minus sign difference from gluino-quark-squark feynman rules
            AL=SFMIX(IFL,1)
            BL=-SFMIX(IFL,3)
            AR=SFMIX(IFL,2)
            BR=-SFMIX(IFL,4)
C...F1 -> F CHI
            IF(ILR.EQ.1) THEN
              CA=AL
              CB=BL
C...F2 -> F CHI
            ELSE
              CA=AR
              CB=BR
            ENDIF
            LKNT=LKNT+1
            XMA2=XMJ**2
            XMB2=XMF**2
            XL=PYLAMF(XMI2,XMA2,XMB2)
            XLAM(LKNT)=4D0/8D0*AS/4D0/XMI3*SQRT(XL)*((XMI2+XMB2-XMA2)*
     &      (CA**2+CB**2)-4D0*CA*CB*XMI*XMF)
            IDLAM(LKNT,1)=ILR*KSUSY1+IFL
            IDLAM(LKNT,2)=-IFL
            IDLAM(LKNT,3)=0
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=0
          ENDIF
  100   CONTINUE
  110 CONTINUE
 
C...3-BODY DECAYS TO GAUGINO FERMION-FERMION
C...GLUINO -> NI Q QBAR
      DO 170 IX=1,4
        XMJ=SMZ(IX)
        AXMJ=ABS(XMJ)
        IF(AXMI.GE.AXMJ) THEN
          DO 120 I=1,4
            ZMIXC(IX,I)=DCMPLX(ZMIX(IX,I),ZMIXI(IX,I))
  120     CONTINUE
          OLPP=DCMPLX(COS(RMSS(32)),SIN(RMSS(32)))/SR2
          ORPP=DCONJG(OLPP)
          XXC(1)=0D0
          XXC(2)=XMJ
          XXC(3)=0D0
          XXC(4)=XMI
          IA=1
          XXC(5)=PMAS(PYCOMP(KSUSY1+IA),1)
          XXC(6)=PMAS(PYCOMP(KSUSY2+IA),1)
          XXC(7)=XXC(5)
          XXC(8)=XXC(6)
          XXC(9)=1D6
          XXC(10)=0D0
          EI=KCHG(IA,1)/3D0
          T3I=SIGN(1D0,EI+1D-6)/2D0
          GLIJ=(T3I*ZMIXC(IX,2)-TANW*(T3I-EI)*ZMIXC(IX,1))*OLPP
          GRIJ=ZMIXC(IX,1)*(EI*TANW)*ORPP
          CXC(1)=0D0
          CXC(2)=-GLIJ
          CXC(3)=0D0
          CXC(4)=DCONJG(GLIJ)
          CXC(5)=0D0
          CXC(6)=GRIJ
          CXC(7)=0D0
          CXC(8)=-DCONJG(GRIJ)
          S12MIN=0D0
          S12MAX=(AXMI-AXMJ)**2
          IF( XXC(5).LT.AXMI .OR. XXC(6).LT.AXMI ) GOTO 130
          IF(AXMI.GE.AXMJ+2D0*PMAS(1,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1*AS/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-2)
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=1
            IDLAM(LKNT,3)=-1
          ENDIF
          IF(AXMI.GE.AXMJ+2D0*PMAS(3,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=3
            IDLAM(LKNT,3)=-3
          ENDIF
  130     CONTINUE
          IF(AXMI.GE.AXMJ+2D0*PMAS(5,1)) THEN
            PMOLD=PMAS(PYCOMP(KSUSY1+5),1)
            IF(AXMI.GT.PMAS(PYCOMP(KSUSY2+5),1)+PMAS(5,1)) THEN
              GOTO 140
            ELSEIF(AXMI.GT.PMAS(PYCOMP(KSUSY1+5),1)+PMAS(5,1)) THEN
              PMAS(PYCOMP(KSUSY1+5),1)=100D0*XMI
            ENDIF
            CALL PYTBBN(IX,100,-1D0/3D0,XMI,GAM)
            LKNT=LKNT+1
            XLAM(LKNT)=GAM
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=5
            IDLAM(LKNT,3)=-5
            PMAS(PYCOMP(KSUSY1+5),1)=PMOLD
          ENDIF
C...U-TYPE QUARKS
  140     CONTINUE
          IA=2
          XXC(5)=PMAS(PYCOMP(KSUSY1+IA),1)
          XXC(6)=PMAS(PYCOMP(KSUSY2+IA),1)
C        IF( XXC(5).LT.AXMI .OR. XXC(6).LT.AXMI ) GOTO 290
          XXC(7)=XXC(5)
          XXC(8)=XXC(6)
          EI=KCHG(IA,1)/3D0
          T3I=SIGN(1D0,EI+1D-6)/2D0
          GLIJ=(T3I*ZMIXC(IX,2)-TANW*(T3I-EI)*ZMIXC(IX,1))*OLPP
          GRIJ=ZMIXC(IX,1)*(EI*TANW)*ORPP
          CXC(2)=-GLIJ
          CXC(4)=DCONJG(GLIJ)
          CXC(6)=GRIJ
          CXC(8)=-DCONJG(GRIJ)
          IF( XXC(5).LT.AXMI .OR. XXC(6).LT.AXMI ) GOTO 150
          IF(AXMI.GE.AXMJ+2D0*PMAS(2,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1*AS/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-2)
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=2
            IDLAM(LKNT,3)=-2
          ENDIF
          IF(AXMI.GE.AXMJ+2D0*PMAS(4,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=4
            IDLAM(LKNT,3)=-4
          ENDIF
  150     CONTINUE
C...INCLUDE THE DECAY GLUINO -> NJ + T + T~
C...IF THE DECAY GLUINO -> ST + T CANNOT OCCUR
          XMF=PMAS(6,1)
          IF(AXMI.GE.AXMJ+2D0*XMF) THEN
            PMOLD=PMAS(PYCOMP(KSUSY1+6),1)
            IF(AXMI.GT.PMAS(PYCOMP(KSUSY2+6),1)+XMF) THEN
              GOTO 160
            ELSEIF(AXMI.GT.PMAS(PYCOMP(KSUSY1+6),1)+XMF) THEN
              PMAS(PYCOMP(KSUSY1+6),1)=100D0*XMI
            ENDIF
            CALL PYTBBN(IX,100,2D0/3D0,XMI,GAM)
            LKNT=LKNT+1
            XLAM(LKNT)=GAM
            IDLAM(LKNT,1)=KFNCHI(IX)
            IDLAM(LKNT,2)=6
            IDLAM(LKNT,3)=-6
            PMAS(PYCOMP(KSUSY1+6),1)=PMOLD
          ENDIF
  160     CONTINUE
        ENDIF
  170 CONTINUE
 
C...GLUINO -> CI Q QBAR'
      DO 210 IX=1,2
        XMJ=SMW(IX)
        AXMJ=ABS(XMJ)
        IF(AXMI.GE.AXMJ) THEN
          DO 180 I=1,2
            VMIXC(IX,I)=DCMPLX(VMIX(IX,I),VMIXI(IX,I))
            UMIXC(IX,I)=DCMPLX(UMIX(IX,I),UMIXI(IX,I))
  180     CONTINUE
          S12MIN=0D0
          S12MAX=(AXMI-AXMJ)**2
          XXC(1)=0D0
          XXC(2)=XMJ
          XXC(3)=0D0
          XXC(4)=XMI
          XXC(5)=PMAS(PYCOMP(KSUSY1+1),1)
          XXC(6)=PMAS(PYCOMP(KSUSY1+2),1)
          XXC(9)=1D6
          XXC(10)=0D0
          OLPP=DCMPLX(COS(RMSS(32)),SIN(RMSS(32)))
          ORPP=DCONJG(OLPP)
          CXC(1)=DCMPLX(0D0,0D0)
          CXC(3)=DCMPLX(0D0,0D0)
          CXC(5)=DCMPLX(0D0,0D0)
          CXC(7)=DCMPLX(0D0,0D0)
          CXC(2)=UMIXC(IX,1)*OLPP/SR2
          CXC(4)=-DCONJG(VMIXC(IX,1))*ORPP/SR2
          CXC(6)=DCMPLX(0D0,0D0)
          CXC(8)=DCMPLX(0D0,0D0)
          IF(XXC(5).LT.AXMI) THEN
            XXC(5)=1D6
          ELSEIF(XXC(6).LT.AXMI) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(6)
          XXC(8)=XXC(5)
          IF( XXC(5).LT.AXMI .OR. XXC(6).LT.AXMI ) GOTO 190
          IF(AXMI.GE.AXMJ+PMAS(1,1)+PMAS(2,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=0.5D0*C1*AS/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
            IDLAM(LKNT,1)=KFCCHI(IX)
            IDLAM(LKNT,2)=1
            IDLAM(LKNT,3)=-2
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
          ENDIF
          IF(AXMI.GE.AXMJ+PMAS(3,1)+PMAS(4,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(IX)
            IDLAM(LKNT,2)=3
            IDLAM(LKNT,3)=-4
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
          ENDIF
  190     CONTINUE
 
          XMF=PMAS(6,1)
          XMFP=PMAS(5,1)
          IF(AXMI.GE.AXMJ+XMF+XMFP) THEN
            IF(XMI.GT.MIN(PMAS(PYCOMP(KSUSY1+5),1)+XMFP,
     $      PMAS(PYCOMP(KSUSY2+6),1)+XMF)) GOTO 200
            PMOLT2=PMAS(PYCOMP(KSUSY2+6),1)
            PMOLB2=PMAS(PYCOMP(KSUSY2+5),1)
            PMOLT1=PMAS(PYCOMP(KSUSY1+6),1)
            PMOLB1=PMAS(PYCOMP(KSUSY1+5),1)
            IF(XMI.GT.PMOLT2+XMF) PMAS(PYCOMP(KSUSY2+6),1)=100D0*AXMI
            IF(XMI.GT.PMOLT1+XMF) PMAS(PYCOMP(KSUSY1+6),1)=100D0*AXMI
            IF(XMI.GT.PMOLB2+XMFP) PMAS(PYCOMP(KSUSY2+5),1)=100D0*AXMI
            IF(XMI.GT.PMOLB1+XMFP) PMAS(PYCOMP(KSUSY1+5),1)=100D0*AXMI
            CALL PYTBBC(IX,100,XMI,GAM)
            LKNT=LKNT+1
            XLAM(LKNT)=GAM
            IDLAM(LKNT,1)=KFCCHI(IX)
            IDLAM(LKNT,2)=5
            IDLAM(LKNT,3)=-6
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
            PMAS(PYCOMP(KSUSY2+6),1)=PMOLT2
            PMAS(PYCOMP(KSUSY2+5),1)=PMOLB2
            PMAS(PYCOMP(KSUSY1+6),1)=PMOLT1
            PMAS(PYCOMP(KSUSY1+5),1)=PMOLB1
          ENDIF
  200     CONTINUE
        ENDIF
  210 CONTINUE
 
C...R-parity violating (3-body) decays.
      CALL PYRVGL(KFIN,XLAM,IDLAM,LKNT)
 
      IKNT=LKNT
      XLAM(0)=0D0
      DO 220 I=1,IKNT
        IF(XLAM(I).LT.0D0) XLAM(I)=0D0
        XLAM(0)=XLAM(0)+XLAM(I)
  220 CONTINUE
      IF(XLAM(0).EQ.0D0) XLAM(0)=1D-6
 
      RETURN
      END
