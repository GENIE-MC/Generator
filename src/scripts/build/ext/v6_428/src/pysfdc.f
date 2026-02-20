 
C*********************************************************************
 
C...PYSFDC
C...Calculates decays of sfermions.
 
      SUBROUTINE PYSFDC(KFIN,XLAM,IDLAM,IKNT)
 
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
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/,/PYSSMT/
 
C...Local variables.
      COMPLEX*16 ZMIXC(4,4),VMIXC(2,2),UMIXC(2,2)
      COMPLEX*16 CAL,CAR,CBL,CBR,CALP,CARP,CBLP,CBRP,CA,CB
      INTEGER KFIN,KCIN
      DOUBLE PRECISION XMI,XMJ,XMF,XMSF1,XMSF2,XMW,XMW2,XMZ,AXMJ
      DOUBLE PRECISION XMI2,XMI3,XMA2,XMB2,XMFP
      DOUBLE PRECISION PYLAMF,XL
      DOUBLE PRECISION TANW,XW,AEM,C1,AS
      DOUBLE PRECISION AL,AR,BL,BR
      DOUBLE PRECISION CH1,CH2,CH3,CH4
      DOUBLE PRECISION XMBOT,XMTOP
      DOUBLE PRECISION XLAM(0:400)
      INTEGER IDLAM(400,3)
      INTEGER LKNT,IX,ILR,IDU,J,I,IKNT,IFL,II
      DOUBLE PRECISION SR2
      DOUBLE PRECISION CBETA,SBETA
      DOUBLE PRECISION CW
      DOUBLE PRECISION BETA,ALFA,XMU,AT,AB,ATRIT,ATRIB,ATRIL
      DOUBLE PRECISION COSA,SINA,TANB
      DOUBLE PRECISION PYALEM,PI,PYALPS,EI
      DOUBLE PRECISION GHRR,GHLL,GHLR,XMB,BLR
      INTEGER IG,KF1,KF2
      INTEGER IGG(4),KFNCHI(4),KFCCHI(2)
      DATA IGG/23,25,35,36/
      DATA PI/3.141592654D0/
      DATA SR2/1.4142136D0/
      DATA KFNCHI/1000022,1000023,1000025,1000035/
      DATA KFCCHI/1000024,1000037/
 
C...COUNT THE NUMBER OF DECAY MODES
      LKNT=0
 
C...NO NU_R DECAYS
      IF(KFIN.EQ.KSUSY2+12.OR.KFIN.EQ.KSUSY2+14.OR.
     &KFIN.EQ.KSUSY2+16) RETURN
 
      XMW=PMAS(24,1)
      XMW2=XMW**2
      XMZ=PMAS(23,1)
      XW=PARU(102)
      TANW = SQRT(XW/(1D0-XW))
      CW=SQRT(1D0-XW)
 
      DO 110 I=1,4
        DO 100 J=1,4
          ZMIXC(J,I)=DCMPLX(ZMIX(J,I),ZMIXI(J,I))
  100   CONTINUE
  110 CONTINUE
      DO 130 I=1,2
        DO 120 J=1,2
           VMIXC(J,I)=DCMPLX(VMIX(J,I),VMIXI(J,I))
           UMIXC(J,I)=DCMPLX(UMIX(J,I),UMIXI(J,I))
  120   CONTINUE
  130 CONTINUE
 
C...KCIN
      KCIN=PYCOMP(KFIN)
C...ILR is 1 for left and 2 for right.
      ILR=KFIN/KSUSY1
C...IFL is matching non-SUSY flavour.
      IFL=MOD(KFIN,KSUSY1)
C...IDU is weak isospin, 1 for down and 2 for up.
      IDU=2-MOD(IFL,2)
 
      XMI=PMAS(KCIN,1)
      XMI2=XMI**2
      AEM=PYALEM(XMI2)
      AS =PYALPS(XMI2)
      C1=AEM/XW
      XMI3=XMI**3
      EI=KCHG(IFL,1)/3D0
 
      XMBOT=PYMRUN(5,XMI2)
      XMTOP=PYMRUN(6,XMI2)
 
      TANB=RMSS(5)
      BETA=ATAN(TANB)
      ALFA=RMSS(18)
      CBETA=COS(BETA)
      SBETA=TANB*CBETA
      SINA=SIN(ALFA)
      COSA=COS(ALFA)
      XMU=-RMSS(4)
      ATRIT=RMSS(16)
      ATRIB=RMSS(15)
      ATRIL=RMSS(17)
 
C...2-BODY DECAYS OF SFERMION -> GRAVITINO + FERMION
 
      IF(IMSS(11).EQ.1) THEN
        XMP=RMSS(29)
        IDG=39+KSUSY1
        XMGR=PMAS(PYCOMP(IDG),1)
        XFAC=(XMI2/(XMP*XMGR))**2*XMI/48D0/PI
        IF(IFL.EQ.5) THEN
          XMF=XMBOT
        ELSEIF(IFL.EQ.6) THEN
          XMF=XMTOP
        ELSE
          XMF=PMAS(IFL,1)
        ENDIF
        IF(XMI.GT.XMGR+XMF) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=IFL
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*(1D0-XMF**2/XMI2)**4
        ENDIF
      ENDIF
 
C...2-BODY DECAYS OF SFERMION -> FERMION + GAUGE/GAUGINO
 
C...CHARGED DECAYS:
      DO 140 IX=1,2
C...DI -> U CHI1-,CHI2-
        IF(IDU.EQ.1) THEN
          XMFP=PMAS(IFL+1,1)
          XMF =PMAS(IFL,1)
C...UI -> D CHI1+,CHI2+
        ELSE
          XMFP=PMAS(IFL-1,1)
          XMF =PMAS(IFL,1)
        ENDIF
        XMJ=SMW(IX)
        AXMJ=ABS(XMJ)
        IF(XMI.GE.AXMJ+XMFP) THEN
          XMA2=XMJ**2
          XMB2=XMFP**2
          IF(IDU.EQ.2) THEN
            IF(IFL.EQ.6) THEN
              XMFP=XMBOT
              XMF =XMTOP
            ELSEIF(IFL.LT.6) THEN
              XMF=0D0
              XMFP=0D0
            ENDIF
            CBL=VMIXC(IX,1)
            CAL=-XMFP*UMIXC(IX,2)/SR2/XMW/CBETA
            CBR=-XMF*VMIXC(IX,2)/SR2/XMW/SBETA
            CAR=0D0
          ELSE
            IF(IFL.EQ.5) THEN
              XMF =XMBOT
              XMFP=XMTOP
            ELSEIF(IFL.LT.5) THEN
              XMF=0D0
              XMFP=0D0
            ENDIF
            CBL=UMIXC(IX,1)
            CAL=-XMFP*VMIXC(IX,2)/SR2/XMW/SBETA
            CBR=-XMF*UMIXC(IX,2)/SR2/XMW/CBETA
            CAR=0D0
          ENDIF
 
          CALP=SFMIX(IFL,1)*CAL + SFMIX(IFL,2)*CAR
          CBLP=SFMIX(IFL,1)*CBL + SFMIX(IFL,2)*CBR
          CARP=SFMIX(IFL,4)*CAR + SFMIX(IFL,3)*CAL
          CBRP=SFMIX(IFL,4)*CBR + SFMIX(IFL,3)*CBL
          CAL=CALP
          CBL=CBLP
          CAR=CARP
          CBR=CBRP
 
C...F1 -> F` CHI
          IF(ILR.EQ.1) THEN
            CA=CAL
            CB=CBL
C...F2 -> F` CHI
          ELSE
            CA=CAR
            CB=CBR
          ENDIF
          LKNT=LKNT+1
          XL=PYLAMF(XMI2,XMA2,XMB2)
C...SPIN AVERAGE = 1/1 NOT 1/2....NO COLOR ENHANCEMENT
          XLAM(LKNT)=2D0*C1/8D0/XMI3*SQRT(XL)*((XMI2-XMB2-XMA2)*
     &    (ABS(CA)**2+ABS(CB)**2)-4D0*DBLE(CA*DCONJG(CB))*XMJ*XMFP)
          IDLAM(LKNT,3)=0
          IF(IDU.EQ.1) THEN
            IDLAM(LKNT,1)=-KFCCHI(IX)
            IDLAM(LKNT,2)=IFL+1
          ELSE
            IDLAM(LKNT,1)=KFCCHI(IX)
            IDLAM(LKNT,2)=IFL-1
          ENDIF
        ENDIF
  140 CONTINUE
 
C...NEUTRAL DECAYS
      DO 150 IX=1,4
C...DI -> D CHI10
        XMF=PMAS(IFL,1)
        XMJ=SMZ(IX)
        AXMJ=ABS(XMJ)
        IF(XMI.GE.AXMJ+XMF) THEN
          XMA2=XMJ**2
          XMB2=XMF**2
          IF(IDU.EQ.1) THEN
            IF(IFL.EQ.5) THEN
              XMF=XMBOT
            ELSEIF(IFL.LT.5) THEN
              XMF=0D0
            ENDIF
            CBL=-ZMIXC(IX,2)+TANW*ZMIXC(IX,1)*(2D0*EI+1)
            CAL=XMF*ZMIXC(IX,3)/XMW/CBETA
            CAR=-2D0*EI*TANW*ZMIXC(IX,1)
            CBR=CAL
          ELSE
            IF(IFL.EQ.6) THEN
              XMF=XMTOP
            ELSEIF(IFL.LT.5) THEN
              XMF=0D0
            ENDIF
            CBL=ZMIXC(IX,2)+TANW*ZMIXC(IX,1)*(2D0*EI-1)
            CAL=XMF*ZMIXC(IX,4)/XMW/SBETA
            CAR=-2D0*EI*TANW*ZMIXC(IX,1)
            CBR=CAL
          ENDIF
 
          CALP=SFMIX(IFL,1)*CAL + SFMIX(IFL,2)*CAR
          CBLP=SFMIX(IFL,1)*CBL + SFMIX(IFL,2)*CBR
          CARP=SFMIX(IFL,4)*CAR + SFMIX(IFL,3)*CAL
          CBRP=SFMIX(IFL,4)*CBR + SFMIX(IFL,3)*CBL
          CAL=CALP
          CBL=CBLP
          CAR=CARP
          CBR=CBRP
 
C...F1 -> F CHI
          IF(ILR.EQ.1) THEN
            CA=CAL
            CB=CBL
C...F2 -> F CHI
          ELSE
            CA=CAR
            CB=CBR
          ENDIF
          LKNT=LKNT+1
          XL=PYLAMF(XMI2,XMA2,XMB2)
C...SPIN AVERAGE = 1/1 NOT 1/2....NO COLOR ENHANCEMENT
          XLAM(LKNT)=C1/8D0/XMI3*SQRT(XL)*((XMI2-XMB2-XMA2)*
     &    (ABS(CA)**2+ABS(CB)**2)-4D0*DBLE(CA*DCONJG(CB))*XMJ*XMF)
          IDLAM(LKNT,1)=KFNCHI(IX)
          IDLAM(LKNT,2)=IFL
          IDLAM(LKNT,3)=0
        ENDIF
  150 CONTINUE
 
C...2-BODY DECAYS TO SM GAUGE AND HIGGS BOSONS
C...IG=23,25,35,36
      DO 160 II=1,4
        IG=IGG(II)
        IF(ILR.EQ.1) GOTO 160
        XMB=PMAS(IG,1)
        XMSF1=PMAS(PYCOMP(KFIN-KSUSY1),1)
        IF(XMI.LT.XMSF1+XMB) GOTO 160
        IF(IG.EQ.23) THEN
          BL=-SIGN(.5D0,EI)/CW+EI*XW/CW
          BR=EI*XW/CW
          BLR=0D0
        ELSEIF(IG.EQ.25) THEN
          IF(IFL.EQ.5) THEN
            XMF=XMBOT
          ELSEIF(IFL.EQ.6) THEN
            XMF=XMTOP
          ELSEIF(IFL.LT.5) THEN
            XMF=0D0
          ELSE
            XMF=PMAS(IFL,1)
          ENDIF
          IF(IDU.EQ.2) THEN
            GHLL=XMZ/CW*(0.5D0-EI*XW)*(-SIN(ALFA+BETA))+
     &      XMF**2/XMW*COSA/SBETA
            GHRR=XMZ/CW*(EI*XW)*(-SIN(ALFA+BETA))+
     &      XMF**2/XMW*COSA/SBETA
          ELSE
            GHLL=XMZ/CW*(0.5D0-EI*XW)*(-SIN(ALFA+BETA))+
     &      XMF**2/XMW*(-SINA)/CBETA
            GHRR=XMZ/CW*(EI*XW)*(-SIN(ALFA+BETA))+
     &      XMF**2/XMW*(-SINA)/CBETA
          ENDIF
          IF(IFL.EQ.5) THEN
            AT=ATRIB
          ELSEIF(IFL.EQ.6) THEN
            AT=ATRIT
          ELSEIF(IFL.EQ.15) THEN
            AT=ATRIL
          ELSE
            AT=0D0
          ENDIF
C.........need to complexify
          IF(IDU.EQ.2) THEN
            GHLR=XMF/2D0/XMW/SBETA*(-XMU*SINA+
     &      AT*COSA)
          ELSE
            GHLR=XMF/2D0/XMW/CBETA*(XMU*COSA-
     &      AT*SINA)
          ENDIF
          BL=GHLL
          BR=GHRR
          BLR=-GHLR
        ELSEIF(IG.EQ.35) THEN
          IF(IFL.EQ.5) THEN
            XMF=XMBOT
          ELSEIF(IFL.EQ.6) THEN
            XMF=XMTOP
          ELSEIF(IFL.LT.5) THEN
            XMF=0D0
          ELSE
            XMF=PMAS(IFL,1)
          ENDIF
          IF(IDU.EQ.2) THEN
            GHLL=XMZ/CW*(0.5D0-EI*XW)*COS(ALFA+BETA)+
     &      XMF**2/XMW*SINA/SBETA
            GHRR=XMZ/CW*(EI*XW)*COS(ALFA+BETA)+
     &      XMF**2/XMW*SINA/SBETA
          ELSE
            GHLL=XMZ/CW*(0.5D0-EI*XW)*COS(ALFA+BETA)+
     &      XMF**2/XMW*COSA/CBETA
            GHRR=XMZ/CW*(EI*XW)*COS(ALFA+BETA)+
     &      XMF**2/XMW*COSA/CBETA
          ENDIF
          IF(IFL.EQ.5) THEN
            AT=ATRIB
          ELSEIF(IFL.EQ.6) THEN
            AT=ATRIT
          ELSEIF(IFL.EQ.15) THEN
            AT=ATRIL
          ELSE
            AT=0D0
          ENDIF
C.........Need to complexify
          IF(IDU.EQ.2) THEN
            GHLR=XMF/2D0/XMW/SBETA*(XMU*COSA+
     &      AT*SINA)
          ELSE
            GHLR=XMF/2D0/XMW/CBETA*(XMU*SINA+
     &      AT*COSA)
          ENDIF
          BL=GHLL
          BR=GHRR
          BLR=GHLR
        ELSEIF(IG.EQ.36) THEN
          GHLL=0D0
          GHRR=0D0
          IF(IFL.EQ.5) THEN
            XMF=XMBOT
          ELSEIF(IFL.EQ.6) THEN
            XMF=XMTOP
          ELSEIF(IFL.LT.5) THEN
            XMF=0D0
          ELSE
            XMF=PMAS(IFL,1)
          ENDIF
          IF(IFL.EQ.5) THEN
            AT=ATRIB
          ELSEIF(IFL.EQ.6) THEN
            AT=ATRIT
          ELSEIF(IFL.EQ.15) THEN
            AT=ATRIL
          ELSE
            AT=0D0
          ENDIF
C.........Need to complexify
          IF(IDU.EQ.2) THEN
            GHLR=XMF/2D0/XMW*(-XMU+AT/TANB)
          ELSE
            GHLR=XMF/2D0/XMW/(-XMU+AT*TANB)
          ENDIF
          BL=GHLL
          BR=GHRR
          BLR=GHLR
        ENDIF
        AL=SFMIX(IFL,1)*SFMIX(IFL,3)*BL+
     &  SFMIX(IFL,2)*SFMIX(IFL,4)*BR+
     &  (SFMIX(IFL,1)*SFMIX(IFL,4)+SFMIX(IFL,3)*SFMIX(IFL,2))*BLR
        XL=PYLAMF(XMI2,XMSF1**2,XMB**2)
        LKNT=LKNT+1
        IF(IG.EQ.23) THEN
          XLAM(LKNT)=C1/4D0/XMI3*XL**1.5D0/XMB**2*AL**2
        ELSE
          XLAM(LKNT)=C1/4D0/XMI3*SQRT(XL)*AL**2
        ENDIF
        IDLAM(LKNT,3)=0
        IDLAM(LKNT,1)=KFIN-KSUSY1
        IDLAM(LKNT,2)=IG
  160 CONTINUE
 
C...SF -> SF' + W
      XMB=PMAS(24,1)
      IF(MOD(IFL,2).EQ.0) THEN
        KF1=KSUSY1+IFL-1
      ELSE
        KF1=KSUSY1+IFL+1
      ENDIF
      KF2=KF1+KSUSY1
      XMSF1=PMAS(PYCOMP(KF1),1)
      XMSF2=PMAS(PYCOMP(KF2),1)
      IF(XMI.GT.XMB+XMSF1) THEN
        IF(MOD(IFL,2).EQ.0) THEN
          IF(ILR.EQ.1) THEN
            AL=1D0/SR2*SFMIX(IFL,1)*SFMIX(IFL-1,1)
          ELSE
            AL=1D0/SR2*SFMIX(IFL,3)*SFMIX(IFL-1,1)
          ENDIF
        ELSE
          IF(ILR.EQ.1) THEN
            AL=1D0/SR2*SFMIX(IFL,1)*SFMIX(IFL+1,1)
          ELSE
            AL=1D0/SR2*SFMIX(IFL,3)*SFMIX(IFL+1,1)
          ENDIF
        ENDIF
        XL=PYLAMF(XMI2,XMSF1**2,XMB**2)
        LKNT=LKNT+1
        XLAM(LKNT)=C1/4D0/XMI3*XL**1.5D0/XMB**2*AL**2
        IDLAM(LKNT,3)=0
        IDLAM(LKNT,1)=KF1
        IDLAM(LKNT,2)=SIGN(24,KCHG(IFL,1))
      ENDIF
      IF(XMI.GT.XMB+XMSF2) THEN
        IF(MOD(IFL,2).EQ.0) THEN
          IF(ILR.EQ.1) THEN
            AL=1D0/SR2*SFMIX(IFL,1)*SFMIX(IFL-1,3)
          ELSE
            AL=1D0/SR2*SFMIX(IFL,3)*SFMIX(IFL-1,3)
          ENDIF
        ELSE
          IF(ILR.EQ.1) THEN
            AL=1D0/SR2*SFMIX(IFL,1)*SFMIX(IFL+1,3)
          ELSE
            AL=1D0/SR2*SFMIX(IFL,3)*SFMIX(IFL+1,3)
          ENDIF
        ENDIF
        XL=PYLAMF(XMI2,XMSF2**2,XMB**2)
        LKNT=LKNT+1
        XLAM(LKNT)=C1/4D0/XMI3*XL**1.5D0/XMB**2*AL**2
        IDLAM(LKNT,3)=0
        IDLAM(LKNT,1)=KF2
        IDLAM(LKNT,2)=SIGN(24,KCHG(IFL,1))
      ENDIF
 
C...SF -> SF' + HC
      XMB=PMAS(37,1)
      IF(MOD(IFL,2).EQ.0) THEN
        KF1=KSUSY1+IFL-1
      ELSE
        KF1=KSUSY1+IFL+1
      ENDIF
      KF2=KF1+KSUSY1
      XMSF1=PMAS(PYCOMP(KF1),1)
      XMSF2=PMAS(PYCOMP(KF2),1)
      IF(XMI.GT.XMB+XMSF1) THEN
        XMF=0D0
        XMFP=0D0
        AT=0D0
        AB=0D0
        IF(MOD(IFL,2).EQ.0) THEN
C...T1-> B1 HC
          IF(ILR.EQ.1) THEN
            CH1=-SFMIX(IFL,1)*SFMIX(IFL-1,1)
            CH2= SFMIX(IFL,2)*SFMIX(IFL-1,2)
            CH3=-SFMIX(IFL,1)*SFMIX(IFL-1,2)
            CH4=-SFMIX(IFL,2)*SFMIX(IFL-1,1)
C...T2-> B1 HC
          ELSE
            CH1= SFMIX(IFL,3)*SFMIX(IFL-1,1)
            CH2=-SFMIX(IFL,4)*SFMIX(IFL-1,2)
            CH3= SFMIX(IFL,3)*SFMIX(IFL-1,2)
            CH4= SFMIX(IFL,4)*SFMIX(IFL-1,1)
          ENDIF
          IF(IFL.EQ.6) THEN
            XMF=XMTOP
            XMFP=XMBOT
            AT=ATRIT
            AB=ATRIB
          ENDIF
        ELSE
C...B1 -> T1 HC
          IF(ILR.EQ.1) THEN
            CH1=-SFMIX(IFL+1,1)*SFMIX(IFL,1)
            CH2= SFMIX(IFL+1,2)*SFMIX(IFL,2)
            CH3=-SFMIX(IFL+1,1)*SFMIX(IFL,2)
            CH4=-SFMIX(IFL+1,2)*SFMIX(IFL,1)
C...B2-> T1 HC
          ELSE
            CH1= SFMIX(IFL,3)*SFMIX(IFL+1,1)
            CH2=-SFMIX(IFL,4)*SFMIX(IFL+1,2)
            CH3= SFMIX(IFL,4)*SFMIX(IFL+1,1)
            CH4= SFMIX(IFL,3)*SFMIX(IFL+1,2)
          ENDIF
          IF(IFL.EQ.5) THEN
            XMF=XMTOP
            XMFP=XMBOT
            AT=ATRIT
            AB=ATRIB
          ENDIF
        ENDIF
        XL=PYLAMF(XMI2,XMSF1**2,XMB**2)
        LKNT=LKNT+1
C.......Need to complexify
        AL=CH1*(XMW2*2D0*CBETA*SBETA-XMFP**2*TANB-XMF**2/TANB)+
     &  CH2*2D0*XMF*XMFP/(2D0*CBETA*SBETA)+
     &  CH3*XMFP*(-XMU+AB*TANB)+CH4*XMF*(-XMU+AT/TANB)
        XLAM(LKNT)=C1/8D0/XMI3*SQRT(XL)/XMW2*AL**2
        IDLAM(LKNT,3)=0
        IDLAM(LKNT,1)=KF1
        IDLAM(LKNT,2)=SIGN(37,KCHG(IFL,1))
      ENDIF
      IF(XMI.GT.XMB+XMSF2) THEN
        XMF=0D0
        XMFP=0D0
        AT=0D0
        AB=0D0
        IF(MOD(IFL,2).EQ.0) THEN
C...T1-> B2 HC
          IF(ILR.EQ.1) THEN
            CH1= SFMIX(IFL-1,3)*SFMIX(IFL,1)
            CH2=-SFMIX(IFL-1,4)*SFMIX(IFL,2)
            CH3= SFMIX(IFL-1,4)*SFMIX(IFL,1)
            CH4= SFMIX(IFL-1,3)*SFMIX(IFL,2)
C...T2-> B2 HC
          ELSE
            CH1= -SFMIX(IFL,3)*SFMIX(IFL-1,3)
            CH2= SFMIX(IFL,4)*SFMIX(IFL-1,4)
            CH3= -SFMIX(IFL,3)*SFMIX(IFL-1,4)
            CH4= -SFMIX(IFL,4)*SFMIX(IFL-1,3)
          ENDIF
          IF(IFL.EQ.6) THEN
            XMF=XMTOP
            XMFP=XMBOT
            AT=ATRIT
            AB=ATRIB
          ENDIF
        ELSE
C...B1 -> T2 HC
          IF(ILR.EQ.1) THEN
            CH1= SFMIX(IFL+1,3)*SFMIX(IFL,1)
            CH2=-SFMIX(IFL+1,4)*SFMIX(IFL,2)
            CH3= SFMIX(IFL+1,3)*SFMIX(IFL,2)
            CH4= SFMIX(IFL+1,4)*SFMIX(IFL,1)
C...B2-> T2 HC
          ELSE
            CH1= -SFMIX(IFL+1,3)*SFMIX(IFL,3)
            CH2= SFMIX(IFL+1,4)*SFMIX(IFL,4)
            CH3= -SFMIX(IFL+1,3)*SFMIX(IFL,4)
            CH4= -SFMIX(IFL+1,4)*SFMIX(IFL,3)
          ENDIF
          IF(IFL.EQ.5) THEN
            XMF=XMTOP
            XMFP=XMBOT
            AT=ATRIT
            AB=ATRIB
          ENDIF
        ENDIF
        XL=PYLAMF(XMI2,XMSF1**2,XMB**2)
        LKNT=LKNT+1
C.......Need to complexify
        AL=CH1*(XMW2*2D0*CBETA*SBETA-XMFP**2*TANB-XMF**2/TANB)+
     &  CH2*2D0*XMF*XMFP/(2D0*CBETA*SBETA)+
     &  CH3*XMFP*(-XMU+AB*TANB)+CH4*XMF*(-XMU+AT/TANB)
        XLAM(LKNT)=C1/8D0/XMI3*SQRT(XL)/XMW2*AL**2
        IDLAM(LKNT,3)=0
        IDLAM(LKNT,1)=KF2
        IDLAM(LKNT,2)=SIGN(37,KCHG(IFL,1))
      ENDIF
 
C...2-BODY DECAYS OF SQUARK -> QUARK GLUINO
 
      IF(IFL.LE.6) THEN
        XMFP=0D0
        XMF=0D0
        IF(IFL.EQ.6) XMF=PMAS(6,1)
        IF(IFL.EQ.5) XMF=PMAS(5,1)
        XMJ=PMAS(PYCOMP(KSUSY1+21),1)
        AXMJ=ABS(XMJ)
        IF(XMI.GE.AXMJ+XMF) THEN
          AL=-SFMIX(IFL,3)
          BL=SFMIX(IFL,1)
          AR=-SFMIX(IFL,4)
          BR=SFMIX(IFL,2)
C...F1 -> F CHI
          IF(ILR.EQ.1) THEN
            XCA=AL
            XCB=BL
C...F2 -> F CHI
          ELSE
            XCA=AR
            XCB=BR
          ENDIF
          LKNT=LKNT+1
          XMA2=XMJ**2
          XMB2=XMF**2
          XL=PYLAMF(XMI2,XMA2,XMB2)
          XLAM(LKNT)=4D0/3D0*AS/2D0/XMI3*SQRT(XL)*((XMI2-XMB2-XMA2)*
     &    (XCA**2+XCB**2)+4D0*XCA*XCB*XMJ*XMF)
          IDLAM(LKNT,1)=KSUSY1+21
          IDLAM(LKNT,2)=IFL
          IDLAM(LKNT,3)=0
        ENDIF
      ENDIF
 
C...IF NOTHING ELSE FOR T1, THEN T1* -> C+CHI0
      IF(KFIN.EQ.KSUSY1+6.AND.PMAS(KCIN,1).GT.
     &PMAS(PYCOMP(KSUSY1+22),1)+PMAS(4,1)) THEN
C...THIS IS A BACK-OF-THE-ENVELOPE ESTIMATE
C...M = 1/(16PI**2)G**3 = G*2/(4PI) G/(4PI) = C1 * G/(4PI)
C...M*M = C1**2 * G**2/(16PI**2)
C...G = 1/(8PI)P/MI**2 * M*M = C1**3/(32PI**2)*LAM/(2*MI**3)
        LKNT=LKNT+1
        XL=PYLAMF(XMI2,0D0,PMAS(PYCOMP(KSUSY1+22),1)**2)
        XLAM(LKNT)=C1**3/64D0/PI**2/XMI3*SQRT(XL)
        IF(XLAM(LKNT).EQ.0) XLAM(LKNT)=1D-3
        IDLAM(LKNT,1)=KSUSY1+22
        IDLAM(LKNT,2)=4
        IDLAM(LKNT,3)=0
      ENDIF
 
C...R-violating sfermion decays (SKANDS).
      CALL PYRVSF(KFIN,XLAM,IDLAM,LKNT)
 
      IKNT=LKNT
      XLAM(0)=0D0
      DO 170 I=1,IKNT
        IF(XLAM(I).LT.0D0) XLAM(I)=0D0
        XLAM(0)=XLAM(0)+XLAM(I)
  170 CONTINUE
      IF(XLAM(0).EQ.0D0) XLAM(0)=1D-3
 
      RETURN
      END
