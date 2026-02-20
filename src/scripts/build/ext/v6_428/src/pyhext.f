 
C*********************************************************************
 
C...PYHEXT
C...Calculates the non-standard decay modes of the Higgs boson.
C...
C...Author:  Stephen Mrenna
C...Last Update:  April 2001
C......Allow complex values for Z,U, and V
 
      SUBROUTINE PYHEXT(KFIN,XLAM,IDLAM,IKNT)
 
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
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
      SAVE /PYDAT1/,/PYDAT2/,/PYPARS/,/PYMSSM/,/PYSSMT/
 
C...Local variables.
      COMPLEX*16 ZMIXC(4,4),VMIXC(2,2),UMIXC(2,2),OLPP,ORPP
      COMPLEX*16 QIJ,RIJ,F21K,F12K
      INTEGER KFIN
      DOUBLE PRECISION XMI,XMJ,XMF,XMW,XMW2,XMZ,AXMJ,AXMI
      DOUBLE PRECISION XMI2,XMI3,XMJ2
      DOUBLE PRECISION PYLAMF,XL,CF,EI
      INTEGER IDU,IFL
      DOUBLE PRECISION TANW,XW,AEM,C1,AS
      DOUBLE PRECISION PYH2XX,GHLL,GHRR,GHLR
      DOUBLE PRECISION XLAM(0:400)
      INTEGER IDLAM(400,3)
      INTEGER LKNT,IH,J,IJ,I,IKNT,IK
      INTEGER ITH(4)
      INTEGER KFNCHI(4),KFCCHI(2)
      DOUBLE PRECISION ETAH(3),CH(3),DH(3),EH(3)
      DOUBLE PRECISION SR2
      DOUBLE PRECISION BETA,ALFA
      DOUBLE PRECISION CBETA,SBETA,GR,GL,TANB
      DOUBLE PRECISION PYALEM
      DOUBLE PRECISION AL,AR,ALR
      DOUBLE PRECISION XMK,AXMK,COSA,SINA,CW,XML
      DOUBLE PRECISION XMUZ,ATRIT,ATRIB,ATRIL
      DOUBLE PRECISION XMJL,XMJR,XM1,XM2
      DATA ITH/25,35,36,37/
      DATA ETAH/1D0,1D0,-1D0/
      DATA SR2/1.4142136D0/
      DATA KFNCHI/1000022,1000023,1000025,1000035/
      DATA KFCCHI/1000024,1000037/
 
C...COUNT THE NUMBER OF DECAY MODES
      LKNT=IKNT
 
      XMW=PMAS(24,1)
      XMW2=XMW**2
      XMZ=PMAS(23,1)
      XW=PARU(102)
      TANW = SQRT(XW/(1D0-XW))
      CW=SQRT(1D0-XW)
 
C...1 - 4 DEPENDING ON Higgs species.
      IH=1
      IF(KFIN.EQ.ITH(2)) IH=2
      IF(KFIN.EQ.ITH(3)) IH=3
      IF(KFIN.EQ.ITH(4)) IH=4
 
      XMI=PMAS(KFIN,1)
      XMI2=XMI**2
      AXMI=ABS(XMI)
      AEM=PYALEM(XMI2)
      C1=AEM/XW
      XMI3=ABS(XMI**3)
 
      TANB=RMSS(5)
      BETA=ATAN(TANB)
      CBETA=COS(BETA)
      SBETA=TANB*CBETA
      ALFA=RMSS(18)
      COSA=COS(ALFA)
      SINA=SIN(ALFA)
      ATRIT=RMSS(16)
      ATRIB=RMSS(15)
      ATRIL=RMSS(17)
      XMUZ=-RMSS(4)
 
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
 
 
      IF(IH.EQ.4) GOTO 220
 
C...CHECK ALL 2-BODY DECAYS TO GAUGE AND HIGGS BOSONS
C...H0_K -> CHI0_I + CHI0_J
      EH(2)=SINA
      EH(1)=COSA
      EH(3)=CBETA
      DH(2)=COSA
      DH(1)=-SINA
      DH(3)=SBETA
      DO 150 IJ=1,4
        XMJ=SMZ(IJ)
        AXMJ=ABS(XMJ)
        DO 140 IK=1,IJ
          XMK=SMZ(IK)
          AXMK=ABS(XMK)
          IF(AXMI.GE.AXMJ+AXMK) THEN
            LKNT=LKNT+1
            QIJ=ZMIXC(IK,3)*ZMIXC(IJ,2)+
     &      ZMIXC(IJ,3)*ZMIXC(IK,2)-
     &      TANW*(ZMIXC(IK,3)*ZMIXC(IJ,1)+
     &      ZMIXC(IJ,3)*ZMIXC(IK,1))
            RIJ=ZMIXC(IK,4)*ZMIXC(IJ,2)+
     &      ZMIXC(IJ,4)*ZMIXC(IK,2)-
     &      TANW*(ZMIXC(IK,4)*ZMIXC(IJ,1)+
     &      ZMIXC(IJ,4)*ZMIXC(IK,1))
            F21K=0.5D0*DCONJG(QIJ*DH(IH)-RIJ*EH(IH))
            F12K=0.5D0*(QIJ*DH(IH)-RIJ*EH(IH))
C...SIGN OF MASSES I,J
            XML=XMK*ETAH(IH)
            GX2=ABS(F12K)**2+ABS(F21K)**2
            GLR=DBLE(F12K*DCONJG(F21K))
            XLAM(LKNT)=PYH2XX(C1,XMI,XMJ,XML,GX2,GLR)
            IF(IJ.EQ.IK) XLAM(LKNT)=XLAM(LKNT)*0.5D0
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=KFNCHI(IK)
            IDLAM(LKNT,3)=0
          ENDIF
  140   CONTINUE
  150 CONTINUE
 
C...H0_K -> CHI+_I CHI-_J
      DO 170 IJ=1,2
        XMJ=SMW(IJ)
        AXMJ=ABS(XMJ)
        DO 160 IK=1,2
          XMK=SMW(IK)
          AXMK=ABS(XMK)
          IF(AXMI.GE.AXMJ+AXMK) THEN
            LKNT=LKNT+1
            OLPP=DCONJG(VMIXC(IJ,1)*UMIXC(IK,2)*DH(IH) +
     &      VMIXC(IJ,2)*UMIXC(IK,1)*EH(IH))/SR2
            ORPP=(VMIXC(IK,1)*UMIXC(IJ,2)*DH(IH) +
     &      VMIXC(IK,2)*UMIXC(IJ,1)*EH(IH))/SR2
            GX2=ABS(OLPP)**2+ABS(ORPP)**2
            GLR=DBLE(OLPP*DCONJG(ORPP))
            XML=XMK*ETAH(IH)
            XLAM(LKNT)=PYH2XX(C1,XMI,XMJ,XML,GX2,GLR)
            IDLAM(LKNT,1)=KFCCHI(IJ)
            IDLAM(LKNT,2)=-KFCCHI(IK)
            IDLAM(LKNT,3)=0
          ENDIF
  160   CONTINUE
  170 CONTINUE
 
C...HIGGS TO SFERMION SFERMION
      DO 200 IFL=1,16
        IF(IFL.GE.7.AND.IFL.LE.10) GOTO 200
        IJ=KSUSY1+IFL
        XMJL=PMAS(PYCOMP(IJ),1)
        XMJR=PMAS(PYCOMP(IJ+KSUSY1),1)
        IF(AXMI.GE.2D0*MIN(XMJL,XMJR)) THEN
          XMJ=XMJL
          XMJ2=XMJ**2
          XL=PYLAMF(XMI2,XMJ2,XMJ2)
          XMF=PMAS(IFL,1)
          EI=KCHG(IFL,1)/3D0
          IDU=2-MOD(IFL,2)
 
          IF(IH.EQ.1) THEN
            IF(IDU.EQ.1) THEN
              GHLL=-XMZ/CW*(0.5D0+EI*XW)*SIN(ALFA+BETA)+
     &        XMF**2/XMW*SINA/CBETA
              GHRR=XMZ/CW*(EI*XW)*SIN(ALFA+BETA)+
     &        XMF**2/XMW*SINA/CBETA
              IF(IFL.EQ.5) THEN
                GHLR=-XMF/2D0/XMW/CBETA*(XMUZ*COSA-
     &          ATRIB*SINA)
              ELSEIF(IFL.EQ.15) THEN
                GHLR=-XMF/2D0/XMW/CBETA*(XMUZ*COSA-
     &          ATRIL*SINA)
              ELSE
                GHLR=0D0
              ENDIF
            ELSE
              GHLL=XMZ/CW*(0.5D0-EI*XW)*SIN(ALFA+BETA)-
     &        XMF**2/XMW*COSA/SBETA
              GHRR=XMZ/CW*(EI*XW)*SIN(ALFA+BETA)-
     &        XMF**2/XMW*COSA/SBETA
              IF(IFL.EQ.6) THEN
                GHLR=XMF/2D0/XMW/SBETA*(XMUZ*SINA-
     &          ATRIT*COSA)
              ELSE
                GHLR=0D0
              ENDIF
            ENDIF
 
          ELSEIF(IH.EQ.2) THEN
            IF(IDU.EQ.1) THEN
              GHLL=XMZ/CW*(0.5D0+EI*XW)*COS(ALFA+BETA)-
     &        XMF**2/XMW*COSA/CBETA
              GHRR=-XMZ/CW*(EI*XW)*COS(ALFA+BETA)-
     &        XMF**2/XMW*COSA/CBETA
              IF(IFL.EQ.5) THEN
                GHLR=-XMF/2D0/XMW/CBETA*(XMUZ*SINA+
     &          ATRIB*COSA)
              ELSEIF(IFL.EQ.15) THEN
                GHLR=-XMF/2D0/XMW/CBETA*(XMUZ*SINA+
     &          ATRIL*COSA)
              ELSE
                GHLR=0D0
              ENDIF
            ELSE
              GHLL=-XMZ/CW*(0.5D0-EI*XW)*COS(ALFA+BETA)-
     &        XMF**2/XMW*SINA/SBETA
              GHRR=-XMZ/CW*(EI*XW)*COS(ALFA+BETA)-
     &        XMF**2/XMW*SINA/SBETA
              IF(IFL.EQ.6) THEN
                GHLR=-XMF/2D0/XMW/SBETA*(XMUZ*COSA+
     &          ATRIT*SINA)
              ELSE
                GHLR=0D0
              ENDIF
            ENDIF
 
          ELSEIF(IH.EQ.3) THEN
            GHLL=0D0
            GHRR=0D0
            GHLR=0D0
            IF(IDU.EQ.1) THEN
              IF(IFL.EQ.5) THEN
                GHLR=XMF/2D0/XMW*(ATRIB*TANB-XMUZ)
              ELSEIF(IFL.EQ.15) THEN
                GHLR=XMF/2D0/XMW*(ATRIL*TANB-XMUZ)
              ENDIF
            ELSE
              IF(IFL.EQ.6) THEN
                GHLR=XMF/2D0/XMW*(ATRIT/TANB-XMUZ)
              ENDIF
            ENDIF
          ENDIF
          IF(IH.EQ.3) GOTO 180
 
          AL=SFMIX(IFL,1)**2
          AR=SFMIX(IFL,2)**2
          ALR=SFMIX(IFL,1)*SFMIX(IFL,2)
          IF(IFL.LE.6) THEN
            CF=3D0
          ELSE
            CF=1D0
          ENDIF
 
          IF(AXMI.GE.2D0*XMJ) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*
     &      (GHLL*AL+GHRR*AR
     &      +2D0*GHLR*ALR)**2
            IDLAM(LKNT,1)=IJ
            IDLAM(LKNT,2)=-IJ
            IDLAM(LKNT,3)=0
          ENDIF
 
          IF(AXMI.GE.2D0*XMJR) THEN
            LKNT=LKNT+1
            AL=SFMIX(IFL,3)**2
            AR=SFMIX(IFL,4)**2
            ALR=SFMIX(IFL,3)*SFMIX(IFL,4)
            XMJ=XMJR
            XMJ2=XMJ**2
            XL=PYLAMF(XMI2,XMJ2,XMJ2)
            XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*
     &      (GHLL*AL+GHRR*AR
     &      +2D0*GHLR*ALR)**2
            IDLAM(LKNT,1)=IJ+KSUSY1
            IDLAM(LKNT,2)=-(IJ+KSUSY1)
            IDLAM(LKNT,3)=0
          ENDIF
  180     CONTINUE
 
          IF(AXMI.GE.XMJL+XMJR) THEN
            LKNT=LKNT+1
            AL=SFMIX(IFL,1)*SFMIX(IFL,3)
            AR=SFMIX(IFL,2)*SFMIX(IFL,4)
            ALR=SFMIX(IFL,1)*SFMIX(IFL,4)+SFMIX(IFL,2)*SFMIX(IFL,3)
            XMJ=XMJR
            XMJ2=XMJ**2
            XL=PYLAMF(XMI2,XMJ2,XMJL**2)
            XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*
     &      (GHLL*AL+GHRR*AR)**2
            IDLAM(LKNT,1)=IJ
            IDLAM(LKNT,2)=-(IJ+KSUSY1)
            IDLAM(LKNT,3)=0
            LKNT=LKNT+1
            IDLAM(LKNT,1)=-IJ
            IDLAM(LKNT,2)=IJ+KSUSY1
            IDLAM(LKNT,3)=0
            XLAM(LKNT)=XLAM(LKNT-1)
          ENDIF
        ENDIF
  190   CONTINUE
  200 CONTINUE
  210 CONTINUE
 
      GOTO 270
  220 CONTINUE
 
C...H+ -> CHI+_I + CHI0_J
      DO 240 IJ=1,4
        XMJ=SMZ(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
        DO 230 IK=1,2
          XMK=SMW(IK)
          AXMK=ABS(XMK)
          IF(AXMI.GE.AXMJ+AXMK) THEN
            LKNT=LKNT+1
            OLPP=CBETA*DCONJG(ZMIXC(IJ,4)*VMIXC(IK,1)+(ZMIXC(IJ,2)+
     &      ZMIXC(IJ,1)*TANW)*VMIXC(IK,2)/SR2)
            ORPP=SBETA*(ZMIXC(IJ,3)*UMIXC(IK,1)-
     &      (ZMIXC(IJ,2)+ZMIXC(IJ,1)*TANW)*UMIXC(IK,2)/SR2)
            GX2=ABS(OLPP)**2+ABS(ORPP)**2
            GLR=DBLE(OLPP*DCONJG(ORPP))
            XLAM(LKNT)=PYH2XX(C1,XMI,XMJ,-XMK,GX2,GLR)
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=KFCCHI(IK)
            IDLAM(LKNT,3)=0
          ENDIF
  230   CONTINUE
  240 CONTINUE
 
      GL=-XMW/SR2*(SIN(2D0*BETA)-PMAS(6,1)**2/TANB/XMW2)
      GR=-PMAS(6,1)/SR2/XMW*(XMUZ-ATRIT/TANB)
      AL=0D0
      AR=0D0
      CF=3D0
 
C...H+ -> T_1 B_1~
      XM1=PMAS(PYCOMP(KSUSY1+6),1)
      XM2=PMAS(PYCOMP(KSUSY1+5),1)
      IF(XMI.GE.XM1+XM2) THEN
        XL=PYLAMF(XMI2,XM1**2,XM2**2)
        LKNT=LKNT+1
        XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*
     &  (GL*SFMIX(6,1)*SFMIX(5,1)+GR*SFMIX(6,2)*SFMIX(5,1))**2
        IDLAM(LKNT,1)=KSUSY1+6
        IDLAM(LKNT,2)=-(KSUSY1+5)
        IDLAM(LKNT,3)=0
      ENDIF
 
C...H+ -> T_2 B_1~
      XM1=PMAS(PYCOMP(KSUSY2+6),1)
      XM2=PMAS(PYCOMP(KSUSY1+5),1)
      IF(XMI.GE.XM1+XM2) THEN
        XL=PYLAMF(XMI2,XM1**2,XM2**2)
        LKNT=LKNT+1
        XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*
     &  (GL*SFMIX(6,3)*SFMIX(5,1)+GR*SFMIX(6,4)*SFMIX(5,1))**2
        IDLAM(LKNT,1)=KSUSY2+6
        IDLAM(LKNT,2)=-(KSUSY1+5)
        IDLAM(LKNT,3)=0
      ENDIF
 
C...H+ -> T_1 B_2~
      XM1=PMAS(PYCOMP(KSUSY1+6),1)
      XM2=PMAS(PYCOMP(KSUSY2+5),1)
      IF(XMI.GE.XM1+XM2) THEN
        XL=PYLAMF(XMI2,XM1**2,XM2**2)
        LKNT=LKNT+1
        XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*
     &  (GL*SFMIX(6,1)*SFMIX(5,3)+GR*SFMIX(6,2)*SFMIX(5,3))**2
        IDLAM(LKNT,1)=KSUSY1+6
        IDLAM(LKNT,2)=-(KSUSY2+5)
        IDLAM(LKNT,3)=0
      ENDIF
 
C...H+ -> T_2 B_2~
      XM1=PMAS(PYCOMP(KSUSY2+6),1)
      XM2=PMAS(PYCOMP(KSUSY2+5),1)
      IF(XMI.GE.XM1+XM2) THEN
        XL=PYLAMF(XMI2,XM1**2,XM2**2)
        LKNT=LKNT+1
        XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*
     &  (GL*SFMIX(6,3)*SFMIX(5,3)+GR*SFMIX(6,4)*SFMIX(5,3))**2
        IDLAM(LKNT,1)=KSUSY2+6
        IDLAM(LKNT,2)=-(KSUSY2+5)
        IDLAM(LKNT,3)=0
      ENDIF
 
C...H+ -> UL DL~
      GL=-XMW/SR2*SIN(2D0*BETA)
      DO 250 IJ=1,3,2
        XM1=PMAS(PYCOMP(KSUSY1+IJ),1)
        XM2=PMAS(PYCOMP(KSUSY1+IJ+1),1)
        IF(XMI.GE.XM1+XM2) THEN
          XL=PYLAMF(XMI2,XM1**2,XM2**2)
          LKNT=LKNT+1
          XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*GL**2
          IDLAM(LKNT,1)=-(KSUSY1+IJ)
          IDLAM(LKNT,2)=KSUSY1+IJ+1
          IDLAM(LKNT,3)=0
        ENDIF
  250 CONTINUE
 
C...H+ -> EL~ NUL
      CF=1D0
      DO 260 IJ=11,13,2
        XM1=PMAS(PYCOMP(KSUSY1+IJ),1)
        XM2=PMAS(PYCOMP(KSUSY1+IJ+1),1)
        IF(XMI.GE.XM1+XM2) THEN
          XL=PYLAMF(XMI2,XM1**2,XM2**2)
          LKNT=LKNT+1
          XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*GL**2
          IDLAM(LKNT,1)=-(KSUSY1+IJ)
          IDLAM(LKNT,2)=KSUSY1+IJ+1
          IDLAM(LKNT,3)=0
        ENDIF
  260 CONTINUE
 
C...H+ -> TAU1 NUTAUL
      XM1=PMAS(PYCOMP(KSUSY1+15),1)
      XM2=PMAS(PYCOMP(KSUSY1+16),1)
      IF(XMI.GE.XM1+XM2) THEN
        XL=PYLAMF(XMI2,XM1**2,XM2**2)
        LKNT=LKNT+1
        XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*GL**2*SFMIX(15,1)**2
        IDLAM(LKNT,1)=-(KSUSY1+15)
        IDLAM(LKNT,2)= KSUSY1+16
        IDLAM(LKNT,3)=0
      ENDIF
 
C...H+ -> TAU2 NUTAUL
      XM1=PMAS(PYCOMP(KSUSY2+15),1)
      XM2=PMAS(PYCOMP(KSUSY1+16),1)
      IF(XMI.GE.XM1+XM2) THEN
        XL=PYLAMF(XMI2,XM1**2,XM2**2)
        LKNT=LKNT+1
        XLAM(LKNT)=CF*SQRT(XL)/4D0*C1/XMI3*GL**2*SFMIX(15,3)**2
        IDLAM(LKNT,1)=-(KSUSY2+15)
        IDLAM(LKNT,2)= KSUSY1+16
        IDLAM(LKNT,3)=0
      ENDIF
 
  270 CONTINUE
      IKNT=LKNT
      XLAM(0)=0D0
      DO 280 I=1,IKNT
        IF(XLAM(I).LE.0D0) XLAM(I)=0D0
        XLAM(0)=XLAM(0)+XLAM(I)
  280 CONTINUE
      IF(XLAM(0).EQ.0D0) XLAM(0)=1D-6
 
      RETURN
      END
