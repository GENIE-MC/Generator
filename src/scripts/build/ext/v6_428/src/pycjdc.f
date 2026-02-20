 
C*********************************************************************
 
C...PYCJDC
C...Calculate decay widths for the charginos (admixtures of
C...charged Wino and charged Higgsino.
 
C...Input:  KCIN = KF code for particle
C...Output: XLAM = widths
C...        IDLAM = KF codes for decay particles
C...        IKNT = number of decay channels defined
C...AUTHOR: STEPHEN MRENNA
C...Last change:
C...10-16-95:  force decay chi^+_1 -> chi^0_1 e+ nu_e
C...when CHIENU .NE. 0
 
      SUBROUTINE PYCJDC(KFIN,XLAM,IDLAM,IKNT)
 
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
      COMPLEX*16 ZMIXC(4,4),VMIXC(2,2),UMIXC(2,2),OLPP,ORPP
      COMPLEX*16 CAL,CBL,CAR,CBR,CA,CB
      INTEGER KFIN,KCIN
      DOUBLE PRECISION XMI,XMJ,XMF,XMSF1,XMSF2,XMW,XMW2,
     &XMZ,XMZ2,AXMJ,AXMI
      DOUBLE PRECISION S12MIN,S12MAX
      DOUBLE PRECISION XMI2,XMI3,XMJ2,XMH,XMH2,XMHP,XMA2,XMB2,XMK
      DOUBLE PRECISION PYLAMF,XL
      DOUBLE PRECISION TANW,XW,AEM,C1,AS,EI,T3I,BETA,ALFA
      DOUBLE PRECISION PYX2XH,PYX2XG
      DOUBLE PRECISION XLAM(0:400)
      INTEGER IDLAM(400,3)
      INTEGER LKNT,IX,IH,J,IJ,I,IKNT
      INTEGER ITH(3)
      INTEGER ITHC
      DOUBLE PRECISION ETAH(3),DH(3),EH(3)
      DOUBLE PRECISION SR2
      DOUBLE PRECISION CBETA,SBETA,TANB
 
      DOUBLE PRECISION PYALEM,PI,PYALPS
      DOUBLE PRECISION FCOL
      INTEGER KF1,KF2,ISF
      INTEGER KFNCHI(4),KFCCHI(2)
 
      DOUBLE PRECISION TEMP
      EXTERNAL PYGAUS,PYXXZ6
      DOUBLE PRECISION PYGAUS,PYXXZ6
      DOUBLE PRECISION PREC
      DATA ITH/25,35,36/
      DATA ITHC/37/
      DATA ETAH/1D0,1D0,-1D0/
      DATA SR2/1.4142136D0/
      DATA PI/3.141592654D0/
      DATA PREC/1D-2/
      DATA KFNCHI/1000022,1000023,1000025,1000035/
      DATA KFCCHI/1000024,1000037/
 
C...COUNT THE NUMBER OF DECAY MODES
      LKNT=0
      XMW=PMAS(24,1)
      XMW2=XMW**2
      XMZ=PMAS(23,1)
      XMZ2=XMZ**2
      XW=1D0-XMW2/XMZ2
      XW1=1D0-XW
      TANW = SQRT(XW/XW1)
 
C...1 OR 2 DEPENDING ON CHARGINO TYPE
      IX=1
      IF(KFIN.EQ.KFCCHI(2)) IX=2
      KCIN=PYCOMP(KFIN)
 
      XMI=SMW(IX)
      XMI2=XMI**2
      AXMI=ABS(XMI)
      AEM=PYALEM(XMI2)
      AS =PYALPS(XMI2)
      C1=AEM/XW
      XMI3=ABS(XMI**3)
      TANB=RMSS(5)
      BETA=ATAN(TANB)
      CBETA=COS(BETA)
      SBETA=TANB*CBETA
      ALFA=RMSS(18)
 
      DO 110 I=1,2
        DO 100 J=1,2
          VMIXC(J,I)=DCMPLX(VMIX(J,I),VMIXI(J,I))
          UMIXC(J,I)=DCMPLX(UMIX(J,I),UMIXI(J,I))
  100   CONTINUE
  110 CONTINUE
 
C...GRAVITINO DECAY MODES
 
      IF(IMSS(11).EQ.1) THEN
        XMP=RMSS(29)
        IDG=39+KSUSY1
        XMGR=PMAS(PYCOMP(IDG),1)
C        SINW=SQRT(XW)
C        COSW=SQRT(1D0-XW)
        XFAC=(XMI2/(XMP*XMGR))**2*AXMI/48D0/PI
        IF(AXMI.GT.XMGR+XMW) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=24
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*(
     &  .5D0*(ABS(VMIXC(IX,1))**2+ABS(UMIXC(IX,1))**2)+
     &  .5D0*((ABS(VMIXC(IX,2))*SBETA)**2+(ABS(UMIXC(IX,2))*CBETA)**2))*
     &  (1D0-XMW2/XMI2)**4
        ENDIF
        IF(AXMI.GT.XMGR+PMAS(37,1)) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=37
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*(.5D0*((ABS(VMIXC(IX,2))*CBETA)**2+
     &   (ABS(UMIXC(IX,2))*SBETA)**2))
     &   *(1D0-PMAS(37,1)**2/XMI2)**4
       ENDIF
      ENDIF
 
C...CHECK ALL 2-BODY DECAYS TO GAUGE AND HIGGS BOSONS
      IF(IX.EQ.1) GOTO 170
      XMJ=SMW(1)
      AXMJ=ABS(XMJ)
      XMJ2=XMJ**2
 
C...CHI_2+ -> CHI_1+ + Z0
      IF(AXMI.GE.AXMJ+XMZ) THEN
        LKNT=LKNT+1
        IJ=1
        OLPP=-VMIXC(IJ,1)*DCONJG(VMIXC(IX,1))-
     &  VMIXC(IJ,2)*DCONJG(VMIXC(IX,2))/2D0
        ORPP=-UMIXC(IX,1)*DCONJG(UMIXC(IJ,1))-
     &  UMIXC(IX,2)*DCONJG(UMIXC(IJ,2))/2D0
        GX2=ABS(OLPP)**2+ABS(ORPP)**2
        GLR=DBLE(OLPP*DCONJG(ORPP))
        XLAM(LKNT)=PYX2XG(C1/XMW2,XMI,XMJ,XMZ,GX2,GLR)
        IDLAM(LKNT,1)=KFCCHI(1)
        IDLAM(LKNT,2)=23
        IDLAM(LKNT,3)=0
 
C...CHARGED LEPTONS
      ELSEIF(AXMI.GE.AXMJ) THEN
        S12MIN=0D0
        S12MAX=(AXMI-AXMJ)**2
        IA=11
        JA=12
        EI=KCHG(IABS(IA),1)/3D0
        T3I=SIGN(1D0,EI+1D-6)/2D0
        XXC(1)=0D0
        XXC(2)=XMJ
        XXC(3)=0D0
        XXC(4)=XMI
        XXC(5)=PMAS(PYCOMP(KSUSY1+JA),1)
        XXC(6)=1D6
        XXC(9)=PMAS(23,1)
        XXC(10)=PMAS(23,2)
        IJ=1
        OLPP=-VMIXC(IJ,1)*DCONJG(VMIXC(IX,1))-
     &  VMIXC(IJ,2)*DCONJG(VMIXC(IX,2))/2D0
        ORPP=-UMIXC(IX,1)*DCONJG(UMIXC(IJ,1))-
     &  UMIXC(IX,2)*DCONJG(UMIXC(IJ,2))/2D0
        CXC(1)=DCMPLX((T3I-XW*EI)/XW/XW1)*ORPP
        CXC(2)=DCMPLX(0D0,0D0)
        CXC(3)=DCMPLX((T3I-XW*EI)/XW/XW1)*OLPP
        CXC(4)=-VMIXC(IJ,1)*DCONJG(VMIXC(IX,1))*DCMPLX(T3I/XW)
        CXC(5)=-DCMPLX(EI/XW1)*ORPP
        CXC(6)=DCMPLX(0D0,0D0)
        CXC(7)=-DCMPLX(EI/XW1)*OLPP
        CXC(8)=DCMPLX(0D0,0D0)
        IF( XXC(5).LT.AXMI ) THEN
          XXC(5)=1D6
        ENDIF
        XXC(7)=XXC(5)
        XXC(8)=XXC(6)
        IF(AXMI.GE.AXMJ+2D0*PMAS(11,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=11
          IDLAM(LKNT,3)=-11
          IF(AXMI.GE.AXMJ+2D0*PMAS(13,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(1)
            IDLAM(LKNT,2)=13
            IDLAM(LKNT,3)=-13
          ENDIF
          IF(AXMI.GE.AXMJ+2D0*PMAS(15,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(1)
            IDLAM(LKNT,2)=15
            IDLAM(LKNT,3)=-15
          ENDIF
        ENDIF
 
C...NEUTRINOS
  120   CONTINUE
        IA=12
        JA=11
        EI=KCHG(IABS(IA),1)/3D0
        T3I=SIGN(1D0,EI+1D-6)/2D0
        XXC(5)=PMAS(PYCOMP(KSUSY1+JA),1)
        XXC(6)=1D6
        CXC(1)=DCMPLX((T3I-XW*EI)/XW/XW1)*ORPP
        CXC(3)=DCMPLX((T3I-XW*EI)/XW/XW1)*OLPP
        CXC(4)=-UMIXC(IJ,1)*DCONJG(UMIXC(IX,1))*DCMPLX(T3I/XW)
        CXC(5)=-DCMPLX(EI/XW1)*ORPP
        CXC(7)=-DCMPLX(EI/XW1)*OLPP
        IF( XXC(5).LT.AXMI ) THEN
          XXC(5)=1D6
        ENDIF
        XXC(7)=XXC(5)
        XXC(8)=XXC(6)
        IF(AXMI.GE.AXMJ+2D0*PMAS(12,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=12
          IDLAM(LKNT,3)=-12
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=14
          IDLAM(LKNT,3)=-14
        ENDIF
        IF(AXMI.GE.AXMJ+2D0*PMAS(16,1)) THEN
          IF(ABS(SFMIX(15,1)).GT.ABS(SFMIX(15,2))) THEN
            XXC(5)=PMAS(PYCOMP(KSUSY1+15),1)
          ELSE
            XXC(5)=PMAS(PYCOMP(KSUSY2+15),1)
          ENDIF
          IF( XXC(5).LT.AXMI ) THEN
            XXC(5)=1D6
          ENDIF
          XXC(7)=XXC(5)
          LKNT=LKNT+1
          XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=16
          IDLAM(LKNT,3)=-16
        ENDIF
 
C...D-TYPE QUARKS
  130   CONTINUE
        IA=1
        JA=2
        EI=KCHG(IABS(IA),1)/3D0
        T3I=SIGN(1D0,EI+1D-6)/2D0
        XXC(5)=PMAS(PYCOMP(KSUSY1+JA),1)
        XXC(6)=1D6
        CXC(1)=DCMPLX((T3I-XW*EI)/XW/XW1)*ORPP
        CXC(2)=DCMPLX(0D0,0D0)
        CXC(3)=DCMPLX((T3I-XW*EI)/XW/XW1)*OLPP
        CXC(4)=-VMIXC(IJ,1)*DCONJG(VMIXC(IX,1))*DCMPLX(T3I/XW)
        CXC(5)=-DCMPLX(EI/XW1)*ORPP
        CXC(6)=DCMPLX(0D0,0D0)
        CXC(7)=-DCMPLX(EI/XW1)*OLPP
        CXC(8)=DCMPLX(0D0,0D0)
        IF( XXC(5).LT.AXMI ) THEN
          XXC(5)=1D6
        ENDIF
        XXC(7)=XXC(5)
        XXC(8)=XXC(6)
        IF(AXMI.GE.AXMJ+2D0*PMAS(1,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=3D0*C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=1
          IDLAM(LKNT,3)=-1
          IF(AXMI.GE.AXMJ+2D0*PMAS(3,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(1)
            IDLAM(LKNT,2)=3
            IDLAM(LKNT,3)=-3
          ENDIF
        ENDIF
        IF(AXMI.GE.AXMJ+2D0*PMAS(5,1)) THEN
          IF(ABS(SFMIX(5,1)).GT.ABS(SFMIX(5,2))) THEN
            XXC(5)=PMAS(PYCOMP(KSUSY1+5),1)
          ELSE
            XXC(5)=PMAS(PYCOMP(KSUSY2+5),1)
          ENDIF
          IF( XXC(5).LT.AXMI ) THEN
            XXC(5)=1D6
          ENDIF
          XXC(7)=XXC(5)
          LKNT=LKNT+1
          XLAM(LKNT)=3D0*C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=5
          IDLAM(LKNT,3)=-5
        ENDIF
 
C...U-TYPE QUARKS
  140   CONTINUE
        IA=2
        JA=1
        EI=KCHG(IABS(IA),1)/3D0
        T3I=SIGN(1D0,EI+1D-6)/2D0
        XXC(5)=PMAS(PYCOMP(KSUSY1+JA),1)
        XXC(6)=1D6
        CXC(1)=DCMPLX((T3I-XW*EI)/XW/XW1)*ORPP
        CXC(2)=DCMPLX(0D0,0D0)
        CXC(3)=DCMPLX((T3I-XW*EI)/XW/XW1)*OLPP
        CXC(4)=-UMIXC(IJ,1)*DCONJG(UMIXC(IX,1))*DCMPLX(T3I/XW)
        CXC(5)=-DCMPLX(EI/XW1)*ORPP
        CXC(6)=DCMPLX(0D0,0D0)
        CXC(7)=-DCMPLX(EI/XW1)*OLPP
        CXC(8)=DCMPLX(0D0,0D0)
        IF( XXC(5).LT.AXMI ) THEN
          XXC(5)=1D6
        ENDIF
        XXC(7)=XXC(5)
        XXC(8)=XXC(6)
        IF(AXMI.GE.AXMJ+2D0*PMAS(2,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=3D0*C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=2
          IDLAM(LKNT,3)=-2
          IF(AXMI.GE.AXMJ+2D0*PMAS(4,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(1)
            IDLAM(LKNT,2)=4
            IDLAM(LKNT,3)=-4
          ENDIF
        ENDIF
  150   CONTINUE
      ENDIF
 
C...CHI_2+ -> CHI_1+ + H0_K
      EH(2)=COS(ALFA)
      EH(1)=SIN(ALFA)
      EH(3)=-SBETA
      DH(2)=-SIN(ALFA)
      DH(1)=COS(ALFA)
      DH(3)=COS(BETA)
      DO 160 IH=1,3
        XMH=PMAS(ITH(IH),1)
        XMH2=XMH**2
C...NO 3-BODY OPTION
        IF(AXMI.GE.AXMJ+XMH) THEN
          LKNT=LKNT+1
          XL=PYLAMF(XMI2,XMJ2,XMH2)
          OLPP=(VMIXC(2,1)*DCONJG(UMIXC(1,2))*EH(IH) -
     &    VMIXC(2,2)*DCONJG(UMIXC(1,1))*DH(IH))/SR2
          ORPP=(DCONJG(VMIXC(1,1))*UMIXC(2,2)*EH(IH) -
     &    DCONJG(VMIXC(1,2))*UMIXC(2,1)*DH(IH))/SR2
          XMK=XMJ*ETAH(IH)
          GX2=ABS(OLPP)**2+ABS(ORPP)**2
          GLR=DBLE(OLPP*DCONJG(ORPP))
          XLAM(LKNT)=PYX2XH(C1,XMI,XMK,XMH,GX2,GLR)
          IDLAM(LKNT,1)=KFCCHI(1)
          IDLAM(LKNT,2)=ITH(IH)
          IDLAM(LKNT,3)=0
        ENDIF
  160 CONTINUE
 
C...CHI1 JUMPS TO HERE
  170 CONTINUE
 
C...CHI+_I -> CHI0_J + W+
      DO 220 IJ=1,4
        XMJ=SMZ(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
        IF(AXMI.GE.AXMJ+XMW) THEN
          LKNT=LKNT+1
          DO 180 I=1,4
            ZMIXC(IJ,I)=DCMPLX(ZMIX(IJ,I),ZMIXI(IJ,I))
  180     CONTINUE
          CXC(1)=(DCONJG(ZMIXC(IJ,2))*VMIXC(IX,1)-
     &    DCONJG(ZMIXC(IJ,4))*VMIXC(IX,2)/SR2)
          CXC(3)=(ZMIXC(IJ,2)*DCONJG(UMIXC(IX,1))+
     &    ZMIXC(IJ,3)*DCONJG(UMIXC(IX,2))/SR2)
          GX2=ABS(CXC(1))**2+ABS(CXC(3))**2
          GLR=DBLE(CXC(1)*DCONJG(CXC(3)))
          XLAM(LKNT)=PYX2XG(C1/XMW2,XMI,XMJ,XMW,GX2,GLR)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=24
          IDLAM(LKNT,3)=0
C...LEPTONS
        ELSEIF(AXMI.GE.AXMJ) THEN
          S12MIN=0D0
          S12MAX=(AXMI-AXMJ)**2
          DO 190 I=1,4
            ZMIXC(IJ,I)=DCMPLX(ZMIX(IJ,I),ZMIXI(IJ,I))
  190     CONTINUE
          CXC(1)=(DCONJG(ZMIXC(IJ,2))*VMIXC(IX,1)-
     &    DCONJG(ZMIXC(IJ,4))*VMIXC(IX,2)/SR2)/SR2
          CXC(3)=(ZMIXC(IJ,2)*DCONJG(UMIXC(IX,1))+
     &    ZMIXC(IJ,3)*DCONJG(UMIXC(IX,2))/SR2)/SR2
          CXC(5)=DCMPLX(0D0,0D0)
          CXC(7)=DCMPLX(0D0,0D0)
          IA=11
          JA=12
          EI=KCHG(IA,1)/3D0
          T3I=SIGN(1D0,EI+1D-6)/2D0
          EJ=KCHG(JA,1)/3D0
          T3J=SIGN(1D0,EJ+1D-6)/2D0
          CXC(2)=VMIXC(IX,1)*DCONJG(ZMIXC(IJ,1)*(EJ-T3J)*
     &    TANW+ZMIXC(IJ,2)*T3J)/SR2
          CXC(4)=-DCONJG(UMIXC(IX,1))*(
     &    ZMIXC(IJ,1)*(EI-T3I)*TANW+ZMIXC(IJ,2)*T3I)/SR2
          CXC(6)=DCMPLX(0D0,0D0)
          CXC(8)=DCMPLX(0D0,0D0)
          XXC(1)=0D0
          XXC(2)=XMJ
          XXC(3)=0D0
          XXC(4)=XMI
          XXC(5)=PMAS(PYCOMP(KSUSY1+JA),1)
          XXC(6)=PMAS(PYCOMP(KSUSY1+IA),1)
          XXC(9)=PMAS(24,1)
          XXC(10)=PMAS(24,2)
CCC          IF( XXC(5).LT.AXMI .AND. XXC(6).LT.AXMI ) GOTO 190
          IF(XXC(5).LT.AXMI) THEN
            XXC(5)=1D6
          ELSEIF(XXC(6).LT.AXMI) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(6)
          XXC(8)=XXC(5)
C...1/(2PI)**3*/(32*M**3)*G^4, G^2/(4*PI)= AEM/XW,
C...--> 1/(16PI)/M**3*(AEM/XW)**2
          IF(AXMI.GE.AXMJ+PMAS(11,1)+PMAS(12,1)) THEN
            LKNT=LKNT+1
            TEMP=PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*TEMP
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=-11
            IDLAM(LKNT,3)=12
C...ONLY DECAY CHI+1 -> E+ NU_E
            IF( IMSS(12).NE. 0 ) GOTO 260
            IF(AXMI.GE.AXMJ+PMAS(13,1)+PMAS(14,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFNCHI(IJ)
              IDLAM(LKNT,2)=-13
              IDLAM(LKNT,3)=14
            ENDIF
          ENDIF
          IF(AXMI.GE.AXMJ+PMAS(15,1)+PMAS(16,1)) THEN
            LKNT=LKNT+1
            IF(ABS(SFMIX(15,1)).GT.ABS(SFMIX(15,2))) THEN
              XXC(6)=PMAS(PYCOMP(KSUSY1+15),1)
            ELSE
              XXC(6)=PMAS(PYCOMP(KSUSY2+15),1)
            ENDIF
            XXC(5)=PMAS(PYCOMP(KSUSY1+16),1)
            IF(XXC(5).LT.AXMI) THEN
              XXC(5)=1D6
            ELSEIF(XXC(6).LT.AXMI) THEN
              XXC(6)=1D6
            ENDIF
            XXC(7)=XXC(6)
            XXC(8)=XXC(5)
            TEMP=PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*TEMP
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=-15
            IDLAM(LKNT,3)=16
          ENDIF
 
C...NOW, DO THE QUARKS
  200     CONTINUE
          IA=1
          JA=2
          EI=KCHG(IA,1)/3D0
          T3I=SIGN(1D0,EI+1D-6)/2D0
          EJ=KCHG(JA,1)/3D0
          T3J=SIGN(1D0,EJ+1D-6)/2D0
          CXC(2)=VMIXC(IX,1)*DCONJG(ZMIXC(IJ,1)*(EJ-T3J)*
     &    TANW+ZMIXC(IJ,2)*T3J)
          CXC(4)=-DCONJG(UMIXC(IX,1))*(
     &    ZMIXC(IJ,1)*(EI-T3I)*TANW+ZMIXC(IJ,2)*T3I)
          XXC(5)=PMAS(PYCOMP(KSUSY1+JA),1)
          XXC(6)=PMAS(PYCOMP(KSUSY1+IA),1)
          IF( XXC(5).LT.AXMI .AND. XXC(6).LT.AXMI ) GOTO 210
          IF(XXC(5).LT.AXMI) THEN
            XXC(5)=1D6
          ENDIF
          IF(XXC(6).LT.AXMI) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(6)
          XXC(8)=XXC(5)
          IF(AXMI.GE.AXMJ+PMAS(1,1)+PMAS(2,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=3D0*C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=-1
            IDLAM(LKNT,3)=2
            IF(AXMI.GE.AXMJ+PMAS(3,1)+PMAS(4,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFNCHI(IJ)
              IDLAM(LKNT,2)=-3
              IDLAM(LKNT,3)=4
            ENDIF
          ENDIF
  210     CONTINUE
        ENDIF
  220 CONTINUE
 
C...CHI+_I -> CHI0_J + H+
      DO 230 IJ=1,4
        XMJ=SMZ(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
        XMHP=PMAS(ITHC,1)
        IF(AXMI.GE.AXMJ+XMHP) THEN
          LKNT=LKNT+1
          OLPP=CBETA*(ZMIXC(IJ,4)*DCONJG(VMIXC(IX,1))+(ZMIXC(IJ,2)+
     &    ZMIXC(IJ,1)*TANW)*DCONJG(VMIXC(IX,2))/SR2)
          ORPP=SBETA*(DCONJG(ZMIXC(IJ,3))*UMIXC(IX,1)-
     &    (DCONJG(ZMIXC(IJ,2))+DCONJG(ZMIXC(IJ,1))*TANW)*
     &    UMIXC(IX,2)/SR2)
          GX2=ABS(OLPP)**2+ABS(ORPP)**2
          GLR=DBLE(OLPP*DCONJG(ORPP))
          XLAM(LKNT)=PYX2XH(C1,XMI,XMJ,XMHP,GX2,GLR)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=ITHC
          IDLAM(LKNT,3)=0
        ELSE
 
        ENDIF
  230 CONTINUE
 
C...2-BODY DECAYS TO FERMION SFERMION
      DO 240 J=1,16
        IF(J.GE.7.AND.J.LE.10) GOTO 240
        IF(MOD(J,2).EQ.0) THEN
          KF1=KSUSY1+J-1
        ELSE
          KF1=KSUSY1+J+1
        ENDIF
        KF2=KF1+KSUSY1
        XMSF1=PMAS(PYCOMP(KF1),1)
        XMSF2=PMAS(PYCOMP(KF2),1)
        XMF=PMAS(J,1)
        IF(J.LE.6) THEN
          FCOL=3D0
        ELSE
          FCOL=1D0
        ENDIF
 
C...U~ D_L
        IF(MOD(J,2).EQ.0) THEN
          XMFP=PMAS(J-1,1)
          CAL=UMIXC(IX,1)
          CBL=-XMF*VMIXC(IX,2)/XMW/SBETA/SR2
          CAR=-XMFP*UMIXC(IX,2)/XMW/CBETA/SR2
          CBR=0D0
          ISF=J-1
        ELSE
          XMFP=PMAS(J+1,1)
          CAL=VMIXC(IX,1)
          CBL=-XMF*UMIXC(IX,2)/XMW/CBETA/SR2
          CBR=0D0
          CAR=-XMFP*VMIXC(IX,2)/XMW/SBETA/SR2
          ISF=J+1
        ENDIF
 
C...~U_L D
        IF(AXMI.GE.XMF+XMSF1) THEN
          LKNT=LKNT+1
          XMA2=XMSF1**2
          XMB2=XMF**2
          XL=PYLAMF(XMI2,XMA2,XMB2)
          CA=CAL*SFMIX(ISF,1)+CAR*SFMIX(ISF,2)
          CB=CBL*SFMIX(ISF,1)+CBR*SFMIX(ISF,2)
          XLAM(LKNT)=FCOL*C1/8D0/XMI3*SQRT(XL)*( (XMI2+XMB2-XMA2)*
     &    (ABS(CA)**2+ABS(CB)**2)+4D0*DBLE(CA*DCONJG(CB))*XMF*XMI)
          IDLAM(LKNT,3)=0
          IF(MOD(J,2).EQ.0) THEN
            IDLAM(LKNT,1)=-KF1
            IDLAM(LKNT,2)=J
          ELSE
            IDLAM(LKNT,1)=KF1
            IDLAM(LKNT,2)=-J
          ENDIF
        ENDIF
 
C...U~ D_R
        IF(AXMI.GE.XMF+XMSF2) THEN
          LKNT=LKNT+1
          XMA2=XMSF2**2
          XMB2=XMF**2
          CA=CAL*SFMIX(ISF,3)+CAR*SFMIX(ISF,4)
          CB=CBL*SFMIX(ISF,3)+CBR*SFMIX(ISF,4)
          XL=PYLAMF(XMI2,XMA2,XMB2)
          XLAM(LKNT)=FCOL*C1/8D0/XMI3*SQRT(XL)*( (XMI2+XMB2-XMA2)*
     &    (ABS(CA)**2+ABS(CB)**2)+4D0*DBLE(CA*DCONJG(CB))*XMF*XMI)
          IDLAM(LKNT,3)=0
          IF(MOD(J,2).EQ.0) THEN
            IDLAM(LKNT,1)=-KF2
            IDLAM(LKNT,2)=J
          ELSE
            IDLAM(LKNT,1)=KF2
            IDLAM(LKNT,2)=-J
          ENDIF
        ENDIF
  240 CONTINUE
 
C...3-BODY DECAY TO Q Q~' GLUINO, ONLY IF IT CANNOT PROCEED THROUGH
C...A 2-BODY -- 2-BODY CHAIN
      XMJ=PMAS(PYCOMP(KSUSY1+21),1)
      IF(AXMI.GE.XMJ) THEN
        AXMJ=ABS(XMJ)
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
        IF( XXC(5).LT.AXMI .OR. XXC(6).LT.AXMI ) GOTO 250
        IF(AXMI.GE.AXMJ+PMAS(1,1)+PMAS(2,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=4D0*C1*AS/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
          IDLAM(LKNT,1)=KSUSY1+21
          IDLAM(LKNT,2)=-1
          IDLAM(LKNT,3)=2
          IF(AXMI.GE.AXMJ+PMAS(3,1)+PMAS(4,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KSUSY1+21
            IDLAM(LKNT,2)=-3
            IDLAM(LKNT,3)=4
          ENDIF
        ENDIF
  250   CONTINUE
      ENDIF
 
C...R-violating decay modes (SKANDS).
      CALL PYRVCH(KFIN,XLAM,IDLAM,LKNT)
 
  260 IKNT=LKNT
      XLAM(0)=0D0
      DO 270 I=1,IKNT
        XLAM(0)=XLAM(0)+XLAM(I)
        IF(XLAM(I).LT.0D0) THEN
          WRITE(MSTU(11),*) ' XLAM(I) = ',XLAM(I),KCIN,
     &    (IDLAM(I,J),J=1,3)
          XLAM(I)=0D0
        ENDIF
  270 CONTINUE
      IF(XLAM(0).EQ.0D0) THEN
        XLAM(0)=1D-6
        WRITE(MSTU(11),*) ' XLAM(0) = ',XLAM(0)
        WRITE(MSTU(11),*) LKNT
        WRITE(MSTU(11),*) (XLAM(J),J=1,LKNT)
      ENDIF
 
      RETURN
      END
