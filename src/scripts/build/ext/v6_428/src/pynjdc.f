 
C*********************************************************************
 
C...PYNJDC
C...Calculates decay widths for the neutralinos (admixtures of
C...Bino, W3-ino, Higgs1-ino, Higgs2-ino)
 
C...Input:  KCIN = KF code for particle
C...Output: XLAM = widths
C...        IDLAM = KF codes for decay particles
C...        IKNT = number of decay channels defined
C...AUTHOR: STEPHEN MRENNA
C...Last change:
C...10-15-95:  force decay chi^0_2 -> chi^0_1 + gamma
C...when CHIGAMMA .NE. 0
C...10 FEB 96:  Calculate this decay for small tan(beta)
 
      SUBROUTINE PYNJDC(KFIN,XLAM,IDLAM,IKNT)
 
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
c      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
c     &SFMIX(16,4)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
C      COMMON/PYINTS/XXM(20)
      COMPLEX*16 CXC
      COMMON/PYINTC/XXC(10),CXC(8)
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/,/PYSSMT/,/PYINTC/
 
C...Local variables.
      COMPLEX*16 ZMIXC(4,4),VMIXC(2,2),UMIXC(2,2),OLPP,ORPP,GLIJ,GRIJ
      COMPLEX*16 QIJ,RIJ,F21K,F12K,CAL,CAR,CBL,CBR,CA,CB
      INTEGER KFIN
      DOUBLE PRECISION XMI,XMJ,XMF,XMSF1,XMSF2,XMW,XMW2,
     &XMZ,XMZ2,AXMJ,AXMI
      DOUBLE PRECISION S12MIN,S12MAX
      DOUBLE PRECISION XMI2,XMI3,XMJ2,XMH,XMH2,XMHP,XMA2,XMB2
      DOUBLE PRECISION PYLAMF,XL
      DOUBLE PRECISION TANW,XW,AEM,C1,AS,EI,T3I
      DOUBLE PRECISION PYX2XH,PYX2XG
      DOUBLE PRECISION XLAM(0:400)
      INTEGER IDLAM(400,3)
      INTEGER LKNT,IX,IH,J,IJ,I,IKNT,FID
      INTEGER ITH(3),KF1,KF2
      INTEGER ITHC
      DOUBLE PRECISION DH(3),EH(3)
      DOUBLE PRECISION SR2
      DOUBLE PRECISION CBETA,SBETA
      DOUBLE PRECISION GAMCON,XMT1,XMT2
      DOUBLE PRECISION PYALEM,PI,PYALPS
      DOUBLE PRECISION RAT1,RAT2
      DOUBLE PRECISION T3T,FCOL
      DOUBLE PRECISION ALFA,BETA,TANB
      DOUBLE PRECISION PYXXGA
      EXTERNAL PYGAUS,PYXXZ6
      DOUBLE PRECISION PYGAUS,PYXXZ6
      DOUBLE PRECISION PREC
      INTEGER KFNCHI(4),KFCCHI(2)
      DATA ITH/25,35,36/
      DATA ITHC/37/
      DATA PREC/1D-2/
      DATA PI/3.141592654D0/
      DATA SR2/1.4142136D0/
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
 
C...IX IS 1 - 4 DEPENDING ON SEQUENCE NUMBER
      IX=1
      IF(KFIN.EQ.KFNCHI(2)) IX=2
      IF(KFIN.EQ.KFNCHI(3)) IX=3
      IF(KFIN.EQ.KFNCHI(4)) IX=4
 
      XMI=SMZ(IX)
      XMI2=XMI**2
      AXMI=ABS(XMI)
      AEM=PYALEM(XMI2)
      AS =PYALPS(XMI2)
      C1=AEM/XW
      XMI3=ABS(XMI**3)
 
      TANB=RMSS(5)
      BETA=ATAN(TANB)
      ALFA=RMSS(18)
      CBETA=COS(BETA)
      SBETA=TANB*CBETA
      CALFA=COS(ALFA)
      SALFA=SIN(ALFA)
 
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
 
C...CHECK ALL 2-BODY DECAYS TO GAUGE AND HIGGS BOSONS
      IF(IX.EQ.1.AND.IMSS(11).EQ.0) GOTO 300
 
C...FORCE CHI0_2 -> CHI0_1 + GAMMA
      IF(IX.EQ.2 .AND. IMSS(10).NE.0 ) THEN
        XMJ=SMZ(1)
        AXMJ=ABS(XMJ)
        LKNT=LKNT+1
        GAMCON=AEM**3/8D0/PI/XMW2/XW
        XMT1=(PMAS(PYCOMP(KSUSY1+6),1)/PMAS(6,1))**2
        XMT2=(PMAS(PYCOMP(KSUSY2+6),1)/PMAS(6,1))**2
        XLAM(LKNT)=PYXXGA(GAMCON,AXMI,AXMJ,XMT1,XMT2)
        IDLAM(LKNT,1)=KSUSY1+22
        IDLAM(LKNT,2)=22
        IDLAM(LKNT,3)=0
        WRITE(MSTU(11),*) 'FORCED N2 -> N1 + GAMMA ',XLAM(LKNT)
        GOTO 340
      ENDIF
 
C...GRAVITINO DECAY MODES
 
      IF(IMSS(11).EQ.1) THEN
        XMP=RMSS(29)
        IDG=39+KSUSY1
        XMGR=PMAS(PYCOMP(IDG),1)
        SINW=SQRT(XW)
        COSW=SQRT(1D0-XW)
        XFAC=(XMI2/(XMP*XMGR))**2*AXMI/48D0/PI
        IF(AXMI.GT.XMGR+PMAS(22,1)) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=22
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*ABS(ZMIXC(IX,1)*COSW+ZMIXC(IX,2)*SINW)**2
        ENDIF
        IF(AXMI.GT.XMGR+XMZ) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=23
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*(ABS(ZMIXC(IX,1)*SINW-ZMIXC(IX,2)*COSW)**2 +
     $  .5D0*ABS(ZMIXC(IX,3)*CBETA-ZMIXC(IX,4)*SBETA)**2)*
     &  (1D0-XMZ2/XMI2)**4
        ENDIF
        IF(AXMI.GT.XMGR+PMAS(25,1)) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=25
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*(ABS(ZMIXC(IX,3)*SALFA-ZMIXC(IX,4)*CALFA)**2)*
     $  .5D0*(1D0-PMAS(25,1)**2/XMI2)**4
        ENDIF
        IF(AXMI.GT.XMGR+PMAS(35,1)) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=35
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*(ABS(ZMIXC(IX,3)*CALFA+ZMIXC(IX,4)*SALFA)**2)*
     $  .5D0*(1D0-PMAS(35,1)**2/XMI2)**4
        ENDIF
        IF(AXMI.GT.XMGR+PMAS(36,1)) THEN
          LKNT=LKNT+1
          IDLAM(LKNT,1)=IDG
          IDLAM(LKNT,2)=36
          IDLAM(LKNT,3)=0
          XLAM(LKNT)=XFAC*(ABS(ZMIXC(IX,3)*SBETA+ZMIXC(IX,4)*CBETA)**2)*
     $  .5D0*(1D0-PMAS(36,1)**2/XMI2)**4
        ENDIF
        IF(IX.EQ.1) GOTO 300
      ENDIF
 
      DO 220 IJ=1,IX-1
        XMJ=SMZ(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
 
C...CHI0_I -> CHI0_J + GAMMA
        IF(AXMI.GE.AXMJ.AND.SBETA/CBETA.LE.2D0) THEN
          RAT1=ABS(ZMIXC(IJ,1))**2+ABS(ZMIXC(IJ,2))**2
          RAT1=RAT1/( 1D-6+ABS(ZMIXC(IX,3))**2+ABS(ZMIXC(IX,4))**2 )
          RAT2=ABS(ZMIXC(IX,1))**2+ABS(ZMIXC(IX,2))**2
          RAT2=RAT2/( 1D-6+ABS(ZMIXC(IJ,3))**2+ABS(ZMIXC(IJ,4))**2 )
          IF((RAT1.GT. 0.90D0 .AND. RAT1.LT. 1.10D0) .OR.
     &    (RAT2.GT. 0.90D0 .AND. RAT2.LT. 1.10D0)) THEN
            LKNT=LKNT+1
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=22
            IDLAM(LKNT,3)=0
            GAMCON=AEM**3/8D0/PI/XMW2/XW
            XMT1=(PMAS(PYCOMP(KSUSY1+6),1)/PMAS(6,1))**2
            XMT2=(PMAS(PYCOMP(KSUSY2+6),1)/PMAS(6,1))**2
            XLAM(LKNT)=PYXXGA(GAMCON,AXMI,AXMJ,XMT1,XMT2)
          ENDIF
        ENDIF
 
C...CHI0_I -> CHI0_J + Z0
        IF(AXMI.GE.AXMJ+XMZ) THEN
          LKNT=LKNT+1
          OLPP=(ZMIXC(IX,3)*DCONJG(ZMIXC(IJ,3))-
     &    ZMIXC(IX,4)*DCONJG(ZMIXC(IJ,4)))/2D0
          ORPP=-DCONJG(OLPP)
          GX2=ABS(OLPP)**2+ABS(ORPP)**2
          GLR=DBLE(OLPP*DCONJG(ORPP))
          XLAM(LKNT)=PYX2XG(C1/XMW2,XMI,XMJ,XMZ,GX2,GLR)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=23
          IDLAM(LKNT,3)=0
        ELSEIF(AXMI.GE.AXMJ) THEN
          XXC(1)=0D0
          XXC(2)=XMJ
          XXC(3)=0D0
          XXC(4)=XMI
          XXC(9)=XMZ
          XXC(10)=PMAS(23,2)
          OLPP=(ZMIXC(IX,3)*DCONJG(ZMIXC(IJ,3))-
     &    ZMIXC(IX,4)*DCONJG(ZMIXC(IJ,4)))/2D0
          ORPP=DCONJG(OLPP)
C...CHARGED LEPTONS
          FID=11
          XXC(5)=PMAS(PYCOMP(KSUSY1+FID),1)
          XXC(6)=PMAS(PYCOMP(KSUSY2+FID),1)
          EI=KCHG(FID,1)/3D0
          T3I=SIGN(1D0,EI+1D-6)/2D0
          GLIJ=(T3I*ZMIXC(IX,2)-TANW*(T3I-EI)*ZMIXC(IX,1))*
     &    DCONJG(T3I*ZMIXC(IJ,2)-TANW*(T3I-EI)*ZMIXC(IJ,1))
          GRIJ=ZMIXC(IX,1)*DCONJG(ZMIXC(IJ,1))*(EI*TANW)**2
          CXC(1)=DCMPLX((T3I-EI*XW)/XW1)*OLPP
          CXC(2)=-GLIJ
          CXC(3)=-DCMPLX((T3I-EI*XW)/XW1)*ORPP
          CXC(4)=DCONJG(GLIJ)
          CXC(5)=-DCMPLX((EI*XW)/XW1)*OLPP
          CXC(6)=GRIJ
          CXC(7)=DCMPLX((EI*XW)/XW1)*ORPP
          CXC(8)=-DCONJG(GRIJ)
          S12MIN=0D0
          S12MAX=(AXMI-AXMJ)**2
          IF( XXC(5).LT.AXMI ) THEN
            XXC(5)=1D6
          ENDIF
          IF(XXC(6).LT.AXMI ) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(5)
          XXC(8)=XXC(6)
 
          IF(AXMI.GE.AXMJ+2D0*PMAS(11,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-3)
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=FID
            IDLAM(LKNT,3)=-FID
            IF(AXMI.GE.AXMJ+2D0*PMAS(13,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFNCHI(IJ)
              IDLAM(LKNT,2)=13
              IDLAM(LKNT,3)=-13
            ENDIF
          ENDIF
  140     CONTINUE
          IF(ABS(SFMIX(15,1)).GT.ABS(SFMIX(15,2))) THEN
            XXC(5)=PMAS(PYCOMP(KSUSY1+15),1)
            XXC(6)=PMAS(PYCOMP(KSUSY2+15),1)
          ELSE
            XXC(6)=PMAS(PYCOMP(KSUSY1+15),1)
            XXC(5)=PMAS(PYCOMP(KSUSY2+15),1)
          ENDIF
          IF( XXC(5).LT.AXMI ) THEN
            XXC(5)=1D6
          ENDIF
          IF(XXC(6).LT.AXMI ) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(5)
          XXC(8)=XXC(6)
 
          IF(AXMI.GE.AXMJ+2D0*PMAS(15,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-3)
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=15
            IDLAM(LKNT,3)=-15
          ENDIF
 
C...NEUTRINOS
  150     CONTINUE
          FID=12
          XXC(5)=PMAS(PYCOMP(KSUSY1+FID),1)
          XXC(6)=PMAS(PYCOMP(KSUSY2+FID),1)
          EI=KCHG(FID,1)/3D0
          T3I=SIGN(1D0,EI+1D-6)/2D0
          GLIJ=(T3I*ZMIXC(IX,2)-TANW*(T3I-EI)*ZMIXC(IX,1))*
     &    DCONJG(T3I*ZMIXC(IJ,2)-TANW*(T3I-EI)*ZMIXC(IJ,1))
          GRIJ=ZMIXC(IX,1)*DCONJG(ZMIXC(IJ,1))*(EI*TANW)**2
          CXC(1)=DCMPLX((T3I-EI*XW)/XW1)*OLPP
          CXC(2)=-GLIJ
          CXC(3)=-DCMPLX((T3I-EI*XW)/XW1)*ORPP
          CXC(4)=DCONJG(GLIJ)
          CXC(5)=-DCMPLX((EI*XW)/XW1)*OLPP
          CXC(6)=GRIJ
          CXC(7)=DCMPLX((EI*XW)/XW1)*ORPP
          CXC(8)=-DCONJG(GRIJ)
          S12MIN=0D0
          S12MAX=(AXMI-AXMJ)**2
          IF( XXC(5).LT.AXMI ) THEN
            XXC(5)=1D6
          ENDIF
          IF( XXC(6).LT.AXMI ) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(5)
          XXC(8)=XXC(6)
 
          LKNT=LKNT+1
          XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-3)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=12
          IDLAM(LKNT,3)=-12
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=14
          IDLAM(LKNT,3)=-14
  160     CONTINUE
 
          IF(PMAS(PYCOMP(KSUSY1+16),1).NE.PMAS(PYCOMP(KSUSY1+12),1))
     &    THEN
            XXC(5)=PMAS(PYCOMP(KSUSY1+16),1)
            IF( XXC(5).LT.AXMI ) THEN
              XXC(5)=1D6
            ENDIF
            XXC(7)=XXC(5)
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-3)
          ELSE
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
          ENDIF
          IDLAM(LKNT,1)=KFNCHI(IJ)
          IDLAM(LKNT,2)=16
          IDLAM(LKNT,3)=-16
C...D-TYPE QUARKS
  170     CONTINUE
          FID=1
          XXC(5)=PMAS(PYCOMP(KSUSY1+FID),1)
          XXC(6)=PMAS(PYCOMP(KSUSY2+FID),1)
          EI=KCHG(FID,1)/3D0
          T3I=SIGN(1D0,EI+1D-6)/2D0
          GLIJ=(T3I*ZMIXC(IX,2)-TANW*(T3I-EI)*ZMIXC(IX,1))*
     &    DCONJG(T3I*ZMIXC(IJ,2)-TANW*(T3I-EI)*ZMIXC(IJ,1))
          GRIJ=ZMIXC(IX,1)*DCONJG(ZMIXC(IJ,1))*(EI*TANW)**2
          CXC(1)=DCMPLX((T3I-EI*XW)/XW1)*OLPP
          CXC(2)=-GLIJ
          CXC(3)=-DCMPLX((T3I-EI*XW)/XW1)*ORPP
          CXC(4)=DCONJG(GLIJ)
          CXC(5)=-DCMPLX((EI*XW)/XW1)*OLPP
          CXC(6)=GRIJ
          CXC(7)=DCMPLX((EI*XW)/XW1)*ORPP
          CXC(8)=-DCONJG(GRIJ)
          S12MIN=0D0
          S12MAX=(AXMI-AXMJ)**2
          IF( XXC(5).LT.AXMI ) THEN
            XXC(5)=1D6
          ENDIF
          IF( XXC(6).LT.AXMI ) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(5)
          XXC(8)=XXC(6)
 
          IF(AXMI.GE.AXMJ+2D0*PMAS(1,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-3)*3D0
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=1
            IDLAM(LKNT,3)=-1
            IF(AXMI.GE.AXMJ+2D0*PMAS(3,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFNCHI(IJ)
              IDLAM(LKNT,2)=3
              IDLAM(LKNT,3)=-3
            ENDIF
          ENDIF
  180     CONTINUE
          IF(ABS(SFMIX(5,1)).GT.ABS(SFMIX(5,2))) THEN
            XXC(5)=PMAS(PYCOMP(KSUSY1+5),1)
            XXC(6)=PMAS(PYCOMP(KSUSY2+5),1)
          ELSE
            XXC(6)=PMAS(PYCOMP(KSUSY1+5),1)
            XXC(5)=PMAS(PYCOMP(KSUSY2+5),1)
          ENDIF
          IF( XXC(5).LT.AXMI .AND. XXC(6).LT.AXMI ) GOTO 190
          IF(XXC(5).LT.AXMI) THEN
            XXC(5)=1D6
          ELSEIF(XXC(6).LT.AXMI) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(5)
          XXC(8)=XXC(6)
          IF(AXMI.GE.AXMJ+2D0*PMAS(5,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-3)*3D0
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=5
            IDLAM(LKNT,3)=-5
          ENDIF
 
C...U-TYPE QUARKS
  190     CONTINUE
          FID=2
          XXC(5)=PMAS(PYCOMP(KSUSY1+FID),1)
          XXC(6)=PMAS(PYCOMP(KSUSY2+FID),1)
          EI=KCHG(FID,1)/3D0
          T3I=SIGN(1D0,EI+1D-6)/2D0
          GLIJ=(T3I*ZMIXC(IX,2)-TANW*(T3I-EI)*ZMIXC(IX,1))*
     &    DCONJG(T3I*ZMIXC(IJ,2)-TANW*(T3I-EI)*ZMIXC(IJ,1))
          GRIJ=ZMIXC(IX,1)*DCONJG(ZMIXC(IJ,1))*(EI*TANW)**2
          CXC(1)=DCMPLX((T3I-EI*XW)/XW1)*OLPP
          CXC(2)=-GLIJ
          CXC(3)=-DCMPLX((T3I-EI*XW)/XW1)*ORPP
          CXC(4)=DCONJG(GLIJ)
          CXC(5)=-DCMPLX((EI*XW)/XW1)*OLPP
          CXC(6)=GRIJ
          CXC(7)=DCMPLX((EI*XW)/XW1)*ORPP
          CXC(8)=-DCONJG(GRIJ)
 
          IF( XXC(5).LT.AXMI .AND. XXC(6).LT.AXMI ) GOTO 200
          IF(XXC(5).LT.AXMI) THEN
            XXC(5)=1D6
          ELSEIF(XXC(6).LT.AXMI) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(5)
          XXC(8)=XXC(6)
 
          IF(AXMI.GE.AXMJ+2D0*PMAS(2,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-3)*3D0
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=2
            IDLAM(LKNT,3)=-2
            IF(AXMI.GE.AXMJ+2D0*PMAS(4,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFNCHI(IJ)
              IDLAM(LKNT,2)=4
              IDLAM(LKNT,3)=-4
            ENDIF
          ENDIF
  200     CONTINUE
        ENDIF
 
C...CHI0_I -> CHI0_J + H0_K
        EH(1)=SIN(ALFA)
        EH(2)=COS(ALFA)
        EH(3)=-SIN(BETA)
        DH(1)=COS(ALFA)
        DH(2)=-SIN(ALFA)
        DH(3)=COS(BETA)
        QIJ=ZMIXC(IX,3)*DCONJG(ZMIXC(IJ,2))+
     &  DCONJG(ZMIXC(IJ,3))*ZMIXC(IX,2)-
     &  TANW*(ZMIXC(IX,3)*DCONJG(ZMIXC(IJ,1))+
     &  DCONJG(ZMIXC(IJ,3))*ZMIXC(IX,1))
        RIJ=DCONJG(ZMIXC(IX,4))*ZMIXC(IJ,2)+
     &  ZMIXC(IJ,4)*DCONJG(ZMIXC(IX,2))-
     &  TANW*(DCONJG(ZMIXC(IX,4))*ZMIXC(IJ,1)+
     &  ZMIXC(IJ,4)*DCONJG(ZMIXC(IX,1)))
        DO 210 IH=1,3
          XMH=PMAS(ITH(IH),1)
          XMH2=XMH**2
          IF(AXMI.GE.AXMJ+XMH) THEN
            LKNT=LKNT+1
            XL=PYLAMF(XMI2,XMJ2,XMH2)
            F21K=0.5D0*(QIJ*EH(IH)+RIJ*DH(IH))
            F12K=F21K
C...SIGN OF MASSES I,J
            XMK=XMJ
            IF(IH.EQ.3) XMK=-XMK
            GX2=ABS(F21K)**2+ABS(F12K)**2
            GLR=DBLE(F21K*DCONJG(F12K))
            XLAM(LKNT)=PYX2XH(C1,XMI,XMK,XMH,GX2,GLR)
            IDLAM(LKNT,1)=KFNCHI(IJ)
            IDLAM(LKNT,2)=ITH(IH)
            IDLAM(LKNT,3)=0
          ENDIF
  210   CONTINUE
  220 CONTINUE
 
C...CHI0_I -> CHI+_J + W-
      DO 260 IJ=1,2
        XMJ=SMW(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
        IF(AXMI.GE.AXMJ+XMW) THEN
          LKNT=LKNT+1
          CXC(1)=(DCONJG(ZMIXC(IX,2))*VMIXC(IJ,1)-
     &    DCONJG(ZMIXC(IX,4))*VMIXC(IJ,2)/SR2)
          CXC(3)=(ZMIXC(IX,2)*DCONJG(UMIXC(IJ,1))+
     &    ZMIXC(IX,3)*DCONJG(UMIXC(IJ,2))/SR2)
          GX2=ABS(CXC(1))**2+ABS(CXC(3))**2
          GLR=DBLE(CXC(1)*DCONJG(CXC(3)))
          XLAM(LKNT)=PYX2XG(C1/XMW2,XMI,XMJ,XMW,GX2,GLR)
          IDLAM(LKNT,1)=KFCCHI(IJ)
          IDLAM(LKNT,2)=-24
          IDLAM(LKNT,3)=0
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=-KFCCHI(IJ)
          IDLAM(LKNT,2)=24
          IDLAM(LKNT,3)=0
        ELSEIF(AXMI.GE.AXMJ) THEN
          S12MIN=0D0
          S12MAX=(AXMI-AXMJ)**2
          RT2I = 1D0/SQRT(2D0)
          CXC(1)=(DCONJG(ZMIXC(IX,2))*VMIXC(IJ,1)-
     &    DCONJG(ZMIXC(IX,4))*VMIXC(IJ,2)*RT2I)*RT2I
          CXC(3)=(ZMIXC(IX,2)*DCONJG(UMIXC(IJ,1))+
     &    ZMIXC(IX,3)*DCONJG(UMIXC(IJ,2))*RT2I)*RT2I
          CXC(5)=DCMPLX(0D0,0D0)
          CXC(7)=DCMPLX(0D0,0D0)
          IA=11
          JA=12
          EI=KCHG(IA,1)/3D0
          T3I=SIGN(1D0,EI+1D-6)/2D0
          EJ=KCHG(JA,1)/3D0
          T3J=SIGN(1D0,EJ+1D-6)/2D0
          CXC(2)=VMIXC(IJ,1)*DCONJG(ZMIXC(IX,1)*(EJ-T3J)*
     &    TANW+ZMIXC(IX,2)*T3J)*RT2I
          CXC(4)=-DCONJG(UMIXC(IJ,1))*(
     &    ZMIXC(IX,1)*(EI-T3I)*TANW+ZMIXC(IX,2)*T3I)*RT2I
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
          IF( XXC(5).LT.AXMI .AND. XXC(6).LT.AXMI ) GOTO 230
          IF(XXC(5).LT.AXMI) THEN
            XXC(5)=1D6
          ELSEIF(XXC(6).LT.AXMI) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(6)
          XXC(8)=XXC(5)
          IF(AXMI.GE.AXMJ+PMAS(11,1)+PMAS(12,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
            IDLAM(LKNT,1)=KFCCHI(IJ)
            IDLAM(LKNT,2)=11
            IDLAM(LKNT,3)=-12
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
            IF(AXMI.GE.AXMJ+PMAS(13,1)+PMAS(14,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFCCHI(IJ)
              IDLAM(LKNT,2)=13
              IDLAM(LKNT,3)=-14
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
              IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
              IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
            ENDIF
          ENDIF
  230     CONTINUE
          IF(ABS(SFMIX(15,1)).GT.ABS(SFMIX(15,2))) THEN
            XXC(5)=PMAS(PYCOMP(KSUSY1+15),1)
            XXC(6)=PMAS(PYCOMP(KSUSY1+16),1)
          ELSE
            XXC(5)=PMAS(PYCOMP(KSUSY2+15),1)
            XXC(6)=PMAS(PYCOMP(KSUSY1+16),1)
          ENDIF
          IF(XXC(5).LT.AXMI) THEN
            XXC(5)=1D6
          ENDIF
          IF(XXC(6).LT.AXMI) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(6)
          XXC(8)=XXC(5)
          IF(AXMI.GE.AXMJ+PMAS(15,1)+PMAS(16,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KFCCHI(IJ)
            IDLAM(LKNT,2)=15
            IDLAM(LKNT,3)=-16
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
          ENDIF
 
C...NOW, DO THE QUARKS
  240     CONTINUE
          IA=1
          JA=2
          EI=KCHG(IA,1)/3D0
          T3I=SIGN(1D0,EI+1D-6)/2D0
          EJ=KCHG(JA,1)/3D0
          T3J=SIGN(1D0,EJ+1D-6)/2D0
          CXC(2)=VMIXC(IJ,1)*DCONJG(ZMIXC(IX,1)*(EJ-T3J)*
     &    TANW+ZMIXC(IX,2)*T3J)
          CXC(4)=-DCONJG(UMIXC(IJ,1))*(
     &    ZMIXC(IX,1)*(EI-T3I)*TANW+ZMIXC(IX,2)*T3I)
          XXC(5)=PMAS(PYCOMP(KSUSY1+IA),1)
          XXC(6)=PMAS(PYCOMP(KSUSY1+JA),1)
          IF(XXC(5).LT.AXMI) THEN
            XXC(5)=1D6
          ENDIF
          IF(XXC(6).LT.AXMI) THEN
            XXC(6)=1D6
          ENDIF
          XXC(7)=XXC(6)
          XXC(8)=XXC(5)
          IF(AXMI.GE.AXMJ+PMAS(2,1)+PMAS(1,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=3D0*C1**2/XMI3/(16D0*PI)*
     &      PYGAUS(PYXXZ6,S12MIN,S12MAX,PREC)
            IDLAM(LKNT,1)=KFCCHI(IJ)
            IDLAM(LKNT,2)=1
            IDLAM(LKNT,3)=-2
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
            IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
            IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
            IF(AXMI.GE.AXMJ+PMAS(3,1)+PMAS(4,1)) THEN
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=KFCCHI(IJ)
              IDLAM(LKNT,2)=3
              IDLAM(LKNT,3)=-4
              LKNT=LKNT+1
              XLAM(LKNT)=XLAM(LKNT-1)
              IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
              IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
              IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
            ENDIF
          ENDIF
  250     CONTINUE
        ENDIF
  260 CONTINUE
  270 CONTINUE
 
C...CHI0_I -> CHI+_I + H-
      DO 280 IJ=1,2
        XMJ=SMW(IJ)
        AXMJ=ABS(XMJ)
        XMJ2=XMJ**2
        XMHP=PMAS(ITHC,1)
        IF(AXMI.GE.AXMJ+XMHP) THEN
          LKNT=LKNT+1
          OLPP=CBETA*(ZMIXC(IX,4)*DCONJG(VMIXC(IJ,1))+(ZMIXC(IX,2)+
     &    ZMIXC(IX,1)*TANW)*DCONJG(VMIXC(IJ,2))/SR2)
          ORPP=SBETA*(DCONJG(ZMIXC(IX,3))*UMIXC(IJ,1)-
     &    (DCONJG(ZMIXC(IX,2))+DCONJG(ZMIXC(IX,1))*TANW)*
     &    UMIXC(IJ,2)/SR2)
          GX2=ABS(OLPP)**2+ABS(ORPP)**2
          GLR=DBLE(OLPP*DCONJG(ORPP))
          XLAM(LKNT)=PYX2XH(C1,XMI,XMJ,XMHP,GX2,GLR)
          IDLAM(LKNT,1)=KFCCHI(IJ)
          IDLAM(LKNT,2)=-ITHC
          IDLAM(LKNT,3)=0
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
          IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
          IDLAM(LKNT,3)=-IDLAM(LKNT-1,3)
        ELSE
 
        ENDIF
  280 CONTINUE
 
C...2-BODY DECAYS TO FERMION SFERMION
      DO 290 J=1,16
        IF(J.GE.7.AND.J.LE.10) GOTO 290
        KF1=KSUSY1+J
        KF2=KSUSY2+J
        XMSF1=PMAS(PYCOMP(KF1),1)
        XMSF2=PMAS(PYCOMP(KF2),1)
        XMF=PMAS(J,1)
        IF(J.LE.6) THEN
          FCOL=3D0
        ELSE
          FCOL=1D0
        ENDIF
 
        EI=KCHG(J,1)/3D0
        T3T=SIGN(1D0,EI)
        IF(J.EQ.12.OR.J.EQ.14.OR.J.EQ.16) T3T=1D0
        IF(MOD(J,2).EQ.0) THEN
          CBL=T3T*ZMIXC(IX,2)+TANW*ZMIXC(IX,1)*(2D0*EI-T3T)
          CAL=XMF*ZMIXC(IX,4)/XMW/SBETA
          CAR=-2D0*EI*TANW*ZMIXC(IX,1)
          CBR=CAL
        ELSE
          CBL=T3T*ZMIXC(IX,2)+TANW*ZMIXC(IX,1)*(2D0*EI-T3T)
          CAL=XMF*ZMIXC(IX,3)/XMW/CBETA
          CAR=-2D0*EI*TANW*ZMIXC(IX,1)
          CBR=CAL
        ENDIF
 
C...D~ D_L
        IF(AXMI.GE.XMF+XMSF1) THEN
          LKNT=LKNT+1
          XMA2=XMSF1**2
          XMB2=XMF**2
          XL=PYLAMF(XMI2,XMA2,XMB2)
          CA=CAL*SFMIX(J,1)+CAR*SFMIX(J,2)
          CB=CBL*SFMIX(J,1)+CBR*SFMIX(J,2)
          XLAM(LKNT)=0.5D0*FCOL*C1/8D0/XMI3*SQRT(XL)*( (XMI2+XMB2-XMA2)*
     &    (ABS(CA)**2+ABS(CB)**2)+4D0*DBLE(CA*DCONJG(CB))*XMF*XMI)
          IDLAM(LKNT,1)=KF1
          IDLAM(LKNT,2)=-J
          IDLAM(LKNT,3)=0
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
          IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
          IDLAM(LKNT,3)=0
        ENDIF
 
C...D~ D_R
        IF(AXMI.GE.XMF+XMSF2) THEN
          LKNT=LKNT+1
          XMA2=XMSF2**2
          XMB2=XMF**2
          CA=CAL*SFMIX(J,3)+CAR*SFMIX(J,4)
          CB=CBL*SFMIX(J,3)+CBR*SFMIX(J,4)
          XL=PYLAMF(XMI2,XMA2,XMB2)
          XLAM(LKNT)=0.5D0*FCOL*C1/8D0/XMI3*SQRT(XL)*( (XMI2+XMB2-XMA2)*
     &    (ABS(CA)**2+ABS(CB)**2)+4D0*DBLE(CA*DCONJG(CB))*XMF*XMI)
          IDLAM(LKNT,1)=KF2
          IDLAM(LKNT,2)=-J
          IDLAM(LKNT,3)=0
          LKNT=LKNT+1
          XLAM(LKNT)=XLAM(LKNT-1)
          IDLAM(LKNT,1)=-IDLAM(LKNT-1,1)
          IDLAM(LKNT,2)=-IDLAM(LKNT-1,2)
          IDLAM(LKNT,3)=0
        ENDIF
  290 CONTINUE
  300 CONTINUE
C...3-BODY DECAY TO Q Q~ GLUINO
      XMJ=PMAS(PYCOMP(KSUSY1+21),1)
      IF(AXMI.GE.XMJ) THEN
        RT2I = 1D0/SQRT(2D0)
        OLPP=DCMPLX(COS(RMSS(32)),SIN(RMSS(32)))*RT2I
        ORPP=DCONJG(OLPP)
        AXMJ=ABS(XMJ)
        XXC(1)=0D0
        XXC(2)=XMJ
        XXC(3)=0D0
        XXC(4)=XMI
        FID=1
        XXC(5)=PMAS(PYCOMP(KSUSY1+FID),1)
        XXC(6)=PMAS(PYCOMP(KSUSY2+FID),1)
        XXC(7)=XXC(5)
        XXC(8)=XXC(6)
        XXC(9)=1D6
        XXC(10)=0D0
        EI=KCHG(FID,1)/3D0
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
CMRENNA.This statement must be here to define S12MAX
        IF( XXC(5).LT.AXMI .OR. XXC(6).LT.AXMI ) GOTO 310
C...ALL QUARKS BUT T
        IF(AXMI.GE.AXMJ+2D0*PMAS(1,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=4D0*C1*AS/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-3)
          IDLAM(LKNT,1)=KSUSY1+21
          IDLAM(LKNT,2)=1
          IDLAM(LKNT,3)=-1
          IF(AXMI.GE.AXMJ+2D0*PMAS(3,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KSUSY1+21
            IDLAM(LKNT,2)=3
            IDLAM(LKNT,3)=-3
          ENDIF
        ENDIF
  310   CONTINUE
        IF(ABS(SFMIX(5,1)).GT.ABS(SFMIX(5,2))) THEN
          XXC(5)=PMAS(PYCOMP(KSUSY1+5),1)
          XXC(6)=PMAS(PYCOMP(KSUSY2+5),1)
        ELSE
          XXC(6)=PMAS(PYCOMP(KSUSY1+5),1)
          XXC(5)=PMAS(PYCOMP(KSUSY2+5),1)
        ENDIF
        IF( XXC(5).LT.AXMI .OR. XXC(6).LT.AXMI ) GOTO 320
        XXC(7)=XXC(5)
        XXC(8)=XXC(6)
        IF(AXMI.GE.AXMJ+2D0*PMAS(5,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=0.5D0*C1*AS/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-3)
          IDLAM(LKNT,1)=KSUSY1+21
          IDLAM(LKNT,2)=5
          IDLAM(LKNT,3)=-5
        ENDIF
C...U-TYPE QUARKS
  320   CONTINUE
        FID=2
        XXC(5)=PMAS(PYCOMP(KSUSY1+FID),1)
        XXC(6)=PMAS(PYCOMP(KSUSY2+FID),1)
        IF( XXC(5).LT.AXMI .OR. XXC(6).LT.AXMI ) GOTO 330
        XXC(7)=XXC(5)
        XXC(8)=XXC(6)
        EI=KCHG(FID,1)/3D0
        T3I=SIGN(1D0,EI+1D-6)/2D0
        GLIJ=(T3I*ZMIXC(IX,2)-TANW*(T3I-EI)*ZMIXC(IX,1))*OLPP
        GRIJ=ZMIXC(IX,1)*(EI*TANW)*ORPP
        CXC(2)=-GLIJ
        CXC(4)=DCONJG(GLIJ)
        CXC(6)=GRIJ
        CXC(8)=-DCONJG(GRIJ)
        IF(AXMI.GE.AXMJ+2D0*PMAS(2,1)) THEN
          LKNT=LKNT+1
          XLAM(LKNT)=0.5D0*C1*AS/XMI3/(16D0*PI)*
     &    PYGAUS(PYXXZ6,S12MIN,S12MAX,1D-3)
          IDLAM(LKNT,1)=KSUSY1+21
          IDLAM(LKNT,2)=2
          IDLAM(LKNT,3)=-2
          IF(AXMI.GE.AXMJ+2D0*PMAS(4,1)) THEN
            LKNT=LKNT+1
            XLAM(LKNT)=XLAM(LKNT-1)
            IDLAM(LKNT,1)=KSUSY1+21
            IDLAM(LKNT,2)=4
            IDLAM(LKNT,3)=-4
          ENDIF
        ENDIF
  330   CONTINUE
      ENDIF
 
C...R-violating decay modes (SKANDS).
      CALL PYRVNE(KFIN,XLAM,IDLAM,LKNT)
 
  340 IKNT=LKNT
      XLAM(0)=0D0
      DO 350 I=1,IKNT
        IF(XLAM(I).LT.0D0) XLAM(I)=0D0
        XLAM(0)=XLAM(0)+XLAM(I)
  350 CONTINUE
      IF(XLAM(0).EQ.0D0) XLAM(0)=1D-6
 
      RETURN
      END
