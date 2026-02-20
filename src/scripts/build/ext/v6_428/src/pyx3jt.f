 
C*********************************************************************
 
C...PYX3JT
C...Selects the kinematical variables of three-jet events.
 
      SUBROUTINE PYX3JT(NJET,CUT,KFL,ECM,X1,X2)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
C...Local array.
      DIMENSION ZHUP(5,12)
 
C...Coefficients of Zhu second order parametrization.
      DATA ((ZHUP(IC1,IC2),IC2=1,12),IC1=1,5)/
     &18.29D0,  89.56D0,  4.541D0,  -52.09D0, -109.8D0,  24.90D0,
     &11.63D0,  3.683D0,  17.50D0,0.002440D0, -1.362D0,-0.3537D0,
     &11.42D0,  6.299D0, -22.55D0,  -8.915D0,  59.25D0, -5.855D0,
     &-32.85D0, -1.054D0, -16.90D0,0.006489D0,-0.8156D0,0.01095D0,
     &7.847D0, -3.964D0, -35.83D0,   1.178D0,  29.39D0, 0.2806D0,
     &47.82D0, -12.36D0, -56.72D0, 0.04054D0,-0.4365D0, 0.6062D0,
     &5.441D0, -56.89D0, -50.27D0,   15.13D0,  114.3D0, -18.19D0,
     &97.05D0, -1.890D0, -139.9D0, 0.08153D0,-0.4984D0, 0.9439D0,
     &-17.65D0,  51.44D0, -58.32D0,   70.95D0, -255.7D0, -78.99D0,
     &476.9D0,  29.65D0, -239.3D0,  0.4745D0, -1.174D0,  6.081D0/
 
C...Dilogarithm of x for x<0.5 (x>0.5 obtained by analytic trick).
      DILOG(X)=X+X**2/4D0+X**3/9D0+X**4/16D0+X**5/25D0+X**6/36D0+
     &X**7/49D0
 
C...Event type. Mass effect factors and other common constants.
      MSTJ(120)=2
      MSTJ(121)=0
      PMQ=PYMASS(KFL)
      QME=(2D0*PMQ/ECM)**2
      IF(MSTJ(109).NE.1) THEN
        CUTL=LOG(CUT)
        CUTD=LOG(1D0/CUT-2D0)
        IF(MSTJ(109).EQ.0) THEN
          CF=4D0/3D0
          CN=3D0
          TR=2D0
          WTMX=MIN(20D0,37D0-6D0*CUTD)
          IF(MSTJ(110).EQ.2) WTMX=2D0*(7.5D0+80D0*CUT)
        ELSE
          CF=1D0
          CN=0D0
          TR=12D0
          WTMX=0D0
        ENDIF
 
C...Alpha_strong and effects of optimized Q^2 scale. Maximum weight.
        ALS2PI=PARU(118)/PARU(2)
        WTOPT=0D0
        IF(MSTJ(111).EQ.1) WTOPT=(33D0-2D0*MSTU(112))/6D0*
     &  LOG(PARJ(169))*ALS2PI
        WTMAX=MAX(0D0,1D0+WTOPT+ALS2PI*WTMX)
 
C...Choose three-jet events in allowed region.
  100   NJET=3
  110   Y13L=CUTL+CUTD*PYR(0)
        Y23L=CUTL+CUTD*PYR(0)
        Y13=EXP(Y13L)
        Y23=EXP(Y23L)
        Y12=1D0-Y13-Y23
        IF(Y12.LE.CUT) GOTO 110
        IF(Y13**2+Y23**2+2D0*Y12.LE.2D0*PYR(0)) GOTO 110
 
C...Second order corrections.
        IF(MSTJ(101).EQ.2.AND.MSTJ(110).LE.1) THEN
          Y12L=LOG(Y12)
          Y13M=LOG(1D0-Y13)
          Y23M=LOG(1D0-Y23)
          Y12M=LOG(1D0-Y12)
          IF(Y13.LE.0.5D0) Y13I=DILOG(Y13)
          IF(Y13.GE.0.5D0) Y13I=1.644934D0-Y13L*Y13M-DILOG(1D0-Y13)
          IF(Y23.LE.0.5D0) Y23I=DILOG(Y23)
          IF(Y23.GE.0.5D0) Y23I=1.644934D0-Y23L*Y23M-DILOG(1D0-Y23)
          IF(Y12.LE.0.5D0) Y12I=DILOG(Y12)
          IF(Y12.GE.0.5D0) Y12I=1.644934D0-Y12L*Y12M-DILOG(1D0-Y12)
          WT1=(Y13**2+Y23**2+2D0*Y12)/(Y13*Y23)
          WT2=CF*(-2D0*(CUTL-Y12L)**2-3D0*CUTL-1D0+3.289868D0+
     &    2D0*(2D0*CUTL-Y12L)*CUT/Y12)+
     &    CN*((CUTL-Y12L)**2-(CUTL-Y13L)**2-(CUTL-Y23L)**2-
     &    11D0*CUTL/6D0+67D0/18D0+1.644934D0-(2D0*CUTL-Y12L)*CUT/Y12+
     &    (2D0*CUTL-Y13L)*CUT/Y13+(2D0*CUTL-Y23L)*CUT/Y23)+
     &    TR*(2D0*CUTL/3D0-10D0/9D0)+
     &    CF*(Y12/(Y12+Y13)+Y12/(Y12+Y23)+(Y12+Y23)/Y13+(Y12+Y13)/Y23+
     &    Y13L*(4D0*Y12**2+2D0*Y12*Y13+4D0*Y12*Y23+Y13*Y23)/
     &    (Y12+Y23)**2+Y23L*(4D0*Y12**2+2D0*Y12*Y23+4D0*Y12*Y13+
     &    Y13*Y23)/(Y12+Y13)**2)/WT1+
     &    CN*(Y13L*Y13/(Y12+Y23)+Y23L*Y23/(Y12+Y13))/WT1+(CN-2D0*CF)*
     &    ((Y12**2+(Y12+Y13)**2)*(Y12L*Y23L-Y12L*Y12M-Y23L*
     &    Y23M+1.644934D0-Y12I-Y23I)/(Y13*Y23)+(Y12**2+(Y12+Y23)**2)*
     &    (Y12L*Y13L-Y12L*Y12M-Y13L*Y13M+1.644934D0-Y12I-Y13I)/
     &    (Y13*Y23)+(Y13**2+Y23**2)/(Y13*Y23*(Y13+Y23))-
     &    2D0*Y12L*Y12**2/(Y13+Y23)**2-4D0*Y12L*Y12/(Y13+Y23))/WT1-
     &    CN*(Y13L*Y23L-Y13L*Y13M-Y23L*Y23M+1.644934D0-Y13I-Y23I)
          IF(1D0+WTOPT+ALS2PI*WT2.LE.0D0) MSTJ(121)=1
          IF(1D0+WTOPT+ALS2PI*WT2.LE.WTMAX*PYR(0)) GOTO 110
          PARJ(156)=(WTOPT+ALS2PI*WT2)/(1D0+WTOPT+ALS2PI*WT2)
 
        ELSEIF(MSTJ(101).EQ.2.AND.MSTJ(110).EQ.2) THEN
C...Second order corrections; Zhu parametrization of ERT.
          ZX=(Y23-Y13)**2
          ZY=1D0-Y12
          IZA=0
          DO 120 IY=1,5
            IF(ABS(CUT-0.01D0*IY).LT.0.0001D0) IZA=IY
  120     CONTINUE
          IF(IZA.NE.0) THEN
            IZ=IZA
            WT2=ZHUP(IZ,1)+ZHUP(IZ,2)*ZX+ZHUP(IZ,3)*ZX**2+(ZHUP(IZ,4)+
     &      ZHUP(IZ,5)*ZX)*ZY+(ZHUP(IZ,6)+ZHUP(IZ,7)*ZX)*ZY**2+
     &      (ZHUP(IZ,8)+ZHUP(IZ,9)*ZX)*ZY**3+ZHUP(IZ,10)/(ZX-ZY**2)+
     &      ZHUP(IZ,11)/(1D0-ZY)+ZHUP(IZ,12)/ZY
          ELSE
            IZ=100D0*CUT
            WTL=ZHUP(IZ,1)+ZHUP(IZ,2)*ZX+ZHUP(IZ,3)*ZX**2+(ZHUP(IZ,4)+
     &      ZHUP(IZ,5)*ZX)*ZY+(ZHUP(IZ,6)+ZHUP(IZ,7)*ZX)*ZY**2+
     &      (ZHUP(IZ,8)+ZHUP(IZ,9)*ZX)*ZY**3+ZHUP(IZ,10)/(ZX-ZY**2)+
     &      ZHUP(IZ,11)/(1D0-ZY)+ZHUP(IZ,12)/ZY
            IZ=IZ+1
            WTU=ZHUP(IZ,1)+ZHUP(IZ,2)*ZX+ZHUP(IZ,3)*ZX**2+(ZHUP(IZ,4)+
     &      ZHUP(IZ,5)*ZX)*ZY+(ZHUP(IZ,6)+ZHUP(IZ,7)*ZX)*ZY**2+
     &      (ZHUP(IZ,8)+ZHUP(IZ,9)*ZX)*ZY**3+ZHUP(IZ,10)/(ZX-ZY**2)+
     &      ZHUP(IZ,11)/(1D0-ZY)+ZHUP(IZ,12)/ZY
            WT2=WTL+(WTU-WTL)*(100D0*CUT+1D0-IZ)
          ENDIF
          IF(1D0+WTOPT+2D0*ALS2PI*WT2.LE.0D0) MSTJ(121)=1
          IF(1D0+WTOPT+2D0*ALS2PI*WT2.LE.WTMAX*PYR(0)) GOTO 110
          PARJ(156)=(WTOPT+2D0*ALS2PI*WT2)/(1D0+WTOPT+2D0*ALS2PI*WT2)
        ENDIF
 
C...Impose mass cuts (gives two jets). For fixed jet number new try.
        X1=1D0-Y23
        X2=1D0-Y13
        X3=1D0-Y12
        IF(4D0*Y23*Y13*Y12/X3**2.LE.QME) NJET=2
        IF(MOD(MSTJ(103),4).GE.2.AND.IABS(MSTJ(101)).LE.1.AND.QME*X3+
     &  0.5D0*QME**2+(0.5D0*QME+0.25D0*QME**2)*((1D0-X2)/(1D0-X1)+
     &  (1D0-X1)/(1D0-X2)).GT.(X1**2+X2**2)*PYR(0)) NJET=2
        IF(MSTJ(101).EQ.-1.AND.NJET.EQ.2) GOTO 100
 
C...Scalar gluon model (first order only, no mass effects).
      ELSE
  130   NJET=3
  140   X3=SQRT(4D0*CUT**2+PYR(0)*((1D0-CUT)**2-4D0*CUT**2))
        IF(LOG((X3-CUT)/CUT).LE.PYR(0)*LOG((1D0-2D0*CUT)/CUT)) GOTO 140
        YD=SIGN(2D0*CUT*((X3-CUT)/CUT)**PYR(0)-X3,PYR(0)-0.5D0)
        X1=1D0-0.5D0*(X3+YD)
        X2=1D0-0.5D0*(X3-YD)
        IF(4D0*(1D0-X1)*(1D0-X2)*(1D0-X3)/X3**2.LE.QME) NJET=2
        IF(MSTJ(102).GE.2) THEN
          IF(X3**2-2D0*(1D0+X3)*(1D0-X1)*(1D0-X2)*PARJ(171).LT.
     &    X3**2*PYR(0)) NJET=2
        ENDIF
        IF(MSTJ(101).EQ.-1.AND.NJET.EQ.2) GOTO 130
      ENDIF
 
      RETURN
      END
