 
C*********************************************************************
 
C...PYPDPI
C...Gives pi+ parton distribution according to two different
C...parametrizations.
 
      SUBROUTINE PYPDPI(X,Q2,XPPI)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYDAT1/,/PYPARS/,/PYINT1/
C...Local arrays.
      DIMENSION XPPI(-6:6),COW(3,5,4,2),XQ(9),TS(6)
 
C...The following data lines are coefficients needed in the
C...Owens pion parton distribution parametrizations, see below.
C...Expansion coefficients for up and down valence quark distributions.
      DATA ((COW(IP,IS,1,1),IS=1,5),IP=1,3)/
     &4.0000D-01,  7.0000D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     &-6.2120D-02,  6.4780D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     &-7.1090D-03,  1.3350D-02,  0.0000D+00,  0.0000D+00,  0.0000D+00/
      DATA ((COW(IP,IS,1,2),IS=1,5),IP=1,3)/
     &4.0000D-01,  6.2800D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     &-5.9090D-02,  6.4360D-01,  0.0000D+00,  0.0000D+00,  0.0000D+00,
     &-6.5240D-03,  1.4510D-02,  0.0000D+00,  0.0000D+00,  0.0000D+00/
C...Expansion coefficients for gluon distribution.
      DATA ((COW(IP,IS,2,1),IS=1,5),IP=1,3)/
     &8.8800D-01,  0.0000D+00,  3.1100D+00,  6.0000D+00,  0.0000D+00,
     &-1.8020D+00, -1.5760D+00, -1.3170D-01,  2.8010D+00, -1.7280D+01,
     &1.8120D+00,  1.2000D+00,  5.0680D-01, -1.2160D+01,  2.0490D+01/
      DATA ((COW(IP,IS,2,2),IS=1,5),IP=1,3)/
     &7.9400D-01,  0.0000D+00,  2.8900D+00,  6.0000D+00,  0.0000D+00,
     &-9.1440D-01, -1.2370D+00,  5.9660D-01, -3.6710D+00, -8.1910D+00,
     &5.9660D-01,  6.5820D-01, -2.5500D-01, -2.3040D+00,  7.7580D+00/
C...Expansion coefficients for (up+down+strange) quark sea distribution.
      DATA ((COW(IP,IS,3,1),IS=1,5),IP=1,3)/
     &9.0000D-01,  0.0000D+00,  5.0000D+00,  0.0000D+00,  0.0000D+00,
     &-2.4280D-01, -2.1200D-01,  8.6730D-01,  1.2660D+00,  2.3820D+00,
     &1.3860D-01,  3.6710D-03,  4.7470D-02, -2.2150D+00,  3.4820D-01/
      DATA ((COW(IP,IS,3,2),IS=1,5),IP=1,3)/
     &9.0000D-01,  0.0000D+00,  5.0000D+00,  0.0000D+00,  0.0000D+00,
     &-1.4170D-01, -1.6970D-01, -2.4740D+00, -2.5340D+00,  5.6210D-01,
     &-1.7400D-01, -9.6230D-02,  1.5750D+00,  1.3780D+00, -2.7010D-01/
C...Expansion coefficients for charm quark sea distribution.
      DATA ((COW(IP,IS,4,1),IS=1,5),IP=1,3)/
     &0.0000D+00, -2.2120D-02,  2.8940D+00,  0.0000D+00,  0.0000D+00,
     &7.9280D-02, -3.7850D-01,  9.4330D+00,  5.2480D+00,  8.3880D+00,
     &-6.1340D-02, -1.0880D-01, -1.0852D+01, -7.1870D+00, -1.1610D+01/
      DATA ((COW(IP,IS,4,2),IS=1,5),IP=1,3)/
     &0.0000D+00, -8.8200D-02,  1.9240D+00,  0.0000D+00,  0.0000D+00,
     &6.2290D-02, -2.8920D-01,  2.4240D-01, -4.4630D+00, -8.3670D-01,
     &-4.0990D-02, -1.0820D-01,  2.0360D+00,  5.2090D+00, -4.8400D-02/
 
C...Euler's beta function, requires ordinary Gamma function
      EULBET(X,Y)=PYGAMM(X)*PYGAMM(Y)/PYGAMM(X+Y)
 
C...Reset output array.
      DO 100 KFL=-6,6
        XPPI(KFL)=0D0
  100 CONTINUE
 
      IF(MSTP(53).LE.2) THEN
C...Pion parton distributions from Owens.
C...Allowed variable range: 4 GeV^2 < Q^2 < approx 2000 GeV^2.
 
C...Determine set, Lambda and s expansion variable.
        NSET=MSTP(53)
        IF(NSET.EQ.1) ALAM=0.2D0
        IF(NSET.EQ.2) ALAM=0.4D0
        VINT(231)=4D0
        IF(MSTP(57).LE.0) THEN
          SD=0D0
        ELSE
          Q2IN=MIN(2D3,MAX(4D0,Q2))
          SD=LOG(LOG(Q2IN/ALAM**2)/LOG(4D0/ALAM**2))
        ENDIF
 
C...Calculate parton distributions.
        DO 120 KFL=1,4
          DO 110 IS=1,5
            TS(IS)=COW(1,IS,KFL,NSET)+COW(2,IS,KFL,NSET)*SD+
     &      COW(3,IS,KFL,NSET)*SD**2
  110     CONTINUE
          IF(KFL.EQ.1) THEN
            XQ(KFL)=X**TS(1)*(1D0-X)**TS(2)/EULBET(TS(1),TS(2)+1D0)
          ELSE
            XQ(KFL)=TS(1)*X**TS(2)*(1D0-X)**TS(3)*(1D0+TS(4)*X+
     &      TS(5)*X**2)
          ENDIF
  120   CONTINUE
 
C...Put into output array.
        XPPI(0)=XQ(2)
        XPPI(1)=XQ(3)/6D0
        XPPI(2)=XQ(1)+XQ(3)/6D0
        XPPI(3)=XQ(3)/6D0
        XPPI(4)=XQ(4)
        XPPI(-1)=XQ(1)+XQ(3)/6D0
        XPPI(-2)=XQ(3)/6D0
        XPPI(-3)=XQ(3)/6D0
        XPPI(-4)=XQ(4)
 
C...Leading order pion parton distributions from Glueck, Reya and Vogt.
C...Allowed variable range: 0.25 GeV^2 < Q^2 < 10^8 GeV^2 and
C...10^-5 < x < 1.
      ELSE
 
C...Determine s expansion variable and some x expressions.
        VINT(231)=0.25D0
        IF(MSTP(57).LE.0) THEN
          SD=0D0
        ELSE
          Q2IN=MIN(1D8,MAX(0.25D0,Q2))
          SD=LOG(LOG(Q2IN/0.232D0**2)/LOG(0.25D0/0.232D0**2))
        ENDIF
        SD2=SD**2
        XL=-LOG(X)
        XS=SQRT(X)
 
C...Evaluate valence, gluon and sea distributions.
        XFVAL=(0.519D0+0.180D0*SD-0.011D0*SD2)*X**(0.499D0-0.027D0*SD)*
     &  (1D0+(0.381D0-0.419D0*SD)*XS)*(1D0-X)**(0.367D0+0.563D0*SD)
        XFGLU=(X**(0.482D0+0.341D0*SQRT(SD))*((0.678D0+0.877D0*
     &  SD-0.175D0*SD2)+
     &  (0.338D0-1.597D0*SD)*XS+(-0.233D0*SD+0.406D0*SD2)*X)+
     &  SD**0.599D0*EXP(-(0.618D0+2.070D0*SD)+SQRT(3.676D0*SD**1.263D0*
     &  XL)))*
     &  (1D0-X)**(0.390D0+1.053D0*SD)
        XFSEA=SD**0.55D0*(1D0-0.748D0*XS+(0.313D0+0.935D0*SD)*X)*(1D0-
     &  X)**3.359D0*
     &  EXP(-(4.433D0+1.301D0*SD)+SQRT((9.30D0-0.887D0*SD)*SD**0.56D0*
     &  XL))/
     &  XL**(2.538D0-0.763D0*SD)
        IF(SD.LE.0.888D0) THEN
          XFCHM=0D0
        ELSE
          XFCHM=(SD-0.888D0)**1.02D0*(1D0+1.008D0*X)*(1D0-X)**(1.208D0+
     &    0.771D0*SD)*
     &    EXP(-(4.40D0+1.493D0*SD)+SQRT((2.032D0+1.901D0*SD)*SD**0.39D0*
     &    XL))
        ENDIF
        IF(SD.LE.1.351D0) THEN
          XFBOT=0D0
        ELSE
          XFBOT=(SD-1.351D0)**1.03D0*(1D0-X)**(0.697D0+0.855D0*SD)*
     &    EXP(-(4.51D0+1.490D0*SD)+SQRT((3.056D0+1.694D0*SD)*SD**0.39D0*
     &    XL))
        ENDIF
 
C...Put into output array.
        XPPI(0)=XFGLU
        XPPI(1)=XFSEA
        XPPI(2)=XFSEA
        XPPI(3)=XFSEA
        XPPI(4)=XFCHM
        XPPI(5)=XFBOT
        DO 130 KFL=1,5
          XPPI(-KFL)=XPPI(KFL)
  130   CONTINUE
        XPPI(2)=XPPI(2)+XFVAL
        XPPI(-1)=XPPI(-1)+XFVAL
      ENDIF
 
      RETURN
      END
