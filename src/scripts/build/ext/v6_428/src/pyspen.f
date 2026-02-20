 
C***********************************************************************
 
C...PYSPEN
C...Calculates real and imaginary part of Spence function; see
C...G. 't Hooft and M. Veltman, Nucl. Phys. B153 (1979) 365.
 
      FUNCTION PYSPEN(XREIN,XIMIN,IREIM)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
C...Local array and data.
      DIMENSION B(0:14)
      DATA B/
     &1.000000D+00,        -5.000000D-01,         1.666667D-01,
     &0.000000D+00,        -3.333333D-02,         0.000000D+00,
     &2.380952D-02,         0.000000D+00,        -3.333333D-02,
     &0.000000D+00,         7.575757D-02,         0.000000D+00,
     &-2.531135D-01,         0.000000D+00,         1.166667D+00/
 
      XRE=XREIN
      XIM=XIMIN
      IF(ABS(1D0-XRE).LT.1D-6.AND.ABS(XIM).LT.1D-6) THEN
        IF(IREIM.EQ.1) PYSPEN=PARU(1)**2/6D0
        IF(IREIM.EQ.2) PYSPEN=0D0
        RETURN
      ENDIF
 
      XMOD=SQRT(XRE**2+XIM**2)
      IF(XMOD.LT.1D-6) THEN
        IF(IREIM.EQ.1) PYSPEN=0D0
        IF(IREIM.EQ.2) PYSPEN=0D0
        RETURN
      ENDIF
 
      XARG=SIGN(ACOS(XRE/XMOD),XIM)
      SP0RE=0D0
      SP0IM=0D0
      SGN=1D0
      IF(XMOD.GT.1D0) THEN
        ALGXRE=LOG(XMOD)
        ALGXIM=XARG-SIGN(PARU(1),XARG)
        SP0RE=-PARU(1)**2/6D0-(ALGXRE**2-ALGXIM**2)/2D0
        SP0IM=-ALGXRE*ALGXIM
        SGN=-1D0
        XMOD=1D0/XMOD
        XARG=-XARG
        XRE=XMOD*COS(XARG)
        XIM=XMOD*SIN(XARG)
      ENDIF
      IF(XRE.GT.0.5D0) THEN
        ALGXRE=LOG(XMOD)
        ALGXIM=XARG
        XRE=1D0-XRE
        XIM=-XIM
        XMOD=SQRT(XRE**2+XIM**2)
        XARG=SIGN(ACOS(XRE/XMOD),XIM)
        ALGYRE=LOG(XMOD)
        ALGYIM=XARG
        SP0RE=SP0RE+SGN*(PARU(1)**2/6D0-(ALGXRE*ALGYRE-ALGXIM*ALGYIM))
        SP0IM=SP0IM-SGN*(ALGXRE*ALGYIM+ALGXIM*ALGYRE)
        SGN=-SGN
      ENDIF
 
      XRE=1D0-XRE
      XIM=-XIM
      XMOD=SQRT(XRE**2+XIM**2)
      XARG=SIGN(ACOS(XRE/XMOD),XIM)
      ZRE=-LOG(XMOD)
      ZIM=-XARG
 
      SPRE=0D0
      SPIM=0D0
      SAVERE=1D0
      SAVEIM=0D0
      DO 100 I=0,14
        IF(MAX(ABS(SAVERE),ABS(SAVEIM)).LT.1D-30) GOTO 110
        TERMRE=(SAVERE*ZRE-SAVEIM*ZIM)/DBLE(I+1)
        TERMIM=(SAVERE*ZIM+SAVEIM*ZRE)/DBLE(I+1)
        SAVERE=TERMRE
        SAVEIM=TERMIM
        SPRE=SPRE+B(I)*TERMRE
        SPIM=SPIM+B(I)*TERMIM
  100 CONTINUE
 
  110 IF(IREIM.EQ.1) PYSPEN=SP0RE+SGN*SPRE
      IF(IREIM.EQ.2) PYSPEN=SP0IM+SGN*SPIM
 
      RETURN
      END
