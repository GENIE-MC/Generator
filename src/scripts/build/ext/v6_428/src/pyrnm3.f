 
C*********************************************************************
 
C...PYRNM3
C...Calculates the running of M3, the SU(3) gluino mass parameter.
 
      FUNCTION PYRNM3(RGUT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
 
C...Local variables.
      DOUBLE PRECISION R
      DOUBLE PRECISION TOL
      EXTERNAL PYALPS
      DOUBLE PRECISION PYALPS
      DATA TOL/0.001D0/
      DATA R/0.61803399D0/
 
      C=1D0-R
 
      BX=RGUT*PYALPS(RGUT**2)
      AX=MIN(50D0,BX*0.5D0)
      CX=MAX(2000D0,2D0*BX)
 
      X0=AX
      X3=CX
      IF(ABS(CX-BX).GT.ABS(BX-AX))THEN
        X1=BX
        X2=BX+C*(CX-BX)
      ELSE
        X2=BX
        X1=BX-C*(BX-AX)
      ENDIF
      AS1=PYALPS(X1**2)
      F1=ABS(X1-RGUT*AS1)
      AS2=PYALPS(X2**2)
      F2=ABS(X2-RGUT*AS2)
  100 IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2))) THEN
        IF(F2.LT.F1) THEN
          X0=X1
          X1=X2
          X2=R*X1+C*X3
          F1=F2
          AS2=PYALPS(X2**2)
          F2=ABS(X2-RGUT*AS2)
        ELSE
          X3=X2
          X2=X1
          X1=R*X2+C*X0
          F2=F1
          AS1=PYALPS(X1**2)
          F1=ABS(X1-RGUT*AS1)
        ENDIF
        GOTO 100
      ENDIF
      IF(F1.LT.F2) THEN
        PYRNM3=X1
        XMIN=X1
      ELSE
        PYRNM3=X2
        XMIN=X2
      ENDIF
 
      RETURN
      END
