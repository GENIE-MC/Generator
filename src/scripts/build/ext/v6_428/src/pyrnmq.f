 
C*********************************************************************
 
C...PYRNMQ
C...Determines the running mass of Squarks.
 
      FUNCTION PYRNMQ(ID,DTERM)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblock.
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      SAVE /PYMSSM/
 
C...Local variables.
      DOUBLE PRECISION PI,R
      DOUBLE PRECISION TOL
      DOUBLE PRECISION CI(3)
      EXTERNAL PYALPS
      DOUBLE PRECISION PYALPS
      DATA TOL/0.001D0/
      DATA PI,R/3.141592654D0,.61803399D0/
      DATA CI/0.47D0,0.07D0,0.02D0/
 
      C=1D0-R
      CA=CI(ID)
      AG=(0.71D0)**2/4D0/PI
      AG=RMSS(20)
      XM0=RMSS(8)
      XMG=RMSS(1)
      XM02=XM0*XM0
      XMG2=XMG*XMG
 
      AS=PYALPS(XM02+6D0*XMG2)
      CG=8D0/9D0*((AS/AG)**2-1D0)
      BX=XM02+(CA+CG)*XMG2+DTERM
      AX=MIN(50D0**2,0.5D0*BX)
      CX=MAX(2000D0**2,2D0*BX)
 
      X0=AX
      X3=CX
      IF(ABS(CX-BX).GT.ABS(BX-AX))THEN
        X1=BX
        X2=BX+C*(CX-BX)
      ELSE
        X2=BX
        X1=BX-C*(BX-AX)
      ENDIF
      AS1=PYALPS(X1)
      CG=8D0/9D0*((AS1/AG)**2-1D0)
      F1=ABS(XM02+(CA+CG)*XMG2+DTERM-X1)
      AS2=PYALPS(X2)
      CG=8D0/9D0*((AS2/AG)**2-1D0)
      F2=ABS(XM02+(CA+CG)*XMG2+DTERM-X2)
  100 IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2))) THEN
        IF(F2.LT.F1) THEN
          X0=X1
          X1=X2
          X2=R*X1+C*X3
          F1=F2
          AS2=PYALPS(X2)
          CG=8D0/9D0*((AS2/AG)**2-1D0)
          F2=ABS(XM02+(CA+CG)*XMG2+DTERM-X2)
        ELSE
          X3=X2
          X2=X1
          X1=R*X2+C*X0
          F2=F1
          AS1=PYALPS(X1)
          CG=8D0/9D0*((AS1/AG)**2-1D0)
          F1=ABS(XM02+(CA+CG)*XMG2+DTERM-X1)
        ENDIF
        GOTO 100
      ENDIF
      IF(F1.LT.F2) THEN
        PYRNMQ=X1
        XMIN=X1
      ELSE
        PYRNMQ=X2
        XMIN=X2
      ENDIF
 
      RETURN
      END
