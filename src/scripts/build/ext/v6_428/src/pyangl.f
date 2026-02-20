 
C*********************************************************************
 
C...PYANGL
C...Reconstructs an angle from given x and y coordinates.
 
      FUNCTION PYANGL(X,Y)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
      PYANGL=0D0
      R=SQRT(X**2+Y**2)
      IF(R.LT.1D-20) RETURN
      IF(ABS(X)/R.LT.0.8D0) THEN
        PYANGL=SIGN(ACOS(X/R),Y)
      ELSE
        PYANGL=ASIN(Y/R)
        IF(X.LT.0D0.AND.PYANGL.GE.0D0) THEN
          PYANGL=PARU(1)-PYANGL
        ELSEIF(X.LT.0D0) THEN
          PYANGL=-PARU(1)-PYANGL
        ENDIF
      ENDIF
 
      RETURN
      END
