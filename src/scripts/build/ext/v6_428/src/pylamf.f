 
C*********************************************************************
 
C...PYLAMF
C...The standard lambda function.
 
      FUNCTION PYLAMF(X,Y,Z)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
 
C...Local variables.
      DOUBLE PRECISION PYLAMF,X,Y,Z
 
      PYLAMF=(X-(Y+Z))**2-4D0*Y*Z
      IF(PYLAMF.LT.0D0) PYLAMF=0D0
 
      RETURN
      END
