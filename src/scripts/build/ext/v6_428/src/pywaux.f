 
C***********************************************************************
 
C...PYWAUX
C...Calculates real and imaginary parts of the auxiliary functions W1
C...and W2; see R. K. Ellis, I. Hinchliffe, M. Soldate and J. J. van
C...der Bij, Nucl. Phys. B297 (1988) 221.
 
      SUBROUTINE PYWAUX(IAUX,EPS,WRE,WIM)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
      ASINH(X)=LOG(X+SQRT(X**2+1D0))
      ACOSH(X)=LOG(X+SQRT(X**2-1D0))
 
      IF(EPS.LT.0D0) THEN
        IF(IAUX.EQ.1) WRE=2D0*SQRT(1D0-EPS)*ASINH(SQRT(-1D0/EPS))
        IF(IAUX.EQ.2) WRE=4D0*(ASINH(SQRT(-1D0/EPS)))**2
        WIM=0D0
      ELSEIF(EPS.LT.1D0) THEN
        IF(IAUX.EQ.1) WRE=2D0*SQRT(1D0-EPS)*ACOSH(SQRT(1D0/EPS))
        IF(IAUX.EQ.2) WRE=4D0*(ACOSH(SQRT(1D0/EPS)))**2-PARU(1)**2
        IF(IAUX.EQ.1) WIM=-PARU(1)*SQRT(1D0-EPS)
        IF(IAUX.EQ.2) WIM=-4D0*PARU(1)*ACOSH(SQRT(1D0/EPS))
      ELSE
        IF(IAUX.EQ.1) WRE=2D0*SQRT(EPS-1D0)*ASIN(SQRT(1D0/EPS))
        IF(IAUX.EQ.2) WRE=-4D0*(ASIN(SQRT(1D0/EPS)))**2
        WIM=0D0
      ENDIF
 
      RETURN
      END
