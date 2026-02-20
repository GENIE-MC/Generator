 
C*********************************************************************
 
C...PYFISB
C...Auxiliary routine to PYFINT for SUSY Higgs calculations.
 
      FUNCTION PYFISB(X)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblock.
      COMMON/PYINTS/XXM(20)
      SAVE/PYINTS/
 
      PYFISB = LOG(ABS(X*XXM(2)+(1-X)*XXM(3)-X*(1-X)*XXM(1))/
     &(X*(XXM(2)-XXM(3))+XXM(3)))
 
      RETURN
      END
