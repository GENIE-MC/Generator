 
 
 
 
 
C*********************************************************************
 
C...PYFINT
C...Auxiliary routine to PYPOLE for SUSY Higgs calculations.
 
      FUNCTION PYFINT(A,B,C)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblock.
      COMMON/PYINTS/XXM(20)
      SAVE/PYINTS/
 
C...Local variables.
      EXTERNAL PYFISB
      DOUBLE PRECISION PYFISB
 
      XXM(1)=A
      XXM(2)=B
      XXM(3)=C
      XLO=0D0
      XHI=1D0
      PYFINT  = PYGAUS(PYFISB,XLO,XHI,1D-3)
 
      RETURN
      END
