 
C*********************************************************************
 
C...PYH2XX
C...Calculates the decay rate for a Higgs to an ino pair.
 
      FUNCTION PYH2XX(C1,XM1,XM2,XM3,GX2,GLR)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
C...Local variables.
      DOUBLE PRECISION PYH2XX,XM1,XM2,XM3,GL,GR
      DOUBLE PRECISION XL,PYLAMF,C1
      DOUBLE PRECISION XMI2,XMJ2,XMK2,XMI3
 
      XMI2=XM1**2
      XMI3=ABS(XM1**3)
      XMJ2=XM2**2
      XMK2=XM3**2
      XL=PYLAMF(XMI2,XMJ2,XMK2)
      PYH2XX=C1/4D0/XMI3*SQRT(XL)
     &*(GX2*(XMI2-XMJ2-XMK2)-
     &4D0*GLR*XM3*XM2)
      IF(PYH2XX.LT.0D0) PYH2XX=0D0
 
      RETURN
      END
