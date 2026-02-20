 
 
C*********************************************************************
 
C...PYPTDI
C...Generates transverse momentum according to a Gaussian.
 
      SUBROUTINE PYPTDI(KFL,PX,PY)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
C...Generate p_T and azimuthal angle, gives p_x and p_y.
      KFLA=IABS(KFL)
      PT=PARJ(21)*SQRT(-LOG(MAX(1D-10,PYR(0))))
      IF(PARJ(23).GT.PYR(0)) PT=PARJ(24)*PT
      IF(MSTJ(91).EQ.1) PT=PARJ(22)*PT
      IF(KFLA.EQ.0.AND.MSTJ(13).LE.0) PT=0D0
      PHI=PARU(2)*PYR(0)
      PX=PT*COS(PHI)
      PY=PT*SIN(PHI)
 
      RETURN
      END
