 
C*********************************************************************
 
C...PYCHGE
C...Gives three times the charge for a particle/parton.
 
      FUNCTION PYCHGE(KF)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT2/
 
C...Read out charge and change sign for antiparticle.
      PYCHGE=0
      KC=PYCOMP(KF)
      IF(KC.NE.0) PYCHGE=KCHG(KC,1)*ISIGN(1,KF)
 
      RETURN
      END
