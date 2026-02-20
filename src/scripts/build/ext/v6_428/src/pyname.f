 
C*********************************************************************
 
C...PYNAME
C...Gives the particle/parton name as a character string.
 
      SUBROUTINE PYNAME(KF,CHAU)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT4/
C...Local character variable.
      CHARACTER CHAU*16
 
C...Read out code with distinction particle/antiparticle.
      CHAU=' '
      KC=PYCOMP(KF)
      IF(KC.NE.0) CHAU=CHAF(KC,(3-ISIGN(1,KF))/2)
 
 
      RETURN
      END
