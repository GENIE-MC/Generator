 
C*********************************************************************
 
C...PYHFTH
C...Gives threshold attractive/repulsive factor for heavy flavour
C...production.
 
      FUNCTION PYHFTH(SH,SQM,FRATT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYDAT1/,/PYPARS/,/PYINT1/
 
C...Value for alpha_strong.
      IF(MSTP(35).LE.1) THEN
        ALSSG=PARP(35)
      ELSE
        MST115=MSTU(115)
        MSTU(115)=MSTP(36)
        Q2BN=SQRT(MAX(1D0,SQM*((SQRT(SH)-2D0*SQRT(SQM))**2+
     &  PARP(36)**2)))
        ALSSG=PYALPS(Q2BN)
        MSTU(115)=MST115
      ENDIF
 
C...Evaluate attractive and repulsive factors.
      XATTR=4D0*PARU(1)*ALSSG/(3D0*SQRT(MAX(1D-20,1D0-4D0*SQM/SH)))
      FATTR=XATTR/(1D0-EXP(-MIN(50D0,XATTR)))
      XREPU=PARU(1)*ALSSG/(6D0*SQRT(MAX(1D-20,1D0-4D0*SQM/SH)))
      FREPU=XREPU/(EXP(MIN(50D0,XREPU))-1D0)
      PYHFTH=FRATT*FATTR+(1D0-FRATT)*FREPU
      VINT(138)=PYHFTH
 
      RETURN
      END
