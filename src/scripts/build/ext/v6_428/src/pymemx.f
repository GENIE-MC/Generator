 
C*********************************************************************
 
C...PYMEMX
C...Generates maximum ME weight in some initial-state showers.
C...Inparameter MECOR: kind of hard scattering process
C...Outparameter WTFF: maximum weight for fermion -> fermion
C...             WTGF: maximum weight for gluon/photon -> fermion
C...             WTFG: maximum weight for fermion -> gluon/photon
C...             WTGG: maximum weight for gluon -> gluon
 
      SUBROUTINE PYMEMX(MECOR,WTFF,WTGF,WTFG,WTGG)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      SAVE /PYJETS/,/PYDAT1/,/PYPARS/,/PYINT1/,/PYINT2/
 
C...Default maximum weight.
      WTFF=1D0
      WTGF=1D0
      WTFG=1D0
      WTGG=1D0
 
C...Select maximum weight by process.
      IF(MECOR.EQ.1) THEN
        WTFF=1D0
        WTGF=3D0
      ELSEIF(MECOR.EQ.2) THEN
        WTFG=1D0
        WTGG=1D0
      ENDIF
 
      RETURN
      END
