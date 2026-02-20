 
C*********************************************************************
 
C...PYMRUN
C...Gives the running, current-algebra mass of a d, u, s, c or b quark,
C...for Higgs couplings. Everything else sent on to PYMASS.
 
      FUNCTION PYMRUN(KF,Q2)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /PYDAT1/,/PYDAT2/,/PYPARS/
 
C...Most masses not handled here.
      KFA=IABS(KF)
      IF(KFA.EQ.0.OR.KFA.GT.6) THEN
        PYMRUN=PYMASS(KF)
 
C...Current-algebra masses, but no Q2 dependence.
      ELSEIF(MSTP(37).NE.1.OR.MSTP(2).LE.0) THEN
        PYMRUN=PARF(90+KFA)
 
C...Running current-algebra masses.
      ELSE
        AS=PYALPS(Q2)
        PYMRUN=PARF(90+KFA)*
     &  (LOG(MAX(4D0,PARP(37)**2*PARF(90+KFA)**2/PARU(117)**2))/
     &  LOG(MAX(4D0,Q2/PARU(117)**2)))**(12D0/(33D0-2D0*MSTU(118)))
      ENDIF
 
      RETURN
      END
