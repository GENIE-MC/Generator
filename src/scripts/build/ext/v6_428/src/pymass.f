 
C*********************************************************************
 
C...PYMASS
C...Gives the mass of a particle/parton.
 
      FUNCTION PYMASS(KF)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 
C...Reset variables. Compressed code. Special case for popcorn diquarks.
      PYMASS=0D0
      KFA=IABS(KF)
      KC=PYCOMP(KF)
      IF(KC.EQ.0) THEN
        MSTJ(93)=0
        RETURN
      ENDIF
 
C...Guarantee use of constituent masses for internal checks.
      IF((MSTJ(93).EQ.1.OR.MSTJ(93).EQ.2).AND.
     &(KFA.LE.10.OR.MOD(KFA/10,10).EQ.0)) THEN
        IF(KFA.LE.5) THEN
          PYMASS=PARF(100+KFA)
          IF(MSTJ(93).EQ.2) PYMASS=MAX(0D0,PYMASS-PARF(121))
        ELSEIF(KFA.LE.10) THEN
          PYMASS=PMAS(KFA,1)
        ELSEIF(MSTJ(93).EQ.1) THEN
          PYMASS=PARF(100+MOD(KFA/1000,10))+PARF(100+MOD(KFA/100,10))
        ELSE
          PYMASS=MAX(0D0,PMAS(KC,1)-PARF(122)-2D0*PARF(112)/3D0)
        ENDIF
 
C...Other masses can be read directly off table.
      ELSE
        PYMASS=PMAS(KC,1)
      ENDIF
 
C...Optional mass broadening according to truncated Breit-Wigner
C...(either in m or in m^2).
      IF(MSTJ(24).GE.1.AND.PMAS(KC,2).GT.1D-4) THEN
        IF(MSTJ(24).EQ.1.OR.(MSTJ(24).EQ.2.AND.KFA.GT.100)) THEN
          PYMASS=PYMASS+0.5D0*PMAS(KC,2)*TAN((2D0*PYR(0)-1D0)*
     &    ATAN(2D0*PMAS(KC,3)/PMAS(KC,2)))
        ELSE
          PM0=PYMASS
          PMLOW=ATAN((MAX(0D0,PM0-PMAS(KC,3))**2-PM0**2)/
     &    (PM0*PMAS(KC,2)))
          PMUPP=ATAN(((PM0+PMAS(KC,3))**2-PM0**2)/(PM0*PMAS(KC,2)))
          PYMASS=SQRT(MAX(0D0,PM0**2+PM0*PMAS(KC,2)*TAN(PMLOW+
     &    (PMUPP-PMLOW)*PYR(0))))
        ENDIF
      ENDIF
      MSTJ(93)=0
 
      RETURN
      END
