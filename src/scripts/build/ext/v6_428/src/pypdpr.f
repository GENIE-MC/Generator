 
C*********************************************************************
 
C...PYPDPR
C...Gives proton parton distributions according to a few different
C...parametrizations.
 
      SUBROUTINE PYPDPR(X,Q2,XPPR)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/
C...Arrays and data.
      DIMENSION XPPR(-6:6),Q2MIN(16)
      DATA Q2MIN/ 2.56D0, 2.56D0, 2.56D0, 0.4D0, 0.4D0, 0.4D0,
     &1.0D0, 1.0D0, 2*0D0, 0.25D0, 5D0, 5D0, 4D0, 4D0, 0D0/
 
C...Reset output array.
      DO 100 KFL=-6,6
        XPPR(KFL)=0D0
  100 CONTINUE
 
C...Common preliminaries.
      NSET=MAX(1,MIN(16,MSTP(51)))
      IF(NSET.EQ.9.OR.NSET.EQ.10) NSET=6
      VINT(231)=Q2MIN(NSET)
      IF(MSTP(57).EQ.0) THEN
        Q2L=Q2MIN(NSET)
      ELSE
        Q2L=MAX(Q2MIN(NSET),Q2)
      ENDIF
 
      IF(NSET.GE.1.AND.NSET.LE.3) THEN
C...Interface to the CTEQ 3 parton distributions.
        QRT=SQRT(MAX(1D0,Q2L))
 
C...Loop over flavours.
        DO 110 I=-6,6
          IF(I.LE.0) THEN
            XPPR(I)=PYCTEQ(NSET,I,X,QRT)
          ELSEIF(I.LE.2) THEN
            XPPR(I)=PYCTEQ(NSET,I,X,QRT)+XPPR(-I)
          ELSE
            XPPR(I)=XPPR(-I)
          ENDIF
  110   CONTINUE
 
      ELSEIF(NSET.GE.4.AND.NSET.LE.6) THEN
C...Interface to the GRV 94 distributions.
        IF(NSET.EQ.4) THEN
          CALL PYGRVL (X, Q2L, UV, DV, DEL, UDB, SB, CHM, BOT, GL)
        ELSEIF(NSET.EQ.5) THEN
          CALL PYGRVM (X, Q2L, UV, DV, DEL, UDB, SB, CHM, BOT, GL)
        ELSE
          CALL PYGRVD (X, Q2L, UV, DV, DEL, UDB, SB, CHM, BOT, GL)
        ENDIF
 
C...Put into output array.
        XPPR(0)=GL
        XPPR(-1)=0.5D0*(UDB+DEL)
        XPPR(-2)=0.5D0*(UDB-DEL)
        XPPR(-3)=SB
        XPPR(-4)=CHM
        XPPR(-5)=BOT
        XPPR(1)=DV+XPPR(-1)
        XPPR(2)=UV+XPPR(-2)
        XPPR(3)=SB
        XPPR(4)=CHM
        XPPR(5)=BOT
 
      ELSEIF(NSET.EQ.7) THEN
C...Interface to the CTEQ 5L parton distributions.
C...Range of validity 10^-6 < x < 1, 1 < Q < 10^4 extended by
C...freezing x*f(x,Q2) at borders.
        QRT=SQRT(MAX(1D0,MIN(1D8,Q2L)))
        XIN=MAX(1D-6,MIN(1D0,X))
 
C...Loop over flavours (with u <-> d notation mismatch).
        SUMUDB=PYCT5L(-1,XIN,QRT)
        RATUDB=PYCT5L(-2,XIN,QRT)
        DO 120 I=-5,2
          IF(I.EQ.1) THEN
            XPPR(I)=XIN*PYCT5L(2,XIN,QRT)
          ELSEIF(I.EQ.2) THEN
            XPPR(I)=XIN*PYCT5L(1,XIN,QRT)
          ELSEIF(I.EQ.-1) THEN
            XPPR(I)=XIN*SUMUDB*RATUDB/(1D0+RATUDB)
          ELSEIF(I.EQ.-2) THEN
            XPPR(I)=XIN*SUMUDB/(1D0+RATUDB)
          ELSE
            XPPR(I)=XIN*PYCT5L(I,XIN,QRT)
            IF(I.LT.0) XPPR(-I)=XPPR(I)
          ENDIF
  120   CONTINUE
 
      ELSEIF(NSET.EQ.8) THEN
C...Interface to the CTEQ 5M1 parton distributions.
        QRT=SQRT(MAX(1D0,MIN(1D8,Q2L)))
        XIN=MAX(1D-6,MIN(1D0,X))
 
C...Loop over flavours (with u <-> d notation mismatch).
        SUMUDB=PYCT5M(-1,XIN,QRT)
        RATUDB=PYCT5M(-2,XIN,QRT)
        DO 130 I=-5,2
          IF(I.EQ.1) THEN
            XPPR(I)=XIN*PYCT5M(2,XIN,QRT)
          ELSEIF(I.EQ.2) THEN
            XPPR(I)=XIN*PYCT5M(1,XIN,QRT)
          ELSEIF(I.EQ.-1) THEN
            XPPR(I)=XIN*SUMUDB*RATUDB/(1D0+RATUDB)
          ELSEIF(I.EQ.-2) THEN
            XPPR(I)=XIN*SUMUDB/(1D0+RATUDB)
          ELSE
            XPPR(I)=XIN*PYCT5M(I,XIN,QRT)
            IF(I.LT.0) XPPR(-I)=XPPR(I)
          ENDIF
  130   CONTINUE
 
      ELSEIF(NSET.GE.11.AND.NSET.LE.15) THEN
C...GRV92LO, EHLQ1, EHLQ2, DO1 AND DO2 distributions:
C...obsolete but offers backwards compatibility.
        CALL PYPDPO(X,Q2L,XPPR)
 
C...Symmetric choice for debugging only
      ELSEIF(NSET.EQ.16) THEN
        XPPR(0)=.5D0/X
        XPPR(1)=.05D0/X
        XPPR(2)=.05D0/X
        XPPR(3)=.05D0/X
        XPPR(4)=.05D0/X
        XPPR(5)=.05D0/X
        XPPR(-1)=.05D0/X
        XPPR(-2)=.05D0/X
        XPPR(-3)=.05D0/X
        XPPR(-4)=.05D0/X
        XPPR(-5)=.05D0/X
 
      ENDIF
 
      RETURN
      END
