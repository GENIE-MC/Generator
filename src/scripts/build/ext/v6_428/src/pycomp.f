 
C*********************************************************************
 
C...PYCOMP
C...Compress the standard KF codes for use in mass and decay arrays;
C...also checks whether a given code actually is defined.
 
      FUNCTION PYCOMP(KF)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
C...Local arrays and saved data.
      DIMENSION KFORD(100:500),KCORD(101:500)
      SAVE KFORD,KCORD,NFORD,KFLAST,KCLAST
 
C...Whenever necessary reorder codes for faster search.
      IF(MSTU(20).EQ.0) THEN
        NFORD=100
        KFORD(100)=0
        DO 120 I=101,500
          KFA=KCHG(I,4)
          IF(KFA.LE.100) GOTO 120
          NFORD=NFORD+1
          DO 100 I1=NFORD-1,0,-1
            IF(KFA.GE.KFORD(I1)) GOTO 110
            KFORD(I1+1)=KFORD(I1)
            KCORD(I1+1)=KCORD(I1)
  100     CONTINUE
  110     KFORD(I1+1)=KFA
          KCORD(I1+1)=I
  120   CONTINUE
        MSTU(20)=1
        KFLAST=0
        KCLAST=0
      ENDIF
 
C...Fast action if same code as in latest call.
      IF(KF.EQ.KFLAST) THEN
        PYCOMP=KCLAST
        RETURN
      ENDIF
 
C...Starting values. Remove internal diquark flags.
      PYCOMP=0
      KFA=IABS(KF)
      IF(MOD(KFA/10,10).EQ.0.AND.KFA.LT.100000
     &     .AND.MOD(KFA/1000,10).GT.0) KFA=MOD(KFA,10000)
 
C...Simple cases: direct translation.
      IF(KFA.GT.KFORD(NFORD)) THEN
      ELSEIF(KFA.LE.100) THEN
        PYCOMP=KFA
 
C...Else binary search.
      ELSE
        IMIN=100
        IMAX=NFORD+1
  130   IAVG=(IMIN+IMAX)/2
        IF(KFORD(IAVG).GT.KFA) THEN
          IMAX=IAVG
          IF(IMAX.GT.IMIN+1) GOTO 130
        ELSEIF(KFORD(IAVG).LT.KFA) THEN
          IMIN=IAVG
          IF(IMAX.GT.IMIN+1) GOTO 130
        ELSE
          PYCOMP=KCORD(IAVG)
        ENDIF
      ENDIF
 
C...Check if antiparticle allowed.
      IF(PYCOMP.NE.0.AND.KF.LT.0) THEN
        IF(KCHG(PYCOMP,3).EQ.0) PYCOMP=0
      ENDIF
 
C...Save codes for possible future fast action.
      KFLAST=KF
      KCLAST=PYCOMP
 
      RETURN
      END
