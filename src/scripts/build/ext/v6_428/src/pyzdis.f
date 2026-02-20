 
C*********************************************************************
 
C...PYZDIS
C...Generates the longitudinal splitting variable z.
 
      SUBROUTINE PYZDIS(KFL1,KFL2,PR,Z)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 
C...Check if heavy flavour fragmentation.
      KFLA=IABS(KFL1)
      KFLB=IABS(KFL2)
      KFLH=KFLA
      IF(KFLA.GE.10) KFLH=MOD(KFLA/1000,10)
 
C...Lund symmetric scaling function: determine parameters of shape.
      IF(MSTJ(11).EQ.1.OR.(MSTJ(11).EQ.3.AND.KFLH.LE.3).OR.
     &MSTJ(11).GE.4) THEN
        FA=PARJ(41)
        IF(MSTJ(91).EQ.1) FA=PARJ(43)
        IF(KFLB.GE.10) FA=FA+PARJ(45)
        FBB=PARJ(42)
        IF(MSTJ(91).EQ.1) FBB=PARJ(44)
        FB=FBB*PR
        FC=1D0
        IF(KFLA.GE.10) FC=FC-PARJ(45)
        IF(KFLB.GE.10) FC=FC+PARJ(45)
        IF(MSTJ(11).GE.4.AND.(KFLH.EQ.4.OR.KFLH.EQ.5)) THEN
          FRED=PARJ(46)
          IF(MSTJ(11).EQ.5.AND.KFLH.EQ.5) FRED=PARJ(47)
          FC=FC+FRED*FBB*PARF(100+KFLH)**2
        ENDIF
        MC=1
        IF(ABS(FC-1D0).GT.0.01D0) MC=2
 
C...Determine position of maximum. Special cases for a = 0 or a = c.
        IF(FA.LT.0.02D0) THEN
          MA=1
          ZMAX=1D0
          IF(FC.GT.FB) ZMAX=FB/FC
        ELSEIF(ABS(FC-FA).LT.0.01D0) THEN
          MA=2
          ZMAX=FB/(FB+FC)
        ELSE
          MA=3
          ZMAX=0.5D0*(FB+FC-SQRT((FB-FC)**2+4D0*FA*FB))/(FC-FA)
          IF(ZMAX.GT.0.9999D0.AND.FB.GT.100D0) ZMAX=MIN(ZMAX,1D0-FA/FB)
        ENDIF
 
C...Subdivide z range if distribution very peaked near endpoint.
        MMAX=2
        IF(ZMAX.LT.0.1D0) THEN
          MMAX=1
          ZDIV=2.75D0*ZMAX
          IF(MC.EQ.1) THEN
            FINT=1D0-LOG(ZDIV)
          ELSE
            ZDIVC=ZDIV**(1D0-FC)
            FINT=1D0+(1D0-1D0/ZDIVC)/(FC-1D0)
          ENDIF
        ELSEIF(ZMAX.GT.0.85D0.AND.FB.GT.1D0) THEN
          MMAX=3
          FSCB=SQRT(4D0+(FC/FB)**2)
          ZDIV=FSCB-1D0/ZMAX-(FC/FB)*LOG(ZMAX*0.5D0*(FSCB+FC/FB))
          IF(MA.GE.2) ZDIV=ZDIV+(FA/FB)*LOG(1D0-ZMAX)
          ZDIV=MIN(ZMAX,MAX(0D0,ZDIV))
          FINT=1D0+FB*(1D0-ZDIV)
        ENDIF
 
C...Choice of z, preweighted for peaks at low or high z.
  100   Z=PYR(0)
        FPRE=1D0
        IF(MMAX.EQ.1) THEN
          IF(FINT*PYR(0).LE.1D0) THEN
            Z=ZDIV*Z
          ELSEIF(MC.EQ.1) THEN
            Z=ZDIV**Z
            FPRE=ZDIV/Z
          ELSE
            Z=(ZDIVC+Z*(1D0-ZDIVC))**(1D0/(1D0-FC))
            FPRE=(ZDIV/Z)**FC
          ENDIF
        ELSEIF(MMAX.EQ.3) THEN
          IF(FINT*PYR(0).LE.1D0) THEN
            Z=ZDIV+LOG(Z)/FB
            FPRE=EXP(FB*(Z-ZDIV))
          ELSE
            Z=ZDIV+Z*(1D0-ZDIV)
          ENDIF
        ENDIF
 
C...Weighting according to correct formula.
        IF(Z.LE.0D0.OR.Z.GE.1D0) GOTO 100
        FEXP=FC*LOG(ZMAX/Z)+FB*(1D0/ZMAX-1D0/Z)
        IF(MA.GE.2) FEXP=FEXP+FA*LOG((1D0-Z)/(1D0-ZMAX))
        FVAL=EXP(MAX(-50D0,MIN(50D0,FEXP)))
        IF(FVAL.LT.PYR(0)*FPRE) GOTO 100
 
C...Generate z according to Field-Feynman, SLAC, (1-z)**c OR z**c.
      ELSE
        FC=PARJ(50+MAX(1,KFLH))
        IF(MSTJ(91).EQ.1) FC=PARJ(59)
  110   Z=PYR(0)
        IF(FC.GE.0D0.AND.FC.LE.1D0) THEN
          IF(FC.GT.PYR(0)) Z=1D0-Z**(1D0/3D0)
        ELSEIF(FC.GT.-1.AND.FC.LT.0D0) THEN
          IF(-4D0*FC*Z*(1D0-Z)**2.LT.PYR(0)*((1D0-Z)**2-FC*Z)**2)
     &    GOTO 110
        ELSE
          IF(FC.GT.0D0) Z=1D0-Z**(1D0/FC)
          IF(FC.LT.0D0) Z=Z**(-1D0/FC)
        ENDIF
      ENDIF
 
      RETURN
      END
