 
C*********************************************************************
 
C...PYEEVT
C...Handles the generation of an e+e- annihilation jet event.
 
      SUBROUTINE PYEEVT(KFL,ECM)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
 
C...Check input parameters.
      IF(MSTU(12).NE.12345) CALL PYLIST(0)
      IF(KFL.LT.0.OR.KFL.GT.8) THEN
        CALL PYERRM(16,'(PYEEVT:) called with unknown flavour code')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      IF(KFL.LE.5) ECMMIN=PARJ(127)+2.02D0*PARF(100+MAX(1,KFL))
      IF(KFL.GE.6) ECMMIN=PARJ(127)+2.02D0*PMAS(KFL,1)
      IF(ECM.LT.ECMMIN) THEN
        CALL PYERRM(16,'(PYEEVT:) called with too small CM energy')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
 
C...Check consistency of MSTJ options set.
      IF(MSTJ(109).EQ.2.AND.MSTJ(110).NE.1) THEN
        CALL PYERRM(6,
     &  '(PYEEVT:) MSTJ(109) value requires MSTJ(110) = 1')
        MSTJ(110)=1
      ENDIF
      IF(MSTJ(109).EQ.2.AND.MSTJ(111).NE.0) THEN
        CALL PYERRM(6,
     &  '(PYEEVT:) MSTJ(109) value requires MSTJ(111) = 0')
        MSTJ(111)=0
      ENDIF
 
C...Initialize alpha_strong and total cross-section.
      MSTU(111)=MSTJ(108)
      IF(MSTJ(108).EQ.2.AND.(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.1))
     &MSTU(111)=1
      PARU(112)=PARJ(121)
      IF(MSTU(111).EQ.2) PARU(112)=PARJ(122)
      IF(MSTJ(116).GT.0.AND.(MSTJ(116).GE.2.OR.ABS(ECM-PARJ(151)).GE.
     &PARJ(139).OR.10*MSTJ(102)+KFL.NE.MSTJ(119))) CALL PYXTEE(KFL,ECM,
     &XTOT)
      IF(MSTJ(116).GE.3) MSTJ(116)=1
      PARJ(171)=0D0
 
C...Add initial e+e- to event record (documentation only).
      NTRY=0
  100 NTRY=NTRY+1
      IF(NTRY.GT.100) THEN
        CALL PYERRM(14,'(PYEEVT:) caught in an infinite loop')
        RETURN
      ENDIF
      MSTU(24)=0
      NC=0
      IF(MSTJ(115).GE.2) THEN
        NC=NC+2
        CALL PY1ENT(NC-1,11,0.5D0*ECM,0D0,0D0)
        K(NC-1,1)=21
        CALL PY1ENT(NC,-11,0.5D0*ECM,PARU(1),0D0)
        K(NC,1)=21
      ENDIF
 
C...Radiative photon (in initial state).
      MK=0
      ECMC=ECM
      IF(MSTJ(107).GE.1.AND.MSTJ(116).GE.1) CALL PYRADK(ECM,MK,PAK,
     &THEK,PHIK,ALPK)
      IF(MK.EQ.1) ECMC=SQRT(ECM*(ECM-2D0*PAK))
      IF(MSTJ(115).GE.1.AND.MK.EQ.1) THEN
        NC=NC+1
        CALL PY1ENT(NC,22,PAK,THEK,PHIK)
        K(NC,3)=MIN(MSTJ(115)/2,1)
      ENDIF
 
C...Virtual exchange boson (gamma or Z0).
      IF(MSTJ(115).GE.3) THEN
        NC=NC+1
        KF=22
        IF(MSTJ(102).EQ.2) KF=23
        MSTU10=MSTU(10)
        MSTU(10)=1
        P(NC,5)=ECMC
        CALL PY1ENT(NC,KF,ECMC,0D0,0D0)
        K(NC,1)=21
        K(NC,3)=1
        MSTU(10)=MSTU10
      ENDIF
 
C...Choice of flavour and jet configuration.
      CALL PYXKFL(KFL,ECM,ECMC,KFLC)
      IF(KFLC.EQ.0) GOTO 100
      CALL PYXJET(ECMC,NJET,CUT)
      KFLN=21
      IF(NJET.EQ.4) CALL PYX4JT(NJET,CUT,KFLC,ECMC,KFLN,X1,X2,X4,
     &X12,X14)
      IF(NJET.EQ.3) CALL PYX3JT(NJET,CUT,KFLC,ECMC,X1,X3)
      IF(NJET.EQ.2) MSTJ(120)=1
 
C...Fill jet configuration and origin.
      IF(NJET.EQ.2.AND.MSTJ(101).NE.5) CALL PY2ENT(NC+1,KFLC,-KFLC,ECMC)
      IF(NJET.EQ.2.AND.MSTJ(101).EQ.5) CALL PY2ENT(-(NC+1),KFLC,-KFLC,
     &ECMC)
      IF(NJET.EQ.3) CALL PY3ENT(NC+1,KFLC,21,-KFLC,ECMC,X1,X3)
      IF(NJET.EQ.4.AND.KFLN.EQ.21) CALL PY4ENT(NC+1,KFLC,KFLN,KFLN,
     &-KFLC,ECMC,X1,X2,X4,X12,X14)
      IF(NJET.EQ.4.AND.KFLN.NE.21) CALL PY4ENT(NC+1,KFLC,-KFLN,KFLN,
     &-KFLC,ECMC,X1,X2,X4,X12,X14)
      IF(MSTU(24).NE.0) GOTO 100
      DO 110 IP=NC+1,N
        K(IP,3)=K(IP,3)+MIN(MSTJ(115)/2,1)+(MSTJ(115)/3)*(NC-1)
  110 CONTINUE
 
C...Angular orientation according to matrix element.
      IF(MSTJ(106).EQ.1) THEN
        CALL PYXDIF(NC,NJET,KFLC,ECMC,CHI,THE,PHI)
        CALL PYROBO(NC+1,N,0D0,CHI,0D0,0D0,0D0)
        CALL PYROBO(NC+1,N,THE,PHI,0D0,0D0,0D0)
      ENDIF
 
C...Rotation and boost from radiative photon.
      IF(MK.EQ.1) THEN
        DBEK=-PAK/(ECM-PAK)
        NMIN=NC+1-MSTJ(115)/3
        CALL PYROBO(NMIN,N,0D0,-PHIK,0D0,0D0,0D0)
        CALL PYROBO(NMIN,N,ALPK,0D0,DBEK*SIN(THEK),0D0,DBEK*COS(THEK))
        CALL PYROBO(NMIN,N,0D0,PHIK,0D0,0D0,0D0)
      ENDIF
 
C...Generate parton shower. Rearrange along strings and check.
      IF(MSTJ(101).EQ.5) THEN
        CALL PYSHOW(N-1,N,ECMC)
        MSTJ14=MSTJ(14)
        IF(MSTJ(105).EQ.-1) MSTJ(14)=-1
        IF(MSTJ(105).GE.0) MSTU(28)=0
        CALL PYPREP(0)
        MSTJ(14)=MSTJ14
        IF(MSTJ(105).GE.0.AND.MSTU(28).NE.0) GOTO 100
      ENDIF
 
C...Fragmentation/decay generation. Information for PYTABU.
      IF(MSTJ(105).EQ.1) CALL PYEXEC
      MSTU(161)=KFLC
      MSTU(162)=-KFLC
 
      RETURN
      END
