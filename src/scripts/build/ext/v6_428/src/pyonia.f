 
C*********************************************************************
 
C...PYONIA
C...Generates Upsilon and toponium decays into three gluons
C...or two gluons and a photon.
 
      SUBROUTINE PYONIA(KFL,ECM)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
 
C...Printout. Check input parameters.
      IF(MSTU(12).NE.12345) CALL PYLIST(0)
      IF(KFL.LT.0.OR.KFL.GT.8) THEN
        CALL PYERRM(16,'(PYONIA:) called with unknown flavour code')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      IF(ECM.LT.PARJ(127)+2.02D0*PARF(101)) THEN
        CALL PYERRM(16,'(PYONIA:) called with too small CM energy')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
 
C...Initial e+e- and onium state (optional).
      NC=0
      IF(MSTJ(115).GE.2) THEN
        NC=NC+2
        CALL PY1ENT(NC-1,11,0.5D0*ECM,0D0,0D0)
        K(NC-1,1)=21
        CALL PY1ENT(NC,-11,0.5D0*ECM,PARU(1),0D0)
        K(NC,1)=21
      ENDIF
      KFLC=IABS(KFL)
      IF(MSTJ(115).GE.3.AND.KFLC.GE.5) THEN
        NC=NC+1
        KF=110*KFLC+3
        MSTU10=MSTU(10)
        MSTU(10)=1
        P(NC,5)=ECM
        CALL PY1ENT(NC,KF,ECM,0D0,0D0)
        K(NC,1)=21
        K(NC,3)=1
        MSTU(10)=MSTU10
      ENDIF
 
C...Choose x1 and x2 according to matrix element.
      NTRY=0
  100 X1=PYR(0)
      X2=PYR(0)
      X3=2D0-X1-X2
      IF(X3.GE.1D0.OR.((1D0-X1)/(X2*X3))**2+((1D0-X2)/(X1*X3))**2+
     &((1D0-X3)/(X1*X2))**2.LE.2D0*PYR(0)) GOTO 100
      NTRY=NTRY+1
      NJET=3
      IF(MSTJ(101).LE.4) CALL PY3ENT(NC+1,21,21,21,ECM,X1,X3)
      IF(MSTJ(101).GE.5) CALL PY3ENT(-(NC+1),21,21,21,ECM,X1,X3)
 
C...Photon-gluon-gluon events. Small system modifications. Jet origin.
      MSTU(111)=MSTJ(108)
      IF(MSTJ(108).EQ.2.AND.(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.1))
     &MSTU(111)=1
      PARU(112)=PARJ(121)
      IF(MSTU(111).EQ.2) PARU(112)=PARJ(122)
      QF=0D0
      IF(KFLC.NE.0) QF=KCHG(KFLC,1)/3D0
      RGAM=7.2D0*QF**2*PARU(101)/PYALPS(ECM**2)
      MK=0
      ECMC=ECM
      IF(PYR(0).GT.RGAM/(1D0+RGAM)) THEN
        IF(1D0-MAX(X1,X2,X3).LE.MAX((PARJ(126)/ECM)**2,PARJ(125)))
     &  NJET=2
        IF(NJET.EQ.2.AND.MSTJ(101).LE.4) CALL PY2ENT(NC+1,21,21,ECM)
        IF(NJET.EQ.2.AND.MSTJ(101).GE.5) CALL PY2ENT(-(NC+1),21,21,ECM)
      ELSE
        MK=1
        ECMC=SQRT(1D0-X1)*ECM
        IF(ECMC.LT.2D0*PARJ(127)) GOTO 100
        K(NC+1,1)=1
        K(NC+1,2)=22
        K(NC+1,4)=0
        K(NC+1,5)=0
        IF(MSTJ(101).GE.5) K(NC+2,4)=MSTU(5)*(NC+3)
        IF(MSTJ(101).GE.5) K(NC+2,5)=MSTU(5)*(NC+3)
        IF(MSTJ(101).GE.5) K(NC+3,4)=MSTU(5)*(NC+2)
        IF(MSTJ(101).GE.5) K(NC+3,5)=MSTU(5)*(NC+2)
        NJET=2
        IF(ECMC.LT.4D0*PARJ(127)) THEN
          MSTU10=MSTU(10)
          MSTU(10)=1
          P(NC+2,5)=ECMC
          CALL PY1ENT(NC+2,83,0.5D0*(X2+X3)*ECM,PARU(1),0D0)
          MSTU(10)=MSTU10
          NJET=0
        ENDIF
      ENDIF
      DO 110 IP=NC+1,N
        K(IP,3)=K(IP,3)+(MSTJ(115)/2)+(KFLC/5)*(MSTJ(115)/3)*(NC-1)
  110 CONTINUE
 
C...Differential cross-sections. Upper limit for cross-section.
      IF(MSTJ(106).EQ.1) THEN
        SQ2=SQRT(2D0)
        HF1=1D0-PARJ(131)*PARJ(132)
        HF3=PARJ(133)**2
        CT13=(X1*X3-2D0*X1-2D0*X3+2D0)/(X1*X3)
        ST13=SQRT(1D0-CT13**2)
        SIGL=0.5D0*X3**2*((1D0-X2)**2+(1D0-X3)**2)*ST13**2
        SIGU=(X1*(1D0-X1))**2+(X2*(1D0-X2))**2+(X3*(1D0-X3))**2-SIGL
        SIGT=0.5D0*SIGL
        SIGI=(SIGL*CT13/ST13+0.5D0*X1*X3*(1D0-X2)**2*ST13)/SQ2
        SIGMAX=(2D0*HF1+HF3)*ABS(SIGU)+2D0*(HF1+HF3)*ABS(SIGL)+2D0*(HF1+
     &  2D0*HF3)*ABS(SIGT)+2D0*SQ2*(HF1+2D0*HF3)*ABS(SIGI)
 
C...Angular orientation of event.
  120   CHI=PARU(2)*PYR(0)
        CTHE=2D0*PYR(0)-1D0
        PHI=PARU(2)*PYR(0)
        CCHI=COS(CHI)
        SCHI=SIN(CHI)
        C2CHI=COS(2D0*CHI)
        S2CHI=SIN(2D0*CHI)
        THE=ACOS(CTHE)
        STHE=SIN(THE)
        C2PHI=COS(2D0*(PHI-PARJ(134)))
        S2PHI=SIN(2D0*(PHI-PARJ(134)))
        SIG=((1D0+CTHE**2)*HF1+STHE**2*C2PHI*HF3)*SIGU+2D0*(STHE**2*HF1-
     &  STHE**2*C2PHI*HF3)*SIGL+2D0*(STHE**2*C2CHI*HF1+((1D0+CTHE**2)*
     &  C2CHI*C2PHI-2D0*CTHE*S2CHI*S2PHI)*HF3)*SIGT-
     &  2D0*SQ2*(2D0*STHE*CTHE*CCHI*HF1-2D0*STHE*
     &  (CTHE*CCHI*C2PHI-SCHI*S2PHI)*HF3)*SIGI
        IF(SIG.LT.SIGMAX*PYR(0)) GOTO 120
        CALL PYROBO(NC+1,N,0D0,CHI,0D0,0D0,0D0)
        CALL PYROBO(NC+1,N,THE,PHI,0D0,0D0,0D0)
      ENDIF
 
C...Generate parton shower. Rearrange along strings and check.
      IF(MSTJ(101).GE.5.AND.NJET.GE.2) THEN
        CALL PYSHOW(NC+MK+1,-NJET,ECMC)
        MSTJ14=MSTJ(14)
        IF(MSTJ(105).EQ.-1) MSTJ(14)=-1
        IF(MSTJ(105).GE.0) MSTU(28)=0
        CALL PYPREP(0)
        MSTJ(14)=MSTJ14
        IF(MSTJ(105).GE.0.AND.MSTU(28).NE.0) GOTO 100
      ENDIF
 
C...Generate fragmentation. Information for PYTABU:
      IF(MSTJ(105).EQ.1) CALL PYEXEC
      MSTU(161)=110*KFLC+3
      MSTU(162)=0
 
      RETURN
      END
