
C*********************************************************************
 
C...PYEXEC
C...Administrates the fragmentation and decay chain.
 
      SUBROUTINE PYEXEC
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYINT1/,/PYINT4/
C...Local array.
      DIMENSION PS(2,6),IJOIN(100)
 
C...Initialize and reset.
      MSTU(24)=0
      IF(MSTU(12).NE.12345) CALL PYLIST(0)
      MSTU(29)=0
      MSTU(31)=MSTU(31)+1
      MSTU(1)=0
      MSTU(2)=0
      MSTU(3)=0
      IF(MSTU(17).LE.0) MSTU(90)=0
      MCONS=1
 
C...Sum up momentum, energy and charge for starting entries.
      NSAV=N
      DO 110 I=1,2
        DO 100 J=1,6
          PS(I,J)=0D0
  100   CONTINUE
  110 CONTINUE
      DO 130 I=1,N
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 130
        DO 120 J=1,4
          PS(1,J)=PS(1,J)+P(I,J)
  120   CONTINUE
        PS(1,6)=PS(1,6)+PYCHGE(K(I,2))
  130 CONTINUE
      PARU(21)=PS(1,4)
 
C...Start by all decays of coloured resonances involved in shower.
      NORIG=N
      DO 140 I=1,NORIG
        IF(K(I,1).EQ.3) THEN
          KC=PYCOMP(K(I,2))
          IF(MWID(KC).NE.0.AND.KCHG(KC,2).NE.0) CALL PYRESD(I)
        ENDIF
  140 CONTINUE
 
C...Prepare system for subsequent fragmentation/decay.
      CALL PYPREP(0)
      IF(MINT(51).NE.0) RETURN
 
C...Loop through jet fragmentation and particle decays.
      MBE=0
  150 MBE=MBE+1
      IP=0
  160 IP=IP+1
      KC=0
      IF(K(IP,1).GT.0.AND.K(IP,1).LE.10) KC=PYCOMP(K(IP,2))
      IF(KC.EQ.0) THEN
 
C...Deal with any remaining undecayed resonance
C...(normally the task of PYEVNT, so seldom used).
      ELSEIF(MWID(KC).NE.0) THEN
        IBEG=IP
        IF(KCHG(KC,2).NE.0.AND.K(I,1).NE.3) THEN
          IBEG=IP+1
  170     IBEG=IBEG-1
          IF(IBEG.GE.2.AND.K(IBEG,1).EQ.2) GOTO 170
          IF(K(IBEG,1).NE.2) IBEG=IBEG+1
          IEND=IP-1
  180     IEND=IEND+1
          IF(IEND.LT.N.AND.K(IEND,1).EQ.2) GOTO 180
          IF(IEND.LT.N.AND.KCHG(PYCOMP(K(IEND,2)),2).EQ.0) GOTO 180
          NJOIN=0
          DO 190 I=IBEG,IEND
            IF(KCHG(PYCOMP(K(IEND,2)),2).NE.0) THEN
              NJOIN=NJOIN+1
              IJOIN(NJOIN)=I
            ENDIF
  190     CONTINUE
        ENDIF
        CALL PYRESD(IP)
        CALL PYPREP(IBEG)
        IF(MINT(51).NE.0) RETURN
 
C...Particle decay if unstable and allowed. Save long-lived particle
C...decays until second pass after Bose-Einstein effects.
      ELSEIF(KCHG(KC,2).EQ.0) THEN
        IF(MSTJ(21).GE.1.AND.MDCY(KC,1).GE.1.AND.(MSTJ(51).LE.0.OR.MBE
     &  .EQ.2.OR.PMAS(KC,2).GE.PARJ(91).OR.IABS(K(IP,2)).EQ.311))
     &  CALL PYDECY(IP)
 
C...Decay products may develop a shower.
        IF(MSTJ(92).GT.0) THEN
          IP1=MSTJ(92)
          QMAX=SQRT(MAX(0D0,(P(IP1,4)+P(IP1+1,4))**2-(P(IP1,1)+P(IP1+1,
     &    1))**2-(P(IP1,2)+P(IP1+1,2))**2-(P(IP1,3)+P(IP1+1,3))**2))
          MINT(33)=0
          CALL PYSHOW(IP1,IP1+1,QMAX)
          CALL PYPREP(IP1)
          IF(MINT(51).NE.0) RETURN
          MSTJ(92)=0
        ELSEIF(MSTJ(92).LT.0) THEN
          IP1=-MSTJ(92)
          MINT(33)=0
          CALL PYSHOW(IP1,-3,P(IP,5))
          CALL PYPREP(IP1)
          IF(MINT(51).NE.0) RETURN
          MSTJ(92)=0
        ENDIF
 
C...Jet fragmentation: string or independent fragmentation.
      ELSEIF(K(IP,1).EQ.1.OR.K(IP,1).EQ.2) THEN
        MFRAG=MSTJ(1)
        IF(MFRAG.GE.1.AND.K(IP,1).EQ.1) MFRAG=2
        IF(MSTJ(21).GE.2.AND.K(IP,1).EQ.2.AND.N.GT.IP) THEN
          IF(K(IP+1,1).EQ.1.AND.K(IP+1,3).EQ.K(IP,3).AND.
     &    K(IP,3).GT.0.AND.K(IP,3).LT.IP) THEN
            IF(KCHG(PYCOMP(K(K(IP,3),2)),2).EQ.0) MFRAG=MIN(1,MFRAG)
          ENDIF
        ENDIF
        IF(MFRAG.EQ.1) CALL PYSTRF(IP)
        IF(MFRAG.EQ.2) CALL PYINDF(IP)
        IF(MFRAG.EQ.2.AND.K(IP,1).EQ.1) MCONS=0
        IF(MFRAG.EQ.2.AND.(MSTJ(3).LE.0.OR.MOD(MSTJ(3),5).EQ.0)) MCONS=0
      ENDIF
 
C...Loop back if enough space left in PYJETS and no error abort.
      IF(MSTU(24).NE.0.AND.MSTU(21).GE.2) THEN
      ELSEIF(IP.LT.N.AND.N.LT.MSTU(4)-20-MSTU(32)) THEN
        GOTO 160
      ELSEIF(IP.LT.N) THEN
        CALL PYERRM(11,'(PYEXEC:) no more memory left in PYJETS')
      ENDIF
 
C...Include simple Bose-Einstein effect parametrization if desired.
      IF(MBE.EQ.1.AND.MSTJ(51).GE.1) THEN
        CALL PYBOEI(NSAV)
        GOTO 150
      ENDIF
 
C...Check that momentum, energy and charge were conserved.
      DO 210 I=1,N
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 210
        DO 200 J=1,4
          PS(2,J)=PS(2,J)+P(I,J)
  200   CONTINUE
        PS(2,6)=PS(2,6)+PYCHGE(K(I,2))
  210 CONTINUE
      PDEV=(ABS(PS(2,1)-PS(1,1))+ABS(PS(2,2)-PS(1,2))+ABS(PS(2,3)-
     &PS(1,3))+ABS(PS(2,4)-PS(1,4)))/(1D0+ABS(PS(2,4))+ABS(PS(1,4)))
      IF(MCONS.EQ.1.AND.PDEV.GT.PARU(11)) CALL PYERRM(15,
     &'(PYEXEC:) four-momentum was not conserved')
      IF(MCONS.EQ.1.AND.ABS(PS(2,6)-PS(1,6)).GT.0.1D0) CALL PYERRM(15,
     &'(PYEXEC:) charge was not conserved')
 
      RETURN
      END
