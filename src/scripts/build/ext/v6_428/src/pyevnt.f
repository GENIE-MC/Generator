 
C*********************************************************************
 
C...PYEVNT
C...Administers the generation of a high-pT event via calls to
C...a number of subroutines.
 
      SUBROUTINE PYEVNT
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (MAXNUR=1000)
C...Commonblocks.
      COMMON/PYPART/NPART,NPARTD,IPART(MAXNUR),PTPART(MAXNUR)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYCTAG/NCT,MCT(4000,2)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      SAVE /PYJETS/,/PYDAT1/,/PYCTAG/,/PYDAT2/,/PYDAT3/,/PYPARS/,
     &/PYINT1/,/PYINT2/,/PYINT4/,/PYINT5/
C...Local array.
      DIMENSION VTX(4)
 
C...Optionally let PYEVNW do the whole job.
      IF(MSTP(81).GE.20) THEN
        CALL PYEVNW
        RETURN
      ENDIF
 
C...Stop if no subprocesses on.
      IF(MINT(121).EQ.1.AND.MSTI(53).EQ.1) THEN
        WRITE(MSTU(11),5100)
        CALL PYSTOP(1)
      ENDIF
 
C...Initial values for some counters.
      MSTU(1)=0
      MSTU(2)=0
      N=0
      MINT(5)=MINT(5)+1
      MINT(7)=0
      MINT(8)=0
      MINT(30)=0
      MINT(83)=0
      MINT(84)=MSTP(126)
      MSTU(24)=0
      MSTU70=0
      MSTJ14=MSTJ(14)
C...Normally, use K(I,4:5) colour info rather than /PYCTAG/.
      NCT=0
      MINT(33)=0
 
C...Let called routines know call is from PYEVNT (not PYEVNW).
      MINT(35)=1
      IF (MSTP(81).GE.10) MINT(35)=2
 
C...If variable energies: redo incoming kinematics and cross-section.
      MSTI(61)=0
      IF(MSTP(171).EQ.1) THEN
        CALL PYINKI(1)
        IF(MSTI(61).EQ.1) THEN
          MINT(5)=MINT(5)-1
          RETURN
        ENDIF
        IF(MINT(121).GT.1) CALL PYSAVE(3,1)
        CALL PYXTOT
      ENDIF
 
C...Loop over number of pileup events; check space left.
      IF(MSTP(131).LE.0) THEN
        NPILE=1
      ELSE
        CALL PYPILE(2)
        NPILE=MINT(81)
      ENDIF
      DO 270 IPILE=1,NPILE
        IF(MINT(84)+100.GE.MSTU(4)) THEN
          CALL PYERRM(11,
     &    '(PYEVNT:) no more space in PYJETS for pileup events')
          IF(MSTU(21).GE.1) GOTO 280
        ENDIF
        MINT(82)=IPILE
 
C...Generate variables of hard scattering.
        MINT(51)=0
        MSTI(52)=0
  100   CONTINUE
        IF(MINT(51).NE.0.OR.MSTU(24).NE.0) MSTI(52)=MSTI(52)+1
        MINT(31)=0
        MINT(39)=0
        MINT(51)=0
        MINT(57)=0
        CALL PYRAND
        IF(MSTI(61).EQ.1) THEN
          MINT(5)=MINT(5)-1
          RETURN
        ENDIF
        IF(MINT(51).EQ.2) RETURN
        ISUB=MINT(1)
        IF(MSTP(111).EQ.-1) GOTO 260
 
C...Loopback point if PYPREP fails, especially for junction topologies.
        NPREP=0
        MNT31S=MINT(31)
  110   NPREP=NPREP+1
        MINT(31)=MNT31S
 
        IF((ISUB.LE.90.OR.ISUB.GE.95).AND.ISUB.NE.99) THEN
C...Hard scattering (including low-pT):
C...reconstruct kinematics and colour flow of hard scattering.
          MINT31=MINT(31)
  120     MINT(31)=MINT31
          MINT(51)=0
          CALL PYSCAT
          IF(MINT(51).EQ.1) GOTO 100
          IPU1=MINT(84)+1
          IPU2=MINT(84)+2
          IF(ISUB.EQ.95) GOTO 140
 
C...Reset statistics on activity in event.
        DO 130 J=351,359
          MINT(J)=0
          VINT(J)=0D0
  130   CONTINUE
 
C...Showering of initial state partons (optional).
          NFIN=N
          ALAMSV=PARJ(81)
          PARJ(81)=PARP(72)
          IF(MSTP(61).GE.1.AND.MINT(47).GE.2.AND.MINT(111).NE.12)
     &    CALL PYSSPA(IPU1,IPU2)
          PARJ(81)=ALAMSV
          IF(MINT(51).EQ.1) GOTO 100

C...pT-ordered FSR off ISR (optional, must have at least 2 partons)
          IF (NPART.GE.2.AND.(MSTJ(41).EQ.11.OR.MSTJ(41).EQ.12)) THEN
            PTMAX=0.5*SQRT(PARP(71))*VINT(55)
            CALL PYPTFS(3,PTMAX,0D0,PTGEN)
          ENDIF
 
C...Showering of final state partons (optional).
          ALAMSV=PARJ(81)
          PARJ(81)=PARP(72)
          IF(MSTP(71).GE.1.AND.ISET(ISUB).GE.2.AND.ISET(ISUB).LE.10)
     &    THEN
            IPU3=MINT(84)+3
            IPU4=MINT(84)+4
            IF(ISET(ISUB).EQ.5) IPU4=-3
            QMAX=VINT(55)
            IF(ISET(ISUB).EQ.2) QMAX=SQRT(PARP(71))*VINT(55)
            CALL PYSHOW(IPU3,IPU4,QMAX)
          ELSEIF(ISET(ISUB).EQ.11) THEN
            CALL PYADSH(NFIN)
          ENDIF
          PARJ(81)=ALAMSV
 
C...Allow possibility for user to abort event generation.
          IVETO=0
          IF(IPILE.EQ.1.AND.MSTP(143).EQ.1) CALL PYVETO(IVETO)
          IF(IVETO.EQ.1) GOTO 100
 
C...Decay of final state resonances.
          MINT(32)=0
          IF(MSTP(41).GE.1.AND.ISET(ISUB).LE.10) CALL PYRESD(0)
          IF(MINT(51).EQ.1) GOTO 100
          MINT(52)=N
 
 
C...Multiple interactions - PYTHIA 6.3 intermediate style.
  140     IF(MSTP(81).GE.10.AND.MINT(50).EQ.1) THEN
            IF(ISUB.EQ.95) MINT(31)=MINT(31)+1
            CALL PYMIGN(6)
            IF(MINT(51).EQ.1) GOTO 100
            MINT(53)=N
 
C...Beam remnant flavour and colour assignments - new scheme.
            CALL PYMIHK
            IF(MINT(51).EQ.1.AND.MINT(57).GE.1.AND.MINT(57).LE.5)
     &      GOTO 120
            IF(MINT(51).EQ.1) GOTO 100
 
C...Primordial kT and beam remnant momentum sharing - new scheme.
            CALL PYMIRM
            IF(MINT(51).EQ.1.AND.MINT(57).GE.1.AND.MINT(57).LE.5)
     &      GOTO 120
            IF(MINT(51).EQ.1) GOTO 100
            IF(ISUB.EQ.95) MINT(31)=MINT(31)-1
 
C...Multiple interactions - PYTHIA 6.2 style.
          ELSEIF(MINT(111).NE.12) THEN
            IF (MSTP(81).GE.1.AND.MINT(50).EQ.1.AND.ISUB.NE.95) THEN
              CALL PYMULT(6)
              MINT(53)=N
            ENDIF
 
C...Hadron remnants and primordial kT.
            CALL PYREMN(IPU1,IPU2)
            IF(MINT(51).EQ.1.AND.MINT(57).GE.1.AND.MINT(57).LE.5) GOTO
     &           110
            IF(MINT(51).EQ.1) GOTO 100
          ENDIF
 
        ELSEIF(ISUB.NE.99) THEN
C...Diffractive and elastic scattering.
          CALL PYDIFF
 
        ELSE
C...DIS scattering (photon flux external).
          CALL PYDISG
          IF(MINT(51).EQ.1) GOTO 100
        ENDIF
 
C...Check that no odd resonance left undecayed.
        MINT(54)=N
        IF(MSTP(111).GE.1) THEN
          NFIX=N
          DO 150 I=MINT(84)+1,NFIX
            IF(K(I,1).GE.1.AND.K(I,1).LE.10.AND.K(I,2).NE.21.AND.
     &      K(I,2).NE.22) THEN
              KCA=PYCOMP(K(I,2))
              IF(MWID(KCA).NE.0.AND.MDCY(KCA,1).GE.1) THEN
                CALL PYRESD(I)
                IF(MINT(51).EQ.1) GOTO 100
              ENDIF
            ENDIF
  150     CONTINUE
        ENDIF
 
C...Boost hadronic subsystem to overall rest frame.
C..(Only relevant when photon inside lepton beam.)
        IF(MINT(141).NE.0.OR.MINT(142).NE.0) CALL PYGAGA(4,WTGAGA)
 
C...Recalculate energies from momenta and masses (if desired).
        IF(MSTP(113).GE.1) THEN
          DO 160 I=MINT(83)+1,N
            IF(K(I,1).GT.0.AND.K(I,1).LE.10) P(I,4)=SQRT(P(I,1)**2+
     &      P(I,2)**2+P(I,3)**2+P(I,5)**2)
  160     CONTINUE
          NRECAL=N
        ENDIF
 
C...Colour reconnection before string formation
        IF (MSTP(95).GE.2) CALL PYFSCR(MINT(84)+1)

C...Rearrange partons along strings, check invariant mass cuts.
        MSTU(28)=0
        IF(MSTP(111).LE.0) MSTJ(14)=-1
        CALL PYPREP(MINT(84)+1)
        MSTJ(14)=MSTJ14
        IF(MINT(51).EQ.1.AND.MSTU(24).EQ.1) THEN
          MSTU(24)=0
          GOTO 100
        ENDIF
        IF (MINT(51).EQ.1.AND.NPREP.LE.5) GOTO 110
        IF (MINT(51).EQ.1) GOTO 100
        IF(MSTP(112).EQ.1.AND.MSTU(28).EQ.3) GOTO 100
        IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) THEN
          DO 190 I=MINT(84)+1,N
            IF(K(I,2).EQ.94) THEN
              DO 180 I1=I+1,MIN(N,I+10)
                IF(K(I1,3).EQ.I) THEN
                  K(I1,3)=MOD(K(I1,4)/MSTU(5),MSTU(5))
                  IF(K(I1,3).EQ.0) THEN
                    DO 170 II=MINT(84)+1,I-1
                        IF(K(II,2).EQ.K(I1,2)) THEN
                          IF(MOD(K(II,4),MSTU(5)).EQ.I1.OR.
     &                    MOD(K(II,5),MSTU(5)).EQ.I1) K(I1,3)=II
                        ENDIF
  170               CONTINUE
                    IF(K(I+1,3).EQ.0) K(I+1,3)=K(I,3)
                  ENDIF
                ENDIF
  180         CONTINUE
            ENDIF
  190     CONTINUE
          CALL PYEDIT(12)
          CALL PYEDIT(14)
          IF(MSTP(125).EQ.0) CALL PYEDIT(15)
          IF(MSTP(125).EQ.0) MINT(4)=0
          DO 210 I=MINT(83)+1,N
            IF(K(I,1).EQ.11.AND.K(I,4).EQ.0.AND.K(I,5).EQ.0) THEN
              DO 200 I1=I+1,N
                IF(K(I1,3).EQ.I.AND.K(I,4).EQ.0) K(I,4)=I1
                IF(K(I1,3).EQ.I) K(I,5)=I1
  200         CONTINUE
            ENDIF
  210     CONTINUE
        ENDIF
 
C...Introduce separators between sections in PYLIST event listing.
        IF(IPILE.EQ.1.AND.MSTP(125).LE.0) THEN
          MSTU70=1
          MSTU(71)=N
        ELSEIF(IPILE.EQ.1) THEN
          MSTU70=3
          MSTU(71)=2
          MSTU(72)=MINT(4)
          MSTU(73)=N
        ENDIF
 
C...Go back to lab frame (needed for vertices, also in fragmentation).
        CALL PYFRAM(1)
 
C...Set nonvanishing production vertex (optional).
        IF(MSTP(151).EQ.1) THEN
          DO 220 J=1,4
            VTX(J)=PARP(150+J)*SQRT(-2D0*LOG(MAX(1D-10,PYR(0))))*
     &      SIN(PARU(2)*PYR(0))
  220     CONTINUE
          DO 240 I=MINT(83)+1,N
            DO 230 J=1,4
              V(I,J)=V(I,J)+VTX(J)
  230       CONTINUE
  240     CONTINUE
        ENDIF
 
C...Perform hadronization (if desired).
        IF(MSTP(111).GE.1) THEN
          CALL PYEXEC
          IF(MSTU(24).NE.0) GOTO 100
        ENDIF
        IF(MSTP(113).GE.1) THEN
          DO 250 I=NRECAL,N
            IF(P(I,5).GT.0D0) P(I,4)=SQRT(P(I,1)**2+
     &      P(I,2)**2+P(I,3)**2+P(I,5)**2)
  250     CONTINUE
        ENDIF
        IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) CALL PYEDIT(14)
 
C...Store event information and calculate Monte Carlo estimates of
C...subprocess cross-sections.
  260   IF(IPILE.EQ.1) CALL PYDOCU
 
C...Set counters for current pileup event and loop to next one.
        MSTI(41)=IPILE
        IF(IPILE.GE.2.AND.IPILE.LE.10) MSTI(40+IPILE)=ISUB
        IF(MSTU70.LT.10) THEN
          MSTU70=MSTU70+1
          MSTU(70+MSTU70)=N
        ENDIF
        MINT(83)=N
        MINT(84)=N+MSTP(126)
        IF(IPILE.LT.NPILE) CALL PYFRAM(2)
  270 CONTINUE
 
C...Generic information on pileup events. Reconstruct missing history.
      IF(MSTP(131).EQ.1.AND.MSTP(133).GE.1) THEN
        PARI(91)=VINT(132)
        PARI(92)=VINT(133)
        PARI(93)=VINT(134)
        IF(MSTP(133).GE.2) PARI(93)=PARI(93)*XSEC(0,3)/VINT(131)
      ENDIF
      CALL PYEDIT(16)
 
C...Transform to the desired coordinate frame.
  280 CALL PYFRAM(MSTP(124))
      MSTU(70)=MSTU70
      PARU(21)=VINT(1)
 
C...Error messages
 5100 FORMAT(1X,'Error: no subprocess switched on.'/
     &1X,'Execution stopped.')
 
      RETURN
      END
