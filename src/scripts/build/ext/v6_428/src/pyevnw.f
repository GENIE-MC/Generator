 
C*********************************************************************
 
C...PYEVNW
C...Administers the generation of a high-pT event via calls to
C...a number of subroutines for the new multiple interactions and
C...showering framework.
 
      SUBROUTINE PYEVNW
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (MAXNUR=1000)
C...Commonblocks.
      COMMON/PYPART/NPART,NPARTD,IPART(MAXNUR),PTPART(MAXNUR)
C...Commonblocks.
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
      COMMON/PYINTM/KFIVAL(2,3),NMI(2),IMI(2,800,2),NVC(2,-6:6),
     &     XASSOC(2,-6:6,240),XPSVC(-6:6,-1:240),PVCTOT(2,-1:1),
     &     XMI(2,240),PT2MI(240),IMISEP(0:240)
      SAVE /PYJETS/,/PYCTAG/,/PYDAT1/,/PYDAT2/,/PYDAT3/,
     &     /PYPARS/,/PYINT1/,/PYINT2/,/PYINT4/,/PYINT5/,/PYINTM/
C...Local arrays.
      DIMENSION VTX(4)
 
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
C...Normally, use K(I,4:5) colour info rather than /PYCT/.
      NCT=0
      MINT(33)=0
C...Zero counters for pT-ordered showers (failsafe)
      NPART=0
      NPARTD=0
 
C...Let called routines know call is from PYEVNW (not PYEVNT).
      MINT(35)=3
 
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
      DO 300 IPILE=1,NPILE
        IF(MINT(84)+100.GE.MSTU(4)) THEN
          CALL PYERRM(11,
     &    '(PYEVNW:) no more space in PYJETS for pileup events')
          IF(MSTU(21).GE.1) GOTO 310
        ENDIF
        MINT(82)=IPILE
 
C...Generate variables of hard scattering.
        MINT(51)=0
        MSTI(52)=0
        LOOPHS  =0
  100   CONTINUE
        LOOPHS  = LOOPHS + 1
        IF(MINT(51).NE.0.OR.MSTU(24).NE.0) MSTI(52)=MSTI(52)+1
        IF(LOOPHS.GE.10) THEN
          CALL PYERRM(19,'(PYEVNW:) failed to evolve shower or '
     &        //'multiple interactions. Returning.')
          MINT(51)=1
          RETURN
        ENDIF
        MINT(31)=0
        MINT(39)=0
        MINT(36)=0
        MINT(51)=0
        MINT(57)=0
        CALL PYRAND
        IF(MSTI(61).EQ.1) THEN
          MINT(5)=MINT(5)-1
          RETURN
        ENDIF
        IF(MINT(51).EQ.2) RETURN
        ISUB=MINT(1)
        IF(MSTP(111).EQ.-1) GOTO 290
 
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
          NPARTD=N
          NFIN=N
 
C...Intertwined initial state showers and multiple interactions.
C...Force no IS showers if no pdfs defined: MSTP(61) -> 0 for PYEVOL.
C...Force no MI if cross section not known: MSTP(81) -> 0 for PYEVOL.
          MSTP61=MSTP(61)
          IF (MINT(47).LT.2) MSTP(61)=0
          MSTP81=MSTP(81)
          IF (MINT(50).EQ.0) MSTP(81)=0
          IF ((MSTP(61).GE.1.OR.MOD(MSTP(81),10).GE.0).AND.
     &    MINT(111).NE.12) THEN
C...Absolute max pT2 scale for evolution: phase space limit.
            PT2MXS=0.25D0*VINT(2)
C...Check if more constrained by ISR and MI max scales:
            PT2MXS=MIN(PT2MXS,MAX(MAX(1D0,PARP(67))*VINT(56),VINT(62)))
C...Loopback point in case of failure in evolution.
            LOOP=0
  130       LOOP=LOOP+1
            MINT(51)=0
            IF(LOOP.GT.100) THEN
              CALL PYERRM(9,'(PYEVNW:) failed to evolve shower or '
     &             //'multiple interactions. Trying new point.')
              MINT(51)=1
              RETURN
            ENDIF
 
C...Pre-initialization of interleaved MI/ISR/JI evolution, only done
C...once per event. (E.g. compute constants and save variables to be
C...restored later in case of failure.)
            IF (LOOP.EQ.1) CALL PYEVOL(-1,DUMMY1,DUMMY2)
 
C...Initialize interleaved MI/ISR/JI evolution.
C...PT2MAX: absolute upper limit for evolution - Initialization may
C...        return a PT2MAX which is lower than this.
C...PT2MIN: absolute lower limit for evolution - Initialization may
C...        return a PT2MIN which is larger than this (e.g. Lambda_QCD).
            PT2MAX=PT2MXS
            PT2MIN=0D0
            CALL PYEVOL(0,PT2MAX,PT2MIN)
C...If failed to initialize evolution, generate a new hard process
            IF (MINT(51).EQ.1) GOTO 100
 
C...Perform interleaved MI/ISR/JI evolution from PT2MAX to PT2MIN.
C...In principle factorized, so can be stopped and restarted.
C...Example: stop/start at pT=10 GeV. (Commented out for now.)
C            PT2MED=MAX(10D0**2,PT2MIN)
C            CALL PYEVOL(1,PT2MAX,PT2MED)
C            IF (MINT(51).EQ.1) GOTO 160
C            PT2MAX=PT2MED
            CALL PYEVOL(1,PT2MAX,PT2MIN)
C...If fatal error (e.g., massive hard-process initiator, but no available 
C...phase space for creation), generate a new hard process
            IF (MINT(51).EQ.2) GOTO 100
C...If smaller error, just try running evolution again
            IF (MINT(51).EQ.1) GOTO 130
 
C...Finalize interleaved MI/ISR/JI evolution.
            CALL PYEVOL(2,PT2MAX,PT2MIN)
            IF (MINT(51).EQ.1) GOTO 130
 
          ENDIF
          MSTP(61)=MSTP61
          MSTP(81)=MSTP81
          IF(MINT(51).EQ.1) GOTO 100
C...(MINT(52) is actually obsolete in this routine. Set anyway
C...to ensure PYDOCU stable.)
          MINT(52)=N
          MINT(53)=N
 
C...Beam remnants - new scheme.
  140     IF(MINT(50).EQ.1) THEN
            IF (ISUB.EQ.95) MINT(31)=1
 
C...Beam remnant flavour and colour assignments - new scheme.
            CALL PYMIHK
            IF(MINT(51).EQ.1.AND.MINT(57).GE.1.AND.MINT(57).LE.5)
     &           GOTO 120
            IF(MINT(51).EQ.1) GOTO 100
 
C...Primordial kT and beam remnant momentum sharing - new scheme.
            CALL PYMIRM
            IF(MINT(51).EQ.1.AND.MINT(57).GE.1.AND.MINT(57).LE.5)
     &      GOTO 120
            IF(MINT(51).EQ.1) GOTO 100
            IF (ISUB.EQ.95) MINT(31)=0
          ELSEIF(MINT(111).NE.12) THEN
C...Hadron remnants and primordial kT - old model.
C...Happens e.g. for direct photon on one side.
            IPU1=IMI(1,1,1)
            IPU2=IMI(2,1,1)
            CALL PYREMN(IPU1,IPU2)
            IF(MINT(51).EQ.1.AND.MINT(57).GE.1.AND.MINT(57).LE.5) GOTO
     &           110
            IF(MINT(51).EQ.1) GOTO 100
C...PYREMN does not set colour tags for BRs, so needs to be done now.
            DO 160 I=MINT(53)+1,N
              DO 150 KCS=4,5
                IDA=MOD(K(I,KCS),MSTU(5))
                IF (IDA.NE.0) THEN
                  MCT(I,KCS-3)=MCT(IDA,6-KCS)
                ELSE
                  MCT(I,KCS-3)=0
                ENDIF
  150         CONTINUE
  160       CONTINUE
C...Instruct PYPREP to use colour tags
            MINT(33)=1

            DO 360 MQGST=1,2
              DO 350 I=MINT(84)+1,N
  
C...Look for coloured string endpoint, or (later) leftover gluon.
                IF (K(I,1).NE.3) GOTO 350
                KC=PYCOMP(K(I,2))
                IF(KC.EQ.0) GOTO 350
                KQ=KCHG(KC,2)
                IF(KQ.EQ.0.OR.(MQGST.EQ.1.AND.KQ.EQ.2)) GOTO 350
  
C...  Pick up loose string end with no previous tag.
                KCS=4
                IF(KQ*ISIGN(1,K(I,2)).LT.0) KCS=5
                IF(MCT(I,KCS-3).NE.0) GOTO 350
                  
                CALL PYCTTR(I,KCS,I)
                IF(MINT(51).NE.0) RETURN
  
 350          CONTINUE
 360        CONTINUE
C...Now delete any colour processing information if set (since partons
C...otherwise not FS showered!)
            DO 170 I=MINT(84)+1,N
              IF (I.LE.N) THEN
                K(I,4)=MOD(K(I,4),MSTU(5)**2)
                K(I,5)=MOD(K(I,5),MSTU(5)**2)
              ENDIF
  170       CONTINUE
          ENDIF
 
C...Showering of final state partons (optional).
          ALAMSV=PARJ(81)
          PARJ(81)=PARP(72)
          IF(MSTP(71).GE.1.AND.ISET(ISUB).GE.1.AND.ISET(ISUB).LE.10)
     &    THEN
            QMAX=VINT(55)
            IF(ISET(ISUB).EQ.2) QMAX=SQRT(PARP(71))*VINT(55)
            CALL PYPTFS(1,QMAX,0D0,PTGEN)
C...External processes: handle successive showers.
          ELSEIF(ISET(ISUB).EQ.11) THEN
            CALL PYADSH(NFIN)
          ENDIF
          PARJ(81)=ALAMSV

C...Allow possibility for user to abort event generation.
          IVETO=0
          IF(IPILE.EQ.1.AND.MSTP(143).EQ.1) CALL PYVETO(IVETO) ! sm
          IF(IVETO.EQ.1) THEN
C...........No reason to count this as an error
            LOOPHS = LOOPHS-1
            GOTO 100
          ENDIF

 
C...Decay of final state resonances.
          MINT(32)=0
          IF(MSTP(41).GE.1.AND.ISET(ISUB).LE.10) THEN
            CALL PYRESD(0)
            IF(MINT(51).NE.0) GOTO 100
          ENDIF
 
          IF(MINT(51).EQ.1) GOTO 100
 
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
          DO 180 I=MINT(84)+1,NFIX
            IF(K(I,1).GE.1.AND.K(I,1).LE.10.AND.K(I,2).NE.21.AND.
     &      K(I,2).NE.22) THEN
              KCA=PYCOMP(K(I,2))
              IF(MWID(KCA).NE.0.AND.MDCY(KCA,1).GE.1) THEN
                CALL PYRESD(I)
                IF(MINT(51).EQ.1) GOTO 100
              ENDIF
            ENDIF
  180     CONTINUE
        ENDIF
 
C...Boost hadronic subsystem to overall rest frame.
C..(Only relevant when photon inside lepton beam.)
        IF(MINT(141).NE.0.OR.MINT(142).NE.0) CALL PYGAGA(4,WTGAGA)
 
C...Recalculate energies from momenta and masses (if desired).
        IF(MSTP(113).GE.1) THEN
          DO 190 I=MINT(83)+1,N
            IF(K(I,1).GT.0.AND.K(I,1).LE.10) P(I,4)=SQRT(P(I,1)**2+
     &      P(I,2)**2+P(I,3)**2+P(I,5)**2)
  190     CONTINUE
          NRECAL=N
        ENDIF
 
C...Colour reconnection before string formation
        CALL PYFSCR(MINT(84)+1)
 
C...Rearrange partons along strings, check invariant mass cuts.
        MSTU(28)=0
        IF(MSTP(111).LE.0) MSTJ(14)=-1
        CALL PYPREP(MINT(84)+1)
        MSTJ(14)=MSTJ14
        IF(MINT(51).EQ.1.AND.MSTU(24).EQ.1) THEN
          MSTU(24)=0
          GOTO 100
        ENDIF
        IF(MINT(51).EQ.1) GOTO 110
        IF(MSTP(112).EQ.1.AND.MSTU(28).EQ.3) GOTO 100
        IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) THEN
          DO 220 I=MINT(84)+1,N
            IF(K(I,2).EQ.94) THEN
              DO 210 I1=I+1,MIN(N,I+10)
                IF(K(I1,3).EQ.I) THEN
                  K(I1,3)=MOD(K(I1,4)/MSTU(5),MSTU(5))
                  IF(K(I1,3).EQ.0) THEN
                    DO 200 II=MINT(84)+1,I-1
                        IF(K(II,2).EQ.K(I1,2)) THEN
                          IF(MOD(K(II,4),MSTU(5)).EQ.I1.OR.
     &                    MOD(K(II,5),MSTU(5)).EQ.I1) K(I1,3)=II
                        ENDIF
  200               CONTINUE
                    IF(K(I+1,3).EQ.0) K(I+1,3)=K(I,3)
                  ENDIF
                ENDIF
  210         CONTINUE
C...Also collapse particles decaying to themselves (if same KS)
C...Sep 22 2009: Commented out by PS following suggestion by TS to fix 
C...problem with history point-backs in new shower, where a particle is
C...copied with a new momentum when it is the recoiler.
C            ELSEIF (K(I,1).GT.0.AND.K(I,4).EQ.K(I,5).AND.K(I,4).GT.0
C     &            .AND.K(I,4).LT.N) THEN
C              IDA=K(I,4)
C              IF (K(IDA,1).EQ.K(I,1).AND.K(IDA,2).EQ.K(I,2)) THEN
C                K(I,1)=0
C              ENDIF
            ENDIF
  220     CONTINUE
          CALL PYEDIT(12)
          CALL PYEDIT(14)
          IF(MSTP(125).EQ.0) CALL PYEDIT(15)
          IF(MSTP(125).EQ.0) MINT(4)=0
          DO 240 I=MINT(83)+1,N
            IF(K(I,1).EQ.11.AND.K(I,4).EQ.0.AND.K(I,5).EQ.0) THEN
              DO 230 I1=I+1,N
                IF(K(I1,3).EQ.I.AND.K(I,4).EQ.0) K(I,4)=I1
                IF(K(I1,3).EQ.I) K(I,5)=I1
  230         CONTINUE
            ENDIF
  240     CONTINUE
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
          DO 250 J=1,4
            VTX(J)=PARP(150+J)*SQRT(-2D0*LOG(MAX(1D-10,PYR(0))))*
     &      SIN(PARU(2)*PYR(0))
  250     CONTINUE
          DO 270 I=MINT(83)+1,N
            DO 260 J=1,4
              V(I,J)=V(I,J)+VTX(J)
  260       CONTINUE
  270     CONTINUE
        ENDIF
 
C...Perform hadronization (if desired).
        IF(MSTP(111).GE.1) THEN
          CALL PYEXEC
          IF(MSTU(24).NE.0) GOTO 100
        ENDIF
        IF(MSTP(113).GE.1) THEN
          DO 280 I=NRECAL,N
            IF(P(I,5).GT.0D0) P(I,4)=SQRT(P(I,1)**2+
     &      P(I,2)**2+P(I,3)**2+P(I,5)**2)
  280     CONTINUE
        ENDIF
        IF(MSTP(125).EQ.0.OR.MSTP(125).EQ.1) CALL PYEDIT(14)
 
C...Store event information and calculate Monte Carlo estimates of
C...subprocess cross-sections.
  290   IF(IPILE.EQ.1) CALL PYDOCU
 
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
  300 CONTINUE
 
C...Generic information on pileup events. Reconstruct missing history.
      IF(MSTP(131).EQ.1.AND.MSTP(133).GE.1) THEN
        PARI(91)=VINT(132)
        PARI(92)=VINT(133)
        PARI(93)=VINT(134)
        IF(MSTP(133).GE.2) PARI(93)=PARI(93)*XSEC(0,3)/VINT(131)
      ENDIF
      CALL PYEDIT(16)
 
C...Transform to the desired coordinate frame.
  310 CALL PYFRAM(MSTP(124))
      MSTU(70)=MSTU70
      PARU(21)=VINT(1)
 
C...Error messages
 5100 FORMAT(1X,'Error: no subprocess switched on.'/
     &1X,'Execution stopped.')
 
      RETURN
      END
