 
C*********************************************************************
 
C...PYFSCR
C...Performs colour annealing.
C...MSTP(95) : CR Type
C...         = 1  : old cut-and-paste reconnections, handled in PYMIHK
C...         = 2  : Type I(no gg loops); hadron-hadron only
C...         = 3  : Type I(no gg loops); all beams
C...         = 4  : Type II(gg loops)  ; hadron-hadron only
C...         = 5  : Type II(gg loops)  ; all beams
C...         = 6  : Type S             ; hadron-hadron only
C...         = 7  : Type S             ; all beams
C...         = 8  : Type P             ; hadron-hadron only
C...         = 9  : Type P             ; all beams
C...Types I and II are described in Sandhoff+Skands, in hep-ph/0604120.
C...Type S is driven by starting only from free triplets, not octets.
C...Type P is also driven by free triplets, but the reconnect probability
C...is computed from the string density per unit rapidity, where the axis
C...with respect to which the rapidity is computed is the Thrust axis of the
C...event.
C...A string piece remains unchanged with probability
C...    PKEEP = (1-PARP(78))**N
C...This scaling corresponds to each string piece having to go through
C...N other ones, each with probability PARP(78) for reconnection.
C...For types I, II, and S, N is chosen simply as the number of multiple
C...interactions, for a rough scaling with the general level of activity.
C...For type P, N is chosen to be the number of string pieces in a given
C...interval of rapidity (minus one, since the string doesn't reconnect
C...with itself), and the reconnect probability is interpreted as the
C...probability per unit rapidity.
C...It also also possible to apply a dampening factor to the CR strength,
C...using PARP(77), which will cause reconnections among high-pT string
C...pieces to be suppressed.
 
      SUBROUTINE PYFSCR(IP)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYINT1/MINT(400),VINT(400)
C...The common block of colour tags.
      COMMON/PYCTAG/NCT,MCT(4000,2)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYINT1/,/PYCTAG/,
     &/PYPARS/
C...MCN: Temporary storage of new colour tags
      INTEGER MCN(4000,2)
C...Arrays for storing color strings
      PARAMETER (NBINY=100)
      INTEGER ICR(4000),MSCR(4000)
      INTEGER IOPT(4000), NSTRY(NBINY)
      DOUBLE PRECISION RLOPTC(4000)
 
C...Function to give four-product.
      FOUR(I,J)=P(I,4)*P(J,4)
     &          -P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)
 
C...Check valid range of MSTP(95), local copy
      IF (MSTP(95).LE.1.OR.MSTP(95).GE.10) RETURN
      MSTP95=MOD(MSTP(95),10)
C...Set whether CR allowed inside resonance systems or not
C...(not implemented yet)
C      MRESCR=1
C      IF (MSTP(95).GE.10) MRESCR=0
 
C...Check whether colour tags already defined
      IF (MINT(33).EQ.0) THEN
C...Erase any existing colour tags for this event
        DO 100 I=1,N
          MCT(I,1)=0
          MCT(I,2)=0
  100   CONTINUE
C...Create colour tags for this event
        DO 120 I=1,N
          IF (K(I,1).EQ.3) THEN
            DO 110 KCS=4,5
              KCSIN=KCS
              IF (MCT(I,KCSIN-3).EQ.0) THEN
                CALL PYCTTR(I,KCSIN,I)
              ENDIF
  110       CONTINUE
          ENDIF
  120   CONTINUE
C...Instruct PYPREP to use colour tags
        MINT(33)=1
      ENDIF
 
C...For MSTP(95) even, only apply to hadron-hadron
      KA1=IABS(MINT(11))
      KA2=IABS(MINT(12))
      IF (MOD(MSTP(95),2).EQ.0.AND.(KA1.LT.100.OR.KA2.LT.100)) GOTO 9999
 
C...Initialize new tag array (but do not delete old yet)
      LCT=NCT
      DO 130 I=MAX(1,IP),N
         MCN(I,1)=0
         MCN(I,2)=0
  130 CONTINUE
 
C...For Paquis type, determine thrust axis (default along Z axis)
      TX=0D0
      TY=0D0
      TZ=1D0
      IF (MSTP95.GE.8) THEN
        CALL PYTHRU(THRDUM,OBLDUM)
        TX = P(N+1,1)
        TY = P(N+1,2)
        TZ = P(N+1,3)
      ENDIF
 
C...For each final-state dipole, check whether string should be
C...preserved.
      NCR=0
      IA=0
      IC=0
      RAPMAX=0.0
 
      ICTMIN=NCT
      DO 150 ICT=1,NCT
        IA=0
        IC=0
        DO 140 I=MAX(1,IP),N
          IF (K(I,1).EQ.3.AND.MCT(I,1).EQ.ICT) IC=I
          IF (K(I,1).EQ.3.AND.MCT(I,2).EQ.ICT) IA=I
  140   CONTINUE
        IF (IC.NE.0.AND.IA.NE.0) THEN
C...Save smallest NCT value so far
          ICTMIN = MIN(ICTMIN,ICT)
C...For Paquis algorithm, just store all string pieces for now
          IF (MSTP95.GE.8) THEN
C...  Add coloured parton
            NCR=NCR+1
            ICR(NCR)=IC
            MSCR(NCR)=1
            IOPT(NCR)=0
C...  Store rapidity (along Thrust axis) in RLOPT for the time being
C...  Add pion mass headroom to energy for this calculation
            EET = P(IC,4)*SQRT(1D0+(0.135D0/P(IC,4))**2)
            PZT = P(IC,1)*TX+P(IC,2)*TY+P(IC,3)*TZ
            RLOPTC(NCR)=LOG((EET+PZT)/(EET-PZT))
C...  Add anti-coloured parton
            NCR       = NCR+1
            ICR(NCR)  = IA
            MSCR(NCR) = 2
            IOPT(NCR) = 0
C...  Store rapidity (along Thrust axis) in RLOPT for the time being
            EET = P(IA,4)*SQRT(1D0+(0.135D0/P(IA,4))**2)
            PZT = P(IA,1)*TX+P(IA,2)*TY+P(IA,3)*TZ
            RLOPTC(NCR)=LOG((EET+PZT)/(EET-PZT))
C...  Keep track of largest endpoint "rapidity"
            RAPMAX = MAX(RAPMAX,ABS(RLOPTC(NCR)))
            RAPMAX = MAX(RAPMAX,ABS(RLOPTC(NCR-1)))
          ELSE
            CRMODF=1D0
C...  Opt: suppress breakup of high-boost string pieces (i.e., let them escape)
C...  (so far ignores the possibility that the whole "muck" may be moving.)
            IF (PARP(77).GT.0D0) THEN
              PT2STR=(P(IA,1)+P(IC,1))**2+(P(IA,2)+P(IC,2))**2
C...  For lepton-lepton, use actual p2/m2, otherwise approximate p2 ~ 3/2 pT2
              IF (KA1.LT.100.AND.KA2.LT.100) THEN
                P2STR = PT2STR + (P(IA,3)+P(IC,3))**2
              ELSE
                P2STR = 3D0/2D0 * PT2STR
              ENDIF
              RM2STR=(P(IA,4)+P(IC,4))**2-(P(IA,3)+P(IC,3))**2-PT2STR
              RM2STR=MAX(RM2STR,PMAS(PYCOMP(111),1)**2)
C...  Estimate number of particles ~ log(M2), cut off at 1.
              RLOGM2=MAX(1D0,LOG(RM2STR))
              P2AVG=P2STR/RLOGM2
C...  Supress reconnection probability by 1/(1+P77*P2AVG)
              CRMODF=1D0/(1D0+PARP(77)**2*P2AVG)
            ENDIF
            PKEEP=(1D0-PARP(78)*CRMODF)**MINT(31)
            IF (PYR(0).LE.PKEEP) THEN
              LCT=LCT+1
              MCN(IC,1)=LCT
              MCN(IA,2)=LCT
            ELSE
C...  Add coloured parton
              NCR=NCR+1
              ICR(NCR)=IC
              MSCR(NCR)=1
              IOPT(NCR)=0
              RLOPTC(NCR)=1D19
C...  Add anti-coloured parton
              NCR=NCR+1
              ICR(NCR)=IA
              MSCR(NCR)=2
              IOPT(NCR)=0
              RLOPTC(NCR)=1D19
            ENDIF
          ENDIF
        ENDIF
  150 CONTINUE
 
C...PAQUIS TYPE
      IF (MSTP95.GE.8) THEN
C...  For Paquis type, make "histogram" of string densities along thrust axis
        RAPMIN = -RAPMAX
        DRAP   = 2*RAPMAX/(1D0*NBINY)
C...  Explicitly zero histogram bin content
        DO 160 IBINY=1,NBINY
          NSTRY(IBINY)=0
  160   CONTINUE
        DO 180 ISTR=1,NCR-1,2
          IC = ICR(ISTR)
          IA = ICR(ISTR+1)
          Y1 = MIN(RLOPTC(ISTR),RLOPTC(ISTR+1))
          Y2 = MAX(RLOPTC(ISTR),RLOPTC(ISTR+1))
          DO 170 IBINY=1,NBINY
            YBINLO = RAPMIN + (IBINY-1)*DRAP
C...  If bin inside string piece, add 1 in this bin
C...  (Strictly speaking: if it starts before midpoint and ends after midpoint)
            IF (Y1.LE.YBINLO+0.5*DRAP.AND.Y2.GE.YBINLO+0.5*DRAP)
     &           NSTRY(IBINY) = NSTRY(IBINY) + 1
  170     CONTINUE
  180   CONTINUE
C...  Loop over pieces to find individual reconnect probability
        DO 200 IS=1,NCR-1,2
          DNSUM  = 0D0
          DNAVG  = 0D0
C...Beginning at Y = RAPMIN = -RAPMAX, ending at Y = RAPMAX
          RBINLO = (MIN(RLOPTC(IS),RLOPTC(IS+1))-RAPMIN)/DRAP + 0.5
          RBINHI = (MAX(RLOPTC(IS),RLOPTC(IS+1))-RAPMIN)/DRAP + 0.5
C...Make sure integer bin numbers lie inside proper range
          IBINLO = MAX(1,MIN(NBINY,NINT(RBINLO)))
          IBINHI = MAX(1,MIN(NBINY,NINT(RBINHI)))
C...Size of rapidity bins (is < DRAP if piece smaller than one bin)
C...(also smaller than DRAP if a one-unit wide piece is stretched
C... over 2 bins, thus making the computation more accurate)
          DRAPAV = (RBINHI-RBINLO)/(IBINHI-IBINLO+1)*DRAP
C...  Decide whether to suppress reconnections in high-pT string pieces
          CRMODF = 1D0
          IF (PARP(77).GT.0D0) THEN
C...  Total string piece energy, momentum squared, and components
            EES  =  P(ICR(IS),4) + P(ICR(IS+1),4)
            PPS2 = (P(ICR(IS),1)+ P(ICR(IS+1),1))**2
     &           + (P(ICR(IS),2)+ P(ICR(IS+1),2))**2
     &           + (P(ICR(IS),3)+ P(ICR(IS+1),3))**2
            PZTS = P(ICR(IS),1)*TX+P(ICR(IS),2)*TY+P(ICR(IS),3)*TZ
     &           + P(ICR(IS+1),1)*TX+P(ICR(IS+1),2)*TY+P(ICR(IS+1),3)*TZ
            PTTS = SQRT(PPS2 - PZTS**2)
C...  Mass of string piece in units of mpi (at least 1)
            RMPI2  = 0.135D0
            RM2STR = MAX(RMPI2,EES**2 - PPS2)
C...  Estimate number of pions ~ log(M2) (at least 1)
            RNPI   = LOG(RM2STR/RMPI2)+1D0
            PT2AVG = (PTTS / RNPI)**2
C...  Supress reconnection probability by 1/(1+P77*P2AVG)
            CRMODF=1D0/(1D0+PARP(77)**2*PT2AVG)
          ENDIF
          PKEEP = 1.0
          DO 190 IBINY=IBINLO,IBINHI
C            DNSUM = DNSUM + 1D0
            DNOVL = MAX(0,NSTRY(IBINY)-1)
            PKEEP = PKEEP * (1D0-CRMODF*PARP(78))**(DRAPAV*DNOVL)
C            DNAVG = DNAVG + MAX(1,NSTRY(IBINY))
  190     CONTINUE
C          DNAVG = DNAVG / DNSUM
C...  If keeping string piece, save
          IF (PYR(0).LE.PKEEP) THEN
            LCT = LCT+1
            MCN(ICR(IS),1)=LCT
            MCN(ICR(IS+1),2)=LCT
          ENDIF
  200   CONTINUE
      ENDIF
 
C...Skip if there is only one possibility
      IF (NCR.LE.2) THEN
        GOTO 9999
      ENDIF
 
C...Reorder, so ordered in I (in order to correspond to old algorithm)
      NLOOP=0
  210 NLOOP=NLOOP+1
      MORD=1
      DO 220 IC1=1,NCR-1
        I1=ICR(IC1)
        I2=ICR(IC1+1)
        IF (I1.GT.I2) THEN
          IT=I1
          MST=MSCR(IC1)
          ICR(IC1)=I2
          MSCR(IC1)=MSCR(IC1+1)
          ICR(IC1+1)=IT
          MSCR(IC1+1)=MST
          MORD=0
        ENDIF
  220 CONTINUE
C...Max do 1000 reordering loops
      IF (MORD.EQ.0.AND.NLOOP.LE.1000) GOTO 210
 
C...PS: 03 May 2010
C...For Seattle and Paquis types, check if there is a dangling tag
C...Needed for special case when entire reconnected state was one or
C...more gluon loops in original topology in which case these CR
C...algorithms need to be told they shouldn't look for a dangling tag.
      M3FREE=0
      IF (MSTP95.GE.6.AND.MSTP95.LE.9) THEN
        DO 230 IC1=1,NCR
          I1=ICR(IC1)
C...Color charge
          MCI=KCHG(PYCOMP(K(I1,2)),2)*ISIGN(1,K(I1,2))
          IF (MCI.EQ.1.AND.MCN(I1,1).EQ.0) M3FREE=1
          IF (MCI.EQ.-1.AND.MCN(I1,2).EQ.0) M3FREE=1
          IF (MCI.EQ.2) THEN
            IF (MCN(I1,1).NE.0.AND.MCN(I1,2).EQ.0) M3FREE=1
            IF (MCN(I1,2).NE.0.AND.MCN(I1,1).EQ.0) M3FREE=1
          ENDIF
  230   CONTINUE
      ENDIF
 
C...Loop over CR partons
C...(Ignore junctions for now.)
      NLOOP=0
  240 NLOOP=NLOOP+1
      RLMAX=0D0
      ICRMAX=0
C...Loop over coloured partons
      DO 260 IC1=1,NCR
C...Retrieve parton Event Record index and Colour Side
        I=ICR(IC1)
        MSI=MSCR(IC1)
C...Skip already connected partons
        IF (MCN(I,MSI).NE.0) GOTO 260
C...Shorthand for colour charge
        MCI=KCHG(PYCOMP(K(I,2)),2)*ISIGN(1,K(I,2))
C...For Seattle algorithm, only start from partons with one dangling
C...colour tag (unless there aren't any, cf. M3FREE above.)
        IF (MSTP(95).GE.6.AND.MSTP(95).LE.9) THEN
          IF (MCI.EQ.2.AND.MCN(I,1).EQ.0.AND.MCN(I,2).EQ.0
     &         .AND.M3FREE.EQ.1) THEN
            GOTO 260
          ENDIF
        ENDIF
C...Retrieve saved optimal partner
        IO=IOPT(IC1)
        IF (IO.NE.0) THEN
C...Reject saved optimal partner if latter is now connected
C...(Also reject if using model S1, since saved partner may
C...now give rise to gg loop.)
          IF (MCN(IO,3-MSI).NE.0.OR.MSTP(95).LE.3) THEN
            IOPT(IC1)=0
            RLOPTC(IC1)=1D19
          ENDIF
        ENDIF
        RLOPT=RLOPTC(IC1)
C...Search for new optimal partner if necessary
        IF (IOPT(IC1).EQ.0) THEN
          MBROPT=0
          MGGOPT=0
          RLOPT=1D19
C...Loop over partons you can connect to
          DO 250 IC2=1,NCR
            J=ICR(IC2)
            MSJ=MSCR(IC2)
C...Skip if already connected
            IF (MCN(J,MSJ).NE.0) GOTO 250
C...Skip if this not colour-anticolour pair
            IF (MSI.EQ.MSJ) GOTO 250
C...And do not let gluons connect to themselves
            IF (I.EQ.J) GOTO 250
C...Suppress direct connections between partons in same Beam Remnant
            MBRSTR=0
            IF (K(I,3).LE.2.AND.K(I,3).GE.1.AND.K(I,3).EQ.K(J,3))
     &          MBRSTR=1
C...Shorthand for colour charge
            MCJ=KCHG(PYCOMP(K(J,2)),2)*ISIGN(1,K(J,2))
C...Check for gluon loops
            MGGSTR=0
            IF (MCJ.EQ.2.AND.MCI.EQ.2) THEN
              IF (MCN(I,2).EQ.MCN(J,1).AND.MSTP(95).LE.3.AND.
     &            MCN(I,2).NE.0) MGGSTR=1
            ENDIF
C...Save connection with smallest lambda measure
            RL=FOUR(I,J)
C...If best so far was a BR string and this is not, also save.
C...If best so far was a gg string and this is not, also save.
C...NB: this is not fool-proof. If the algorithm finds a BR or gg
C...string with a small Lambda measure as the last step, this connection
C...will be saved regardless of whether other possibilities existed.
C...I.e., there should really be a check whether another possibility has
C...already been found, but since these models are now actively in use
C...and uncertainties are anyway large, the algorithm is left as it is.
C...(correction --> Pythia 8 ?)
            IF (RL.LT.RLOPT.OR.(RL.EQ.RLOPT.AND.PYR(0).LE.0.5D0)
     &          .OR.(MBROPT.EQ.1.AND.MBRSTR.EQ.0)
     &          .OR.(MGGOPT.EQ.1.AND.MGGSTR.EQ.0)) THEN
C...Paquis type: fix problem above
              MPAQ = 0
              IF (MSTP95.GE.8.AND.RLOPT.LE.1D18) THEN
                IF (MBRSTR.EQ.1.AND.MBROPT.EQ.0) MPAQ=1
                IF (MGGSTR.EQ.1.AND.MGGOPT.EQ.0) MPAQ=1
              ENDIF
              IF (MPAQ.EQ.0) THEN
                RLOPT=RL
                RLOPTC(IC1)=RLOPT
                IOPT(IC1)=J
                MBROPT=MBRSTR
                MGGOPT=MGGSTR
              ENDIF
            ENDIF
  250     CONTINUE
        ENDIF
        IF (IOPT(IC1).NE.0) THEN
C...Save pair with largest RLOPT so far
          IF (RLOPT.GE.RLMAX) THEN
            ICRMAX=IC1
            RLMAX=RLOPT
          ENDIF
        ENDIF
  260 CONTINUE
C...Save and iterate
      ICMAX=0
      IF (ICRMAX.GT.0) THEN
        LCT=LCT+1
        ILMAX=ICR(ICRMAX)
        JLMAX=IOPT(ICRMAX)
        ICMAX=MSCR(ICRMAX)
        JCMAX=3-ICMAX
        MCN(ILMAX,ICMAX)=LCT
        MCN(JLMAX,JCMAX)=LCT
        IF (NLOOP.LE.2*(N-IP)) THEN
          GOTO 240
        ELSE
          CALL PYERRM(31,' PYFSCR: infinite loop in color annealing')
          CALL PYSTOP(11)
        ENDIF
      ELSE
C...Save and exit. First check for leftover gluon(s)
        DO 290 I=MAX(1,IP),N
C...Check colour charge
          MCI=KCHG(PYCOMP(K(I,2)),2)*ISIGN(1,K(I,2))
          IF (K(I,1).NE.3.OR.MCI.NE.2) GOTO 290
          IF(MCN(I,1).EQ.0.AND.MCN(I,2).EQ.0) THEN
C...Decide where to put left-over gluon (minimal insertion)
            ICMAX=0
            RLMAX=1D19
C...PS: Bug fix 30 Apr 2010: try all lines, not just reconnected ones
            DO 280 KCT=ICTMIN,LCT
              IC=0
              IA=0
              DO 270 IT=MAX(1,IP),N
                IF (IT.EQ.I.OR.K(IT,1).NE.3) GOTO 270
                IF (MCN(IT,1).EQ.KCT) IC=IT
                IF (MCN(IT,2).EQ.KCT) IA=IT
  270         CONTINUE
C...Skip if this color tag no longer present in event record
              IF (IC.EQ.0.OR.IA.EQ.0) GOTO 280
              RL=FOUR(IC,I)*FOUR(IA,I)
              IF (RL.LT.RLMAX) THEN
                RLMAX=RL
                ICMAX=IC
                IAMAX=IA
              ENDIF
  280       CONTINUE
            LCT=LCT+1
            MCN(I,1)=MCN(ICMAX,1)
            MCN(I,2)=LCT
            MCN(ICMAX,1)=LCT
          ENDIF
  290   CONTINUE
C...Here we need to loop over entire event.
        DO 300 IZ=MAX(1,IP),N
C...Do not erase parton shower colour history
          IF (K(IZ,1).NE.3) GOTO 300
C...Check colour charge
          MCI=KCHG(PYCOMP(K(IZ,2)),2)*ISIGN(1,K(IZ,2))
          IF (MCI.EQ.0) GOTO 300
          IF (MCN(IZ,1).NE.0) MCT(IZ,1)=MCN(IZ,1)
          IF (MCN(IZ,2).NE.0) MCT(IZ,2)=MCN(IZ,2)
  300   CONTINUE
      ENDIF
 
 9999 RETURN
      END
