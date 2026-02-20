 
C*********************************************************************
 
C...PYMIHK
C...Finds left-behind remnant flavour content and hooks up
C...the colour flow between the hard scattering and remnants
 
      SUBROUTINE PYMIHK
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...The event record
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Parameters
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
C...The common block of dangling ends
      COMMON/PYINTM/KFIVAL(2,3),NMI(2),IMI(2,800,2),NVC(2,-6:6),
     &     XASSOC(2,-6:6,240),XPSVC(-6:6,-1:240),PVCTOT(2,-1:1),
     &     XMI(2,240),PT2MI(240),IMISEP(0:240)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/,/PYINTM/
C...Local variables
      PARAMETER (NERSIZ=4000)
      COMMON /PYCBLS/MCO(NERSIZ,2),NCC,JCCO(NERSIZ,2),JCCN(NERSIZ,2)
     &     ,MACCPT
      COMMON /PYCTAG/NCT,MCT(NERSIZ,2)
      SAVE /PYCBLS/,/PYCTAG/
      DIMENSION JST(2,3),IV(2,3),IDQ(3),NVSUM(2),NBRTOT(2),NG(2)
     &     ,ITJUNC(2),MOUT(2),INSR(1000,3),ISTR(6),YMI(240)
      DATA NERRPR/0/
      SAVE NERRPR
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,3)*P(J,3)-P(I,2)*P(J,2)-P(I,1)*P(J,1)
 
C...Set up error checkers
      IBOOST=0
 
C...Initialize colour arrays: MCO (Original) and MCT (New)
      DO 110 I=MINT(84)+1,NERSIZ
        DO 100 JC=1,2
          MCT(I,JC)=0
          MCO(I,JC)=0
  100   CONTINUE
C...Also zero colour tracing information, if existed.
        IF (I.LE.N) THEN
          K(I,4)=MOD(K(I,4),MSTU(5)**2)
          K(I,5)=MOD(K(I,5),MSTU(5)**2)
        ENDIF
  110 CONTINUE
 
C...Initialize colour tag collapse arrays:
C...JCCO (Original) and JCCN (New).
      DO 130 MG=MINT(84)+1,NERSIZ
        DO 120 JC=1,2
          JCCO(MG,JC)=0
          JCCN(MG,JC)=0
  120   CONTINUE
  130 CONTINUE
 
C...Zero gluon insertion array
      DO 150 IM=1,1000
        DO 140 J=1,3
          INSR(IM,J)=0
  140   CONTINUE
  150 CONTINUE
 
C...Compute hard scattering system rapidities
      IF (MSTP(89).EQ.1) THEN
        DO 160 IM=1,240
          IF (IM.LE.MINT(31)) THEN
            YMI(IM)=LOG(XMI(1,IM)/XMI(2,IM))
          ELSE
C...Set (unsigned) rapidity = 100 for beam remnant systems.
            YMI(IM)=100D0
          ENDIF
  160   CONTINUE
      ENDIF
 
C...Treat each side separately
      DO 290 JS=1,2
 
C...Initialize side.
        NG(JS)=0
        JV=0
        KFS=ISIGN(1,MINT(10+JS))
 
C...Set valence content of pi0, gamma, K0S, K0L if not yet done.
        IF(KFIVAL(JS,1).EQ.0) THEN
          IF(MINT(10+JS).EQ.111) THEN
            KFIVAL(JS,1)=INT(1.5D0+PYR(0))
            KFIVAL(JS,2)=-KFIVAL(JS,1)
          ELSEIF(MINT(10+JS).EQ.22) THEN
            PYRKF=PYR(0)
            KFIVAL(JS,1)=1
            IF(PYRKF.GT.0.1D0) KFIVAL(JS,1)=2
            IF(PYRKF.GT.0.5D0) KFIVAL(JS,1)=3
            IF(PYRKF.GT.0.6D0) KFIVAL(JS,1)=4
            KFIVAL(JS,2)=-KFIVAL(JS,1)
          ELSEIF(MINT(10+JS).EQ.130.OR.MINT(10+JS).EQ.310) THEN
            IF(PYR(0).GT.0.5D0) THEN
              KFIVAL(JS,1)=1
              KFIVAL(JS,2)=-3
            ELSE
              KFIVAL(JS,1)=3
              KFIVAL(JS,2)=-1
            ENDIF
          ENDIF
        ENDIF
 
C...Initialize beam remnant sea and valence content flavour by flavour.
        NVSUM(JS)=0
        NBRTOT(JS)=0
        DO 210 JFA=1,6
C...Count up original number of JFA valence quarks and antiquarks.
          NVALQ=0
          NVALQB=0
          NSEA=0
          DO 170 J=1,3
            IF(KFIVAL(JS,J).EQ.JFA) NVALQ=NVALQ+1
            IF(KFIVAL(JS,J).EQ.-JFA) NVALQB=NVALQB+1
  170     CONTINUE
          NVSUM(JS)=NVSUM(JS)+NVALQ+NVALQB
C...Subtract kicked out valence and determine sea from flavour cons.
          DO 180 IM=1,NMI(JS)
            IFL = K(IMI(JS,IM,1),2)
            IFA = IABS(IFL)
            IFS = ISIGN(1,IFL)
            IF (IFL.EQ.JFA.AND.IMI(JS,IM,2).EQ.0) THEN
C...Subtract K.O. valence quark from remainder.
              NVALQ=NVALQ-1
              JV=NVSUM(JS)-NVALQ-NVALQB
              IV(JS,JV)=IMI(JS,IM,1)
            ELSEIF (IFL.EQ.-JFA.AND.IMI(JS,IM,2).EQ.0) THEN
C...Subtract K.O. valence antiquark from remainder.
              NVALQB=NVALQB-1
              JV=NVSUM(JS)-NVALQ-NVALQB
              IV(JS,JV)=IMI(JS,IM,1)
            ELSEIF (IFA.EQ.JFA) THEN
C...Outside sea without companion: add opposite sea flavour inside.
              IF (IMI(JS,IM,2).LT.0) NSEA=NSEA-IFS
            ENDIF
  180     CONTINUE
C...Check if space left in PYJETS for additional BR flavours
          NFLSUM=IABS(NSEA)+NVALQ+NVALQB
          NBRTOT(JS)=NBRTOT(JS)+NFLSUM
          IF (N+NFLSUM+1.GT.MSTU(4)) THEN
            CALL PYERRM(11,'(PYMIHK:) no more memory left in PYJETS')
            MINT(51)=1
            RETURN
          ENDIF
C...Add required val+sea content to beam remnant.
          IF (NFLSUM.GT.0) THEN
            DO 200 IA=1,NFLSUM
C...Insert beam remnant quark as p.t. symbolic parton in ER.
              N=N+1
              DO 190 IX=1,5
                K(N,IX)=0
                P(N,IX)=0D0
                V(N,IX)=0D0
  190         CONTINUE
              K(N,1)=3
              K(N,2)=ISIGN(JFA,NSEA)
              IF (IA.LE.NVALQ) K(N,2)=JFA
              IF (IA.GT.NVALQ.AND.IA.LE.NVALQ+NVALQB) K(N,2)=-JFA
              K(N,3)=MINT(83)+JS
C...Also update NMI, IMI, and IV arrays.
              NMI(JS)=NMI(JS)+1
              IMI(JS,NMI(JS),1)=N
              IMI(JS,NMI(JS),2)=-1
              IF (IA.LE.NVALQ+NVALQB) THEN
                IMI(JS,NMI(JS),2)=0
                JV=JV+1
                IV(JS,JV)=IMI(JS,NMI(JS),1)
              ENDIF
  200       CONTINUE
          ENDIF
  210   CONTINUE
 
        IM=0
  220   IM=IM+1
        IF (IM.LE.NMI(JS)) THEN
          IF (K(IMI(JS,IM,1),2).EQ.21) THEN
            NG(JS)=NG(JS)+1
C...Add fictitious parent gluons for companion pairs.
          ELSEIF (IMI(JS,IM,2).NE.0.AND.K(IMI(JS,IM,1),2).GT.0) THEN
C...Randomly assign companions to sea quarks which have none.
            IF (IMI(JS,IM,2).LT.0) THEN
              IMC=PYR(0)*NMI(JS)
  230         IMC=MOD(IMC,NMI(JS))+1
              IF (K(IMI(JS,IMC,1),2).NE.-K(IMI(JS,IM,1),2)) GOTO 230
              IF (IMI(JS,IMC,2).GE.0) GOTO 230
              IMI(JS, IM,2) = IMI(JS,IMC,1)
              IMI(JS,IMC,2) = IMI(JS, IM,1)
            ENDIF
C...Add fictitious parent gluon
            N=N+1
            DO 240 IX=1,5
              K(N,IX)=0
              P(N,IX)=0D0
              V(N,IX)=0D0
  240       CONTINUE
            K(N,1)=14
            K(N,2)=21
            K(N,3)=MINT(83)+JS
C...Set gluon (anti-)colour daughter pointers
            K(N,4)=IMI(JS, IM,1)
            K(N,5)=IMI(JS, IM,2)
C...Set quark (anti-)colour parent pointers
            K(IMI(JS, IM,2),5)=K(IMI(JS, IM,2),5)+MSTU(5)*N
            K(IMI(JS, IM,1),4)=K(IMI(JS, IM,1),4)+MSTU(5)*N
C...Add gluon to IMI
            NMI(JS)=NMI(JS)+1
            IMI(JS,NMI(JS),1)=N
            IMI(JS,NMI(JS),2)=0
          ENDIF
          GOTO 220
        ENDIF
 
C...If incoming (anti-)baryon, insert inside (anti-)junction.
C...Set up initial v-v-j-v configuration. Otherwise set up
C...mesonic v-vbar configuration
        IF (IABS(MINT(10+JS)).GT.1000) THEN
C...Determine junction type (1: B=1 2: B=-1)
          ITJUNC(JS) = (3-KFS)/2
C...Insert junction.
          N=N+1
          DO 250 IX=1,5
            K(N,IX)=0
            P(N,IX)=0D0
            V(N,IX)=0D0
  250     CONTINUE
C...Set special junction codes:
          K(N,1)=42
          K(N,2)=88
C...Set parent to side.
          K(N,3)=MINT(83)+JS
          K(N,4)=ITJUNC(JS)*MSTU(5)
          K(N,5)=0
C...Connect valence quarks to junction.
          MOUT(JS)=0
          MANTI=ITJUNC(JS)-1
C...Set (anti)colour mother = junction.
          DO 260 JV=1,3
            K(IV(JS,JV),4+MANTI)=MOD(K(IV(JS,JV),4+MANTI),MSTU(5))
     &           +MSTU(5)*N
C...Keep track of partons adjacent to junction:
            JST(JS,JV)=IV(JS,JV)
  260     CONTINUE
        ELSE
C...Mesons: set up initial q-qbar topology
          ITJUNC(JS)=0
          IF (K(IV(JS,1),2).GT.0) THEN
            IQ=IV(JS,1)
            IQBAR=IV(JS,2)
          ELSE
            IQ=IV(JS,2)
            IQBAR=IV(JS,1)
          ENDIF
          IV(JS,3)=0
          JST(JS,1)=IQ
          JST(JS,2)=IQBAR
          JST(JS,3)=0
          K(IQ,4)=MOD(K(IQ,4),MSTU(5))+MSTU(5)*IQBAR
          K(IQBAR,5)=MOD(K(IQBAR,5),MSTU(5))+MSTU(5)*IQ
C...Special for mesons. Insert gluon if BR empty.
          IF (NBRTOT(JS).EQ.0) THEN
            N=N+1
            DO 270 IX=1,5
              K(N,IX)=0
              P(N,IX)=0D0
              V(N,IX)=0D0
  270       CONTINUE
            K(N,1)=3
            K(N,2)=21
            K(N,3)=MINT(83)+JS
            K(N,4)=0
            K(N,5)=0
            NBRTOT(JS)=1
            NG(JS)=NG(JS)+1
C...Add gluon to IMI
            NMI(JS)=NMI(JS)+1
            IMI(JS,NMI(JS),1)=N
            IMI(JS,NMI(JS),2)=0
          ENDIF
          MOUT(JS)=0
        ENDIF
 
C...Count up number of valence quarks outside BR.
        DO 280 JV=1,3
          IF (JST(JS,JV).LE.MINT(53).AND.JST(JS,JV).GT.0)
     &         MOUT(JS)=MOUT(JS)+1
  280   CONTINUE
 
  290 CONTINUE
 
C...Now both sides have been prepared in an initial vvjv (baryonic) or
C...v(g)vbar (mesonic) configuration.
 
C...Create colour line tags starting from initiators.
      NCT=0
      DO 320 IM=1,MINT(31)
C...Consider each side in turn.
        DO 310 JS=1,2
          I1=IMI(JS,IM,1)
          I2=IMI(3-JS,IM,1)
          DO 300 JCS=4,5
            IF (K(I1,2).NE.21.AND.(9-2*JCS).NE.ISIGN(1,K(I1,2)))
     &           GOTO 300
            IF (K(I1,JCS)/MSTU(5)**2.NE.0) GOTO 300
 
            KCS=JCS
            CALL PYCTTR(I1,KCS,I2)
            IF(MINT(51).NE.0) RETURN
 
  300     CONTINUE
  310   CONTINUE
  320 CONTINUE
 
      DO 340 JS=1,2
C...Create colour tags for beam remnant partons.
        DO 330 IM=MINT(31)+1,NMI(JS)
          IP=IMI(JS,IM,1)
          IF (K(IP,2).NE.21) THEN
            JC=(3-ISIGN(1,K(IP,2)))/2
            IF (MCT(IP,JC).EQ.0) THEN
              NCT=NCT+1
              MCT(IP,JC)=NCT
            ENDIF
          ELSE
C...Gluons
            ICD=K(IP,4)
            IAD=K(IP,5)
            IF (ICD.NE.0) THEN
C...Fictituous gluons just inherit from their quark daughters.
              ICC=MCT(ICD,1)
              IAC=MCT(IAD,2)
            ELSE
C...Real beam remnant gluons get their own colours
              ICC=NCT+1
              IAC=NCT+2
              NCT=NCT+2
            ENDIF
            MCT(IP,1)=ICC
            MCT(IP,2)=IAC
          ENDIF
  330   CONTINUE
  340 CONTINUE
 
C...Create colour tags for colour lines which are detached from the
C...initial state.
 
      DO 360 MQGST=1,2
        DO 350 I=MINT(84)+1,N
 
C...Look for coloured string endpoint, or (later) leftover gluon.
          IF (K(I,1).NE.3) GOTO 350
          KC=PYCOMP(K(I,2))
          IF(KC.EQ.0) GOTO 350
          KQ=KCHG(KC,2)
          IF(KQ.EQ.0.OR.(MQGST.EQ.1.AND.KQ.EQ.2)) GOTO 350
 
C...Pick up loose string end with no previous tag.
          KCS=4
          IF(KQ*ISIGN(1,K(I,2)).LT.0) KCS=5
          IF(MCT(I,KCS-3).NE.0) GOTO 350
 
          CALL PYCTTR(I,KCS,I)
          IF(MINT(51).NE.0) RETURN
 
  350   CONTINUE
  360 CONTINUE
 
C...Store original colour tags
      DO 370 I=MINT(84)+1,N
        MCO(I,1)=MCT(I,1)
        MCO(I,2)=MCT(I,2)
  370 CONTINUE
 
C...Iteratively add gluons to already existing string pieces, enforcing
C...various possible orderings, and rejecting insertions that would give
C...rise to singlet gluons.
C...<kappa tau> normalization.
      RM0=1.5D0
      MRETRY=0
      PARP80=PARP(80)
 
C...Set up simplified kinematics.
C...Boost hard interaction systems.
      IBOOST=IBOOST+1
      DO 380 IM=1,MINT(31)
        BETA=(XMI(1,IM)-XMI(2,IM))/(XMI(1,IM)+XMI(2,IM))
        CALL PYROBO(IMISEP(IM-1)+1,IMISEP(IM),0D0,0D0,0D0,0D0,BETA)
  380 CONTINUE
C...Assign preliminary beam remnant momenta.
      DO 390 I=MINT(53)+1,N
        JS=K(I,3)
        P(I,1)=0D0
        P(I,2)=0D0
        IF (K(I,2).NE.88) THEN
          P(I,4)=0.5D0*VINT(142+JS)*VINT(1)/MAX(1,NMI(JS)-MINT(31))
          P(I,3)=P(I,4)
          IF (JS.EQ.2) P(I,3)=-P(I,3)
        ELSE
C...Junctions are wildcards for the present.
          P(I,4)=0D0
          P(I,3)=0D0
        ENDIF
  390 CONTINUE
 
C...Reset colour processing information.
  400 DO 410 I=MINT(84)+1,N
        K(I,4)=MOD(K(I,4),MSTU(5)**2)
        K(I,5)=MOD(K(I,5),MSTU(5)**2)
  410 CONTINUE
 
      NCC=0
      DO 430 JS=1,2
C...If meson,  without gluon in BR, collapse q-qbar colour tags:
        IF (ITJUNC(JS).EQ.0) THEN
          JC1=MCT(JST(JS,1),1)
          JC2=MCT(JST(JS,2),2)
          NCC=NCC+1
          JCCO(NCC,1)=MAX(JC1,JC2)
          JCCO(NCC,2)=MIN(JC1,JC2)
C...Collapse colour tags in event record
          DO 420 I=MINT(84)+1,N
            IF (MCT(I,1).EQ.JCCO(NCC,1)) MCT(I,1)=JCCO(NCC,2)
            IF (MCT(I,2).EQ.JCCO(NCC,1)) MCT(I,2)=JCCO(NCC,2)
  420     CONTINUE
        ENDIF
  430 CONTINUE
 
  440 JS=1
      IF (PYR(0).GT.0.5D0.OR.NG(1).EQ.0) JS=2
      IF (NG(JS).GT.0) THEN
        NOPT=0
        RLOPT=1D9
C...Start at random gluon (optimizes speed for random attachments)
        NMGL=0
        IMGL=PYR(0)*NMI(JS)+1
  450   IMGL=MOD(IMGL,NMI(JS))+1
        NMGL=NMGL+1
C...Only loop through NMI once (with upper limit to save time)
        IF (NMGL.LE.NMI(JS).AND.NOPT.LE.3) THEN
          IGL  = IMI(JS,IMGL,1)
C...If not gluon or if already connected, try next.
          IF (K(IGL,2).NE.21.OR.K(IGL,4)/MSTU(5).NE.0
     &         .OR.K(IGL,5)/MSTU(5).NE.0) GOTO 450
C...Now loop through all possible insertions of this gluon.
          NMP1=0
          IMP1=PYR(0)*NMI(JS)+1
  460     IMP1=MOD(IMP1,NMI(JS))+1
          NMP1=NMP1+1
          IF (IMP1.EQ.IMGL) GOTO 460
C...Only loop through NMI once (with upper limit to save time).
          IF (NMP1.LE.NMI(JS).AND.NOPT.LE.3) THEN
            IP1  = IMI(JS,IMP1,1)
C...Try both colour mother and colour anti-mother.
C...Randomly select which one to try first.
            NANTI=0
            MANTI=PYR(0)*2
  470       MANTI=MOD(MANTI+1,2)
            NANTI=NANTI+1
            IF (NANTI.LE.2) THEN
              IP2 =MOD(K(IP1,4+MANTI)/MSTU(5),MSTU(5))
C...Reject if no appropriate mother (or if mother is fictitious
C...parent gluon.)
              IF (IP2.LE.0) GOTO 470
              IF (K(IP2,2).EQ.21.AND.IP2.GT.MINT(53)) GOTO 470
C...Also reject if this link has already been tried.
              IF (K(IP1,4+MANTI)/MSTU(5)**2.EQ.2) GOTO 470
              IF (K(IP2,5-MANTI)/MSTU(5)**2.EQ.2) GOTO 470
C...Set flag to indicate that this link has now been tried for this
C...gluon. IP2 may be junction, which has several mothers.
              K(IP1,4+MANTI)=K(IP1,4+MANTI)+2*MSTU(5)**2
              IF (K(IP2,2).NE.88) THEN
                K(IP2,5-MANTI)=K(IP2,5-MANTI)+2*MSTU(5)**2
              ENDIF
 
C...JCG1: Original colour tag of gluon on IP1 side
C...JCG2: Original colour tag of gluon on IP2 side
C...JCP1: Original colour tag of IP1 on gluon side
C...JCP2: Original colour tag of IP2 on gluon side.
              JCG1=MCO(IGL,2-MANTI)
              JCG2=MCO(IGL,1+MANTI)
              JCP1=MCO(IP1,1+MANTI)
              JCP2=MCO(IP2,2-MANTI)
 
              CALL PYMIHG(JCP1,JCG1,JCP2,JCG2)
C...Reject gluon attachments that give rise to singlet gluons.
              IF (MACCPT.EQ.0) GOTO 470
 
C...Update colours
              JCG1=MCT(IGL,2-MANTI)
              JCG2=MCT(IGL,1+MANTI)
              JCP1=MCT(IP1,1+MANTI)
              JCP2=MCT(IP2,2-MANTI)
 
C...Select whether to accept this insertion
              IF (MSTP(89).EQ.0) THEN
C...Random insertions: no measure.
                RL=1D0
C...For random ordering, we want to suppress beam remnant breakups
C...already at this point.
                IF (IP1.GT.MINT(53).AND.IP2.GT.MINT(53)
     &               .AND.MOUT(JS).NE.0.AND.PYR(0).GT.PARP80) THEN
                  NMP1=0
                  NMGL=0
                  GOTO 470
                ENDIF
              ELSEIF (MSTP(89).EQ.1) THEN
C...Rapidity ordering:
C...YGL = Rapidity of gluon.
                YGL=YMI(IMGL)
C...If fictitious gluon
                IF (YGL.EQ.100D0) THEN
                  YGL=(3-2*JS)*100D0
                  IDA1=MOD(K(IGL,4),MSTU(5))
                  IDA2=MOD(K(IGL,5),MSTU(5))
                  DO 480 IMT=1,NMI(JS)
C...Select (arbitrarily) the most central daughter.
                    IF (IMI(JS,IMT,1).EQ.IDA1.OR.IMI(JS,IMT,1).EQ.IDA2)
     &                   THEN
                      IF (ABS(YGL).GT.ABS(YMI(IMT))) YGL=YMI(IMT)
                    ENDIF
  480             CONTINUE
                ENDIF
C...YP1 = Rapidity IP1
                YP1=YMI(IMP1)
C...If fictitious gluon
                IF (YP1.EQ.100D0) THEN
                  YP1=(3-2*JS)*YP1
                  IDA1=MOD(K(IP1,4),MSTU(5))
                  IDA2=MOD(K(IP1,5),MSTU(5))
                  DO 490 IMT=1,NMI(JS)
C...Select (arbitrarily) the most central daughter.
                    IF (IMI(JS,IMT,1).EQ.IDA1.OR.IMI(JS,IMT,1).EQ.IDA2)
     &                   THEN
                      IF (ABS(YP1).GT.ABS(YMI(IMT))) YP1=YMI(IMT)
                    ENDIF
  490             CONTINUE
                ENDIF
C...YP2 = Rapidity of mother system
                IF (K(IP2,2).NE.88) THEN
                  DO 500 IMT=1,NMI(JS)
                    IF (IMI(JS,IMT,1).EQ.IP2) YP2=YMI(IMT)
  500             CONTINUE
C...If fictitious gluon
                  IF (YP2.EQ.100D0) THEN
                    YP2=(3-2*JS)*YP2
                    IDA1=MOD(K(IP2,4),MSTU(5))
                    IDA2=MOD(K(IP2,5),MSTU(5))
                    DO 510 IMT=1,NMI(JS)
C...Select (arbitrarily) the most central daughter.
                      IF (IMI(JS,IMT,1).EQ.IDA1.OR.IMI(JS,IMT,1).EQ.IDA2
     &                     ) THEN
                        IF (ABS(YP2).GT.ABS(YMI(IMT))) YP2=YMI(IMT)
                      ENDIF
  510               CONTINUE
                  ENDIF
C...Assign (arbitrarily) 100D0 to junction also
                ELSE
                  YP2=(3-2*JS)*100D0
                ENDIF
                RL=ABS(YGL-YP1)+ABS(YGL-YP2)
              ELSEIF (MSTP(89).EQ.2) THEN
C...Lambda ordering:
C...Compute lambda measure for this insertion.
                RL=1D0
                DO 520 IST=1,6
                  ISTR(IST)=0
  520           CONTINUE
C...If IP2 is junction, not caught below.
                IF (JCP2.EQ.0) THEN
                  ITJU=MOD(K(IP2,4)/MSTU(5),MSTU(5))
C...Anti-junction is colour endpoint et vv., always on JCG2.
                  ISTR(5-ITJU)=IP2
                ENDIF
                DO 530 I=MINT(84)+1,N
                  IF (K(I,1).LT.10) THEN
C...The new string pieces
                    IF (MCT(I,1).EQ.JCG1) ISTR(1)=I
                    IF (MCT(I,2).EQ.JCG1) ISTR(2)=I
                    IF (MCT(I,1).EQ.JCG2) ISTR(3)=I
                    IF (MCT(I,2).EQ.JCG2) ISTR(4)=I
                  ENDIF
  530           CONTINUE
C...Also identify junctions as string endpoints.
                DO 540 I=MINT(84)+1,N
                  ICMO=MOD(K(I,4)/MSTU(5),MSTU(5))
                  IAMO=MOD(K(I,5)/MSTU(5),MSTU(5))
C...Find partons adjacent to junctions.
                  IF (ICMO.GT.0.AND.ICMO.LE.N) THEN
                    IF (K(ICMO,1).EQ.42.AND.MCT(I,1).EQ.JCG1.AND.ISTR(2)
     &                  .EQ.0) ISTR(2) = ICMO
                    IF (K(ICMO,1).EQ.42.AND.MCT(I,1).EQ.JCG2.AND.ISTR(4)
     &                  .EQ.0) ISTR(4) = ICMO
                  ENDIF
                  IF (IAMO.GT.0.AND.IAMO.LE.N) THEN
                    IF (K(IAMO,1).EQ.42.AND.MCT(I,2).EQ.JCG1.AND.ISTR(1)
     &                  .EQ.0) ISTR(1) = IAMO
                    IF (K(IAMO,1).EQ.42.AND.MCT(I,2).EQ.JCG2.AND.ISTR(3)
     &                  .EQ.0) ISTR(3) = IAMO
                  ENDIF
  540           CONTINUE
C...The old string piece
                ISTR(5)=ISTR(1+2*MANTI)
                ISTR(6)=ISTR(4-2*MANTI)
                IF (ISTR(1).EQ.0.OR.ISTR(2).EQ.0.OR.ISTR(3).EQ.0.OR.
     &              ISTR(4).EQ.0.OR.ISTR(5).EQ.0.OR.ISTR(6).EQ.0) THEN
C...If one or more of the colour tags for this connection is/are still
C...dangling, skip this attempt for the time being. 
                  RL=1D6
                ELSE
                  RL=MAX(1D0,FOUR(ISTR(1),ISTR(2)))*MAX(1D0,FOUR(ISTR(3)
     &                ,ISTR(4)))/MAX(1D0,FOUR(ISTR(5),ISTR(6)))
                  RL=LOG(RL)
                ENDIF
              ENDIF
C...Allow some breadth to speed things up.
              IF (ABS(1D0-RL/RLOPT).LT.0.05D0) THEN
                NOPT=NOPT+1
              ELSEIF (RL.GT.RLOPT) THEN
                GOTO 470
              ELSE
                NOPT=1
                RLOPT=RL
              ENDIF
C...INSR(NOPT,1)=Gluon colour mother
C...INSR(NOPT,2)=Gluon
C...INSR(NOPT,3)=Gluon anticolour mother
              IF (NOPT.GT.1000) GOTO 470
              INSR(NOPT,1+2*MANTI)=IP2
              INSR(NOPT,2)=IGL
              INSR(NOPT,3-2*MANTI)=IP1
              IF (MSTP(89).GT.0.OR.NOPT.EQ.0) GOTO 470
            ENDIF
            IF (MSTP(89).GT.0.OR.NOPT.EQ.0) GOTO 460
          ENDIF
C...Reset link test information.
          DO 550 I=MINT(84)+1,N
            K(I,4)=MOD(K(I,4),MSTU(5)**2)
            K(I,5)=MOD(K(I,5),MSTU(5)**2)
  550     CONTINUE
          IF (MSTP(89).GT.0.OR.NOPT.EQ.0) GOTO 450
        ENDIF
C...Now we have a list of best gluon insertions, none of which cause
C...singlets to arise. If list is empty, try again a few times. Note:
C...this should never happen if we have a meson with a gluon inserted
C...in the beam remnant, since that breaks up the colour line.
        IF (NOPT.EQ.0) THEN
C...Abandon BR-g-BR suppression for retries. This is not serious, it
C...just means we happened to start with trying a bad sequence.
          PARP80=1D0
          IF (MRETRY.LE.10.AND.(ITJUNC(1).NE.0.OR.JST(1,3).EQ.0).AND
     &         .(ITJUNC(2).NE.0.OR.JST(2,3).EQ.0)) THEN
            MRETRY=MRETRY+1
            DO 590 JS=1,2
              IF (ITJUNC(JS).NE.0) THEN
                JST(JS,1)=IV(JS,1)
                JST(JS,2)=IV(JS,2)
                JST(JS,3)=IV(JS,3)
C...Reset valence quark parent pointers
                DO 560 I=MINT(53)+1,N
                  IF (K(I,2).EQ.88.AND.K(I,3).EQ.JS) IJU=I
  560           CONTINUE
                MANTI=ITJUNC(JS)-1
C...Set (anti)colour mother = junction.
                DO 570 JV=1,3
                  K(IV(JS,JV),4+MANTI)=MOD(K(IV(JS,JV),4+MANTI),MSTU(5))
     &                 +MSTU(5)*IJU
  570           CONTINUE
              ELSE
C...Same for mesons. JST unchanged, so needn't be restored.
                IQ=JST(JS,1)
                IQBAR=JST(JS,2)
                K(IQ,4)=MOD(K(IQ,4),MSTU(5))+MSTU(5)*IQBAR
                K(IQBAR,5)=MOD(K(IQBAR,5),MSTU(5))+MSTU(5)*IQ
              ENDIF
C...Also reset gluon parent pointers.
              NG(JS)=0
              DO 580 IM=1,NMI(JS)
                I=IMI(JS,IM,1)
                IF (K(I,2).EQ.21) THEN
                  K(I,4)=MOD(K(I,4),MSTU(5))
                  K(I,5)=MOD(K(I,5),MSTU(5))
                  NG(JS)=NG(JS)+1
                ENDIF
  580         CONTINUE
  590       CONTINUE
C...Reset colour tags
            DO 600 I=MINT(84)+1,N
              MCT(I,1)=MCO(I,1)
              MCT(I,2)=MCO(I,2)
  600       CONTINUE
            GOTO 400
          ELSE
            IF(NERRPR.LT.5) THEN
              NERRPR=NERRPR+1
              CALL PYLIST(4)
              CALL PYERRM(19,'(PYMIHK:) No physical colour flow found!')
              WRITE(MSTU(11),*) 'NG:', NG,'   MOUT:', MOUT(JS)
            ENDIF
C...Kill event and start another.
            MINT(51)=1
            RETURN
          ENDIF
        ELSE
C...Select between insertions, suppressing insertions wholly in the BR.
          IIN=PYR(0)*NOPT+1
  610     IIN=MOD(IIN,NOPT)+1
          IF (INSR(IIN,1).GT.MINT(53).AND.INSR(IIN,3).GT.MINT(53)
     &         .AND.MOUT(JS).NE.0.AND.PYR(0).GT.PARP80) GOTO 610
        ENDIF
 
C...Now we know which gluon to insert where. Colour tags in JCCO and
C...colour connection information should be updated, NG(JS) should be
C...counted down, and a new loop performed if there are still gluons
C...left on any side.
        ICM=INSR(IIN,1)
        IACM=INSR(IIN,3)
        IGL=INSR(IIN,2)
C...JCG : Original gluon colour tag
C...JCAG: Original gluon anticolour tag.
C...JCM : Original anticolour tag of gluon colour mother
C...JACM: Original colour tag of gluon anticolour mother
        JCG=MCO(IGL,1)
        JCM=MCO(ICM,2)
        JACG=MCO(IGL,2)
        JACM=MCO(IACM,1)
 
        CALL PYMIHG(JACM,JACG,JCM,JCG)
        IF (MACCPT.EQ.0) THEN
          IF(NERRPR.LT.5) THEN
            NERRPR=NERRPR+1
            CALL PYLIST(4)
            CALL PYERRM(11,'(PYMIHK:) Unphysical colour flow!')
            WRITE(MSTU(11),*) 'attaching', IGL,' between', ICM, IACM
          ENDIF
C...Kill event and start another.
          MINT(51)=1
          RETURN
        ELSE
C...If everything went fine, store new JCCN in JCCO.
          NCC=NCC+1
          DO 620 ICC=1,NCC
            JCCO(ICC,1)=JCCN(ICC,1)
            JCCO(ICC,2)=JCCN(ICC,2)
  620     CONTINUE
        ENDIF
 
C...One gluon attached is counted as equivalent to one end outside.
        MOUT(JS)=1
C...Set IGL colour mother = ICM.
        K(IGL,4)=MOD(K(IGL,4),MSTU(5))+MSTU(5)*ICM
C...Set ICM anticolour mother = IGL colour.
        IF (K(ICM,2).NE.88) THEN
          K(ICM,5)=MOD(K(ICM,5),MSTU(5))+MSTU(5)*IGL
        ELSE
C...If ICM is junction, just update JST array for now.
          DO 630 MSJ=1,3
            IF (JST(JS,MSJ).EQ.IACM) JST(JS,MSJ)=IGL
  630     CONTINUE
        ENDIF
C...Set IGL anticolour mother = IACM.
        K(IGL,5)=MOD(K(IGL,5),MSTU(5))+MSTU(5)*IACM
C...Set IACM anticolour mother = IGL anticolour.
        IF (K(IACM,2).NE.88) THEN
          K(IACM,4)=MOD(K(IACM,4),MSTU(5))+MSTU(5)*IGL
        ELSE
C...If IACM is junction, just update JST array for now.
          DO 640 MSJ=1,3
            IF (JST(JS,MSJ).EQ.ICM) JST(JS,MSJ)=IGL
  640     CONTINUE
        ENDIF
C...Count down # unconnected gluons.
        NG(JS)=NG(JS)-1
      ENDIF
      IF (NG(1).GT.0.OR.NG(2).GT.0) GOTO 440
 
      DO 840 JS=1,2
C...Collapse fictitious gluons.
        DO 670 IGL=MINT(53)+1,N
          IF (K(IGL,2).EQ.21.AND.K(IGL,3).EQ.MINT(83)+JS.AND.
     &         K(IGL,1).EQ.14) THEN
            ICM=K(IGL,4)/MSTU(5)
            IAM=K(IGL,5)/MSTU(5)
            ICD=MOD(K(IGL,4),MSTU(5))
            IAD=MOD(K(IGL,5),MSTU(5))
C...Set gluon daughters pointing to gluon mothers
            K(IAD,5)=MOD(K(IAD,5),MSTU(5))+MSTU(5)*IAM
            K(ICD,4)=MOD(K(ICD,4),MSTU(5))+MSTU(5)*ICM
C...Set gluon mothers pointing to gluon daughters.
            IF (K(ICM,2).NE.88) THEN
              K(ICM,5)=MOD(K(ICM,5),MSTU(5))+MSTU(5)*ICD
            ELSE
C...Special case: mother=junction. Just update JST array for now.
              DO 650 MSJ=1,3
                IF (JST(JS,MSJ).EQ.IGL) JST(JS,MSJ)=ICD
  650         CONTINUE
            ENDIF
            IF (K(IAM,2).NE.88) THEN
              K(IAM,4)=MOD(K(IAM,4),MSTU(5))+MSTU(5)*IAD
            ELSE
              DO 660 MSJ=1,3
                IF (JST(JS,MSJ).EQ.IGL) JST(JS,MSJ)=IAD
  660         CONTINUE
            ENDIF
          ENDIF
  670   CONTINUE
 
C...Erase collapsed gluons from NMI and IMI (but keep them in ER)
        IM=NMI(JS)+1
  680   IM=IM-1
        IF (IM.GT.MINT(31).AND.K(IMI(JS,IM,1),2).NE.21) GOTO 680
        IF (IM.GT.MINT(31)) THEN
          NMI(JS)=NMI(JS)-1
          DO 690 IMR=IM,NMI(JS)
            IMI(JS,IMR,1)=IMI(JS,IMR+1,1)
            IMI(JS,IMR,2)=IMI(JS,IMR+1,2)
  690     CONTINUE
          GOTO 680
        ENDIF
 
C...Finally, connect junction.
        IF (ITJUNC(JS).NE.0) THEN
          DO 700 I=MINT(53)+1,N
            IF (K(I,2).EQ.88.AND.K(I,3).EQ.MINT(83)+JS) IJU=I
  700     CONTINUE
C...NBRJQ counts # of jq, NBRVQ # of jv, inside BR.
          NBRJQ =0
          NBRVQ =0
          DO 720 MSJ=1,3
            IDQ(MSJ)=0
C...Find jq with no glue inbetween inside beam remnant.
            IF (JST(JS,MSJ).GT.MINT(53).AND.IABS(K(JST(JS,MSJ),2)).LE.5)
     &           THEN
              NBRJQ=NBRJQ+1
C...Set IDQ = -I if q non-valence and = +I if q valence.
              IDQ(NBRJQ)=-JST(JS,MSJ)
              DO 710 JV=1,3
                IF (IV(JS,JV).EQ.JST(JS,MSJ)) THEN
                  IDQ(NBRJQ)=JST(JS,MSJ)
                  NBRVQ=NBRVQ+1
                ENDIF
  710         CONTINUE
            ENDIF
            I12=MOD(MSJ+1,2)
            I45=5
            IF (MSJ.EQ.3) I45=4
            K(IJU,I45)=K(IJU,I45)+(MSTU(5)**I12)*JST(JS,MSJ)
  720     CONTINUE
 
C...Check if diquark can be formed.
          IF ((MSTP(88).GE.0.AND.NBRVQ.GE.2).OR.(NBRJQ.GE.2.AND.MSTP(88)
     &         .GE.1)) THEN
C...If there is less than 2 valence quarks connected to junction
C...and MSTP(88)>1, use random non-valence quarks to fill up.
            IF (NBRVQ.LE.1) THEN
              NDIQ=NBRVQ
  730         JFLIP=NBRJQ*PYR(0)+1
              IF (IDQ(JFLIP).LT.0) THEN
                IDQ(JFLIP)=-IDQ(JFLIP)
                NDIQ=NDIQ+1
              ENDIF
              IF (NDIQ.LE.1) GOTO 730
            ENDIF
C...Place selected quarks first in IDQ, ordered in flavour.
            DO 740 JDQ=1,3
              IF (IDQ(JDQ).LE.0) THEN
                ITEMP1  = IDQ(JDQ)
                IDQ(JDQ)= IDQ(3)
                IDQ(3)  = -ITEMP1
                IF (IABS(K(IDQ(1),2)).LT.IABS(K(IDQ(2),2))) THEN
                  ITEMP1  = IDQ(1)
                  IDQ(1)  = IDQ(2)
                  IDQ(2)  = ITEMP1
                ENDIF
              ENDIF
  740       CONTINUE
C...Choose diquark spin.
            IF (NBRVQ.EQ.2) THEN
C...If the selected quarks are both valence, we may use SU(6) rules
C...to figure out which spin the diquark has, by a subdivision of the
C...original beam hadron into the selected diquark system plus a kicked
C...out quark, IKO.
              JKO=6
              DO 760 JDQ=1,2
                DO 750 JV=1,3
                  IF (IDQ(JDQ).EQ.IV(JS,JV)) JKO=JKO-JV
  750           CONTINUE
  760         CONTINUE
              IKO=IV(JS,JKO)
              CALL PYSPLI(MINT(10+JS),K(IKO,2),KFDUM,KFDQ)
            ELSE
C...If one or more of the selected quarks are not valence, we cannot use
C...SU(6) subdivisions of the original beam hadron. Instead, with the
C...flavours of the diquark already selected, we assume for now
C...50:50 spin-1:spin-0 (where spin-0 possible).
              KFDQ=1000*K(IDQ(1),2)+100*K(IDQ(2),2)
              IS=3
              IF (K(IDQ(1),2).NE.K(IDQ(2),2).AND.
     &           (1D0+3D0*PARJ(4))*PYR(0).LT.1D0) IS=1
              KFDQ=KFDQ+ISIGN(IS,KFDQ)
            ENDIF
 
C...Collapse diquark-j-quark system to baryon, if allowed and possible.
C...Note: third quark can per definition not also be valence,
C...therefore we can only do this if we are allowed to use sea quarks.
  770       IF (IDQ(3).NE.0.AND.MSTP(88).GE.2) THEN
              NTRY=0
  780         NTRY=NTRY+1
              CALL PYKFDI(KFDQ,K(IABS(IDQ(3)),2),KFDUM,KFBAR)
              IF (KFBAR.EQ.0.AND.NTRY.LE.100) THEN
                GOTO 780
              ELSEIF(NTRY.GT.100) THEN
C...If no baryon can be found, give up and form diquark.
                IDQ(3)=0
                GOTO 770
              ELSE
C...Replace junction by baryon.
                K(IJU,1)=1
                K(IJU,2)=KFBAR
                K(IJU,3)=MINT(83)+JS
                K(IJU,4)=0
                K(IJU,5)=0
                P(IJU,5)=PYMASS(KFBAR)
                DO 790 MSJ=1,3
C...Prepare removal of participating quarks from ER.
                  K(JST(JS,MSJ),1)=-1
  790           CONTINUE
              ENDIF
            ELSE
C...If collapse to baryon not possible or not allowed, replace junction
C...by diquark. This way, collapsed gluons that were pointing at the
C...junction will now point (correctly) at diquark.
              MANTI=ITJUNC(JS)-1
              K(IJU,1)=3
              K(IJU,2)=KFDQ
              K(IJU,3)=MINT(83)+JS
              K(IJU,4)=0
              K(IJU,5)=0
              DO 800 MSJ=1,3
                IP=JST(JS,MSJ)
                IF (IP.NE.IDQ(1).AND.IP.NE.IDQ(2)) THEN
                  K(IJU,4+MANTI)=0
                  K(IJU,5-MANTI)=IP*MSTU(5)
                  K(IP,4+MANTI)=MOD(K(IP,4+MANTI),MSTU(5))+
     &                 MSTU(5)*IJU
                  MCT(IJU,2-MANTI)=MCT(IP,1+MANTI)
                ELSE
C...Prepare removal of participating quarks from ER.
                  K(IP,1)=-1
                ENDIF
  800         CONTINUE
            ENDIF
 
C...Update so ER pointers to collapsed quarks
C...now go to collapsed object.
            DO 820 I=MINT(84)+1,N
              IF ((K(I,3).EQ.MINT(83)+JS.OR.K(I,3).EQ.MINT(83)+2+JS).AND
     &             .K(I,1).GT.0) THEN
                DO 810 ISID=4,5
                  IMO=K(I,ISID)/MSTU(5)
                  IDA=MOD(K(I,ISID),MSTU(5))
                  IF (IMO.GT.0) THEN
                    IF (K(IMO,1).EQ.-1) IMO=IJU
                  ENDIF
                  IF (IDA.GT.0) THEN
                    IF (K(IDA,1).EQ.-1) IDA=IJU
                  ENDIF
                  K(I,ISID)=IDA+MSTU(5)*IMO
  810           CONTINUE
              ENDIF
  820       CONTINUE
          ENDIF
        ENDIF
 
C...Finally, if beam remnant is empty, insert a gluon in beam remnant.
C...(this only happens for baryons, where we want to force the gluon
C...to sit next to the junction. Mesons handled above.)
        IF (NBRTOT(JS).EQ.0) THEN
          N=N+1
          DO 830 IX=1,5
            K(N,IX)=0
            P(N,IX)=0D0
            V(N,IX)=0D0
  830     CONTINUE
          IGL=N
          K(IGL,1)=3
          K(IGL,2)=21
          K(IGL,3)=MINT(83)+JS
          IF (ITJUNC(JS).NE.0) THEN
C...Incoming baryons. Pick random leg in JST (NVSUM = 3 for baryons)
            JLEG=PYR(0)*NVSUM(JS)+1
            I1=JST(JS,JLEG)
            JST(JS,JLEG)=IGL
            JCT=MCT(I1,ITJUNC(JS))
            MCT(IGL,3-ITJUNC(JS))=JCT
            NCT=NCT+1
            MCT(IGL,ITJUNC(JS))=NCT
            MANTI=ITJUNC(JS)-1
          ELSE
C...Meson. Should not happen.
            CALL PYERRM(19,'(PYMIHK:) Empty meson beam remnant')
            IF(NERRPR.LT.5) THEN
              WRITE(MSTU(11),*) 'This should not have been possible!'
              CALL PYLIST(4)
              NERRPR=NERRPR+1
            ENDIF
            MINT(51)=1
            RETURN
          ENDIF
          I2=MOD(K(I1,4+MANTI)/MSTU(5),MSTU(5))
          K(I1,4+MANTI)=MOD(K(I1,4+MANTI),MSTU(5))+MSTU(5)*IGL
          K(IGL,5-MANTI)=MOD(K(IGL,5-MANTI),MSTU(5))+MSTU(5)*I1
          K(IGL,4+MANTI)=MOD(K(IGL,4+MANTI),MSTU(5))+MSTU(5)*I2
          IF (K(I2,2).NE.88) THEN
            K(I2,5-MANTI)=MOD(K(I2,5-MANTI),MSTU(5))+MSTU(5)*IGL
          ELSE
            IF (MOD(K(I2,4),MSTU(5)).EQ.I1) THEN
              K(I2,4)=(K(I2,4)/MSTU(5))*MSTU(5)+IGL
            ELSEIF(MOD(K(I2,5)/MSTU(5),MSTU(5)).EQ.I1) THEN
              K(I2,5)=MOD(K(I2,5),MSTU(5))+MSTU(5)*IGL
            ELSE
              K(I2,5)=(K(I2,5)/MSTU(5))*MSTU(5)+IGL
            ENDIF
          ENDIF
        ENDIF
  840 CONTINUE
 
C...Remove collapsed quarks and junctions from ER and update IMI.
      CALL PYEDIT(11)
 
C...Also update beam remnant part of IMI.
      NMI(1)=MINT(31)
      NMI(2)=MINT(31)
      DO 850 I=MINT(53)+1,N
        IF (K(I,1).LE.0) GOTO 850
C...Restore BR quark/diquark/baryon pointers in IMI.
        IF ((K(I,2).NE.21.OR.K(I,1).NE.14).AND.K(I,2).NE.88) THEN
          JS=K(I,3)-MINT(83)
          NMI(JS)=NMI(JS)+1
          IMI(JS,NMI(JS),1)=I
          IMI(JS,NMI(JS),2)=0
        ENDIF
  850 CONTINUE
 
C...Restore companion information from collapsed gluons.
      DO 870 I=MINT(53)+1,N
        IF (K(I,2).EQ.21.AND.K(I,1).EQ.14) THEN
          JS=K(I,3)-MINT(83)
          JCD=MOD(K(I,4),MSTU(5))
          JAD=MOD(K(I,5),MSTU(5))
          DO 860 IM=1,NMI(JS)
            IF (IMI(JS,IM,1).EQ.JCD) IMC=IM
            IF (IMI(JS,IM,1).EQ.JAD) IMA=IM
  860     CONTINUE
          IMI(JS,IMC,2)=IMI(JS,IMA,1)
          IMI(JS,IMA,2)=IMI(JS,IMC,1)
        ENDIF
  870 CONTINUE
 
C...Renumber colour lines (since some have disappeared)
      JCT=0
      JCD=0
  880 JCT=JCT+1
      MFOUND=0
      I=MINT(84)
  890 I=I+1
      IF (I.EQ.N+1) THEN
        IF (MFOUND.EQ.0) JCD=JCD+1
      ELSEIF (MCT(I,1).EQ.JCT.AND.K(I,1).GE.1) THEN
        MCT(I,1)=JCT-JCD
        MFOUND=1
      ELSEIF (MCT(I,2).EQ.JCT.AND.K(I,1).GE.1) THEN
        MCT(I,2)=JCT-JCD
        MFOUND=1
      ENDIF
      IF (I.LE.N) GOTO 890
      IF (JCT.LT.NCT) GOTO 880
      NCT=JCT-JCD
 
C...Reset hard interaction subsystems to their CM frames.
      IF (IBOOST.EQ.1) THEN
        DO 900 IM=1,MINT(31)
          BETA=-(XMI(1,IM)-XMI(2,IM))/(XMI(1,IM)+XMI(2,IM))
          CALL PYROBO(IMISEP(IM-1)+1,IMISEP(IM),0D0,0D0,0D0,0D0,BETA)
  900   CONTINUE
C...Zero beam remnant longitudinal momenta and energies
        DO 910 I=MINT(53)+1,N
          P(I,3)=0D0
          P(I,4)=0D0
  910   CONTINUE
      ELSE
        CALL PYERRM(9
     &       ,'(PYMIHK:) Inconsistent kinematics. Too many boosts.')
C...Kill event and start another.
        MINT(51)=1
        RETURN
      ENDIF
 
 9999 RETURN
      END
