 
C***********************************************************************
 
C...PYRECO
C...Handles the possibility of colour reconnection in W+W- events,
C...Based on the main scenarios of the Sjostrand and Khoze study:
C...I, II, II', intermediate and instantaneous; plus one model
C...along the lines of the Gustafson and Hakkinen: GH.
C...Note: also handles Z0 Z0 and W-W+ events, but notation below
C...is as if first resonance is W+ and second W-.
 
      SUBROUTINE PYRECO(IW1,IW2,NSD1,NAFT1)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter value; number of points in MC integration.
      PARAMETER (NPT=100)
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/
C...Local arrays.
      DIMENSION NBEG(2),NEND(2),INP(50),INM(50),BEWW(3),XP(3),XM(3),
     &V1(3),V2(3),BETP(50,4),DIRP(50,3),BETM(50,4),DIRM(50,3),
     &XD(4),XB(4),IAP(NPT),IAM(NPT),WTA(NPT),V1P(3),V2P(3),V1M(3),
     &V2M(3),Q(4,3),XPP(3),XMM(3),IPC(20),IMC(20),TC(0:20),TPC(20),
     &TMC(20),IJOIN(100)
 
C...Functions to give four-product and to do determinants.
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)
      DETER(I,J,L)=Q(I,1)*Q(J,2)*Q(L,3)-Q(I,1)*Q(L,2)*Q(J,3)+
     &Q(J,1)*Q(L,2)*Q(I,3)-Q(J,1)*Q(I,2)*Q(L,3)+
     &Q(L,1)*Q(I,2)*Q(J,3)-Q(L,1)*Q(J,2)*Q(I,3)
 
C...Only allow fraction of recoupling for GH, intermediate and
C...instantaneous.
      IF(MSTP(115).EQ.5.OR.MSTP(115).EQ.11.OR.MSTP(115).EQ.12) THEN
        IF(PYR(0).GT.PARP(120)) RETURN
      ENDIF
      ISUB=MINT(1)
 
C...Common part for scenarios I, II, II', and GH.
      IF(MSTP(115).EQ.1.OR.MSTP(115).EQ.2.OR.MSTP(115).EQ.3.OR.
     &MSTP(115).EQ.5) THEN
 
C...Read out frequently-used parameters.
        PI=PARU(1)
        HBAR=PARU(3)
        PMW=PMAS(24,1)
        IF(ISUB.EQ.22) PMW=PMAS(23,1)
        PGW=PMAS(24,2)
        IF(ISUB.EQ.22) PGW=PMAS(23,2)
        TFRAG=PARP(115)
        RHAD=PARP(116)
        FACT=PARP(117)
        BLOWR=PARP(118)
        BLOWT=PARP(119)
 
C...Find range of decay products of the W's.
C...Background: the W's are stored in IW1 and IW2.
C...Their direct decay products in NSD1+1 through NSD1+4.
C...Products after shower (if any) in NSD1+5 through NAFT1
C...for first W and in NAFT1+1 through N for the second.
        IF(NAFT1.GT.NSD1+4) THEN
          NBEG(1)=NSD1+5
          NEND(1)=NAFT1
        ELSE
          NBEG(1)=NSD1+1
          NEND(1)=NSD1+2
        ENDIF
        IF(N.GT.NAFT1) THEN
          NBEG(2)=NAFT1+1
          NEND(2)=N
        ELSE
          NBEG(2)=NSD1+3
          NEND(2)=NSD1+4
        ENDIF
 
C...Rearrange parton shower products along strings.
        NOLD=N
        CALL PYPREP(NSD1+1)
        IF(MINT(51).NE.0) RETURN
 
C...Find partons pointing back to W+ and W-; store them with quark
C...end of string first.
        NNP=0
        NNM=0
        ISGP=0
        ISGM=0
        DO 120 I=NOLD+1,N
          IF(K(I,1).NE.1.AND.K(I,1).NE.2) GOTO 120
          IF(IABS(K(I,2)).GE.22) GOTO 120
          IF(K(I,3).GE.NBEG(1).AND.K(I,3).LE.NEND(1)) THEN
            IF(ISGP.EQ.0) ISGP=ISIGN(1,K(I,2))
            NNP=NNP+1
            IF(ISGP.EQ.1) THEN
              INP(NNP)=I
            ELSE
              DO 100 I1=NNP,2,-1
                INP(I1)=INP(I1-1)
  100         CONTINUE
              INP(1)=I
            ENDIF
            IF(K(I,1).EQ.1) ISGP=0
          ELSEIF(K(I,3).GE.NBEG(2).AND.K(I,3).LE.NEND(2)) THEN
            IF(ISGM.EQ.0) ISGM=ISIGN(1,K(I,2))
            NNM=NNM+1
            IF(ISGM.EQ.1) THEN
              INM(NNM)=I
            ELSE
              DO 110 I1=NNM,2,-1
                INM(I1)=INM(I1-1)
  110         CONTINUE
              INM(1)=I
            ENDIF
            IF(K(I,1).EQ.1) ISGM=0
          ENDIF
  120   CONTINUE
 
C...Boost to W+W- rest frame (not strictly needed).
        DO 130 J=1,3
          BEWW(J)=(P(IW1,J)+P(IW2,J))/(P(IW1,4)+P(IW2,4))
  130   CONTINUE
        CALL PYROBO(IW1,IW1,0D0,0D0,-BEWW(1),-BEWW(2),-BEWW(3))
        CALL PYROBO(IW2,IW2,0D0,0D0,-BEWW(1),-BEWW(2),-BEWW(3))
        CALL PYROBO(NOLD+1,N,0D0,0D0,-BEWW(1),-BEWW(2),-BEWW(3))
 
C...Select decay vertices of W+ and W-.
        TP=HBAR*(-LOG(PYR(0)))*P(IW1,4)/
     &  SQRT((P(IW1,5)**2-PMW**2)**2+(P(IW1,5)**2*PGW/PMW)**2)
        TM=HBAR*(-LOG(PYR(0)))*P(IW2,4)/
     &  SQRT((P(IW2,5)**2-PMW**2)**2+(P(IW2,5)**2*PGW/PMW)**2)
        GTMAX=MAX(TP,TM)
        DO 140 J=1,3
          XP(J)=TP*P(IW1,J)/P(IW1,4)
          XM(J)=TM*P(IW2,J)/P(IW2,4)
  140   CONTINUE
 
C...Begin scenario I specifics.
        IF(MSTP(115).EQ.1) THEN
 
C...Reconstruct velocity and direction of W+ string pieces.
          DO 170 IIP=1,NNP-1
            IF(K(INP(IIP),2).LT.0) GOTO 170
            I1=INP(IIP)
            I2=INP(IIP+1)
            P1A=SQRT(P(I1,1)**2+P(I1,2)**2+P(I1,3)**2)
            P2A=SQRT(P(I2,1)**2+P(I2,2)**2+P(I2,3)**2)
            DO 150 J=1,3
              V1(J)=P(I1,J)/P1A
              V2(J)=P(I2,J)/P2A
              BETP(IIP,J)=0.5D0*(V1(J)+V2(J))
              DIRP(IIP,J)=V1(J)-V2(J)
  150       CONTINUE
            BETP(IIP,4)=1D0/SQRT(1D0-BETP(IIP,1)**2-BETP(IIP,2)**2-
     &      BETP(IIP,3)**2)
            DIRL=SQRT(DIRP(IIP,1)**2+DIRP(IIP,2)**2+DIRP(IIP,3)**2)
            DO 160 J=1,3
              DIRP(IIP,J)=DIRP(IIP,J)/DIRL
  160       CONTINUE
  170     CONTINUE
 
C...Reconstruct velocity and direction of W- string pieces.
          DO 200 IIM=1,NNM-1
            IF(K(INM(IIM),2).LT.0) GOTO 200
            I1=INM(IIM)
            I2=INM(IIM+1)
            P1A=SQRT(P(I1,1)**2+P(I1,2)**2+P(I1,3)**2)
            P2A=SQRT(P(I2,1)**2+P(I2,2)**2+P(I2,3)**2)
            DO 180 J=1,3
              V1(J)=P(I1,J)/P1A
              V2(J)=P(I2,J)/P2A
              BETM(IIM,J)=0.5D0*(V1(J)+V2(J))
              DIRM(IIM,J)=V1(J)-V2(J)
  180       CONTINUE
            BETM(IIM,4)=1D0/SQRT(1D0-BETM(IIM,1)**2-BETM(IIM,2)**2-
     &      BETM(IIM,3)**2)
            DIRL=SQRT(DIRM(IIM,1)**2+DIRM(IIM,2)**2+DIRM(IIM,3)**2)
            DO 190 J=1,3
              DIRM(IIM,J)=DIRM(IIM,J)/DIRL
  190       CONTINUE
  200     CONTINUE
 
C...Loop over number of space-time points.
          NACC=0
          SUM=0D0
          DO 250 IPT=1,NPT
 
C...Pick x,y,z,t Gaussian (width RHAD and TFRAG, respectively).
            R=SQRT(-LOG(PYR(0)))
            PHI=2D0*PI*PYR(0)
            X=BLOWR*RHAD*R*COS(PHI)
            Y=BLOWR*RHAD*R*SIN(PHI)
            R=SQRT(-LOG(PYR(0)))
            PHI=2D0*PI*PYR(0)
            Z=BLOWR*RHAD*R*COS(PHI)
            T=GTMAX+BLOWT*SQRT(0.5D0)*TFRAG*R*ABS(SIN(PHI))
 
C...Reject impossible points. Weight for sample distribution.
            IF(T**2-X**2-Y**2-Z**2.LT.0D0) GOTO 250
            WTSMP=EXP(-(X**2+Y**2+Z**2)/(BLOWR*RHAD)**2)*
     &      EXP(-2D0*(T-GTMAX)**2/(BLOWT*TFRAG)**2)
 
C...Loop over W+ string pieces and find one with largest weight.
            IMAXP=0
            WTMAXP=1D-10
            XD(1)=X-XP(1)
            XD(2)=Y-XP(2)
            XD(3)=Z-XP(3)
            XD(4)=T-TP
            DO 220 IIP=1,NNP-1
              IF(K(INP(IIP),2).LT.0) GOTO 220
              BED=BETP(IIP,1)*XD(1)+BETP(IIP,2)*XD(2)+BETP(IIP,3)*XD(3)
              BEDG=BETP(IIP,4)*(BETP(IIP,4)*BED/(1D0+BETP(IIP,4))-XD(4))
              DO 210 J=1,3
                XB(J)=XD(J)+BEDG*BETP(IIP,J)
  210         CONTINUE
              XB(4)=BETP(IIP,4)*(XD(4)-BED)
              SR2=XB(1)**2+XB(2)**2+XB(3)**2
              SZ2=(DIRP(IIP,1)*XB(1)+DIRP(IIP,2)*XB(2)+
     &        DIRP(IIP,3)*XB(3))**2
              WTP=EXP(-(SR2-SZ2)/(2D0*RHAD**2))*EXP(-(XB(4)**2-SZ2)/
     &        TFRAG**2)
              IF(XB(4)-SQRT(SR2).LT.0D0) WTP=0D0
              IF(WTP.GT.WTMAXP) THEN
                IMAXP=IIP
                WTMAXP=WTP
              ENDIF
  220       CONTINUE
 
C...Loop over W- string pieces and find one with largest weight.
            IMAXM=0
            WTMAXM=1D-10
            XD(1)=X-XM(1)
            XD(2)=Y-XM(2)
            XD(3)=Z-XM(3)
            XD(4)=T-TM
            DO 240 IIM=1,NNM-1
              IF(K(INM(IIM),2).LT.0) GOTO 240
              BED=BETM(IIM,1)*XD(1)+BETM(IIM,2)*XD(2)+BETM(IIM,3)*XD(3)
              BEDG=BETM(IIM,4)*(BETM(IIM,4)*BED/(1D0+BETM(IIM,4))-XD(4))
              DO 230 J=1,3
                XB(J)=XD(J)+BEDG*BETM(IIM,J)
  230         CONTINUE
              XB(4)=BETM(IIM,4)*(XD(4)-BED)
              SR2=XB(1)**2+XB(2)**2+XB(3)**2
              SZ2=(DIRM(IIM,1)*XB(1)+DIRM(IIM,2)*XB(2)+
     &        DIRM(IIM,3)*XB(3))**2
              WTM=EXP(-(SR2-SZ2)/(2D0*RHAD**2))*EXP(-(XB(4)**2-SZ2)/
     &        TFRAG**2)
              IF(XB(4)-SQRT(SR2).LT.0D0) WTM=0D0
              IF(WTM.GT.WTMAXM) THEN
                IMAXM=IIM
                WTMAXM=WTM
              ENDIF
  240       CONTINUE
 
C...Result of integration.
            WT=0D0
            IF(IMAXP.NE.0.AND.IMAXM.NE.0) THEN
              WT=WTMAXP*WTMAXM/WTSMP
              SUM=SUM+WT
              NACC=NACC+1
              IAP(NACC)=IMAXP
              IAM(NACC)=IMAXM
              WTA(NACC)=WT
            ENDIF
  250     CONTINUE
          RES=BLOWR**3*BLOWT*SUM/NPT
 
C...Decide whether to reconnect and, if so, where.
          IACC=0
          PREC=1D0-EXP(-FACT*RES)
          IF(PREC.GT.PYR(0)) THEN
            RSUM=PYR(0)*SUM
            DO 260 IA=1,NACC
              IACC=IA
              RSUM=RSUM-WTA(IA)
              IF(RSUM.LE.0D0) GOTO 270
  260       CONTINUE
  270       IIP=IAP(IACC)
            IIM=IAM(IACC)
          ENDIF
 
C...Begin scenario II and II' specifics.
        ELSEIF(MSTP(115).EQ.2.OR.MSTP(115).EQ.3) THEN
 
C...Loop through all string pieces, one from W+ and one from W-.
          NCROSS=0
          TC(0)=0D0
          DO 340 IIP=1,NNP-1
            IF(K(INP(IIP),2).LT.0) GOTO 340
            I1P=INP(IIP)
            I2P=INP(IIP+1)
            DO 330 IIM=1,NNM-1
              IF(K(INM(IIM),2).LT.0) GOTO 330
              I1M=INM(IIM)
              I2M=INM(IIM+1)
 
C...Find endpoint velocity vectors.
              DO 280 J=1,3
                V1P(J)=P(I1P,J)/P(I1P,4)
                V2P(J)=P(I2P,J)/P(I2P,4)
                V1M(J)=P(I1M,J)/P(I1M,4)
                V2M(J)=P(I2M,J)/P(I2M,4)
  280         CONTINUE
 
C...Define q matrix and find t.
              DO 290 J=1,3
                Q(1,J)=V2P(J)-V1P(J)
                Q(2,J)=-(V2M(J)-V1M(J))
                Q(3,J)=XP(J)-XM(J)-TP*V1P(J)+TM*V1M(J)
                Q(4,J)=V1P(J)-V1M(J)
  290         CONTINUE
              T=-DETER(1,2,3)/DETER(1,2,4)
 
C...Find alpha and beta; i.e. coordinates of crossing point.
              S11=Q(1,1)*(T-TP)
              S12=Q(2,1)*(T-TM)
              S13=Q(3,1)+Q(4,1)*T
              S21=Q(1,2)*(T-TP)
              S22=Q(2,2)*(T-TM)
              S23=Q(3,2)+Q(4,2)*T
              DEN=S11*S22-S12*S21
              ALP=(S12*S23-S22*S13)/DEN
              BET=(S21*S13-S11*S23)/DEN
 
C...Check if solution acceptable.
              IANSW=1
              IF(T.LT.GTMAX) IANSW=0
              IF(ALP.LT.0D0.OR.ALP.GT.1D0) IANSW=0
              IF(BET.LT.0D0.OR.BET.GT.1D0) IANSW=0
 
C...Find point of crossing and check that not inconsistent.
              DO 300 J=1,3
                XPP(J)=XP(J)+(V1P(J)+ALP*(V2P(J)-V1P(J)))*(T-TP)
                XMM(J)=XM(J)+(V1M(J)+BET*(V2M(J)-V1M(J)))*(T-TM)
  300         CONTINUE
              D2PM=(XPP(1)-XMM(1))**2+(XPP(2)-XMM(2))**2+
     &        (XPP(3)-XMM(3))**2
              D2P=XPP(1)**2+XPP(2)**2+XPP(3)**2
              D2M=XMM(1)**2+XMM(2)**2+XMM(3)**2
              IF(D2PM.GT.1D-4*(D2P+D2M)) IANSW=-1
 
C...Find string eigentimes at crossing.
              IF(IANSW.EQ.1) THEN
                TAUP=SQRT(MAX(0D0,(T-TP)**2-(XPP(1)-XP(1))**2-
     &          (XPP(2)-XP(2))**2-(XPP(3)-XP(3))**2))
                TAUM=SQRT(MAX(0D0,(T-TM)**2-(XMM(1)-XM(1))**2-
     &          (XMM(2)-XM(2))**2-(XMM(3)-XM(3))**2))
              ELSE
                TAUP=0D0
                TAUM=0D0
              ENDIF
 
C...Order crossings by time. End loop over crossings.
              IF(IANSW.EQ.1.AND.NCROSS.LT.20) THEN
                NCROSS=NCROSS+1
                DO 310 I1=NCROSS,1,-1
                  IF(T.GT.TC(I1-1).OR.I1.EQ.1) THEN
                    IPC(I1)=IIP
                    IMC(I1)=IIM
                    TC(I1)=T
                    TPC(I1)=TAUP
                    TMC(I1)=TAUM
                    GOTO 320
                  ELSE
                    IPC(I1)=IPC(I1-1)
                    IMC(I1)=IMC(I1-1)
                    TC(I1)=TC(I1-1)
                    TPC(I1)=TPC(I1-1)
                    TMC(I1)=TMC(I1-1)
                  ENDIF
  310           CONTINUE
  320           CONTINUE
              ENDIF
  330       CONTINUE
  340     CONTINUE
 
C...Loop over crossings; find first (if any) acceptable one.
          IACC=0
          IF(NCROSS.GE.1) THEN
            DO 350 IC=1,NCROSS
              PNFRAG=EXP(-(TPC(IC)**2+TMC(IC)**2)/TFRAG**2)
              IF(PNFRAG.GT.PYR(0)) THEN
C...Scenario II: only compare with fragmentation time.
                IF(MSTP(115).EQ.2) THEN
                  IACC=IC
                  IIP=IPC(IACC)
                  IIM=IMC(IACC)
                  GOTO 360
C...Scenario II': also require that string length decreases.
                ELSE
                  IIP=IPC(IC)
                  IIM=IMC(IC)
                  I1P=INP(IIP)
                  I2P=INP(IIP+1)
                  I1M=INM(IIM)
                  I2M=INM(IIM+1)
                  ELOLD=FOUR(I1P,I2P)*FOUR(I1M,I2M)
                  ELNEW=FOUR(I1P,I2M)*FOUR(I1M,I2P)
                  IF(ELNEW.LT.ELOLD) THEN
                    IACC=IC
                    IIP=IPC(IACC)
                    IIM=IMC(IACC)
                    GOTO 360
                  ENDIF
                ENDIF
              ENDIF
  350       CONTINUE
  360       CONTINUE
          ENDIF
 
C...Begin scenario GH specifics.
        ELSEIF(MSTP(115).EQ.5) THEN
 
C...Loop through all string pieces, one from W+ and one from W-.
          IACC=0
          ELMIN=1D0
          DO 380 IIP=1,NNP-1
            IF(K(INP(IIP),2).LT.0) GOTO 380
            I1P=INP(IIP)
            I2P=INP(IIP+1)
            DO 370 IIM=1,NNM-1
              IF(K(INM(IIM),2).LT.0) GOTO 370
              I1M=INM(IIM)
              I2M=INM(IIM+1)
 
C...Look for largest decrease of (exponent of) Lambda measure.
              ELOLD=FOUR(I1P,I2P)*FOUR(I1M,I2M)
              ELNEW=FOUR(I1P,I2M)*FOUR(I1M,I2P)
              ELDIF=ELNEW/MAX(1D-10,ELOLD)
              IF(ELDIF.LT.ELMIN) THEN
                IACC=IIP+IIM
                ELMIN=ELDIF
                IPC(1)=IIP
                IMC(1)=IIM
              ENDIF
  370       CONTINUE
  380     CONTINUE
          IIP=IPC(1)
          IIM=IMC(1)
        ENDIF
 
C...Common for scenarios I, II, II' and GH: reconnect strings.
        IF(IACC.NE.0) THEN
          MINT(32)=1
          NJOIN=0
          DO 390 IS=1,NNP+NNM
            NJOIN=NJOIN+1
            IF(IS.LE.IIP) THEN
              I=INP(IS)
            ELSEIF(IS.LE.IIP+NNM-IIM) THEN
              I=INM(IS-IIP+IIM)
            ELSEIF(IS.LE.IIP+NNM) THEN
              I=INM(IS-IIP-NNM+IIM)
            ELSE
              I=INP(IS-NNM)
            ENDIF
            IJOIN(NJOIN)=I
            IF(K(I,2).LT.0) THEN
              CALL PYJOIN(NJOIN,IJOIN)
              NJOIN=0
            ENDIF
  390     CONTINUE
 
C...Restore original event record if no reconnection.
        ELSE
          DO 400 I=NSD1+1,NOLD
            IF(K(I,1).EQ.13.OR.K(I,1).EQ.14) THEN
              K(I,4)=MOD(K(I,4),MSTU(5)**2)
              K(I,5)=MOD(K(I,5),MSTU(5)**2)
            ENDIF
  400     CONTINUE
          DO 410 I=NOLD+1,N
            K(K(I,3),1)=3
  410     CONTINUE
          N=NOLD
        ENDIF
 
C...Boost back system.
        CALL PYROBO(IW1,IW1,0D0,0D0,BEWW(1),BEWW(2),BEWW(3))
        CALL PYROBO(IW2,IW2,0D0,0D0,BEWW(1),BEWW(2),BEWW(3))
        IF(N.GT.NOLD) CALL PYROBO(NOLD+1,N,0D0,0D0,
     &  BEWW(1),BEWW(2),BEWW(3))
 
C...Common part for intermediate and instantaneous scenarios.
      ELSEIF(MSTP(115).EQ.11.OR.MSTP(115).EQ.12) THEN
        MINT(32)=1
 
C...Remove old shower products and reset showering ones.
        N=NSD1+4
        DO 420 I=NSD1+1,NSD1+4
          K(I,1)=3
          K(I,4)=MOD(K(I,4),MSTU(5)**2)
          K(I,5)=MOD(K(I,5),MSTU(5)**2)
  420   CONTINUE
 
C...Identify quark-antiquark pairs.
        IQ1=NSD1+1
        IQ2=NSD1+2
        IQ3=NSD1+3
        IF(K(IQ1,2)*K(IQ3,2).LT.0) IQ3=NSD1+4
        IQ4=2*NSD1+7-IQ3
 
C...Reconnect strings.
        IJOIN(1)=IQ1
        IJOIN(2)=IQ4
        CALL PYJOIN(2,IJOIN)
        IJOIN(1)=IQ3
        IJOIN(2)=IQ2
        CALL PYJOIN(2,IJOIN)
 
C...Do new parton showers in intermediate scenario.
        IF(MSTP(71).GE.1.AND.MSTP(115).EQ.11) THEN
          MSTJ50=MSTJ(50)
          MSTJ(50)=0
          CALL PYSHOW(IQ1,IQ2,P(IW1,5))
          CALL PYSHOW(IQ3,IQ4,P(IW2,5))
          MSTJ(50)=MSTJ50
 
C...Do new parton showers in instantaneous scenario.
        ELSEIF(MSTP(71).GE.1.AND.MSTP(115).EQ.12) THEN
          PPM2=(P(IQ1,4)+P(IQ4,4))**2-(P(IQ1,1)+P(IQ4,1))**2-
     &    (P(IQ1,2)+P(IQ4,2))**2-(P(IQ1,3)+P(IQ4,3))**2
          PPM=SQRT(MAX(0D0,PPM2))
          CALL PYSHOW(IQ1,IQ4,PPM)
          PPM2=(P(IQ3,4)+P(IQ2,4))**2-(P(IQ3,1)+P(IQ2,1))**2-
     &    (P(IQ3,2)+P(IQ2,2))**2-(P(IQ3,3)+P(IQ2,3))**2
          PPM=SQRT(MAX(0D0,PPM2))
          CALL PYSHOW(IQ3,IQ2,PPM)
        ENDIF
      ENDIF
 
      RETURN
      END
