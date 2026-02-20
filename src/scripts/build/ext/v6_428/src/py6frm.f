 
C*********************************************************************
 
C...PY6FRM
C...An interface from a six-fermion generator to include
C...parton showers and hadronization.
 
      SUBROUTINE PY6FRM(P12,P13,P21,P23,P31,P32,PTOP,IRAD,ITAU,ICOM)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYJETS/,/PYDAT1/
C...Local arrays.
      DIMENSION IJOIN(2),INTAU(6),BETA(3),BETAO(3),BETAN(3)
 
C...Call PYHEPC to convert input from HEPEVT to PYJETS common.
      IF(ICOM.EQ.0) THEN
        MSTU(28)=0
        CALL PYHEPC(2)
      ENDIF
 
C...Loop through entries and pick up all final fermions/antifermions.
      I1=0
      I2=0
      I3=0
      I4=0
      I5=0
      I6=0
      DO 100 I=1,N
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 100
      KFA=IABS(K(I,2))
      IF((KFA.GE.1.AND.KFA.LE.6).OR.(KFA.GE.11.AND.KFA.LE.16)) THEN
        IF(K(I,2).GT.0) THEN
          IF(I1.EQ.0) THEN
            I1=I
          ELSEIF(I3.EQ.0) THEN
            I3=I
          ELSEIF(I5.EQ.0) THEN
            I5=I
          ELSE
            CALL PYERRM(16,'(PY6FRM:) more than three fermions')
          ENDIF
        ELSE
          IF(I2.EQ.0) THEN
            I2=I
          ELSEIF(I4.EQ.0) THEN
            I4=I
          ELSEIF(I6.EQ.0) THEN
            I6=I
          ELSE
            CALL PYERRM(16,'(PY6FRM:) more than three antifermions')
          ENDIF
        ENDIF
      ENDIF
  100 CONTINUE
 
C...Check that event is arranged according to conventions.
      IF(I5.EQ.0.OR.I6.EQ.0) THEN
        CALL PYERRM(16,'(PY6FRM:) event contains too few fermions')
      ENDIF
      IF(I2.LT.I1.OR.I3.LT.I2.OR.I4.LT.I3.OR.I5.LT.I4.OR.I6.LT.I5) THEN
        CALL PYERRM(6,'(PY6FRM:) fermions arranged in wrong order')
      ENDIF
 
C...Check which fermion pairs are quarks and which leptons.
      IF(IABS(K(I1,2)).LT.10.AND.IABS(K(I2,2)).LT.10) THEN
        IQL12=1
      ELSEIF(IABS(K(I1,2)).GT.10.AND.IABS(K(I2,2)).GT.10) THEN
        IQL12=2
      ELSE
        CALL PYERRM(16,'(PY6FRM:) first fermion pair inconsistent')
      ENDIF
      IF(IABS(K(I3,2)).LT.10.AND.IABS(K(I4,2)).LT.10) THEN
        IQL34=1
      ELSEIF(IABS(K(I3,2)).GT.10.AND.IABS(K(I4,2)).GT.10) THEN
        IQL34=2
      ELSE
        CALL PYERRM(16,'(PY6FRM:) second fermion pair inconsistent')
      ENDIF
      IF(IABS(K(I5,2)).LT.10.AND.IABS(K(I6,2)).LT.10) THEN
        IQL56=1
      ELSEIF(IABS(K(I5,2)).GT.10.AND.IABS(K(I6,2)).GT.10) THEN
        IQL56=2
      ELSE
        CALL PYERRM(16,'(PY6FRM:) third fermion pair inconsistent')
      ENDIF
 
C...Decide whether to allow or not photon radiation in showers.
      MSTJ(41)=2
      IF(IRAD.EQ.0) MSTJ(41)=1
 
C...Allow dipole pairings only among leptons and quarks separately.
      P12D=P12
      P13D=0D0
      IF(IQL34.EQ.IQL56) P13D=P13
      P21D=0D0
      IF(IQL12.EQ.IQL34) P21D=P21
      P23D=0D0
      IF(IQL12.EQ.IQL34.AND.IQL12.EQ.IQL56) P23D=P23
      P31D=0D0
      IF(IQL12.EQ.IQL34.AND.IQL12.EQ.IQL56) P31D=P31
      P32D=0D0
      IF(IQL12.EQ.IQL56) P32D=P32
 
C...Decide whether t+tbar.
      ITOP=0
      IF(PYR(0).LT.PTOP) THEN
        ITOP=1
 
C...If t+tbar: reconstruct t's.
        IT=N+1
        ITB=N+2
        DO 110 J=1,5
          K(IT,J)=0
          K(ITB,J)=0
          P(IT,J)=P(I1,J)+P(I3,J)+P(I4,J)
          P(ITB,J)=P(I2,J)+P(I5,J)+P(I6,J)
          V(IT,J)=0D0
          V(ITB,J)=0D0
  110   CONTINUE
        K(IT,1)=1
        K(ITB,1)=1
        K(IT,2)=6
        K(ITB,2)=-6
        P(IT,5)=SQRT(MAX(0D0,P(IT,4)**2-P(IT,1)**2-P(IT,2)**2-
     &  P(IT,3)**2))
        P(ITB,5)=SQRT(MAX(0D0,P(ITB,4)**2-P(ITB,1)**2-P(ITB,2)**2-
     &  P(ITB,3)**2))
        N=N+2
 
C...If t+tbar: colour join t's and let them shower.
        IJOIN(1)=IT
        IJOIN(2)=ITB
        CALL PYJOIN(2,IJOIN)
        PMTTS=(P(IT,4)+P(ITB,4))**2-(P(IT,1)+P(ITB,1))**2-
     &  (P(IT,2)+P(ITB,2))**2-(P(IT,3)+P(ITB,3))**2
        CALL PYSHOW(IT,ITB,SQRT(MAX(0D0,PMTTS)))
 
C...If t+tbar: pick up the t's after shower.
        ITNEW=IT
        ITBNEW=ITB
        DO 120 I=ITB+1,N
          IF(K(I,2).EQ.6) ITNEW=I
          IF(K(I,2).EQ.-6) ITBNEW=I
  120   CONTINUE
 
C...If t+tbar: loop over two top systems.
        DO 200 IT1=1,2
          IF(IT1.EQ.1) THEN
            ITO=IT
            ITN=ITNEW
            IBO=I1
            IW1=I3
            IW2=I4
          ELSE
            ITO=ITB
            ITN=ITBNEW
            IBO=I2
            IW1=I5
            IW2=I6
          ENDIF
          IF(IABS(K(IBO,2)).NE.5) CALL PYERRM(6,
     &    '(PY6FRM:) not b in t decay')
 
C...If t+tbar: find boost from original to new top frame.
          DO 130 J=1,3
            BETAO(J)=P(ITO,J)/P(ITO,4)
            BETAN(J)=P(ITN,J)/P(ITN,4)
  130     CONTINUE
 
C...If t+tbar: boost copy of b by t shower and connect it in colour.
          N=N+1
          IB=N
          K(IB,1)=3
          K(IB,2)=K(IBO,2)
          K(IB,3)=ITN
          DO 140 J=1,5
            P(IB,J)=P(IBO,J)
            V(IB,J)=0D0
  140     CONTINUE
          CALL PYROBO(IB,IB,0D0,0D0,-BETAO(1),-BETAO(2),-BETAO(3))
          CALL PYROBO(IB,IB,0D0,0D0,BETAN(1),BETAN(2),BETAN(3))
          K(IB,4)=MSTU(5)*ITN
          K(IB,5)=MSTU(5)*ITN
          K(ITN,4)=K(ITN,4)+IB
          K(ITN,5)=K(ITN,5)+IB
          K(ITN,1)=K(ITN,1)+10
          K(IBO,1)=K(IBO,1)+10
 
C...If t+tbar: construct W recoiling against b.
          N=N+1
          IW=N
          DO 150 J=1,5
            K(IW,J)=0
            V(IW,J)=0D0
  150     CONTINUE
          K(IW,1)=1
          KCHW=PYCHGE(K(IW1,2))+PYCHGE(K(IW2,2))
          IF(IABS(KCHW).EQ.3) THEN
            K(IW,2)=ISIGN(24,KCHW)
          ELSE
            CALL PYERRM(16,'(PY6FRM:) fermion pair inconsistent with W')
          ENDIF
          K(IW,3)=IW1
 
C...If t+tbar: construct W momentum, including boost by t shower.
          DO 160 J=1,4
            P(IW,J)=P(IW1,J)+P(IW2,J)
  160     CONTINUE
          P(IW,5)=SQRT(MAX(0D0,P(IW,4)**2-P(IW,1)**2-P(IW,2)**2-
     &    P(IW,3)**2))
          CALL PYROBO(IW,IW,0D0,0D0,-BETAO(1),-BETAO(2),-BETAO(3))
          CALL PYROBO(IW,IW,0D0,0D0,BETAN(1),BETAN(2),BETAN(3))
 
C...If t+tbar: boost b and W to top rest frame.
          DO 170 J=1,3
            BETA(J)=(P(IB,J)+P(IW,J))/(P(IB,4)+P(IW,4))
  170     CONTINUE
          CALL PYROBO(IB,IB,0D0,0D0,-BETA(1),-BETA(2),-BETA(3))
          CALL PYROBO(IW,IW,0D0,0D0,-BETA(1),-BETA(2),-BETA(3))
 
C...If t+tbar: let b shower and pick up modified W.
          PMTS=(P(IB,4)+P(IW,4))**2-(P(IB,1)+P(IW,1))**2-
     &    (P(IB,2)+P(IW,2))**2-(P(IB,3)+P(IW,3))**2
          CALL PYSHOW(IB,IW,SQRT(MAX(0D0,PMTS)))
          DO 180 I=IW,N
            IF(IABS(K(I,2)).EQ.24) IWM=I
  180     CONTINUE
 
C...If t+tbar: take copy of W decay products.
          DO 190 J=1,5
            K(N+1,J)=K(IW1,J)
            P(N+1,J)=P(IW1,J)
            V(N+1,J)=V(IW1,J)
            K(N+2,J)=K(IW2,J)
            P(N+2,J)=P(IW2,J)
            V(N+2,J)=V(IW2,J)
  190     CONTINUE
          K(IW1,1)=K(IW1,1)+10
          K(IW2,1)=K(IW2,1)+10
          K(IWM,1)=K(IWM,1)+10
          K(IWM,4)=N+1
          K(IWM,5)=N+2
          K(N+1,3)=IWM
          K(N+2,3)=IWM
          IF(IT1.EQ.1) THEN
            I3=N+1
            I4=N+2
          ELSE
            I5=N+1
            I6=N+2
          ENDIF
          N=N+2
 
C...If t+tbar: boost W decay products, first by effects of t shower,
C...then by those of b shower. b and its shower simple boost back.
          CALL PYROBO(N-1,N,0D0,0D0,-BETAO(1),-BETAO(2),-BETAO(3))
          CALL PYROBO(N-1,N,0D0,0D0,BETAN(1),BETAN(2),BETAN(3))
          CALL PYROBO(N-1,N,0D0,0D0,-BETA(1),-BETA(2),-BETA(3))
          CALL PYROBO(N-1,N,0D0,0D0,-P(IW,1)/P(IW,4),
     &    -P(IW,2)/P(IW,4),-P(IW,3)/P(IW,4))
          CALL PYROBO(N-1,N,0D0,0D0,P(IWM,1)/P(IWM,4),
     &    P(IWM,2)/P(IWM,4),P(IWM,3)/P(IWM,4))
          CALL PYROBO(IB,IB,0D0,0D0,BETA(1),BETA(2),BETA(3))
          CALL PYROBO(IW,N,0D0,0D0,BETA(1),BETA(2),BETA(3))
  200   CONTINUE
      ENDIF
 
C...Decide on dipole pairing.
      IP1=I1
      IP3=I3
      IP5=I5
      PRN=PYR(0)*(P12D+P13D+P21D+P23D+P31D+P32D)
      IF(ITOP.EQ.1.OR.PRN.LT.P12D) THEN
        IP2=I2
        IP4=I4
        IP6=I6
      ELSEIF(PRN.LT.P12D+P13D) THEN
        IP2=I2
        IP4=I6
        IP6=I4
      ELSEIF(PRN.LT.P12D+P13D+P21D) THEN
        IP2=I4
        IP4=I2
        IP6=I6
      ELSEIF(PRN.LT.P12D+P13D+P21D+P23D) THEN
        IP2=I4
        IP4=I6
        IP6=I2
      ELSEIF(PRN.LT.P12D+P13D+P21D+P23D+P31D) THEN
        IP2=I6
        IP4=I2
        IP6=I4
      ELSE
        IP2=I6
        IP4=I4
        IP6=I2
      ENDIF
 
C...Do colour joinings and parton showers
C...(except ones already made for t+tbar).
      IF(ITOP.EQ.0) THEN
        IF(IQL12.EQ.1) THEN
          IJOIN(1)=IP1
          IJOIN(2)=IP2
          CALL PYJOIN(2,IJOIN)
        ENDIF
        IF(IQL12.EQ.1.OR.IRAD.EQ.1) THEN
          PM12S=(P(IP1,4)+P(IP2,4))**2-(P(IP1,1)+P(IP2,1))**2-
     &    (P(IP1,2)+P(IP2,2))**2-(P(IP1,3)+P(IP2,3))**2
          CALL PYSHOW(IP1,IP2,SQRT(MAX(0D0,PM12S)))
        ENDIF
      ENDIF
      IF(IQL34.EQ.1) THEN
        IJOIN(1)=IP3
        IJOIN(2)=IP4
        CALL PYJOIN(2,IJOIN)
      ENDIF
      IF(IQL34.EQ.1.OR.IRAD.EQ.1) THEN
        PM34S=(P(IP3,4)+P(IP4,4))**2-(P(IP3,1)+P(IP4,1))**2-
     &  (P(IP3,2)+P(IP4,2))**2-(P(IP3,3)+P(IP4,3))**2
        CALL PYSHOW(IP3,IP4,SQRT(MAX(0D0,PM34S)))
      ENDIF
      IF(IQL56.EQ.1) THEN
        IJOIN(1)=IP5
        IJOIN(2)=IP6
        CALL PYJOIN(2,IJOIN)
      ENDIF
      IF(IQL56.EQ.1.OR.IRAD.EQ.1) THEN
        PM56S=(P(IP5,4)+P(IP6,4))**2-(P(IP5,1)+P(IP6,1))**2-
     &  (P(IP5,2)+P(IP6,2))**2-(P(IP5,3)+P(IP6,3))**2
        CALL PYSHOW(IP5,IP6,SQRT(MAX(0D0,PM56S)))
      ENDIF
 
C...Do fragmentation and decays. Possibly except tau decay.
      IF(ITAU.EQ.0) THEN
        NTAU=0
        DO 210 I=1,N
        IF(IABS(K(I,2)).EQ.15.AND.K(I,1).EQ.1) THEN
          NTAU=NTAU+1
          INTAU(NTAU)=I
          K(I,1)=11
        ENDIF
  210   CONTINUE
      ENDIF
      CALL PYEXEC
      IF(ITAU.EQ.0) THEN
        DO 220 I=1,NTAU
        K(INTAU(I),1)=1
  220   CONTINUE
      ENDIF
 
C...Call PYHEPC to convert output from PYJETS to HEPEVT common.
      IF(ICOM.EQ.0) THEN
        MSTU(28)=0
        CALL PYHEPC(1)
      ENDIF
 
      END
