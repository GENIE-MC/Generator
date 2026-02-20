 
C*********************************************************************
 
C...PY4JET
C...An interface from a four-parton generator to include
C...parton showers and hadronization.
 
      SUBROUTINE PY4JET(PMAX,IRAD,ICOM)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYJETS/,/PYDAT1/
C...Local arrays.
      DIMENSION IJOIN(2),PTOT(4),BETA(3)
 
C...Call PYHEPC to convert input from HEPEVT to PYJETS common.
      IF(ICOM.EQ.0) THEN
        MSTU(28)=0
        CALL PYHEPC(2)
      ENDIF
 
C...Loop through entries and pick up all final partons.
      I1=0
      I2=0
      I3=0
      I4=0
      DO 100 I=1,N
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 100
      KFA=IABS(K(I,2))
      IF((KFA.GE.1.AND.KFA.LE.6).OR.KFA.EQ.21) THEN
        IF(K(I,2).GT.0.AND.K(I,2).LE.6) THEN
          IF(I1.EQ.0) THEN
            I1=I
          ELSEIF(I3.EQ.0) THEN
            I3=I
          ELSE
            CALL PYERRM(16,'(PY4JET:) more than two quarks')
          ENDIF
        ELSEIF(K(I,2).LT.0) THEN
          IF(I2.EQ.0) THEN
            I2=I
          ELSEIF(I4.EQ.0) THEN
            I4=I
          ELSE
            CALL PYERRM(16,'(PY4JET:) more than two antiquarks')
          ENDIF
        ELSE
          IF(I3.EQ.0) THEN
            I3=I
          ELSEIF(I4.EQ.0) THEN
            I4=I
          ELSE
            CALL PYERRM(16,'(PY4JET:) more than two gluons')
          ENDIF
        ENDIF
      ENDIF
  100 CONTINUE
 
C...Check that event is arranged according to conventions.
      IF(I1.EQ.0.OR.I2.EQ.0.OR.I3.EQ.0.OR.I4.EQ.0) THEN
        CALL PYERRM(16,'(PY4JET:) event contains too few partons')
      ENDIF
      IF(I2.LT.I1.OR.I3.LT.I2.OR.I4.LT.I3) THEN
        CALL PYERRM(6,'(PY4JET:) partons arranged in wrong order')
      ENDIF
 
C...Check whether second pair are quarks or gluons.
      IF(IABS(K(I3,2)).LT.10.AND.IABS(K(I4,2)).LT.10) THEN
        IQG34=1
      ELSEIF(K(I3,2).EQ.21.AND.K(I4,2).EQ.21) THEN
        IQG34=2
      ELSE
        CALL PYERRM(16,'(PY4JET:) second parton pair inconsistent')
      ENDIF
 
C...Boost partons to their cm frame.
      DO 110 J=1,4
        PTOT(J)=P(I1,J)+P(I2,J)+P(I3,J)+P(I4,J)
  110 CONTINUE
      ECM=SQRT(MAX(0D0,PTOT(4)**2-PTOT(1)**2-PTOT(2)**2-PTOT(3)**2))
      DO 120 J=1,3
        BETA(J)=PTOT(J)/PTOT(4)
  120 CONTINUE
      CALL PYROBO(I1,I1,0D0,0D0,-BETA(1),-BETA(2),-BETA(3))
      CALL PYROBO(I2,I2,0D0,0D0,-BETA(1),-BETA(2),-BETA(3))
      CALL PYROBO(I3,I3,0D0,0D0,-BETA(1),-BETA(2),-BETA(3))
      CALL PYROBO(I4,I4,0D0,0D0,-BETA(1),-BETA(2),-BETA(3))
      NSAV=N
 
C...Decide and set up shower history for q qbar q' qbar' events.
      IF(IQG34.EQ.1) THEN
        W1=PY4JTW(0,I1,I3,I4)
        W2=PY4JTW(0,I2,I3,I4)
        IF(W1.GT.PYR(0)*(W1+W2)) THEN
          CALL PY4JTS(0,I1,I3,I4,I2,QMAX)
        ELSE
          CALL PY4JTS(0,I2,I3,I4,I1,QMAX)
        ENDIF
 
C...Decide and set up shower history for q qbar g g events.
      ELSE
        W1=PY4JTW(I1,I3,I2,I4)
        W2=PY4JTW(I1,I4,I2,I3)
        W3=PY4JTW(0,I3,I1,I4)
        W4=PY4JTW(0,I4,I1,I3)
        W5=PY4JTW(0,I3,I2,I4)
        W6=PY4JTW(0,I4,I2,I3)
        W7=PY4JTW(0,I1,I3,I4)
        W8=PY4JTW(0,I2,I3,I4)
        WR=(W1+W2+W3+W4+W5+W6+W7+W8)*PYR(0)
        IF(W1.GT.WR) THEN
          CALL PY4JTS(I1,I3,I2,I4,0,QMAX)
        ELSEIF(W1+W2.GT.WR) THEN
          CALL PY4JTS(I1,I4,I2,I3,0,QMAX)
        ELSEIF(W1+W2+W3.GT.WR) THEN
          CALL PY4JTS(0,I3,I1,I4,I2,QMAX)
        ELSEIF(W1+W2+W3+W4.GT.WR) THEN
          CALL PY4JTS(0,I4,I1,I3,I2,QMAX)
        ELSEIF(W1+W2+W3+W4+W5.GT.WR) THEN
          CALL PY4JTS(0,I3,I2,I4,I1,QMAX)
        ELSEIF(W1+W2+W3+W4+W5+W6.GT.WR) THEN
          CALL PY4JTS(0,I4,I2,I3,I1,QMAX)
        ELSEIF(W1+W2+W3+W4+W5+W6+W7.GT.WR) THEN
          CALL PY4JTS(0,I1,I3,I4,I2,QMAX)
        ELSE
          CALL PY4JTS(0,I2,I3,I4,I1,QMAX)
        ENDIF
      ENDIF
 
C...Boost back original partons and mark them as deleted.
      CALL PYROBO(I1,I1,0D0,0D0,BETA(1),BETA(2),BETA(3))
      CALL PYROBO(I2,I2,0D0,0D0,BETA(1),BETA(2),BETA(3))
      CALL PYROBO(I3,I3,0D0,0D0,BETA(1),BETA(2),BETA(3))
      CALL PYROBO(I4,I4,0D0,0D0,BETA(1),BETA(2),BETA(3))
      K(I1,1)=K(I1,1)+10
      K(I2,1)=K(I2,1)+10
      K(I3,1)=K(I3,1)+10
      K(I4,1)=K(I4,1)+10
 
C...Rotate shower initiating partons to be along z axis.
      PHI=PYANGL(P(NSAV+1,1),P(NSAV+1,2))
      CALL PYROBO(NSAV+1,NSAV+6,0D0,-PHI,0D0,0D0,0D0)
      THE=PYANGL(P(NSAV+1,3),P(NSAV+1,1))
      CALL PYROBO(NSAV+1,NSAV+6,-THE,0D0,0D0,0D0,0D0)
 
C...Set up copy of shower initiating partons as on mass shell.
      DO 140 I=N+1,N+2
        DO 130 J=1,5
          K(I,J)=0
          P(I,J)=0D0
          V(I,J)=V(I1,J)
  130   CONTINUE
        K(I,1)=1
        K(I,2)=K(I-6,2)
  140 CONTINUE
      IF(K(NSAV+1,2).EQ.K(I1,2)) THEN
        K(N+1,3)=I1
        P(N+1,5)=P(I1,5)
        K(N+2,3)=I2
        P(N+2,5)=P(I2,5)
      ELSE
        K(N+1,3)=I2
        P(N+1,5)=P(I2,5)
        K(N+2,3)=I1
        P(N+2,5)=P(I1,5)
      ENDIF
      PABS=SQRT(MAX(0D0,(ECM**2-P(N+1,5)**2-P(N+2,5)**2)**2-
     &(2D0*P(N+1,5)*P(N+2,5))**2))/(2D0*ECM)
      P(N+1,3)=PABS
      P(N+1,4)=SQRT(PABS**2+P(N+1,5)**2)
      P(N+2,3)=-PABS
      P(N+2,4)=SQRT(PABS**2+P(N+2,5)**2)
      N=N+2
 
C...Decide whether to allow or not photon radiation in showers.
C...Connect up colours.
      MSTJ(41)=2
      IF(IRAD.EQ.0) MSTJ(41)=1
      IJOIN(1)=N-1
      IJOIN(2)=N
      CALL PYJOIN(2,IJOIN)
 
C...Decide on maximum virtuality and do parton shower.
      IF(PMAX.LT.PARJ(82)) THEN
        PQMAX=QMAX
      ELSE
        PQMAX=PMAX
      ENDIF
      CALL PYSHOW(NSAV+1,-100,PQMAX)
 
C...Rotate and boost back system.
      CALL PYROBO(NSAV+1,N,THE,PHI,BETA(1),BETA(2),BETA(3))
 
C...Do fragmentation and decays.
      CALL PYEXEC
 
C...Call PYHEPC to convert output from PYJETS to HEPEVT common.
      IF(ICOM.EQ.0) THEN
        MSTU(28)=0
        CALL PYHEPC(1)
      ENDIF
 
      RETURN
      END
