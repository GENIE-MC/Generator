 
C*********************************************************************
 
C...PYGAGA
C...For lepton beams it gives photon-hadron or photon-photon systems
C...to be treated with the ordinary machinery and combines this with a
C...description of the lepton -> lepton + photon branching.
 
      SUBROUTINE PYGAGA(IGAGA,WTGAGA)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYSUBS/,/PYPARS/,/PYINT1/,
     &/PYINT5/
C...Local variables and data statement.
      DIMENSION PMS(2),XMIN(2),XMAX(2),Q2MIN(2),Q2MAX(2),PMC(3),
     &X(2),Q2(2),Y(2),THETA(2),PHI(2),PT(2),BETA(3)
      SAVE PMS,XMIN,XMAX,Q2MIN,Q2MAX,PMC,X,Q2,THETA,PHI,PT,W2MIN
      DATA EPS/1D-4/
 
C...Initialize generation of photons inside leptons.
      IF(IGAGA.EQ.1) THEN
 
C...Save quantities on incoming lepton system.
        VINT(301)=VINT(1)
        VINT(302)=VINT(2)
        PMS(1)=VINT(303)**2
        IF(MINT(141).EQ.0) PMS(1)=SIGN(VINT(3)**2,VINT(3))
        PMS(2)=VINT(304)**2
        IF(MINT(142).EQ.0) PMS(2)=SIGN(VINT(4)**2,VINT(4))
        PMC(3)=VINT(302)-PMS(1)-PMS(2)
        W2MIN=MAX(CKIN(77),2D0*CKIN(3),2D0*CKIN(5))**2
 
C...Calculate range of x and Q2 values allowed in generation.
        DO 100 I=1,2
          PMC(I)=VINT(302)+PMS(I)-PMS(3-I)
          IF(MINT(140+I).NE.0) THEN
            XMIN(I)=MAX(CKIN(59+2*I),EPS)
            XMAX(I)=MIN(CKIN(60+2*I),1D0-2D0*VINT(301)*SQRT(PMS(I))/
     &      PMC(I),1D0-EPS)
            YMIN=MAX(CKIN(71+2*I),EPS)
            YMAX=MIN(CKIN(72+2*I),1D0-EPS)
            IF(CKIN(64+2*I).GT.0D0) XMIN(I)=MAX(XMIN(I),
     &      (YMIN*PMC(3)-CKIN(64+2*I))/PMC(I))
            XMAX(I)=MIN(XMAX(I),(YMAX*PMC(3)-CKIN(63+2*I))/PMC(I))
            THEMIN=MAX(CKIN(67+2*I),0D0)
            THEMAX=MIN(CKIN(68+2*I),PARU(1))
            IF(CKIN(68+2*I).LT.0D0) THEMAX=PARU(1)
            Q2MIN(I)=MAX(CKIN(63+2*I),XMIN(I)**2*PMS(I)/(1D0-XMIN(I))+
     &      ((1D0-XMAX(I))*(VINT(302)-2D0*PMS(3-I))-
     &      2D0*PMS(I)/(1D0-XMAX(I)))*SIN(THEMIN/2D0)**2,0D0)
            Q2MAX(I)=XMAX(I)**2*PMS(I)/(1D0-XMAX(I))+
     &      ((1D0-XMIN(I))*(VINT(302)-2D0*PMS(3-I))-
     &      2D0*PMS(I)/(1D0-XMIN(I)))*SIN(THEMAX/2D0)**2
            IF(CKIN(64+2*I).GT.0D0) Q2MAX(I)=MIN(CKIN(64+2*I),Q2MAX(I))
C...W limits when lepton on one side only.
            IF(MINT(143-I).EQ.0) THEN
              XMIN(I)=MAX(XMIN(I),(W2MIN-PMS(3-I))/PMC(I))
              IF(CKIN(78).GT.0D0) XMAX(I)=MIN(XMAX(I),
     &        (CKIN(78)**2-PMS(3-I))/PMC(I))
            ENDIF
          ENDIF
  100   CONTINUE
 
C...W limits when lepton on both sides.
        IF(MINT(141).NE.0.AND.MINT(142).NE.0) THEN
          IF(CKIN(78).GT.0D0) XMAX(1)=MIN(XMAX(1),
     &    (CKIN(78)**2+PMC(3)-PMC(2)*XMIN(2))/PMC(1))
          IF(CKIN(78).GT.0D0) XMAX(2)=MIN(XMAX(2),
     &    (CKIN(78)**2+PMC(3)-PMC(1)*XMIN(1))/PMC(2))
          IF(IABS(MINT(141)).NE.IABS(MINT(142))) THEN
            XMIN(1)=MAX(XMIN(1),(PMS(1)-PMS(2)+VINT(302)*(W2MIN-
     &      PMS(1)-PMS(2))/(PMC(2)*XMAX(2)+PMS(1)-PMS(2)))/PMC(1))
            XMIN(2)=MAX(XMIN(2),(PMS(2)-PMS(1)+VINT(302)*(W2MIN-
     &      PMS(1)-PMS(2))/(PMC(1)*XMAX(1)+PMS(2)-PMS(1)))/PMC(2))
          ELSE
            XMIN(1)=MAX(XMIN(1),W2MIN/(VINT(302)*XMAX(2)))
            XMIN(2)=MAX(XMIN(2),W2MIN/(VINT(302)*XMAX(1)))
          ENDIF
        ENDIF
 
C...Q2 and W values and photon flux weight factors for initialization.
      ELSEIF(IGAGA.EQ.2) THEN
        ISUB=MINT(1)
        MINT(15)=0
        MINT(16)=0
 
C...W value for photon on one or both sides, and for processes
C...with gamma-gamma cross section peaked at small shat.
        IF(MINT(141).NE.0.AND.MINT(142).EQ.0) THEN
          VINT(2)=VINT(302)+PMS(1)-PMC(1)*(1D0-XMAX(1))
        ELSEIF(MINT(141).EQ.0.AND.MINT(142).NE.0) THEN
          VINT(2)=VINT(302)+PMS(2)-PMC(2)*(1D0-XMAX(2))
        ELSEIF(ISUB.GE.137.AND.ISUB.LE.140) THEN
          VINT(2)=MAX(CKIN(77)**2,12D0*MAX(CKIN(3),CKIN(5))**2)
          IF(CKIN(78).GT.0D0) VINT(2)=MIN(VINT(2),CKIN(78)**2)
        ELSE
          VINT(2)=XMAX(1)*XMAX(2)*VINT(302)
          IF(CKIN(78).GT.0D0) VINT(2)=MIN(VINT(2),CKIN(78)**2)
        ENDIF
        VINT(1)=SQRT(MAX(0D0,VINT(2)))
 
C...Upper estimate of photon flux weight factor.
C...Initialization Q2 scale. Flag incoming unresolved photon.
        WTGAGA=1D0
        DO 110 I=1,2
          IF(MINT(140+I).NE.0) THEN
            WTGAGA=WTGAGA*2D0*(PARU(101)/PARU(2))*
     &      LOG(XMAX(I)/XMIN(I))*LOG(Q2MAX(I)/Q2MIN(I))
            IF(ISUB.EQ.99.AND.MINT(106+I).EQ.4.AND.MINT(109-I).EQ.3)
     &      THEN
              Q2INIT=5D0+Q2MIN(3-I)
            ELSEIF(ISUB.EQ.99.AND.MINT(106+I).EQ.4) THEN
              Q2INIT=PMAS(PYCOMP(113),1)**2+Q2MIN(3-I)
            ELSEIF(ISUB.EQ.132.OR.ISUB.EQ.134.OR.ISUB.EQ.136) THEN
              Q2INIT=MAX(CKIN(1),2D0*CKIN(3),2D0*CKIN(5))**2/3D0
            ELSEIF((ISUB.EQ.138.AND.I.EQ.2).OR.
     &      (ISUB.EQ.139.AND.I.EQ.1)) THEN
              Q2INIT=VINT(2)/3D0
            ELSEIF(ISUB.EQ.140) THEN
              Q2INIT=VINT(2)/2D0
            ELSE
              Q2INIT=Q2MIN(I)
            ENDIF
            VINT(2+I)=-SQRT(MAX(Q2MIN(I),MIN(Q2MAX(I),Q2INIT)))
            IF(MSTP(14).EQ.0.OR.(ISUB.GE.131.AND.ISUB.LE.140))
     &      MINT(14+I)=22
            VINT(306+I)=VINT(2+I)**2
          ENDIF
  110   CONTINUE
        VINT(320)=WTGAGA
 
C...Update pTmin and cross section information.
        IF(MSTP(82).LE.1) THEN
          PTMN=PARP(81)*(VINT(1)/PARP(89))**PARP(90)
        ELSE
          PTMN=PARP(82)*(VINT(1)/PARP(89))**PARP(90)
        ENDIF
        VINT(149)=4D0*PTMN**2/VINT(2)
        VINT(154)=PTMN
        CALL PYXTOT
        VINT(318)=VINT(317)
 
C...Generate photons inside leptons and
C...calculate photon flux weight factors.
      ELSEIF(IGAGA.EQ.3) THEN
        ISUB=MINT(1)
        MINT(15)=0
        MINT(16)=0
 
C...Generate phase space point and check against cuts.
        LOOP=0
  120   LOOP=LOOP+1
        DO 130 I=1,2
          IF(MINT(140+I).NE.0) THEN
C...Pick x and Q2
            X(I)=XMIN(I)*(XMAX(I)/XMIN(I))**PYR(0)
            Q2(I)=Q2MIN(I)*(Q2MAX(I)/Q2MIN(I))**PYR(0)
C...Cuts on internal consistency in x and Q2.
            IF(Q2(I).LT.X(I)**2*PMS(I)/(1D0-X(I))) GOTO 120
            IF(Q2(I).GT.(1D0-X(I))*(VINT(302)-2D0*PMS(3-I))-
     &      (2D0-X(I)**2)*PMS(I)/(1D0-X(I))) GOTO 120
C...Cuts on y and theta.
            Y(I)=(PMC(I)*X(I)+Q2(I))/PMC(3)
            IF(Y(I).LT.CKIN(71+2*I).OR.Y(I).GT.CKIN(72+2*I)) GOTO 120
            RAT=((1D0-X(I))*Q2(I)-X(I)**2*PMS(I))/
     &      ((1D0-X(I))**2*(VINT(302)-2D0*PMS(3-I)-2D0*PMS(I)))
            THETA(I)=2D0*ASIN(SQRT(MAX(0D0,MIN(1D0,RAT))))
            IF(THETA(I).LT.CKIN(67+2*I)) GOTO 120
            IF(CKIN(68+2*I).GT.0D0.AND.THETA(I).GT.CKIN(68+2*I))
     &      GOTO 120
 
C...Phi angle isotropic. Reconstruct pT.
            PHI(I)=PARU(2)*PYR(0)
            PT(I)=SQRT(((1D0-X(I))*PMC(I))**2/(4D0*VINT(302))-
     &      PMS(I))*SIN(THETA(I))
 
C...Store info on variables selected, for documentation purposes.
            VINT(2+I)=-SQRT(Q2(I))
            VINT(304+I)=X(I)
            VINT(306+I)=Q2(I)
            VINT(308+I)=Y(I)
            VINT(310+I)=THETA(I)
            VINT(312+I)=PHI(I)
          ELSE
            VINT(304+I)=1D0
            VINT(306+I)=0D0
            VINT(308+I)=1D0
            VINT(310+I)=0D0
            VINT(312+I)=0D0
          ENDIF
  130   CONTINUE
 
C...Cut on W combines info from two sides.
        IF(MINT(141).NE.0.AND.MINT(142).NE.0) THEN
          W2=-Q2(1)-Q2(2)+0.5D0*X(1)*PMC(1)*X(2)*PMC(2)/VINT(302)-
     &    2D0*PT(1)*PT(2)*COS(PHI(1)-PHI(2))+2D0*
     &    SQRT((0.5D0*X(1)*PMC(1)/VINT(301))**2+Q2(1)-PT(1)**2)*
     &    SQRT((0.5D0*X(2)*PMC(2)/VINT(301))**2+Q2(2)-PT(2)**2)
          IF(W2.LT.W2MIN) GOTO 120
          IF(CKIN(78).GT.0D0.AND.W2.GT.CKIN(78)**2) GOTO 120
          PMS1=-Q2(1)
          PMS2=-Q2(2)
        ELSEIF(MINT(141).NE.0) THEN
          W2=(VINT(302)+PMS(1))*X(1)+PMS(2)*(1D0-X(1))
          PMS1=-Q2(1)
          PMS2=PMS(2)
        ELSEIF(MINT(142).NE.0) THEN
          W2=(VINT(302)+PMS(2))*X(2)+PMS(1)*(1D0-X(2))
          PMS1=PMS(1)
          PMS2=-Q2(2)
        ENDIF
 
C...Store kinematics info for photon(s) in subsystem cm frame.
        VINT(2)=W2
        VINT(1)=SQRT(W2)
        VINT(291)=0D0
        VINT(292)=0D0
        VINT(293)=0.5D0*SQRT((W2-PMS1-PMS2)**2-4D0*PMS1*PMS2)/VINT(1)
        VINT(294)=0.5D0*(W2+PMS1-PMS2)/VINT(1)
        VINT(295)=SIGN(SQRT(ABS(PMS1)),PMS1)
        VINT(296)=0D0
        VINT(297)=0D0
        VINT(298)=-VINT(293)
        VINT(299)=0.5D0*(W2+PMS2-PMS1)/VINT(1)
        VINT(300)=SIGN(SQRT(ABS(PMS2)),PMS2)
 
C...Assign weight for photon flux; different for transverse and
C...longitudinal photons. Flag incoming unresolved photon.
        WTGAGA=1D0
        DO 140 I=1,2
          IF(MINT(140+I).NE.0) THEN
            WTGAGA=WTGAGA*2D0*(PARU(101)/PARU(2))*
     &      LOG(XMAX(I)/XMIN(I))*LOG(Q2MAX(I)/Q2MIN(I))
            IF(MSTP(16).EQ.0) THEN
              XY=X(I)
            ELSE
              WTGAGA=WTGAGA*X(I)/Y(I)
              XY=Y(I)
            ENDIF
            IF(ISUB.EQ.132.OR.ISUB.EQ.134.OR.ISUB.EQ.136) THEN
              WTGAGA=WTGAGA*(1D0-XY)
            ELSEIF(I.EQ.1.AND.(ISUB.EQ.139.OR.ISUB.EQ.140)) THEN
              WTGAGA=WTGAGA*(1D0-XY)
            ELSEIF(I.EQ.2.AND.(ISUB.EQ.138.OR.ISUB.EQ.140)) THEN
              WTGAGA=WTGAGA*(1D0-XY)
            ELSE
              WTGAGA=WTGAGA*(0.5D0*(1D0+(1D0-XY)**2)-
     &        PMS(I)*XY**2/Q2(I))
            ENDIF
            IF(MINT(106+I).EQ.0) MINT(14+I)=22
          ENDIF
  140   CONTINUE
        VINT(319)=WTGAGA
        MINT(143)=LOOP
 
C...Update pTmin and cross section information.
        IF(MSTP(82).LE.1) THEN
          PTMN=PARP(81)*(VINT(1)/PARP(89))**PARP(90)
        ELSE
          PTMN=PARP(82)*(VINT(1)/PARP(89))**PARP(90)
        ENDIF
        VINT(149)=4D0*PTMN**2/VINT(2)
        VINT(154)=PTMN
        CALL PYXTOT
 
C...Reconstruct kinematics of photons inside leptons.
      ELSEIF(IGAGA.EQ.4) THEN
 
C...Make place for incoming particles and scattered leptons.
        MOVE=3
        IF(MINT(141).NE.0.AND.MINT(142).NE.0) MOVE=4
        MINT(4)=MINT(4)+MOVE
        DO 160 I=MINT(84)-MOVE,MINT(83)+1,-1
          IF(K(I,1).EQ.21) THEN
            DO 150 J=1,5
              K(I+MOVE,J)=K(I,J)
              P(I+MOVE,J)=P(I,J)
              V(I+MOVE,J)=V(I,J)
  150       CONTINUE
            IF(K(I,3).GT.MINT(83).AND.K(I,3).LE.MINT(84))
     &      K(I+MOVE,3)=K(I,3)+MOVE
            IF(K(I,4).GT.MINT(83).AND.K(I,4).LE.MINT(84))
     &      K(I+MOVE,4)=K(I,4)+MOVE
            IF(K(I,5).GT.MINT(83).AND.K(I,5).LE.MINT(84))
     &      K(I+MOVE,5)=K(I,5)+MOVE
          ENDIF
  160   CONTINUE
        DO 170 I=MINT(84)+1,N
          IF(K(I,3).GT.MINT(83).AND.K(I,3).LE.MINT(84))
     &    K(I,3)=K(I,3)+MOVE
  170   CONTINUE
 
C...Fill in incoming particles.
        DO 190 I=MINT(83)+1,MINT(83)+MOVE
          DO 180 J=1,5
            K(I,J)=0
            P(I,J)=0D0
            V(I,J)=0D0
  180     CONTINUE
  190   CONTINUE
        DO 200 I=1,2
          K(MINT(83)+I,1)=21
          IF(MINT(140+I).NE.0) THEN
            K(MINT(83)+I,2)=MINT(140+I)
            P(MINT(83)+I,5)=VINT(302+I)
          ELSE
            K(MINT(83)+I,2)=MINT(10+I)
            P(MINT(83)+I,5)=VINT(2+I)
          ENDIF
          P(MINT(83)+I,3)=0.5D0*SQRT((PMC(3)**2-4D0*PMS(1)*PMS(2))/
     &    VINT(302))*(-1D0)**(I+1)
          P(MINT(83)+I,4)=0.5D0*PMC(I)/VINT(301)
  200   CONTINUE
 
C...New mother-daughter relations in documentation section.
        IF(MINT(141).NE.0.AND.MINT(142).NE.0) THEN
          K(MINT(83)+1,4)=MINT(83)+3
          K(MINT(83)+1,5)=MINT(83)+5
          K(MINT(83)+2,4)=MINT(83)+4
          K(MINT(83)+2,5)=MINT(83)+6
          K(MINT(83)+3,3)=MINT(83)+1
          K(MINT(83)+5,3)=MINT(83)+1
          K(MINT(83)+4,3)=MINT(83)+2
          K(MINT(83)+6,3)=MINT(83)+2
        ELSEIF(MINT(141).NE.0) THEN
          K(MINT(83)+1,4)=MINT(83)+3
          K(MINT(83)+1,5)=MINT(83)+4
          K(MINT(83)+2,4)=MINT(83)+5
          K(MINT(83)+3,3)=MINT(83)+1
          K(MINT(83)+4,3)=MINT(83)+1
          K(MINT(83)+5,3)=MINT(83)+2
        ELSEIF(MINT(142).NE.0) THEN
          K(MINT(83)+1,4)=MINT(83)+4
          K(MINT(83)+2,4)=MINT(83)+3
          K(MINT(83)+2,5)=MINT(83)+5
          K(MINT(83)+3,3)=MINT(83)+2
          K(MINT(83)+4,3)=MINT(83)+1
          K(MINT(83)+5,3)=MINT(83)+2
        ENDIF
 
C...Fill scattered lepton(s).
        DO 210 I=1,2
          IF(MINT(140+I).NE.0) THEN
            LSC=MINT(83)+MIN(I+2,MOVE)
            K(LSC,1)=21
            K(LSC,2)=MINT(140+I)
            P(LSC,1)=PT(I)*COS(PHI(I))
            P(LSC,2)=PT(I)*SIN(PHI(I))
            P(LSC,4)=(1D0-X(I))*P(MINT(83)+I,4)
            P(LSC,3)=SQRT(P(LSC,4)**2-PMS(I))*COS(THETA(I))*
     &      (-1D0)**(I-1)
            P(LSC,5)=VINT(302+I)
          ENDIF
  210   CONTINUE
 
C...Find incoming four-vectors to subprocess.
        K(N+1,1)=21
        IF(MINT(141).NE.0) THEN
          DO 220 J=1,4
            P(N+1,J)=P(MINT(83)+1,J)-P(MINT(83)+3,J)
  220     CONTINUE
        ELSE
          DO 230 J=1,4
            P(N+1,J)=P(MINT(83)+1,J)
  230     CONTINUE
        ENDIF
        K(N+2,1)=21
        IF(MINT(142).NE.0) THEN
          DO 240 J=1,4
            P(N+2,J)=P(MINT(83)+2,J)-P(MINT(83)+MOVE,J)
  240     CONTINUE
        ELSE
          DO 250 J=1,4
            P(N+2,J)=P(MINT(83)+2,J)
  250     CONTINUE
        ENDIF
 
C...Define boost and rotation between hadronic subsystem and
C...collision rest frame; boost hadronic subsystem to this frame.
        DO 260 J=1,3
          BETA(J)=(P(N+1,J)+P(N+2,J))/(P(N+1,4)+P(N+2,4))
  260   CONTINUE
        CALL PYROBO(N+1,N+2,0D0,0D0,-BETA(1),-BETA(2),-BETA(3))
        BPHI=PYANGL(P(N+1,1),P(N+1,2))
        CALL PYROBO(N+1,N+2,0D0,-BPHI,0D0,0D0,0D0)
        BTHETA=PYANGL(P(N+1,3),P(N+1,1))
        CALL PYROBO(MINT(83)+MOVE+1,N,BTHETA,BPHI,BETA(1),BETA(2),
     &  BETA(3))
 
C...Add on scattered leptons to final state.
        DO 280 I=1,2
          IF(MINT(140+I).NE.0) THEN
            LSC=MINT(83)+MIN(I+2,MOVE)
            N=N+1
            DO 270 J=1,5
              K(N,J)=K(LSC,J)
              P(N,J)=P(LSC,J)
              V(N,J)=V(LSC,J)
  270       CONTINUE
            K(N,1)=1
            K(N,3)=LSC
          ENDIF
  280   CONTINUE
      ENDIF
 
      RETURN
      END
