 
C*********************************************************************
 
C...PYREMN
C...Adds on target remnants (one or two from each side) and
C...includes primordial kT for hadron beams.
 
      SUBROUTINE PYREMN(IPU1,IPU2)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/
C...Local arrays.
      DIMENSION KFLCH(2),KFLSP(2),CHI(2),PMS(0:6),IS(2),ISN(2),ROBO(5),
     &PSYS(0:2,5),PMIN(0:2),QOLD(4),QNEW(4),DBE(3),PSUM(4)
 
C...Find event type and remaining energy.
      ISUB=MINT(1)
      NS=N
      IF(MINT(50).EQ.0.OR.MOD(MSTP(81),10).LE.0) THEN
        VINT(143)=1D0-VINT(141)
        VINT(144)=1D0-VINT(142)
      ENDIF
 
C...Define initial partons.
      NTRY=0
  100 NTRY=NTRY+1
      DO 130 JT=1,2
        I=MINT(83)+JT+2
        IF(JT.EQ.1) IPU=IPU1
        IF(JT.EQ.2) IPU=IPU2
        K(I,1)=21
        K(I,2)=K(IPU,2)
        K(I,3)=I-2
        PMS(JT)=0D0
        VINT(156+JT)=0D0
        VINT(158+JT)=0D0
        IF(MINT(47).EQ.1) THEN
          DO 110 J=1,5
            P(I,J)=P(I-2,J)
  110     CONTINUE
        ELSEIF(ISUB.EQ.95) THEN
          K(I,2)=21
        ELSE
          P(I,5)=P(IPU,5)
 
C...No primordial kT, or chosen according to truncated Gaussian or
C...exponential, or (for photon) predetermined or power law.
  120     IF(MINT(40+JT).EQ.2.AND.MINT(10+JT).NE.22) THEN
            IF(MSTP(91).LE.0) THEN
              PT=0D0
            ELSEIF(MSTP(91).EQ.1) THEN
              PT=PARP(91)*SQRT(-LOG(PYR(0)))
            ELSE
              RPT1=PYR(0)
              RPT2=PYR(0)
              PT=-PARP(92)*LOG(RPT1*RPT2)
            ENDIF
            IF(PT.GT.PARP(93)) GOTO 120
          ELSEIF(MINT(106+JT).EQ.3) THEN
            PTA=SQRT(VINT(282+JT))
            PTB=0D0
            IF(MSTP(66).EQ.5.AND.MSTP(93).EQ.1) THEN
              PTB=PARP(99)*SQRT(-LOG(PYR(0)))
            ELSEIF(MSTP(66).EQ.5.AND.MSTP(93).EQ.2) THEN
              RPT1=PYR(0)
              RPT2=PYR(0)
              PTB=-PARP(99)*LOG(RPT1*RPT2)
            ENDIF
            IF(PTB.GT.PARP(100)) GOTO 120
            PT=SQRT(PTA**2+PTB**2+2D0*PTA*PTB*COS(PARU(2)*PYR(0)))
            PT=PT*0.8D0**MINT(57)
            IF(NTRY.GT.10) PT=PT*0.8D0**(NTRY-10)
          ELSEIF(IABS(MINT(14+JT)).LE.8.OR.MINT(14+JT).EQ.21) THEN
            IF(MSTP(93).LE.0) THEN
              PT=0D0
            ELSEIF(MSTP(93).EQ.1) THEN
              PT=PARP(99)*SQRT(-LOG(PYR(0)))
            ELSEIF(MSTP(93).EQ.2) THEN
              RPT1=PYR(0)
              RPT2=PYR(0)
              PT=-PARP(99)*LOG(RPT1*RPT2)
            ELSEIF(MSTP(93).EQ.3) THEN
              HA=PARP(99)**2
              HB=PARP(100)**2
              PT=SQRT(MAX(0D0,HA*(HA+HB)/(HA+HB-PYR(0)*HB)-HA))
            ELSE
              HA=PARP(99)**2
              HB=PARP(100)**2
              IF(MSTP(93).EQ.5) HB=MIN(VINT(48),PARP(100)**2)
              PT=SQRT(MAX(0D0,HA*((HA+HB)/HA)**PYR(0)-HA))
            ENDIF
            IF(PT.GT.PARP(100)) GOTO 120
          ELSE
            PT=0D0
          ENDIF
          VINT(156+JT)=PT
          PHI=PARU(2)*PYR(0)
          P(I,1)=PT*COS(PHI)
          P(I,2)=PT*SIN(PHI)
          PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2
        ENDIF
  130 CONTINUE
      IF(MINT(47).EQ.1) RETURN
 
C...Kinematics construction for initial partons.
      I1=MINT(83)+3
      I2=MINT(83)+4
      IF(ISUB.EQ.95) THEN
        SHS=0D0
        SHR=0D0
      ELSE
        SHS=VINT(141)*VINT(142)*VINT(2)+(P(I1,1)+P(I2,1))**2+
     &  (P(I1,2)+P(I2,2))**2
        SHR=SQRT(MAX(0D0,SHS))
        IF((SHS-PMS(1)-PMS(2))**2-4D0*PMS(1)*PMS(2).LE.0D0) GOTO 100
        P(I1,4)=0.5D0*(SHR+(PMS(1)-PMS(2))/SHR)
        P(I1,3)=SQRT(MAX(0D0,P(I1,4)**2-PMS(1)))
        P(I2,4)=SHR-P(I1,4)
        P(I2,3)=-P(I1,3)
 
C...Transform partons to overall CM-frame.
        ROBO(3)=(P(I1,1)+P(I2,1))/SHR
        ROBO(4)=(P(I1,2)+P(I2,2))/SHR
        CALL PYROBO(I1,I2,0D0,0D0,-ROBO(3),-ROBO(4),0D0)
        ROBO(2)=PYANGL(P(I1,1),P(I1,2))
        CALL PYROBO(I1,I2,0D0,-ROBO(2),0D0,0D0,0D0)
        ROBO(1)=PYANGL(P(I1,3),P(I1,1))
        CALL PYROBO(I1,I2,-ROBO(1),0D0,0D0,0D0,0D0)
        CALL PYROBO(I2+1,MINT(52),0D0,-ROBO(2),0D0,0D0,0D0)
        CALL PYROBO(I1,MINT(52),ROBO(1),ROBO(2),ROBO(3),ROBO(4),0D0)
        ROBO(5)=(VINT(141)-VINT(142))/(VINT(141)+VINT(142))
        CALL PYROBO(I1,MINT(52),0D0,0D0,0D0,0D0,ROBO(5))
      ENDIF
 
C...Optionally fix up x and Q2 definitions for leptoproduction.
      IDISXQ=0
      IF((MINT(43).EQ.2.OR.MINT(43).EQ.3).AND.((ISUB.EQ.10.AND.
     &MSTP(23).GE.1).OR.(ISUB.EQ.83.AND.MSTP(23).GE.2))) IDISXQ=1
      IF(IDISXQ.EQ.1) THEN
 
C...Find where incoming and outgoing leptons/partons are sitting.
        LESD=1
        IF(MINT(42).EQ.1) LESD=2
        LPIN=MINT(83)+3-LESD
        LEIN=MINT(84)+LESD
        LQIN=MINT(84)+3-LESD
        LEOUT=MINT(84)+2+LESD
        LQOUT=MINT(84)+5-LESD
        IF(K(LEIN,3).GT.LEIN) LEIN=K(LEIN,3)
        IF(K(LQIN,3).GT.LQIN) LQIN=K(LQIN,3)
        LSCMS=0
        DO 140 I=MINT(84)+5,N
          IF(K(I,2).EQ.94) THEN
            LSCMS=I
            LEOUT=I+LESD
            LQOUT=I+3-LESD
          ENDIF
  140   CONTINUE
        LQBG=IPU1
        IF(LESD.EQ.1) LQBG=IPU2
 
C...Calculate actual and wanted momentum transfer.
        XNOM=VINT(43-LESD)
        Q2NOM=-VINT(45)
        HPK=2D0*(P(LPIN,4)*P(LEIN,4)-P(LPIN,1)*P(LEIN,1)-
     &  P(LPIN,2)*P(LEIN,2)-P(LPIN,3)*P(LEIN,3))*
     &  (P(MINT(83)+LESD,4)*VINT(40+LESD)/P(LEIN,4))
        HPT2=MAX(0D0,Q2NOM*(1D0-Q2NOM/(XNOM*HPK)))
        FAC=SQRT(HPT2/(P(LEOUT,1)**2+P(LEOUT,2)**2))
        P(N+1,1)=FAC*P(LEOUT,1)
        P(N+1,2)=FAC*P(LEOUT,2)
        P(N+1,3)=0.25D0*((HPK-Q2NOM/XNOM)/P(LPIN,4)-
     &  Q2NOM/(P(MINT(83)+LESD,4)*VINT(40+LESD)))*(-1)**(LESD+1)
        P(N+1,4)=SQRT(P(LEOUT,5)**2+P(N+1,1)**2+P(N+1,2)**2+
     &  P(N+1,3)**2)
        DO 150 J=1,4
          QOLD(J)=P(LEIN,J)-P(LEOUT,J)
          QNEW(J)=P(LEIN,J)-P(N+1,J)
  150   CONTINUE
 
C...Boost outgoing electron and daughters.
        IF(LSCMS.EQ.0) THEN
          DO 160 J=1,4
            P(LEOUT,J)=P(N+1,J)
  160     CONTINUE
        ELSE
          DO 170 J=1,3
            P(N+2,J)=(P(N+1,J)-P(LEOUT,J))/(P(N+1,4)+P(LEOUT,4))
  170     CONTINUE
          PINV=2D0/(1D0+P(N+2,1)**2+P(N+2,2)**2+P(N+2,3)**2)
          DO 180 J=1,3
            DBE(J)=PINV*P(N+2,J)
  180     CONTINUE
          DO 200 I=LSCMS+1,N
            IORIG=I
  190       IORIG=K(IORIG,3)
            IF(IORIG.GT.LEOUT) GOTO 190
            IF(I.EQ.LEOUT.OR.IORIG.EQ.LEOUT)
     &      CALL PYROBO(I,I,0D0,0D0,DBE(1),DBE(2),DBE(3))
  200     CONTINUE
        ENDIF
 
C...Copy shower initiator and all outgoing partons.
        NCOP=N+1
        K(NCOP,3)=LQBG
        DO 210 J=1,5
          P(NCOP,J)=P(LQBG,J)
  210   CONTINUE
        DO 240 I=MINT(84)+1,N
          ICOP=0
          IF(K(I,1).GT.10) GOTO 240
          IF(I.EQ.LQBG.OR.I.EQ.LQOUT) THEN
            ICOP=I
          ELSE
            IORIG=I
  220       IORIG=K(IORIG,3)
            IF(IORIG.EQ.LQBG.OR.IORIG.EQ.LQOUT) THEN
              ICOP=IORIG
            ELSEIF(IORIG.GT.MINT(84).AND.IORIG.LE.N) THEN
              GOTO 220
            ENDIF
          ENDIF
          IF(ICOP.NE.0) THEN
            NCOP=NCOP+1
            K(NCOP,3)=I
            DO 230 J=1,5
              P(NCOP,J)=P(I,J)
  230       CONTINUE
          ENDIF
  240   CONTINUE
 
C...Calculate relative rescaling factors.
        SLC=3-2*LESD
        PLCSUM=0D0
        DO 250 I=N+2,NCOP
          PLCSUM=PLCSUM+(P(I,4)+SLC*P(I,3))
  250   CONTINUE
        DO 260 I=N+2,NCOP
          V(I,1)=(P(I,4)+SLC*P(I,3))/PLCSUM
  260   CONTINUE
 
C...Transfer extra three-momentum of current.
        DO 280 I=N+2,NCOP
          DO 270 J=1,3
            P(I,J)=P(I,J)+V(I,1)*(QNEW(J)-QOLD(J))
  270     CONTINUE
          P(I,4)=SQRT(P(I,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2)
  280   CONTINUE
 
C...Iterate change of initiator momentum to get energy right.
        ITER=0
  290   ITER=ITER+1
        PEEX=-P(N+1,4)-QNEW(4)
        PEMV=-P(N+1,3)/P(N+1,4)
        DO 300 I=N+2,NCOP
          PEEX=PEEX+P(I,4)
          PEMV=PEMV+V(I,1)*P(I,3)/P(I,4)
  300   CONTINUE
        IF(ABS(PEMV).LT.1D-10) THEN
          MINT(51)=1
          MINT(57)=MINT(57)+1
          RETURN
        ENDIF
        PZCH=-PEEX/PEMV
        P(N+1,3)=P(N+1,3)+PZCH
        P(N+1,4)=SQRT(P(N+1,5)**2+P(N+1,1)**2+P(N+1,2)**2+P(N+1,3)**2)
        DO 310 I=N+2,NCOP
          P(I,3)=P(I,3)+V(I,1)*PZCH
          P(I,4)=SQRT(P(I,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2)
  310   CONTINUE
        IF(ITER.LT.10.AND.ABS(PEEX).GT.1D-6*P(N+1,4)) GOTO 290
 
C...Modify momenta in event record.
        HBE=2D0*(P(N+1,4)+P(LQBG,4))*(P(N+1,3)-P(LQBG,3))/
     &  ((P(N+1,4)+P(LQBG,4))**2+(P(N+1,3)-P(LQBG,3))**2)
        IF(ABS(HBE).GE.1D0) THEN
          MINT(51)=1
          MINT(57)=MINT(57)+1
          RETURN
        ENDIF
        I=MINT(83)+5-LESD
        CALL PYROBO(I,I,0D0,0D0,0D0,0D0,HBE)
        DO 330 I=N+1,NCOP
          ICOP=K(I,3)
          DO 320 J=1,4
            P(ICOP,J)=P(I,J)
  320     CONTINUE
  330   CONTINUE
      ENDIF
 
C...Check minimum invariant mass of remnant system(s).
      PSYS(0,4)=P(I1,4)+P(I2,4)+0.5D0*VINT(1)*(VINT(151)+VINT(152))
      PSYS(0,3)=P(I1,3)+P(I2,3)+0.5D0*VINT(1)*(VINT(151)-VINT(152))
      PMS(0)=MAX(0D0,PSYS(0,4)**2-PSYS(0,3)**2)
      PMIN(0)=SQRT(PMS(0))
      DO 340 JT=1,2
        PSYS(JT,4)=0.5D0*VINT(1)*VINT(142+JT)
        PSYS(JT,3)=PSYS(JT,4)*(-1)**(JT-1)
        PMIN(JT)=0D0
        IF(MINT(44+JT).EQ.1) GOTO 340
        MINT(105)=MINT(102+JT)
        MINT(109)=MINT(106+JT)
        CALL PYSPLI(MINT(10+JT),MINT(12+JT),KFLCH(JT),KFLSP(JT))
        IF(MINT(51).NE.0) THEN
          MINT(57)=MINT(57)+1
          RETURN
        ENDIF
        IF(KFLCH(JT).NE.0) PMIN(JT)=PMIN(JT)+PYMASS(KFLCH(JT))
        IF(KFLSP(JT).NE.0) PMIN(JT)=PMIN(JT)+PYMASS(KFLSP(JT))
        IF(KFLCH(JT)*KFLSP(JT).NE.0) PMIN(JT)=PMIN(JT)+0.5D0*PARP(111)
        PMIN(JT)=SQRT(PMIN(JT)**2+P(MINT(83)+JT+2,1)**2+
     &  P(MINT(83)+JT+2,2)**2)
  340 CONTINUE
      IF(PMIN(0)+PMIN(1)+PMIN(2).GT.VINT(1).OR.(MINT(45).GE.2.AND.
     &PMIN(1).GT.PSYS(1,4)).OR.(MINT(46).GE.2.AND.PMIN(2).GT.
     &PSYS(2,4))) THEN
        MINT(51)=1
        MINT(57)=MINT(57)+1
        RETURN
      ENDIF
 
C...Loop over two remnants; skip if none there.
      I=NS
      DO 410 JT=1,2
        ISN(JT)=0
        IF(MINT(44+JT).EQ.1) GOTO 410
        IF(JT.EQ.1) IPU=IPU1
        IF(JT.EQ.2) IPU=IPU2
 
C...Store first remnant parton.
        I=I+1
        IS(JT)=I
        ISN(JT)=1
        DO 350 J=1,5
          K(I,J)=0
          P(I,J)=0D0
          V(I,J)=0D0
  350   CONTINUE
        K(I,1)=1
        K(I,2)=KFLSP(JT)
        K(I,3)=MINT(83)+JT
        P(I,5)=PYMASS(K(I,2))
 
C...First parton colour connections and kinematics.
        KCOL=KCHG(PYCOMP(KFLSP(JT)),2)
        IF(KCOL.EQ.2) THEN
          K(I,1)=3
          K(I,4)=MSTU(5)*IPU+IPU
          K(I,5)=MSTU(5)*IPU+IPU
          K(IPU,4)=MOD(K(IPU,4),MSTU(5))+MSTU(5)*I
          K(IPU,5)=MOD(K(IPU,5),MSTU(5))+MSTU(5)*I
        ELSEIF(KCOL.NE.0) THEN
          K(I,1)=3
          KFLS=(3-KCOL*ISIGN(1,KFLSP(JT)))/2
          K(I,KFLS+3)=IPU
          K(IPU,6-KFLS)=MOD(K(IPU,6-KFLS),MSTU(5))+MSTU(5)*I
        ENDIF
        IF(KFLCH(JT).EQ.0) THEN
          P(I,1)=-P(MINT(83)+JT+2,1)
          P(I,2)=-P(MINT(83)+JT+2,2)
          PMS(JT)=P(I,5)**2+P(I,1)**2+P(I,2)**2
          PSYS(JT,3)=SQRT(MAX(0D0,PSYS(JT,4)**2-PMS(JT)))*(-1)**(JT-1)
          P(I,3)=PSYS(JT,3)
          P(I,4)=PSYS(JT,4)
 
C...When extra remnant parton or hadron: store extra remnant.
        ELSE
          I=I+1
          ISN(JT)=2
          DO 360 J=1,5
            K(I,J)=0
            P(I,J)=0D0
            V(I,J)=0D0
  360     CONTINUE
          K(I,1)=1
          K(I,2)=KFLCH(JT)
          K(I,3)=MINT(83)+JT
          P(I,5)=PYMASS(K(I,2))
 
C...Find parton colour connections of extra remnant.
          KCOL=KCHG(PYCOMP(KFLCH(JT)),2)
          IF(KCOL.EQ.2) THEN
            K(I,1)=3
            K(I,4)=MSTU(5)*IPU+IPU
            K(I,5)=MSTU(5)*IPU+IPU
            K(IPU,4)=MOD(K(IPU,4),MSTU(5))+MSTU(5)*I
            K(IPU,5)=MOD(K(IPU,5),MSTU(5))+MSTU(5)*I
          ELSEIF(KCOL.NE.0) THEN
            K(I,1)=3
            KFLS=(3-KCOL*ISIGN(1,KFLCH(JT)))/2
            K(I,KFLS+3)=IPU
            K(IPU,6-KFLS)=MOD(K(IPU,6-KFLS),MSTU(5))+MSTU(5)*I
          ENDIF
 
C...Relative transverse momentum when two remnants.
          LOOP=0
  370     LOOP=LOOP+1
          CALL PYPTDI(1,P(I-1,1),P(I-1,2))
          IF(IABS(MINT(10+JT)).LT.20) THEN
            P(I-1,1)=0D0
            P(I-1,2)=0D0
          ELSE
            P(I-1,1)=P(I-1,1)-0.5D0*P(MINT(83)+JT+2,1)
            P(I-1,2)=P(I-1,2)-0.5D0*P(MINT(83)+JT+2,2)
          ENDIF
          PMS(JT+2)=P(I-1,5)**2+P(I-1,1)**2+P(I-1,2)**2
          P(I,1)=-P(MINT(83)+JT+2,1)-P(I-1,1)
          P(I,2)=-P(MINT(83)+JT+2,2)-P(I-1,2)
          PMS(JT+4)=P(I,5)**2+P(I,1)**2+P(I,2)**2
 
C...Meson or baryon; photon as meson. For splitup below.
          IMB=1
          IF(MOD(MINT(10+JT)/1000,10).NE.0) IMB=2
 
C***Relative distribution for electron into two electrons. Temporary!
          IF(IABS(MINT(10+JT)).LT.20.AND.MINT(14+JT).EQ.-MINT(10+JT))
     &    THEN
            CHI(JT)=PYR(0)
 
C...Relative distribution of electron energy into electron plus parton.
          ELSEIF(IABS(MINT(10+JT)).LT.20) THEN
            XHRD=VINT(140+JT)
            XE=VINT(154+JT)
            CHI(JT)=(XE-XHRD)/(1D0-XHRD)
 
C...Relative distribution of energy for particle into two jets.
          ELSEIF(IABS(KFLCH(JT)).LE.10.OR.KFLCH(JT).EQ.21) THEN
            CHIK=PARP(92+2*IMB)
            IF(MSTP(92).LE.1) THEN
              IF(IMB.EQ.1) CHI(JT)=PYR(0)
              IF(IMB.EQ.2) CHI(JT)=1D0-SQRT(PYR(0))
            ELSEIF(MSTP(92).EQ.2) THEN
              CHI(JT)=1D0-PYR(0)**(1D0/(1D0+CHIK))
            ELSEIF(MSTP(92).EQ.3) THEN
              CUT=2D0*0.3D0/VINT(1)
  380         CHI(JT)=PYR(0)**2
              IF((CHI(JT)**2/(CHI(JT)**2+CUT**2))**0.25D0*
     &        (1D0-CHI(JT))**CHIK.LT.PYR(0)) GOTO 380
            ELSEIF(MSTP(92).EQ.4) THEN
              CUT=2D0*0.3D0/VINT(1)
              CUTR=(1D0+SQRT(1D0+CUT**2))/CUT
  390         CHIR=CUT*CUTR**PYR(0)
              CHI(JT)=(CHIR**2-CUT**2)/(2D0*CHIR)
              IF((1D0-CHI(JT))**CHIK.LT.PYR(0)) GOTO 390
            ELSE
              CUT=2D0*0.3D0/VINT(1)
              CUTA=CUT**(1D0-PARP(98))
              CUTB=(1D0+CUT)**(1D0-PARP(98))
  400         CHI(JT)=(CUTA+PYR(0)*(CUTB-CUTA))**(1D0/(1D0-PARP(98)))
              IF(((CHI(JT)+CUT)**2/(2D0*(CHI(JT)**2+CUT**2)))**
     &        (0.5D0*PARP(98))*(1D0-CHI(JT))**CHIK.LT.PYR(0)) GOTO 400
            ENDIF
 
C...Relative distribution of energy for particle into jet plus particle.
          ELSE
            IF(MSTP(94).LE.1) THEN
              IF(IMB.EQ.1) CHI(JT)=PYR(0)
              IF(IMB.EQ.2) CHI(JT)=1D0-SQRT(PYR(0))
              IF(MOD(KFLCH(JT)/1000,10).NE.0) CHI(JT)=1D0-CHI(JT)
            ELSEIF(MSTP(94).EQ.2) THEN
              CHI(JT)=1D0-PYR(0)**(1D0/(1D0+PARP(93+2*IMB)))
              IF(MOD(KFLCH(JT)/1000,10).NE.0) CHI(JT)=1D0-CHI(JT)
            ELSEIF(MSTP(94).EQ.3) THEN
              CALL PYZDIS(1,0,PMS(JT+4),ZZ)
              CHI(JT)=ZZ
            ELSE
              CALL PYZDIS(1000,0,PMS(JT+4),ZZ)
              CHI(JT)=ZZ
            ENDIF
          ENDIF
 
C...Construct total transverse mass; reject if too large.
          CHI(JT)=MAX(1D-8,MIN(1D0-1D-8,CHI(JT)))
          PMS(JT)=PMS(JT+4)/CHI(JT)+PMS(JT+2)/(1D0-CHI(JT))
          IF(PMS(JT).GT.PSYS(JT,4)**2) THEN
            IF(LOOP.LT.100) THEN
              GOTO 370
            ELSE
              MINT(51)=1
              MINT(57)=MINT(57)+1
              RETURN
            ENDIF
          ENDIF
          PSYS(JT,3)=SQRT(MAX(0D0,PSYS(JT,4)**2-PMS(JT)))*(-1)**(JT-1)
          VINT(158+JT)=CHI(JT)
 
C...Subdivide longitudinal momentum according to value selected above.
          PW1=CHI(JT)*(PSYS(JT,4)+ABS(PSYS(JT,3)))
          P(IS(JT)+1,4)=0.5D0*(PW1+PMS(JT+4)/PW1)
          P(IS(JT)+1,3)=0.5D0*(PW1-PMS(JT+4)/PW1)*(-1)**(JT-1)
          P(IS(JT),4)=PSYS(JT,4)-P(IS(JT)+1,4)
          P(IS(JT),3)=PSYS(JT,3)-P(IS(JT)+1,3)
        ENDIF
  410 CONTINUE
      N=I
 
C...Check if longitudinal boosts needed - if so pick two systems.
      PDEV=ABS(PSYS(0,4)+PSYS(1,4)+PSYS(2,4)-VINT(1))+
     &ABS(PSYS(0,3)+PSYS(1,3)+PSYS(2,3))
      IF(PDEV.LE.1D-6*VINT(1)) RETURN
      IF(ISN(1).EQ.0) THEN
        IR=0
        IL=2
      ELSEIF(ISN(2).EQ.0) THEN
        IR=1
        IL=0
      ELSEIF(VINT(143).GT.0.2D0.AND.VINT(144).GT.0.2D0) THEN
        IR=1
        IL=2
      ELSEIF(VINT(143).GT.0.2D0) THEN
        IR=1
        IL=0
      ELSEIF(VINT(144).GT.0.2D0) THEN
        IR=0
        IL=2
      ELSEIF(PMS(1)/PSYS(1,4)**2.GT.PMS(2)/PSYS(2,4)**2) THEN
        IR=1
        IL=0
      ELSE
        IR=0
        IL=2
      ENDIF
      IG=3-IR-IL
 
C...E+-pL wanted for system to be modified.
      IF((IG.EQ.1.AND.ISN(1).EQ.0).OR.(IG.EQ.2.AND.ISN(2).EQ.0)) THEN
        PPB=VINT(1)
        PNB=VINT(1)
      ELSE
        PPB=VINT(1)-(PSYS(IG,4)+PSYS(IG,3))
        PNB=VINT(1)-(PSYS(IG,4)-PSYS(IG,3))
      ENDIF
 
C...To keep x and Q2 in leptoproduction: do not count scattered lepton.
      IF(IDISXQ.EQ.1.AND.IG.NE.0) THEN
        PPB=PPB-(PSYS(0,4)+PSYS(0,3))
        PNB=PNB-(PSYS(0,4)-PSYS(0,3))
        DO 420 J=1,4
          PSYS(0,J)=0D0
  420   CONTINUE
        DO 450 I=MINT(84)+1,NS
          IF(K(I,1).GT.10) GOTO 450
          INCL=0
          IORIG=I
  430     IF(IORIG.EQ.LQOUT.OR.IORIG.EQ.LPIN+2) INCL=1
          IORIG=K(IORIG,3)
          IF(IORIG.GT.LPIN) GOTO 430
          IF(INCL.EQ.0) GOTO 450
          DO 440 J=1,4
            PSYS(0,J)=PSYS(0,J)+P(I,J)
  440     CONTINUE
  450   CONTINUE
        PMS(0)=MAX(0D0,PSYS(0,4)**2-PSYS(0,3)**2)
        PPB=PPB+(PSYS(0,4)+PSYS(0,3))
        PNB=PNB+(PSYS(0,4)-PSYS(0,3))
      ENDIF
 
C...Construct longitudinal boosts.
      DPMTB=PPB*PNB
      DPMTR=PMS(IR)
      DPMTL=PMS(IL)
      DSQLAM=SQRT(MAX(0D0,(DPMTB-DPMTR-DPMTL)**2-4D0*DPMTR*DPMTL))
      IF(DSQLAM.LE.1D-6*DPMTB) THEN
        MINT(51)=1
        MINT(57)=MINT(57)+1
        RETURN
      ENDIF
      DSQSGN=SIGN(1D0,PSYS(IR,3)*PSYS(IL,4)-PSYS(IL,3)*PSYS(IR,4))
      DRKR=(DPMTB+DPMTR-DPMTL+DSQLAM*DSQSGN)/
     &(2D0*(PSYS(IR,4)+PSYS(IR,3))*PNB)
      DRKL=(DPMTB+DPMTL-DPMTR+DSQLAM*DSQSGN)/
     &(2D0*(PSYS(IL,4)-PSYS(IL,3))*PPB)
      DBER=(DRKR**2-1D0)/(DRKR**2+1D0)
      DBEL=-(DRKL**2-1D0)/(DRKL**2+1D0)
 
C...Perform longitudinal boosts.
      IF(IR.EQ.1.AND.ISN(1).EQ.1.AND.DBER.LE.-0.99999999D0) THEN
        P(IS(1),3)=0D0
        P(IS(1),4)=SQRT(P(IS(1),5)**2+P(IS(1),1)**2+P(IS(1),2)**2)
      ELSEIF(IR.EQ.1) THEN
        CALL PYROBO(IS(1),IS(1)+ISN(1)-1,0D0,0D0,0D0,0D0,DBER)
      ELSEIF(IDISXQ.EQ.1) THEN
        DO 470 I=I1,NS
          INCL=0
          IORIG=I
  460     IF(IORIG.EQ.LQOUT.OR.IORIG.EQ.LPIN+2) INCL=1
          IORIG=K(IORIG,3)
          IF(IORIG.GT.LPIN) GOTO 460
          IF(INCL.EQ.1) CALL PYROBO(I,I,0D0,0D0,0D0,0D0,DBER)
  470   CONTINUE
      ELSE
        CALL PYROBO(I1,NS,0D0,0D0,0D0,0D0,DBER)
      ENDIF
      IF(IL.EQ.2.AND.ISN(2).EQ.1.AND.DBEL.GE.0.99999999D0) THEN
        P(IS(2),3)=0D0
        P(IS(2),4)=SQRT(P(IS(2),5)**2+P(IS(2),1)**2+P(IS(2),2)**2)
      ELSEIF(IL.EQ.2) THEN
        CALL PYROBO(IS(2),IS(2)+ISN(2)-1,0D0,0D0,0D0,0D0,DBEL)
      ELSEIF(IDISXQ.EQ.1) THEN
        DO 490 I=I1,NS
          INCL=0
          IORIG=I
  480     IF(IORIG.EQ.LQOUT.OR.IORIG.EQ.LPIN+2) INCL=1
          IORIG=K(IORIG,3)
          IF(IORIG.GT.LPIN) GOTO 480
          IF(INCL.EQ.1) CALL PYROBO(I,I,0D0,0D0,0D0,0D0,DBEL)
  490   CONTINUE
      ELSE
        CALL PYROBO(I1,NS,0D0,0D0,0D0,0D0,DBEL)
      ENDIF
 
C...Final check that energy-momentum conservation worked.
      PESUM=0D0
      PZSUM=0D0
      DO 500 I=MINT(84)+1,N
        IF(K(I,1).GT.10) GOTO 500
        PESUM=PESUM+P(I,4)
        PZSUM=PZSUM+P(I,3)
  500 CONTINUE
      PDEV=ABS(PESUM-VINT(1))+ABS(PZSUM)
      IF(PDEV.GT.1D-4*VINT(1)) THEN
        MINT(51)=1
        MINT(57)=MINT(57)+1
        RETURN
      ENDIF
 
C...Calculate rotation and boost from overall CM frame to
C...hadronic CM frame in leptoproduction.
      MINT(91)=0
      IF(MINT(82).EQ.1.AND.(MINT(43).EQ.2.OR.MINT(43).EQ.3)) THEN
        MINT(91)=1
        LESD=1
        IF(MINT(42).EQ.1) LESD=2
        LPIN=MINT(83)+3-LESD
 
C...Sum upp momenta of everything not lepton or photon to define boost.
        DO 510 J=1,4
          PSUM(J)=0D0
  510   CONTINUE
        DO 530 I=1,N
          IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 530
          IF(IABS(K(I,2)).GE.11.AND.IABS(K(I,2)).LE.20) GOTO 530
          IF(K(I,2).EQ.22) GOTO 530
          DO 520 J=1,4
            PSUM(J)=PSUM(J)+P(I,J)
  520     CONTINUE
  530   CONTINUE
        VINT(223)=-PSUM(1)/PSUM(4)
        VINT(224)=-PSUM(2)/PSUM(4)
        VINT(225)=-PSUM(3)/PSUM(4)
 
C...Boost incoming hadron to hadronic CM frame to determine rotations.
        K(N+1,1)=1
        DO 540 J=1,5
          P(N+1,J)=P(LPIN,J)
          V(N+1,J)=V(LPIN,J)
  540   CONTINUE
        CALL PYROBO(N+1,N+1,0D0,0D0,VINT(223),VINT(224),VINT(225))
        VINT(222)=-PYANGL(P(N+1,1),P(N+1,2))
        CALL PYROBO(N+1,N+1,0D0,VINT(222),0D0,0D0,0D0)
        IF(LESD.EQ.2) THEN
          VINT(221)=-PYANGL(P(N+1,3),P(N+1,1))
        ELSE
          VINT(221)=PYANGL(-P(N+1,3),P(N+1,1))
        ENDIF
      ENDIF
 
      RETURN
      END
