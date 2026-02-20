 
C*********************************************************************
 
C...PYPDFU
C...Gives electron, muon, tau, photon, pi+, neutron, proton and hyperon
C...parton distributions according to a few different parametrizations.
C...Note that what is coded is x times the probability distribution,
C...i.e. xq(x,Q2) etc.
 
      SUBROUTINE PYPDFU(KF,X,Q2,XPQ)
 
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
      COMMON/PYINT8/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
     &XPDIR(-6:6)
      COMMON/PYINT9/VXPVMD(-6:6),VXPANL(-6:6),VXPANH(-6:6),VXPDGM(-6:6)
      COMMON/PYINTM/KFIVAL(2,3),NMI(2),IMI(2,800,2),NVC(2,-6:6),
     &     XASSOC(2,-6:6,240),XPSVC(-6:6,-1:240),PVCTOT(2,-1:1),
     &     XMI(2,240),PT2MI(240),IMISEP(0:240)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/,/PYINT8/,
     &/PYINT9/,/PYINTM/
C...Local arrays.
      DIMENSION XPQ(-25:25),XPEL(-25:25),XPGA(-6:6),VXPGA(-6:6),
     &XPPI(-6:6),XPPR(-6:6),XPVAL(-6:6),PPAR(6,2)
      SAVE PPAR
 
C...Interface to PDFLIB.
      COMMON/W50513/XMIN,XMAX,Q2MIN,Q2MAX
      SAVE /W50513/
      DOUBLE PRECISION XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     &VALUE(20),XMIN,XMAX,Q2MIN,Q2MAX
      CHARACTER*20 PARM(20)
      DATA VALUE/20*0D0/,PARM/20*' '/
 
C...Data related to Schuler-Sjostrand photon distributions.
      DATA ALAMGA/0.2D0/, PMCGA/1.3D0/, PMBGA/4.6D0/
 
C...Valence PDF momentum integral parametrizations PER PARTON!
      DATA (PPAR(1,IPAR),IPAR=1,2) /0.385D0,1.60D0/
      DATA (PPAR(2,IPAR),IPAR=1,2) /0.480D0,1.56D0/
      PAVG(IFL,Q2)=PPAR(IFL,1)/(1D0+PPAR(IFL,2)*
     &LOG(LOG(MAX(Q2,1D0)/0.04D0)))
 
C...Reset parton distributions.
      MINT(92)=0
      DO 100 KFL=-25,25
        XPQ(KFL)=0D0
  100 CONTINUE
      DO 110 KFL=-6,6
        XPVAL(KFL)=0D0
  110 CONTINUE
 
C...Check x and particle species.
      IF(X.LE.0D0.OR.X.GE.1D0) THEN
        WRITE(MSTU(11),5000) X
        GOTO 9999
      ENDIF
      KFA=IABS(KF)
      IF(KFA.NE.11.AND.KFA.NE.13.AND.KFA.NE.15.AND.KFA.NE.22.AND.
     &KFA.NE.211.AND.KFA.NE.2112.AND.KFA.NE.2212.AND.KFA.NE.3122.AND.
     &KFA.NE.3112.AND.KFA.NE.3212.AND.KFA.NE.3222.AND.KFA.NE.3312.AND.
     &KFA.NE.3322.AND.KFA.NE.3334.AND.KFA.NE.111.AND.KFA.NE.321.AND.
     &KFA.NE.310.AND.KFA.NE.130) THEN
        WRITE(MSTU(11),5100) KF
        GOTO 9999
      ENDIF
 
C...Electron (or muon or tau) parton distribution call.
      IF(KFA.EQ.11.OR.KFA.EQ.13.OR.KFA.EQ.15) THEN
        CALL PYPDEL(KFA,X,Q2,XPEL)
        DO 120 KFL=-25,25
          XPQ(KFL)=XPEL(KFL)
  120   CONTINUE
 
C...Photon parton distribution call (VDM+anomalous).
      ELSEIF(KFA.EQ.22.AND.MINT(109).LE.1) THEN
        IF(MSTP(56).EQ.1.AND.MSTP(55).EQ.1) THEN
          CALL PYPDGA(X,Q2,XPGA)
          DO 130 KFL=-6,6
            XPQ(KFL)=XPGA(KFL)
  130     CONTINUE
          XPVU=4D0*(XPQ(2)-XPQ(1))/3D0
          XPVAL(1)=XPVU/4D0
          XPVAL(2)=XPVU
          XPVAL(3)=MIN(XPQ(3),XPVU/4D0)
          XPVAL(4)=MIN(XPQ(4),XPVU)
          XPVAL(5)=MIN(XPQ(5),XPVU/4D0)
          XPVAL(-1)=XPVAL(1)
          XPVAL(-2)=XPVAL(2)
          XPVAL(-3)=XPVAL(3)
          XPVAL(-4)=XPVAL(4)
          XPVAL(-5)=XPVAL(5)
        ELSEIF(MSTP(56).EQ.1.AND.MSTP(55).GE.5.AND.MSTP(55).LE.8) THEN
          Q2MX=Q2
          P2MX=0.36D0
          IF(MSTP(55).GE.7) P2MX=4.0D0
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          CALL PYGGAM(MSTP(55)-4,X,Q2MX,P2,MSTP(60),F2GAM,XPGA)
          DO 140 KFL=-6,6
            XPQ(KFL)=XPGA(KFL)
            XPVAL(KFL)=VXPDGM(KFL)
  140     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(56).EQ.1.AND.MSTP(55).GE.9.AND.MSTP(55).LE.12) THEN
          Q2MX=Q2
          P2MX=0.36D0
          IF(MSTP(55).GE.11) P2MX=4.0D0
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          CALL PYGGAM(MSTP(55)-8,X,Q2MX,P2,MSTP(60),F2GAM,XPGA)
          DO 150 KFL=-6,6
            XPQ(KFL)=XPVMD(KFL)+XPANL(KFL)+XPBEH(KFL)+XPDIR(KFL)
            XPVAL(KFL)=VXPVMD(KFL)+VXPANL(KFL)+XPBEH(KFL)+XPDIR(KFL)
  150     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(56).EQ.2) THEN
C...Call PDFLIB parton distributions.
          PARM(1)='NPTYPE'
          VALUE(1)=3
          PARM(2)='NGROUP'
          VALUE(2)=MSTP(55)/1000
          PARM(3)='NSET'
          VALUE(3)=MOD(MSTP(55),1000)
          IF(MINT(93).NE.3000000+MSTP(55)) THEN
            CALL PDFSET(PARM,VALUE)
            MINT(93)=3000000+MSTP(55)
          ENDIF
          XX=X
          QQ2=MAX(0D0,Q2MIN,Q2)
          IF(MSTP(57).EQ.0) QQ2=Q2MIN
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          IP2=MSTP(60)
          IF(MSTP(55).EQ.5004) THEN
            IF(5D0*P2.LT.QQ2.AND.
     &      QQ2.GT.0.6D0.AND.QQ2.LT.5D4.AND.
     &      P2.GE.0D0.AND.P2.LT.10D0.AND.
     &      XX.GT.1D-4.AND.XX.LT.1D0) THEN
              CALL STRUCTP(XX,QQ2,P2,IP2,UPV,DNV,USEA,DSEA,STR,CHM,
     &        BOT,TOP,GLU)
            ELSE
              UPV=0D0
              DNV=0D0
              USEA=0D0
              DSEA=0D0
              STR=0D0
              CHM=0D0
              BOT=0D0
              TOP=0D0
              GLU=0D0
            ENDIF
          ELSE
            IF(P2.LT.QQ2) THEN
              CALL STRUCTP(XX,QQ2,P2,IP2,UPV,DNV,USEA,DSEA,STR,CHM,
     &        BOT,TOP,GLU)
            ELSE
              UPV=0D0
              DNV=0D0
              USEA=0D0
              DSEA=0D0
              STR=0D0
              CHM=0D0
              BOT=0D0
              TOP=0D0
              GLU=0D0
            ENDIF
          ENDIF
          VINT(231)=Q2MIN
          XPQ(0)=GLU
          XPQ(1)=DNV
          XPQ(-1)=DNV
          XPQ(2)=UPV
          XPQ(-2)=UPV
          XPQ(3)=STR
          XPQ(-3)=STR
          XPQ(4)=CHM
          XPQ(-4)=CHM
          XPQ(5)=BOT
          XPQ(-5)=BOT
          XPQ(6)=TOP
          XPQ(-6)=TOP
          XPVU=4D0*(XPQ(2)-XPQ(1))/3D0
          XPVAL(1)=XPVU/4D0
          XPVAL(2)=XPVU
          XPVAL(3)=MIN(XPQ(3),XPVU/4D0)
          XPVAL(4)=MIN(XPQ(4),XPVU)
          XPVAL(5)=MIN(XPQ(5),XPVU/4D0)
          XPVAL(-1)=XPVAL(1)
          XPVAL(-2)=XPVAL(2)
          XPVAL(-3)=XPVAL(3)
          XPVAL(-4)=XPVAL(4)
          XPVAL(-5)=XPVAL(5)
        ELSE
          WRITE(MSTU(11),5200) KF,MSTP(56),MSTP(55)
        ENDIF
 
C...Pion/gammaVDM parton distribution call.
      ELSEIF(KFA.EQ.211.OR.KFA.EQ.111.OR.KFA.EQ.321.OR.KFA.EQ.130.OR.
     &KFA.EQ.310.OR.(KFA.EQ.22.AND.MINT(109).EQ.2)) THEN
        IF(KFA.EQ.22.AND.MSTP(56).EQ.1.AND.MSTP(55).GE.5.AND.
     &  MSTP(55).LE.12) THEN
          ISET=1+MOD(MSTP(55)-1,4)
          Q2MX=Q2
          P2MX=0.36D0
          IF(ISET.GE.3) P2MX=4.0D0
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          CALL PYGGAM(ISET,X,Q2MX,P2,MSTP(60),F2GAM,XPGA)
          DO 160 KFL=-6,6
            XPQ(KFL)=XPVMD(KFL)
            XPVAL(KFL)=VXPVMD(KFL)
  160     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(54).EQ.1.AND.MSTP(53).GE.1.AND.MSTP(53).LE.3) THEN
          CALL PYPDPI(X,Q2,XPPI)
          DO 170 KFL=-6,6
            XPQ(KFL)=XPPI(KFL)
  170     CONTINUE
          XPVAL(2)=XPQ(2)-XPQ(-2)
          XPVAL(-1)=XPQ(-1)-XPQ(1)
        ELSEIF(MSTP(54).EQ.2) THEN
C...Call PDFLIB parton distributions.
          PARM(1)='NPTYPE'
          VALUE(1)=2
          PARM(2)='NGROUP'
          VALUE(2)=MSTP(53)/1000
          PARM(3)='NSET'
          VALUE(3)=MOD(MSTP(53),1000)
          IF(MINT(93).NE.2000000+MSTP(53)) THEN
            CALL PDFSET(PARM,VALUE)
            MINT(93)=2000000+MSTP(53)
          ENDIF
          XX=X
          QQ=SQRT(MAX(0D0,Q2MIN,Q2))
          IF(MSTP(57).EQ.0) QQ=SQRT(Q2MIN)
          CALL STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
          VINT(231)=Q2MIN
          XPQ(0)=GLU
          XPQ(1)=DSEA
          XPQ(-1)=UPV+DSEA
          XPQ(2)=UPV+USEA
          XPQ(-2)=USEA
          XPQ(3)=STR
          XPQ(-3)=STR
          XPQ(4)=CHM
          XPQ(-4)=CHM
          XPQ(5)=BOT
          XPQ(-5)=BOT
          XPQ(6)=TOP
          XPQ(-6)=TOP
          XPVAL(2)=UPV
          XPVAL(-1)=UPV
        ELSE
          WRITE(MSTU(11),5200) KF,MSTP(54),MSTP(53)
        ENDIF
 
C...Anomalous photon parton distribution call.
      ELSEIF(KFA.EQ.22.AND.MINT(109).EQ.3) THEN
        Q2MX=Q2
        P2MX=PARP(15)**2
        IF(MSTP(56).EQ.1.AND.MSTP(55).LE.8) THEN
          IF(MSTP(55).EQ.5.OR.MSTP(55).EQ.6) P2MX=0.36D0
          IF(MSTP(55).EQ.7.OR.MSTP(55).EQ.8) P2MX=4.0D0
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          CALL PYGGAM(MSTP(55)-4,X,Q2MX,P2,MSTP(60),F2GM,XPGA)
          DO 180 KFL=-6,6
            XPQ(KFL)=XPANL(KFL)+XPANH(KFL)
            XPVAL(KFL)=VXPANL(KFL)+VXPANH(KFL)
  180     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(56).EQ.1) THEN
          IF(MSTP(55).EQ.9.OR.MSTP(55).EQ.10) P2MX=0.36D0
          IF(MSTP(55).EQ.11.OR.MSTP(55).EQ.12) P2MX=4.0D0
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          P2=0D0
          IF(VINT(120).LT.0D0) P2=VINT(120)**2
          CALL PYGGAM(MSTP(55)-8,X,Q2MX,P2,MSTP(60),F2GM,XPGA)
          DO 190 KFL=-6,6
            XPQ(KFL)=MAX(0D0,XPANL(KFL)+XPBEH(KFL)+XPDIR(KFL))
            XPVAL(KFL)=MAX(0D0,VXPANL(KFL)+XPBEH(KFL)+XPDIR(KFL))
  190     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(56).EQ.2) THEN
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          CALL PYGANO(0,X,Q2MX,P2MX,ALAMGA,XPGA,VXPGA)
          DO 200 KFL=-6,6
            XPQ(KFL)=XPGA(KFL)
            XPVAL(KFL)=VXPGA(KFL)
  200     CONTINUE
          VINT(231)=P2MX
        ELSEIF(MSTP(55).GE.1.AND.MSTP(55).LE.5) THEN
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          CALL PYGVMD(0,MSTP(55),X,Q2MX,P2MX,PARP(1),XPGA,VXPGA)
          DO 210 KFL=-6,6
            XPQ(KFL)=XPGA(KFL)
            XPVAL(KFL)=VXPGA(KFL)
  210     CONTINUE
          VINT(231)=P2MX
        ELSE
  220     RKF=11D0*PYR(0)
          KFR=1
          IF(RKF.GT.1D0) KFR=2
          IF(RKF.GT.5D0) KFR=3
          IF(RKF.GT.6D0) KFR=4
          IF(RKF.GT.10D0) KFR=5
          IF(KFR.EQ.4.AND.Q2.LT.PMCGA**2) GOTO 220
          IF(KFR.EQ.5.AND.Q2.LT.PMBGA**2) GOTO 220
          IF(MSTP(57).EQ.0) Q2MX=P2MX
          CALL PYGVMD(0,KFR,X,Q2MX,P2MX,PARP(1),XPGA,VXPGA)
          DO 230 KFL=-6,6
            XPQ(KFL)=XPGA(KFL)
            XPVAL(KFL)=VXPGA(KFL)
  230     CONTINUE
          VINT(231)=P2MX
        ENDIF
 
C...Proton parton distribution call.
      ELSE
        IF(MSTP(52).EQ.1.AND.MSTP(51).GE.1.AND.MSTP(51).LE.20) THEN
          CALL PYPDPR(X,Q2,XPPR)
          DO 240 KFL=-6,6
            XPQ(KFL)=XPPR(KFL)
  240     CONTINUE
C...Force VAL > 0 (can be < 0 at very small Q2 and small x apparently)
          XPVAL(1)=MAX(0D0,XPQ(1)-XPQ(-1))
          XPVAL(2)=MAX(0D0,XPQ(2)-XPQ(-2))
        ELSEIF(MSTP(52).EQ.2) THEN
C...Call PDFLIB parton distributions.
          PARM(1)='NPTYPE'
          VALUE(1)=1
          PARM(2)='NGROUP'
          VALUE(2)=MSTP(51)/1000
          PARM(3)='NSET'
          VALUE(3)=MOD(MSTP(51),1000)
          IF(MINT(93).NE.1000000+MSTP(51)) THEN
            CALL PDFSET(PARM,VALUE)
            MINT(93)=1000000+MSTP(51)
          ENDIF
          XX=X
          QQ=SQRT(MAX(0D0,Q2MIN,Q2))
          IF(MSTP(57).EQ.0) QQ=SQRT(Q2MIN)
          CALL STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
          VINT(231)=Q2MIN
          XPQ(0)=GLU
          XPQ(1)=DNV+DSEA
          XPQ(-1)=DSEA
          XPQ(2)=UPV+USEA
          XPQ(-2)=USEA
          XPQ(3)=STR
          XPQ(-3)=STR
          XPQ(4)=CHM
          XPQ(-4)=CHM
          XPQ(5)=BOT
          XPQ(-5)=BOT
          XPQ(6)=TOP
          XPQ(-6)=TOP
          XPVAL(1)=DNV
          XPVAL(2)=UPV
        ELSE
          WRITE(MSTU(11),5200) KF,MSTP(52),MSTP(51)
        ENDIF
      ENDIF
 
C...Isospin average for pi0/gammaVDM.
      IF(KFA.EQ.111.OR.(KFA.EQ.22.AND.MINT(109).EQ.2)) THEN
        IF(KFA.EQ.22.AND.MSTP(55).GE.5.AND.MSTP(55).LE.12) THEN
          XPV=XPQ(2)-XPQ(1)
          XPQ(2)=XPQ(1)
          XPQ(-2)=XPQ(-1)
        ELSE
          XPS=0.5D0*(XPQ(1)+XPQ(-2))
          XPV=0.5D0*(XPQ(2)+XPQ(-1))-XPS
          XPQ(2)=XPS
          XPQ(-1)=XPS
        ENDIF
        XPVL=0.5D0*(XPVAL(1)+XPVAL(2)+XPVAL(-1)+XPVAL(-2))+
     &  XPVAL(3)+XPVAL(4)+XPVAL(5)
        DO 250 KFL=-6,6
          XPVAL(KFL)=0D0
  250   CONTINUE
        IF(KFA.EQ.22.AND.MINT(105).LE.223) THEN
          XPQ(1)=XPQ(1)+0.2D0*XPV
          XPQ(2)=XPQ(2)+0.8D0*XPV
          XPVAL(1)=0.2D0*XPVL
          XPVAL(2)=0.8D0*XPVL
        ELSEIF(KFA.EQ.22.AND.MINT(105).EQ.333) THEN
          XPQ(3)=XPQ(3)+XPV
          XPVAL(3)=XPVL
        ELSEIF(KFA.EQ.22.AND.MINT(105).EQ.443) THEN
          XPQ(4)=XPQ(4)+XPV
          XPVAL(4)=XPVL
          IF(MSTP(55).GE.9) THEN
            DO 260 KFL=-6,6
              XPQ(KFL)=0D0
  260       CONTINUE
          ENDIF
        ELSE
          XPQ(1)=XPQ(1)+0.5D0*XPV
          XPQ(2)=XPQ(2)+0.5D0*XPV
          XPVAL(1)=0.5D0*XPVL
          XPVAL(2)=0.5D0*XPVL
        ENDIF
        DO 270 KFL=1,6
          XPQ(-KFL)=XPQ(KFL)
          XPVAL(-KFL)=XPVAL(KFL)
  270   CONTINUE
 
C...Rescale for gammaVDM by effective gamma -> rho coupling.
C+++Do not rescale?
        IF(KFA.EQ.22.AND.MINT(109).EQ.2.AND..NOT.(MSTP(56).EQ.1
     &  .AND.MSTP(55).GE.5.AND.MSTP(55).LE.12)) THEN
          DO 280 KFL=-6,6
            XPQ(KFL)=VINT(281)*XPQ(KFL)
            XPVAL(KFL)=VINT(281)*XPVAL(KFL)
  280     CONTINUE
          VINT(232)=VINT(281)*XPV
        ENDIF
 
C...Simple recipes for kaons.
      ELSEIF(KFA.EQ.321) THEN
        XPQ(-3)=XPQ(-3)+XPQ(-1)-XPQ(1)
        XPQ(-1)=XPQ(1)
        XPVAL(-3)=XPVAL(-1)
        XPVAL(-1)=0D0
      ELSEIF(KFA.EQ.130.OR.KFA.EQ.310) THEN
        XPS=0.5D0*(XPQ(1)+XPQ(-2))
        XPV=0.5D0*(XPQ(2)+XPQ(-1))-XPS
        XPQ(2)=XPS
        XPQ(-1)=XPS
        XPQ(1)=XPQ(1)+0.5D0*XPV
        XPQ(-1)=XPQ(-1)+0.5D0*XPV
        XPQ(3)=XPQ(3)+0.5D0*XPV
        XPQ(-3)=XPQ(-3)+0.5D0*XPV
        XPV=0.5D0*(XPVAL(2)+XPVAL(-1))
        XPVAL(2)=0D0
        XPVAL(-1)=0D0
        XPVAL(1)=0.5D0*XPV
        XPVAL(-1)=0.5D0*XPV
        XPVAL(3)=0.5D0*XPV
        XPVAL(-3)=0.5D0*XPV
 
C...Isospin conjugation for neutron.
      ELSEIF(KFA.EQ.2112) THEN
        XPSV=XPQ(1)
        XPQ(1)=XPQ(2)
        XPQ(2)=XPSV
        XPSV=XPQ(-1)
        XPQ(-1)=XPQ(-2)
        XPQ(-2)=XPSV
        XPSV=XPVAL(1)
        XPVAL(1)=XPVAL(2)
        XPVAL(2)=XPSV
 
C...Simple recipes for hyperon (average valence parton distribution).
      ELSEIF(KFA.EQ.3122.OR.KFA.EQ.3112.OR.KFA.EQ.3212.OR.KFA.EQ.3222
     &  .OR.KFA.EQ.3312.OR.KFA.EQ.3322.OR.KFA.EQ.3334) THEN
        XPV=(XPQ(1)+XPQ(2)-XPQ(-1)-XPQ(-2))/3D0
        XPS=0.5D0*(XPQ(-1)+XPQ(-2))
        XPQ(1)=XPS
        XPQ(2)=XPS
        XPQ(-1)=XPS
        XPQ(-2)=XPS
        XPQ(KFA/1000)=XPQ(KFA/1000)+XPV
        XPQ(MOD(KFA/100,10))=XPQ(MOD(KFA/100,10))+XPV
        XPQ(MOD(KFA/10,10))=XPQ(MOD(KFA/10,10))+XPV
        XPV=(XPVAL(1)+XPVAL(2))/3D0
        XPVAL(1)=0D0
        XPVAL(2)=0D0
        XPVAL(KFA/1000)=XPVAL(KFA/1000)+XPV
        XPVAL(MOD(KFA/100,10))=XPVAL(MOD(KFA/100,10))+XPV
        XPVAL(MOD(KFA/10,10))=XPVAL(MOD(KFA/10,10))+XPV
      ENDIF
 
C...Charge conjugation for antiparticle.
      IF(KF.LT.0) THEN
        DO 290 KFL=1,25
          IF(KFL.EQ.21.OR.KFL.EQ.22.OR.KFL.EQ.23.OR.KFL.EQ.25) GOTO 290
          XPSV=XPQ(KFL)
          XPQ(KFL)=XPQ(-KFL)
          XPQ(-KFL)=XPSV
  290   CONTINUE
        DO 300 KFL=1,6
          XPSV=XPVAL(KFL)
          XPVAL(KFL)=XPVAL(-KFL)
          XPVAL(-KFL)=XPSV
  300  CONTINUE
      ENDIF
 
C...MULTIPLE INTERACTIONS - PDF RESHAPING.
C...Set side.
      JS=MINT(30)
C...Only reshape PDFs for the non-first interactions;
C...But need valence/sea separation already from first interaction.
      IF ((JS.EQ.1.OR.JS.EQ.2).AND.MINT(35).GE.2) THEN
        KFVSEL=KFIVAL(JS,1)
C...If valence quark kicked out of pi0 or gamma then that decides
C...whether we should consider state as d dbar, u ubar, s sbar, etc.
        IF(KFVSEL.NE.0.AND.(KFA.EQ.111.OR.KFA.EQ.22)) THEN
          XPVL=0D0
          DO 310 KFL=1,6
            XPVL=XPVL+XPVAL(KFL)
            XPQ(KFL)=MAX(0D0,XPQ(KFL)-XPVAL(KFL))
            XPVAL(KFL)=0D0
  310     CONTINUE
          XPQ(IABS(KFVSEL))=XPQ(IABS(KFVSEL))+XPVL
          XPVAL(IABS(KFVSEL))=XPVL
          DO 320 KFL=1,6
            XPQ(-KFL)=XPQ(KFL)
            XPVAL(-KFL)=XPVAL(KFL)
  320     CONTINUE
 
C...If valence quark kicked out of K0S or K0S then that decides whether
C...we should consider state as d sbar or s dbar.
        ELSEIF(KFVSEL.NE.0.AND.(KFA.EQ.130.OR.KFA.EQ.310)) THEN
          KFS=1
          IF(KFVSEL.EQ.-1.OR.KFVSEL.EQ.3) KFS=-1
          XPQ(KFS)=XPQ(KFS)+XPVAL(-KFS)
          XPVAL(KFS)=XPVAL(KFS)+XPVAL(-KFS)
          XPQ(-KFS)=MAX(0D0,XPQ(-KFS)-XPVAL(-KFS))
          XPVAL(-KFS)=0D0
          KFS=-3*KFS
          XPQ(KFS)=XPQ(KFS)+XPVAL(-KFS)
          XPVAL(KFS)=XPVAL(KFS)+XPVAL(-KFS)
          XPQ(-KFS)=MAX(0D0,XPQ(-KFS)-XPVAL(-KFS))
          XPVAL(-KFS)=0D0
        ENDIF
 
C...XPQ distributions are nominal for a (signed) beam particle
C...of KF type, with 1-Sum(x_prev) rescaled to 1.
        CMPFAC=1D0
        NRESC=0
 345    NRESC=NRESC+1
        PVCTOT(JS,-1)=0D0
        PVCTOT(JS, 0)=0D0
        PVCTOT(JS, 1)=0D0
        DO 350 IFL=-6,6
          IF(IFL.EQ.0) GOTO 350
 
C...Count up number of original IFL valence quarks.
          IVORG=0
          IF(KFIVAL(JS,1).EQ.IFL) IVORG=IVORG+1
          IF(KFIVAL(JS,2).EQ.IFL) IVORG=IVORG+1
          IF(KFIVAL(JS,3).EQ.IFL) IVORG=IVORG+1
C...For pi0/gamma/K0S/K0L without valence flavour decided yet, here
C...bookkeep as if d dbar (for total momentum sum in valence sector).
          IF(KFIVAL(JS,1).EQ.0.AND.IABS(IFL).EQ.1) IVORG=1
C...Count down number of remaining IFL valence quarks. Skip current
C...interaction initiator.
          IVREM=IVORG
          DO 330 I1=1,NMI(JS)
            IF (I1.EQ.MINT(36)) GOTO 330
            IF (K(IMI(JS,I1,1),2).EQ.IFL.AND.IMI(JS,I1,2).EQ.0)
     &           IVREM=IVREM-1
  330     CONTINUE
 
C...Separate out original VALENCE and SEA content.
          VAL=XPVAL(IFL)
          SEA=MAX(0D0,XPQ(IFL)-VAL)
          XPSVC(IFL,0)=VAL
          XPSVC(IFL,-1)=SEA
 
C...Rescale valence content if changed.
          IF (IVORG.NE.0.AND.IVREM.NE.IVORG) XPSVC(IFL,0)=
     &    (VAL*IVREM)/IVORG
 
C...Momentum integrals of original and removed valence quarks.
          IF(IVORG.NE.0) THEN
C...For p/n/pbar/nbar beams can split into d_val and u_val.
C...Isospin conjugation for neutrons
            IF(KFA.EQ.2212.OR.KFA.EQ.2112) THEN
              IAFLP=IABS(IFL)
              IF (KFA.EQ.2112) IAFLP=3-IAFLP
              VPAVG=PAVG(IAFLP,Q2)
C...For other baryons average d_val and u_val, like for PDFs.
            ELSEIF(KFA.GT.1000) THEN
              VPAVG=(PAVG(1,Q2)+2D0*PAVG(2,Q2))/3D0
C...For mesons and photon average d_val and u_val and scale by 3/2.
C...Very crude, especially for photon.
            ELSE
              VPAVG=0.5D0*(PAVG(1,Q2)+2D0*PAVG(2,Q2))
            ENDIF
            PVCTOT(JS,-1)=PVCTOT(JS,-1)+IVORG*VPAVG
            PVCTOT(JS, 0)=PVCTOT(JS, 0)+(IVORG-IVREM)*VPAVG
          ENDIF
 
C...Now add companions (at X with partner having been at Z=XASSOC).
C...NOTE: due to the assumed simple x scaling, the partner was at what
C...corresponds to a higher Z than XASSOC, if there were intermediate
C...scatterings. Nothing done about that for the moment.
          DO 340 IVC=1,NVC(JS,IFL)
C...Skip companions that have been kicked out
            IF (XASSOC(JS,IFL,IVC).LE.0D0) THEN
              XPSVC(IFL,IVC)=0D0
              GOTO 340
            ELSE
C...Momentum fraction of the partner quark.
C...Use rescaled YS = XS/(1-Sum_rest) where X and XS are not in "rest".
              XS=XASSOC(JS,IFL,IVC)
              XREM=VINT(142+JS)
              YS=XS/(XREM+XS)
C...Momentum fraction of the companion quark.
C...Rescale from X = x/XREM to Y = x/(1-Sum_rest) -> factor (1-YS).
              Y=X*(1D0-YS)
              XPSVC(IFL,IVC)=PYFCMP(Y/CMPFAC,YS/CMPFAC,MSTP(87))
C...Add to momentum sum, with rescaling compensation factor.
              XCFAC=(XREM+XS)/XREM*CMPFAC
              PVCTOT(JS,1)=PVCTOT(JS,1)+XCFAC*PYPCMP(YS/CMPFAC,MSTP(87))
            ENDIF
  340     CONTINUE
  350   CONTINUE
 
C...Wait until all flavours treated, then rescale seas and gluon.
        XPSVC(0,-1)=XPQ(0)
        XPSVC(0,0)=0D0
        RSFAC=1D0+(PVCTOT(JS,0)-PVCTOT(JS,1))/(1D0-PVCTOT(JS,-1))
        IF (RSFAC.LE.0D0) THEN
C...First calculate factor needed to exactly restore pz cons.
          IF (NRESC.EQ.1) CMPFAC =
     &         (1D0-(PVCTOT(JS,-1)-PVCTOT(JS,0)))/PVCTOT(JS,1)
C...Add a bit of headroom
          CMPFAC=0.99*CMPFAC
C...Try a few times if more headroom is needed, then print error message.
          IF (NRESC.LE.10) GOTO 345
          CALL PYERRM(15,
     &         '(PYPDFU:) Negative reshaping factor persists!')
          WRITE(MSTU(11),5300) (PVCTOT(JS,ITMP),ITMP=-1,1), RSFAC
          RSFAC=0D0
        ENDIF
        DO 370 IFL=-6,6
          XPSVC(IFL,-1)=RSFAC*XPSVC(IFL,-1)
C...Also store resulting distributions in XPQ
          XPQ(IFL)=0D0
          DO 360 ISVC=-1,NVC(JS,IFL)
            XPQ(IFL)=XPQ(IFL)+XPSVC(IFL,ISVC)
  360     CONTINUE
  370   CONTINUE
C...Save companion reweighting factor for PYPTIS.
        VINT(140)=CMPFAC
      ENDIF
 
 
C...Allow gluon also in position 21.
      XPQ(21)=XPQ(0)
 
C...Check positivity and reset above maximum allowed flavour.
      DO 380 KFL=-25,25
        XPQ(KFL)=MAX(0D0,XPQ(KFL))
        IF(IABS(KFL).GT.MSTP(58).AND.IABS(KFL).LE.8) XPQ(KFL)=0D0
  380 CONTINUE
 
C...Formats for error printouts.
 5000 FORMAT(' Error: x value outside physical range; x =',1P,D12.3)
 5100 FORMAT(' Error: illegal particle code for parton distribution;',
     &' KF =',I5)
 5200 FORMAT(' Error: unknown parton distribution; KF, library, set =',
     &3I5)
 5300 FORMAT(' Original valence momentum fraction : ',F6.3/
     &       ' Removed valence momentum fraction  : ',F6.3/
     &       ' Added companion momentum fraction  : ',F6.3/
     &       ' Resulting rescale factor           : ',F6.3)
 
C...Reset side pointer and return
 9999 MINT(30)=0
 
      RETURN
      END
