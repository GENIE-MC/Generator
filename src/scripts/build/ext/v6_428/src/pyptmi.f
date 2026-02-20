 
C*********************************************************************
 
C...PYPTMI
C...Handles the generation of additional interactions in the new
C...multiple interactions framework.
C...MODE=-1 : Initalize MI from scratch.
C...MODE= 0 : Generate trial interaction. Start at PT2NOW, solve
C...         Sudakov for PT2, abort if below PT2CUT.
C...MODE= 1 : Accept interaction at PT2NOW and store variables.
C...MODE= 2 : Decide sea/val/cmp for kicked-out quark at PT2NOW
C...PT2NOW  : Starting (max) PT2 scale for evolution.
C...PT2CUT  : Lower limit for evolution.
C...PT2     : Result of evolution. Generated PT2 for trial interaction.
C...IFAIL   : Status return code.
C...         = 0: All is well.
C...         < 0: Phase space exhausted, generation to be terminated.
C...         > 0: Additional interaction vetoed, but continue evolution.
 
      SUBROUTINE PYPTMI(MODE,PT2NOW,PT2CUT,PT2,IFAIL)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement for maximum size of showers.
      PARAMETER (MAXNUR=1000)
C...Commonblocks.
      COMMON/PYPART/NPART,NPARTD,IPART(MAXNUR),PTPART(MAXNUR)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      COMMON/PYINTM/KFIVAL(2,3),NMI(2),IMI(2,800,2),NVC(2,-6:6),
     &     XASSOC(2,-6:6,240),XPSVC(-6:6,-1:240),PVCTOT(2,-1:1),
     &     XMI(2,240),PT2MI(240),IMISEP(0:240)
      COMMON/PYISMX/MIMX,JSMX,KFLAMX,KFLCMX,KFBEAM(2),NISGEN(2,240),
     &     PT2MX,PT2AMX,ZMX,RM2CMX,Q2BMX,PHIMX
      COMMON/PYCTAG/NCT,MCT(4000,2)
C...Local arrays and saved variables.
      DIMENSION WDTP(0:400),WDTE(0:400,0:5),XPQ(-25:25)
 
      SAVE /PYPART/,/PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYPARS/,
     &     /PYINT1/,/PYINT2/,/PYINT3/,/PYINT5/,/PYINT7/,/PYINTM/,
     &     /PYISMX/,/PYCTAG/
      SAVE NCHN,XT2FAC,SIGS
 
      IFAIL=0
C...Set MI subprocess = QCD 2 -> 2.
      ISUB=96
 
C----------------------------------------------------------------------
C...MODE=-1: Initialize from scratch
      IF (MODE.EQ.-1) THEN
C...Initialize PT2 array.
        PT2MI(1)=VINT(54)
C...Initialize list of incoming beams and partons from two sides.
        DO 110 JS=1,2
          DO 100 MI=1,240
            IMI(JS,MI,1)=0
            IMI(JS,MI,2)=0
  100     CONTINUE
          NMI(JS)=1
          IMI(JS,1,1)=MINT(84)+JS
          IMI(JS,1,2)=0
          XMI(JS,1)=VINT(40+JS)
C...Rescale x values to fractions of photon energy.
          IF(MINT(18+JS).EQ.1) XMI(JS,1)=VINT(40+JS)/VINT(154+JS)
C...Hard reset: hard interaction initiators motherless by definition.
          K(MINT(84)+JS,3)=2+JS
          K(MINT(84)+JS,4)=MOD(K(MINT(84)+JS,4),MSTU(5))
          K(MINT(84)+JS,5)=MOD(K(MINT(84)+JS,5),MSTU(5))
  110   CONTINUE
        IMISEP(0)=MINT(84)
        IMISEP(1)=N
        IF (MOD(MSTP(81),10).GE.1) THEN
          IF(MSTP(82).LE.1) THEN
            SIGRAT=XSEC(ISUB,1)/MAX(1D-10,VINT(315)*VINT(316)*SIGT(0,0
     &           ,5))
            IF(MINT(141).NE.0.OR.MINT(142).NE.0) SIGRAT=SIGRAT*
     &           VINT(317)/(VINT(318)*VINT(320))
            XT2FAC=SIGRAT*VINT(149)/(1D0-VINT(149))
          ELSE
            XT2FAC=VINT(146)*VINT(148)*XSEC(ISUB,1)/
     &           MAX(1D-10,SIGT(0,0,5))*VINT(149)*(1D0+VINT(149))
          ENDIF
        ENDIF
C...Zero entries relating to scatterings beyond the first.
        DO 120 MI=2,240
          IMI(1,MI,1)=0
          IMI(2,MI,1)=0
          IMI(1,MI,2)=0
          IMI(2,MI,2)=0
          IMISEP(MI)=IMISEP(1)
          PT2MI(MI)=0D0
          XMI(1,MI)=0D0
          XMI(2,MI)=0D0
  120   CONTINUE
C...Initialize factors for PDF reshaping.
        DO 140 JS=1,2
          KFBEAM(JS)=MINT(10+JS)
          IF(MINT(18+JS).EQ.1) KFBEAM(JS)=22
          KFABM=IABS(KFBEAM(JS))
          KFSBM=ISIGN(1,KFBEAM(JS))
 
C...Zero flavour content of incoming beam particle.
          KFIVAL(JS,1)=0
          KFIVAL(JS,2)=0
          KFIVAL(JS,3)=0
C...  Flavour content of baryon.
          IF(KFABM.GT.1000) THEN
            KFIVAL(JS,1)=KFSBM*MOD(KFABM/1000,10)
            KFIVAL(JS,2)=KFSBM*MOD(KFABM/100,10)
            KFIVAL(JS,3)=KFSBM*MOD(KFABM/10,10)
C...  Flavour content of pi+-, K+-.
          ELSEIF(KFABM.EQ.211) THEN
            KFIVAL(JS,1)=KFSBM*2
            KFIVAL(JS,2)=-KFSBM
          ELSEIF(KFABM.EQ.321) THEN
            KFIVAL(JS,1)=-KFSBM*3
            KFIVAL(JS,2)=KFSBM*2
C...  Flavour content of pi0, gamma, K0S, K0L not defined yet.
          ENDIF
 
C...Zero initial valence and companion content.
          DO 130 IFL=-6,6
            NVC(JS,IFL)=0
  130     CONTINUE
  140   CONTINUE
C...Set up colour line tags starting from hard interaction initiators.
        NCT=0
C...Reset colour tag array and colour processing flags.
        DO 150 I=IMISEP(0)+1,N
          MCT(I,1)=0
          MCT(I,2)=0
          K(I,4)=MOD(K(I,4),MSTU(5)**2)
          K(I,5)=MOD(K(I,5),MSTU(5)**2)
  150   CONTINUE
C...  Consider each side in turn.
        DO 170 JS=1,2
          I1=IMI(JS,1,1)
          I2=IMI(3-JS,1,1)
          DO 160 JCS=4,5
            IF (K(I1,2).NE.21.AND.(9-2*JCS).NE.ISIGN(1,K(I1,2)))
     &           GOTO 160
            IF (K(I1,JCS)/MSTU(5)**2.NE.0) GOTO 160
            KCS=JCS
            CALL PYCTTR(I1,KCS,I2)
            IF(MINT(51).NE.0) RETURN
  160     CONTINUE
  170   CONTINUE
 
C...Range checking for companion quark pdf large-x param.
        IF (MSTP(87).LT.0) THEN
          CALL PYERRM(19,'(PYPTMI:) MSTP(87) out of range. Forced'//
     &         ' MSTP(87)=0')
          MSTP(87)=0
        ELSEIF (MSTP(87).GT.4) THEN
          CALL PYERRM(19,'(PYPTMI:) MSTP(87) out of range. Forced'//
     &         ' MSTP(87)=4')
          MSTP(87)=4
        ENDIF
 
C----------------------------------------------------------------------
C...MODE=0: Generate trial interaction. Return codes:
C...IFAIL < 0: Phase space exhausted, generation to be terminated.
C...IFAIL = 0: Additional interaction generated at PT2.
C...IFAIL > 0: Additional interaction vetoed, but continue evolution.
      ELSEIF (MODE.EQ.0) THEN
C...Abolute MI max scale = VINT(62)
        XT2=4D0*MIN(PT2NOW,VINT(62))/VINT(2)
  180   IF(MSTP(82).LE.1) THEN
          XT2=XT2FAC*XT2/(XT2FAC-XT2*LOG(PYR(0)))
          IF(XT2.LT.VINT(149)) IFAIL=-2
        ELSE
          IF(XT2.LE.0.01001D0*VINT(149)) THEN
            IFAIL=-3
          ELSE
            XT2=XT2FAC*(XT2+VINT(149))/(XT2FAC-(XT2+VINT(149))*
     &           LOG(PYR(0)))-VINT(149)
          ENDIF
        ENDIF
C...Also exit if below lower limit or if higher trial branching
C...already found.
        PT2=0.25D0*VINT(2)*XT2
        IF (PT2.LE.PT2CUT) IFAIL=-4
        IF (PT2.LE.PT2MX) IFAIL=-5
        IF (IFAIL.NE.0) THEN
          PT2=0D0
          RETURN
        ENDIF
        IF(MSTP(82).GE.2) PT2=MAX(0.25D0*VINT(2)*0.01D0*VINT(149),PT2)
        VINT(25)=4D0*PT2/VINT(2)
        XT2=VINT(25)
 
C...Choose tau and y*. Calculate cos(theta-hat).
        IF(PYR(0).LE.COEF(ISUB,1)) THEN
          TAUT=(2D0*(1D0+SQRT(1D0-XT2))/XT2-1D0)**PYR(0)
          TAU=XT2*(1D0+TAUT)**2/(4D0*TAUT)
        ELSE
          TAU=XT2*(1D0+TAN(PYR(0)*ATAN(SQRT(1D0/XT2-1D0)))**2)
        ENDIF
        VINT(21)=TAU
C...New: require shat > 1.
        IF(TAU*VINT(2).LT.1D0) GOTO 180
        CALL PYKLIM(2)
        RYST=PYR(0)
        MYST=1
        IF(RYST.GT.COEF(ISUB,8)) MYST=2
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)) MYST=3
        CALL PYKMAP(2,MYST,PYR(0))
        VINT(23)=SQRT(MAX(0D0,1D0-XT2/TAU))*(-1)**INT(1.5D0+PYR(0))
 
C...Check that x not used up. Accept or reject kinematical variables.
        X1M=SQRT(TAU)*EXP(VINT(22))
        X2M=SQRT(TAU)*EXP(-VINT(22))
        IF(VINT(143)-X1M.LT.0.01D0.OR.VINT(144)-X2M.LT.0.01D0) GOTO 180
        VINT(71)=0.5D0*VINT(1)*SQRT(XT2)
        NCHN=0
        CALL PYSIGH(NCHN,SIGS)
        IF(MINT(141).NE.0.OR.MINT(142).NE.0) SIGS=SIGS*VINT(320)
        IF(SIGS.LT.XSEC(ISUB,1)*PYR(0)) GOTO 180
        IF(MINT(141).NE.0.OR.MINT(142).NE.0) SIGS=SIGS/VINT(320)
 
C...Save if highest PT so far.
        IF (PT2.GT.PT2MX) THEN
          JSMX=0
          MIMX=MINT(31)+1
          PT2MX=PT2
        ENDIF
 
C----------------------------------------------------------------------
C...MODE=1: Generate and save accepted scattering.
      ELSEIF (MODE.EQ.1) THEN
        PT2=PT2NOW
C...Reset K, P, V, and MCT vectors.
        DO 200 I=N+1,N+4
          DO 190 J=1,5
            K(I,J)=0
            P(I,J)=0D0
            V(I,J)=0D0
  190     CONTINUE
          MCT(I,1)=0
          MCT(I,2)=0
  200   CONTINUE
 
        NTRY=0
C...Choose flavour of reacting partons (and subprocess).
  210   NTRY=NTRY+1
        IF (NTRY.GT.50) THEN
          CALL PYERRM(9,'(PYPTMI:) Unable to generate additional '
     &               //'interaction. Giving up!')
          MINT(51)=1
          RETURN
        ENDIF
        RSIGS=SIGS*PYR(0)
        DO 220 ICHN=1,NCHN
          KFL1=ISIG(ICHN,1)
          KFL2=ISIG(ICHN,2)
          ICONMI=ISIG(ICHN,3)
          RSIGS=RSIGS-SIGH(ICHN)
          IF(RSIGS.LE.0D0) GOTO 230
  220   CONTINUE
 
C...Reassign to appropriate process codes.
  230   ISUBMI=ICONMI/10
        ICONMI=MOD(ICONMI,10)
 
C...Choose new quark flavour for annihilation graphs
        IF(ISUBMI.EQ.12.OR.ISUBMI.EQ.53) THEN
          SH=VINT(21)*VINT(2)
          CALL PYWIDT(21,SH,WDTP,WDTE)
  240     RKFL=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))*PYR(0)
          DO 250 I=1,MDCY(21,3)
            KFLF=KFDP(I+MDCY(21,2)-1,1)
            RKFL=RKFL-(WDTE(I,1)+WDTE(I,2)+WDTE(I,4))
            IF(RKFL.LE.0D0) GOTO 260
  250     CONTINUE
  260     IF(ISUBMI.EQ.53.AND.ICONMI.LE.2) THEN
            IF(KFLF.GE.4) GOTO 240
          ELSEIF(ISUBMI.EQ.53.AND.ICONMI.LE.4) THEN
            KFLF=4
            ICONMI=ICONMI-2
          ELSEIF(ISUBMI.EQ.53) THEN
            KFLF=5
            ICONMI=ICONMI-4
          ENDIF
        ENDIF
 
C...Final state flavours and colour flow: default values
        JS=1
        KFL3=KFL1
        KFL4=KFL2
        KCC=20
        KCS=ISIGN(1,KFL1)
 
        IF(ISUBMI.EQ.11) THEN
C...f + f' -> f + f' (g exchange); th = (p(f)-p(f))**2
          KCC=ICONMI
          IF(KFL1*KFL2.LT.0) KCC=KCC+2
 
        ELSEIF(ISUBMI.EQ.12) THEN
C...f + fbar -> f' + fbar'; th = (p(f)-p(f'))**2
          KFL3=ISIGN(KFLF,KFL1)
          KFL4=-KFL3
          KCC=4
 
        ELSEIF(ISUBMI.EQ.13) THEN
C...f + fbar -> g + g; th arbitrary
          KFL3=21
          KFL4=21
          KCC=ICONMI+4
 
        ELSEIF(ISUBMI.EQ.28) THEN
C...f + g -> f + g; th = (p(f)-p(f))**2
          IF(KFL1.EQ.21) JS=2
          KCC=ICONMI+6
          IF(KFL1.EQ.21) KCC=KCC+2
          IF(KFL1.NE.21) KCS=ISIGN(1,KFL1)
          IF(KFL2.NE.21) KCS=ISIGN(1,KFL2)
 
        ELSEIF(ISUBMI.EQ.53) THEN
C...g + g -> f + fbar; th arbitrary
          KCS=(-1)**INT(1.5D0+PYR(0))
          KFL3=ISIGN(KFLF,KCS)
          KFL4=-KFL3
          KCC=ICONMI+10
 
        ELSEIF(ISUBMI.EQ.68) THEN
C...g + g -> g + g; th arbitrary
          KCC=ICONMI+12
          KCS=(-1)**INT(1.5D0+PYR(0))
        ENDIF
 
C...Check that massive sea quarks have non-zero phase space for g -> Q Q
        IF (IABS(KFL3).EQ.4.OR.IABS(KFL4).EQ.4.OR.IABS(KFL3).EQ.5
     &       .OR.IABS(KFL4).EQ.5) THEN
          RMMAX2=MAX(PMAS(PYCOMP(KFL3),1),PMAS(PYCOMP(KFL4),1))**2
          IF (PT2.LE.1.05*RMMAX2) THEN
            IF (NTRY.EQ.2) CALL PYERRM(9,'(PYPTMI:) Heavy quarks'
     &           //' too close to threshold (2nd try).')
            GOTO 210
          ENDIF
        ENDIF
 
C...Store flavours of scattering.
        MINT(13)=KFL1
        MINT(14)=KFL2
        MINT(15)=KFL1
        MINT(16)=KFL2
        MINT(21)=KFL3
        MINT(22)=KFL4
 
C...Set flavours and mothers of scattering partons.
        K(N+1,1)=14
        K(N+2,1)=14
        K(N+3,1)=3
        K(N+4,1)=3
        K(N+1,2)=KFL1
        K(N+2,2)=KFL2
        K(N+3,2)=KFL3
        K(N+4,2)=KFL4
        K(N+1,3)=MINT(83)+1
        K(N+2,3)=MINT(83)+2
        K(N+3,3)=N+1
        K(N+4,3)=N+2
 
C...Store colour connection indices.
        DO 270 J=1,2
          JC=J
          IF(KCS.EQ.-1) JC=3-J
          IF(ICOL(KCC,1,JC).NE.0) K(N+1,J+3)=N+ICOL(KCC,1,JC)
          IF(ICOL(KCC,2,JC).NE.0) K(N+2,J+3)=N+ICOL(KCC,2,JC)
          IF(ICOL(KCC,3,JC).NE.0) K(N+3,J+3)=MSTU(5)*(N+ICOL(KCC,3,JC))
          IF(ICOL(KCC,4,JC).NE.0) K(N+4,J+3)=MSTU(5)*(N+ICOL(KCC,4,JC))
  270   CONTINUE
 
C...Store incoming and outgoing partons in their CM-frame.
        SHR=SQRT(VINT(21))*VINT(1)
        P(N+1,3)=0.5D0*SHR
        P(N+1,4)=0.5D0*SHR
        P(N+2,3)=-0.5D0*SHR
        P(N+2,4)=0.5D0*SHR
        P(N+3,5)=PYMASS(K(N+3,2))
        P(N+4,5)=PYMASS(K(N+4,2))
        IF(P(N+3,5)+P(N+4,5).GE.SHR) THEN
          IFAIL=1
          RETURN
        ENDIF
        P(N+3,4)=0.5D0*(SHR+(P(N+3,5)**2-P(N+4,5)**2)/SHR)
        P(N+3,3)=SQRT(MAX(0D0,P(N+3,4)**2-P(N+3,5)**2))
        P(N+4,4)=SHR-P(N+3,4)
        P(N+4,3)=-P(N+3,3)
 
C...Rotate outgoing partons using cos(theta)=(th-uh)/lam(sh,sqm3,sqm4)
        PHI=PARU(2)*PYR(0)
        CALL PYROBO(N+3,N+4,ACOS(VINT(23)),PHI,0D0,0D0,0D0)
 
C...Global statistics.
        MINT(351)=MINT(351)+1
        VINT(351)=VINT(351)+SQRT(P(N+3,1)**2+P(N+3,2)**2)
        IF (MINT(351).EQ.1) VINT(356)=SQRT(P(N+3,1)**2+P(N+3,2)**2)
 
C...Keep track of loose colour ends and information on scattering.
        MINT(31)=MINT(31)+1
        MINT(36)=MINT(31)
        PT2MI(MINT(36))=PT2
        IMISEP(MINT(31))=N+4
        DO 280 JS=1,2
          IMI(JS,MINT(31),1)=N+JS
          IMI(JS,MINT(31),2)=0
          XMI(JS,MINT(31))=VINT(40+JS)
          NMI(JS)=NMI(JS)+1
C...Update cumulative counters
          VINT(142+JS)=VINT(142+JS)-VINT(40+JS)
          VINT(150+JS)=VINT(150+JS)+VINT(40+JS)
  280   CONTINUE
 
C...Add to list of final state partons
        IPART(NPART+1)=N+3
        IPART(NPART+2)=N+4
        PTPART(NPART+1)=SQRT(PT2)
        PTPART(NPART+2)=SQRT(PT2)
        NPART=NPART+2
 
C...Initialize ISR
        NISGEN(1,MINT(31))=0
        NISGEN(2,MINT(31))=0
 
C...Update ER
        N=N+4
        IF(N.GT.MSTU(4)-MSTU(32)-10) THEN
          CALL PYERRM(11,'(PYMIGN:) no more memory left in PYJETS')
          MINT(51)=1
          RETURN
        ENDIF
 
C...Finally, assign colour tags to new partons
        DO 300 JS=1,2
          I1=IMI(JS,MINT(31),1)
          I2=IMI(3-JS,MINT(31),1)
          DO 290 JCS=4,5
            IF (K(I1,2).NE.21.AND.(9-2*JCS).NE.ISIGN(1,K(I1,2)))
     &           GOTO 290
            IF (K(I1,JCS)/MSTU(5)**2.NE.0) GOTO 290
            KCS=JCS
            CALL PYCTTR(I1,KCS,I2)
            IF(MINT(51).NE.0) RETURN
  290     CONTINUE
  300   CONTINUE
 
C----------------------------------------------------------------------
C...MODE=2: Decide whether quarks in last scattering were valence,
C...companion, or sea.
      ELSEIF (MODE.EQ.2) THEN
        JS=MINT(30)
        MI=MINT(36)
        PT2=PT2NOW
        KFSBM=ISIGN(1,MINT(10+JS))
        IFL=K(IMI(JS,MI,1),2)
        IMI(JS,MI,2)=0
        IF (IABS(IFL).GE.6) THEN
          IF (IABS(IFL).EQ.6) THEN
            CALL PYERRM(29,'(PYPTMI:) top in initial state!')
          ENDIF
          RETURN
        ENDIF
C...Get PDFs at X(rescaled) and PT2 of the current initiator.
C...(Do not include the parton itself in the X rescaling.)
        X=XMI(JS,MI)
        XRSC=X/(VINT(142+JS)+X)
C...Note: XPSVC = x*pdf.
        MINT(30)=JS
        CALL PYPDFU(KFBEAM(JS),XRSC,PT2,XPQ)
        SEA=XPSVC(IFL,-1)
        VAL=XPSVC(IFL,0) 
C...Ensure that pdfs are positive definite   
        IF (SEA.LT.0D0) THEN
          CALL PYERRM(9,'(PYPTMI:) Sea distribution negative.')
          SEA=MAX(0D0,SEA)
        ELSEIF (VAL.LT.0D0) THEN
          CALL PYERRM(9,'(PYPTMI:) Val distribution negative.')
          VAL=MAX(0D0,VAL)          
        ENDIF
        CMP=0D0
        DO 310 IVC=1,NVC(JS,IFL)
          CMP=CMP+XPSVC(IFL,IVC)
  310   CONTINUE
C...PS 05 Aug 2012: bug fix to prevent heavy companion quarks from being
C...picked up by MPI (necessary since intertwining not implemented)
C...Here simply reclassify companions as ordinary SEA. Will give 
C...additional spurious companions, but is simplest solution.
        IF (IABS(IFL).EQ.4.OR.IABS(IFL).EQ.5) THEN
          SEA = SEA + CMP
          CMP = 0D0
        ENDIF
 
        NTRY=0
C...Decide (Extra factor x cancels in the dvision).
  320   RVCS=PYR(0)*(SEA+VAL+CMP)
        IVNOW=1
        NTRY=NTRY+1
  330   IF (RVCS.LE.VAL.AND.IVNOW.GE.1) THEN
C...Safety check that valence present; pi0/gamma/K0S/K0L special cases.
          IVNOW=0
          IF(KFIVAL(JS,1).EQ.IFL) IVNOW=IVNOW+1
          IF(KFIVAL(JS,2).EQ.IFL) IVNOW=IVNOW+1
          IF(KFIVAL(JS,3).EQ.IFL) IVNOW=IVNOW+1
          IF(KFIVAL(JS,1).EQ.0) THEN
            IF(KFBEAM(JS).EQ.111.AND.IABS(IFL).LE.2) IVNOW=1
            IF(KFBEAM(JS).EQ.22.AND.IABS(IFL).LE.5) IVNOW=1
            IF((KFBEAM(JS).EQ.130.OR.KFBEAM(JS).EQ.310).AND.
     &           (IABS(IFL).EQ.1.OR.IABS(IFL).EQ.3)) IVNOW=1
          ELSE
C...Count down valence remaining. Do not count current scattering.
            DO 340 I1=1,NMI(JS)
              IF (I1.EQ.MINT(36)) GOTO 340
              IF (K(IMI(JS,I1,1),2).EQ.IFL.AND.IMI(JS,I1,2).EQ.0)
     &             IVNOW=IVNOW-1
  340       CONTINUE
          ENDIF
          IF(IVNOW.EQ.0) GOTO 330
C...Mark valence.
          IMI(JS,MI,2)=0
C...Sets valence content of gamma, pi0, K0S, K0L if not done.
          IF(KFIVAL(JS,1).EQ.0) THEN
            IF(KFBEAM(JS).EQ.111.OR.KFBEAM(JS).EQ.22) THEN
              KFIVAL(JS,1)=IFL
              KFIVAL(JS,2)=-IFL
            ELSEIF(KFBEAM(JS).EQ.130.OR.KFBEAM(JS).EQ.310) THEN
              KFIVAL(JS,1)=IFL
              IF(IABS(IFL).EQ.1) KFIVAL(JS,2)=ISIGN(3,-IFL)
              IF(IABS(IFL).NE.1) KFIVAL(JS,2)=ISIGN(1,-IFL)
            ENDIF
          ENDIF
 
        ELSEIF (RVCS.LE.VAL+SEA) THEN
C...If sea, add opposite sign companion parton. Store X and I.
          NVC(JS,-IFL)=NVC(JS,-IFL)+1
          XASSOC(JS,-IFL,NVC(JS,-IFL))=XMI(JS,MI)
C...Set pointer to companion
          IMI(JS,MI,2)=-NVC(JS,-IFL)
 
        ELSE
C...If companion, check whether we've got any in the books
          IF (NVC(JS,IFL).EQ.0) THEN
            CMP=0D0
C...Only report error first time for this event
            IF (NTRY.EQ.1) 
     &           CALL PYERRM(9,'(PYPTMI:) No cmp quark, but pdf != 0!')
C...Try a few times
            IF (NTRY.LE.10) THEN
              GOTO 320
C... But if it stil fails, abort this event
            ELSE
              MINT(51)=1
              RETURN
            ENDIF
          ENDIF
C...If several possibilities, decide which one
          CMPSUM=VAL+SEA
          ISEL=0
  350     ISEL=ISEL+1
          CMPSUM=CMPSUM+XPSVC(IFL,ISEL)
          IF (RVCS.GT.CMPSUM.AND.ISEL.LT.NVC(JS,IFL)) GOTO 350
C...Find original sea (anti-)quark. Do not consider current scattering.
          IASSOC=0
          DO 360 I1=1,NMI(JS)
            IF (I1.EQ.MINT(36)) GOTO 360
            IF (K(IMI(JS,I1,1),2).NE.-IFL) GOTO 360
            IF (-IMI(JS,I1,2).EQ.ISEL) THEN
              IMI(JS,MI,2)=IMI(JS,I1,1)
              IMI(JS,I1,2)=IMI(JS,MI,1)
            ENDIF
  360     CONTINUE
C...Mark companion "out-kicked".
          XASSOC(JS,IFL,ISEL)=-XASSOC(JS,IFL,ISEL)
        ENDIF
 
      ENDIF
      RETURN
      END
