
C*********************************************************************
 
C...PYPTIS
C...Generates pT-ordered spacelike initial-state parton showers and
C...trial joinings.
C...MODE=-1: Initialize ISR from scratch, starting from the hardest
C...         interaction initiators at PT2NOW.
C...MODE= 0: Generate a trial branching on interaction MINT(36), side
C...         MINT(30). Start evolution at PT2NOW, solve Sudakov for PT2.
C...         Store in /PYISMX/ if PT2 is largest so far. Abort if PT2
C...         is below PT2CUT.
C...         (Also generate test joinings if MSTP(96)=1.)
C...MODE= 1: Accept stored shower branching. Update event record etc.
C...PT2NOW : Starting (max) PT2 scale for evolution.
C...PT2CUT : Lower limit for evolution.
C...PT2    : Result of evolution. Generated PT2 for trial emission.
C...IFAIL  : Status return code. IFAIL=0 when all is well.
 
      SUBROUTINE PYPTIS(MODE,PT2NOW,PT2CUT,PT2,IFAIL)
 
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
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINTM/KFIVAL(2,3),NMI(2),IMI(2,800,2),NVC(2,-6:6),
     &     XASSOC(2,-6:6,240),XPSVC(-6:6,-1:240),PVCTOT(2,-1:1),
     &     XMI(2,240),PT2MI(240),IMISEP(0:240)
      COMMON/PYISMX/MIMX,JSMX,KFLAMX,KFLCMX,KFBEAM(2),NISGEN(2,240),
     &     PT2MX,PT2AMX,ZMX,RM2CMX,Q2BMX,PHIMX
      COMMON/PYCTAG/NCT,MCT(4000,2)
      COMMON/PYISJN/MJN1MX,MJN2MX,MJOIND(2,240)
      SAVE /PYPART/,/PYJETS/,/PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/,
     &     /PYINT2/,/PYINTM/,/PYISMX/,/PYCTAG/,/PYISJN/
C...Local variables
      DIMENSION ZSAV(2,240),PT2SAV(2,240),
     &     XFB(-25:25),XFA(-25:25),XFN(-25:25),XFJ(-25:25),
     &     WTAP(-25:25),WTPDF(-25:25),SHTNOW(240),
     &     WTAPJ(240),WTPDFJ(240),X1(240),Y(240)
      SAVE ZSAV,PT2SAV,XFB,XFA,XFN,WTAP,WTPDF,XMXC,SHTNOW,
     &     RMB2,RMC2,ALAM3,ALAM4,ALAM5,TMIN,PTEMAX,WTEMAX,AEM2PI
C...For check on excessive weights.
      CHARACTER CHWT*12
 
C...Only give errors for very large weights, otherwise just warnings
      DATA WTEMAX /1.5D0/
C...Only give errors for large pT, otherwise just warnings
      DATA PTEMAX /5D0/
 
      IFAIL=-1
 
C----------------------------------------------------------------------
C...MODE=-1: Initialize initial state showers from scratch, i.e.
C...starting from the hardest interaction initiators.
      IF (MODE.EQ.-1) THEN
C...Set hard scattering SHAT.
        SHTNOW(1)=VINT(44)
C...Mass thresholds and Lambda for QCD evolution.
        AEM2PI=PARU(101)/PARU(2)
        RMB=PMAS(5,1)
        RMC=PMAS(4,1)
        ALAM4=PARP(61)
        IF(MSTU(112).LT.4) ALAM4=PARP(61)*(PARP(61)/RMC)**(2D0/25D0)
        IF(MSTU(112).GT.4) ALAM4=PARP(61)*(RMB/PARP(61))**(2D0/25D0)
        ALAM5=ALAM4*(ALAM4/RMB)**(2D0/23D0)
        ALAM3=ALAM4*(RMC/ALAM4)**(2D0/27D0)
C...Optionally use Lambda_MC = Lambda_CMW 
        IF (MSTP(64).EQ.3) THEN
          ALAM5 = ALAM5 * 1.569 
          ALAM4 = ALAM4 * 1.618 
          ALAM3 = ALAM3 * 1.661 
        ENDIF
        RMB2=RMB**2
        RMC2=RMC**2
C...Massive quark forced creation threshold (in M**2).
        TMIN=1.01D0
C...Set upper limit for X (ensures some X left for beam remnant).
        XMXC=1D0-2D0*PARP(111)/VINT(1)
 
        IF (MSTP(61).GE.1) THEN
C...Initial values: flavours, momenta, virtualities.
          DO 100 JS=1,2
            NISGEN(JS,1)=0
 
C...Special kinematics check for c/b quarks (that g -> c cbar or
C...b bbar kinematically possible).
            KFLB=K(IMI(JS,1,1),2)
            KFLCB=IABS(KFLB)
            IF(KFBEAM(JS).NE.22.AND.(KFLCB.EQ.4.OR.KFLCB.EQ.5)) THEN
C...Check PT2MAX > mQ^2
              IF (VINT(56).LT.1.05D0*PMAS(PYCOMP(KFLCB),1)**2) THEN
                CALL PYERRM(9,'(PYPTIS:) PT2MAX < 1.05 * MQ**2. '//
     &               'No Q creation possible.')
                MINT(51)=1
                RETURN
              ELSE
C...Check for physical z values (m == MQ / sqrt(s))
C...For creation diagram, x < z < (1-m)/(1+m(1-m))
                FMQ=PMAS(KFLCB,1)/SQRT(SHTNOW(1))
                ZMXCR=(1D0-FMQ)/(1D0+FMQ*(1D0-FMQ))
                IF (XMI(JS,1).GT.0.9D0*ZMXCR) THEN
                  CALL PYERRM(9,'(PYPTIS:) No physical z value for '//
     &                 'Q creation.')
                  MINT(51)=1
                  RETURN
                ENDIF
              ENDIF
            ENDIF
  100     CONTINUE
        ENDIF
 
        MINT(354)=0
C...Zero joining array
        DO 110 MJ=1,240
          MJOIND(1,MJ)=0
          MJOIND(2,MJ)=0
  110   CONTINUE
 
C----------------------------------------------------------------------
C...MODE= 0: Generate a trial branching on interaction MINT(36) side
C...MINT(30). Store if emission PT2 scale is largest so far.
C...Also generate test joinings if MSTP(96)=1.
      ELSEIF(MODE.EQ.0) THEN
        IFAIL=-1
        MECOR=0
        ISUB=MINT(1)
        JS=MINT(30)
C...No shower for structureless beam
        IF (MINT(44+JS).EQ.1) RETURN
        MI=MINT(36)
        SHAT=VINT(44)
C...Absolute shower max scale = VINT(56)
        IF (MSTP(67).NE.0) THEN
          PT2 = MIN(PT2NOW,VINT(56))
        ELSE
C...For MSTP(67)=0, adjust starting scale by PARP(67)
          PT2=MIN(PT2NOW,PARP(67)*VINT(56))
        ENDIF
        IF (NISGEN(1,MI).EQ.0.AND.NISGEN(2,MI).EQ.0) SHTNOW(MI)=SHAT
C...Define for which processes ME corrections have been implemented.
        IF(MSTP(68).EQ.1.OR.MSTP(68).EQ.3) THEN
          IF(ISUB.EQ.1.OR.ISUB.EQ.2.OR.ISUB.EQ.141.OR.ISUB.EQ
     &         .142.OR.ISUB.EQ.144) MECOR=1
          IF(ISUB.EQ.102.OR.ISUB.EQ.152.OR.ISUB.EQ.157) MECOR=2
          IF(ISUB.EQ.3.OR.ISUB.EQ.151.OR.ISUB.EQ.156) MECOR=3
C...Calculate preweighting factor for ME-corrected processes.
          IF(MECOR.GE.1) CALL PYMEMX(MECOR,WTFF,WTGF,WTFG,WTGG)
        ENDIF
C...Basic info on daughter for which to find mother.
        KFLB=K(IMI(JS,MI,1),2)
        KFLBA=IABS(KFLB)
C...KSVCB: -1 for sea or first companion, 0 for valence or gluon, >1 for
C...second companion.
        KSVCB=MAX(-1,IMI(JS,MI,2))
C...Treat "first" companion of a pair like an ordinary sea quark
C...(except that creation diagram is not allowed)
        IF(IMI(JS,MI,2).GT.IMISEP(MI)) KSVCB=-1
C...X (rescaled to [0,1])
        XB=XMI(JS,MI)/VINT(142+JS)
C...Massive quarks (use physical masses.)
        RMQ2=0D0
        MQMASS=0
        IF (KFLBA.EQ.4.OR.KFLBA.EQ.5) THEN
          RMQ2=RMC2
          IF (KFLBA.EQ.5) RMQ2=RMB2
C...Special threshold treatment for non-photon beams
          IF (KFBEAM(JS).NE.22) MQMASS=KFLBA
C...Check that not below mass threshold.
          IF(MQMASS.GT.0.AND.PT2.LT.TMIN*RMQ2) THEN
            CALL PYERRM(9,'(PYPTIS:) PT2 < 1.01 * MQ**2. '//
     &        'No Q creation possible.')
            MINT(51)=1
C...Special return code if failing before any evolution at all: bad event
            IF (NISGEN(1,MI).EQ.0.AND.NISGEN(2,MI).EQ.0) MINT(51)=2
            RETURN
          ENDIF

        ENDIF
 
C...Flags for parton distribution calls.
        MINT(105)=MINT(102+JS)
        MINT(109)=MINT(106+JS)
        VINT(120)=VINT(2+JS)
 
C...Calculate initial parton distribution weights.
        IF(XB.GE.XMXC) THEN
          RETURN
        ELSEIF(MQMASS.EQ.0) THEN
          CALL PYPDFU(KFBEAM(JS),XB,PT2,XFB)
        ELSE
C...Initialize massive quark PT2 dependent pdf underestimate.
          PT20=PT2
          CALL PYPDFU(KFBEAM(JS),XB,PT20,XFB)
C.!.Tentative treatment of massive valence quarks.
          XQ0=MAX(1D-10,XPSVC(KFLB,KSVCB))
          XG0=XFB(21)
          TPM0=LOG(PT20/RMQ2)
          WPDF0=TPM0*XG0/XQ0
        ENDIF
        IF (KFLBA.LE.6) THEN
C...For quarks, only include respective sea, val, or cmp part.
          IF (KSVCB.LE.0) THEN
            XFB(KFLB)=XPSVC(KFLB,KSVCB)
          ELSE
C...Find companion's companion
            MISEA=0
  120       MISEA=MISEA+1
            IF (IMI(JS,MISEA,2).NE.IMI(JS,MI,1)) GOTO 120
            XS=XMI(JS,MISEA)
            XREM=VINT(142+JS)
            YS=XS/(XREM+XS)
C...Momentum fraction of the companion quark.
C...Rescale from XB = x/XREM to YB = x/(1-Sum_rest) -> factor (1-YS).
            YB=XB*(1D0-YS)
            XFB(KFLB)=PYFCMP(YB/VINT(140),YS/VINT(140),MSTP(87))
          ENDIF
        ENDIF
 
C...Determine overestimated z range: switch at c and b masses.
  130   IF (PT2.GT.TMIN*RMB2) THEN
          IZRG=3
          PT2MNE=MAX(TMIN*RMB2,PT2CUT)
          B0=23D0/6D0
          ALAM2=ALAM5**2
        ELSEIF(PT2.GT.TMIN*RMC2) THEN
          IZRG=2
          PT2MNE=MAX(TMIN*RMC2,PT2CUT)
          B0=25D0/6D0
          ALAM2=ALAM4**2
        ELSE
          IZRG=1
          PT2MNE=PT2CUT
          B0=27D0/6D0
          ALAM2=ALAM3**2
        ENDIF
C...Divide Lambda by PARP(64) (equivalent to mult pT2 by PARP(64))
        ALAM2=ALAM2/PARP(64)
C...Overestimated ZMAX:
        IF (MQMASS.EQ.0) THEN
C...Massless
          ZMAX=1D0-0.5D0*(PT2MNE/SHTNOW(MI))*(SQRT(1D0+4D0*SHTNOW(MI)
     &         /PT2MNE)-1D0)
        ELSE
C...Massive (limit for bremsstrahlung diagram > creation)
          FMQ=SQRT(RMQ2/SHTNOW(MI))
          ZMAX=1D0/(1D0+FMQ)
        ENDIF
        ZMIN=XB/XMXC
 
C...If kinematically impossible then do not evolve.
        IF(PT2.LT.PT2CUT.OR.ZMAX.LE.ZMIN) RETURN
 
C...Reset Altarelli-Parisi and PDF weights.
        DO 140 KFL=-5,5
          WTAP(KFL)=0D0
          WTPDF(KFL)=0D0
  140   CONTINUE
        WTAP(21)=0D0
        WTPDF(21)=0D0
C...Zero joining weights and compute X(partner) and X(mother) values.
        NJN=0
        IF (MSTP(96).NE.0) THEN
          DO 150 MJ=1,MINT(31)
            WTAPJ(MJ)=0D0
            WTPDFJ(MJ)=0D0
            X1(MJ)=XMI(JS,MJ)/(VINT(142+JS)+XMI(JS,MJ))
            Y(MJ)=(XMI(JS,MI)+XMI(JS,MJ))/(VINT(142+JS)+XMI(JS,MJ)
     &           +XMI(JS,MI))
  150     CONTINUE
        ENDIF
 
C...Approximate Altarelli-Parisi weights (integrated AP dz).
C...q -> q, g -> q or q -> q + gamma (already set which).
        IF(KFLBA.LE.5) THEN
C...Val and cmp quarks get an extra sqrt(z) to smooth their bumps.
          IF (KSVCB.LT.0) THEN
            WTAP(KFLB)=(8D0/3D0)*LOG((1D0-ZMIN)/(1D0-ZMAX))
          ELSE
            RMIN=(1+SQRT(ZMIN))/(1-SQRT(ZMIN))
            RMAX=(1+SQRT(ZMAX))/(1-SQRT(ZMAX))
            WTAP(KFLB)=(8D0/3D0)*LOG(RMAX/RMIN)
          ENDIF
          WTAP(21)=0.5D0*(ZMAX-ZMIN)
          WTAPE=(2D0/9D0)*LOG((1D0-ZMIN)/(1D0-ZMAX))
          IF(MOD(KFLBA,2).EQ.0) WTAPE=4D0*WTAPE
          IF(MECOR.GE.1.AND.NISGEN(JS,MI).EQ.0) THEN
            WTAP(KFLB)=WTFF*WTAP(KFLB)
            WTAP(21)=WTGF*WTAP(21)
            WTAPE=WTFF*WTAPE
          ENDIF
          IF(MSTP(61).EQ.1) WTAPE=0D0
          IF (KSVCB.GE.1) THEN
C...Kill normal creation but add joining diagrams for cmp quark.
            WTAP(21)=0D0
            IF (KFLBA.EQ.4.OR.KFLBA.EQ.5) THEN
              CALL PYERRM(9,'(PYPTIS:) Sorry, I got a heavy companion'//
     &             " quark here. Not handled yet, giving up!")
              PT2=0D0
              MINT(51)=1
              RETURN
            ENDIF
C...Check for possible joinings
            IF (MSTP(96).NE.0.AND.MJOIND(JS,MI).EQ.0) THEN
C...Find companion's companion.
              MJ=0
  160         MJ=MJ+1
              IF (IMI(JS,MJ,2).NE.IMI(JS,MI,1)) GOTO 160
              IF (MJOIND(JS,MJ).EQ.0) THEN
                Y(MI)=YB+YS
                Z=YB/Y(MI)
                WTAPJ(MJ)=Z*(1D0-Z)*0.5D0*(Z**2+(1D0-Z)**2)
                IF (WTAPJ(MJ).GT.1D-6) THEN
                  NJN=1
                ELSE
                  WTAPJ(MJ)=0D0
                ENDIF
              ENDIF
C...Add trial gluon joinings.
              DO 170 MJ=1,MINT(31)
                KFLC=K(IMI(JS,MJ,1),2)
                IF (KFLC.NE.21.OR.MJOIND(JS,MJ).NE.0) GOTO 170
                Z=XMI(JS,MJ)/(XMI(JS,MI)+XMI(JS,MJ))
                WTAPJ(MJ)=6D0*(Z**2+(1D0-Z)**2)
                IF (WTAPJ(MJ).GT.1D-6) THEN
                  NJN=NJN+1
                ELSE
                  WTAPJ(MJ)=0D0
                ENDIF
  170         CONTINUE
            ENDIF
          ELSEIF (IMI(JS,MI,2).GE.0) THEN
C...Kill creation diagram for val quarks and sea quarks with companions.
            WTAP(21)=0D0
          ELSEIF (MQMASS.EQ.0) THEN
C...Extra safety factor for massless sea quark creation.
            WTAP(21)=WTAP(21)*1.25D0
          ENDIF
 
C...  q -> g, g -> g.
        ELSEIF(KFLB.EQ.21) THEN
C...Here we decide later whether a quark picked up is valence or
C...sea, so we maintain the extra factor sqrt(z) since we deal
C...with the *sum* of sea and valence in this context.
          WTAPQ=(16D0/3D0)*(SQRT(1D0/ZMIN)-SQRT(1D0/ZMAX))
C...new: do not allow backwards evol to pick up heavy flavour.
          DO 180 KFL=1,MIN(3,MSTP(58))
            WTAP(KFL)=WTAPQ
            WTAP(-KFL)=WTAPQ
  180     CONTINUE
          WTAP(21)=6D0*LOG(ZMAX*(1D0-ZMIN)/(ZMIN*(1D0-ZMAX)))
          IF(MECOR.GE.1.AND.NISGEN(JS,MI).EQ.0) THEN
            WTAPQ=WTFG*WTAPQ
            WTAP(21)=WTGG*WTAP(21)
          ENDIF
C...Check for possible joinings (companions handled separately above)
          IF (MSTP(96).NE.0.AND.MINT(31).GE.2.AND.MJOIND(JS,MI).EQ.0)
     &         THEN
            DO 190 MJ=1,MINT(31)
              IF (MJ.EQ.MI.OR.MJOIND(JS,MJ).NE.0) GOTO 190
              KSVCC=IMI(JS,MJ,2)
              IF (IMI(JS,MJ,2).GT.IMISEP(MJ)) KSVCC=-1
              IF (KSVCC.GE.1) GOTO 190
              KFLC=K(IMI(JS,MJ,1),2)
C...Only try g -> g + g once.
              IF (MJ.GT.MI.AND.KFLC.EQ.21) GOTO 190
              Z=XMI(JS,MJ)/(XMI(JS,MI)+XMI(JS,MJ))
              IF (KFLC.EQ.21) THEN
                WTAPJ(MJ)=6D0*(Z**2+(1D0-Z)**2)
              ELSE
                WTAPJ(MJ)=Z*4D0/3D0*(1D0+Z**2)
              ENDIF
              IF (WTAPJ(MJ).GT.1D-6) THEN
                NJN=NJN+1
              ELSE
                WTAPJ(MJ)=0D0
              ENDIF
  190       CONTINUE
          ENDIF
        ENDIF
 
C...Initialize massive quark evolution
        IF (MQMASS.NE.0) THEN
          RML=(RMQ2+VINT(18))/ALAM2
          TML=LOG(RML)
          TPL=LOG((PT2+VINT(18))/ALAM2)
          TPM=LOG((PT2+VINT(18))/RMQ2)
          WN=WTAP(21)*WPDF0/B0
        ENDIF
 
 
C...Loopback point for iteration
        NTRY=0
        NTHRES=0
  200   NTRY=NTRY+1
        IF(NTRY.GT.500) THEN
          CALL PYERRM(9,'(PYPTIS:) failed to evolve shower.')
          MINT(51)=1
          RETURN
        ENDIF
 
C...  Calculate PDF weights and sum for evolution rate.
        WTSUM=0D0
        XFBO=MAX(1D-10,XFB(KFLB))
        DO 210 KFL=-5,5
          WTPDF(KFL)=XFB(KFL)/XFBO
          WTSUM=WTSUM+WTAP(KFL)*WTPDF(KFL)
  210   CONTINUE
C...Only add gluon mother diagram for massless KFLB.
        IF(MQMASS.EQ.0) THEN
          WTPDF(21)=XFB(21)/XFBO
          WTSUM=WTSUM+WTAP(21)*WTPDF(21)
        ENDIF
        WTSUM=MAX(0.0001D0,WTSUM)
        WTSUMS=WTSUM
C...Add joining diagrams where applicable.
        WTJOIN=0D0
        IF (MSTP(96).NE.0.AND.NJN.NE.0) THEN
          DO 220 MJ=1,MINT(31)
            IF (WTAPJ(MJ).LT.1D-3) GOTO 220
            WTPDFJ(MJ)=1D0/XFBO
C...x and x*pdf (+ sea/val) for parton C.
            KFLC=K(IMI(JS,MJ,1),2)
            KFLCA=IABS(KFLC)
            KSVCC=MAX(-1,IMI(JS,MJ,2))
            IF (IMI(JS,MJ,2).GT.IMISEP(MJ)) KSVCC=-1
            MINT(30)=JS
            MINT(36)=MJ
            CALL PYPDFU(KFBEAM(JS),X1(MJ),PT2,XFJ)
            MINT(36)=MI
            IF (KFLCA.LE.6.AND.KSVCC.LE.0) THEN
              XFJ(KFLC)=XPSVC(KFLC,KSVCC)
            ELSEIF (KSVCC.GE.1) THEN
              print*, 'error! parton C is companion!'
            ENDIF
            WTPDFJ(MJ)=WTPDFJ(MJ)/XFJ(KFLC)
C...x and x*pdf (+ sea/val) for parton A.
            KFLA=21
            KSVCA=0
            IF (KFLCA.EQ.21.AND.KFLBA.LE.5) THEN
              KFLA=KFLB
              KSVCA=KSVCB
            ELSEIF (KFLBA.EQ.21.AND.KFLCA.LE.5) THEN
              KFLA=KFLC
              KSVCA=KSVCC
            ENDIF
            MINT(30)=JS
            IF (KSVCA.LE.0) THEN
C...Consider C the "evolved" parton if B is gluon. Val/sea
C...counting will then be done correctly in PYPDFU.
              IF (KFLBA.EQ.21) MINT(36)=MJ
              CALL PYPDFU(KFBEAM(JS),Y(MJ),PT2,XFJ)
              MINT(36)=MI
              IF (IABS(KFLA).LE.6) XFJ(KFLA)=XPSVC(KFLA,KSVCA)
            ELSE
C...If parton A is companion, use Y(MI) and YS in call to PYFCMP.
              XFJ(KFLA)=PYFCMP(Y(MI)/VINT(140),YS/VINT(140),MSTP(87))
            ENDIF
            WTPDFJ(MJ)=XFJ(KFLA)*WTPDFJ(MJ)
            WTJOIN=WTJOIN+WTAPJ(MJ)*WTPDFJ(MJ)
  220     CONTINUE
        ENDIF
 
C...Pick normal pT2 (in overestimated z range).
  230   PT2OLD=PT2
        WTSUM=WTSUMS
        PT2=ALAM2*((PT2+VINT(18))/ALAM2)**(PYR(0)**(B0/WTSUM))-VINT(18)
        KFLC=21
 
C...Evolve q -> q gamma separately, pick it if larger pT.
        IF(KFLBA.LE.5.AND.MSTP(61).GE.2) THEN
          PT2QED=(PT2OLD+VINT(18))*PYR(0)**(1D0/(AEM2PI*WTAPE))-VINT(18)
          IF(PT2QED.GT.PT2) THEN
            PT2=PT2QED
            KFLC=22
            KFLA=KFLB
          ENDIF
        ENDIF
 
C...  Evolve massive quark creation separately.
        MCRQQ=0
        IF (MQMASS.NE.0) THEN
          PT2CR=(RMQ2+VINT(18))*(RML**(TPM/(TPL*PYR(0)**(-TML/WN)-TPM)))
     &         -VINT(18)
C...If massive quark also on opposite side, ensure sufficient remaining 
C...phase space also for creation of that quark
          TMINQQ = TMIN
          KFLOPP = K(IMI(3-JS,MI,1),2)
          IF (ABS(KFLOPP).EQ.4.OR.ABS(KFLOPP).EQ.5) TMINQQ = 1.05
C...Ensure mininimum PT2CR and force creation near threshold.
          IF (PT2CR.LT.TMINQQ*RMQ2) THEN
            NTHRES=NTHRES+1
            IF (NTHRES.GT.50) THEN
              CALL PYERRM(9,'(PYPTIS:) no phase space left for '//
     &             'massive quark creation. Gave up trying.')
              MINT(51)=1
C...Special return code if failing before any evolution at all: bad event
              IF (NISGEN(1,MI).EQ.0.AND.NISGEN(2,MI).EQ.0) MINT(51)=2
              RETURN
            ENDIF
            PT2=0D0
            PT2CR=TMINQQ*RMQ2
C...Signal that massive quark creation is being forced
            MCRQQ=2
          ENDIF
C...  Select largest PT2 (brems or creation):
          IF (PT2CR.GT.PT2) THEN
            MCRQQ=MAX(MCRQQ,1)
            WTSUM=0D0
            PT2=PT2CR
            KFLA=21
          ELSE
            MCRQQ=0
            KFLA=KFLB
          ENDIF
C...  Compute logarithms for this PT2
          TPL=LOG((PT2+VINT(18))/ALAM2)
          TPM=LOG((PT2+VINT(18))/(RMQ2+VINT(18)))
          WTCRQQ=TPM/LOG(PT2/RMQ2)
        ENDIF
 
C...Evolve joining separately
        MJOIN=0
        IF (MSTP(96).NE.0.AND.NJN.NE.0) THEN
          PT2JN=ALAM2*((PT2OLD+VINT(18))/ALAM2)**(PYR(0)**(B0/WTJOIN))
     &         -VINT(18)
          IF (PT2JN.GE.PT2) THEN
            MJOIN=1
            PT2=PT2JN
          ENDIF
        ENDIF
 
C...Loopback if crossed c/b mass thresholds.
        IF(IZRG.EQ.3.AND.PT2.LT.RMB2) THEN
          PT2=RMB2
         GOTO 130
        ELSEIF(IZRG.EQ.2.AND.PT2.LT.RMC2) THEN
          PT2=RMC2
          GOTO 130
        ENDIF
 
C...Speed up shower. Skip if higher-PT acceptable branching
C...already found somewhere else.
C...Also finish if below lower cutoff.
        IF ((PT2-PT2MX).LT.-0.001.OR.PT2.LT.PT2CUT) RETURN

C...Select parton A flavour (massive Q handled above.)
        IF (MQMASS.EQ.0.AND.KFLC.NE.22.AND.MJOIN.EQ.0) THEN
          WTRAN=PYR(0)*WTSUM
          KFLA=-6
  240     KFLA=KFLA+1
          WTRAN=WTRAN-WTAP(KFLA)*WTPDF(KFLA)
          IF(KFLA.LE.5.AND.WTRAN.GT.0D0) GOTO 240
          IF(KFLA.EQ.6) KFLA=21
        ELSEIF (MJOIN.EQ.1) THEN
C...Tentative joining accept/reject.
          WTRAN=PYR(0)*WTJOIN
          MJ=0
  250     MJ=MJ+1
          WTRAN=WTRAN-WTAPJ(MJ)*WTPDFJ(MJ)
          IF(MJ.LE.MINT(31)-1.AND.WTRAN.GT.0D0) GOTO 250
          IF(MJOIND(JS,MJ).NE.0.OR.MJOIND(JS,MI).NE.0) THEN
            CALL PYERRM(9,'(PYPTIS:) Attempted double joining.'//
     &           ' Rejected.')
            GOTO 230
          ENDIF
C...x*pdf (+ sea/val) at new pT2 for parton B.
          IF (KSVCB.LE.0) THEN
            MINT(30)=JS
            CALL PYPDFU(KFBEAM(JS),XB,PT2,XFB)
            IF (KFLBA.LE.6) XFB(KFLB)=XPSVC(KFLB,KSVCB)
          ELSE
C...Companion distributions do not evolve.
            XFB(KFLB)=XFBO
          ENDIF
          WTVETO=1D0/WTPDFJ(MJ)/XFB(KFLB)
          KFLC=K(IMI(JS,MJ,1),2)
          KFLCA=IABS(KFLC)
          KSVCC=MAX(-1,IMI(JS,MJ,2))
          IF (KSVCB.GE.1) KSVCC=-1
C...x*pdf (+ sea/val) at new pT2 for parton C.
          MINT(30)=JS
          MINT(36)=MJ
          CALL PYPDFU(KFBEAM(JS),X1(MJ),PT2,XFJ)
          MINT(36)=MI
          IF (KFLCA.LE.6.AND.KSVCC.LE.0) XFJ(KFLC)=XPSVC(KFLC,KSVCC)
          WTVETO=WTVETO/XFJ(KFLC)
C...x and x*pdf (+ sea/val) at new pT2 for parton A.
          KFLA=21
          KSVCA=0
          IF (KFLCA.EQ.21.AND.KFLBA.LE.5) THEN
            KFLA=KFLB
            KSVCA=KSVCB
          ELSEIF (KFLBA.EQ.21.AND.KFLCA.LE.5) THEN
            KFLA=KFLC
            KSVCA=KSVCC
          ENDIF
          IF (KSVCA.LE.0) THEN
            MINT(30)=JS
            IF (KFLB.EQ.21) MINT(36)=MJ
            CALL PYPDFU(KFBEAM(JS),Y(MJ),PT2,XFJ)
            MINT(36)=MI
            IF (IABS(KFLA).LE.6) XFJ(KFLA)=XPSVC(KFLA,KSVCA)
          ELSE
            XFJ(KFLA)=PYFCMP(Y(MJ)/VINT(140),YS/VINT(140),MSTP(87))
          ENDIF
C...PS 05 Aug 2012: bug fix to prevent heavy companion quarks from being
C...picked up by ISR (necessary since intertwining not implemented)
C...Here simply kill backwards-evolution probability.
          IF (KFLB.EQ.21.AND.(IABS(KFLA).EQ.4.OR.IABS(KFLA).EQ.5)) THEN
            IF (KSVCA.GE.1) WTVETO = 0D0
          ENDIF
          WTVETO=WTVETO*XFJ(KFLA)
C...Monte Carlo veto to accept trial joining
          IF (WTVETO.LT.PYR(0)) GOTO 200
C...If accept, save PT2 of this joining.
          IF (PT2.GT.PT2MX) THEN
            PT2MX=PT2
            JSMX=2+JS
            MJN1MX=MJ
            MJN2MX=MI
            WTAPJ(MJ)=0D0
            NJN=0
          ENDIF
C...Exit and continue evolution.
          GOTO 390
        ENDIF
        KFLAA=IABS(KFLA)
 
C...Choose z value (still in overestimated range) and corrective weight.
C...Unphysical z will be rejected below when Q2 has is computed.
        WTZ=0D0
 
C...Note: ME and MQ>0 give corrections to overall weights, not shapes.
C...q -> q + g or q -> q + gamma (already set which).
        IF (KFLAA.LE.5.AND.KFLBA.LE.5) THEN
          IF (KSVCB.LT.0) THEN
            Z=1D0-(1D0-ZMIN)*((1D0-ZMAX)/(1D0-ZMIN))**PYR(0)
          ELSE
            ZFAC=RMIN*(RMAX/RMIN)**PYR(0)
            Z=((1-ZFAC)/(1+ZFAC))**2
          ENDIF
          WTZ=0.5D0*(1D0+Z**2)
C...Massive weight correction.
          IF (KFLBA.GE.4) WTZ=WTZ-Z*(1D0-Z)**2*RMQ2/PT2
C...Valence quark weight correction (extra sqrt)
          IF (KSVCB.GE.0) WTZ=WTZ*SQRT(Z)
 
C...q -> g + q.
C...NB: MQ>0 not yet implemented. Forced absent above.
        ELSEIF (KFLAA.LE.5.AND.KFLB.EQ.21) THEN
          KFLC=KFLA
          Z=ZMAX/(1D0+PYR(0)*(SQRT(ZMAX/ZMIN)-1D0))**2
          WTZ=0.5D0*(1D0+(1D0-Z)**2)*SQRT(Z)
 
C...g -> q + qbar.
        ELSEIF (KFLA.EQ.21.AND.KFLBA.LE.5) THEN
          KFLC=-KFLB
          Z=ZMIN+PYR(0)*(ZMAX-ZMIN)
          WTZ=Z**2+(1D0-Z)**2
C...Massive correction
          IF (MQMASS.NE.0) THEN
            WTZ=WTZ+2D0*Z*(1D0-Z)*RMQ2/PT2
C...Extra safety margin for light sea quark creation
          ELSEIF (KSVCB.LT.0) THEN
            WTZ=WTZ/1.25D0
          ENDIF
 
C...g -> g + g.
        ELSEIF (KFLA.EQ.21.AND.KFLB.EQ.21) THEN
          KFLC=21
          Z=1D0/(1D0+((1D0-ZMIN)/ZMIN)*((1D0-ZMAX)*ZMIN/
     &         (ZMAX*(1D0-ZMIN)))**PYR(0))
          WTZ=(1D0-Z*(1D0-Z))**2
        ENDIF
 
C...Derive Q2 from pT2.
        Q2B=PT2/(1D0-Z)
        IF (KFLBA.GE.4) Q2B=Q2B-RMQ2
 
C...Loopback if outside allowed z range for given pT2.
        RM2C=PYMASS(KFLC)**2
        PT2ADJ=Q2B-Z*(SHTNOW(MI)+Q2B)*(Q2B+RM2C)/SHTNOW(MI)
        IF (PT2ADJ.LT.1D-6) GOTO 230
 
C...Size of phase space and coherence suppression: MSTP(67) and MSTP(62)
C...No modification for very first emission if using ME correction
        MSTP67 = MSTP(67)
        IF (MECOR.GE.1.AND.NISGEN(1,MI).EQ.0.AND.NISGEN(2,MI).EQ.0) THEN
          MSTP67 = 0
        ENDIF
 
C...For 1st branching, limit phase space by s-hat with color-partner
C...(prevent infinite loop by limiting number of NTRY)
        IF (MSTP67.GE.1.AND.NISGEN(JS,MI).EQ.0.AND.NTRY.LE.200) THEN
          MSIDE=1
          IDIP=IMI(JS,MI,1)
C...Use anticolor tag for antiquark, or for gluon half the time
          IF ((KFLB.LT.0.AND.KFLBA.LT.10).OR.
     &         (KFLB.EQ.21.AND.PYR(0).GT.0.5)) MSIDE=2
C...Tag
          MCTAG=MCT(IDIP,MSIDE)
C...Default is to set up phase space using the opposite incoming parton
          JDIP=IMI(3-JS,MI,1)
          NDIP=0

C...Alternatively, look for final-state color partner (pick last if several)
          DO 260 IFS=1,NPART
            MCJ = MCT(IPART(IFS),MSIDE)
            IF (MCJ.NE.MCTAG) GOTO 260
C...Pick last matching final-state partner if several
C...(if no matching final-state partner, defaults back to annihilation)
            KSJ = K(IPART(IFS),1)
            IF (KSJ.GE.1.AND.KSJ.LT.10) THEN
              JDIP=IPART(IFS)
              NDIP=NDIP+1
            ENDIF
  260     CONTINUE

C...Compute momentum transfer: sdip = -t = - (p1 - p2)^2
C...(also works for annihilation since incoming massless, so shat = -(p1 - p2)^2)
          SDIP=ABS(((P(IDIP,4)-P(JDIP,4))**2-(P(IDIP,3)-P(JDIP,3))**2
     &        -(P(IDIP,2)-P(JDIP,2))**2-(P(IDIP,1)-P(JDIP,1))**2))

          IF (MSTP67.EQ.1) THEN
C...1 Option to completely kill radiation above s_dip * PARP(67)
            IF (4D0*PT2.GT.PARP(67)*SDIP) GOTO 230
          ELSE IF (MSTP67.EQ.2) THEN
C...2 Option to allow suppressed unordered radiation above s_dip * PARP(67)
C...  (-> improved power showers?)
            IF (4D0*PT2*PYR(0).GT.PARP(67)*SDIP) GOTO 230
          ENDIF
          
C...For subsequent branchings, loopback if nonordered in angle/rapidity
        ELSE IF (MSTP(62).GE.3.AND.NISGEN(JS,MI).GE.1) THEN
          IF(PT2.GT.((1D0-Z)/(Z*(1D0-ZSAV(JS,MI))))**2*PT2SAV(JS,MI))
     &         GOTO 230
        ENDIF
 
C...Select phi angle of branching at random.
        PHI=PARU(2)*PYR(0)
 
C...Matrix-element corrections for some processes.
        IF (MECOR.GE.1.AND.NISGEN(JS,MI).EQ.0) THEN
          IF (KFLAA.LE.20.AND.KFLBA.LE.20) THEN
            CALL PYMEWT(MECOR,1,Q2B*SHAT/SHTNOW(MI),Z,PHI,WTME)
            WTZ=WTZ*WTME/WTFF
          ELSEIF((KFLA.EQ.21.OR.KFLA.EQ.22).AND.KFLBA.LE.20) THEN
            CALL PYMEWT(MECOR,2,Q2B*SHAT/SHTNOW(MI),Z,PHI,WTME)
            WTZ=WTZ*WTME/WTGF
          ELSEIF(KFLAA.LE.20.AND.(KFLB.EQ.21.OR.KFLB.EQ.22)) THEN
            CALL PYMEWT(MECOR,3,Q2B*SHAT/SHTNOW(MI),Z,PHI,WTME)
            WTZ=WTZ*WTME/WTFG
          ELSEIF(KFLA.EQ.21.AND.KFLB.EQ.21) THEN
            CALL PYMEWT(MECOR,4,Q2B*SHAT/SHTNOW(MI),Z,PHI,WTME)
            WTZ=WTZ*WTME/WTGG
          ENDIF
        ENDIF
 
C...Parton distributions at new pT2 but old x.
        MINT(30)=JS
        CALL PYPDFU(KFBEAM(JS),XB,PT2,XFN)
C...Treat val and cmp separately
        IF (KFLBA.LE.6.AND.KSVCB.LE.0) XFN(KFLB)=XPSVC(KFLB,KSVCB)
        IF (KSVCB.GE.1)
     &       XFN(KFLB)=PYFCMP(YB/VINT(140),YS/VINT(140),MSTP(87))
        XFBN=XFN(KFLB)
        IF(XFBN.LT.1D-20) THEN
          IF(KFLA.EQ.KFLB) THEN
            WTAP(KFLB)=0D0
            GOTO 200
          ELSE
            XFBN=1D-10
            XFN(KFLB)=XFBN
          ENDIF
        ENDIF
        DO 270 KFL=-5,5
          XFB(KFL)=XFN(KFL)
  270   CONTINUE
        XFB(21)=XFN(21)
 
C...Parton distributions at new pT2 and new x.
        XA=XB/Z
        MINT(30)=JS
        CALL PYPDFU(KFBEAM(JS),XA,PT2,XFA)
        IF (KFLBA.LE.5.AND.KFLAA.LE.5) THEN
C...q -> q + g: only consider respective sea, val, or cmp content.
          IF (KSVCB.LE.0) THEN
            XFA(KFLA)=XPSVC(KFLA,KSVCB)
          ELSE
            YA=XA*(1D0-YS)
            XFA(KFLB)=PYFCMP(YA/VINT(140),YS/VINT(140),MSTP(87))
          ENDIF
        ENDIF
        XFAN=XFA(KFLA)
        IF(XFAN.LT.1D-20) THEN
          GOTO 200
        ENDIF
 
C...If weighting fails continue evolution.
        WTTOT=0D0
        IF (MCRQQ.EQ.0) THEN
          WTPDFA=1D0/WTPDF(KFLA)
          WTTOT=WTZ*XFAN/XFBN*WTPDFA
        ELSEIF(MCRQQ.EQ.1) THEN
          WTPDFA=TPM/WPDF0
          WTTOT=WTCRQQ*WTZ*XFAN/XFBN*WTPDFA
          XBEST=TPM/TPM0*XQ0
        ELSEIF(MCRQQ.EQ.2) THEN
C...Force massive quark creation.
          WTTOT=1D0
        ENDIF
 
C...Loop back if trial emission fails.
        IF(WTTOT.GE.0D0.AND.WTTOT.LT.PYR(0)) GOTO 200
        WTACC=((1D0+PT2)/(0.25D0+PT2))**2
        IF(WTTOT.LT.0D0) THEN
          WRITE(CHWT,'(1P,E12.4)') WTTOT
          CALL PYERRM(19,'(PYPTIS:) Weight '//CHWT//' negative')
        ELSEIF(WTTOT.GT.WTACC) THEN
          WRITE(CHWT,'(1P,E12.4)') WTTOT
          IF (PT2.GT.PTEMAX.OR.WTTOT.GE.WTEMAX) THEN
C...Too high weight: write out as error, but do not update error counter
            IF(MSTU(29).EQ.0) MSTU(23)=MSTU(23)-1
            CALL PYERRM(19,
     &         '(PYPTIS:) Weight '//CHWT//' above unity')
            IF (PT2.GT.PTEMAX) PTEMAX=PT2
            IF (WTTOT.GT.WTEMAX) WTEMAX=WTTOT
          ELSE
            CALL PYERRM(9,
     &         '(PYPTIS:) Weight '//CHWT//' above unity')
          ENDIF
C...Useful for debugging but commented out for distribution:
C          print*, 'JS, MI',JS, MI
C          print*, 'PT:',SQRT(PT2), ' MCRQQ',MCRQQ
C          print*, 'A -> B C',KFLA, KFLB, KFLC
C          XFAO=XFBO/WTPDFA
C          print*, 'WT(Z,XFA,XFB)',WTZ, XFAN/XFAO, XFBO/XFBN
        ENDIF
 
C...Special for PT2 = PT2MX (e.g., if two incoming massive quarks 
C...simultaneously reached their creation thresholds) 
        IF (ABS(PT2-PT2MX).LT.0.001) THEN
          IF (PYR(0).GT.0.5) PT2=1.0001*PT2MX
        ENDIF

C...Save acceptable branching.
        IF(PT2.GT.PT2MX) THEN
          MIMX=MINT(36)
          JSMX=JS
          PT2MX=PT2
          KFLAMX=KFLA
          KFLCMX=KFLC
          RM2CMX=RM2C
          Q2BMX=Q2B
          ZMX=Z
          PT2AMX=PT2ADJ
          PHIMX=PHI
        ENDIF
 
C----------------------------------------------------------------------
C...MODE= 1: Accept stored shower branching. Update event record etc.
      ELSEIF (MODE.EQ.1) THEN
        MI=MIMX
        JS=JSMX
        SHAT=SHTNOW(MI)
        SIDE=3D0-2D0*JS
C...Shift down rest of event record to make room for insertion.
        IT=IMISEP(MI)+1
        IM=IT+1
        IS=IMI(JS,MI,1)
        DO 290 I=N,IT,-1
          IF (K(I,3).GE.IT) K(I,3)=K(I,3)+2
          KT1=K(I,4)/MSTU(5)**2
          KT2=K(I,5)/MSTU(5)**2
          ID1=MOD(K(I,4),MSTU(5))
          ID2=MOD(K(I,5),MSTU(5))
          IM1=MOD(K(I,4)/MSTU(5),MSTU(5))
          IM2=MOD(K(I,5)/MSTU(5),MSTU(5))
          IF (ID1.GE.IT) ID1=ID1+2
          IF (ID2.GE.IT) ID2=ID2+2
          IF (IM1.GE.IT) IM1=IM1+2
          IF (IM2.GE.IT) IM2=IM2+2
          K(I,4)=KT1*MSTU(5)**2+IM1*MSTU(5)+ID1
          K(I,5)=KT2*MSTU(5)**2+IM2*MSTU(5)+ID2
          DO 280 IX=1,5
            K(I+2,IX)=K(I,IX)
            P(I+2,IX)=P(I,IX)
            V(I+2,IX)=V(I,IX)
  280     CONTINUE
          MCT(I+2,1)=MCT(I,1)
          MCT(I+2,2)=MCT(I,2)
  290   CONTINUE
        N=N+2
C...Also update shifted-down pointers in IMI, IMISEP, and IPART.
        DO 300 JI=1,MINT(31)
          IF (IMI(1,JI,1).GE.IT) IMI(1,JI,1)=IMI(1,JI,1)+2
          IF (IMI(1,JI,2).GE.IT) IMI(1,JI,2)=IMI(1,JI,2)+2
          IF (IMI(2,JI,1).GE.IT) IMI(2,JI,1)=IMI(2,JI,1)+2
          IF (IMI(2,JI,2).GE.IT) IMI(2,JI,2)=IMI(2,JI,2)+2
          IF (JI.GE.MI) IMISEP(JI)=IMISEP(JI)+2
C...Also update companion pointers to the present mother.
          IF (IMI(JS,JI,2).EQ.IS) IMI(JS,JI,2)=IM
  300   CONTINUE
        DO 310 IFS=1,NPART
          IF (IPART(IFS).GE.IT) IPART(IFS)=IPART(IFS)+2
  310   CONTINUE
C...Zero entries dedicated for new timelike and mother partons.
        DO 330 I=IT,IT+1
          DO 320 J=1,5
            K(I,J)=0
            P(I,J)=0D0
            V(I,J)=0D0
  320     CONTINUE
          MCT(I,1)=0
          MCT(I,2)=0
  330   CONTINUE
 
C...Define timelike and new mother partons. History.
        K(IT,1)=3
        K(IT,2)=KFLCMX
        K(IM,1)=14
        K(IM,2)=KFLAMX
        K(IS,3)=IM
        K(IT,3)=IM
C...Set mother origin = side.
        K(IM,3)=MINT(83)+JS+2
        IF(MI.GE.2) K(IM,3)=MINT(83)+JS
 
C...Define colour flow of branching.
        IM1=IM
        IM2=IM
C...q -> q + gamma.
        IF(K(IT,2).EQ.22) THEN
          K(IT,1)=1
          ID1=IS
          ID2=IS
C...q -> q + g.
        ELSEIF(K(IM,2).GT.0.AND.K(IM,2).LE.5.AND.K(IT,2).EQ.21) THEN
          ID1=IT
          ID2=IS
C...q -> g + q.
        ELSEIF(K(IM,2).GT.0.AND.K(IM,2).LE.5) THEN
          ID1=IS
          ID2=IT
C...qbar -> qbar + g.
        ELSEIF(K(IM,2).LT.0.AND.K(IM,2).GE.-5.AND.K(IT,2).EQ.21) THEN
          ID1=IS
          ID2=IT
C...qbar -> g + qbar.
        ELSEIF(K(IM,2).LT.0.AND.K(IM,2).GE.-5) THEN
          ID1=IT
          ID2=IS
C...g -> g + g; g -> q + qbar..
        ELSEIF((K(IT,2).EQ.21.AND.PYR(0).GT.0.5D0).OR.K(IT,2).LT.0) THEN
          ID1=IS
          ID2=IT
        ELSE
          ID1=IT
          ID2=IS
        ENDIF
        IF(IM1.EQ.IM) K(IM1,4)=K(IM1,4)+ID1
        IF(IM2.EQ.IM) K(IM2,5)=K(IM2,5)+ID2
        K(ID1,4)=K(ID1,4)+MSTU(5)*IM1
        K(ID2,5)=K(ID2,5)+MSTU(5)*IM2
        IF(ID1.NE.ID2) THEN
          K(ID1,5)=K(ID1,5)+MSTU(5)*ID2
          K(ID2,4)=K(ID2,4)+MSTU(5)*ID1
        ENDIF
        IF(K(IT,1).EQ.1) THEN
          K(IT,4)=0
          K(IT,5)=0
        ENDIF
C...Update IMI and colour tag arrays.
        IMI(JS,MI,1)=IM
        DO 340 MC=1,2
          MCT(IT,MC)=0
          MCT(IM,MC)=0
  340   CONTINUE
        DO 350 JCS=4,5
          KCS=JCS
C...If mother flag not yet set for spacelike parton, trace it.
          IF (K(IS,KCS)/MSTU(5)**2.LE.1) CALL PYCTTR(IS,-KCS,IM)
          IF(MINT(51).NE.0) RETURN
  350   CONTINUE
        DO 360 JCS=4,5
          KCS=JCS
C...If mother flag not yet set for timelike parton, trace it.
          IF (K(IT,KCS)/MSTU(5)**2.LE.1) CALL PYCTTR(IT,KCS,IM)
          IF(MINT(51).NE.0) RETURN
  360   CONTINUE
 
C...Boost recoiling parton to compensate for Q2 scale.
        BETAZ=SIDE*(1D0-(1D0+Q2BMX/SHAT)**2)/
     &  (1D0+(1D0+Q2BMX/SHAT)**2)
        IR=IMI(3-JS,MI,1)
        CALL PYROBO(IR,IR,0D0,0D0,0D0,0D0,BETAZ)
 
C...Define system to be rotated and boosted
C...(not including the 2 just added partons)
C...(but including the docu lines for first interaction)
        IMIN=IMISEP(MI-1)+1
        IF (MI.EQ.1) IMIN=MINT(83)+5
        IMAX=IMISEP(MI)-2
 
C...Rotate back system in phi to compensate for subsequent rotation.
        CALL PYROBO(IMIN,IMAX,0D0,-PHIMX,0D0,0D0,0D0)
 
C...Define kinematics of new partons in old frame.
        IMAX=IMISEP(MI)
        P(IM,1)=SQRT(PT2AMX)*SHAT/(ZMX*(SHAT+Q2BMX))
        P(IM,3)=0.5D0*SQRT(SHAT)*((SHAT-Q2BMX)/((SHAT
     &       +Q2BMX)*ZMX)+(Q2BMX+RM2CMX)/SHAT)*SIDE
        P(IM,4)=SQRT(P(IM,1)**2+P(IM,3)**2)
        P(IT,1)=P(IM,1)
        P(IT,3)=P(IM,3)-0.5D0*(SHAT+Q2BMX)/SQRT(SHAT)*SIDE
        P(IT,4)=SQRT(P(IT,1)**2+P(IT,3)**2+RM2CMX)
        P(IT,5)=SQRT(RM2CMX)
 
C...Update internal line, now spacelike
        P(IS,1)=P(IM,1)-P(IT,1)
        P(IS,2)=P(IM,2)-P(IT,2)
        P(IS,3)=P(IM,3)-P(IT,3)
        P(IS,4)=P(IM,4)-P(IT,4)
        P(IS,5)=P(IS,4)**2-P(IS,1)**2-P(IS,2)**2-P(IS,3)**2
C...Represent spacelike virtualities as -sqrt(abs(Q2)) .
        IF (P(IS,5).LT.0D0) THEN
          P(IS,5)=-SQRT(ABS(P(IS,5)))
        ELSE
          P(IS,5)=SQRT(P(IS,5))
        ENDIF
 
C...Boost entire system and rotate to new frame.
C...(including docu lines)
        BETAX=(P(IM,1)+P(IR,1))/(P(IM,4)+P(IR,4))
        BETAZ=(P(IM,3)+P(IR,3))/(P(IM,4)+P(IR,4))
        IF(BETAX**2+BETAZ**2.GE.1D0) THEN
          CALL PYERRM(1,'(PYPTIS:) boost bigger than unity')
          MINT(51)=1
          IFAIL=-1
          RETURN
        ENDIF
        CALL PYROBO(IMIN,IMAX,0D0,0D0,-BETAX,0D0,-BETAZ)
        I1=IMI(1,MI,1)
        THETA=PYANGL(P(I1,3),P(I1,1))
        CALL PYROBO(IMIN,IMAX,-THETA,PHIMX,0D0,0D0,0D0)
 
C...Global statistics.
        MINT(352)=MINT(352)+1
        VINT(352)=VINT(352)+SQRT(P(IT,1)**2+P(IT,2)**2)
        IF (MINT(352).EQ.1) VINT(357)=SQRT(P(IT,1)**2+P(IT,2)**2)
 
C...Add parton with relevant pT scale for timelike shower.
        IF (K(IT,2).NE.22) THEN
          NPART=NPART+1
          IPART(NPART)=IT
          PTPART(NPART)=SQRT(PT2AMX)
        ENDIF
 
C...Update saved variables.
        SHTNOW(MIMX)=SHTNOW(MIMX)/ZMX
        NISGEN(JSMX,MIMX)=NISGEN(JSMX,MIMX)+1
        XMI(JSMX,MIMX)=XMI(JSMX,MIMX)/ZMX
        PT2SAV(JSMX,MIMX)=PT2MX
        ZSAV(JS,MIMX)=ZMX
 
        KSA=IABS(K(IS,2))
        KMA=IABS(K(IM,2))
        IF (KSA.EQ.21.AND.KMA.GE.1.AND.KMA.LE.5) THEN
C...Gluon reconstructs to quark.
C...Decide whether newly created quark is valence or sea:
          MINT(30)=JS
          CALL PYPTMI(2,PT2NOW,PTDUM1,PTDUM2,IFAIL)
          IF(MINT(51).NE.0) RETURN
        ENDIF
        IF(KSA.GE.1.AND.KSA.LE.5.AND.KMA.EQ.21) THEN
C...Quark reconstructs to gluon.
C...Now some guy may have lost his companion. Check.
          ICMP=IMI(JS,MI,2)
          IF (ICMP.GT.0) THEN
            CALL PYERRM(9,'(PYPTIS:) Sorry, companion quark radiated'
     &           //' away. Cannot handle that yet. Giving up.')
            MINT(51)=1
            RETURN
          ELSEIF(ICMP.LT.0) THEN
C...A sea quark with companion still in BR was reconstructed to a gluon.
C...Companion should now be removed from the beam remnant.
C...(Momentum integral is automatically updated in next call to PYPDFU.)
            ICMP=-ICMP
            IFL=-K(IS,2)
            DO 380 JCMP=ICMP,NVC(JS,IFL)-1
              XASSOC(JS,IFL,JCMP)=XASSOC(JS,IFL,JCMP+1)
              DO 370 JI=1,MINT(31)
                KMI=-IMI(JS,JI,2)
                JFL=-K(IMI(JS,JI,1),2)
                IF (KMI.EQ.JCMP+1.AND.JFL.EQ.IFL) IMI(JS,JI,2)=IMI(JS,JI
     &               ,2)+1
  370         CONTINUE
  380       CONTINUE
            NVC(JS,IFL)=NVC(JS,IFL)-1
          ENDIF
C...Set gluon IMI(JS,MI,2) = 0.
          IMI(JS,MI,2)=0
        ELSEIF(KSA.GE.1.AND.KSA.LE.5.AND.KMA.NE.21) THEN
C...Quark reconstructing to quark. If sea with companion still in BR
C...then update associated x value.
C...(Momentum integral is automatically updated in next call to PYPDFU.)
          IF (IMI(JS,MI,2).LT.0) THEN
            ICMP=-IMI(JS,MI,2)
            IFL=-K(IS,2)
            XASSOC(JS,IFL,ICMP)=XMI(JSMX,MIMX)
          ENDIF
        ENDIF
 
      ENDIF
 
C...If reached this point, normal exit.
  390 IFAIL=0
 
      RETURN
      END
