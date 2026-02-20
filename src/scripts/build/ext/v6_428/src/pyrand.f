 
C*********************************************************************
 
C...PYRAND
C...Generates quantities characterizing the high-pT scattering at the
C...parton level according to the matrix elements. Chooses incoming,
C...reacting partons, their momentum fractions and one of the possible
C...subprocesses.
 
      SUBROUTINE PYRAND
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
 
C...User process initialization and event commonblocks.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPRUP/,/HEPEUP/
 
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYTCCO/COEFX(194:380,2)
      COMMON/TCPARA/IRES,JRES,XMAS(3),XWID(3),YMAS(2),YWID(2)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,/PYINT1/,
     &/PYINT2/,/PYINT3/,/PYINT4/,/PYINT5/,/PYINT7/,/PYMSSM/,/PYTCCO/,
     &/TCPARA/
C...Local arrays.
      DIMENSION XPQ(-25:25),PMM(2),PDIF(4),BHAD(4),PMMN(2)
 
C...Parameters and data used in elastic/diffractive treatment.
      DATA EPS/0.0808D0/, ALP/0.25D0/, CRES/2D0/, PMRC/1.062D0/,
     &SMP/0.880D0/, BHAD/2.3D0,1.4D0,1.4D0,0.23D0/
 
C...Initial values, specifically for (first) semihard interaction.
      MINT(10)=0
      MINT(17)=0
      MINT(18)=0
      VINT(143)=1D0
      VINT(144)=1D0
      VINT(157)=0D0
      VINT(158)=0D0
      MFAIL=0
      IF(MSTP(171).EQ.1.AND.MSTP(172).EQ.2) MFAIL=1
      ISUB=0
      ISTSB=0
      LOOP=0
  100 LOOP=LOOP+1
      MINT(51)=0
      MINT(143)=1
      VINT(97)=1D0
 
C...Start by assuming incoming photon is entering subprocess.
      IF(MINT(11).EQ.22) THEN
         MINT(15)=22
         VINT(307)=VINT(3)**2
      ENDIF
      IF(MINT(12).EQ.22) THEN
         MINT(16)=22
         VINT(308)=VINT(4)**2
      ENDIF
      MINT(103)=MINT(11)
      MINT(104)=MINT(12)
 
C...Choice of process type - first event of pileup.
      INMULT=0
      IF(MINT(82).EQ.1.AND.ISUB.GE.91.AND.ISUB.LE.96) THEN
      ELSEIF(MINT(82).EQ.1) THEN
 
C...For gamma-p or gamma-gamma first pick between alternatives.
        IGA=0
        IF(MINT(121).GT.1) CALL PYSAVE(4,IGA)
        MINT(122)=IGA
 
C...For real gamma + gamma with different nature, flip at random.
        IF(MINT(11).EQ.22.AND.MINT(12).EQ.22.AND.MINT(123).GE.4.AND.
     &  MSTP(14).LE.10.AND.PYR(0).GT.0.5D0) THEN
          MINTSV=MINT(41)
          MINT(41)=MINT(42)
          MINT(42)=MINTSV
          MINTSV=MINT(45)
          MINT(45)=MINT(46)
          MINT(46)=MINTSV
          MINTSV=MINT(107)
          MINT(107)=MINT(108)
          MINT(108)=MINTSV
          IF(MINT(47).EQ.2.OR.MINT(47).EQ.3) MINT(47)=5-MINT(47)
        ENDIF
 
C...Pick process type, possibly by user process machinery.
C...(If the latter, also event will be picked here.)
        IF(MINT(111).GE.11.AND.IABS(IDWTUP).EQ.2.AND.LOOP.GE.2) THEN
          CALL UPEVNT
          CALL PYUPRE
        ELSEIF(MINT(111).GE.11.AND.IABS(IDWTUP).GE.3) THEN
          CALL UPEVNT
          CALL PYUPRE
          ISUB=0
  110     ISUB=ISUB+1
          IF((ISET(ISUB).NE.11.OR.KFPR(ISUB,2).NE.IDPRUP).AND.
     &    ISUB.LT.500) GOTO 110
        ELSE
          RSUB=XSEC(0,1)*PYR(0)
          DO 120 I=1,500
            IF(MSUB(I).NE.1.OR.I.EQ.96) GOTO 120
            ISUB=I
            RSUB=RSUB-XSEC(I,1)
            IF(RSUB.LE.0D0) GOTO 130
  120     CONTINUE
  130     IF(ISUB.EQ.95) ISUB=96
          IF(ISUB.EQ.96) INMULT=1
          IF(ISET(ISUB).EQ.11) THEN
            IDPRUP=KFPR(ISUB,2)
            CALL UPEVNT
            CALL PYUPRE
          ENDIF
        ENDIF
 
C...Choice of inclusive process type - pileup events.
      ELSEIF(MINT(82).GE.2.AND.ISUB.EQ.0) THEN
        RSUB=VINT(131)*PYR(0)
        ISUB=96
        IF(RSUB.GT.SIGT(0,0,5)) ISUB=94
        IF(RSUB.GT.SIGT(0,0,5)+SIGT(0,0,4)) ISUB=93
        IF(RSUB.GT.SIGT(0,0,5)+SIGT(0,0,4)+SIGT(0,0,3)) ISUB=92
        IF(RSUB.GT.SIGT(0,0,5)+SIGT(0,0,4)+SIGT(0,0,3)+SIGT(0,0,2))
     &  ISUB=91
        IF(ISUB.EQ.96) INMULT=1
      ENDIF
 
C...Choice of photon energy and flux factor inside lepton.
      IF(MINT(141).NE.0.OR.MINT(142).NE.0) THEN
        CALL PYGAGA(3,WTGAGA)
        IF(ISUB.GE.131.AND.ISUB.LE.140) THEN
          CKIN(3)=MAX(VINT(285),VINT(154))
          CKIN(1)=2D0*CKIN(3)
        ENDIF
C...When necessary set direct/resolved photon by hand.
      ELSEIF(MINT(15).EQ.22.OR.MINT(16).EQ.22) THEN
        IF(MINT(15).EQ.22.AND.MINT(41).EQ.2) MINT(15)=0
        IF(MINT(16).EQ.22.AND.MINT(42).EQ.2) MINT(16)=0
      ENDIF
 
C...Restrict direct*resolved processes to pTmin >= Q,
C...to avoid doublecounting  with DIS.
      IF(MSTP(18).EQ.3.AND.ISUB.GE.131.AND.ISUB.LE.136) THEN
        IF(MINT(15).EQ.22) THEN
          CKIN(3)=MAX(VINT(285),VINT(154),ABS(VINT(3)))
        ELSE
          CKIN(3)=MAX(VINT(285),VINT(154),ABS(VINT(4)))
        ENDIF
        CKIN(1)=2D0*CKIN(3)
      ENDIF
 
C...Set up for multiple interactions (may include impact parameter).
      IF(INMULT.EQ.1) THEN
        IF(MINT(35).LE.1) CALL PYMULT(2)
        IF(MINT(35).GE.2) CALL PYMIGN(2)
      ENDIF
 
C...Loopback point for minimum bias in photon physics.
      LOOP2=0
  140 LOOP2=LOOP2+1
      IF(MINT(82).EQ.1) NGEN(0,1)=NGEN(0,1)+MINT(143)
      IF(MINT(82).EQ.1) NGEN(ISUB,1)=NGEN(ISUB,1)+MINT(143)
      IF(ISUB.EQ.96.AND.LOOP2.EQ.1.AND.MINT(82).EQ.1)
     &NGEN(97,1)=NGEN(97,1)+MINT(143)
      MINT(1)=ISUB
      ISTSB=ISET(ISUB)
 
C...Random choice of flavour for some SUSY processes.
      IF(ISUB.GE.201.AND.ISUB.LE.301) THEN
C...~e_L ~nu_e or ~mu_L ~nu_mu.
        IF(ISUB.EQ.210) THEN
          KFPR(ISUB,1)=KSUSY1+11+2*INT(0.5D0+PYR(0))
          KFPR(ISUB,2)=KFPR(ISUB,1)+1
C...~nu_e ~nu_e(bar) or ~nu_mu ~nu_mu(bar).
        ELSEIF(ISUB.EQ.213) THEN
          KFPR(ISUB,1)=KSUSY1+12+2*INT(0.5D0+PYR(0))
          KFPR(ISUB,2)=KFPR(ISUB,1)
C...~q ~chi/~g; ~q = ~d, ~u, ~s, ~c or ~b.
        ELSEIF(ISUB.GE.246.AND.ISUB.LE.259.AND.ISUB.NE.255.AND.
     &  ISUB.NE.257) THEN
          IF(ISUB.GE.258) THEN
            RKF=4D0
          ELSE
            RKF=5D0
          ENDIF
          IF(MOD(ISUB,2).EQ.0) THEN
            KFPR(ISUB,1)=KSUSY1+1+INT(RKF*PYR(0))
          ELSE
            KFPR(ISUB,1)=KSUSY2+1+INT(RKF*PYR(0))
          ENDIF
C...~q1 ~q2; ~q = ~d, ~u, ~s, or ~c.
        ELSEIF(ISUB.GE.271.AND.ISUB.LE.276) THEN
          IF(ISUB.EQ.271.OR.ISUB.EQ.274) THEN
            KSU1=KSUSY1
            KSU2=KSUSY1
          ELSEIF(ISUB.EQ.272.OR.ISUB.EQ.275) THEN
            KSU1=KSUSY2
            KSU2=KSUSY2
          ELSEIF(PYR(0).LT.0.5D0) THEN
            KSU1=KSUSY1
            KSU2=KSUSY2
          ELSE
            KSU1=KSUSY2
            KSU2=KSUSY1
          ENDIF
          KFPR(ISUB,1)=KSU1+1+INT(4D0*PYR(0))
          KFPR(ISUB,2)=KSU2+1+INT(4D0*PYR(0))
C...~q ~q(bar);  ~q = ~d, ~u, ~s, or ~c.
        ELSEIF(ISUB.EQ.277.OR.ISUB.EQ.279) THEN
          KFPR(ISUB,1)=KSUSY1+1+INT(4D0*PYR(0))
          KFPR(ISUB,2)=KFPR(ISUB,1)
        ELSEIF(ISUB.EQ.278.OR.ISUB.EQ.280) THEN
          KFPR(ISUB,1)=KSUSY2+1+INT(4D0*PYR(0))
          KFPR(ISUB,2)=KFPR(ISUB,1)
C...~q1 ~q2; ~q = ~d, ~u, ~s, or ~c.
        ELSEIF(ISUB.GE.281.AND.ISUB.LE.286) THEN
          IF(ISUB.EQ.281.OR.ISUB.EQ.284) THEN
            KSU1=KSUSY1
            KSU2=KSUSY1
          ELSEIF(ISUB.EQ.282.OR.ISUB.EQ.285) THEN
            KSU1=KSUSY2
            KSU2=KSUSY2
          ELSEIF(PYR(0).LT.0.5D0) THEN
            KSU1=KSUSY1
            KSU2=KSUSY2
          ELSE
            KSU1=KSUSY2
            KSU2=KSUSY1
          ENDIF
          IF(ISUB.EQ.281.OR.ISUB.LE.283) THEN
            RKF=5D0
          ELSE
            RKF=4D0
          ENDIF
          KFPR(ISUB,2)=KSU2+1+INT(RKF*PYR(0))
        ENDIF
      ENDIF
 
C...Random choice of flavours for some UED processes
c...The production processes can generate a doublet pair,
c...a singlet pair, or a doublet + singlet.
      IF(ISUB.EQ.313)THEN
C...q + q -> q*_Di + q*_Dj, q*_Si + q*_Sj
         IF(PYR(0).LE.0.1)THEN
            KFPR(ISUB,1)=5100001
         ELSE
            KFPR(ISUB,1)=5100002
         ENDIF
         KFPR(ISUB,2)=KFPR(ISUB,1)
      ELSEIF(ISUB.EQ.314.OR.ISUB.EQ.315)THEN
C...g + g -> q*_D + q*_Dbar, q*_S + q*_Sbar
C...q + qbar -> q*_D + q*_Dbar, q*_S + q*_Sbar
         IF(PYR(0).LE.0.1)THEN
            KFPR(ISUB,1)=5100001
         ELSE
            KFPR(ISUB,1)=5100002
         ENDIF
         KFPR(ISUB,2)=-KFPR(ISUB,1)
      ELSEIF(ISUB.EQ.316)THEN
C...qi + qbarj -> q*_Di + q*_Sbarj
         IF(PYR(0).LE.0.5)THEN
            KFPR(ISUB,1)=5100001
c Changed from private pythia6410_ued code
c            KFPR(ISUB,2)=-5010001
            KFPR(ISUB,2)=-6100002
         ELSE
            KFPR(ISUB,1)=5100002
c Changed from private pythia6410_ued code
c            KFPR(ISUB,2)=-5010002
            KFPR(ISUB,2)=-6100001
         ENDIF
      ELSEIF(ISUB.EQ.317)THEN
C...qi + qbarj -> q*_Di + q*_Dbarj, q*_Si + q*_Dbarj
         IF(PYR(0).LE.0.5)THEN
            KFPR(ISUB,1)=5100001
            KFPR(ISUB,2)=-5100002
         ELSE
            KFPR(ISUB,1)=5100002
            KFPR(ISUB,2)=-5100001
         ENDIF
      ELSEIF(ISUB.EQ.318)THEN
C...qi + qj -> q*_Di + q*_Sj
         IF(PYR(0).LE.0.5)THEN
            KFPR(ISUB,1)=5100001
            KFPR(ISUB,2)=6100002
         ELSE
            KFPR(ISUB,1)=5100002
            KFPR(ISUB,2)=6100001
         ENDIF
      ENDIF

C...Find resonances (explicit or implicit in cross-section).
      MINT(72)=0
      KFR1=0
      IF(ISTSB.EQ.1.OR.ISTSB.EQ.3.OR.ISTSB.EQ.5) THEN
        KFR1=KFPR(ISUB,1)
      ELSEIF(ISUB.EQ.24.OR.ISUB.EQ.25.OR.ISUB.EQ.110.OR.ISUB.EQ.165.OR.
     &  ISUB.EQ.171.OR.ISUB.EQ.176) THEN
        KFR1=23
      ELSEIF(ISUB.EQ.23.OR.ISUB.EQ.26.OR.ISUB.EQ.166.OR.ISUB.EQ.172.OR.
     &  ISUB.EQ.177) THEN
        KFR1=24
      ELSEIF(ISUB.GE.71.AND.ISUB.LE.77) THEN
        KFR1=25
        IF(MSTP(46).EQ.5) THEN
          KFR1=89
          PMAS(89,1)=PARP(45)
          PMAS(89,2)=PARP(45)**3/(96D0*PARU(1)*PARP(47)**2)
        ENDIF
      ELSEIF(ISUB.EQ.481) THEN
        KFR1=9900001
      ENDIF
      CKMX=CKIN(2)
      IF(CKMX.LE.0D0) CKMX=VINT(1)
      KCR1=PYCOMP(KFR1)
      IF(KCR1.EQ.0) KFR1=0
      IF(KFR1.NE.0) THEN
        IF(CKIN(1).GT.PMAS(KCR1,1)+20D0*PMAS(KCR1,2).OR.
     &  CKMX.LT.PMAS(KCR1,1)-20D0*PMAS(KCR1,2)) KFR1=0
      ENDIF
      IF(KFR1.NE.0) THEN
        TAUR1=PMAS(KCR1,1)**2/VINT(2)
        GAMR1=PMAS(KCR1,1)*PMAS(KCR1,2)/VINT(2)
        MINT(72)=1
        MINT(73)=KFR1
        VINT(73)=TAUR1
        VINT(74)=GAMR1
      ENDIF
      KFR2=0
      KFR3=0
      IF(ISUB.EQ.141.OR.ISUB.EQ.194.OR.ISUB.EQ.195.OR.
     $(ISUB.GE.361.AND.ISUB.LE.380))
     $THEN
        KFR2=23
        IF(ISUB.EQ.141) THEN
          KCR2=PYCOMP(KFR2)
          IF(CKIN(1).GT.PMAS(KCR2,1)+20D0*PMAS(KCR2,2).OR.
     &     CKMX.LT.PMAS(KCR2,1)-20D0*PMAS(KCR2,2)) THEN
            KFR2=0
          ELSE
            TAUR2=PMAS(KCR2,1)**2/VINT(2)            
            GAMR2=PMAS(KCR2,1)*PMAS(KCR2,2)/VINT(2)
            MINT(72)=2
            MINT(74)=KFR2
            VINT(75)=TAUR2
            VINT(76)=GAMR2
          ENDIF
C...3 resonances at work:   rho, omega, a
        ELSEIF(ISUB.EQ.194.OR.(ISUB.GE.361.AND.ISUB.LE.368)
     &     .OR.ISUB.EQ.379.OR.ISUB.EQ.380) THEN
          MINT(72)=IRES
          IF(IRES.GE.1) THEN
            VINT(73)=XMAS(1)**2/VINT(2)
            VINT(74)=XMAS(1)*XWID(1)/VINT(2)
            TAUR1=VINT(73)
            GAMR1=VINT(74)
            KFR1=1
          ENDIF
          IF(IRES.GE.2) THEN
            VINT(75)=XMAS(2)**2/VINT(2)
            VINT(76)=XMAS(2)*XWID(2)/VINT(2)
            TAUR2=VINT(75)
            GAMR2=VINT(76)
            KFR2=2
          ENDIF
          IF(IRES.EQ.3) THEN
            VINT(77)=XMAS(3)**2/VINT(2)
            VINT(78)=XMAS(3)*XWID(3)/VINT(2)
            TAUR3=VINT(77)
            GAMR3=VINT(78)
            KFR3=3
          ENDIF
C...Charged current:  rho+- and a+-
        ELSEIF(ISUB.EQ.195.OR.ISUB.GE.370.AND.ISUB.LE.378) THEN
          MINT(72)=IRES
          IF(JRES.GE.1) THEN
            VINT(73)=YMAS(1)**2/VINT(2)
            VINT(74)=YMAS(1)*YWID(1)/VINT(2)
            KFR1=1
            TAUR1=VINT(73)
            GAMR1=VINT(74)
          ENDIF
          IF(JRES.GE.2) THEN
            VINT(75)=YMAS(2)**2/VINT(2)
            VINT(76)=YMAS(2)*YWID(2)/VINT(2)
            KFR2=2
            TAUR2=VINT(73)
            GAMR2=VINT(74)
          ENDIF
          KFR3=0
        ENDIF
        IF(ISUB.NE.141) THEN
          IF(KFR3.NE.0.AND.KFR2.NE.0.AND.KFR1.NE.0) THEN

          ELSEIF(KFR1.NE.0.AND.KFR2.NE.0) THEN
            MINT(72)=2
          ELSEIF(KFR1.NE.0.AND.KFR3.NE.0) THEN
            MINT(72)=2
            MINT(74)=KFR3
            VINT(75)=TAUR3
            VINT(76)=GAMR3
          ELSEIF(KFR2.NE.0.AND.KFR3.NE.0) THEN
            MINT(72)=2
            MINT(73)=KFR2
            VINT(73)=TAUR2
            VINT(74)=GAMR2
            MINT(74)=KFR3
            VINT(75)=TAUR3
            VINT(76)=GAMR3
          ELSEIF(KFR1.NE.0) THEN
            MINT(72)=1
          ELSEIF(KFR2.NE.0) THEN
            MINT(72)=1
            MINT(73)=KFR2
            VINT(73)=TAUR2
            VINT(74)=GAMR2
          ELSEIF(KFR3.NE.0) THEN
            MINT(72)=1
            MINT(73)=KFR3
            VINT(73)=TAUR3
            VINT(74)=GAMR3
          ELSE
            MINT(72)=0
          ENDIF
        ELSE
          IF(KFR2.NE.0.AND.KFR1.NE.0) THEN

          ELSEIF(KFR2.NE.0) THEN
            KFR1=KFR2
            TAUR1=TAUR2
            GAMR1=GAMR2
            MINT(72)=1
            MINT(73)=KFR1
            VINT(73)=TAUR1
            VINT(74)=GAMR1
            KFR2=0
          ELSE
            MINT(72)=0
          ENDIF
        ENDIF
      ENDIF
 
C...Find product masses and minimum pT of process,
C...optionally with broadening according to a truncated Breit-Wigner.
      VINT(63)=0D0
      VINT(64)=0D0
      MINT(71)=0
      VINT(71)=CKIN(3)
      IF(MINT(82).GE.2) VINT(71)=0D0
      VINT(80)=1D0
      IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN
        NBW=0
        DO 160 I=1,2
          PMMN(I)=0D0
          IF(KFPR(ISUB,I).EQ.0) THEN
          ELSEIF(MSTP(42).LE.0.OR.PMAS(PYCOMP(KFPR(ISUB,I)),2).LT.
     &      PARP(41)) THEN
            VINT(62+I)=PMAS(PYCOMP(KFPR(ISUB,I)),1)**2
          ELSE
            NBW=NBW+1
C...This prevents SUSY/t particles from becoming too light.
            KFLW=KFPR(ISUB,I)
            IF(KFLW/KSUSY1.EQ.1.OR.KFLW/KSUSY1.EQ.2) THEN
              KCW=PYCOMP(KFLW)
              PMMN(I)=PMAS(KCW,1)
              DO 150 IDC=MDCY(KCW,2),MDCY(KCW,2)+MDCY(KCW,3)-1
                IF(MDME(IDC,1).GT.0.AND.BRAT(IDC).GT.1E-4) THEN
                  PMSUM=PMAS(PYCOMP(KFDP(IDC,1)),1)+
     &            PMAS(PYCOMP(KFDP(IDC,2)),1)
                  IF(KFDP(IDC,3).NE.0) PMSUM=PMSUM+
     &            PMAS(PYCOMP(KFDP(IDC,3)),1)
                  PMMN(I)=MIN(PMMN(I),PMSUM)
                ENDIF
  150         CONTINUE
            ELSEIF(KFLW.EQ.6) THEN
              PMMN(I)=PMAS(24,1)+PMAS(5,1)
            ENDIF
          ENDIF
  160   CONTINUE
        IF(NBW.GE.1) THEN
          CKIN41=CKIN(41)
          CKIN43=CKIN(43)
          CKIN(41)=MAX(PMMN(1),CKIN(41))
          CKIN(43)=MAX(PMMN(2),CKIN(43))
          CALL PYOFSH(4,0,KFPR(ISUB,1),KFPR(ISUB,2),0D0,PQM3,PQM4)
          CKIN(41)=CKIN41
          CKIN(43)=CKIN43
          IF(MINT(51).EQ.1) THEN
            IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
            IF(MFAIL.EQ.1) THEN
              MSTI(61)=1
              RETURN
            ENDIF
            GOTO 100
          ENDIF
          VINT(63)=PQM3**2
          VINT(64)=PQM4**2
        ENDIF
        IF(MIN(VINT(63),VINT(64)).LT.CKIN(6)**2) MINT(71)=1
        IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5))
      ENDIF
 
C...Prepare for additional variable choices in 2 -> 3.
      IF(ISTSB.EQ.5) THEN
        VINT(201)=0D0
        IF(KFPR(ISUB,2).GT.0) VINT(201)=PMAS(PYCOMP(KFPR(ISUB,2)),1)
        VINT(206)=VINT(201)
        IF(ISUB.EQ.401.OR.ISUB.EQ.402) VINT(206)=PMAS(5,1)
        VINT(204)=PMAS(23,1)
        IF(ISUB.EQ.124.OR.ISUB.EQ.174.OR.ISUB.EQ.179.OR.ISUB.EQ.351)
     &   VINT(204)=PMAS(24,1) 
        IF(ISUB.EQ.352) VINT(204)=PMAS(PYCOMP(9900024),1)
        IF(ISUB.EQ.121.OR.ISUB.EQ.122.OR.ISUB.EQ.181.OR.ISUB.EQ.182.OR.
     &    ISUB.EQ.186.OR.ISUB.EQ.187.OR.ISUB.EQ.401.OR.ISUB.EQ.402)
     &         VINT(204)=VINT(201)
        VINT(209)=VINT(204)
          IF(ISUB.EQ.401.OR.ISUB.EQ.402) VINT(209)=VINT(206)
      ENDIF
 
C...Select incoming VDM particle (rho/omega/phi/J/psi).
      IF(ISTSB.NE.0.AND.(MINT(101).GE.2.OR.MINT(102).GE.2).AND.
     &(MINT(123).EQ.2.OR.MINT(123).EQ.3.OR.MINT(123).EQ.7)) THEN
        VRN=PYR(0)*SIGT(0,0,5)
        IF(MINT(101).LE.1) THEN
          I1MN=0
          I1MX=0
        ELSE
          I1MN=1
          I1MX=MINT(101)
        ENDIF
        IF(MINT(102).LE.1) THEN
          I2MN=0
          I2MX=0
        ELSE
          I2MN=1
          I2MX=MINT(102)
        ENDIF
        DO 180 I1=I1MN,I1MX
          KFV1=110*I1+3
          DO 170 I2=I2MN,I2MX
            KFV2=110*I2+3
            VRN=VRN-SIGT(I1,I2,5)
            IF(VRN.LE.0D0) GOTO 190
  170     CONTINUE
  180   CONTINUE
  190   IF(MINT(101).GE.2) MINT(103)=KFV1
        IF(MINT(102).GE.2) MINT(104)=KFV2
      ENDIF
 
      IF(ISTSB.EQ.0) THEN
C...Elastic scattering or single or double diffractive scattering.
 
C...Select incoming particle (rho/omega/phi/J/psi for VDM) and mass.
        MINT(103)=MINT(11)
        MINT(104)=MINT(12)
        PMM(1)=VINT(3)
        PMM(2)=VINT(4)
        IF(MINT(101).GE.2.OR.MINT(102).GE.2) THEN
          JJ=ISUB-90
          VRN=PYR(0)*SIGT(0,0,JJ)
          IF(MINT(101).LE.1) THEN
            I1MN=0
            I1MX=0
          ELSE
            I1MN=1
            I1MX=MINT(101)
          ENDIF
          IF(MINT(102).LE.1) THEN
            I2MN=0
            I2MX=0
          ELSE
            I2MN=1
            I2MX=MINT(102)
          ENDIF
          DO 210 I1=I1MN,I1MX
            KFV1=110*I1+3
            DO 200 I2=I2MN,I2MX
              KFV2=110*I2+3
              VRN=VRN-SIGT(I1,I2,JJ)
              IF(VRN.LE.0D0) GOTO 220
  200       CONTINUE
  210     CONTINUE
  220     IF(MINT(101).GE.2) THEN
            MINT(103)=KFV1
            PMM(1)=PYMASS(KFV1)
          ENDIF
          IF(MINT(102).GE.2) THEN
            MINT(104)=KFV2
            PMM(2)=PYMASS(KFV2)
          ENDIF
        ENDIF
        VINT(67)=PMM(1)
        VINT(68)=PMM(2)
 
C...Select mass for GVMD states (rejecting previous assignment).
        Q0S=4D0*PARP(15)**2
        Q1S=4D0*VINT(154)**2
        LOOP3=0
  230   LOOP3=LOOP3+1
        DO 240 JT=1,2
          IF(MINT(106+JT).EQ.3) THEN
            PS=VINT(2+JT)**2
            PMM(JT)=SQRT((Q0S+PS)*(Q1S+PS)/
     &      (Q0S+PYR(0)*(Q1S-Q0S)+PS)-PS)
            IF(MINT(102+JT).GE.333) PMM(JT)=PMM(JT)-
     &      PMAS(PYCOMP(113),1)+PMAS(PYCOMP(MINT(102+JT)),1)
          ENDIF
  240   CONTINUE
        IF(PMM(1)+PMM(2)+PARP(104).GE.VINT(1)) THEN
          IF(LOOP3.LT.100.AND.(MINT(107).EQ.3.OR.MINT(108).EQ.3))
     &    GOTO 230
          GOTO 100
        ENDIF
 
C...Side/sides of diffractive system.
        MINT(17)=0
        MINT(18)=0
        IF(ISUB.EQ.92.OR.ISUB.EQ.94) MINT(17)=1
        IF(ISUB.EQ.93.OR.ISUB.EQ.94) MINT(18)=1
 
C...Find masses of particles and minimal masses of diffractive states.
        DO 250 JT=1,2
          PDIF(JT)=PMM(JT)
          VINT(68+JT)=PDIF(JT)
          IF(MINT(16+JT).EQ.1) PDIF(JT)=PDIF(JT)+PARP(102)
  250   CONTINUE
        SH=VINT(2)
        SQM1=PMM(1)**2
        SQM2=PMM(2)**2
        SQM3=PDIF(1)**2
        SQM4=PDIF(2)**2
        SMRES1=(PMM(1)+PMRC)**2
        SMRES2=(PMM(2)+PMRC)**2
 
C...Find elastic slope and lower limit diffractive slope.
        IHA=MAX(2,IABS(MINT(103))/110)
        IF(IHA.GE.5) IHA=1
        IHB=MAX(2,IABS(MINT(104))/110)
        IF(IHB.GE.5) IHB=1
        IF(ISUB.EQ.91) THEN
          BMN=2D0*BHAD(IHA)+2D0*BHAD(IHB)+4D0*SH**EPS-4.2D0
        ELSEIF(ISUB.EQ.92) THEN
          BMN=MAX(2D0,2D0*BHAD(IHB))
        ELSEIF(ISUB.EQ.93) THEN
          BMN=MAX(2D0,2D0*BHAD(IHA))
        ELSEIF(ISUB.EQ.94) THEN
          BMN=2D0*ALP*4D0
        ENDIF
 
C...Determine maximum possible t range and coefficient of generation.
        SQLA12=(SH-SQM1-SQM2)**2-4D0*SQM1*SQM2
        SQLA34=(SH-SQM3-SQM4)**2-4D0*SQM3*SQM4
        THA=SH-(SQM1+SQM2+SQM3+SQM4)+(SQM1-SQM2)*(SQM3-SQM4)/SH
        THB=SQRT(MAX(0D0,SQLA12))*SQRT(MAX(0D0,SQLA34))/SH
        THC=(SQM3-SQM1)*(SQM4-SQM2)+(SQM1+SQM4-SQM2-SQM3)*
     &  (SQM1*SQM4-SQM2*SQM3)/SH
        THL=-0.5D0*(THA+THB)
        THU=THC/THL
        THRND=EXP(MAX(-50D0,BMN*(THL-THU)))-1D0
 
C...Select diffractive mass/masses according to dm^2/m^2.
        LOOP3=0
  260   LOOP3=LOOP3+1
        DO 270 JT=1,2
          IF(MINT(16+JT).EQ.0) THEN
            PDIF(2+JT)=PDIF(JT)
          ELSE
            PMMIN=PDIF(JT)
            PMMAX=MAX(VINT(2+JT),VINT(1)-PDIF(3-JT))
            PDIF(2+JT)=PMMIN*(PMMAX/PMMIN)**PYR(0)
          ENDIF
  270   CONTINUE
        SQM3=PDIF(3)**2
        SQM4=PDIF(4)**2
 
C..Additional mass factors, including resonance enhancement.
        IF(PDIF(3)+PDIF(4).GE.VINT(1)) THEN
          IF(LOOP3.LT.100) GOTO 260
          GOTO 100
        ENDIF
        IF(ISUB.EQ.92) THEN
          FSD=(1D0-SQM3/SH)*(1D0+CRES*SMRES1/(SMRES1+SQM3))
          IF(FSD.LT.PYR(0)*(1D0+CRES)) GOTO 260
        ELSEIF(ISUB.EQ.93) THEN
          FSD=(1D0-SQM4/SH)*(1D0+CRES*SMRES2/(SMRES2+SQM4))
          IF(FSD.LT.PYR(0)*(1D0+CRES)) GOTO 260
        ELSEIF(ISUB.EQ.94) THEN
          FDD=(1D0-(PDIF(3)+PDIF(4))**2/SH)*(SH*SMP/
     &    (SH*SMP+SQM3*SQM4))*(1D0+CRES*SMRES1/(SMRES1+SQM3))*
     &    (1D0+CRES*SMRES2/(SMRES2+SQM4))
          IF(FDD.LT.PYR(0)*(1D0+CRES)**2) GOTO 260
        ENDIF
 
C...Select t according to exp(Bmn*t) and correct to right slope.
        TH=THU+LOG(1D0+THRND*PYR(0))/BMN
        IF(ISUB.GE.92) THEN
          IF(ISUB.EQ.92) THEN
            BADD=2D0*ALP*LOG(SH/SQM3)
            IF(BHAD(IHB).LT.1D0) BADD=MAX(0D0,BADD+2D0*BHAD(IHB)-2D0)
          ELSEIF(ISUB.EQ.93) THEN
            BADD=2D0*ALP*LOG(SH/SQM4)
            IF(BHAD(IHA).LT.1D0) BADD=MAX(0D0,BADD+2D0*BHAD(IHA)-2D0)
          ELSEIF(ISUB.EQ.94) THEN
            BADD=2D0*ALP*(LOG(EXP(4D0)+SH/(ALP*SQM3*SQM4))-4D0)
          ENDIF
          IF(EXP(MAX(-50D0,BADD*(TH-THU))).LT.PYR(0)) GOTO 260
        ENDIF
 
C...Check whether m^2 and t choices are consistent.
        SQLA34=(SH-SQM3-SQM4)**2-4D0*SQM3*SQM4
        THA=SH-(SQM1+SQM2+SQM3+SQM4)+(SQM1-SQM2)*(SQM3-SQM4)/SH
        THB=SQRT(MAX(0D0,SQLA12))*SQRT(MAX(0D0,SQLA34))/SH
        IF(THB.LE.1D-8) GOTO 260
        THC=(SQM3-SQM1)*(SQM4-SQM2)+(SQM1+SQM4-SQM2-SQM3)*
     &  (SQM1*SQM4-SQM2*SQM3)/SH
        THLM=-0.5D0*(THA+THB)
        THUM=THC/THLM
        IF(TH.LT.THLM.OR.TH.GT.THUM) GOTO 260
 
C...Information to output.
        VINT(21)=1D0
        VINT(22)=0D0
        VINT(23)=MIN(1D0,MAX(-1D0,(THA+2D0*TH)/THB))
        VINT(45)=TH
        VINT(59)=2D0*SQRT(MAX(0D0,-(THC+THA*TH+TH**2)))/THB
        VINT(63)=PDIF(3)**2
        VINT(64)=PDIF(4)**2
        VINT(283)=PMM(1)**2/4D0
        VINT(284)=PMM(2)**2/4D0
 
C...Note: in the following, by In is meant the integral over the
C...quantity multiplying coefficient cn.
C...Choose tau according to h1(tau)/tau, where
C...h1(tau) = c1 + I1/I2*c2*1/tau + I1/I3*c3*1/(tau+tau_R) +
C...I1/I4*c4*tau/((s*tau-m^2)^2+(m*Gamma)^2) +
C...I1/I5*c5*1/(tau+tau_R') +
C...I1/I6*c6*tau/((s*tau-m'^2)^2+(m'*Gamma')^2) +
C...I1/I7*c7*tau/(1.-tau), and
C...c1 + c2 + c3 + c4 + c5 + c6 + c7 = 1.
      ELSEIF(ISTSB.GE.1.AND.ISTSB.LE.5) THEN
        CALL PYKLIM(1)
        IF(MINT(51).NE.0) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(MFAIL.EQ.1) THEN
            MSTI(61)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
        RTAU=PYR(0)
        MTAU=1
        IF(RTAU.GT.COEF(ISUB,1)) MTAU=2
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)) MTAU=3
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)) MTAU=4
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4))
     &  MTAU=5
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)+
     &  COEF(ISUB,5)) MTAU=6
        IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)+COEF(ISUB,4)+
     &  COEF(ISUB,5)+COEF(ISUB,6)) MTAU=7
C...Additional check to handle techni-processes with extra resonance
C....Only modify tau treatment
        IF(ISUB.EQ.194.OR.ISUB.EQ.195.OR.(ISUB.GE.361.AND.ISUB.LE.380))
     &   THEN
          IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)
     &     +COEF(ISUB,4)+COEF(ISUB,5)+COEF(ISUB,6)+COEF(ISUB,7)) MTAU=8
          IF(RTAU.GT.COEF(ISUB,1)+COEF(ISUB,2)+COEF(ISUB,3)
     &     +COEF(ISUB,4)+COEF(ISUB,5)+COEF(ISUB,6)+COEF(ISUB,7)
     &     +COEFX(ISUB,1)) MTAU=9
        ENDIF
        CALL PYKMAP(1,MTAU,PYR(0))
 
C...2 -> 3, 4 processes:
C...Choose tau' according to h4(tau,tau')/tau', where
C...h4(tau,tau') = c1 + I1/I2*c2*(1 - tau/tau')^3/tau' +
C...I1/I3*c3*1/(1 - tau'), and c1 + c2 + c3 = 1.
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
          CALL PYKLIM(4)
          IF(MINT(51).NE.0) THEN
            IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
            IF(MFAIL.EQ.1) THEN
              MSTI(61)=1
              RETURN
            ENDIF
            GOTO 100
          ENDIF
          RTAUP=PYR(0)
          MTAUP=1
          IF(RTAUP.GT.COEF(ISUB,18)) MTAUP=2
          IF(RTAUP.GT.COEF(ISUB,18)+COEF(ISUB,19)) MTAUP=3
          CALL PYKMAP(4,MTAUP,PYR(0))
        ENDIF
 
C...Choose y* according to h2(y*), where
C...h2(y*) = I0/I1*c1*(y*-y*min) + I0/I2*c2*(y*max-y*) +
C...I0/I3*c3*1/cosh(y*) + I0/I4*c4*1/(1-exp(y*-y*max)) +
C...I0/I5*c5*1/(1-exp(-y*-y*min)), I0 = y*max-y*min,
C...and c1 + c2 + c3 + c4 + c5 = 1.
        CALL PYKLIM(2)
        IF(MINT(51).NE.0) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(MFAIL.EQ.1) THEN
            MSTI(61)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
        RYST=PYR(0)
        MYST=1
        IF(RYST.GT.COEF(ISUB,8)) MYST=2
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)) MYST=3
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)+COEF(ISUB,10)) MYST=4
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)+COEF(ISUB,10)+
     &  COEF(ISUB,11)) MYST=5
        CALL PYKMAP(2,MYST,PYR(0))
 
C...2 -> 2 processes:
C...Choose cos(theta-hat) (cth) according to h3(cth), where
C...h3(cth) = c0 + I0/I1*c1*1/(A - cth) + I0/I2*c2*1/(A + cth) +
C...I0/I3*c3*1/(A - cth)^2 + I0/I4*c4*1/(A + cth)^2,
C...A = 1 + 2*(m3*m4/sh)^2 (= 1 for massless products),
C...and c0 + c1 + c2 + c3 + c4 = 1.
        CALL PYKLIM(3)
        IF(MINT(51).NE.0) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(MFAIL.EQ.1) THEN
            MSTI(61)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
        IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN
          RCTH=PYR(0)
          MCTH=1
          IF(RCTH.GT.COEF(ISUB,13)) MCTH=2
          IF(RCTH.GT.COEF(ISUB,13)+COEF(ISUB,14)) MCTH=3
          IF(RCTH.GT.COEF(ISUB,13)+COEF(ISUB,14)+COEF(ISUB,15)) MCTH=4
          IF(RCTH.GT.COEF(ISUB,13)+COEF(ISUB,14)+COEF(ISUB,15)+
     &    COEF(ISUB,16)) MCTH=5
          CALL PYKMAP(3,MCTH,PYR(0))
        ENDIF
 
C...2 -> 3 : select pT1, phi1, pT2, phi2, y3 for 3 outgoing.
        IF(ISTSB.EQ.5) THEN
          CALL PYKMAP(5,0,0D0)
          IF(MINT(51).NE.0) THEN
            IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
            IF(MFAIL.EQ.1) THEN
              MSTI(61)=1
              RETURN
            ENDIF
            GOTO 100
          ENDIF
        ENDIF
 
C...DIS as f + gamma* -> f process: set dummy values.
      ELSEIF(ISTSB.EQ.8) THEN
        VINT(21)=0.9D0
        VINT(22)=0D0
        VINT(23)=0D0
        VINT(47)=0D0
        VINT(48)=0D0
 
C...Low-pT or multiple interactions (first semihard interaction).
      ELSEIF(ISTSB.EQ.9) THEN
        IF(MINT(35).LE.1) CALL PYMULT(3)
        IF(MINT(35).GE.2) CALL PYMIGN(3)
        ISUB=MINT(1)
 
C...Study user-defined process: kinematics plus weight.
      ELSEIF(ISTSB.EQ.11) THEN
        IF(IDWTUP.GT.0.AND.XWGTUP.LT.0D0) CALL
     &  PYERRM(26,'(PYRAND:) Negative XWGTUP for user process')
        MSTI(51)=0
        IF(NUP.LE.0) THEN
          MINT(51)=2
          MSTI(51)=1
          IF(MINT(82).EQ.1) THEN
            NGEN(0,1)=NGEN(0,1)-1
            NGEN(ISUB,1)=NGEN(ISUB,1)-1
          ENDIF
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          RETURN
        ENDIF
 
C...Extract cross section event weight.
        IF(IABS(IDWTUP).EQ.1.OR.IABS(IDWTUP).EQ.4) THEN
          SIGS=1D-9*XWGTUP
        ELSE
          SIGS=1D-9*XSECUP(KFPR(ISUB,1))
        ENDIF
        IF(IABS(IDWTUP).GE.1.AND.IABS(IDWTUP).LE.3) THEN
          VINT(97)=SIGN(1D0,XWGTUP)
        ELSE
          VINT(97)=1D-9*XWGTUP
        ENDIF
 
C...Construct 'trivial' kinematical variables needed.
        KFL1=IDUP(1)
        KFL2=IDUP(2)
        VINT(41)=PUP(4,1)/EBMUP(1)
        VINT(42)=PUP(4,2)/EBMUP(2)
        IF (VINT(41).GT.1.000001.OR.VINT(42).GT.1.000001) THEN
          CALL PYERRM(9,'(PYRAND:) x > 1 in external event '//
     &        '(listing follows):') 
          CALL PYLIST(7)
        ENDIF
        VINT(21)=VINT(41)*VINT(42)
        VINT(22)=0.5D0*LOG(VINT(41)/VINT(42))
        VINT(44)=VINT(21)*VINT(2)
        VINT(43)=SQRT(MAX(0D0,VINT(44)))
        VINT(55)=SCALUP
        IF(SCALUP.LE.0D0) VINT(55)=VINT(43)
        VINT(56)=VINT(55)**2
        VINT(57)=AQEDUP
        VINT(58)=AQCDUP
 
C...Construct other kinematical variables needed (approximately).
        VINT(23)=0D0
        VINT(26)=VINT(21)
        VINT(45)=-0.5D0*VINT(44)
        VINT(46)=-0.5D0*VINT(44)
        VINT(49)=VINT(43)
        VINT(50)=VINT(44)
        VINT(51)=VINT(55)
        VINT(52)=VINT(56)
        VINT(53)=VINT(55)
        VINT(54)=VINT(56)
        VINT(25)=0D0
        VINT(48)=0D0
        IF(ISTUP(1).NE.-1.OR.ISTUP(2).NE.-1) CALL PYERRM(26,
     &  '(PYRAND:) unacceptable ISTUP code for incoming particles')
        DO 280 IUP=3,NUP
          IF(ISTUP(IUP).LT.1.OR.ISTUP(IUP).GT.3) CALL PYERRM(26,
     &    '(PYRAND:) unacceptable ISTUP code for particles')
          IF(ISTUP(IUP).EQ.1) VINT(25)=VINT(25)+2D0*(PUP(5,IUP)**2+
     &    PUP(1,IUP)**2+PUP(2,IUP)**2)/VINT(2)
          IF(ISTUP(IUP).EQ.1) VINT(48)=VINT(48)+0.5D0*(PUP(1,IUP)**2+
     &    PUP(2,IUP)**2)
  280   CONTINUE
        VINT(47)=SQRT(VINT(48))
      ENDIF
 
C...Choose azimuthal angle.
      VINT(24)=0D0
      IF(ISTSB.NE.11) VINT(24)=PARU(2)*PYR(0)
 
C...Check against user cuts on kinematics at parton level.
      MINT(51)=0
      IF((ISUB.LE.90.OR.ISUB.GT.100).AND.ISTSB.LE.10) CALL PYKLIM(0)
      IF(MINT(51).NE.0) THEN
        IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
        IF(MFAIL.EQ.1) THEN
          MSTI(61)=1
          RETURN
        ENDIF
        GOTO 100
      ENDIF
      IF(MINT(82).EQ.1.AND.MSTP(141).GE.1.AND.ISTSB.LE.10) THEN
        MCUT=0
        IF(MSUB(91)+MSUB(92)+MSUB(93)+MSUB(94)+MSUB(95).EQ.0)
     &  CALL PYKCUT(MCUT)
        IF(MCUT.NE.0) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(MFAIL.EQ.1) THEN
            MSTI(61)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
      ENDIF
 
      IF(ISTSB.LE.10) THEN
C...  If internal process, call PYSIGH
        CALL PYSIGH(NCHN,SIGS)
      ELSE
C...  If external process, still have to set MI starting scale 
        IF (MSTP(86).EQ.1) THEN
C...  Limit phase space by xT2 of hard interaction
C...  (gives undercounting of MI when ext proc != dijets)
          XT2GMX = VINT(25)
        ELSE
C...  All accessible phase space allowed
C...  (gives double counting of MI when ext proc = dijets)
          XT2GMX = (1D0-VINT(41))*(1D0-VINT(42))
        ENDIF
        VINT(62)=0.25D0*XT2GMX*VINT(2)
        VINT(61)=SQRT(MAX(0D0,VINT(62)))
      ENDIF
      
      SIGSOR=SIGS
      SIGLPT=SIGT(0,0,5)*VINT(315)*VINT(316)
 
C...Multiply cross section by lepton -> photon flux factor.
      IF(MINT(141).NE.0.OR.MINT(142).NE.0) THEN
        SIGS=WTGAGA*SIGS
        DO 290 ICHN=1,NCHN
          SIGH(ICHN)=WTGAGA*SIGH(ICHN)
  290   CONTINUE
        SIGLPT=WTGAGA*SIGLPT
      ENDIF
 
C...Multiply cross-section by user-defined weights.
      IF(MSTP(173).EQ.1) THEN
        SIGS=PARP(173)*SIGS
        DO 300 ICHN=1,NCHN
          SIGH(ICHN)=PARP(173)*SIGH(ICHN)
  300   CONTINUE
        SIGLPT=PARP(173)*SIGLPT
      ENDIF
      WTXS=1D0
      SIGSWT=SIGS
      VINT(99)=1D0
      VINT(100)=1D0
      IF(MINT(82).EQ.1.AND.MSTP(142).GE.1) THEN
        IF(ISUB.NE.96.AND.MSUB(91)+MSUB(92)+MSUB(93)+MSUB(94)+
     &  MSUB(95).EQ.0) CALL PYEVWT(WTXS)
        SIGSWT=WTXS*SIGS
        VINT(99)=WTXS
        IF(MSTP(142).EQ.1) VINT(100)=1D0/WTXS
      ENDIF
 
C...Calculations for Monte Carlo estimate of all cross-sections.
      IF(MINT(82).EQ.1.AND.ISUB.LE.90.OR.ISUB.GE.96) THEN
        IF(MSTP(142).LE.1) THEN
          XSEC(ISUB,2)=XSEC(ISUB,2)+SIGS
        ELSE
          XSEC(ISUB,2)=XSEC(ISUB,2)+SIGSWT
        ENDIF
      ELSEIF(MINT(82).EQ.1) THEN
        XSEC(ISUB,2)=XSEC(ISUB,2)+SIGS
      ENDIF
      IF((ISUB.EQ.95.OR.ISUB.EQ.96).AND.LOOP2.EQ.1.AND.
     &MINT(82).EQ.1) XSEC(97,2)=XSEC(97,2)+SIGLPT
 
C...Multiple interactions: store results of cross-section calculation.
      IF(MINT(50).EQ.1.AND.MSTP(82).GE.3) THEN
        VINT(153)=SIGSOR
        IF(MINT(35).LE.1) CALL PYMULT(4)
        IF(MINT(35).GE.2) CALL PYMIGN(4)
      ENDIF
 
C...Ratio of actual to maximum cross section.
      IF(ISTSB.NE.11) THEN
        VIOL=SIGSWT/XSEC(ISUB,1)
        IF(ISUB.EQ.96.AND.MSTP(173).EQ.1) VIOL=VIOL/PARP(174)
      ELSEIF(IDWTUP.EQ.1.OR.IDWTUP.EQ.2) THEN
        VIOL=XWGTUP/XMAXUP(KFPR(ISUB,1))
      ELSEIF(IDWTUP.EQ.-1.OR.IDWTUP.EQ.-2) THEN
        VIOL=ABS(XWGTUP)/ABS(XMAXUP(KFPR(ISUB,1)))
      ELSE
        VIOL=1D0
      ENDIF
 
C...Check that weight not negative.
      IF(MSTP(123).LE.0) THEN
        IF(VIOL.LT.-1D-3) THEN
          WRITE(MSTU(11),5000) VIOL,NGEN(0,3)+1
          IF(MSTP(122).GE.1) WRITE(MSTU(11),5100) ISUB,VINT(21),
     &    VINT(22),VINT(23),VINT(26)
          CALL PYSTOP(2)
        ENDIF
      ELSE
        IF(VIOL.LT.MIN(-1D-3,VINT(109))) THEN
          VINT(109)=VIOL
          IF(MSTP(123).LE.2) WRITE(MSTU(11),5200) VIOL,NGEN(0,3)+1
          IF(MSTP(122).GE.1) WRITE(MSTU(11),5100) ISUB,VINT(21),
     &    VINT(22),VINT(23),VINT(26)
        ENDIF
      ENDIF
 
C...Weighting using estimate of maximum of differential cross-section.
      RATND=1D0
      IF(MFAIL.EQ.0.AND.ISUB.NE.95.AND.ISUB.NE.96) THEN
        IF(VIOL.LT.PYR(0)) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(ISUB.GE.91.AND.ISUB.LE.94) ISUB=0
          GOTO 100
        ENDIF
      ELSEIF(MFAIL.EQ.0) THEN
        RATND=SIGLPT/XSEC(95,1)
        VIOL=VIOL/RATND
        IF(LOOP2.EQ.1.AND.RATND.LT.PYR(0)) THEN
          IF(VIOL.GT.PYR(0).AND.MINT(82).EQ.1.AND.MSUB(95).EQ.1.AND.
     &    (ISUB.LE.90.OR.ISUB.GE.95)) NGEN(95,1)=NGEN(95,1)+MINT(143)
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          ISUB=0
          GOTO 100
        ENDIF
        IF(VIOL.LT.PYR(0)) THEN
          GOTO 140
        ENDIF
      ELSEIF(ISUB.NE.95.AND.ISUB.NE.96) THEN
        IF(VIOL.LT.PYR(0)) THEN
          MSTI(61)=1
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          RETURN
        ENDIF
      ELSE
        RATND=SIGLPT/XSEC(95,1)
        IF(LOOP.EQ.1.AND.RATND.LT.PYR(0)) THEN
          MSTI(61)=1
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          RETURN
        ENDIF
        VIOL=VIOL/RATND
        IF(VIOL.LT.PYR(0)) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          GOTO 100
        ENDIF
      ENDIF
 
C...Check for possible violation of estimated maximum of differential
C...cross-section used in weighting.
      IF(MSTP(123).LE.0) THEN
        IF(VIOL.GT.1D0) THEN
          WRITE(MSTU(11),5300) VIOL,NGEN(0,3)+1
          IF(MSTP(122).GE.2) WRITE(MSTU(11),5100) ISUB,VINT(21),
     &    VINT(22),VINT(23),VINT(26)
          CALL PYSTOP(2)
        ENDIF
      ELSEIF(MSTP(123).EQ.1) THEN
        IF(VIOL.GT.VINT(108)) THEN
          VINT(108)=VIOL
          IF(VIOL.GT.1.0001D0) THEN
            MINT(10)=1
            WRITE(MSTU(11),5400) VIOL,NGEN(0,3)+1
            IF(MSTP(122).GE.2) WRITE(MSTU(11),5100) ISUB,VINT(21),
     &      VINT(22),VINT(23),VINT(26)
          ENDIF
        ENDIF
      ELSEIF(VIOL.GT.VINT(108)) THEN
        VINT(108)=VIOL
        IF(VIOL.GT.1D0) THEN
          MINT(10)=1
          IF(MSTP(123).EQ.2) WRITE(MSTU(11),5400) VIOL,NGEN(0,3)+1
          IF(ISTSB.EQ.11.AND.(IABS(IDWTUP).EQ.1.OR.IABS(IDWTUP).EQ.2))
     &    THEN
            XMAXUP(KFPR(ISUB,1))=VIOL*XMAXUP(KFPR(ISUB,1))
            IF(KFPR(ISUB,1).LE.9) THEN
              IF(MSTP(123).EQ.2) WRITE(MSTU(11),5800) KFPR(ISUB,1),
     &        XMAXUP(KFPR(ISUB,1))
            ELSEIF(KFPR(ISUB,1).LE.99) THEN
              IF(MSTP(123).EQ.2) WRITE(MSTU(11),5900) KFPR(ISUB,1),
     &        XMAXUP(KFPR(ISUB,1))
            ELSE
              IF(MSTP(123).EQ.2) WRITE(MSTU(11),6000) KFPR(ISUB,1),
     &        XMAXUP(KFPR(ISUB,1))
            ENDIF
          ENDIF
          IF(ISTSB.NE.11.OR.IABS(IDWTUP).EQ.1) THEN
            XDIF=XSEC(ISUB,1)*(VIOL-1D0)
            XSEC(ISUB,1)=XSEC(ISUB,1)+XDIF
            IF(MSUB(ISUB).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GT.96))
     &      XSEC(0,1)=XSEC(0,1)+XDIF
            IF(MSTP(122).GE.2) WRITE(MSTU(11),5100) ISUB,VINT(21),
     &      VINT(22),VINT(23),VINT(26)
            IF(ISUB.LE.9) THEN
              IF(MSTP(123).EQ.2) WRITE(MSTU(11),5500) ISUB,XSEC(ISUB,1)
            ELSEIF(ISUB.LE.99) THEN
              IF(MSTP(123).EQ.2) WRITE(MSTU(11),5600) ISUB,XSEC(ISUB,1)
            ELSE
              IF(MSTP(123).EQ.2) WRITE(MSTU(11),5700) ISUB,XSEC(ISUB,1)
            ENDIF
          ENDIF
          VINT(108)=1D0
        ENDIF
      ENDIF
 
C...Multiple interactions: choose impact parameter (if not already done).
      IF(MINT(39).EQ.0) VINT(148)=1D0
      IF(MINT(50).EQ.1.AND.(ISUB.LE.90.OR.ISUB.GE.96).AND.
     &MSTP(82).GE.3) THEN
        IF(MINT(35).LE.1) CALL PYMULT(5)
        IF(MINT(35).GE.2) CALL PYMIGN(5)
        IF(VINT(150).LT.PYR(0)) THEN
          IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
          IF(MFAIL.EQ.1) THEN
            MSTI(61)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
      ENDIF
      IF(MINT(82).EQ.1) NGEN(0,2)=NGEN(0,2)+1
      IF(MINT(82).EQ.1.AND.MSUB(95).EQ.1) THEN
        IF(ISUB.LE.90.OR.ISUB.GE.95) NGEN(95,1)=NGEN(95,1)+MINT(143)
        IF(ISUB.LE.90.OR.ISUB.GE.96) NGEN(96,2)=NGEN(96,2)+1
      ENDIF
      IF(ISUB.LE.90.OR.ISUB.GE.96) MINT(31)=MINT(31)+1
 
C...Choose flavour of reacting partons (and subprocess).
      IF(ISTSB.GE.11) GOTO 320
      RSIGS=SIGS*PYR(0)
      QT2=VINT(48)
      RQQBAR=PARP(87)*(1D0-(QT2/(QT2+(PARP(88)*PARP(82)*
     &(VINT(1)/PARP(89))**PARP(90))**2))**2)
      IF(ISUB.NE.95.AND.(ISUB.NE.96.OR.MSTP(82).LE.1.OR.
     &PYR(0).GT.RQQBAR)) THEN
        DO 310 ICHN=1,NCHN
          KFL1=ISIG(ICHN,1)
          KFL2=ISIG(ICHN,2)
          MINT(2)=ISIG(ICHN,3)
          RSIGS=RSIGS-SIGH(ICHN)
          IF(RSIGS.LE.0D0) GOTO 320
  310   CONTINUE
 
C...Multiple interactions: choose qqbar preferentially at small pT.
      ELSEIF(ISUB.EQ.96) THEN
        MINT(105)=MINT(103)
        MINT(109)=MINT(107)
        CALL PYSPLI(MINT(11),21,KFL1,KFLDUM)
        MINT(105)=MINT(104)
        MINT(109)=MINT(108)
        CALL PYSPLI(MINT(12),21,KFL2,KFLDUM)
        MINT(1)=11
        MINT(2)=1
        IF(KFL1.EQ.KFL2.AND.PYR(0).LT.0.5D0) MINT(2)=2
 
C...Low-pT: choose string drawing configuration.
      ELSE
        KFL1=21
        KFL2=21
        RSIGS=6D0*PYR(0)
        MINT(2)=1
        IF(RSIGS.GT.1D0) MINT(2)=2
        IF(RSIGS.GT.2D0) MINT(2)=3
      ENDIF
 
C...Reassign QCD process. Partons before initial state radiation.
  320 IF(MINT(2).GT.10) THEN
        MINT(1)=MINT(2)/10
        MINT(2)=MOD(MINT(2),10)
      ENDIF
      IF(MINT(82).EQ.1.AND.MSTP(111).GE.0) NGEN(MINT(1),2)=
     &NGEN(MINT(1),2)+1
      MINT(15)=KFL1
      MINT(16)=KFL2
      MINT(13)=MINT(15)
      MINT(14)=MINT(16)
      VINT(141)=VINT(41)
      VINT(142)=VINT(42)
      VINT(151)=0D0
      VINT(152)=0D0
 
C...Calculate x value of photon for parton inside photon inside e.
      DO 350 JT=1,2
        MINT(18+JT)=0
        VINT(154+JT)=0D0
        MSPLI=0
        IF(JT.EQ.1.AND.MINT(43).LE.2) MSPLI=1
        IF(JT.EQ.2.AND.MOD(MINT(43),2).EQ.1) MSPLI=1
        IF(IABS(MINT(14+JT)).LE.8.OR.MINT(14+JT).EQ.21) MSPLI=MSPLI+1
        IF(MSPLI.EQ.2) THEN
          KFLH=MINT(14+JT)
          XHRD=VINT(140+JT)
          Q2HRD=VINT(54)
          MINT(105)=MINT(102+JT)
          MINT(109)=MINT(106+JT)
          VINT(120)=VINT(2+JT)
          IF(MSTP(57).LE.1) THEN
            CALL PYPDFU(22,XHRD,Q2HRD,XPQ)
          ELSE
            CALL PYPDFL(22,XHRD,Q2HRD,XPQ)
          ENDIF
          WTMX=4D0*XPQ(KFLH)
          IF(MSTP(13).EQ.2) THEN
            Q2PMS=Q2HRD/PMAS(11,1)**2
            WTMX=WTMX*LOG(MAX(2D0,Q2PMS*(1D0-XHRD)/XHRD**2))
          ENDIF
  330     XE=XHRD**PYR(0)
          XG=MIN(1D0-1D-10,XHRD/XE)
          IF(MSTP(57).LE.1) THEN
            CALL PYPDFU(22,XG,Q2HRD,XPQ)
          ELSE
            CALL PYPDFL(22,XG,Q2HRD,XPQ)
          ENDIF
          WT=(1D0+(1D0-XE)**2)*XPQ(KFLH)
          IF(MSTP(13).EQ.2) WT=WT*LOG(MAX(2D0,Q2PMS*(1D0-XE)/XE**2))
          IF(WT.LT.PYR(0)*WTMX) GOTO 330
          MINT(18+JT)=1
          VINT(154+JT)=XE
          DO 340 KFLS=-25,25
            XSFX(JT,KFLS)=XPQ(KFLS)
  340     CONTINUE
        ENDIF
  350 CONTINUE
 
C...Pick scale where photon is resolved.
      Q0S=PARP(15)**2
      Q1S=VINT(154)**2
      VINT(283)=0D0
      IF(MINT(107).EQ.3) THEN
        IF(MSTP(66).EQ.1) THEN
          VINT(283)=Q0S*(VINT(54)/Q0S)**PYR(0)
        ELSEIF(MSTP(66).EQ.2) THEN
          PS=VINT(3)**2
          Q2EFF=VINT(54)*((Q0S+PS)/(VINT(54)+PS))*
     &    EXP(PS*(VINT(54)-Q0S)/((VINT(54)+PS)*(Q0S+PS)))
          Q2INT=SQRT(Q0S*Q2EFF)
          VINT(283)=Q2INT*(VINT(54)/Q2INT)**PYR(0)
        ELSEIF(MSTP(66).EQ.3) THEN
          VINT(283)=Q0S*(Q1S/Q0S)**PYR(0)
        ELSEIF(MSTP(66).GE.4) THEN
          PS=0.25D0*VINT(3)**2
          VINT(283)=(Q0S+PS)*(Q1S+PS)/
     &    (Q0S+PYR(0)*(Q1S-Q0S)+PS)-PS
        ENDIF
      ENDIF
      VINT(284)=0D0
      IF(MINT(108).EQ.3) THEN
        IF(MSTP(66).EQ.1) THEN
          VINT(284)=Q0S*(VINT(54)/Q0S)**PYR(0)
        ELSEIF(MSTP(66).EQ.2) THEN
          PS=VINT(4)**2
          Q2EFF=VINT(54)*((Q0S+PS)/(VINT(54)+PS))*
     &    EXP(PS*(VINT(54)-Q0S)/((VINT(54)+PS)*(Q0S+PS)))
          Q2INT=SQRT(Q0S*Q2EFF)
          VINT(284)=Q2INT*(VINT(54)/Q2INT)**PYR(0)
        ELSEIF(MSTP(66).EQ.3) THEN
          VINT(284)=Q0S*(Q1S/Q0S)**PYR(0)
        ELSEIF(MSTP(66).GE.4) THEN
          PS=0.25D0*VINT(4)**2
          VINT(284)=(Q0S+PS)*(Q1S+PS)/
     &    (Q0S+PYR(0)*(Q1S-Q0S)+PS)-PS
        ENDIF
      ENDIF
      IF(MINT(121).GT.1) CALL PYSAVE(2,IGA)
 
C...Format statements for differential cross-section maximum violations.
 5000 FORMAT(/1X,'Error: negative cross-section fraction',1P,D11.3,1X,
     &'in event',1X,I7,'D0'/1X,'Execution stopped!')
 5100 FORMAT(1X,'ISUB = ',I3,'; Point of violation:'/1X,'tau =',1P,
     &D11.3,', y* =',D11.3,', cthe = ',0P,F11.7,', tau'' =',1P,D11.3)
 5200 FORMAT(/1X,'Warning: negative cross-section fraction',1P,D11.3,1X,
     &'in event',1X,I7)
 5300 FORMAT(/1X,'Error: maximum violated by',1P,D11.3,1X,
     &'in event',1X,I7,'D0'/1X,'Execution stopped!')
 5400 FORMAT(/1X,'Advisory warning: maximum violated by',1P,D11.3,1X,
     &'in event',1X,I7)
 5500 FORMAT(1X,'XSEC(',I1,',1) increased to',1P,D11.3)
 5600 FORMAT(1X,'XSEC(',I2,',1) increased to',1P,D11.3)
 5700 FORMAT(1X,'XSEC(',I3,',1) increased to',1P,D11.3)
 5800 FORMAT(1X,'XMAXUP(',I1,') increased to',1P,D11.3)
 5900 FORMAT(1X,'XMAXUP(',I2,') increased to',1P,D11.3)
 6000 FORMAT(1X,'XMAXUP(',I3,') increased to',1P,D11.3)

      RETURN
      END
