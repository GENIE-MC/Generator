 
C*********************************************************************
 
C...PYINBM
C...Identifies the two incoming particles and the choice of frame.
 
       SUBROUTINE PYINBM(CHFRAM,CHBEAM,CHTARG,WIN)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
 
C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/
 
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYSUBS/,/PYPARS/,/PYINT1/
 
C...Local arrays, character variables and data.
      CHARACTER CHFRAM*12,CHBEAM*12,CHTARG*12,CHCOM(3)*12,CHALP(2)*26,
     &CHIDNT(3)*12,CHTEMP*12,CHCDE(39)*12,CHINIT*76,CHNAME*16
      DIMENSION LEN(3),KCDE(39),PM(2)
      DATA CHALP/'abcdefghijklmnopqrstuvwxyz',
     &'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA CHCDE/    'e-          ','e+          ','nu_e        ',
     &'nu_ebar     ','mu-         ','mu+         ','nu_mu       ',
     &'nu_mubar    ','tau-        ','tau+        ','nu_tau      ',
     &'nu_taubar   ','pi+         ','pi-         ','n0          ',
     &'nbar0       ','p+          ','pbar-       ','gamma       ',
     &'lambda0     ','sigma-      ','sigma0      ','sigma+      ',
     &'xi-         ','xi0         ','omega-      ','pi0         ',
     &'reggeon     ','pomeron     ','gamma/e-    ','gamma/e+    ',
     &'gamma/mu-   ','gamma/mu+   ','gamma/tau-  ','gamma/tau+  ',
     &'k+          ','k-          ','ks0         ','kl0         '/
      DATA KCDE/11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,
     &211,-211,2112,-2112,2212,-2212,22,3122,3112,3212,3222,
     &3312,3322,3334,111,110,990,6*22,321,-321,310,130/
 
C...Store initial energy. Default frame.
      VINT(290)=WIN
      MINT(111)=0
 
C...Special user process initialization; convert to normal input.
      IF(CHFRAM(1:1).EQ.'u'.OR.CHFRAM(1:1).EQ.'U') THEN
        MINT(111)=11
        IF(PDFGUP(1).EQ.-9.OR.PDFGUP(2).EQ.-9) MINT(111)=12
        CALL PYNAME(IDBMUP(1),CHNAME)
        CHBEAM=CHNAME(1:12)
        CALL PYNAME(IDBMUP(2),CHNAME)
        CHTARG=CHNAME(1:12)
      ENDIF
 
C...Convert character variables to lowercase and find their length.
      CHCOM(1)=CHFRAM
      CHCOM(2)=CHBEAM
      CHCOM(3)=CHTARG
      DO 130 I=1,3
        LEN(I)=12
        DO 110 LL=12,1,-1
          IF(LEN(I).EQ.LL.AND.CHCOM(I)(LL:LL).EQ.' ') LEN(I)=LL-1
          DO 100 LA=1,26
            IF(CHCOM(I)(LL:LL).EQ.CHALP(2)(LA:LA)) CHCOM(I)(LL:LL)=
     &      CHALP(1)(LA:LA)
  100     CONTINUE
  110   CONTINUE
        CHIDNT(I)=CHCOM(I)
 
C...Fix up bar, underscore and charge in particle name (if needed).
        DO 120 LL=1,10
          IF(CHIDNT(I)(LL:LL).EQ.'~') THEN
            CHTEMP=CHIDNT(I)
            CHIDNT(I)=CHTEMP(1:LL-1)//'bar'//CHTEMP(LL+1:10)//'  '
          ENDIF
  120   CONTINUE
        IF(CHIDNT(I)(1:2).EQ.'nu'.AND.CHIDNT(I)(3:3).NE.'_') THEN
          CHTEMP=CHIDNT(I)
          CHIDNT(I)='nu_'//CHTEMP(3:7)
        ELSEIF(CHIDNT(I)(1:2).EQ.'n ') THEN
          CHIDNT(I)(1:3)='n0 '
        ELSEIF(CHIDNT(I)(1:4).EQ.'nbar') THEN
          CHIDNT(I)(1:5)='nbar0'
        ELSEIF(CHIDNT(I)(1:2).EQ.'p ') THEN
          CHIDNT(I)(1:3)='p+ '
        ELSEIF(CHIDNT(I)(1:4).EQ.'pbar'.OR.
     &    CHIDNT(I)(1:2).EQ.'p-') THEN
          CHIDNT(I)(1:5)='pbar-'
        ELSEIF(CHIDNT(I)(1:6).EQ.'lambda') THEN
          CHIDNT(I)(7:7)='0'
        ELSEIF(CHIDNT(I)(1:3).EQ.'reg') THEN
          CHIDNT(I)(1:7)='reggeon'
        ELSEIF(CHIDNT(I)(1:3).EQ.'pom') THEN
          CHIDNT(I)(1:7)='pomeron'
        ENDIF
  130 CONTINUE
 
C...Identify free initialization.
      IF(CHCOM(1)(1:2).EQ.'no') THEN
        MINT(65)=1
        RETURN
      ENDIF
 
C...Identify incoming beam and target particles.
      DO 160 I=1,2
        DO 140 J=1,39
          IF(CHIDNT(I+1).EQ.CHCDE(J)) MINT(10+I)=KCDE(J)
  140   CONTINUE
        PM(I)=PYMASS(MINT(10+I))
        VINT(2+I)=PM(I)
        MINT(140+I)=0
        IF(MINT(10+I).EQ.22.AND.CHIDNT(I+1)(6:6).EQ.'/') THEN
          CHTEMP=CHIDNT(I+1)(7:12)//' '
          DO 150 J=1,12
            IF(CHTEMP.EQ.CHCDE(J)) MINT(140+I)=KCDE(J)
  150	  CONTINUE
          PM(I)=PYMASS(MINT(140+I))
          VINT(302+I)=PM(I)
        ENDIF
  160 CONTINUE
      IF(MINT(11).EQ.0) WRITE(MSTU(11),5000) CHBEAM(1:LEN(2))
      IF(MINT(12).EQ.0) WRITE(MSTU(11),5100) CHTARG(1:LEN(3))
      IF(MINT(11).EQ.0.OR.MINT(12).EQ.0) CALL PYSTOP(7)
 
C...Identify choice of frame and input energies.
      CHINIT=' '
 
C...Events defined in the CM frame.
      IF(CHCOM(1)(1:2).EQ.'cm') THEN
        MINT(111)=1
        S=WIN**2
        IF(MSTP(122).GE.1) THEN
          IF(CHCOM(2)(1:1).NE.'e') THEN
            LOFFS=(31-(LEN(2)+LEN(3)))/2
            CHINIT(LOFFS+1:76)='PYTHIA will be initialized for a '//
     &      CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//
     &      ' collider'//' '
          ELSE
            LOFFS=(30-(LEN(2)+LEN(3)))/2
            CHINIT(LOFFS+1:76)='PYTHIA will be initialized for an '//
     &      CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//
     &      ' collider'//' '
          ENDIF
          WRITE(MSTU(11),5200) CHINIT
          WRITE(MSTU(11),5300) WIN
        ENDIF
 
C...Events defined in fixed target frame.
      ELSEIF(CHCOM(1)(1:3).EQ.'fix') THEN
        MINT(111)=2
        S=PM(1)**2+PM(2)**2+2D0*PM(2)*SQRT(PM(1)**2+WIN**2)
        IF(MSTP(122).GE.1) THEN
          LOFFS=(29-(LEN(2)+LEN(3)))/2
          CHINIT(LOFFS+1:76)='PYTHIA will be initialized for '//
     &    CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//
     &    ' fixed target'//' '
          WRITE(MSTU(11),5200) CHINIT
          WRITE(MSTU(11),5400) WIN
          WRITE(MSTU(11),5500) SQRT(S)
        ENDIF
 
C...Frame defined by user three-vectors.
      ELSEIF(CHCOM(1)(1:1).EQ.'3') THEN
        MINT(111)=3
        P(1,5)=PM(1)
        P(2,5)=PM(2)
        P(1,4)=SQRT(P(1,1)**2+P(1,2)**2+P(1,3)**2+P(1,5)**2)
        P(2,4)=SQRT(P(2,1)**2+P(2,2)**2+P(2,3)**2+P(2,5)**2)
        S=(P(1,4)+P(2,4))**2-(P(1,1)+P(2,1))**2-(P(1,2)+P(2,2))**2-
     &  (P(1,3)+P(2,3))**2
        IF(MSTP(122).GE.1) THEN
          LOFFS=(22-(LEN(2)+LEN(3)))/2
          CHINIT(LOFFS+1:76)='PYTHIA will be initialized for '//
     &    CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//
     &    ' user configuration'//' '
          WRITE(MSTU(11),5200) CHINIT
          WRITE(MSTU(11),5600)
          WRITE(MSTU(11),5700) CHCOM(2),P(1,1),P(1,2),P(1,3),P(1,4)
          WRITE(MSTU(11),5700) CHCOM(3),P(2,1),P(2,2),P(2,3),P(2,4)
          WRITE(MSTU(11),5500) SQRT(MAX(0D0,S))
        ENDIF
 
C...Frame defined by user four-vectors.
      ELSEIF(CHCOM(1)(1:1).EQ.'4') THEN
        MINT(111)=4
        PMS1=P(1,4)**2-P(1,1)**2-P(1,2)**2-P(1,3)**2
        P(1,5)=SIGN(SQRT(ABS(PMS1)),PMS1)
        PMS2=P(2,4)**2-P(2,1)**2-P(2,2)**2-P(2,3)**2
        P(2,5)=SIGN(SQRT(ABS(PMS2)),PMS2)
        S=(P(1,4)+P(2,4))**2-(P(1,1)+P(2,1))**2-(P(1,2)+P(2,2))**2-
     &  (P(1,3)+P(2,3))**2
        IF(MSTP(122).GE.1) THEN
          LOFFS=(22-(LEN(2)+LEN(3)))/2
          CHINIT(LOFFS+1:76)='PYTHIA will be initialized for '//
     &    CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//
     &    ' user configuration'//' '
          WRITE(MSTU(11),5200) CHINIT
          WRITE(MSTU(11),5600)
          WRITE(MSTU(11),5700) CHCOM(2),P(1,1),P(1,2),P(1,3),P(1,4)
          WRITE(MSTU(11),5700) CHCOM(3),P(2,1),P(2,2),P(2,3),P(2,4)
          WRITE(MSTU(11),5500) SQRT(MAX(0D0,S))
        ENDIF
 
C...Frame defined by user five-vectors.
      ELSEIF(CHCOM(1)(1:1).EQ.'5') THEN
        MINT(111)=5
        S=(P(1,4)+P(2,4))**2-(P(1,1)+P(2,1))**2-(P(1,2)+P(2,2))**2-
     &  (P(1,3)+P(2,3))**2
        IF(MSTP(122).GE.1) THEN
          LOFFS=(22-(LEN(2)+LEN(3)))/2
          CHINIT(LOFFS+1:76)='PYTHIA will be initialized for '//
     &    CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//
     &    ' user configuration'//' '
          WRITE(MSTU(11),5200) CHINIT
          WRITE(MSTU(11),5600)
          WRITE(MSTU(11),5700) CHCOM(2),P(1,1),P(1,2),P(1,3),P(1,4)
          WRITE(MSTU(11),5700) CHCOM(3),P(2,1),P(2,2),P(2,3),P(2,4)
          WRITE(MSTU(11),5500) SQRT(MAX(0D0,S))
        ENDIF
 
C...Frame defined by HEPRUP common block.
      ELSEIF(MINT(111).GE.11) THEN
        S=(EBMUP(1)+EBMUP(2))**2-(SQRT(MAX(0D0,EBMUP(1)**2-PM(1)**2))-
     &  SQRT(MAX(0D0,EBMUP(2)**2-PM(2)**2)))**2
        IF(MSTP(122).GE.1) THEN
          LOFFS=(22-(LEN(2)+LEN(3)))/2
          CHINIT(LOFFS+1:76)='PYTHIA will be initialized for '//
     &    CHCOM(2)(1:LEN(2))//' on '//CHCOM(3)(1:LEN(3))//
     &    ' user configuration'//' '
          WRITE(MSTU(11),5200) CHINIT
          WRITE(MSTU(11),6000) EBMUP(1),EBMUP(2)
          WRITE(MSTU(11),5500) SQRT(MAX(0D0,S))
        ENDIF
 
C...Unknown frame. Error for too low CM energy.
      ELSE
        WRITE(MSTU(11),5800) CHFRAM(1:LEN(1))
        CALL PYSTOP(7)
      ENDIF
      IF(S.LT.PARP(2)**2) THEN
        WRITE(MSTU(11),5900) SQRT(S)
        CALL PYSTOP(7)
      ENDIF
 
C...Formats for initialization and error information.
 5000 FORMAT(1X,'Error: unrecognized beam particle ''',A,'''D0'/
     &1X,'Execution stopped!')
 5100 FORMAT(1X,'Error: unrecognized target particle ''',A,'''D0'/
     &1X,'Execution stopped!')
 5200 FORMAT(/1X,78('=')/1X,'I',76X,'I'/1X,'I',A76,'I')
 5300 FORMAT(1X,'I',18X,'at',1X,F10.3,1X,'GeV center-of-mass energy',
     &19X,'I'/1X,'I',76X,'I'/1X,78('='))
 5400 FORMAT(1X,'I',22X,'at',1X,F10.3,1X,'GeV/c lab-momentum',22X,'I')
 5500 FORMAT(1X,'I',76X,'I'/1X,'I',11X,'corresponding to',1X,F10.3,1X,
     &'GeV center-of-mass energy',12X,'I'/1X,'I',76X,'I'/1X,78('='))
 5600 FORMAT(1X,'I',76X,'I'/1X,'I',18X,'px (GeV/c)',3X,'py (GeV/c)',3X,
     &'pz (GeV/c)',6X,'E (GeV)',9X,'I')
 5700 FORMAT(1X,'I',8X,A8,4(2X,F10.3,1X),8X,'I')
 5800 FORMAT(1X,'Error: unrecognized coordinate frame ''',A,'''D0'/
     &1X,'Execution stopped!')
 5900 FORMAT(1X,'Error: too low CM energy,',F8.3,' GeV for event ',
     &'generation.'/1X,'Execution stopped!')
 6000 FORMAT(1X,'I',12X,'with',1X,F10.3,1X,'GeV on',1X,F10.3,1X,
     &'GeV beam energies',13X,'I')
 
      RETURN
      END
