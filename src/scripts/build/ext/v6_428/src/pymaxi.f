 
C*********************************************************************
 
C...PYMAXI
C...Finds optimal set of coefficients for kinematical variable selection
C...and the maximum of the part of the differential cross-section used
C...in the event weighting.
 
      SUBROUTINE PYMAXI
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
 
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
      COMMON/PYINT6/PROC(0:500)
      CHARACTER PROC*28
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      COMMON/PYTCSM/ITCM(0:99),RTCM(0:99)
      COMMON/PYTCCO/COEFX(194:380,2)
      COMMON/TCPARA/IRES,JRES,XMAS(3),XWID(3),YMAS(2),YWID(2)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,/PYINT1/,
     &/PYINT2/,/PYINT3/,/PYINT4/,/PYINT5/,/PYINT6/,/PYINT7/,/PYTCCO/,
     &/PYTCSM/,/TCPARA/
C...Local arrays, character variables and data.
      LOGICAL IOK
      CHARACTER CVAR(4)*4
      DIMENSION NPTS(4),MVARPT(500,4),VINTPT(500,30),SIGSPT(500),
     &NAREL(9),WTREL(9),WTMAT(9,9),WTRELN(9),COEFU(9),COEFO(9),
     &IACCMX(4),SIGSMX(4),SIGSSM(3),PMMN(2),WTRSAV(9),TEMPC(9),
     &IQ(9),IP(9)
      DATA CVAR/'tau ','tau''','y*  ','cth '/
      DATA SIGSSM/3*0D0/
 
C...Initial values and loop over subprocesses.
      NPOSI=0
      VINT(143)=1D0
      VINT(144)=1D0
      XSEC(0,1)=0D0
      ITECH=0
      DO 460 ISUB=1,500
        MINT(1)=ISUB
        MINT(51)=0
 
C...Find maximum weight factors for photon flux.
        IF(MSUB(ISUB).EQ.1.OR.(ISUB.GE.91.AND.ISUB.LE.100)) THEN
          IF(MINT(141).NE.0.OR.MINT(142).NE.0) CALL PYGAGA(2,WTGAGA)
        ENDIF
 
C...Select subprocess to study: skip cases not applicable.
        IF(ISET(ISUB).EQ.11) THEN
          IF(MSUB(ISUB).NE.1) GOTO 460
C...User process intialization: cross section model dependent.
          IF(IABS(IDWTUP).EQ.1) THEN
            IF(IDWTUP.GT.0.AND.XMAXUP(KFPR(ISUB,1)).LT.0D0) CALL
     &      PYERRM(26,'(PYMAXI:) Negative XMAXUP for user process')
            XSEC(ISUB,1)=1.00000001D-9*ABS(XMAXUP(KFPR(ISUB,1)))
          ELSE
            IF((IDWTUP.EQ.2.OR.IDWTUP.EQ.3).AND.
     &      XSECUP(KFPR(ISUB,1)).LT.0D0) CALL
     &      PYERRM(26,'(PYMAXI:) Negative XSECUP for user process')
            IF(IDWTUP.EQ.2.AND.XMAXUP(KFPR(ISUB,1)).LT.0D0) CALL
     &      PYERRM(26,'(PYMAXI:) Negative XMAXUP for user process')
            XSEC(ISUB,1)=1.00000001D-9*ABS(XSECUP(KFPR(ISUB,1)))
          ENDIF
          IF(MINT(141).NE.0.OR.MINT(142).NE.0) XSEC(ISUB,1)=
     &    WTGAGA*XSEC(ISUB,1)
          NPOSI=NPOSI+1
          GOTO 450
        ELSEIF(ISUB.GE.91.AND.ISUB.LE.95) THEN
          CALL PYSIGH(NCHN,SIGS)
          XSEC(ISUB,1)=SIGS
          IF(MINT(141).NE.0.OR.MINT(142).NE.0) XSEC(ISUB,1)=
     &    WTGAGA*XSEC(ISUB,1)
          IF(MSUB(ISUB).NE.1) GOTO 460
          NPOSI=NPOSI+1
          GOTO 450
        ELSEIF(ISUB.EQ.99.AND.MSUB(ISUB).EQ.1) THEN
          CALL PYSIGH(NCHN,SIGS)
          XSEC(ISUB,1)=SIGS
          IF(MINT(141).NE.0.OR.MINT(142).NE.0) XSEC(ISUB,1)=
     &    WTGAGA*XSEC(ISUB,1)
          IF(XSEC(ISUB,1).EQ.0D0) THEN
            MSUB(ISUB)=0
          ELSE
            NPOSI=NPOSI+1
          ENDIF
          GOTO 450
        ELSEIF(ISUB.EQ.96) THEN
          IF(MINT(50).EQ.0) GOTO 460
          IF(MSUB(95).NE.1.AND.MOD(MSTP(81),10).LE.0.AND.MSTP(131).LE.0)
     &    GOTO 460
          IF(MINT(49).EQ.0.AND.MSTP(131).EQ.0) GOTO 460
        ELSEIF(ISUB.EQ.11.OR.ISUB.EQ.12.OR.ISUB.EQ.13.OR.ISUB.EQ.28.OR.
     &    ISUB.EQ.53.OR.ISUB.EQ.68) THEN
          IF(MSUB(ISUB).NE.1.OR.MSUB(95).EQ.1) GOTO 460
        ELSEIF(ISUB.GE.381.AND.ISUB.LE.386) THEN
          IF(MSUB(ISUB).NE.1.OR.MSUB(95).EQ.1) GOTO 460
        ELSE
          IF(MSUB(ISUB).NE.1) GOTO 460
        ENDIF
        ISTSB=ISET(ISUB)
        IF(ISUB.EQ.96) ISTSB=2
        IF(MSTP(122).GE.2) WRITE(MSTU(11),5000) ISUB
        MWTXS=0
        IF(MSTP(142).GE.1.AND.ISUB.NE.96.AND.MSUB(91)+MSUB(92)+MSUB(93)+
     &  MSUB(94)+MSUB(95).EQ.0) MWTXS=1
 
C...Find resonances (explicit or implicit in cross-section).
        MINT(72)=0
        KFR1=0
        IF(ISTSB.EQ.1.OR.ISTSB.EQ.3.OR.ISTSB.EQ.5) THEN
          KFR1=KFPR(ISUB,1)
        ELSEIF(ISUB.EQ.24.OR.ISUB.EQ.25.OR.ISUB.EQ.110.OR.ISUB.EQ.165
     &    .OR.ISUB.EQ.171.OR.ISUB.EQ.176) THEN
          KFR1=23
        ELSEIF(ISUB.EQ.23.OR.ISUB.EQ.26.OR.ISUB.EQ.166.OR.ISUB.EQ.172
     &    .OR.ISUB.EQ.177) THEN
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
     &    CKMX.LT.PMAS(KCR1,1)-20D0*PMAS(KCR1,2)) KFR1=0
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
     $  (ISUB.GE.361.AND.ISUB.LE.380))
     $  THEN
          KFR2=23
          IF(ISUB.EQ.141) THEN
            KCR2=PYCOMP(KFR2)
            IF(CKIN(1).GT.PMAS(KCR2,1)+20D0*PMAS(KCR2,2).OR.
     &       CKMX.LT.PMAS(KCR2,1)-20D0*PMAS(KCR2,2)) THEN
              KFR2=0
            ELSE
              TAUR2=PMAS(KCR2,1)**2/VINT(2)            
              GAMR2=PMAS(KCR2,1)*PMAS(KCR2,2)/VINT(2)
              MINT(72)=2
              MINT(74)=KFR2
              VINT(75)=TAUR2
              VINT(76)=GAMR2
            ENDIF
          ELSEIF(ITECH.EQ.0) THEN
            ALPRHT=2.16D0*(3D0/DBLE(ITCM(1)))
            ITECH=1
            KFR1=KTECHN+113              
            KCR1=PYCOMP(KFR1)
            KFR2=KTECHN+223
            KCR2=PYCOMP(KFR2)
            KFR3=KTECHN+115
            KCR3=PYCOMP(KFR3)
            IRES=0
C...Order the resonances
            IF(PMAS(KCR3,1).LT.PMAS(KCR2,1)) THEN
              KCT=KCR3
              KCR3=KCR2
              KCR2=KCT
            ENDIF
            IF(PMAS(KCR3,1).LT.PMAS(KCR1,1)) THEN
              KCT=KCR3
              KCR3=KCR1
              KCR1=KCT
            ENDIF
            IF(PMAS(KCR2,1).LT.PMAS(KCR1,1)) THEN
              KCT=KCR2
              KCR2=KCR1
              KCR1=KCT
            ENDIF
            DO 101 I=1,3
              IF(I.EQ.1) THEN
                SHN0=PMAS(KCR1,1)**2
              ELSEIF(I.EQ.2) THEN
                IF(ABS(PMAS(KCR2,1)-PMAS(KCR1,1)).LE.1D-6) GOTO 101
                SHN0=PMAS(KCR2,1)**2
              ELSEIF(I.EQ.3) THEN
                IF(ABS(PMAS(KCR3,1)-PMAS(KCR3,1)).LE.1D-6) GOTO 101
                SHN0=PMAS(KCR3,1)**2
              ENDIF
              AEM=PYALEM(SHN0)
              FAR=SQRT(AEM/ALPRHT)              
              SHN=SHN0*(1D0-FAR)
              CALL PYTECM(SHN,S1,WIDO,1)
              RES=SHN-S1
              SHN=S1*.99D0
              SHSTEP=2D0
 102          SHN=SHN+SHSTEP
              CALL PYTECM(SHN,S1,WIDO,1)
              IF(RES.LT.0D0.AND.SHN-S1.GE.0D0) THEN
                IOK=.FALSE.
                IF(IRES.GT.0) THEN
                  IF(ABS(SQRT(S1)-XMAS(IRES)).GT.1D-6) IOK=.TRUE.
                ELSEIF(IRES.EQ.0) THEN
                  IOK=.TRUE.
                ENDIF
                IF(IOK) THEN
                  IRES=IRES+1
                  XMAS(IRES)=SQRT(S1)
                  XWID(IRES)=WIDO
                ENDIF
              ENDIF
              RES=SHN-S1
              IF(IRES.LT.3.AND.SHN.LT.SHN0*(1D0+FAR)) GOTO 102
 101        CONTINUE
            JRES=0
            KFR1=KTECHN+213              
            KCR1=PYCOMP(KFR1)
            KFR2=KTECHN+215
            KCR2=PYCOMP(KFR2)
            IF(PMAS(KCR2,1).LT.PMAS(KCR1,1)) THEN
              KCT=KCR2
              KCR2=KCR1
              KCR1=KCT
            ENDIF
            DO 103 I=1,2
              IF(I.EQ.1) THEN
                SHN0=PMAS(KCR1,1)**2
              ELSEIF(I.EQ.2) THEN
                IF(ABS(PMAS(KCR2,1)-PMAS(KCR1,1)).LE.1D-6) GOTO 103
                SHN0=PMAS(KCR2,1)**2
              ENDIF
              AEM=PYALEM(SHN0)
              FAR=SQRT(AEM/ALPRHT)              
              SHN=SHN0*(1D0-FAR)
              CALL PYTECM(SHN,S1,WIDO,2)
              RES=SHN-S1
              SHN=S1*.99D0
              SHSTEP=2D0
 104          SHN=SHN+SHSTEP
              CALL PYTECM(SHN,S1,WIDO,2)
              IF(RES.LT.0D0.AND.SHN-S1.GE.0D0) THEN
                IOK=.FALSE.
                IF(JRES.GT.0) THEN
                  IF(ABS(SQRT(S1)-XMAS(IRES)).GT.1D-6) IOK=.TRUE.
                ELSEIF(JRES.EQ.0) THEN
                  IOK=.TRUE.
                ENDIF
                IF(IOK) THEN
                  JRES=JRES+1
                  YMAS(JRES)=SQRT(S1)
                  YWID(JRES)=WIDO
                ENDIF
              ENDIF
              RES=SHN-S1
              IF(JRES.LT.2.AND.SHN.LT.SHN0*(1D0+FAR)) GOTO 104
 103        CONTINUE
          ENDIF
          IF(ISUB.EQ.194.OR.(ISUB.GE.361.AND.ISUB.LE.368).OR.
     &     ISUB.EQ.379.OR.ISUB.EQ.380) THEN
            MINT(72)=IRES
            IF(IRES.GE.1) THEN
              VINT(73)=XMAS(1)**2/VINT(2)
              VINT(74)=XMAS(1)*XWID(1)/VINT(2)
              TAUR1=VINT(73)
              GAMR1=VINT(74)
              XM1=XMAS(1)
              XG1=XWID(1)
              KFR1=1
            ENDIF
            IF(IRES.GE.2) THEN
              VINT(75)=XMAS(2)**2/VINT(2)
              VINT(76)=XMAS(2)*XWID(2)/VINT(2)
              TAUR2=VINT(75)
              GAMR2=VINT(76)
              XM2=XMAS(2)
              XG2=XWID(2)
              KFR2=2
            ENDIF
            IF(IRES.EQ.3) THEN
              VINT(77)=XMAS(3)**2/VINT(2)
              VINT(78)=XMAS(3)*XWID(3)/VINT(2)
              TAUR3=VINT(77)
              GAMR3=VINT(78)
              XM3=XMAS(3)
              XG3=XWID(3)
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
              XM1=YMAS(1)
              XG1=YWID(1)
            ENDIF
            IF(JRES.GE.2) THEN
              VINT(75)=YMAS(2)**2/VINT(2)
              VINT(76)=YMAS(2)*YWID(2)/VINT(2)
              KFR2=2
              TAUR2=VINT(73)
              GAMR2=VINT(74)
              XM2=YMAS(2)
              XG2=YWID(2)
            ENDIF
            KFR3=0
          ENDIF
          IF(ISUB.NE.141) THEN
            IF(KFR1.NE.0.AND.(CKIN(1).GT.(XM1+20D0*XG1)
     &       .OR.CKMX.LT.(XM1-20D0*XG1))) KFR1=0
            IF(KFR2.NE.0.AND.(CKIN(1).GT.(XM2+20D0*XG2)
     &       .OR.CKMX.LT.(XM2-20D0*XG2))) KFR2=0
            IF(KFR3.NE.0.AND.(CKIN(1).GT.(XM3+20D0*XG3)
     &       .OR.CKMX.LT.(XM3-20D0*XG3))) KFR3=0
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
 
C...Find product masses and minimum pT of process.
        SQM3=0D0
        SQM4=0D0
        MINT(71)=0
        VINT(71)=CKIN(3)
        VINT(80)=1D0
        IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN
          NBW=0
          DO 110 I=1,2
            PMMN(I)=0D0
            IF(KFPR(ISUB,I).EQ.0) THEN
            ELSEIF(MSTP(42).LE.0.OR.PMAS(PYCOMP(KFPR(ISUB,I)),2).LT.
     &        PARP(41)) THEN
              IF(I.EQ.1) SQM3=PMAS(PYCOMP(KFPR(ISUB,I)),1)**2
              IF(I.EQ.2) SQM4=PMAS(PYCOMP(KFPR(ISUB,I)),1)**2
            ELSE
              NBW=NBW+1
C...This prevents SUSY/t particles from becoming too light.
              KFLW=KFPR(ISUB,I)
              IF(KFLW/KSUSY1.EQ.1.OR.KFLW/KSUSY1.EQ.2) THEN
                KCW=PYCOMP(KFLW)
                PMMN(I)=PMAS(KCW,1)
                DO 100 IDC=MDCY(KCW,2),MDCY(KCW,2)+MDCY(KCW,3)-1
                  IF(MDME(IDC,1).GT.0.AND.BRAT(IDC).GT.1E-4) THEN
                    PMSUM=PMAS(PYCOMP(KFDP(IDC,1)),1)+
     &              PMAS(PYCOMP(KFDP(IDC,2)),1)
                    IF(KFDP(IDC,3).NE.0) PMSUM=PMSUM+
     &              PMAS(PYCOMP(KFDP(IDC,3)),1)
                    PMMN(I)=MIN(PMMN(I),PMSUM)
                  ENDIF
  100           CONTINUE
              ELSEIF(KFLW.EQ.6) THEN
                PMMN(I)=PMAS(24,1)+PMAS(5,1)
              ENDIF
            ENDIF
  110     CONTINUE
          IF(NBW.GE.1) THEN
            CKIN41=CKIN(41)
            CKIN43=CKIN(43)
            CKIN(41)=MAX(PMMN(1),CKIN(41))
            CKIN(43)=MAX(PMMN(2),CKIN(43))
            CALL PYOFSH(3,0,KFPR(ISUB,1),KFPR(ISUB,2),0D0,PQM3,PQM4)
            CKIN(41)=CKIN41
            CKIN(43)=CKIN43
            IF(MINT(51).EQ.1) THEN
              WRITE(MSTU(11),5100) ISUB
              MSUB(ISUB)=0
              GOTO 460
            ENDIF
            SQM3=PQM3**2
            SQM4=PQM4**2
          ENDIF
          IF(MIN(SQM3,SQM4).LT.CKIN(6)**2) MINT(71)=1
          IF(MINT(71).EQ.1) VINT(71)=MAX(CKIN(3),CKIN(5))
          IF(ISUB.EQ.96.AND.MSTP(82).LE.1) THEN
            VINT(71)=PARP(81)*(VINT(1)/PARP(89))**PARP(90)
          ELSEIF(ISUB.EQ.96) THEN
            VINT(71)=0.08D0*PARP(82)*(VINT(1)/PARP(89))**PARP(90)
          ENDIF
        ENDIF
        VINT(63)=SQM3
        VINT(64)=SQM4
 
C...Prepare for additional variable choices in 2 -> 3.
        IF(ISTSB.EQ.5) THEN
          VINT(201)=0D0
          IF(KFPR(ISUB,2).GT.0) VINT(201)=PMAS(PYCOMP(KFPR(ISUB,2)),1)
          VINT(206)=VINT(201)
          IF(ISUB.EQ.401.OR.ISUB.EQ.402) VINT(206)=PMAS(5,1)
          VINT(204)=PMAS(23,1)
          IF(ISUB.EQ.124.OR.ISUB.EQ.351) VINT(204)=PMAS(24,1)
          IF(ISUB.EQ.352) VINT(204)=PMAS(PYCOMP(9900024),1)
          IF(ISUB.EQ.121.OR.ISUB.EQ.122.OR.ISUB.EQ.181.OR.ISUB.EQ.182
     &    .OR.ISUB.EQ.186.OR.ISUB.EQ.187.OR.ISUB.EQ.401.OR.ISUB.EQ.402)
     &         VINT(204)=VINT(201)
          VINT(209)=VINT(204)
          IF(ISUB.EQ.401.OR.ISUB.EQ.402) VINT(209)=VINT(206)
        ENDIF
 
C...Number of points for each variable: tau, tau', y*, cos(theta-hat).
        IPEAK7=0
        NPTS(1)=2+2*MINT(72)
        IF(MINT(47).EQ.1) THEN
          IF(ISTSB.EQ.1.OR.ISTSB.EQ.2) NPTS(1)=1
        ELSEIF(MINT(47).GE.5) THEN
          IF(ISTSB.LE.2.OR.ISTSB.GT.5) THEN
            NPTS(1)=NPTS(1)+1
            IPEAK7=1
          ENDIF
        ENDIF
        NPTS(2)=1
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
          IF(MINT(47).GE.2) NPTS(2)=2
          IF(MINT(47).GE.5) NPTS(2)=3
        ENDIF
        NPTS(3)=1
        IF(MINT(47).EQ.4.OR.MINT(47).EQ.5) THEN
          NPTS(3)=3
          IF(MINT(45).EQ.3) NPTS(3)=NPTS(3)+1
          IF(MINT(46).EQ.3) NPTS(3)=NPTS(3)+1
        ENDIF
        NPTS(4)=1
        IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) NPTS(4)=5
        NTRY=NPTS(1)*NPTS(2)*NPTS(3)*NPTS(4)
 
C...Reset coefficients of cross-section weighting.
        DO 120 J=1,20
          COEF(ISUB,J)=0D0
  120   CONTINUE
        IF(ISUB.EQ.194.OR.ISUB.EQ.195.OR.(ISUB.GE.361
     &   .AND.ISUB.LE.380)) THEN
          DO 125 J=1,2
            COEFX(ISUB,J)=0D0
 125      CONTINUE
        ENDIF
        COEF(ISUB,1)=1D0
        COEF(ISUB,8)=0.5D0
        COEF(ISUB,9)=0.5D0
        COEF(ISUB,13)=1D0
        COEF(ISUB,18)=1D0
        MCTH=0
        MTAUP=0
        METAUP=0
        VINT(23)=0D0
        VINT(26)=0D0
        SIGSAM=0D0
 
C...Find limits and select tau, y*, cos(theta-hat) and tau' values,
C...in grid of phase space points.
        CALL PYKLIM(1)
        METAU=MINT(51)
        NACC=0
        DO 150 ITRY=1,NTRY
          MINT(51)=0
          IF(METAU.EQ.1) GOTO 150
          IF(MOD(ITRY-1,NPTS(2)*NPTS(3)*NPTS(4)).EQ.0) THEN
            MTAU=1+(ITRY-1)/(NPTS(2)*NPTS(3)*NPTS(4))
            IF(MINT(72).LE.2.AND.MTAU.GT.2+2*MINT(72)) THEN
              MTAU=7
            ELSEIF(MINT(72).EQ.3.AND.IPEAK7.EQ.0.AND.MTAU.GE.7) THEN
              MTAU=MTAU+1              
            ENDIF
            RTAU=0.5D0
C...Special case when both resonances have same mass,
C...as is often the case in process 194.
c           IF(MINT(72).GE.2) THEN
c             IF(ABS(PMAS(KCR2,1)-PMAS(KCR1,1)).LT.
c    &        0.01D0*(PMAS(KCR2,1)+PMAS(KCR1,1))) THEN
c               IF(MTAU.EQ.3.OR.MTAU.EQ.4) THEN
c                 RTAU=0.4D0
c               ELSEIF(MTAU.EQ.5.OR.MTAU.EQ.6) THEN
c                 RTAU=0.6D0
c               ENDIF
c             ENDIF
c           ENDIF
            CALL PYKMAP(1,MTAU,RTAU)
            IF(ISTSB.GE.3.AND.ISTSB.LE.5) CALL PYKLIM(4)
            METAUP=MINT(51)
          ENDIF
          IF(METAUP.EQ.1) GOTO 150
          IF(ISTSB.GE.3.AND.ISTSB.LE.5.AND.MOD(ITRY-1,NPTS(3)*NPTS(4))
     &    .EQ.0) THEN
            MTAUP=1+MOD((ITRY-1)/(NPTS(3)*NPTS(4)),NPTS(2))
            CALL PYKMAP(4,MTAUP,0.5D0)
          ENDIF
          IF(MOD(ITRY-1,NPTS(3)*NPTS(4)).EQ.0) THEN
            CALL PYKLIM(2)
            MEYST=MINT(51)
          ENDIF
          IF(MEYST.EQ.1) GOTO 150
          IF(MOD(ITRY-1,NPTS(4)).EQ.0) THEN
            MYST=1+MOD((ITRY-1)/NPTS(4),NPTS(3))
            IF(MYST.EQ.4.AND.MINT(45).NE.3) MYST=5
            CALL PYKMAP(2,MYST,0.5D0)
            CALL PYKLIM(3)
            MECTH=MINT(51)
          ENDIF
          IF(MECTH.EQ.1) GOTO 150
          IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN
            MCTH=1+MOD(ITRY-1,NPTS(4))
            CALL PYKMAP(3,MCTH,0.5D0)
          ENDIF
          IF(ISUB.EQ.96) VINT(25)=VINT(21)*(1D0-VINT(23)**2)
 
C...Store position and limits.
          MINT(51)=0
          CALL PYKLIM(0)
          IF(MINT(51).EQ.1) GOTO 150
          NACC=NACC+1
          MVARPT(NACC,1)=MTAU
          MVARPT(NACC,2)=MTAUP
          MVARPT(NACC,3)=MYST
          MVARPT(NACC,4)=MCTH
          DO 130 J=1,30
            VINTPT(NACC,J)=VINT(10+J)
  130     CONTINUE
 
C...Normal case: calculate cross-section.
          IF(ISTSB.NE.5) THEN
            CALL PYSIGH(NCHN,SIGS)
            IF(MWTXS.EQ.1) THEN
              CALL PYEVWT(WTXS)
              SIGS=WTXS*SIGS
            ENDIF
 
C..2 -> 3: find highest value out of a number of tries.
          ELSE
            SIGS=0D0
            DO 140 IKIN3=1,MSTP(129)
              CALL PYKMAP(5,0,0D0)
              IF(MINT(51).EQ.1) GOTO 140
              CALL PYSIGH(NCHN,SIGTMP)
              IF(MWTXS.EQ.1) THEN
                CALL PYEVWT(WTXS)
                SIGTMP=WTXS*SIGTMP
              ENDIF
              IF(SIGTMP.GT.SIGS) SIGS=SIGTMP
  140       CONTINUE
          ENDIF
 
C...Store cross-section.
          SIGSPT(NACC)=SIGS
          IF(SIGS.GT.SIGSAM) SIGSAM=SIGS
          IF(MSTP(122).GE.2) WRITE(MSTU(11),5200) MTAU,MYST,MCTH,MTAUP,
     &    VINT(21),VINT(22),VINT(23),VINT(26),SIGS
  150   CONTINUE
        IF(NACC.EQ.0) THEN
          WRITE(MSTU(11),5100) ISUB
          MSUB(ISUB)=0
          GOTO 460
        ELSEIF(SIGSAM.EQ.0D0) THEN
          WRITE(MSTU(11),5300) ISUB
          MSUB(ISUB)=0
          GOTO 460
        ENDIF
        IF(ISUB.NE.96) NPOSI=NPOSI+1
 
C...Calculate integrals in tau over maximal phase space limits.
        TAUMIN=VINT(11)
        TAUMAX=VINT(31)
        ATAU1=LOG(TAUMAX/TAUMIN)
        IF(NPTS(1).GE.2) THEN
          ATAU2=(TAUMAX-TAUMIN)/(TAUMAX*TAUMIN)
        ENDIF
        IF(NPTS(1).GE.4) THEN
          ATAU3=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR1)/(TAUMAX+TAUR1))/TAUR1
          ATAU4=(ATAN((TAUMAX-TAUR1)/GAMR1)-ATAN((TAUMIN-TAUR1)/GAMR1))/
     &    GAMR1
        ENDIF
        IF(NPTS(1).GE.6) THEN
          ATAU5=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR2)/(TAUMAX+TAUR2))/TAUR2
          ATAU6=(ATAN((TAUMAX-TAUR2)/GAMR2)-ATAN((TAUMIN-TAUR2)/GAMR2))/
     &    GAMR2
        ENDIF
        IF(NPTS(1).GE.8) THEN
          ATAU8=LOG(TAUMAX/TAUMIN*(TAUMIN+TAUR3)/(TAUMAX+TAUR3))/TAUR3
          ATAU9=(ATAN((TAUMAX-TAUR3)/GAMR3)-ATAN((TAUMIN-TAUR3)/GAMR3))/
     &    GAMR3
        ENDIF
        IF(IPEAK7.EQ.1) THEN
          ATAU7=LOG(MAX(2D-10,1D0-TAUMIN)/MAX(2D-10,1D0-TAUMAX))
        ENDIF
 
C...Reset. Sum up cross-sections in points calculated.
        DO 320 IVAR=1,4
          IF(NPTS(IVAR).EQ.1) GOTO 320
          IF(ISUB.EQ.96.AND.IVAR.EQ.4) GOTO 320
          NBIN=NPTS(IVAR)
          DO 170 J1=1,NBIN
            NAREL(J1)=0
            WTREL(J1)=0D0
            COEFU(J1)=0D0
            DO 160 J2=1,NBIN
              WTMAT(J1,J2)=0D0
  160       CONTINUE
  170     CONTINUE
          DO 180 IACC=1,NACC
            IBIN=MVARPT(IACC,IVAR)
            IF(IVAR.EQ.1) THEN
              IF(IBIN.GT.7.AND.IPEAK7.EQ.0) THEN
                IBIN=IBIN-1
              ELSEIF(IBIN.EQ.7.AND.IPEAK7.EQ.1.AND.MSTP(72).LT.3) THEN
                IBIN=3+2*MINT(72)
              ENDIF
            ENDIF
            IF(IVAR.EQ.3.AND.IBIN.EQ.5.AND.MINT(45).NE.3) IBIN=4
            NAREL(IBIN)=NAREL(IBIN)+1
            WTREL(IBIN)=WTREL(IBIN)+SIGSPT(IACC)
 
C...Sum up tau cross-section pieces in points used.
            IF(IVAR.EQ.1) THEN
              TAU=VINTPT(IACC,11)
              WTMAT(IBIN,1)=WTMAT(IBIN,1)+1D0
              WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ATAU1/ATAU2)/TAU
              IF(NBIN.GE.4) THEN
                WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ATAU1/ATAU3)/(TAU+TAUR1)
                WTMAT(IBIN,4)=WTMAT(IBIN,4)+(ATAU1/ATAU4)*TAU/
     &          ((TAU-TAUR1)**2+GAMR1**2)
              ENDIF
              IF(NBIN.GE.6) THEN
                WTMAT(IBIN,5)=WTMAT(IBIN,5)+(ATAU1/ATAU5)/(TAU+TAUR2)
                WTMAT(IBIN,6)=WTMAT(IBIN,6)+(ATAU1/ATAU6)*TAU/
     &          ((TAU-TAUR2)**2+GAMR2**2)
              ENDIF
              IF(MINT(72).LE.2.AND.IPEAK7.EQ.1) THEN
                WTMAT(IBIN,3+2*MINT(72))=WTMAT(IBIN,3+2*MINT(72))
     &           +(ATAU1/ATAU7)*TAU/MAX(2D-10,1D0-TAU)
              ELSEIF(MINT(72).EQ.3.AND.IPEAK7.EQ.1) THEN
                WTMAT(IBIN,7)=WTMAT(IBIN,7)
     &           +(ATAU1/ATAU7)*TAU/MAX(2D-10,1D0-TAU)
              ENDIF
              IF(MINT(72).EQ.3) THEN
                WTMAT(IBIN,7+IPEAK7)=WTMAT(IBIN,7+IPEAK7)
     &           +(ATAU1/ATAU8)/(TAU+TAUR3)
                WTMAT(IBIN,8+IPEAK7)=WTMAT(IBIN,8+IPEAK7)
     &           +(ATAU1/ATAU9)*TAU/((TAU-TAUR3)**2+GAMR3**2)
              ENDIF
C...Sum up tau' cross-section pieces in points used.
            ELSEIF(IVAR.EQ.2) THEN
              TAU=VINTPT(IACC,11)
              TAUP=VINTPT(IACC,16)
              TAUPMN=VINTPT(IACC,6)
              TAUPMX=VINTPT(IACC,26)
              ATAUP1=LOG(TAUPMX/TAUPMN)
              ATAUP2=((1D0-TAU/TAUPMX)**4-(1D0-TAU/TAUPMN)**4)/(4D0*TAU)
              WTMAT(IBIN,1)=WTMAT(IBIN,1)+1D0
              WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ATAUP1/ATAUP2)*
     &        (1D0-TAU/TAUP)**3/TAUP
              IF(NBIN.GE.3) THEN
                ATAUP3=LOG(MAX(2D-10,1D0-TAUPMN)/MAX(2D-10,1D0-TAUPMX))
                WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ATAUP1/ATAUP3)*
     &          TAUP/MAX(2D-10,1D0-TAUP)
              ENDIF
 
C...Sum up y* cross-section pieces in points used.
            ELSEIF(IVAR.EQ.3) THEN
              YST=VINTPT(IACC,12)
              YSTMIN=VINTPT(IACC,2)
              YSTMAX=VINTPT(IACC,22)
              AYST0=YSTMAX-YSTMIN
              AYST1=0.5D0*(YSTMAX-YSTMIN)**2
              AYST2=AYST1
              AYST3=2D0*(ATAN(EXP(YSTMAX))-ATAN(EXP(YSTMIN)))
              WTMAT(IBIN,1)=WTMAT(IBIN,1)+(AYST0/AYST1)*(YST-YSTMIN)
              WTMAT(IBIN,2)=WTMAT(IBIN,2)+(AYST0/AYST2)*(YSTMAX-YST)
              WTMAT(IBIN,3)=WTMAT(IBIN,3)+(AYST0/AYST3)/COSH(YST)
              IF(MINT(45).EQ.3) THEN
                TAUE=VINTPT(IACC,11)
                IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=VINTPT(IACC,16)
                YST0=-0.5D0*LOG(TAUE)
                AYST4=LOG(MAX(1D-10,EXP(YST0-YSTMIN)-1D0)/
     &          MAX(1D-10,EXP(YST0-YSTMAX)-1D0))
                WTMAT(IBIN,4)=WTMAT(IBIN,4)+(AYST0/AYST4)/
     &          MAX(1D-10,1D0-EXP(YST-YST0))
              ENDIF
              IF(MINT(46).EQ.3) THEN
                TAUE=VINTPT(IACC,11)
                IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=VINTPT(IACC,16)
                YST0=-0.5D0*LOG(TAUE)
                AYST5=LOG(MAX(1D-10,EXP(YST0+YSTMAX)-1D0)/
     &          MAX(1D-10,EXP(YST0+YSTMIN)-1D0))
                WTMAT(IBIN,NBIN)=WTMAT(IBIN,NBIN)+(AYST0/AYST5)/
     &          MAX(1D-10,1D0-EXP(-YST-YST0))
              ENDIF
 
C...Sum up cos(theta-hat) cross-section pieces in points used.
            ELSE
              RM34=MAX(1D-20,2D0*SQM3*SQM4/(VINTPT(IACC,11)*VINT(2))**2)
              RSQM=1D0+RM34
              CTHMAX=SQRT(1D0-4D0*VINT(71)**2/(TAUMAX*VINT(2)))
              CTHMIN=-CTHMAX
              IF(CTHMAX.GT.0.9999D0) RM34=MAX(RM34,2D0*VINT(71)**2/
     &        (TAUMAX*VINT(2)))
              ACTH1=CTHMAX-CTHMIN
              ACTH2=LOG(MAX(RM34,RSQM-CTHMIN)/MAX(RM34,RSQM-CTHMAX))
              ACTH3=LOG(MAX(RM34,RSQM+CTHMAX)/MAX(RM34,RSQM+CTHMIN))
              ACTH4=1D0/MAX(RM34,RSQM-CTHMAX)-1D0/MAX(RM34,RSQM-CTHMIN)
              ACTH5=1D0/MAX(RM34,RSQM+CTHMIN)-1D0/MAX(RM34,RSQM+CTHMAX)
              CTH=VINTPT(IACC,13)
              WTMAT(IBIN,1)=WTMAT(IBIN,1)+1D0
              WTMAT(IBIN,2)=WTMAT(IBIN,2)+(ACTH1/ACTH2)/
     &        MAX(RM34,RSQM-CTH)
              WTMAT(IBIN,3)=WTMAT(IBIN,3)+(ACTH1/ACTH3)/
     &        MAX(RM34,RSQM+CTH)
              WTMAT(IBIN,4)=WTMAT(IBIN,4)+(ACTH1/ACTH4)/
     &        MAX(RM34,RSQM-CTH)**2
              WTMAT(IBIN,5)=WTMAT(IBIN,5)+(ACTH1/ACTH5)/
     &        MAX(RM34,RSQM+CTH)**2
            ENDIF
  180     CONTINUE
 
C...Check that equation system solvable.
          IF(MSTP(122).GE.2) WRITE(MSTU(11),5400) CVAR(IVAR)
          MSOLV=1
          WTRELS=0D0
          DO 190 IBIN=1,NBIN
            IF(MSTP(122).GE.2) WRITE(MSTU(11),5500) (WTMAT(IBIN,IRED),
     &      IRED=1,NBIN),WTREL(IBIN)
            IF(NAREL(IBIN).EQ.0) MSOLV=0
            WTRELS=WTRELS+WTREL(IBIN)
  190     CONTINUE
          IF(ABS(WTRELS).LT.1D-20) MSOLV=0
 
C...Solve to find relative importance of cross-section pieces.
          IF(MSOLV.EQ.1) THEN
            DO 200 IBIN=1,NBIN
              WTRELN(IBIN)=MAX(0.1D0,WTREL(IBIN)/WTRELS)
              WTRSAV(IBIN)=WTREL(IBIN)
  200       CONTINUE
C...Auxiliary vectors to record order of permutations
            DO I=1,NBIN
              IP(I) = I
              IQ(I) = I
            ENDDO
            DO 230 IRED=1,NBIN-1
              MROW=IRED
              RESMAX=ABS(WTREL(MROW))
C...Find row with largest residual
              DO JBIN=IRED+1,NBIN
                IF(RESMAX.LT.ABS(WTREL(JBIN))) THEN
                  MROW=JBIN
                  RESMAX=ABS(WTREL(MROW))
                ENDIF
              ENDDO
              IF(RESMAX.LT.1D-20) THEN
                MSOLV=0
                GOTO 260
              ENDIF
              MCOL = IRED
              AMAX = ABS(WTMAT(MROW,MCOL))
C...Find column with largest entry
              DO JBIN=IRED+1,NBIN
                IF (AMAX.LT.ABS(WTMAT(MROW,JBIN))) THEN
                  MCOL = JBIN
                  AMAX = ABS(WTMAT(MROW,MCOL))
                ENDIF
              ENDDO
C...Swap rows if necessary
              IF(MROW.NE.IRED) THEN
                DO JBIN=1,NBIN
                  TMPE=WTMAT(IRED,JBIN)
                  WTMAT(IRED,JBIN)=WTMAT(MROW,JBIN)
                  WTMAT(MROW,JBIN)=TMPE
                ENDDO
                TMPE=WTREL(IRED)
                WTREL(IRED)=WTREL(MROW)
                WTREL(MROW)=TMPE
                MTMP=IQ(IRED)
                IQ(IRED)=IQ(MROW)
                IQ(MROW)=MTMP
              ENDIF
C...Swap columns if necessary
              IF(MCOL.NE.IRED) THEN
                DO JBIN=1,NBIN
                  TMPE=WTMAT(JBIN,IRED)
                  WTMAT(JBIN,IRED)=WTMAT(JBIN,MCOL)
                  WTMAT(JBIN,MCOL)=TMPE
                ENDDO
                MTMP=IP(IRED)
                IP(IRED)=IP(MCOL)
                IP(MCOL)=MTMP
              ENDIF
C...Begin eliminating equations
              DO 220 IBIN=IRED+1,NBIN
                IF(ABS(WTMAT(IRED,IRED)).LT.1D-20) THEN
                  MSOLV=0
                  GOTO 260
                ENDIF
C                RQT=WTMAT(IBIN,IRED)/WTMAT(IRED,IRED)
                RQTU=WTMAT(IBIN,IRED)
                RQTL=WTMAT(IRED,IRED)
C...Switch order of operations
                WTREL(IBIN)=WTREL(IBIN)-RQTU*
     $            (WTREL(IRED)/RQTL)
                DO 210 ICOE=IRED,NBIN
                   WTMAT(IBIN,ICOE)=WTMAT(IBIN,ICOE)-
     $                RQTU*(WTMAT(IRED,ICOE)/RQTL)
  210           CONTINUE
  220         CONTINUE
  230       CONTINUE
            DO 250 IRED=NBIN,1,-1
              DO 240 ICOE=IRED+1,NBIN
                WTREL(IRED)=WTREL(IRED)-WTMAT(IRED,ICOE)*COEFU(ICOE)
  240         CONTINUE
              IF(ABS(WTMAT(IRED,IRED)).LT.1D-20) THEN
                MSOLV=0
                GOTO 260
              ENDIF
              COEFU(IRED)=WTREL(IRED)/WTMAT(IRED,IRED)
              TEMPC(IRED)=COEFU(IRED)
  250       CONTINUE
C...Return to original order
            DO IBIN=1,NBIN
              MTMP=IP(IBIN)
              COEFU(MTMP)=TEMPC(IBIN)
            ENDDO
          ENDIF
 
C...Share evenly if failure.
  260     IF(MSOLV.EQ.0) THEN
            DO 270 IBIN=1,NBIN
              COEFU(IBIN)=1D0
              WTRELN(IBIN)=0.1D0
              IF(WTRELS.GT.0D0) WTRELN(IBIN)=MAX(0.1D0,
     &        WTRSAV(IBIN)/WTRELS)
  270       CONTINUE
          ENDIF
 
C...Normalize coefficients, with piece shared democratically.
          COEFSU=0D0
          WTRELS=0D0
          DO 280 IBIN=1,NBIN
            COEFU(IBIN)=MAX(0D0,COEFU(IBIN))
            COEFSU=COEFSU+COEFU(IBIN)
            WTRELS=WTRELS+WTRELN(IBIN)
  280     CONTINUE
          IF(COEFSU.GT.0D0) THEN
            DO 290 IBIN=1,NBIN
              COEFO(IBIN)=PARP(122)/NBIN+(1D0-PARP(122))*0.5D0*
     &        (COEFU(IBIN)/COEFSU+WTRELN(IBIN)/WTRELS)
  290       CONTINUE
          ELSE
            DO 300 IBIN=1,NBIN
              COEFO(IBIN)=1D0/NBIN
  300       CONTINUE
          ENDIF
          IF(IVAR.EQ.1) IOFF=0
          IF(IVAR.EQ.2) IOFF=17
          IF(IVAR.EQ.3) IOFF=7
          IF(IVAR.EQ.4) IOFF=12
          DO 310 IBIN=1,NBIN
            ICOF=IOFF+IBIN
            IF(IVAR.EQ.1) THEN
              IF(IBIN.EQ.NBIN.AND.(MINT(72).LE.2.AND.IPEAK7.EQ.1)) THEN
                ICOF=7
              ENDIF
            ENDIF
            IF(IVAR.EQ.3.AND.IBIN.EQ.4.AND.MINT(45).NE.3) ICOF=ICOF+1
            IF(IVAR.EQ.1.AND.IBIN.GE.7+IPEAK7.AND.MINT(72).EQ.3) THEN
              COEFX(ISUB,IBIN-6-IPEAK7)=COEFO(IBIN)
            ELSE
              COEF(ISUB,ICOF)=COEFO(IBIN)
            ENDIF
  310     CONTINUE
          
          IF(MSTP(122).GE.2) WRITE(MSTU(11),5600) CVAR(IVAR),
     &       (COEFO(IBIN),IBIN=1,NBIN)

  320   CONTINUE
 
C...Find two most promising maxima among points previously determined.
        DO 330 J=1,4
          IACCMX(J)=0
          SIGSMX(J)=0D0
  330   CONTINUE
        NMAX=0
        DO 390 IACC=1,NACC
          DO 340 J=1,30
            VINT(10+J)=VINTPT(IACC,J)
  340     CONTINUE
          IF(ISTSB.NE.5) THEN
            CALL PYSIGH(NCHN,SIGS)
            IF(MWTXS.EQ.1) THEN
              CALL PYEVWT(WTXS)
              SIGS=WTXS*SIGS
            ENDIF
          ELSE
            SIGS=0D0
            DO 350 IKIN3=1,MSTP(129)
              CALL PYKMAP(5,0,0D0)
              IF(MINT(51).EQ.1) GOTO 350
              CALL PYSIGH(NCHN,SIGTMP)
              IF(MWTXS.EQ.1) THEN
                CALL PYEVWT(WTXS)
                SIGTMP=WTXS*SIGTMP
              ENDIF
              IF(SIGTMP.GT.SIGS) SIGS=SIGTMP
  350       CONTINUE
          ENDIF
          IEQ=0
          DO 360 IMV=1,NMAX
            IF(ABS(SIGS-SIGSMX(IMV)).LT.1D-4*(SIGS+SIGSMX(IMV))) IEQ=IMV
  360     CONTINUE
          IF(IEQ.EQ.0) THEN
            DO 370 IMV=NMAX,1,-1
              IIN=IMV+1
              IF(SIGS.LE.SIGSMX(IMV)) GOTO 380
              IACCMX(IMV+1)=IACCMX(IMV)
              SIGSMX(IMV+1)=SIGSMX(IMV)
  370       CONTINUE
            IIN=1
  380       IACCMX(IIN)=IACC
            SIGSMX(IIN)=SIGS
            IF(NMAX.LE.1) NMAX=NMAX+1
          ENDIF
  390   CONTINUE
 
C...Read out starting position for search.
        IF(MSTP(122).GE.2) WRITE(MSTU(11),5700)
        SIGSAM=SIGSMX(1)
        DO 440 IMAX=1,NMAX
          IACC=IACCMX(IMAX)
          MTAU=MVARPT(IACC,1)
          MTAUP=MVARPT(IACC,2)
          MYST=MVARPT(IACC,3)
          MCTH=MVARPT(IACC,4)
          VTAU=0.5D0
          VYST=0.5D0
          VCTH=0.5D0
          VTAUP=0.5D0
 
C...Starting point and step size in parameter space.
          DO 430 IRPT=1,2
            DO 420 IVAR=1,4
              IF(NPTS(IVAR).EQ.1) GOTO 420
              IF(IVAR.EQ.1) VVAR=VTAU
              IF(IVAR.EQ.2) VVAR=VTAUP
              IF(IVAR.EQ.3) VVAR=VYST
              IF(IVAR.EQ.4) VVAR=VCTH
              IF(IVAR.EQ.1) MVAR=MTAU
              IF(IVAR.EQ.2) MVAR=MTAUP
              IF(IVAR.EQ.3) MVAR=MYST
              IF(IVAR.EQ.4) MVAR=MCTH
              IF(IRPT.EQ.1) VDEL=0.1D0
              IF(IRPT.EQ.2) VDEL=MAX(0.01D0,MIN(0.05D0,VVAR-0.02D0,
     &        0.98D0-VVAR))
              IF(IRPT.EQ.1) VMAR=0.02D0
              IF(IRPT.EQ.2) VMAR=0.002D0
              IMOV0=1
              IF(IRPT.EQ.1.AND.IVAR.EQ.1) IMOV0=0
              DO 410 IMOV=IMOV0,8
 
C...Define new point in parameter space.
                IF(IMOV.EQ.0) THEN
                  INEW=2
                  VNEW=VVAR
                ELSEIF(IMOV.EQ.1) THEN
                  INEW=3
                  VNEW=VVAR+VDEL
                ELSEIF(IMOV.EQ.2) THEN
                  INEW=1
                  VNEW=VVAR-VDEL
                ELSEIF(SIGSSM(3).GE.MAX(SIGSSM(1),SIGSSM(2)).AND.
     &            VVAR+2D0*VDEL.LT.1D0-VMAR) THEN
                  VVAR=VVAR+VDEL
                  SIGSSM(1)=SIGSSM(2)
                  SIGSSM(2)=SIGSSM(3)
                  INEW=3
                  VNEW=VVAR+VDEL
                ELSEIF(SIGSSM(1).GE.MAX(SIGSSM(2),SIGSSM(3)).AND.
     &            VVAR-2D0*VDEL.GT.VMAR) THEN
                  VVAR=VVAR-VDEL
                  SIGSSM(3)=SIGSSM(2)
                  SIGSSM(2)=SIGSSM(1)
                  INEW=1
                  VNEW=VVAR-VDEL
                ELSEIF(SIGSSM(3).GE.SIGSSM(1)) THEN
                  VDEL=0.5D0*VDEL
                  VVAR=VVAR+VDEL
                  SIGSSM(1)=SIGSSM(2)
                  INEW=2
                  VNEW=VVAR
                ELSE
                  VDEL=0.5D0*VDEL
                  VVAR=VVAR-VDEL
                  SIGSSM(3)=SIGSSM(2)
                  INEW=2
                  VNEW=VVAR
                ENDIF
 
C...Convert to relevant variables and find derived new limits.
                ILERR=0
                IF(IVAR.EQ.1) THEN
                  VTAU=VNEW
                  CALL PYKMAP(1,MTAU,VTAU)
                  IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
                    CALL PYKLIM(4)
                    IF(MINT(51).EQ.1) ILERR=1
                  ENDIF
                ENDIF
                IF(IVAR.LE.2.AND.ISTSB.GE.3.AND.ISTSB.LE.5.AND.
     &          ILERR.EQ.0) THEN
                  IF(IVAR.EQ.2) VTAUP=VNEW
                  CALL PYKMAP(4,MTAUP,VTAUP)
                ENDIF
                IF(IVAR.LE.2.AND.ILERR.EQ.0) THEN
                  CALL PYKLIM(2)
                  IF(MINT(51).EQ.1) ILERR=1
                ENDIF
                IF(IVAR.LE.3.AND.ILERR.EQ.0) THEN
                  IF(IVAR.EQ.3) VYST=VNEW
                  CALL PYKMAP(2,MYST,VYST)
                  CALL PYKLIM(3)
                  IF(MINT(51).EQ.1) ILERR=1
                ENDIF
                IF((ISTSB.EQ.2.OR.ISTSB.EQ.4.OR.ISTSB.EQ.6).AND.
     &          ILERR.EQ.0) THEN
                  IF(IVAR.EQ.4) VCTH=VNEW
                  CALL PYKMAP(3,MCTH,VCTH)
                ENDIF
                IF(ISUB.EQ.96) VINT(25)=VINT(21)*(1.-VINT(23)**2)
 
C...Evaluate cross-section. Save new maximum. Final maximum.
                IF(ILERR.NE.0) THEN
                   SIGS=0.
                ELSEIF(ISTSB.NE.5) THEN
                  CALL PYSIGH(NCHN,SIGS)
                  IF(MWTXS.EQ.1) THEN
                    CALL PYEVWT(WTXS)
                    SIGS=WTXS*SIGS
                  ENDIF
                ELSE
                  SIGS=0D0
                  DO 400 IKIN3=1,MSTP(129)
                    CALL PYKMAP(5,0,0D0)
                    IF(MINT(51).EQ.1) GOTO 400
                    CALL PYSIGH(NCHN,SIGTMP)
                    IF(MWTXS.EQ.1) THEN
                        CALL PYEVWT(WTXS)
                        SIGTMP=WTXS*SIGTMP
                    ENDIF
                    IF(SIGTMP.GT.SIGS) SIGS=SIGTMP
  400             CONTINUE
                ENDIF
                SIGSSM(INEW)=SIGS
                IF(SIGS.GT.SIGSAM) SIGSAM=SIGS
                IF(MSTP(122).GE.2) WRITE(MSTU(11),5800) IMAX,IVAR,MVAR,
     &          IMOV,VNEW,VINT(21),VINT(22),VINT(23),VINT(26),SIGS
  410         CONTINUE
  420       CONTINUE
  430     CONTINUE
  440   CONTINUE
        IF(MSTP(121).EQ.1) SIGSAM=PARP(121)*SIGSAM
        XSEC(ISUB,1)=1.05D0*SIGSAM
C...Add extra headroom for UED
        IF(ISUB.GT.310.AND.ISUB.LT.320) XSEC(ISUB,1)=XSEC(ISUB,1)*1.1D0
        IF(MINT(141).NE.0.OR.MINT(142).NE.0) XSEC(ISUB,1)=
     &  WTGAGA*XSEC(ISUB,1)
  450   CONTINUE
        IF(MSTP(173).EQ.1.AND.ISUB.NE.96) XSEC(ISUB,1)=
     &  PARP(174)*XSEC(ISUB,1)
        IF(ISUB.NE.96) XSEC(0,1)=XSEC(0,1)+XSEC(ISUB,1)
  460 CONTINUE
      MINT(51)=0
 
C...Print summary table.
      IF(MINT(121).EQ.1.AND.NPOSI.EQ.0) THEN
        IF(MSTP(127).NE.1) THEN
          WRITE(MSTU(11),5900)
          CALL PYSTOP(1)
        ELSE
          WRITE(MSTU(11),6400)
          MSTI(53)=1
        ENDIF
      ENDIF
      IF(MSTP(122).GE.1) THEN
        WRITE(MSTU(11),6000)
        WRITE(MSTU(11),6100)
        DO 470 ISUB=1,500
          IF(MSUB(ISUB).NE.1.AND.ISUB.NE.96) GOTO 470
          IF(ISUB.EQ.96.AND.MINT(50).EQ.0) GOTO 470
          IF(ISUB.EQ.96.AND.MSUB(95).NE.1.AND.MOD(MSTP(81),10).LE.0)
     &    GOTO 470
          IF(ISUB.EQ.96.AND.MINT(49).EQ.0.AND.MSTP(131).EQ.0) GOTO 470
          IF(MSUB(95).EQ.1.AND.(ISUB.EQ.11.OR.ISUB.EQ.12.OR.ISUB.EQ.13
     &    .OR.ISUB.EQ.28.OR.ISUB.EQ.53.OR.ISUB.EQ.68)) GOTO 470
          IF(MSUB(95).EQ.1.AND.ISUB.GE.381.AND.ISUB.LE.386) GOTO 470
          WRITE(MSTU(11),6200) ISUB,PROC(ISUB),XSEC(ISUB,1)
  470   CONTINUE
        WRITE(MSTU(11),6300)
      ENDIF
 
C...Format statements for maximization results.
 5000 FORMAT(/1X,'Coefficient optimization and maximum search for ',
     &'subprocess no',I4/1X,'Coefficient modes     tau',10X,'y*',9X,
     &'cth',9X,'tau''',7X,'sigma')
 5100 FORMAT(1X,'Warning: requested subprocess ',I3,' has no allowed ',
     &'phase space.'/1X,'Process switched off!')
 5200 FORMAT(1X,4I4,F12.8,F12.6,F12.7,F12.8,1P,D12.4)
 5300 FORMAT(1X,'Warning: requested subprocess ',I3,' has vanishing ',
     &'cross-section.'/1X,'Process switched off!')
 5400 FORMAT(1X,'Coefficients of equation system to be solved for ',A4)
 5500 FORMAT(1X,1P,10D11.3)
 5600 FORMAT(1X,'Result for ',A4,':',9F9.4)
 5700 FORMAT(1X,'Maximum search for given coefficients'/2X,'MAX VAR ',
     &'MOD MOV   VNEW',7X,'tau',7X,'y*',8X,'cth',7X,'tau''',7X,'sigma')
 5800 FORMAT(1X,4I4,F8.4,F11.7,F9.3,F11.6,F11.7,1P,D12.4)
 5900 FORMAT(1X,'Error: no requested process has non-vanishing ',
     &'cross-section.'/1X,'Execution stopped!')
 6000 FORMAT(/1X,8('*'),1X,'PYMAXI: summary of differential ',
     &'cross-section maximum search',1X,8('*'))
 6100 FORMAT(/11X,58('=')/11X,'I',38X,'I',17X,'I'/11X,'I  ISUB  ',
     &'Subprocess name',15X,'I  Maximum value  I'/11X,'I',38X,'I',
     &17X,'I'/11X,58('=')/11X,'I',38X,'I',17X,'I')
 6200 FORMAT(11X,'I',2X,I3,3X,A28,2X,'I',2X,1P,D12.4,3X,'I')
 6300 FORMAT(11X,'I',38X,'I',17X,'I'/11X,58('='))
 6400 FORMAT(1X,'Error: no requested process has non-vanishing ',
     &'cross-section.'/
     &1X,'Execution will stop if you try to generate events.')
 
      RETURN
      END
