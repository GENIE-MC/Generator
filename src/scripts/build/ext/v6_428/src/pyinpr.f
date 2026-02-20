 
C*********************************************************************
 
C...PYINPR
C...Selects partonic subprocesses to be included in the simulation.
 
      SUBROUTINE PYINPR
 
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
 
C...Commonblocks and character variables.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT6/PROC(0:500)
      CHARACTER PROC*28
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,/PYINT1/,
     &/PYINT2/,/PYINT6/
      CHARACTER CHIPR*10

 
C...Reset processes to be included.
      IF(MSEL.NE.0) THEN
        DO 100 I=1,500
          MSUB(I)=0
  100   CONTINUE
      ENDIF
 
C...Set running pTmin scale.
      IF(MSTP(82).LE.1) THEN
        PTMRUN=PARP(81)*(VINT(1)/PARP(89))**PARP(90)
      ELSE
        PTMRUN=PARP(82)*(VINT(1)/PARP(89))**PARP(90)
      ENDIF
 
C...Begin by assuming incoming photon to enter subprocess.
      IF(MINT(11).EQ.22) MINT(15)=22
      IF(MINT(12).EQ.22) MINT(16)=22
 
C...For e-gamma with MSTP(14)=10 allow mixture of VMD and anomalous.
      IF(MINT(121).EQ.2.AND.MSTP(14).EQ.10) THEN
        MSUB(10)=1
        MINT(123)=MINT(122)+1
 
C...For gamma-p or gamma-gamma with MSTP(14) = 10, 20, 25 or 30
C...allow mixture.
C...Here also set a few parameters otherwise normally not touched.
      ELSEIF(MINT(121).GT.1) THEN
 
C...Parton distributions dampened at small Q2; go to low energies,
C...alpha_s <1; no minimum pT cut-off a priori.
        IF(MSTP(18).EQ.2) THEN
          MSTP(57)=3
          PARP(2)=2D0
          PARU(115)=1D0
          CKIN(5)=0.2D0
          CKIN(6)=0.2D0
        ENDIF
 
C...Define pT cut-off parameters and whether run involves low-pT.
        PTMVMD=PTMRUN
        VINT(154)=PTMVMD
        PTMDIR=PTMVMD
        IF(MSTP(18).EQ.2) PTMDIR=PARP(15)
        PTMANO=PTMVMD
        IF(MSTP(15).EQ.5) PTMANO=0.60D0+
     &  0.125D0*LOG(1D0+0.10D0*VINT(1))**2
        IPTL=1
        IF(VINT(285).GT.MAX(PTMVMD,PTMDIR,PTMANO)) IPTL=0
        IF(MSEL.EQ.2) IPTL=1
 
C...Set up for p/gamma * gamma; real or virtual photons.
        IF(MINT(121).EQ.3.OR.MINT(121).EQ.6.OR.(MINT(121).EQ.4.AND.
     &  MSTP(14).EQ.30)) THEN
 
C...Set up for p/VMD * VMD.
        IF(MINT(122).EQ.1) THEN
          MINT(123)=2
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(IPTL.EQ.1) MSUB(95)=1
          IF(MSEL.EQ.2) THEN
            MSUB(91)=1
            MSUB(92)=1
            MSUB(93)=1
            MSUB(94)=1
          ENDIF
          IF(IPTL.EQ.1) CKIN(3)=0D0
 
C...Set up for p/VMD * direct gamma.
        ELSEIF(MINT(122).EQ.2) THEN
          MINT(123)=0
          IF(MINT(121).EQ.6) MINT(123)=5
          MSUB(131)=1
          MSUB(132)=1
          MSUB(135)=1
          MSUB(136)=1
          IF(IPTL.EQ.1) CKIN(3)=PTMDIR
 
C...Set up for p/VMD * anomalous gamma.
        ELSEIF(MINT(122).EQ.3) THEN
          MINT(123)=3
          IF(MINT(121).EQ.6) MINT(123)=7
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(IPTL.EQ.1) MSUB(95)=1
          IF(MSEL.EQ.2) THEN
            MSUB(91)=1
            MSUB(92)=1
            MSUB(93)=1
            MSUB(94)=1
          ENDIF
          IF(IPTL.EQ.1) CKIN(3)=0D0
 
C...Set up for DIS * p.
        ELSEIF(MINT(122).EQ.4.AND.(IABS(MINT(11)).GT.100.OR.
     &  IABS(MINT(12)).GT.100)) THEN
          MINT(123)=8
          IF(IPTL.EQ.1) MSUB(99)=1
 
C...Set up for direct * direct gamma (switch off leptons).
        ELSEIF(MINT(122).EQ.4) THEN
          MINT(123)=0
          MSUB(137)=1
          MSUB(138)=1
          MSUB(139)=1
          MSUB(140)=1
          DO 110 II=MDCY(22,2),MDCY(22,2)+MDCY(22,3)-1
            IF(IABS(KFDP(II,1)).GE.10) MDME(II,1)=MIN(0,MDME(II,1))
  110     CONTINUE
          IF(IPTL.EQ.1) CKIN(3)=PTMDIR
 
C...Set up for direct * anomalous gamma.
        ELSEIF(MINT(122).EQ.5) THEN
          MINT(123)=6
          MSUB(131)=1
          MSUB(132)=1
          MSUB(135)=1
          MSUB(136)=1
          IF(IPTL.EQ.1) CKIN(3)=PTMANO
 
C...Set up for anomalous * anomalous gamma.
        ELSEIF(MINT(122).EQ.6) THEN
          MINT(123)=3
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(IPTL.EQ.1) MSUB(95)=1
          IF(MSEL.EQ.2) THEN
            MSUB(91)=1
            MSUB(92)=1
            MSUB(93)=1
            MSUB(94)=1
          ENDIF
          IF(IPTL.EQ.1) CKIN(3)=0D0
        ENDIF
 
C...Set up for gamma* * gamma*; virtual photons = dir, VMD, anom.
        ELSEIF(MINT(121).EQ.9.OR.MINT(121).EQ.13) THEN
 
C...Set up for direct * direct gamma (switch off leptons).
        IF(MINT(122).EQ.1) THEN
          MINT(123)=0
          MSUB(137)=1
          MSUB(138)=1
          MSUB(139)=1
          MSUB(140)=1
          DO 120 II=MDCY(22,2),MDCY(22,2)+MDCY(22,3)-1
            IF(IABS(KFDP(II,1)).GE.10) MDME(II,1)=MIN(0,MDME(II,1))
  120     CONTINUE
          IF(IPTL.EQ.1) CKIN(3)=PTMDIR
 
C...Set up for direct * VMD and VMD * direct gamma.
        ELSEIF(MINT(122).EQ.2.OR.MINT(122).EQ.4) THEN
          MINT(123)=5
          MSUB(131)=1
          MSUB(132)=1
          MSUB(135)=1
          MSUB(136)=1
          IF(IPTL.EQ.1) CKIN(3)=PTMDIR
 
C...Set up for direct * anomalous and anomalous * direct gamma.
        ELSEIF(MINT(122).EQ.3.OR.MINT(122).EQ.7) THEN
          MINT(123)=6
          MSUB(131)=1
          MSUB(132)=1
          MSUB(135)=1
          MSUB(136)=1
          IF(IPTL.EQ.1) CKIN(3)=PTMANO
 
C...Set up for VMD*VMD.
        ELSEIF(MINT(122).EQ.5) THEN
          MINT(123)=2
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(IPTL.EQ.1) MSUB(95)=1
          IF(MSEL.EQ.2) THEN
            MSUB(91)=1
            MSUB(92)=1
            MSUB(93)=1
            MSUB(94)=1
          ENDIF
          IF(IPTL.EQ.1) CKIN(3)=0D0
 
C...Set up for VMD * anomalous and anomalous * VMD gamma.
        ELSEIF(MINT(122).EQ.6.OR.MINT(122).EQ.8) THEN
          MINT(123)=7
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(IPTL.EQ.1) MSUB(95)=1
          IF(MSEL.EQ.2) THEN
            MSUB(91)=1
            MSUB(92)=1
            MSUB(93)=1
            MSUB(94)=1
          ENDIF
          IF(IPTL.EQ.1) CKIN(3)=0D0
 
C...Set up for anomalous * anomalous gamma.
        ELSEIF(MINT(122).EQ.9) THEN
          MINT(123)=3
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(IPTL.EQ.1) MSUB(95)=1
          IF(MSEL.EQ.2) THEN
            MSUB(91)=1
            MSUB(92)=1
            MSUB(93)=1
            MSUB(94)=1
          ENDIF
          IF(IPTL.EQ.1) CKIN(3)=0D0
 
C...Set up for DIS * VMD and VMD * DIS gamma.
        ELSEIF(MINT(122).EQ.10.OR.MINT(122).EQ.12) THEN
          MINT(123)=8
          IF(IPTL.EQ.1) MSUB(99)=1
 
C...Set up for DIS * anomalous and anomalous * DIS gamma.
        ELSEIF(MINT(122).EQ.11.OR.MINT(122).EQ.13) THEN
          MINT(123)=9
          IF(IPTL.EQ.1) MSUB(99)=1
        ENDIF
 
C...Set up for gamma* * p; virtual photons = dir, res.
        ELSEIF(MINT(121).EQ.2) THEN
 
C...Set up for direct * p.
        IF(MINT(122).EQ.1) THEN
          MINT(123)=0
          MSUB(131)=1
          MSUB(132)=1
          MSUB(135)=1
          MSUB(136)=1
          IF(IPTL.EQ.1) CKIN(3)=PTMDIR
 
C...Set up for resolved * p.
        ELSEIF(MINT(122).EQ.2) THEN
          MINT(123)=1
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(IPTL.EQ.1) MSUB(95)=1
          IF(MSEL.EQ.2) THEN
            MSUB(91)=1
            MSUB(92)=1
            MSUB(93)=1
            MSUB(94)=1
          ENDIF
          IF(IPTL.EQ.1) CKIN(3)=0D0
        ENDIF
 
C...Set up for gamma* * gamma*; virtual photons = dir, res.
        ELSEIF(MINT(121).EQ.4) THEN
 
C...Set up for direct * direct gamma (switch off leptons).
        IF(MINT(122).EQ.1) THEN
          MINT(123)=0
          MSUB(137)=1
          MSUB(138)=1
          MSUB(139)=1
          MSUB(140)=1
          DO 130 II=MDCY(22,2),MDCY(22,2)+MDCY(22,3)-1
            IF(IABS(KFDP(II,1)).GE.10) MDME(II,1)=MIN(0,MDME(II,1))
  130     CONTINUE
          IF(IPTL.EQ.1) CKIN(3)=PTMDIR
 
C...Set up for direct * resolved and resolved * direct gamma.
        ELSEIF(MINT(122).EQ.2.OR.MINT(122).EQ.3) THEN
          MINT(123)=5
          MSUB(131)=1
          MSUB(132)=1
          MSUB(135)=1
          MSUB(136)=1
          IF(IPTL.EQ.1) CKIN(3)=PTMDIR
 
C...Set up for resolved * resolved gamma.
        ELSEIF(MINT(122).EQ.4) THEN
          MINT(123)=2
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          IF(IPTL.EQ.1) MSUB(95)=1
          IF(MSEL.EQ.2) THEN
            MSUB(91)=1
            MSUB(92)=1
            MSUB(93)=1
            MSUB(94)=1
          ENDIF
          IF(IPTL.EQ.1) CKIN(3)=0D0
        ENDIF
 
C...End of special set up for gamma-p and gamma-gamma.
        ENDIF
        CKIN(1)=2D0*CKIN(3)
      ENDIF
 
C...Flavour information for individual beams.
      DO 140 I=1,2
        MINT(40+I)=1
        IF(MINT(123).GE.1.AND.MINT(10+I).EQ.22) MINT(40+I)=2
        IF(IABS(MINT(10+I)).GT.100) MINT(40+I)=2
        MINT(44+I)=MINT(40+I)
        IF(MSTP(11).GE.1.AND.(IABS(MINT(10+I)).EQ.11.OR.
     &  IABS(MINT(10+I)).EQ.13.OR.IABS(MINT(10+I)).EQ.15)) MINT(44+I)=3
  140 CONTINUE
 
C...If two real gammas, whereof one direct, pick the first.
C...For two virtual photons, keep requested order.
      IF(MINT(11).EQ.22.AND.MINT(12).EQ.22) THEN
        IF(MSTP(14).LE.10.AND.MINT(123).GE.4.AND.MINT(123).LE.6) THEN
          MINT(41)=1
          MINT(45)=1
        ELSEIF(MSTP(14).EQ.12.OR.MSTP(14).EQ.13.OR.MSTP(14).EQ.22.OR.
     &  MSTP(14).EQ.26.OR.MSTP(14).EQ.27) THEN
          MINT(41)=1
          MINT(45)=1
        ELSEIF(MSTP(14).EQ.14.OR.MSTP(14).EQ.17.OR.MSTP(14).EQ.23.OR.
     &  MSTP(14).EQ.28.OR.MSTP(14).EQ.29) THEN
          MINT(42)=1
          MINT(46)=1
        ELSEIF((MSTP(14).EQ.20.OR.MSTP(14).EQ.30).AND.(MINT(122).EQ.2
     &  .OR.MINT(122).EQ.3.OR.MINT(122).EQ.10.OR.MINT(122).EQ.11)) THEN
          MINT(41)=1
          MINT(45)=1
        ELSEIF((MSTP(14).EQ.20.OR.MSTP(14).EQ.30).AND.(MINT(122).EQ.4
     &  .OR.MINT(122).EQ.7.OR.MINT(122).EQ.12.OR.MINT(122).EQ.13)) THEN
          MINT(42)=1
          MINT(46)=1
        ELSEIF(MSTP(14).EQ.25.AND.MINT(122).EQ.2) THEN
          MINT(41)=1
          MINT(45)=1
        ELSEIF(MSTP(14).EQ.25.AND.MINT(122).EQ.3) THEN
          MINT(42)=1
          MINT(46)=1
        ENDIF
      ELSEIF(MINT(11).EQ.22.OR.MINT(12).EQ.22) THEN
        IF(MSTP(14).EQ.26.OR.MSTP(14).EQ.28.OR.MINT(122).EQ.4) THEN
          IF(MINT(11).EQ.22) THEN
            MINT(41)=1
            MINT(45)=1
          ELSE
            MINT(42)=1
            MINT(46)=1
          ENDIF
        ENDIF
        IF(MINT(123).GE.4.AND.MINT(123).LE.7) CALL PYERRM(26,
     &  '(PYINPR:) unallowed MSTP(14) code for single photon')
      ENDIF
 
C...Flavour information on combination of incoming particles.
      MINT(43)=2*MINT(41)+MINT(42)-2
      MINT(44)=MINT(43)
      IF(MINT(123).LE.0) THEN
        IF(MINT(11).EQ.22) MINT(43)=MINT(43)+2
        IF(MINT(12).EQ.22) MINT(43)=MINT(43)+1
      ELSEIF(MINT(123).LE.3) THEN
        IF(MINT(11).EQ.22) MINT(44)=MINT(44)-2
        IF(MINT(12).EQ.22) MINT(44)=MINT(44)-1
      ELSEIF(MINT(11).EQ.22.AND.MINT(12).EQ.22) THEN
        MINT(43)=4
        MINT(44)=1
      ENDIF
      MINT(47)=2*MIN(2,MINT(45))+MIN(2,MINT(46))-2
      IF(MIN(MINT(45),MINT(46)).EQ.3) MINT(47)=5
      IF(MINT(45).EQ.1.AND.MINT(46).EQ.3) MINT(47)=6
      IF(MINT(45).EQ.3.AND.MINT(46).EQ.1) MINT(47)=7
      MINT(50)=0
      IF(MINT(41).EQ.2.AND.MINT(42).EQ.2.AND.MINT(111).NE.12) MINT(50)=1
      MINT(107)=0
      MINT(108)=0
      IF(MINT(121).EQ.9.OR.MINT(121).EQ.13) THEN
        IF((MINT(122).GE.4.AND.MINT(122).LE.6).OR.MINT(122).EQ.12)
     &  MINT(107)=2
        IF((MINT(122).GE.7.AND.MINT(122).LE.9).OR.MINT(122).EQ.13)
     &  MINT(107)=3
        IF(MINT(122).EQ.10.OR.MINT(122).EQ.11) MINT(107)=4
        IF(MINT(122).EQ.2.OR.MINT(122).EQ.5.OR.MINT(122).EQ.8.OR.
     &  MINT(122).EQ.10) MINT(108)=2
        IF(MINT(122).EQ.3.OR.MINT(122).EQ.6.OR.MINT(122).EQ.9.OR.
     &  MINT(122).EQ.11) MINT(108)=3
        IF(MINT(122).EQ.12.OR.MINT(122).EQ.13) MINT(108)=4
      ELSEIF(MINT(121).EQ.4.AND.MSTP(14).EQ.25) THEN
        IF(MINT(122).GE.3) MINT(107)=1
        IF(MINT(122).EQ.2.OR.MINT(122).EQ.4) MINT(108)=1
      ELSEIF(MINT(121).EQ.2) THEN
        IF(MINT(122).EQ.2.AND.MINT(11).EQ.22) MINT(107)=1
        IF(MINT(122).EQ.2.AND.MINT(12).EQ.22) MINT(108)=1
      ELSE
        IF(MINT(11).EQ.22) THEN
          MINT(107)=MINT(123)
          IF(MINT(123).GE.4) MINT(107)=0
          IF(MINT(123).EQ.7) MINT(107)=2
          IF(MSTP(14).EQ.26.OR.MSTP(14).EQ.27) MINT(107)=4
          IF(MSTP(14).EQ.28) MINT(107)=2
          IF(MSTP(14).EQ.29) MINT(107)=3
          IF(MSTP(14).EQ.30.AND.MINT(121).EQ.4.AND.MINT(122).EQ.4)
     &    MINT(107)=4
        ENDIF
        IF(MINT(12).EQ.22) THEN
          MINT(108)=MINT(123)
          IF(MINT(123).GE.4) MINT(108)=MINT(123)-3
          IF(MINT(123).EQ.7) MINT(108)=3
          IF(MSTP(14).EQ.26) MINT(108)=2
          IF(MSTP(14).EQ.27) MINT(108)=3
          IF(MSTP(14).EQ.28.OR.MSTP(14).EQ.29) MINT(108)=4
          IF(MSTP(14).EQ.30.AND.MINT(121).EQ.4.AND.MINT(122).EQ.4)
     &    MINT(108)=4
        ENDIF
        IF(MINT(11).EQ.22.AND.MINT(12).EQ.22.AND.(MSTP(14).EQ.14.OR.
     &  MSTP(14).EQ.17.OR.MSTP(14).EQ.18.OR.MSTP(14).EQ.23)) THEN
          MINTTP=MINT(107)
          MINT(107)=MINT(108)
          MINT(108)=MINTTP
        ENDIF
      ENDIF
      IF(MINT(15).EQ.22.AND.MINT(41).EQ.2) MINT(15)=0
      IF(MINT(16).EQ.22.AND.MINT(42).EQ.2) MINT(16)=0
 
C...Select default processes according to incoming beams
C...(already done for gamma-p and gamma-gamma with
C...MSTP(14) = 10, 20, 25 or 30).
      IF(MINT(121).GT.1) THEN
      ELSEIF(MSEL.EQ.1.OR.MSEL.EQ.2) THEN
 
        IF(MINT(43).EQ.1) THEN
C...Lepton + lepton -> gamma/Z0 or W.
          IF(MINT(11)+MINT(12).EQ.0) MSUB(1)=1
          IF(MINT(11)+MINT(12).NE.0) MSUB(2)=1
 
        ELSEIF(MINT(43).LE.3.AND.MINT(123).EQ.0.AND.
     &    (MINT(11).EQ.22.OR.MINT(12).EQ.22)) THEN
C...Unresolved photon + lepton: Compton scattering.
          MSUB(133)=1
          MSUB(134)=1
 
        ELSEIF((MINT(123).EQ.8.OR.MINT(123).EQ.9).AND.(MINT(11).EQ.22
     &  .OR.MINT(12).EQ.22)) THEN
C...DIS as pure gamma* + f -> f process.
          MSUB(99)=1
 
        ELSEIF(MINT(43).LE.3) THEN
C...Lepton + hadron: deep inelastic scattering.
          MSUB(10)=1
 
        ELSEIF(MINT(123).EQ.0.AND.MINT(11).EQ.22.AND.
     &    MINT(12).EQ.22) THEN
C...Two unresolved photons: fermion pair production,
C...exclude lepton pairs.
          DO 150 ISUB=137,140
            MSUB(ISUB)=1
  150     CONTINUE
          DO 160 II=MDCY(22,2),MDCY(22,2)+MDCY(22,3)-1
            IF(IABS(KFDP(II,1)).GE.10) MDME(II,1)=MIN(0,MDME(II,1))
  160     CONTINUE
          PTMDIR=PTMRUN
          IF(MSTP(18).EQ.2) PTMDIR=PARP(15)
          IF(CKIN(3).LT.PTMRUN.OR.MSEL.EQ.2) CKIN(3)=PTMDIR
          CKIN(1)=MAX(CKIN(1),2D0*CKIN(3))
 
        ELSEIF((MINT(123).EQ.0.AND.(MINT(11).EQ.22.OR.MINT(12).EQ.22))
     &    .OR.(MINT(123).GE.4.AND.MINT(123).LE.6.AND.MINT(11).EQ.22.AND.
     &    MINT(12).EQ.22)) THEN
C...Unresolved photon + hadron: photon-parton scattering.
          DO 170 ISUB=131,136
            MSUB(ISUB)=1
  170     CONTINUE
 
        ELSEIF(MSEL.EQ.1) THEN
C...High-pT QCD processes:
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          PTMN=PTMRUN
          VINT(154)=PTMN
          IF(CKIN(3).LT.PTMN) MSUB(95)=1
          IF(MSUB(95).EQ.1.AND.MINT(50).EQ.0) MSUB(95)=0
 
        ELSE
C...All QCD processes:
          MSUB(11)=1
          MSUB(12)=1
          MSUB(13)=1
          MSUB(28)=1
          MSUB(53)=1
          MSUB(68)=1
          MSUB(91)=1
          MSUB(92)=1
          MSUB(93)=1
          MSUB(94)=1
          MSUB(95)=1
        ENDIF
 
      ELSEIF(MSEL.GE.4.AND.MSEL.LE.8) THEN
C...Heavy quark production.
        MSUB(81)=1
        MSUB(82)=1
        MSUB(84)=1
        DO 180 J=1,MIN(8,MDCY(21,3))
          MDME(MDCY(21,2)+J-1,1)=0
  180   CONTINUE
        MDME(MDCY(21,2)+MSEL-1,1)=1
        MSUB(85)=1
        DO 190 J=1,MIN(12,MDCY(22,3))
          MDME(MDCY(22,2)+J-1,1)=0
  190   CONTINUE
        MDME(MDCY(22,2)+MSEL-1,1)=1
 
      ELSEIF(MSEL.EQ.10) THEN
C...Prompt photon production:
        MSUB(14)=1
        MSUB(18)=1
        MSUB(29)=1
 
      ELSEIF(MSEL.EQ.11) THEN
C...Z0/gamma* production:
        MSUB(1)=1
 
      ELSEIF(MSEL.EQ.12) THEN
C...W+/- production:
        MSUB(2)=1
 
      ELSEIF(MSEL.EQ.13) THEN
C...Z0 + jet:
        MSUB(15)=1
        MSUB(30)=1
 
      ELSEIF(MSEL.EQ.14) THEN
C...W+/- + jet:
        MSUB(16)=1
        MSUB(31)=1
 
      ELSEIF(MSEL.EQ.15) THEN
C...Z0 & W+/- pair production:
        MSUB(19)=1
        MSUB(20)=1
        MSUB(22)=1
        MSUB(23)=1
        MSUB(25)=1
 
      ELSEIF(MSEL.EQ.16) THEN
C...h0 production:
        MSUB(3)=1
        MSUB(102)=1
        MSUB(103)=1
        MSUB(123)=1
        MSUB(124)=1
 
      ELSEIF(MSEL.EQ.17) THEN
C...h0 & Z0 or W+/- pair production:
        MSUB(24)=1
        MSUB(26)=1
 
      ELSEIF(MSEL.EQ.18) THEN
C...h0 production; interesting processes in e+e-.
        MSUB(24)=1
        MSUB(103)=1
        MSUB(123)=1
        MSUB(124)=1
 
      ELSEIF(MSEL.EQ.19) THEN
C...h0, H0 and A0 production; interesting processes in e+e-.
        MSUB(24)=1
        MSUB(103)=1
        MSUB(123)=1
        MSUB(124)=1
        MSUB(153)=1
        MSUB(171)=1
        MSUB(173)=1
        MSUB(174)=1
        MSUB(158)=1
        MSUB(176)=1
        MSUB(178)=1
        MSUB(179)=1
 
      ELSEIF(MSEL.EQ.21) THEN
C...Z'0 production:
        MSUB(141)=1
 
      ELSEIF(MSEL.EQ.22) THEN
C...W'+/- production:
        MSUB(142)=1
 
      ELSEIF(MSEL.EQ.23) THEN
C...H+/- production:
        MSUB(143)=1
 
      ELSEIF(MSEL.EQ.24) THEN
C...R production:
        MSUB(144)=1
 
      ELSEIF(MSEL.EQ.25) THEN
C...LQ (leptoquark) production.
        MSUB(145)=1
        MSUB(162)=1
        MSUB(163)=1
        MSUB(164)=1
 
      ELSEIF(MSEL.GE.35.AND.MSEL.LE.38) THEN
C...Production of one heavy quark (W exchange):
        MSUB(83)=1
        DO 200 J=1,MIN(8,MDCY(21,3))
          MDME(MDCY(21,2)+J-1,1)=0
  200   CONTINUE
        MDME(MDCY(21,2)+MSEL-31,1)=1
 
CMRENNA++Define SUSY alternatives.
      ELSEIF(MSEL.EQ.39) THEN
C...Turn on all SUSY processes.
        IF(MINT(43).EQ.4) THEN
C...Hadron-hadron processes.
          DO 210 I=201,296
            IF(ISET(I).GE.0) MSUB(I)=1
  210     CONTINUE
        ELSEIF(MINT(43).EQ.1) THEN
C...Lepton-lepton processes: QED production of squarks.
          DO 220 I=201,214
            MSUB(I)=1
  220     CONTINUE
          MSUB(210)=0
          MSUB(211)=0
          MSUB(212)=0
          DO 230 I=216,228
            MSUB(I)=1
  230     CONTINUE
          DO 240 I=261,263
            MSUB(I)=1
  240     CONTINUE
          MSUB(277)=1
          MSUB(278)=1
        ENDIF
 
      ELSEIF(MSEL.EQ.40) THEN
C...Gluinos and squarks.
        IF(MINT(43).EQ.4) THEN
          MSUB(243)=1
          MSUB(244)=1
          MSUB(258)=1
          MSUB(259)=1
          MSUB(261)=1
          MSUB(262)=1
          MSUB(264)=1
          MSUB(265)=1
          DO 250 I=271,296
            MSUB(I)=1
  250     CONTINUE
        ELSEIF(MINT(43).EQ.1) THEN
          MSUB(277)=1
          MSUB(278)=1
        ENDIF
 
      ELSEIF(MSEL.EQ.41) THEN
C...Stop production.
        MSUB(261)=1
        MSUB(262)=1
        MSUB(263)=1
        IF(MINT(43).EQ.4) THEN
          MSUB(264)=1
          MSUB(265)=1
        ENDIF
 
      ELSEIF(MSEL.EQ.42) THEN
C...Slepton production.
        DO 260 I=201,214
          MSUB(I)=1
  260   CONTINUE
        IF(MINT(43).NE.4) THEN
          MSUB(210)=0
          MSUB(211)=0
          MSUB(212)=0
        ENDIF
 
      ELSEIF(MSEL.EQ.43) THEN
C...Neutralino/Chargino + Gluino/Squark.
        IF(MINT(43).EQ.4) THEN
          DO 270 I=237,242
            MSUB(I)=1
  270     CONTINUE
          DO 280 I=246,254
            MSUB(I)=1
  280     CONTINUE
          MSUB(256)=1
        ENDIF
 
      ELSEIF(MSEL.EQ.44) THEN
C...Neutralino/Chargino pair production.
        IF(MINT(43).EQ.4) THEN
          DO 290 I=216,236
            MSUB(I)=1
  290     CONTINUE
        ELSEIF(MINT(43).EQ.1) THEN
          DO 300 I=216,228
            MSUB(I)=1
  300     CONTINUE
        ENDIF
 
      ELSEIF(MSEL.EQ.45) THEN
C...Sbottom production.
        MSUB(287)=1
        MSUB(288)=1
        IF(MINT(43).EQ.4) THEN
          DO 310 I=281,296
            MSUB(I)=1
  310     CONTINUE
        ENDIF
 
      ELSEIF(MSEL.EQ.50) THEN
C...Pair production of technipions and gauge bosons.
        DO 320 I=361,368
          MSUB(I)=1
  320   CONTINUE
        IF(MINT(43).EQ.4) THEN
          DO 330 I=370,377
            MSUB(I)=1
  330     CONTINUE
        ENDIF
 
      ELSEIF(MSEL.EQ.51) THEN
C...QCD 2 -> 2 processes with compositeness/technicolor modifications.
        DO 340 I=381,386
          MSUB(I)=1
  340   CONTINUE
 
      ELSEIF(MSEL.EQ.61) THEN
C...Charmonium production in colour octet model, with recoiling parton.
        DO 342 I=421,439
          MSUB(I)=1
 342   CONTINUE
 
      ELSEIF(MSEL.EQ.62) THEN
C...Bottomonium production in colour octet model, with recoiling parton.
        DO 344 I=461,479
          MSUB(I)=1
 344   CONTINUE
 
      ELSEIF(MSEL.EQ.63) THEN
C...Charmonium and bottomonium production in colour octet model.
        DO 346 I=421,439
          MSUB(I)=1
          MSUB(I+40)=1
 346   CONTINUE
      ENDIF
 
C...Find heaviest new quark flavour allowed in processes 81-84.
      KFLQM=1
      DO 350 I=1,MIN(8,MDCY(21,3))
        IDC=I+MDCY(21,2)-1
        IF(MDME(IDC,1).LE.0) GOTO 350
        KFLQM=I
  350 CONTINUE
      IF(MSTP(7).GE.1.AND.MSTP(7).LE.8.AND.(MSEL.LE.3.OR.MSEL.GE.9))
     &KFLQM=MSTP(7)
      MINT(55)=KFLQM
      KFPR(81,1)=KFLQM
      KFPR(81,2)=KFLQM
      KFPR(82,1)=KFLQM
      KFPR(82,2)=KFLQM
      KFPR(83,1)=KFLQM
      KFPR(84,1)=KFLQM
      KFPR(84,2)=KFLQM
 
C...Find heaviest new fermion flavour allowed in process 85.
      KFLFM=1
      DO 360 I=1,MIN(12,MDCY(22,3))
        IDC=I+MDCY(22,2)-1
        IF(MDME(IDC,1).LE.0) GOTO 360
        KFLFM=KFDP(IDC,1)
  360 CONTINUE
      IF(((MSTP(7).GE.1.AND.MSTP(7).LE.8).OR.(MSTP(7).GE.11.AND.
     &MSTP(7).LE.18)).AND.(MSEL.LE.3.OR.MSEL.GE.9)) KFLFM=MSTP(7)
      MINT(56)=KFLFM
      KFPR(85,1)=KFLFM
      KFPR(85,2)=KFLFM

C...Initialize Generic Processes
      KFGEN=9900001
      KCGEN=PYCOMP(KFGEN)
      IF(KCGEN.GT.0) THEN
        IDCY=MDCY(KCGEN,2)
        IF(IDCY.GT.0) THEN
          KFF1=KFDP(IDCY+1,1)
          KFF2=KFDP(IDCY+1,2)
          KCF1=PYCOMP(KFF1)
          KCF2=PYCOMP(KFF2)
          JCOL1=IABS(KCHG(KCF1,2))
          IF(JCOL1.EQ.1) THEN
            KF1=KFF1
            KF2=KFF2
          ELSE
            KF1=KFF2
            KF2=KFF1
          ENDIF
          KFPR(481,1)=KF1
          KFPR(481,2)=KF2
          KFPR(482,1)=KF1
          KFPR(482,2)=KF2
        ENDIF
        IF(KFDP(IDCY,1).EQ.21.OR.KFDP(IDCY,2).EQ.21) THEN
          KFIN(1,0)=1
          KFIN(2,0)=1
        ENDIF
      ENDIF
 
C...Import relevant information on external user processes.
      IF(MINT(111).GE.11) THEN
        IPYPR=0
        DO 390 IUP=1,NPRUP
C...Find next empty PYTHIA process number slot and enable it.
  370     IPYPR=IPYPR+1
          IF(IPYPR.GT.500) CALL PYERRM(26,
     &    '(PYINPR.) no more empty slots for user processes')
          IF(ISET(IPYPR).GE.0.AND.ISET(IPYPR).LE.9) GOTO 370
          IF(IPYPR.GE.91.AND.IPYPR.LE.100) GOTO 370
          ISET(IPYPR)=11
C...Overwrite KFPR with references back to process number and ID.
          KFPR(IPYPR,1)=IUP
          KFPR(IPYPR,2)=LPRUP(IUP)
C...Process title.
          WRITE(CHIPR,'(I10)') LPRUP(IUP)
          ICHIN=1
          DO 380 ICH=1,9
            IF(CHIPR(ICH:ICH).EQ.' ') ICHIN=ICH+1
  380     CONTINUE
          PROC(IPYPR)='User process '//CHIPR(ICHIN:10)//' '
C...Switch on process.
          MSUB(IPYPR)=1
  390   CONTINUE
      ENDIF

      RETURN
      END
