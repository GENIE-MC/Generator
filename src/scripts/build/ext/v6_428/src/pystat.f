 
 
C***********************************************************************
 
C...PYSTAT
C...Prints out information about cross-sections, decay widths, branching
C...ratios, kinematical limits, status codes and parameter values.
 
      SUBROUTINE PYSTAT(MSTAT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
      PARAMETER (EPS=1D-3)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT6/PROC(0:500)
      CHARACTER PROC*28, CHTMP*16
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYMSRV/RVLAM(3,3,3), RVLAMP(3,3,3), RVLAMB(3,3,3)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,/PYINT1/,
     &/PYINT2/,/PYINT4/,/PYINT5/,/PYINT6/,/PYMSSM/,/PYMSRV/
C...Local arrays, character variables and data.
      DIMENSION WDTP(0:400),WDTE(0:400,0:5),NMODES(0:20),PBRAT(10)
      CHARACTER PROGA(6)*28,CHAU*16,CHKF*16,CHD1*16,CHD2*16,CHD3*16,
     &CHIN(2)*12,STATE(-1:5)*4,CHKIN(21)*18,DISGA(2)*28,
     &PROGG9(13)*28,PROGG4(4)*28,PROGG2(2)*28,PROGP4(4)*28
      CHARACTER*24 CHD0, CHDC(10)
      CHARACTER*6 DNAME(3)
      DATA PROGA/
     &'VMD/hadron * VMD            ','VMD/hadron * direct         ',
     &'VMD/hadron * anomalous      ','direct * direct             ',
     &'direct * anomalous          ','anomalous * anomalous       '/
      DATA DISGA/'e * VMD','e * anomalous'/
      DATA PROGG9/
     &'direct * direct             ','direct * VMD                ',
     &'direct * anomalous          ','VMD * direct                ',
     &'VMD * VMD                   ','VMD * anomalous             ',
     &'anomalous * direct          ','anomalous * VMD             ',
     &'anomalous * anomalous       ','DIS * VMD                   ',
     &'DIS * anomalous             ','VMD * DIS                   ',
     &'anomalous * DIS             '/
      DATA PROGG4/
     &'direct * direct             ','direct * resolved           ',
     &'resolved * direct           ','resolved * resolved         '/
      DATA PROGG2/
     &'direct * hadron             ','resolved * hadron           '/
      DATA PROGP4/
     &'VMD * hadron                ','direct * hadron             ',
     &'anomalous * hadron          ','DIS * hadron                '/
      DATA STATE/'----','off ','on  ','on/+','on/-','on/1','on/2'/,
     &CHKIN/' m_hard (GeV/c^2) ',' p_T_hard (GeV/c) ',
     &'m_finite (GeV/c^2)','   y*_subsystem   ','     y*_large     ',
     &'     y*_small     ','    eta*_large    ','    eta*_small    ',
     &'cos(theta*)_large ','cos(theta*)_small ','       x_1        ',
     &'       x_2        ','       x_F        ',' cos(theta_hard)  ',
     &'m''_hard (GeV/c^2) ','       tau        ','        y*        ',
     &'cos(theta_hard^-) ','cos(theta_hard^+) ','      x_T^2       ',
     &'       tau''       '/
      DATA DNAME /'q     ','lepton','nu    '/
 
C...Cross-sections.
      IF(MSTAT.LE.1) THEN
        IF(MINT(121).GT.1) CALL PYSAVE(5,0)
        WRITE(MSTU(11),5000)
        WRITE(MSTU(11),5100)
        WRITE(MSTU(11),5200) 0,PROC(0),NGEN(0,3),NGEN(0,1),XSEC(0,3)
        DO 100 I=1,500
          IF(MSUB(I).NE.1) GOTO 100
          WRITE(MSTU(11),5200) I,PROC(I),NGEN(I,3),NGEN(I,1),XSEC(I,3)
  100   CONTINUE
        IF(MINT(121).GT.1) THEN
          WRITE(MSTU(11),5300)
          DO 110 IGA=1,MINT(121)
            CALL PYSAVE(3,IGA)
            IF(MINT(121).EQ.2.AND.MSTP(14).EQ.10) THEN
              WRITE(MSTU(11),5200) IGA,DISGA(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ELSEIF(MINT(121).EQ.9.OR.MINT(121).EQ.13) THEN
              WRITE(MSTU(11),5200) IGA,PROGG9(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ELSEIF(MINT(121).EQ.4.AND.MSTP(14).EQ.30) THEN
              WRITE(MSTU(11),5200) IGA,PROGP4(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ELSEIF(MINT(121).EQ.4) THEN
              WRITE(MSTU(11),5200) IGA,PROGG4(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ELSEIF(MINT(121).EQ.2) THEN
              WRITE(MSTU(11),5200) IGA,PROGG2(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ELSE
              WRITE(MSTU(11),5200) IGA,PROGA(IGA),NGEN(0,3),NGEN(0,1),
     &        XSEC(0,3)
            ENDIF
  110     CONTINUE
          CALL PYSAVE(5,0)
        ENDIF
        WRITE(MSTU(11),5400) MSTU(23),MSTU(30),MSTU(27),
     &  1D0-DBLE(NGEN(0,3))/MAX(1D0,DBLE(NGEN(0,2)))
 
C...Decay widths and branching ratios.
      ELSEIF(MSTAT.EQ.2) THEN
        WRITE(MSTU(11),5500)
        WRITE(MSTU(11),5600)
        DO 140 KC=1,500
          KF=KCHG(KC,4)
          CALL PYNAME(KF,CHKF)
          IOFF=0
          IF(KC.LE.22) THEN
            IF(KC.GT.2*MSTP(1).AND.KC.LE.10) GOTO 140
            IF(KC.GT.10+2*MSTP(1).AND.KC.LE.20) GOTO 140
            IF(KC.LE.5.OR.(KC.GE.11.AND.KC.LE.16)) IOFF=1
            IF(KC.EQ.18.AND.PMAS(18,1).LT.1D0) IOFF=1
            IF(KC.EQ.21.OR.KC.EQ.22) IOFF=1
          ELSE
            IF(MWID(KC).LE.0) GOTO 140
            IF(IMSS(1).LE.0.AND.(KF/KSUSY1.EQ.1.OR.
     &      KF/KSUSY1.EQ.2)) GOTO 140
          ENDIF
C...Off-shell branchings.
          IF(IOFF.EQ.1) THEN
            NGP=0
            IF(KC.LE.20) NGP=(MOD(KC,10)+1)/2
            IF(NGP.LE.MSTP(1)) WRITE(MSTU(11),5700) KF,CHKF(1:10),
     &      PMAS(KC,1),0D0,0D0,STATE(MDCY(KC,1)),0D0
            DO 120 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              NGP1=0
              IF(IABS(KFDP(IDC,1)).LE.20) NGP1=
     &        (MOD(IABS(KFDP(IDC,1)),10)+1)/2
              NGP2=0
              IF(IABS(KFDP(IDC,2)).LE.20) NGP2=
     &        (MOD(IABS(KFDP(IDC,2)),10)+1)/2
              CALL PYNAME(KFDP(IDC,1),CHD1)
              CALL PYNAME(KFDP(IDC,2),CHD2)
              IF(KFDP(IDC,3).EQ.0) THEN
                IF(MDME(IDC,2).EQ.102.AND.NGP1.LE.MSTP(1).AND.
     &          NGP2.LE.MSTP(1)) WRITE(MSTU(11),5800) IDC,CHD1(1:10),
     &          CHD2(1:10),0D0,0D0,STATE(MDME(IDC,1)),0D0
              ELSE
                CALL PYNAME(KFDP(IDC,3),CHD3)
                IF(MDME(IDC,2).EQ.102.AND.NGP1.LE.MSTP(1).AND.
     &          NGP2.LE.MSTP(1)) WRITE(MSTU(11),5900) IDC,CHD1(1:10),
     &          CHD2(1:10),CHD3(1:10),0D0,0D0,STATE(MDME(IDC,1)),0D0
              ENDIF
  120       CONTINUE
C...On-shell decays.
          ELSE
            CALL PYWIDT(KF,PMAS(KC,1)**2,WDTP,WDTE)
            BRFIN=1D0
            IF(WDTE(0,0).LE.0D0) BRFIN=0D0
            WRITE(MSTU(11),5700) KF,CHKF(1:10),PMAS(KC,1),WDTP(0),1D0,
     &      STATE(MDCY(KC,1)),BRFIN
            DO 130 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              NGP1=0
              IF(IABS(KFDP(IDC,1)).LE.20) NGP1=
     &        (MOD(IABS(KFDP(IDC,1)),10)+1)/2
              NGP2=0
              IF(IABS(KFDP(IDC,2)).LE.20) NGP2=
     &        (MOD(IABS(KFDP(IDC,2)),10)+1)/2
              BRPRI=0D0
              IF(WDTP(0).GT.0D0) BRPRI=WDTP(J)/WDTP(0)
              BRFIN=0D0
              IF(WDTE(0,0).GT.0D0) BRFIN=WDTE(J,0)/WDTE(0,0)
              CALL PYNAME(KFDP(IDC,1),CHD1)
              CALL PYNAME(KFDP(IDC,2),CHD2)
              IF(KFDP(IDC,3).EQ.0) THEN
                IF(NGP1.LE.MSTP(1).AND.NGP2.LE.MSTP(1))
     &          WRITE(MSTU(11),5800) IDC,CHD1(1:10),
     &          CHD2(1:10),WDTP(J),BRPRI,
     &          STATE(MDME(IDC,1)),BRFIN
              ELSE
                CALL PYNAME(KFDP(IDC,3),CHD3)
                IF(NGP1.LE.MSTP(1).AND.NGP2.LE.MSTP(1))
     &          WRITE(MSTU(11),5900) IDC,CHD1(1:10),
     &          CHD2(1:10),CHD3(1:10),WDTP(J),BRPRI,
     &          STATE(MDME(IDC,1)),BRFIN
              ENDIF
  130       CONTINUE
          ENDIF
  140   CONTINUE
        WRITE(MSTU(11),6000)
 
C...Allowed incoming partons/particles at hard interaction.
      ELSEIF(MSTAT.EQ.3) THEN
        WRITE(MSTU(11),6100)
        CALL PYNAME(MINT(11),CHAU)
        CHIN(1)=CHAU(1:12)
        CALL PYNAME(MINT(12),CHAU)
        CHIN(2)=CHAU(1:12)
        WRITE(MSTU(11),6200) CHIN(1),CHIN(2)
        DO 150 I=-20,22
          IF(I.EQ.0) GOTO 150
          IA=IABS(I)
          IF(IA.GT.MSTP(58).AND.IA.LE.10) GOTO 150
          IF(IA.GT.10+2*MSTP(1).AND.IA.LE.20) GOTO 150
          CALL PYNAME(I,CHAU)
          WRITE(MSTU(11),6300) CHAU,STATE(KFIN(1,I)),CHAU,
     &    STATE(KFIN(2,I))
  150   CONTINUE
        WRITE(MSTU(11),6400)
 
C...User-defined limits on kinematical variables.
      ELSEIF(MSTAT.EQ.4) THEN
        WRITE(MSTU(11),6500)
        WRITE(MSTU(11),6600)
        SHRMAX=CKIN(2)
        IF(SHRMAX.LT.0D0) SHRMAX=VINT(1)
        WRITE(MSTU(11),6700) CKIN(1),CHKIN(1),SHRMAX
        PTHMIN=MAX(CKIN(3),CKIN(5))
        PTHMAX=CKIN(4)
        IF(PTHMAX.LT.0D0) PTHMAX=0.5D0*SHRMAX
        WRITE(MSTU(11),6800) CKIN(3),PTHMIN,CHKIN(2),PTHMAX
        WRITE(MSTU(11),6900) CHKIN(3),CKIN(6)
        DO 160 I=4,14
          WRITE(MSTU(11),6700) CKIN(2*I-1),CHKIN(I),CKIN(2*I)
  160   CONTINUE
        SPRMAX=CKIN(32)
        IF(SPRMAX.LT.0D0) SPRMAX=VINT(1)
        WRITE(MSTU(11),6700) CKIN(31),CHKIN(15),SPRMAX
        WRITE(MSTU(11),7000)
 
C...Status codes and parameter values.
      ELSEIF(MSTAT.EQ.5) THEN
        WRITE(MSTU(11),7100)
        WRITE(MSTU(11),7200)
        DO 170 I=1,100
          WRITE(MSTU(11),7300) I,MSTP(I),PARP(I),100+I,MSTP(100+I),
     &    PARP(100+I)
  170   CONTINUE
 
C...List of all processes implemented in the program.
      ELSEIF(MSTAT.EQ.6) THEN
        WRITE(MSTU(11),7400)
        WRITE(MSTU(11),7500)
        DO 180 I=1,500
          IF(ISET(I).LT.0) GOTO 180
          WRITE(MSTU(11),7600) I,PROC(I),ISET(I),KFPR(I,1),KFPR(I,2)
  180   CONTINUE
        WRITE(MSTU(11),7700)
 
      ELSEIF(MSTAT.EQ.7) THEN
      WRITE (MSTU(11),8000)
      NMODES(0)=0
      NMODES(10)=0
      NMODES(9)=0
      DO 290 ILR=1,2
        DO 280 KFSM=1,16
          KFSUSY=ILR*KSUSY1+KFSM
          NRVDC=0
C...SDOWN DECAYS
          IF (KFSM.EQ.1.OR.KFSM.EQ.3.OR.KFSM.EQ.5) THEN
            NRVDC=3
            DO 190 I=1,NRVDC
              PBRAT(I)=0D0
              NMODES(I)=0
  190       CONTINUE
            CALL PYNAME(KFSUSY,CHTMP)
            CHD0=CHTMP//' '
            CHDC(1)=DNAME(3) // ' + ' // DNAME(1)
            CHDC(2)=DNAME(2) // ' + ' // DNAME(1)
            CHDC(3)=DNAME(1) // ' + ' // DNAME(1)
            KC=PYCOMP(KFSUSY)
            DO 200 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              ID1=IABS(KFDP(IDC,1))
              ID2=IABS(KFDP(IDC,2))
              IF (KFDP(IDC,3).EQ.0) THEN
                IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND.(ID2
     &               .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  PBRAT(1)=PBRAT(1)+BRAT(IDC)
                  NMODES(1)=NMODES(1)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &                 .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6)) THEN
                  PBRAT(2)=PBRAT(2)+BRAT(IDC)
                  NMODES(2)=NMODES(2)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ELSE IF ((ID1.EQ.2.OR.ID1.EQ.4.OR.ID1.EQ.6).AND
     &                 .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  PBRAT(3)=PBRAT(3)+BRAT(IDC)
                  NMODES(3)=NMODES(3)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
              ENDIF
  200       CONTINUE
          ENDIF
C...SUP DECAYS
          IF (KFSM.EQ.2.OR.KFSM.EQ.4.OR.KFSM.EQ.6) THEN
            NRVDC=2
            DO 210 I=1,NRVDC
              NMODES(I)=0
              PBRAT(I)=0D0
  210       CONTINUE
            CALL PYNAME(KFSUSY,CHTMP)
            CHD0=CHTMP//' '
            CHDC(1)=DNAME(2) // ' + ' // DNAME(1)
            CHDC(2)=DNAME(1) // ' + ' // DNAME(1)
            KC=PYCOMP(KFSUSY)
            DO 220 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              ID1=IABS(KFDP(IDC,1))
              ID2=IABS(KFDP(IDC,2))
              IF (KFDP(IDC,3).EQ.0) THEN
                IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND.(ID2
     &               .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  PBRAT(1)=PBRAT(1)+BRAT(IDC)
                  NMODES(1)=NMODES(1)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ELSE IF ((ID1.EQ.1.OR.ID1.EQ.3.OR.ID1.EQ.5).AND.(ID2
     &               .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  PBRAT(2)=PBRAT(2)+BRAT(IDC)
                  NMODES(2)=NMODES(2)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
              ENDIF
  220       CONTINUE
          ENDIF
C...SLEPTON DECAYS
          IF (KFSM.EQ.11.OR.KFSM.EQ.13.OR.KFSM.EQ.15) THEN
            NRVDC=2
            DO 230 I=1,NRVDC
              PBRAT(I)=0D0
              NMODES(I)=0
  230       CONTINUE
            CALL PYNAME(KFSUSY,CHTMP)
            CHD0=CHTMP//' '
            CHDC(1)=DNAME(3) // ' + ' // DNAME(2)
            CHDC(2)=DNAME(1) // ' + ' // DNAME(1)
            KC=PYCOMP(KFSUSY)
            DO 240 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              ID1=IABS(KFDP(IDC,1))
              ID2=IABS(KFDP(IDC,2))
              IF (KFDP(IDC,3).EQ.0) THEN
                IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND.(ID2
     &               .EQ.11.OR.ID2.EQ.13.OR.ID2.EQ.15)) THEN
                  PBRAT(1)=PBRAT(1)+BRAT(IDC)
                  NMODES(1)=NMODES(1)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
                IF ((ID1.EQ.2.OR.ID1.EQ.4.OR.ID1.EQ.6).AND.(ID2
     &               .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  PBRAT(2)=PBRAT(2)+BRAT(IDC)
                  NMODES(2)=NMODES(2)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
              ENDIF
  240       CONTINUE
          ENDIF
C...SNEUTRINO DECAYS
          IF ((KFSM.EQ.12.OR.KFSM.EQ.14.OR.KFSM.EQ.16).AND.ILR.EQ.1)
     &         THEN
            NRVDC=2
            DO 250 I=1,NRVDC
              PBRAT(I)=0D0
              NMODES(I)=0
  250       CONTINUE
            CALL PYNAME(KFSUSY,CHTMP)
            CHD0=CHTMP//' '
            CHDC(1)=DNAME(2) // ' + ' // DNAME(2)
            CHDC(2)=DNAME(1) // ' + ' // DNAME(1)
            KC=PYCOMP(KFSUSY)
            DO 260 J=1,MDCY(KC,3)
              IDC=J+MDCY(KC,2)-1
              ID1=IABS(KFDP(IDC,1))
              ID2=IABS(KFDP(IDC,2))
              IF (KFDP(IDC,3).EQ.0) THEN
                IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND.(ID2
     &               .EQ.11.OR.ID2.EQ.13.OR.ID2.EQ.15)) THEN
                  PBRAT(1)=PBRAT(1)+BRAT(IDC)
                  NMODES(1)=NMODES(1)+1
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
                IF ((ID1.EQ.1.OR.ID1.EQ.3.OR.ID1.EQ.5).AND.(ID2
     &               .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5)) THEN
                  NMODES(2)=NMODES(2)+1
                  PBRAT(2)=PBRAT(2)+BRAT(IDC)
                  IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
                  IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
                ENDIF
              ENDIF
  260       CONTINUE
          ENDIF
          IF (NRVDC.NE.0) THEN
            DO 270 I=1,NRVDC
              WRITE (MSTU(11),8200) CHD0, CHDC(I), PBRAT(I), NMODES(I)
              NMODES(0)=NMODES(0)+NMODES(I)
  270       CONTINUE
          ENDIF
  280   CONTINUE
  290 CONTINUE
      DO 370 KFSM=21,37
        KFSUSY=KSUSY1+KFSM
        NRVDC=0
C...NEUTRALINO DECAYS
        IF (KFSM.EQ.22.OR.KFSM.EQ.23.OR.KFSM.EQ.25.OR.KFSM.EQ.35) THEN
          NRVDC=4
          DO 300 I=1,NRVDC
            PBRAT(I)=0D0
            NMODES(I)=0
  300     CONTINUE
          CALL PYNAME(KFSUSY,CHTMP)
          CHD0=CHTMP//' '
          CHDC(1)=DNAME(3) // ' + ' // DNAME(2) // ' + ' // DNAME(2)
          CHDC(2)=DNAME(3) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(3)=DNAME(2) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(4)=DNAME(1) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          KC=PYCOMP(KFSUSY)
          DO 310 J=1,MDCY(KC,3)
            IDC=J+MDCY(KC,2)-1
            ID1=IABS(KFDP(IDC,1))
            ID2=IABS(KFDP(IDC,2))
            ID3=IABS(KFDP(IDC,3))
            IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND.(ID2
     &           .EQ.11.OR.ID2.EQ.13.OR.ID2.EQ.15).AND.(ID3.EQ.11.OR
     &           .ID3.EQ.13.OR.ID3.EQ.15)) THEN
              PBRAT(1)=PBRAT(1)+BRAT(IDC)
              NMODES(1)=NMODES(1)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ.1
     &             .OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(2)=PBRAT(2)+BRAT(IDC)
              NMODES(2)=NMODES(2)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &             .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6).AND.(ID3.EQ.1
     &             .OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(3)=PBRAT(3)+BRAT(IDC)
              NMODES(3)=NMODES(3)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.2.OR.ID1.EQ.4.OR.ID1.EQ.6).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ.1
     &             .OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(4)=PBRAT(4)+BRAT(IDC)
              NMODES(4)=NMODES(4)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ENDIF
  310     CONTINUE
        ENDIF
C...CHARGINO DECAYS
        IF (KFSM.EQ.24.OR.KFSM.EQ.37) THEN
          NRVDC=5
          DO 320 I=1,NRVDC
            PBRAT(I)=0D0
            NMODES(I)=0
  320     CONTINUE
          CALL PYNAME(KFSUSY,CHTMP)
          CHD0=CHTMP//' '
          CHDC(1)=DNAME(3) // ' + ' // DNAME(3) // ' + ' // DNAME(2)
          CHDC(2)=DNAME(2) // ' + ' // DNAME(2) // ' + ' // DNAME(2)
          CHDC(3)=DNAME(3) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(4)=DNAME(2) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(5)=DNAME(1) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          KC=PYCOMP(KFSUSY)
          DO 330 J=1,MDCY(KC,3)
            IDC=J+MDCY(KC,2)-1
            ID1=IABS(KFDP(IDC,1))
            ID2=IABS(KFDP(IDC,2))
            ID3=IABS(KFDP(IDC,3))
            IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND.(ID2
     &           .EQ.11.OR.ID2.EQ.13.OR.ID2.EQ.15).AND.(ID3.EQ.12.OR
     &           .ID3.EQ.14.OR.ID3.EQ.16)) THEN
              PBRAT(1)=PBRAT(1)+BRAT(IDC)
              NMODES(1)=NMODES(1)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND
     &             .(ID2.EQ.12.OR.ID2.EQ.14.OR.ID2.EQ.16).AND.(ID3.EQ
     &             .11.OR.ID3.EQ.13.OR.ID3.EQ.15)) THEN
              PBRAT(1)=PBRAT(1)+BRAT(IDC)
              NMODES(1)=NMODES(1)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &             .(ID2.EQ.11.OR.ID2.EQ.13.OR.ID2.EQ.15).AND.(ID3.EQ
     &             .11.OR.ID3.EQ.13.OR.ID3.EQ.15)) THEN
              PBRAT(2)=PBRAT(2)+BRAT(IDC)
              NMODES(2)=NMODES(2)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ
     &             .2.OR.ID3.EQ.4.OR.ID3.EQ.6)) THEN
              PBRAT(3)=PBRAT(3)+BRAT(IDC)
              NMODES(3)=NMODES(3)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND
     &             .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6).AND.(ID3.EQ
     &             .1.OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(3)=PBRAT(3)+BRAT(IDC)
              NMODES(3)=NMODES(3)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &             .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6).AND.(ID3.EQ
     &             .2.OR.ID3.EQ.4.OR.ID3.EQ.6)) THEN
              PBRAT(4)=PBRAT(4)+BRAT(IDC)
              NMODES(4)=NMODES(4)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ
     &             .1.OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(4)=PBRAT(4)+BRAT(IDC)
              NMODES(4)=NMODES(4)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.2.OR.ID1.EQ.4.OR.ID1.EQ.6).AND
     &             .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6).AND.(ID3.EQ
     &             .1.OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(5)=PBRAT(5)+BRAT(IDC)
              NMODES(5)=NMODES(5)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.1.OR.ID1.EQ.3.OR.ID1.EQ.5).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ
     &             .1.OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(5)=PBRAT(5)+BRAT(IDC)
              NMODES(5)=NMODES(5)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ENDIF
  330     CONTINUE
        ENDIF
C...GLUINO DECAYS
        IF (KFSM.EQ.21) THEN
          NRVDC=3
          DO 340 I=1,NRVDC
            PBRAT(I)=0D0
            NMODES(I)=0
  340     CONTINUE
          CALL PYNAME(KFSUSY,CHTMP)
          CHD0=CHTMP//' '
          CHDC(1)=DNAME(3) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(2)=DNAME(2) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          CHDC(3)=DNAME(1) // ' + ' // DNAME(1) // ' + ' // DNAME(1)
          KC=PYCOMP(KFSUSY)
          DO 350 J=1,MDCY(KC,3)
            IDC=J+MDCY(KC,2)-1
            ID1=IABS(KFDP(IDC,1))
            ID2=IABS(KFDP(IDC,2))
            ID3=IABS(KFDP(IDC,3))
            IF ((ID1.EQ.12.OR.ID1.EQ.14.OR.ID1.EQ.16).AND.(ID2
     &           .EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ.1.OR
     &           .ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(1)=PBRAT(1)+BRAT(IDC)
              NMODES(1)=NMODES(1)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.11.OR.ID1.EQ.13.OR.ID1.EQ.15).AND
     &             .(ID2.EQ.2.OR.ID2.EQ.4.OR.ID2.EQ.6).AND.(ID3.EQ.1
     &             .OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(2)=PBRAT(2)+BRAT(IDC)
              NMODES(2)=NMODES(2)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ELSE IF ((ID1.EQ.2.OR.ID1.EQ.4.OR.ID1.EQ.6).AND
     &             .(ID2.EQ.1.OR.ID2.EQ.3.OR.ID2.EQ.5).AND.(ID3.EQ.1
     &             .OR.ID3.EQ.3.OR.ID3.EQ.5)) THEN
              PBRAT(3)=PBRAT(3)+BRAT(IDC)
              NMODES(3)=NMODES(3)+1
              IF (BRAT(IDC).GT.0D0) NMODES(10)=NMODES(10)+1
              IF (BRAT(IDC).GT.EPS) NMODES(9)=NMODES(9)+1
            ENDIF
  350     CONTINUE
        ENDIF
 
        IF (NRVDC.NE.0) THEN
          DO 360 I=1,NRVDC
            WRITE (MSTU(11),8200) CHD0, CHDC(I), PBRAT(I), NMODES(I)
            NMODES(0)=NMODES(0)+NMODES(I)
  360     CONTINUE
        ENDIF
  370 CONTINUE
      WRITE (MSTU(11),8100) NMODES(0), NMODES(10), NMODES(9)
 
      IF (IMSS(51).GE.1.OR.IMSS(52).GE.1.OR.IMSS(53).GE.1) THEN
        WRITE (MSTU(11),8500)
        DO 400 IRV=1,3
          DO 390 JRV=1,3
            DO 380 KRV=1,3
              WRITE (MSTU(11),8700) IRV,JRV,KRV,RVLAM(IRV,JRV,KRV)
     &             ,RVLAMP(IRV,JRV,KRV),RVLAMB(IRV,JRV,KRV)
  380       CONTINUE
  390     CONTINUE
  400   CONTINUE
        WRITE (MSTU(11),8600)
      ENDIF
      ENDIF
 
C...Formats for printouts.
 5000 FORMAT('1',9('*'),1X,'PYSTAT:  Statistics on Number of ',
     &'Events and Cross-sections',1X,9('*'))
 5100 FORMAT(/1X,78('=')/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',12X,
     &'Subprocess',12X,'I',6X,'Number of points',6X,'I',4X,'Sigma',3X,
     &'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',34('-'),'I',28('-'),
     &'I',4X,'(mb)',4X,'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',1X,
     &'N:o',1X,'Type',25X,'I',4X,'Generated',9X,'Tried',1X,'I',12X,
     &'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')/1X,'I',34X,'I',28X,
     &'I',12X,'I')
 5200 FORMAT(1X,'I',1X,I3,1X,A28,1X,'I',1X,I12,1X,I13,1X,'I',1X,1P,
     &D10.3,1X,'I')
 5300 FORMAT(1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')/
     &1X,'I',34X,'I',28X,'I',12X,'I')
 5400 FORMAT(1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')//
     &1X,'********* Total number of errors, excluding junctions =',
     &1X,I8,' *************'/
     &1X,'********* Total number of errors, including junctions =',
     &1X,I8,' *************'/
     &1X,'********* Total number of warnings =                   ',
     &1X,I8,' *************'/
     &1X,'********* Fraction of events that fail fragmentation ',
     &'cuts =',1X,F8.5,' *********'/)
 5500 FORMAT('1',27('*'),1X,'PYSTAT:  Decay Widths and Branching ',
     &'Ratios',1X,27('*'))
 5600 FORMAT(/1X,98('=')/1X,'I',49X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/
     &1X,'I',5X,'Mother  -->  Branching/Decay Channel',8X,'I',1X,
     &'Width (GeV)',1X,'I',7X,'B.R.',1X,'I',1X,'Stat',1X,'I',2X,
     &'Eff. B.R.',1X,'I'/1X,'I',49X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/
     &1X,98('='))
 5700 FORMAT(1X,'I',49X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,'I',1X,
     &I8,2X,A10,3X,'(m =',F10.3,')',2X,'-->',5X,'I',2X,1P,D10.3,0P,1X,
     &'I',1X,1P,D10.3,0P,1X,'I',1X,A4,1X,'I',1X,1P,D10.3,0P,1X,'I')
 5800 FORMAT(1X,'I',1X,I8,2X,A10,1X,'+',1X,A10,15X,'I',2X,
     &1P,D10.3,0P,1X,'I',1X,1P,D10.3,0P,1X,'I',1X,A4,1X,'I',1X,
     &1P,D10.3,0P,1X,'I')
 5900 FORMAT(1X,'I',1X,I8,2X,A10,1X,'+',1X,A10,1X,'+',1X,A10,2X,'I',2X,
     &1P,D10.3,0P,1X,'I',1X,1P,D10.3,0P,1X,'I',1X,A4,1X,'I',1X,
     &1P,D10.3,0P,1X,'I')
 6000 FORMAT(1X,'I',49X,'I',13X,'I',12X,'I',6X,'I',12X,'I'/1X,98('='))
 6100 FORMAT('1',7('*'),1X,'PYSTAT: Allowed Incoming Partons/',
     &'Particles at Hard Interaction',1X,7('*'))
 6200 FORMAT(/1X,78('=')/1X,'I',38X,'I',37X,'I'/1X,'I',1X,
     &'Beam particle:',1X,A12,10X,'I',1X,'Target particle:',1X,A12,7X,
     &'I'/1X,'I',38X,'I',37X,'I'/1X,'I',1X,'Content',6X,'State',19X,
     &'I',1X,'Content',6X,'State',18X,'I'/1X,'I',38X,'I',37X,'I'/1X,
     &78('=')/1X,'I',38X,'I',37X,'I')
 6300 FORMAT(1X,'I',1X,A9,5X,A4,19X,'I',1X,A9,5X,A4,18X,'I')
 6400 FORMAT(1X,'I',38X,'I',37X,'I'/1X,78('='))
 6500 FORMAT('1',12('*'),1X,'PYSTAT: User-Defined Limits on ',
     &'Kinematical Variables',1X,12('*'))
 6600 FORMAT(/1X,78('=')/1X,'I',76X,'I')
 6700 FORMAT(1X,'I',16X,1P,D10.3,0P,1X,'<',1X,A,1X,'<',1X,1P,D10.3,0P,
     &16X,'I')
 6800 FORMAT(1X,'I',3X,1P,D10.3,0P,1X,'(',1P,D10.3,0P,')',1X,'<',1X,A,
     &1X,'<',1X,1P,D10.3,0P,16X,'I')
 6900 FORMAT(1X,'I',29X,A,1X,'=',1X,1P,D10.3,0P,16X,'I')
 7000 FORMAT(1X,'I',76X,'I'/1X,78('='))
 7100 FORMAT('1',12('*'),1X,'PYSTAT: Summary of Status Codes and ',
     &'Parameter Values',1X,12('*'))
 7200 FORMAT(/3X,'I',4X,'MSTP(I)',9X,'PARP(I)',20X,'I',4X,'MSTP(I)',9X,
     &'PARP(I)'/)
 7300 FORMAT(1X,I3,5X,I6,6X,1P,D10.3,0P,18X,I3,5X,I6,6X,1P,D10.3)
 7400 FORMAT('1',13('*'),1X,'PYSTAT: List of implemented processes',
     &1X,13('*'))
 7500 FORMAT(/1X,65('=')/1X,'I',34X,'I',28X,'I'/1X,'I',12X,
     &'Subprocess',12X,'I',1X,'ISET',2X,'KFPR(I,1)',2X,'KFPR(I,2)',1X,
     &'I'/1X,'I',34X,'I',28X,'I'/1X,65('=')/1X,'I',34X,'I',28X,'I')
 7600 FORMAT(1X,'I',1X,I3,1X,A28,1X,'I',1X,I4,1X,I10,1X,I10,1X,'I')
 7700 FORMAT(1X,'I',34X,'I',28X,'I'/1X,65('='))
 8000 FORMAT(1X/ 1X/
     &     17X,'Sums over R-Violating branching ratios',1X/ 1X
     &     /1X,70('=')/1X,'I',50X,'I',11X,'I',5X,'I'/1X,'I',4X
     &     ,'Mother  -->  Sum over final state flavours',4X,'I',2X
     &     ,'BR(sum)',2X,'I',2X,'N',2X,'I'/1X,'I',50X,'I',11X,'I',5X,'I'
     &     /1X,70('=')/1X,'I',50X,'I',11X,'I',5X,'I')
 8100 FORMAT(1X,'I',50X,'I',11X,'I',5X,'I'/1X,70('=')/1X,'I',1X
     &     ,'Total number of R-Violating modes :',3X,I5,24X,'I'/
     &     1X,'I',1X,'Total number with non-vanishing BR :',2X,I5,24X
     &     ,'I'/1X,'I',1X,'Total number with BR > 0.001 :',8X,I5,24X,'I'
     &     /1X,70('='))
 8200 FORMAT(1X,'I',1X,A9,1X,'-->',1X,A24,11X,
     &     'I',2X,1P,D8.2,0P,1X,'I',2X,I2,1X,'I')
 8300 FORMAT(1X,'I',50X,'I',11X,'I',5X,'I')
 8500 FORMAT(1X/ 1X/
     &     1X,'R-Violating couplings',1X/ 1X /
     &     1X,55('=')/
     &     1X,'I',1X,'IJK',1X,'I',2X,'LAMBDA(IJK)',2X,'I',2X
     &     ,'LAMBDA''(IJK)',1X,'I',1X,"LAMBDA''(IJK)",1X,'I'/1X,'I',5X
     &     ,'I',15X,'I',15X,'I',15X,'I')
 8600 FORMAT(1X,55('='))
 8700 FORMAT(1X,'I',1X,I1,I1,I1,1X,'I',1X,1P,D13.3,0P,1X,'I',1X,1P
     &     ,D13.3,0P,1X,'I',1X,1P,D13.3,0P,1X,'I')
 
      RETURN
      END
