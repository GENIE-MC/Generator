 
C*********************************************************************
 
C...PYONOF
C...Switches on and off decay channel by search for match.
 
      SUBROUTINE PYONOF(CHIN)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      SAVE /PYDAT1/,/PYDAT3/
C...Local arrays and character variables.
      INTEGER KFCMP(10),KFTMP(10)
      CHARACTER CHIN*(*),CHTMP*104,CHFIX*104,CHMODE*10,CHCODE*8,
     &CHALP(2)*26
      DATA CHALP/'abcdefghijklmnopqrstuvwxyz',
     &'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

C...Determine length of character variable.
      CHTMP=CHIN//' '
      LBEG=0
  100 LBEG=LBEG+1
      IF(CHTMP(LBEG:LBEG).EQ.' ') GOTO 100
      LEND=LBEG-1
  105 LEND=LEND+1
      IF(LEND.LE.100.AND.CHTMP(LEND:LEND).NE.'!') GOTO 105
  110 LEND=LEND-1
      IF(CHTMP(LEND:LEND).EQ.' ') GOTO 110
      LEN=1+LEND-LBEG
      CHFIX(1:LEN)=CHTMP(LBEG:LEND)

C...Find colon separator and particle code.
      LCOLON=0
  120 LCOLON=LCOLON+1
      IF(CHFIX(LCOLON:LCOLON).NE.':') GOTO 120
      CHCODE=' '
      CHCODE(10-LCOLON:8)=CHFIX(1:LCOLON-1)
      READ(CHCODE,'(I8)',ERR=300) KF
      KC=PYCOMP(KF)

C...Done if unknown code or no decay channels.
      IF(KC.EQ.0) THEN
        CALL PYERRM(18,'(PYONOF:) unrecognized particle '//CHCODE)
        RETURN
      ENDIF
      IDCBEG=MDCY(KC,2)
      IDCLEN=MDCY(KC,3)
      IF(IDCBEG.EQ.0.OR.IDCLEN.EQ.0) THEN
        CALL PYERRM(18,'(PYONOF:) no decay channels for '//CHCODE)
        RETURN
      ENDIF

C...Find command name up to blank or equal sign.
      LSEP=LCOLON
  130 LSEP=LSEP+1
      IF(LSEP.LE.LEN.AND.CHFIX(LSEP:LSEP).NE.' '.AND.
     &CHFIX(LSEP:LSEP).NE.'=') GOTO 130
      CHMODE=' '
      LMODE=LSEP-LCOLON-1
      CHMODE(1:LMODE)=CHFIX(LCOLON+1:LSEP-1)

C...Convert to uppercase.
      DO 150 LCOM=1,LMODE
        DO 140 LALP=1,26
          IF(CHMODE(LCOM:LCOM).EQ.CHALP(1)(LALP:LALP)) 
     &    CHMODE(LCOM:LCOM)=CHALP(2)(LALP:LALP)
  140   CONTINUE
  150 CONTINUE

C...Identify command. Failed if not identified.
      MODE=0
      IF(CHMODE.EQ.'ALLOFF') MODE=1
      IF(CHMODE.EQ.'ALLON') MODE=2
      IF(CHMODE.EQ.'OFFIFANY') MODE=3
      IF(CHMODE.EQ.'ONIFANY') MODE=4
      IF(CHMODE.EQ.'OFFIFALL') MODE=5
      IF(CHMODE.EQ.'ONIFALL') MODE=6
      IF(CHMODE.EQ.'OFFIFMATCH') MODE=7
      IF(CHMODE.EQ.'ONIFMATCH') MODE=8
      IF(MODE.EQ.0) THEN
        CALL PYERRM(18,'(PYONOF:) unknown command '//CHMODE)
        RETURN
      ENDIF

C...Simple cases when all on or all off.
      IF(MODE.EQ.1.OR.MODE.EQ.2) THEN
        WRITE(MSTU(11),1000) KF,CHMODE
        DO 160 IDC=IDCBEG,IDCBEG+IDCLEN-1
          IF(MDME(IDC,1).LT.0) GOTO 160
          MDME(IDC,1)=MODE-1
  160   CONTINUE
        RETURN
      ENDIF

C...Identify matching list.
      NCMP=0
      LBEG=LSEP
  170 LBEG=LBEG+1
      IF(LBEG.GT.LEN) GOTO 190
      IF(LBEG.LT.LEN.AND.(CHFIX(LBEG:LBEG).EQ.' '.OR.
     &CHFIX(LBEG:LBEG).EQ.'='.OR.CHFIX(LBEG:LBEG).EQ.',')) GOTO 170
      LEND=LBEG-1
  180 LEND=LEND+1
      IF(LEND.LT.LEN.AND.CHFIX(LEND:LEND).NE.' '.AND.
     &CHFIX(LEND:LEND).NE.'='.AND.CHFIX(LEND:LEND).NE.',') GOTO 180
      IF(LEND.LT.LEN) LEND=LEND-1
      CHCODE=' '
      CHCODE(8-LEND+LBEG:8)=CHFIX(LBEG:LEND)
      READ(CHCODE,'(I8)',ERR=300) KFREAD
      NCMP=NCMP+1
      KFCMP(NCMP)=IABS(KFREAD)
      LBEG=LEND
      IF(NCMP.LT.10) GOTO 170
  190 CONTINUE
      WRITE(MSTU(11),1100) KF,CHMODE,(KFCMP(ICMP),ICMP=1,NCMP)

C...Only one matching required.
      IF(MODE.EQ.3.OR.MODE.EQ.4) THEN
        DO 220 IDC=IDCBEG,IDCBEG+IDCLEN-1
          IF(MDME(IDC,1).LT.0) GOTO 220
          DO 210 IKF=1,5
            KFNOW=IABS(KFDP(IDC,IKF))
            IF(KFNOW.EQ.0) GOTO 210
            DO 200 ICMP=1,NCMP
              IF(KFCMP(ICMP).EQ.KFNOW) THEN
                MDME(IDC,1)=MODE-3
                GOTO 220
              ENDIF
  200      CONTINUE
  210     CONTINUE
  220   CONTINUE
        RETURN
      ENDIF

C...Multiple matchings required.
      DO 260 IDC=IDCBEG,IDCBEG+IDCLEN-1
        IF(MDME(IDC,1).LT.0) GOTO 260
        NTMP=NCMP
        DO 230 ITMP=1,NTMP
          KFTMP(ITMP)=KFCMP(ITMP)
  230   CONTINUE  
        NFIN=0 
        DO 250 IKF=1,5
          KFNOW=IABS(KFDP(IDC,IKF))
          IF(KFNOW.EQ.0) GOTO 250
          NFIN=NFIN+1
          DO 240 ITMP=1,NTMP
            IF(KFTMP(ITMP).EQ.KFNOW) THEN
              KFTMP(ITMP)=KFTMP(NTMP) 
              NTMP=NTMP-1
              GOTO 250
            ENDIF
  240     CONTINUE
  250   CONTINUE
        IF(NTMP.EQ.0.AND.MODE.LE.6) MDME(IDC,1)=MODE-5
        IF(NTMP.EQ.0.AND.NFIN.EQ.NCMP.AND.MODE.GE.7) 
     &  MDME(IDC,1)=MODE-7
  260 CONTINUE
      RETURN

C...Error exit for impossible read of particle code.
  300 CALL PYERRM(18,'(PYONOF:) could not interpret particle code '
     &//CHCODE)

C...Formats for output.
 1000 FORMAT(' Decays for',I8,' set ',A10)
 1100 FORMAT(' Decays for',I8,' set ',A10,' if match',10I8)

      RETURN
      END
