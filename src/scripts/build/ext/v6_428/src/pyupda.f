 
C*********************************************************************
 
C...PYUPDA
C...Facilitates the updating of particle and decay data
C...by allowing it to be done in an external file.
 
      SUBROUTINE PYUPDA(MUPDA,LFN)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYDAT4/,/PYINT4/
C...Local arrays, character variables and data.
      CHARACTER CHINL*120,CHKF*9,CHVAR(22)*9,CHLIN*72,
     &CHBLK(20)*72,CHOLD*16,CHTMP*16,CHNEW*16,CHCOM*24
      DATA CHVAR/ 'KCHG(I,1)','KCHG(I,2)','KCHG(I,3)','KCHG(I,4)',
     &'PMAS(I,1)','PMAS(I,2)','PMAS(I,3)','PMAS(I,4)','MDCY(I,1)',
     &'MDCY(I,2)','MDCY(I,3)','MDME(I,1)','MDME(I,2)','BRAT(I)  ',
     &'KFDP(I,1)','KFDP(I,2)','KFDP(I,3)','KFDP(I,4)','KFDP(I,5)',
     &'CHAF(I,1)','CHAF(I,2)','MWID(I)  '/
 
C...Write header if not yet done.
      IF(MSTU(12).NE.12345) CALL PYLIST(0)
 
C...Write information on file for editing.
      IF(MUPDA.EQ.1) THEN
        DO 110 KC=1,500
          WRITE(LFN,5000) KCHG(KC,4),(CHAF(KC,J1),J1=1,2),
     &    (KCHG(KC,J2),J2=1,3),(PMAS(KC,J3),J3=1,4),
     &    MWID(KC),MDCY(KC,1)
          DO 100 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1
            WRITE(LFN,5100) MDME(IDC,1),MDME(IDC,2),BRAT(IDC),
     &      (KFDP(IDC,J),J=1,5)
  100     CONTINUE
  110   CONTINUE
 
C...Read complete set of information from edited file or
C...read partial set of new or updated information from edited file.
      ELSEIF(MUPDA.EQ.2.OR.MUPDA.EQ.3) THEN
 
C...Reset counters.
        KCC=100
        NDC=0
        CHKF='         '
        IF(MUPDA.EQ.2) THEN
          DO 120 I=1,MSTU(6)
            KCHG(I,4)=0
  120     CONTINUE
        ELSE
          DO 130 KC=1,MSTU(6)
            IF(KC.GT.100.AND.KCHG(KC,4).GT.100) KCC=KC
            NDC=MAX(NDC,MDCY(KC,2)+MDCY(KC,3)-1)
  130     CONTINUE
        ENDIF
 
C...Begin of loop: read new line; unknown whether particle or
C...decay data.
  140   READ(LFN,5200,END=190) CHINL
 
C...Identify particle code and whether already defined  (for MUPDA=3).
        IF(CHINL(2:10).NE.'         ') THEN
          CHKF=CHINL(2:10)
          READ(CHKF,5300) KF
          IF(MUPDA.EQ.2) THEN
            IF(KF.LE.100) THEN
              KC=KF
            ELSE
              KCC=KCC+1
              KC=KCC
            ENDIF
          ELSE
            KCREP=0
            IF(KF.LE.100) THEN
              KCREP=KF
            ELSE
              DO 150 KCR=101,KCC
                IF(KCHG(KCR,4).EQ.KF) KCREP=KCR
  150         CONTINUE
            ENDIF
C...Remove duplicate old decay data.
            IF(KCREP.NE.0.AND.MDCY(KCREP,3).GT.0) THEN
              IDCREP=MDCY(KCREP,2)
              NDCREP=MDCY(KCREP,3)
              DO 160 I=1,KCC
                IF(MDCY(I,2).GT.IDCREP) MDCY(I,2)=MDCY(I,2)-NDCREP
  160         CONTINUE
              DO 180 I=IDCREP,NDC-NDCREP
                MDME(I,1)=MDME(I+NDCREP,1)
                MDME(I,2)=MDME(I+NDCREP,2)
                BRAT(I)=BRAT(I+NDCREP)
                DO 170 J=1,5
                  KFDP(I,J)=KFDP(I+NDCREP,J)
  170           CONTINUE
  180         CONTINUE
              NDC=NDC-NDCREP
              KC=KCREP
            ELSEIF(KCREP.NE.0) THEN
              KC=KCREP
            ELSE
              KCC=KCC+1
              KC=KCC
            ENDIF
          ENDIF
 
C...Study line with particle data.
          IF(KC.GT.MSTU(6)) CALL PYERRM(27,
     &    '(PYUPDA:) Particle arrays full by KF ='//CHKF)
          READ(CHINL,5000) KCHG(KC,4),(CHAF(KC,J1),J1=1,2),
     &    (KCHG(KC,J2),J2=1,3),(PMAS(KC,J3),J3=1,4),
     &    MWID(KC),MDCY(KC,1)
          MDCY(KC,2)=0
          MDCY(KC,3)=0
 
C...Study line with decay data.
        ELSE
          NDC=NDC+1
          IF(NDC.GT.MSTU(7)) CALL PYERRM(27,
     &    '(PYUPDA:) Decay data arrays full by KF ='//CHKF)
          IF(MDCY(KC,2).EQ.0) MDCY(KC,2)=NDC
          MDCY(KC,3)=MDCY(KC,3)+1
          READ(CHINL,5100) MDME(NDC,1),MDME(NDC,2),BRAT(NDC),
     &    (KFDP(NDC,J),J=1,5)
        ENDIF
 
C...End of loop; ensure that PYCOMP tables are updated.
        GOTO 140
  190   CONTINUE
        MSTU(20)=0
 
C...Perform possible tests that new information is consistent.
        DO 220 KC=1,MSTU(6)
          KF=KCHG(KC,4)
          IF(KF.EQ.0) GOTO 220
          WRITE(CHKF,5300) KF
          IF(MIN(PMAS(KC,1),PMAS(KC,2),PMAS(KC,3),PMAS(KC,1)-PMAS(KC,3),
     &    PMAS(KC,4)).LT.0D0.OR.MDCY(KC,3).LT.0) CALL PYERRM(17,
     &    '(PYUPDA:) Mass/width/life/(# channels) wrong for KF ='//CHKF)
          BRSUM=0D0
          DO 210 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1
            IF(MDME(IDC,2).GT.80) GOTO 210
            KQ=KCHG(KC,1)
            PMS=PMAS(KC,1)-PMAS(KC,3)-PARJ(64)
            MERR=0
            DO 200 J=1,5
              KP=KFDP(IDC,J)
              IF(KP.EQ.0.OR.KP.EQ.81.OR.IABS(KP).EQ.82) THEN
                IF(KP.EQ.81) KQ=0
              ELSEIF(PYCOMP(KP).EQ.0) THEN
                MERR=3
              ELSE
                KQ=KQ-PYCHGE(KP)
                KPC=PYCOMP(KP)
                PMS=PMS-PMAS(KPC,1)
                IF(MSTJ(24).GT.0) PMS=PMS+0.5D0*MIN(PMAS(KPC,2),
     &          PMAS(KPC,3))
              ENDIF
  200       CONTINUE
            IF(KQ.NE.0) MERR=MAX(2,MERR)
            IF(MWID(KC).EQ.0.AND.KF.NE.311.AND.PMS.LT.0D0)
     &      MERR=MAX(1,MERR)
            IF(MERR.EQ.3) CALL PYERRM(17,
     &      '(PYUPDA:) Unknown particle code in decay of KF ='//CHKF)
            IF(MERR.EQ.2) CALL PYERRM(17,
     &      '(PYUPDA:) Charge not conserved in decay of KF ='//CHKF)
            IF(MERR.EQ.1) CALL PYERRM(7,
     &      '(PYUPDA:) Kinematically unallowed decay of KF ='//CHKF)
            BRSUM=BRSUM+BRAT(IDC)
  210     CONTINUE
          WRITE(CHTMP,5500) BRSUM
          IF(ABS(BRSUM).GT.0.0005D0.AND.ABS(BRSUM-1D0).GT.0.0005D0)
     &    CALL PYERRM(7,'(PYUPDA:) Sum of branching ratios is '//
     &    CHTMP(9:16)//' for KF ='//CHKF)
  220   CONTINUE
 
C...Write DATA statements for inclusion in program.
      ELSEIF(MUPDA.EQ.4) THEN
 
C...Find out how many codes and decay channels are actually used.
        KCC=0
        NDC=0
        DO 230 I=1,MSTU(6)
          IF(KCHG(I,4).NE.0) THEN
            KCC=I
            NDC=MAX(NDC,MDCY(I,2)+MDCY(I,3)-1)
          ENDIF
  230   CONTINUE
 
C...Initialize writing of DATA statements for inclusion in program.
        DO 300 IVAR=1,22
          NDIM=MSTU(6)
          IF(IVAR.GE.12.AND.IVAR.LE.19) NDIM=MSTU(7)
          NLIN=1
          CHLIN=' '
          CHLIN(7:35)='DATA ('//CHVAR(IVAR)//',I=   1,    )/'
          LLIN=35
          CHOLD='START'
 
C...Loop through variables for conversion to characters.
          DO 280 IDIM=1,NDIM
            IF(IVAR.EQ.1) WRITE(CHTMP,5400) KCHG(IDIM,1)
            IF(IVAR.EQ.2) WRITE(CHTMP,5400) KCHG(IDIM,2)
            IF(IVAR.EQ.3) WRITE(CHTMP,5400) KCHG(IDIM,3)
            IF(IVAR.EQ.4) WRITE(CHTMP,5400) KCHG(IDIM,4)
            IF(IVAR.EQ.5) WRITE(CHTMP,5500) PMAS(IDIM,1)
            IF(IVAR.EQ.6) WRITE(CHTMP,5500) PMAS(IDIM,2)
            IF(IVAR.EQ.7) WRITE(CHTMP,5500) PMAS(IDIM,3)
            IF(IVAR.EQ.8) WRITE(CHTMP,5500) PMAS(IDIM,4)
            IF(IVAR.EQ.9) WRITE(CHTMP,5400) MDCY(IDIM,1)
            IF(IVAR.EQ.10) WRITE(CHTMP,5400) MDCY(IDIM,2)
            IF(IVAR.EQ.11) WRITE(CHTMP,5400) MDCY(IDIM,3)
            IF(IVAR.EQ.12) WRITE(CHTMP,5400) MDME(IDIM,1)
            IF(IVAR.EQ.13) WRITE(CHTMP,5400) MDME(IDIM,2)
            IF(IVAR.EQ.14) WRITE(CHTMP,5600) BRAT(IDIM)
            IF(IVAR.EQ.15) WRITE(CHTMP,5400) KFDP(IDIM,1)
            IF(IVAR.EQ.16) WRITE(CHTMP,5400) KFDP(IDIM,2)
            IF(IVAR.EQ.17) WRITE(CHTMP,5400) KFDP(IDIM,3)
            IF(IVAR.EQ.18) WRITE(CHTMP,5400) KFDP(IDIM,4)
            IF(IVAR.EQ.19) WRITE(CHTMP,5400) KFDP(IDIM,5)
            IF(IVAR.EQ.20) CHTMP=CHAF(IDIM,1)
            IF(IVAR.EQ.21) CHTMP=CHAF(IDIM,2)
            IF(IVAR.EQ.22) WRITE(CHTMP,5400) MWID(IDIM)
 
C...Replace variables beyond what is properly defined.
            IF(IVAR.LE.4) THEN
              IF(IDIM.GT.KCC) CHTMP='               0'
            ELSEIF(IVAR.LE.8) THEN
              IF(IDIM.GT.KCC) CHTMP='             0.0'
            ELSEIF(IVAR.LE.11) THEN
              IF(IDIM.GT.KCC) CHTMP='               0'
            ELSEIF(IVAR.LE.13) THEN
              IF(IDIM.GT.NDC) CHTMP='               0'
            ELSEIF(IVAR.LE.14) THEN
              IF(IDIM.GT.NDC) CHTMP='             0.0'
            ELSEIF(IVAR.LE.19) THEN
              IF(IDIM.GT.NDC) CHTMP='               0'
            ELSEIF(IVAR.LE.21) THEN
              IF(IDIM.GT.KCC) CHTMP='                '
            ELSE
              IF(IDIM.GT.KCC) CHTMP='               0'
            ENDIF
 
C...Length of variable, trailing decimal zeros, quotation marks.
            LLOW=1
            LHIG=1
            DO 240 LL=1,16
              IF(CHTMP(17-LL:17-LL).NE.' ') LLOW=17-LL
              IF(CHTMP(LL:LL).NE.' ') LHIG=LL
  240       CONTINUE
            CHNEW=CHTMP(LLOW:LHIG)//' '
            LNEW=1+LHIG-LLOW
            IF((IVAR.GE.5.AND.IVAR.LE.8).OR.IVAR.EQ.14) THEN
              LNEW=LNEW+1
  250         LNEW=LNEW-1
              IF(LNEW.GE.2.AND.CHNEW(LNEW:LNEW).EQ.'0') GOTO 250
              IF(CHNEW(LNEW:LNEW).EQ.'.') LNEW=LNEW-1
              IF(LNEW.EQ.0) THEN
                CHNEW(1:3)='0D0'
                LNEW=3
              ELSE
                CHNEW(LNEW+1:LNEW+2)='D0'
                LNEW=LNEW+2
              ENDIF
            ELSEIF(IVAR.EQ.20.OR.IVAR.EQ.21) THEN
              DO 260 LL=LNEW,1,-1
                IF(CHNEW(LL:LL).EQ.'''') THEN
                  CHTMP=CHNEW
                  CHNEW=CHTMP(1:LL)//''''//CHTMP(LL+1:11)
                  LNEW=LNEW+1
                ENDIF
  260         CONTINUE
              LNEW=MIN(14,LNEW)
              CHTMP=CHNEW
              CHNEW(1:LNEW+2)=''''//CHTMP(1:LNEW)//''''
              LNEW=LNEW+2
            ENDIF
 
C...Form composite character string, often including repetition counter.
            IF(CHNEW.NE.CHOLD) THEN
              NRPT=1
              CHOLD=CHNEW
              CHCOM=CHNEW
              LCOM=LNEW
            ELSE
              LRPT=LNEW+1
              IF(NRPT.GE.2) LRPT=LNEW+3
              IF(NRPT.GE.10) LRPT=LNEW+4
              IF(NRPT.GE.100) LRPT=LNEW+5
              IF(NRPT.GE.1000) LRPT=LNEW+6
              LLIN=LLIN-LRPT
              NRPT=NRPT+1
              WRITE(CHTMP,5400) NRPT
              LRPT=1
              IF(NRPT.GE.10) LRPT=2
              IF(NRPT.GE.100) LRPT=3
              IF(NRPT.GE.1000) LRPT=4
              CHCOM(1:LRPT+1+LNEW)=CHTMP(17-LRPT:16)//'*'//CHNEW(1:LNEW)
              LCOM=LRPT+1+LNEW
            ENDIF
 
C...Add characters to end of line, to new line (after storing old line),
C...or to new block of lines (after writing old block).
            IF(LLIN+LCOM.LE.70) THEN
              CHLIN(LLIN+1:LLIN+LCOM+1)=CHCOM(1:LCOM)//','
              LLIN=LLIN+LCOM+1
            ELSEIF(NLIN.LE.19) THEN
              CHLIN(LLIN+1:72)=' '
              CHBLK(NLIN)=CHLIN
              NLIN=NLIN+1
              CHLIN(6:6+LCOM+1)='&'//CHCOM(1:LCOM)//','
              LLIN=6+LCOM+1
            ELSE
              CHLIN(LLIN:72)='/'//' '
              CHBLK(NLIN)=CHLIN
              WRITE(CHTMP,5400) IDIM-NRPT
              CHBLK(1)(30:33)=CHTMP(13:16)
              DO 270 ILIN=1,NLIN
                WRITE(LFN,5700) CHBLK(ILIN)
  270         CONTINUE
              NLIN=1
              CHLIN=' '
              CHLIN(7:35+LCOM+1)='DATA ('//CHVAR(IVAR)//
     &        ',I=    ,    )/'//CHCOM(1:LCOM)//','
              WRITE(CHTMP,5400) IDIM-NRPT+1
              CHLIN(25:28)=CHTMP(13:16)
              LLIN=35+LCOM+1
            ENDIF
  280     CONTINUE
 
C...Write final block of lines.
          CHLIN(LLIN:72)='/'//' '
          CHBLK(NLIN)=CHLIN
          WRITE(CHTMP,5400) NDIM
          CHBLK(1)(30:33)=CHTMP(13:16)
          DO 290 ILIN=1,NLIN
            WRITE(LFN,5700) CHBLK(ILIN)
  290     CONTINUE
  300   CONTINUE
      ENDIF
 
C...Formats for reading and writing particle data.
 5000 FORMAT(1X,I9,2X,A16,2X,A16,3I3,3F12.5,1P,E13.5,2I3)
 5100 FORMAT(10X,2I5,F12.6,5I10)
 5200 FORMAT(A120)
 5300 FORMAT(I9)
 5400 FORMAT(I16)
 5500 FORMAT(F16.5)
 5600 FORMAT(F16.6)
 5700 FORMAT(A72)
 
      RETURN
      END
