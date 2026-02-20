 
C*********************************************************************
 
C...PYSGHF
C...Subprocess cross sections for heavy flavour production,
C...open and closed.
C...Auxiliary to PYSIGH.
 
      SUBROUTINE PYSGHF(NCHN,SIGS)
 
C...Double precision and integer declarations
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYSGCM/ISUB,ISUBSV,MMIN1,MMAX1,MMIN2,MMAX2,MMINA,MMAXA,
     &KFAC(2,-40:40),COMFAC,FACK,FACA,SH,TH,UH,SH2,TH2,UH2,SQM3,SQM4,
     &SHR,SQPTH,TAUP,BE34,CTH,X(2),SQMZ,SQMW,GMMZ,GMMW,
     &AEM,AS,XW,XW1,XWC,XWV,POLL,POLR,POLLL,POLRR
      SAVE /PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT3/,
     &/PYINT4/,/PYSGCM/
C...Local arrays
      DIMENSION WDTP(0:400),WDTE(0:400,0:5)
 
C...Determine where are charmonium/bottomonium wave function parameters.
      IONIUM=140
      IF(ISUB.GE.461.AND.ISUB.LE.479) IONIUM=145
 
C...Convert bottomonium process into equivalent charmonium ones.
      IF(ISUB.GE.461.AND.ISUB.LE.479) ISUB=ISUB-40
 
C...Differential cross section expressions.
 
      IF(ISUB.LE.100) THEN
        IF(ISUB.EQ.81) THEN
C...q + qbar -> Q + Qbar
          SQMAVG=0.5D0*(SQM3+SQM4)-0.25D0*(SQM3-SQM4)**2/SH
          THQ=-0.5D0*SH*(1D0-BE34*CTH)
          UHQ=-0.5D0*SH*(1D0+BE34*CTH)
          FACQQB=COMFAC*AS**2*4D0/9D0*((THQ**2+UHQ**2)/SH2+
     &    2D0*SQMAVG/SH)
          IF(MSTP(35).GE.1) FACQQB=FACQQB*PYHFTH(SH,SQMAVG,0D0)
          WID2=1D0
          IF(MINT(55).EQ.6) WID2=WIDS(6,1)
          IF(MINT(55).EQ.7.OR.MINT(55).EQ.8) WID2=WIDS(MINT(55),1)
          FACQQB=FACQQB*WID2
          DO 100 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 100
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQB
  100     CONTINUE
 
        ELSEIF(ISUB.EQ.82) THEN
C...g + g -> Q + Qbar
          SQMAVG=0.5D0*(SQM3+SQM4)-0.25D0*(SQM3-SQM4)**2/SH
          THQ=-0.5D0*SH*(1D0-BE34*CTH)
          UHQ=-0.5D0*SH*(1D0+BE34*CTH)
          THUHQ=THQ*UHQ-SQMAVG*SH
          IF(MSTP(34).EQ.0) THEN
            FACQQ1=UHQ/THQ-2D0*UHQ**2/SH2+4D0*(SQMAVG/SH)*THUHQ/THQ**2
            FACQQ2=THQ/UHQ-2D0*THQ**2/SH2+4D0*(SQMAVG/SH)*THUHQ/UHQ**2
          ELSE
            FACQQ1=UHQ/THQ-2.25D0*UHQ**2/SH2+4.5D0*(SQMAVG/SH)*THUHQ/
     &      THQ**2+0.5D0*SQMAVG*(THQ+SQMAVG)/THQ**2-SQMAVG**2/(SH*THQ)
            FACQQ2=THQ/UHQ-2.25D0*THQ**2/SH2+4.5D0*(SQMAVG/SH)*THUHQ/
     &      UHQ**2+0.5D0*SQMAVG*(UHQ+SQMAVG)/UHQ**2-SQMAVG**2/(SH*UHQ)
          ENDIF
          FACQQ1=COMFAC*FACA*AS**2*(1D0/6D0)*FACQQ1
          FACQQ2=COMFAC*FACA*AS**2*(1D0/6D0)*FACQQ2
          IF(MSTP(35).GE.1) THEN
            FATRE=PYHFTH(SH,SQMAVG,2D0/7D0)
            FACQQ1=FACQQ1*FATRE
            FACQQ2=FACQQ2*FATRE
          ENDIF
          WID2=1D0
          IF(MINT(55).EQ.6) WID2=WIDS(6,1)
          IF(MINT(55).EQ.7.OR.MINT(55).EQ.8) WID2=WIDS(MINT(55),1)
          FACQQ1=FACQQ1*WID2
          FACQQ2=FACQQ2*WID2
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 110
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ1
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=2
          SIGH(NCHN)=FACQQ2
  110     CONTINUE
 
        ELSEIF(ISUB.EQ.83) THEN
C...f + q -> f' + Q
          FACQQS=COMFAC*(0.5D0*AEM/XW)**2*SH*(SH-SQM3)/(SQMW-TH)**2
          FACQQU=COMFAC*(0.5D0*AEM/XW)**2*UH*(UH-SQM3)/(SQMW-TH)**2
          DO 130 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 130
            DO 120 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 120
              IF(I*J.GT.0.AND.MOD(IABS(I+J),2).EQ.0) GOTO 120
              IF(I*J.LT.0.AND.MOD(IABS(I+J),2).EQ.1) GOTO 120
              IF(IABS(I).LT.MINT(55).AND.MOD(IABS(I+MINT(55)),2).EQ.1)
     &        THEN
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=1
                IF(MOD(MINT(55),2).EQ.0) FACCKM=VCKM(MINT(55)/2,
     &          (IABS(I)+1)/2)*VINT(180+J)
                IF(MOD(MINT(55),2).EQ.1) FACCKM=VCKM(IABS(I)/2,
     &          (MINT(55)+1)/2)*VINT(180+J)
                WID2=1D0
                IF(I.GT.0) THEN
                  IF(MINT(55).EQ.6) WID2=WIDS(6,2)
                  IF(MINT(55).EQ.7.OR.MINT(55).EQ.8) WID2=
     &            WIDS(MINT(55),2)
                ELSE
                  IF(MINT(55).EQ.6) WID2=WIDS(6,3)
                  IF(MINT(55).EQ.7.OR.MINT(55).EQ.8) WID2=
     &            WIDS(MINT(55),3)
                ENDIF
                IF(I*J.GT.0) SIGH(NCHN)=FACQQS*FACCKM*WID2
                IF(I*J.LT.0) SIGH(NCHN)=FACQQU*FACCKM*WID2
              ENDIF
              IF(IABS(J).LT.MINT(55).AND.MOD(IABS(J+MINT(55)),2).EQ.1)
     &        THEN
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=2
                IF(MOD(MINT(55),2).EQ.0) FACCKM=VCKM(MINT(55)/2,
     &          (IABS(J)+1)/2)*VINT(180+I)
                IF(MOD(MINT(55),2).EQ.1) FACCKM=VCKM(IABS(J)/2,
     &          (MINT(55)+1)/2)*VINT(180+I)
                WID2=1D0
                IF(J.GT.0) THEN
                  IF(MINT(55).EQ.6) WID2=WIDS(6,2)
                  IF(MINT(55).EQ.7.OR.MINT(55).EQ.8) WID2=
     &            WIDS(MINT(55),2)
                ELSE
                  IF(MINT(55).EQ.6) WID2=WIDS(6,3)
                  IF(MINT(55).EQ.7.OR.MINT(55).EQ.8) WID2=
     &            WIDS(MINT(55),3)
                ENDIF
                IF(I*J.GT.0) SIGH(NCHN)=FACQQS*FACCKM*WID2
                IF(I*J.LT.0) SIGH(NCHN)=FACQQU*FACCKM*WID2
              ENDIF
  120       CONTINUE
  130     CONTINUE
 
        ELSEIF(ISUB.EQ.84) THEN
C...g + gamma -> Q + Qbar
          SQMAVG=0.5D0*(SQM3+SQM4)-0.25D0*(SQM3-SQM4)**2/SH
          THQ=-0.5D0*SH*(1D0-BE34*CTH)
          UHQ=-0.5D0*SH*(1D0+BE34*CTH)
          FACQQ=COMFAC*AS*AEM*(KCHG(IABS(MINT(55)),1)/3D0)**2*
     &    (THQ**2+UHQ**2+4D0*SQMAVG*SH*(1D0-SQMAVG*SH/(THQ*UHQ)))/
     &    (THQ*UHQ)
          IF(MSTP(35).GE.1) FACQQ=FACQQ*PYHFTH(SH,SQMAVG,0D0)
          WID2=1D0
          IF(MINT(55).EQ.6) WID2=WIDS(6,1)
          IF(MINT(55).EQ.7.OR.MINT(55).EQ.8) WID2=WIDS(MINT(55),1)
          FACQQ=FACQQ*WID2
          IF(KFAC(1,21)*KFAC(2,22).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=22
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQ
          ENDIF
          IF(KFAC(1,22)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=22
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQ
          ENDIF
 
        ELSEIF(ISUB.EQ.85) THEN
C...gamma + gamma -> F + Fbar (heavy fermion, quark or lepton)
          SQMAVG=0.5D0*(SQM3+SQM4)-0.25D0*(SQM3-SQM4)**2/SH
          THQ=-0.5D0*SH*(1D0-BE34*CTH)
          UHQ=-0.5D0*SH*(1D0+BE34*CTH)
          FACFF=COMFAC*AEM**2*(KCHG(IABS(MINT(56)),1)/3D0)**4*2D0*
     &    ((1D0-PARJ(131)*PARJ(132))*(THQ*UHQ-SQMAVG*SH)*
     &    (UHQ**2+THQ**2+2D0*SQMAVG*SH)+(1D0+PARJ(131)*PARJ(132))*
     &    SQMAVG*SH**2*(SH-2D0*SQMAVG))/(THQ*UHQ)**2
          IF(IABS(MINT(56)).LT.10) FACFF=3D0*FACFF
          IF(IABS(MINT(56)).LT.10.AND.MSTP(35).GE.1)
     &    FACFF=FACFF*PYHFTH(SH,SQMAVG,1D0)
          WID2=1D0
          IF(MINT(56).EQ.6) WID2=WIDS(6,1)
          IF(MINT(56).EQ.7.OR.MINT(56).EQ.8) WID2=WIDS(MINT(56),1)
          IF(MINT(56).EQ.17) WID2=WIDS(17,1)
          FACFF=FACFF*WID2
          IF(KFAC(1,22)*KFAC(2,22).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=22
            ISIG(NCHN,2)=22
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACFF
          ENDIF
 
        ELSEIF(ISUB.EQ.86) THEN
C...g + g -> J/Psi + g
          FACQQG=COMFAC*AS**3*(5D0/9D0)*PARP(38)*SQRT(SQM3)*
     &    (((SH*(SH-SQM3))**2+(TH*(TH-SQM3))**2+(UH*(UH-SQM3))**2)/
     &    ((TH-SQM3)*(UH-SQM3))**2)/(SH-SQM3)**2
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG
          ENDIF
 
        ELSEIF(ISUB.EQ.87) THEN
C...g + g -> chi_0c + g
          PGTW=(SH*TH+TH*UH+UH*SH)/SH2
          QGTW=(SH*TH*UH)/SH**3
          RGTW=SQM3/SH
          FACQQG=COMFAC*AS**3*4D0*(PARP(39)/SQRT(SQM3))*(1D0/SH)*
     &    (9D0*RGTW**2*PGTW**4*(RGTW**4-2D0*RGTW**2*PGTW+PGTW**2)-
     &    6D0*RGTW*PGTW**3*QGTW*(2D0*RGTW**4-5D0*RGTW**2*PGTW+PGTW**2)-
     &    PGTW**2*QGTW**2*(RGTW**4+2D0*RGTW**2*PGTW-PGTW**2)+
     &    2D0*RGTW*PGTW*QGTW**3*(RGTW**2-PGTW)+6D0*RGTW**2*QGTW**4)/
     &    (QGTW*(QGTW-RGTW*PGTW)**4)
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG
          ENDIF
 
        ELSEIF(ISUB.EQ.88) THEN
C...g + g -> chi_1c + g
          PGTW=(SH*TH+TH*UH+UH*SH)/SH2
          QGTW=(SH*TH*UH)/SH**3
          RGTW=SQM3/SH
          FACQQG=COMFAC*AS**3*12D0*(PARP(39)/SQRT(SQM3))*(1D0/SH)*
     &    PGTW**2*(RGTW*PGTW**2*(RGTW**2-4D0*PGTW)+2D0*QGTW*(-RGTW**4+
     &    5D0*RGTW**2*PGTW+PGTW**2)-15D0*RGTW*QGTW**2)/
     &    (QGTW-RGTW*PGTW)**4
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG
          ENDIF
 
        ELSEIF(ISUB.EQ.89) THEN
C...g + g -> chi_2c + g
          PGTW=(SH*TH+TH*UH+UH*SH)/SH2
          QGTW=(SH*TH*UH)/SH**3
          RGTW=SQM3/SH
          FACQQG=COMFAC*AS**3*4D0*(PARP(39)/SQRT(SQM3))*(1D0/SH)*
     &    (12D0*RGTW**2*PGTW**4*(RGTW**4-2D0*RGTW**2*PGTW+PGTW**2)-
     &    3D0*RGTW*PGTW**3*QGTW*(8D0*RGTW**4-RGTW**2*PGTW+4D0*PGTW**2)+
     &    2D0*PGTW**2*QGTW**2*(-7D0*RGTW**4+43D0*RGTW**2*PGTW+PGTW**2)+
     &    RGTW*PGTW*QGTW**3*(16D0*RGTW**2-61D0*PGTW)+12D0*RGTW**2*
     &    QGTW**4)/(QGTW*(QGTW-RGTW*PGTW)**4)
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG
          ENDIF
        ENDIF
 
      ELSEIF(ISUB.LE.200) THEN
        IF(ISUB.EQ.104) THEN
C...g + g -> chi_c0.
          KC=PYCOMP(10441)
          FACBW=COMFAC*12D0*AS**2*PARP(39)*PMAS(KC,2)/
     &    ((SH-PMAS(KC,1)**2)**2+(PMAS(KC,1)*PMAS(KC,2))**2)
          IF(ABS(SQRT(SH)-PMAS(KC,1)).GT.50D0*PMAS(KC,2)) FACBW=0D0
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACBW
          ENDIF
 
        ELSEIF(ISUB.EQ.105) THEN
C...g + g -> chi_c2.
          KC=PYCOMP(445)
          FACBW=COMFAC*16D0*AS**2*PARP(39)*PMAS(KC,2)/
     &    ((SH-PMAS(KC,1)**2)**2+(PMAS(KC,1)*PMAS(KC,2))**2)
          IF(ABS(SQRT(SH)-PMAS(KC,1)).GT.50D0*PMAS(KC,2)) FACBW=0D0
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACBW
          ENDIF
 
        ELSEIF(ISUB.EQ.106) THEN
C...g + g -> J/Psi + gamma.
          EQ=KCHG(MOD(KFPR(ISUB,1)/10,10),1)/3D0
          FACQQG=COMFAC*AEM*EQ**2*AS**2*(4D0/3D0)*PARP(38)*SQRT(SQM3)*
     &    (((SH*(SH-SQM3))**2+(TH*(TH-SQM3))**2+(UH*(UH-SQM3))**2)/
     &    ((TH-SQM3)*(UH-SQM3))**2)/(SH-SQM3)**2
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG
          ENDIF
 
        ELSEIF(ISUB.EQ.107) THEN
C...g + gamma -> J/Psi + g.
          EQ=KCHG(MOD(KFPR(ISUB,1)/10,10),1)/3D0
          FACQQG=COMFAC*AEM*EQ**2*AS**2*(32D0/3D0)*PARP(38)*SQRT(SQM3)*
     &    (((SH*(SH-SQM3))**2+(TH*(TH-SQM3))**2+(UH*(UH-SQM3))**2)/
     &    ((TH-SQM3)*(UH-SQM3))**2)/(SH-SQM3)**2
          IF(KFAC(1,21)*KFAC(2,22).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=22
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG
          ENDIF
          IF(KFAC(1,22)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=22
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG
          ENDIF
 
        ELSEIF(ISUB.EQ.108) THEN
C...gamma + gamma -> J/Psi + gamma.
          EQ=KCHG(MOD(KFPR(ISUB,1)/10,10),1)/3D0
          FACQQG=COMFAC*AEM**3*EQ**6*384D0*PARP(38)*SQRT(SQM3)*
     &    (((SH*(SH-SQM3))**2+(TH*(TH-SQM3))**2+(UH*(UH-SQM3))**2)/
     &    ((TH-SQM3)*(UH-SQM3))**2)/(SH-SQM3)**2
          IF(KFAC(1,22)*KFAC(2,22).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=22
            ISIG(NCHN,2)=22
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG
          ENDIF
        ENDIF
 
C...QUARKONIA+++
C...Additional code by Stefan Wolf
      ELSE
 
C...Common code for quarkonium production.
        SHTH=SH+TH
        THUH=TH+UH
        UHSH=UH+SH
        SHTH2=SHTH**2
        THUH2=THUH**2
        UHSH2=UHSH**2
        IF ( (ISUB.GE.421.AND.ISUB.LE.424).OR.
     &       (ISUB.GE.431.AND.ISUB.LE.433)) THEN
          SQMQQ=SQM3
        ELSEIF((ISUB.GE.425.AND.ISUB.LE.430).OR.
     &         (ISUB.GE.434.AND.ISUB.LE.439)) THEN
          SQMQQ=SQM4
        ENDIF
        SQMQQR=SQRT(SQMQQ)
        IF(MSTP(145).EQ.1) THEN
           IF ( (ISUB.GE.421.AND.ISUB.LE.427).OR.
     &          (ISUB.GE.431.AND.ISUB.LE.436)) THEN
              AQ=UHSH/(2D0*X(1)) + SHTH/(2D0*X(2))
              BQ=UHSH/(2D0*X(1)) - SHTH/(2D0*X(2))
              ATILK1=X(1)*VINT(2)/2D0-UHSH/(2D0*SQMQQ)*AQ
              ATILK2=X(2)*VINT(2)/2D0-SHTH/(2D0*SQMQQ)*AQ
              BTILK1=-X(1)*VINT(2)/2D0-UHSH/(2D0*SQMQQ)*BQ
              BTILK2=X(2)*VINT(2)/2D0-SHTH/(2D0*SQMQQ)*BQ
           ELSEIF( (ISUB.GE.428.AND.ISUB.LE.430).OR.
     &             ISUB.GE.437) THEN
              AQ=SHTH/(2D0*X(1)) + UHSH/(2D0*X(2))
              BQ=SHTH/(2D0*X(1)) - UHSH/(2D0*X(2))
              ATILK1=X(1)*VINT(2)/2D0-SHTH/(2D0*SQMQQ)*AQ
              ATILK2=X(2)*VINT(2)/2D0-UHSH/(2D0*SQMQQ)*AQ
              BTILK1=-X(1)*VINT(2)/2D0-SHTH/(2D0*SQMQQ)*BQ
              BTILK2=X(2)*VINT(2)/2D0-UHSH/(2D0*SQMQQ)*BQ
           ENDIF
           AQ2=AQ**2
           BQ2=BQ**2
           SMQQ2=SQMQQ*VINT(2)
C...Polarisation frames
           IF(MSTP(146).EQ.1) THEN
C...Recoil frame
              POLH1=SQRT(AQ2-SMQQ2)
              POLH2=SQRT(VINT(2)*(AQ2-BQ2-SMQQ2))
              AZ=-SQMQQR/POLH1
              BZ=0D0
              AX=AQ*BQ/(POLH1*POLH2)
              BX=-POLH1/POLH2
           ELSEIF(MSTP(146).EQ.2) THEN
C...Gottfried Jackson frame
              POLH1=AQ+BQ
              POLH2=POLH1*SQRT(VINT(2)*(AQ2-BQ2-SMQQ2))
              AZ=SQMQQR/POLH1
              BZ=AZ
              AX=-(BQ2+AQ*BQ+SMQQ2)/POLH2
              BX=(AQ2+AQ*BQ-SMQQ2)/POLH2
           ELSEIF(MSTP(146).EQ.3) THEN
C...Target frame
              POLH1=AQ-BQ
              POLH2=POLH1*SQRT(VINT(2)*(AQ2-BQ2-SMQQ2))
              AZ=-SQMQQR/POLH1
              BZ=-AZ
              AX=-(BQ2-AQ*BQ+SMQQ2)/POLH2
              BX=-(AQ2-AQ*BQ-SMQQ2)/POLH2
           ELSEIF(MSTP(146).EQ.4) THEN
C...Collins Soper frame
              POLH1=AQ2-BQ2
              POLH2=SQRT(VINT(2)*POLH1)
              AZ=-BQ/POLH2
              BZ=AQ/POLH2
              AX=-SQMQQR*AQ/SQRT(POLH1*(POLH1-SMQQ2))
              BX=SQMQQR*BQ/SQRT(POLH1*(POLH1-SMQQ2))
           ENDIF
C...Contract EL1(lam) EL2(lam') with K1 and K2 (initial parton momenta)
           EL1K10=AZ*ATILK1+BZ*BTILK1
           EL1K20=AZ*ATILK2+BZ*BTILK2
           EL2K10=EL1K10
           EL2K20=EL1K20
           EL1K11=1D0/SQRT(2D0)*(AX*ATILK1+BX*BTILK1)
           EL1K21=1D0/SQRT(2D0)*(AX*ATILK2+BX*BTILK2)
           EL2K11=EL1K11
           EL2K21=EL1K21
        ENDIF
 
        IF(ISUB.EQ.421) THEN
C...g + g -> QQ~[3S11] + g
          IF(MSTP(145).EQ.0) THEN
*            FACQQG=COMFAC*PARU(1)*AS**3*(10D0/81D0)*SQMQQR*
*     &            (SH2*THUH2+TH2*UHSH2+UH2*SHTH2)/(SHTH2*THUH2*UHSH2)
            FACQQG=COMFAC*PARU(1)*AS**3*(10D0/81D0)*SQMQQR*
     &            (SH2*THUH2+TH2*UHSH2+UH2*SHTH2)/SHTH2/THUH2/UHSH2
*            FACQQG=COMFAC*PARU(1)*AS**3*(10D0/81D0)*SQMQQR*
*     &           (SH2/(SHTH2*UHSH2)+TH2/(SHTH2*THUH2)+UH2/(THUH2*UHSH2))
          ELSE
            FF=-PARU(1)*AS**3*(10D0/81D0)*SQMQQR/THUH2/SHTH2/UHSH2
            AA=(SHTH2*UH2+UHSH2*TH2+THUH2*SH2)/2D0
            BB=2D0*(SH2+TH2)
            CC=2D0*(SH2+UH2)
            DD=2D0*SH2
            IF(MSTP(147).EQ.0) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=2D0*(-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11)))
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K10+CC*EL1K21*EL2K20
     &              +DD*(EL1K11*EL2K20+EL1K21*EL2K10))
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG*PARP(IONIUM+1)
          ENDIF
 
        ELSEIF(ISUB.EQ.422) THEN
C...g + g -> QQ~[3S18] + g
          IF(MSTP(145).EQ.0) THEN
            FACQQG=-COMFAC*PARU(1)*AS**3*(1D0/72D0)*
     &            (16D0*SQMQQ**2-27D0*(SHTH2+THUH2+UHSH2))/
     &            (SQMQQ*SQMQQR)*
     &            ((SH2*THUH2+TH2*UHSH2+UH2*SHTH2)/SHTH2/THUH2/UHSH2)
          ELSE
            FF=PARU(1)*AS**3*(16D0*SQMQQ**2-27D0*(SHTH2+THUH2+UHSH2))/
     &            (72D0*SQMQQ*SQMQQR*SHTH2*THUH2*UHSH2)
            AA=(SHTH2*UH2+UHSH2*TH2+THUH2*SH2)/2D0
            BB=2D0*(SH2+TH2)
            CC=2D0*(SH2+UH2)
            DD=2D0*SH2
            IF(MSTP(147).EQ.0) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=2D0*(-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11)))
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K10+CC*EL1K21*EL2K20
     &              +DD*(EL1K11*EL2K20+EL1K21*EL2K10))
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
C...Split total contribution into different colour flows just like
C...in g g -> g g (recalculate kinematics for massless partons).
          THP=-0.5D0*SH*(1D0-CTH)
          UHP=-0.5D0*SH*(1D0+CTH)
          FACGG1=(SH/THP)**2+2D0*SH/THP+3D0+2D0*THP/SH+(THP/SH)**2
          FACGG2=(UHP/SH)**2+2D0*UHP/SH+3D0+2D0*SH/UHP+(SH/UHP)**2
          FACGG3=(THP/UHP)**2+2D0*THP/UHP+3D0+2D0*UHP/THP+(UHP/THP)**2
          FACGGS=FACGG1+FACGG2+FACGG3
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
             NCHN=NCHN+1
             ISIG(NCHN,1)=21
             ISIG(NCHN,2)=21
             ISIG(NCHN,3)=1
             SIGH(NCHN)=FACQQG*PARP(IONIUM+2)*FACGG1/FACGGS
             NCHN=NCHN+1
             ISIG(NCHN,1)=21
             ISIG(NCHN,2)=21
             ISIG(NCHN,3)=2
             SIGH(NCHN)=FACQQG*PARP(IONIUM+2)*FACGG2/FACGGS
             NCHN=NCHN+1
             ISIG(NCHN,1)=21
             ISIG(NCHN,2)=21
             ISIG(NCHN,3)=3
             SIGH(NCHN)=FACQQG*PARP(IONIUM+2)*FACGG3/FACGGS
          ENDIF
 
        ELSEIF(ISUB.EQ.423) THEN
C...g + g -> QQ~[1S08] + g
          IF(MSTP(145).EQ.0) THEN
*            FACQQG=COMFAC*PARU(1)*AS**3*(5D0/16D0)*
*     &           (SHTH2*UH2+THUH2*SH2+UHSH2*TH2)/(SQMQQR*SH*TH*UH)*
*     &           (12D0*SQMQQ*SH*TH*UH+SHTH2**2+THUH2**2+UHSH2**2)/
*     &           (SHTH2*THUH2*UHSH2)
            FACQQG=COMFAC*PARU(1)*AS**3*(5D0/16D0)*SQMQQR*
     &            (UH2/(THUH2*UHSH2)+SH2/(SHTH2*UHSH2)+
     &            TH2/(SHTH2*THUH2))*
     &            (12D0+(SHTH2**2+THUH2**2+UHSH2**2)/(SQMQQ*SH*TH*UH))
          ELSE
            FA=PARU(1)*AS**3*(5D0/48D0)*SQMQQR*
     &            (UH2/(THUH2*UHSH2)+SH2/(SHTH2*UHSH2)+
     &            TH2/(SHTH2*THUH2))*
     &            (12D0+(SHTH2**2+THUH2**2+UHSH2**2)/(SQMQQ*SH*TH*UH))
            IF(MSTP(147).EQ.0) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=COMFAC*2D0*FA
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=0D0
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=0D0
            ENDIF
          ENDIF
C...Split total contribution into different colour flows just like
C...in g g -> g g (recalculate kinematics for massless partons).
          THP=-0.5D0*SH*(1D0-CTH)
          UHP=-0.5D0*SH*(1D0+CTH)
          FACGG1=(SH/THP)**2+2D0*SH/THP+3D0+2D0*THP/SH+(THP/SH)**2
          FACGG2=(UHP/SH)**2+2D0*UHP/SH+3D0+2D0*SH/UHP+(SH/UHP)**2
          FACGG3=(THP/UHP)**2+2D0*THP/UHP+3D0+2D0*UHP/THP+(UHP/THP)**2
          FACGGS=FACGG1+FACGG2+FACGG3
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
             NCHN=NCHN+1
             ISIG(NCHN,1)=21
             ISIG(NCHN,2)=21
             ISIG(NCHN,3)=1
             SIGH(NCHN)=FACQQG*PARP(IONIUM+3)*FACGG1/FACGGS
             NCHN=NCHN+1
             ISIG(NCHN,1)=21
             ISIG(NCHN,2)=21
             ISIG(NCHN,3)=2
             SIGH(NCHN)=FACQQG*PARP(IONIUM+3)*FACGG2/FACGGS
             NCHN=NCHN+1
             ISIG(NCHN,1)=21
             ISIG(NCHN,2)=21
             ISIG(NCHN,3)=3
             SIGH(NCHN)=FACQQG*PARP(IONIUM+3)*FACGG3/FACGGS
          ENDIF
 
        ELSEIF(ISUB.EQ.424) THEN
C...g + g -> QQ~[3PJ8] + g
          POLY=SH2+SH*TH+TH2
          IF(MSTP(145).EQ.0) THEN
            FACQQG=COMFAC*5D0*PARU(1)*AS**3*(3D0*SH*TH*SHTH*POLY**4
     &            -SQMQQ*POLY**2*(7D0*SH**6+36D0*SH**5*TH+45D0*SH**4*TH2
     &            +28D0*SH**3*TH**3+45D0*SH2*TH**4+36D0*SH*TH**5
     &            +7D0*TH**6)
     &            +SQMQQ**2*SHTH*(35D0*SH**8+169D0*SH**7*TH
     &            +299D0*SH**6*TH2+401D0*SH**5*TH**3+418D0*SH**4*TH**4
     &            +401D0*SH**3*TH**5+299D0*SH2*TH**6+169D0*SH*TH**7
     &            +35D0*TH**8)
     &            -SQMQQ**3*(84D0*SH**8+432D0*SH**7*TH+905D0*SH**6*TH2
     &            +1287D0*SH**5*TH**3+1436D0*SH**4*TH**4
     &            +1287D0*SH**3*TH**5+905D0*SH2*TH**6+432D0*SH*TH**7
     &            +84D0*TH**8)
     &            +SQMQQ**4*SHTH*(126D0*SH**6+451D0*SH**5*TH
     &            +677D0*SH**4*TH2+836D0*SH**3*TH**3+677D0*SH2*TH**4
     &            +451D0*SH*TH**5+126D0*TH**6)
     &            -3D0*SQMQQ**5*(42D0*SH**6+171D0*SH**5*TH
     &            +304D0*SH**4*TH2+362D0*SH**3*TH**3+304D0*SH2*TH**4
     &            +171D0*SH*TH**5+42D0*TH**6)
     &            +2D0*SQMQQ**6*SHTH*(42D0*SH**4+106D0*SH**3*TH
     &            +119D0*SH2*TH2+106D0*SH*TH**3+42D0*TH**4)
     &            -SQMQQ**7*(35D0*SH**4+99D0*SH**3*TH+120D0*SH2*TH2
     &            +99D0*SH*TH**3+35D0*TH**4)
     &            +7D0*SQMQQ**8*SHTH*POLY)/
     &            (SH*TH*UH*SQMQQR*SQMQQ*
     &            SHTH*SHTH2*THUH*THUH2*UHSH*UHSH2)
          ELSE
            FF=-5D0*PARU(1)*AS**3/(SH2*TH2*UH2
     &            *SQMQQR*SQMQQ*SHTH*SHTH2*THUH*THUH2*UHSH*UHSH2)
            AA=SH*TH*UH*(SH*TH*SHTH*POLY**4
     &           -SQMQQ*SHTH2*POLY**2*
     &           (SH**4+6D0*SH**3*TH-6D0*SH2*TH2+6D0*SH*TH**3+TH**4)
     &           +SQMQQ**2*SHTH*(5D0*SH**8+35D0*SH**7*TH+49D0*SH**6*TH2
     &           +57D0*SH**5*TH**3+46D0*SH**4*TH**4+57D0*SH**3*TH**5
     &           +49D0*SH2*TH**6+35D0*SH*TH**7+5D0*TH**8)
     &           -SQMQQ**3*(16D0*SH**8+104D0*SH**7*TH+215D0*SH**6*TH2
     &           +291D0*SH**5*TH**3+316D0*SH**4*TH**4+291D0*SH**3*TH**5
     &           +215D0*SH2*TH**6+104D0*SH*TH**7+16D0*TH**8)
     &           +SQMQQ**4*SHTH*(34D0*SH**6+145D0*SH**5*TH
     &           +211D0*SH**4*TH2+262D0*SH**3*TH**3+211D0*SH2*TH**4
     &           +145D0*SH*TH**5+34D0*TH**6)
     &           -SQMQQ**5*(44D0*SH**6+193D0*SH**5*TH+346D0*SH**4*TH2
     &           +410D0*SH**3*TH**3+346D0*SH2*TH**4+193D0*SH*TH**5
     &           +44D0*TH**6)
     &           +2D0*SQMQQ**6*SHTH*(17D0*SH**4+45D0*SH**3*TH
     &           +49D0*SH2*TH2+45D0*SH*TH**3+17D0*TH**4)
     &           -SQMQQ**7*(3D0*SH2+2D0*SH*TH+3D0*TH2)
     &           *(5D0*SH2+11D0*SH*TH+5D0*TH2)
     &           +3D0*SQMQQ**8*SHTH*POLY)
            BB=4D0*SHTH2*POLY**3
     &           *(SH**4+SH**3*TH-SH2*TH2+SH*TH**3+TH**4)
     &           -SQMQQ*SHTH*(20D0*SH**10+84D0*SH**9*TH+166D0*SH**8*TH2
     &           +231D0*SH**7*TH**3+250D0*SH**6*TH**4+250D0*SH**5*TH**5
     &           +250D0*SH**4*TH**6+231D0*SH**3*TH**7+166D0*SH2*TH**8
     &           +84D0*SH*TH**9+20D0*TH**10)
     &           +SQMQQ**2*SHTH2*(40D0*SH**8+86D0*SH**7*TH
     &           +66D0*SH**6*TH2+67D0*SH**5*TH**3+6D0*SH**4*TH**4
     &           +67D0*SH**3*TH**5+66D0*SH2*TH**6+86D0*SH*TH**7
     &           +40D0*TH**8)
     &           -SQMQQ**3*SHTH*(40D0*SH**8+57D0*SH**7*TH
     &           -110D0*SH**6*TH2-263D0*SH**5*TH**3-384D0*SH**4*TH**4
     &           -263D0*SH**3*TH**5-110D0*SH2*TH**6+57D0*SH*TH**7
     &           +40D0*TH**8)
     &           +SQMQQ**4*(20D0*SH**8-33D0*SH**7*TH-368D0*SH**6*TH2
     &           -751D0*SH**5*TH**3-920D0*SH**4*TH**4-751D0*SH**3*TH**5
     &           -368D0*SH2*TH**6-33D0*SH*TH**7+20D0*TH**8)
     &           -SQMQQ**5*SHTH*(4D0*SH**6-81D0*SH**5*TH-242D0*SH**4*TH2
     &           -250D0*SH**3*TH**3-242D0*SH2*TH**4-81D0*SH*TH**5
     &           +4D0*TH**6)
     &           -SQMQQ**6*SH*TH*(41D0*SH**4+120D0*SH**3*TH
     &           +142D0*SH2*TH2+120D0*SH*TH**3+41D0*TH**4)
     &           +8D0*SQMQQ**7*SH*TH*SHTH*POLY
            CC=4D0*TH2*POLY**3
     &           *(-SH**4-2D0*SH**3*TH+2D0*SH2*TH2+3D0*SH*TH**3+TH**4)
     &           -SQMQQ*TH2*(-20D0*SH**9-56D0*SH**8*TH-24D0*SH**7*TH2
     &           +147D0*SH**6*TH**3+409D0*SH**5*TH**4+599D0*SH**4*TH**5
     &           +571D0*SH**3*TH**6+370D0*SH2*TH**7+148D0*SH*TH**8
     &           +28D0*TH**9)
     &           +SQMQQ**2*(4D0*SH**10+20D0*SH**9*TH-16D0*SH**8*TH2
     &           -48D0*SH**7*TH**3+150D0*SH**6*TH**4+611D0*SH**5*TH**5
     &           +1060D0*SH**4*TH**6+1155D0*SH**3*TH**7+854D0*SH2*TH**8
     &           +394D0*SH*TH**9+84D0*TH**10)
     &           -SQMQQ**3*SHTH*(20D0*SH**8+68D0*SH**7*TH-20D0*SH**6*TH2
     &           +32D0*SH**5*TH**3+286D0*SH**4*TH**4+577D0*SH**3*TH**5
     &           +618D0*SH2*TH**6+443D0*SH*TH**7+140D0*TH**8)
     &           +SQMQQ**4*(40D0*SH**8+152D0*SH**7*TH+94D0*SH**6*TH2
     &           +38D0*SH**5*TH**3+290D0*SH**4*TH**4+631D0*SH**3*TH**5
     &           +738D0*SH2*TH**6+513D0*SH*TH**7+140D0*TH**8)
     &           -SQMQQ**5*(40D0*SH**7+129D0*SH**6*TH+53D0*SH**5*TH2
     &           +7D0*SH**4*TH**3+129D0*SH**3*TH**4+264D0*SH2*TH**5
     &           +266D0*SH*TH**6+84D0*TH**7)
     &           +SQMQQ**6*(20D0*SH**6+55D0*SH**5*TH+2D0*SH**4*TH2
     &           -15D0*SH**3*TH**3+30D0*SH2*TH**4+76D0*SH*TH**5
     &           +28D0*TH**6)
     &           -SQMQQ**7*SHTH*(4D0*SH**4+7D0*SH**3*TH-14D0*SH2*TH2
     &           +7D0*SH*TH**3+4*TH**4)
     &           +SQMQQ**8*SH*(SH-TH)**2*TH
            DD=2D0*TH2*SHTH2*POLY**3
     &           *(-SH2+2*SH*TH+2*TH2)
     &           +SQMQQ*(4D0*SH**11+22D0*SH**10*TH+70D0*SH**9*TH2
     &           +115D0*SH**8*TH**3+71D0*SH**7*TH**4-119D0*SH**6*TH**5
     &           -381D0*SH**5*TH**6-552D0*SH**4*TH**7-512D0*SH**3*TH**8
     &           -320D0*SH2*TH**9-126D0*SH*TH**10-24D0*TH**11)
     &           -SQMQQ**2*SHTH*(20D0*SH**9+84D0*SH**8*TH
     &           +212D0*SH**7*TH2+247D0*SH**6*TH**3+105D0*SH**5*TH**4
     &           -178D0*SH**4*TH**5-380D0*SH**3*TH**6-364D0*SH2*TH**7
     &           -210D0*SH*TH**8-60D0*TH**9)
     &           +SQMQQ**3*SHTH*(40D0*SH**8+159D0*SH**7*TH
     &           +374D0*SH**6*TH2+404D0*SH**5*TH**3+192D0*SH**4*TH**4
     &           -141D0*SH**3*TH**5-264D0*SH2*TH**6-216D0*SH*TH**7
     &           -80D0*TH**8)
     &           -SQMQQ**4*(40D0*SH**8+197D0*SH**7*TH+506D0*SH**6*TH2
     &           +672D0*SH**5*TH**3+460D0*SH**4*TH**4+79D0*SH**3*TH**5
     &           -138D0*SH2*TH**6-164D0*SH*TH**7-60D0*TH**8)
     &           +SQMQQ**5*(20D0*SH**7+107D0*SH**6*TH+267D0*SH**5*TH2
     &           +307D0*SH**4*TH**3+185D0*SH**3*TH**4+56D0*SH2*TH**5
     &           -30D0*SH*TH**6-24D0*TH**7)
     &           -SQMQQ**6*(4D0*SH**6+31D0*SH**5*TH+74D0*SH**4*TH2
     &           +71D0*SH**3*TH**3+46D0*SH2*TH**4+10D0*SH*TH**5
     &           -4D0*TH**6)
     &           +4D0*SQMQQ**7*SH*TH*SHTH*POLY
            IF(MSTP(147).EQ.0) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=2D0*(-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11)))
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K10+CC*EL1K21*EL2K20
     &              +DD*(EL1K11*EL2K20+EL1K21*EL2K10))
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
C...Split total contribution into different colour flows just like
C...in g g -> g g (recalculate kinematics for massless partons).
          THP=-0.5D0*SH*(1D0-CTH)
          UHP=-0.5D0*SH*(1D0+CTH)
          FACGG1=(SH/THP)**2+2D0*SH/THP+3D0+2D0*THP/SH+(THP/SH)**2
          FACGG2=(UHP/SH)**2+2D0*UHP/SH+3D0+2D0*SH/UHP+(SH/UHP)**2
          FACGG3=(THP/UHP)**2+2D0*THP/UHP+3D0+2D0*UHP/THP+(UHP/THP)**2
          FACGGS=FACGG1+FACGG2+FACGG3
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
             NCHN=NCHN+1
             ISIG(NCHN,1)=21
             ISIG(NCHN,2)=21
             ISIG(NCHN,3)=1
             SIGH(NCHN)=FACQQG*PARP(IONIUM+4)*FACGG1/FACGGS
             NCHN=NCHN+1
             ISIG(NCHN,1)=21
             ISIG(NCHN,2)=21
             ISIG(NCHN,3)=2
             SIGH(NCHN)=FACQQG*PARP(IONIUM+4)*FACGG2/FACGGS
             NCHN=NCHN+1
             ISIG(NCHN,1)=21
             ISIG(NCHN,2)=21
             ISIG(NCHN,3)=3
             SIGH(NCHN)=FACQQG*PARP(IONIUM+4)*FACGG3/FACGGS
          ENDIF
 
        ELSEIF(ISUB.EQ.425) THEN
C...q + g -> q + QQ~[3S18]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=-COMFAC*PARU(1)*AS**3*(1D0/27D0)*
     &            (4D0*(SH2+UH2)-SH*UH)*(SHTH2+THUH2)/
     &            (SQMQQ*SQMQQR*SH*UH*UHSH2)
          ELSE
            FF=PARU(1)*AS**3*(4D0*(SH2+UH2)-SH*UH)/
     &            (54D0*SQMQQ*SQMQQR*SH*UH*UHSH2)
            AA=SHTH2+THUH2
            BB=4D0
            CC=8D0
            DD=4D0
            IF(MSTP(147).EQ.0) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=2D0*(-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11)))
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K10+CC*EL1K21*EL2K20
     &              +DD*(EL1K11*EL2K20+EL1K21*EL2K10))
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
C...Split total contribution into different colour flows just like
C...in ISUB.EQ.28 [f + g -> f + g (q + g -> q + g only)]
C...(recalculate kinematics for massless partons).
          THP=-0.5D0*SH*(1D0-CTH)
          UHP=-0.5D0*SH*(1D0+CTH)
          FACQG1=9D0/4D0*(UHP/THP)**2-UHP/SH
          FACQG2=9D0/4D0*(SH/THP)**2-SH/UHP
          FACQGS=FACQG1+FACQG2
          DO 2442 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 2442
            DO 2441 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 2441
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 2441
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQQG*PARP(IONIUM+2)*FACQG1/FACQGS
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=2
              SIGH(NCHN)=FACQQG*PARP(IONIUM+2)*FACQG2/FACQGS
 2441       CONTINUE
 2442     CONTINUE
 
        ELSEIF(ISUB.EQ.426) THEN
C...q + g -> q + QQ~[1S08]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=-COMFAC*PARU(1)*AS**3*(5D0/18D0)*
     &            (SH2+UH2)/(SQMQQR*TH*UHSH2)
          ELSE
            FA=-PARU(1)*AS**3*(5D0/54D0)*(SH2+UH2)/(SQMQQR*TH*UHSH2)
            IF(MSTP(147).EQ.0) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=COMFAC*2D0*FA
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=0D0
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=0D0
            ENDIF
          ENDIF
C...Split total contribution into different colour flows just like
C...in ISUB.EQ.28 [f + g -> f + g (q + g -> q + g only)]
C...(recalculate kinematics for massless partons).
          THP=-0.5D0*SH*(1D0-CTH)
          UHP=-0.5D0*SH*(1D0+CTH)
          FACQG1=9D0/4D0*(UHP/THP)**2-UHP/SH
          FACQG2=9D0/4D0*(SH/THP)**2-SH/UHP
          FACQGS=FACQG1+FACQG2
          DO 2444 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 2444
            DO 2443 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 2443
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 2443
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQQG*PARP(IONIUM+3)*FACQG1/FACQGS
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=2
              SIGH(NCHN)=FACQQG*PARP(IONIUM+3)*FACQG2/FACQGS
 2443       CONTINUE
 2444     CONTINUE
 
        ELSEIF(ISUB.EQ.427) THEN
C...q + g -> q + QQ~[3PJ8]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=-COMFAC*PARU(1)*AS**3*(10D0/9D0)*
     &            ((7D0*UHSH+8D0*TH)*(SH2+UH2)
     &            +4D0*TH*(2D0*SQMQQ**2-SHTH2-THUH2))/
     &            (SQMQQ*SQMQQR*TH*UHSH2*UHSH)
          ELSE
            FF=10D0*PARU(1)*AS**3/
     &            (9D0*SQMQQ*SQMQQR*TH2*UHSH2*UHSH)
            AA=TH*UHSH*(2D0*SQMQQ**2+SHTH2+THUH2)
            BB=8D0*(SHTH2+TH*UH)
            CC=8D0*UHSH*(SHTH+THUH)
            DD=4D0*(2D0*SQMQQ*SH+TH*UHSH)
            IF(MSTP(147).EQ.0) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=2D0*(-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11)))
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K10+CC*EL1K21*EL2K20
     &              +DD*(EL1K11*EL2K20+EL1K21*EL2K10))
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
C...Split total contribution into different colour flows just like
C...in ISUB.EQ.28 [f + g -> f + g (q + g -> q + g only)]
C...(recalculate kinematics for massless partons).
          THP=-0.5D0*SH*(1D0-CTH)
          UHP=-0.5D0*SH*(1D0+CTH)
          FACQG1=9D0/4D0*(UHP/THP)**2-UHP/SH
          FACQG2=9D0/4D0*(SH/THP)**2-SH/UHP
          FACQGS=FACQG1+FACQG2
          DO 2446 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 2446
            DO 2445 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 2445
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 2445
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQQG*PARP(IONIUM+4)*FACQG1/FACQGS
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=2
              SIGH(NCHN)=FACQQG*PARP(IONIUM+4)*FACQG2/FACQGS
 2445       CONTINUE
 2446     CONTINUE
 
        ELSEIF(ISUB.EQ.428) THEN
C...q + q~ -> g + QQ~[3S18]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=COMFAC*PARU(1)*AS**3*(8D0/81D0)*
     &            (4D0*(TH2+UH2)-TH*UH)*(SHTH2+UHSH2)/
     &            (SQMQQ*SQMQQR*TH*UH*THUH2)
          ELSE
            FF=-4D0*PARU(1)*AS**3*(4D0*(TH2+UH2)-TH*UH)/
     &            (81D0*SQMQQ*SQMQQR*TH*UH*THUH2)
            AA=SHTH2+UHSH2
            BB=4D0
            CC=4D0
            DD=0D0
            IF(MSTP(147).EQ.0) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=2D0*(-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11)))
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K10+CC*EL1K21*EL2K20
     &              +DD*(EL1K11*EL2K20+EL1K21*EL2K10))
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
C...Split total contribution into different colour flows just like
C...in ISUB.EQ.13 [f + fbar -> g + g (q + qbar -> g + g only)]
C...(recalculate kinematics for massless partons).
          THP=-0.5D0*SH*(1D0-CTH)
          UHP=-0.5D0*SH*(1D0+CTH)
          FACGG1=UH/TH-9D0/4D0*UH2/SH2
          FACGG2=TH/UH-9D0/4D0*TH2/SH2
          FACGGS=FACGG1+FACGG2
          DO 2447 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &            KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 2447
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG*PARP(IONIUM+2)*FACGG1/FACGGS
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=2
            SIGH(NCHN)=FACQQG*PARP(IONIUM+2)*FACGG2/FACGGS
 2447     CONTINUE
 
        ELSEIF(ISUB.EQ.429) THEN
C...q + q~ -> g + QQ~[1S08]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=COMFAC*PARU(1)*AS**3*(20D0/27D0)*
     &            (TH2+UH2)/(SQMQQR*SH*THUH2)
          ELSE
            FA=PARU(1)*AS**3*(20D0/81D0)*(TH2+UH2)/(SQMQQR*SH*THUH2)
            IF(MSTP(147).EQ.0) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=COMFAC*2D0*FA
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=0D0
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=0D0
            ENDIF
          ENDIF
C...Split total contribution into different colour flows just like
C...in ISUB.EQ.13 [f + fbar -> g + g (q + qbar -> g + g only)]
C...(recalculate kinematics for massless partons).
          THP=-0.5D0*SH*(1D0-CTH)
          UHP=-0.5D0*SH*(1D0+CTH)
          FACGG1=UH/TH-9D0/4D0*UH2/SH2
          FACGG2=TH/UH-9D0/4D0*TH2/SH2
          FACGGS=FACGG1+FACGG2
          DO 2448 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &            KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 2448
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG*PARP(IONIUM+3)*FACGG1/FACGGS
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=2
            SIGH(NCHN)=FACQQG*PARP(IONIUM+3)*FACGG2/FACGGS
 2448     CONTINUE
 
        ELSEIF(ISUB.EQ.430) THEN
C...q + q~ -> g + QQ~[3PJ8]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=COMFAC*PARU(1)*AS**3*(80D0/27D0)*
     &            ((7D0*THUH+8D0*SH)*(TH2+UH2)
     &            +4D0*SH*(2D0*SQMQQ**2-SHTH2-UHSH2))/
     &            (SQMQQ*SQMQQR*SH*THUH2*THUH)
          ELSE
            FF=-80D0*PARU(1)*AS**3/(27D0*SQMQQ*SQMQQR*SH2*THUH2*THUH)
            AA=SH*THUH*(2D0*SQMQQ**2+SHTH2+UHSH2)
            BB=8D0*(UHSH2+SH*TH)
            CC=8D0*(SHTH2+SH*UH)
            DD=4D0*(SHTH2+UHSH2+SH*SQMQQ-SQMQQ**2)
            IF(MSTP(147).EQ.0) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=2D0*(-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11)))
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K10*EL2K10+CC*EL1K20*EL2K20
     &              +DD*(EL1K10*EL2K20+EL1K20*EL2K10))
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=-AA+SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K10+CC*EL1K21*EL2K20
     &              +DD*(EL1K11*EL2K20+EL1K21*EL2K10))
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=SQMQQ*(BB*EL1K11*EL2K11+CC*EL1K21*EL2K21
     &              +DD*(EL1K11*EL2K21+EL1K21*EL2K11))
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
C...Split total contribution into different colour flows just like
C...in ISUB.EQ.13 [f + fbar -> g + g (q + qbar -> g + g only)]
C...(recalculate kinematics for massless partons).
          THP=-0.5D0*SH*(1D0-CTH)
          UHP=-0.5D0*SH*(1D0+CTH)
          FACGG1=UH/TH-9D0/4D0*UH2/SH2
          FACGG2=TH/UH-9D0/4D0*TH2/SH2
          FACGGS=FACGG1+FACGG2
          DO 2449 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &            KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 2449
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG*PARP(IONIUM+4)*FACGG1/FACGGS
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=2
            SIGH(NCHN)=FACQQG*PARP(IONIUM+4)*FACGG2/FACGGS
 2449     CONTINUE
 
        ELSEIF(ISUB.EQ.431) THEN
C...g + g -> QQ~[3P01] + g
          PGTW=(SH*TH+TH*UH+UH*SH)/SH2
          QGTW=(SH*TH*UH)/SH**3
          RGTW=SQMQQ/SH
          IF(MSTP(145).EQ.0) THEN
            FACQQG=COMFAC*PARU(1)*AS**3*8D0/(9D0*SQMQQR*SH)*
     &            (9D0*RGTW**2*PGTW**4*
     &            (RGTW**4-2D0*RGTW**2*PGTW+PGTW**2)
     &            -6D0*RGTW*PGTW**3*QGTW*
     &            (2D0*RGTW**4-5D0*RGTW**2*PGTW+PGTW**2)
     &            -PGTW**2*QGTW**2*(RGTW**4+2D0*RGTW**2*PGTW-PGTW**2)
     &            +2D0*RGTW*PGTW*QGTW**3*(RGTW**2-PGTW)
     &            +6D0*RGTW**2*QGTW**4)/(QGTW*(QGTW-RGTW*PGTW)**4)
          ELSE
            FC1=PARU(1)*AS**3*8D0/(27D0*SQMQQR*SH)*
     &            (9D0*RGTW**2*PGTW**4*
     &            (RGTW**4-2D0*RGTW**2*PGTW+PGTW**2)
     &            -6D0*RGTW*PGTW**3*QGTW*
     &            (2D0*RGTW**4-5D0*RGTW**2*PGTW+PGTW**2)
     &            -PGTW**2*QGTW**2*(RGTW**4+2D0*RGTW**2*PGTW-PGTW**2)
     &            +2D0*RGTW*PGTW*QGTW**3*(RGTW**2-PGTW)
     &            +6D0*RGTW**2*QGTW**4)/(QGTW*(QGTW-RGTW*PGTW)**4)
            IF(MSTP(147).EQ.0) THEN
               FACQQG=COMFAC*FC1
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=COMFAC*2D0*FC1
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=COMFAC*FC1
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=COMFAC*FC1
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=0D0
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=0D0
            ENDIF
          ENDIF
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG*PARP(IONIUM+5)
          ENDIF
 
        ELSEIF(ISUB.EQ.432) THEN
C...g + g -> QQ~[3P11] + g
          PGTW=(SH*TH+TH*UH+UH*SH)/SH2
          QGTW=(SH*TH*UH)/SH**3
          RGTW=SQMQQ/SH
          IF(MSTP(145).EQ.0) THEN
            FACQQG=COMFAC*PARU(1)*AS**3*8D0/(3D0*SQMQQR*SH)*
     &            PGTW**2*(RGTW*PGTW**2*(RGTW**2-4D0*PGTW)
     &            +2D0*QGTW*(-RGTW**4+5D0*RGTW**2*PGTW+PGTW**2)
     &            -15D0*RGTW*QGTW**2)/(QGTW-RGTW*PGTW)**4
          ELSE
            FF=4D0/3D0*PARU(1)*AS**3*SQMQQR/SHTH2**2/THUH2**2/UHSH2**2
            C1=(4D0*PGTW**5+23D0*PGTW**2*QGTW**2
     &            +(-14D0*PGTW**3*QGTW+3D0*QGTW**3)*RGTW
     &            -(PGTW**4+2D0*PGTW*QGTW**2)*RGTW**2
     &            +3D0*PGTW**2*QGTW*RGTW**3)*SH2**5
            C2=2D0*SHTH2*(SH2*THUH*(SH*THUH*(SH-TH)*(SH-UH)
     &            -TH*UH*(TH-UH)**2)+SH2**2*(TH-UH)*(TH2+UH2-SH*THUH)
     &            *(PGTW**2-QGTW*(SH+2D0*UH)/SH))
            C3=2D0*UHSH2*(SH2*THUH*(SH*THUH*(SH-TH)*(SH-UH)
     &            -TH*UH*(TH-UH)**2)-SH2**2*(TH-UH)*(TH2+UH2-SH*THUH)
     &            *(PGTW**2-QGTW*(SH+2D0*TH)/SH))
            C4=-4D0*THUH*(TH-UH)**2*
     &            (TH**3*UH**3+SH2**2*(2D0*TH+UH)*(TH+2D0*UH)
     &            -SH2*TH*UH*(TH2+UH2))
     &            +4D0*THUH2*(SH**3*(SH2**2+TH2**2+UH2**2)
     &            -SH*TH*UH*(SH2**2+TH*UH*(TH2-3D0*TH*UH+UH2)
     &            +SH2*(5D0*THUH2-17D0*TH*UH)))
            IF(MSTP(147).EQ.0) THEN
               FACQQG=-C1+C2*EL1K10*EL2K10+C3*EL1K20*EL2K20
     &              +C4*(EL1K10*EL2K20+EL1K20*EL2K10)/2D0
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=2D0*(-C1+C2*EL1K11*EL2K11+C3*EL1K21*EL2K21
     &              +C4*(EL1K11*EL2K21+EL1K21*EL2K11)/2D0)
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=-C1+C2*EL1K10*EL2K10+C3*EL1K20*EL2K20
     &              +C4*(EL1K10*EL2K20+EL1K20*EL2K10)/2D0
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=-C1+C2*EL1K11*EL2K11+C3*EL1K21*EL2K21
     &              +C4*(EL1K11*EL2K21+EL1K21*EL2K11)/2D0
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=C2*EL1K11*EL2K10+C3*EL1K21*EL2K20
     &              +C4*(EL1K11*EL2K20+EL1K21*EL2K10)/2D0
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=C2*EL1K11*EL2K11+C3*EL1K21*EL2K21
     &              +C4*(EL1K11*EL2K21+EL1K21*EL2K11)/2D0
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG*PARP(IONIUM+5)
          ENDIF
 
        ELSEIF(ISUB.EQ.433) THEN
C...g + g -> QQ~[3P21] + g
          PGTW=(SH*TH+TH*UH+UH*SH)/SH2
          QGTW=(SH*TH*UH)/SH**3
          RGTW=SQMQQ/SH
          IF(MSTP(145).EQ.0) THEN
            FACQQG=COMFAC*PARU(1)*AS**3*8D0/(9D0*SQMQQR*SH)*
     &            (12D0*RGTW**2*PGTW**4*
     &            (RGTW**4-2D0*RGTW**2*PGTW+PGTW**2)
     &            -3D0*RGTW*PGTW**3*QGTW*
     &            (8D0*RGTW**4-RGTW**2*PGTW+4D0*PGTW**2)
     &            +2D0*PGTW**2*QGTW**2*
     &            (-7D0*RGTW**4+43D0*RGTW**2*PGTW+PGTW**2)
     &            +RGTW*PGTW*QGTW**3*(16D0*RGTW**2-61D0*PGTW)
     &            +12D0*RGTW**2*QGTW**4)/(QGTW*(QGTW-RGTW*PGTW)**4)
          ELSE
            FF=(16D0*PARU(1)*AS**3*SQMQQ*SQMQQR)/
     &            (3D0*SH2*TH2*UH2*SHTH2**2*THUH2**2*UHSH2**2)
            C1=PGTW**2*QGTW*(PGTW*RGTW-QGTW)**2*(RGTW**2-2D0*PGTW)
     &            *SH*SH2**7
            C2=2D0*SHTH2*(-SH2**3*TH2**3-SH**5*TH**5*UH*SHTH
     &            +SH2**2*TH2**2*UH2*(8D0*SHTH2-5D0*SH*TH)
     &            +SH**3*TH**3*UH**3*SHTH*(17D0*SHTH2-2D0*SH*TH)
     &            +SH2*TH2*UH2**2*(105D0*SH2*TH2+64D0*SH*TH*(SH2+TH2)
     &            +10D0*(SH2**2+TH2**2))
     &            +SH2*TH2*UH**5*SHTH*(32D0*SHTH2+7D0*SH*TH)
     &            -UH2**3*(SH2**3-87D0*SH**3*TH**3+TH2**3
     &            -45D0*SH2*TH2*(SH2+TH2)-5D0*SH*TH*(SH2**2+TH2**2))
     &            +SH*TH*UH**7*SHTH*(7D0*SHTH2+12D0*SH*TH)
     &            +4D0*SH*TH*UH2**4*SHTH2)
            C3=2D0*UHSH2*(-SH2**3*UH2**3-SH**5*UH**5*TH*UHSH
     &            +SH2**2*UH2**2*TH2*(8D0*UHSH2-5D0*SH*UH)
     &            +SH**3*UH**3*TH**3*UHSH*(17D0*UHSH2-2D0*SH*UH)
     &            +SH2*UH2*TH2**2*(105D0*SH2*UH2+64D0*SH*UH*(SH2+UH2)
     &            +10D0*(SH2**2+UH2**2))
     &            +SH2*UH2*TH**5*UHSH*(32D0*UHSH2+7D0*SH*UH)
     &            -TH2**3*(SH2**3-87D0*SH**3*UH**3+UH2**3
     &            -45D0*SH2*UH2*(SH2+UH2)-5D0*SH*UH*(SH2**2+UH2**2))
     &            +SH*UH*TH**7*UHSH*(7D0*UHSH2+12D0*SH*UH)
     &            +4D0*SH*UH*TH2**4*UHSH2)
            C4=-2D0*SHTH*UHSH*(-2D0*TH2**3*UH2**3
     &            -SH**5*TH2*UH2*THUH*(5D0*TH+3D0*UH)*(3D0*TH+5D0*UH)
     &            +SH2**3*(2D0*TH+UH)*(TH+2D0*UH)*(TH2-UH2)**2
     &            -SH*TH2**2*UH2**2*THUH*(5D0*THUH2-4D0*TH*UH)
     &            -SH2*TH**3*UH**3*THUH2*(13D0*THUH2-16D0*TH*UH)
     &            -SH**3*TH2*UH2*(92D0*TH2*UH2*THUH
     &            +53D0*TH*UH*(TH**3+UH**3)+11D0*(TH**5+UH**5))
     &            -SH2**2*TH*UH*(114D0*TH**3*UH**3
     &            +83D0*TH2*UH2*(TH2+UH2)+28D0*TH*UH*(TH2**2+UH2**2)
     &            +3D0*(TH2**3+UH2**3)))
            C5=4D0*SH*TH*UH2*SHTH2*(2D0*SH*TH+SH*UH+TH*UH)**2
     &            *(2D0*UH*SQMQQ**2+SHTH*(SH*TH-UH2))
            C6=4D0*SH*UH*TH2*UHSH2*(2D0*SH*UH+SH*TH+TH*UH)**2
     &            *(2D0*TH*SQMQQ**2+UHSH*(SH*UH-TH2))
            C7=4D0*SH*TH*UH2*SHTH*(SH2**2*TH**3*(11D0*SH+16D0*TH)
     &            +SH**3*TH2*UH*(31D0*SH2+83D0*SH*TH+61D0*TH2)
     &            +SH2*TH*UH2*(19D0*SH**3+110D0*SH2*TH+156D0*SH*TH2+
     &            82D0*TH**3)
     &            +SH*TH*UH**3*(43D0*SH**3+132D0*SH2*TH+124D0*SH*TH2
     &            +45D0*TH**3)
     &            +TH*UH2**2*(37D0*SH**3+68D0*SH2*TH+43D0*SH*TH2+
     &            8D0*TH**3)
     &            +TH*UH**5*(11D0*SH2+13D0*SH*TH+5D0*TH2)
     &            +SH**3*UH**3*(3D0*UHSH2-2D0*SH*UH)
     &            +TH**5*UHSH*(5D0*UHSH2+2D0*SH*UH))
            C8=4D0*SH*UH*TH2*UHSH*(SH2**2*UH**3*(11D0*SH+16D0*UH)
     &            +SH**3*UH2*TH*(31D0*SH2+83D0*SH*UH+61D0*UH2)
     &            +SH2*UH*TH2*(19D0*SH**3+110D0*SH2*UH+156D0*SH*UH2+
     &            82D0*UH**3)
     &            +SH*UH*TH**3*(43D0*SH**3+132D0*SH2*UH+124D0*SH*UH2
     &            +45D0*UH**3)
     &            +UH*TH2**2*(37D0*SH**3+68D0*SH2*UH+43D0*SH*UH2+
     &            8D0*UH**3)
     &            +UH*TH**5*(11D0*SH2+13D0*SH*UH+5D0*UH2)
     &            +SH**3*TH**3*(3D0*SHTH2-2D0*SH*TH)
     &            +UH**5*SHTH*(5D0*SHTH2+2D0*SH*TH))
            C9=4D0*SHTH*UHSH*(2D0*TH**5*UH**5*THUH
     &            +4D0*SH*TH2**2*UH2**2*THUH2
     &            -SH2*TH**3*UH**3*THUH*(TH2+UH2)
     &            -2D0*SH**3*TH2*UH2*(THUH2**2+2D0*TH*UH*THUH2-TH2*UH2)
     &            +SH2**2*TH*UH*THUH*(-TH*UH*THUH2+3D0*(TH2**2+UH2**2))
     &            +SH**5*(4D0*TH2*UH2*(THUH2-TH*UH)
     &            +5D0*TH*UH*(TH2**2+UH2**2)+2D0*(TH2**3+UH2**3)))
            C0=-4D0*(2D0*TH2**3*UH2**3*SQMQQ
     &            -SH2*TH2**2*UH2**2*THUH*(19D0*THUH2-4D0*TH*UH)
     &            -SH**3*TH**3*UH**3*THUH2*(32D0*THUH2+29D0*TH*UH)
     &            -SH2**2*TH2*UH2*THUH*(264D0*TH2*UH2
     &            +136D0*TH*UH*(TH2+UH2)+15D0*(TH2**2+UH2**2))
     &            +SH**5*TH*UH*(-428D0*TH**3*UH**3
     &            -256D0*TH2*UH2*(TH2+UH2)-43D0*TH*UH*(TH2**2+UH2**2)
     &            +2D0*(TH2**3+UH2**3))
     &            +SH**7*(-46D0*TH**3*UH**3-21D0*TH2*UH2*(TH2+UH2)
     &            +2D0*TH*UH*(TH2**2+UH2**2)+2D0*(TH2**3+UH2**3))
     &            +SH2**3*THUH*(-134*TH**3*UH**3-53D0*TH2*UH2*(TH2+UH2)
     &            +4D0*TH*UH*(TH2**2+UH2**2)+2D0*(TH2**3+UH2**3)))
            IF(MSTP(147).EQ.0) THEN
               FACQQG=1D0/3D0*(C1*3D0
     &              -C2*(2D0*EL1K10*EL2K10+EL1K11*EL2K11)
     &              -C3*(2D0*EL1K20*EL2K20+EL1K21*EL2K21)
     &              -C4*(2D0*EL1K10*EL2K20+EL1K11*EL2K21)
     &              +C5*2D0*(EL1K10*EL2K10-EL1K11*EL2K11)**2
     &              +C6*2D0*(EL1K20*EL2K20-EL1K21*EL2K21)**2
     &              +C7*2D0*(EL1K10*EL2K10-EL1K11*EL2K11)
     &                     *(EL1K10*EL2K20-EL1K11*EL2K21)
     &              +C8*2D0*(EL1K20*EL2K20-EL1K21*EL2K21)
     &                     *(EL1K10*EL2K20-EL1K11*EL2K21)
     &              +C9*2D0*(EL1K10*EL2K10-EL1K11*EL2K11)
     &                     *(EL1K20*EL2K20-EL1K21*EL2K21)
     &              +C0*2D0*(EL1K10*EL2K20-EL1K11*EL2K21)**2)
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=C1*2D0
     &              -C2*(EL1K10*EL2K10+EL1K11*EL2K11)
     &              -C3*(EL1K20*EL2K20+EL1K21*EL2K21)
     &              -C4*(EL1K10*EL2K20+EL1K11*EL2K21)
     &              +C5*4D0*EL1K10*EL2K10*EL1K11*EL2K11
     &              +C6*4D0*EL1K20*EL2K20*EL1K21*EL2K21
     &              +C7*2D0*(EL1K10*EL2K10*EL1K11*EL2K21
     &                      +EL1K10*EL2K20*EL1K11*EL2K11)
     &              +C8*2D0*(EL1K20*EL2K20*EL1K11*EL2K21
     &                      +EL1K10*EL2K20*EL1K21*EL2K21)
     &              +C9*4D0*EL1K10*EL2K20*EL1K11*EL2K21
     &              +C0*(EL1K10*EL2K10*EL1K21*EL2K21
     &              +2D0*EL1K10*EL2K20*EL1K11*EL2K21
     &                  +EL1K20*EL2K20*EL1K11*EL2K11)
            ELSEIF(MSTP(147).EQ.2) THEN
               FACQQG=2D0*(C1
     &              -C2*EL1K11*EL2K11
     &              -C3*EL1K21*EL2K21
     &              -C4*EL1K11*EL2K21
     &              +C5*(EL1K11*EL2K11)**2
     &              +C6*(EL1K21*EL2K21)**2
     &              +C7*EL1K11*EL2K11*EL1K11*EL2K21
     &              +C8*EL1K21*EL2K21*EL1K11*EL2K21
     &              +(C9+C0)*(EL1K11*EL2K21)**2)
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG*PARP(IONIUM+5)
          ENDIF
 
        ELSEIF(ISUB.EQ.434) THEN
C...q + g -> q + QQ~[3P01]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=-COMFAC*PARU(1)*AS**3*(16D0/81D0)*
     &            (TH-3D0*SQMQQ)**2*(SH2+UH2)/(SQMQQR*TH*UHSH2**2)
          ELSE
            FA=-PARU(1)*AS**3*(16D0/243D0)*
     &            (TH-3D0*SQMQQ)**2*(SH2+UH2)/(SQMQQR*TH*UHSH2**2)
            IF(MSTP(147).EQ.0) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=COMFAC*2D0*FA
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=0D0
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=0D0
            ENDIF
          ENDIF
          DO 2452 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 2452
            DO 2451 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 2451
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 2451
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQQG*PARP(IONIUM+5)
 2451       CONTINUE
 2452     CONTINUE
 
        ELSEIF(ISUB.EQ.435) THEN
C...q + g -> q + QQ~[3P11]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=-COMFAC*PARU(1)*AS**3*(32D0/27D0)*
     &            (4D0*SQMQQ*SH*UH+TH*(SH2+UH2))/(SQMQQR*UHSH2**2)
          ELSE
            FF=(64D0*PARU(1)*AS**3*SQMQQR)/(27D0*UHSH2**2)
            C1=SH*UH
            C2=2D0*SH
            C3=0D0
            C4=2D0*(SH-UH)
            IF(MSTP(147).EQ.0) THEN
               FACQQG=-C1+C2*EL1K10*EL2K10+C3*EL1K20*EL2K20
     &              +C4*(EL1K10*EL2K20+EL1K20*EL2K10)/2D0
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=2D0*(-C1+C2*EL1K11*EL2K11+C3*EL1K21*EL2K21
     &              +C4*(EL1K11*EL2K21+EL1K21*EL2K11)/2D0)
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=-C1+C2*EL1K10*EL2K10+C3*EL1K20*EL2K20
     &              +C4*(EL1K10*EL2K20+EL1K20*EL2K10)/2D0
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=-C1+C2*EL1K11*EL2K11+C3*EL1K21*EL2K21
     &              +C4*(EL1K11*EL2K21+EL1K21*EL2K11)/2D0
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=C2*EL1K11*EL2K10+C3*EL1K21*EL2K20
     &              +C4*(EL1K11*EL2K20+EL1K21*EL2K10)/2D0
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=C2*EL1K11*EL2K11+C3*EL1K21*EL2K21
     &              +C4*(EL1K11*EL2K21+EL1K21*EL2K11)/2D0
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
          DO 2454 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 2454
            DO 2453 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 2453
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 2453
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQQG*PARP(IONIUM+5)
 2453       CONTINUE
 2454     CONTINUE
 
        ELSEIF(ISUB.EQ.436) THEN
C...q + g -> q + QQ~[3P21]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=-COMFAC*PARU(1)*AS**3*(32D0/81D0)*
     &            ((6D0*SQMQQ**2+TH2)*UHSH2
     &            -2D0*SH*UH*(TH2+6D0*SQMQQ*UHSH))/
     &            (SQMQQR*TH*UHSH2**2)
          ELSE
            FF=-(32D0*PARU(1)*AS**3*SQMQQ*SQMQQR)/(27D0*TH2*UHSH2**2)
            C1=TH*UHSH2
            C2=4D0*(SH2+TH2+2D0*TH*UHSH)
            C3=4D0*UHSH2
            C4=8D0*SH*UHSH
            C5=8D0*TH
            C6=0D0
            C7=16D0*TH
            C8=0D0
            C9=-16D0*UHSH
            C0=16D0*SQMQQ
            IF(MSTP(147).EQ.0) THEN
               FACQQG=1D0/3D0*(C1*3D0
     &              -C2*(2D0*EL1K10*EL2K10+EL1K11*EL2K11)
     &              -C3*(2D0*EL1K20*EL2K20+EL1K21*EL2K21)
     &              -C4*(2D0*EL1K10*EL2K20+EL1K11*EL2K21)
     &              +C5*2D0*(EL1K10*EL2K10-EL1K11*EL2K11)**2
     &              +C6*2D0*(EL1K20*EL2K20-EL1K21*EL2K21)**2
     &              +C7*2D0*(EL1K10*EL2K10-EL1K11*EL2K11)
     &                     *(EL1K10*EL2K20-EL1K11*EL2K21)
     &              +C8*2D0*(EL1K20*EL2K20-EL1K21*EL2K21)
     &                     *(EL1K10*EL2K20-EL1K11*EL2K21)
     &              +C9*2D0*(EL1K10*EL2K10-EL1K11*EL2K11)
     &                     *(EL1K20*EL2K20-EL1K21*EL2K21)
     &              +C0*2D0*(EL1K10*EL2K20-EL1K11*EL2K21)**2)
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=C1*2D0
     &              -C2*(EL1K10*EL2K10+EL1K11*EL2K11)
     &              -C3*(EL1K20*EL2K20+EL1K21*EL2K21)
     &              -C4*(EL1K10*EL2K20+EL1K11*EL2K21)
     &              +C5*4D0*EL1K10*EL2K10*EL1K11*EL2K11
     &              +C6*4D0*EL1K20*EL2K20*EL1K21*EL2K21
     &              +C7*2D0*(EL1K10*EL2K10*EL1K11*EL2K21
     &                      +EL1K10*EL2K20*EL1K11*EL2K11)
     &              +C8*2D0*(EL1K20*EL2K20*EL1K11*EL2K21
     &                      +EL1K10*EL2K20*EL1K21*EL2K21)
     &              +C9*4D0*EL1K10*EL2K20*EL1K11*EL2K21
     &              +C0*(EL1K10*EL2K10*EL1K21*EL2K21
     &              +2D0*EL1K10*EL2K20*EL1K11*EL2K21
     &                  +EL1K20*EL2K20*EL1K11*EL2K11)
            ELSEIF(MSTP(147).EQ.2) THEN
               FACQQG=2D0*(C1
     &              -C2*EL1K11*EL2K11
     &              -C3*EL1K21*EL2K21
     &              -C4*EL1K11*EL2K21
     &              +C5*(EL1K11*EL2K11)**2
     &              +C6*(EL1K21*EL2K21)**2
     &              +C7*EL1K11*EL2K11*EL1K11*EL2K21
     &              +C8*EL1K21*EL2K21*EL1K11*EL2K21
     &              +(C9+C0)*(EL1K11*EL2K21)**2)
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
          DO 2456 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 2456
            DO 2455 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 2455
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 2455
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQQG*PARP(IONIUM+5)
 2455       CONTINUE
 2456     CONTINUE
 
        ELSEIF(ISUB.EQ.437) THEN
C...q + q~ -> g + QQ~[3P01]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=COMFAC*PARU(1)*AS**3*(128D0/243D0)*
     &            (SH-3D0*SQMQQ)**2*(TH2+UH2)/(SQMQQR*SH*THUH2**2)
          ELSE
            FA=PARU(1)*AS**3*(128D0/729D0)*
     &            (SH-3D0*SQMQQ)**2*(TH2+UH2)/(SQMQQR*SH*THUH2**2)
            IF(MSTP(147).EQ.0) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=COMFAC*2D0*FA
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=COMFAC*FA
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=0D0
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=0D0
            ENDIF
          ENDIF
          DO 2457 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 2457
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG*PARP(IONIUM+5)
 2457     CONTINUE
 
        ELSEIF(ISUB.EQ.438) THEN
C...q + q~ -> g + QQ~[3P11]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=COMFAC*PARU(1)*AS**3*256D0/81D0*
     &            (4D0*SQMQQ*TH*UH+SH*(TH2+UH2))/(SQMQQR*THUH2**2)
          ELSE
            FF=-(512D0*PARU(1)*AS**3*SQMQQR)/(81D0*THUH2**2)
            C1=TH*UH
            C2=2D0*UH
            C3=2D0*TH
            C4=2D0*THUH
            IF(MSTP(147).EQ.0) THEN
               FACQQG=-C1+C2*EL1K10*EL2K10+C3*EL1K20*EL2K20
     &              +C4*(EL1K10*EL2K20+EL1K20*EL2K10)/2D0
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=2D0*(-C1+C2*EL1K11*EL2K11+C3*EL1K21*EL2K21
     &              +C4*(EL1K11*EL2K21+EL1K21*EL2K11)/2D0)
            ELSEIF(MSTP(147).EQ.3) THEN
               FACQQG=-C1+C2*EL1K10*EL2K10+C3*EL1K20*EL2K20
     &              +C4*(EL1K10*EL2K20+EL1K20*EL2K10)/2D0
            ELSEIF(MSTP(147).EQ.4) THEN
               FACQQG=-C1+C2*EL1K11*EL2K11+C3*EL1K21*EL2K21
     &              +C4*(EL1K11*EL2K21+EL1K21*EL2K11)/2D0
            ELSEIF(MSTP(147).EQ.5) THEN
               FACQQG=C2*EL1K11*EL2K10+C3*EL1K21*EL2K20
     &              +C4*(EL1K11*EL2K20+EL1K21*EL2K10)/2D0
            ELSEIF(MSTP(147).EQ.6) THEN
               FACQQG=C2*EL1K11*EL2K11+C3*EL1K21*EL2K21
     &              +C4*(EL1K11*EL2K21+EL1K21*EL2K11)/2D0
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
          DO 2458 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 2458
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG*PARP(IONIUM+5)
 2458     CONTINUE
 
        ELSEIF(ISUB.EQ.439) THEN
C...q + q~ -> g + QQ~[3P21]
          IF(MSTP(145).EQ.0) THEN
            FACQQG=COMFAC*PARU(1)*AS**3*(256D0/243D0)*
     &            ((6D0*SQMQQ**2+SH2)*THUH2
     &            -2D0*TH*UH*(SH2+6D0*SQMQQ*THUH))/
     &            (SQMQQR*SH*THUH2**2)
          ELSE
            FF=(256D0*PARU(1)*AS**3*SQMQQ*SQMQQR)/(81D0*SH2*THUH2**2)
            C1=SH*THUH2
            C2=4D0*(SH2+UH2+2D0*SH*THUH)
            C3=4D0*(SH2+TH2+2D0*SH*THUH)
            C4=8D0*(SH2-TH*UH+2D0*SH*THUH)
            C5=8D0*SH
            C6=C5
            C7=16D0*SH
            C8=C7
            C9=-16D0*THUH
            C0=16D0*SQMQQ
            IF(MSTP(147).EQ.0) THEN
               FACQQG=1D0/3D0*(C1*3D0
     &              -C2*(2D0*EL1K10*EL2K10+EL1K11*EL2K11)
     &              -C3*(2D0*EL1K20*EL2K20+EL1K21*EL2K21)
     &              -C4*(2D0*EL1K10*EL2K20+EL1K11*EL2K21)
     &              +C5*2D0*(EL1K10*EL2K10-EL1K11*EL2K11)**2
     &              +C6*2D0*(EL1K20*EL2K20-EL1K21*EL2K21)**2
     &              +C7*2D0*(EL1K10*EL2K10-EL1K11*EL2K11)
     &                     *(EL1K10*EL2K20-EL1K11*EL2K21)
     &              +C8*2D0*(EL1K20*EL2K20-EL1K21*EL2K21)
     &                     *(EL1K10*EL2K20-EL1K11*EL2K21)
     &              +C9*2D0*(EL1K10*EL2K10-EL1K11*EL2K11)
     &                     *(EL1K20*EL2K20-EL1K21*EL2K21)
     &              +C0*2D0*(EL1K10*EL2K20-EL1K11*EL2K21)**2)
            ELSEIF(MSTP(147).EQ.1) THEN
               FACQQG=C1*2D0
     &              -C2*(EL1K10*EL2K10+EL1K11*EL2K11)
     &              -C3*(EL1K20*EL2K20+EL1K21*EL2K21)
     &              -C4*(EL1K10*EL2K20+EL1K11*EL2K21)
     &              +C5*4D0*EL1K10*EL2K10*EL1K11*EL2K11
     &              +C6*4D0*EL1K20*EL2K20*EL1K21*EL2K21
     &              +C7*2D0*(EL1K10*EL2K10*EL1K11*EL2K21
     &                      +EL1K10*EL2K20*EL1K11*EL2K11)
     &              +C8*2D0*(EL1K20*EL2K20*EL1K11*EL2K21
     &                      +EL1K10*EL2K20*EL1K21*EL2K21)
     &              +C9*4D0*EL1K10*EL2K20*EL1K11*EL2K21
     &              +C0*(EL1K10*EL2K10*EL1K21*EL2K21
     &              +2D0*EL1K10*EL2K20*EL1K11*EL2K21
     &                  +EL1K20*EL2K20*EL1K11*EL2K11)
            ELSEIF(MSTP(147).EQ.2) THEN
               FACQQG=2D0*(C1
     &              -C2*EL1K11*EL2K11
     &              -C3*EL1K21*EL2K21
     &              -C4*EL1K11*EL2K21
     &              +C5*(EL1K11*EL2K11)**2
     &              +C6*(EL1K21*EL2K21)**2
     &              +C7*EL1K11*EL2K11*EL1K11*EL2K21
     &              +C8*EL1K21*EL2K21*EL1K11*EL2K21
     &              +(C9+C0)*(EL1K11*EL2K21)**2)
            ENDIF
            FACQQG=COMFAC*FF*FACQQG
          ENDIF
          DO 2459 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 2459
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQG*PARP(IONIUM+5)
 2459     CONTINUE
        ENDIF
C...QUARKONIA---
 
      ENDIF
 
      RETURN
      END
