 
C*********************************************************************
 
C...PYSGQC
C...Subprocess cross sections for QCD processes,
C...including photons.
C...Auxiliary to PYSIGH.
 
      SUBROUTINE PYSGQC(NCHN,SIGS)
 
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
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      COMMON/PYSGCM/ISUB,ISUBSV,MMIN1,MMAX1,MMIN2,MMAX2,MMINA,MMAXA,
     &KFAC(2,-40:40),COMFAC,FACK,FACA,SH,TH,UH,SH2,TH2,UH2,SQM3,SQM4,
     &SHR,SQPTH,TAUP,BE34,CTH,X(2),SQMZ,SQMW,GMMZ,GMMW,
     &AEM,AS,XW,XW1,XWC,XWV,POLL,POLR,POLLL,POLRR
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYPARS/,/PYINT1/,/PYINT2/,
     &/PYINT3/,/PYINT4/,/PYINT7/,/PYSGCM/
C...Local arrays
      DIMENSION WDTP(0:400),WDTE(0:400,0:5)
 
C...Differential cross section expressions.
 
      IF(ISUB.LE.20) THEN
        IF(ISUB.EQ.10) THEN
C...f + f' -> f + f' (gamma/Z/W exchange)
          FACGGF=COMFAC*AEM**2*2D0*(SH2+UH2)/TH2
          FACGZF=COMFAC*AEM**2*XWC*4D0*SH2/(TH*(TH-SQMZ))
          FACZZF=COMFAC*(AEM*XWC)**2*2D0*SH2/(TH-SQMZ)**2
          FACWWF=COMFAC*(0.5D0*AEM/XW)**2*SH2/(TH-SQMW)**2
          DO 110 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 110
            IA=IABS(I)
            DO 100 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 100
              JA=IABS(J)
C...Electroweak couplings
              EI=KCHG(IA,1)*ISIGN(1,I)/3D0
              AI=SIGN(1D0,KCHG(IA,1)+0.5D0)*ISIGN(1,I)
              VI=AI-4D0*EI*XWV
              EJ=KCHG(JA,1)*ISIGN(1,J)/3D0
              AJ=SIGN(1D0,KCHG(JA,1)+0.5D0)*ISIGN(1,J)
              VJ=AJ-4D0*EJ*XWV
              EPSIJ=ISIGN(1,I*J)
C...gamma/Z exchange, only gamma exchange, or only Z exchange
              IF(MSTP(21).GE.1.AND.MSTP(21).LE.4) THEN
                IF(MSTP(21).EQ.1.OR.MSTP(21).EQ.4) THEN
                  FACNCF=FACGGF*EI**2*EJ**2+FACGZF*EI*EJ*
     &            (VI*VJ*(1D0+UH2/SH2)+AI*AJ*EPSIJ*(1D0-UH2/SH2))+
     &            FACZZF*((VI**2+AI**2)*(VJ**2+AJ**2)*(1D0+UH2/SH2)+
     &            4D0*VI*VJ*AI*AJ*EPSIJ*(1D0-UH2/SH2))
                ELSEIF(MSTP(21).EQ.2) THEN
                  FACNCF=FACGGF*EI**2*EJ**2
                ELSE
                  FACNCF=FACZZF*((VI**2+AI**2)*(VJ**2+AJ**2)*
     &            (1D0+UH2/SH2)+4D0*VI*VJ*AI*AJ*EPSIJ*(1D0-UH2/SH2))
                ENDIF
C...Extrafactor 2 for only one incoming neutrino spin state.
                IF(IA.GT.10.AND.MOD(IA,2).EQ.0) FACNCF=2D0*FACNCF
                IF(JA.GT.10.AND.MOD(JA,2).EQ.0) FACNCF=2D0*FACNCF
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=1
                SIGH(NCHN)=FACNCF
              ENDIF
C...W exchange
              IF((MSTP(21).EQ.1.OR.MSTP(21).EQ.5).AND.AI*AJ.LT.0D0) THEN
                FACCCF=FACWWF*VINT(180+I)*VINT(180+J)
                IF(EPSIJ.LT.0D0) FACCCF=FACCCF*UH2/SH2
                IF(IA.GT.10.AND.MOD(IA,2).EQ.0) FACCCF=2D0*FACCCF
                IF(JA.GT.10.AND.MOD(JA,2).EQ.0) FACCCF=2D0*FACCCF
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=2
                SIGH(NCHN)=FACCCF
              ENDIF
  100       CONTINUE
  110     CONTINUE
 
        ELSEIF(ISUB.EQ.11) THEN
C...f + f' -> f + f' (g exchange)
          FACQQ1=COMFAC*AS**2*4D0/9D0*(SH2+UH2)/TH2
          FACQQB=COMFAC*AS**2*4D0/9D0*((SH2+UH2)/TH2*FACA-
     &    MSTP(34)*2D0/3D0*UH2/(SH*TH))
          FACQQ2=COMFAC*AS**2*4D0/9D0*((SH2+TH2)/UH2-
     &    MSTP(34)*2D0/3D0*SH2/(TH*UH))
          DO 130 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.MSTP(58).OR.KFAC(1,I).EQ.0) GOTO 130
            DO 120 J=MMIN2,MMAX2
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.MSTP(58).OR.KFAC(2,J).EQ.0) GOTO 120
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQQ1
              IF(I.EQ.-J) SIGH(NCHN)=FACQQB
              IF(I.EQ.J) THEN
                SIGH(NCHN)=0.5D0*SIGH(NCHN)
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=2
                SIGH(NCHN)=0.5D0*FACQQ2
              ENDIF
  120       CONTINUE
  130     CONTINUE
 
        ELSEIF(ISUB.EQ.12) THEN
C...f + fbar -> f' + fbar' (q + qbar -> q' + qbar' only)
          CALL PYWIDT(21,SH,WDTP,WDTE)
          FACQQB=COMFAC*AS**2*4D0/9D0*(TH2+UH2)/SH2*
     &    (WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          DO 140 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 140
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQB
  140     CONTINUE
 
        ELSEIF(ISUB.EQ.13) THEN
C...f + fbar -> g + g (q + qbar -> g + g only)
          FACGG1=COMFAC*AS**2*32D0/27D0*(UH/TH-(2D0+MSTP(34)*1D0/4D0)*
     &    UH2/SH2)
          FACGG2=COMFAC*AS**2*32D0/27D0*(TH/UH-(2D0+MSTP(34)*1D0/4D0)*
     &    TH2/SH2)
          DO 150 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 150
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=0.5D0*FACGG1
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=2
            SIGH(NCHN)=0.5D0*FACGG2
  150     CONTINUE
 
        ELSEIF(ISUB.EQ.14) THEN
C...f + fbar -> g + gamma (q + qbar -> g + gamma only)
          FACGG=COMFAC*AS*AEM*8D0/9D0*(TH2+UH2)/(TH*UH)
          DO 160 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 160
            EI=KCHG(IABS(I),1)/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACGG*EI**2
  160     CONTINUE
 
        ELSEIF(ISUB.EQ.18) THEN
C...f + fbar -> gamma + gamma
          FACGG=COMFAC*AEM**2*2D0*(TH2+UH2)/(TH*UH)
          DO 170 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 170
            EI=KCHG(IABS(I),1)/3D0
            FCOI=1D0
            IF(IABS(I).LE.10) FCOI=FACA/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=0.5D0*FACGG*FCOI*EI**4
  170     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.40) THEN
        IF(ISUB.EQ.28) THEN
C...f + g -> f + g (q + g -> q + g only)
          FACQG1=COMFAC*AS**2*4D0/9D0*((2D0+MSTP(34)*1D0/4D0)*UH2/TH2-
     &    UH/SH)*FACA
          FACQG2=COMFAC*AS**2*4D0/9D0*((2D0+MSTP(34)*1D0/4D0)*SH2/TH2-
     &    SH/UH)
          DO 190 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.10) GOTO 190
            DO 180 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 180
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 180
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQG1
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=2
              SIGH(NCHN)=FACQG2
  180       CONTINUE
  190     CONTINUE
 
        ELSEIF(ISUB.EQ.29) THEN
C...f + g -> f + gamma (q + g -> q + gamma only)
          FGQ=COMFAC*FACA*AS*AEM*1D0/3D0*(SH2+UH2)/(-SH*UH)
          DO 210 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 210
            EI=KCHG(IABS(I),1)/3D0
            FACGQ=FGQ*EI**2
            DO 200 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 200
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 200
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACGQ
  200       CONTINUE
  210     CONTINUE
 
        ELSEIF(ISUB.EQ.33) THEN
C...f + gamma -> f + g (q + gamma -> q + g only)
          FGQ=COMFAC*AS*AEM*8D0/3D0*(SH2+UH2)/(-SH*UH)
          DO 230 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 230
            EI=KCHG(IABS(I),1)/3D0
            FACGQ=FGQ*EI**2
            DO 220 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 220
              IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 220
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=22
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACGQ
  220       CONTINUE
  230     CONTINUE
 
        ELSEIF(ISUB.EQ.34) THEN
C...f + gamma -> f + gamma
          FGQ=COMFAC*AEM**2*2D0*(SH2+UH2)/(-SH*UH)
          DO 250 I=MMINA,MMAXA
            IF(I.EQ.0) GOTO 250
            EI=KCHG(IABS(I),1)/3D0
            FACGQ=FGQ*EI**4
            DO 240 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 240
              IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 240
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=22
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACGQ
  240       CONTINUE
  250     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.80) THEN
        IF(ISUB.EQ.53) THEN
C...g + g -> f + fbar (g + g -> q + qbar only)
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 270
          IDC0=MDCY(21,2)-1
C...Begin by d, u, s flavours.
          FLAVWT=0D0
          IF(MDME(IDC0+1,1).GE.1) FLAVWT=FLAVWT+
     &    SQRT(MAX(0D0,1D0-4D0*PMAS(1,1)**2/SH))
          IF(MDME(IDC0+2,1).GE.1) FLAVWT=FLAVWT+
     &    SQRT(MAX(0D0,1D0-4D0*PMAS(2,1)**2/SH))
          IF(MDME(IDC0+3,1).GE.1) FLAVWT=FLAVWT+
     &    SQRT(MAX(0D0,1D0-4D0*PMAS(3,1)**2/SH))
          FACQQ1=COMFAC*AS**2*1D0/6D0*(UH/TH-(2D0+MSTP(34)*1D0/4D0)*
     &    UH2/SH2)*FLAVWT*FACA
          FACQQ2=COMFAC*AS**2*1D0/6D0*(TH/UH-(2D0+MSTP(34)*1D0/4D0)*
     &    TH2/SH2)*FLAVWT*FACA
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
C...Next c and b flavours: modified that and uhat for fixed
C...cos(theta-hat).
          DO 260 IFL=4,5
          SQMAVG=PMAS(IFL,1)**2
          IF(MDME(IDC0+IFL,1).GE.1.AND.SH.GT.4.04D0*SQMAVG) THEN
            BE34=SQRT(1D0-4D0*SQMAVG/SH)
            THQ=-0.5D0*SH*(1D0-BE34*CTH)
            UHQ=-0.5D0*SH*(1D0+BE34*CTH)
            THUHQ=THQ*UHQ-SQMAVG*SH
            IF(MSTP(34).EQ.0) THEN
              FACQQ1=UHQ/THQ-2D0*UHQ**2/SH2+4D0*(SQMAVG/SH)*THUHQ/THQ**2
              FACQQ2=THQ/UHQ-2D0*THQ**2/SH2+4D0*(SQMAVG/SH)*THUHQ/UHQ**2
            ELSE
              FACQQ1=UHQ/THQ-2.25D0*UHQ**2/SH2+4.5D0*(SQMAVG/SH)*THUHQ/
     &        THQ**2+0.5D0*SQMAVG*(THQ+SQMAVG)/THQ**2-SQMAVG**2/(SH*THQ)
              FACQQ2=THQ/UHQ-2.25D0*THQ**2/SH2+4.5D0*(SQMAVG/SH)*THUHQ/
     &        UHQ**2+0.5D0*SQMAVG*(UHQ+SQMAVG)/UHQ**2-SQMAVG**2/(SH*UHQ)
            ENDIF
            FACQQ1=COMFAC*FACA*AS**2*(1D0/6D0)*FACQQ1*BE34
            FACQQ2=COMFAC*FACA*AS**2*(1D0/6D0)*FACQQ2*BE34
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1+2*(IFL-3)
            SIGH(NCHN)=FACQQ1
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=2+2*(IFL-3)
            SIGH(NCHN)=FACQQ2
          ENDIF
  260     CONTINUE
  270     CONTINUE
 
        ELSEIF(ISUB.EQ.54) THEN
C...g + gamma -> f + fbar (g + gamma -> q + qbar only)
          CALL PYWIDT(21,SH,WDTP,WDTE)
          WDTESU=0D0
          DO 280 I=1,MIN(8,MDCY(21,3))
            EF=KCHG(I,1)/3D0
            WDTESU=WDTESU+EF**2*(WDTE(I,1)+WDTE(I,2)+WDTE(I,3)+
     &      WDTE(I,4))
  280     CONTINUE
          FACQQ=COMFAC*AEM*AS*WDTESU*(TH2+UH2)/(TH*UH)
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
 
        ELSEIF(ISUB.EQ.58) THEN
C...gamma + gamma -> f + fbar
          CALL PYWIDT(22,SH,WDTP,WDTE)
          WDTESU=0D0
          DO 290 I=1,MIN(12,MDCY(22,3))
            IF(I.LE.8) EF= KCHG(I,1)/3D0
            IF(I.GE.9) EF= KCHG(9+2*(I-8),1)/3D0
            WDTESU=WDTESU+EF**2*(WDTE(I,1)+WDTE(I,2)+WDTE(I,3)+
     &      WDTE(I,4))
  290     CONTINUE
          FACFF=COMFAC*AEM**2*WDTESU*2D0*(TH2+UH2)/(TH*UH)
          IF(KFAC(1,22)*KFAC(2,22).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=22
            ISIG(NCHN,2)=22
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACFF
          ENDIF
 
        ELSEIF(ISUB.EQ.68) THEN
C...g + g -> g + g
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 300
          FACGG1=COMFAC*AS**2*9D0/4D0*(SH2/TH2+2D0*SH/TH+3D0+2D0*TH/SH+
     &    TH2/SH2)*FACA
          FACGG2=COMFAC*AS**2*9D0/4D0*(UH2/SH2+2D0*UH/SH+3D0+2D0*SH/UH+
     &    SH2/UH2)*FACA
          FACGG3=COMFAC*AS**2*9D0/4D0*(TH2/UH2+2D0*TH/UH+3D0+2D0*UH/TH+
     &    UH2/TH2)
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=0.5D0*FACGG1
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=2
          SIGH(NCHN)=0.5D0*FACGG2
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=3
          SIGH(NCHN)=0.5D0*FACGG3
  300     CONTINUE
 
        ELSEIF(ISUB.EQ.80) THEN
C...q + gamma -> q' + pi+/-
          FQPI=COMFAC*(2D0*AEM/9D0)*(-SH/TH)*(1D0/SH2+1D0/TH2)
          ASSH=PYALPS(MAX(0.5D0,0.5D0*SH))
          Q2FPSH=0.55D0/LOG(MAX(2D0,2D0*SH))
          DELSH=UH*SQRT(ASSH*Q2FPSH)
          ASUH=PYALPS(MAX(0.5D0,-0.5D0*UH))
          Q2FPUH=0.55D0/LOG(MAX(2D0,-2D0*UH))
          DELUH=SH*SQRT(ASUH*Q2FPUH)
          DO 320 I=MAX(-2,MMINA),MIN(2,MMAXA)
            IF(I.EQ.0) GOTO 320
            EI=KCHG(IABS(I),1)/3D0
            EJ=SIGN(1D0-ABS(EI),EI)
            DO 310 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 310
              IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 310
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=22
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FQPI*(EI*DELSH+EJ*DELUH)**2
  310       CONTINUE
  320     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.100) THEN
        IF(ISUB.EQ.91) THEN
C...Elastic scattering
          SIGS=VINT(315)*VINT(316)*SIGT(0,0,1)
 
        ELSEIF(ISUB.EQ.92) THEN
C...Single diffractive scattering (first side, i.e. XB)
          SIGS=VINT(315)*VINT(316)*SIGT(0,0,2)
 
        ELSEIF(ISUB.EQ.93) THEN
C...Single diffractive scattering (second side, i.e. AX)
          SIGS=VINT(315)*VINT(316)*SIGT(0,0,3)
 
        ELSEIF(ISUB.EQ.94) THEN
C...Double diffractive scattering
          SIGS=VINT(315)*VINT(316)*SIGT(0,0,4)
 
        ELSEIF(ISUB.EQ.95) THEN
C...Low-pT scattering
          SIGS=VINT(315)*VINT(316)*SIGT(0,0,5)
 
        ELSEIF(ISUB.EQ.96) THEN
C...Multiple interactions: sum of QCD processes
          CALL PYWIDT(21,SH,WDTP,WDTE)
 
C...q + q' -> q + q'
          FACQQ1=COMFAC*AS**2*4D0/9D0*(SH2+UH2)/TH2
          FACQQB=COMFAC*AS**2*4D0/9D0*((SH2+UH2)/TH2*FACA-
     &    MSTP(34)*2D0/3D0*UH2/(SH*TH))
          FACQQ2=COMFAC*AS**2*4D0/9D0*(SH2+TH2)/UH2
          FACQQI=-COMFAC*AS**2*4D0/9D0*MSTP(34)*2D0/3D0*SH2/(TH*UH)
          RATQQI=(FACQQ1+FACQQ2+FACQQI)/(FACQQ1+FACQQ2)
          DO 340 I=-5,5
            IF(I.EQ.0) GOTO 340
            DO 330 J=-5,5
              IF(J.EQ.0) GOTO 330
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=111
              SIGH(NCHN)=FACQQ1
              IF(I.EQ.-J) SIGH(NCHN)=FACQQB
              IF(I.EQ.J) THEN
                SIGH(NCHN)=0.5D0*FACQQ1*RATQQI
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=112
                SIGH(NCHN)=0.5D0*FACQQ2*RATQQI
              ENDIF
  330       CONTINUE
  340     CONTINUE
 
C...q + qbar -> q' + qbar' or g + g
          FACQQB=COMFAC*AS**2*4D0/9D0*(TH2+UH2)/SH2*
     &    (WDTE(0,1)+WDTE(0,2)+WDTE(0,3)+WDTE(0,4))
          FACGG1=COMFAC*AS**2*32D0/27D0*(UH/TH-(2D0+MSTP(34)*1D0/4D0)*
     &    UH2/SH2)
          FACGG2=COMFAC*AS**2*32D0/27D0*(TH/UH-(2D0+MSTP(34)*1D0/4D0)*
     &    TH2/SH2)
          DO 350 I=-5,5
            IF(I.EQ.0) GOTO 350
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=121
            SIGH(NCHN)=FACQQB
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=131
            SIGH(NCHN)=0.5D0*FACGG1
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=132
            SIGH(NCHN)=0.5D0*FACGG2
  350     CONTINUE
 
C...q + g -> q + g
          FACQG1=COMFAC*AS**2*4D0/9D0*((2D0+MSTP(34)*1D0/4D0)*UH2/TH2-
     &    UH/SH)*FACA
          FACQG2=COMFAC*AS**2*4D0/9D0*((2D0+MSTP(34)*1D0/4D0)*SH2/TH2-
     &    SH/UH)
          DO 370 I=-5,5
            IF(I.EQ.0) GOTO 370
            DO 360 ISDE=1,2
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=281
              SIGH(NCHN)=FACQG1
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=282
              SIGH(NCHN)=FACQG2
  360       CONTINUE
  370     CONTINUE
 
C...g + g -> q + qbar (only d, u, s)
          IDC0=MDCY(21,2)-1
          FLAVWT=0D0
          IF(MDME(IDC0+1,1).GE.1) FLAVWT=FLAVWT+
     &    SQRT(MAX(0D0,1D0-4D0*PMAS(1,1)**2/SH))
          IF(MDME(IDC0+2,1).GE.1) FLAVWT=FLAVWT+
     &    SQRT(MAX(0D0,1D0-4D0*PMAS(2,1)**2/SH))
          IF(MDME(IDC0+3,1).GE.1) FLAVWT=FLAVWT+
     &    SQRT(MAX(0D0,1D0-4D0*PMAS(3,1)**2/SH))
          FACQQ1=COMFAC*AS**2*1D0/6D0*(UH/TH-(2D0+MSTP(34)*1D0/4D0)*
     &    UH2/SH2)*FLAVWT*FACA
          FACQQ2=COMFAC*AS**2*1D0/6D0*(TH/UH-(2D0+MSTP(34)*1D0/4D0)*
     &    TH2/SH2)*FLAVWT*FACA
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=531
          SIGH(NCHN)=FACQQ1
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=532
          SIGH(NCHN)=FACQQ2
 
C...g + g -> c + cbar, b + bbar: modified that/uhat for fixed
C...cos(theta-hat)
          DO 380 IFL=4,5
          SQMAVG=PMAS(IFL,1)**2
          IF(MDME(IDC0+IFL,1).GE.1.AND.SH.GT.4.04D0*SQMAVG) THEN
            BE34=SQRT(1D0-4D0*SQMAVG/SH)
            THQ=-0.5D0*SH*(1D0-BE34*CTH)
            UHQ=-0.5D0*SH*(1D0+BE34*CTH)
            THUHQ=THQ*UHQ-SQMAVG*SH
            IF(MSTP(34).EQ.0) THEN
              FACQQ1=UHQ/THQ-2D0*UHQ**2/SH2+4D0*(SQMAVG/SH)*THUHQ/THQ**2
              FACQQ2=THQ/UHQ-2D0*THQ**2/SH2+4D0*(SQMAVG/SH)*THUHQ/UHQ**2
            ELSE
              FACQQ1=UHQ/THQ-2.25D0*UHQ**2/SH2+4.5D0*(SQMAVG/SH)*THUHQ/
     &        THQ**2+0.5D0*SQMAVG*(THQ+SQMAVG)/THQ**2-SQMAVG**2/(SH*THQ)
              FACQQ2=THQ/UHQ-2.25D0*THQ**2/SH2+4.5D0*(SQMAVG/SH)*THUHQ/
     &        UHQ**2+0.5D0*SQMAVG*(UHQ+SQMAVG)/UHQ**2-SQMAVG**2/(SH*UHQ)
            ENDIF
            FACQQ1=COMFAC*FACA*AS**2*(1D0/6D0)*FACQQ1*BE34
            FACQQ2=COMFAC*FACA*AS**2*(1D0/6D0)*FACQQ2*BE34
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=531+2*(IFL-3)
            SIGH(NCHN)=FACQQ1
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=532+2*(IFL-3)
            SIGH(NCHN)=FACQQ2
          ENDIF
  380     CONTINUE
 
C...g + g -> g + g
          FACGG1=COMFAC*AS**2*9D0/4D0*(SH2/TH2+2D0*SH/TH+3D0+
     &    2D0*TH/SH+TH2/SH2)*FACA
          FACGG2=COMFAC*AS**2*9D0/4D0*(UH2/SH2+2D0*UH/SH+3D0+
     &    2D0*SH/UH+SH2/UH2)*FACA
          FACGG3=COMFAC*AS**2*9D0/4D0*(TH2/UH2+2D0*TH/UH+3+
     &    2D0*UH/TH+UH2/TH2)
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=681
          SIGH(NCHN)=0.5D0*FACGG1
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=682
          SIGH(NCHN)=0.5D0*FACGG2
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=683
          SIGH(NCHN)=0.5D0*FACGG3
 
        ELSEIF(ISUB.EQ.99) THEN
C...f + gamma* -> f.
          IF(MINT(107).EQ.4) THEN
            Q2GA=VINT(307)
            P2GA=VINT(308)
            ISDE=2
          ELSE
            Q2GA=VINT(308)
            P2GA=VINT(307)
            ISDE=1
          ENDIF
          COMFAC=PARU(5)*4D0*PARU(1)**2*PARU(101)*VINT(315)*VINT(316)
          PM2RHO=PMAS(PYCOMP(113),1)**2
          IF(MSTP(19).EQ.0) THEN
            COMFAC=COMFAC/Q2GA
          ELSEIF(MSTP(19).EQ.1) THEN
            COMFAC=COMFAC/(Q2GA+PM2RHO)
          ELSEIF(MSTP(19).EQ.2) THEN
            COMFAC=COMFAC*Q2GA/(Q2GA+PM2RHO)**2
          ELSE
            COMFAC=COMFAC*Q2GA/(Q2GA+PM2RHO)**2
            W2GA=VINT(2)
            IF(MINT(11).EQ.22.AND.MINT(12).EQ.22) THEN
              RDRDS=4.1D-3*W2GA**2.167D0/((Q2GA+0.15D0*W2GA)**2*
     &        Q2GA**0.75D0)*(1D0+0.11D0*Q2GA*P2GA/(1D0+0.02D0*P2GA**2))
              XGA=Q2GA/(W2GA+VINT(307)+VINT(308))
            ELSE
              RDRDS=1.5D-4*W2GA**2.167D0/((Q2GA+0.041D0*W2GA)**2*
     &        Q2GA**0.57D0)
              XGA=Q2GA/(W2GA+Q2GA-PMAS(PYCOMP(MINT(10+ISDE)),1)**2)
            ENDIF
            COMFAC=COMFAC*EXP(-MAX(1D-10,RDRDS))
            IF(MSTP(19).EQ.4) COMFAC=COMFAC/MAX(1D-2,1D0-XGA)
          ENDIF
          DO 390 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(ISDE,I).EQ.0) GOTO 390
            IF(IABS(I).LT.10.AND.IABS(I).GT.MSTP(58)) GOTO 390
            EI=KCHG(IABS(I),1)/3D0
            NCHN=NCHN+1
            ISIG(NCHN,ISDE)=I
            ISIG(NCHN,3-ISDE)=22
            ISIG(NCHN,3)=1
            SIGH(NCHN)=COMFAC*EI**2
  390     CONTINUE
        ENDIF
 
      ELSE
        IF(ISUB.EQ.114.OR.ISUB.EQ.115) THEN
C...g + g -> gamma + gamma or g + g -> g + gamma
          A0STUR=0D0
          A0STUI=0D0
          A0TSUR=0D0
          A0TSUI=0D0
          A0UTSR=0D0
          A0UTSI=0D0
          A1STUR=0D0
          A1STUI=0D0
          A2STUR=0D0
          A2STUI=0D0
          ALST=LOG(-SH/TH)
          ALSU=LOG(-SH/UH)
          ALTU=LOG(TH/UH)
          IMAX=2*MSTP(1)
          IF(MSTP(38).GE.1.AND.MSTP(38).LE.8) IMAX=MSTP(38)
          DO 400 I=1,IMAX
            EI=KCHG(IABS(I),1)/3D0
            EIWT=EI**2
            IF(ISUB.EQ.115) EIWT=EI
            SQMQ=PMAS(I,1)**2
            EPSS=4D0*SQMQ/SH
            EPST=4D0*SQMQ/TH
            EPSU=4D0*SQMQ/UH
            IF((MSTP(38).GE.1.AND.MSTP(38).LE.8).OR.EPSS.LT.1D-4) THEN
              B0STUR=1D0+(TH-UH)/SH*ALTU+0.5D0*(TH2+UH2)/SH2*(ALTU**2+
     &        PARU(1)**2)
              B0STUI=0D0
              B0TSUR=1D0+(SH-UH)/TH*ALSU+0.5D0*(SH2+UH2)/TH2*ALSU**2
              B0TSUI=-PARU(1)*((SH-UH)/TH+(SH2+UH2)/TH2*ALSU)
              B0UTSR=1D0+(SH-TH)/UH*ALST+0.5D0*(SH2+TH2)/UH2*ALST**2
              B0UTSI=-PARU(1)*((SH-TH)/UH+(SH2+TH2)/UH2*ALST)
              B1STUR=-1D0
              B1STUI=0D0
              B2STUR=-1D0
              B2STUI=0D0
            ELSE
              CALL PYWAUX(1,EPSS,W1SR,W1SI)
              CALL PYWAUX(1,EPST,W1TR,W1TI)
              CALL PYWAUX(1,EPSU,W1UR,W1UI)
              CALL PYWAUX(2,EPSS,W2SR,W2SI)
              CALL PYWAUX(2,EPST,W2TR,W2TI)
              CALL PYWAUX(2,EPSU,W2UR,W2UI)
              CALL PYI3AU(EPSS,TH/UH,Y3STUR,Y3STUI)
              CALL PYI3AU(EPSS,UH/TH,Y3SUTR,Y3SUTI)
              CALL PYI3AU(EPST,SH/UH,Y3TSUR,Y3TSUI)
              CALL PYI3AU(EPST,UH/SH,Y3TUSR,Y3TUSI)
              CALL PYI3AU(EPSU,SH/TH,Y3USTR,Y3USTI)
              CALL PYI3AU(EPSU,TH/SH,Y3UTSR,Y3UTSI)
              B0STUR=1D0+(1D0+2D0*TH/SH)*W1TR+(1D0+2D0*UH/SH)*W1UR+
     &        0.5D0*((TH2+UH2)/SH2-EPSS)*(W2TR+W2UR)-
     &        0.25D0*EPST*(1D0-0.5D0*EPSS)*(Y3SUTR+Y3TUSR)-
     &        0.25D0*EPSU*(1D0-0.5D0*EPSS)*(Y3STUR+Y3UTSR)+
     &        0.25D0*(-2D0*(TH2+UH2)/SH2+4D0*EPSS+EPST+EPSU+
     &        0.5D0*EPST*EPSU)*(Y3TSUR+Y3USTR)
              B0STUI=(1D0+2D0*TH/SH)*W1TI+(1D0+2D0*UH/SH)*W1UI+
     &        0.5D0*((TH2+UH2)/SH2-EPSS)*(W2TI+W2UI)-
     &        0.25D0*EPST*(1D0-0.5D0*EPSS)*(Y3SUTI+Y3TUSI)-
     &        0.25D0*EPSU*(1D0-0.5D0*EPSS)*(Y3STUI+Y3UTSI)+
     &        0.25D0*(-2D0*(TH2+UH2)/SH2+4D0*EPSS+EPST+EPSU+
     &        0.5D0*EPST*EPSU)*(Y3TSUI+Y3USTI)
              B0TSUR=1D0+(1D0+2D0*SH/TH)*W1SR+(1D0+2D0*UH/TH)*W1UR+
     &        0.5D0*((SH2+UH2)/TH2-EPST)*(W2SR+W2UR)-
     &        0.25D0*EPSS*(1D0-0.5D0*EPST)*(Y3TUSR+Y3SUTR)-
     &        0.25D0*EPSU*(1D0-0.5D0*EPST)*(Y3TSUR+Y3USTR)+
     &        0.25D0*(-2D0*(SH2+UH2)/TH2+4D0*EPST+EPSS+EPSU+
     &        0.5D0*EPSS*EPSU)*(Y3STUR+Y3UTSR)
              B0TSUI=(1D0+2D0*SH/TH)*W1SI+(1D0+2D0*UH/TH)*W1UI+
     &        0.5D0*((SH2+UH2)/TH2-EPST)*(W2SI+W2UI)-
     &        0.25D0*EPSS*(1D0-0.5D0*EPST)*(Y3TUSI+Y3SUTI)-
     &        0.25D0*EPSU*(1D0-0.5D0*EPST)*(Y3TSUI+Y3USTI)+
     &        0.25D0*(-2D0*(SH2+UH2)/TH2+4D0*EPST+EPSS+EPSU+
     &        0.5D0*EPSS*EPSU)*(Y3STUI+Y3UTSI)
              B0UTSR=1D0+(1D0+2D0*TH/UH)*W1TR+(1D0+2D0*SH/UH)*W1SR+
     &        0.5D0*((TH2+SH2)/UH2-EPSU)*(W2TR+W2SR)-
     &        0.25D0*EPST*(1D0-0.5D0*EPSU)*(Y3USTR+Y3TSUR)-
     &        0.25D0*EPSS*(1D0-0.5D0*EPSU)*(Y3UTSR+Y3STUR)+
     &        0.25D0*(-2D0*(TH2+SH2)/UH2+4D0*EPSU+EPST+EPSS+
     &        0.5D0*EPST*EPSS)*(Y3TUSR+Y3SUTR)
              B0UTSI=(1D0+2D0*TH/UH)*W1TI+(1D0+2D0*SH/UH)*W1SI+
     &        0.5D0*((TH2+SH2)/UH2-EPSU)*(W2TI+W2SI)-
     &        0.25D0*EPST*(1D0-0.5D0*EPSU)*(Y3USTI+Y3TSUI)-
     &        0.25D0*EPSS*(1D0-0.5D0*EPSU)*(Y3UTSI+Y3STUI)+
     &        0.25D0*(-2D0*(TH2+SH2)/UH2+4D0*EPSU+EPST+EPSS+
     &        0.5D0*EPST*EPSS)*(Y3TUSI+Y3SUTI)
              B1STUR=-1D0-0.25D0*(EPSS+EPST+EPSU)*(W2SR+W2TR+W2UR)+
     &        0.25D0*(EPSU+0.5D0*EPSS*EPST)*(Y3SUTR+Y3TUSR)+
     &        0.25D0*(EPST+0.5D0*EPSS*EPSU)*(Y3STUR+Y3UTSR)+
     &        0.25D0*(EPSS+0.5D0*EPST*EPSU)*(Y3TSUR+Y3USTR)
              B1STUI=-0.25D0*(EPSS+EPST+EPSU)*(W2SI+W2TI+W2UI)+
     &        0.25D0*(EPSU+0.5D0*EPSS*EPST)*(Y3SUTI+Y3TUSI)+
     &        0.25D0*(EPST+0.5D0*EPSS*EPSU)*(Y3STUI+Y3UTSI)+
     &        0.25D0*(EPSS+0.5D0*EPST*EPSU)*(Y3TSUI+Y3USTI)
              B2STUR=-1D0+0.125D0*EPSS*EPST*(Y3SUTR+Y3TUSR)+
     &        0.125D0*EPSS*EPSU*(Y3STUR+Y3UTSR)+
     &        0.125D0*EPST*EPSU*(Y3TSUR+Y3USTR)
              B2STUI=0.125D0*EPSS*EPST*(Y3SUTI+Y3TUSI)+
     &        0.125D0*EPSS*EPSU*(Y3STUI+Y3UTSI)+
     &        0.125D0*EPST*EPSU*(Y3TSUI+Y3USTI)
            ENDIF
            A0STUR=A0STUR+EIWT*B0STUR
            A0STUI=A0STUI+EIWT*B0STUI
            A0TSUR=A0TSUR+EIWT*B0TSUR
            A0TSUI=A0TSUI+EIWT*B0TSUI
            A0UTSR=A0UTSR+EIWT*B0UTSR
            A0UTSI=A0UTSI+EIWT*B0UTSI
            A1STUR=A1STUR+EIWT*B1STUR
            A1STUI=A1STUI+EIWT*B1STUI
            A2STUR=A2STUR+EIWT*B2STUR
            A2STUI=A2STUI+EIWT*B2STUI
  400     CONTINUE
          ASQSUM=A0STUR**2+A0STUI**2+A0TSUR**2+A0TSUI**2+A0UTSR**2+
     &    A0UTSI**2+4D0*A1STUR**2+4D0*A1STUI**2+A2STUR**2+A2STUI**2
          FACGG=COMFAC*FACA/(16D0*PARU(1)**2)*AS**2*AEM**2*ASQSUM
          FACGP=COMFAC*FACA*5D0/(192D0*PARU(1)**2)*AS**3*AEM*ASQSUM
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 410
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          IF(ISUB.EQ.114) SIGH(NCHN)=0.5D0*FACGG
          IF(ISUB.EQ.115) SIGH(NCHN)=FACGP
  410     CONTINUE
 
        ELSEIF(ISUB.EQ.131.OR.ISUB.EQ.132) THEN
C...f + gamma*_(T,L) -> f + g (q + gamma*_(T,L) -> q + g only)
          PH=0D0
          IF(MINT(15).EQ.22.AND.MINT(107).EQ.0.AND.VINT(3).LT.0D0)
     &    PH=VINT(3)**2
          IF(MINT(16).EQ.22.AND.MINT(108).EQ.0.AND.VINT(4).LT.0D0)
     &    PH=VINT(4)**2
          IF(ISUB.EQ.131) THEN
            FGQ=COMFAC*AS*AEM*8D0/3D0*SH**2/(SH+PH)**2*
     &      ((SH2+UH2-2D0*PH*TH)/(-SH*UH)-2D0*PH*TH/(SH+PH)**2)
          ELSE
            FGQ=COMFAC*AS*AEM*8D0/3D0*SH**2/(SH+PH)**4*(-4D0*PH*TH)
          ENDIF
          DO 430 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 430
            EI=KCHG(IABS(I),1)/3D0
            FACGQ=FGQ*EI**2
            DO 420 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 420
              IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 420
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=22
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACGQ
  420       CONTINUE
  430     CONTINUE
 
        ELSEIF(ISUB.EQ.133.OR.ISUB.EQ.134) THEN
C...f + gamma*_(T,L) -> f + gamma
          PH=0D0
          IF(MINT(15).EQ.22.AND.MINT(107).EQ.0.AND.VINT(3).LT.0D0)
     &    PH=VINT(3)**2
          IF(MINT(16).EQ.22.AND.MINT(108).EQ.0.AND.VINT(4).LT.0D0)
     &    PH=VINT(4)**2
          IF(ISUB.EQ.133) THEN
            FGQ=COMFAC*AEM**2*2D0*SH**2/(SH+PH)**2*
     &      ((SH2+UH2-2D0*PH*TH)/(-SH*UH)-2D0*PH*TH/(SH+PH)**2)
          ELSE
            FGQ=COMFAC*AEM**2*2D0*SH**2/(SH+PH)**4*(-4D0*PH*TH)
          ENDIF
          DO 450 I=MMINA,MMAXA
            IF(I.EQ.0) GOTO 450
            EI=KCHG(IABS(I),1)/3D0
            FACGQ=FGQ*EI**4
            DO 440 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 440
              IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 440
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=22
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACGQ
  440       CONTINUE
  450     CONTINUE
 
        ELSEIF(ISUB.EQ.135.OR.ISUB.EQ.136) THEN
C...g + gamma*_(T,L) -> f + fbar (g + gamma*_(T,L) -> q + qbar only)
          PH=0D0
          IF(MINT(15).EQ.22.AND.MINT(107).EQ.0.AND.VINT(3).LT.0D0)
     &    PH=VINT(3)**2
          IF(MINT(16).EQ.22.AND.MINT(108).EQ.0.AND.VINT(4).LT.0D0)
     &    PH=VINT(4)**2
          CALL PYWIDT(21,SH,WDTP,WDTE)
          WDTESU=0D0
          DO 460 I=1,MIN(8,MDCY(21,3))
            EF=KCHG(I,1)/3D0
            WDTESU=WDTESU+EF**2*(WDTE(I,1)+WDTE(I,2)+WDTE(I,3)+
     &      WDTE(I,4))
  460     CONTINUE
          IF(ISUB.EQ.135) THEN
            FACQQ=COMFAC*AEM*AS*WDTESU*SH**2/(SH+PH)**2*
     &      ((TH2+UH2-2D0*PH*SH)/(TH*UH)+4D0*PH*SH/(SH+PH)**2)
          ELSE
            FACQQ=COMFAC*AEM*AS*WDTESU*SH**2/(SH+PH)**4*8D0*PH*SH
          ENDIF
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
 
        ELSEIF(ISUB.GE.137.AND.ISUB.LE.140) THEN
C...gamma*_(T,L) + gamma*_(T,L) -> f + fbar
          PH1=0D0
          IF(VINT(3).LT.0D0) PH1=VINT(3)**2
          PH2=0D0
          IF(VINT(4).LT.0D0) PH2=VINT(4)**2
          CALL PYWIDT(22,SH,WDTP,WDTE)
          WDTESU=0D0
          DO 470 I=1,MIN(12,MDCY(22,3))
            IF(I.LE.8) EF= KCHG(I,1)/3D0
            IF(I.GE.9) EF= KCHG(9+2*(I-8),1)/3D0
            WDTESU=WDTESU+EF**2*(WDTE(I,1)+WDTE(I,2)+WDTE(I,3)+
     &      WDTE(I,4))
  470     CONTINUE
          DLAMB2=(TH+UH)**2-4D0*PH1*PH2
          IF(ISUB.EQ.137) THEN
            FPARAM=-SH*(TH+UH)/DLAMB2
            FACFF=COMFAC*AEM**2*WDTESU*2D0*SH2/(DLAMB2*TH2*UH2)*
     &      (TH*UH-PH1*PH2)*((TH2+UH2)*(1D0-2D0*FPARAM*(1D0-FPARAM))-
     &      2D0*PH1*PH2*FPARAM**2)
          ELSEIF(ISUB.EQ.138) THEN
            FACFF=COMFAC*AEM**2*WDTESU*4D0*SH2*SH/(DLAMB2**2*TH2*UH2)*
     &      PH2*(4D0*(TH*UH-PH1*PH2)*(TH*UH+PH1*SH*(TH-UH)**2/DLAMB2)+
     &      2D0*PH1**2*(TH-UH)**2)
          ELSEIF(ISUB.EQ.139) THEN
            FACFF=COMFAC*AEM**2*WDTESU*4D0*SH2*SH/(DLAMB2**2*TH2*UH2)*
     &      PH1*(4D0*(TH*UH-PH1*PH2)*(TH*UH+PH2*SH*(TH-UH)**2/DLAMB2)+
     &      2D0*PH2**2*(TH-UH)**2)
          ELSE
            FACFF=COMFAC*AEM**2*WDTESU*32D0*SH2**2/(DLAMB2**3*TH2*UH2)*
     &      PH1*PH2*(TH*UH-PH1*PH2)*(TH-UH)**2
          ENDIF
          IF(KFAC(1,22)*KFAC(2,22).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=22
            ISIG(NCHN,2)=22
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACFF
          ENDIF
 
        ENDIF
      ENDIF
 
      RETURN
      END
