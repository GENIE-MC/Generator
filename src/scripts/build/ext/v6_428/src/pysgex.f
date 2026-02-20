 
C*********************************************************************
 
C...PYSGEX
C...Subprocess cross sections for assorted exotic processes,
C...including Z'/W'/LQ/R/f*/H++/Z_R/W_R/G*.
C...Auxiliary to PYSIGH.
 
      SUBROUTINE PYSGEX(NCHN,SIGS)
 
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
      COMMON/PYTCSM/ITCM(0:99),RTCM(0:99)
      COMMON/PYSGCM/ISUB,ISUBSV,MMIN1,MMAX1,MMIN2,MMAX2,MMINA,MMAXA,
     &KFAC(2,-40:40),COMFAC,FACK,FACA,SH,TH,UH,SH2,TH2,UH2,SQM3,SQM4,
     &SHR,SQPTH,TAUP,BE34,CTH,X(2),SQMZ,SQMW,GMMZ,GMMW,
     &AEM,AS,XW,XW1,XWC,XWV,POLL,POLR,POLLL,POLRR
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYPARS/,/PYINT1/,/PYINT2/,
     &/PYINT3/,/PYINT4/,/PYTCSM/,/PYSGCM/
C...Local arrays
      DIMENSION WDTP(0:400),WDTE(0:400,0:5)
 
C...Differential cross section expressions.
 
      IF(ISUB.LE.160) THEN
        IF(ISUB.EQ.141) THEN
C...f + fbar -> gamma*/Z0/Z'0
          SQMZP=PMAS(32,1)**2
          MINT(61)=2
          CALL PYWIDT(32,SH,WDTP,WDTE)
          HP0=AEM/3D0*SH
          HP1=AEM/3D0*XWC*SH
          HP2=HP1
          HS=SHR*VINT(117)
          HSP=SHR*WDTP(0)
          FACZP=4D0*COMFAC*3D0
          DO 100 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 100
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            IA=IABS(I)
            IF(IA.LT.10) THEN
              IF(IA.LE.2) THEN
                VPI=PARU(123-2*MOD(IABS(I),2))
                API=PARU(124-2*MOD(IABS(I),2))
              ELSEIF(IA.LE.4) THEN
                VPI=PARJ(182-2*MOD(IABS(I),2))
                API=PARJ(183-2*MOD(IABS(I),2))
              ELSE
                VPI=PARJ(190-2*MOD(IABS(I),2))
                API=PARJ(191-2*MOD(IABS(I),2))
              ENDIF
            ELSE
              IF(IA.LE.12) THEN
                VPI=PARU(127-2*MOD(IABS(I),2))
                API=PARU(128-2*MOD(IABS(I),2))
              ELSEIF(IA.LE.14) THEN
                VPI=PARJ(186-2*MOD(IABS(I),2))
                API=PARJ(187-2*MOD(IABS(I),2))
              ELSE
                VPI=PARJ(194-2*MOD(IABS(I),2))
                API=PARJ(195-2*MOD(IABS(I),2))
              ENDIF
            ENDIF
            HI0=HP0
            IF(IABS(I).LE.10) HI0=HI0*FACA/3D0
            HI1=HP1
            IF(IABS(I).LE.10) HI1=HI1*FACA/3D0
            HI2=HP2
            IF(IABS(I).LE.10) HI2=HI2*FACA/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
C...Special case: if only branching ratios known then use them.
            IF(MWID(32).EQ.2.AND.MSTP(44).EQ.3) THEN
              HI=0D0
              IF(IA.LT.10) THEN
                HI=SHR*WDTP(IA)*FACA/9D0
              ELSEIF(IA.LT.20) THEN
                HI=SHR*WDTP(IA-2)
              ENDIF
              HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
              SIGH(NCHN)=HI*FACZP*HF/((SH-SQMZP)**2+HSP**2)
            ELSE
C...Normal cross section.
              SIGH(NCHN)=FACZP*(EI**2/SH2*HI0*HP0*VINT(111)+EI*VI*
     &        (1D0-SQMZ/SH)/((SH-SQMZ)**2+HS**2)*(HI0*HP1+HI1*HP0)*
     &        VINT(112)+EI*VPI*(1D0-SQMZP/SH)/((SH-SQMZP)**2+HSP**2)*
     &        (HI0*HP2+HI2*HP0)*VINT(113)+(VI**2+AI**2)/
     &        ((SH-SQMZ)**2+HS**2)*HI1*HP1*VINT(114)+(VI*VPI+AI*API)*
     &        ((SH-SQMZ)*(SH-SQMZP)+HS*HSP)/(((SH-SQMZ)**2+HS**2)*
     &        ((SH-SQMZP)**2+HSP**2))*(HI1*HP2+HI2*HP1)*VINT(115)+
     &        (VPI**2+API**2)/((SH-SQMZP)**2+HSP**2)*HI2*HP2*VINT(116))
            ENDIF
  100     CONTINUE
 
        ELSEIF(ISUB.EQ.142) THEN
C...f + fbar' -> W'+/-
          SQMWP=PMAS(34,1)**2
          CALL PYWIDT(34,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=4D0*COMFAC/((SH-SQMWP)**2+HS**2)*3D0
          HP=AEM/(24D0*XW)*SH
          DO 120 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 120
            IA=IABS(I)
            DO 110 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 110
              JA=IABS(J)
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 110
              IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10))
     &        GOTO 110
              KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
C...Special case: if only branching ratios known then use them.
              IF(MWID(34).EQ.2) THEN
                HI=0D0
                DO 105 IDC=MDCY(34,2),MDCY(34,2)+MDCY(34,3)-1
                  IF((IA.EQ.IABS(KFDP(IDC,1)).AND.JA.EQ.
     &            IABS(KFDP(IDC,2))).OR.(IA.EQ.IABS(KFDP(IDC,2))
     &            .AND.JA.EQ.IABS(KFDP(IDC,1))))
     &             HI=SHR*WDTP(IDC+1-MDCY(34,2))
  105           CONTINUE
                IF(IA.LT.10) HI=HI*FACA/9D0
              ELSE
C...Normal cross section.
                HI=HP*(PARU(133)**2+PARU(134)**2)
                IF(IA.LE.10) HI=HP*(PARU(131)**2+PARU(132)**2)*
     &          VCKM((IA+1)/2,(JA+1)/2)*FACA/3D0
              ENDIF 
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              HF=SHR*(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))
              SIGH(NCHN)=HI*FACBW*HF
  110       CONTINUE
  120     CONTINUE
 
        ELSEIF(ISUB.EQ.144) THEN
C...f + fbar' -> R
          SQMR=PMAS(41,1)**2
          CALL PYWIDT(41,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=4D0*COMFAC/((SH-SQMR)**2+HS**2)*3D0
          HP=AEM/(12D0*XW)*SH
          DO 140 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 140
            IA=IABS(I)
            DO 130 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 130
              JA=IABS(J)
              IF(I*J.GT.0.OR.IABS(IA-JA).NE.2) GOTO 130
              HI=HP
              IF(IA.LE.10) HI=HI*FACA/3D0
              HF=SHR*(WDTE(0,1)+WDTE(0,(10-(I+J))/4)+WDTE(0,4))
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=HI*FACBW*HF
  130       CONTINUE
  140     CONTINUE
 
        ELSEIF(ISUB.EQ.145) THEN
C...q + l -> LQ (leptoquark)
          SQMLQ=PMAS(42,1)**2
          CALL PYWIDT(42,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=4D0*COMFAC/((SH-SQMLQ)**2+HS**2)
          IF(ABS(SHR-PMAS(42,1)).GT.PARP(48)*PMAS(42,2)) FACBW=0D0
          HP=AEM/4D0*SH
          KFLQQ=KFDP(MDCY(42,2),1)
          KFLQL=KFDP(MDCY(42,2),2)
          DO 160 I=MMIN1,MMAX1
            IF(KFAC(1,I).EQ.0) GOTO 160
            IA=IABS(I)
            IF(IA.NE.KFLQQ.AND.IA.NE.IABS(KFLQL)) GOTO 160
            DO 150 J=MMIN2,MMAX2
              IF(KFAC(2,J).EQ.0) GOTO 150
              JA=IABS(J)
              IF(JA.NE.KFLQQ.AND.JA.NE.IABS(KFLQL)) GOTO 150
              IF(I*J.NE.KFLQQ*KFLQL) GOTO 150
              IF(JA.EQ.IA) GOTO 150
              IF(IA.EQ.KFLQQ) KCHLQ=ISIGN(1,I)
              IF(JA.EQ.KFLQQ) KCHLQ=ISIGN(1,J)
              HI=HP*PARU(151)
              HF=SHR*(WDTE(0,1)+WDTE(0,(5-KCHLQ)/2)+WDTE(0,4))
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=HI*FACBW*HF
  150       CONTINUE
  160     CONTINUE
 
        ELSEIF(ISUB.EQ.146) THEN
C...e + gamma* -> e* (excited lepton)
          KFQSTR=KFPR(ISUB,1)
          KCQSTR=PYCOMP(KFQSTR)
          KFQEXC=MOD(KFQSTR,KEXCIT)
          CALL PYWIDT(KFQSTR,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=COMFAC/((SH-PMAS(KCQSTR,1)**2)**2+HS**2)
          QF=-RTCM(43)/2D0-RTCM(44)/2D0
          FACBW=FACBW*AEM*QF**2*SH/RTCM(41)**2
          IF(ABS(SHR-PMAS(KCQSTR,1)).GT.PARP(48)*PMAS(KCQSTR,2))
     &    FACBW=0D0
          HP=SH
          DO 180 I=-KFQEXC,KFQEXC,2*KFQEXC
            DO 170 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 170
              IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 170
              HI=HP
              IF(I.GT.0) HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
              IF(I.LT.0) HF=SHR*(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=22
              ISIG(NCHN,3)=1
              SIGH(NCHN)=HI*FACBW*HF
  170       CONTINUE
  180     CONTINUE
 
        ELSEIF(ISUB.EQ.147.OR.ISUB.EQ.148) THEN
C...d + g -> d* and u + g -> u* (excited quarks)
          KFQSTR=KFPR(ISUB,1)
          KCQSTR=PYCOMP(KFQSTR)
          KFQEXC=MOD(KFQSTR,KEXCIT)
          CALL PYWIDT(KFQSTR,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=COMFAC/((SH-PMAS(KCQSTR,1)**2)**2+HS**2)
          FACBW=FACBW*AS*RTCM(45)**2*SH/(3D0*RTCM(41)**2)
          IF(ABS(SHR-PMAS(KCQSTR,1)).GT.PARP(48)*PMAS(KCQSTR,2))
     &    FACBW=0D0
          HP=SH
          DO 200 I=-KFQEXC,KFQEXC,2*KFQEXC
            DO 190 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 190
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 190
              HI=HP
              IF(I.GT.0) HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
              IF(I.LT.0) HF=SHR*(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=HI*FACBW*HF
  190       CONTINUE
  200     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.190) THEN
        IF(ISUB.EQ.162) THEN
C...q + g -> LQ + lbar; LQ=leptoquark
          SQMLQ=PMAS(42,1)**2
          FACLQ=COMFAC*FACA*PARU(151)*(AS*AEM/6D0)*(-TH/SH)*
     &    (UH2+SQMLQ**2)/(UH-SQMLQ)**2
          KFLQQ=KFDP(MDCY(42,2),1)
          DO 220 I=MMINA,MMAXA
            IF(IABS(I).NE.KFLQQ) GOTO 220
            KCHLQ=ISIGN(1,I)
            DO 210 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 210
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 210
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACLQ*WIDS(42,(5-KCHLQ)/2)
  210       CONTINUE
  220     CONTINUE
 
        ELSEIF(ISUB.EQ.163) THEN
C...g + g -> LQ + LQbar; LQ=leptoquark
          SQMLQ=PMAS(42,1)**2
          FACLQ=COMFAC*FACA*WIDS(42,1)*(AS**2/2D0)*
     &    (7D0/48D0+3D0*(UH-TH)**2/(16D0*SH2))*(1D0+2D0*SQMLQ*TH/
     &    (TH-SQMLQ)**2+2D0*SQMLQ*UH/(UH-SQMLQ)**2+4D0*SQMLQ**2/
     &    ((TH-SQMLQ)*(UH-SQMLQ)))
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 230
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
C...Since don't know proper colour flow, randomize between alternatives
          ISIG(NCHN,3)=INT(1.5D0+PYR(0))
          SIGH(NCHN)=FACLQ
  230     CONTINUE
 
        ELSEIF(ISUB.EQ.164) THEN
C...q + qbar -> LQ + LQbar; LQ=leptoquark
          DELTA=0.25D0*(SQM3-SQM4)**2/SH
          SQMLQ=0.5D0*(SQM3+SQM4)-DELTA
          TH=TH-DELTA
          UH=UH-DELTA
C          SQMLQ=PMAS(42,1)**2
          FACLQA=COMFAC*WIDS(42,1)*(AS**2/9D0)*
     &    (SH*(SH-4D0*SQMLQ)-(UH-TH)**2)/SH2
          FACLQS=COMFAC*WIDS(42,1)*((PARU(151)**2*AEM**2/8D0)*
     &    (-SH*TH-(SQMLQ-TH)**2)/TH2+(PARU(151)*AEM*AS/18D0)*
     &    ((SQMLQ-TH)*(UH-TH)+SH*(SQMLQ+TH))/(SH*TH))
          KFLQQ=KFDP(MDCY(42,2),1)
          DO 240 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 240
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACLQA
            IF(IABS(I).EQ.KFLQQ) SIGH(NCHN)=FACLQA+FACLQS
  240     CONTINUE
 
        ELSEIF(ISUB.EQ.167.OR.ISUB.EQ.168) THEN
C...q + q' -> q" + d* and q + q' -> q" + u* (excited quarks)
          KFQSTR=KFPR(ISUB,2)
          KCQSTR=PYCOMP(KFQSTR)
          KFQEXC=MOD(KFQSTR,KEXCIT)
          FACQSA=COMFAC*(SH/RTCM(41)**2)**2*(1D0-SQM4/SH)
          FACQSB=COMFAC*0.25D0*(SH/RTCM(41)**2)**2*(1D0-SQM4/SH)*
     &    (1D0+SQM4/SH)*(1D0+CTH)*(1D0+((SH-SQM4)/(SH+SQM4))*CTH)
C...Propagators: as simulated in PYOFSH and as desired
          GMMQ=PMAS(KCQSTR,1)*PMAS(KCQSTR,2)
          HBW4=GMMQ/((SQM4-PMAS(KCQSTR,1)**2)**2+GMMQ**2)
          CALL PYWIDT(KFQSTR,SQM4,WDTP,WDTE)
          GMMQC=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMQC/((SQM4-PMAS(KCQSTR,1)**2)**2+GMMQC**2)
          FACQSA=FACQSA*HBW4C/HBW4
          FACQSB=FACQSB*HBW4C/HBW4
C...Branching ratios.
          BRPOS=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)
          BRNEG=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)
          DO 260 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.6.OR.KFAC(1,I).EQ.0) GOTO 260
            DO 250 J=MMIN2,MMAX2
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.6.OR.KFAC(2,J).EQ.0) GOTO 250
              IF(IA.EQ.KFQEXC.AND.I.EQ.J) THEN
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=1
                IF(I.GT.0) SIGH(NCHN)=(4D0/3D0)*FACQSA*BRPOS
                IF(I.LT.0) SIGH(NCHN)=(4D0/3D0)*FACQSA*BRNEG
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=2
                IF(J.GT.0) SIGH(NCHN)=(4D0/3D0)*FACQSA*BRPOS
                IF(J.LT.0) SIGH(NCHN)=(4D0/3D0)*FACQSA*BRNEG
              ELSEIF((IA.EQ.KFQEXC.OR.JA.EQ.KFQEXC).AND.I*J.GT.0) THEN
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=1
                IF(JA.EQ.KFQEXC) ISIG(NCHN,3)=2
                IF(ISIG(NCHN,ISIG(NCHN,3)).GT.0) SIGH(NCHN)=FACQSA*BRPOS
                IF(ISIG(NCHN,ISIG(NCHN,3)).LT.0) SIGH(NCHN)=FACQSA*BRNEG
              ELSEIF(IA.EQ.KFQEXC.AND.I.EQ.-J) THEN
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=1
                IF(I.GT.0) SIGH(NCHN)=(8D0/3D0)*FACQSB*BRPOS
                IF(I.LT.0) SIGH(NCHN)=(8D0/3D0)*FACQSB*BRNEG
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=2
                IF(J.GT.0) SIGH(NCHN)=(8D0/3D0)*FACQSB*BRPOS
                IF(J.LT.0) SIGH(NCHN)=(8D0/3D0)*FACQSB*BRNEG
              ELSEIF(I.EQ.-J) THEN
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=1
                IF(I.GT.0) SIGH(NCHN)=FACQSB*BRPOS
                IF(I.LT.0) SIGH(NCHN)=FACQSB*BRNEG
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=2
                IF(J.GT.0) SIGH(NCHN)=FACQSB*BRPOS
                IF(J.LT.0) SIGH(NCHN)=FACQSB*BRNEG
              ELSEIF(IA.EQ.KFQEXC.OR.JA.EQ.KFQEXC) THEN
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=1
                IF(JA.EQ.KFQEXC) ISIG(NCHN,3)=2
                IF(ISIG(NCHN,ISIG(NCHN,3)).GT.0) SIGH(NCHN)=FACQSB*BRPOS
                IF(ISIG(NCHN,ISIG(NCHN,3)).LT.0) SIGH(NCHN)=FACQSB*BRNEG
              ENDIF
  250       CONTINUE
  260     CONTINUE
 
        ELSEIF(ISUB.EQ.169) THEN
C...q + qbar -> e + e* (excited lepton)
          KFQSTR=KFPR(ISUB,2)
          KCQSTR=PYCOMP(KFQSTR)
          KFQEXC=MOD(KFQSTR,KEXCIT)
          FACQSB=(COMFAC/12D0)*(SH/RTCM(41)**2)**2*(1D0-SQM4/SH)*
     &    (1D0+SQM4/SH)*(1D0+CTH)*(1D0+((SH-SQM4)/(SH+SQM4))*CTH)
C...Propagators: as simulated in PYOFSH and as desired
          GMMQ=PMAS(KCQSTR,1)*PMAS(KCQSTR,2)
          HBW4=GMMQ/((SQM4-PMAS(KCQSTR,1)**2)**2+GMMQ**2)
          CALL PYWIDT(KFQSTR,SQM4,WDTP,WDTE)
          GMMQC=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMQC/((SQM4-PMAS(KCQSTR,1)**2)**2+GMMQC**2)
          FACQSB=FACQSB*HBW4C/HBW4
C...Branching ratios.
          BRPOS=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/WDTP(0)
          BRNEG=(WDTE(0,1)+WDTE(0,3)+WDTE(0,4))/WDTP(0)
          DO 270 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.6.OR.KFAC(1,I).EQ.0) GOTO 270
            J=-I
            JA=IABS(J)
            IF(J.EQ.0.OR.JA.GT.6.OR.KFAC(2,J).EQ.0) GOTO 270
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=J
            ISIG(NCHN,3)=1
            IF(I.GT.0) SIGH(NCHN)=FACQSB*BRPOS
            IF(I.LT.0) SIGH(NCHN)=FACQSB*BRNEG
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=J
            ISIG(NCHN,3)=2
            IF(J.GT.0) SIGH(NCHN)=FACQSB*BRPOS
            IF(J.LT.0) SIGH(NCHN)=FACQSB*BRNEG
  270     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.360) THEN
        IF(ISUB.EQ.341.OR.ISUB.EQ.342) THEN
C...l + l -> H_L++/-- or H_R++/--.
          KFRES=KFPR(ISUB,1)
          KFREC=PYCOMP(KFRES)
          CALL PYWIDT(KFRES,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=8D0*COMFAC/((SH-PMAS(KFREC,1)**2)**2+HS**2)
          DO 290 I=MMIN1,MMAX1
            IA=IABS(I)
            IF((IA.NE.11.AND.IA.NE.13.AND.IA.NE.15).OR.KFAC(1,I).EQ.0)
     &      GOTO 290
            DO 280 J=MMIN2,MMAX2
              JA=IABS(J)
              IF((JA.NE.11.AND.JA.NE.13.AND.JA.NE.15).OR.KFAC(2,J).EQ.0)
     &        GOTO 280
              IF(I*J.LT.0) GOTO 280
              KCHH=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              HI=SH*PARP(181+3*((IA-11)/2)+(JA-11)/2)**2/(8D0*PARU(1))
              HF=SHR*(WDTE(0,1)+WDTE(0,(5-KCHH/2)/2)+WDTE(0,4))
              SIGH(NCHN)=HI*FACBW*HF
  280       CONTINUE
  290     CONTINUE
 
        ELSEIF(ISUB.GE.343.AND.ISUB.LE.348) THEN
C...l + gamma -> H_L++/-- l' or l + gamma -> H_R++/-- l'.
          KFRES=KFPR(ISUB,1)
          KFREC=PYCOMP(KFRES)
C...Propagators: as simulated in PYOFSH and as desired
          HBW3=PMAS(KFREC,1)*PMAS(KFREC,2)/((SQM3-PMAS(KFREC,1)**2)**2+
     &    (PMAS(KFREC,1)*PMAS(KFREC,2))**2)
          CALL PYWIDT(KFRES,SQM3,WDTP,WDTE)
          GMMC=SQRT(SQM3)*WDTP(0)
          HBW3C=GMMC/((SQM3-PMAS(KFREC,1)**2)**2+GMMC**2)
          FHCC=COMFAC*AEM*HBW3C/HBW3
          DO 310 I=MMINA,MMAXA
            IA=IABS(I)
            IF(IA.NE.11.AND.IA.NE.13.AND.IA.NE.15) GOTO 310
            SQML=PMAS(IA,1)**2
            J=ISIGN(KFPR(ISUB,2),-I)
            KCHH=ISIGN(2,KCHG(IA,1)*ISIGN(1,I))
            WIDSC=(WDTE(0,1)+WDTE(0,(5-KCHH/2)/2)+WDTE(0,4))/WDTP(0)
            SMM1=8D0*(SH+TH-SQM3)*(SH+TH-2D0*SQM3-SQML-SQM4)/
     &      (UH-SQM3)**2
            SMM2=2D0*((2D0*SQM3-3D0*SQML)*SQM4+(SQML-2D0*SQM4)*TH-
     &      (TH-SQM4)*SH)/(TH-SQM4)**2
            SMM3=2D0*((2D0*SQM3-3D0*SQM4+TH)*SQML-(2D0*SQML-SQM4+TH)*
     &      SH)/(SH-SQML)**2
            SMM12=4D0*((2D0*SQML-SQM4-2D0*SQM3+TH)*SH+(TH-3D0*SQM3-
     &      3D0*SQM4)*TH+(2D0*SQM3-2D0*SQML+3D0*SQM4)*SQM3)/
     &      ((UH-SQM3)*(TH-SQM4))
            SMM13=-4D0*((TH+SQML-2D0*SQM4)*TH-(SQM3+3D0*SQML-2D0*SQM4)*
     &      SQM3+(SQM3+3D0*SQML+TH)*SH-(TH-SQM3+SH)**2)/
     &      ((UH-SQM3)*(SH-SQML))
            SMM23=-4D0*((SQML-SQM4+SQM3)*TH-SQM3**2+SQM3*(SQML+SQM4)-
     &      3D0*SQML*SQM4-(SQML-SQM4-SQM3+TH)*SH)/
     &      ((SH-SQML)*(TH-SQM4))
            SMM=(SH/(SH-SQML))**2*(SMM1+SMM2+SMM3+SMM12+SMM13+SMM23)*
     &      PARP(181+3*((IA-11)/2)+(IABS(J)-11)/2)**2/(4D0*PARU(1))
            DO 300 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 300
              IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 300
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=22
              ISIG(NCHN,3)=0
              SIGH(NCHN)=FHCC*SMM*WIDSC
  300       CONTINUE
  310     CONTINUE
 
        ELSEIF(ISUB.EQ.349.OR.ISUB.EQ.350) THEN
C...f + fbar -> H_L++ + H_L-- or H_R++ + H_R--
          KFRES=KFPR(ISUB,1)
          KFREC=PYCOMP(KFRES)
          SQMH=PMAS(KFREC,1)**2
          GMMH=PMAS(KFREC,1)*PMAS(KFREC,2)
C...Propagators: H++/-- as simulated in PYOFSH and as desired
          HBW3=GMMH/((SQM3-SQMH)**2+GMMH**2)
          CALL PYWIDT(KFRES,SQM3,WDTP,WDTE)
          GMMH3=SQRT(SQM3)*WDTP(0)
          HBW3C=GMMH3/((SQM3-SQMH)**2+GMMH3**2)
          HBW4=GMMH/((SQM4-SQMH)**2+GMMH**2)
          CALL PYWIDT(KFRES,SQM4,WDTP,WDTE)
          GMMH4=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMH4/((SQM4-SQMH)**2+GMMH4**2)
C...Kinematical and coupling functions
          FACHH=COMFAC*(HBW3C/HBW3)*(HBW4C/HBW4)*(TH*UH-SQM3*SQM4)
          XWHH=(1D0-2D0*XWV)/(8D0*XWV*(1D0-XWV))
C...Loop over allowed flavours
          DO 320 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 320
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI+0.1D0)
            VI=AI-4D0*EI*XWV
            FCOI=1D0
            IF(IABS(I).LE.10) FCOI=FACA/3D0
            IF(ISUB.EQ.349) THEN
              HBWZ=1D0/((SH-SQMZ)**2+GMMZ**2)
              IF(IABS(I).LT.10) THEN
                DSIGHH=8D0*AEM**2*(EI**2/SH2+
     &          2D0*EI*VI*XWHH*(SH-SQMZ)*HBWZ/SH+
     &          (VI**2+AI**2)*XWHH**2*HBWZ)
              ELSE
                IAOFF=181+3*((IABS(I)-11)/2)
                HSUM=(PARP(IAOFF)**2+PARP(IAOFF+1)**2+PARP(IAOFF+2)**2)/
     &          (4D0*PARU(1))
                DSIGHH=8D0*AEM**2*(EI**2/SH2+
     &          2D0*EI*VI*XWHH*(SH-SQMZ)*HBWZ/SH+
     &          (VI**2+AI**2)*XWHH**2*HBWZ)+
     &          8D0*AEM*(EI*HSUM/(SH*TH)+
     &          (VI+AI)*XWHH*HSUM*(SH-SQMZ)*HBWZ/TH)+
     &          4D0*HSUM**2/TH2
              ENDIF
            ELSE
              IF(IABS(I).LT.10) THEN
                DSIGHH=8D0*AEM**2*EI**2/SH2
              ELSE
                IAOFF=181+3*((IABS(I)-11)/2)
                HSUM=(PARP(IAOFF)**2+PARP(IAOFF+1)**2+PARP(IAOFF+2)**2)/
     &          (4D0*PARU(1))
                DSIGHH=8D0*AEM**2*EI**2/SH2+8D0*AEM*EI*HSUM/(SH*TH)+
     &          4D0*HSUM**2/TH2
              ENDIF
            ENDIF
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACHH*FCOI*DSIGHH
  320     CONTINUE
 
        ELSEIF(ISUB.EQ.351.OR.ISUB.EQ.352) THEN
C...f + f' -> f" + f"' + H++/-- (W+/- + W+/- -> H++/-- as inner process)
          KFRES=KFPR(ISUB,1)
          KFREC=PYCOMP(KFRES)
          SQMH=PMAS(KFREC,1)**2
          IF(ISUB.EQ.351) FACNOR=PARP(190)**8*PARP(192)**2
          IF(ISUB.EQ.352) FACNOR=PARP(191)**6*2D0*
     &    PMAS(PYCOMP(9900024),1)**2
          FACWW=COMFAC*FACNOR*TAUP*VINT(2)*VINT(219)
          FACPRT=1D0/((VINT(204)**2-VINT(215))*
     &    (VINT(209)**2-VINT(216)))
          FACPRU=1D0/((VINT(204)**2+2D0*VINT(217))*
     &    (VINT(209)**2+2D0*VINT(218)))
          CALL PYWIDT(KFRES,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=(1D0/PARU(1))*VINT(2)/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(KFREC,1)).GT.PARP(48)*PMAS(KFREC,2))
     &    FACBW=0D0
          DO 340 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 340
            IF(ISUB.EQ.352.AND.IABS(I).GT.10) GOTO 340
            KCHWI=(1-2*MOD(IABS(I),2))*ISIGN(1,I)
            DO 330 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 330
              IF(ISUB.EQ.352.AND.IABS(J).GT.10) GOTO 330
              KCHWJ=(1-2*MOD(IABS(J),2))*ISIGN(1,J)
              KCHH=KCHWI+KCHWJ
              IF(IABS(KCHH).NE.2) GOTO 330
              FACLR=VINT(180+I)*VINT(180+J)
              HF=SHR*(WDTE(0,1)+WDTE(0,(5-KCHH/2)/2)+WDTE(0,4))
              IF(I.EQ.J.AND.IABS(I).GT.10) THEN
                FACPRP=0.5D0*(FACPRT+FACPRU)**2
              ELSE
                FACPRP=FACPRT**2
              ENDIF
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACLR*FACWW*FACPRP*FACBW*HF
  330       CONTINUE
  340     CONTINUE
 
        ELSEIF(ISUB.EQ.353) THEN
C...f + fbar -> Z_R0
          SQMZR=PMAS(PYCOMP(KFPR(ISUB,1)),1)**2
          CALL PYWIDT(KFPR(ISUB,1),SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=4D0*COMFAC/((SH-SQMZR)**2+HS**2)*3D0
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          HP=(AEM/(3D0*(1D0-2D0*XW)))*XWC*SH
          DO 350 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 350
            IF(IABS(I).LE.8) THEN
              EI=KCHG(IABS(I),1)/3D0
              AI=SIGN(1D0,EI+0.1D0)*(1D0-2D0*XW)
              VI=SIGN(1D0,EI+0.1D0)-4D0*EI*XW
            ELSE
              AI=-(1D0-2D0*XW)
              VI=-1D0+4D0*XW
            ENDIF
            HI=HP*(VI**2+AI**2)
            IF(IABS(I).LE.10) HI=HI*FACA/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=HI*FACBW*HF
  350     CONTINUE
 
        ELSEIF(ISUB.EQ.354) THEN
C...f + fbar' -> W_R+/-
          SQMWR=PMAS(PYCOMP(KFPR(ISUB,1)),1)**2
          CALL PYWIDT(KFPR(ISUB,1),SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=4D0*COMFAC/((SH-SQMWR)**2+HS**2)*3D0
          HP=AEM/(24D0*XW)*SH
          DO 370 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 370
            IA=IABS(I)
            DO 360 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 360
              JA=IABS(J)
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 360
              IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10))
     &        GOTO 360
              KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
              HI=HP*2D0
              IF(IA.LE.10) HI=HI*VCKM((IA+1)/2,(JA+1)/2)*FACA/3D0
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              HF=SHR*(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))
              SIGH(NCHN)=HI*FACBW*HF
  360       CONTINUE
  370     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.400) THEN
        IF(ISUB.EQ.391) THEN
C...f + fbar -> G*.
          KFGSTR=KFPR(ISUB,1)
          KCGSTR=PYCOMP(KFGSTR)
          CALL PYWIDT(KFGSTR,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          FACG=COMFAC*PARP(50)**2/(16D0*PARU(1))*SH*HF/
     &    ((SH-PMAS(KCGSTR,1)**2)**2+HS**2)
C...Modify cross section in wings of peak.
          FACG = FACG * SH**2 / PMAS(KCGSTR,1)**4
          DO 380 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 380
            HI=1D0
            IF(IABS(I).LE.10) HI=HI*FACA/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACG*HI
  380     CONTINUE
 
        ELSEIF(ISUB.EQ.392) THEN
C...g + g -> G*.
          KFGSTR=KFPR(ISUB,1)
          KCGSTR=PYCOMP(KFGSTR)
          CALL PYWIDT(KFGSTR,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          FACG=COMFAC*PARP(50)**2/(32D0*PARU(1))*SH*HF/
     &    ((SH-PMAS(KCGSTR,1)**2)**2+HS**2)
C...Modify cross section in wings of peak.
          FACG = FACG * SH**2 / PMAS(KCGSTR,1)**4
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 390
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACG
  390     CONTINUE
 
        ELSEIF(ISUB.EQ.393) THEN
C...q + qbar -> g + G*.
          KFGSTR=KFPR(ISUB,2)
          KCGSTR=PYCOMP(KFGSTR)
          FACG=COMFAC*PARP(50)**2*AS*SH/(72D0*PARU(1)*SQM4)*
     &    (4D0*(TH2+UH2)/SH2+9D0*(TH+UH)/SH+(TH2/UH+UH2/TH)/SH+
     &    3D0*(4D0+TH/UH+UH/TH)+4D0*(SH/UH+SH/TH)+
     &    2D0*SH2/(TH*UH))
C...Propagators: as simulated in PYOFSH and as desired
          GMMG=PMAS(KCGSTR,1)*PMAS(KCGSTR,2)
          HBW4=GMMG/((SQM4-PMAS(KCGSTR,1)**2)**2+GMMG**2)
          CALL PYWIDT(KFGSTR,SQM4,WDTP,WDTE)
          HS=SQRT(SQM4)*WDTP(0)
          HF=SQRT(SQM4)*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          HBW4C=HF/((SQM4-PMAS(KCGSTR,1)**2)**2+HS**2)
          FACG=FACG*HBW4C/HBW4
          DO 400 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 400
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACG
  400     CONTINUE
 
        ELSEIF(ISUB.EQ.394) THEN
C...q + g -> q + G*.
          KFGSTR=KFPR(ISUB,2)
          KCGSTR=PYCOMP(KFGSTR)
          FACG=-COMFAC*PARP(50)**2*AS*SH/(192D0*PARU(1)*SQM4)*
     &    (4D0*(SH2+UH2)/(TH*SH)+9D0*(SH+UH)/SH+SH/UH+UH2/SH2+
     &    3D0*TH*(4D0+SH/UH+UH/SH)/SH+4D0*TH2*(1D0/UH+1D0/SH)/SH+
     &    2D0*TH2*TH/(UH*SH2))
C...Propagators: as simulated in PYOFSH and as desired
          GMMG=PMAS(KCGSTR,1)*PMAS(KCGSTR,2)
          HBW4=GMMG/((SQM4-PMAS(KCGSTR,1)**2)**2+GMMG**2)
          CALL PYWIDT(KFGSTR,SQM4,WDTP,WDTE)
          HS=SQRT(SQM4)*WDTP(0)
          HF=SQRT(SQM4)*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          HBW4C=HF/((SQM4-PMAS(KCGSTR,1)**2)**2+HS**2)
          FACG=FACG*HBW4C/HBW4
          DO 420 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 420
            DO 410 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 410
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 410
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACG
  410       CONTINUE
  420     CONTINUE
 
        ELSEIF(ISUB.EQ.395) THEN
C...g + g -> g + G*.
          KFGSTR=KFPR(ISUB,2)
          KCGSTR=PYCOMP(KFGSTR)
          FACG=COMFAC*3D0*PARP(50)**2*AS*SH/(32D0*PARU(1)*SQM4)*
     &    ((TH2+TH*UH+UH2)**2/(SH2*TH*UH)+2D0*(TH2/UH+UH2/TH)/SH+
     &    3D0*(TH/UH+UH/TH)+2D0*(SH/UH+SH/TH)+SH2/(TH*UH))
C...Propagators: as simulated in PYOFSH and as desired
          GMMG=PMAS(KCGSTR,1)*PMAS(KCGSTR,2)
          HBW4=GMMG/((SQM4-PMAS(KCGSTR,1)**2)**2+GMMG**2)
          CALL PYWIDT(KFGSTR,SQM4,WDTP,WDTE)
          HS=SQRT(SQM4)*WDTP(0)
          HF=SQRT(SQM4)*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          HBW4C=HF/((SQM4-PMAS(KCGSTR,1)**2)**2+HS**2)
          FACG=FACG*HBW4C/HBW4
          IF(KFAC(1,21)*KFAC(2,21).NE.0) THEN
            NCHN=NCHN+1
            ISIG(NCHN,1)=21
            ISIG(NCHN,2)=21
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACG
          ENDIF
        ENDIF
      ELSEIF(ISUB.LE.500) THEN
        IF(ISUBSV.EQ.481) ISUB=482
c...  GENERIC 2->(1)->2
        IF(ISUB.EQ.482) THEN
          KFRES=9900001
          KCRES=PYCOMP(KFRES)
          IF(KCRES.EQ.0) RETURN
          IDCY=MDCY(KCRES,2)
          KCOL=KCHG(KCRES,2)
          KCEM=KCHG(KCRES,1)
          FACT=COMFAC
          KCF1=PYCOMP(KFPR(ISUB,1))
          KCF2=PYCOMP(KFPR(ISUB,2))
          IF(ISUBSV.EQ.481) THEN
            SQMZR=PMAS(KCRES,1)**2
            CALL PYWIDT(KFRES,SH,WDTP,WDTE)
            HS=SHR*WDTP(0)
            FACBW=SH2/((SH-SQMZR)**2+HS**2)
            FACT=FACT*FACBW
          ELSE
            SQMH=PMAS(KCF1,1)**2
            GMMH=PMAS(KCF1,1)*PMAS(KCF1,2)
C...Propagators: as simulated in PYOFSH and as desired
            HBW3=GMMH/((SQM3-SQMH)**2+GMMH**2)
            CALL PYWIDT(KFPR(ISUB,1),SQM3,WDTP,WDTE)
            GMMH3=SQRT(SQM3)*WDTP(0)
            HBW3C=GMMH3/((SQM3-SQMH)**2+GMMH3**2)
            SQMH=PMAS(KCF2,1)**2
            GMMH=PMAS(KCF2,1)*PMAS(KCF2,2)
            HBW4=GMMH/((SQM4-SQMH)**2+GMMH**2)
            CALL PYWIDT(KFPR(ISUB,2),SQM4,WDTP,WDTE)
            GMMH4=SQRT(SQM4)*WDTP(0)
            HBW4C=GMMH4/((SQM4-SQMH)**2+GMMH4**2)
            FACT=FACT*(HBW3C/HBW3)*(HBW4C/HBW4)
          ENDIF

          KCI1=ABS(PYCOMP(KFDP(IDCY,1)))
          KCI2=ABS(PYCOMP(KFDP(IDCY,2)))
          JCOL1=SIGN(KCHG(KCF1,2),KFPR(ISUB,1))
          JCOL2=SIGN(KCHG(KCF2,2),KFPR(ISUB,2))
          IF(KCOL.EQ.0) THEN
            NCOL=1
          ELSEIF(KCI1.EQ.21.AND.KCI2.EQ.21.AND.KCOL.EQ.2) THEN
            IF(JCOL1.EQ.2.AND.JCOL2.EQ.2) THEN
              NCOL=3
            ELSE
              NCOL=2
            ENDIF
          ELSEIF(KCOL.EQ.-1.OR.KCOL.EQ.1) THEN
            NCOL=2
          ELSEIF(KCI1.EQ.21.AND.KCI2.EQ.21.AND.JCOL1.EQ.0.AND.
     $      JCOL2.EQ.0) THEN
            NCOL=1
          ELSEIF(KCOL.EQ.2.AND.((JCOL1.EQ.0.AND.JCOL2.EQ.2).OR.
     $      (JCOL1.EQ.2.AND.JCOL2.EQ.0))) THEN
            NCOL=1
          ELSE
            NCOL=2
          ENDIF
          DO 440 I=MMIN1,MMAX1
            IF(KFAC(1,I).EQ.0) GOTO 440
            IP=I
            IF(IP.EQ.0) IP=21
            IA=ABS(IP)
            DO 430 J=MMIN2,MMAX2
              IF(KFAC(2,J).EQ.0) GOTO 430
              JP=J
              IF(JP.EQ.0) JP=21
              JA=ABS(JP)
              IF((IA.EQ.KCI1.AND.JA.EQ.KCI2).OR.
     $          (JA.EQ.KCI1.AND.IA.EQ.KCI2)) THEN
                KCHW=KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J)
                IF(ABS(KCHW).EQ.ABS(KCEM)) THEN
                  DO II=1,NCOL
                    NCHN=NCHN+1
                    ISIG(NCHN,1)=IP
                    ISIG(NCHN,2)=JP
                    ISIG(NCHN,3)=II
                    SIGH(NCHN)=FACT/NCOL
                  ENDDO
                ENDIF
              ENDIF
 430        CONTINUE
 440      CONTINUE
        ENDIF
      ENDIF
 
      RETURN
      END
