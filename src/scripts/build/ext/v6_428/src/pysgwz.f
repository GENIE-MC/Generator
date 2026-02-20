 
C*********************************************************************
 
C...PYSGWZ
C...Subprocess cross sections for W/Z processes,
C...except that longitudinal WW scattering is in Higgs sector.
C...Auxiliary to PYSIGH.
 
      SUBROUTINE PYSGWZ(NCHN,SIGS)
 
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
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
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
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,/PYINT1/,
     &/PYINT2/,/PYINT3/,/PYINT4/,/PYTCSM/,/PYSGCM/
C...Local arrays and complex numbers
      DIMENSION WDTP(0:400),WDTE(0:400,0:5),HGZ(6,3),HL3(3),HR3(3),
     &HL4(3),HR4(3)
      COMPLEX*16 COULCK,COULCP,COULCD,COULCR,COULCS
 
C...Differential cross section expressions.
 
      IF(ISUB.LE.20) THEN
        IF(ISUB.EQ.1) THEN
C...f + fbar -> gamma*/Z0
          MINT(61)=2
          CALL PYWIDT(23,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACZ=4D0*COMFAC*3D0
          HP0=AEM/3D0*SH
          HP1=AEM/3D0*XWC*SH
          DO 100 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 100
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            HI0=HP0
            IF(IABS(I).LE.10) HI0=HI0*FACA/3D0
            HI1=HP1
            IF(IABS(I).LE.10) HI1=HI1*FACA/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACZ*(EI**2/SH2*HI0*HP0*VINT(111)+
     &      EI*VI*(1D0-SQMZ/SH)/((SH-SQMZ)**2+HS**2)*
     &      (HI0*HP1+HI1*HP0)*VINT(112)+(VI**2+AI**2)/
     &      ((SH-SQMZ)**2+HS**2)*HI1*HP1*VINT(114))
  100     CONTINUE
 
        ELSEIF(ISUB.EQ.2) THEN
C...f + fbar' -> W+/-
          CALL PYWIDT(24,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=4D0*COMFAC/((SH-SQMW)**2+HS**2)*3D0
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
              HI=HP*2D0
              IF(IA.LE.10) HI=HI*VCKM((IA+1)/2,(JA+1)/2)*FACA/3D0
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              HF=SHR*(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))
              SIGH(NCHN)=HI*FACBW*HF
  110       CONTINUE
  120     CONTINUE
 
        ELSEIF(ISUB.EQ.15) THEN
C...f + fbar -> g + (gamma*/Z0) (q + qbar -> g + (gamma*/Z0) only)
          FACZG=COMFAC*AS*AEM*(8D0/9D0)*(TH2+UH2+2D0*SQM4*SH)/(TH*UH)
C...gamma, gamma/Z interference and Z couplings to final fermion pairs
          HFGG=0D0
          HFGZ=0D0
          HFZZ=0D0
          RADC4=1D0+PYALPS(SQM4)/PARU(1)
          DO 130 I=1,MIN(16,MDCY(23,3))
            IDC=I+MDCY(23,2)-1
            IF(MDME(IDC,1).LT.0) GOTO 130
            IMDM=0
            IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2.OR.MDME(IDC,1).EQ.4)
     &      IMDM=1
            IF(I.LE.8) THEN
              EF=KCHG(I,1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*XWV
            ELSEIF(I.LE.16) THEN
              EF=KCHG(I+2,1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*XWV
            ENDIF
            RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
            IF(4D0*RM1.LT.1D0) THEN
              FCOF=1D0
              IF(I.LE.8) FCOF=3D0*RADC4
              BE34=SQRT(MAX(0D0,1D0-4D0*RM1))
              IF(IMDM.EQ.1) THEN
                HFGG=HFGG+FCOF*EF**2*(1D0+2D0*RM1)*BE34
                HFGZ=HFGZ+FCOF*EF*VF*(1D0+2D0*RM1)*BE34
                HFZZ=HFZZ+FCOF*(VF**2*(1D0+2D0*RM1)+
     &          AF**2*(1D0-4D0*RM1))*BE34
              ENDIF
            ENDIF
  130     CONTINUE
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=(1D0/PARU(1))*GMMZ/((SQM4-SQMZ)**2+GMMZ**2)
          MINT15=MINT(15)
          MINT(15)=1
          MINT(61)=1
          CALL PYWIDT(23,SQM4,WDTP,WDTE)
          MINT(15)=MINT15
          HFAEM=(PARU(108)/PARU(2))*(2D0/3D0)
          HFGG=HFGG*HFAEM*VINT(111)/SQM4
          HFGZ=HFGZ*HFAEM*VINT(112)/SQM4
          HFZZ=HFZZ*HFAEM*VINT(114)/SQM4
C...Loop over flavours; consider full gamma/Z structure
          DO 140 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 140
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACZG*(EI**2*HFGG+EI*VI*HFGZ+
     &      (VI**2+AI**2)*HFZZ)/HBW4
  140     CONTINUE
 
        ELSEIF(ISUB.EQ.16) THEN
C...f + fbar' -> g + W+/- (q + qbar' -> g + W+/- only)
          FACWG=COMFAC*AS*AEM/XW*2D0/9D0*(TH2+UH2+2D0*SQM4*SH)/(TH*UH)
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=GMMW/((SQM4-SQMW)**2+GMMW**2)
          CALL PYWIDT(24,SQM4,WDTP,WDTE)
          GMMWC=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMWC/((SQM4-SQMW)**2+GMMWC**2)
          FACWG=FACWG*HBW4C/HBW4
          DO 160 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.10.OR.KFAC(1,I).EQ.0) GOTO 160
            DO 150 J=MMIN2,MMAX2
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.10.OR.KFAC(2,J).EQ.0) GOTO 150
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 150
              KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
              WIDSC=(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))/WDTP(0)
              FCKM=VCKM((IA+1)/2,(JA+1)/2)
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACWG*FCKM*WIDSC
  150       CONTINUE
  160     CONTINUE
 
        ELSEIF(ISUB.EQ.19) THEN
C...f + fbar -> gamma + (gamma*/Z0)
          FACGZ=COMFAC*2D0*AEM**2*(TH2+UH2+2D0*SQM4*SH)/(TH*UH)
C...gamma, gamma/Z interference and Z couplings to final fermion pairs
          HFGG=0D0
          HFGZ=0D0
          HFZZ=0D0
          RADC4=1D0+PYALPS(SQM4)/PARU(1)
          DO 170 I=1,MIN(16,MDCY(23,3))
            IDC=I+MDCY(23,2)-1
            IF(MDME(IDC,1).LT.0) GOTO 170
            IMDM=0
            IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2.OR.MDME(IDC,1).EQ.4)
     &      IMDM=1
            IF(I.LE.8) THEN
              EF=KCHG(I,1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*XWV
            ELSEIF(I.LE.16) THEN
              EF=KCHG(I+2,1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*XWV
            ENDIF
            RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
            IF(4D0*RM1.LT.1D0) THEN
              FCOF=1D0
              IF(I.LE.8) FCOF=3D0*RADC4
              BE34=SQRT(MAX(0D0,1D0-4D0*RM1))
              IF(IMDM.EQ.1) THEN
                HFGG=HFGG+FCOF*EF**2*(1D0+2D0*RM1)*BE34
                HFGZ=HFGZ+FCOF*EF*VF*(1D0+2D0*RM1)*BE34
                HFZZ=HFZZ+FCOF*(VF**2*(1D0+2D0*RM1)+
     &          AF**2*(1D0-4D0*RM1))*BE34
              ENDIF
            ENDIF
  170     CONTINUE
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=(1D0/PARU(1))*GMMZ/((SQM4-SQMZ)**2+GMMZ**2)
          MINT15=MINT(15)
          MINT(15)=1
          MINT(61)=1
          CALL PYWIDT(23,SQM4,WDTP,WDTE)
          MINT(15)=MINT15
          HFAEM=(PARU(108)/PARU(2))*(2D0/3D0)
          HFGG=HFGG*HFAEM*VINT(111)/SQM4
          HFGZ=HFGZ*HFAEM*VINT(112)/SQM4
          HFZZ=HFZZ*HFAEM*VINT(114)/SQM4
C...Loop over flavours; consider full gamma/Z structure
          DO 180 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 180
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            FCOI=1D0
            IF(IABS(I).LE.10) FCOI=FACA/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACGZ*FCOI*EI**2*(EI**2*HFGG+EI*VI*HFGZ+
     &      (VI**2+AI**2)*HFZZ)/HBW4
  180     CONTINUE
 
        ELSEIF(ISUB.EQ.20) THEN
C...f + fbar' -> gamma + W+/-
          FACGW=COMFAC*0.5D0*AEM**2/XW
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=GMMW/((SQM4-SQMW)**2+GMMW**2)
          CALL PYWIDT(24,SQM4,WDTP,WDTE)
          GMMWC=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMWC/((SQM4-SQMW)**2+GMMWC**2)
          FACGW=FACGW*HBW4C/HBW4
C...Anomalous couplings
          TERM1=(TH2+UH2+2D0*SQM4*SH)/(TH*UH)
          TERM2=0D0
          TERM3=0D0
          IF(ITCM(5).GE.1.AND.ITCM(5).LE.4) THEN
            TERM2=RTCM(46)*(TH-UH)/(TH+UH)
            TERM3=0.5D0*RTCM(46)**2*(TH*UH+(TH2+UH2)*SH/
     &      (4D0*SQMW))/(TH+UH)**2
          ENDIF
          DO 200 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.20.OR.KFAC(1,I).EQ.0) GOTO 200
            DO 190 J=MMIN2,MMAX2
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.20.OR.KFAC(2,J).EQ.0) GOTO 190
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 190
              IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10))
     &        GOTO 190
              KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
              WIDSC=(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))/WDTP(0)
              IF(IA.LE.10) THEN
                FACWR=UH/(TH+UH)-1D0/3D0
                FCKM=VCKM((IA+1)/2,(JA+1)/2)
                FCOI=FACA/3D0
              ELSE
                FACWR=-TH/(TH+UH)
                FCKM=1D0
                FCOI=1D0
              ENDIF
              FACWK=TERM1*FACWR**2+TERM2*FACWR+TERM3
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACGW*FACWK*FCOI*FCKM*WIDSC
  190       CONTINUE
  200     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.40) THEN
        IF(ISUB.EQ.22) THEN
C...f + fbar -> (gamma*/Z0) + (gamma*/Z0)
C...Kinematics dependence
          FACZZ=COMFAC*AEM**2*((TH2+UH2+2D0*(SQM3+SQM4)*SH)/(TH*UH)-
     &    SQM3*SQM4*(1D0/TH2+1D0/UH2))
C...gamma, gamma/Z interference and Z couplings to final fermion pairs
          DO 220 I=1,6
            DO 210 J=1,3
              HGZ(I,J)=0D0
  210       CONTINUE
  220     CONTINUE
          RADC3=1D0+PYALPS(SQM3)/PARU(1)
          RADC4=1D0+PYALPS(SQM4)/PARU(1)
          DO 230 I=1,MIN(16,MDCY(23,3))
            IDC=I+MDCY(23,2)-1
            IF(MDME(IDC,1).LT.0) GOTO 230
            IMDM=0
            IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2) IMDM=1
            IF(MDME(IDC,1).EQ.4.OR.MDME(IDC,1).EQ.5) IMDM=MDME(IDC,1)-2
            IF(I.LE.8) THEN
              EF=KCHG(I,1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*XWV
            ELSEIF(I.LE.16) THEN
              EF=KCHG(I+2,1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*XWV
            ENDIF
            RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM3
            IF(4D0*RM1.LT.1D0) THEN
              FCOF=1D0
              IF(I.LE.8) FCOF=3D0*RADC3
              BE34=SQRT(MAX(0D0,1D0-4D0*RM1))
              IF(IMDM.GE.1) THEN
                HGZ(1,IMDM)=HGZ(1,IMDM)+FCOF*EF**2*(1D0+2D0*RM1)*BE34
                HGZ(2,IMDM)=HGZ(2,IMDM)+FCOF*EF*VF*(1D0+2D0*RM1)*BE34
                HGZ(3,IMDM)=HGZ(3,IMDM)+FCOF*(VF**2*(1D0+2D0*RM1)+
     &          AF**2*(1D0-4D0*RM1))*BE34
              ENDIF
            ENDIF
            RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
            IF(4D0*RM1.LT.1D0) THEN
              FCOF=1D0
              IF(I.LE.8) FCOF=3D0*RADC4
              BE34=SQRT(MAX(0D0,1D0-4D0*RM1))
              IF(IMDM.GE.1) THEN
                HGZ(4,IMDM)=HGZ(4,IMDM)+FCOF*EF**2*(1D0+2D0*RM1)*BE34
                HGZ(5,IMDM)=HGZ(5,IMDM)+FCOF*EF*VF*(1D0+2D0*RM1)*BE34
                HGZ(6,IMDM)=HGZ(6,IMDM)+FCOF*(VF**2*(1D0+2D0*RM1)+
     &          AF**2*(1D0-4D0*RM1))*BE34
              ENDIF
            ENDIF
  230     CONTINUE
C...Propagators: as simulated in PYOFSH and as desired
          HBW3=(1D0/PARU(1))*GMMZ/((SQM3-SQMZ)**2+GMMZ**2)
          HBW4=(1D0/PARU(1))*GMMZ/((SQM4-SQMZ)**2+GMMZ**2)
          MINT15=MINT(15)
          MINT(15)=1
          MINT(61)=1
          CALL PYWIDT(23,SQM3,WDTP,WDTE)
          MINT(15)=MINT15
          HFAEM=(PARU(108)/PARU(2))*(2D0/3D0)
          DO 240 J=1,3
            HGZ(1,J)=HGZ(1,J)*HFAEM*VINT(111)/SQM3
            HGZ(2,J)=HGZ(2,J)*HFAEM*VINT(112)/SQM3
            HGZ(3,J)=HGZ(3,J)*HFAEM*VINT(114)/SQM3
  240     CONTINUE
          MINT15=MINT(15)
          MINT(15)=1
          MINT(61)=1
          CALL PYWIDT(23,SQM4,WDTP,WDTE)
          MINT(15)=MINT15
          HFAEM=(PARU(108)/PARU(2))*(2D0/3D0)
          DO 250 J=1,3
            HGZ(4,J)=HGZ(4,J)*HFAEM*VINT(111)/SQM4
            HGZ(5,J)=HGZ(5,J)*HFAEM*VINT(112)/SQM4
            HGZ(6,J)=HGZ(6,J)*HFAEM*VINT(114)/SQM4
  250     CONTINUE
C...Loop over flavours; separate left- and right-handed couplings
          DO 270 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 270
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            VALI=VI-AI
            VARI=VI+AI
            FCOI=1D0
            IF(IABS(I).LE.10) FCOI=FACA/3D0
            DO 260 J=1,3
              HL3(J)=EI**2*HGZ(1,J)+EI*VALI*HGZ(2,J)+VALI**2*HGZ(3,J)
              HR3(J)=EI**2*HGZ(1,J)+EI*VARI*HGZ(2,J)+VARI**2*HGZ(3,J)
              HL4(J)=EI**2*HGZ(4,J)+EI*VALI*HGZ(5,J)+VALI**2*HGZ(6,J)
              HR4(J)=EI**2*HGZ(4,J)+EI*VARI*HGZ(5,J)+VARI**2*HGZ(6,J)
  260       CONTINUE
            FACLR=HL3(1)*HL4(1)+HL3(1)*(HL4(2)+HL4(3))+
     &      HL4(1)*(HL3(2)+HL3(3))+HL3(2)*HL4(3)+HL4(2)*HL3(3)+
     &      HR3(1)*HR4(1)+HR3(1)*(HR4(2)+HR4(3))+
     &      HR4(1)*(HR3(2)+HR3(3))+HR3(2)*HR4(3)+HR4(2)*HR3(3)
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=0.5D0*FACZZ*FCOI*FACLR/(HBW3*HBW4)
  270     CONTINUE
 
        ELSEIF(ISUB.EQ.23) THEN
C...f + fbar' -> Z0 + W+/- (Z0 only, i.e. no gamma* admixture.)
          FACZW=COMFAC*0.5D0*(AEM/XW)**2
          FACZW=FACZW*WIDS(23,2)
          THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
          FACBW=1D0/((SH-SQMW)**2+GMMW**2)
          DO 290 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.20.OR.KFAC(1,I).EQ.0) GOTO 290
            DO 280 J=MMIN2,MMAX2
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.20.OR.KFAC(2,J).EQ.0) GOTO 280
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 280
              IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10))
     &        GOTO 280
              KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
              EI=KCHG(IA,1)/3D0
              AI=SIGN(1D0,EI+0.1D0)
              VI=AI-4D0*EI*XWV
              EJ=KCHG(JA,1)/3D0
              AJ=SIGN(1D0,EJ+0.1D0)
              VJ=AJ-4D0*EJ*XWV
              IF(VI+AI.GT.0) THEN
                VISAV=VI
                AISAV=AI
                VI=VJ
                AI=AJ
                VJ=VISAV
                AJ=AISAV
              ENDIF
              FCKM=1D0
              IF(IA.LE.10) FCKM=VCKM((IA+1)/2,(JA+1)/2)
              FCOI=1D0
              IF(IA.LE.10) FCOI=FACA/3D0
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACZW*FCOI*FCKM*(FACBW*((9D0-8D0*XW)/4D0*THUH+
     &        (8D0*XW-6D0)/4D0*SH*(SQM3+SQM4))+(THUH-SH*(SQM3+SQM4))*
     &        (SH-SQMW)*FACBW*0.5D0*((VJ+AJ)/TH-(VI+AI)/UH)+
     &        THUH/(16D0*XW1)*((VJ+AJ)**2/TH2+(VI+AI)**2/UH2)+
     &        SH*(SQM3+SQM4)/(8D0*XW1)*(VI+AI)*(VJ+AJ)/(TH*UH))*
     &        WIDS(24,(5-KCHW)/2)
C***Protect against slightly negative cross sections. (Reason yet to be
C***sorted out. One possibility: addition of width to the W propagator.)
              SIGH(NCHN)=MAX(0D0,SIGH(NCHN))
  280       CONTINUE
  290     CONTINUE
 
        ELSEIF(ISUB.EQ.25) THEN
C...f + fbar -> W+ + W-
C...Propagators: Z0, W+- as simulated in PYOFSH and as desired
          GMMZC=GMMZ
          HBWZC=SH**2/((SH-SQMZ)**2+GMMZC**2)
          HBW3=GMMW/((SQM3-SQMW)**2+GMMW**2)
          CALL PYWIDT(24,SQM3,WDTP,WDTE)
          GMMW3=SQRT(SQM3)*WDTP(0)
          HBW3C=GMMW3/((SQM3-SQMW)**2+GMMW3**2)
          HBW4=GMMW/((SQM4-SQMW)**2+GMMW**2)
          CALL PYWIDT(24,SQM4,WDTP,WDTE)
          GMMW4=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMW4/((SQM4-SQMW)**2+GMMW4**2)
C...Kinematical functions
          THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
          THUH34=(2D0*SH*(SQM3+SQM4)+THUH)/(SQM3*SQM4)
          GS=(((SH-SQM3-SQM4)**2-4D0*SQM3*SQM4)*THUH34+12D0*THUH)/SH2
          GT=THUH34+4D0*THUH/TH2
          GST=((SH-SQM3-SQM4)*THUH34+4D0*(SH*(SQM3+SQM4)-THUH)/TH)/SH
          GU=THUH34+4D0*THUH/UH2
          GSU=((SH-SQM3-SQM4)*THUH34+4D0*(SH*(SQM3+SQM4)-THUH)/UH)/SH
C...Common factors and couplings
          FACWW=COMFAC*(HBW3C/HBW3)*(HBW4C/HBW4)
          FACWW=FACWW*WIDS(24,1)
          CGG=AEM**2/2D0
          CGZ=AEM**2/(4D0*XW)*HBWZC*(1D0-SQMZ/SH)
          CZZ=AEM**2/(32D0*XW**2)*HBWZC
          CNG=AEM**2/(4D0*XW)
          CNZ=AEM**2/(16D0*XW**2)*HBWZC*(1D0-SQMZ/SH)
          CNN=AEM**2/(16D0*XW**2)
C...Coulomb factor for W+W- pair
          IF(MSTP(40).GE.1.AND.MSTP(40).LE.3) THEN
            COULE=(SH-4D0*SQMW)/(4D0*PMAS(24,1))
            COULP=MAX(1D-10,0.5D0*BE34*SQRT(SH))
            IF(COULE.LT.100D0*PMAS(24,2)) THEN
              COULP1=SQRT(0.5D0*PMAS(24,1)*(SQRT(COULE**2+
     &        PMAS(24,2)**2)-COULE))
            ELSE
              COULP1=SQRT(0.5D0*PMAS(24,1)*(0.5D0*PMAS(24,2)**2/COULE))
            ENDIF
            IF(COULE.GT.-100D0*PMAS(24,2)) THEN
              COULP2=SQRT(0.5D0*PMAS(24,1)*(SQRT(COULE**2+
     &        PMAS(24,2)**2)+COULE))
            ELSE
              COULP2=SQRT(0.5D0*PMAS(24,1)*(0.5D0*PMAS(24,2)**2/
     &        ABS(COULE)))
            ENDIF
            IF(MSTP(40).EQ.1) THEN
              COULDC=PARU(1)-2D0*ATAN((COULP1**2+COULP2**2-COULP**2)/
     &        MAX(1D-10,2D0*COULP*COULP1))
              FACCOU=1D0+0.5D0*PARU(101)*COULDC/MAX(1D-5,BE34)
            ELSEIF(MSTP(40).EQ.2) THEN
              COULCK=DCMPLX(DBLE(COULP1),DBLE(COULP2))
              COULCP=DCMPLX(0D0,DBLE(COULP))
              COULCD=(COULCK+COULCP)/(COULCK-COULCP)
              COULCR=1D0+DBLE(PARU(101)*SQRT(SH))/
     &        (4D0*COULCP)*LOG(COULCD)
              COULCS=DCMPLX(0D0,0D0)
              NSTP=100
              DO 300 ISTP=1,NSTP
                COULXX=(ISTP-0.5)/NSTP
                COULCS=COULCS+(1D0/COULXX)*LOG((1D0+COULXX*COULCD)/
     &          (1D0+COULXX/COULCD))
  300         CONTINUE
              COULCR=COULCR+DBLE(PARU(101)**2*SH)/(16D0*COULCP*COULCK)*
     &        (COULCS/NSTP)
              FACCOU=ABS(COULCR)**2
            ELSEIF(MSTP(40).EQ.3) THEN
              COULDC=PARU(1)-2D0*(1D0-BE34)**2*ATAN((COULP1**2+
     &        COULP2**2-COULP**2)/MAX(1D-10,2D0*COULP*COULP1))
              FACCOU=1D0+0.5D0*PARU(101)*COULDC/MAX(1D-5,BE34)
            ENDIF
          ELSEIF(MSTP(40).EQ.4) THEN
            FACCOU=1D0+0.5D0*PARU(101)*PARU(1)/MAX(1D-5,BE34)
          ELSE
            FACCOU=1D0
          ENDIF
          VINT(95)=FACCOU
          FACWW=FACWW*FACCOU
C...Loop over allowed flavours
          DO 310 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 310
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI+0.1D0)
            VI=AI-4D0*EI*XWV
            FCOI=1D0
            IF(IABS(I).LE.10) FCOI=FACA/3D0
            IF(MSTP(50).LE.0.OR.IABS(I).LE.10) THEN
              IF(AI.LT.0D0) THEN
                DSIGWW=(CGG*EI**2+CGZ*VI*EI+CZZ*(VI**2+AI**2))*GS+
     &          (CNG*EI+CNZ*(VI+AI))*GST+CNN*GT
              ELSE
                DSIGWW=(CGG*EI**2+CGZ*VI*EI+CZZ*(VI**2+AI**2))*GS-
     &          (CNG*EI+CNZ*(VI+AI))*GSU+CNN*GU
              ENDIF
            ELSE
              XMW02=0.5D0*(SQM3+SQM4)-0.25D0*(SQM3-SQM4)**2/SH
              BET=SQRT(1D0-4D0*XMW02/SH)
              GAT=1D0/SQRT(1D0-BET**2)
              STHE2=1D0-CTH**2
              AMPZG=BET**3*(16D0+(4D0*BET**2*GAT**2+3D0/GAT**2)*STHE2)
              AMPNU=BET*(2D0+BET**2*GAT**2*STHE2/2D0+
     &        2D0*BET**2*(1D0-BET**2)*STHE2/(1D0-2D0*BET*CTH+BET**2)**2)
              AMPNG=BET*((1D0+BET**2)*(4D0+BET**2*GAT**2*STHE2)+
     &        2D0*(1D0-BET**2)*(BET**2*STHE2-2D0*(1D0-BET**2))/
     &        (1D0-2D0*BET*CTH+BET**2))
              PROPI1=(0.25D0*SQMZ/XMW02)*HBWZC*(1D0-SQMZ/SH)
              PROPI2=(0.25D0*SQMZ/XMW02)**2*HBWZC
              A0=(2D0*(XMW02/SQMZ)-(1D0-BET**2)*XW)*POLL
              A1=(2D0*(XMW02/SQMZ)**2-2*XMW02/SQMZ*(1D0-BET**2)*XW)*POLL
              A2=(1D0-BET**2)**2*XW**2*(POLR+POLL)/2D0
              ATOT=AMPNU*POLL+(A1+A2)*PROPI2*AMPZG-A0*PROPI1*AMPNG
              ATOT=ATOT*CNN/SQMW*SH/BET*2D0
              DSIGWW=ATOT
            ENDIF
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACWW*FCOI*DSIGWW
  310     CONTINUE
 
        ELSEIF(ISUB.EQ.30) THEN
C...f + g -> f + (gamma*/Z0) (q + g -> q + (gamma*/Z0) only)
          FZQ=COMFAC*FACA*AS*AEM*(1D0/3D0)*(SH2+UH2+2D0*SQM4*TH)/
     &    (-SH*UH)
C...gamma, gamma/Z interference and Z couplings to final fermion pairs
          HFGG=0D0
          HFGZ=0D0
          HFZZ=0D0
          RADC4=1D0+PYALPS(SQM4)/PARU(1)
          DO 320 I=1,MIN(16,MDCY(23,3))
            IDC=I+MDCY(23,2)-1
            IF(MDME(IDC,1).LT.0) GOTO 320
            IMDM=0
            IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2.OR.MDME(IDC,1).EQ.4)
     &      IMDM=1
            IF(I.LE.8) THEN
              EF=KCHG(I,1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*XWV
            ELSEIF(I.LE.16) THEN
              EF=KCHG(I+2,1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*XWV
            ENDIF
            RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
            IF(4D0*RM1.LT.1D0) THEN
              FCOF=1D0
              IF(I.LE.8) FCOF=3D0*RADC4
              BE34=SQRT(MAX(0D0,1D0-4D0*RM1))
              IF(IMDM.EQ.1) THEN
                HFGG=HFGG+FCOF*EF**2*(1D0+2D0*RM1)*BE34
                HFGZ=HFGZ+FCOF*EF*VF*(1D0+2D0*RM1)*BE34
                HFZZ=HFZZ+FCOF*(VF**2*(1D0+2D0*RM1)+
     &          AF**2*(1D0-4D0*RM1))*BE34
              ENDIF
            ENDIF
  320     CONTINUE
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=(1D0/PARU(1))*GMMZ/((SQM4-SQMZ)**2+GMMZ**2)
          MINT15=MINT(15)
          MINT(15)=1
          MINT(61)=1
          CALL PYWIDT(23,SQM4,WDTP,WDTE)
          MINT(15)=MINT15
          HFAEM=(PARU(108)/PARU(2))*(2D0/3D0)
          HFGG=HFGG*HFAEM*VINT(111)/SQM4
          HFGZ=HFGZ*HFAEM*VINT(112)/SQM4
          HFZZ=HFZZ*HFAEM*VINT(114)/SQM4
C...Loop over flavours; consider full gamma/Z structure
          DO 340 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 340
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            FACZQ=FZQ*(EI**2*HFGG+EI*VI*HFGZ+
     &      (VI**2+AI**2)*HFZZ)/HBW4
            DO 330 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 330
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 330
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACZQ
  330       CONTINUE
  340     CONTINUE
 
        ELSEIF(ISUB.EQ.31) THEN
C...f + g -> f' + W+/- (q + g -> q' + W+/- only)
          FACWQ=COMFAC*FACA*AS*AEM/XW*1D0/12D0*
     &    (SH2+UH2+2D0*SQM4*TH)/(-SH*UH)
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=GMMW/((SQM4-SQMW)**2+GMMW**2)
          CALL PYWIDT(24,SQM4,WDTP,WDTE)
          GMMWC=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMWC/((SQM4-SQMW)**2+GMMWC**2)
          FACWQ=FACWQ*HBW4C/HBW4
          DO 360 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 360
            IA=IABS(I)
            KCHW=ISIGN(1,KCHG(IA,1)*ISIGN(1,I))
            WIDSC=(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))/WDTP(0)
            DO 350 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 350
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 350
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACWQ*VINT(180+I)*WIDSC
  350       CONTINUE
  360     CONTINUE
 
        ELSEIF(ISUB.EQ.35) THEN
C...f + gamma -> f + (gamma*/Z0)
          IF(MINT(15).EQ.22.AND.VINT(3).LT.0D0) THEN
            FZQN=SH2+UH2+2D0*(SQM4-VINT(3)**2)*TH
            FZQDTM=VINT(3)**2*SQM4-SH*(UH-VINT(4)**2)
          ELSEIF(MINT(16).EQ.22.AND.VINT(4).LT.0D0) THEN
            FZQN=SH2+UH2+2D0*(SQM4-VINT(4)**2)*TH
            FZQDTM=VINT(4)**2*SQM4-SH*(UH-VINT(3)**2)
          ELSE
            FZQN=SH2+UH2+2D0*SQM4*TH
            FZQDTM=-SH*UH
          ENDIF
          FZQN=COMFAC*2D0*AEM**2*MAX(0D0,FZQN)
C...gamma, gamma/Z interference and Z couplings to final fermion pairs
          HFGG=0D0
          HFGZ=0D0
          HFZZ=0D0
          RADC4=1D0+PYALPS(SQM4)/PARU(1)
          DO 370 I=1,MIN(16,MDCY(23,3))
            IDC=I+MDCY(23,2)-1
            IF(MDME(IDC,1).LT.0) GOTO 370
            IMDM=0
            IF(MDME(IDC,1).EQ.1.OR.MDME(IDC,1).EQ.2.OR.MDME(IDC,1).EQ.4)
     &      IMDM=1
            IF(I.LE.8) THEN
              EF=KCHG(I,1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*XWV
            ELSEIF(I.LE.16) THEN
              EF=KCHG(I+2,1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*XWV
            ENDIF
            RM1=PMAS(IABS(KFDP(IDC,1)),1)**2/SQM4
            IF(4D0*RM1.LT.1D0) THEN
              FCOF=1D0
              IF(I.LE.8) FCOF=3D0*RADC4
              BE34=SQRT(MAX(0D0,1D0-4D0*RM1))
              IF(IMDM.EQ.1) THEN
                HFGG=HFGG+FCOF*EF**2*(1D0+2D0*RM1)*BE34
                HFGZ=HFGZ+FCOF*EF*VF*(1D0+2D0*RM1)*BE34
                HFZZ=HFZZ+FCOF*(VF**2*(1D0+2D0*RM1)+
     &          AF**2*(1D0-4D0*RM1))*BE34
              ENDIF
            ENDIF
  370     CONTINUE
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=(1D0/PARU(1))*GMMZ/((SQM4-SQMZ)**2+GMMZ**2)
          MINT15=MINT(15)
          MINT(15)=1
          MINT(61)=1
          CALL PYWIDT(23,SQM4,WDTP,WDTE)
          MINT(15)=MINT15
          HFAEM=(PARU(108)/PARU(2))*(2D0/3D0)
          HFGG=HFGG*HFAEM*VINT(111)/SQM4
          HFGZ=HFGZ*HFAEM*VINT(112)/SQM4
          HFZZ=HFZZ*HFAEM*VINT(114)/SQM4
C...Loop over flavours; consider full gamma/Z structure
          DO 390 I=MMINA,MMAXA
            IF(I.EQ.0) GOTO 390
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            FACZQ=EI**2*(EI**2*HFGG+EI*VI*HFGZ+
     &      (VI**2+AI**2)*HFZZ)/HBW4
            FZQD=MAX(PMAS(IABS(I),1)**2*SQM4,FZQDTM)
            DO 380 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 380
              IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 380
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=22
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACZQ*FZQN/FZQD
  380       CONTINUE
  390     CONTINUE
 
        ELSEIF(ISUB.EQ.36) THEN
C...f + gamma -> f' + W+/-
          FWQ=COMFAC*AEM**2/(2D0*XW)*
     &    (SH2+UH2+2D0*SQM4*TH)/(SQPTH*SQM4-SH*UH)
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=GMMW/((SQM4-SQMW)**2+GMMW**2)
          CALL PYWIDT(24,SQM4,WDTP,WDTE)
          GMMWC=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMWC/((SQM4-SQMW)**2+GMMWC**2)
          FWQ=FWQ*HBW4C/HBW4
          DO 410 I=MMINA,MMAXA
            IF(I.EQ.0) GOTO 410
            IA=IABS(I)
            EIA=ABS(KCHG(IABS(I),1)/3D0)
            FACWQ=FWQ*(EIA-SH/(SH+UH))**2
            KCHW=ISIGN(1,KCHG(IA,1)*ISIGN(1,I))
            WIDSC=(WDTE(0,1)+WDTE(0,(5-KCHW)/2)+WDTE(0,4))/WDTP(0)
            DO 400 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,22).EQ.0) GOTO 400
              IF(ISDE.EQ.2.AND.KFAC(1,22)*KFAC(2,I).EQ.0) GOTO 400
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=22
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACWQ*VINT(180+I)*WIDSC
  400       CONTINUE
  410     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.100) THEN
        IF(ISUB.EQ.69) THEN
C...gamma + gamma -> W+ + W-
          SQMWE=MAX(0.5D0*SQMW,SQRT(SQM3*SQM4))
          FPROP=SH2/((SQMWE-TH)*(SQMWE-UH))
          FACWW=COMFAC*6D0*AEM**2*(1D0-FPROP*(4D0/3D0+2D0*SQMWE/SH)+
     &    FPROP**2*(2D0/3D0+2D0*(SQMWE/SH)**2))*WIDS(24,1)
          IF(KFAC(1,22)*KFAC(2,22).EQ.0) GOTO 420
          NCHN=NCHN+1
          ISIG(NCHN,1)=22
          ISIG(NCHN,2)=22
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACWW
  420     CONTINUE
 
        ELSEIF(ISUB.EQ.70) THEN
C...gamma + W+/- -> Z0 + W+/-
          SQMWE=MAX(0.5D0*SQMW,SQRT(SQM3*SQM4))
          FPROP=(TH-SQMWE)**2/(-SH*(SQMWE-UH))
          FACZW=COMFAC*6D0*AEM**2*(XW1/XW)*
     &    (1D0-FPROP*(4D0/3D0+2D0*SQMWE/(TH-SQMWE))+
     &    FPROP**2*(2D0/3D0+2D0*(SQMWE/(TH-SQMWE))**2))*WIDS(23,2)
          DO 440 KCHW=1,-1,-2
            DO 430 ISDE=1,2
              IF(KFAC(ISDE,22)*KFAC(3-ISDE,24*KCHW).EQ.0) GOTO 430
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=22
              ISIG(NCHN,3-ISDE)=24*KCHW
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACZW*WIDS(24,(5-KCHW)/2)
  430       CONTINUE
  440     CONTINUE
        ENDIF
      ENDIF
 
      RETURN
      END
