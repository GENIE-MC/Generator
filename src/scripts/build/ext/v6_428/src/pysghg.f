 
C*********************************************************************
 
C...PYSGHG
C...Subprocess cross sections for Higgs processes,
C...except Higgs pairs in PYSGSU, but including WW scattering.
C...Auxiliary to PYSIGH.
 
      SUBROUTINE PYSGHG(NCHN,SIGS)
 
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
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSGCM/ISUB,ISUBSV,MMIN1,MMAX1,MMIN2,MMAX2,MMINA,MMAXA,
     &KFAC(2,-40:40),COMFAC,FACK,FACA,SH,TH,UH,SH2,TH2,UH2,SQM3,SQM4,
     &SHR,SQPTH,TAUP,BE34,CTH,X(2),SQMZ,SQMW,GMMZ,GMMW,
     &AEM,AS,XW,XW1,XWC,XWV,POLL,POLR,POLLL,POLRR
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYPARS/,/PYINT1/,/PYINT2/,
     &/PYINT3/,/PYINT4/,/PYSUBS/,/PYMSSM/,/PYSGCM/
C...Local arrays and complex variables
      DIMENSION WDTP(0:400),WDTE(0:400,0:5)
      COMPLEX*16 A004,A204,A114,A00U,A20U,A11U
      COMPLEX*16 CIGTOT,CIZTOT,F0ALP,F1ALP,F2ALP,F0BET,F1BET,F2BET,FIF
 
C...Convert H or A process into equivalent h one
      IHIGG=1
      KFHIGG=25
      IF(ISUB.EQ.401.OR.ISUB.EQ.402) THEN
         KFHIGG=KFPR(ISUB,1)
      END IF
      IF((ISUB.GE.151.AND.ISUB.LE.160).OR.(ISUB.GE.171.AND.
     &ISUB.LE.190)) THEN
        IHIGG=2
        IF(MOD(ISUB-1,10).GE.5) IHIGG=3
        KFHIGG=33+IHIGG
        IF(ISUB.EQ.151.OR.ISUB.EQ.156) ISUB=3
        IF(ISUB.EQ.152.OR.ISUB.EQ.157) ISUB=102
        IF(ISUB.EQ.153.OR.ISUB.EQ.158) ISUB=103
        IF(ISUB.EQ.171.OR.ISUB.EQ.176) ISUB=24
        IF(ISUB.EQ.172.OR.ISUB.EQ.177) ISUB=26
        IF(ISUB.EQ.173.OR.ISUB.EQ.178) ISUB=123
        IF(ISUB.EQ.174.OR.ISUB.EQ.179) ISUB=124
        IF(ISUB.EQ.181.OR.ISUB.EQ.186) ISUB=121
        IF(ISUB.EQ.182.OR.ISUB.EQ.187) ISUB=122
        IF(ISUB.EQ.183.OR.ISUB.EQ.188) ISUB=111
        IF(ISUB.EQ.184.OR.ISUB.EQ.189) ISUB=112
        IF(ISUB.EQ.185.OR.ISUB.EQ.190) ISUB=113
      ENDIF
      SQMH=PMAS(KFHIGG,1)**2
      GMMH=PMAS(KFHIGG,1)*PMAS(KFHIGG,2)
 
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
      IF((MSTP(46).GE.3.AND.MSTP(46).LE.6).AND.(ISUB.EQ.71.OR.ISUB.EQ.
     &72.OR.ISUB.EQ.73.OR.ISUB.EQ.76.OR.ISUB.EQ.77)) THEN
C...Calculate M_R and N_R functions for Higgs-like and QCD-like models
        IF(MSTP(46).LE.4) THEN
          HDTLH=LOG(PMAS(25,1)/PARP(44))
          HDTMR=(4.5D0*PARU(1)/SQRT(3D0)-74D0/9D0)/8D0+HDTLH/12D0
          HDTNR=-1D0/18D0+HDTLH/6D0
        ELSE
          HDTNM=0.125D0*(1D0/(288D0*PARU(1)**2)+(PARP(47)/PARP(45))**2)
          HDTLQ=LOG(PARP(45)/PARP(44))
          HDTMR=-(4D0*PARU(1))**2*0.5D0*HDTNM+HDTLQ/12D0
          HDTNR=(4D0*PARU(1))**2*HDTNM+HDTLQ/6D0
        ENDIF
 
C...Calculate lowest and next-to-lowest order partial wave amplitudes
        HDTV=1D0/(16D0*PARU(1)*PARP(47)**2)
        A00L=DBLE(HDTV*SH)
        A20L=-0.5D0*A00L
        A11L=A00L/6D0
        HDTLS=LOG(SH/PARP(44)**2)
        A004=DBLE((HDTV*SH)**2/(4D0*PARU(1)))*
     &  CMPLX(DBLE((176D0*HDTMR+112D0*HDTNR)/3D0+11D0/27D0-
     &  (50D0/9D0)*HDTLS),DBLE(4D0*PARU(1)))
        A204=DBLE((HDTV*SH)**2/(4D0*PARU(1)))*
     &  CMPLX(DBLE(32D0*(HDTMR+2D0*HDTNR)/3D0+25D0/54D0-
     &  (20D0/9D0)*HDTLS),DBLE(PARU(1)))
        A114=DBLE((HDTV*SH)**2/(6D0*PARU(1)))*
     &  CMPLX(DBLE(4D0*(-2D0*HDTMR+HDTNR)-1D0/18D0),DBLE(PARU(1)/6D0))
 
C...Unitarize partial wave amplitudes with Pade or K-matrix method
        IF(MSTP(46).EQ.3.OR.MSTP(46).EQ.5) THEN
          A00U=A00L/(1D0-A004/A00L)
          A20U=A20L/(1D0-A204/A20L)
          A11U=A11L/(1D0-A114/A11L)
        ELSE
          A00U=(A00L+DBLE(A004))/(1D0-DCMPLX(0.D0,A00L+DBLE(A004)))
          A20U=(A20L+DBLE(A204))/(1D0-DCMPLX(0.D0,A20L+DBLE(A204)))
          A11U=(A11L+DBLE(A114))/(1D0-DCMPLX(0.D0,A11L+DBLE(A114)))
        ENDIF
      ENDIF
 
C...Differential cross section expressions.
 
      IF(ISUB.LE.60) THEN
        IF(ISUB.EQ.3) THEN
C...f + fbar -> h0 (or H0, or A0)
          CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=4D0*COMFAC/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(KFHIGG,1)).GT.PARP(48)*PMAS(KFHIGG,2))
     &    FACBW=0D0
          HP=AEM/(8D0*XW)*SH/SQMW*SH
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          DO 100 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 100
            IA=IABS(I)
            RMQ=PYMRUN(IA,SH)**2/SH
            HI=HP*RMQ
            IF(IA.LE.10) HI=HP*RMQ*FACA/3D0
            IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
              IKFI=1
              IF(IA.LE.10.AND.MOD(IA,2).EQ.0) IKFI=2
              IF(IA.GT.10) IKFI=3
              HI=HI*PARU(150+10*IHIGG+IKFI)**2
              IF(IMSS(1).NE.0.AND.IA.EQ.5) THEN
                HI=HI/(1D0+RMSS(41))**2
                IF(IHIGG.NE.3) THEN
                  HI=HI*(1D0+RMSS(41)*PARU(152+10*IHIGG)/
     &            PARU(151+10*IHIGG))**2
                ENDIF
              ENDIF
            ENDIF
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=HI*FACBW*HF
  100     CONTINUE
 
        ELSEIF(ISUB.EQ.5) THEN
C...Z0 + Z0 -> h0
          CALL PYWIDT(25,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=4D0*COMFAC/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(25,1)).GT.PARP(48)*PMAS(25,2)) FACBW=0D0
          HP=AEM/(8D0*XW)*SH/SQMW*SH
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          HI=HP/4D0
          FACI=8D0/(PARU(1)**2*XW1)*(AEM*XWC)**2
          DO 120 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 120
            DO 110 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 110
              EI=KCHG(IABS(I),1)/3D0
              AI=SIGN(1D0,EI)
              VI=AI-4D0*EI*XWV
              EJ=KCHG(IABS(J),1)/3D0
              AJ=SIGN(1D0,EJ)
              VJ=AJ-4D0*EJ*XWV
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACI*(VI**2+AI**2)*(VJ**2+AJ**2)*HI*FACBW*HF
  110       CONTINUE
  120     CONTINUE
 
        ELSEIF(ISUB.EQ.8) THEN
C...W+ + W- -> h0
          CALL PYWIDT(25,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=4D0*COMFAC/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(25,1)).GT.PARP(48)*PMAS(25,2)) FACBW=0D0
          HP=AEM/(8D0*XW)*SH/SQMW*SH
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          HI=HP/2D0
          FACI=1D0/(4D0*PARU(1)**2)*(AEM/XW)**2
          DO 140 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 140
            EI=SIGN(1D0,DBLE(I))*KCHG(IABS(I),1)
            DO 130 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 130
              EJ=SIGN(1D0,DBLE(J))*KCHG(IABS(J),1)
              IF(EI*EJ.GT.0D0) GOTO 130
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACI*VINT(180+I)*VINT(180+J)*HI*FACBW*HF
  130       CONTINUE
  140     CONTINUE
 
        ELSEIF(ISUB.EQ.24) THEN
C...f + fbar -> Z0 + h0 (or H0, or A0)
C...Propagators: Z0, h0 as simulated in PYOFSH and as desired
          HBW3=GMMZ/((SQM3-SQMZ)**2+GMMZ**2)
          CALL PYWIDT(23,SQM3,WDTP,WDTE)
          GMMZ3=SQRT(SQM3)*WDTP(0)
          HBW3C=GMMZ3/((SQM3-SQMZ)**2+GMMZ3**2)
          HBW4=GMMH/((SQM4-SQMH)**2+GMMH**2)
          CALL PYWIDT(KFHIGG,SQM4,WDTP,WDTE)
          GMMH4=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMH4/((SQM4-SQMH)**2+GMMH4**2)
          THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
          FACHZ=COMFAC*(HBW3C/HBW3)*(HBW4C/HBW4)*8D0*(AEM*XWC)**2*
     &    (THUH+2D0*SH*SQM3)/((SH-SQMZ)**2+GMMZ**2)
          FACHZ=FACHZ*WIDS(23,2)*WIDS(KFHIGG,2)
          IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACHZ=FACHZ*
     &    PARU(154+10*IHIGG)**2
          DO 150 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 150
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            FCOI=1D0
            IF(IABS(I).LE.10) FCOI=FACA/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACHZ*FCOI*(VI**2+AI**2)
  150     CONTINUE
 
        ELSEIF(ISUB.EQ.26) THEN
C...f + fbar' -> W+/- + h0 (or H0, or A0)
C...Propagators: W+-, h0 as simulated in PYOFSH and as desired
          HBW3=GMMW/((SQM3-SQMW)**2+GMMW**2)
          CALL PYWIDT(24,SQM3,WDTP,WDTE)
          GMMW3=SQRT(SQM3)*WDTP(0)
          HBW3C=GMMW3/((SQM3-SQMW)**2+GMMW3**2)
          HBW4=GMMH/((SQM4-SQMH)**2+GMMH**2)
          CALL PYWIDT(KFHIGG,SQM4,WDTP,WDTE)
          GMMH4=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMH4/((SQM4-SQMH)**2+GMMH4**2)
          THUH=MAX(TH*UH-SQM3*SQM4,SH*CKIN(3)**2)
          FACHW=COMFAC*0.125D0*(AEM/XW)**2*(THUH+2D0*SH*SQM3)/
     &    ((SH-SQMW)**2+GMMW**2)*(HBW3C/HBW3)*(HBW4C/HBW4)
          FACHW=FACHW*WIDS(KFHIGG,2)
          IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACHW=FACHW*
     &    PARU(155+10*IHIGG)**2
          DO 170 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.20.OR.KFAC(1,I).EQ.0) GOTO 170
            DO 160 J=MMIN2,MMAX2
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.20.OR.KFAC(1,J).EQ.0) GOTO 160
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 160
              IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10))
     &        GOTO 160
              KCHW=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
              FCKM=1D0
              IF(IA.LE.10) FCKM=VCKM((IA+1)/2,(JA+1)/2)
              FCOI=1D0
              IF(IA.LE.10) FCOI=FACA/3D0
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACHW*FCOI*FCKM*WIDS(24,(5-KCHW)/2)
  160       CONTINUE
  170     CONTINUE
 
        ELSEIF(ISUB.EQ.32) THEN
C...f + g -> f + h0 (q + g -> q + h0 only)
          FHCQ=COMFAC*FACA*AS*AEM/XW*1D0/24D0
C...H propagator: as simulated in PYOFSH and as desired
          SQMHC=PMAS(25,1)**2
          GMMHC=PMAS(25,1)*PMAS(25,2)
          HBW4=GMMHC/((SQM4-SQMHC)**2+GMMHC**2)
          CALL PYWIDT(25,SQM4,WDTP,WDTE)
          GMMHCC=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMHCC/((SQM4-SQMHC)**2+GMMHCC**2)
          FHCQ=FHCQ*HBW4C/HBW4
          DO 190 I=MMINA,MMAXA
            IA=IABS(I)
            IF(IA.NE.5) GOTO 190
            SQML=PYMRUN(IA,SH)**2
            SQMQ=PMAS(IA,1)**2
            FACHCQ=FHCQ*SQML/SQMW*
     &      (SH/(SQMQ-UH)+2D0*SQMQ*(SQM4-UH)/(SQMQ-UH)**2+(SQMQ-UH)/SH-
     &      2D0*SQMQ/(SQMQ-UH)+2D0*(SQM4-UH)/(SQMQ-UH)*
     &      (SQM4-SQMQ-SH)/SH)
            DO 180 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 180
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 180
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACHCQ*WIDS(25,2)
  180       CONTINUE
  190     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.80) THEN
        IF(ISUB.EQ.71) THEN
C...Z0 + Z0 -> Z0 + Z0
          IF(SH.LE.4.01D0*SQMZ) GOTO 220
 
          IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons
            BE2=1D0-4D0*SQMZ/SH
            TH=-0.5D0*SH*BE2*(1D0-CTH)
            UH=-0.5D0*SH*BE2*(1D0+CTH)
            IF(MAX(TH,UH).GT.-1D0) GOTO 220
            SHANG=1D0/XW1*SQMW/SQMZ*(1D0+BE2)**2
            ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
            ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
            THANG=1D0/XW1*SQMW/SQMZ*(BE2-CTH)**2
            ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG
            ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG
            UHANG=1D0/XW1*SQMW/SQMZ*(BE2+CTH)**2
            AUHRE=(UH-SQMH)/((UH-SQMH)**2+GMMH**2)*UHANG
            AUHIM=-GMMH/((UH-SQMH)**2+GMMH**2)*UHANG
            FACZZ=COMFAC*1D0/(4096D0*PARU(1)**2*16D0*XW1**2)*
     &      (AEM/XW)**4*(SH/SQMW)**2*(SQMZ/SQMW)*SH2
            IF(MSTP(46).LE.0) FACZZ=FACZZ*(ASHRE**2+ASHIM**2)
            IF(MSTP(46).EQ.1) FACZZ=FACZZ*((ASHRE+ATHRE+AUHRE)**2+
     &      (ASHIM+ATHIM+AUHIM)**2)
            IF(MSTP(46).EQ.2) FACZZ=0D0
 
          ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
            FACZZ=COMFAC*(AEM/(16D0*PARU(1)*XW*XW1))**2*(64D0/9D0)*
     &      ABS(A00U+2D0*A20U)**2
          ENDIF
          FACZZ=FACZZ*WIDS(23,1)
 
          DO 210 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 210
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            AVI=AI**2+VI**2
            DO 200 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 200
              EJ=KCHG(IABS(J),1)/3D0
              AJ=SIGN(1D0,EJ)
              VJ=AJ-4D0*EJ*XWV
              AVJ=AJ**2+VJ**2
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=0.5D0*FACZZ*AVI*AVJ
  200       CONTINUE
  210     CONTINUE
  220     CONTINUE
 
        ELSEIF(ISUB.EQ.72) THEN
C...Z0 + Z0 -> W+ + W-
          IF(SH.LE.4.01D0*SQMZ) GOTO 250
 
          IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons
            BE2=SQRT((1D0-4D0*SQMW/SH)*(1D0-4D0*SQMZ/SH))
            CTH2=CTH**2
            TH=-0.5D0*SH*(1D0-2D0*(SQMW+SQMZ)/SH-BE2*CTH)
            UH=-0.5D0*SH*(1D0-2D0*(SQMW+SQMZ)/SH+BE2*CTH)
            IF(MAX(TH,UH).GT.-1D0) GOTO 250
            SHANG=4D0*SQRT(SQMW/(SQMZ*XW1))*(1D0-2D0*SQMW/SH)*
     &      (1D0-2D0*SQMZ/SH)
            ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
            ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
            ATWRE=XW1/SQMZ*SH/(TH-SQMW)*((CTH-BE2)**2*(3D0/2D0+BE2/2D0*
     &      CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4D0*
     &      ((SQMW+SQMZ)/SH*(1D0-3D0*CTH2)+8D0*SQMW*SQMZ/SH2*
     &      (2D0*CTH2-1D0)+4D0*(SQMW**2+SQMZ**2)/SH2*CTH2+
     &      2D0*(SQMW+SQMZ)/SH*BE2*CTH))
            ATWIM=0D0
            AUWRE=XW1/SQMZ*SH/(UH-SQMW)*((CTH+BE2)**2*(3D0/2D0-BE2/2D0*
     &      CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4D0*
     &      ((SQMW+SQMZ)/SH*(1D0-3D0*CTH2)+8D0*SQMW*SQMZ/SH2*
     &      (2D0*CTH2-1D0)+4D0*(SQMW**2+SQMZ**2)/SH2*CTH2-
     &      2D0*(SQMW+SQMZ)/SH*BE2*CTH))
            AUWIM=0D0
            A4RE=2D0*XW1/SQMZ*(3D0-CTH2-4D0*(SQMW+SQMZ)/SH)
            A4IM=0D0
            FACWW=COMFAC*1D0/(4096D0*PARU(1)**2*16D0*XW1**2)*
     &      (AEM/XW)**4*(SH/SQMW)**2*(SQMZ/SQMW)*SH2
            IF(MSTP(46).LE.0) FACWW=FACWW*(ASHRE**2+ASHIM**2)
            IF(MSTP(46).EQ.1) FACWW=FACWW*((ASHRE+ATWRE+AUWRE+A4RE)**2+
     &      (ASHIM+ATWIM+AUWIM+A4IM)**2)
            IF(MSTP(46).EQ.2) FACWW=FACWW*((ATWRE+AUWRE+A4RE)**2+
     &      (ATWIM+AUWIM+A4IM)**2)
 
          ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
            FACWW=COMFAC*(AEM/(16D0*PARU(1)*XW*XW1))**2*(64D0/9D0)*
     &      ABS(A00U-A20U)**2
          ENDIF
          FACWW=FACWW*WIDS(24,1)
 
          DO 240 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 240
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            AVI=AI**2+VI**2
            DO 230 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 230
              EJ=KCHG(IABS(J),1)/3D0
              AJ=SIGN(1D0,EJ)
              VJ=AJ-4D0*EJ*XWV
              AVJ=AJ**2+VJ**2
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACWW*AVI*AVJ
  230       CONTINUE
  240     CONTINUE
  250     CONTINUE
 
        ELSEIF(ISUB.EQ.73) THEN
C...Z0 + W+/- -> Z0 + W+/-
          IF(SH.LE.2D0*SQMZ+2D0*SQMW) GOTO 280
 
          IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons
            BE2=1D0-2D0*(SQMZ+SQMW)/SH+((SQMZ-SQMW)/SH)**2
            EP1=1D0-(SQMZ-SQMW)/SH
            EP2=1D0+(SQMZ-SQMW)/SH
            TH=-0.5D0*SH*BE2*(1D0-CTH)
            UH=(SQMZ-SQMW)**2/SH-0.5D0*SH*BE2*(1D0+CTH)
            IF(MAX(TH,UH).GT.-1D0) GOTO 280
            THANG=(BE2-EP1*CTH)*(BE2-EP2*CTH)
            ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG
            ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG
            ASWRE=-XW1/SQMZ*SH/(SH-SQMW)*(-BE2*(EP1+EP2)**4*CTH+
     &      1D0/4D0*(BE2+EP1*EP2)**2*((EP1-EP2)**2-4D0*BE2*CTH)+
     &      2D0*BE2*(BE2+EP1*EP2)*(EP1+EP2)**2*CTH-
     &      1D0/16D0*SH/SQMW*(EP1**2-EP2**2)**2*(BE2+EP1*EP2)**2)
            ASWIM=0D0
            AUWRE=XW1/SQMZ*SH/(UH-SQMW)*(-BE2*(EP2+EP1*CTH)*
     &      (EP1+EP2*CTH)*(BE2+EP1*EP2)+BE2*(EP2+EP1*CTH)*
     &      (BE2+EP1*EP2*CTH)*(2D0*EP2-EP2*CTH+EP1)-
     &      BE2*(EP2+EP1*CTH)**2*(BE2-EP2**2*CTH)-1D0/8D0*
     &      (BE2+EP1*EP2*CTH)**2*((EP1+EP2)**2+2D0*BE2*(1D0-CTH))+
     &      1D0/32D0*SH/SQMW*(BE2+EP1*EP2*CTH)**2*
     &      (EP1**2-EP2**2)**2-BE2*(EP1+EP2*CTH)*(EP2+EP1*CTH)*
     &      (BE2+EP1*EP2)+BE2*(EP1+EP2*CTH)*(BE2+EP1*EP2*CTH)*
     &      (2D0*EP1-EP1*CTH+EP2)-BE2*(EP1+EP2*CTH)**2*
     &      (BE2-EP1**2*CTH)-1D0/8D0*(BE2+EP1*EP2*CTH)**2*
     &      ((EP1+EP2)**2+2D0*BE2*(1D0-CTH))+1D0/32D0*SH/SQMW*
     &      (BE2+EP1*EP2*CTH)**2*(EP1**2-EP2**2)**2)
            AUWIM=0D0
            A4RE=XW1/SQMZ*(EP1**2*EP2**2*(CTH**2-1D0)-
     &      2D0*BE2*(EP1**2+EP2**2+EP1*EP2)*CTH-2D0*BE2*EP1*EP2)
            A4IM=0D0
            FACZW=COMFAC*1D0/(4096D0*PARU(1)**2*4D0*XW1)*(AEM/XW)**4*
     &      (SH/SQMW)**2*SQRT(SQMZ/SQMW)*SH2
            IF(MSTP(46).LE.0) FACZW=0D0
            IF(MSTP(46).EQ.1) FACZW=FACZW*((ATHRE+ASWRE+AUWRE+A4RE)**2+
     &      (ATHIM+ASWIM+AUWIM+A4IM)**2)
            IF(MSTP(46).EQ.2) FACZW=FACZW*((ASWRE+AUWRE+A4RE)**2+
     &      (ASWIM+AUWIM+A4IM)**2)
 
          ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
            FACZW=COMFAC*AEM**2/(64D0*PARU(1)**2*XW**2*XW1)*16D0*
     &      ABS(A20U+3D0*A11U*DBLE(CTH))**2
          ENDIF
          FACZW=FACZW*WIDS(23,2)
 
          DO 270 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 270
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            AVI=AI**2+VI**2
            KCHWI=ISIGN(1,KCHG(IABS(I),1)*ISIGN(1,I))
            DO 260 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 260
              EJ=KCHG(IABS(J),1)/3D0
              AJ=SIGN(1D0,EJ)
              VJ=AI-4D0*EJ*XWV
              AVJ=AJ**2+VJ**2
              KCHWJ=ISIGN(1,KCHG(IABS(J),1)*ISIGN(1,J))
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACZW*AVI*VINT(180+J)*WIDS(24,(5-KCHWJ)/2)
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=2
              SIGH(NCHN)=FACZW*VINT(180+I)*WIDS(24,(5-KCHWI)/2)*AVJ
  260       CONTINUE
  270     CONTINUE
  280     CONTINUE
 
        ELSEIF(ISUB.EQ.75) THEN
C...W+ + W- -> gamma + gamma
 
        ELSEIF(ISUB.EQ.76) THEN
C...W+ + W- -> Z0 + Z0
          IF(SH.LE.4.01D0*SQMZ) GOTO 310
 
          IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons
            BE2=SQRT((1D0-4D0*SQMW/SH)*(1D0-4D0*SQMZ/SH))
            CTH2=CTH**2
            TH=-0.5D0*SH*(1D0-2D0*(SQMW+SQMZ)/SH-BE2*CTH)
            UH=-0.5D0*SH*(1D0-2D0*(SQMW+SQMZ)/SH+BE2*CTH)
            IF(MAX(TH,UH).GT.-1D0) GOTO 310
            SHANG=4D0*SQRT(SQMW/(SQMZ*XW1))*(1D0-2D0*SQMW/SH)*
     &      (1D0-2D0*SQMZ/SH)
            ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
            ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
            ATWRE=XW1/SQMZ*SH/(TH-SQMW)*((CTH-BE2)**2*(3D0/2D0+BE2/2D0*
     &      CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4D0*
     &      ((SQMW+SQMZ)/SH*(1D0-3D0*CTH2)+8D0*SQMW*SQMZ/SH2*
     &      (2D0*CTH2-1D0)+4D0*(SQMW**2+SQMZ**2)/SH2*CTH2+
     &      2D0*(SQMW+SQMZ)/SH*BE2*CTH))
            ATWIM=0D0
            AUWRE=XW1/SQMZ*SH/(UH-SQMW)*((CTH+BE2)**2*(3D0/2D0-BE2/2D0*
     &      CTH-(SQMW+SQMZ)/SH+(SQMW-SQMZ)**2/(SH*SQMW))+4D0*
     &      ((SQMW+SQMZ)/SH*(1D0-3D0*CTH2)+8D0*SQMW*SQMZ/SH2*
     &      (2D0*CTH2-1D0)+4D0*(SQMW**2+SQMZ**2)/SH2*CTH2-
     &      2D0*(SQMW+SQMZ)/SH*BE2*CTH))
            AUWIM=0D0
            A4RE=2D0*XW1/SQMZ*(3D0-CTH2-4D0*(SQMW+SQMZ)/SH)
            A4IM=0D0
            FACZZ=COMFAC*1D0/(4096D0*PARU(1)**2)*(AEM/XW)**4*
     &      (SH/SQMW)**2*SH2
            IF(MSTP(46).LE.0) FACZZ=FACZZ*(ASHRE**2+ASHIM**2)
            IF(MSTP(46).EQ.1) FACZZ=FACZZ*((ASHRE+ATWRE+AUWRE+A4RE)**2+
     &      (ASHIM+ATWIM+AUWIM+A4IM)**2)
            IF(MSTP(46).EQ.2) FACZZ=FACZZ*((ATWRE+AUWRE+A4RE)**2+
     &      (ATWIM+AUWIM+A4IM)**2)
 
          ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
            FACZZ=COMFAC*(AEM/(4D0*PARU(1)*XW))**2*(64D0/9D0)*
     &      ABS(A00U-A20U)**2
          ENDIF
          FACZZ=FACZZ*WIDS(23,1)
 
          DO 300 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 300
            EI=SIGN(1D0,DBLE(I))*KCHG(IABS(I),1)
            DO 290 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 290
              EJ=SIGN(1D0,DBLE(J))*KCHG(IABS(J),1)
              IF(EI*EJ.GT.0D0) GOTO 290
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=0.5D0*FACZZ*VINT(180+I)*VINT(180+J)
  290       CONTINUE
  300     CONTINUE
  310     CONTINUE
 
        ELSEIF(ISUB.EQ.77) THEN
C...W+/- + W+/- -> W+/- + W+/-
          IF(SH.LE.4.01D0*SQMW) GOTO 340
 
          IF(MSTP(46).LE.2) THEN
C...Exact scattering ME:s for on-mass-shell gauge bosons
            BE2=1D0-4D0*SQMW/SH
            BE4=BE2**2
            CTH2=CTH**2
            CTH3=CTH**3
            TH=-0.5D0*SH*BE2*(1D0-CTH)
            UH=-0.5D0*SH*BE2*(1D0+CTH)
            IF(MAX(TH,UH).GT.-1D0) GOTO 340
            SHANG=(1D0+BE2)**2
            ASHRE=(SH-SQMH)/((SH-SQMH)**2+GMMH**2)*SHANG
            ASHIM=-GMMH/((SH-SQMH)**2+GMMH**2)*SHANG
            THANG=(BE2-CTH)**2
            ATHRE=(TH-SQMH)/((TH-SQMH)**2+GMMH**2)*THANG
            ATHIM=-GMMH/((TH-SQMH)**2+GMMH**2)*THANG
            UHANG=(BE2+CTH)**2
            AUHRE=(UH-SQMH)/((UH-SQMH)**2+GMMH**2)*UHANG
            AUHIM=-GMMH/((UH-SQMH)**2+GMMH**2)*UHANG
            SGZANG=1D0/SQMW*BE2*(3D0-BE2)**2*CTH
            ASGRE=XW*SGZANG
            ASGIM=0D0
            ASZRE=XW1*SH/(SH-SQMZ)*SGZANG
            ASZIM=0D0
            TGZANG=1D0/SQMW*(BE2*(4D0-2D0*BE2+BE4)+BE2*(4D0-10D0*BE2+
     &      BE4)*CTH+(2D0-11D0*BE2+10D0*BE4)*CTH2+BE2*CTH3)
            ATGRE=0.5D0*XW*SH/TH*TGZANG
            ATGIM=0D0
            ATZRE=0.5D0*XW1*SH/(TH-SQMZ)*TGZANG
            ATZIM=0D0
            UGZANG=1D0/SQMW*(BE2*(4D0-2D0*BE2+BE4)-BE2*(4D0-10D0*BE2+
     &      BE4)*CTH+(2D0-11D0*BE2+10D0*BE4)*CTH2-BE2*CTH3)
            AUGRE=0.5D0*XW*SH/UH*UGZANG
            AUGIM=0D0
            AUZRE=0.5D0*XW1*SH/(UH-SQMZ)*UGZANG
            AUZIM=0D0
            A4ARE=1D0/SQMW*(1D0+2D0*BE2-6D0*BE2*CTH-CTH2)
            A4AIM=0D0
            A4SRE=2D0/SQMW*(1D0+2D0*BE2-CTH2)
            A4SIM=0D0
            FWW=COMFAC*1D0/(4096D0*PARU(1)**2)*(AEM/XW)**4*
     &      (SH/SQMW)**2*SH2
            IF(MSTP(46).LE.0) THEN
              AWWARE=ASHRE
              AWWAIM=ASHIM
              AWWSRE=0D0
              AWWSIM=0D0
            ELSEIF(MSTP(46).EQ.1) THEN
              AWWARE=ASHRE+ATHRE+ASGRE+ASZRE+ATGRE+ATZRE+A4ARE
              AWWAIM=ASHIM+ATHIM+ASGIM+ASZIM+ATGIM+ATZIM+A4AIM
              AWWSRE=-ATHRE-AUHRE+ATGRE+ATZRE+AUGRE+AUZRE+A4SRE
              AWWSIM=-ATHIM-AUHIM+ATGIM+ATZIM+AUGIM+AUZIM+A4SIM
            ELSE
              AWWARE=ASGRE+ASZRE+ATGRE+ATZRE+A4ARE
              AWWAIM=ASGIM+ASZIM+ATGIM+ATZIM+A4AIM
              AWWSRE=ATGRE+ATZRE+AUGRE+AUZRE+A4SRE
              AWWSIM=ATGIM+ATZIM+AUGIM+AUZIM+A4SIM
            ENDIF
            AWWA2=AWWARE**2+AWWAIM**2
            AWWS2=AWWSRE**2+AWWSIM**2
 
          ELSE
C...Strongly interacting Z_L/W_L model of Dobado, Herrero, Terron
            FWWA=COMFAC*(AEM/(4D0*PARU(1)*XW))**2*(64D0/9D0)*
     &      ABS(A00U+0.5D0*A20U+4.5D0*A11U*DBLE(CTH))**2
            FWWS=COMFAC*(AEM/(4D0*PARU(1)*XW))**2*64D0*ABS(A20U)**2
          ENDIF
 
          DO 330 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 330
            EI=SIGN(1D0,DBLE(I))*KCHG(IABS(I),1)
            DO 320 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 320
              EJ=SIGN(1D0,DBLE(J))*KCHG(IABS(J),1)
              IF(EI*EJ.LT.0D0) THEN
C...W+W-
                IF(MSTP(45).EQ.1) GOTO 320
                IF(MSTP(46).LE.2) FACWW=FWW*AWWA2*WIDS(24,1)
                IF(MSTP(46).GE.3) FACWW=FWWA*WIDS(24,1)
              ELSE
C...W+W+/W-W-
                IF(MSTP(45).EQ.2) GOTO 320
                IF(MSTP(46).LE.2) FACWW=FWW*AWWS2
                IF(MSTP(46).GE.3) FACWW=FWWS
                IF(EI.GT.0D0) FACWW=FACWW*WIDS(24,4)
                IF(EI.LT.0D0) FACWW=FACWW*WIDS(24,5)
              ENDIF
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACWW*VINT(180+I)*VINT(180+J)
              IF(EI*EJ.GT.0D0) SIGH(NCHN)=0.5D0*SIGH(NCHN)
  320       CONTINUE
  330     CONTINUE
  340     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.120) THEN
        IF(ISUB.EQ.102) THEN
C...g + g -> h0 (or H0, or A0)
          CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          FACBW=4D0*COMFAC/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(KFHIGG,1)).GT.PARP(48)*PMAS(KFHIGG,2))
     &    FACBW=0D0
C...PS: Only use fixed-width when using SLHA decay table for this Higgs
          IF (IMSS(22).GE.1.AND.MWID(KFHIGG).EQ.2) THEN
            WDTP13=0D0
            DO 345 IDC=MDCY(KFHIGG,2),MDCY(KFHIGG,2)+MDCY(KFHIGG,3)-1
              IF(KFDP(IDC,1).EQ.21.AND.KFDP(IDC,2).EQ.21.AND.
     &            KFDP(IDC,3).EQ.0) WDTP13=PMAS(KFHIGG,2)*BRAT(IDC)
 345        CONTINUE
            IF(WDTP13.EQ.0D0) CALL PYERRM(26,
     &          '(PYSGHG:) did not find Higgs -> g g channel')  
            HI=SHR*WDTP13/32D0
          ELSE
            HI=SHR*WDTP(13)/32D0 
          ENDIF
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 350
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=HI*FACBW*HF
  350     CONTINUE
 
        ELSEIF(ISUB.EQ.103) THEN
C...gamma + gamma -> h0 (or H0, or A0)
          CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          FACBW=4D0*COMFAC/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(KFHIGG,1)).GT.PARP(48)*PMAS(KFHIGG,2))
     &    FACBW=0D0
C...PS: Only use fixed-width when using SLHA decay table for this Higgs
          IF (IMSS(22).GE.1.AND.MWID(KFHIGG).EQ.2) THEN
            WDTP14=0D0
            DO 355 IDC=MDCY(KFHIGG,2),MDCY(KFHIGG,2)+MDCY(KFHIGG,3)-1
              IF(KFDP(IDC,1).EQ.22.AND.KFDP(IDC,2).EQ.22.AND.
     &            KFDP(IDC,3).EQ.0) WDTP14=PMAS(KFHIGG,2)*BRAT(IDC)
 355        CONTINUE
            IF(WDTP14.EQ.0D0) CALL PYERRM(26,
     &          '(PYSGHG:) did not find Higgs -> gamma gamma channel') 
            HI=SHR*WDTP14*2D0
          ELSE
            HI=SHR*WDTP(14)*2D0
          ENDIF
          IF(KFAC(1,22)*KFAC(2,22).EQ.0) GOTO 360
          NCHN=NCHN+1
          ISIG(NCHN,1)=22
          ISIG(NCHN,2)=22
          ISIG(NCHN,3)=1
          SIGH(NCHN)=HI*FACBW*HF
  360     CONTINUE
 
        ELSEIF(ISUB.EQ.110) THEN
C...f + fbar -> gamma + h0
          THUH=MAX(TH*UH,SH*CKIN(3)**2)
          FACHG=COMFAC*(3D0*AEM**4)/(2D0*PARU(1)**2*XW*SQMW)*SH*THUH
          FACHG=FACHG*WIDS(KFHIGG,2)
C...Calculate loop contributions for intermediate gamma* and Z0
          CIGTOT=DCMPLX(0D0,0D0)
          CIZTOT=DCMPLX(0D0,0D0)
          JMAX=3*MSTP(1)+1
          DO 370 J=1,JMAX
            IF(J.LE.2*MSTP(1)) THEN
              FNC=1D0
              EJ=KCHG(J,1)/3D0
              AJ=SIGN(1D0,EJ+0.1D0)
              VJ=AJ-4D0*EJ*XWV
              BALP=SQM4/(2D0*PMAS(J,1))**2
              BBET=SH/(2D0*PMAS(J,1))**2
            ELSEIF(J.LE.3*MSTP(1)) THEN
              FNC=3D0
              JL=2*(J-2*MSTP(1))-1
              EJ=KCHG(10+JL,1)/3D0
              AJ=SIGN(1D0,EJ+0.1D0)
              VJ=AJ-4D0*EJ*XWV
              BALP=SQM4/(2D0*PMAS(10+JL,1))**2
              BBET=SH/(2D0*PMAS(10+JL,1))**2
            ELSE
              BALP=SQM4/(2D0*PMAS(24,1))**2
              BBET=SH/(2D0*PMAS(24,1))**2
            ENDIF
            BABI=1D0/(BALP-BBET)
            IF(BALP.LT.1D0) THEN
              F0ALP=DCMPLX(DBLE(ASIN(SQRT(BALP))),0D0)
              F1ALP=F0ALP**2
            ELSE
              F0ALP=DCMPLX(DBLE(LOG(SQRT(BALP)+SQRT(BALP-1D0))),
     &        -DBLE(0.5D0*PARU(1)))
              F1ALP=-F0ALP**2
            ENDIF
            F2ALP=DBLE(SQRT(ABS(BALP-1D0)/BALP))*F0ALP
            IF(BBET.LT.1D0) THEN
              F0BET=DCMPLX(DBLE(ASIN(SQRT(BBET))),0D0)
              F1BET=F0BET**2
            ELSE
              F0BET=DCMPLX(DBLE(LOG(SQRT(BBET)+SQRT(BBET-1D0))),
     &        -DBLE(0.5D0*PARU(1)))
              F1BET=-F0BET**2
            ENDIF
            F2BET=DBLE(SQRT(ABS(BBET-1D0)/BBET))*F0BET
            IF(J.LE.3*MSTP(1)) THEN
              FIF=DBLE(0.5D0*BABI)+DBLE(BABI**2)*(DBLE(0.5D0*(1D0-BALP+
     &        BBET))*(F1BET-F1ALP)+DBLE(BBET)*(F2BET-F2ALP))
              CIGTOT=CIGTOT+DBLE(FNC*EJ**2)*FIF
              CIZTOT=CIZTOT+DBLE(FNC*EJ*VJ)*FIF
            ELSE
              TXW=XW/XW1
              CIGTOT=CIGTOT-0.5*(DBLE(BABI*(1.5D0+BALP))+DBLE(BABI**2)*
     &        (DBLE(1.5D0-3D0*BALP+4D0*BBET)*(F1BET-F1ALP)+
     &        DBLE(BBET*(2D0*BALP+3D0))*(F2BET-F2ALP)))
              CIZTOT=CIZTOT-DBLE(0.5D0*BABI*XW1)*(DBLE(5D0-TXW+2D0*BALP*
     &        (1D0-TXW))*(1D0+DBLE(2D0*BABI*BBET)*(F2BET-F2ALP))+
     &        DBLE(BABI*(4D0*BBET*(3D0-TXW)-(2D0*BALP-1D0)*(5D0-TXW)))*
     &        (F1BET-F1ALP))
            ENDIF
  370     CONTINUE
          CIGTOT=CIGTOT/DBLE(SH)
          CIZTOT=CIZTOT*DBLE(XWC)/DCMPLX(DBLE(SH-SQMZ),DBLE(GMMZ))
C...Loop over initial flavours
          DO 380 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 380
            EI=KCHG(IABS(I),1)/3D0
            AI=SIGN(1D0,EI)
            VI=AI-4D0*EI*XWV
            FCOI=1D0
            IF(IABS(I).LE.10) FCOI=FACA/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACHG*FCOI*(ABS(DBLE(EI)*CIGTOT+DBLE(VI)*
     &      CIZTOT)**2+AI**2*ABS(CIZTOT)**2)
  380     CONTINUE
 
        ELSEIF(ISUB.EQ.111) THEN
C...f + fbar -> g + h0 (q + qbar -> g + h0 only)
          IF(MSTP(38).NE.0) THEN
C...Simple case: only do gg <-> h exactly.
          CALL PYWIDT(KFHIGG,SQM4,WDTP,WDTE)
C...PS: Only use fixed-width when using SLHA decay table for this Higgs
          IF (IMSS(22).GE.1.AND.MWID(KFHIGG).EQ.2) THEN
            WDTP13=0D0
            DO 385 IDC=MDCY(KFHIGG,2),MDCY(KFHIGG,2)+MDCY(KFHIGG,3)-1
              IF(KFDP(IDC,1).EQ.21.AND.KFDP(IDC,2).EQ.21.AND.
     &            KFDP(IDC,3).EQ.0) WDTP13=PMAS(KFHIGG,2)*BRAT(IDC)
 385        CONTINUE
            IF(WDTP13.EQ.0D0) CALL PYERRM(26,
     &          '(PYSGHG:) did not find Higgs -> g g channel')  
            FACGH=COMFAC*FACA*(2D0/9D0)*AS*(WDTP13/SQRT(SQM4))*
     &          (TH**2+UH**2)/(SH*SQM4)
          ELSE
            FACGH=COMFAC*FACA*(2D0/9D0)*AS*(WDTP(13)/SQRT(SQM4))*
     &          (TH**2+UH**2)/(SH*SQM4)
          ENDIF
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=GMMH/((SQM4-SQMH)**2+GMMH**2)
          GMMHC=SQRT(SQM4)*WDTP(0)
          HBW4C=SQRT(SQM4)*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/
     &    ((SQM4-SQMH)**2+GMMHC**2)
          FACGH=FACGH*HBW4C/HBW4
          ELSE
C...Messy case: do full loop integrals
          A5STUR=0D0
          A5STUI=0D0
          DO 390 I=1,2*MSTP(1)
            SQMQ=PMAS(I,1)**2
            EPSS=4D0*SQMQ/SH
            EPSH=4D0*SQMQ/SQMH
            CALL PYWAUX(1,EPSS,W1SR,W1SI)
            CALL PYWAUX(1,EPSH,W1HR,W1HI)
            CALL PYWAUX(2,EPSS,W2SR,W2SI)
            CALL PYWAUX(2,EPSH,W2HR,W2HI)
            A5STUR=A5STUR+EPSH*(1D0+SH/(TH+UH)*(W1SR-W1HR)+
     &      (0.25D0-SQMQ/(TH+UH))*(W2SR-W2HR))
            A5STUI=A5STUI+EPSH*(SH/(TH+UH)*(W1SI-W1HI)+
     &      (0.25D0-SQMQ/(TH+UH))*(W2SI-W2HI))
  390     CONTINUE
          FACGH=COMFAC*FACA/(144D0*PARU(1)**2)*AEM/XW*AS**3*SQMH/SQMW*
     &    SQMH/SH*(UH**2+TH**2)/(UH+TH)**2*(A5STUR**2+A5STUI**2)
          FACGH=FACGH*WIDS(25,2)
          ENDIF
          DO 400 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 400
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACGH
  400     CONTINUE
 
        ELSEIF(ISUB.EQ.112) THEN
C...f + g -> f + h0 (q + g -> q + h0 only)
          IF(MSTP(38).NE.0) THEN
C...Simple case: only do gg <-> h exactly.
          CALL PYWIDT(KFHIGG,SQM4,WDTP,WDTE)
C...PS: Only use fixed-width when using SLHA decay table for this Higgs
          IF (IMSS(22).GE.1.AND.MWID(KFHIGG).EQ.2) THEN
            WDTP13=0D0
            DO 405 IDC=MDCY(KFHIGG,2),MDCY(KFHIGG,2)+MDCY(KFHIGG,3)-1
              IF(KFDP(IDC,1).EQ.21.AND.KFDP(IDC,2).EQ.21.AND.
     &            KFDP(IDC,3).EQ.0) WDTP13=PMAS(KFHIGG,2)*BRAT(IDC)
 405        CONTINUE
            IF(WDTP13.EQ.0D0) CALL PYERRM(26,
     &          '(PYSGHG:) did not find Higgs -> g g channel')  
            FACQH=COMFAC*FACA*(1D0/12D0)*AS*(WDTP13/SQRT(SQM4))*
     &          (SH**2+UH**2)/(-TH*SQM4)
          ELSE
            FACQH=COMFAC*FACA*(1D0/12D0)*AS*(WDTP(13)/SQRT(SQM4))*
     &          (SH**2+UH**2)/(-TH*SQM4)
          ENDIF
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=GMMH/((SQM4-SQMH)**2+GMMH**2)
          GMMHC=SQRT(SQM4)*WDTP(0)
          HBW4C=SQRT(SQM4)*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/
     &    ((SQM4-SQMH)**2+GMMHC**2)
          FACQH=FACQH*HBW4C/HBW4
          ELSE
C...Messy case: do full loop integrals
          A5TSUR=0D0
          A5TSUI=0D0
          DO 410 I=1,2*MSTP(1)
            SQMQ=PMAS(I,1)**2
            EPST=4D0*SQMQ/TH
            EPSH=4D0*SQMQ/SQMH
            CALL PYWAUX(1,EPST,W1TR,W1TI)
            CALL PYWAUX(1,EPSH,W1HR,W1HI)
            CALL PYWAUX(2,EPST,W2TR,W2TI)
            CALL PYWAUX(2,EPSH,W2HR,W2HI)
            A5TSUR=A5TSUR+EPSH*(1D0+TH/(SH+UH)*(W1TR-W1HR)+
     &      (0.25D0-SQMQ/(SH+UH))*(W2TR-W2HR))
            A5TSUI=A5TSUI+EPSH*(TH/(SH+UH)*(W1TI-W1HI)+
     &      (0.25D0-SQMQ/(SH+UH))*(W2TI-W2HI))
  410     CONTINUE
          FACQH=COMFAC*FACA/(384D0*PARU(1)**2)*AEM/XW*AS**3*SQMH/SQMW*
     &    SQMH/(-TH)*(UH**2+SH**2)/(UH+SH)**2*(A5TSUR**2+A5TSUI**2)
          FACQH=FACQH*WIDS(25,2)
          ENDIF
          DO 430 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 430
            DO 420 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 420
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 420
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQH
  420       CONTINUE
  430     CONTINUE
 
        ELSEIF(ISUB.EQ.113) THEN
C...g + g -> g + h0
          IF(MSTP(38).NE.0) THEN
C...Simple case: only do gg <-> h exactly.
          CALL PYWIDT(KFHIGG,SQM4,WDTP,WDTE)
C...PS: Only use fixed-width when using SLHA decay table for this Higgs
          IF (IMSS(22).GE.1.AND.MWID(KFHIGG).EQ.2) THEN
            WDTP13=0D0
            DO 435 IDC=MDCY(KFHIGG,2),MDCY(KFHIGG,2)+MDCY(KFHIGG,3)-1
              IF(KFDP(IDC,1).EQ.21.AND.KFDP(IDC,2).EQ.21.AND.
     &            KFDP(IDC,3).EQ.0) WDTP13=PMAS(KFHIGG,2)*BRAT(IDC)
 435        CONTINUE
            IF(WDTP13.EQ.0D0) CALL PYERRM(26,
     &          '(PYSGHG:) did not find Higgs -> g g channel')  
            FACGH=COMFAC*FACA*(3D0/16D0)*AS*(WDTP13/SQRT(SQM4))*
     &          (SH**4+TH**4+UH**4+SQM4**4)/(SH*TH*UH*SQM4)
          ELSE
            FACGH=COMFAC*FACA*(3D0/16D0)*AS*(WDTP(13)/SQRT(SQM4))*
     &          (SH**4+TH**4+UH**4+SQM4**4)/(SH*TH*UH*SQM4)
          ENDIF
C...Propagators: as simulated in PYOFSH and as desired
          HBW4=GMMH/((SQM4-SQMH)**2+GMMH**2)
          GMMHC=SQRT(SQM4)*WDTP(0)
          HBW4C=SQRT(SQM4)*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))/
     &    ((SQM4-SQMH)**2+GMMHC**2)
          FACGH=FACGH*HBW4C/HBW4
          ELSE
C...Messy case: do full loop integrals
          A2STUR=0D0
          A2STUI=0D0
          A2USTR=0D0
          A2USTI=0D0
          A2TUSR=0D0
          A2TUSI=0D0
          A4STUR=0D0
          A4STUI=0D0
          DO 440 I=1,2*MSTP(1)
            SQMQ=PMAS(I,1)**2
            EPSS=4D0*SQMQ/SH
            EPST=4D0*SQMQ/TH
            EPSU=4D0*SQMQ/UH
            EPSH=4D0*SQMQ/SQMH
            IF(EPSH.LT.1D-6) GOTO 440
            CALL PYWAUX(1,EPSS,W1SR,W1SI)
            CALL PYWAUX(1,EPST,W1TR,W1TI)
            CALL PYWAUX(1,EPSU,W1UR,W1UI)
            CALL PYWAUX(1,EPSH,W1HR,W1HI)
            CALL PYWAUX(2,EPSS,W2SR,W2SI)
            CALL PYWAUX(2,EPST,W2TR,W2TI)
            CALL PYWAUX(2,EPSU,W2UR,W2UI)
            CALL PYWAUX(2,EPSH,W2HR,W2HI)
            CALL PYI3AU(EPSS,TH/UH,Y3STUR,Y3STUI)
            CALL PYI3AU(EPSS,UH/TH,Y3SUTR,Y3SUTI)
            CALL PYI3AU(EPST,SH/UH,Y3TSUR,Y3TSUI)
            CALL PYI3AU(EPST,UH/SH,Y3TUSR,Y3TUSI)
            CALL PYI3AU(EPSU,SH/TH,Y3USTR,Y3USTI)
            CALL PYI3AU(EPSU,TH/SH,Y3UTSR,Y3UTSI)
            CALL PYI3AU(EPSH,SQMH/SH*TH/UH,YHSTUR,YHSTUI)
            CALL PYI3AU(EPSH,SQMH/SH*UH/TH,YHSUTR,YHSUTI)
            CALL PYI3AU(EPSH,SQMH/TH*SH/UH,YHTSUR,YHTSUI)
            CALL PYI3AU(EPSH,SQMH/TH*UH/SH,YHTUSR,YHTUSI)
            CALL PYI3AU(EPSH,SQMH/UH*SH/TH,YHUSTR,YHUSTI)
            CALL PYI3AU(EPSH,SQMH/UH*TH/SH,YHUTSR,YHUTSI)
            W3STUR=YHSTUR-Y3STUR-Y3UTSR
            W3STUI=YHSTUI-Y3STUI-Y3UTSI
            W3SUTR=YHSUTR-Y3SUTR-Y3TUSR
            W3SUTI=YHSUTI-Y3SUTI-Y3TUSI
            W3TSUR=YHTSUR-Y3TSUR-Y3USTR
            W3TSUI=YHTSUI-Y3TSUI-Y3USTI
            W3TUSR=YHTUSR-Y3TUSR-Y3SUTR
            W3TUSI=YHTUSI-Y3TUSI-Y3SUTI
            W3USTR=YHUSTR-Y3USTR-Y3TSUR
            W3USTI=YHUSTI-Y3USTI-Y3TSUI
            W3UTSR=YHUTSR-Y3UTSR-Y3STUR
            W3UTSI=YHUTSI-Y3UTSI-Y3STUI
            B2STUR=SQMQ/SQMH**2*(SH*(UH-SH)/(SH+UH)+2D0*TH*UH*
     &      (UH+2D0*SH)/(SH+UH)**2*(W1TR-W1HR)+(SQMQ-SH/4D0)*
     &      (0.5D0*W2SR+0.5D0*W2HR-W2TR+W3STUR)+SH2*(2D0*SQMQ/
     &      (SH+UH)**2-0.5D0/(SH+UH))*(W2TR-W2HR)+0.5D0*TH*UH/SH*
     &      (W2HR-2D0*W2TR)+0.125D0*(SH-12D0*SQMQ-4D0*TH*UH/SH)*W3TSUR)
            B2STUI=SQMQ/SQMH**2*(2D0*TH*UH*(UH+2D0*SH)/(SH+UH)**2*
     &      (W1TI-W1HI)+(SQMQ-SH/4D0)*(0.5D0*W2SI+0.5D0*W2HI-W2TI+
     &      W3STUI)+SH2*(2D0*SQMQ/(SH+UH)**2-0.5D0/(SH+UH))*
     &      (W2TI-W2HI)+0.5D0*TH*UH/SH*(W2HI-2D0*W2TI)+0.125D0*
     &      (SH-12D0*SQMQ-4D0*TH*UH/SH)*W3TSUI)
            B2SUTR=SQMQ/SQMH**2*(SH*(TH-SH)/(SH+TH)+2D0*UH*TH*
     &      (TH+2D0*SH)/(SH+TH)**2*(W1UR-W1HR)+(SQMQ-SH/4D0)*
     &      (0.5D0*W2SR+0.5D0*W2HR-W2UR+W3SUTR)+SH2*(2D0*SQMQ/
     &      (SH+TH)**2-0.5D0/(SH+TH))*(W2UR-W2HR)+0.5D0*UH*TH/SH*
     &      (W2HR-2D0*W2UR)+0.125D0*(SH-12D0*SQMQ-4D0*UH*TH/SH)*W3USTR)
            B2SUTI=SQMQ/SQMH**2*(2D0*UH*TH*(TH+2D0*SH)/(SH+TH)**2*
     &      (W1UI-W1HI)+(SQMQ-SH/4D0)*(0.5D0*W2SI+0.5D0*W2HI-W2UI+
     &      W3SUTI)+SH2*(2D0*SQMQ/(SH+TH)**2-0.5D0/(SH+TH))*
     &      (W2UI-W2HI)+0.5D0*UH*TH/SH*(W2HI-2D0*W2UI)+0.125D0*
     &      (SH-12D0*SQMQ-4D0*UH*TH/SH)*W3USTI)
            B2TSUR=SQMQ/SQMH**2*(TH*(UH-TH)/(TH+UH)+2D0*SH*UH*
     &      (UH+2D0*TH)/(TH+UH)**2*(W1SR-W1HR)+(SQMQ-TH/4D0)*
     &      (0.5D0*W2TR+0.5D0*W2HR-W2SR+W3TSUR)+TH2*(2D0*SQMQ/
     &      (TH+UH)**2-0.5D0/(TH+UH))*(W2SR-W2HR)+0.5D0*SH*UH/TH*
     &      (W2HR-2D0*W2SR)+0.125D0*(TH-12D0*SQMQ-4D0*SH*UH/TH)*W3STUR)
            B2TSUI=SQMQ/SQMH**2*(2D0*SH*UH*(UH+2D0*TH)/(TH+UH)**2*
     &      (W1SI-W1HI)+(SQMQ-TH/4D0)*(0.5D0*W2TI+0.5D0*W2HI-W2SI+
     &      W3TSUI)+TH2*(2D0*SQMQ/(TH+UH)**2-0.5D0/(TH+UH))*
     &      (W2SI-W2HI)+0.5D0*SH*UH/TH*(W2HI-2D0*W2SI)+0.125D0*
     &      (TH-12D0*SQMQ-4D0*SH*UH/TH)*W3STUI)
            B2TUSR=SQMQ/SQMH**2*(TH*(SH-TH)/(TH+SH)+2D0*UH*SH*
     &      (SH+2D0*TH)/(TH+SH)**2*(W1UR-W1HR)+(SQMQ-TH/4D0)*
     &      (0.5D0*W2TR+0.5D0*W2HR-W2UR+W3TUSR)+TH2*(2D0*SQMQ/
     &      (TH+SH)**2-0.5D0/(TH+SH))*(W2UR-W2HR)+0.5D0*UH*SH/TH*
     &      (W2HR-2D0*W2UR)+0.125D0*(TH-12D0*SQMQ-4D0*UH*SH/TH)*W3UTSR)
            B2TUSI=SQMQ/SQMH**2*(2D0*UH*SH*(SH+2D0*TH)/(TH+SH)**2*
     &      (W1UI-W1HI)+(SQMQ-TH/4D0)*(0.5D0*W2TI+0.5D0*W2HI-W2UI+
     &      W3TUSI)+TH2*(2D0*SQMQ/(TH+SH)**2-0.5D0/(TH+SH))*
     &      (W2UI-W2HI)+0.5D0*UH*SH/TH*(W2HI-2D0*W2UI)+0.125D0*
     &      (TH-12D0*SQMQ-4D0*UH*SH/TH)*W3UTSI)
            B2USTR=SQMQ/SQMH**2*(UH*(TH-UH)/(UH+TH)+2D0*SH*TH*
     &      (TH+2D0*UH)/(UH+TH)**2*(W1SR-W1HR)+(SQMQ-UH/4D0)*
     &      (0.5D0*W2UR+0.5D0*W2HR-W2SR+W3USTR)+UH2*(2D0*SQMQ/
     &      (UH+TH)**2-0.5D0/(UH+TH))*(W2SR-W2HR)+0.5D0*SH*TH/UH*
     &      (W2HR-2D0*W2SR)+0.125D0*(UH-12D0*SQMQ-4D0*SH*TH/UH)*W3SUTR)
            B2USTI=SQMQ/SQMH**2*(2D0*SH*TH*(TH+2D0*UH)/(UH+TH)**2*
     &      (W1SI-W1HI)+(SQMQ-UH/4D0)*(0.5D0*W2UI+0.5D0*W2HI-W2SI+
     &      W3USTI)+UH2*(2D0*SQMQ/(UH+TH)**2-0.5D0/(UH+TH))*
     &      (W2SI-W2HI)+0.5D0*SH*TH/UH*(W2HI-2D0*W2SI)+0.125D0*
     &      (UH-12D0*SQMQ-4D0*SH*TH/UH)*W3SUTI)
            B2UTSR=SQMQ/SQMH**2*(UH*(SH-UH)/(UH+SH)+2D0*TH*SH*
     &      (SH+2D0*UH)/(UH+SH)**2*(W1TR-W1HR)+(SQMQ-UH/4D0)*
     &      (0.5D0*W2UR+0.5D0*W2HR-W2TR+W3UTSR)+UH2*(2D0*SQMQ/
     &      (UH+SH)**2-0.5D0/(UH+SH))*(W2TR-W2HR)+0.5D0*TH*SH/UH*
     &      (W2HR-2D0*W2TR)+0.125D0*(UH-12D0*SQMQ-4D0*TH*SH/UH)*W3TUSR)
            B2UTSI=SQMQ/SQMH**2*(2D0*TH*SH*(SH+2D0*UH)/(UH+SH)**2*
     &      (W1TI-W1HI)+(SQMQ-UH/4D0)*(0.5D0*W2UI+0.5D0*W2HI-W2TI+
     &      W3UTSI)+UH2*(2D0*SQMQ/(UH+SH)**2-0.5D0/(UH+SH))*
     &      (W2TI-W2HI)+0.5D0*TH*SH/UH*(W2HI-2D0*W2TI)+0.125D0*
     &      (UH-12D0*SQMQ-4D0*TH*SH/UH)*W3TUSI)
            B4STUR=0.25D0*EPSH*(-2D0/3D0+0.25D0*(EPSH-1D0)*
     &      (W2SR-W2HR+W3STUR))
            B4STUI=0.25D0*EPSH*0.25D0*(EPSH-1D0)*(W2SI-W2HI+W3STUI)
            B4TUSR=0.25D0*EPSH*(-2D0/3D0+0.25D0*(EPSH-1D0)*
     &      (W2TR-W2HR+W3TUSR))
            B4TUSI=0.25D0*EPSH*0.25D0*(EPSH-1D0)*(W2TI-W2HI+W3TUSI)
            B4USTR=0.25D0*EPSH*(-2D0/3D0+0.25D0*(EPSH-1D0)*
     &      (W2UR-W2HR+W3USTR))
            B4USTI=0.25D0*EPSH*0.25D0*(EPSH-1D0)*(W2UI-W2HI+W3USTI)
            A2STUR=A2STUR+B2STUR+B2SUTR
            A2STUI=A2STUI+B2STUI+B2SUTI
            A2USTR=A2USTR+B2USTR+B2UTSR
            A2USTI=A2USTI+B2USTI+B2UTSI
            A2TUSR=A2TUSR+B2TUSR+B2TSUR
            A2TUSI=A2TUSI+B2TUSI+B2TSUI
            A4STUR=A4STUR+B4STUR+B4USTR+B4TUSR
            A4STUI=A4STUI+B4STUI+B4USTI+B4TUSI
  440     CONTINUE
          FACGH=COMFAC*FACA*3D0/(128D0*PARU(1)**2)*AEM/XW*AS**3*
     &    SQMH/SQMW*SQMH**3/(SH*TH*UH)*(A2STUR**2+A2STUI**2+A2USTR**2+
     &    A2USTI**2+A2TUSR**2+A2TUSI**2+A4STUR**2+A4STUI**2)
          FACGH=FACGH*WIDS(25,2)
          ENDIF
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 450
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACGH
  450     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.170) THEN
        IF(ISUB.EQ.121) THEN
C...g + g -> Q + Qbar + h0
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 460
          IA=KFPR(ISUBSV,2)
          PMF=PYMRUN(IA,SH)
          FACQQH=COMFAC*(4D0*PARU(1)*AEM/XW)*(4D0*PARU(1)*AS)**2*
     &    (0.5D0*PMF/PMAS(24,1))**2
          WID2=1D0
          IF(IA.EQ.6.OR.IA.EQ.7.OR.IA.EQ.8) WID2=WIDS(IA,1)
          FACQQH=FACQQH*WID2
          IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
            IKFI=1
            IF(IA.LE.10.AND.MOD(IA,2).EQ.0) IKFI=2
            IF(IA.GT.10) IKFI=3
            FACQQH=FACQQH*PARU(150+10*IHIGG+IKFI)**2
            IF(IMSS(1).NE.0.AND.IA.EQ.5) THEN
              FACQQH=FACQQH/(1D0+RMSS(41))**2
              IF(IHIGG.NE.3) THEN
                FACQQH=FACQQH*(1D0+RMSS(41)*PARU(152+10*IHIGG)/
     &          PARU(151+10*IHIGG))**2
              ENDIF
            ENDIF
          ENDIF
          CALL PYQQBH(WTQQBH)
          CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          FACBW=(1D0/PARU(1))*VINT(2)*HF/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(KFHIGG,1)).GT.PARP(48)*PMAS(KFHIGG,2))
     &    FACBW=0D0
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQH*WTQQBH*FACBW
  460     CONTINUE
 
        ELSEIF(ISUB.EQ.122) THEN
C...q + qbar -> Q + Qbar + h0
          IA=KFPR(ISUBSV,2)
          PMF=PYMRUN(IA,SH)
          FACQQH=COMFAC*(4D0*PARU(1)*AEM/XW)*(4D0*PARU(1)*AS)**2*
     &    (0.5D0*PMF/PMAS(24,1))**2
          WID2=1D0
          IF(IA.EQ.6.OR.IA.EQ.7.OR.IA.EQ.8) WID2=WIDS(IA,1)
          FACQQH=FACQQH*WID2
          IF(MSTP(4).GE.1.OR.IHIGG.GE.2) THEN
            IKFI=1
            IF(IA.LE.10.AND.MOD(IA,2).EQ.0) IKFI=2
            IF(IA.GT.10) IKFI=3
            FACQQH=FACQQH*PARU(150+10*IHIGG+IKFI)**2
            IF(IMSS(1).NE.0.AND.IA.EQ.5) THEN
              FACQQH=FACQQH/(1D0+RMSS(41))**2
              IF(IHIGG.NE.3) THEN
                FACQQH=FACQQH*(1D0+RMSS(41)*PARU(152+10*IHIGG)/
     &          PARU(151+10*IHIGG))**2
              ENDIF
            ENDIF
          ENDIF
          CALL PYQQBH(WTQQBH)
          CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          FACBW=(1D0/PARU(1))*VINT(2)*HF/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(KFHIGG,1)).GT.PARP(48)*PMAS(KFHIGG,2))
     &    FACBW=0D0
          DO 470 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 470
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQH*WTQQBH*FACBW
  470     CONTINUE
 
        ELSEIF(ISUB.EQ.123) THEN
C...f + f' -> f + f' + h0 (or H0, or A0) (Z0 + Z0 -> h0 as
C...inner process)
          FACNOR=COMFAC*(4D0*PARU(1)*AEM/(XW*XW1))**3*SQMZ/32D0
          IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACNOR=FACNOR*
     &    PARU(154+10*IHIGG)**2
          FACPRP=1D0/((VINT(215)-VINT(204)**2)*
     &    (VINT(216)-VINT(209)**2))**2
          FACZZ1=FACNOR*FACPRP*(0.5D0*TAUP*VINT(2))*VINT(219)
          FACZZ2=FACNOR*FACPRP*VINT(217)*VINT(218)
          CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          FACBW=(1D0/PARU(1))*VINT(2)*HF/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(KFHIGG,1)).GT.PARP(48)*PMAS(KFHIGG,2))
     &    FACBW=0D0
          DO 490 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 490
            IA=IABS(I)
            DO 480 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 480
              JA=IABS(J)
              EI=KCHG(IA,1)*ISIGN(1,I)/3D0
              AI=SIGN(1D0,KCHG(IA,1)+0.5D0)*ISIGN(1,I)
              VI=AI-4D0*EI*XWV
              EJ=KCHG(JA,1)*ISIGN(1,J)/3D0
              AJ=SIGN(1D0,KCHG(JA,1)+0.5D0)*ISIGN(1,J)
              VJ=AJ-4D0*EJ*XWV
              FACLR1=(VI**2+AI**2)*(VJ**2+AJ**2)+4D0*VI*AI*VJ*AJ
              FACLR2=(VI**2+AI**2)*(VJ**2+AJ**2)-4D0*VI*AI*VJ*AJ
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=(FACLR1*FACZZ1+FACLR2*FACZZ2)*FACBW
  480       CONTINUE
  490     CONTINUE
 
        ELSEIF(ISUB.EQ.124) THEN
C...f + f' -> f" + f"' + h0 (or H0, or A0) (W+ + W- -> h0 as
C...inner process)
          FACNOR=COMFAC*(4D0*PARU(1)*AEM/XW)**3*SQMW
          IF(MSTP(4).GE.1.OR.IHIGG.GE.2) FACNOR=FACNOR*
     &    PARU(155+10*IHIGG)**2
          FACPRP=1D0/((VINT(215)-VINT(204)**2)*
     &    (VINT(216)-VINT(209)**2))**2
          FACWW=FACNOR*FACPRP*(0.5D0*TAUP*VINT(2))*VINT(219)
          CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          HF=SHR*(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))
          FACBW=(1D0/PARU(1))*VINT(2)*HF/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(KFHIGG,1)).GT.PARP(48)*PMAS(KFHIGG,2))
     &    FACBW=0D0
          DO 510 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 510
            EI=SIGN(1D0,DBLE(I))*KCHG(IABS(I),1)
            DO 500 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 500
              EJ=SIGN(1D0,DBLE(J))*KCHG(IABS(J),1)
              IF(EI*EJ.GT.0D0) GOTO 500
              FACLR=VINT(180+I)*VINT(180+J)
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACLR*FACWW*FACBW
  500       CONTINUE
  510     CONTINUE
 
        ELSEIF(ISUB.EQ.143) THEN
C...f + fbar' -> H+/-
          SQMHC=PMAS(37,1)**2
          CALL PYWIDT(37,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=4D0*COMFAC/((SH-SQMHC)**2+HS**2)
          HP=AEM/(8D0*XW)*SH/SQMW*SH
          DO 530 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.KFAC(1,I).EQ.0) GOTO 530
            IA=IABS(I)
            IM=(MOD(IA,10)+1)/2
            DO 520 J=MMIN2,MMAX2
              IF(J.EQ.0.OR.KFAC(2,J).EQ.0) GOTO 520
              JA=IABS(J)
              JM=(MOD(JA,10)+1)/2
              IF(I*J.GT.0.OR.IA.EQ.JA.OR.IM.NE.JM) GOTO 520
              IF((IA.LE.10.AND.JA.GT.10).OR.(IA.GT.10.AND.JA.LE.10))
     &        GOTO 520
              IF(MOD(IA,2).EQ.0) THEN
                IU=IA
                IL=JA
              ELSE
                IU=JA
                IL=IA
              ENDIF
              RML=PYMRUN(IL,SH)**2/SH
              RMU=PYMRUN(IU,SH)**2/SH
              HI=HP*(RML*PARU(141)**2+RMU/PARU(141)**2)
              IF(IA.LE.10) HI=HI*FACA/3D0
              KCHHC=(KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J))/3
              HF=SHR*(WDTE(0,1)+WDTE(0,(5-KCHHC)/2)+WDTE(0,4))
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=HI*FACBW*HF
  520       CONTINUE
  530     CONTINUE
 
        ELSEIF(ISUB.EQ.161) THEN
C...f + g -> f' + H+/- (b + g -> t + H+/- only)
C...(choice of only b and t to avoid kinematics problems)
          FHCQ=COMFAC*FACA*AS*AEM/XW*1D0/24
C...H propagator: as simulated in PYOFSH and as desired
          SQMHC=PMAS(37,1)**2
          GMMHC=PMAS(37,1)*PMAS(37,2)
          HBW4=GMMHC/((SQM4-SQMHC)**2+GMMHC**2)
          CALL PYWIDT(37,SQM4,WDTP,WDTE)
          GMMHCC=SQRT(SQM4)*WDTP(0)
          HBW4C=GMMHCC/((SQM4-SQMHC)**2+GMMHCC**2)
          FHCQ=FHCQ*HBW4C/HBW4
          Q2RM=SH
          IF(MSTP(32).EQ.12) Q2RM=PARP(194)
          DO 550 I=MMINA,MMAXA
            IA=IABS(I)
            IF(IA.NE.5) GOTO 550
            SQML=PYMRUN(IA,Q2RM)**2
            IUA=IA+MOD(IA,2)
            SQMQ=PYMRUN(IUA,Q2RM)**2
            FACHCQ=FHCQ*(SQML*PARU(141)**2+SQMQ/PARU(141)**2)/SQMW*
     &      (SH/(SQMQ-UH)+2D0*SQMQ*(SQMHC-UH)/(SQMQ-UH)**2+(SQMQ-UH)/SH-
     &      2D0*SQMQ/(SQMQ-UH)+2D0*(SQMHC-UH)/(SQMQ-UH)*
     &      (SQMHC-SQMQ-SH)/SH)
            KCHHC=ISIGN(1,KCHG(IA,1)*ISIGN(1,I))
            DO 540 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 540
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 540
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACHCQ*WIDS(37,(5-KCHHC)/2)
              IF(IUA.EQ.6) SIGH(NCHN)=SIGH(NCHN)*WIDS(6,(5+KCHHC)/2)
  540       CONTINUE
  550     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.402) THEN
        IF(ISUB.EQ.401) THEN
C...  g + g -> t + bbar + H-
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 560
          IA=KFPR(ISUBSV,2)
          CALL PYSTBH(WTTBH)
          CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=(1D0/PARU(1))*VINT(2)*HS/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(KFHIGG,1)).GT.PARP(48)*PMAS(KFHIGG,2))
     &       FACBW=0D0
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=2d0*COMFAC*WTTBH*FACBW
c     Since we don't know yet if H+ or H-, assume H+
c     when calculating suppression due to closed channels.
          SIGH(NCHN)=SIGH(NCHN)*WIDS(37,2)*WIDS(6,3)
          IF(ABS(WIDS(37,2)-WIDS(37,3))
     &       .GE.1D-6*(WIDS(37,2)+WIDS(37,3)).OR.
     &       ABS(WIDS(6,2)-WIDS(6,3))
     &       .GE.1D-6*(WIDS(6,2)+WIDS(6,3))) THEN
            WRITE(*,*)'Error: Process 401 cannot handle different'
            WRITE(*,*)'decays for H+ and H- or t and tbar.'
            WRITE(*,*)'Execution stopped.'
            CALL PYSTOP(108)
          END IF
 560      CONTINUE
 
        ELSEIF(ISUB.EQ.402) THEN
C...  q + qbar -> t + bbar + H-
          IA=KFPR(ISUBSV,2)
          CALL PYSTBH(WTTBH)
          CALL PYWIDT(KFHIGG,SH,WDTP,WDTE)
          HS=SHR*WDTP(0)
          FACBW=(1D0/PARU(1))*VINT(2)*HS/((SH-SQMH)**2+HS**2)
          IF(ABS(SHR-PMAS(KFHIGG,1)).GT.PARP(48)*PMAS(KFHIGG,2))
     &       FACBW=0D0
          DO 570 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &         KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 570
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=2d0*COMFAC*WTTBH*FACBW
c     Since we don't know yet if H+ or H-, assume H+
c     when calculating suppression due to closed channels.
            SIGH(NCHN)=SIGH(NCHN)*WIDS(37,2)*WIDS(6,3)
            IF(ABS(WIDS(37,2)-WIDS(37,3))/(WIDS(37,2)+WIDS(37,3))
     &         .GE.1D-6.OR.
     &         ABS(WIDS(6,2)-WIDS(6,3))/(WIDS(6,2)+WIDS(6,3))
     &         .GE.1D-6) THEN
              WRITE(*,*)'Error: Process 402 cannot handle different'
              WRITE(*,*)'decays for H+ and H- or t and tbar.'
              WRITE(*,*)'Execution stopped.'
              CALL PYSTOP(108)
            END IF
 570      CONTINUE
        ENDIF
      ENDIF
 
      RETURN
      END
