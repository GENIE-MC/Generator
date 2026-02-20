 
C*********************************************************************
 
C...PYSGSU
C...Subprocess cross sections for SUSY processes,
C...including Higgs pair production.
C...Auxiliary to PYSIGH.
 
      SUBROUTINE PYSGSU(NCHN,SIGS)
 
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
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
      COMMON/PYSGCM/ISUB,ISUBSV,MMIN1,MMAX1,MMIN2,MMAX2,MMINA,MMAXA,
     &KFAC(2,-40:40),COMFAC,FACK,FACA,SH,TH,UH,SH2,TH2,UH2,SQM3,SQM4,
     &SHR,SQPTH,TAUP,BE34,CTH,X(2),SQMZ,SQMW,GMMZ,GMMW,
     &AEM,AS,XW,XW1,XWC,XWV,POLL,POLR,POLLL,POLRR
      SAVE /PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT3/,
     &/PYINT4/,/PYMSSM/,/PYSSMT/,/PYSGCM/
C...Local arrays and complex variables
      DIMENSION WDTP(0:400),WDTE(0:400,0:5)
      COMPLEX*16 OLPP,ORPP,OLP,ORP,OL,OR,QLL,QLR
      COMPLEX*16 QRR,QRL,GLIJ,GRIJ,PROPW,PROPZ
      COMPLEX*16 ZMIXC(4,4),UMIXC(2,2),VMIXC(2,2)
 
CMRENNA++
C...Z and W width, combinations of weak mixing angle
      ZWID=PMAS(23,2)
      WWID=PMAS(24,2)
      TANW=SQRT(XW/XW1)
      CT2W=(1D0-2D0*XW)/(2D0*XW/TANW)
 
C...Convert almost equivalent SUSY processes into each other
C...Extract differences in flavours and couplings
 
C...Sleptons and sneutrinos
      IF(ISUB.EQ.201.OR.ISUB.EQ.204.OR.ISUB.EQ.207) THEN
        KFID=MOD(KFPR(ISUB,1),KSUSY1)
        ISUB=201
        ILR=0
      ELSEIF(ISUB.EQ.202.OR.ISUB.EQ.205.OR.ISUB.EQ.208) THEN
        KFID=MOD(KFPR(ISUB,1),KSUSY1)
        ISUB=201
        ILR=1
      ELSEIF(ISUB.EQ.203.OR.ISUB.EQ.206.OR.ISUB.EQ.209) THEN
        KFID=MOD(KFPR(ISUB,1),KSUSY1)
        ISUB=203
      ELSEIF(ISUB.GE.210.AND.ISUB.LE.212) THEN
        IF(ISUB.EQ.210) THEN
          RKF=2.0D0
        ELSEIF(ISUB.EQ.211) THEN
          RKF=SFMIX(15,1)**2
        ELSEIF(ISUB.EQ.212) THEN
          RKF=SFMIX(15,2)**2
        ENDIF
          ISUB=210
      ELSEIF(ISUB.EQ.213.OR.ISUB.EQ.214) THEN
        IF(ISUB.EQ.213) THEN
          KFID=MOD(KFPR(ISUB,1),KSUSY1)
          RKF=2.0D0
        ELSEIF(ISUB.EQ.214) THEN
          KFID=16
          RKF=1.0D0
        ENDIF
        ISUB=213
 
C...Neutralinos
      ELSEIF(ISUB.GE.216.AND.ISUB.LE.225) THEN
        IF(ISUB.EQ.216) THEN
          IZID1=1
          IZID2=1
        ELSEIF(ISUB.EQ.217) THEN
          IZID1=2
          IZID2=2
        ELSEIF(ISUB.EQ.218) THEN
          IZID1=3
          IZID2=3
        ELSEIF(ISUB.EQ.219) THEN
          IZID1=4
          IZID2=4
        ELSEIF(ISUB.EQ.220) THEN
          IZID1=1
          IZID2=2
        ELSEIF(ISUB.EQ.221) THEN
          IZID1=1
          IZID2=3
        ELSEIF(ISUB.EQ.222) THEN
          IZID1=1
          IZID2=4
        ELSEIF(ISUB.EQ.223) THEN
          IZID1=2
          IZID2=3
        ELSEIF(ISUB.EQ.224) THEN
          IZID1=2
          IZID2=4
        ELSEIF(ISUB.EQ.225) THEN
          IZID1=3
          IZID2=4
        ENDIF
        ISUB=216
 
C...Charginos
      ELSEIF(ISUB.GE.226.AND.ISUB.LE.228) THEN
        IF(ISUB.EQ.226) THEN
          IZID1=1
          IZID2=1
        ELSEIF(ISUB.EQ.227) THEN
          IZID1=2
          IZID2=2
        ELSEIF(ISUB.EQ.228) THEN
          IZID1=1
          IZID2=2
        ENDIF
        ISUB=226
 
C...Neutralino + chargino
      ELSEIF(ISUB.GE.229.AND.ISUB.LE.236) THEN
        IF(ISUB.EQ.229) THEN
          IZID1=1
          IZID2=1
        ELSEIF(ISUB.EQ.230) THEN
          IZID1=1
          IZID2=2
        ELSEIF(ISUB.EQ.231) THEN
          IZID1=1
          IZID2=3
        ELSEIF(ISUB.EQ.232) THEN
          IZID1=1
          IZID2=4
        ELSEIF(ISUB.EQ.233) THEN
          IZID1=2
          IZID2=1
        ELSEIF(ISUB.EQ.234) THEN
          IZID1=2
          IZID2=2
        ELSEIF(ISUB.EQ.235) THEN
          IZID1=2
          IZID2=3
        ELSEIF(ISUB.EQ.236) THEN
          IZID1=2
          IZID2=4
        ENDIF
        ISUB=229
 
C...Gluino + neutralino
      ELSEIF(ISUB.GE.237.AND.ISUB.LE.240) THEN
        IF(ISUB.EQ.237) THEN
          IZID=1
        ELSEIF(ISUB.EQ.238) THEN
          IZID=2
        ELSEIF(ISUB.EQ.239) THEN
          IZID=3
        ELSEIF(ISUB.EQ.240) THEN
          IZID=4
        ENDIF
        ISUB=237
 
C...Gluino + chargino
      ELSEIF(ISUB.GE.241.AND.ISUB.LE.242) THEN
        IF(ISUB.EQ.241) THEN
          IZID=1
        ELSEIF(ISUB.EQ.242) THEN
          IZID=2
        ENDIF
        ISUB=241
 
C...Squark + neutralino
      ELSEIF(ISUB.GE.246.AND.ISUB.LE.253) THEN
        ILR=0
        IF(MOD(ISUB,2).NE.0) ILR=1
        IF(ISUB.LE.247) THEN
          IZID=1
        ELSEIF(ISUB.LE.249) THEN
          IZID=2
        ELSEIF(ISUB.LE.251) THEN
          IZID=3
        ELSEIF(ISUB.LE.253) THEN
          IZID=4
        ENDIF
        ISUB=246
        RKF=5D0
 
C...Squark + chargino
      ELSEIF(ISUB.GE.254.AND.ISUB.LE.257) THEN
        IF(ISUB.LE.255) THEN
          IZID=1
        ELSEIF(ISUB.LE.257) THEN
          IZID=2
        ENDIF
        IF(MOD(ISUB,2).EQ.0) THEN
          ILR=0
        ELSE
          ILR=1
        ENDIF
        ISUB=254
        RKF=5D0
 
C...Squark + gluino
      ELSEIF(ISUB.EQ.258.OR.ISUB.EQ.259) THEN
        ISUB=258
        RKF=4D0
 
C...Stops
      ELSEIF(ISUB.EQ.261.OR.ISUB.EQ.262) THEN
        ILR=0
        IF(ISUB.EQ.262) ILR=1
        ISUB=261
      ELSEIF(ISUB.EQ.265) THEN
        ISUB=264
 
C...Squarks
      ELSEIF(ISUB.GE.271.AND.ISUB.LE.280) THEN
        ILR=0
        IF(ISUB.LE.273) THEN
          IF(ISUB.EQ.273) ILR=1
          ISUB=271
          RKF=16D0
        ELSEIF(ISUB.LE.276) THEN
          IF(ISUB.EQ.276) ILR=1
          ISUB=274
          RKF=16D0
        ELSEIF(ISUB.LE.278) THEN
          IF(ISUB.EQ.278) ILR=1
          ISUB=277
          RKF=4D0
        ELSE
          IF(ISUB.EQ.280) ILR=1
          ISUB=279
          RKF=4D0
        ENDIF
C...Sbottoms
      ELSEIF(ISUB.GE.281.AND.ISUB.LE.296) THEN
        ILR=0
        IF(ISUB.LE.283) THEN
          IF(ISUB.EQ.283) ILR=1
          ISUB=271
          RKF=4D0
        ELSEIF(ISUB.LE.286) THEN
          IF(ISUB.EQ.286) ILR=1
          ISUB=274
          RKF=4D0
        ELSEIF(ISUB.LE.288) THEN
          IF(ISUB.EQ.288) ILR=1
          ISUB=277
          RKF=1D0
        ELSEIF(ISUB.LE.290) THEN
          IF(ISUB.EQ.290) ILR=1
          ISUB=279
          RKF=1D0
        ELSEIF(ISUB.LE.293) THEN
          IF(ISUB.EQ.293) ILR=1
          ISUB=271
          RKF=1D0
        ELSEIF(ISUB.EQ.296) THEN
          ILR=1
          ISUB=274
          RKF=1D0
C...Squark + gluino
        ELSEIF(ISUB.EQ.294.OR.ISUB.EQ.295) THEN
          ISUB=258
          RKF=1D0
        ENDIF
C...H+/- + H0
      ELSEIF(ISUB.EQ.297.OR.ISUB.EQ.298) THEN
        IF(ISUB.EQ.297) THEN
          RKF=.5D0*PARU(195)**2
        ELSEIF(ISUB.EQ.298) THEN
          RKF=.5D0*(1D0-PARU(195)**2)
        ENDIF
        ISUB=210
C...A0 + H0
      ELSEIF(ISUB.EQ.299.OR.ISUB.EQ.300) THEN
        IF(ISUB.EQ.299) THEN
          RKF=PARU(186)**2
          KFID=25
        ELSEIF(ISUB.EQ.300) THEN
          RKF=PARU(187)**2
          KFID=35
        ENDIF
        ISUB=213
C...H+ + H-
      ELSEIF(ISUB.EQ.301) THEN
        KFID=37
        RKF=1D0
        ISUB=201
      ENDIF
 
C...Supersymmetric processes - all of type 2 -> 2 :
C...correct final-state Breit-Wigners from fixed to running width.
      IF(MSTP(42).GT.0) THEN
        DO 100 I=1,2
        KFLW=KFPR(ISUBSV,I)
        KCW=PYCOMP(KFLW)
        IF(PMAS(KCW,2).LT.PARP(41)) GOTO 100
        IF(I.EQ.1) SQMI=SQM3
        IF(I.EQ.2) SQMI=SQM4
        SQMS=PMAS(KCW,1)**2
        GMMS=PMAS(KCW,1)*PMAS(KCW,2)
        HBWS=GMMS/((SQMI-SQMS)**2+GMMS**2)
        CALL PYWIDT(KFLW,SQMI,WDTP,WDTE)
        GMMI=SQRT(SQMI)*WDTP(0)
        HBWI=GMMI/((SQMI-SQMS)**2+GMMI**2)
        COMFAC=COMFAC*(HBWI/HBWS)
  100   CONTINUE
      ENDIF
 
C...Differential cross section expressions.
 
      IF(ISUB.LE.210) THEN
        IF(ISUB.EQ.201) THEN
C...f + fbar -> e_L + e_Lbar
          COMFAC=COMFAC*WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
          DO 130 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 130
            EI=KCHG(IA,1)/3D0
            TT3I=SIGN(1D0,EI+1D-6)/2D0
            EJ=-1D0
            TT3J=-1D0/2D0
            FCOL=1D0
C...Color factor for e+ e-
            IF(IA.GE.11) FCOL=3D0
            IF(ISUBSV.EQ.301) THEN
              A1=1D0
              A2=0D0
            ELSEIF(ILR.EQ.1) THEN
              A1=SFMIX(KFID,3)**2
              A2=SFMIX(KFID,4)**2
            ELSEIF(ILR.EQ.0) THEN
              A1=SFMIX(KFID,1)**2
              A2=SFMIX(KFID,2)**2
            ENDIF
            XLQ=(TT3J-EJ*XW)*A1
            XRQ=(-EJ*XW)*A2
            XLF=(TT3I-EI*XW)
            XRF=(-EI*XW)
            TAA=(EI*EJ)**2*(POLL+POLR)
            TZZ=(XLF**2*POLL+XRF**2*POLR)*(XLQ+XRQ)**2/XW**2/XW1**2
            TZZ=TZZ/((1D0-SQMZ/SH)**2+SQMZ*ZWID/SH**2)
            TAZ=2D0*EI*EJ*(XLQ+XRQ)*(XLF*POLL+XRF*POLR)/XW/XW1
            TAZ=TAZ/((1D0-SQMZ/SH)**2+SQMZ*(ZWID/SH)**2)*(1D0-SQMZ/SH)
            TNN=0.0D0
            TAN=0.0D0
            TZN=0.0D0
            IF(IA.GE.11.AND.IA.LE.18.AND.KFID.EQ.IA) THEN
              FAC2=SQRT(2D0)
              TNN1=0D0
              TNN2=0D0
              TNN3=0D0
              DO 120 II=1,4
                DK=1D0/(TH-SMZ(II)**2)
                FLEK=-FAC2*(TT3I*ZMIX(II,2)-TANW*(TT3I-EI)*
     &          ZMIX(II,1))
                FREK=FAC2*TANW*EI*ZMIX(II,1)
                TNN1=TNN1+FLEK**2*DK
                TNN2=TNN2+FREK**2*DK
                DO 110 JJ=1,4
                  DL=1D0/(TH-SMZ(JJ)**2)
                  FLEL=-FAC2*(TT3J*ZMIX(JJ,2)-TANW*(TT3J-EJ)*
     &            ZMIX(JJ,1))
                  FREL=FAC2*TANW*EJ*ZMIX(JJ,1)
                  TNN3=TNN3+FLEK*FREK*FLEL*FREL*DK*DL*SMZ(II)*SMZ(JJ)
  110           CONTINUE
  120         CONTINUE
              TNN=(UH*TH-SQM3*SQM4)*(A1**2*TNN1**2*POLL+
     &        A2**2*TNN2**2*POLR)
              TNN=(TNN+SH*A1*A2*TNN3*((1D0-PARJ(131))*(1D0-PARJ(132))+
     &        (1D0+PARJ(131))*(1D0+PARJ(132))))/4D0/XW**2
              TZN=(UH*TH-SQM3*SQM4)*(XLQ+XRQ)*
     &        (TNN1*XLF*A1*POLL+TNN2*XRF*A2*POLR)
              TZN=TZN/((1D0-SQMZ/SH)**2+SQMZ*(ZWID/SH)**2)*
     &        (1D0-SQMZ/SH)/SH
              TZN=TZN/XW**2/XW1
              TAN=EI*EJ*(UH*TH-SQM3*SQM4)/SH*(A1*TNN1*POLL+
     &        A2*TNN2*POLR)/XW
            ENDIF
            FACQQ1=COMFAC*AEM**2*(TAA+TZZ+TAZ)*FCOL/3D0
            FACQQ1=FACQQ1*( UH*TH-SQM3*SQM4 )/SH**2
            FACQQ2=COMFAC*AEM**2*(TNN+TZN+TAN)*FCOL/3D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQ1+FACQQ2
  130     CONTINUE
 
        ELSEIF(ISUB.EQ.203) THEN
C...f + fbar -> e_L + e_Rbar
          DO 160 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 160
            EI=KCHG(IABS(I),1)/3D0
            TT3I=SIGN(1D0,EI)/2D0
            EJ=-1
            TT3J=-1D0/2D0
            FCOL=1D0
C...Color factor for e+ e-
            IF(IA.GE.11) FCOL=3D0
            A1=SFMIX(KFID,1)**2
            A2=SFMIX(KFID,2)**2
            XLQ=(TT3J-EJ*XW)
            XRQ=(-EJ*XW)
            XLF=(TT3I-EI*XW)
            XRF=(-EI*XW)
            TZZ=(XLF**2*POLL+XRF**2*POLR)*(XLQ-XRQ)**2
     &      /XW**2/XW1**2*A1*A2
            TZZ=TZZ/((1D0-SQMZ/SH)**2+SQMZ*(ZWID/SH)**2)
            TNN=0.0D0
            TZN=0.0D0
            TNNA=0D0
            TNNB=0D0
            IF(IA.GE.11.AND.IA.LE.18.AND.KFID.EQ.IA) THEN
              FAC2=SQRT(2D0)
              TNN1=0D0
              TNN2=0D0
              TNN3=0D0
              DO 150 II=1,4
                DK=1D0/(TH-SMZ(II)**2)
                FLEK=-FAC2*(TT3I*ZMIX(II,2)-TANW*(TT3I-EI)*
     &          ZMIX(II,1))
                FREK=FAC2*TANW*EI*ZMIX(II,1)
                TNN1=TNN1+FLEK**2*DK
                TNN2=TNN2+FREK**2*DK
                DO 140 JJ=1,4
                  DL=1D0/(TH-SMZ(JJ)**2)
                  FLEL=-FAC2*(TT3J*ZMIX(JJ,2)-TANW*(TT3J-EJ)*
     &            ZMIX(JJ,1))
                  FREL=FAC2*TANW*EJ*ZMIX(JJ,1)
                  TNN3=TNN3+FLEK*FREK*FLEL*FREL*DK*DL*SMZ(II)*SMZ(JJ)
  140           CONTINUE
  150         CONTINUE
              TNN=(UH*TH-SQM3*SQM4)*A1*A2*(TNN2**2*POLR+TNN1**2*POLL)
              TNNA=(TNN+SH*(A1**2*POLLL+A2**2*POLRR)*TNN3)/4D0
              TNNB=(TNN+SH*(A1**2*POLRR+A2**2*POLLL)*TNN3)/4D0
              TZN=(UH*TH-SQM3*SQM4)*A1*A2
              TZN=TZN*(XLQ-XRQ)*(XLF*TNN1*POLL-XRF*TNN2*POLR)/XW1
              TZN=TZN/((1D0-SQMZ/SH)**2+SQMZ*(ZWID/SH)**2)*
     &        (1D0-SQMZ/SH)/SH
            ENDIF
            FACQQ0=COMFAC*AEM**2*TZZ*FCOL/3D0*(UH*TH-SQM3*SQM4)/SH2
            FACQQ2=COMFAC*AEM**2/XW**2*(TNNA+TZN)*FCOL/3D0
            FACQQ1=COMFAC*AEM**2/XW**2*(TNNB+TZN)*FCOL/3D0
C%%%%%%%%%%%
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=(FACQQ0+FACQQ1)*WIDS(PYCOMP(KFPR(ISUBSV,1)),2)*
     &      WIDS(PYCOMP(KFPR(ISUBSV,2)),3)
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=2
            SIGH(NCHN)=(FACQQ0+FACQQ2)*WIDS(PYCOMP(KFPR(ISUBSV,1)),3)*
     &      WIDS(PYCOMP(KFPR(ISUBSV,2)),2)
  160     CONTINUE
 
        ELSEIF(ISUB.EQ.210) THEN
C...q + qbar' -> W*- > ~l_L + ~nu_L
          FAC0=RKF*COMFAC*AEM**2/XW**2/12D0
          FAC1=(TH*UH-SQM3*SQM4)/((SH-SQMW)**2+WWID**2*SQMW)
          DO 180 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.10.OR.KFAC(1,I).EQ.0) GOTO 180
            DO 170 J=MMIN2,MMAX2
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.10.OR.KFAC(2,J).EQ.0) GOTO 170
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 170
              FCKM=3D0
              IF(IA.LE.10) FCKM=VCKM((IA+1)/2,(JA+1)/2)
              KCHSUM=KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J)
              KCHW=2
              IF(KCHSUM.LT.0) KCHW=3
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              IF(ISUBSV.EQ.297.OR.ISUBSV.EQ.298) THEN
                FACR=WIDS(PYCOMP(KFPR(ISUBSV,1)),5-KCHW)*
     &          WIDS(PYCOMP(KFPR(ISUBSV,2)),2)
              ELSE
                FACR=WIDS(PYCOMP(KFPR(ISUBSV,1)),5-KCHW)*
     &          WIDS(PYCOMP(KFPR(ISUBSV,2)),KCHW)
              ENDIF
              SIGH(NCHN)=FAC0*FAC1*FCKM*FACR
  170       CONTINUE
  180     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.220) THEN
        IF(ISUB.EQ.213) THEN
C...f + fbar -> ~nu_L + ~nu_Lbar
          IF(ISUBSV.EQ.299.OR.ISUBSV.EQ.300) THEN
            FACR=WIDS(PYCOMP(KFPR(ISUBSV,1)),2)*
     &      WIDS(PYCOMP(KFPR(ISUBSV,2)),2)
          ELSE
            FACR=WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
          ENDIF
          COMFAC=COMFAC*FACR
          PROPZ2=(SH-SQMZ)**2+ZWID**2*SQMZ
          XLL=0.5D0
          XLR=0.0D0
          DO 190 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 190
            EI=KCHG(IA,1)/3D0
            FCOL=1D0
C...Color factor for e+ e-
            IF(IA.GE.11) FCOL=3D0
            XLQ=(SIGN(1D0,EI)-2D0*EI*XW)/2D0
            XRQ=-EI*XW
            TZC=0.0D0
            TCC=0.0D0
            IF(IA.GE.11.AND.KFID.EQ.IA+1) THEN
              TZC=VMIX(1,1)**2/(TH-SMW(1)**2)+VMIX(2,1)**2/
     &        (TH-SMW(2)**2)
              TCC=TZC**2
              TZC=TZC/XW1*(SH-SQMZ)/PROPZ2*XLQ*XLL
            ENDIF
            FACQQ1=(XLQ**2+XRQ**2)*(XLL+XLR)**2/XW1**2/PROPZ2
            FACQQ2=TZC+TCC/4D0
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=(FACQQ1+FACQQ2)*RKF*(UH*TH-SQM3*SQM4)*COMFAC
     &      *AEM**2*FCOL/3D0/XW**2
  190     CONTINUE
 
        ELSEIF(ISUB.EQ.216) THEN
C...q + qbar -> ~chi0_1 + ~chi0_1
          IF(IZID1.EQ.IZID2) THEN
            COMFAC=COMFAC*WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
          ELSE
            COMFAC=COMFAC*WIDS(PYCOMP(KFPR(ISUBSV,1)),2)*
     &      WIDS(PYCOMP(KFPR(ISUBSV,2)),2)
          ENDIF
          FACXX=COMFAC*AEM**2/3D0/XW**2
          IF(IZID1.EQ.IZID2) FACXX=FACXX/2D0
          ZM12=SQM3
          ZM22=SQM4
          WU2 = (UH-ZM12)*(UH-ZM22)
          WT2 = (TH-ZM12)*(TH-ZM22)
          WS2 = SMZ(IZID1)*SMZ(IZID2)*SH
          PROPZ2 = (SH-SQMZ)**2 + SQMZ*ZWID**2
          PROPZ=DCMPLX(SH-SQMZ,-ZWID*PMAS(23,1))/DCMPLX(PROPZ2)
          DO 200 I=1,4
            ZMIXC(IZID1,I)=DCMPLX(ZMIX(IZID1,I),ZMIXI(IZID1,I))
            IF(IZID2.NE.IZID1) THEN
              ZMIXC(IZID2,I)=DCMPLX(ZMIX(IZID2,I),ZMIXI(IZID2,I))
            ENDIF
  200     CONTINUE
          OLPP=(ZMIXC(IZID1,3)*DCONJG(ZMIXC(IZID2,3))-
     &    ZMIXC(IZID1,4)*DCONJG(ZMIXC(IZID2,4)))/2D0
          ORPP=DCONJG(OLPP)
          DO 210 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 210
            EI=KCHG(IABS(I),1)/3D0
            T3I=SIGN(1D0,EI+1D-6)/2D0
            XML2=PMAS(PYCOMP(KSUSY1+IABS(I)),1)**2
            XMR2=PMAS(PYCOMP(KSUSY2+IABS(I)),1)**2
            GLIJ=(T3I*ZMIXC(IZID1,2)-TANW*(T3I-EI)*ZMIXC(IZID1,1))*
     &      DCONJG(T3I*ZMIXC(IZID2,2)-TANW*(T3I-EI)*ZMIXC(IZID2,1))
            GRIJ=ZMIXC(IZID1,1)*DCONJG(ZMIXC(IZID2,1))*(EI*TANW)**2
            QLL=DCMPLX((T3I-EI*XW)/XW1)*OLPP*PROPZ-GLIJ/DCMPLX(UH-XML2)
            QLR=-DCMPLX((T3I-EI*XW)/XW1)*ORPP*PROPZ+DCONJG(GLIJ)
     &      /DCMPLX(TH-XML2)
            QRL=-DCMPLX((EI*XW)/XW1)*OLPP*PROPZ+GRIJ/DCMPLX(TH-XMR2)
            QRR=DCMPLX((EI*XW)/XW1)*ORPP*PROPZ
     &      -DCONJG(GRIJ)/DCMPLX(UH-XMR2)
            FCOL=1D0
            IF(IABS(I).GE.11) FCOL=3D0
            FACGG1=(ABS(QLL)**2*POLL+ABS(QRR)**2*POLR)*WU2+
     &      (ABS(QRL)**2*POLR+ABS(QLR)**2*POLL)*WT2+
     &      2D0*DBLE(QLR*DCONJG(QLL)*POLL+
     &      QRL*DCONJG(QRR)*POLR)*WS2
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACXX*FACGG1*FCOL
  210     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.230) THEN
        IF(ISUB.EQ.226) THEN
C...f + fbar -> ~chi+_1 + ~chi-_1
          FACXX=COMFAC*AEM**2/3D0
          ZM12=SQM3
          ZM22=SQM4
          WU2 = (UH-ZM12)*(UH-ZM22)
          WT2 = (TH-ZM12)*(TH-ZM22)
          WS2 = SMW(IZID1)*SMW(IZID2)*SH
          PROPZ2 = (SH-SQMZ)**2 + SQMZ*ZWID**2
          PROPZ=DCMPLX(SH-SQMZ,-ZWID*PMAS(23,1))/DCMPLX(PROPZ2)
          DIFF=0D0
          IF(IZID1.EQ.IZID2) DIFF=1D0
          DO 220 I=1,2
            VMIXC(IZID1,I)=DCMPLX(VMIX(IZID1,I),VMIXI(IZID1,I))
            UMIXC(IZID1,I)=DCMPLX(UMIX(IZID1,I),UMIXI(IZID1,I))
            IF(IZID2.NE.IZID1) THEN
              VMIXC(IZID2,I)=DCMPLX(VMIX(IZID2,I),VMIXI(IZID2,I))
              UMIXC(IZID2,I)=DCMPLX(UMIX(IZID2,I),UMIXI(IZID2,I))
            ENDIF
  220     CONTINUE
          OLP=-VMIXC(IZID2,1)*DCONJG(VMIXC(IZID1,1))-
     &    VMIXC(IZID2,2)*DCONJG(VMIXC(IZID1,2))/2D0+DCMPLX(XW*DIFF)
          ORP=-UMIXC(IZID1,1)*DCONJG(UMIXC(IZID2,1))-
     &    UMIXC(IZID1,2)*DCONJG(UMIXC(IZID2,2))/2D0+DCMPLX(XW*DIFF)
          DO 230 I=MMINA,MMAXA
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 230
            EI=KCHG(IABS(I),1)/3D0
            T3I=SIGN(1D0,EI+1D-6)/2D0
            QRL=DCMPLX(-EI/SH*DIFF)-DCMPLX(EI/XW1)*PROPZ*ORP
            QLL=DCMPLX(-EI/SH*DIFF)+DCMPLX((T3I-XW*EI)/XW/XW1)*PROPZ*ORP
            QRR=DCMPLX(-EI/SH*DIFF)-DCMPLX(EI/XW1)*PROPZ*OLP
            IF(MOD(I,2).EQ.0) THEN
              XML2=PMAS(PYCOMP(KSUSY1+IABS(I)-1),1)**2
              QLR=DCMPLX(-EI/SH*DIFF)+DCMPLX((T3I-XW*EI)/XW/XW1)*
     &        PROPZ*OLP-UMIXC(IZID2,1)*DCONJG(UMIXC(IZID1,1))*
     &        DCMPLX(T3I/XW/(TH-XML2))
            ELSE
              XML2=PMAS(PYCOMP(KSUSY1+IABS(I)+1),1)**2
              QLR=DCMPLX(-EI/SH*DIFF)+DCMPLX((T3I-XW*EI)/XW/XW1)*
     &        PROPZ*OLP-VMIXC(IZID2,1)*DCONJG(VMIXC(IZID1,1))*
     &        DCMPLX(T3I/XW/(TH-XML2))
            ENDIF
            FCOL=1D0
            IF(IABS(I).GE.11) FCOL=3D0
            FACSUM=((ABS(QLL)**2*POLL+ABS(QRR)**2*POLR)*WU2+
     &      (ABS(QRL)**2*POLR+ABS(QLR)**2*POLL)*WT2+
     &      2D0*DBLE(QLR*DCONJG(QLL)*POLL+
     &      QRL*DCONJG(QRR)*POLR)*WS2)*FACXX*FCOL
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            IF(IZID1.EQ.IZID2) THEN
              SIGH(NCHN)=FACSUM*WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
            ELSE
              SIGH(NCHN)=FACSUM*WIDS(PYCOMP(KFPR(ISUBSV,1)),3)*
     &        WIDS(PYCOMP(KFPR(ISUBSV,2)),2)
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=-I
              ISIG(NCHN,3)=2
              SIGH(NCHN)=FACSUM*WIDS(PYCOMP(KFPR(ISUBSV,1)),2)*
     &        WIDS(PYCOMP(KFPR(ISUBSV,2)),3)
            ENDIF
  230     CONTINUE
 
        ELSEIF(ISUB.EQ.229) THEN
C...q + qbar' -> ~chi0_1 + ~chi+-_1
          FACXX=COMFAC*AEM**2/6D0/XW**2
          ZM12=SQM3
          ZM22=SQM4
          WU2 = (UH-ZM12)*(UH-ZM22)
          WT2 = (TH-ZM12)*(TH-ZM22)
          WS2 = SMW(IZID1)*SMZ(IZID2)*SH
          RT2I = 1D0/SQRT(2D0)
          PROPW = DCMPLX(SH-SQMW,-WWID*PMAS(24,1))/
     &    DCMPLX((SH-SQMW)**2+WWID**2*SQMW,0D0)
          DO 240 I=1,2
            VMIXC(IZID1,I)=DCMPLX(VMIX(IZID1,I),VMIXI(IZID1,I))
            UMIXC(IZID1,I)=DCMPLX(UMIX(IZID1,I),UMIXI(IZID1,I))
  240     CONTINUE
          DO 250 I=1,4
            ZMIXC(IZID2,I)=DCMPLX(ZMIX(IZID2,I),ZMIXI(IZID2,I))
  250     CONTINUE
          OL=(DCONJG(ZMIXC(IZID2,2))*VMIXC(IZID1,1)-
     &    DCONJG(ZMIXC(IZID2,4))*VMIXC(IZID1,2)*RT2I)*PROPW
          OR=(ZMIXC(IZID2,2)*DCONJG(UMIXC(IZID1,1))+
     &    ZMIXC(IZID2,3)*DCONJG(UMIXC(IZID1,2))*RT2I)*PROPW
 
          DO 270 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.20.OR.KFAC(1,I).EQ.0) GOTO 270
            EI=KCHG(IA,1)/3D0
            T3I=SIGN(1D0,EI+1D-6)/2D0
            DO 260 J=MMIN2,MMAX2
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.20.OR.KFAC(2,J).EQ.0) GOTO 260
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 260
              EJ=KCHG(JA,1)/3D0
              T3J=SIGN(1D0,EJ+1D-6)/2D0
              FCKM=3D0
              IF(IA.LE.10) FCKM=VCKM((IA+1)/2,(JA+1)/2)
              KCHSUM=KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J)
              KCHW=2
              IF(KCHSUM.LT.0) KCHW=3
              IF(MOD(IA,2).EQ.0) THEN
                ZMI2  = PMAS(PYCOMP(KSUSY1+IA),1)**2
                ZMJ2  = PMAS(PYCOMP(KSUSY1+JA),1)**2
                QLL=OL+VMIXC(IZID1,1)*DCONJG(ZMIXC(IZID2,1)*(EI-T3I)*
     &          TANW+ZMIXC(IZID2,2)*T3I)/DCMPLX(UH-ZMI2)
                QLR=OR-DCONJG(UMIXC(IZID1,1))*(
     &          ZMIXC(IZID2,1)*(EJ-T3J)*TANW+ZMIXC(IZID2,2)*T3J)
     &          /DCMPLX(TH-ZMJ2)
              ELSE
                ZMI2  = PMAS(PYCOMP(KSUSY1+JA),1)**2
                ZMJ2  = PMAS(PYCOMP(KSUSY1+IA),1)**2
                QLL=OL+VMIXC(IZID1,1)*DCONJG(ZMIXC(IZID2,1)*(EJ-T3J)*
     &          TANW+ZMIXC(IZID2,2)*T3J)/DCMPLX(UH-ZMJ2)
                QLR=OR-DCONJG(UMIXC(IZID1,1))*(
     &          ZMIXC(IZID2,1)*(EI-T3I)*TANW+ZMIXC(IZID2,2)*T3I)
     &          /DCMPLX(TH-ZMI2)
              ENDIF
              ZINTR=DBLE(QLR*DCONJG(QLL))
              FACGG1=FACXX*(ABS(QLL)**2*WU2+ABS(QLR)**2*WT2+
     &        2D0*ZINTR*WS2)
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACGG1*FCKM*WIDS(PYCOMP(KFPR(ISUBSV,1)),2)*
     &        WIDS(PYCOMP(KFPR(ISUBSV,2)),KCHW)
  260       CONTINUE
  270     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.240) THEN
        IF(ISUB.EQ.237) THEN
C...q + qbar -> gluino + ~chi0_1
          COMFAC=COMFAC*WIDS(PYCOMP(KFPR(ISUBSV,1)),2)*
     &    WIDS(PYCOMP(KFPR(ISUBSV,2)),2)
          ASYUK=RMSS(42)*AS
          FAC0=COMFAC*ASYUK*AEM*4D0/9D0/XW
          GM2=SQM3
          ZM2=SQM4
          DO 280 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 280
            EI=KCHG(IABS(I),1)/3D0
            IA=IABS(I)
            XLQC = -TANW*EI*ZMIX(IZID,1)
            XRQC =(SIGN(1D0,EI)*ZMIX(IZID,2)-TANW*
     &      (SIGN(1D0,EI)-2D0*EI)*ZMIX(IZID,1))/2D0
            XLQ2=XLQC**2
            XRQ2=XRQC**2
            XML2=PMAS(PYCOMP(KSUSY1+IA),1)**2
            XMR2=PMAS(PYCOMP(KSUSY2+IA),1)**2
            ATKIN=(TH-GM2)*(TH-ZM2)/(TH-XML2)**2
            AUKIN=(UH-GM2)*(UH-ZM2)/(UH-XML2)**2
            ATUKIN=SMZ(IZID)*SQRT(GM2)*SH/(TH-XML2)/(UH-XML2)
            SGCHIL=XLQ2*(ATKIN+AUKIN-2D0*ATUKIN)
            ATKIN=(TH-GM2)*(TH-ZM2)/(TH-XMR2)**2
            AUKIN=(UH-GM2)*(UH-ZM2)/(UH-XMR2)**2
            ATUKIN=SMZ(IZID)*SQRT(GM2)*SH/(TH-XMR2)/(UH-XMR2)
            SGCHIR=XRQ2*(ATKIN+AUKIN-2D0*ATUKIN)
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FAC0*(SGCHIL+SGCHIR)
  280     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.250) THEN
        IF(ISUB.EQ.241) THEN
C...q + qbar' -> ~chi+-_1 + gluino
          FACWG=COMFAC*AS*AEM/XW*2D0/9D0
          GM2=SQM3
          ZM2=SQM4
          FAC01=2D0*UMIX(IZID,1)*VMIX(IZID,1)
          FAC0=UMIX(IZID,1)**2
          FAC1=VMIX(IZID,1)**2
          DO 300 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.10.OR.KFAC(1,I).EQ.0) GOTO 300
            DO 290 J=MMIN2,MMAX2
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.10.OR.KFAC(2,J).EQ.0) GOTO 290
              IF(I*J.GT.0.OR.MOD(IA+JA,2).EQ.0) GOTO 290
              FCKM=1D0
              IF(IA.LE.10) FCKM=VCKM((IA+1)/2,(JA+1)/2)
              KCHSUM=KCHG(IA,1)*ISIGN(1,I)+KCHG(JA,1)*ISIGN(1,J)
              KCHW=2
              IF(KCHSUM.LT.0) KCHW=3
              XMU2=PMAS(PYCOMP(KSUSY1+2),1)**2
              XMD2=PMAS(PYCOMP(KSUSY1+1),1)**2
              ATKIN=(TH-GM2)*(TH-ZM2)/(TH-XMU2)**2
              AUKIN=(UH-GM2)*(UH-ZM2)/(UH-XMD2)**2
              ATUKIN=SMW(IZID)*SQRT(GM2)*SH/(TH-XMU2)/(UH-XMD2)
              XMU2=PMAS(PYCOMP(KSUSY2+2),1)**2
              XMD2=PMAS(PYCOMP(KSUSY2+1),1)**2
              ATKIN=(ATKIN+(TH-GM2)*(TH-ZM2)/(TH-XMU2)**2)/2D0
              AUKIN=(AUKIN+(UH-GM2)*(UH-ZM2)/(UH-XMD2)**2)/2D0
              ATUKIN=(ATUKIN+SMW(IZID)*SQRT(GM2)*
     &        SH/(TH-XMU2)/(UH-XMD2))/2D0
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACWG*FCKM*(FAC0*ATKIN+FAC1*AUKIN-
     &        FAC01*ATUKIN)*WIDS(PYCOMP(KFPR(ISUBSV,1)),2)*
     &        WIDS(PYCOMP(KFPR(ISUBSV,2)),KCHW)
  290       CONTINUE
  300     CONTINUE
 
        ELSEIF(ISUB.EQ.243) THEN
C...q + qbar -> gluino + gluino
          COMFAC=COMFAC*WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
          XMT=SQM3-TH
          XMU=SQM3-UH
          DO 310 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 310
            NCHN=NCHN+1
            XSU=PMAS(PYCOMP(KSUSY1+IABS(I)),1)**2-UH
            XST=PMAS(PYCOMP(KSUSY1+IABS(I)),1)**2-TH
            FACGG1=COMFAC*AS**2*8D0/3D0*( (XMT**2+XMU**2+
     &      2D0*SQM3*SH)/SH2 + RMSS(42)**2*(4D0/9D0*(XMT**2/XST**2+
     &      XMU**2/XSU**2) + SQM3*SH/XST/XSU/9D0) - RMSS(42)*(
     &      (XMT**2+SH*SQM3)/SH/XST + (XMU**2+SH*SQM3)/SH/XSU ))
            XSU=PMAS(PYCOMP(KSUSY2+IABS(I)),1)**2-UH
            XST=PMAS(PYCOMP(KSUSY2+IABS(I)),1)**2-TH
            FACGG2=COMFAC*AS**2*8D0/3D0*( (XMT**2+XMU**2+
     &      2D0*SQM3*SH)/SH2 + RMSS(42)**2*(4D0/9D0*(XMT**2/XST**2+
     &      XMU**2/XSU**2) + SQM3*SH/XST/XSU/9D0) - RMSS(42)*(
     &      (XMT**2+SH*SQM3)/SH/XST + (XMU**2+SH*SQM3)/SH/XSU ))
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
C...1/2 for identical particles
            SIGH(NCHN)=0.25D0*(FACGG1+FACGG2)
  310     CONTINUE
 
        ELSEIF(ISUB.EQ.244) THEN
C...g + g -> gluino + gluino
          COMFAC=COMFAC*WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
          XMT=SQM3-TH
          XMU=SQM3-UH
          FACQQ1=COMFAC*AS**2*9D0/4D0*(
     &    (XMT*XMU-2D0*SQM3*(TH+SQM3))/XMT**2 -
     &    (XMT*XMU+SQM3*(UH-TH))/SH/XMT )
          FACQQ2=COMFAC*AS**2*9D0/4D0*(
     &    (XMU*XMT-2D0*SQM3*(UH+SQM3))/XMU**2 -
     &    (XMU*XMT+SQM3*(TH-UH))/SH/XMU )
          FACQQ3=COMFAC*AS**2*9D0/4D0*(2D0*XMT*XMU/SH2 +
     &    SQM3*(SH-4D0*SQM3)/XMT/XMU)
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 320
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ1/2D0
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=2
          SIGH(NCHN)=FACQQ2/2D0
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=3
          SIGH(NCHN)=FACQQ3/2D0
  320     CONTINUE
 
        ELSEIF(ISUB.EQ.246) THEN
C...g + q_j -> ~chi0_1 + ~q_j
          FAC0=COMFAC*AS*AEM/6D0/XW
          ZM2=SQM4
          QM2=SQM3
          FACZQ0=FAC0*( (ZM2-TH)/SH +
     &    (UH-ZM2)*(UH+QM2)/(UH-QM2)**2 -
     &    (SH*(UH+ZM2)+2D0*(QM2-ZM2)*(ZM2-UH))/SH/(UH-QM2) )
          KFNSQ=MOD(KFPR(ISUBSV,1),KSUSY1)
          DO 340 I=-KFNSQ,KFNSQ,2*KFNSQ
            IF(I.LT.MMINA.OR.I.GT.MMAXA) GOTO 340
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 340
            EI=KCHG(IABS(I),1)/3D0
            IA=IABS(I)
            XRQZ = -TANW*EI*ZMIX(IZID,1)
            XLQZ =(SIGN(1D0,EI)*ZMIX(IZID,2)-TANW*
     &      (SIGN(1D0,EI)-2D0*EI)*ZMIX(IZID,1))/2D0
            IF(ILR.EQ.0) THEN
              BS=XLQZ**2*SFMIX(IA,1)**2+XRQZ**2*SFMIX(IA,2)**2
            ELSE
              BS=XLQZ**2*SFMIX(IA,3)**2+XRQZ**2*SFMIX(IA,4)**2
            ENDIF
            FACZQ=FACZQ0*BS
            KCHQ=2
            IF(I.LT.0) KCHQ=3
            DO 330 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 330
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 330
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACZQ*RKF*WIDS(PYCOMP(KFPR(ISUBSV,1)),KCHQ)*
     &        WIDS(PYCOMP(KFPR(ISUBSV,2)),2)
  330       CONTINUE
  340     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.260) THEN
        IF(ISUB.EQ.254) THEN
C...g + q_j -> ~chi1_1 + ~q_i
          FAC0=COMFAC*AS*AEM/12D0/XW
          ZM2=SQM4
          QM2=SQM3
          AU=UMIX(IZID,1)**2
          AD=VMIX(IZID,1)**2
          FACZQ0=FAC0*( (ZM2-TH)/SH +
     &    (UH-ZM2)*(UH+QM2)/(UH-QM2)**2 -
     &    (SH*(UH+ZM2)+2D0*(QM2-ZM2)*(ZM2-UH))/SH/(UH-QM2) )
          KFNSQ1=MOD(KFPR(ISUBSV,1),KSUSY1)
          IF(MOD(KFNSQ1,2).EQ.0) THEN
            KFNSQ=KFNSQ1-1
            KCHW=2
          ELSE
            KFNSQ=KFNSQ1+1
            KCHW=3
          ENDIF
          DO 360 I=-KFNSQ,KFNSQ,2*KFNSQ
            IF(I.LT.MMINA.OR.I.GT.MMAXA) GOTO 360
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58)) GOTO 360
            IA=IABS(I)
            IF(MOD(IA,2).EQ.0) THEN
              FACZQ=FACZQ0*AU
            ELSE
              FACZQ=FACZQ0*AD
            ENDIF
            FACZQ=FACZQ*SFMIX(KFNSQ1,1+2*ILR)**2
            KCHQ=2
            IF(I.LT.0) KCHQ=3
            KCHWQ=KCHW
            IF(I.LT.0) KCHWQ=5-KCHW
            DO 350 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 350
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 350
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACZQ*RKF*WIDS(PYCOMP(KFPR(ISUBSV,1)),KCHQ)*
     &        WIDS(PYCOMP(KFPR(ISUBSV,2)),KCHWQ)
  350       CONTINUE
  360     CONTINUE
 
        ELSEIF(ISUB.EQ.258) THEN
C...g + q_j -> gluino + ~q_i
          XG2=SQM4
          XQ2=SQM3
          XMT=XG2-TH
          XMU=XG2-UH
          XST=XQ2-TH
          XSU=XQ2-UH
          FACQG1=0.5D0*4D0/9D0*XMT/SH + (XMT*SH+2D0*XG2*XST)/XMT**2 -
     &    ( (SH-XQ2+XG2)*(-XST)-SH*XG2 )/SH/(-XMT) +
     &    0.5D0*1D0/2D0*( XST*(TH+2D0*UH+XG2)-XMT*(SH-2D0*XST) +
     &    (-XMU)*(TH+XG2+2D0*XQ2) )/2D0/XMT/XSU
          FACQG2= 4D0/9D0*(-XMU)*(UH+XQ2)/XSU**2 + 1D0/18D0*
     &    (SH*(UH+XG2)
     &    +2D0*(XQ2-XG2)*XMU)/SH/(-XSU) + 0.5D0*4D0/9D0*XMT/SH +
     &    0.5D0*1D0/2D0*(XST*(TH+2D0*UH+XG2)-XMT*(SH-2D0*XST)+
     &    (-XMU)*(TH+XG2+2D0*XQ2))/2D0/XMT/XSU
          ASYUK=RMSS(42)*AS
          FACQG1=COMFAC*AS*ASYUK*FACQG1/2D0
          FACQG2=COMFAC*AS*ASYUK*FACQG2/2D0
          KFNSQ=MOD(KFPR(ISUBSV,1),KSUSY1)
          DO 380 I=-KFNSQ,KFNSQ,2*KFNSQ
            IF(I.LT.MMINA.OR.I.GT.MMAXA) GOTO 380
            IF(I.EQ.0.OR.IABS(I).GT.10) GOTO 380
            KCHQ=2
            IF(I.LT.0) KCHQ=3
            FACSEL=RKF*WIDS(PYCOMP(KFPR(ISUBSV,1)),KCHQ)*
     &      WIDS(PYCOMP(KFPR(ISUBSV,2)),2)
            DO 370 ISDE=1,2
              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 370
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 370
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQG1*FACSEL
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=2
              SIGH(NCHN)=FACQG2*FACSEL
  370       CONTINUE
  380     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.270) THEN
        IF(ISUB.EQ.261) THEN
C...q_i + q_ibar -> ~t_1 + ~t_1bar
          FACQQ1=COMFAC*( (UH*TH-SQM3*SQM4)/ SH**2 )*
     &    WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
          KFNSQ=MOD(KFPR(ISUBSV,1),KSUSY1)
          FAC0=AS**2*4D0/9D0
          DO 390 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 390
            IF(IA.GE.11.AND.IA.LE.18) THEN
              EI=KCHG(IA,1)/3D0
              EJ=KCHG(KFNSQ,1)/3D0
              T3I=SIGN(1D0,EI)/2D0
              T3J=SIGN(1D0,EJ)/2D0
              XLQ=2D0*(T3J-EJ*XW)*SFMIX(KFNSQ,2*ILR+1)**2
              XRQ=2D0*(-EJ*XW)*SFMIX(KFNSQ,2*ILR+2)**2
              XLF=2D0*(T3I-EI*XW)
              XRF=2D0*(-EI*XW)
              TAA=0.5D0*(EI*EJ)**2
              TZZ=(XLF**2+XRF**2)*(XLQ+XRQ)**2/64D0/XW**2/XW1**2
              TZZ=TZZ/((1D0-SQMZ/SH)**2+SQMZ*(ZWID/SH)**2)
              TAZ=EI*EJ*(XLQ+XRQ)*(XLF+XRF)/8D0/XW/XW1
              TAZ=TAZ/((1D0-SQMZ/SH)**2+SQMZ*(ZWID/SH)**2)*(1D0-SQMZ/SH)
              FAC0=AEM**2*12D0*(TAA+TZZ+TAZ)
            ENDIF
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQ1*FAC0
  390     CONTINUE
 
        ELSEIF(ISUB.EQ.263) THEN
C...f + fbar -> ~t1 + ~t2bar
          DO 400 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 400
            EI=KCHG(IABS(I),1)/3D0
            TT3I=SIGN(1D0,EI)/2D0
            EJ=2D0/3D0
            TT3J=1D0/2D0
            FCOL=1D0
C...Color factor for e+ e-
            IF(IA.GE.11) FCOL=3D0
            XLQ=2D0*(TT3J-EJ*XW)
            XRQ=2D0*(-EJ*XW)
            XLF=2D0*(TT3I-EI*XW)
            XRF=2D0*(-EI*XW)
            TZZ=(XLF**2+XRF**2)*(XLQ-XRQ)**2/64D0/XW**2/XW1**2
            TZZ=TZZ*(SFMIX(6,1)*SFMIX(6,2))**2
            TZZ=TZZ/((1D0-SQMZ/SH)**2+SQMZ*(ZWID/SH)**2)
C...Factor of 2 for t1 t2bar + t2 t1bar
C...PS: bug fix 24 Aug 2010. Factor 2 accounted for by the 2 channels.
            FACQQ1=COMFAC*AEM**2*TZZ*FCOL*4D0
            FACQQ1=FACQQ1*( UH*TH-SQM3*SQM4 )/SH2
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQ1*WIDS(PYCOMP(KFPR(ISUBSV,1)),2)*
     &      WIDS(PYCOMP(KFPR(ISUBSV,2)),3)
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=2
            SIGH(NCHN)=FACQQ1*WIDS(PYCOMP(KFPR(ISUBSV,1)),3)*
     &      WIDS(PYCOMP(KFPR(ISUBSV,2)),2)
  400     CONTINUE
 
        ELSEIF(ISUB.EQ.264) THEN
C...g + g -> ~t_1 + ~t_1bar
          XSU=SQM3-UH
          XST=SQM3-TH
          FAC0=COMFAC*AS**2*(7D0/48D0+3D0*(UH-TH)**2/16D0/SH2 )*0.5D0*
     &    WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
          FACQQ1=FAC0*(0.5D0+2D0*SQM3*TH/XST**2 + 2D0*SQM3**2/XSU/XST)
          FACQQ2=FAC0*(0.5D0+2D0*SQM3*UH/XSU**2 + 2D0*SQM3**2/XSU/XST)
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 410
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
  410     CONTINUE
        ENDIF
 
      ELSEIF(ISUB.LE.280) THEN
        IF(ISUB.EQ.271) THEN
C...q + q' -> ~q + ~q' (~g exchange)
          XMG2=PMAS(PYCOMP(KSUSY1+21),1)**2
          XMT=XMG2-TH
          XMU=XMG2-UH
          XSU1=SQM3-UH
          XSU2=SQM4-UH
          XST1=SQM3-TH
          XST2=SQM4-TH
          ASYUK=RMSS(42)*AS
          IF(ILR.EQ.1) THEN
            FACQQ1=COMFAC*ASYUK**2*4D0/9D0*( -(XST1*XST2+SH*TH)/XMT**2 )
            FACQQ2=COMFAC*ASYUK**2*4D0/9D0*( -(XSU1*XSU2+SH*UH)/XMU**2 )
            FACQQB=0.0D0
          ELSE
            FACQQ1=0.5D0*COMFAC*ASYUK**2*4D0/9D0*( SH*XMG2/XMT**2 )
            FACQQ2=0.5D0*COMFAC*ASYUK**2*4D0/9D0*( SH*XMG2/XMU**2 )
            FACQQB=0.5D0*COMFAC*ASYUK**2*4D0/9D0*( -2D0*SH*XMG2/3D0/
     &      XMT/XMU )
          ENDIF
          KFNSQI=MOD(KFPR(ISUBSV,1),KSUSY1)
          KFNSQJ=MOD(KFPR(ISUBSV,2),KSUSY1)
          DO 430 I=-KFNSQI,KFNSQI,2*KFNSQI
            IF(I.LT.MMIN1.OR.I.GT.MMAX1) GOTO 430
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.MSTP(58).OR.KFAC(1,I).EQ.0) GOTO 430
            KCHQ=2
            IF(I.LT.0) KCHQ=3
            DO 420 J=-KFNSQJ,KFNSQJ,2*KFNSQJ
              IF(J.LT.MMIN2.OR.J.GT.MMAX2) GOTO 420
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.MSTP(58).OR.KFAC(2,J).EQ.0) GOTO 420
              IF(I*J.LT.0) GOTO 420
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQQ1*RKF*WIDS(PYCOMP(KFPR(ISUBSV,1)),KCHQ)*
     &        WIDS(PYCOMP(KFPR(ISUBSV,2)),KCHQ)
              IF(I.EQ.J) THEN
                IF(ILR.EQ.0) THEN
                  SIGH(NCHN)=0.5D0*(FACQQ1+0.5D0*FACQQB)*RKF*
     &            WIDS(PYCOMP(KFPR(ISUBSV,1)),KCHQ+2)
                ELSE
                  SIGH(NCHN)=0.5D0*FACQQ1*RKF*
     &            WIDS(PYCOMP(KFPR(ISUBSV,1)),KCHQ)*
     &            WIDS(PYCOMP(KFPR(ISUBSV,2)),KCHQ)
                ENDIF
                NCHN=NCHN+1
                ISIG(NCHN,1)=I
                ISIG(NCHN,2)=J
                ISIG(NCHN,3)=2
                IF(ILR.EQ.0) THEN
                  SIGH(NCHN)=0.5D0*(FACQQ2+0.5D0*FACQQB)*RKF*
     &            WIDS(PYCOMP(KFPR(ISUBSV,1)),KCHQ+2)
                ELSE
                  SIGH(NCHN)=0.5D0*FACQQ2*RKF*
     &            WIDS(PYCOMP(KFPR(ISUBSV,1)),KCHQ)*
     &            WIDS(PYCOMP(KFPR(ISUBSV,2)),KCHQ)
                ENDIF
              ENDIF
  420       CONTINUE
  430     CONTINUE
 
        ELSEIF(ISUB.EQ.274) THEN
C...q + qbar' -> ~q + ~qbar'
          XMG2=PMAS(PYCOMP(KSUSY1+21),1)**2
          XMT=XMG2-TH
          XMU=XMG2-UH
          IF(ILR.EQ.0) THEN
C...Mrenna...Normalization.and.1/XMT
            FACQQ1=COMFAC*AS**2*2D0/9D0*(
     &      (UH*TH-SQM3*SQM4)/XMT**2 )*RMSS(42)**2
            FACQQB=COMFAC*AS**2*4D0/9D0*(
     &      (UH*TH-SQM3*SQM4)/SH2 )
C...Mrenna..Switched sign to agree with Eichten, Dawson, etc.
            FACQQI=COMFAC*AS**2*4D0/27D0*(
     &      (UH*TH-SQM3*SQM4)/SH/XMT )*RMSS(42)
            FACQQB=FACQQB+FACQQ1+FACQQI
          ELSE
            FACQQ1=COMFAC*AS**2*4D0/9D0*( XMG2*SH/XMT**2 )*RMSS(42)**2
            FACQQB=FACQQ1
          ENDIF
          KFNSQI=MOD(KFPR(ISUBSV,1),KSUSY1)
          KFNSQJ=MOD(KFPR(ISUBSV,2),KSUSY1)
          DO 450 I=-KFNSQI,KFNSQI,2*KFNSQI
            IF(I.LT.MMIN1.OR.I.GT.MMAX1) GOTO 450
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.MSTP(58).OR.KFAC(1,I).EQ.0) GOTO 450
            KCHQ=2
            IF(I.LT.0) KCHQ=3
            DO 440 J=-KFNSQJ,KFNSQJ,2*KFNSQJ
              IF(J.LT.MMIN2.OR.J.GT.MMAX2) GOTO 440
              JA=IABS(J)
              IF(J.EQ.0.OR.JA.GT.MSTP(58).OR.KFAC(2,J).EQ.0) GOTO 440
              IF(I*J.GT.0) GOTO 440
              NCHN=NCHN+1
              ISIG(NCHN,1)=I
              ISIG(NCHN,2)=J
              ISIG(NCHN,3)=1
              SIGH(NCHN)=FACQQ1*RKF*WIDS(PYCOMP(KFPR(ISUBSV,1)),KCHQ)*
     &        WIDS(PYCOMP(KFPR(ISUBSV,2)),5-KCHQ)
              IF(ILR.EQ.0.AND.I.EQ.-J) SIGH(NCHN)=FACQQB*RKF*
     &        WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
  440       CONTINUE
  450     CONTINUE
 
        ELSEIF(ISUB.EQ.277) THEN
C...q_i + q_ibar -> ~q_j + ~q_jbar ,i .ne. j
C...if i .eq. j covered in 274
          FACQQ1=COMFAC*( (UH*TH-SQM3*SQM4)/ SH**2 )
          KFNSQ=MOD(KFPR(ISUBSV,1),KSUSY1)
          FAC0=0D0
          DO 460 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.(IA.GT.MSTP(58).AND.IA.LE.10).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 460
            IF(IA.EQ.KFNSQ) GOTO 460
            IF(IA.EQ.11.OR.IA.EQ.13.OR.IA.EQ.15) THEN
              EI=KCHG(IA,1)/3D0
              EJ=KCHG(KFNSQ,1)/3D0
              T3J=SIGN(0.5D0,EJ)
              T3I=SIGN(1D0,EI)/2D0
              IF(ILR.EQ.0) THEN
                XLQ=2D0*(T3J-EJ*XW)*SFMIX(KFNSQ,1)
                XRQ=2D0*(-EJ*XW)*SFMIX(KFNSQ,2)
              ELSE
                XLQ=2D0*(T3J-EJ*XW)*SFMIX(KFNSQ,3)
                XRQ=2D0*(-EJ*XW)*SFMIX(KFNSQ,4)
              ENDIF
              XLF=2D0*(T3I-EI*XW)
              XRF=2D0*(-EI*XW)
              IF(ILR.EQ.0) THEN
                XRQ=0D0
              ELSE
                XLQ=0D0
              ENDIF
              TAA=0.5D0*(EI*EJ)**2
              TZZ=(XLF**2+XRF**2)*(XLQ+XRQ)**2/64D0/XW**2/XW1**2
              TZZ=TZZ/((1D0-SQMZ/SH)**2+SQMZ*(ZWID/SH)**2)
              TAZ=EI*EJ*(XLQ+XRQ)*(XLF+XRF)/8D0/XW/XW1
              TAZ=TAZ/((1D0-SQMZ/SH)**2+SQMZ*(ZWID/SH)**2)*(1D0-SQMZ/SH)
              FAC0=AEM**2*12D0*(TAA+TZZ+TAZ)
            ELSEIF(IA.LE.6) THEN
              FAC0=AS**2*8D0/9D0/2D0
            ENDIF
            NCHN=NCHN+1
            ISIG(NCHN,1)=I
            ISIG(NCHN,2)=-I
            ISIG(NCHN,3)=1
            SIGH(NCHN)=FACQQ1*FAC0*RKF*WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
  460     CONTINUE
 
        ELSEIF(ISUB.EQ.279) THEN
C...g + g -> ~q_j + ~q_jbar
          XSU=SQM3-UH
          XST=SQM3-TH
C...4=RKF because ~t ~tbar and ~b ~bbar treated separately
          FAC0=RKF*COMFAC*AS**2*( 7D0/48D0+3D0*(UH-TH)**2/16D0/SH2 )
          FACQQ1=FAC0*(0.5D0+2D0*SQM3*TH/XST**2 + 2D0*SQM3**2/XSU/XST)
          FACQQ2=FAC0*(0.5D0+2D0*SQM3*UH/XSU**2 + 2D0*SQM3**2/XSU/XST)
          IF(KFAC(1,21)*KFAC(2,21).EQ.0) GOTO 470
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=1
          SIGH(NCHN)=FACQQ1/2D0*WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
          NCHN=NCHN+1
          ISIG(NCHN,1)=21
          ISIG(NCHN,2)=21
          ISIG(NCHN,3)=2
          SIGH(NCHN)=FACQQ2/2D0*WIDS(PYCOMP(KFPR(ISUBSV,1)),1)
  470     CONTINUE
 
        ENDIF
      ENDIF
CMRENNA--
 
      RETURN
      END
