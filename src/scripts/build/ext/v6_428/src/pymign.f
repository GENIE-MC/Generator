 
C*********************************************************************
 
C...PYMIGN
C...Initializes treatment of new multiple interactions scenario,
C...selects kinematics of hardest interaction if low-pT physics
C...included in run, and generates all non-hardest interactions.
 
      SUBROUTINE PYMIGN(MMUL)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      EXTERNAL PYALPS
      DOUBLE PRECISION PYALPS
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      COMMON/PYINTM/KFIVAL(2,3),NMI(2),IMI(2,800,2),NVC(2,-6:6),
     &     XASSOC(2,-6:6,240),XPSVC(-6:6,-1:240),PVCTOT(2,-1:1),
     &     XMI(2,240),PT2MI(240),IMISEP(0:240)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,
     &/PYINT1/,/PYINT2/,/PYINT3/,/PYINT5/,/PYINT7/,/PYINTM/
C...Local arrays and saved variables.
      DIMENSION NMUL(20),SIGM(20),KSTR(500,2),VINTSV(80),
     &WDTP(0:400),WDTE(0:400,0:5),XPQ(-25:25),KSAV(4,5),PSAV(4,5)
      SAVE XT2,XT2FAC,XC2,XTS,IRBIN,RBIN,NMUL,SIGM,P83A,P83B,P83C,
     &CQ2I,CQ2R,PIK,BDIV,B,PLOWB,PHIGHB,PALLB,S4A,S4B,S4C,POWIP,
     &RPWIP,B2RPDV,B2RPMX,BAVG,VNT145,VNT146,VNT147
 
C...Initialization of multiple interaction treatment.
      IF(MMUL.EQ.1) THEN
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5000) MSTP(82)
        ISUB=96
        MINT(1)=96
        VINT(63)=0D0
        VINT(64)=0D0
        VINT(143)=1D0
        VINT(144)=1D0
 
C...Loop over phase space points: xT2 choice in 20 bins.
  100   SIGSUM=0D0
        DO 120 IXT2=1,20
          NMUL(IXT2)=MSTP(83)
          SIGM(IXT2)=0D0
          DO 110 ITRY=1,MSTP(83)
            RSCA=0.05D0*((21-IXT2)-PYR(0))
            XT2=VINT(149)*(1D0+VINT(149))/(VINT(149)+RSCA)-VINT(149)
            XT2=MAX(0.01D0*VINT(149),XT2)
            VINT(25)=XT2
 
C...Choose tau and y*. Calculate cos(theta-hat).
            IF(PYR(0).LE.COEF(ISUB,1)) THEN
              TAUT=(2D0*(1D0+SQRT(1D0-XT2))/XT2-1D0)**PYR(0)
              TAU=XT2*(1D0+TAUT)**2/(4D0*TAUT)
            ELSE
              TAU=XT2*(1D0+TAN(PYR(0)*ATAN(SQRT(1D0/XT2-1D0)))**2)
            ENDIF
            VINT(21)=TAU
            CALL PYKLIM(2)
            RYST=PYR(0)
            MYST=1
            IF(RYST.GT.COEF(ISUB,8)) MYST=2
            IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)) MYST=3
            CALL PYKMAP(2,MYST,PYR(0))
            VINT(23)=SQRT(MAX(0D0,1D0-XT2/TAU))*(-1)**INT(1.5D0+PYR(0))
 
C...Calculate differential cross-section.
            VINT(71)=0.5D0*VINT(1)*SQRT(XT2)
            CALL PYSIGH(NCHN,SIGS)
            SIGM(IXT2)=SIGM(IXT2)+SIGS
  110     CONTINUE
          SIGSUM=SIGSUM+SIGM(IXT2)
  120   CONTINUE
        SIGSUM=SIGSUM/(20D0*MSTP(83))
 
C...Reject result if sigma(parton-parton) is smaller than hadronic one.
        IF(SIGSUM.LT.1.1D0*SIGT(0,0,5)) THEN
          IF(MSTP(122).GE.1) WRITE(MSTU(11),5100)
     &    PARP(82)*(VINT(1)/PARP(89))**PARP(90),SIGSUM
          PARP(82)=0.9D0*PARP(82)
          VINT(149)=4D0*(PARP(82)*(VINT(1)/PARP(89))**PARP(90))**2/
     &    VINT(2)
          GOTO 100
        ENDIF
        IF(MSTP(122).GE.1) WRITE(MSTU(11),5200)
     &  PARP(82)*(VINT(1)/PARP(89))**PARP(90), SIGSUM
 
C...Start iteration to find k factor.
        YKE=SIGSUM/MAX(1D-10,SIGT(0,0,5))
        P83A=(1D0-PARP(83))**2
        P83B=2D0*PARP(83)*(1D0-PARP(83))
        P83C=PARP(83)**2
        CQ2I=1D0/PARP(84)**2
        CQ2R=2D0/(1D0+PARP(84)**2)
        SO=0.5D0
        XI=0D0
        YI=0D0
        XF=0D0
        YF=0D0
        XK=0.5D0
        IIT=0
  130   IF(IIT.EQ.0) THEN
          XK=2D0*XK
        ELSEIF(IIT.EQ.1) THEN
          XK=0.5D0*XK
        ELSE
          XK=XI+(YKE-YI)*(XF-XI)/(YF-YI)
        ENDIF
 
C...Evaluate overlap integrals. Find where to divide the b range.
        IF(MSTP(82).EQ.2) THEN
          SP=0.5D0*PARU(1)*(1D0-EXP(-XK))
          SOP=SP/PARU(1)
        ELSE
          IF(MSTP(82).EQ.3) THEN
            DELTAB=0.02D0
          ELSEIF(MSTP(82).EQ.4) THEN
            DELTAB=MIN(0.01D0,0.05D0*PARP(84))
          ELSE
            POWIP=MAX(0.4D0,PARP(83))
            RPWIP=2D0/POWIP-1D0
            DELTAB=MAX(0.02D0,0.02D0*(2D0/POWIP)**(1D0/POWIP))
            SO=0D0
          ENDIF
          SP=0D0
          SOP=0D0
          BSP=0D0
          SOHIGH=0D0
          IBDIV=0
          B=-0.5D0*DELTAB
  140     B=B+DELTAB
          IF(MSTP(82).EQ.3) THEN
            OV=EXP(-B**2)/PARU(2)
          ELSEIF(MSTP(82).EQ.4) THEN
            OV=(P83A*EXP(-MIN(50D0,B**2))+
     &      P83B*CQ2R*EXP(-MIN(50D0,B**2*CQ2R))+
     &      P83C*CQ2I*EXP(-MIN(50D0,B**2*CQ2I)))/PARU(2)
          ELSE
            OV=EXP(-B**POWIP)/PARU(2)
            SO=SO+PARU(2)*B*DELTAB*OV
          ENDIF
          IF(IBDIV.EQ.1) SOHIGH=SOHIGH+PARU(2)*B*DELTAB*OV
          PACC=1D0-EXP(-MIN(50D0,PARU(1)*XK*OV))
          SP=SP+PARU(2)*B*DELTAB*PACC
          SOP=SOP+PARU(2)*B*DELTAB*OV*PACC
          BSP=BSP+B*PARU(2)*B*DELTAB*PACC
          IF(IBDIV.EQ.0.AND.PARU(1)*XK*OV.LT.1D0) THEN
            IBDIV=1 
            BDIV=B+0.5D0*DELTAB
          ENDIF
          IF(B.LT.1D0.OR.B*PACC.GT.1D-6) GOTO 140
        ENDIF
        YK=PARU(1)*XK*SO/SP
 
C...Continue iteration until convergence.
        IF(YK.LT.YKE) THEN
          XI=XK
          YI=YK
          IF(IIT.EQ.1) IIT=2
        ELSE
          XF=XK
          YF=YK
          IF(IIT.EQ.0) IIT=1
        ENDIF
        IF(ABS(YK-YKE).GE.1D-5*YKE) GOTO 130
 
C...Store some results for subsequent use.
        BAVG=BSP/SP
        VINT(145)=SIGSUM
        VINT(146)=SOP/SO
        VINT(147)=SOP/SP
        VNT145=VINT(145)
        VNT146=VINT(146)
        VNT147=VINT(147)
C...PIK = PARU(1)*XK = (VINT(146)/VINT(147))*sigma_jet/sigma_nondiffr.
        PIK=(VNT146/VNT147)*YKE

C...Find relative weight for low and high impact parameter..
      PLOWB=PARU(1)*BDIV**2
      IF(MSTP(82).EQ.3) THEN
        PHIGHB=PIK*0.5*EXP(-BDIV**2)
      ELSEIF(MSTP(82).EQ.4) THEN
        S4A=P83A*EXP(-BDIV**2)
        S4B=P83B*EXP(-BDIV**2*CQ2R)
        S4C=P83C*EXP(-BDIV**2*CQ2I)
        PHIGHB=PIK*0.5*(S4A+S4B+S4C)
      ELSEIF(PARP(83).GE.1.999D0) THEN
        PHIGHB=PIK*SOHIGH
        B2RPDV=BDIV**POWIP
      ELSE
        PHIGHB=PIK*SOHIGH
        B2RPDV=BDIV**POWIP
        B2RPMX=MAX(2D0*RPWIP,B2RPDV)
      ENDIF 
      PALLB=PLOWB+PHIGHB
 
C...Initialize iteration in xT2 for hardest interaction.
      ELSEIF(MMUL.EQ.2) THEN
        VINT(145)=VNT145
        VINT(146)=VNT146
        VINT(147)=VNT147
        IF(MSTP(82).LE.0) THEN
        ELSEIF(MSTP(82).EQ.1) THEN
          XT2=1D0
          SIGRAT=XSEC(96,1)/MAX(1D-10,VINT(315)*VINT(316)*SIGT(0,0,5))
          IF(MINT(141).NE.0.OR.MINT(142).NE.0) SIGRAT=SIGRAT*
     &    VINT(317)/(VINT(318)*VINT(320))
          XT2FAC=SIGRAT*VINT(149)/(1D0-VINT(149))
        ELSEIF(MSTP(82).EQ.2) THEN
          XT2=1D0
          XT2FAC=VNT146*XSEC(96,1)/MAX(1D-10,SIGT(0,0,5))*
     &    VINT(149)*(1D0+VINT(149))
        ELSE
          XC2=4D0*CKIN(3)**2/VINT(2)
          IF(CKIN(3).LE.CKIN(5).OR.MINT(82).GE.2) XC2=0D0
        ENDIF

C...Select impact parameter for hardest interaction.
        IF(MSTP(82).LE.2) RETURN
  142   IF(PYR(0)*PALLB.LT.PLOWB) THEN
C...Treatment in low b region.
          MINT(39)=1
          B=BDIV*SQRT(PYR(0)) 
          IF(MSTP(82).EQ.3) THEN
            OV=EXP(-B**2)/PARU(2)
          ELSEIF(MSTP(82).EQ.4) THEN
            OV=(P83A*EXP(-MIN(50D0,B**2))+
     &      P83B*CQ2R*EXP(-MIN(50D0,B**2*CQ2R))+
     &      P83C*CQ2I*EXP(-MIN(50D0,B**2*CQ2I)))/PARU(2)
          ELSE
            OV=EXP(-B**POWIP)/PARU(2)
          ENDIF  
          VINT(148)=OV/VNT147
          PACC=1D0-EXP(-MIN(50D0,PIK*OV))
          XT2=1D0
          XT2FAC=VNT146*VINT(148)*XSEC(96,1)/MAX(1D-10,SIGT(0,0,5))*
     &    VINT(149)*(1D0+VINT(149))
        ELSE
C...Treatment in high b region.
          MINT(39)=2
          IF(MSTP(82).EQ.3) THEN
            B=SQRT(BDIV**2-LOG(PYR(0)))
            OV=EXP(-B**2)/PARU(2)
          ELSEIF(MSTP(82).EQ.4) THEN
            S4RNDM=PYR(0)*(S4A+S4B+S4C)
            IF(S4RNDM.LT.S4A) THEN
              B=SQRT(BDIV**2-LOG(PYR(0)))
            ELSEIF(S4RNDM.LT.S4A+S4B) THEN
              B=SQRT(BDIV**2-LOG(PYR(0))/CQ2R)
            ELSE
              B=SQRT(BDIV**2-LOG(PYR(0))/CQ2I)
            ENDIF    
            OV=(P83A*EXP(-MIN(50D0,B**2))+
     &      P83B*CQ2R*EXP(-MIN(50D0,B**2*CQ2R))+
     &      P83C*CQ2I*EXP(-MIN(50D0,B**2*CQ2I)))/PARU(2)
          ELSEIF(PARP(83).GE.1.999D0) THEN
  144       B2RPW=B2RPDV-LOG(PYR(0))
            ACCIP=(B2RPW/B2RPDV)**RPWIP
            IF(ACCIP.LT.PYR(0)) GOTO 144
            OV=EXP(-B2RPW)/PARU(2)
            B=B2RPW**(1D0/POWIP)
          ELSE
  146       B2RPW=B2RPDV-2D0*LOG(PYR(0))
            ACCIP=(B2RPW/B2RPMX)**RPWIP*EXP(-0.5D0*(B2RPW-B2RPMX))
            IF(ACCIP.LT.PYR(0)) GOTO 146
            OV=EXP(-B2RPW)/PARU(2)
            B=B2RPW**(1D0/POWIP)
          ENDIF  
          VINT(148)=OV/VNT147
          PACC=(1D0-EXP(-MIN(50D0,PIK*OV)))/(PIK*OV)
        ENDIF
        IF(PACC.LT.PYR(0)) GOTO 142
        VINT(139)=B/BAVG
 
      ELSEIF(MMUL.EQ.3) THEN
C...Low-pT or multiple interactions (first semihard interaction):
C...choose xT2 according to dpT2/pT2**2*exp(-(sigma above pT2)/norm)
C...or (MSTP(82)>=2) dpT2/(pT2+pT0**2)**2*exp(-....).
        ISUB=MINT(1)
        VINT(145)=VNT145
        VINT(146)=VNT146
        VINT(147)=VNT147
        IF(MSTP(82).LE.0) THEN
          XT2=0D0
        ELSEIF(MSTP(82).EQ.1) THEN
          XT2=XT2FAC*XT2/(XT2FAC-XT2*LOG(PYR(0)))
C...Use with "Sudakov" for low b values when impact parameter dependence.
        ELSEIF(MSTP(82).EQ.2.OR.MINT(39).EQ.1) THEN
          IF(XT2.LT.1D0.AND.EXP(-XT2FAC*XT2/(VINT(149)*(XT2+
     &    VINT(149)))).GT.PYR(0)) XT2=1D0
          IF(XT2.GE.1D0) THEN
            XT2=(1D0+VINT(149))*XT2FAC/(XT2FAC-(1D0+VINT(149))*LOG(1D0-
     &      PYR(0)*(1D0-EXP(-XT2FAC/(VINT(149)*(1D0+VINT(149)))))))-
     &      VINT(149)
          ELSE
            XT2=-XT2FAC/LOG(EXP(-XT2FAC/(XT2+VINT(149)))+PYR(0)*
     &      (EXP(-XT2FAC/VINT(149))-EXP(-XT2FAC/(XT2+VINT(149)))))-
     &      VINT(149)
          ENDIF
          XT2=MAX(0.01D0*VINT(149),XT2)
C...Use without "Sudakov" for high b values when impact parameter dep.
        ELSE
          XT2=(XC2+VINT(149))*(1D0+VINT(149))/(1D0+VINT(149)-
     &    PYR(0)*(1D0-XC2))-VINT(149)
          XT2=MAX(0.01D0*VINT(149),XT2)
        ENDIF
        VINT(25)=XT2
 
C...Low-pT: choose xT2, tau, y* and cos(theta-hat) fixed.
        IF(MSTP(82).LE.1.AND.XT2.LT.VINT(149)) THEN
          IF(MINT(82).EQ.1) NGEN(0,1)=NGEN(0,1)-MINT(143)
          IF(MINT(82).EQ.1) NGEN(ISUB,1)=NGEN(ISUB,1)-MINT(143)
          ISUB=95
          MINT(1)=ISUB
          VINT(21)=1D-12*VINT(149)
          VINT(22)=0D0
          VINT(23)=0D0
          VINT(25)=1D-12*VINT(149)
 
        ELSE
C...Multiple interactions (first semihard interaction).
C...Choose tau and y*. Calculate cos(theta-hat).
          IF(PYR(0).LE.COEF(ISUB,1)) THEN
            TAUT=(2D0*(1D0+SQRT(1D0-XT2))/XT2-1D0)**PYR(0)
            TAU=XT2*(1D0+TAUT)**2/(4D0*TAUT)
          ELSE
            TAU=XT2*(1D0+TAN(PYR(0)*ATAN(SQRT(1D0/XT2-1D0)))**2)
          ENDIF
          VINT(21)=TAU
          CALL PYKLIM(2)
          RYST=PYR(0)
          MYST=1
          IF(RYST.GT.COEF(ISUB,8)) MYST=2
          IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)) MYST=3
          CALL PYKMAP(2,MYST,PYR(0))
          VINT(23)=SQRT(MAX(0D0,1D0-XT2/TAU))*(-1)**INT(1.5D0+PYR(0))
        ENDIF
        VINT(71)=0.5D0*VINT(1)*SQRT(VINT(25))
 
C...Store results of cross-section calculation.
      ELSEIF(MMUL.EQ.4) THEN
        ISUB=MINT(1)
        VINT(145)=VNT145
        VINT(146)=VNT146
        VINT(147)=VNT147
        XTS=VINT(25)
        IF(ISET(ISUB).EQ.1) XTS=VINT(21)
        IF(ISET(ISUB).EQ.2)
     &  XTS=(4D0*VINT(48)+2D0*VINT(63)+2D0*VINT(64))/VINT(2)
        IF(ISET(ISUB).GE.3.AND.ISET(ISUB).LE.5) XTS=VINT(26)
        RBIN=MAX(0.000001D0,MIN(0.999999D0,XTS*(1D0+VINT(149))/
     &  (XTS+VINT(149))))
        IRBIN=INT(1D0+20D0*RBIN)
        IF(ISUB.EQ.96.AND.MSTP(171).EQ.0) THEN
          NMUL(IRBIN)=NMUL(IRBIN)+1
          SIGM(IRBIN)=SIGM(IRBIN)+VINT(153)
        ENDIF
 
C...Choose impact parameter if not already done.
      ELSEIF(MMUL.EQ.5) THEN
        ISUB=MINT(1)
        VINT(145)=VNT145
        VINT(146)=VNT146
        VINT(147)=VNT147
  150   IF(MINT(39).GT.0) THEN
        ELSEIF(MSTP(82).EQ.3) THEN
          EXPB2=PYR(0)
          B2=-LOG(PYR(0))
          VINT(148)=EXPB2/(PARU(2)*VNT147)
          VINT(139)=SQRT(B2)/BAVG
        ELSEIF(MSTP(82).EQ.4) THEN
          RTYPE=PYR(0)
          IF(RTYPE.LT.P83A) THEN
            B2=-LOG(PYR(0))
          ELSEIF(RTYPE.LT.P83A+P83B) THEN
            B2=-LOG(PYR(0))/CQ2R
          ELSE
            B2=-LOG(PYR(0))/CQ2I
          ENDIF
          VINT(148)=(P83A*EXP(-MIN(50D0,B2))+
     &    P83B*CQ2R*EXP(-MIN(50D0,B2*CQ2R))+
     &    P83C*CQ2I*EXP(-MIN(50D0,B2*CQ2I)))/(PARU(2)*VNT147)
          VINT(139)=SQRT(B2)/BAVG
        ELSEIF(PARP(83).GE.1.999D0) THEN
          POWIP=MAX(2D0,PARP(83))
          RPWIP=2D0/POWIP-1D0
          PROB1=POWIP/(2D0*EXP(-1D0)+POWIP)
  160     IF(PYR(0).LT.PROB1) THEN
            B2RPW=PYR(0)**(0.5D0*POWIP)
            ACCIP=EXP(-B2RPW)
          ELSE
            B2RPW=1D0-LOG(PYR(0))
            ACCIP=B2RPW**RPWIP
          ENDIF
          IF(ACCIP.LT.PYR(0)) GOTO 160
          VINT(148)=EXP(-B2RPW)/(PARU(2)*VNT147)
          VINT(139)=B2RPW**(1D0/POWIP)/BAVG
        ELSE
          POWIP=MAX(0.4D0,PARP(83))
          RPWIP=2D0/POWIP-1D0
          PROB1=RPWIP/(RPWIP+2D0**RPWIP*EXP(-RPWIP))
  170     IF(PYR(0).LT.PROB1) THEN
            B2RPW=2D0*RPWIP*PYR(0)
            ACCIP=(B2RPW/RPWIP)**RPWIP*EXP(RPWIP-B2RPW)
          ELSE
            B2RPW=2D0*(RPWIP-LOG(PYR(0)))
            ACCIP=(0.5D0*B2RPW/RPWIP)**RPWIP*EXP(RPWIP-0.5D0*B2RPW)
          ENDIF
          IF(ACCIP.LT .PYR(0)) GOTO 170
          VINT(148)=EXP(-B2RPW)/(PARU(2)*VNT147)
          VINT(139)=B2RPW**(1D0/POWIP)/BAVG
        ENDIF
 
C...Multiple interactions (variable impact parameter) : reject with
C...probability exp(-overlap*cross-section above pT/normalization).
C...Does not apply to low-b region, where "Sudakov" already included.
        VINT(150)=1D0 
        IF(MINT(39).NE.1) THEN
          RNCOR=(IRBIN-20D0*RBIN)*NMUL(IRBIN)
          SIGCOR=(IRBIN-20D0*RBIN)*SIGM(IRBIN)
          DO 180 IBIN=IRBIN+1,20
            RNCOR=RNCOR+NMUL(IBIN)
            SIGCOR=SIGCOR+SIGM(IBIN)
  180     CONTINUE
          SIGABV=(SIGCOR/RNCOR)*VINT(149)*(1D0-XTS)/(XTS+VINT(149))
          IF(MSTP(171).EQ.1) SIGABV=SIGABV*VINT(2)/VINT(289)
          VINT(150)=EXP(-MIN(50D0,VNT146*VINT(148)*
     &    SIGABV/MAX(1D-10,SIGT(0,0,5))))
        ENDIF
        IF(MSTP(86).EQ.3.OR.(MSTP(86).EQ.2.AND.ISUB.NE.11.AND.
     &  ISUB.NE.12.AND.ISUB.NE.13.AND.ISUB.NE.28.AND.ISUB.NE.53
     &  .AND.ISUB.NE.68.AND.ISUB.NE.95.AND.ISUB.NE.96)) THEN
          IF(VINT(150).LT.PYR(0)) GOTO 150
          VINT(150)=1D0
        ENDIF
 
C...Generate additional multiple semihard interactions.
      ELSEIF(MMUL.EQ.6) THEN
 
C...Save data for hardest initeraction, to be restored.
        ISUBSV=MINT(1)
        VINT(145)=VNT145
        VINT(146)=VNT146
        VINT(147)=VNT147
        M13SV=MINT(13)
        M14SV=MINT(14)
        M15SV=MINT(15)
        M16SV=MINT(16)
        M21SV=MINT(21)
        M22SV=MINT(22)
        DO 190 J=11,80
          VINTSV(J)=VINT(J)
  190   CONTINUE
        V141SV=VINT(141)
        V142SV=VINT(142)
 
C...Store data on hardest interaction.
        XMI(1,1)=VINT(141)
        XMI(2,1)=VINT(142)
        PT2MI(1)=VINT(54)
        IMISEP(0)=MINT(84)
        IMISEP(1)=N
 
C...Change process to generate; sum of x values so far.
        ISUB=96
        MINT(1)=96
        VINT(143)=1D0-VINT(141)
        VINT(144)=1D0-VINT(142)
        VINT(151)=0D0
        VINT(152)=0D0
 
C...Initialize factors for PDF reshaping.
        DO 230 JS=1,2
          KFBEAM=MINT(10+JS)
          KFABM=IABS(KFBEAM)
          KFSBM=ISIGN(1,KFBEAM)
 
C...Zero flavour content of incoming beam particle.
          KFIVAL(JS,1)=0
          KFIVAL(JS,2)=0
          KFIVAL(JS,3)=0
C...Flavour content of baryon.
          IF(KFABM.GT.1000) THEN
            KFIVAL(JS,1)=KFSBM*MOD(KFABM/1000,10)
            KFIVAL(JS,2)=KFSBM*MOD(KFABM/100,10)
            KFIVAL(JS,3)=KFSBM*MOD(KFABM/10,10)
C...Flavour content of pi+-, K+-.
          ELSEIF(KFABM.EQ.211) THEN
            KFIVAL(JS,1)=KFSBM*2
            KFIVAL(JS,2)=-KFSBM
          ELSEIF(KFABM.EQ.321) THEN
            KFIVAL(JS,1)=-KFSBM*3
            KFIVAL(JS,2)=KFSBM*2
C...Flavour content of pi0, gamma, K0S, K0L not defined yet.
          ENDIF
 
C...Zero initial valence and companion content.
          DO 200 IFL=-6,6
            NVC(JS,IFL)=0
  200     CONTINUE
 
C...Initiate listing of all incoming partons from two sides.
          NMI(JS)=0
          DO 210 I=MINT(84)+1,N
            IF(K(I,3).EQ.MINT(83)+2+JS) THEN
              IMI(JS,1,1)=I
              IMI(JS,1,2)=0
            ENDIF
  210     CONTINUE
 
C...Decide whether quarks in hard scattering were valence or sea.
          IFL=K(IMI(JS,1,1),2)
          IF (IABS(IFL).GT.6) GOTO 230
 
C...Get PDFs at X and Q2 of the parton shower initiator for the
C...hard scattering.
          X=VINT(140+JS)
          IF(MSTP(61).GE.1) THEN
            Q2=PARP(62)**2
          ELSE
            Q2=VINT(54)
          ENDIF
C...Note: XPSVC = x*pdf.
          MINT(30)=JS
          CALL PYPDFU(KFBEAM,X,Q2,XPQ)
          SEA=XPSVC(IFL,-1)
          VAL=XPSVC(IFL,0)
 
C...Decide (Extra factor x cancels in the division).
          RVCS=PYR(0)*(SEA+VAL)
          IVNOW=1
  220     IF (RVCS.LE.VAL.AND.IVNOW.GE.1) THEN
C...Safety check that valence present; pi0/gamma/K0S/K0L special cases.
            IVNOW=0
            IF(KFIVAL(JS,1).EQ.IFL) IVNOW=IVNOW+1
            IF(KFIVAL(JS,2).EQ.IFL) IVNOW=IVNOW+1
            IF(KFIVAL(JS,3).EQ.IFL) IVNOW=IVNOW+1
            IF(KFIVAL(JS,1).EQ.0) THEN
              IF(KFBEAM.EQ.111.AND.IABS(IFL).LE.2) IVNOW=1
              IF(KFBEAM.EQ.22.AND.IABS(IFL).LE.5) IVNOW=1
              IF((KFBEAM.EQ.130.OR.KFBEAM.EQ.310).AND.
     &        (IABS(IFL).EQ.1.OR.IABS(IFL).EQ.3)) IVNOW=1
            ENDIF
            IF(IVNOW.EQ.0) GOTO 220
C...Mark valence.
            IMI(JS,1,2)=0
C...Sets valence content of gamma, pi0, K0S, K0L if not done.
            IF(KFIVAL(JS,1).EQ.0) THEN
              IF(KFBEAM.EQ.111.OR.KFBEAM.EQ.22) THEN
                KFIVAL(JS,1)=IFL
                KFIVAL(JS,2)=-IFL
              ELSEIF(KFBEAM.EQ.130.OR.KFBEAM.EQ.310) THEN
                KFIVAL(JS,1)=IFL
                IF(IABS(IFL).EQ.1) KFIVAL(JS,2)=ISIGN(3,-IFL)
                IF(IABS(IFL).NE.1) KFIVAL(JS,2)=ISIGN(1,-IFL)
              ENDIF
            ENDIF
 
C...If sea, add opposite sign companion parton. Store X and I.
          ELSE
            NVC(JS,-IFL)=NVC(JS,-IFL)+1
            XASSOC(JS,-IFL,NVC(JS,-IFL))=X
C...Set pointer to companion
            IMI(JS,1,2)=-NVC(JS,-IFL)
          ENDIF
  230   CONTINUE
 
C...Update counter number of multiple interactions.
        NMI(1)=1
        NMI(2)=1
 
C...Set up starting values for iteration in xT2.
        IF(MSTP(86).EQ.3.OR.(MSTP(86).EQ.2.AND.ISUBSV.NE.11.AND.
     &  ISUBSV.NE.12.AND.ISUBSV.NE.13.AND.ISUBSV.NE.28.AND.
     &  ISUBSV.NE.53.AND.ISUBSV.NE.68.AND.ISUBSV.NE.95.AND.
     &  ISUBSV.NE.96)) THEN
          XT2=(1D0-VINT(141))*(1D0-VINT(142))
        ELSE
          XT2=VINT(25)
          IF(ISET(ISUBSV).EQ.1) XT2=VINT(21)
          IF(ISET(ISUBSV).EQ.2)
     &    XT2=(4D0*VINT(48)+2D0*VINT(63)+2D0*VINT(64))/VINT(2)
          IF(ISET(ISUBSV).GE.3.AND.ISET(ISUBSV).LE.5) XT2=VINT(26)
        ENDIF
        IF(MSTP(82).LE.1) THEN
          SIGRAT=XSEC(ISUB,1)/MAX(1D-10,VINT(315)*VINT(316)*SIGT(0,0,5))
          IF(MINT(141).NE.0.OR.MINT(142).NE.0) SIGRAT=SIGRAT*
     &    VINT(317)/(VINT(318)*VINT(320))
          XT2FAC=SIGRAT*VINT(149)/(1D0-VINT(149))
        ELSE
          XT2FAC=VNT146*VINT(148)*XSEC(ISUB,1)/
     &    MAX(1D-10,SIGT(0,0,5))*VINT(149)*(1D0+VINT(149))
        ENDIF
        VINT(63)=0D0
        VINT(64)=0D0
 
C...Iterate downwards in xT2.
  240   IF((MINT(35).EQ.2.AND.MSTP(81).EQ.10).OR.ISUBSV.EQ.95) THEN
          XT2=0D0
          GOTO 440
        ELSEIF(MSTP(82).LE.1) THEN
          XT2=XT2FAC*XT2/(XT2FAC-XT2*LOG(PYR(0)))
          IF(XT2.LT.VINT(149)) GOTO 440
        ELSE
          IF(XT2.LE.0.01001D0*VINT(149)) GOTO 440
          XT2=XT2FAC*(XT2+VINT(149))/(XT2FAC-(XT2+VINT(149))*
     &    LOG(PYR(0)))-VINT(149)
          IF(XT2.LE.0D0) GOTO 440
          XT2=MAX(0.01D0*VINT(149),XT2)
        ENDIF
        VINT(25)=XT2
 
C...Choose tau and y*. Calculate cos(theta-hat).
        IF(PYR(0).LE.COEF(ISUB,1)) THEN
          TAUT=(2D0*(1D0+SQRT(1D0-XT2))/XT2-1D0)**PYR(0)
          TAU=XT2*(1D0+TAUT)**2/(4D0*TAUT)
        ELSE
          TAU=XT2*(1D0+TAN(PYR(0)*ATAN(SQRT(1D0/XT2-1D0)))**2)
        ENDIF
        VINT(21)=TAU
C...New: require shat > 1.
        IF(TAU*VINT(2).LT.1D0) GOTO 240
        CALL PYKLIM(2)
        RYST=PYR(0)
        MYST=1
        IF(RYST.GT.COEF(ISUB,8)) MYST=2
        IF(RYST.GT.COEF(ISUB,8)+COEF(ISUB,9)) MYST=3
        CALL PYKMAP(2,MYST,PYR(0))
        VINT(23)=SQRT(MAX(0D0,1D0-XT2/TAU))*(-1)**INT(1.5D0+PYR(0))
 
C...Check that x not used up. Accept or reject kinematical variables.
        X1M=SQRT(TAU)*EXP(VINT(22))
        X2M=SQRT(TAU)*EXP(-VINT(22))
        IF(VINT(143)-X1M.LT.0.01D0.OR.VINT(144)-X2M.LT.0.01D0) GOTO 240
        VINT(71)=0.5D0*VINT(1)*SQRT(XT2)
        CALL PYSIGH(NCHN,SIGS)
        IF(MINT(141).NE.0.OR.MINT(142).NE.0) SIGS=SIGS*VINT(320)
        IF(SIGS.LT.XSEC(ISUB,1)*PYR(0)) GOTO 240
        IF(MINT(141).NE.0.OR.MINT(142).NE.0) SIGS=SIGS/VINT(320)
 
C...Reset K, P and V vectors.
        DO 260 I=N+1,N+4
          DO 250 J=1,5
            K(I,J)=0
            P(I,J)=0D0
            V(I,J)=0D0
  250     CONTINUE
  260   CONTINUE
        PT=0.5D0*VINT(1)*SQRT(XT2)
 
C...Choose flavour of reacting partons (and subprocess).
        RSIGS=SIGS*PYR(0)
        DO 270 ICHN=1,NCHN
          KFL1=ISIG(ICHN,1)
          KFL2=ISIG(ICHN,2)
          ICONMI=ISIG(ICHN,3)
          RSIGS=RSIGS-SIGH(ICHN)
          IF(RSIGS.LE.0D0) GOTO 280
  270   CONTINUE
 
C...Reassign to appropriate process codes.
  280   ISUBMI=ICONMI/10
        ICONMI=MOD(ICONMI,10)
 
C...Choose new quark flavour for annihilation graphs
        IF(ISUBMI.EQ.12.OR.ISUBMI.EQ.53) THEN
          SH=TAU*VINT(2)
          CALL PYWIDT(21,SH,WDTP,WDTE)
  290     RKFL=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))*PYR(0)
          DO 300 I=1,MDCY(21,3)
            KFLF=KFDP(I+MDCY(21,2)-1,1)
            RKFL=RKFL-(WDTE(I,1)+WDTE(I,2)+WDTE(I,4))
            IF(RKFL.LE.0D0) GOTO 310
  300     CONTINUE
  310     IF(ISUBMI.EQ.53.AND.ICONMI.LE.2) THEN
            IF(KFLF.GE.4) GOTO 290
          ELSEIF(ISUBMI.EQ.53.AND.ICONMI.LE.4) THEN
            KFLF=4
            ICONMI=ICONMI-2
          ELSEIF(ISUBMI.EQ.53) THEN
            KFLF=5
            ICONMI=ICONMI-4
          ENDIF
        ENDIF
 
C...Final state flavours and colour flow: default values
        JS=1
        KFL3=KFL1
        KFL4=KFL2
        KCC=20
        KCS=ISIGN(1,KFL1)
 
        IF(ISUBMI.EQ.11) THEN
C...f + f' -> f + f' (g exchange); th = (p(f)-p(f))**2
          KCC=ICONMI
          IF(KFL1*KFL2.LT.0) KCC=KCC+2
 
        ELSEIF(ISUBMI.EQ.12) THEN
C...f + fbar -> f' + fbar'; th = (p(f)-p(f'))**2
          KFL3=ISIGN(KFLF,KFL1)
          KFL4=-KFL3
          KCC=4
 
        ELSEIF(ISUBMI.EQ.13) THEN
C...f + fbar -> g + g; th arbitrary
          KFL3=21
          KFL4=21
          KCC=ICONMI+4
 
        ELSEIF(ISUBMI.EQ.28) THEN
C...f + g -> f + g; th = (p(f)-p(f))**2
          IF(KFL1.EQ.21) JS=2
          KCC=ICONMI+6
          IF(KFL1.EQ.21) KCC=KCC+2
          IF(KFL1.NE.21) KCS=ISIGN(1,KFL1)
          IF(KFL2.NE.21) KCS=ISIGN(1,KFL2)
 
        ELSEIF(ISUBMI.EQ.53) THEN
C...g + g -> f + fbar; th arbitrary
          KCS=(-1)**INT(1.5D0+PYR(0))
          KFL3=ISIGN(KFLF,KCS)
          KFL4=-KFL3
          KCC=ICONMI+10
 
        ELSEIF(ISUBMI.EQ.68) THEN
C...g + g -> g + g; th arbitrary
          KCC=ICONMI+12
          KCS=(-1)**INT(1.5D0+PYR(0))
        ENDIF
 
C...Store flavours of scattering.
        MINT(13)=KFL1
        MINT(14)=KFL2
        MINT(15)=KFL1
        MINT(16)=KFL2
        MINT(21)=KFL3
        MINT(22)=KFL4
 
C...Set flavours and mothers of scattering partons.
        K(N+1,1)=14
        K(N+2,1)=14
        K(N+3,1)=3
        K(N+4,1)=3
        K(N+1,2)=KFL1
        K(N+2,2)=KFL2
        K(N+3,2)=KFL3
        K(N+4,2)=KFL4
        K(N+1,3)=MINT(83)+1
        K(N+2,3)=MINT(83)+2
        K(N+3,3)=N+1
        K(N+4,3)=N+2
 
C...Store colour connection indices.
        DO 320 J=1,2
          JC=J
          IF(KCS.EQ.-1) JC=3-J
          IF(ICOL(KCC,1,JC).NE.0) K(N+1,J+3)=N+ICOL(KCC,1,JC)
          IF(ICOL(KCC,2,JC).NE.0) K(N+2,J+3)=N+ICOL(KCC,2,JC)
          IF(ICOL(KCC,3,JC).NE.0) K(N+3,J+3)=MSTU(5)*(N+ICOL(KCC,3,JC))
          IF(ICOL(KCC,4,JC).NE.0) K(N+4,J+3)=MSTU(5)*(N+ICOL(KCC,4,JC))
  320   CONTINUE
 
C...Store incoming and outgoing partons in their CM-frame.
        SHR=SQRT(TAU)*VINT(1)
        P(N+1,3)=0.5D0*SHR
        P(N+1,4)=0.5D0*SHR
        P(N+2,3)=-0.5D0*SHR
        P(N+2,4)=0.5D0*SHR
        P(N+3,5)=PYMASS(K(N+3,2))
        P(N+4,5)=PYMASS(K(N+4,2))
        IF(P(N+3,5)+P(N+4,5).GE.SHR) GOTO 240
        P(N+3,4)=0.5D0*(SHR+(P(N+3,5)**2-P(N+4,5)**2)/SHR)
        P(N+3,3)=SQRT(MAX(0D0,P(N+3,4)**2-P(N+3,5)**2))
        P(N+4,4)=SHR-P(N+3,4)
        P(N+4,3)=-P(N+3,3)
 
C...Rotate outgoing partons using cos(theta)=(th-uh)/lam(sh,sqm3,sqm4)
        PHI=PARU(2)*PYR(0)
        CALL PYROBO(N+3,N+4,ACOS(VINT(23)),PHI,0D0,0D0,0D0)
 
C...Set up default values before showers.
        MINT(31)=MINT(31)+1
        IPU1=N+1
        IPU2=N+2
        IPU3=N+3
        IPU4=N+4
        VINT(141)=VINT(41)
        VINT(142)=VINT(42)
        N=N+4
 
C...Showering of initial state partons (optional).
C...Note: no showering of final state partons here; it comes later.
        IF(MSTP(84).GE.1.AND.MSTP(61).GE.1) THEN
          MINT(51)=0
          ALAMSV=PARJ(81)
          PARJ(81)=PARP(72)
          NSAV=N
          DO 340 I=1,4
            DO 330 J=1,5
              KSAV(I,J)=K(N-4+I,J)
              PSAV(I,J)=P(N-4+I,J)
  330       CONTINUE
  340     CONTINUE
          CALL PYSSPA(IPU1,IPU2)
          PARJ(81)=ALAMSV
C...If shower failed then restore to situation before shower.
          IF(MINT(51).GE.1) THEN
            N=NSAV
            DO 360 I=1,4
              DO 350 J=1,5
                K(N-4+I,J)=KSAV(I,J)
                P(N-4+I,J)=PSAV(I,J)
  350         CONTINUE
  360       CONTINUE
            IPU1=N-3
            IPU2=N-2
            VINT(141)=VINT(41)
            VINT(142)=VINT(42)
          ENDIF
        ENDIF
 
C...Keep track of loose colour ends and information on scattering.
  370   IMI(1,MINT(31),1)=IPU1
        IMI(2,MINT(31),1)=IPU2
        IMI(1,MINT(31),2)=0
        IMI(2,MINT(31),2)=0
        XMI(1,MINT(31))=VINT(141)
        XMI(2,MINT(31))=VINT(142)
        PT2MI(MINT(31))=VINT(54)
        IMISEP(MINT(31))=N
 
C...Decide whether quarks in last scattering were valence, companion or
C...sea.
        DO 430 JS=1,2
          KFBEAM=MINT(10+JS)
          KFSBM=ISIGN(1,MINT(10+JS))
          IFL=K(IMI(JS,MINT(31),1),2)
          IMI(JS,MINT(31),2)=0
          IF (IABS(IFL).GT.6) GOTO 430
 
C...Get PDFs at X and Q2 of the parton shower initiator for the
C...last scattering. At this point VINT(143:144) do not yet
C...include the scattered x values VINT(141:142).
          X=VINT(140+JS)/VINT(142+JS)
          IF(MSTP(84).GE.1.AND.MSTP(61).GE.1) THEN
            Q2=PARP(62)**2
          ELSE
            Q2=VINT(54)
          ENDIF
C...Note: XPSVC = x*pdf.
          MINT(30)=JS
          CALL PYPDFU(KFBEAM,X,Q2,XPQ)
          SEA=XPSVC(IFL,-1)
          VAL=XPSVC(IFL,0)
          CMP=0D0
          DO 380 IVC=1,NVC(JS,IFL)
            CMP=CMP+XPSVC(IFL,IVC)
  380     CONTINUE
 
C...Decide (Extra factor x cancels in the dvision).
          RVCS=PYR(0)*(SEA+VAL+CMP)
          IVNOW=1
  390     IF (RVCS.LE.VAL.AND.IVNOW.GE.1) THEN
C...Safety check that valence present; pi0/gamma/K0S/K0L special cases.
            IVNOW=0
            IF(KFIVAL(JS,1).EQ.IFL) IVNOW=IVNOW+1
            IF(KFIVAL(JS,2).EQ.IFL) IVNOW=IVNOW+1
            IF(KFIVAL(JS,3).EQ.IFL) IVNOW=IVNOW+1
            IF(KFIVAL(JS,1).EQ.0) THEN
              IF(KFBEAM.EQ.111.AND.IABS(IFL).LE.2) IVNOW=1
              IF(KFBEAM.EQ.22.AND.IABS(IFL).LE.5) IVNOW=1
              IF((KFBEAM.EQ.130.OR.KFBEAM.EQ.310).AND.
     &        (IABS(IFL).EQ.1.OR.IABS(IFL).EQ.3)) IVNOW=1
            ELSE
              DO 400 I1=1,NMI(JS)
                IF (K(IMI(JS,I1,1),2).EQ.IFL.AND.IMI(JS,I1,2).EQ.0)
     &            IVNOW=IVNOW-1
  400         CONTINUE
            ENDIF
            IF(IVNOW.EQ.0) GOTO 390
C...Mark valence.
            IMI(JS,MINT(31),2)=0
C...Sets valence content of gamma, pi0, K0S, K0L if not done.
            IF(KFIVAL(JS,1).EQ.0) THEN
              IF(KFBEAM.EQ.111.OR.KFBEAM.EQ.22) THEN
                KFIVAL(JS,1)=IFL
                KFIVAL(JS,2)=-IFL
              ELSEIF(KFBEAM.EQ.130.OR.KFBEAM.EQ.310) THEN
                KFIVAL(JS,1)=IFL
                IF(IABS(IFL).EQ.1) KFIVAL(JS,2)=ISIGN(3,-IFL)
                IF(IABS(IFL).NE.1) KFIVAL(JS,2)=ISIGN(1,-IFL)
              ENDIF
            ENDIF
 
          ELSEIF (RVCS.LE.VAL+SEA.OR.NVC(JS,IFL).EQ.0) THEN
C...If sea, add opposite sign companion parton. Store X and I.
            NVC(JS,-IFL)=NVC(JS,-IFL)+1
            XASSOC(JS,-IFL,NVC(JS,-IFL))=X
C...Set pointer to companion
            IMI(JS,MINT(31),2)=-NVC(JS,-IFL)
          ELSE
C...If companion, decide which one.
            CMPSUM=VAL+SEA
            ISEL=0
  410       ISEL=ISEL+1
            CMPSUM=CMPSUM+XPSVC(IFL,ISEL)
            IF (RVCS.GT.CMPSUM.AND.ISEL.LT.NVC(JS,IFL)) GOTO 410
C...Find original sea (anti-)quark:
            IASSOC=0
            DO 420 I1=1,NMI(JS)
              IF (K(IMI(JS,I1,1),2).NE.-IFL) GOTO 420
              IF (-IMI(JS,I1,2).EQ.ISEL) THEN
                IMI(JS,MINT(31),2)=IMI(JS,I1,1)
                IMI(JS,I1,2)=IMI(JS,MINT(31),1)
              ENDIF
  420       CONTINUE
C...Change X to what associated companion had, so that the correct
C...amount of momentum can be subtracted from the companion sum below.
            X=XASSOC(JS,IFL,ISEL)
C...Mark companion read.
            XASSOC(JS,IFL,ISEL)=0D0
          ENDIF
 430    CONTINUE
 
C...Global statistics.
        MINT(351)=MINT(351)+1
        VINT(351)=VINT(351)+PT
        IF (MINT(351).EQ.1) VINT(356)=PT
 
C...Update remaining energy and other counters.
        IF(N.GT.MSTU(4)-MSTU(32)-10) THEN
          CALL PYERRM(11,'(PYMIGN:) no more memory left in PYJETS')
          MINT(51)=1
          RETURN
        ENDIF
        NMI(1)=NMI(1)+1
        NMI(2)=NMI(2)+1
        VINT(151)=VINT(151)+VINT(41)
        VINT(152)=VINT(152)+VINT(42)
        VINT(143)=VINT(143)-VINT(141)
        VINT(144)=VINT(144)-VINT(142)
 
C...Iterate, with more interactions allowed.
        IF(MINT(31).LT.240) GOTO 240
 440    CONTINUE
 
C...Restore saved quantities for hardest interaction.
        MINT(1)=ISUBSV
        MINT(13)=M13SV
        MINT(14)=M14SV
        MINT(15)=M15SV
        MINT(16)=M16SV
        MINT(21)=M21SV
        MINT(22)=M22SV
        DO 450 J=11,80
          VINT(J)=VINTSV(J)
  450   CONTINUE
        VINT(141)=V141SV
        VINT(142)=V142SV
 
      ENDIF
 
C...Format statements for printout.
 5000 FORMAT(/1X,'****** PYMIGN: initialization of multiple inter',
     &'actions for MSTP(82) =',I2,' ******')
 5100 FORMAT(8X,'pT0 =',F5.2,' GeV gives sigma(parton-parton) =',1P,
     &D9.2,' mb: rejected')
 5200 FORMAT(8X,'pT0 =',F5.2,' GeV gives sigma(parton-parton) =',1P,
     &D9.2,' mb: accepted')
 
      RETURN
      END
