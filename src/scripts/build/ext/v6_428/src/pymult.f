 
C*********************************************************************
 
C...PYMULT
C...Initializes treatment of multiple interactions, selects kinematics
C...of hardest interaction if low-pT physics included in run, and
C...generates all non-hardest interactions.
 
      SUBROUTINE PYMULT(MMUL)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYSUBS/,/PYPARS/,/PYINT1/,
     &/PYINT2/,/PYINT3/,/PYINT5/,/PYINT7/
C...Local arrays and saved variables.
      DIMENSION NMUL(20),SIGM(20),KSTR(500,2),VINTSV(80)
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

C...Find relative weight for low and high impact parameter.
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
          VINT(21)=0.01D0*VINT(149)
          VINT(22)=0D0
          VINT(23)=0D0
          VINT(25)=0.01D0*VINT(149)
 
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
        ISUBSV=MINT(1)
        VINT(145)=VNT145
        VINT(146)=VNT146
        VINT(147)=VNT147
        DO 190 J=11,80
          VINTSV(J)=VINT(J)
  190   CONTINUE
        ISUB=96
        MINT(1)=96
        VINT(151)=0D0
        VINT(152)=0D0
 
C...Reconstruct strings in hard scattering.
        NMAX=MINT(84)+4
        IF(ISET(ISUBSV).EQ.1) NMAX=MINT(84)+2
        IF(ISET(ISUBSV).EQ.11) NMAX=MINT(84)+2+MINT(3)
        NSTR=0
        DO 210 I=MINT(84)+1,NMAX
          KCS=KCHG(PYCOMP(K(I,2)),2)*ISIGN(1,K(I,2))
          IF(KCS.EQ.0) GOTO 210
          DO 200 J=1,4
            IF(KCS.EQ.1.AND.(J.EQ.2.OR.J.EQ.4)) GOTO 200
            IF(KCS.EQ.-1.AND.(J.EQ.1.OR.J.EQ.3)) GOTO 200
            IF(J.LE.2) THEN
              IST=MOD(K(I,J+3)/MSTU(5),MSTU(5))
            ELSE
              IST=MOD(K(I,J+1),MSTU(5))
            ENDIF
            IF(IST.LT.MINT(84).OR.IST.GT.I) GOTO 200
            IF(KCHG(PYCOMP(K(IST,2)),2).EQ.0) GOTO 200
            NSTR=NSTR+1
            IF(J.EQ.1.OR.J.EQ.4) THEN
              KSTR(NSTR,1)=I
              KSTR(NSTR,2)=IST
            ELSE
              KSTR(NSTR,1)=IST
              KSTR(NSTR,2)=I
            ENDIF
  200     CONTINUE
  210   CONTINUE
 
C...Set up starting values for iteration in xT2.
        XT2=4D0*VINT(62)/VINT(2)
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
        VINT(143)=1D0-VINT(141)
        VINT(144)=1D0-VINT(142)
 
C...Iterate downwards in xT2.
  220   IF(MSTP(82).LE.1) THEN
          XT2=XT2FAC*XT2/(XT2FAC-XT2*LOG(PYR(0)))
          IF(XT2.LT.VINT(149)) GOTO 270
        ELSE
          IF(XT2.LE.0.01001D0*VINT(149)) GOTO 270
          XT2=XT2FAC*(XT2+VINT(149))/(XT2FAC-(XT2+VINT(149))*
     &    LOG(PYR(0)))-VINT(149)
          IF(XT2.LE.0D0) GOTO 270
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
        IF(VINT(143)-X1M.LT.0.01D0.OR.VINT(144)-X2M.LT.0.01D0) GOTO 220
        VINT(71)=0.5D0*VINT(1)*SQRT(XT2)
        CALL PYSIGH(NCHN,SIGS)
        IF(MINT(141).NE.0.OR.MINT(142).NE.0) SIGS=SIGS*VINT(320)
        IF(SIGS.LT.XSEC(ISUB,1)*PYR(0)) GOTO 220
 
C...Reset K, P and V vectors. Select some variables.
        DO 240 I=N+1,N+2
          DO 230 J=1,5
            K(I,J)=0
            P(I,J)=0D0
            V(I,J)=0D0
  230     CONTINUE
  240   CONTINUE
        RFLAV=PYR(0)
        PT=0.5D0*VINT(1)*SQRT(XT2)
        PHI=PARU(2)*PYR(0)
        CTH=VINT(23)
 
C...Add first parton to event record.
        K(N+1,1)=3
        K(N+1,2)=21
        IF(RFLAV.GE.MAX(PARP(85),PARP(86))) K(N+1,2)=
     &  1+INT((2D0+PARJ(2))*PYR(0))
        P(N+1,1)=PT*COS(PHI)
        P(N+1,2)=PT*SIN(PHI)
        P(N+1,3)=0.25D0*VINT(1)*(VINT(41)*(1D0+CTH)-VINT(42)*(1D0-CTH))
        P(N+1,4)=0.25D0*VINT(1)*(VINT(41)*(1D0+CTH)+VINT(42)*(1D0-CTH))
        P(N+1,5)=0D0
 
C...Add second parton to event record.
        K(N+2,1)=3
        K(N+2,2)=21
        IF(K(N+1,2).NE.21) K(N+2,2)=-K(N+1,2)
        P(N+2,1)=-P(N+1,1)
        P(N+2,2)=-P(N+1,2)
        P(N+2,3)=0.25D0*VINT(1)*(VINT(41)*(1D0-CTH)-VINT(42)*(1D0+CTH))
        P(N+2,4)=0.25D0*VINT(1)*(VINT(41)*(1D0-CTH)+VINT(42)*(1D0+CTH))
        P(N+2,5)=0D0
 
        IF(RFLAV.LT.PARP(85).AND.NSTR.GE.1) THEN
C....Choose relevant string pieces to place gluons on.
          DO 260 I=N+1,N+2
            DMIN=1D8
            DO 250 ISTR=1,NSTR
              I1=KSTR(ISTR,1)
              I2=KSTR(ISTR,2)
              DIST=(P(I,4)*P(I1,4)-P(I,1)*P(I1,1)-P(I,2)*P(I1,2)-
     &        P(I,3)*P(I1,3))*(P(I,4)*P(I2,4)-P(I,1)*P(I2,1)-
     &        P(I,2)*P(I2,2)-P(I,3)*P(I2,3))/MAX(1D0,P(I1,4)*P(I2,4)-
     &        P(I1,1)*P(I2,1)-P(I1,2)*P(I2,2)-P(I1,3)*P(I2,3))
              IF(ISTR.EQ.1.OR.DIST.LT.DMIN) THEN
                DMIN=DIST
                IST1=I1
                IST2=I2
                ISTM=ISTR
              ENDIF
  250       CONTINUE
 
C....Colour flow adjustments, new string pieces.
            IF(K(IST1,4)/MSTU(5).EQ.IST2) K(IST1,4)=MSTU(5)*I+
     &      MOD(K(IST1,4),MSTU(5))
            IF(MOD(K(IST1,5),MSTU(5)).EQ.IST2) K(IST1,5)=
     &      MSTU(5)*(K(IST1,5)/MSTU(5))+I
            K(I,5)=MSTU(5)*IST1
            K(I,4)=MSTU(5)*IST2
            IF(K(IST2,5)/MSTU(5).EQ.IST1) K(IST2,5)=MSTU(5)*I+
     &      MOD(K(IST2,5),MSTU(5))
            IF(MOD(K(IST2,4),MSTU(5)).EQ.IST1) K(IST2,4)=
     &      MSTU(5)*(K(IST2,4)/MSTU(5))+I
            KSTR(ISTM,2)=I
            KSTR(NSTR+1,1)=I
            KSTR(NSTR+1,2)=IST2
            NSTR=NSTR+1
  260     CONTINUE
 
C...String drawing and colour flow for gluon loop.
        ELSEIF(K(N+1,2).EQ.21) THEN
          K(N+1,4)=MSTU(5)*(N+2)
          K(N+1,5)=MSTU(5)*(N+2)
          K(N+2,4)=MSTU(5)*(N+1)
          K(N+2,5)=MSTU(5)*(N+1)
          KSTR(NSTR+1,1)=N+1
          KSTR(NSTR+1,2)=N+2
          KSTR(NSTR+2,1)=N+2
          KSTR(NSTR+2,2)=N+1
          NSTR=NSTR+2
 
C...String drawing and colour flow for qqbar pair.
        ELSE
          K(N+1,4)=MSTU(5)*(N+2)
          K(N+2,5)=MSTU(5)*(N+1)
          KSTR(NSTR+1,1)=N+1
          KSTR(NSTR+1,2)=N+2
          NSTR=NSTR+1
        ENDIF
 
C...Global statistics.
        MINT(351)=MINT(351)+1
        VINT(351)=VINT(351)+PT
        IF (MINT(351).EQ.1) VINT(356)=PT
 
C...Update remaining energy; iterate.
        N=N+2
        IF(N.GT.MSTU(4)-MSTU(32)-10) THEN
          CALL PYERRM(11,'(PYMULT:) no more memory left in PYJETS')
          MINT(51)=1
          RETURN
        ENDIF
        MINT(31)=MINT(31)+1
        VINT(151)=VINT(151)+VINT(41)
        VINT(152)=VINT(152)+VINT(42)
        VINT(143)=VINT(143)-VINT(41)
        VINT(144)=VINT(144)-VINT(42)
C...Allow FSR for UE (always handle with old showers)
        IF(MSTP(152).EQ.1) THEN
          M41SAV=MSTJ(41)
          IF (MSTJ(41).EQ.10) MSTJ(41)=2
          MSTJ(41)=MOD(MSTJ(41),10)
          CALL PYSHOW(N-1,N,SQRT(PARP(71))*PT)
          MSTJ(41)=M41SAV
        ENDIF
        IF(MINT(31).LT.240) GOTO 220
  270   CONTINUE
        MINT(1)=ISUBSV
        DO 280 J=11,80
          VINT(J)=VINTSV(J)
  280   CONTINUE
      ENDIF
 
C...Format statements for printout.
 5000 FORMAT(/1X,'****** PYMULT: initialization of multiple inter',
     &'actions for MSTP(82) =',I2,' ******')
 5100 FORMAT(8X,'pT0 =',F5.2,' GeV gives sigma(parton-parton) =',1P,
     &D9.2,' mb: rejected')
 5200 FORMAT(8X,'pT0 =',F5.2,' GeV gives sigma(parton-parton) =',1P,
     &D9.2,' mb: accepted')
 
      RETURN
      END
