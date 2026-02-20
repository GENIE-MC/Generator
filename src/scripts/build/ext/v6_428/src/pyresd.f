C*********************************************************************
 
C...PYRESD
C...Allows resonances to decay (including parton showers for hadronic
C...channels).
 
      SUBROUTINE PYRESD(IRES)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Parameter statement for maximum size of showers.
      PARAMETER (MAXNUR=1000)
C...Commonblocks.
      COMMON/PYPART/NPART,NPARTD,IPART(MAXNUR),PTPART(MAXNUR)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYCTAG/NCT,MCT(4000,2)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYPUED/IUED(0:99),RUED(0:99)
      SAVE /PYPART/,/PYJETS/,/PYCTAG/,/PYDAT1/,/PYDAT2/,/PYDAT3/,
     &/PYSUBS/,/PYPARS/,/PYINT1/,/PYINT2/,/PYINT4/,/PYPUED/
C...Local arrays and complex and character variables.
      DIMENSION IREF(50,8),KDCY(3),KFL1(3),KFL2(3),KFL3(3),KEQL(3),
     &KCQM(3),KCQ1(3),KCQ2(3),KCQ3(3),NSD(3),PMMN(4),ILIN(6),
     &HGZ(3,3),COUP(6,4),CORL(2,2,2),PK(6,4),PKK(6,6),CTHE(3),
     &PHI(3),WDTP(0:400),WDTE(0:400,0:5),DPMO(5),VDCY(4),
     &ITJUNC(3),CTM2(3),KCQ(0:10),IANT(4),ITRI(4),IOCT(4),KCQ4(3),
     &KFL4(3)
      COMPLEX FGK,HA(6,6),HC(6,6)
      REAL TIR,UIR
      CHARACTER CODE*9,MASS*9
C...Local arrays.
      DIMENSION PV(10,5),RORD(10),UE(3),BE(3),WTCOR(10)
      DATA WTCOR/2D0,5D0,15D0,60D0,250D0,1500D0,1.2D4,1.2D5,150D0,16D0/
  
C...Functions: momentum in two-particle decays and four-product.
      PAWT(A,B,C)=SQRT((A**2-(B+C)**2)*(A**2-(B-C)**2))/(2D0*A)
 
C...The F, Xi and Xj functions of Gunion and Kunszt
C...(Phys. Rev. D33, 665, plus errata from the authors).
      FGK(I1,I2,I3,I4,I5,I6)=4.*HA(I1,I3)*HC(I2,I6)*(HA(I1,I5)*
     &HC(I1,I4)+HA(I3,I5)*HC(I3,I4))
      DIGK(DT,DU)=-4D0*D34*D56+DT*(3D0*DT+4D0*DU)+DT**2*(DT*DU/
     &(D34*D56)-2D0*(1D0/D34+1D0/D56)*(DT+DU)+2D0*(D34/D56+D56/D34))
      DJGK(DT,DU)=8D0*(D34+D56)**2-8D0*(D34+D56)*(DT+DU)-6D0*DT*DU-
     &2D0*DT*DU*(DT*DU/(D34*D56)-2D0*(1D0/D34+1D0/D56)*(DT+DU)+
     &2D0*(D34/D56+D56/D34))
 
C...Some general constants.
      XW=PARU(102)
      XWV=XW
      IF(MSTP(8).GE.2) XW=1D0-(PMAS(24,1)/PMAS(23,1))**2
      XW1=1D0-XW
      SQMZ=PMAS(23,1)**2
 
      GMMZ=PMAS(23,1)*PMAS(23,2)
      SQMW=PMAS(24,1)**2
      GMMW=PMAS(24,1)*PMAS(24,2)
      SH=VINT(44)
 
C...Boost and rotate to rest frame of incoming partons, 
C...to get proper amount of smearing of decay angles.
      IBST=0
      IF(IRES.EQ.0) THEN
        IBST=1
        IIN1=MINT(84)+1
        IIN2=MINT(84)+2
C...Bug fix 09 OCT 2008 (PS) at 6.4.18: in new shower, the incoming partons 
C...(101,102) are off shell and can have inconsistent momenta, resulting 
C...in boosts larger than unity. However, the corresponding docu partons 
C...(5,6) are kept on shell, and have consistent momenta that can be used 
C...to derive this boost instead. Ultimately, should change the way the new 
C...shower stores intermediate partons, but just using partons (5,6) for now 
C...does define the boost and furnishes a quick and much needed solution.
        IF (MINT(35).EQ.3) THEN
          IIN1=MINT(83)+5
          IIN2=MINT(83)+6
        ENDIF
        ETOTIN=P(IIN1,4)+P(IIN2,4)
        BEXIN=(P(IIN1,1)+P(IIN2,1))/ETOTIN
        BEYIN=(P(IIN1,2)+P(IIN2,2))/ETOTIN
        BEZIN=(P(IIN1,3)+P(IIN2,3))/ETOTIN
        CALL PYROBO(MINT(83)+7,N,0D0,0D0,-BEXIN,-BEYIN,-BEZIN)
        PHIIN=PYANGL(P(MINT(84)+1,1),P(MINT(84)+1,2))
        CALL PYROBO(MINT(83)+7,N,0D0,-PHIIN,0D0,0D0,0D0)
        THEIN=PYANGL(P(MINT(84)+1,3),P(MINT(84)+1,1))
        CALL PYROBO(MINT(83)+7,N,-THEIN,0D0,0D0,0D0,0D0)
      ENDIF
 
C...Reset original resonance configuration.
      DO 100 JT=1,8
        IREF(1,JT)=0
  100 CONTINUE
 
C...Define initial one, two or three objects for subprocess.
      IHDEC=0
      IF(IRES.EQ.0) THEN
        ISUB=MINT(1)
        IF(ISET(ISUB).EQ.1.OR.ISET(ISUB).EQ.3) THEN
          IREF(1,1)=MINT(84)+2+ISET(ISUB)
          IREF(1,4)=MINT(83)+6+ISET(ISUB)
          JTMAX=1
        ELSEIF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.4) THEN
          IREF(1,1)=MINT(84)+1+ISET(ISUB)
          IREF(1,2)=MINT(84)+2+ISET(ISUB)
          IREF(1,4)=MINT(83)+5+ISET(ISUB)
          IREF(1,5)=MINT(83)+6+ISET(ISUB)
          JTMAX=2
        ELSEIF(ISET(ISUB).EQ.5) THEN
          IREF(1,1)=MINT(84)+3
          IREF(1,2)=MINT(84)+4
          IREF(1,3)=MINT(84)+5
          IREF(1,4)=MINT(83)+7
          IREF(1,5)=MINT(83)+8
          IREF(1,6)=MINT(83)+9
          JTMAX=3
        ENDIF
 
C...Define original resonance for odd cases.
      ELSE
        ISUB=0
        IF(K(IRES,2).EQ.25.OR.K(IRES,2).EQ.35.OR.K(IRES,2).EQ.36)
     &  IHDEC=1
        IF(IHDEC.EQ.1) ISUB=3
        IREF(1,1)=IRES
        IREF(1,4)=K(IRES,3)
        IRESTM=IRES
        IF(IREF(1,4).GT.MINT(84)) THEN
  110     ITMPMO=IREF(1,4)
          IF(K(ITMPMO,2).EQ.94) THEN
            IREF(1,4)=K(ITMPMO,3)+(IRESTM-ITMPMO-1)
            IF(K(IREF(1,4),3).LE.MINT(84)) IREF(1,4)=K(IREF(1,4),3)
          ELSEIF(K(ITMPMO,2).EQ.K(IRES,2)) THEN
            IRESTM=ITMPMO
C...Explicitly check that reference particle exists, otherwise stop recursion
            IF(ITMPMO.GT.0.AND.K(ITMPMO,3).GT.0) THEN
              IREF(1,4)=K(ITMPMO,3)
              GOTO 110
            ENDIF
          ENDIF
        ENDIF
        IF(IREF(1,4).GT.MINT(84)) THEN
          EMATCH=1D10
          IREF14=IREF(1,4)
          DO 120 II=MINT(83)+7,MINT(83)+MINT(4)
            IF(K(II,2).EQ.K(IRES,2).AND.ABS(P(II,4)-P(IREF14,4)).LT.
     &      EMATCH) THEN
              IREF(1,4)=II
              EMATCH=ABS(P(II,4)-P(IREF14,4))
            ENDIF
  120     CONTINUE
        ENDIF
        JTMAX=1
      ENDIF
 
C...Check if initial resonance has been moved (in resonance + jet).
      DO 140 JT=1,3
        IF(IREF(1,JT).GT.0) THEN
          IF(K(IREF(1,JT),1).GT.10) THEN
            KFA=IABS(K(IREF(1,JT),2))
            IF(KFA.GE.6.AND.KCHG(PYCOMP(KFA),2).NE.0) THEN
              KDA1=MOD(K(IREF(1,JT),4),MSTU(5))
              KDA2=MOD(K(IREF(1,JT),5),MSTU(5))
              IF(KDA1.GT.IREF(1,JT).AND.KDA1.LE.N) THEN
                IF(K(KDA1,2).EQ.21) KDA1=K(KDA1,5)/MSTU(5)
              ENDIF
              IF(KDA2.GT.IREF(1,JT).AND.KDA2.LE.N) THEN
                IF(K(KDA2,2).EQ.21) KDA2=K(KDA2,4)/MSTU(5)
              ENDIF
              DO 130 I=IREF(1,JT)+1,N
                IF(K(I,2).EQ.K(IREF(1,JT),2).AND.(I.EQ.KDA1.OR.
     &          I.EQ.KDA2)) THEN
                  IREF(1,JT)=I
                  KDA1=MOD(K(IREF(1,JT),4),MSTU(5))
                  KDA2=MOD(K(IREF(1,JT),5),MSTU(5))
                  IF(KDA1.GT.IREF(1,JT).AND.KDA1.LE.N) THEN
                    IF(K(KDA1,2).EQ.21) KDA1=K(KDA1,5)/MSTU(5)
                  ENDIF
                  IF(KDA2.GT.IREF(1,JT).AND.KDA2.LE.N) THEN
                    IF(K(KDA2,2).EQ.21) KDA2=K(KDA2,4)/MSTU(5)
                  ENDIF
                ENDIF
  130         CONTINUE
            ELSE
              KDA=MOD(K(IREF(1,JT),4),MSTU(5))
              IF(MWID(PYCOMP(KFA)).NE.0.AND.KDA.GT.1) IREF(1,JT)=KDA
            ENDIF
          ENDIF
        ENDIF
  140 CONTINUE
 
C...Set decay vertex for initial resonances
      DO 160 JT=1,JTMAX
        DO 150 I=1,4
          V(IREF(1,JT),I)=0D0
  150   CONTINUE
  160 CONTINUE
 
C...Loop over decay history.
      NP=1
      IP=0
  170 IP=IP+1
      NINH=0
      JTMAX=2
      IF(IREF(IP,2).EQ.0) JTMAX=1
      IF(IREF(IP,3).NE.0) JTMAX=3
      IT4=0
      NSAV=N
 
C...Check for Higgs which appears as decay product of user-process.
      IF(ISUB.EQ.0) THEN
        IHDEC=0
        IF(IREF(IP,7).EQ.25.OR.IREF(IP,7).EQ.35.OR.IREF(IP,7)
     &  .EQ.36) IHDEC=1
        IF(IHDEC.EQ.1) ISUB=3
      ENDIF
 
C...Start treatment of one, two or three resonances in parallel.
  180 N=NSAV
      DO 340 JT=1,JTMAX
        ID=IREF(IP,JT)
        KDCY(JT)=0
        KFL1(JT)=0
        KFL2(JT)=0
        KFL3(JT)=0
        KFL4(JT)=0
        KEQL(JT)=0
        NSD(JT)=ID
        ITJUNC(JT)=0
 
C...Check whether particle can/is allowed to decay.
        IF(ID.EQ.0) GOTO 330
        KFA=IABS(K(ID,2))
        KCA=PYCOMP(KFA)
        IF(MWID(KCA).EQ.0) GOTO 330
        IF(K(ID,1).GT.10.OR.MDCY(KCA,1).EQ.0) GOTO 330
        IF(KFA.EQ.6.OR.KFA.EQ.7.OR.KFA.EQ.8.OR.KFA.EQ.17.OR.
     &  KFA.EQ.18) IT4=IT4+1
        K(ID,4)=MSTU(5)*(K(ID,4)/MSTU(5))
        K(ID,5)=MSTU(5)*(K(ID,5)/MSTU(5))
 
C...Choose lifetime and determine decay vertex.
        IF(K(ID,1).EQ.5) THEN
          V(ID,5)=0D0
        ELSEIF(K(ID,1).NE.4) THEN
          V(ID,5)=-PMAS(KCA,4)*LOG(PYR(0))
        ENDIF
        DO 190 J=1,4
          VDCY(J)=V(ID,J)+V(ID,5)*P(ID,J)/P(ID,5)
  190   CONTINUE
 
C...Determine whether decay allowed or not.
        MOUT=0
        IF(MSTJ(22).EQ.2) THEN
          IF(PMAS(KCA,4).GT.PARJ(71)) MOUT=1
        ELSEIF(MSTJ(22).EQ.3) THEN
          IF(VDCY(1)**2+VDCY(2)**2+VDCY(3)**2.GT.PARJ(72)**2) MOUT=1
        ELSEIF(MSTJ(22).EQ.4) THEN
          IF(VDCY(1)**2+VDCY(2)**2.GT.PARJ(73)**2) MOUT=1
          IF(ABS(VDCY(3)).GT.PARJ(74)) MOUT=1
        ENDIF
        IF(MOUT.EQ.1.AND.K(ID,1).NE.5) THEN
          K(ID,1)=4
          GOTO 330
        ENDIF
 
C...Info for selection of decay channel: sign, pairings.
        IF(KCHG(KCA,3).EQ.0) THEN
          IPM=2
        ELSE
          IPM=(5-ISIGN(1,K(ID,2)))/2
        ENDIF
        KFB=0
        IF(JTMAX.EQ.2) THEN
          KFB=IABS(K(IREF(IP,3-JT),2))
        ELSEIF(JTMAX.EQ.3) THEN
          JT2=JT+1-3*(JT/3)
          KFB=IABS(K(IREF(IP,JT2),2))
          IF(KFB.NE.KFA) THEN
            JT2=JT+2-3*((JT+1)/3)
            KFB=IABS(K(IREF(IP,JT2),2))
          ENDIF
        ENDIF
 
C...Select decay channel.
        IF(ISUB.EQ.1.OR.ISUB.EQ.15.OR.ISUB.EQ.19.OR.ISUB.EQ.22.OR.
     &  ISUB.EQ.30.OR.ISUB.EQ.35.OR.ISUB.EQ.141) MINT(61)=1
        CALL PYWIDT(KFA,P(ID,5)**2,WDTP,WDTE)
        WDTE0S=WDTE(0,1)+WDTE(0,IPM)+WDTE(0,4)
        IF(KFB.EQ.KFA) WDTE0S=WDTE0S+WDTE(0,5)
        IF(WDTE0S.LE.0D0) GOTO 330
        RKFL=WDTE0S*PYR(0)
        IDL=0
  200   IDL=IDL+1
        IDC=IDL+MDCY(KCA,2)-1
        RKFL=RKFL-(WDTE(IDL,1)+WDTE(IDL,IPM)+WDTE(IDL,4))
        IF(KFB.EQ.KFA) RKFL=RKFL-WDTE(IDL,5)
        IF(IDL.LT.MDCY(KCA,3).AND.RKFL.GT.0D0) GOTO 200
 
        NPROD=0
C...Read out flavours and colour charges of decay channel chosen.
        KCQM(JT)=KCHG(KCA,2)*ISIGN(1,K(ID,2))
        IF(KCQM(JT).EQ.-2) KCQM(JT)=2
        KFL1(JT)=KFDP(IDC,1)*ISIGN(1,K(ID,2))
        KFC1A=PYCOMP(IABS(KFL1(JT)))
        IF(KCHG(KFC1A,3).EQ.0) KFL1(JT)=IABS(KFL1(JT))
        NPROD=NPROD+1
        KCQ1(JT)=KCHG(KFC1A,2)*ISIGN(1,KFL1(JT))
        IF(KCQ1(JT).EQ.-2) KCQ1(JT)=2
        KFL2(JT)=KFDP(IDC,2)*ISIGN(1,K(ID,2))
        KFC2A=PYCOMP(IABS(KFL2(JT)))
        IF(KCHG(KFC2A,3).EQ.0) KFL2(JT)=IABS(KFL2(JT))
        KCQ2(JT)=KCHG(KFC2A,2)*ISIGN(1,KFL2(JT))
        IF(KCQ2(JT).EQ.-2) KCQ2(JT)=2
        NPROD=NPROD+1
        KFL3(JT)=KFDP(IDC,3)*ISIGN(1,K(ID,2))
        KCQ3(JT)=0
        KFL4(JT)=KFDP(IDC,4)*ISIGN(1,K(ID,2))
        KCQ4(JT)=0        
        IF(KFL3(JT).NE.0) THEN
          KFC3A=PYCOMP(IABS(KFL3(JT)))
          IF(KCHG(KFC3A,3).EQ.0) KFL3(JT)=IABS(KFL3(JT))
          KCQ3(JT)=KCHG(KFC3A,2)*ISIGN(1,KFL3(JT))
          IF(KCQ3(JT).EQ.-2) KCQ3(JT)=2
          NPROD=NPROD+1
          IF(KFL4(JT).NE.0) THEN
            KFC4A=PYCOMP(IABS(KFL4(JT)))
            IF(KCHG(KFC4A,3).EQ.0) KFL4(JT)=IABS(KFL4(JT))
            KCQ4(JT)=KCHG(KFC4A,2)*ISIGN(1,KFL4(JT))
            IF(KCQ4(JT).EQ.-2) KCQ4(JT)=2
            NPROD=NPROD+1
          ENDIF
        ENDIF
 
C...Set/save further info on channel.
        KDCY(JT)=1
        IF(KFB.EQ.KFA) KEQL(JT)=MDME(IDC,1)
        NSD(JT)=N
        HGZ(JT,1)=VINT(111)
        HGZ(JT,2)=VINT(112)
        HGZ(JT,3)=VINT(114)
        JTZ=JT
 
        PXSUM=0D0
C...Select masses; to begin with assume resonances narrow.
        DO 220 I=1,4
          P(N+I,5)=0D0
          PMMN(I)=0D0
          IF(I.EQ.1) THEN
            KFLW=IABS(KFL1(JT))
            KCW=KFC1A
          ELSEIF(I.EQ.2) THEN
            KFLW=IABS(KFL2(JT))
            KCW=KFC2A
          ELSEIF(I.EQ.3) THEN
            IF(KFL3(JT).EQ.0) GOTO 220
            KFLW=IABS(KFL3(JT))
            KCW=KFC3A
          ELSEIF(I.EQ.4) THEN
            IF(KFL4(JT).EQ.0) GOTO 220
            KFLW=IABS(KFL4(JT))
            KCW=KFC4A
          ENDIF
          P(N+I,5)=PMAS(KCW,1)
          PXSUM=PXSUM+P(N+I,5)
CMRENNA++
C...This prevents SUSY/t particles from becoming too light.
          IF(KFLW/KSUSY1.EQ.1.OR.KFLW/KSUSY1.EQ.2) THEN
            PMMN(I)=PMAS(KCW,1)
            DO 210 IDC=MDCY(KCW,2),MDCY(KCW,2)+MDCY(KCW,3)-1
              IF(MDME(IDC,1).GT.0.AND.BRAT(IDC).GT.1E-4) THEN
                PMSUM=PMAS(PYCOMP(KFDP(IDC,1)),1)+
     &              PMAS(PYCOMP(KFDP(IDC,2)),1)
                IF(KFDP(IDC,3).NE.0) PMSUM=PMSUM+
     &              PMAS(PYCOMP(KFDP(IDC,3)),1)
                IF(KFDP(IDC,4).NE.0) PMSUM=PMSUM+
     &              PMAS(PYCOMP(KFDP(IDC,4)),1)
                PMMN(I)=MIN(PMMN(I),PMSUM)
              ENDIF
 210        CONTINUE
C   MRENNA--
          ELSEIF(KFLW.EQ.6) THEN
            PMMN(I)=PMAS(24,1)+PMAS(5,1)
          ENDIF
C...UED: select a graviton mass from continuous distribution
C...(stored in PMAS(39,1) so no value returned)
          IF (IUED(1).EQ.1.AND.IUED(2).EQ.1.AND.KFLW.EQ.39) 
     &         CALL PYGRAM(1)
 220    CONTINUE
        
C...Check which two out of three are widest.
        IWID1=1
        IWID2=2
        PWID1=PMAS(KFC1A,2)
        PWID2=PMAS(KFC2A,2)
        KFLW1=IABS(KFL1(JT))
        KFLW2=IABS(KFL2(JT))
        IF(KFL3(JT).NE.0) THEN
          PWID3=PMAS(KFC3A,2)
          IF(PWID3.GT.PWID1.AND.PWID2.GE.PWID1) THEN
            IWID1=3
            PWID1=PWID3
            KFLW1=IABS(KFL3(JT))
          ELSEIF(PWID3.GT.PWID2) THEN
            IWID2=3
            PWID2=PWID3
            KFLW2=IABS(KFL3(JT))
          ENDIF
        ENDIF
        IF(KFL4(JT).NE.0) THEN
          PWID4=PMAS(KFC4A,2)
          IF(PWID4.GT.PWID1.AND.PWID2.GE.PWID1) THEN
            IWID1=4
            PWID1=PWID4
            KFLW1=IABS(KFL4(JT))
          ELSEIF(PWID4.GT.PWID2) THEN
            IWID2=4
            PWID2=PWID4
            KFLW2=IABS(KFL4(JT))
          ENDIF
        ENDIF
 
C...If all narrow then only check that masses consistent.
        IF(MSTP(42).LE.0.OR.(PWID1.LT.PARP(41).AND.
     &  PWID2.LT.PARP(41))) THEN
CMRENNA++
C....Handle near degeneracy cases.
          IF(KFA/KSUSY1.EQ.1.OR.KFA/KSUSY1.EQ.2) THEN
            IF(P(N+1,5)+P(N+2,5)+P(N+3,5).GT.P(ID,5)) THEN
              P(N+1,5)=P(ID,5)-P(N+2,5)-0.5D0
              IF(P(N+1,5).LT.0D0) P(N+1,5)=0D0
            ENDIF
          ENDIF
CMRENNA--
          IF(PXSUM.GT.P(ID,5)) THEN
            CALL PYERRM(13,'(PYRESD:) daughter masses too large')
            MINT(51)=1
            GOTO 720
          ELSEIF(PXSUM+PARJ(64).GT.P(ID,5)) THEN
            CALL PYERRM(3,'(PYRESD:) masses+PARJ(64) too large')
            MINT(51)=1
            GOTO 720
          ENDIF
 
C...For three wide resonances select narrower of three
C...according to BW decoupled from rest.
        ELSE
          PMTOT=P(ID,5)
          IF(KFL3(JT).NE.0) THEN
            IWID3=6-IWID1-IWID2
            KFLW3=IABS(KFL1(JT))+IABS(KFL2(JT))+IABS(KFL3(JT))-
     &      KFLW1-KFLW2
            LOOP=0
  230       LOOP=LOOP+1
            P(N+IWID3,5)=PYMASS(KFLW3)
            IF(LOOP.LE.10.AND. P(N+IWID3,5).LE.PMMN(IWID3)) GOTO 230
            PMTOT=PMTOT-P(N+IWID3,5)
          ENDIF
C...Select other two correlated within remaining phase space.
          IF(IP.EQ.1) THEN
            CKIN45=CKIN(45)
            CKIN47=CKIN(47)
            CKIN(45)=MAX(PMMN(IWID1),CKIN(45))
            CKIN(47)=MAX(PMMN(IWID2),CKIN(47))
            CALL PYOFSH(2,KFA,KFLW1,KFLW2,PMTOT,P(N+IWID1,5),
     &      P(N+IWID2,5))
            CKIN(45)=CKIN45
            CKIN(47)=CKIN47
          ELSE
            CKIN(49)=PMMN(IWID1)
            CKIN(50)=PMMN(IWID2)
            CALL PYOFSH(5,KFA,KFLW1,KFLW2,PMTOT,P(N+IWID1,5),
     &      P(N+IWID2,5))
            CKIN(49)=0D0
            CKIN(50)=0D0
          ENDIF
          IF(MINT(51).EQ.1) GOTO 720
        ENDIF
 
C...Begin fill decay products, with colour flow for coloured objects.
        MSTU10=MSTU(10)
        MSTU(10)=1
        MSTU(19)=1


C...Three-body decays 
        IF(KFL3(JT).NE.0.OR.KFL4(JT).NE.0) THEN
          DO 250 I=N+1,N+NPROD
            DO 240 J=1,5
              K(I,J)=0
              V(I,J)=0D0
  240       CONTINUE
            MCT(I,1)=0
            MCT(I,2)=0
  250     CONTINUE
          K(N+1,1)=1
          K(N+1,2)=KFL1(JT)
          K(N+2,1)=1
          K(N+2,2)=KFL2(JT)
          K(N+3,1)=1
          K(N+3,2)=KFL3(JT)
          IF(KFL4(JT).NE.0) THEN
            K(N+4,1)=1
            K(N+4,2)=KFL4(JT)
          ENDIF
          IDIN=ID

C...Generate kinematics (default is flat)
          IF(KFL4(JT).EQ.0) THEN
            CALL PYTBDY(IDIN)
          ELSE
            PS=P(N+1,5)+P(N+2,5)+P(N+3,5)+P(N+4,5)
            ND=4
            PV(1,1)=0D0
            PV(1,2)=0D0
            PV(1,3)=0D0
            PV(1,4)=P(IDIN,5)
            PV(1,5)=P(IDIN,5)
C...Calculate maximum weight ND-particle decay.
            PV(ND,5)=P(N+ND,5)
            WTMAX=1D0/WTCOR(ND-2)
            PMAX=PV(1,5)-PS+P(N+ND,5)
            PMIN=0D0
            DO 381 IL=ND-1,1,-1
              PMAX=PMAX+P(N+IL,5)
              PMIN=PMIN+P(N+IL+1,5)
              WTMAX=WTMAX*PAWT(PMAX,PMIN,P(N+IL,5))
 381        CONTINUE

C...M-generator gives weight. If rejected, try again.

 411        RORD(1)=1D0
            DO 441 IL1=2,ND-1
              RSAV=PYR(0)
              DO 421 IL2=IL1-1,1,-1
                IF(RSAV.LE.RORD(IL2)) GOTO 431
                RORD(IL2+1)=RORD(IL2)
 421          CONTINUE
 431          RORD(IL2+1)=RSAV
 441        CONTINUE
            RORD(ND)=0D0
            WT=1D0
            DO 451 IL=ND-1,1,-1
              PV(IL,5)=PV(IL+1,5)+P(N+IL,5)+(RORD(IL)-RORD(IL+1))*
     &             (PV(1,5)-PS)
              WT=WT*PAWT(PV(IL,5),PV(IL+1,5),P(N+IL,5))
 451        CONTINUE
            IF(WT.LT.PYR(0)*WTMAX) GOTO 411

C...Perform two-particle decays in respective CM frame.
            DO 481 IL=1,ND-1
              PA=PAWT(PV(IL,5),PV(IL+1,5),P(N+IL,5))
              UE(3)=2D0*PYR(0)-1D0
              PHIX=PARU(2)*PYR(0)
              UE(1)=SQRT(1D0-UE(3)**2)*COS(PHIX)
              UE(2)=SQRT(1D0-UE(3)**2)*SIN(PHIX)
              DO 471 J=1,3
                P(N+IL,J)=PA*UE(J)
                PV(IL+1,J)=-PA*UE(J)
 471          CONTINUE
              P(N+IL,4)=SQRT(PA**2+P(N+IL,5)**2)
              PV(IL+1,4)=SQRT(PA**2+PV(IL+1,5)**2)
 481        CONTINUE

C...Lorentz transform decay products to lab frame.
            DO 491 J=1,4
              P(N+ND,J)=PV(ND,J)
 491        CONTINUE
            DO 531 IL=ND-1,1,-1
              DO 501 J=1,3
                BE(J)=PV(IL,J)/PV(IL,4)
 501          CONTINUE
              GA=PV(IL,4)/PV(IL,5)
              DO 521 I=N+IL,N+ND
                BEP=BE(1)*P(I,1)+BE(2)*P(I,2)+BE(3)*P(I,3)
                DO 511 J=1,3
                  P(I,J)=P(I,J)+GA*(GA*BEP/(1D0+GA)+P(I,4))*BE(J)
 511            CONTINUE
                P(I,4)=GA*(P(I,4)+BEP)
 521          CONTINUE
 531        CONTINUE

          ENDIF

C...Set generic colour flows whenever unambiguous,
C...(independently of the order of the decay products)
C...Sum up total colour content
          NANT=0
          NTRI=0
          NOCT=0
          KCQ(0)=KCQM(JT)
          KCQ(1)=KCQ1(JT)
          KCQ(2)=KCQ2(JT)
          KCQ(3)=KCQ3(JT)
          KCQ(4)=KCQ4(JT)
          DO 255 J=0,NPROD
            IF (KCQ(J).EQ.-1) THEN
              NANT=NANT+1
              IANT(NANT)=N+J
            ELSEIF (KCQ(J).EQ.1) THEN
              NTRI=NTRI+1              
              ITRI(NTRI)=N+J
            ELSEIF (KCQ(J).EQ.2) THEN 
              NOCT=NOCT+1
              IOCT(NOCT)=N+J
            ENDIF
 255      CONTINUE
          
C...Set color flow for generic 1 -> N processes (N arbitrary)
          IF (NTRI.EQ.0.AND.NANT.EQ.0.AND.NOCT.EQ.0) THEN
C...All singlets: do nothing
            
          ELSEIF (NOCT.EQ.2.AND.NTRI.EQ.0.AND.NANT.EQ.0) THEN
C...Two octets, zero triplets, n singlets:
            IF (KCQ(0).EQ.2) THEN
C...8 -> 8 + n(1) 
              K(ID,4)=K(ID,4)+IOCT(2)
              K(ID,5)=K(ID,5)+IOCT(2)
              K(IOCT(2),1)=3
              K(IOCT(2),4)=MSTU(5)*ID
              K(IOCT(2),5)=MSTU(5)*ID
              MCT(IOCT(2),1)=MCT(ID,1)
              MCT(IOCT(2),2)=MCT(ID,2)
            ELSE
C...1 -> 8 + 8 + n(1)
              K(IOCT(1),1)=3
              K(IOCT(1),4)=MSTU(5)*IOCT(2)
              K(IOCT(1),5)=MSTU(5)*IOCT(2)
              K(IOCT(2),1)=3
              K(IOCT(2),4)=MSTU(5)*IOCT(1)
              K(IOCT(2),5)=MSTU(5)*IOCT(1)
              NCT=NCT+1
              MCT(IOCT(1),1)=NCT
              MCT(IOCT(2),2)=NCT
              NCT=NCT+1
              MCT(IOCT(2),1)=NCT
              MCT(IOCT(1),2)=NCT
            ENDIF
            
          ELSEIF (NTRI+NANT.EQ.2.AND.NOCT.EQ.0) THEN
C...Two triplets, zero octets, n singlets.            
            IF (KCQ(0).EQ.1) THEN
C...3 -> 3 + n(1)
              K(ID,4)=K(ID,4)+ITRI(2)
              K(ITRI(2),1)=3
              K(ITRI(2),4)=MSTU(5)*ID
              MCT(ITRI(2),1)=MCT(ID,1)
            ELSEIF (KCQ(0).EQ.-1) THEN
C...3bar -> 3bar + n(1)              
              K(ID,5)=K(ID,5)+IANT(2)
              K(IANT(2),1)=3
              K(IANT(2),5)=MSTU(5)*ID
              MCT(IANT(2),2)=MCT(ID,2)
            ELSE
C...1 -> 3 + 3bar + n(1)
              K(ITRI(1),1)=3
              K(ITRI(1),4)=MSTU(5)*IANT(1)
              K(IANT(1),1)=3
              K(IANT(1),5)=MSTU(5)*ITRI(1)
              NCT=NCT+1
              MCT(ITRI(1),1)=NCT
              MCT(IANT(1),2)=NCT
            ENDIF
            
          ELSEIF(NTRI+NANT.EQ.2.AND.NOCT.EQ.1) THEN
C...Two triplets, one octet, n singlets.            
            IF (KCQ(0).EQ.2) THEN
C...8 -> 3 + 3bar + n(1)
              K(ID,4)=K(ID,4)+ITRI(1)
              K(ID,5)=K(ID,5)+IANT(1)
              K(ITRI(1),1)=3
              K(ITRI(1),4)=MSTU(5)*ID
              K(IANT(1),1)=3
              K(IANT(1),5)=MSTU(5)*ID
              MCT(ITRI(1),1)=MCT(ID,1)
              MCT(IANT(1),2)=MCT(ID,2)
            ELSEIF (KCQ(0).EQ.1) THEN
C...3 -> 8 + 3 + n(1)
              K(ID,4)=K(ID,4)+IOCT(1)
              K(IOCT(1),1)=3
              K(IOCT(1),4)=MSTU(5)*ID
              K(IOCT(1),5)=MSTU(5)*ITRI(2)
              K(ITRI(2),1)=3
              K(ITRI(2),4)=MSTU(5)*IOCT(1)
              MCT(IOCT(1),1)=MCT(ID,1)
              NCT=NCT+1
              MCT(IOCT(1),2)=NCT
              MCT(ITRI(2),1)=NCT
            ELSEIF (KCQ(0).EQ.-1) THEN
C...3bar -> 8 + 3bar + n(1)
              K(ID,5)=K(ID,5)+IOCT(1)
              K(IOCT(1),1)=3
              K(IOCT(1),5)=MSTU(5)*ID
              K(IOCT(1),4)=MSTU(5)*IANT(2)
              K(IANT(2),1)=3
              K(IANT(2),5)=MSTU(5)*IOCT(1)
              MCT(IOCT(1),2)=MCT(ID,2)
              NCT=NCT+1
              MCT(IOCT(1),1)=NCT
              MCT(IANT(2),2)=NCT
            ELSE
C...1 -> 3 + 3bar + 8 + n(1)
              K(ITRI(1),1)=3
              K(ITRI(1),4)=MSTU(5)*IOCT(1)
              K(IOCT(1),1)=3
              K(IOCT(1),5)=MSTU(5)*ITRI(1)
              K(IOCT(1),4)=MSTU(5)*IANT(1)
              K(IANT(1),1)=3
              K(IANT(1),5)=MSTU(5)*IOCT(1)
              NCT=NCT+1
              MCT(ITRI(1),1)=NCT
              MCT(IOCT(1),2)=NCT
              NCT=NCT+1
              MCT(IOCT(1),1)=NCT
              MCT(IANT(1),2)=NCT
            ENDIF
         ELSEIF(NTRI+NANT.EQ.4) THEN
C...
            IF (KCQ(0).EQ.1) THEN
C...3 -> 3 + n(1) -> 3 + 3bar
              K(ID,4)=K(ID,4)+ITRI(2)
              K(ITRI(2),1)=3
              K(ITRI(2),4)=MSTU(5)*ID
              MCT(ITRI(2),1)=MCT(ID,1)
              K(ITRI(3),1)=3
              K(ITRI(3),4)=MSTU(5)*IANT(1)
              K(IANT(1),1)=3
              K(IANT(1),5)=MSTU(5)*ITRI(3)
              NCT=NCT+1
              MCT(ITRI(3),1)=NCT
              MCT(IANT(1),2)=NCT
            ELSEIF (KCQ(0).EQ.-1) THEN
C...3bar -> 3bar + n(1) -> 3 + 3bar             
              K(ID,5)=K(ID,5)+IANT(2)
              K(IANT(2),1)=3
              K(IANT(2),5)=MSTU(5)*ID
              MCT(IANT(2),2)=MCT(ID,2)
              K(ITRI(1),1)=3
              K(ITRI(1),4)=MSTU(5)*IANT(3)
              K(IANT(3),1)=3
              K(IANT(3),5)=MSTU(5)*ITRI(1)
              NCT=NCT+1
              MCT(ITRI(1),1)=NCT
              MCT(IANT(3),2)=NCT
            ENDIF
          ELSEIF(KFL4(JT).NE.0) THEN
            CALL PYERRM(21,'(PYRESD:) unknown 4-bdy decay')
CPS-- End of generic cases 
C...(could three octets also be handled?)
C...(could (some of) the RPV cases be made generic as well?)

C...Special cases (= old treatment)
C...Set colour flow for t -> W + b + Z.
          ELSEIF(KFA.EQ.6) THEN
            K(N+2,1)=3
            ISID=4
            IF(KCQM(JT).EQ.-1) ISID=5
            IDAU=N+2
            K(ID,ISID)=K(ID,ISID)+IDAU
            K(IDAU,ISID)=MSTU(5)*ID
 
C...Set colour flow in three-body decays - programmed as special cases.
 
          ELSEIF(KFC2A.LE.6) THEN
            K(N+2,1)=3
            K(N+3,1)=3
            ISID=4
            IF(KFL2(JT).LT.0) ISID=5
            K(N+2,ISID)=MSTU(5)*(N+3)
            K(N+3,9-ISID)=MSTU(5)*(N+2)
C...PS++: Bugfix 16 MAR 2006 for 3-body squark decays (e.g. via SLHA)
          ELSEIF(KFA.GT.KSUSY1.AND.MOD(KFA,KSUSY1).LT.10
     &          .AND.KFL3(JT).NE.0) THEN
            KQSUMA=IABS(KCQ1(JT))+IABS(KCQ2(JT))+IABS(KCQ3(JT))
C...3-body decays of squarks to colour singlets plus one quark
            IF (KQSUMA.EQ.1) THEN
C...Find quark
              IQ=0
              IF (KCQ1(JT).NE.0) IQ=1
              IF (KCQ2(JT).NE.0) IQ=2
              IF (KCQ3(JT).NE.0) IQ=3
              ISID=4
              IF (K(N+IQ,2).LT.0) ISID=5
              K(N+IQ,1)=3
              K(ID,ISID)=K(ID,ISID)+(N+IQ)
              K(N+IQ,ISID)=MSTU(5)*ID
            ENDIF
C...PS--
          ELSEIF(KFL1(JT).EQ.KSUSY1+21) THEN
            K(N+1,1)=3
            K(N+2,1)=3
            K(N+3,1)=3
            ISID=4
            IF(KFL2(JT).LT.0) ISID=5
            K(N+1,ISID)=MSTU(5)*(N+2)
            K(N+1,9-ISID)=MSTU(5)*(N+3)
            K(N+2,ISID)=MSTU(5)*(N+1)
            K(N+3,9-ISID)=MSTU(5)*(N+1)
          ELSEIF(KFA.EQ.KSUSY1+21) THEN
            K(N+2,1)=3
            K(N+3,1)=3
            ISID=4
            IF(KFL2(JT).LT.0) ISID=5
            K(ID,ISID)=K(ID,ISID)+(N+2)
            K(ID,9-ISID)=K(ID,9-ISID)+(N+3)
            K(N+2,ISID)=MSTU(5)*ID
            K(N+3,9-ISID)=MSTU(5)*ID
CMRENNA--
 
          ELSEIF(KFA.GE.KSUSY1+22.AND.KFA.LE.KSUSY1+37.AND.
     &    IABS(KCQ2(JT)).EQ.1) THEN
            K(N+2,1)=3
            K(N+3,1)=3
            ISID=4
            IF(KFL2(JT).LT.0) ISID=5
            K(N+2,ISID)=MSTU(5)*(N+3)
            K(N+3,9-ISID)=MSTU(5)*(N+2)
          ENDIF
           
CXXX      NSAV=N
          
C...Set colour flow in three-body decays with baryon number violation.
C...Neutralino and chargino decays first.
          KCQSUM=KCQ1(JT)+KCQ2(JT)+KCQ3(JT)
          IF(KCQM(JT).EQ.0.AND.IABS(KCQSUM).EQ.3) THEN
            ITJUNC(JT)=(1+(1-KCQ1(JT))/2)
            K(N+4,4)=ITJUNC(JT)*MSTU(5)
C...Insert junction to keep track of colours.
            IF(KCQ1(JT).NE.0) K(N+1,1)=3
            IF(KCQ2(JT).NE.0) K(N+2,1)=3
            IF(KCQ3(JT).NE.0) K(N+3,1)=3
C...Set special junction codes:
            K(N+4,1)=42
            K(N+4,2)=88
 
C...Order decay products by invariant mass. (will be used in PYSTRF).
            PM12=P(N+1,4)*P(N+2,4)-P(N+1,1)*P(N+2,1)-P(N+1,2)*P(N+2,2)-
     &      P(N+1,3)*P(N+2,3)
            PM13=P(N+1,4)*P(N+3,4)-P(N+1,1)*P(N+3,1)-P(N+1,2)*P(N+3,2)-
     &      P(N+1,3)*P(N+3,3)
            PM23=P(N+2,4)*P(N+3,4)-P(N+2,1)*P(N+3,1)-P(N+2,2)*P(N+3,2)-
     &      P(N+2,3)*P(N+3,3)
            IF(PM12.LT.PM13.AND.PM12.LT.PM23) THEN
              K(N+4,4)=N+3+K(N+4,4)
              K(N+4,5)=N+1+MSTU(5)*(N+2)
            ELSEIF(PM13.LT.PM23) THEN
              K(N+4,4)=N+2+K(N+4,4)
              K(N+4,5)=N+1+MSTU(5)*(N+3)
            ELSE
              K(N+4,4)=N+1+K(N+4,4)
              K(N+4,5)=N+2+MSTU(5)*(N+3)
            ENDIF
            DO 260 J=1,5
              P(N+4,J)=0D0
              V(N+4,J)=0D0
  260       CONTINUE
C...Connect daughters to junction.
            DO 270 II=N+1,N+3
              K(II,4)=0
              K(II,5)=0
              K(II,ITJUNC(JT)+3)=MSTU(5)*(N+4)
  270       CONTINUE
C...Particle counter should be stepped up one extra for junction.
            N=N+1
 
C...Gluino decays.
          ELSEIF (KCQM(JT).EQ.2.AND.IABS(KCQSUM).EQ.3) THEN
            ITJUNC(JT)=(5+(1-KCQ1(JT))/2)
            K(N+4,4)=ITJUNC(JT)*MSTU(5)
C...Insert junction to keep track of colours.
            IF(KCQ1(JT).NE.0) K(N+1,1)=3
            IF(KCQ2(JT).NE.0) K(N+2,1)=3
            IF(KCQ3(JT).NE.0) K(N+3,1)=3
            K(N+4,1)=42
            K(N+4,2)=88
            DO 280 J=1,5
              P(N+4,J)=0D0
              V(N+4,J)=0D0
  280       CONTINUE
            CTMSUM=0D0
            DO 290 II=N+1,N+3
              K(II,4)=0
              K(II,5)=0
C...Start by connecting all daughters to junction.
              K(II,ITJUNC(JT)-1)=MSTU(5)*(N+4)
C...Only consider colour topologies with off shell resonances.
              RMQ1=PMAS(PYCOMP(K(II,2)),1)
              RMRES=PMAS(PYCOMP(KSUSY1+IABS(K(II,2))),1)
              RMGLU=PMAS(PYCOMP(KSUSY1+21),1)
              IF (RMGLU-RMQ1.LT.RMRES) THEN
C...Calculate propagators for each colour topology.
                RM2Q23=RMGLU**2+RMQ1**2-2D0*(P(II,4)*P(ID,4)+P(II,1)
     &               *P(ID,1)+P(II,2)*P(ID,2)+P(II,3)*P(ID,3))
                CTM2(II-N)=1D0/(RM2Q23-RMRES**2)**2
              ELSE
                CTM2(II-N)=0D0
              ENDIF
              CTMSUM=CTMSUM+CTM2(II-N)
  290       CONTINUE
            CTMSUM=PYR(0)*CTMSUM
C...Select colour topology J, with most off shell least likely.
            J=0
  300       J=J+1
            CTMSUM=CTMSUM-CTM2(J)
            IF (CTMSUM.GT.0D0) GOTO 300
C...The lucky winner gets its colour (anti-colour) directly from gluino.
            K(N+J,ITJUNC(JT)-1)=MSTU(5)*ID
            K(ID,ITJUNC(JT)-1)=N+J+(K(ID,ITJUNC(JT)-1)/MSTU(5))*MSTU(5)
C...The other gluino colour is connected to junction
            K(ID,10-ITJUNC(JT))=N+4+(K(ID,10-ITJUNC(JT))/MSTU(5))*
     &      MSTU(5)
            K(N+4,4)=K(N+4,4)+ID
C...Lastly, connect junction to remaining daughters.
            K(N+4,5)=N+1+MOD(J,3)+MSTU(5)*(N+1+MOD(J+1,3))
C...Particle counter should be stepped up one extra for junction.
            N=N+1
          ENDIF
 
C...Update particle counter.
          N=N+NPROD

C...2) Everything else two-body decay.
        ELSE
          CALL PY2ENT(N+1,KFL1(JT),KFL2(JT),P(ID,5))
          MCT(N-1,1)=0
          MCT(N-1,2)=0
          MCT(N,1)=0
          MCT(N,2)=0
C...First set colour flow as if mother colour singlet.
          IF(KCQ1(JT).NE.0) THEN
            K(N-1,1)=3
            IF(KCQ1(JT).NE.-1) K(N-1,4)=MSTU(5)*N
            IF(KCQ1(JT).NE.1) K(N-1,5)=MSTU(5)*N
          ENDIF
          IF(KCQ2(JT).NE.0) THEN
            K(N,1)=3
            IF(KCQ2(JT).NE.-1) K(N,4)=MSTU(5)*(N-1)
            IF(KCQ2(JT).NE.1) K(N,5)=MSTU(5)*(N-1)
          ENDIF
C...Then redirect colour flow if mother (anti)triplet.
          IF(KCQM(JT).EQ.0) THEN
          ELSEIF(KCQM(JT).NE.2) THEN
            ISID=4
            IF(KCQM(JT).EQ.-1) ISID=5
            IDAU=N-1
            IF(KCQ1(JT).EQ.0.OR.KCQ2(JT).EQ.2) IDAU=N
            K(ID,ISID)=K(ID,ISID)+IDAU
            K(IDAU,ISID)=MSTU(5)*ID
C...Then redirect colour flow if mother octet.
          ELSEIF(KCQ1(JT).EQ.0.OR.KCQ2(JT).EQ.0) THEN
            IDAU=N-1
            IF(KCQ1(JT).EQ.0) IDAU=N
            K(ID,4)=K(ID,4)+IDAU
            K(ID,5)=K(ID,5)+IDAU
            K(IDAU,4)=MSTU(5)*ID
            K(IDAU,5)=MSTU(5)*ID
          ELSE
            ISID=4
            IF(KCQ1(JT).EQ.-1) ISID=5
            IF(KCQ1(JT).EQ.2) ISID=INT(4.5D0+PYR(0))
            K(ID,ISID)=K(ID,ISID)+(N-1)
            K(ID,9-ISID)=K(ID,9-ISID)+N
            K(N-1,ISID)=MSTU(5)*ID
            K(N,9-ISID)=MSTU(5)*ID
          ENDIF
 
C...Insert junction
          IF(IABS(KCQ1(JT)+KCQ2(JT)-KCQM(JT)).EQ.3) THEN
            N=N+1
C...~q* mother: type 3 junction. ~q mother: type 4.
            ITJUNC(JT)=(7+KCQM(JT))/2
C...Specify junction KF and set colour flow from junction
            K(N,1)=42
            K(N,2)=88
            K(N,3)=ID
C...Junction type encoded together with mother:
            K(N,4)=ID+ITJUNC(JT)*MSTU(5)
            K(N,5)=N-1+MSTU(5)*(N-2)
C...Zero P and V for junction (V filled later)
            DO 310 J=1,5
              P(N,J)=0D0
              V(N,J)=0D0
  310       CONTINUE
C...Set colour flow from mother to junction
            K(ID,8-ITJUNC(JT))= N + MSTU(5)*(K(ID,8-ITJUNC(JT))/MSTU(5))
C...Set colour flow from daughters to junction
            DO 320 II=N-2,N-1
              K(II,4) = 0
              K(II,5) = 0
C...(Anti-)colour mother is junction.
              K(II,1+ITJUNC(JT)) = MSTU(5)*N
  320       CONTINUE
          ENDIF
        ENDIF
 
C...End loop over resonances for daughter flavour and mass selection.
        MSTU(10)=MSTU10
  330   IF(MWID(KCA).NE.0.AND.(KFL1(JT).EQ.0.OR.KFL3(JT).NE.0))
     &  NINH=NINH+1
        IF(IRES.GT.0.AND.MWID(KCA).NE.0.AND.MDCY(KCA,1).NE.0.AND.
     &  KFL1(JT).EQ.0) THEN
          WRITE(CODE,'(I9)') K(ID,2)
          WRITE(MASS,'(F9.3)') P(ID,5)
          CALL PYERRM(3,'(PYRESD:) Failed to decay particle'//
     &    CODE//' with mass'//MASS)
          MINT(51)=1
          GOTO 720
        ENDIF
  340 CONTINUE
 
C...Check for allowed combinations. Skip if no decays.
      IF(JTMAX.EQ.1) THEN
        IF(KDCY(1).EQ.0) GOTO 710
      ELSEIF(JTMAX.EQ.2) THEN
        IF(KDCY(1).EQ.0.AND.KDCY(2).EQ.0) GOTO 710
        IF(KEQL(1).EQ.4.AND.KEQL(2).EQ.4) GOTO 180
        IF(KEQL(1).EQ.5.AND.KEQL(2).EQ.5) GOTO 180
      ELSEIF(JTMAX.EQ.3) THEN
        IF(KDCY(1).EQ.0.AND.KDCY(2).EQ.0.AND.KDCY(3).EQ.0) GOTO 710
        IF(KEQL(1).EQ.4.AND.KEQL(2).EQ.4) GOTO 180
        IF(KEQL(1).EQ.4.AND.KEQL(3).EQ.4) GOTO 180
        IF(KEQL(2).EQ.4.AND.KEQL(3).EQ.4) GOTO 180
        IF(KEQL(1).EQ.5.AND.KEQL(2).EQ.5) GOTO 180
        IF(KEQL(1).EQ.5.AND.KEQL(3).EQ.5) GOTO 180
        IF(KEQL(2).EQ.5.AND.KEQL(3).EQ.5) GOTO 180
      ENDIF
 
C...Special case: matrix element option for Z0 decay to quarks.
      IF(MSTP(48).EQ.1.AND.ISUB.EQ.1.AND.JTMAX.EQ.1.AND.
     &IABS(MINT(11)).EQ.11.AND.IABS(KFL1(1)).LE.5) THEN
 
C...Check consistency of MSTJ options set.
        IF(MSTJ(109).EQ.2.AND.MSTJ(110).NE.1) THEN
          CALL PYERRM(6,
     &    '(PYRESD:) MSTJ(109) value requires MSTJ(110) = 1')
          MSTJ(110)=1
        ENDIF
        IF(MSTJ(109).EQ.2.AND.MSTJ(111).NE.0) THEN
          CALL PYERRM(6,
     &    '(PYRESD:) MSTJ(109) value requires MSTJ(111) = 0')
 
          MSTJ(111)=0
        ENDIF
 
C...Select alpha_strong behaviour.
        MST111=MSTU(111)
        PAR112=PARU(112)
        MSTU(111)=MSTJ(108)
        IF(MSTJ(108).EQ.2.AND.(MSTJ(101).EQ.0.OR.MSTJ(101).EQ.1))
     &  MSTU(111)=1
        PARU(112)=PARJ(121)
        IF(MSTU(111).EQ.2) PARU(112)=PARJ(122)
 
C...Find axial fraction in total cross section for scalar gluon model.
        PARJ(171)=0D0
        IF((IABS(MSTJ(101)).EQ.1.AND.MSTJ(109).EQ.1).OR.
     &  (MSTJ(101).EQ.5.AND.MSTJ(49).EQ.1)) THEN
          POLL=1D0-PARJ(131)*PARJ(132)
          SFF=1D0/(16D0*XW*XW1)
          SFW=P(ID,5)**4/((P(ID,5)**2-PARJ(123)**2)**2+
     &    (PARJ(123)*PARJ(124))**2)
          SFI=SFW*(1D0-(PARJ(123)/P(ID,5))**2)
          VE=4D0*XW-1D0
          HF1I=SFI*SFF*(VE*POLL+PARJ(132)-PARJ(131))
          HF1W=SFW*SFF**2*((VE**2+1D0)*POLL+2D0*VE*
     &    (PARJ(132)-PARJ(131)))
          KFLC=IABS(KFL1(1))
          PMQ=PYMASS(KFLC)
          QF=KCHG(KFLC,1)/3D0
          VQ=1D0
          IF(MOD(MSTJ(103),2).EQ.1) VQ=SQRT(MAX(0D0,
     &    1D0-(2D0*PMQ/P(ID,5))**2))
          VF=SIGN(1D0,QF)-4D0*QF*XW
          RFV=0.5D0*VQ*(3D0-VQ**2)*(QF**2*POLL-2D0*QF*VF*HF1I+
     &    VF**2*HF1W)+VQ**3*HF1W
          IF(RFV.GT.0D0) PARJ(171)=MIN(1D0,VQ**3*HF1W/RFV)
        ENDIF
 
C...Choice of jet configuration.
        CALL PYXJET(P(ID,5),NJET,CUT)
        KFLC=IABS(KFL1(1))
        KFLN=21
        IF(NJET.EQ.4) THEN
          CALL PYX4JT(NJET,CUT,KFLC,P(ID,5),KFLN,X1,X2,X4,X12,X14)
        ELSEIF(NJET.EQ.3) THEN
          CALL PYX3JT(NJET,CUT,KFLC,P(ID,5),X1,X3)
        ELSE
          MSTJ(120)=1
        ENDIF
 
C...Fill jet configuration; return if incorrect kinematics.
        NC=N-2
        IF(NJET.EQ.2.AND.MSTJ(101).NE.5) THEN
          CALL PY2ENT(NC+1,KFLC,-KFLC,P(ID,5))
        ELSEIF(NJET.EQ.2) THEN
          CALL PY2ENT(-(NC+1),KFLC,-KFLC,P(ID,5))
        ELSEIF(NJET.EQ.3) THEN
          CALL PY3ENT(NC+1,KFLC,21,-KFLC,P(ID,5),X1,X3)
        ELSEIF(KFLN.EQ.21) THEN
          CALL PY4ENT(NC+1,KFLC,KFLN,KFLN,-KFLC,P(ID,5),X1,X2,X4,
     &    X12,X14)
        ELSE
          CALL PY4ENT(NC+1,KFLC,-KFLN,KFLN,-KFLC,P(ID,5),X1,X2,X4,
     &    X12,X14)
        ENDIF
        IF(MSTU(24).NE.0) THEN
          MINT(51)=1
          MSTU(111)=MST111
          PARU(112)=PAR112
          GOTO 720
        ENDIF
 
C...Angular orientation according to matrix element.
        IF(MSTJ(106).EQ.1) THEN
          CALL PYXDIF(NC,NJET,KFLC,P(ID,5),CHIZ,THEZ,PHIZ)
          IF(MINT(11).LT.0) THEZ=PARU(1)-THEZ
          CTHE(1)=COS(THEZ)
          CALL PYROBO(NC+1,N,0D0,CHIZ,0D0,0D0,0D0)
          CALL PYROBO(NC+1,N,THEZ,PHIZ,0D0,0D0,0D0)
        ENDIF
 
C...Boost partons to Z0 rest frame.
        CALL PYROBO(NC+1,N,0D0,0D0,P(ID,1)/P(ID,4),
     &  P(ID,2)/P(ID,4),P(ID,3)/P(ID,4))
 
C...Mark decayed resonance and add documentation lines,
        K(ID,1)=K(ID,1)+10
        IDOC=MINT(83)+MINT(4)
        DO 360 I=NC+1,N
          I1=MINT(83)+MINT(4)+1
          K(I,3)=I1
          IF(MSTP(128).GE.1) K(I,3)=ID
          IF(MSTP(128).LE.1.AND.MINT(4).LT.MSTP(126)) THEN
            MINT(4)=MINT(4)+1
            K(I1,1)=21
            K(I1,2)=K(I,2)
            K(I1,3)=IREF(IP,4)
            DO 350 J=1,5
              P(I1,J)=P(I,J)
  350       CONTINUE
          ENDIF
  360   CONTINUE
 
C...Generate parton shower.
        IF(MSTJ(101).EQ.5.AND.MINT(35).LE.1) THEN
          CALL PYSHOW(N-1,N,P(ID,5))
        ELSEIF(MSTJ(101).EQ.5.AND.MINT(35).GE.2) THEN
          NPART=2
          IPART(1)=N-1
          IPART(2)=N
          PTPART(1)=0.5D0*P(ID,5)
          PTPART(2)=PTPART(1)
          NCT=NCT+1
          IF(K(N-1,2).GT.0) THEN
            MCT(N-1,1)=NCT
            MCT(N,2)=NCT
          ELSE
            MCT(N-1,2)=NCT
            MCT(N,1)=NCT
          ENDIF
          CALL PYPTFS(2,0.5D0*P(ID,5),0D0,PTGEN)
        ENDIF
 
C... End special case for Z0: skip ahead.
        MSTU(111)=MST111
        PARU(112)=PAR112
        GOTO 700
      ENDIF
 
C...Order incoming partons and outgoing resonances.
      IF(JTMAX.EQ.2.AND.ISUB.NE.0.AND.MSTP(47).GE.1.AND.
     &NINH.EQ.0) THEN
        ILIN(1)=MINT(84)+1
        IF(K(MINT(84)+1,2).GT.0) ILIN(1)=MINT(84)+2
        IF(K(ILIN(1),2).EQ.21.OR.K(ILIN(1),2).EQ.22)
     &  ILIN(1)=2*MINT(84)+3-ILIN(1)
        ILIN(2)=2*MINT(84)+3-ILIN(1)
        IMIN=1
        IF(IREF(IP,7).EQ.25.OR.IREF(IP,7).EQ.35.OR.IREF(IP,7)
     &  .EQ.36) IMIN=3
        IMAX=2
        IORD=1
        IF(K(IREF(IP,1),2).EQ.23) IORD=2
        IF(K(IREF(IP,1),2).EQ.24.AND.K(IREF(IP,2),2).EQ.-24) IORD=2
        IAKIPD=IABS(K(IREF(IP,IORD),2))
        IF(IAKIPD.EQ.25.OR.IAKIPD.EQ.35.OR.IAKIPD.EQ.36) IORD=3-IORD
        IF(KDCY(IORD).EQ.0) IORD=3-IORD
 
C...Order decay products of resonances.
        DO 370 JT=IORD,3-IORD,3-2*IORD
          IF(KDCY(JT).EQ.0) THEN
            ILIN(IMAX+1)=NSD(JT)
            IMAX=IMAX+1
          ELSEIF(K(NSD(JT)+1,2).GT.0) THEN
            ILIN(IMAX+1)=N+2*JT-1
            ILIN(IMAX+2)=N+2*JT
            IMAX=IMAX+2
            K(N+2*JT-1,2)=K(NSD(JT)+1,2)
            K(N+2*JT,2)=K(NSD(JT)+2,2)
          ELSE
            ILIN(IMAX+1)=N+2*JT
 
            ILIN(IMAX+2)=N+2*JT-1
            IMAX=IMAX+2
            K(N+2*JT-1,2)=K(NSD(JT)+1,2)
            K(N+2*JT,2)=K(NSD(JT)+2,2)
          ENDIF
  370   CONTINUE
 
C...Find charge, isospin, left- and righthanded couplings.
        DO 390 I=IMIN,IMAX
          DO 380 J=1,4
            COUP(I,J)=0D0
  380     CONTINUE
          KFA=IABS(K(ILIN(I),2))
          IF(KFA.EQ.0.OR.KFA.GT.20) GOTO 390
          COUP(I,1)=KCHG(KFA,1)/3D0
          COUP(I,2)=(-1)**MOD(KFA,2)
          COUP(I,4)=-2D0*COUP(I,1)*XWV
          COUP(I,3)=COUP(I,2)+COUP(I,4)
  390   CONTINUE
 
C...Full propagator dependence and flavour correlations for 2 gamma*/Z.
        IF(ISUB.EQ.22) THEN
          DO 420 I=3,5,2
            I1=IORD
            IF(I.EQ.5) I1=3-IORD
            DO 410 J1=1,2
              DO 400 J2=1,2
                CORL(I/2,J1,J2)=COUP(1,1)**2*HGZ(I1,1)*COUP(I,1)**2/
     &          16D0+COUP(1,1)*COUP(1,J1+2)*HGZ(I1,2)*COUP(I,1)*
     &          COUP(I,J2+2)/4D0+COUP(1,J1+2)**2*HGZ(I1,3)*
     &          COUP(I,J2+2)**2
  400         CONTINUE
  410       CONTINUE
  420     CONTINUE
          COWT12=(CORL(1,1,1)+CORL(1,1,2))*(CORL(2,1,1)+CORL(2,1,2))+
     &    (CORL(1,2,1)+CORL(1,2,2))*(CORL(2,2,1)+CORL(2,2,2))
          COMX12=(CORL(1,1,1)+CORL(1,1,2)+CORL(1,2,1)+CORL(1,2,2))*
     &    (CORL(2,1,1)+CORL(2,1,2)+CORL(2,2,1)+CORL(2,2,2))
 
          IF(COWT12.LT.PYR(0)*COMX12) GOTO 180
        ENDIF
      ENDIF
 
C...Select angular orientation type - Z'/W' only.
      MZPWP=0
      IF(ISUB.EQ.141) THEN
        IF(PYR(0).LT.PARU(130)) MZPWP=1
        IF(IP.EQ.2) THEN
          IF(IABS(K(IREF(2,1),2)).EQ.37) MZPWP=2
          IAKIR=IABS(K(IREF(2,2),2))
          IF(IAKIR.EQ.25.OR.IAKIR.EQ.35.OR.IAKIR.EQ.36) MZPWP=2
          IF(IAKIR.LE.20) MZPWP=2
        ENDIF
        IF(IP.GE.3) MZPWP=2
      ELSEIF(ISUB.EQ.142) THEN
        IF(PYR(0).LT.PARU(136)) MZPWP=1
        IF(IP.EQ.2) THEN
          IAKIR=IABS(K(IREF(2,2),2))
          IF(IAKIR.EQ.25.OR.IAKIR.EQ.35.OR.IAKIR.EQ.36) MZPWP=2
          IF(IAKIR.LE.20) MZPWP=2
        ENDIF
        IF(IP.GE.3) MZPWP=2
      ENDIF
 
C...Select random angles (begin of weighting procedure).
  430 DO 440 JT=1,JTMAX
        IF(KDCY(JT).EQ.0) GOTO 440
        IF(JTMAX.EQ.1.AND.ISUB.NE.0.AND.IHDEC.EQ.0) THEN
          CTHE(JT)=VINT(13)+(VINT(33)-VINT(13)+VINT(34)-VINT(14))*PYR(0)
          IF(CTHE(JT).GT.VINT(33)) CTHE(JT)=CTHE(JT)+VINT(14)-VINT(33)
          PHI(JT)=VINT(24)
        ELSE
          CTHE(JT)=2D0*PYR(0)-1D0
          PHI(JT)=PARU(2)*PYR(0)
        ENDIF
  440 CONTINUE
 
      IF(JTMAX.EQ.2.AND.MSTP(47).GE.1.AND.NINH.EQ.0) THEN
C...Construct massless four-vectors.
        DO 460 I=N+1,N+4
          K(I,1)=1
          DO 450 J=1,5
            P(I,J)=0D0
            V(I,J)=0D0
  450     CONTINUE
  460   CONTINUE
        DO 470 JT=1,JTMAX
          IF(KDCY(JT).EQ.0) GOTO 470
          ID=IREF(IP,JT)
          P(N+2*JT-1,3)=0.5D0*P(ID,5)
          P(N+2*JT-1,4)=0.5D0*P(ID,5)
          P(N+2*JT,3)=-0.5D0*P(ID,5)
          P(N+2*JT,4)=0.5D0*P(ID,5)
          CALL PYROBO(N+2*JT-1,N+2*JT,ACOS(CTHE(JT)),PHI(JT),
     &    P(ID,1)/P(ID,4),P(ID,2)/P(ID,4),P(ID,3)/P(ID,4))
  470   CONTINUE
 
C...Store incoming and outgoing momenta, with random rotation to
C...avoid accidental zeroes in HA expressions.
        IF(ISUB.NE.0) THEN
          DO 490 I=IMIN,IMAX
            K(N+4+I,1)=1
            P(N+4+I,4)=SQRT(P(ILIN(I),1)**2+P(ILIN(I),2)**2+
     &      P(ILIN(I),3)**2+P(ILIN(I),5)**2)
            P(N+4+I,5)=P(ILIN(I),5)
            DO 480 J=1,3
              P(N+4+I,J)=P(ILIN(I),J)
  480       CONTINUE
  490     CONTINUE
  500     THERR=ACOS(2D0*PYR(0)-1D0)
          PHIRR=PARU(2)*PYR(0)
          CALL PYROBO(N+4+IMIN,N+4+IMAX,THERR,PHIRR,0D0,0D0,0D0)
          DO 520 I=IMIN,IMAX
            IF(P(N+4+I,1)**2+P(N+4+I,2)**2.LT.1D-4*(P(N+4+I,1)**2+
     &      P(N+4+I,2)**2+P(N+4+I,3)**2)) GOTO 500
            DO 510 J=1,4
              PK(I,J)=P(N+4+I,J)
  510       CONTINUE
  520     CONTINUE
        ENDIF
 
C...Calculate internal products.
        IF(ISUB.EQ.22.OR.ISUB.EQ.23.OR.ISUB.EQ.25.OR.ISUB.EQ.141.OR.
     &  ISUB.EQ.142) THEN
          DO 540 I1=IMIN,IMAX-1
            DO 530 I2=I1+1,IMAX
              HA(I1,I2)=SNGL(SQRT((PK(I1,4)-PK(I1,3))*(PK(I2,4)+
     &        PK(I2,3))/(1D-20+PK(I1,1)**2+PK(I1,2)**2)))*
     &        CMPLX(SNGL(PK(I1,1)),SNGL(PK(I1,2)))-
     &        SNGL(SQRT((PK(I1,4)+PK(I1,3))*(PK(I2,4)-PK(I2,3))/
     &        (1D-20+PK(I2,1)**2+PK(I2,2)**2)))*
     &        CMPLX(SNGL(PK(I2,1)),SNGL(PK(I2,2)))
              HC(I1,I2)=CONJG(HA(I1,I2))
              IF(I1.LE.2) HA(I1,I2)=CMPLX(0.,1.)*HA(I1,I2)
              IF(I1.LE.2) HC(I1,I2)=CMPLX(0.,1.)*HC(I1,I2)
              HA(I2,I1)=-HA(I1,I2)
              HC(I2,I1)=-HC(I1,I2)
  530       CONTINUE
  540     CONTINUE
        ENDIF
 
C...Calculate four-products.
        IF(ISUB.NE.0) THEN
          DO 560 I=1,2
            DO 550 J=1,4
              PK(I,J)=-PK(I,J)
  550       CONTINUE
  560     CONTINUE
          DO 580 I1=IMIN,IMAX-1
            DO 570 I2=I1+1,IMAX
              PKK(I1,I2)=2D0*(PK(I1,4)*PK(I2,4)-PK(I1,1)*PK(I2,1)-
     &        PK(I1,2)*PK(I2,2)-PK(I1,3)*PK(I2,3))
              PKK(I2,I1)=PKK(I1,I2)
  570       CONTINUE
  580     CONTINUE
        ENDIF
      ENDIF
 
      KFAGM=IABS(IREF(IP,7))
      IF(MSTP(47).LE.0.OR.NINH.NE.0) THEN
C...Isotropic decay selected by user.
        WT=1D0
        WTMAX=1D0
 
      ELSEIF(JTMAX.EQ.3) THEN
C...Isotropic decay when three mother particles.
        WT=1D0
        WTMAX=1D0
 
      ELSEIF(IT4.GE.1) THEN
C... Isotropic decay t -> b + W etc for 4th generation q and l.
        WT=1D0
        WTMAX=1D0
 
      ELSEIF(IREF(IP,7).EQ.25.OR.IREF(IP,7).EQ.35.OR.
     &  IREF(IP,7).EQ.36) THEN
C...Angular weight for h0/A0 -> Z0 + Z0 or W+ + W- -> 4 quarks/leptons.
C...CP-odd case added by Kari Ertresvag Myklevoll.
C...Now also with mixed Higgs CP-states
        ETA=PARP(25)
        IF(IP.EQ.1) WTMAX=SH**2
        IF(IP.GE.2) WTMAX=P(IREF(IP,8),5)**4
        KFA=IABS(K(IREF(IP,1),2))
        KFT=IABS(K(IREF(IP,2),2))
        
        IF((KFA.EQ.KFT).AND.(KFA.EQ.23.OR.KFA.EQ.24).AND.
     &  MSTP(25).GE.3) THEN
C...For mixed CP states need epsilon product.
          P10=PK(3,4)
          P20=PK(4,4)
          P30=PK(5,4)
          P40=PK(6,4)
          P11=PK(3,1)
          P21=PK(4,1)
          P31=PK(5,1)
          P41=PK(6,1)
          P12=PK(3,2)
          P22=PK(4,2)
          P32=PK(5,2)
          P42=PK(6,2)
          P13=PK(3,3)
          P23=PK(4,3)
          P33=PK(5,3)
          P43=PK(6,3)
          EPSI=P10*P21*P32*P43-P10*P21*P33*P42-P10*P22*P31*P43+P10*P22*
     &      P33*P41+P10*P23*P31*P42-P10*P23*P32*P41-P11*P20*P32*P43+P11*
     &      P20*P33*P42+P11*P22*P30*P43-P11*P22*P33*P40-P11*P23*P30*P42+
     &      P11*P23*P32*P40+P12*P20*P31*P43-P12*P20*P33*P41-P12*P21*P30*
     &      P43+P12*P21*P33*P40+P12*P23*P30*P41-P12*P23*P31*P40-P13*P20*
     &      P31*P42+P13*P20*P32*P41+P13*P21*P30*P42-P13*P21*P32*P40-P13*
     &      P22*P30*P41+P13*P22*P31*P40
C...For mixed CP states need gauge boson masses.
          XMA=SQRT(MAX(0D0,(PK(3,4)+PK(4,4))**2-(PK(3,1)+PK(4,1))**2-
     &      (PK(3,2)+PK(4,2))**2-(PK(3,3)+PK(4,3))**2))
          XMB=SQRT(MAX(0D0,(PK(5,4)+PK(6,4))**2-(PK(5,1)+PK(6,1))**2-
     &      (PK(5,2)+PK(6,2))**2-(PK(5,3)+PK(6,3))**2))
          XMV=PMAS(KFA,1)
        ENDIF
 
C...Z decay
        IF(KFA.EQ.23.AND.KFA.EQ.KFT) THEN
          KFLF1A=IABS(KFL1(1))
          EF1=KCHG(KFLF1A,1)/3D0
          AF1=SIGN(1D0,EF1+0.1D0)
          VF1=AF1-4D0*EF1*XWV
          KFLF2A=IABS(KFL1(2))
          EF2=KCHG(KFLF2A,1)/3D0
          AF2=SIGN(1D0,EF2+0.1D0)
          VF2=AF2-4D0*EF2*XWV
          VA12AS=4D0*VF1*AF1*VF2*AF2/((VF1**2+AF1**2)*(VF2**2+AF2**2))
          IF((MSTP(25).EQ.0.AND.IREF(IP,7).NE.36).OR.MSTP(25).EQ.1)
     &      THEN
C...CP-even decay
            WT=8D0*(1D0+VA12AS)*PKK(3,5)*PKK(4,6)+
     &        8D0*(1D0-VA12AS)*PKK(3,6)*PKK(4,5)
          ELSEIF(MSTP(25).LE.2) THEN
C...CP-odd decay
            WT=((PKK(3,5)+PKK(4,6))**2 +(PKK(3,6)+PKK(4,5))**2
     &        -2*PKK(3,4)*PKK(5,6)
     &        -2*(PKK(3,5)*PKK(4,6)-PKK(3,6)*PKK(4,5))**2/
     &        (PKK(3,4)*PKK(5,6))
     &        +VA12AS*(PKK(3,5)+PKK(3,6)-PKK(4,5)-PKK(4,6))*
     &        (PKK(3,5)+PKK(4,5)-PKK(3,6)-PKK(4,6)))/(1+VA12AS)
          ELSE
C...Mixed CP states.
            WT=32D0*(0.25D0*((1D0+VA12AS)*PKK(3,5)*PKK(4,6)
     &        +(1D0-VA12AS)*PKK(3,6)*PKK(4,5))
     &        -0.5D0*ETA/XMV**2*EPSI*((1D0+VA12AS)*(PKK(3,5)+PKK(4,6))
     &        -(1D0-VA12AS)*(PKK(3,6)+PKK(4,5)))
     &        +6.25D-2*ETA**2/XMV**4*(-2D0*PKK(3,4)**2*PKK(5,6)**2
     &        -2D0*(PKK(3,5)*PKK(4,6)-PKK(3,6)*PKK(4,5))**2
     &        +PKK(3,4)*PKK(5,6)
     &        *((PKK(3,5)+PKK(4,6))**2+(PKK(3,6)+PKK(4,5))**2)
     &        +VA12AS*PKK(3,4)*PKK(5,6)
     &        *(PKK(3,5)+PKK(3,6)-PKK(4,5)-PKK(4,6))
     &        *(PKK(3,5)-PKK(3,6)+PKK(4,5)-PKK(4,6))))
     &        /(1D0 +2D0*ETA*XMA*XMB/XMV**2
     &          +2D0*(ETA*XMA*XMB/XMV**2)**2*(1D0+VA12AS))
          ENDIF
 
C...W decay
        ELSEIF(KFA.EQ.24.AND.KFA.EQ.KFT) THEN
          IF((MSTP(25).EQ.0.AND.IREF(IP,7).NE.36).OR.MSTP(25).EQ.1)
     &      THEN
C...CP-even decay
            WT=16D0*PKK(3,5)*PKK(4,6)
          ELSEIF(MSTP(25).LE.2) THEN
C...CP-odd decay
            WT=0.5D0*((PKK(3,5)+PKK(4,6))**2 +(PKK(3,6)+PKK(4,5))**2
     &        -2*PKK(3,4)*PKK(5,6)
     &        -2*(PKK(3,5)*PKK(4,6)-PKK(3,6)*PKK(4,5))**2/
     &        (PKK(3,4)*PKK(5,6))
     &        +(PKK(3,5)+PKK(3,6)-PKK(4,5)-PKK(4,6))*
     &        (PKK(3,5)+PKK(4,5)-PKK(3,6)-PKK(4,6)))
          ELSE
C...Mixed CP states.
            WT=32D0*(0.25D0*2D0*PKK(3,5)*PKK(4,6)
     &        -0.5D0*ETA/XMV**2*EPSI*2D0*(PKK(3,5)+PKK(4,6))
     &        +6.25D-2*ETA**2/XMV**4*(-2D0*PKK(3,4)**2*PKK(5,6)**2
     &        -2D0*(PKK(3,5)*PKK(4,6)-PKK(3,6)*PKK(4,5))**2
     &        +PKK(3,4)*PKK(5,6)
     &        *((PKK(3,5)+PKK(4,6))**2+(PKK(3,6)+PKK(4,5))**2)
     &        +PKK(3,4)*PKK(5,6)
     &        *(PKK(3,5)+PKK(3,6)-PKK(4,5)-PKK(4,6))
     &        *(PKK(3,5)-PKK(3,6)+PKK(4,5)-PKK(4,6))))
     &        /(1D0 +2D0*ETA*XMA*XMB/XMV**2
     &          +(2D0*ETA*XMA*XMB/XMV**2)**2)
          ENDIF
 
C...No angular correlations in other Higgs decays.
        ELSE
          WT=WTMAX
        ENDIF
 
      ELSEIF((KFAGM.EQ.6.OR.KFAGM.EQ.7.OR.KFAGM.EQ.8.OR.
     &  KFAGM.EQ.17.OR.KFAGM.EQ.18).AND.IABS(K(IREF(IP,1),2)).EQ.24)
     &  THEN
C...Angular correlation in f -> f' + W -> f' + 2 quarks/leptons.
        I1=IREF(IP,8)
        IF(MOD(KFAGM,2).EQ.0) THEN
          I2=N+1
          I3=N+2
        ELSE
          I2=N+2
          I3=N+1
        ENDIF
        I4=IREF(IP,2)
        WT=(P(I1,4)*P(I2,4)-P(I1,1)*P(I2,1)-P(I1,2)*P(I2,2)-
     &  P(I1,3)*P(I2,3))*(P(I3,4)*P(I4,4)-P(I3,1)*P(I4,1)-
     &  P(I3,2)*P(I4,2)-P(I3,3)*P(I4,3))
        WTMAX=(P(I1,5)**4-P(IREF(IP,1),5)**4)/8D0
 
      ELSEIF(ISUB.EQ.1) THEN
C...Angular weight for gamma*/Z0 -> 2 quarks/leptons.
        EI=KCHG(IABS(MINT(15)),1)/3D0
        AI=SIGN(1D0,EI+0.1D0)
        VI=AI-4D0*EI*XWV
        EF=KCHG(IABS(KFL1(1)),1)/3D0
        AF=SIGN(1D0,EF+0.1D0)
 
        VF=AF-4D0*EF*XWV
        RMF=MIN(1D0,4D0*PMAS(IABS(KFL1(1)),1)**2/SH)
        WT1=EI**2*VINT(111)*EF**2+EI*VI*VINT(112)*EF*VF+
     &  (VI**2+AI**2)*VINT(114)*(VF**2+(1D0-RMF)*AF**2)
        WT2=RMF*(EI**2*VINT(111)*EF**2+EI*VI*VINT(112)*EF*VF+
     &  (VI**2+AI**2)*VINT(114)*VF**2)
        WT3=SQRT(1D0-RMF)*(EI*AI*VINT(112)*EF*AF+
     &  4D0*VI*AI*VINT(114)*VF*AF)
        WT=WT1*(1D0+CTHE(1)**2)+WT2*(1D0-CTHE(1)**2)+
     &  2D0*WT3*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1))
        WTMAX=2D0*(WT1+ABS(WT3))
 
      ELSEIF(ISUB.EQ.2) THEN
C...Angular weight for W+/- -> 2 quarks/leptons.
        RM3=PMAS(IABS(KFL1(1)),1)**2/SH
        RM4=PMAS(IABS(KFL2(1)),1)**2/SH
        BE34=SQRT(MAX(0D0,(1D0-RM3-RM4)**2-4D0*RM3*RM4))
        WT=(1D0+BE34*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1)))**2-(RM3-RM4)**2
        WTMAX=4D0
 
      ELSEIF(ISUB.EQ.15.OR.ISUB.EQ.19) THEN
C...Angular weight for f + fbar -> gluon/gamma + (gamma*/Z0) ->
C...-> gluon/gamma + 2 quarks/leptons.
        CLILF=COUP(1,1)**2*HGZ(JTZ,1)*COUP(3,1)**2/16D0+
     &  COUP(1,1)*COUP(1,3)*HGZ(JTZ,2)*COUP(3,1)*COUP(3,3)/4D0+
     &  COUP(1,3)**2*HGZ(JTZ,3)*COUP(3,3)**2
        CLIRF=COUP(1,1)**2*HGZ(JTZ,1)*COUP(3,1)**2/16D0+
     &  COUP(1,1)*COUP(1,3)*HGZ(JTZ,2)*COUP(3,1)*COUP(3,4)/4D0+
     &  COUP(1,3)**2*HGZ(JTZ,3)*COUP(3,4)**2
        CRILF=COUP(1,1)**2*HGZ(JTZ,1)*COUP(3,1)**2/16D0+
     &  COUP(1,1)*COUP(1,4)*HGZ(JTZ,2)*COUP(3,1)*COUP(3,3)/4D0+
     &  COUP(1,4)**2*HGZ(JTZ,3)*COUP(3,3)**2
        CRIRF=COUP(1,1)**2*HGZ(JTZ,1)*COUP(3,1)**2/16D0+
     &  COUP(1,1)*COUP(1,4)*HGZ(JTZ,2)*COUP(3,1)*COUP(3,4)/4D0+
     &  COUP(1,4)**2*HGZ(JTZ,3)*COUP(3,4)**2
        WT=(CLILF+CRIRF)*(PKK(1,3)**2+PKK(2,4)**2)+
     &  (CLIRF+CRILF)*(PKK(1,4)**2+PKK(2,3)**2)
        WTMAX=(CLILF+CLIRF+CRILF+CRIRF)*
     &  ((PKK(1,3)+PKK(1,4))**2+(PKK(2,3)+PKK(2,4))**2)
 
      ELSEIF(ISUB.EQ.16.OR.ISUB.EQ.20) THEN
C...Angular weight for f + fbar' -> gluon/gamma + W+/- ->
C...-> gluon/gamma + 2 quarks/leptons.
        WT=PKK(1,3)**2+PKK(2,4)**2
        WTMAX=(PKK(1,3)+PKK(1,4))**2+(PKK(2,3)+PKK(2,4))**2
 
      ELSEIF(ISUB.EQ.22) THEN
C...Angular weight for f + fbar -> Z0 + Z0 -> 4 quarks/leptons.
        S34=P(IREF(IP,IORD),5)**2
        S56=P(IREF(IP,3-IORD),5)**2
        TI=PKK(1,3)+PKK(1,4)+S34
        UI=PKK(1,5)+PKK(1,6)+S56
        TIR=REAL(TI)
        UIR=REAL(UI)
        FGK135=ABS(FGK(1,2,3,4,5,6)/TIR+FGK(1,2,5,6,3,4)/UIR)**2
        FGK145=ABS(FGK(1,2,4,3,5,6)/TIR+FGK(1,2,5,6,4,3)/UIR)**2
        FGK136=ABS(FGK(1,2,3,4,6,5)/TIR+FGK(1,2,6,5,3,4)/UIR)**2
        FGK146=ABS(FGK(1,2,4,3,6,5)/TIR+FGK(1,2,6,5,4,3)/UIR)**2
        FGK253=ABS(FGK(2,1,5,6,3,4)/TIR+FGK(2,1,3,4,5,6)/UIR)**2
        FGK263=ABS(FGK(2,1,6,5,3,4)/TIR+FGK(2,1,3,4,6,5)/UIR)**2
        FGK254=ABS(FGK(2,1,5,6,4,3)/TIR+FGK(2,1,4,3,5,6)/UIR)**2
        FGK264=ABS(FGK(2,1,6,5,4,3)/TIR+FGK(2,1,4,3,6,5)/UIR)**2
 
        WT=
     &  CORL(1,1,1)*CORL(2,1,1)*FGK135+CORL(1,1,2)*CORL(2,1,1)*FGK145+
     &  CORL(1,1,1)*CORL(2,1,2)*FGK136+CORL(1,1,2)*CORL(2,1,2)*FGK146+
     &  CORL(1,2,1)*CORL(2,2,1)*FGK253+CORL(1,2,2)*CORL(2,2,1)*FGK263+
     &  CORL(1,2,1)*CORL(2,2,2)*FGK254+CORL(1,2,2)*CORL(2,2,2)*FGK264
        WTMAX=16D0*((CORL(1,1,1)+CORL(1,1,2))*(CORL(2,1,1)+CORL(2,1,2))+
     &  (CORL(1,2,1)+CORL(1,2,2))*(CORL(2,2,1)+CORL(2,2,2)))*S34*S56*
     &  ((TI**2+UI**2+2D0*SH*(S34+S56))/(TI*UI)-S34*S56*(1D0/TI**2+
     &  1D0/UI**2))
 
      ELSEIF(ISUB.EQ.23) THEN
C...Angular weight for f + fbar' -> Z0 + W+/- -> 4 quarks/leptons.
        D34=P(IREF(IP,IORD),5)**2
        D56=P(IREF(IP,3-IORD),5)**2
        DT=PKK(1,3)+PKK(1,4)+D34
        DU=PKK(1,5)+PKK(1,6)+D56
        FACBW=1D0/((SH-SQMW)**2+GMMW**2)
        CAWZ=COUP(2,3)/DT-2D0*XW1*COUP(1,2)*(SH-SQMW)*FACBW
        CBWZ=COUP(1,3)/DU+2D0*XW1*COUP(1,2)*(SH-SQMW)*FACBW
        FGK135=ABS(REAL(CAWZ)*FGK(1,2,3,4,5,6)+
 
     &  REAL(CBWZ)*FGK(1,2,5,6,3,4))
        FGK136=ABS(REAL(CAWZ)*FGK(1,2,3,4,6,5)+
     &  REAL(CBWZ)*FGK(1,2,6,5,3,4))
        WT=(COUP(5,3)*FGK135)**2+(COUP(5,4)*FGK136)**2
        WTMAX=4D0*D34*D56*(COUP(5,3)**2+COUP(5,4)**2)*(CAWZ**2*
     &  DIGK(DT,DU)+CBWZ**2*DIGK(DU,DT)+CAWZ*CBWZ*DJGK(DT,DU))
 
      ELSEIF(ISUB.EQ.24.OR.ISUB.EQ.171.OR.ISUB.EQ.176) THEN
C...Angular weight for f + fbar -> Z0 + h0 -> 2 quarks/leptons + h0
C...(or H0, or A0).
        WT=((COUP(1,3)*COUP(3,3))**2+(COUP(1,4)*COUP(3,4))**2)*
     &  PKK(1,3)*PKK(2,4)+((COUP(1,3)*COUP(3,4))**2+(COUP(1,4)*
     &  COUP(3,3))**2)*PKK(1,4)*PKK(2,3)
        WTMAX=(COUP(1,3)**2+COUP(1,4)**2)*(COUP(3,3)**2+COUP(3,4)**2)*
     &  (PKK(1,3)+PKK(1,4))*(PKK(2,3)+PKK(2,4))
 
      ELSEIF(ISUB.EQ.25) THEN
C...Angular weight for f + fbar -> W+ + W- -> 4 quarks/leptons.
        POLR=(1D0+PARJ(132))*(1D0-PARJ(131))
        POLL=(1D0-PARJ(132))*(1D0+PARJ(131))
        D34=P(IREF(IP,IORD),5)**2
        D56=P(IREF(IP,3-IORD),5)**2
        DT=PKK(1,3)+PKK(1,4)+D34
        DU=PKK(1,5)+PKK(1,6)+D56
        FACBW=1D0/((SH-SQMZ)**2+SQMZ*PMAS(23,2)**2)
        CDWW=(COUP(1,3)*SQMZ*(SH-SQMZ)*FACBW+COUP(1,2))/SH
        CAWW=CDWW+0.5D0*(COUP(1,2)+1D0)/DT
        CBWW=CDWW+0.5D0*(COUP(1,2)-1D0)/DU
        CCWW=COUP(1,4)*SQMZ*(SH-SQMZ)*FACBW/SH
        FGK135=ABS(REAL(CAWW)*FGK(1,2,3,4,5,6)-
     &  REAL(CBWW)*FGK(1,2,5,6,3,4))
        FGK253=ABS(FGK(2,1,5,6,3,4)-FGK(2,1,3,4,5,6))
        IF(MSTP(50).LE.0) THEN
          WT=FGK135**2+(CCWW*FGK253)**2
          WTMAX=4D0*D34*D56*(CAWW**2*DIGK(DT,DU)+CBWW**2*DIGK(DU,DT)-
     &    CAWW*CBWW*DJGK(DT,DU)+CCWW**2*(DIGK(DT,DU)+DIGK(DU,DT)-
     &    DJGK(DT,DU)))
        ELSE
          WT=POLL*FGK135**2+POLR*(CCWW*FGK253)**2
          WTMAX=4D0*D34*D56*(POLL*(CAWW**2*DIGK(DT,DU)+
     &    CBWW**2*DIGK(DU,DT)-CAWW*CBWW*DJGK(DT,DU))+
     &    POLR*CCWW**2*(DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU)))
        ENDIF
 
      ELSEIF(ISUB.EQ.26.OR.ISUB.EQ.172.OR.ISUB.EQ.177) THEN
C...Angular weight for f + fbar' -> W+/- + h0 -> 2 quarks/leptons + h0
C...(or H0, or A0).
        WT=PKK(1,3)*PKK(2,4)
        WTMAX=(PKK(1,3)+PKK(1,4))*(PKK(2,3)+PKK(2,4))
 
      ELSEIF(ISUB.EQ.30.OR.ISUB.EQ.35) THEN
C...Angular weight for f + g/gamma -> f + (gamma*/Z0)
C...-> f + 2 quarks/leptons.
        CLILF=COUP(1,1)**2*HGZ(JTZ,1)*COUP(3,1)**2/16D0+
     &  COUP(1,1)*COUP(1,3)*HGZ(JTZ,2)*COUP(3,1)*COUP(3,3)/4D0+
     &  COUP(1,3)**2*HGZ(JTZ,3)*COUP(3,3)**2
        CLIRF=COUP(1,1)**2*HGZ(JTZ,1)*COUP(3,1)**2/16D0+
     &  COUP(1,1)*COUP(1,3)*HGZ(JTZ,2)*COUP(3,1)*COUP(3,4)/4D0+
     &  COUP(1,3)**2*HGZ(JTZ,3)*COUP(3,4)**2
        CRILF=COUP(1,1)**2*HGZ(JTZ,1)*COUP(3,1)**2/16D0+
     &  COUP(1,1)*COUP(1,4)*HGZ(JTZ,2)*COUP(3,1)*COUP(3,3)/4D0+
     &  COUP(1,4)**2*HGZ(JTZ,3)*COUP(3,3)**2
        CRIRF=COUP(1,1)**2*HGZ(JTZ,1)*COUP(3,1)**2/16D0+
     &  COUP(1,1)*COUP(1,4)*HGZ(JTZ,2)*COUP(3,1)*COUP(3,4)/4D0+
     &  COUP(1,4)**2*HGZ(JTZ,3)*COUP(3,4)**2
        IF(K(ILIN(1),2).GT.0) WT=(CLILF+CRIRF)*(PKK(1,4)**2+
     &  PKK(3,5)**2)+(CLIRF+CRILF)*(PKK(1,3)**2+PKK(4,5)**2)
        IF(K(ILIN(1),2).LT.0) WT=(CLILF+CRIRF)*(PKK(1,3)**2+
     &  PKK(4,5)**2)+(CLIRF+CRILF)*(PKK(1,4)**2+PKK(3,5)**2)
        WTMAX=(CLILF+CLIRF+CRILF+CRIRF)*
     &  ((PKK(1,3)+PKK(1,4))**2+(PKK(3,5)+PKK(4,5))**2)
 
      ELSEIF(ISUB.EQ.31.OR.ISUB.EQ.36) THEN
C...Angular weight for f + g/gamma -> f' + W+/- -> f' + 2 fermions.
        IF(K(ILIN(1),2).GT.0) WT=PKK(1,4)**2+PKK(3,5)**2
        IF(K(ILIN(1),2).LT.0) WT=PKK(1,3)**2+PKK(4,5)**2
        WTMAX=(PKK(1,3)+PKK(1,4))**2+(PKK(3,5)+PKK(4,5))**2
 
      ELSEIF(ISUB.EQ.71.OR.ISUB.EQ.72.OR.ISUB.EQ.73.OR.ISUB.EQ.76.OR.
     &  ISUB.EQ.77) THEN
C...Angular weight for V_L1 + V_L2 -> V_L3 + V_L4 (V = Z/W).
        WT=16D0*PKK(3,5)*PKK(4,6)
        WTMAX=SH**2
 
      ELSEIF(ISUB.EQ.110) THEN
C...Angular weight for f + fbar -> gamma + h0 -> gamma + X is isotropic.
        WT=1D0
        WTMAX=1D0
 
      ELSEIF(ISUB.EQ.141) THEN
C...Special case: if only branching ratios known then isotropic decay.
        IF(MWID(32).EQ.2) THEN
          WT=1D0
          WTMAX=1D0
        ELSEIF(IP.EQ.1.AND.IABS(KFL1(1)).LT.20) THEN 
C...Angular weight for f + fbar -> gamma*/Z0/Z'0 -> 2 quarks/leptons.
C...Couplings of incoming flavour.
          KFAI=IABS(MINT(15))
          EI=KCHG(KFAI,1)/3D0
          AI=SIGN(1D0,EI+0.1D0)
          VI=AI-4D0*EI*XWV
          KFAIC=1
          IF(KFAI.LE.10.AND.MOD(KFAI,2).EQ.0) KFAIC=2
          IF(KFAI.GT.10.AND.MOD(KFAI,2).NE.0) KFAIC=3
          IF(KFAI.GT.10.AND.MOD(KFAI,2).EQ.0) KFAIC=4
          IF(KFAI.LE.2.OR.KFAI.EQ.11.OR.KFAI.EQ.12) THEN
            VPI=PARU(119+2*KFAIC)
            API=PARU(120+2*KFAIC)
          ELSEIF(KFAI.LE.4.OR.KFAI.EQ.13.OR.KFAI.EQ.14) THEN
            VPI=PARJ(178+2*KFAIC)
            API=PARJ(179+2*KFAIC)
          ELSE
            VPI=PARJ(186+2*KFAIC)
            API=PARJ(187+2*KFAIC)
          ENDIF
C...Couplings of final flavour.
          KFAF=IABS(KFL1(1))
          EF=KCHG(KFAF,1)/3D0
          AF=SIGN(1D0,EF+0.1D0)
          VF=AF-4D0*EF*XWV
          KFAFC=1
          IF(KFAF.LE.10.AND.MOD(KFAF,2).EQ.0) KFAFC=2
          IF(KFAF.GT.10.AND.MOD(KFAF,2).NE.0) KFAFC=3
          IF(KFAF.GT.10.AND.MOD(KFAF,2).EQ.0) KFAFC=4
          IF(KFAF.LE.2.OR.KFAF.EQ.11.OR.KFAF.EQ.12) THEN
            VPF=PARU(119+2*KFAFC)
            APF=PARU(120+2*KFAFC)
          ELSEIF(KFAF.LE.4.OR.KFAF.EQ.13.OR.KFAF.EQ.14) THEN
            VPF=PARJ(178+2*KFAFC)
            APF=PARJ(179+2*KFAFC)
          ELSE
            VPF=PARJ(186+2*KFAFC)
            APF=PARJ(187+2*KFAFC)
          ENDIF
C...Asymmetry and weight.
          ASYM=2D0*(EI*AI*VINT(112)*EF*AF+EI*API*VINT(113)*EF*APF+
     &    4D0*VI*AI*VINT(114)*VF*AF+(VI*API+VPI*AI)*VINT(115)*
     &    (VF*APF+VPF*AF)+4D0*VPI*API*VINT(116)*VPF*APF)/
     &    (EI**2*VINT(111)*EF**2+EI*VI*VINT(112)*EF*VF+
     &    EI*VPI*VINT(113)*EF*VPF+(VI**2+AI**2)*VINT(114)*
     &    (VF**2+AF**2)+(VI*VPI+AI*API)*VINT(115)*(VF*VPF+AF*APF)+
     &    (VPI**2+API**2)*VINT(116)*(VPF**2+APF**2))
          WT=1D0+ASYM*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1))+CTHE(1)**2
          WTMAX=2D0+ABS(ASYM)
        ELSEIF(IP.EQ.1.AND.IABS(KFL1(1)).EQ.24) THEN
C...Angular weight for f + fbar -> Z' -> W+ + W-.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          CCOS2=-(1D0/16D0)*((1D0-RM1-RM2)**2-4D0*RM1*RM2)*
     &    (1D0-2D0*RM1-2D0*RM2+RM1**2+RM2**2+10D0*RM1*RM2)
          CFLAT=-CCOS2+0.5D0*(RM1+RM2)*(1D0-2D0*RM1-2D0*RM2+
     &    (RM2-RM1)**2)
          WT=CFLAT+CCOS2*CTHE(1)**2
          WTMAX=CFLAT+MAX(0D0,CCOS2)
        ELSEIF(IP.EQ.1.AND.(KFL1(1).EQ.25.OR.KFL1(1).EQ.35.OR.
     &    IABS(KFL1(1)).EQ.37)) THEN
C...Angular weight for f + fbar -> Z' -> h0 + A0, H0 + A0, H+ + H-.
          WT=1D0-CTHE(1)**2
          WTMAX=1D0
        ELSEIF(IP.EQ.1.AND.KFL2(1).EQ.25) THEN
C...Angular weight for f + fbar -> Z' -> Z0 + h0.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          FLAM2=MAX(0D0,(1D0-RM1-RM2)**2-4D0*RM1*RM2)
          WT=1D0+FLAM2*(1D0-CTHE(1)**2)/(8D0*RM1)
          WTMAX=1D0+FLAM2/(8D0*RM1)
        ELSEIF(MZPWP.EQ.0) THEN
C...Angular weight for f + fbar -> Z' -> W+ + W- -> 4 quarks/leptons
C...(W:s like if intermediate Z).
          D34=P(IREF(IP,IORD),5)**2
          D56=P(IREF(IP,3-IORD),5)**2
          DT=PKK(1,3)+PKK(1,4)+D34
          DU=PKK(1,5)+PKK(1,6)+D56
          FGK135=ABS(FGK(1,2,3,4,5,6)-FGK(1,2,5,6,3,4))
          FGK253=ABS(FGK(2,1,5,6,3,4)-FGK(2,1,3,4,5,6))
          WT=(COUP(1,3)*FGK135)**2+(COUP(1,4)*FGK253)**2
          WTMAX=4D0*D34*D56*(COUP(1,3)**2+COUP(1,4)**2)*
     &    (DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU))
        ELSEIF(MZPWP.EQ.1) THEN
C...Angular weight for f + fbar -> Z' -> W+ + W- -> 4 quarks/leptons
C...(W:s approximately longitudinal, like if intermediate H).
          WT=16D0*PKK(3,5)*PKK(4,6)
          WTMAX=SH**2
        ELSE
C...Angular weight for f + fbar -> Z' -> H+ + H-, Z0 + h0, h0 + A0,
C...H0 + A0 -> 4 quarks/leptons, t + tbar -> b + W+ + bbar + W- .
          WT=1D0
          WTMAX=1D0
        ENDIF
 
      ELSEIF(ISUB.EQ.142) THEN
C...Special case: if only branching ratios known then isotropic decay.
        IF(MWID(34).EQ.2) THEN
          WT=1D0
          WTMAX=1D0
        ELSEIF(IP.EQ.1.AND.IABS(KFL1(1)).LT.20) THEN 
C...Angular weight for f + fbar' -> W'+/- -> 2 quarks/leptons.
          KFAI=IABS(MINT(15))
          KFAIC=1
          IF(KFAI.GT.10) KFAIC=2
          VI=PARU(129+2*KFAIC)
          AI=PARU(130+2*KFAIC)
          KFAF=IABS(KFL1(1))
          KFAFC=1
          IF(KFAF.GT.10) KFAFC=2
          VF=PARU(129+2*KFAFC)
          AF=PARU(130+2*KFAFC)
          ASYM=8D0*VI*AI*VF*AF/((VI**2+AI**2)*(VF**2+AF**2))
          WT=1D0+ASYM*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1))+CTHE(1)**2
          WTMAX=2D0+ABS(ASYM)
        ELSEIF(IP.EQ.1.AND.IABS(KFL2(1)).EQ.23) THEN
C...Angular weight for f + fbar' -> W'+/- -> W+/- + Z0.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          CCOS2=-(1D0/16D0)*((1D0-RM1-RM2)**2-4D0*RM1*RM2)*
     &    (1D0-2D0*RM1-2D0*RM2+RM1**2+RM2**2+10D0*RM1*RM2)
          CFLAT=-CCOS2+0.5D0*(RM1+RM2)*(1D0-2D0*RM1-2D0*RM2+
     &    (RM2-RM1)**2)
          WT=CFLAT+CCOS2*CTHE(1)**2
          WTMAX=CFLAT+MAX(0D0,CCOS2)
        ELSEIF(IP.EQ.1.AND.KFL2(1).EQ.25) THEN
C...Angular weight for f + fbar -> W'+/- -> W+/- + h0.
          RM1=P(NSD(1)+1,5)**2/SH
          RM2=P(NSD(1)+2,5)**2/SH
          FLAM2=MAX(0D0,(1D0-RM1-RM2)**2-4D0*RM1*RM2)
          WT=1D0+FLAM2*(1D0-CTHE(1)**2)/(8D0*RM1)
          WTMAX=1D0+FLAM2/(8D0*RM1)
        ELSEIF(MZPWP.EQ.0) THEN
C...Angular weight for f + fbar' -> W' -> W + Z0 -> 4 quarks/leptons
C...(W/Z like if intermediate W).
          D34=P(IREF(IP,IORD),5)**2
          D56=P(IREF(IP,3-IORD),5)**2
          DT=PKK(1,3)+PKK(1,4)+D34
          DU=PKK(1,5)+PKK(1,6)+D56
          FGK135=ABS(FGK(1,2,3,4,5,6)-FGK(1,2,5,6,3,4))
          FGK136=ABS(FGK(1,2,3,4,6,5)-FGK(1,2,6,5,3,4))
          WT=(COUP(5,3)*FGK135)**2+(COUP(5,4)*FGK136)**2
          WTMAX=4D0*D34*D56*(COUP(5,3)**2+COUP(5,4)**2)*
     &    (DIGK(DT,DU)+DIGK(DU,DT)-DJGK(DT,DU))
        ELSEIF(MZPWP.EQ.1) THEN
C...Angular weight for f + fbar' -> W' -> W + Z0 -> 4 quarks/leptons
C...(W/Z approximately longitudinal, like if intermediate H).
          WT=16D0*PKK(3,5)*PKK(4,6)
          WTMAX=SH**2
        ELSE
C...Angular weight for f + fbar -> W' -> W + h0 -> whatever,
C...t + bbar -> t + W + bbar.
          WT=1D0
          WTMAX=1D0
        ENDIF
 
      ELSEIF(ISUB.EQ.145.OR.ISUB.EQ.162.OR.ISUB.EQ.163.OR.ISUB.EQ.164)
     &  THEN
C...Isotropic decay of leptoquarks (assumed spin 0).
        WT=1D0
        WTMAX=1D0
 
      ELSEIF(ISUB.GE.146.AND.ISUB.LE.148) THEN
C...Decays of (spin 1/2) q*/e* -> q/e + (g,gamma) or (Z0,W+-).
        SIDE=1D0
        IF(MINT(16).EQ.21.OR.MINT(16).EQ.22) SIDE=-1D0
        IF(IP.EQ.1.AND.(KFL1(1).EQ.21.OR.KFL1(1).EQ.22)) THEN
          WT=1D0+SIDE*CTHE(1)
          WTMAX=2D0
        ELSEIF(IP.EQ.1) THEN
 
          RM1=P(NSD(1)+1,5)**2/SH
          WT=1D0+SIDE*CTHE(1)*(1D0-0.5D0*RM1)/(1D0+0.5D0*RM1)
          WTMAX=1D0+(1D0-0.5D0*RM1)/(1D0+0.5D0*RM1)
        ELSE
C...W/Z decay assumed isotropic, since not known.
          WT=1D0
          WTMAX=1D0
        ENDIF
 
      ELSEIF(ISUB.EQ.149) THEN
C...Isotropic decay of techni-eta.
        WT=1D0
        WTMAX=1D0
 
      ELSEIF(ISUB.EQ.191) THEN
        IF(IP.EQ.1.AND.IABS(KFL1(1)).GT.21) THEN
C...Angular weight for f + fbar -> rho_tc0 -> W+ W-,
C...W+ pi_tc-, pi_tc+ W- or pi_tc+ pi_tc-.
          WT=1D0-CTHE(1)**2
          WTMAX=1D0
        ELSEIF(IP.EQ.1) THEN
C...Angular weight for f + fbar -> rho_tc0 -> f fbar.
          CTHESG=CTHE(1)*ISIGN(1,MINT(15))
          XWRHT=(1D0-2D0*XW)/(4D0*XW*(1D0-XW))
          BWZR=XWRHT*SH*(SH-SQMZ)/((SH-SQMZ)**2+GMMZ**2)
          BWZI=XWRHT*SH*GMMZ/((SH-SQMZ)**2+GMMZ**2)
          KFAI=IABS(MINT(15))
          EI=KCHG(KFAI,1)/3D0
          AI=SIGN(1D0,EI+0.1D0)
          VI=AI-4D0*EI*XWV
          VALI=0.5D0*(VI+AI)
          VARI=0.5D0*(VI-AI)
          ALEFTI=(EI+VALI*BWZR)**2+(VALI*BWZI)**2
          ARIGHI=(EI+VARI*BWZR)**2+(VARI*BWZI)**2
          KFAF=IABS(KFL1(1))
          EF=KCHG(KFAF,1)/3D0
          AF=SIGN(1D0,EF+0.1D0)
          VF=AF-4D0*EF*XWV
          VALF=0.5D0*(VF+AF)
          VARF=0.5D0*(VF-AF)
          ALEFTF=(EF+VALF*BWZR)**2+(VALF*BWZI)**2
          ARIGHF=(EF+VARF*BWZR)**2+(VARF*BWZI)**2
          ASAME=ALEFTI*ALEFTF+ARIGHI*ARIGHF
          AFLIP=ALEFTI*ARIGHF+ARIGHI*ALEFTF
          WT=ASAME*(1D0+CTHESG)**2+AFLIP*(1D0-CTHESG)**2
          WTMAX=4D0*MAX(ASAME,AFLIP)
        ELSE
C...Isotropic decay of W/pi_tc produced in rho_tc decay.
          WT=1D0
          WTMAX=1D0
        ENDIF
 
      ELSEIF(ISUB.EQ.192) THEN
        IF(IP.EQ.1.AND.IABS(KFL1(1)).GT.21) THEN
C...Angular weight for f + fbar' -> rho_tc+ -> W+ Z0,
C...W+ pi_tc0, pi_tc+ Z0 or pi_tc+ pi_tc0.
          WT=1D0-CTHE(1)**2
          WTMAX=1D0
        ELSEIF(IP.EQ.1) THEN
C...Angular weight for f + fbar' -> rho_tc+ -> f fbar'.
          CTHESG=CTHE(1)*ISIGN(1,MINT(15))
          WT=(1D0+CTHESG)**2
          WTMAX=4D0
        ELSE
C...Isotropic decay of W/Z/pi_tc produced in rho_tc+ decay.
          WT=1D0
          WTMAX=1D0
        ENDIF
 
      ELSEIF(ISUB.EQ.193) THEN
        IF(IP.EQ.1.AND.IABS(KFL1(1)).GT.21) THEN
C...Angular weight for f + fbar -> omega_tc0 ->
C...gamma pi_tc0 or Z0 pi_tc0.
          WT=1D0+CTHE(1)**2
          WTMAX=2D0
        ELSEIF(IP.EQ.1) THEN
C...Angular weight for f + fbar -> omega_tc0 -> f fbar.
          CTHESG=CTHE(1)*ISIGN(1,MINT(15))
          BWZR=(0.5D0/(1D0-XW))*SH*(SH-SQMZ)/((SH-SQMZ)**2+GMMZ**2)
          BWZI=(0.5D0/(1D0-XW))*SH*GMMZ/((SH-SQMZ)**2+GMMZ**2)
          KFAI=IABS(MINT(15))
          EI=KCHG(KFAI,1)/3D0
          AI=SIGN(1D0,EI+0.1D0)
          VI=AI-4D0*EI*XWV
          VALI=0.5D0*(VI+AI)
          VARI=0.5D0*(VI-AI)
          BLEFTI=(EI-VALI*BWZR)**2+(VALI*BWZI)**2
          BRIGHI=(EI-VARI*BWZR)**2+(VARI*BWZI)**2
          KFAF=IABS(KFL1(1))
          EF=KCHG(KFAF,1)/3D0
          AF=SIGN(1D0,EF+0.1D0)
          VF=AF-4D0*EF*XWV
          VALF=0.5D0*(VF+AF)
          VARF=0.5D0*(VF-AF)
          BLEFTF=(EF-VALF*BWZR)**2+(VALF*BWZI)**2
          BRIGHF=(EF-VARF*BWZR)**2+(VARF*BWZI)**2
          BSAME=BLEFTI*BLEFTF+BRIGHI*BRIGHF
          BFLIP=BLEFTI*BRIGHF+BRIGHI*BLEFTF
          WT=BSAME*(1D0+CTHESG)**2+BFLIP*(1D0-CTHESG)**2
          WTMAX=4D0*MAX(BSAME,BFLIP)
        ELSE
C...Isotropic decay of Z/pi_tc produced in omega_tc decay.
          WT=1D0
          WTMAX=1D0
        ENDIF
 
      ELSEIF(ISUB.EQ.353) THEN
C...Angular weight for Z_R0 -> 2 quarks/leptons.
        EI=KCHG(IABS(MINT(15)),1)/3D0
        AI=SIGN(1D0,EI+0.1D0)
        VI=AI-4D0*EI*XWV
        EF=KCHG(PYCOMP(KFL1(1)),1)/3D0
        AF=SIGN(1D0,EF+0.1D0)
        VF=AF-4D0*EF*XWV
        RMF=MIN(1D0,4D0*PMAS(PYCOMP(KFL1(1)),1)**2/SH)
        WT1=(VI**2+AI**2)*(VF**2+(1D0-RMF)*AF**2)
        WT2=RMF*(VI**2+AI**2)*VF**2
        WT3=SQRT(1D0-RMF)*4D0*VI*AI*VF*AF
        WT=WT1*(1D0+CTHE(1)**2)+WT2*(1D0-CTHE(1)**2)+
     &  2D0*WT3*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1))
        WTMAX=2D0*(WT1+ABS(WT3))
 
      ELSEIF(ISUB.EQ.354) THEN
C...Angular weight for W_R+/- -> 2 quarks/leptons.
        RM3=PMAS(PYCOMP(KFL1(1)),1)**2/SH
        RM4=PMAS(PYCOMP(KFL2(1)),1)**2/SH
        BE34=SQRT(MAX(0D0,(1D0-RM3-RM4)**2-4D0*RM3*RM4))
        WT=(1D0+BE34*CTHE(1)*ISIGN(1,MINT(15)*KFL1(1)))**2-(RM3-RM4)**2
        WTMAX=4D0
 
      ELSEIF(ISUB.EQ.391) THEN
C...Angular weight for f + fbar -> G* -> f + fbar
        IF(IP.EQ.1.AND.IABS(KFL1(1)).LE.18) THEN
          WT=1D0-3D0*CTHE(1)**2+4D0*CTHE(1)**4
          WTMAX=2D0
C...Angular weight for f + fbar -> G* -> gamma + gamma or g + g
C...implemented by M.-C. Lemaire
        ELSEIF(IP.EQ.1.AND.(IABS(KFL1(1)).EQ.21.OR.
     &  IABS(KFL1(1)).EQ.22)) THEN
          WT=1D0-CTHE(1)**4
          WTMAX=1D0
C...Other G* decays not yet implemented angular distributions.
        ELSE
          WT=1D0
          WTMAX=1D0
        ENDIF
 
      ELSEIF(ISUB.EQ.392) THEN
C...Angular weight for g + g -> G* -> f + fbar
        IF(IP.EQ.1.AND.IABS(KFL1(1)).LE.18) THEN
          WT=1D0-CTHE(1)**4
          WTMAX=1D0
C...Angular weight for g + g -> G* -> gamma +gamma or g + g
C...implemented by M.-C. Lemaire
        ELSEIF(IP.EQ.1.AND.(IABS(KFL1(1)).EQ.21.OR.
     &  IABS(KFL1(1)).EQ.22)) THEN
         WT=1D0+6D0*CTHE(1)**2+CTHE(1)**4
          WTMAX=8D0
C...Other G* decays not yet implemented angular distributions.
        ELSE
          WT=1D0
          WTMAX=1D0
        ENDIF
 
C...Obtain correct angular distribution by rejection techniques.
      ELSE
        WT=1D0
        WTMAX=1D0
      ENDIF
      IF(WT.LT.PYR(0)*WTMAX) GOTO 430
  
C...Construct massive four-vectors using angles chosen.
  590 DO 690 JT=1,JTMAX
        IF(KDCY(JT).EQ.0) GOTO 690
        ID=IREF(IP,JT)
        DO 600 J=1,5
          DPMO(J)=P(ID,J)
  600   CONTINUE
        DPMO(4)=SQRT(DPMO(1)**2+DPMO(2)**2+DPMO(3)**2+DPMO(5)**2)
CMRENNA++
        NPROD=2
        IF(KFL3(JT).NE.0) NPROD=3
        IF(KFL4(JT).NE.0) NPROD=4
        CALL PYROBO(NSD(JT)+1,NSD(JT)+NPROD,ACOS(CTHE(JT)),PHI(JT),
     &       DPMO(1)/DPMO(4),DPMO(2)/DPMO(4),DPMO(3)/DPMO(4))
        N0=NSD(JT)+NPROD
 
        DO 610 J=1,4
          VDCY(J)=V(ID,J)+V(ID,5)*P(ID,J)/P(ID,5)
  610   CONTINUE
C...Fill in position of decay vertex.
        DO 630 I=NSD(JT)+1,N0
          DO 620 J=1,4
            V(I,J)=VDCY(J)
  620     CONTINUE
          V(I,5)=0D0
 
  630   CONTINUE
CMRENNA--
 
C...Mark decayed resonances; trace history.
        K(ID,1)=K(ID,1)+10
        KFA=IABS(K(ID,2))
        KCA=PYCOMP(KFA)
        IF(KCQM(JT).NE.0) THEN
C...Do not kill colour flow through coloured resonance!
        ELSE
          K(ID,4)=NSD(JT)+1
          K(ID,5)=NSD(JT)+NPROD
          IF(ITJUNC(JT).NE.0) K(ID,5)=K(ID,5)+1
C...If 3-body or 2-body with junction:
c          IF(KFL3(JT).NE.0.OR.ITJUNC(JT).NE.0) K(ID,5)=NSD(JT)+3
C...If 3-body with junction:
c          IF(ITJUNC(JT).NE.0.AND.KFL3(JT).NE.0) K(ID,5)=NSD(JT)+4
        ENDIF
 
C...Add documentation lines.
        ISUBRG=MAX(1,MIN(500,MINT(1)))
        IF(IRES.EQ.0.OR.ISET(ISUBRG).EQ.11) THEN
          IDOC=MINT(83)+MINT(4)
CMRENNA+++
          IHI=NSD(JT)+NPROD
c          IF(KFL3(JT).NE.0) IHI=IHI+1
          DO 650 I=NSD(JT)+1,IHI
CMRENNA---
            I1=MINT(83)+MINT(4)+1
            K(I,3)=I1
            IF(MSTP(128).GE.1) K(I,3)=ID
            IF(MSTP(128).LE.1.AND.MINT(4).LT.MSTP(126)) THEN
              MINT(4)=MINT(4)+1
              K(I1,1)=21
              K(I1,2)=K(I,2)
              K(I1,3)=IREF(IP,JT+3)
              DO 640 J=1,5
                P(I1,J)=P(I,J)
  640         CONTINUE
            ENDIF
  650     CONTINUE
        ELSE
          K(NSD(JT)+1,3)=ID
          K(NSD(JT)+2,3)=ID
C...If 3-body or 2-body with junction:
          IF(KFL3(JT).NE.0.OR.ITJUNC(JT).GT.0) K(NSD(JT)+3,3)=ID
C...If 3-body with junction:
          IF(KFL3(JT).NE.0.AND.ITJUNC(JT).GT.0) K(NSD(JT)+4,3)=ID
C...If 4-body or 3-body with junction:
          IF(KFL4(JT).NE.0.OR.ITJUNC(JT).GT.0) K(NSD(JT)+4,3)=ID
C...If 4-body with junction:
          IF(KFL4(JT).NE.0.AND.ITJUNC(JT).GT.0) K(NSD(JT)+5,3)=ID
        ENDIF
 
C...Do showering of two or three objects.
        NSHBEF=N
        IF(MSTP(71).GE.1.AND.MINT(35).LE.1) THEN
          IF(KFL3(JT).EQ.0) THEN
            CALL PYSHOW(NSD(JT)+1,NSD(JT)+2,P(ID,5))
          ELSE
            CALL PYSHOW(NSD(JT)+1,-NPROD,P(ID,5))
          ENDIF
 
c...For pT-ordered shower need set up first, especially colour tags.
C...(Need to set up colour tags even if MSTP(71) = 0)
        ELSEIF(MINT(35).GE.2) THEN
          NPART=NPROD
c          IF(KFL3(JT).NE.0) NPART=3
          IPART(1)=NSD(JT)+1
          IPART(2)=NSD(JT)+2
          IPART(3)=NSD(JT)+3
          IPART(4)=NSD(JT)+4
          PTPART(1)=0.5D0*P(ID,5)
          PTPART(2)=PTPART(1)
          PTPART(3)=PTPART(1)
          PTPART(4)=PTPART(1)
          IF(KCQ1(JT).EQ.1.OR.KCQ1(JT).EQ.2) THEN
            MOTHER=K(NSD(JT)+1,4)/MSTU(5)
            IF(MOTHER.LE.NSD(JT)) THEN
              MCT(NSD(JT)+1,1)=MCT(MOTHER,1)
            ELSE
              NCT=NCT+1
              MCT(NSD(JT)+1,1)=NCT
              MCT(MOTHER,2)=NCT
            ENDIF
          ENDIF
          IF(KCQ1(JT).EQ.-1.OR.KCQ1(JT).EQ.2) THEN
            MOTHER=K(NSD(JT)+1,5)/MSTU(5)
            IF(MOTHER.LE.NSD(JT)) THEN
              MCT(NSD(JT)+1,2)=MCT(MOTHER,2)
            ELSE
              NCT=NCT+1
              MCT(NSD(JT)+1,2)=NCT
              MCT(MOTHER,1)=NCT
            ENDIF
          ENDIF
          IF(MCT(NSD(JT)+2,1).EQ.0.AND.(KCQ2(JT).EQ.1.OR.
     &    KCQ2(JT).EQ.2)) THEN
            MOTHER=K(NSD(JT)+2,4)/MSTU(5)
            IF(MOTHER.LE.NSD(JT)) THEN
              MCT(NSD(JT)+2,1)=MCT(MOTHER,1)
            ELSE
              NCT=NCT+1
              MCT(NSD(JT)+2,1)=NCT
              MCT(MOTHER,2)=NCT
            ENDIF
          ENDIF
          IF(MCT(NSD(JT)+2,2).EQ.0.AND.(KCQ2(JT).EQ.-1.OR.
     &    KCQ2(JT).EQ.2)) THEN
            MOTHER=K(NSD(JT)+2,5)/MSTU(5)
            IF(MOTHER.LE.NSD(JT)) THEN
              MCT(NSD(JT)+2,2)=MCT(MOTHER,2)
            ELSE
              NCT=NCT+1
              MCT(NSD(JT)+2,2)=NCT
              MCT(MOTHER,1)=NCT
            ENDIF
          ENDIF
          IF(NPART.EQ.3.AND.MCT(NSD(JT)+3,1).EQ.0.AND.
     &    (KCQ3(JT).EQ.1.OR. KCQ3(JT).EQ.2)) THEN
            MOTHER=K(NSD(JT)+3,4)/MSTU(5)
            MCT(NSD(JT)+3,1)=MCT(MOTHER,1)
          ENDIF
          IF(NPART.EQ.3.AND.MCT(NSD(JT)+3,2).EQ.0.AND.
     &    (KCQ3(JT).EQ.-1.OR.KCQ3(JT).EQ.2)) THEN
            MOTHER=K(NSD(JT)+3,5)/MSTU(5)
            MCT(NSD(JT)+2,2)=MCT(MOTHER,2)
          ENDIF
          IF(NPART.EQ.4.AND.MCT(NSD(JT)+4,1).EQ.0.AND.
     &    (KCQ4(JT).EQ.1.OR. KCQ4(JT).EQ.2)) THEN
            MOTHER=K(NSD(JT)+4,4)/MSTU(5)
            MCT(NSD(JT)+4,1)=MCT(MOTHER,1)
          ENDIF
          IF(NPART.EQ.4.AND.MCT(NSD(JT)+4,2).EQ.0.AND.
     &    (KCQ4(JT).EQ.-1.OR.KCQ4(JT).EQ.2)) THEN
            MOTHER=K(NSD(JT)+4,5)/MSTU(5)
            MCT(NSD(JT)+4,2)=MCT(MOTHER,2)
          ENDIF

          IF (MSTP(71).GE.1) CALL PYPTFS(2,0.5D0*P(ID,5),0D0,PTGEN)
        ENDIF
        NSHAFT=N
        IF(JT.EQ.1) NAFT1=N
 
C...Check if decay products moved by shower.
        NSD1=NSD(JT)+1
        NSD2=NSD(JT)+2
        NSD3=NSD(JT)+3
        NSD4=NSD(JT)+4
C...4-body decays will only work if one of the products is "inert"
        IF(NSHAFT.GT.NSHBEF) THEN
          IF(K(NSD1,1).GT.10) THEN
            DO 660 I=NSHBEF+1,NSHAFT
              IF(K(I,1).LT.10.AND.K(I,2).EQ.K(NSD1,2)) NSD1=I
  660       CONTINUE
          ENDIF
          IF(K(NSD2,1).GT.10) THEN
            DO 670 I=NSHBEF+1,NSHAFT
              IF(K(I,1).LT.10.AND.K(I,2).EQ.K(NSD2,2).AND.
     &        I.NE.NSD1) NSD2=I
  670       CONTINUE
          ENDIF
          IF(KFL3(JT).NE.0.AND.K(NSD3,1).GT.10) THEN
            DO 680 I=NSHBEF+1,NSHAFT
              IF(K(I,1).LT.10.AND.K(I,2).EQ.K(NSD3,2).AND.
     &        I.NE.NSD1.AND.I.NE.NSD2) NSD3=I
  680       CONTINUE
          ENDIF
          IF(KFL4(JT).NE.0.AND.K(NSD4,1).GT.10) THEN
            DO 685 I=NSHBEF+1,NSHAFT
              IF(K(I,1).LT.10.AND.K(I,2).EQ.K(NSD4,2).AND.
     &        I.NE.NSD1.AND.I.NE.NSD2.AND.I.NE.NSD3) NSD4=I
  685       CONTINUE
          ENDIF
        ENDIF
 
C...Store decay products for further treatment.
        IF(KFL4(JT).EQ.0) THEN
          NP=NP+1
          IREF(NP,1)=NSD1
          IREF(NP,2)=NSD2
          IREF(NP,3)=0
          IF(KFL3(JT).NE.0) IREF(NP,3)=NSD3
          IREF(NP,4)=IDOC+1
          IREF(NP,5)=IDOC+2
          IREF(NP,6)=0
          IF(KFL3(JT).NE.0) IREF(NP,6)=IDOC+3
          IREF(NP,7)=K(IREF(IP,JT),2)
          IREF(NP,8)=IREF(IP,JT)
        ELSE
          NSDA=NSD1
          NSDB=NSD2
          NSDC=NSD3
          NP=NP+1
          IREF(NP,4)=IDOC+1
          IREF(NP,5)=IDOC+2
          IREF(NP,6)=IDOC+3
          IF(K(NSD1,1).EQ.1) THEN
            NSDA=NSD4
            IREF(NP,4)=IDOC+4
          ELSEIF(K(NSD2,1).EQ.1) THEN
            NSDB=NSD4
            IREF(NP,5)=IDOC+4
          ELSEIF(K(NSD3,1).EQ.1) THEN
            NSDC=NSD4
            IREF(NP,6)=IDOC+4
          ENDIF
          IREF(NP,1)=NSDA
          IREF(NP,2)=NSDB
          IREF(NP,3)=NSDC
          IREF(NP,7)=K(IREF(IP,JT),2)
          IREF(NP,8)=IREF(IP,JT)
        ENDIF
  690 CONTINUE
 
 
C...Fill information for 2 -> 1 -> 2.
  700 IF(JTMAX.EQ.1.AND.KDCY(1).NE.0.AND.ISUB.NE.0) THEN
        MINT(7)=MINT(83)+6+2*ISET(ISUB)
        MINT(8)=MINT(83)+7+2*ISET(ISUB)
        MINT(25)=KFL1(1)
        MINT(26)=KFL2(1)
        VINT(23)=CTHE(1)
        RM3=P(N-1,5)**2/SH
        RM4=P(N,5)**2/SH
        BE34=SQRT(MAX(0D0,(1D0-RM3-RM4)**2-4D0*RM3*RM4))
        VINT(45)=-0.5D0*SH*(1D0-RM3-RM4-BE34*CTHE(1))
        VINT(46)=-0.5D0*SH*(1D0-RM3-RM4+BE34*CTHE(1))
        VINT(48)=0.25D0*SH*BE34**2*MAX(0D0,1D0-CTHE(1)**2)
        VINT(47)=SQRT(VINT(48))
      ENDIF
 
C...Possibility of colour rearrangement in W+W- events.
      IF((ISUB.EQ.25.OR.ISUB.EQ.22).AND.MSTP(115).GE.1) THEN
        IAKF1=IABS(KFL1(1))
        IAKF2=IABS(KFL1(2))
        IAKF3=IABS(KFL2(1))
        IAKF4=IABS(KFL2(2))
        IF(MIN(IAKF1,IAKF2,IAKF3,IAKF4).GE.1.AND.
     &  MAX(IAKF1,IAKF2,IAKF3,IAKF4).LE.5) CALL
     &  PYRECO(IREF(1,1),IREF(1,2),NSD(1),NAFT1)
        IF(MINT(51).NE.0) RETURN
      ENDIF

C...Loop back if needed.
  710 IF(IP.LT.NP) GOTO 170

C...Boost back to standard frame.
  720 IF(IBST.EQ.1) CALL PYROBO(MINT(83)+7,N,THEIN,PHIIN,BEXIN,BEYIN,
     &BEZIN)

 
      RETURN
      END
