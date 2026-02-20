 
C*********************************************************************
 
C...PYINRE
C...Calculates full and effective widths of gauge bosons, stores
C...masses and widths, rescales coefficients to be used for
C...resonance production generation.
 
      SUBROUTINE PYINRE
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYINT6/PROC(0:500)
      CHARACTER PROC*28
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYDAT4/,/PYSUBS/,/PYPARS/,
     &/PYINT1/,/PYINT2/,/PYINT4/,/PYINT6/,/PYMSSM/
C...Local arrays and data.
      CHARACTER PRTMP*9
      DIMENSION WDTP(0:400),WDTE(0:400,0:5),WDTPM(0:400),
     &WDTEM(0:400,0:5),KCORD(500),PMORD(500)
 
C...Born level couplings in MSSM Higgs doublet sector.
      XW=PARU(102)
      XWV=XW
      IF(MSTP(8).GE.2) XW=1D0-(PMAS(24,1)/PMAS(23,1))**2
      XW1=1D0-XW
      IF(MSTP(4).EQ.2) THEN
        TANBE=PARU(141)
        RATBE=((1D0-TANBE**2)/(1D0+TANBE**2))**2
        SQMZ=PMAS(23,1)**2
        SQMW=PMAS(24,1)**2
        SQMH=PMAS(25,1)**2
        SQMA=SQMH*(SQMZ-SQMH)/(SQMZ*RATBE-SQMH)
        SQMHP=0.5D0*(SQMA+SQMZ+SQRT((SQMA+SQMZ)**2-4D0*SQMA*SQMZ*RATBE))
        SQMHC=SQMA+SQMW
        IF(SQMH.GE.SQMZ.OR.MIN(SQMA,SQMHP,SQMHC).LE.0D0) THEN
          WRITE(MSTU(11),5000)
          CALL PYSTOP(101)
        ENDIF
        PMAS(35,1)=SQRT(SQMHP)
        PMAS(36,1)=SQRT(SQMA)
        PMAS(37,1)=SQRT(SQMHC)
        ALSU=0.5D0*ATAN(2D0*TANBE*(SQMA+SQMZ)/((1D0-TANBE**2)*
     &  (SQMA-SQMZ)))
        BESU=ATAN(TANBE)
        PARU(142)=1D0
        PARU(143)=1D0
        PARU(161)=-SIN(ALSU)/COS(BESU)
        PARU(162)=COS(ALSU)/SIN(BESU)
        PARU(163)=PARU(161)
        PARU(164)=SIN(BESU-ALSU)
        PARU(165)=PARU(164)
        PARU(168)=SIN(BESU-ALSU)+0.5D0*COS(2D0*BESU)*SIN(BESU+ALSU)/XW
        PARU(171)=COS(ALSU)/COS(BESU)
        PARU(172)=SIN(ALSU)/SIN(BESU)
        PARU(173)=PARU(171)
        PARU(174)=COS(BESU-ALSU)
        PARU(175)=PARU(174)
        PARU(176)=COS(2D0*ALSU)*COS(BESU+ALSU)-2D0*SIN(2D0*ALSU)*
     &  SIN(BESU+ALSU)
        PARU(177)=COS(2D0*BESU)*COS(BESU+ALSU)
        PARU(178)=COS(BESU-ALSU)-0.5D0*COS(2D0*BESU)*COS(BESU+ALSU)/XW
        PARU(181)=TANBE
        PARU(182)=1D0/TANBE
        PARU(183)=PARU(181)
        PARU(184)=0D0
        PARU(185)=PARU(184)
        PARU(186)=COS(BESU-ALSU)
        PARU(187)=SIN(BESU-ALSU)
        PARU(188)=PARU(186)
        PARU(189)=PARU(187)
        PARU(190)=0D0
        PARU(195)=COS(BESU-ALSU)
      ENDIF
 
C...Reset effective widths of gauge bosons.
      DO 110 I=1,500
        DO 100 J=1,5
          WIDS(I,J)=1D0
  100   CONTINUE
  110 CONTINUE
 
C...Order resonances by increasing mass (except Z0 and W+/-).
      NRES=0
      DO 140 KC=1,500
        KF=KCHG(KC,4)
        IF(KF.EQ.0) GOTO 140
        IF(MWID(KC).EQ.0) GOTO 140
        IF(KC.EQ.7.OR.KC.EQ.8.OR.KC.EQ.17.OR.KC.EQ.18) THEN
          IF(MSTP(1).LE.3) GOTO 140
        ENDIF
        IF(KF/KSUSY1.EQ.1.OR.KF/KSUSY1.EQ.2) THEN
          IF(IMSS(1).LE.0) GOTO 140
        ENDIF
        NRES=NRES+1
        PMRES=PMAS(KC,1)
        IF(KC.EQ.23.OR.KC.EQ.24) PMRES=0D0
        DO 120 I1=NRES-1,1,-1
          IF(PMRES.GE.PMORD(I1)) GOTO 130
          KCORD(I1+1)=KCORD(I1)
          PMORD(I1+1)=PMORD(I1)
  120   CONTINUE
  130   KCORD(I1+1)=KC
        PMORD(I1+1)=PMRES
  140 CONTINUE
 
C...Loop over possible resonances.
      DO 180 I=1,NRES
        KC=KCORD(I)
        KF=KCHG(KC,4)
 
C...Check that no fourth generation channels on by mistake.
        IF(MSTP(1).LE.3) THEN
          DO 150 J=1,MDCY(KC,3)
            IDC=J+MDCY(KC,2)-1
            KFA1=IABS(KFDP(IDC,1))
            KFA2=IABS(KFDP(IDC,2))
            IF(KFA1.EQ.7.OR.KFA1.EQ.8.OR.KFA1.EQ.17.OR.KFA1.EQ.18.OR.
     &      KFA2.EQ.7.OR.KFA2.EQ.8.OR.KFA2.EQ.17.OR.KFA2.EQ.18)
     &      MDME(IDC,1)=-1
  150     CONTINUE
        ENDIF
 
C...Check that no supersymmetric channels on by mistake.
        IF(IMSS(1).LE.0) THEN
          DO 160 J=1,MDCY(KC,3)
            IDC=J+MDCY(KC,2)-1
            KFA1S=IABS(KFDP(IDC,1))/KSUSY1
            KFA2S=IABS(KFDP(IDC,2))/KSUSY1
            IF(KFA1S.EQ.1.OR.KFA1S.EQ.2.OR.KFA2S.EQ.1.OR.KFA2S.EQ.2)
     &      MDME(IDC,1)=-1
  160     CONTINUE
        ENDIF
 
C...Find mass and evaluate width.
        PMR=PMAS(KC,1)
        IF(KF.EQ.25.OR.KF.EQ.35.OR.KF.EQ.36) MINT(62)=1
        IF(MWID(KC).EQ.3) MINT(63)=1
        CALL PYWIDT(KF,PMR**2,WDTP,WDTE)
        MINT(51)=0
 
C...Evaluate suppression factors due to non-simulated channels.
        IF(KCHG(KC,3).EQ.0) THEN
          WDTP0I=0D0
          IF(WDTP(0).GT.0D0) WDTP0I=1D0/WDTP(0)
          WIDS(KC,1)=((WDTE(0,1)+WDTE(0,2))**2+
     &    2D0*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+
     &    2D0*WDTE(0,4)*WDTE(0,5))*WDTP0I**2
          WIDS(KC,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))*WDTP0I
          WIDS(KC,3)=0D0
          WIDS(KC,4)=0D0
          WIDS(KC,5)=0D0
        ELSE
          IF(MWID(KC).EQ.3) MINT(63)=1
          CALL PYWIDT(-KF,PMR**2,WDTPM,WDTEM)
          MINT(51)=0
          WDTP0I=0D0
          IF(WDTP(0).GT.0D0) WDTP0I=1D0/WDTP(0)
          WIDS(KC,1)=((WDTE(0,1)+WDTE(0,2))*(WDTEM(0,1)+WDTEM(0,3))+
     &    (WDTE(0,1)+WDTE(0,2))*(WDTEM(0,4)+WDTEM(0,5))+
     &    (WDTE(0,4)+WDTE(0,5))*(WDTEM(0,1)+WDTEM(0,3))+
     &    WDTE(0,4)*WDTEM(0,5)+WDTE(0,5)*WDTEM(0,4))*WDTP0I**2
          WIDS(KC,2)=(WDTE(0,1)+WDTE(0,2)+WDTE(0,4))*WDTP0I
          WIDS(KC,3)=(WDTEM(0,1)+WDTEM(0,3)+WDTEM(0,4))*WDTP0I
          WIDS(KC,4)=((WDTE(0,1)+WDTE(0,2))**2+
     &    2D0*(WDTE(0,1)+WDTE(0,2))*(WDTE(0,4)+WDTE(0,5))+
     &    2D0*WDTE(0,4)*WDTE(0,5))*WDTP0I**2
          WIDS(KC,5)=((WDTEM(0,1)+WDTEM(0,3))**2+
     &    2D0*(WDTEM(0,1)+WDTEM(0,3))*(WDTEM(0,4)+WDTEM(0,5))+
     &    2D0*WDTEM(0,4)*WDTEM(0,5))*WDTP0I**2
        ENDIF
 
C...Set resonance widths and branching ratios;
C...also on/off switch for decays.
        IF(MWID(KC).EQ.1.OR.MWID(KC).EQ.3) THEN
          PMAS(KC,2)=WDTP(0)
          PMAS(KC,3)=MIN(0.9D0*PMAS(KC,1),10D0*PMAS(KC,2))
          IF(MSTP(41).EQ.0.OR.MSTP(41).EQ.1) MDCY(KC,1)=MSTP(41)
          DO 170 J=1,MDCY(KC,3)
            IDC=J+MDCY(KC,2)-1
            BRAT(IDC)=0D0
            IF(WDTP(0).GT.0D0) BRAT(IDC)=WDTP(J)/WDTP(0)
  170     CONTINUE
        ENDIF
  180 CONTINUE
 
C...Flavours of leptoquark: redefine charge and name.
      KFLQQ=KFDP(MDCY(42,2),1)
      KFLQL=KFDP(MDCY(42,2),2)
      KCHG(42,1)=KCHG(PYCOMP(KFLQQ),1)*ISIGN(1,KFLQQ)+
     &KCHG(PYCOMP(KFLQL),1)*ISIGN(1,KFLQL)
      LL=1
      IF(IABS(KFLQL).EQ.13) LL=2
      IF(IABS(KFLQL).EQ.15) LL=3
      CHAF(42,1)='LQ_'//CHAF(IABS(KFLQQ),1)(1:1)//
     &CHAF(IABS(KFLQL),1)(1:LL)//' '
      CHAF(42,2)=CHAF(42,2)(1:4+LL)//'bar '
 
C...Special cases in treatment of gamma*/Z0: redefine process name.
      IF(MSTP(43).EQ.1) THEN
        PROC(1)='f + fbar -> gamma*'
        PROC(15)='f + fbar -> g + gamma*'
        PROC(19)='f + fbar -> gamma + gamma*'
        PROC(30)='f + g -> f + gamma*'
        PROC(35)='f + gamma -> f + gamma*'
      ELSEIF(MSTP(43).EQ.2) THEN
        PROC(1)='f + fbar -> Z0'
        PROC(15)='f + fbar -> g + Z0'
        PROC(19)='f + fbar -> gamma + Z0'
        PROC(30)='f + g -> f + Z0'
        PROC(35)='f + gamma -> f + Z0'
      ELSEIF(MSTP(43).EQ.3) THEN
        PROC(1)='f + fbar -> gamma*/Z0'
        PROC(15)='f + fbar -> g + gamma*/Z0'
        PROC(19)='f+ fbar -> gamma + gamma*/Z0'
        PROC(30)='f + g -> f + gamma*/Z0'
        PROC(35)='f + gamma -> f + gamma*/Z0'
      ENDIF
 
C...Special cases in treatment of gamma*/Z0/Z'0: redefine process name.
      IF(MSTP(44).EQ.1) THEN
        PROC(141)='f + fbar -> gamma*'
      ELSEIF(MSTP(44).EQ.2) THEN
        PROC(141)='f + fbar -> Z0'
      ELSEIF(MSTP(44).EQ.3) THEN
        PROC(141)='f + fbar -> Z''0'
      ELSEIF(MSTP(44).EQ.4) THEN
        PROC(141)='f + fbar -> gamma*/Z0'
      ELSEIF(MSTP(44).EQ.5) THEN
        PROC(141)='f + fbar -> gamma*/Z''0'
      ELSEIF(MSTP(44).EQ.6) THEN
        PROC(141)='f + fbar -> Z0/Z''0'
      ELSEIF(MSTP(44).EQ.7) THEN
        PROC(141)='f + fbar -> gamma*/Z0/Z''0'
      ENDIF
 
C...Special cases in treatment of WW -> WW: redefine process name.
      IF(MSTP(45).EQ.1) THEN
        PROC(77)='W+ + W+ -> W+ + W+'
      ELSEIF(MSTP(45).EQ.2) THEN
        PROC(77)='W+ + W- -> W+ + W-'
      ELSEIF(MSTP(45).EQ.3) THEN
        PROC(77)='W+/- + W+/- -> W+/- + W+/-'
      ENDIF

C...Initialize Generic Processes
      KFGEN=9900001
      KCGEN=PYCOMP(KFGEN)
      IF(KCGEN.GT.0) THEN
        IDCY=MDCY(KCGEN,2)
        IF(IDCY.GT.0) THEN
          KFF1=KFDP(IDCY+1,1)
          KFF2=KFDP(IDCY+1,2)
          KCF1=PYCOMP(KFF1)
          KCF2=PYCOMP(KFF2)
          IJ1=1
          IJ2=1
          KCI1=PYCOMP(KFDP(IDCY,1))
          IF(KFDP(IDCY,1).LT.0) IJ1=2
          KCI2=PYCOMP(KFDP(IDCY,2))
          IF(KFDP(IDCY,2).LT.0) IJ2=2
          ITMP1=0
 190      ITMP1=ITMP1+1
          IF(CHAF(KCI1,IJ1)(ITMP1+1:ITMP1+1).NE.' '.AND.ITMP1.LT.4)
     &    GOTO 190
          ITMP2=0
 200      ITMP2=ITMP2+1
          IF(CHAF(KCI2,IJ2)(ITMP2+1:ITMP2+1).NE.' '.AND.ITMP2.LT.4)
     &    GOTO 200          
          PRTMP=CHAF(KCI1,IJ1)(1:ITMP1)//'+'//CHAF(KCI2,IJ2)(1:ITMP2)
          ITMP3=0
 205      ITMP3=ITMP3+1
          IF(PRTMP(ITMP3+1:ITMP3+1).NE.' '.AND.ITMP3.LT.9)
     &    GOTO 205
          PROC(481)=PRTMP(1:ITMP3)//' -> '//CHAF(KCGEN,1)
          IJ1=1
          IJ2=1
          IF(KFF1.LT.0) IJ1=2
          IF(KFF2.LT.0) IJ2=2
          ITMP1=0
 210      ITMP1=ITMP1+1
          IF(CHAF(KCF1,IJ1)(ITMP1+1:ITMP1+1).NE.' '.AND.ITMP1.LT.8)
     &    GOTO 210
          ITMP2=0
 220      ITMP2=ITMP2+1
          IF(CHAF(KCF2,IJ2)(ITMP2+1:ITMP2+1).NE.' '.AND.ITMP2.LT.8)
     &    GOTO 220          
          PROC(482)=PRTMP(1:ITMP3)//' -> '//CHAF(KCF1,IJ1)(1:ITMP1)//
     &    '+'//CHAF(KCF2,IJ2)(1:ITMP2)
        ENDIF
      ENDIF


 
C...Format for error information.
 5000 FORMAT(1X,'Error: unphysical input tan^2(beta) and m_H ',
     &'combination'/1X,'Execution stopped!')
 
      RETURN
      END
