
C*********************************************************************
 
C...PYSSPA
C...Generates spacelike parton showers.
 
      SUBROUTINE PYSSPA(IPU1,IPU2)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (MAXNUR=1000)
C...Commonblocks.
      COMMON/PYPART/NPART,NPARTD,IPART(MAXNUR),PTPART(MAXNUR)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYCTAG/NCT,MCT(4000,2)
      SAVE /PYPART/,/PYJETS/,/PYDAT1/,/PYDAT2/,/PYSUBS/,/PYPARS/,
     &/PYINT1/,/PYINT2/,/PYINT3/,/PYCTAG/
C...Local arrays and data.
      DIMENSION KFLS(4),IS(2),XS(2),ZS(2),Q2S(2),TEVCSV(2),TEVESV(2),
     &XFS(2,-25:25),XFA(-25:25),XFB(-25:25),XFN(-25:25),WTAPC(-25:25),
     &WTAPE(-25:25),WTSF(-25:25),THE2(2),ALAM(2),DQ2(3),DPC(3),DPD(4),
     &DPB(4),ROBO(5),MORE(2),KFBEAM(2),Q2MNCS(2),KCFI(2),NFIS(2),
     &THEFIS(2,2),ISFI(2),DPHI(2),MCESV(2)
      DATA IS/2*0/
 
C...Read out basic information; set global Q^2 scale.
      IPUS1=IPU1
      IPUS2=IPU2
      ISUB=MINT(1)
      Q2MX=VINT(56)
      VINT2R=VINT(2)*VINT(143)*VINT(144)
      IF(ISET(ISUB).EQ.2.OR.ISET(ISUB).EQ.9.OR.ISET(ISUB).EQ.11) Q2MX=
     &MIN(VINT2R,PARP(67)*VINT(56))
      FCQ2MX=1D0
 
C...Define which processes ME corrections have been implemented for.
      MECOR=0
      IF(MSTP(68).EQ.1.OR.MSTP(68).EQ.3) THEN
        IF(ISUB.EQ.1.OR.ISUB.EQ.2.OR.ISUB.EQ.141.OR.ISUB.EQ.142.OR.
     &  ISUB.EQ.144) MECOR=1
        IF(ISUB.EQ.102.OR.ISUB.EQ.152.OR.ISUB.EQ.157) MECOR=2
        IF(ISUB.EQ.3.OR.ISUB.EQ.151.OR.ISUB.EQ.156) MECOR=3
      ENDIF
 
C...Initialize QCD evolution and check phase space.
      Q2MNC=PARP(62)**2
      Q2MNCS(1)=Q2MNC
      Q2MNCS(2)=Q2MNC
      IF(MINT(107).EQ.2.AND.MSTP(66).EQ.2) THEN
        Q0S=PARP(15)**2
        PS=VINT(3)**2
        Q2EFF=VINT(54)*((Q0S+PS)/(VINT(54)+PS))*
     &  EXP(PS*(VINT(54)-Q0S)/((VINT(54)+PS)*(Q0S+PS)))
        Q2INT=SQRT(Q0S*Q2EFF)
        Q2MNCS(1)=MAX(Q2MNC,Q2INT)
      ELSEIF(MINT(107).EQ.3.AND.MSTP(66).GE.1) THEN
        Q2MNCS(1)=MAX(Q2MNC,VINT(283))
      ENDIF
      IF(MINT(108).EQ.2.AND.MSTP(66).EQ.2) THEN
        Q0S=PARP(15)**2
        PS=VINT(4)**2
        Q2EFF=VINT(54)*((Q0S+PS)/(VINT(54)+PS))*
     &  EXP(PS*(VINT(54)-Q0S)/((VINT(54)+PS)*(Q0S+PS)))
        Q2INT=SQRT(Q0S*Q2EFF)
        Q2MNCS(2)=MAX(Q2MNC,Q2INT)
      ELSEIF(MINT(108).EQ.3.AND.MSTP(66).GE.1) THEN
        Q2MNCS(2)=MAX(Q2MNC,VINT(284))
      ENDIF
      MCEV=0
      ALAMS=PARU(112)
      PARU(112)=PARP(61)
      FQ2C=1D0
      TCMX=0D0
      IF(MINT(47).GE.2.AND.(MINT(47).LT.5.OR.MSTP(12).GE.1)) THEN
        MCEV=1
        IF(MSTP(64).EQ.1) FQ2C=PARP(63)
        IF(MSTP(64).EQ.2) FQ2C=PARP(64)
        TCMX=LOG(FQ2C*Q2MX/PARP(61)**2)
        IF(Q2MX.LT.MAX(Q2MNC,2D0*PARP(61)**2).OR.TCMX.LT.0.2D0)
     &  MCEV=0
      ENDIF
 
C...Initialize QED evolution and check phase space.
      MEEV=0
      XEE=1D-10
      SPME=PMAS(11,1)**2
      IF(IABS(MINT(11)).EQ.13.OR.IABS(MINT(12)).EQ.13)
     &SPME=PMAS(13,1)**2
      IF(IABS(MINT(11)).EQ.15.OR.IABS(MINT(12)).EQ.15)
     &SPME=PMAS(15,1)**2
      Q2MNE=MAX(PARP(68)**2,2D0*SPME)
      TEMX=0D0
      FWTE=10D0
      IF(MINT(45).EQ.3.OR.MINT(46).EQ.3) THEN
        MEEV=1
        TEMX=LOG(Q2MX/SPME)
        IF(Q2MX.LE.Q2MNE.OR.TEMX.LT.0.2D0) MEEV=0
      ENDIF
      IF(MSTP(61).GE.2.AND.MCEV.EQ.1.AND.MEEV.EQ.0) THEN
        MEEV=2
        TEMX=TCMX
        FWTE=1D0
      ENDIF
      IF(MCEV.EQ.0.AND.MEEV.EQ.0) RETURN
 
C...Loopback point in case of failure to reconstruct kinematics.
      NS=N
      NPARTS=NPART
      LOOP=0      
      MNT352=MINT(352)
      MNT353=MINT(353)
      VNT352=VINT(352)
      VNT353=VINT(353)
  100 LOOP=LOOP+1
      IF(LOOP.GT.100) THEN
        MINT(51)=1
        RETURN
      ENDIF
      N=NS
      NPART=NPARTS
      MINT(352)=MNT352
      MINT(353)=MNT353
      VINT(352)=VNT352
      VINT(353)=VNT353
 
C...Initial values: flavours, momenta, virtualities.
      DO 120 JT=1,2
        MORE(JT)=1
        KFBEAM(JT)=MINT(10+JT)
        IF(MINT(18+JT).EQ.1)KFBEAM(JT)=22
        KFLS(JT)=MINT(14+JT)
        KFLS(JT+2)=KFLS(JT)
        XS(JT)=VINT(40+JT)
        IF(MINT(18+JT).EQ.1) XS(JT)=VINT(40+JT)/VINT(154+JT)
        IF(MINT(31).GE.2) XS(JT)=XS(JT)/VINT(142+JT)
        ZS(JT)=1D0
        Q2S(JT)=FCQ2MX*Q2MX
        DQ2(JT)=0D0
        TEVCSV(JT)=TCMX
        ALAM(JT)=PARP(61)
        THE2(JT)=1D0
        TEVESV(JT)=TEMX
        MCESV(JT)=0
C...Calculate initial parton distribution weights.
        MINT(105)=MINT(102+JT)
        MINT(109)=MINT(106+JT)
        VINT(120)=VINT(2+JT)
        IF(XS(JT).LT.1D0-XEE) THEN
          IF(MINT(31).GE.2) MINT(30)=JT
          IF(MSTP(57).LE.1) THEN
            CALL PYPDFU(KFBEAM(JT),XS(JT),Q2S(JT),XFB)
          ELSE
            CALL PYPDFL(KFBEAM(JT),XS(JT),Q2S(JT),XFB)
          ENDIF
        ENDIF
        DO 110 KFL=-25,25
          XFS(JT,KFL)=XFB(KFL)
  110   CONTINUE
C...Special kinematics check for c/b quarks (that g -> c cbar or
C...b bbar kinematically possible).
      KFLCB=IABS(KFLS(JT))
      IF(KFBEAM(JT).NE.22.AND.(KFLCB.EQ.4.OR.KFLCB.EQ.5)) THEN
        IF(XS(JT).GT.0.9D0*Q2S(JT)/(PMAS(KFLCB,1)**2+Q2S(JT))) THEN
          MINT(51)=1
          RETURN
        ENDIF
      ENDIF
  120 CONTINUE
      DSH=VINT(44)
      IF(ISET(ISUB).GE.3.AND.ISET(ISUB).LE.5) DSH=VINT(26)*VINT(2)
 
C...Find if interference with final state partons.
      MFIS=0
      IF(MSTP(67).GE.1.AND.MSTP(67).LE.3) MFIS=MSTP(67)
      IF(MFIS.NE.0) THEN
        DO 140 I=1,2
          KCFI(I)=0
          KCA=PYCOMP(IABS(KFLS(I)))
          IF(KCA.NE.0) KCFI(I)=KCHG(KCA,2)*ISIGN(1,KFLS(I))
          NFIS(I)=0
          IF(KCFI(I).NE.0) THEN
            IF(I.EQ.1) IPFS=IPUS1
            IF(I.EQ.2) IPFS=IPUS2
            DO 130 J=1,2
              ICSI=MOD(K(IPFS,3+J),MSTU(5))
              IF(ICSI.GT.0.AND.ICSI.NE.IPUS1.AND.ICSI.NE.IPUS2.AND.
     &        (KCFI(I).EQ.(-1)**(J+1).OR.KCFI(I).EQ.2)) THEN
                NFIS(I)=NFIS(I)+1
                THEFIS(I,NFIS(I))=PYANGL(P(ICSI,3),SQRT(P(ICSI,1)**2+
     &          P(ICSI,2)**2))
                IF(I.EQ.2) THEFIS(I,NFIS(I))=PARU(1)-THEFIS(I,NFIS(I))
              ENDIF
  130       CONTINUE
          ENDIF
  140   CONTINUE
        IF(NFIS(1)+NFIS(2).EQ.0) MFIS=0
      ENDIF
 
C...Pick up leg with highest virtuality.
      JTOLD=1
  150 N=N+1
      JT=1
      IF(N.GT.NS+1.AND.Q2S(2).GT.Q2S(1)) JT=2
      IF(N.EQ.NS+2.AND.JT.EQ.JTOLD) JT=3-JT
      IF(MORE(JT).EQ.0) JT=3-JT
      JTOLD=JT
      KFLB=KFLS(JT)
      XB=XS(JT)
      DO 160 KFL=-25,25
        XFB(KFL)=XFS(JT,KFL)
  160 CONTINUE
      DSHR=2D0*SQRT(DSH)
      DSHZ=DSH/ZS(JT)
 
C...Check if allowed to branch.
      MCEV=0
      IF(IABS(KFLB).LE.10.OR.KFLB.EQ.21) THEN
        MCEV=1
        XEC=MAX(PARP(65)*DSHR/VINT2R,XB*(1D0/(1D0-PARP(66))-1D0))
        IF(XB.GE.1D0-2D0*XEC) MCEV=0
      ENDIF
      MEEV=0
      IF(MINT(44+JT).EQ.3) THEN
        MEEV=1
        IF(XB.GE.1D0-2D0*XEE) MEEV=0
        IF((IABS(KFLB).LE.10.OR.KFLB.EQ.21).AND.XB.GE.1D0-2D0*XEC)
     &  MEEV=0
C***Currently kill QED shower for resolved photoproduction.
        IF(MINT(18+JT).EQ.1) MEEV=0
C***Currently kill shower for W inside electron.
        IF(IABS(KFLB).EQ.24) THEN
          MCEV=0
          MEEV=0
        ENDIF
      ENDIF
      IF(MSTP(61).GE.2.AND.MCEV.EQ.1.AND.MEEV.EQ.0.AND.IABS(KFLB).LE.10)
     &MEEV=2
      IF(MCEV.EQ.0.AND.MEEV.EQ.0) THEN
        Q2B=0D0
        GOTO 260
      ENDIF
 
C...Maximum Q2 with or without Q2 ordering. Effective Lambda and n_f.
      Q2B=Q2S(JT)
      TEVCB=TEVCSV(JT)
      TEVEB=TEVESV(JT)
      IF(MSTP(62).LE.1) THEN
        IF(ZS(JT).GT.0.99999D0) THEN
          Q2B=Q2S(JT)
        ELSE
          Q2B=0.5D0*(1D0/ZS(JT)+1D0)*Q2S(JT)+0.5D0*(1D0/ZS(JT)-1D0)*
     &    (Q2S(3-JT)-DSH+SQRT((DSH+Q2S(1)+Q2S(2))**2+
     &    8D0*Q2S(1)*Q2S(2)*ZS(JT)/(1D0-ZS(JT))))
        ENDIF
        IF(MCEV.EQ.1) TEVCB=LOG(FQ2C*Q2B/ALAM(JT)**2)
        IF(MEEV.EQ.1) TEVEB=LOG(Q2B/SPME)
      ENDIF
      IF(MCEV.EQ.1) THEN
        ALSDUM=PYALPS(FQ2C*Q2B)
        TEVCB=TEVCB+2D0*LOG(ALAM(JT)/PARU(117))
        ALAM(JT)=PARU(117)
        B0=(33D0-2D0*MSTU(118))/6D0
      ENDIF
      IF(MEEV.EQ.2) TEVEB=TEVCB
      TEVCBS=TEVCB
      TEVEBS=TEVEB
 
C...Select side for interference with final state partons.
      IF(MFIS.GE.1.AND.N.LE.NS+2) THEN
        IFI=N-NS
        ISFI(IFI)=0
        IF(IABS(KCFI(IFI)).EQ.1.AND.NFIS(IFI).EQ.1) THEN
          ISFI(IFI)=1
        ELSEIF(KCFI(IFI).EQ.2.AND.NFIS(IFI).EQ.1) THEN
          IF(PYR(0).GT.0.5D0) ISFI(IFI)=1
        ELSEIF(KCFI(IFI).EQ.2.AND.NFIS(IFI).EQ.2) THEN
          ISFI(IFI)=1
          IF(PYR(0).GT.0.5D0) ISFI(IFI)=2
        ENDIF
      ENDIF
 
C...Calculate preweighting factor for ME-corrected processes.
      IF(MECOR.GE.1) CALL PYMEMX(MECOR,WTFF,WTGF,WTFG,WTGG)
 
C...Calculate Altarelli-Parisi weights.
      DO 170 KFL=-25,25
        WTAPC(KFL)=0D0
        WTAPE(KFL)=0D0
        WTSF(KFL)=0D0
  170 CONTINUE
C...q -> q (g or gamma emission), g -> q.
      IF(IABS(KFLB).LE.10) THEN
        WTAPC(KFLB)=(8D0/3D0)*LOG((1D0-XEC-XB)*(XB+XEC)/(XEC*(1D0-XEC)))
        WTAPC(21)=0.5D0*(XB/(XB+XEC)-XB/(1D0-XEC))
        EQ2=1D0/9D0
        IF(MOD(IABS(KFLB),2).EQ.0) EQ2=4D0*EQ2
        IF(MEEV.EQ.2) WTAPE(KFLB)=2.*EQ2*LOG((1D0-XEC-XB)*(XB+XEC)/
     &  (XEC*(1D0-XEC)))
        IF(MECOR.GE.1.AND.(N.EQ.NS+1.OR.N.EQ.NS+2)) THEN
          WTAPC(KFLB)=WTFF*WTAPC(KFLB)
          WTAPC(21)=WTGF*WTAPC(21)
          WTAPE(KFLB)=WTFF*WTAPE(KFLB)
        ENDIF
C...f -> f, gamma -> f.
      ELSEIF(IABS(KFLB).LE.20) THEN
        WTAPF1=LOG((1D0-XEE-XB)*(XB+XEE)/(XEE*(1D0-XEE)))
        WTAPF2=LOG((1D0-XEE-XB)*(1D0-XEE)/(XEE*(XB+XEE)))
        WTAPE(KFLB)=2D0*(WTAPF1+WTAPF2)
        IF(MSTP(12).GE.1) WTAPE(22)=XB/(XB+XEE)-XB/(1D0-XEE)
        IF(MECOR.GE.1.AND.(N.EQ.NS+1.OR.N.EQ.NS+2)) THEN
          WTAPE(KFLB)=WTFF*WTAPE(KFLB)
          WTAPE(22)=WTGF*WTAPE(22)
        ENDIF
C...f -> g, g -> g.
      ELSEIF(KFLB.EQ.21) THEN
        WTAPQ=(16D0/3D0)*(SQRT((1D0-XEC)/XB)-SQRT((XB+XEC)/XB))
        DO 180 KFL=1,MSTP(58)
          WTAPC(KFL)=WTAPQ
          WTAPC(-KFL)=WTAPQ
  180   CONTINUE
        WTAPC(21)=6D0*LOG((1D0-XEC-XB)/XEC)
        IF(MECOR.GE.1.AND.(N.EQ.NS+1.OR.N.EQ.NS+2)) THEN
          DO 190 KFL=1,MSTP(58)
            WTAPC(KFL)=WTFG*WTAPC(KFL)
            WTAPC(-KFL)=WTFG*WTAPC(-KFL)
  190     CONTINUE
          WTAPC(21)=WTGG*WTAPC(21)
        ENDIF
C...f -> gamma, W+, W-.
      ELSEIF(KFLB.EQ.22) THEN
        WTAPF=LOG((1D0-XEE-XB)*(1D0-XEE)/(XEE*(XB+XEE)))/XB
        WTAPE(11)=WTAPF
        WTAPE(-11)=WTAPF
        IF(MECOR.GE.1.AND.(N.EQ.NS+1.OR.N.EQ.NS+2)) THEN
          WTAPE(11)=WTFG*WTAPE(11)
          WTAPE(-11)=WTFG*WTAPE(-11)
        ENDIF
      ELSEIF(KFLB.EQ.24) THEN
        WTAPE(-11)=1D0/(4D0*PARU(102))*LOG((1D0-XEE-XB)*(1D0-XEE)/
     &  (XEE*(XB+XEE)))/XB
      ELSEIF(KFLB.EQ.-24) THEN
        WTAPE(11)=1D0/(4D0*PARU(102))*LOG((1D0-XEE-XB)*(1D0-XEE)/
     &  (XEE*(XB+XEE)))/XB
      ENDIF
 
C...Calculate parton distribution weights and sum.
      NTRY=0
  200 NTRY=NTRY+1
      IF(NTRY.GT.500) THEN
        MINT(51)=1
        RETURN
      ENDIF
      WTSUMC=0D0
      WTSUME=0D0
      XFBO=MAX(1D-10,XFB(KFLB))
      DO 210 KFL=-25,25
        WTSF(KFL)=XFB(KFL)/XFBO
        WTSUMC=WTSUMC+WTAPC(KFL)*WTSF(KFL)
        WTSUME=WTSUME+WTAPE(KFL)*WTSF(KFL)
  210 CONTINUE
      WTSUMC=MAX(0.0001D0,WTSUMC)
      WTSUME=MAX(0.0001D0/FWTE,WTSUME)
 
C...Choose new t: fix alpha_s, alpha_s(Q^2), alpha_s(k_T^2).
      NTRY2=0
  220 NTRY2=NTRY2+1
      IF(NTRY2.GT.500) THEN
        MINT(51)=1
        RETURN
      ENDIF
      IF(MCEV.EQ.1) THEN
        IF(MSTP(64).LE.0) THEN
          TEVCB=TEVCB+LOG(PYR(0))*PARU(2)/(PARU(111)*WTSUMC)
        ELSEIF(MSTP(64).EQ.1) THEN
          TEVCB=TEVCB*EXP(MAX(-50D0,LOG(PYR(0))*B0/WTSUMC))
        ELSE
          TEVCB=TEVCB*EXP(MAX(-50D0,LOG(PYR(0))*B0/(5D0*WTSUMC)))
        ENDIF
      ENDIF
      IF(MEEV.EQ.1) THEN
        TEVEB=TEVEB*EXP(MAX(-50D0,LOG(PYR(0))*PARU(2)/
     &  (PARU(101)*FWTE*WTSUME*TEMX)))
      ELSEIF(MEEV.EQ.2) THEN
        TEVEB=TEVEB+LOG(PYR(0))*PARU(2)/(PARU(101)*WTSUME)
      ENDIF
 
C...Translate t into Q2 scale; choose between QCD and QED evolution.
  230 IF(MCEV.EQ.1) Q2CB=ALAM(JT)**2*EXP(MAX(-50D0,TEVCB))/FQ2C
      IF(MEEV.EQ.1) Q2EB=SPME*EXP(MAX(-50D0,TEVEB))
      IF(MEEV.EQ.2) Q2EB=ALAM(JT)**2*EXP(MAX(-50D0,TEVEB))/FQ2C
C...Ensure that Q2 is above threshold for charm/bottom.
      KFLCB=IABS(KFLB)
      IF(KFBEAM(JT).NE.22.AND.(KFLCB.EQ.4.OR.KFLCB.EQ.5).AND.
     &MCEV.EQ.1) THEN
        IF(Q2CB.LT.PMAS(KFLCB,1)**2) THEN
          Q2CB=1.1D0*PMAS(KFLCB,1)**2
          TEVCB=LOG(FQ2C*Q2B/ALAM(JT)**2)
          FCQ2MX=MIN(2D0,1.05D0*FCQ2MX)
        ENDIF
      ENDIF
      IF(KFBEAM(JT).NE.22.AND.(KFLCB.EQ.4.OR.KFLCB.EQ.5).AND.
     &MEEV.EQ.2) THEN
        IF(Q2EB.LT.PMAS(KFLCB,1)**2) MEEV=0
      ENDIF
      MCE=0
      IF(MCEV.EQ.0.AND.MEEV.EQ.0) THEN
      ELSEIF(MCEV.EQ.1.AND.MEEV.EQ.0) THEN
        IF(Q2CB.GT.Q2MNCS(JT)) MCE=1
      ELSEIF(MCEV.EQ.0.AND.MEEV.EQ.1) THEN
        IF(Q2EB.GT.Q2MNE) MCE=2
      ELSEIF(MCEV.EQ.0.AND.MEEV.EQ.2) THEN
        IF(Q2EB.GT.Q2MNCS(JT)) MCE=2
      ELSEIF(MCEV.EQ.1.AND.MEEV.EQ.2) THEN
        IF(Q2CB.GT.Q2EB.AND.Q2CB.GT.Q2MNCS(JT)) MCE=1
        IF(Q2EB.GT.Q2CB.AND.Q2EB.GT.Q2MNCS(JT)) MCE=2
      ELSEIF(Q2MNCS(JT).GT.Q2MNE) THEN
        MCE=1
        IF(Q2EB.GT.Q2CB.OR.Q2CB.LE.Q2MNCS(JT)) MCE=2
        IF(MCE.EQ.2.AND.Q2EB.LE.Q2MNE) MCE=0
      ELSE
        MCE=2
        IF(Q2CB.GT.Q2EB.OR.Q2EB.LE.Q2MNE) MCE=1
        IF(MCE.EQ.1.AND.Q2CB.LE.Q2MNCS(JT)) MCE=0
      ENDIF
 
C...Evolution possibly ended. Update t values.
      IF(MCE.EQ.0) THEN
        Q2B=0D0
        GOTO 260
      ELSEIF(MCE.EQ.1) THEN
        Q2B=Q2CB
        Q2REF=FQ2C*Q2B
        IF(MEEV.EQ.1) TEVEB=LOG(Q2B/SPME)
        IF(MEEV.EQ.2) TEVEB=LOG(FQ2C*Q2B/ALAM(JT)**2)
      ELSE
        Q2B=Q2EB
        Q2REF=Q2B
        IF(MCEV.EQ.1) TEVCB=LOG(FQ2C*Q2B/ALAM(JT)**2)
      ENDIF
 
C...Select flavour for branching parton.
      IF(MCE.EQ.1) WTRAN=PYR(0)*WTSUMC
      IF(MCE.EQ.2) WTRAN=PYR(0)*WTSUME
      KFLA=-25
  240 KFLA=KFLA+1
      IF(MCE.EQ.1) WTRAN=WTRAN-WTAPC(KFLA)*WTSF(KFLA)
      IF(MCE.EQ.2) WTRAN=WTRAN-WTAPE(KFLA)*WTSF(KFLA)
      IF(KFLA.LE.24.AND.WTRAN.GT.0D0) GOTO 240
      IF(KFLA.EQ.25) THEN
        Q2B=0D0
        GOTO 260
      ENDIF
 
C...Choose z value and corrective weight.
      WTZ=0D0
C...q -> q + g or q -> q + gamma.
      IF(IABS(KFLA).LE.10.AND.IABS(KFLB).LE.10) THEN
        Z=1D0-((1D0-XB-XEC)/(1D0-XEC))*
     &  (XEC*(1D0-XEC)/((XB+XEC)*(1D0-XB-XEC)))**PYR(0)
        WTZ=0.5D0*(1D0+Z**2)
C...q -> g + q.
      ELSEIF(IABS(KFLA).LE.10.AND.KFLB.EQ.21) THEN
        Z=XB/(SQRT(XB+XEC)+PYR(0)*(SQRT(1D0-XEC)-SQRT(XB+XEC)))**2
        WTZ=0.5D0*(1D0+(1D0-Z)**2)*SQRT(Z)
C...f -> f + gamma.
      ELSEIF(IABS(KFLA).LE.20.AND.IABS(KFLB).LE.20) THEN
        IF(WTAPF1.GT.PYR(0)*(WTAPF1+WTAPF2)) THEN
          Z=1D0-((1D0-XB-XEE)/(1D0-XEE))*
     &    (XEE*(1D0-XEE)/((XB+XEE)*(1D0-XB-XEE)))**PYR(0)
        ELSE
          Z=XB+XB*(XEE/(1D0-XEE))*
     &    ((1D0-XB-XEE)*(1D0-XEE)/(XEE*(XB+XEE)))**PYR(0)
        ENDIF
        WTZ=0.5D0*(1D0+Z**2)*(Z-XB)/(1D0-XB)
C...f -> gamma + f.
      ELSEIF(IABS(KFLA).LE.20.AND.KFLB.EQ.22) THEN
        Z=XB+XB*(XEE/(1D0-XEE))*
     &  ((1D0-XB-XEE)*(1D0-XEE)/(XEE*(XB+XEE)))**PYR(0)
        WTZ=0.5D0*(1D0+(1D0-Z)**2)*XB*(Z-XB)/Z
C...f -> W+- + f.
      ELSEIF(IABS(KFLA).LE.20.AND.IABS(KFLB).EQ.24) THEN
        Z=XB+XB*(XEE/(1D0-XEE))*
     &  ((1D0-XB-XEE)*(1D0-XEE)/(XEE*(XB+XEE)))**PYR(0)
        WTZ=0.5D0*(1D0+(1D0-Z)**2)*(XB*(Z-XB)/Z)*
     &  (Q2B/(Q2B+PMAS(24,1)**2))
C...g -> q + qbar.
      ELSEIF(KFLA.EQ.21.AND.IABS(KFLB).LE.10) THEN
        Z=XB/(1D0-XEC)+PYR(0)*(XB/(XB+XEC)-XB/(1D0-XEC))
        WTZ=1D0-2D0*Z*(1D0-Z)
C...g -> g + g.
      ELSEIF(KFLA.EQ.21.AND.KFLB.EQ.21) THEN
        Z=1D0/(1D0+((1D0-XEC-XB)/XB)*(XEC/(1D0-XEC-XB))**PYR(0))
        WTZ=(1D0-Z*(1D0-Z))**2
C...gamma -> f + fbar.
      ELSEIF(KFLA.EQ.22.AND.IABS(KFLB).LE.20) THEN
        Z=XB/(1D0-XEE)+PYR(0)*(XB/(XB+XEE)-XB/(1D0-XEE))
        WTZ=1D0-2D0*Z*(1D0-Z)
      ENDIF
      IF(MCE.EQ.2.AND.MEEV.EQ.1) WTZ=(WTZ/FWTE)*(TEVEB/TEMX)
 
C...Option with resummation of soft gluon emission as effective z shift.
      IF(MCE.EQ.1) THEN
        IF(MSTP(65).GE.1) THEN
          RSOFT=6D0
          IF(KFLB.NE.21) RSOFT=8D0/3D0
          Z=Z*(TEVCB/TEVCSV(JT))**(RSOFT*XEC/((XB+XEC)*B0))
          IF(Z.LE.XB) GOTO 220
        ENDIF
 
C...Option with alpha_s(k_T^2): demand k_T^2 > cutoff, reweight.
        IF(MSTP(64).GE.2) THEN
          IF((1D0-Z)*Q2B.LT.Q2MNCS(JT)) GOTO 220
          ALPRAT=TEVCB/(TEVCB+LOG(1D0-Z))
          IF(ALPRAT.LT.5D0*PYR(0)) GOTO 220
          IF(ALPRAT.GT.5D0) WTZ=WTZ*ALPRAT/5D0
        ENDIF
      ENDIF
 
C...Remove kinematically impossible branchings.
      UHAT=Q2B-DSH*(1D0-Z)/Z
      IF(MSTP(68).GE.0.AND.UHAT.GT.0D0) GOTO 220
 
C...Select phi angle of branching at random.
      PHIBR=PARU(2)*PYR(0)
 
C...Matrix-element corrections for some processes.
      IF(MECOR.GE.1.AND.(N.EQ.NS+1.OR.N.EQ.NS+2)) THEN
        IF(IABS(KFLA).LE.20.AND.IABS(KFLB).LE.20) THEN
          CALL PYMEWT(MECOR,1,Q2B,Z,PHIBR,WTME)
          WTZ=WTZ*WTME/WTFF
        ELSEIF((KFLA.EQ.21.OR.KFLA.EQ.22).AND.IABS(KFLB).LE.20) THEN
          CALL PYMEWT(MECOR,2,Q2B,Z,PHIBR,WTME)
          WTZ=WTZ*WTME/WTGF
        ELSEIF(IABS(KFLA).LE.20.AND.(KFLB.EQ.21.OR.KFLB.EQ.22)) THEN
          CALL PYMEWT(MECOR,3,Q2B,Z,PHIBR,WTME)
          WTZ=WTZ*WTME/WTFG
        ELSEIF(KFLA.EQ.21.AND.KFLB.EQ.21) THEN
          CALL PYMEWT(MECOR,4,Q2B,Z,PHIBR,WTME)
          WTZ=WTZ*WTME/WTGG
        ENDIF
      ENDIF
 
C...Impose angular constraint in first branching from interference
C...with final state partons.
      IF(MCE.EQ.1) THEN
        IF(MFIS.GE.1.AND.N.LE.NS+2.AND.NTRY2.LT.200) THEN
          THE2D=(4D0*Q2B)/(DSH*(1D0-Z))
          IF(N.EQ.NS+1.AND.ISFI(1).GE.1) THEN
            IF(THE2D.GT.THEFIS(1,ISFI(1))**2) GOTO 220
          ELSEIF(N.EQ.NS+2.AND.ISFI(2).GE.1) THEN
            IF(THE2D.GT.THEFIS(2,ISFI(2))**2) GOTO 220
          ENDIF
        ENDIF
 
C...Option with angular ordering requirement.
        IF(MSTP(62).GE.3.AND.NTRY2.LT.200) THEN
          THE2T=(4D0*Z**2*Q2B)/(4D0*Z**2*Q2B+(1D0-Z)*XB**2*VINT2R)
          IF(THE2T.GT.THE2(JT)) GOTO 220
        ENDIF
      ENDIF
 
C...Weighting with new parton distributions.
      MINT(105)=MINT(102+JT)
      MINT(109)=MINT(106+JT)
      VINT(120)=VINT(2+JT)
      IF(MINT(31).GE.2) MINT(30)=JT
      IF(MSTP(57).LE.1) THEN
        CALL PYPDFU(KFBEAM(JT),XB,Q2REF,XFN)
      ELSE
        CALL PYPDFL(KFBEAM(JT),XB,Q2REF,XFN)
      ENDIF
      XFBN=XFN(KFLB)
      IF(XFBN.LT.1D-20) THEN
        IF(KFLA.EQ.KFLB) THEN
          TEVCB=TEVCBS
          TEVEB=TEVEBS
          WTAPC(KFLB)=0D0
          WTAPE(KFLB)=0D0
          GOTO 200
        ELSEIF(MCE.EQ.1.AND.TEVCBS-TEVCB.GT.0.2D0) THEN
          TEVCB=0.5D0*(TEVCBS+TEVCB)
          GOTO 230
        ELSEIF(MCE.EQ.2.AND.TEVEBS-TEVEB.GT.0.2D0) THEN
          TEVEB=0.5D0*(TEVEBS+TEVEB)
          GOTO 230
        ELSE
          XFBN=1D-10
          XFN(KFLB)=XFBN
        ENDIF
      ENDIF
      DO 250 KFL=-25,25
        XFB(KFL)=XFN(KFL)
  250 CONTINUE
      XA=XB/Z
      IF(MINT(31).GE.2) MINT(30)=JT
      IF(MSTP(57).LE.1) THEN
        CALL PYPDFU(KFBEAM(JT),XA,Q2REF,XFA)
      ELSE
        CALL PYPDFL(KFBEAM(JT),XA,Q2REF,XFA)
      ENDIF
      XFAN=XFA(KFLA)
      IF(XFAN.LT.1D-20) GOTO 200
      WTSFA=WTSF(KFLA)
      IF(WTZ*XFAN/XFBN.LT.PYR(0)*WTSFA) GOTO 200
 
C...Define two hard scatterers in their CM-frame.
  260 IF(N.EQ.NS+2) THEN
        DQ2(JT)=Q2B
        DPLCM=SQRT((DSH+DQ2(1)+DQ2(2))**2-4D0*DQ2(1)*DQ2(2))/DSHR
        DO 280 JR=1,2
          I=NS+JR
          IF(JR.EQ.1) IPO=IPUS1
          IF(JR.EQ.2) IPO=IPUS2
          DO 270 J=1,5
            K(I,J)=0
            P(I,J)=0D0
            V(I,J)=0D0
  270     CONTINUE
          K(I,1)=14
          K(I,2)=KFLS(JR+2)
          K(I,4)=IPO
          K(I,5)=IPO
          P(I,3)=DPLCM*(-1)**(JR+1)
          P(I,4)=(DSH+DQ2(3-JR)-DQ2(JR))/DSHR
          P(I,5)=-SQRT(DQ2(JR))
          K(IPO,1)=14
          K(IPO,3)=I
          K(IPO,4)=MOD(K(IPO,4),MSTU(5))+MSTU(5)*I
          K(IPO,5)=MOD(K(IPO,5),MSTU(5))+MSTU(5)*I
          MCT(I,1)=MCT(IPO,1)
          MCT(I,2)=MCT(IPO,2)
  280   CONTINUE
 
C...Find maximum allowed mass of timelike parton.
      ELSEIF(N.GT.NS+2) THEN
        JR=3-JT
        DQ2(3)=Q2B
        DPC(1)=P(IS(1),4)
        DPC(2)=P(IS(2),4)
        DPC(3)=0.5D0*(ABS(P(IS(1),3))+ABS(P(IS(2),3)))
        DPD(1)=DSH+DQ2(JR)+DQ2(JT)
        DPD(2)=DSHZ+DQ2(JR)+DQ2(3)
        DPD(3)=SQRT(DPD(1)**2-4D0*DQ2(JR)*DQ2(JT))
        DPD(4)=SQRT(DPD(2)**2-4D0*DQ2(JR)*DQ2(3))
        IKIN=0
        IF(Q2S(JR).GE.0.25D0*Q2MNC.AND.DPD(1)-DPD(3).GE.
     &  1D-10*DPD(1)) IKIN=1
        IF(IKIN.EQ.0) DMSMA=(DQ2(JT)/ZS(JT)-DQ2(3))*
     &  (DSH/(DSH+DQ2(JT))-DSH/(DSHZ+DQ2(3)))
        IF(IKIN.EQ.1) DMSMA=(DPD(1)*DPD(2)-DPD(3)*DPD(4))/
     &  (2D0*DQ2(JR))-DQ2(JT)-DQ2(3)
 
C...Generate timelike parton shower (if required).
        IT=N
        DO 290 J=1,5
          K(IT,J)=0
          P(IT,J)=0D0
          V(IT,J)=0D0
  290   CONTINUE
C...f -> f + g (gamma).
        IF(IABS(KFLB).LE.20.AND.IABS(KFLS(JT+2)).LE.20) THEN
          K(IT,2)=21
          IF(MCESV(JT).EQ.2.OR.IABS(KFLB).GE.11) K(IT,2)=22
C...f -> g (gamma, W+-) + f.
        ELSEIF(IABS(KFLB).LE.20.AND.IABS(KFLS(JT+2)).GT.20) THEN
          K(IT,2)=KFLB
          IF(KFLS(JT+2).EQ.24) THEN
            K(IT,2)=-12
          ELSEIF(KFLS(JT+2).EQ.-24) THEN
            K(IT,2)=12
          ENDIF
C...g (gamma) -> f + fbar, g + g.
        ELSE
          K(IT,2)=-KFLS(JT+2)
          IF(KFLS(JT+2).GT.20) K(IT,2)=KFLS(JT+2)
        ENDIF
        K(IT,1)=3
        IF((IABS(K(IT,2)).GE.11.AND.IABS(K(IT,2)).LE.18).OR.
     &  IABS(K(IT,2)).EQ.22) K(IT,1)=1
        P(IT,5)=PYMASS(K(IT,2))
        IF(DMSMA.LE.P(IT,5)**2) GOTO 100
        IF(MSTP(63).GE.1.AND.MCESV(JT).EQ.1) THEN
          MSTJ48=MSTJ(48)
          PARJ85=PARJ(85)
          P(IT,4)=(DSHZ-DSH-P(IT,5)**2)/DSHR
          P(IT,3)=SQRT(P(IT,4)**2-P(IT,5)**2)
          IF(MSTP(63).EQ.1) THEN
            Q2TIM=DMSMA
          ELSEIF(MSTP(63).EQ.2) THEN
            Q2TIM=MIN(DMSMA,PARP(71)*Q2S(JT))
          ELSE
            Q2TIM=DMSMA
            MSTJ(48)=1
            IF(IKIN.EQ.0) DPT2=DMSMA*(DSHZ+DQ2(3))/(DSH+DQ2(JT))
            IF(IKIN.EQ.1) DPT2=DMSMA*(0.5D0*DPD(1)*DPD(2)+0.5D0*DPD(3)*
     &      DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)))/(4D0*DSH*DPC(3)**2)
            PARJ(85)=SQRT(MAX(0D0,DPT2))*
     &      (1D0/P(IT,4)+1D0/P(IS(JT),4))
          ENDIF
C...Only do timelike shower here if using PYSHOW
          IF (MSTJ(41).NE.11.AND.MSTJ(41).NE.12) THEN
            CALL PYSHOW(IT,0,SQRT(Q2TIM))
          ENDIF
          MSTJ(48)=MSTJ48
          PARJ(85)=PARJ85
          IF(N.GE.IT+1) P(IT,5)=P(IT+1,5)
        ENDIF
 
C...Reconstruct kinematics of branching: timelike parton shower.
        DMS=P(IT,5)**2
        IF(IKIN.EQ.0) DPT2=(DMSMA-DMS)*(DSHZ+DQ2(3))/(DSH+DQ2(JT))
        IF(IKIN.EQ.1) DPT2=(DMSMA-DMS)*(0.5D0*DPD(1)*DPD(2)+
     &  0.5D0*DPD(3)*DPD(4)-DQ2(JR)*(DQ2(JT)+DQ2(3)+DMS))/
     &  (4D0*DSH*DPC(3)**2)
        IF(DPT2.LT.0D0) GOTO 100
        DPB(1)=(0.5D0*DPD(2)-DPC(JR)*(DSHZ+DQ2(JR)-DQ2(JT)-DMS)/
     &  DSHR)/DPC(3)-DPC(3)
        P(IT,1)=SQRT(DPT2)
        P(IT,3)=DPB(1)*(-1)**(JT+1)
        P(IT,4)=SQRT(DPT2+DPB(1)**2+DMS)
        IF(N.GE.IT+1) THEN
          DPB(1)=SQRT(DPB(1)**2+DPT2)
          DPB(2)=SQRT(DPB(1)**2+DMS)
          DPB(3)=P(IT+1,3)
          DPB(4)=SQRT(DPB(3)**2+DMS)
          DBEZ=(DPB(4)*DPB(1)-DPB(3)*DPB(2))/(DPB(4)*DPB(2)-DPB(3)*
     &    DPB(1))
          CALL PYROBO(IT+1,N,0D0,0D0,0D0,0D0,DBEZ)
          THE=PYANGL(P(IT,3),P(IT,1))
          CALL PYROBO(IT+1,N,THE,0D0,0D0,0D0,0D0)
        ENDIF
 
C...Reconstruct kinematics of branching: spacelike parton.
        DO 300 J=1,5
          K(N+1,J)=0
          P(N+1,J)=0D0
          V(N+1,J)=0D0
  300   CONTINUE
        K(N+1,1)=14
        K(N+1,2)=KFLB
        P(N+1,1)=P(IT,1)
        P(N+1,3)=P(IT,3)+P(IS(JT),3)
        P(N+1,4)=P(IT,4)+P(IS(JT),4)
        P(N+1,5)=-SQRT(DQ2(3))
        MCT(N+1,1)=0
        MCT(N+1,2)=0
 
C...Define colour flow of branching.
        K(IS(JT),3)=N+1
        K(IT,3)=N+1
        IM1=N+1
        IM2=N+1
C...f -> f + gamma (Z, W).
        IF(IABS(K(IT,2)).GE.22) THEN
          K(IT,1)=1
          ID1=IS(JT)
          ID2=IS(JT)
C...f -> gamma (Z, W) + f.
        ELSEIF(IABS(K(IS(JT),2)).GE.22) THEN
          ID1=IT
          ID2=IT
C...gamma -> q + qbar, g + g.
        ELSEIF(K(N+1,2).EQ.22) THEN
          ID1=IS(JT)
          ID2=IT
          IM1=ID2
          IM2=ID1
C...q -> q + g.
        ELSEIF(K(N+1,2).GT.0.AND.K(N+1,2).NE.21.AND.K(IT,2).EQ.21) THEN
          ID1=IT
          ID2=IS(JT)
C...q -> g + q.
        ELSEIF(K(N+1,2).GT.0.AND.K(N+1,2).NE.21) THEN
          ID1=IS(JT)
          ID2=IT
C...qbar -> qbar + g.
        ELSEIF(K(N+1,2).LT.0.AND.K(IT,2).EQ.21) THEN
          ID1=IS(JT)
          ID2=IT
C...qbar -> g + qbar.
        ELSEIF(K(N+1,2).LT.0) THEN
          ID1=IT
          ID2=IS(JT)
C...g -> g + g; g -> q + qbar.
        ELSEIF((K(IT,2).EQ.21.AND.PYR(0).GT.0.5D0).OR.K(IT,2).LT.0) THEN
          ID1=IS(JT)
          ID2=IT
        ELSE
          ID1=IT
          ID2=IS(JT)
        ENDIF
        IF(IM1.EQ.N+1) K(IM1,4)=K(IM1,4)+ID1
        IF(IM2.EQ.N+1) K(IM2,5)=K(IM2,5)+ID2
        K(ID1,4)=K(ID1,4)+MSTU(5)*IM1
        K(ID2,5)=K(ID2,5)+MSTU(5)*IM2
        IF(ID1.NE.ID2) THEN
          K(ID1,5)=K(ID1,5)+MSTU(5)*ID2
          K(ID2,4)=K(ID2,4)+MSTU(5)*ID1
        ENDIF
        N=N+1
        IF(K(IT,1).EQ.1) THEN
          K(IT,4)=0
          K(IT,5)=0
        ENDIF
 
C...Boost to new CM-frame.
        DBSVX=(P(N,1)+P(IS(JR),1))/(P(N,4)+P(IS(JR),4))
        DBSVZ=(P(N,3)+P(IS(JR),3))/(P(N,4)+P(IS(JR),4))
        IF(DBSVX**2+DBSVZ**2.GE.1D0) GOTO 100
        CALL PYROBO(NS+1,N,0D0,0D0,-DBSVX,0D0,-DBSVZ)
        IR=N+(JT-1)*(IS(1)-N)
        CALL PYROBO(NS+1,N,-PYANGL(P(IR,3),P(IR,1)),DPHI(JT),
     &  0D0,0D0,0D0)
 
C...Save timelike parton in PYPART if doing pT-ordered FSR off ISR
        IF (MSTJ(41).EQ.11.OR.MSTJ(41).EQ.12) THEN
          NPART=NPART+1
          IPART(NPART)=IT
          PTPART(NPART)=SQRT(PARP(71)*DPT2)
        ENDIF

C...Global statistics.
        MINT(352)=MINT(352)+1
        VINT(352)=VINT(352)+SQRT(P(IT,1)**2+P(IT,2)**2)
        IF (MINT(352).EQ.1) VINT(357)=SQRT(P(IT,1)**2+P(IT,2)**2)

      ENDIF
 
C...Update kinematics variables.
      IS(JT)=N
      DQ2(JT)=Q2B
      IF(MSTP(62).GE.3.AND.NTRY2.LT.200.AND.MCE.EQ.1) THE2(JT)=THE2T
      DSH=DSHZ
 
C...Save quantities; loop back.
      Q2S(JT)=Q2B
      DPHI(JT)=PHIBR
      MCESV(JT)=MCE
      IF((MCEV.EQ.1.AND.Q2B.GE.0.25D0*Q2MNC).OR.
     &(MEEV.EQ.1.AND.Q2B.GE.Q2MNE)) THEN
        KFLS(JT+2)=KFLS(JT)
        KFLS(JT)=KFLA
        XS(JT)=XA
        ZS(JT)=Z
        DO 310 KFL=-25,25
          XFS(JT,KFL)=XFA(KFL)
  310   CONTINUE
        TEVCSV(JT)=TEVCB
        TEVESV(JT)=TEVEB
      ELSE
        MORE(JT)=0
        IF(JT.EQ.1) IPU1=N
        IF(JT.EQ.2) IPU2=N
      ENDIF
      IF(N.GT.MSTU(4)-MSTU(32)-10) THEN
        CALL PYERRM(11,'(PYSSPA:) no more memory left in PYJETS')
        IF(MSTU(21).GE.1) N=NS
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      IF(MORE(1).EQ.1.OR.MORE(2).EQ.1) GOTO 150
 
C...Boost hard scattering partons to frame of shower initiators.
      DO 320 J=1,3
        ROBO(J+2)=(P(NS+1,J)+P(NS+2,J))/(P(NS+1,4)+P(NS+2,4))
  320 CONTINUE
      K(N+2,1)=1
      DO 330 J=1,5
        P(N+2,J)=P(NS+1,J)
  330 CONTINUE
      CALL PYROBO(N+2,N+2,0D0,0D0,-ROBO(3),-ROBO(4),-ROBO(5))
      ROBO(2)=PYANGL(P(N+2,1),P(N+2,2))
      ROBO(1)=PYANGL(P(N+2,3),SQRT(P(N+2,1)**2+P(N+2,2)**2))
      IMIN=MINT(83)+5
      IF(MINT(31).GE.2) IMIN=MIN(IPUS1,IPUS2)
      CALL PYROBO(IMIN,NS,0D0,-ROBO(2),0D0,0D0,0D0)
      CALL PYROBO(IMIN,NS,ROBO(1),ROBO(2),ROBO(3),ROBO(4),ROBO(5))
 
C...Store user information. Reset Lambda value.
      IF(MINT(31).LE.1) THEN
        K(IPU1,3)=MINT(83)+3
        K(IPU2,3)=MINT(83)+4
      ELSE
        K(IPU1,3)=MINT(83)+1
        K(IPU2,3)=MINT(83)+2
      ENDIF
      DO 340 JT=1,2
        MINT(12+JT)=KFLS(JT)
        VINT(140+JT)=XS(JT)
        IF(MINT(18+JT).EQ.1) VINT(140+JT)=VINT(154+JT)*XS(JT)
        IF(MINT(31).GE.2) VINT(140+JT)=VINT(140+JT)*VINT(142+JT)
  340 CONTINUE
      PARU(112)=ALAMS
 
      RETURN
      END
