 
C*********************************************************************
 
C...PYSHOW
C...Generates timelike parton showers from given partons.
 
      SUBROUTINE PYSHOW(IP1,IP2,QMAX)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
      PARAMETER (MAXNUR=1000)
C...Commonblocks.
      COMMON/PYPART/NPART,NPARTD,IPART(MAXNUR),PTPART(MAXNUR)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYPART/,/PYJETS/,/PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/
C...Local arrays.
      DIMENSION PMTH(5,140),PS(5),PMA(100),PMSD(100),IEP(100),IPA(100),
     &KFLA(100),KFLD(100),KFL(100),ITRY(100),ISI(100),ISL(100),DP(100),
     &DPT(5,4),KSH(0:140),KCII(2),NIIS(2),IIIS(2,2),THEIIS(2,2),
     &PHIIIS(2,2),ISII(2),ISSET(2),ISCOL(0:140),ISCHG(0:140),
     &IREF(1000)
 
C...Check that QMAX not too low.
      IF(MSTJ(41).LE.0) THEN
        RETURN
      ELSEIF(MSTJ(41).EQ.1.OR.MSTJ(41).EQ.11) THEN
        IF(QMAX.LE.PARJ(82).AND.IP2.GE.-80) RETURN
      ELSE
        IF(QMAX.LE.MIN(PARJ(82),PARJ(83),PARJ(90)).AND.IP2.GE.-80)
     &  RETURN
      ENDIF
 
C...Store positions of shower initiating partons.
      MPSPD=0
      IF(IP1.GT.0.AND.IP1.LE.MIN(N,MSTU(4)-MSTU(32)).AND.IP2.EQ.0) THEN
        NPA=1
        IPA(1)=IP1
      ELSEIF(MIN(IP1,IP2).GT.0.AND.MAX(IP1,IP2).LE.MIN(N,MSTU(4)-
     &  MSTU(32))) THEN
        NPA=2
        IPA(1)=IP1
        IPA(2)=IP2
      ELSEIF(IP1.GT.0.AND.IP1.LE.MIN(N,MSTU(4)-MSTU(32)).AND.IP2.LT.0
     &  .AND.IP2.GE.-80) THEN
        NPA=IABS(IP2)
        DO 100 I=1,NPA
          IPA(I)=IP1+I-1
  100   CONTINUE
      ELSEIF(IP1.GT.0.AND.IP1.LE.MIN(N,MSTU(4)-MSTU(32)).AND.
     &IP2.EQ.-100) THEN
        MPSPD=1
        NPA=2
        IPA(1)=IP1+6
        IPA(2)=IP1+7
      ELSE
        CALL PYERRM(12,
     &  '(PYSHOW:) failed to reconstruct showering system')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
 
C...Send off to PYPTFS for pT-ordered evolution if requested,
C...if at least 2 partons, and without predefined shower branchings.
      IF((MSTJ(41).EQ.11.OR.MSTJ(41).EQ.12).AND.NPA.GE.2.AND.
     &MPSPD.EQ.0) THEN
        NPART=NPA
        DO 110 II=1,NPART
          IPART(II)=IPA(II)
          PTPART(II)=0.5D0*QMAX
  110   CONTINUE
        CALL PYPTFS(2,0.5D0*QMAX,0D0,PTGEN)
        RETURN
      ENDIF
 
C...Initialization of cutoff masses etc.
      DO 120 IFL=0,40
        ISCOL(IFL)=0
        ISCHG(IFL)=0
        KSH(IFL)=0
  120 CONTINUE
      ISCOL(21)=1
      KSH(21)=1
      PMTH(1,21)=PYMASS(21)
      PMTH(2,21)=SQRT(PMTH(1,21)**2+0.25D0*PARJ(82)**2)
      PMTH(3,21)=2D0*PMTH(2,21)
      PMTH(4,21)=PMTH(3,21)
      PMTH(5,21)=PMTH(3,21)
      PMTH(1,22)=PYMASS(22)
      PMTH(2,22)=SQRT(PMTH(1,22)**2+0.25D0*PARJ(83)**2)
      PMTH(3,22)=2D0*PMTH(2,22)
      PMTH(4,22)=PMTH(3,22)
      PMTH(5,22)=PMTH(3,22)
      PMQTH1=PARJ(82)
      IF(MSTJ(41).GE.2) PMQTH1=MIN(PARJ(82),PARJ(83))
      PMQT1E=MIN(PMQTH1,PARJ(90))
      PMQTH2=PMTH(2,21)
      IF(MSTJ(41).GE.2) PMQTH2=MIN(PMTH(2,21),PMTH(2,22))
      PMQT2E=MIN(PMQTH2,0.5D0*PARJ(90))
      DO 130 IFL=1,5
        ISCOL(IFL)=1
        IF(MSTJ(41).GE.2) ISCHG(IFL)=1
        KSH(IFL)=1
        PMTH(1,IFL)=PYMASS(IFL)
        PMTH(2,IFL)=SQRT(PMTH(1,IFL)**2+0.25D0*PMQTH1**2)
        PMTH(3,IFL)=PMTH(2,IFL)+PMQTH2
        PMTH(4,IFL)=SQRT(PMTH(1,IFL)**2+0.25D0*PARJ(82)**2)+PMTH(2,21)
        PMTH(5,IFL)=SQRT(PMTH(1,IFL)**2+0.25D0*PARJ(83)**2)+PMTH(2,22)
  130 CONTINUE
      DO 140 IFL=11,15,2
        IF(MSTJ(41).EQ.2.OR.MSTJ(41).GE.4) ISCHG(IFL)=1
        IF(MSTJ(41).EQ.2.OR.MSTJ(41).GE.4) KSH(IFL)=1
        PMTH(1,IFL)=PYMASS(IFL)
        PMTH(2,IFL)=SQRT(PMTH(1,IFL)**2+0.25D0*PARJ(90)**2)
        PMTH(3,IFL)=PMTH(2,IFL)+0.5D0*PARJ(90)
        PMTH(4,IFL)=PMTH(3,IFL)
        PMTH(5,IFL)=PMTH(3,IFL)
  140 CONTINUE
      PT2MIN=MAX(0.5D0*PARJ(82),1.1D0*PARJ(81))**2
      ALAMS=PARJ(81)**2
      ALFM=LOG(PT2MIN/ALAMS)
 
C...Check on phase space available for emission.
      IREJ=0
      DO 150 J=1,5
        PS(J)=0D0
  150 CONTINUE
      PM=0D0
      KFLA(2)=0
      DO 170 I=1,NPA
        KFLA(I)=IABS(K(IPA(I),2))
        PMA(I)=P(IPA(I),5)
C...Special cutoff masses for initial partons (may be a heavy quark,
C...squark, ..., and need not be on the mass shell).
        IR=30+I
        IF(NPA.LE.1) IREF(I)=IR
        IF(NPA.GE.2) IREF(I+1)=IR
        ISCOL(IR)=0
        ISCHG(IR)=0
        KSH(IR)=0
        IF(KFLA(I).LE.8) THEN
          ISCOL(IR)=1
          IF(MSTJ(41).GE.2) ISCHG(IR)=1
        ELSEIF(KFLA(I).EQ.11.OR.KFLA(I).EQ.13.OR.KFLA(I).EQ.15.OR.
     &  KFLA(I).EQ.17) THEN
          IF(MSTJ(41).EQ.2.OR.MSTJ(41).GE.4) ISCHG(IR)=1
        ELSEIF(KFLA(I).EQ.21) THEN
          ISCOL(IR)=1
        ELSEIF((KFLA(I).GE.KSUSY1+1.AND.KFLA(I).LE.KSUSY1+8).OR.
     &  (KFLA(I).GE.KSUSY2+1.AND.KFLA(I).LE.KSUSY2+8)) THEN
          ISCOL(IR)=1
        ELSEIF(KFLA(I).EQ.KSUSY1+21) THEN
          ISCOL(IR)=1
C...QUARKONIA+++
C...same for QQ~[3S18]
        ELSEIF(MSTP(148).GE.1.AND.(KFLA(I).EQ.9900443.OR.
     &  KFLA(I).EQ.9900553)) THEN
          ISCOL(IR)=1
C...QUARKONIA---
        ENDIF

C...Option to switch off radiation from particle KF = MSTJ(39) entirely
C...(only intended for studying the effects of switching such rad on/off)
        IF (MSTJ(39).GT.0.AND.KFLA(I).EQ.MSTJ(39)) THEN
          ISCOL(IR)=0
          ISCHG(IR)=0
        ENDIF

        IF(ISCOL(IR).EQ.1.OR.ISCHG(IR).EQ.1) KSH(IR)=1
        PMTH(1,IR)=PMA(I)
        IF(ISCOL(IR).EQ.1.AND.ISCHG(IR).EQ.1) THEN
          PMTH(2,IR)=SQRT(PMTH(1,IR)**2+0.25D0*PMQTH1**2)
          PMTH(3,IR)=PMTH(2,IR)+PMQTH2
          PMTH(4,IR)=SQRT(PMTH(1,IR)**2+0.25D0*PARJ(82)**2)+PMTH(2,21)
          PMTH(5,IR)=SQRT(PMTH(1,IR)**2+0.25D0*PARJ(83)**2)+PMTH(2,22)
        ELSEIF(ISCOL(IR).EQ.1) THEN
          PMTH(2,IR)=SQRT(PMTH(1,IR)**2+0.25D0*PARJ(82)**2)
          PMTH(3,IR)=PMTH(2,IR)+0.5D0*PARJ(82)
          PMTH(4,IR)=PMTH(3,IR)
          PMTH(5,IR)=PMTH(3,IR)
        ELSEIF(ISCHG(IR).EQ.1) THEN
          PMTH(2,IR)=SQRT(PMTH(1,IR)**2+0.25D0*PARJ(90)**2)
          PMTH(3,IR)=PMTH(2,IR)+0.5D0*PARJ(90)
          PMTH(4,IR)=PMTH(3,IR)
          PMTH(5,IR)=PMTH(3,IR)
        ENDIF
        IF(KSH(IR).EQ.1) PMA(I)=PMTH(3,IR)
        PM=PM+PMA(I)
        IF(KSH(IR).EQ.0.OR.PMA(I).GT.10D0*QMAX) IREJ=IREJ+1
        DO 160 J=1,4
          PS(J)=PS(J)+P(IPA(I),J)
  160   CONTINUE
  170 CONTINUE
      IF(IREJ.EQ.NPA.AND.IP2.GE.-7) RETURN
      PS(5)=SQRT(MAX(0D0,PS(4)**2-PS(1)**2-PS(2)**2-PS(3)**2))
      IF(NPA.EQ.1) PS(5)=PS(4)
      IF(PS(5).LE.PM+PMQT1E) RETURN
 
C...Identify source: q(1), ~q(2), V(3), S(4), chi(5), ~g(6), unknown(0).
      KFSRCE=0
      IF(IP2.LE.0) THEN
      ELSEIF(K(IP1,3).EQ.K(IP2,3).AND.K(IP1,3).GT.0) THEN
        KFSRCE=IABS(K(K(IP1,3),2))
      ELSE
        IPAR1=MAX(1,K(IP1,3))
        IPAR2=MAX(1,K(IP2,3))
        IF(K(IPAR1,3).EQ.K(IPAR2,3).AND.K(IPAR1,3).GT.0)
     &       KFSRCE=IABS(K(K(IPAR1,3),2))
      ENDIF
      ITYPES=0
      IF(KFSRCE.GE.1.AND.KFSRCE.LE.8) ITYPES=1
      IF(KFSRCE.GE.KSUSY1+1.AND.KFSRCE.LE.KSUSY1+8) ITYPES=2
      IF(KFSRCE.GE.KSUSY2+1.AND.KFSRCE.LE.KSUSY2+8) ITYPES=2
      IF(KFSRCE.GE.21.AND.KFSRCE.LE.24) ITYPES=3
      IF(KFSRCE.GE.32.AND.KFSRCE.LE.34) ITYPES=3
      IF(KFSRCE.EQ.25.OR.(KFSRCE.GE.35.AND.KFSRCE.LE.37)) ITYPES=4
      IF(KFSRCE.GE.KSUSY1+22.AND.KFSRCE.LE.KSUSY1+37) ITYPES=5
      IF(KFSRCE.EQ.KSUSY1+21) ITYPES=6
 
C...Identify two primary showerers.
      ITYPE1=0
      IF(KFLA(1).GE.1.AND.KFLA(1).LE.8) ITYPE1=1
      IF(KFLA(1).GE.KSUSY1+1.AND.KFLA(1).LE.KSUSY1+8) ITYPE1=2
      IF(KFLA(1).GE.KSUSY2+1.AND.KFLA(1).LE.KSUSY2+8) ITYPE1=2
      IF(KFLA(1).GE.21.AND.KFLA(1).LE.24) ITYPE1=3
      IF(KFLA(1).GE.32.AND.KFLA(1).LE.34) ITYPE1=3
      IF(KFLA(1).EQ.25.OR.(KFLA(1).GE.35.AND.KFLA(1).LE.37)) ITYPE1=4
      IF(KFLA(1).GE.KSUSY1+22.AND.KFLA(1).LE.KSUSY1+37) ITYPE1=5
      IF(KFLA(1).EQ.KSUSY1+21) ITYPE1=6
      ITYPE2=0
      IF(KFLA(2).GE.1.AND.KFLA(2).LE.8) ITYPE2=1
      IF(KFLA(2).GE.KSUSY1+1.AND.KFLA(2).LE.KSUSY1+8) ITYPE2=2
      IF(KFLA(2).GE.KSUSY2+1.AND.KFLA(2).LE.KSUSY2+8) ITYPE2=2
      IF(KFLA(2).GE.21.AND.KFLA(2).LE.24) ITYPE2=3
      IF(KFLA(2).GE.32.AND.KFLA(2).LE.34) ITYPE2=3
      IF(KFLA(2).EQ.25.OR.(KFLA(2).GE.35.AND.KFLA(2).LE.37)) ITYPE2=4
      IF(KFLA(2).GE.KSUSY1+22.AND.KFLA(2).LE.KSUSY1+37) ITYPE2=5
      IF(KFLA(2).EQ.KSUSY1+21) ITYPE2=6
 
C...Order of showerers. Presence of gluino.
      ITYPMN=MIN(ITYPE1,ITYPE2)
      ITYPMX=MAX(ITYPE1,ITYPE2)
      IORD=1
      IF(ITYPE1.GT.ITYPE2) IORD=2
      IGLUI=0
      IF(ITYPE1.EQ.6.OR.ITYPE2.EQ.6) IGLUI=1
 
C...Check if 3-jet matrix elements to be used.
      M3JC=0
      ALPHA=0.5D0
      IF(NPA.EQ.2.AND.MSTJ(47).GE.1.AND.MPSPD.EQ.0) THEN
        IF(MSTJ(38).NE.0) THEN
          M3JC=MSTJ(38)
          ALPHA=PARJ(80)
          MSTJ(38)=0
        ELSEIF(MSTJ(47).GE.6) THEN
          M3JC=MSTJ(47)
        ELSE
          ICLASS=1
          ICOMBI=4
 
C...Vector/axial vector -> q + qbar; q -> q + V.
          IF(ITYPMN.EQ.1.AND.ITYPMX.EQ.1.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.3)) THEN
            ICLASS=2
            IF(KFSRCE.EQ.21.OR.KFSRCE.EQ.22) THEN
              ICOMBI=1
            ELSEIF(KFSRCE.EQ.23.OR.(KFSRCE.EQ.0.AND.
     &      K(IPA(1),2)+K(IPA(2),2).EQ.0)) THEN
C...gamma*/Z0: assume e+e- initial state if unknown.
              EI=-1D0
              IF(KFSRCE.EQ.23) THEN
                IANNFL=K(K(IP1,3),3)
                IF(IANNFL.NE.0) THEN
                  KANNFL=IABS(K(IANNFL,2))
                  IF(KANNFL.GE.1.AND.KANNFL.LE.18) EI=KCHG(KANNFL,1)/3D0
                ENDIF
              ENDIF
              AI=SIGN(1D0,EI+0.1D0)
              VI=AI-4D0*EI*PARU(102)
              EF=KCHG(KFLA(1),1)/3D0
              AF=SIGN(1D0,EF+0.1D0)
              VF=AF-4D0*EF*PARU(102)
              XWC=1D0/(16D0*PARU(102)*(1D0-PARU(102)))
              SH=PS(5)**2
              SQMZ=PMAS(23,1)**2
              SQWZ=PS(5)*PMAS(23,2)
              SBWZ=1D0/((SH-SQMZ)**2+SQWZ**2)
              VECT=EI**2*EF**2+2D0*EI*VI*EF*VF*XWC*SH*(SH-SQMZ)*SBWZ+
     &        (VI**2+AI**2)*VF**2*XWC**2*SH**2*SBWZ
              AXIV=(VI**2+AI**2)*AF**2*XWC**2*SH**2*SBWZ
              ICOMBI=3
              ALPHA=VECT/(VECT+AXIV)
            ELSEIF(KFSRCE.EQ.24.OR.KFSRCE.EQ.0) THEN
              ICOMBI=4
            ENDIF
C...For chi -> chi q qbar, use V/A -> q qbar as first approximation.
          ELSEIF(ITYPMN.EQ.1.AND.ITYPMX.EQ.1.AND.ITYPES.EQ.5) THEN
            ICLASS=2
          ELSEIF(ITYPMN.EQ.1.AND.ITYPMX.EQ.3.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.1)) THEN
            ICLASS=3
 
C...Scalar/pseudoscalar -> q + qbar; q -> q + S.
          ELSEIF(ITYPMN.EQ.1.AND.ITYPMX.EQ.1.AND.ITYPES.EQ.4) THEN
            ICLASS=4
            IF(KFSRCE.EQ.25.OR.KFSRCE.EQ.35.OR.KFSRCE.EQ.37) THEN
              ICOMBI=1
            ELSEIF(KFSRCE.EQ.36) THEN
              ICOMBI=2
            ENDIF
          ELSEIF(ITYPMN.EQ.1.AND.ITYPMX.EQ.4.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.1)) THEN
            ICLASS=5
 
C...V -> ~q + ~qbar; ~q -> ~q + V; S -> ~q + ~qbar; ~q -> ~q + S.
          ELSEIF(ITYPMN.EQ.2.AND.ITYPMX.EQ.2.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.3)) THEN
            ICLASS=6
          ELSEIF(ITYPMN.EQ.2.AND.ITYPMX.EQ.3.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.2)) THEN
            ICLASS=7
          ELSEIF(ITYPMN.EQ.2.AND.ITYPMX.EQ.2.AND.ITYPES.EQ.4) THEN
            ICLASS=8
          ELSEIF(ITYPMN.EQ.2.AND.ITYPMX.EQ.4.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.2)) THEN
            ICLASS=9
 
C...chi -> q + ~qbar; ~q -> q + chi; q -> ~q + chi.
          ELSEIF(ITYPMN.EQ.1.AND.ITYPMX.EQ.2.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.5)) THEN
            ICLASS=10
          ELSEIF(ITYPMN.EQ.1.AND.ITYPMX.EQ.5.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.2)) THEN
            ICLASS=11
          ELSEIF(ITYPMN.EQ.2.AND.ITYPMX.EQ.5.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.1)) THEN
            ICLASS=12
 
C...~g -> q + ~qbar; ~q -> q + ~g; q -> ~q + ~g.
          ELSEIF(ITYPMN.EQ.1.AND.ITYPMX.EQ.2.AND.ITYPES.EQ.6) THEN
            ICLASS=13
          ELSEIF(ITYPMN.EQ.1.AND.ITYPMX.EQ.6.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.2)) THEN
            ICLASS=14
          ELSEIF(ITYPMN.EQ.2.AND.ITYPMX.EQ.6.AND.(ITYPES.EQ.0.OR.
     &    ITYPES.EQ.1)) THEN
            ICLASS=15
 
C...g -> ~g + ~g (eikonal approximation).
          ELSEIF(ITYPMN.EQ.6.AND.ITYPMX.EQ.6.AND.ITYPES.EQ.0) THEN
            ICLASS=16
          ENDIF

C...Revert to eikonal approximation for gluon in final state.
          IF(KFLA1.EQ.21.OR.KFLA2.EQ.21) ICLASS=1 

          M3JC=5*ICLASS+ICOMBI
        ENDIF
      ENDIF
 
C...Find if interference with initial state partons.
      MIIS=0
      IF(MSTJ(50).GE.1.AND.MSTJ(50).LE.3.AND.NPA.EQ.2.AND.KFSRCE.EQ.0
     &.AND.MPSPD.EQ.0) MIIS=MSTJ(50)
      IF(MSTJ(50).GE.4.AND.MSTJ(50).LE.6.AND.NPA.EQ.2.AND.MPSPD.EQ.0)
     &MIIS=MSTJ(50)-3
      IF(MIIS.NE.0) THEN
        DO 190 I=1,2
          KCII(I)=0
          KCA=PYCOMP(KFLA(I))
          IF(KCA.NE.0) KCII(I)=KCHG(KCA,2)*ISIGN(1,K(IPA(I),2))
          NIIS(I)=0
          IF(KCII(I).NE.0) THEN
            DO 180 J=1,2
              ICSI=MOD(K(IPA(I),3+J)/MSTU(5),MSTU(5))
              IF(ICSI.GT.0.AND.ICSI.NE.IPA(1).AND.ICSI.NE.IPA(2).AND.
     &        (KCII(I).EQ.(-1)**(J+1).OR.KCII(I).EQ.2)) THEN
                NIIS(I)=NIIS(I)+1
                IIIS(I,NIIS(I))=ICSI
              ENDIF
  180       CONTINUE
          ENDIF
  190   CONTINUE
        IF(NIIS(1)+NIIS(2).EQ.0) MIIS=0
      ENDIF
 
C...Boost interfering initial partons to rest frame
C...and reconstruct their polar and azimuthal angles.
      IF(MIIS.NE.0) THEN
        DO 210 I=1,2
          DO 200 J=1,5
            K(N+I,J)=K(IPA(I),J)
            P(N+I,J)=P(IPA(I),J)
            V(N+I,J)=0D0
  200     CONTINUE
  210   CONTINUE
        DO 230 I=3,2+NIIS(1)
          DO 220 J=1,5
            K(N+I,J)=K(IIIS(1,I-2),J)
            P(N+I,J)=P(IIIS(1,I-2),J)
            V(N+I,J)=0D0
  220     CONTINUE
  230   CONTINUE
        DO 250 I=3+NIIS(1),2+NIIS(1)+NIIS(2)
          DO 240 J=1,5
            K(N+I,J)=K(IIIS(2,I-2-NIIS(1)),J)
            P(N+I,J)=P(IIIS(2,I-2-NIIS(1)),J)
            V(N+I,J)=0D0
  240     CONTINUE
  250   CONTINUE
        CALL PYROBO(N+1,N+2+NIIS(1)+NIIS(2),0D0,0D0,-PS(1)/PS(4),
     &  -PS(2)/PS(4),-PS(3)/PS(4))
        PHI=PYANGL(P(N+1,1),P(N+1,2))
        CALL PYROBO(N+1,N+2+NIIS(1)+NIIS(2),0D0,-PHI,0D0,0D0,0D0)
        THE=PYANGL(P(N+1,3),P(N+1,1))
        CALL PYROBO(N+1,N+2+NIIS(1)+NIIS(2),-THE,0D0,0D0,0D0,0D0)
        DO 260 I=3,2+NIIS(1)
          THEIIS(1,I-2)=PYANGL(P(N+I,3),SQRT(P(N+I,1)**2+P(N+I,2)**2))
          PHIIIS(1,I-2)=PYANGL(P(N+I,1),P(N+I,2))
  260   CONTINUE
        DO 270 I=3+NIIS(1),2+NIIS(1)+NIIS(2)
          THEIIS(2,I-2-NIIS(1))=PARU(1)-PYANGL(P(N+I,3),
     &    SQRT(P(N+I,1)**2+P(N+I,2)**2))
          PHIIIS(2,I-2-NIIS(1))=PYANGL(P(N+I,1),P(N+I,2))
  270   CONTINUE
      ENDIF
 
C...Boost 3 or more partons to their rest frame.
      IF(NPA.GE.3) CALL PYROBO(IPA(1),IPA(NPA),0D0,0D0,-PS(1)/PS(4),
     &-PS(2)/PS(4),-PS(3)/PS(4))
 
C...Define imagined single initiator of shower for parton system.
      NS=N
      IF(N.GT.MSTU(4)-MSTU(32)-10) THEN
        CALL PYERRM(11,'(PYSHOW:) no more memory left in PYJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
  280 N=NS
      IF(NPA.GE.2) THEN
        K(N+1,1)=11
        K(N+1,2)=21
        K(N+1,3)=0
        K(N+1,4)=0
        K(N+1,5)=0
        P(N+1,1)=0D0
        P(N+1,2)=0D0
        P(N+1,3)=0D0
        P(N+1,4)=PS(5)
        P(N+1,5)=PS(5)
        V(N+1,5)=PS(5)**2
        N=N+1
        IREF(1)=21
      ENDIF
 
C...Loop over partons that may branch.
      NEP=NPA
      IM=NS
      IF(NPA.EQ.1) IM=NS-1
  290 IM=IM+1
      IF(N.GT.NS) THEN
        IF(IM.GT.N) GOTO 600
        KFLM=IABS(K(IM,2))
        IR=IREF(IM-NS)
        IF(KSH(IR).EQ.0) GOTO 290
        IF(P(IM,5).LT.PMTH(2,IR)) GOTO 290
        IGM=K(IM,3)
      ELSE
        IGM=-1
      ENDIF
      IF(N+NEP.GT.MSTU(4)-MSTU(32)-10) THEN
        CALL PYERRM(11,'(PYSHOW:) no more memory left in PYJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
 
C...Position of aunt (sister to branching parton).
C...Origin and flavour of daughters.
      IAU=0
      IF(IGM.GT.0) THEN
        IF(K(IM-1,3).EQ.IGM) IAU=IM-1
        IF(N.GE.IM+1.AND.K(IM+1,3).EQ.IGM) IAU=IM+1
      ENDIF
      IF(IGM.GE.0) THEN
        K(IM,4)=N+1
        DO 300 I=1,NEP
          K(N+I,3)=IM
  300   CONTINUE
      ELSE
        K(N+1,3)=IPA(1)
      ENDIF
      IF(IGM.LE.0) THEN
        DO 310 I=1,NEP
          K(N+I,2)=K(IPA(I),2)
  310   CONTINUE
      ELSEIF(KFLM.NE.21) THEN
        K(N+1,2)=K(IM,2)
        K(N+2,2)=K(IM,5)
        IREF(N+1-NS)=IREF(IM-NS)
        IREF(N+2-NS)=IABS(K(N+2,2))
      ELSEIF(K(IM,5).EQ.21) THEN
        K(N+1,2)=21
        K(N+2,2)=21
        IREF(N+1-NS)=21
        IREF(N+2-NS)=21
      ELSE
        K(N+1,2)=K(IM,5)
        K(N+2,2)=-K(IM,5)
        IREF(N+1-NS)=IABS(K(N+1,2))
        IREF(N+2-NS)=IABS(K(N+2,2))
      ENDIF
 
C...Reset flags on daughters and tries made.
      DO 320 IP=1,NEP
        K(N+IP,1)=3
        K(N+IP,4)=0
        K(N+IP,5)=0
        KFLD(IP)=IABS(K(N+IP,2))
        IF(KCHG(PYCOMP(KFLD(IP)),2).EQ.0) K(N+IP,1)=1
        ITRY(IP)=0
        ISL(IP)=0
        ISI(IP)=0
        IF(KSH(IREF(N+IP-NS)).EQ.1) ISI(IP)=1
  320 CONTINUE
      ISLM=0
 
C...Maximum virtuality of daughters.
      IF(IGM.LE.0) THEN
        DO 330 I=1,NPA
          IF(NPA.GE.3) P(N+I,4)=P(IPA(I),4)
          P(N+I,5)=MIN(QMAX,PS(5))
          IR=IREF(N+I-NS)
          IF(IP2.LE.-8) P(N+I,5)=MAX(P(N+I,5),2D0*PMTH(3,IR))
          IF(ISI(I).EQ.0) P(N+I,5)=P(IPA(I),5)
  330   CONTINUE
      ELSE
        IF(MSTJ(43).LE.2) PEM=V(IM,2)
        IF(MSTJ(43).GE.3) PEM=P(IM,4)
        P(N+1,5)=MIN(P(IM,5),V(IM,1)*PEM)
        P(N+2,5)=MIN(P(IM,5),(1D0-V(IM,1))*PEM)
        IF(K(N+2,2).EQ.22) P(N+2,5)=PMTH(1,22)
      ENDIF
      DO 340 I=1,NEP
        PMSD(I)=P(N+I,5)
        IF(ISI(I).EQ.1) THEN
          IR=IREF(N+I-NS)
          IF(P(N+I,5).LE.PMTH(3,IR)) P(N+I,5)=PMTH(1,IR)
        ENDIF
        V(N+I,5)=P(N+I,5)**2
  340 CONTINUE
 
C...Choose one of the daughters for evolution.
  350 INUM=0
      IF(NEP.EQ.1) INUM=1
      DO 360 I=1,NEP
        IF(INUM.EQ.0.AND.ISL(I).EQ.1) INUM=I
  360 CONTINUE
      DO 370 I=1,NEP
        IF(INUM.EQ.0.AND.ITRY(I).EQ.0.AND.ISI(I).EQ.1) THEN
          IR=IREF(N+I-NS)
          IF(P(N+I,5).GE.PMTH(2,IR)) INUM=I
        ENDIF
  370 CONTINUE
      IF(INUM.EQ.0) THEN
        RMAX=0D0
        DO 380 I=1,NEP
          IF(ISI(I).EQ.1.AND.PMSD(I).GE.PMQT2E) THEN
            RPM=P(N+I,5)/PMSD(I)
            IR=IREF(N+I-NS)
            IF(RPM.GT.RMAX.AND.P(N+I,5).GE.PMTH(2,IR)) THEN
              RMAX=RPM
              INUM=I
            ENDIF
          ENDIF
  380   CONTINUE
      ENDIF
 
C...Cancel choice of predetermined daughter already treated.
      INUM=MAX(1,INUM)
      INUMT=INUM
      IF(MPSPD.EQ.1.AND.IGM.EQ.0.AND.ITRY(INUMT).GE.1) THEN
        IF(K(IP1-1+INUM,4).GT.0) INUM=3-INUM
      ELSEIF(MPSPD.EQ.1.AND.IM.EQ.NS+2.AND.ITRY(INUMT).GE.1) THEN
        IF(KFLD(INUMT).NE.21.AND.K(IP1+2,4).GT.0) INUM=3-INUM
        IF(KFLD(INUMT).EQ.21.AND.K(IP1+3,4).GT.0) INUM=3-INUM
      ENDIF
 
C...Store information on choice of evolving daughter.
      IEP(1)=N+INUM
      DO 390 I=2,NEP
        IEP(I)=IEP(I-1)+1
        IF(IEP(I).GT.N+NEP) IEP(I)=N+1
  390 CONTINUE
      DO 400 I=1,NEP
        KFL(I)=IABS(K(IEP(I),2))
  400 CONTINUE
      ITRY(INUM)=ITRY(INUM)+1
      IF(ITRY(INUM).GT.200) THEN
        CALL PYERRM(14,'(PYSHOW:) caught in infinite loop')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      Z=0.5D0
      IR=IREF(IEP(1)-NS)
      IF(KSH(IR).EQ.0) GOTO 450
      IF(P(IEP(1),5).LT.PMTH(2,IR)) GOTO 450
 
C...Check if evolution already predetermined for daughter.
      IPSPD=0
      IF(MPSPD.EQ.1.AND.IGM.EQ.0) THEN
        IF(K(IP1-1+INUM,4).GT.0) IPSPD=IP1-1+INUM
      ELSEIF(MPSPD.EQ.1.AND.IM.EQ.NS+2) THEN
        IF(KFL(1).NE.21.AND.K(IP1+2,4).GT.0) IPSPD=IP1+2
        IF(KFL(1).EQ.21.AND.K(IP1+3,4).GT.0) IPSPD=IP1+3
      ENDIF
      IF(INUM.EQ.1.OR.INUM.EQ.2) THEN
        ISSET(INUM)=0
        IF(IPSPD.NE.0) ISSET(INUM)=1
      ENDIF
 
C...Select side for interference with initial state partons.
      IF(MIIS.GE.1.AND.IEP(1).LE.NS+3) THEN
        III=IEP(1)-NS-1
        ISII(III)=0
        IF(IABS(KCII(III)).EQ.1.AND.NIIS(III).EQ.1) THEN
          ISII(III)=1
        ELSEIF(KCII(III).EQ.2.AND.NIIS(III).EQ.1) THEN
          IF(PYR(0).GT.0.5D0) ISII(III)=1
        ELSEIF(KCII(III).EQ.2.AND.NIIS(III).EQ.2) THEN
          ISII(III)=1
          IF(PYR(0).GT.0.5D0) ISII(III)=2
        ENDIF
      ENDIF
 
C...Calculate allowed z range.
      IF(NEP.EQ.1) THEN
        PMED=PS(4)
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN
        PMED=P(IM,5)
      ELSE
        IF(INUM.EQ.1) PMED=V(IM,1)*PEM
        IF(INUM.EQ.2) PMED=(1D0-V(IM,1))*PEM
      ENDIF
      IF(MOD(MSTJ(43),2).EQ.1) THEN
        ZC=PMTH(2,21)/PMED
        ZCE=PMTH(2,22)/PMED
        IF(ISCOL(IR).EQ.0) ZCE=0.5D0*PARJ(90)/PMED
      ELSE
        ZC=0.5D0*(1D0-SQRT(MAX(0D0,1D0-(2D0*PMTH(2,21)/PMED)**2)))
        IF(ZC.LT.1D-6) ZC=(PMTH(2,21)/PMED)**2
        PMTMPE=PMTH(2,22)
        IF(ISCOL(IR).EQ.0) PMTMPE=0.5D0*PARJ(90)
        ZCE=0.5D0*(1D0-SQRT(MAX(0D0,1D0-(2D0*PMTMPE/PMED)**2)))
        IF(ZCE.LT.1D-6) ZCE=(PMTMPE/PMED)**2
      ENDIF
      ZC=MIN(ZC,0.491D0)
      ZCE=MIN(ZCE,0.49991D0)
      IF(((MSTJ(41).EQ.1.AND.ZC.GT.0.49D0).OR.(MSTJ(41).GE.2.AND.
     &MIN(ZC,ZCE).GT.0.4999D0)).AND.IPSPD.EQ.0) THEN
        P(IEP(1),5)=PMTH(1,IR)
        V(IEP(1),5)=P(IEP(1),5)**2
        GOTO 450
      ENDIF
 
C...Integral of Altarelli-Parisi z kernel for QCD.
C...(Includes squark and gluino; with factor N_C/C_F extra for latter).
      IF(MSTJ(49).EQ.0.AND.KFL(1).EQ.21) THEN
        FBR=6D0*LOG((1D0-ZC)/ZC)+MSTJ(45)*0.5D0
C...QUARKONIA+++
C...Evolution of QQ~[3S18] state if MSTP(148)=1.
      ELSEIF(MSTJ(49).EQ.0.AND.MSTP(149).GE.0.AND.
     &       (KFL(1).EQ.9900443.OR.KFL(1).EQ.9900553)) THEN
        FBR=6D0*LOG((1D0-ZC)/ZC)
C...QUARKONIA---
      ELSEIF(MSTJ(49).EQ.0) THEN
        FBR=(8D0/3D0)*LOG((1D0-ZC)/ZC)
        IF(IGLUI.EQ.1.AND.IR.GE.31) FBR=FBR*(9D0/4D0)
 
C...Integral of Altarelli-Parisi z kernel for scalar gluon.
      ELSEIF(MSTJ(49).EQ.1.AND.KFL(1).EQ.21) THEN
        FBR=(PARJ(87)+MSTJ(45)*PARJ(88))*(1D0-2D0*ZC)
      ELSEIF(MSTJ(49).EQ.1) THEN
        FBR=(1D0-2D0*ZC)/3D0
        IF(IGM.EQ.0.AND.M3JC.GE.1) FBR=4D0*FBR
 
C...Integral of Altarelli-Parisi z kernel for Abelian vector gluon.
      ELSEIF(KFL(1).EQ.21) THEN
        FBR=6D0*MSTJ(45)*(0.5D0-ZC)
      ELSE
        FBR=2D0*LOG((1D0-ZC)/ZC)
      ENDIF
 
C...Reset QCD probability for colourless.
      IF(ISCOL(IR).EQ.0) FBR=0D0
 
C...Integral of Altarelli-Parisi kernel for photon emission.
      FBRE=0D0
      IF(MSTJ(41).GE.2.AND.ISCHG(IR).EQ.1) THEN
        IF(KFL(1).LE.18) THEN
          FBRE=(KCHG(KFL(1),1)/3D0)**2*2D0*LOG((1D0-ZCE)/ZCE)
        ENDIF
        IF(MSTJ(41).EQ.10) FBRE=PARJ(84)*FBRE
      ENDIF
 
C...Inner veto algorithm starts. Find maximum mass for evolution.
  410 PMS=V(IEP(1),5)
      IF(IGM.GE.0) THEN
        PM2=0D0
        DO 420 I=2,NEP
          PM=P(IEP(I),5)
          IRI=IREF(IEP(I)-NS)
          IF(KSH(IRI).EQ.1) PM=PMTH(2,IRI)
          PM2=PM2+PM
  420   CONTINUE
        PMS=MIN(PMS,(P(IM,5)-PM2)**2)
      ENDIF
 
C...Select mass for daughter in QCD evolution.
      B0=27D0/6D0
      DO 430 IFF=4,MSTJ(45)
        IF(PMS.GT.4D0*PMTH(2,IFF)**2) B0=(33D0-2D0*IFF)/6D0
  430 CONTINUE
C...Shift m^2 for evolution in Q^2 = m^2 - m(onshell)^2.
      PMSC=MAX(0.5D0*PARJ(82),PMS-PMTH(1,IR)**2)
C...Already predetermined choice.
      IF(IPSPD.NE.0) THEN
        PMSQCD=P(IPSPD,5)**2
      ELSEIF(FBR.LT.1D-3) THEN
        PMSQCD=0D0
      ELSEIF(MSTJ(44).LE.0) THEN
        PMSQCD=PMSC*EXP(MAX(-50D0,LOG(PYR(0))*PARU(2)/(PARU(111)*FBR)))
      ELSEIF(MSTJ(44).EQ.1) THEN
        PMSQCD=4D0*ALAMS*(0.25D0*PMSC/ALAMS)**(PYR(0)**(B0/FBR))
      ELSE
        PMSQCD=PMSC*EXP(MAX(-50D0,ALFM*B0*LOG(PYR(0))/FBR))
      ENDIF
C...Shift back m^2 from evolution in Q^2 = m^2 - m(onshell)^2.
      IF(IPSPD.EQ.0) PMSQCD=PMSQCD+PMTH(1,IR)**2
      IF(ZC.GT.0.49D0.OR.PMSQCD.LE.PMTH(4,IR)**2) PMSQCD=PMTH(2,IR)**2
      V(IEP(1),5)=PMSQCD
      MCE=1
 
C...Select mass for daughter in QED evolution.
      IF(MSTJ(41).GE.2.AND.ISCHG(IR).EQ.1.AND.IPSPD.EQ.0) THEN
C...Shift m^2 for evolution in Q^2 = m^2 - m(onshell)^2.
        PMSE=MAX(0.5D0*PARJ(83),PMS-PMTH(1,IR)**2)
        IF(FBRE.LT.1D-3) THEN
          PMSQED=0D0
        ELSE
          PMSQED=PMSE*EXP(MAX(-50D0,LOG(PYR(0))*PARU(2)/
     &    (PARU(101)*FBRE)))
        ENDIF
C...Shift back m^2 from evolution in Q^2 = m^2 - m(onshell)^2.
        PMSQED=PMSQED+PMTH(1,IR)**2
        IF(ZCE.GT.0.4999D0.OR.PMSQED.LE.PMTH(5,IR)**2) PMSQED=
     &  PMTH(2,IR)**2
        IF(PMSQED.GT.PMSQCD) THEN
          V(IEP(1),5)=PMSQED
          MCE=2
        ENDIF
      ENDIF
 
C...Check whether daughter mass below cutoff.
      P(IEP(1),5)=SQRT(V(IEP(1),5))
      IF(P(IEP(1),5).LE.PMTH(3,IR)) THEN
        P(IEP(1),5)=PMTH(1,IR)
        V(IEP(1),5)=P(IEP(1),5)**2
        GOTO 450
      ENDIF
 
C...Already predetermined choice of z, and flavour in g -> qqbar.
      IF(IPSPD.NE.0) THEN
        IPSGD1=K(IPSPD,4)
        IPSGD2=K(IPSPD,5)
        PMSGD1=P(IPSGD1,5)**2
        PMSGD2=P(IPSGD2,5)**2
        ALAMPS=SQRT(MAX(1D-10,(PMSQCD-PMSGD1-PMSGD2)**2-
     &  4D0*PMSGD1*PMSGD2))
        Z=0.5D0*(PMSQCD*(2D0*P(IPSGD1,4)/P(IPSPD,4)-1D0)+ALAMPS-
     &  PMSGD1+PMSGD2)/ALAMPS
        Z=MAX(0.00001D0,MIN(0.99999D0,Z))
        IF(KFL(1).NE.21) THEN
          K(IEP(1),5)=21
        ELSE
          K(IEP(1),5)=IABS(K(IPSGD1,2))
        ENDIF
 
C...Select z value of branching: q -> qgamma.
      ELSEIF(MCE.EQ.2) THEN
        Z=1D0-(1D0-ZCE)*(ZCE/(1D0-ZCE))**PYR(0)
        IF(1D0+Z**2.LT.2D0*PYR(0)) GOTO 410
        K(IEP(1),5)=22
 
C...QUARKONIA+++
C...Select z value of branching: QQ~[3S18] -> QQ~[3S18]g.
      ELSEIF(MSTJ(49).EQ.0.AND.
     &       (KFL(1).EQ.9900443.OR.KFL(1).EQ.9900553)) THEN
        Z=(1D0-ZC)*(ZC/(1D0-ZC))**PYR(0)
C...Select always the harder 'gluon' if the switch MSTP(149)<=0.
        IF(MSTP(149).LE.0.OR.PYR(0).GT.0.5D0) Z=1D0-Z
        IF((1D0-Z*(1D0-Z))**2.LT.PYR(0)) GOTO 410
        K(IEP(1),5)=21
C...QUARKONIA---
 
C...Select z value of branching: q -> qg, g -> gg, g -> qqbar.
      ELSEIF(MSTJ(49).NE.1.AND.KFL(1).NE.21) THEN
        Z=1D0-(1D0-ZC)*(ZC/(1D0-ZC))**PYR(0)
C...Only do z weighting when no ME correction afterwards.
        IF(M3JC.EQ.0.AND.1D0+Z**2.LT.2D0*PYR(0)) GOTO 410
        K(IEP(1),5)=21
      ELSEIF(MSTJ(49).EQ.0.AND.MSTJ(45)*0.5D0.LT.PYR(0)*FBR) THEN
        Z=(1D0-ZC)*(ZC/(1D0-ZC))**PYR(0)
        IF(PYR(0).GT.0.5D0) Z=1D0-Z
        IF((1D0-Z*(1D0-Z))**2.LT.PYR(0)) GOTO 410
        K(IEP(1),5)=21
      ELSEIF(MSTJ(49).NE.1) THEN
        Z=PYR(0)
        IF(Z**2+(1D0-Z)**2.LT.PYR(0)) GOTO 410
        KFLB=1+INT(MSTJ(45)*PYR(0))
        PMQ=4D0*PMTH(2,KFLB)**2/V(IEP(1),5)
        IF(PMQ.GE.1D0) GOTO 410
        IF(MSTJ(44).LE.2.OR.MSTJ(44).EQ.4) THEN
          IF(Z.LT.ZC.OR.Z.GT.1D0-ZC) GOTO 410
          PMQ0=4D0*PMTH(2,21)**2/V(IEP(1),5)
          IF(MOD(MSTJ(43),2).EQ.0.AND.(1D0+0.5D0*PMQ)*SQRT(1D0-PMQ)
     &    .LT.PYR(0)*(1D0+0.5D0*PMQ0)*SQRT(1D0-PMQ0)) GOTO 410
        ELSE
          IF((1D0+0.5D0*PMQ)*SQRT(1D0-PMQ).LT.PYR(0)) GOTO 410
        ENDIF
        K(IEP(1),5)=KFLB
 
C...Ditto for scalar gluon model.
      ELSEIF(KFL(1).NE.21) THEN
        Z=1D0-SQRT(ZC**2+PYR(0)*(1D0-2D0*ZC))
        K(IEP(1),5)=21
      ELSEIF(PYR(0)*(PARJ(87)+MSTJ(45)*PARJ(88)).LE.PARJ(87)) THEN
        Z=ZC+(1D0-2D0*ZC)*PYR(0)
        K(IEP(1),5)=21
      ELSE
        Z=ZC+(1D0-2D0*ZC)*PYR(0)
        KFLB=1+INT(MSTJ(45)*PYR(0))
        PMQ=4D0*PMTH(2,KFLB)**2/V(IEP(1),5)
        IF(PMQ.GE.1D0) GOTO 410
        K(IEP(1),5)=KFLB
      ENDIF
 
C...Correct to alpha_s(pT^2) (optionally m^2/4 for g -> q qbar).
      IF(MCE.EQ.1.AND.MSTJ(44).GE.2.AND.IPSPD.EQ.0) THEN
        IF(KFL(1).EQ.21.AND.K(IEP(1),5).LT.10.AND.
     &  (MSTJ(44).EQ.3.OR.MSTJ(44).EQ.5)) THEN
          IF(ALFM/LOG(V(IEP(1),5)*0.25D0/ALAMS).LT.PYR(0)) GOTO 410
        ELSE
          PT2APP=Z*(1D0-Z)*V(IEP(1),5)
          IF(MSTJ(44).GE.4) PT2APP=PT2APP*
     &    (1D0-PMTH(1,IR)**2/V(IEP(1),5))**2
          IF(PT2APP.LT.PT2MIN) GOTO 410
          IF(ALFM/LOG(PT2APP/ALAMS).LT.PYR(0)) GOTO 410
        ENDIF
      ENDIF
 
C...Check if z consistent with chosen m.
      IF(KFL(1).EQ.21) THEN
        IRGD1=IABS(K(IEP(1),5))
        IRGD2=IRGD1
      ELSE
        IRGD1=IR
        IRGD2=IABS(K(IEP(1),5))
      ENDIF
      IF(NEP.EQ.1) THEN
        PED=PS(4)
      ELSEIF(NEP.GE.3) THEN
        PED=P(IEP(1),4)
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN
        PED=0.5D0*(V(IM,5)+V(IEP(1),5)-PM2**2)/P(IM,5)
      ELSE
        IF(IEP(1).EQ.N+1) PED=V(IM,1)*PEM
        IF(IEP(1).EQ.N+2) PED=(1D0-V(IM,1))*PEM
      ENDIF
      IF(MOD(MSTJ(43),2).EQ.1) THEN
        PMQTH3=0.5D0*PARJ(82)
        IF(IRGD2.EQ.22) PMQTH3=0.5D0*PARJ(83)
        IF(IRGD2.EQ.22.AND.ISCOL(IR).EQ.0) PMQTH3=0.5D0*PARJ(90)
        PMQ1=(PMTH(1,IRGD1)**2+PMQTH3**2)/V(IEP(1),5)
        PMQ2=(PMTH(1,IRGD2)**2+PMQTH3**2)/V(IEP(1),5)
        ZD=SQRT(MAX(0D0,(1D0-V(IEP(1),5)/PED**2)*((1D0-PMQ1-PMQ2)**2-
     &  4D0*PMQ1*PMQ2)))
        ZH=1D0+PMQ1-PMQ2
      ELSE
        ZD=SQRT(MAX(0D0,1D0-V(IEP(1),5)/PED**2))
        ZH=1D0
      ENDIF
      IF(KFL(1).EQ.21.AND.K(IEP(1),5).LT.10.AND.
     &(MSTJ(44).EQ.3.OR.MSTJ(44).EQ.5)) THEN
      ELSEIF(IPSPD.NE.0) THEN
      ELSE
        ZL=0.5D0*(ZH-ZD)
        ZU=0.5D0*(ZH+ZD)
        IF(Z.LT.ZL.OR.Z.GT.ZU) GOTO 410
      ENDIF
      IF(KFL(1).EQ.21) V(IEP(1),3)=LOG(ZU*(1D0-ZL)/MAX(1D-20,ZL*
     &(1D0-ZU)))
      IF(KFL(1).NE.21) V(IEP(1),3)=LOG((1D0-ZL)/MAX(1D-10,1D0-ZU))
 
C...Width suppression for q -> q + g.
      IF(MSTJ(40).NE.0.AND.KFL(1).NE.21.AND.IPSPD.EQ.0) THEN
        IF(IGM.EQ.0) THEN
          EGLU=0.5D0*PS(5)*(1D0-Z)*(1D0+V(IEP(1),5)/V(NS+1,5))
        ELSE
          EGLU=PMED*(1D0-Z)
        ENDIF
        CHI=PARJ(89)**2/(PARJ(89)**2+EGLU**2)
        IF(MSTJ(40).EQ.1) THEN
          IF(CHI.LT.PYR(0)) GOTO 410
        ELSEIF(MSTJ(40).EQ.2) THEN
          IF(1D0-CHI.LT.PYR(0)) GOTO 410
        ENDIF
      ENDIF
 
C...Three-jet matrix element correction.
      IF(M3JC.GE.1) THEN
        WME=1D0
        WSHOW=1D0
 
C...QED matrix elements: only for massless case so far.
        IF(MCE.EQ.2.AND.IGM.EQ.0) THEN
          X1=Z*(1D0+V(IEP(1),5)/V(NS+1,5))
          X2=1D0-V(IEP(1),5)/V(NS+1,5)
          X3=(1D0-X1)+(1D0-X2)
          KI1=K(IPA(INUM),2)
          KI2=K(IPA(3-INUM),2)
          QF1=KCHG(PYCOMP(KI1),1)*ISIGN(1,KI1)/3D0
          QF2=KCHG(PYCOMP(KI2),1)*ISIGN(1,KI2)/3D0
          WSHOW=QF1**2*(1D0-X1)/X3*(1D0+(X1/(2D0-X2))**2)+
     &    QF2**2*(1D0-X2)/X3*(1D0+(X2/(2D0-X1))**2)
          WME=(QF1*(1D0-X1)/X3-QF2*(1D0-X2)/X3)**2*(X1**2+X2**2)
        ELSEIF(MCE.EQ.2) THEN
 
C...QCD matrix elements, including mass effects.
        ELSEIF(MSTJ(49).NE.1.AND.K(IEP(1),2).NE.21) THEN
          PS1ME=V(IEP(1),5)
          PM1ME=PMTH(1,IR)
          M3JCC=M3JC
          IF(IR.GE.31.AND.IGM.EQ.0) THEN
C...QCD ME: original parton, first branching.
            PM2ME=PMTH(1,63-IR)
            ECMME=PS(5)
          ELSEIF(IR.GE.31) THEN
C...QCD ME: original parton, subsequent branchings.
            PM2ME=PMTH(1,63-IR)
            PEDME=PEM*(V(IM,1)+(1D0-V(IM,1))*PS1ME/V(IM,5))
            ECMME=PEDME+SQRT(MAX(0D0,PEDME**2-PS1ME+PM2ME**2))
          ELSEIF(K(IM,2).EQ.21) THEN
C...QCD ME: secondary partons, first branching.
            PM2ME=PM1ME
            ZMME=V(IM,1)
            IF(IEP(1).GT.IEP(2)) ZMME=1D0-ZMME
            PMLME=SQRT(MAX(0D0,(V(IM,5)-PS1ME-PM2ME**2)**2-
     &      4D0*PS1ME*PM2ME**2))
            PEDME=PEM*(0.5D0*(V(IM,5)-PMLME+PS1ME-PM2ME**2)+PMLME*ZMME)/
     &      V(IM,5)
            ECMME=PEDME+SQRT(MAX(0D0,PEDME**2-PS1ME+PM2ME**2))
            M3JCC=66
          ELSE
C...QCD ME: secondary partons, subsequent branchings.
            PM2ME=PM1ME
            PEDME=PEM*(V(IM,1)+(1D0-V(IM,1))*PS1ME/V(IM,5))
            ECMME=PEDME+SQRT(MAX(0D0,PEDME**2-PS1ME+PM2ME**2))
            M3JCC=66
          ENDIF
C...Construct ME variables.
          R1ME=PM1ME/ECMME
          R2ME=PM2ME/ECMME
          X1=(1D0+PS1ME/ECMME**2-R2ME**2)*(Z+(1D0-Z)*PM1ME**2/PS1ME)
          X2=1D0+R2ME**2-PS1ME/ECMME**2
C...Call ME, with right order important for two inequivalent showerers.
          IF(IR.EQ.IORD+30) THEN
            WME=PYMAEL(M3JCC,X1,X2,R1ME,R2ME,ALPHA)
          ELSE
            WME=PYMAEL(M3JCC,X2,X1,R2ME,R1ME,ALPHA)
          ENDIF
C...Split up total ME when two radiating partons.
          ISPRAD=1
          IF((M3JCC.GE.16.AND.M3JCC.LE.19).OR.
     &    (M3JCC.GE.26.AND.M3JCC.LE.29).OR.
     &    (M3JCC.GE.36.AND.M3JCC.LE.39).OR.
     &    (M3JCC.GE.46.AND.M3JCC.LE.49).OR.
     &    (M3JCC.GE.56.AND.M3JCC.LE.64)) ISPRAD=0
          IF(ISPRAD.EQ.1) WME=WME*MAX(1D-10,1D0+R1ME**2-R2ME**2-X1)/
     &    MAX(1D-10,2D0-X1-X2)
C...Evaluate shower rate to be compared with.
          WSHOW=2D0/(MAX(1D-10,2D0-X1-X2)*
     &    MAX(1D-10,1D0+R2ME**2-R1ME**2-X2))
          IF(IGLUI.EQ.1.AND.IR.GE.31) WSHOW=(9D0/4D0)*WSHOW
        ELSEIF(MSTJ(49).NE.1) THEN
 
C...Toy model scalar theory matrix elements; no mass effects.
        ELSE
          X1=Z*(1D0+V(IEP(1),5)/V(NS+1,5))
          X2=1D0-V(IEP(1),5)/V(NS+1,5)
          X3=(1D0-X1)+(1D0-X2)
          WSHOW=4D0*X3*((1D0-X1)/(2D0-X2)**2+(1D0-X2)/(2D0-X1)**2)
          WME=X3**2
          IF(MSTJ(102).GE.2) WME=X3**2-2D0*(1D0+X3)*(1D0-X1)*(1D0-X2)*
     &    PARJ(171)
        ENDIF
 
        IF(WME.LT.PYR(0)*WSHOW) GOTO 410
      ENDIF
 
C...Impose angular ordering by rejection of nonordered emission.
      IF(MCE.EQ.1.AND.IGM.GT.0.AND.MSTJ(42).GE.2.AND.IPSPD.EQ.0) THEN
        PEMAO=V(IM,1)*P(IM,4)
        IF(IEP(1).EQ.N+2) PEMAO=(1D0-V(IM,1))*P(IM,4)
        IF(IR.GE.31.AND.MSTJ(42).GE.5) THEN
          MAOD=0
        ELSEIF(KFL(1).EQ.21.AND.K(IEP(1),5).LE.10.AND.(MSTJ(42).EQ.4
     &  .OR.MSTJ(42).EQ.7)) THEN
          MAOD=0
        ELSEIF(KFL(1).EQ.21.AND.K(IEP(1),5).LE.10.AND.(MSTJ(42).EQ.3
     &  .OR.MSTJ(42).EQ.6)) THEN
          MAOD=1
          PMDAO=PMTH(2,K(IEP(1),5))
          THE2ID=Z*(1D0-Z)*PEMAO**2/(V(IEP(1),5)-4D0*PMDAO**2)
        ELSE
          MAOD=1
          THE2ID=Z*(1D0-Z)*PEMAO**2/V(IEP(1),5)
          IF(MSTJ(42).GE.3.AND.MSTJ(42).NE.5) THE2ID=THE2ID*
     &    (1D0+PMTH(1,IR)**2*(1D0-Z)/(V(IEP(1),5)*Z))**2
        ENDIF
        MAOM=1
        IAOM=IM
  440   IF(K(IAOM,5).EQ.22) THEN
          IAOM=K(IAOM,3)
          IF(K(IAOM,3).LE.NS) MAOM=0
          IF(MAOM.EQ.1) GOTO 440
        ENDIF
        IF(MAOM.EQ.1.AND.MAOD.EQ.1) THEN
          THE2IM=V(IAOM,1)*(1D0-V(IAOM,1))*P(IAOM,4)**2/V(IAOM,5)
          IF(THE2ID.LT.THE2IM) GOTO 410
        ENDIF
      ENDIF
 
C...Impose user-defined maximum angle at first branching.
      IF(MSTJ(48).EQ.1.AND.IPSPD.EQ.0) THEN
        IF(NEP.EQ.1.AND.IM.EQ.NS) THEN
          THE2ID=Z*(1D0-Z)*PS(4)**2/V(IEP(1),5)
          IF(PARJ(85)**2*THE2ID.LT.1D0) GOTO 410
        ELSEIF(NEP.EQ.2.AND.IEP(1).EQ.NS+2) THEN
          THE2ID=Z*(1D0-Z)*(0.5D0*P(IM,4))**2/V(IEP(1),5)
          IF(PARJ(85)**2*THE2ID.LT.1D0) GOTO 410
        ELSEIF(NEP.EQ.2.AND.IEP(1).EQ.NS+3) THEN
          THE2ID=Z*(1D0-Z)*(0.5D0*P(IM,4))**2/V(IEP(1),5)
          IF(PARJ(86)**2*THE2ID.LT.1D0) GOTO 410
        ENDIF
      ENDIF
 
C...Impose angular constraint in first branching from interference
C...with initial state partons.
      IF(MIIS.GE.2.AND.IEP(1).LE.NS+3) THEN
        THE2D=MAX((1D0-Z)/Z,Z/(1D0-Z))*V(IEP(1),5)/(0.5D0*P(IM,4))**2
        IF(IEP(1).EQ.NS+2.AND.ISII(1).GE.1) THEN
          IF(THE2D.GT.THEIIS(1,ISII(1))**2) GOTO 410
        ELSEIF(IEP(1).EQ.NS+3.AND.ISII(2).GE.1) THEN
          IF(THE2D.GT.THEIIS(2,ISII(2))**2) GOTO 410
        ENDIF
      ENDIF
 
C...End of inner veto algorithm. Check if only one leg evolved so far.
  450 V(IEP(1),1)=Z
      ISL(1)=0
      ISL(2)=0
      IF(NEP.EQ.1) GOTO 490
      IF(NEP.EQ.2.AND.P(IEP(1),5)+P(IEP(2),5).GE.P(IM,5)) GOTO 350
      DO 460 I=1,NEP
        IR=IREF(N+I-NS)
        IF(ITRY(I).EQ.0.AND.KSH(IR).EQ.1) THEN
          IF(P(N+I,5).GE.PMTH(2,IR)) GOTO 350
        ENDIF
  460 CONTINUE
 
C...Check if chosen multiplet m1,m2,z1,z2 is physical.
      IF(NEP.GE.3) THEN
        PMSUM=0D0
        DO 470 I=1,NEP
          PMSUM=PMSUM+P(N+I,5)
  470   CONTINUE
        IF(PMSUM.GE.PS(5)) GOTO 350
      ELSEIF(IGM.EQ.0.OR.MSTJ(43).LE.2.OR.MOD(MSTJ(43),2).EQ.0) THEN
        DO 480 I1=N+1,N+2
          IRDA=IREF(I1-NS)
          IF(KSH(IRDA).EQ.0) GOTO 480
          IF(P(I1,5).LT.PMTH(2,IRDA)) GOTO 480
          IF(IRDA.EQ.21) THEN
            IRGD1=IABS(K(I1,5))
            IRGD2=IRGD1
          ELSE
            IRGD1=IRDA
            IRGD2=IABS(K(I1,5))
          ENDIF
          I2=2*N+3-I1
          IF(IGM.EQ.0.OR.MSTJ(43).LE.2) THEN
            PED=0.5D0*(V(IM,5)+V(I1,5)-V(I2,5))/P(IM,5)
          ELSE
            IF(I1.EQ.N+1) ZM=V(IM,1)
            IF(I1.EQ.N+2) ZM=1D0-V(IM,1)
            PML=SQRT((V(IM,5)-V(N+1,5)-V(N+2,5))**2-
     &      4D0*V(N+1,5)*V(N+2,5))
            PED=PEM*(0.5D0*(V(IM,5)-PML+V(I1,5)-V(I2,5))+PML*ZM)/
     &      V(IM,5)
          ENDIF
          IF(MOD(MSTJ(43),2).EQ.1) THEN
            PMQTH3=0.5D0*PARJ(82)
            IF(IRGD2.EQ.22) PMQTH3=0.5D0*PARJ(83)
            IF(IRGD2.EQ.22.AND.ISCOL(IRDA).EQ.0) PMQTH3=0.5D0*PARJ(90)
            PMQ1=(PMTH(1,IRGD1)**2+PMQTH3**2)/V(I1,5)
            PMQ2=(PMTH(1,IRGD2)**2+PMQTH3**2)/V(I1,5)
            ZD=SQRT(MAX(0D0,(1D0-V(I1,5)/PED**2)*((1D0-PMQ1-PMQ2)**2-
     &      4D0*PMQ1*PMQ2)))
            ZH=1D0+PMQ1-PMQ2
          ELSE
            ZD=SQRT(MAX(0D0,1D0-V(I1,5)/PED**2))
            ZH=1D0
          ENDIF
          IF(IRDA.EQ.21.AND.IRGD1.LT.10.AND.
     &    (MSTJ(44).EQ.3.OR.MSTJ(44).EQ.5)) THEN
          ELSE
            ZL=0.5D0*(ZH-ZD)
            ZU=0.5D0*(ZH+ZD)
            IF(I1.EQ.N+1.AND.(V(I1,1).LT.ZL.OR.V(I1,1).GT.ZU).AND.
     &      ISSET(1).EQ.0) THEN
              ISL(1)=1
            ELSEIF(I1.EQ.N+2.AND.(V(I1,1).LT.ZL.OR.V(I1,1).GT.ZU).AND.
     &      ISSET(2).EQ.0) THEN
              ISL(2)=1
            ENDIF
          ENDIF
          IF(IRDA.EQ.21) V(I1,4)=LOG(ZU*(1D0-ZL)/MAX(1D-20,
     &    ZL*(1D0-ZU)))
          IF(IRDA.NE.21) V(I1,4)=LOG((1D0-ZL)/MAX(1D-10,1D0-ZU))
  480   CONTINUE
        IF(ISL(1).EQ.1.AND.ISL(2).EQ.1.AND.ISLM.NE.0) THEN
          ISL(3-ISLM)=0
          ISLM=3-ISLM
        ELSEIF(ISL(1).EQ.1.AND.ISL(2).EQ.1) THEN
          ZDR1=MAX(0D0,V(N+1,3)/MAX(1D-6,V(N+1,4))-1D0)
          ZDR2=MAX(0D0,V(N+2,3)/MAX(1D-6,V(N+2,4))-1D0)
          IF(ZDR2.GT.PYR(0)*(ZDR1+ZDR2)) ISL(1)=0
          IF(ISL(1).EQ.1) ISL(2)=0
          IF(ISL(1).EQ.0) ISLM=1
          IF(ISL(2).EQ.0) ISLM=2
        ENDIF
        IF(ISL(1).EQ.1.OR.ISL(2).EQ.1) GOTO 350
      ENDIF
      IRD1=IREF(N+1-NS)
      IRD2=IREF(N+2-NS)
      IF(IGM.GT.0) THEN
        IF(MOD(MSTJ(43),2).EQ.1.AND.(P(N+1,5).GE.
     &  PMTH(2,IRD1).OR.P(N+2,5).GE.PMTH(2,IRD2))) THEN
          PMQ1=V(N+1,5)/V(IM,5)
          PMQ2=V(N+2,5)/V(IM,5)
          ZD=SQRT(MAX(0D0,(1D0-V(IM,5)/PEM**2)*((1D0-PMQ1-PMQ2)**2-
     &    4D0*PMQ1*PMQ2)))
          ZH=1D0+PMQ1-PMQ2
          ZL=0.5D0*(ZH-ZD)
          ZU=0.5D0*(ZH+ZD)
          IF(V(IM,1).LT.ZL.OR.V(IM,1).GT.ZU) GOTO 350
        ENDIF
      ENDIF
 
C...Accepted branch. Construct four-momentum for initial partons.
  490 MAZIP=0
      MAZIC=0
      IF(NEP.EQ.1) THEN
        P(N+1,1)=0D0
        P(N+1,2)=0D0
        P(N+1,3)=SQRT(MAX(0D0,(P(IPA(1),4)+P(N+1,5))*(P(IPA(1),4)-
     &  P(N+1,5))))
        P(N+1,4)=P(IPA(1),4)
        V(N+1,2)=P(N+1,4)
      ELSEIF(IGM.EQ.0.AND.NEP.EQ.2) THEN
        PED1=0.5D0*(V(IM,5)+V(N+1,5)-V(N+2,5))/P(IM,5)
        P(N+1,1)=0D0
        P(N+1,2)=0D0
        P(N+1,3)=SQRT(MAX(0D0,(PED1+P(N+1,5))*(PED1-P(N+1,5))))
        P(N+1,4)=PED1
        P(N+2,1)=0D0
        P(N+2,2)=0D0
        P(N+2,3)=-P(N+1,3)
        P(N+2,4)=P(IM,5)-PED1
        V(N+1,2)=P(N+1,4)
        V(N+2,2)=P(N+2,4)
      ELSEIF(NEP.GE.3) THEN
C...Rescale all momenta for energy conservation.
        LOOP=0
        PES=0D0
        PQS=0D0
        DO 510 I=1,NEP
          DO 500 J=1,4
            P(N+I,J)=P(IPA(I),J)
  500     CONTINUE
          PES=PES+P(N+I,4)
          PQS=PQS+P(N+I,5)**2/P(N+I,4)
  510   CONTINUE
  520   LOOP=LOOP+1
        FAC=(PS(5)-PQS)/(PES-PQS)
        PES=0D0
        PQS=0D0
        DO 540 I=1,NEP
          DO 530 J=1,3
            P(N+I,J)=FAC*P(N+I,J)
  530     CONTINUE
          P(N+I,4)=SQRT(P(N+I,5)**2+P(N+I,1)**2+P(N+I,2)**2+P(N+I,3)**2)
          V(N+I,2)=P(N+I,4)
          PES=PES+P(N+I,4)
          PQS=PQS+P(N+I,5)**2/P(N+I,4)
  540   CONTINUE
        IF(LOOP.LT.10.AND.ABS(PES-PS(5)).GT.1D-12*PS(5)) GOTO 520
 
C...Construct transverse momentum for ordinary branching in shower.
      ELSE
        ZM=V(IM,1)
        LOOPPT=0
  550   LOOPPT=LOOPPT+1
        PZM=SQRT(MAX(0D0,(PEM+P(IM,5))*(PEM-P(IM,5))))
        PMLS=(V(IM,5)-V(N+1,5)-V(N+2,5))**2-4D0*V(N+1,5)*V(N+2,5)
        IF(PZM.LE.0D0) THEN
          PTS=0D0
        ELSEIF(K(IM,2).EQ.21.AND.IABS(K(N+1,2)).LE.10.AND.
     &  (MSTJ(44).EQ.3.OR.MSTJ(44).EQ.5)) THEN
          PTS=PMLS*ZM*(1D0-ZM)/V(IM,5)
        ELSEIF(MOD(MSTJ(43),2).EQ.1) THEN
          PTS=(PEM**2*(ZM*(1D0-ZM)*V(IM,5)-(1D0-ZM)*V(N+1,5)-
     &    ZM*V(N+2,5))-0.25D0*PMLS)/PZM**2
        ELSE
          PTS=PMLS*(ZM*(1D0-ZM)*PEM**2/V(IM,5)-0.25D0)/PZM**2
        ENDIF
        IF(PTS.LT.0D0.AND.LOOPPT.LT.10) THEN
          ZM=0.05D0+0.9D0*ZM
          GOTO 550
        ELSEIF(PTS.LT.0D0) THEN
          GOTO 280
        ENDIF
        PT=SQRT(MAX(0D0,PTS))
 
C...Global statistics.
        MINT(353)=MINT(353)+1
        VINT(353)=VINT(353)+PT
        IF (MINT(353).EQ.1) VINT(358)=PT
 
C...Find coefficient of azimuthal asymmetry due to gluon polarization.
        HAZIP=0D0
        IF(MSTJ(49).NE.1.AND.MOD(MSTJ(46),2).EQ.1.AND.K(IM,2).EQ.21
     &  .AND.IAU.NE.0) THEN
          IF(K(IGM,3).NE.0) MAZIP=1
          ZAU=V(IGM,1)
          IF(IAU.EQ.IM+1) ZAU=1D0-V(IGM,1)
          IF(MAZIP.EQ.0) ZAU=0D0
          IF(K(IGM,2).NE.21) THEN
            HAZIP=2D0*ZAU/(1D0+ZAU**2)
          ELSE
            HAZIP=(ZAU/(1D0-ZAU*(1D0-ZAU)))**2
          ENDIF
          IF(K(N+1,2).NE.21) THEN
            HAZIP=HAZIP*(-2D0*ZM*(1D0-ZM))/(1D0-2D0*ZM*(1D0-ZM))
          ELSE
            HAZIP=HAZIP*(ZM*(1D0-ZM)/(1D0-ZM*(1D0-ZM)))**2
          ENDIF
        ENDIF
 
C...Find coefficient of azimuthal asymmetry due to soft gluon
C...interference.
        HAZIC=0D0
        IF(MSTJ(49).NE.2.AND.MSTJ(46).GE.2.AND.(K(N+1,2).EQ.21.OR.
     &  K(N+2,2).EQ.21).AND.IAU.NE.0) THEN
          IF(K(IGM,3).NE.0) MAZIC=N+1
          IF(K(IGM,3).NE.0.AND.K(N+1,2).NE.21) MAZIC=N+2
          IF(K(IGM,3).NE.0.AND.K(N+1,2).EQ.21.AND.K(N+2,2).EQ.21.AND.
     &    ZM.GT.0.5D0) MAZIC=N+2
          IF(K(IAU,2).EQ.22) MAZIC=0
          ZS=ZM
          IF(MAZIC.EQ.N+2) ZS=1D0-ZM
          ZGM=V(IGM,1)
          IF(IAU.EQ.IM-1) ZGM=1D0-V(IGM,1)
          IF(MAZIC.EQ.0) ZGM=1D0
          IF(MAZIC.NE.0) HAZIC=(P(IM,5)/P(IGM,5))*
     &    SQRT((1D0-ZS)*(1D0-ZGM)/(ZS*ZGM))
          HAZIC=MIN(0.95D0,HAZIC)
        ENDIF
      ENDIF
 
C...Construct energies for ordinary branching in shower.
  560 IF(NEP.EQ.2.AND.IGM.GT.0) THEN
        IF(K(IM,2).EQ.21.AND.IABS(K(N+1,2)).LE.10.AND.
     &  (MSTJ(44).EQ.3.OR.MSTJ(44).EQ.5)) THEN
          P(N+1,4)=0.5D0*(PEM*(V(IM,5)+V(N+1,5)-V(N+2,5))+
     &    PZM*SQRT(MAX(0D0,PMLS))*(2D0*ZM-1D0))/V(IM,5)
        ELSEIF(MOD(MSTJ(43),2).EQ.1) THEN
          P(N+1,4)=PEM*V(IM,1)
        ELSE
          P(N+1,4)=PEM*(0.5D0*(V(IM,5)-SQRT(PMLS)+V(N+1,5)-V(N+2,5))+
     &    SQRT(PMLS)*ZM)/V(IM,5)
        ENDIF
 
C...Already predetermined choice of phi angle or not
        PHI=PARU(2)*PYR(0)
        IF(MPSPD.EQ.1.AND.IGM.EQ.NS+1) THEN
          IPSPD=IP1+IM-NS-2
          IF(K(IPSPD,4).GT.0) THEN
            IPSGD1=K(IPSPD,4)
            IF(IM.EQ.NS+2) THEN
              PHI=PYANGL(P(IPSGD1,1),P(IPSGD1,2))
            ELSE
              PHI=PYANGL(-P(IPSGD1,1),P(IPSGD1,2))
            ENDIF
          ENDIF
        ELSEIF(MPSPD.EQ.1.AND.IGM.EQ.NS+2) THEN
          IPSPD=IP1+IM-NS-2
          IF(K(IPSPD,4).GT.0) THEN
            IPSGD1=K(IPSPD,4)
            PHIPSM=PYANGL(P(IPSPD,1),P(IPSPD,2))
            THEPSM=PYANGL(P(IPSPD,3),SQRT(P(IPSPD,1)**2+P(IPSPD,2)**2))
            CALL PYROBO(IPSGD1,IPSGD1,0D0,-PHIPSM,0D0,0D0,0D0)
            CALL PYROBO(IPSGD1,IPSGD1,-THEPSM,0D0,0D0,0D0,0D0)
            PHI=PYANGL(P(IPSGD1,1),P(IPSGD1,2))
            CALL PYROBO(IPSGD1,IPSGD1,THEPSM,PHIPSM,0D0,0D0,0D0)
          ENDIF
        ENDIF
 
C...Construct momenta for ordinary branching in shower.
        P(N+1,1)=PT*COS(PHI)
        P(N+1,2)=PT*SIN(PHI)
        IF(K(IM,2).EQ.21.AND.IABS(K(N+1,2)).LE.10.AND.
     &  (MSTJ(44).EQ.3.OR.MSTJ(44).EQ.5)) THEN
          P(N+1,3)=0.5D0*(PZM*(V(IM,5)+V(N+1,5)-V(N+2,5))+
     &    PEM*SQRT(MAX(0D0,PMLS))*(2D0*ZM-1D0))/V(IM,5)
        ELSEIF(PZM.GT.0D0) THEN
          P(N+1,3)=0.5D0*(V(N+2,5)-V(N+1,5)-V(IM,5)+
     &    2D0*PEM*P(N+1,4))/PZM
        ELSE
          P(N+1,3)=0D0
        ENDIF
        P(N+2,1)=-P(N+1,1)
        P(N+2,2)=-P(N+1,2)
        P(N+2,3)=PZM-P(N+1,3)
        P(N+2,4)=PEM-P(N+1,4)
        IF(MSTJ(43).LE.2) THEN
          V(N+1,2)=(PEM*P(N+1,4)-PZM*P(N+1,3))/P(IM,5)
          V(N+2,2)=(PEM*P(N+2,4)-PZM*P(N+2,3))/P(IM,5)
        ENDIF
      ENDIF
 
C...Rotate and boost daughters.
      IF(IGM.GT.0) THEN
        IF(MSTJ(43).LE.2) THEN
          BEX=P(IGM,1)/P(IGM,4)
          BEY=P(IGM,2)/P(IGM,4)
          BEZ=P(IGM,3)/P(IGM,4)
          GA=P(IGM,4)/P(IGM,5)
          GABEP=GA*(GA*(BEX*P(IM,1)+BEY*P(IM,2)+BEZ*P(IM,3))/(1D0+GA)-
     &    P(IM,4))
        ELSE
          BEX=0D0
          BEY=0D0
          BEZ=0D0
          GA=1D0
          GABEP=0D0
        ENDIF
        PTIMB=SQRT((P(IM,1)+GABEP*BEX)**2+(P(IM,2)+GABEP*BEY)**2)
        THE=PYANGL(P(IM,3)+GABEP*BEZ,PTIMB)
        IF(PTIMB.GT.1D-4) THEN
          PHI=PYANGL(P(IM,1)+GABEP*BEX,P(IM,2)+GABEP*BEY)
        ELSE
          PHI=0D0
        ENDIF
        DO 570 I=N+1,N+2
          DP(1)=COS(THE)*COS(PHI)*P(I,1)-SIN(PHI)*P(I,2)+
     &    SIN(THE)*COS(PHI)*P(I,3)
          DP(2)=COS(THE)*SIN(PHI)*P(I,1)+COS(PHI)*P(I,2)+
     &    SIN(THE)*SIN(PHI)*P(I,3)
          DP(3)=-SIN(THE)*P(I,1)+COS(THE)*P(I,3)
          DP(4)=P(I,4)
          DBP=BEX*DP(1)+BEY*DP(2)+BEZ*DP(3)
          DGABP=GA*(GA*DBP/(1D0+GA)+DP(4))
          P(I,1)=DP(1)+DGABP*BEX
          P(I,2)=DP(2)+DGABP*BEY
          P(I,3)=DP(3)+DGABP*BEZ
          P(I,4)=GA*(DP(4)+DBP)
  570   CONTINUE
      ENDIF
 
C...Weight with azimuthal distribution, if required.
      IF(MAZIP.NE.0.OR.MAZIC.NE.0) THEN
        DO 580 J=1,3
          DPT(1,J)=P(IM,J)
          DPT(2,J)=P(IAU,J)
          DPT(3,J)=P(N+1,J)
  580   CONTINUE
        DPMA=DPT(1,1)*DPT(2,1)+DPT(1,2)*DPT(2,2)+DPT(1,3)*DPT(2,3)
        DPMD=DPT(1,1)*DPT(3,1)+DPT(1,2)*DPT(3,2)+DPT(1,3)*DPT(3,3)
        DPMM=DPT(1,1)**2+DPT(1,2)**2+DPT(1,3)**2
        DO 590 J=1,3
          DPT(4,J)=DPT(2,J)-DPMA*DPT(1,J)/MAX(1D-10,DPMM)
          DPT(5,J)=DPT(3,J)-DPMD*DPT(1,J)/MAX(1D-10,DPMM)
  590   CONTINUE
        DPT(4,4)=SQRT(DPT(4,1)**2+DPT(4,2)**2+DPT(4,3)**2)
        DPT(5,4)=SQRT(DPT(5,1)**2+DPT(5,2)**2+DPT(5,3)**2)
        IF(MIN(DPT(4,4),DPT(5,4)).GT.0.1D0*PARJ(82)) THEN
          CAD=(DPT(4,1)*DPT(5,1)+DPT(4,2)*DPT(5,2)+
     &    DPT(4,3)*DPT(5,3))/(DPT(4,4)*DPT(5,4))
          IF(MAZIP.NE.0) THEN
            IF(1D0+HAZIP*(2D0*CAD**2-1D0).LT.PYR(0)*(1D0+ABS(HAZIP)))
     &      GOTO 560
          ENDIF
          IF(MAZIC.NE.0) THEN
            IF(MAZIC.EQ.N+2) CAD=-CAD
            IF((1D0-HAZIC)*(1D0-HAZIC*CAD)/(1D0+HAZIC**2-2D0*HAZIC*CAD)
     &      .LT.PYR(0)) GOTO 560
          ENDIF
        ENDIF
      ENDIF
 
C...Azimuthal anisotropy due to interference with initial state partons.
      IF(MOD(MIIS,2).EQ.1.AND.IGM.EQ.NS+1.AND.(K(N+1,2).EQ.21.OR.
     &K(N+2,2).EQ.21)) THEN
        III=IM-NS-1
        IF(ISII(III).GE.1) THEN
          IAZIID=N+1
          IF(K(N+1,2).NE.21) IAZIID=N+2
          IF(K(N+1,2).EQ.21.AND.K(N+2,2).EQ.21.AND.
     &    P(N+1,4).GT.P(N+2,4)) IAZIID=N+2
          THEIID=PYANGL(P(IAZIID,3),SQRT(P(IAZIID,1)**2+P(IAZIID,2)**2))
          IF(III.EQ.2) THEIID=PARU(1)-THEIID
          PHIIID=PYANGL(P(IAZIID,1),P(IAZIID,2))
          HAZII=MIN(0.95D0,THEIID/THEIIS(III,ISII(III)))
          CAD=COS(PHIIID-PHIIIS(III,ISII(III)))
          PHIREL=ABS(PHIIID-PHIIIS(III,ISII(III)))
          IF(PHIREL.GT.PARU(1)) PHIREL=PARU(2)-PHIREL
          IF((1D0-HAZII)*(1D0-HAZII*CAD)/(1D0+HAZII**2-2D0*HAZII*CAD)
     &    .LT.PYR(0)) GOTO 560
        ENDIF
      ENDIF
 
C...Continue loop over partons that may branch, until none left.
      IF(IGM.GE.0) K(IM,1)=14
      N=N+NEP
      NEP=2
      IF(N.GT.MSTU(4)-MSTU(32)-10) THEN
        CALL PYERRM(11,'(PYSHOW:) no more memory left in PYJETS')
        IF(MSTU(21).GE.1) N=NS
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      GOTO 290
 
C...Set information on imagined shower initiator.
  600 IF(NPA.GE.2) THEN
        K(NS+1,1)=11
        K(NS+1,2)=94
        K(NS+1,3)=IP1
        IF(IP2.GT.0.AND.IP2.LT.IP1) K(NS+1,3)=IP2
        K(NS+1,4)=NS+2
        K(NS+1,5)=NS+1+NPA
        IIM=1
      ELSE
        IIM=0
      ENDIF
 
C...Reconstruct string drawing information.
      DO 610 I=NS+1+IIM,N
        KQ=KCHG(PYCOMP(K(I,2)),2)
        IF(K(I,1).LE.10.AND.K(I,2).EQ.22) THEN
          K(I,1)=1
        ELSEIF(K(I,1).LE.10.AND.IABS(K(I,2)).GE.11.AND.
     &    IABS(K(I,2)).LE.18) THEN
          K(I,1)=1
        ELSEIF(K(I,1).LE.10) THEN
          K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))
          K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))
        ELSEIF(K(MOD(K(I,4),MSTU(5))+1,2).NE.22) THEN
          ID1=MOD(K(I,4),MSTU(5))
          IF(KQ.EQ.1.AND.K(I,2).GT.0) ID1=MOD(K(I,4),MSTU(5))+1
          IF(KQ.EQ.2.AND.(K(ID1,2).EQ.21.OR.K(ID1+1,2).EQ.21).AND.
     &    PYR(0).GT.0.5D0) ID1=MOD(K(I,4),MSTU(5))+1
          ID2=2*MOD(K(I,4),MSTU(5))+1-ID1
          K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))+ID1
          K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))+ID2
          K(ID1,4)=K(ID1,4)+MSTU(5)*I
          K(ID1,5)=K(ID1,5)+MSTU(5)*ID2
          K(ID2,4)=K(ID2,4)+MSTU(5)*ID1
          K(ID2,5)=K(ID2,5)+MSTU(5)*I
        ELSE
          ID1=MOD(K(I,4),MSTU(5))
          ID2=ID1+1
          K(I,4)=MSTU(5)*(K(I,4)/MSTU(5))+ID1
          K(I,5)=MSTU(5)*(K(I,5)/MSTU(5))+ID1
          IF(KQ.EQ.1.OR.K(ID1,1).GE.11) THEN
            K(ID1,4)=K(ID1,4)+MSTU(5)*I
            K(ID1,5)=K(ID1,5)+MSTU(5)*I
          ELSE
            K(ID1,4)=0
            K(ID1,5)=0
          ENDIF
          K(ID2,4)=0
          K(ID2,5)=0
        ENDIF
  610 CONTINUE
 
C...Transformation from CM frame.
      IF(NPA.EQ.1) THEN
        THE=PYANGL(P(IPA(1),3),SQRT(P(IPA(1),1)**2+P(IPA(1),2)**2))
        PHI=PYANGL(P(IPA(1),1),P(IPA(1),2))
        MSTU(33)=1
        CALL PYROBO(NS+1,N,THE,PHI,0D0,0D0,0D0)
      ELSEIF(NPA.EQ.2) THEN
        BEX=PS(1)/PS(4)
        BEY=PS(2)/PS(4)
        BEZ=PS(3)/PS(4)
        GA=PS(4)/PS(5)
        GABEP=GA*(GA*(BEX*P(IPA(1),1)+BEY*P(IPA(1),2)+BEZ*P(IPA(1),3))
     &  /(1D0+GA)-P(IPA(1),4))
        THE=PYANGL(P(IPA(1),3)+GABEP*BEZ,SQRT((P(IPA(1),1)
     &  +GABEP*BEX)**2+(P(IPA(1),2)+GABEP*BEY)**2))
        PHI=PYANGL(P(IPA(1),1)+GABEP*BEX,P(IPA(1),2)+GABEP*BEY)
        MSTU(33)=1
        CALL PYROBO(NS+1,N,THE,PHI,BEX,BEY,BEZ)
      ELSE
        CALL PYROBO(IPA(1),IPA(NPA),0D0,0D0,PS(1)/PS(4),PS(2)/PS(4),
     &  PS(3)/PS(4))
        MSTU(33)=1
        CALL PYROBO(NS+1,N,0D0,0D0,PS(1)/PS(4),PS(2)/PS(4),PS(3)/PS(4))
      ENDIF
 
C...Decay vertex of shower.
      DO 630 I=NS+1,N
        DO 620 J=1,5
          V(I,J)=V(IP1,J)
  620   CONTINUE
  630 CONTINUE
 
C...Delete trivial shower, else connect initiators.
      IF(N.LE.NS+NPA+IIM) THEN
        N=NS
      ELSE
        DO 640 IP=1,NPA
          K(IPA(IP),1)=14
          K(IPA(IP),4)=K(IPA(IP),4)+NS+IIM+IP
          K(IPA(IP),5)=K(IPA(IP),5)+NS+IIM+IP
          K(NS+IIM+IP,3)=IPA(IP)
          IF(IIM.EQ.1.AND.MSTU(16).NE.2) K(NS+IIM+IP,3)=NS+1
          IF(K(NS+IIM+IP,1).NE.1) THEN
            K(NS+IIM+IP,4)=MSTU(5)*IPA(IP)+K(NS+IIM+IP,4)
            K(NS+IIM+IP,5)=MSTU(5)*IPA(IP)+K(NS+IIM+IP,5)
          ENDIF
  640   CONTINUE
      ENDIF
 
      RETURN
      END
