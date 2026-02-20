 
C*********************************************************************
 
C...PYPDEL
C...Gives electron (or muon, or tau) parton distribution.
 
      SUBROUTINE PYPDEL(KFA,X,Q2,XPEL)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/
C...Local arrays.
      DIMENSION XPEL(-25:25),XPGA(-6:6),SXP(0:6)
 
C...Interface to PDFLIB.
      COMMON/W50513/XMIN,XMAX,Q2MIN,Q2MAX
      SAVE /W50513/
      DOUBLE PRECISION XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU,
     &VALUE(20),XMIN,XMAX,Q2MIN,Q2MAX
      CHARACTER*20 PARM(20)
      DATA VALUE/20*0D0/,PARM/20*' '/
 
C...Some common constants.
      DO 100 KFL=-25,25
        XPEL(KFL)=0D0
  100 CONTINUE
      AEM=PARU(101)
      PME=PMAS(11,1)
      IF(KFA.EQ.13) PME=PMAS(13,1)
      IF(KFA.EQ.15) PME=PMAS(15,1)
      XL=LOG(MAX(1D-10,X))
      X1L=LOG(MAX(1D-10,1D0-X))
      HLE=LOG(MAX(3D0,Q2/PME**2))
      HBE2=(AEM/PARU(1))*(HLE-1D0)
 
C...Electron inside electron, see R. Kleiss et al., in Z physics at
C...LEP 1, CERN 89-08, p. 34
      IF(MSTP(59).LE.1) THEN
        HDE=1D0+(AEM/PARU(1))*(1.5D0*HLE+1.289868D0)+(AEM/PARU(1))**2*
     &  (-2.164868D0*HLE**2+9.840808D0*HLE-10.130464D0)
        HEE=HBE2*(1D0-X)**(HBE2-1D0)*SQRT(MAX(0D0,HDE))-
     &  0.5D0*HBE2*(1D0+X)+HBE2**2/8D0*((1D0+X)*(-4D0*X1L+3D0*XL)-
     &  4D0*XL/(1D0-X)-5D0-X)
      ELSE
        HEE=HBE2*(1D0-X)**(HBE2-1D0)*EXP(0.172784D0*HBE2)/
     &  PYGAMM(1D0+HBE2)-0.5D0*HBE2*(1D0+X)+HBE2**2/8D0*((1D0+X)*
     &  (-4D0*X1L+3D0*XL)-4D0*XL/(1D0-X)-5D0-X)
      ENDIF
C...Zero distribution for very large x and rescale it for intermediate.
      IF(X.GT.1D0-1D-10) THEN
        HEE=0D0
      ELSEIF(X.GT.1D0-1D-7) THEN
        HEE=HEE*1000D0**HBE2/(1000D0**HBE2-1D0)
      ENDIF
      XPEL(KFA)=X*HEE
 
C...Photon and (transverse) W- inside electron.
      AEMP=PYALEM(PME*SQRT(MAX(0D0,Q2)))/PARU(2)
      IF(MSTP(13).LE.1) THEN
        HLG=HLE
      ELSE
        HLG=LOG(MAX(1D0,(PARP(13)/PME**2)*(1D0-X)/X**2))
      ENDIF
      XPEL(22)=AEMP*HLG*(1D0+(1D0-X)**2)
      HLW=LOG(1D0+Q2/PMAS(24,1)**2)/(4D0*PARU(102))
      XPEL(-24)=AEMP*HLW*(1D0+(1D0-X)**2)
 
C...Electron or positron inside photon inside electron.
      IF(KFA.EQ.11.AND.MSTP(12).EQ.1) THEN
        XFSEA=0.5D0*(AEMP*(HLE-1D0))**2*(4D0/3D0+X-X**2-4D0*X**3/3D0+
     &  2D0*X*(1D0+X)*XL)
        XPEL(11)=XPEL(11)+XFSEA
        XPEL(-11)=XFSEA
 
C...Initialize PDFLIB photon parton distributions.
        IF(MSTP(56).EQ.2) THEN
          PARM(1)='NPTYPE'
          VALUE(1)=3
          PARM(2)='NGROUP'
          VALUE(2)=MSTP(55)/1000
          PARM(3)='NSET'
          VALUE(3)=MOD(MSTP(55),1000)
          IF(MINT(93).NE.3000000+MSTP(55)) THEN
            CALL PDFSET(PARM,VALUE)
            MINT(93)=3000000+MSTP(55)
          ENDIF
        ENDIF
 
C...Quarks and gluons inside photon inside electron:
C...numerical convolution required.
        DO 110 KFL=0,6
          SXP(KFL)=0D0
  110   CONTINUE
        SUMXPP=0D0
        ITER=-1
  120   ITER=ITER+1
        SUMXP=SUMXPP
        NSTP=2**(ITER-1)
        IF(ITER.EQ.0) NSTP=2
        DO 130 KFL=0,6
          SXP(KFL)=0.5D0*SXP(KFL)
  130   CONTINUE
        WTSTP=0.5D0/NSTP
        IF(ITER.EQ.0) WTSTP=0.5D0
C...Pick grid of x_{gamma} values logarithmically even.
        DO 150 ISTP=1,NSTP
          IF(ITER.EQ.0) THEN
            XLE=XL*(ISTP-1)
          ELSE
            XLE=XL*(ISTP-0.5D0)/NSTP
          ENDIF
          XE=MIN(1D0-1D-10,EXP(XLE))
          XG=MIN(1D0-1D-10,X/XE)
C...Evaluate photon inside electron parton distribution for convolution.
          XPGP=1D0+(1D0-XE)**2
          IF(MSTP(13).LE.1) THEN
            XPGP=XPGP*HLE
          ELSE
            XPGP=XPGP*LOG(MAX(1D0,(PARP(13)/PME**2)*(1D0-XE)/XE**2))
          ENDIF
C...Evaluate photon parton distributions for convolution.
          IF(MSTP(56).EQ.1) THEN
            IF(MSTP(55).EQ.1) THEN
              CALL PYPDGA(XG,Q2,XPGA)
            ELSEIF(MSTP(55).GE.5.AND.MSTP(55).LE.8) THEN
              Q2MX=Q2
              P2MX=0.36D0
              IF(MSTP(55).GE.7) P2MX=4.0D0
              IF(MSTP(57).EQ.0) Q2MX=P2MX
              P2=0D0
              IF(VINT(120).LT.0D0) P2=VINT(120)**2
              CALL PYGGAM(MSTP(55)-4,XG,Q2MX,P2,MSTP(60),F2GAM,XPGA)
              VINT(231)=P2MX
            ELSEIF(MSTP(55).GE.9.AND.MSTP(55).LE.12) THEN
              Q2MX=Q2
              P2MX=0.36D0
              IF(MSTP(55).GE.11) P2MX=4.0D0
              IF(MSTP(57).EQ.0) Q2MX=P2MX
              P2=0D0
              IF(VINT(120).LT.0D0) P2=VINT(120)**2
              CALL PYGGAM(MSTP(55)-8,XG,Q2MX,P2,MSTP(60),F2GAM,XPGA)
              VINT(231)=P2MX
            ENDIF
            DO 140 KFL=0,5
              SXP(KFL)=SXP(KFL)+WTSTP*XPGP*XPGA(KFL)
  140       CONTINUE
          ELSEIF(MSTP(56).EQ.2) THEN
C...Call PDFLIB parton distributions.
            XX=XG
            QQ=SQRT(MAX(0D0,Q2MIN,Q2))
            IF(MSTP(57).EQ.0) QQ=SQRT(Q2MIN)
            CALL STRUCTM(XX,QQ,UPV,DNV,USEA,DSEA,STR,CHM,BOT,TOP,GLU)
            SXP(0)=SXP(0)+WTSTP*XPGP*GLU
            SXP(1)=SXP(1)+WTSTP*XPGP*DNV
            SXP(2)=SXP(2)+WTSTP*XPGP*UPV
            SXP(3)=SXP(3)+WTSTP*XPGP*STR
            SXP(4)=SXP(4)+WTSTP*XPGP*CHM
            SXP(5)=SXP(5)+WTSTP*XPGP*BOT
            SXP(6)=SXP(6)+WTSTP*XPGP*TOP
          ENDIF
  150   CONTINUE
        SUMXPP=SXP(0)+2D0*SXP(1)+2D0*SXP(2)
        IF(ITER.LE.2.OR.(ITER.LE.7.AND.ABS(SUMXPP-SUMXP).GT.
     &  PARP(14)*(SUMXPP+SUMXP))) GOTO 120
 
C...Put convolution into output arrays.
        FCONV=AEMP*(-XL)
        XPEL(0)=FCONV*SXP(0)
        DO 160 KFL=1,6
          XPEL(KFL)=FCONV*SXP(KFL)
          XPEL(-KFL)=XPEL(KFL)
  160   CONTINUE
      ENDIF
 
      RETURN
      END
