 
C*********************************************************************
 
C...PYDISG
C...Set up a DIS process as gamma* + f -> f, with beam remnant
C...and showering added consecutively. Photon flux by the PYGAGA
C...routine (if at all).
 
      SUBROUTINE PYDISG
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYSUBS/,/PYPARS/,/PYINT1/
C...Local arrays.
      DIMENSION PMS(4)
 
C...Choice of subprocess, number of documentation lines
      IDOC=7
      MINT(3)=IDOC-6
      MINT(4)=IDOC
      IPU1=MINT(84)+1
      IPU2=MINT(84)+2
      IPU3=MINT(84)+3
      ISIDE=1
      IF(MINT(107).EQ.4) ISIDE=2
 
C...Reset K, P and V vectors. Store incoming particles
      DO 110 JT=1,MSTP(126)+20
        I=MINT(83)+JT
        DO 100 J=1,5
          K(I,J)=0
          P(I,J)=0D0
          V(I,J)=0D0
  100   CONTINUE
  110 CONTINUE
      DO 130 JT=1,2
        I=MINT(83)+JT
        K(I,1)=21
        K(I,2)=MINT(10+JT)
        DO 120 J=1,5
          P(I,J)=VINT(285+5*JT+J)
  120   CONTINUE
  130 CONTINUE
      MINT(6)=2
 
C...Store incoming partons in hadronic CM-frame
      DO 140 JT=1,2
        I=MINT(84)+JT
        K(I,1)=14
        K(I,2)=MINT(14+JT)
        K(I,3)=MINT(83)+2+JT
  140 CONTINUE
      IF(MINT(15).EQ.22) THEN
        P(MINT(84)+1,3)=0.5D0*(VINT(1)+VINT(307)/VINT(1))
        P(MINT(84)+1,4)=0.5D0*(VINT(1)-VINT(307)/VINT(1))
        P(MINT(84)+1,5)=-SQRT(VINT(307))
        P(MINT(84)+2,3)=-0.5D0*VINT(307)/VINT(1)
        P(MINT(84)+2,4)=0.5D0*VINT(307)/VINT(1)
        KFRES=MINT(16)
        ISIDE=2
      ELSE
        P(MINT(84)+1,3)=0.5D0*VINT(308)/VINT(1)
        P(MINT(84)+1,4)=0.5D0*VINT(308)/VINT(1)
        P(MINT(84)+2,3)=-0.5D0*(VINT(1)+VINT(308)/VINT(1))
        P(MINT(84)+2,4)=0.5D0*(VINT(1)-VINT(308)/VINT(1))
        P(MINT(84)+1,5)=-SQRT(VINT(308))
        KFRES=MINT(15)
        ISIDE=1
      ENDIF
      SIDESG=(-1D0)**(ISIDE-1)
 
C...Copy incoming partons to documentation lines.
      DO 170 JT=1,2
        I1=MINT(83)+4+JT
        I2=MINT(84)+JT
        K(I1,1)=21
        K(I1,2)=K(I2,2)
        K(I1,3)=I1-2
        DO 150 J=1,5
          P(I1,J)=P(I2,J)
  150   CONTINUE
 
C...Second copy for partons before ISR shower, since no such.
        I1=MINT(83)+2+JT
        K(I1,1)=21
        K(I1,2)=K(I2,2)
        K(I1,3)=I1-2
        DO 160 J=1,5
          P(I1,J)=P(I2,J)
  160   CONTINUE
  170 CONTINUE
 
C...Define initial partons.
      NTRY=0
  180 NTRY=NTRY+1
      IF(NTRY.GT.100) THEN
        MINT(51)=1
        RETURN
      ENDIF
 
C...Scattered quark in hadronic CM frame.
      I=MINT(83)+7
      K(IPU3,1)=3
      K(IPU3,2)=KFRES
      K(IPU3,3)=I
      P(IPU3,5)=PYMASS(KFRES)
      P(IPU3,3)=P(IPU1,3)+P(IPU2,3)
      P(IPU3,4)=P(IPU1,4)+P(IPU2,4)
      P(IPU3,5)=0D0
      K(I,1)=21
      K(I,2)=KFRES
      K(I,3)=MINT(83)+4+ISIDE
      P(I,3)=P(IPU3,3)
      P(I,4)=P(IPU3,4)
      P(I,5)=P(IPU3,5)
      N=IPU3
      MINT(21)=KFRES
      MINT(22)=0
 
C...No primordial kT, or chosen according to truncated Gaussian or
C...exponential, or (for photon) predetermined or power law.
  190 IF(MINT(40+ISIDE).EQ.2.AND.MINT(10+ISIDE).NE.22) THEN
        IF(MSTP(91).LE.0) THEN
          PT=0D0
        ELSEIF(MSTP(91).EQ.1) THEN
          PT=PARP(91)*SQRT(-LOG(PYR(0)))
        ELSE
          RPT1=PYR(0)
          RPT2=PYR(0)
          PT=-PARP(92)*LOG(RPT1*RPT2)
        ENDIF
        IF(PT.GT.PARP(93)) GOTO 190
      ELSEIF(MINT(106+ISIDE).EQ.3) THEN
        PTA=SQRT(VINT(282+ISIDE))
        PTB=0D0
        IF(MSTP(66).EQ.5.AND.MSTP(93).EQ.1) THEN
          PTB=PARP(99)*SQRT(-LOG(PYR(0)))
        ELSEIF(MSTP(66).EQ.5.AND.MSTP(93).EQ.2) THEN
          RPT1=PYR(0)
          RPT2=PYR(0)
          PTB=-PARP(99)*LOG(RPT1*RPT2)
        ENDIF
        IF(PTB.GT.PARP(100)) GOTO 190
        PT=SQRT(PTA**2+PTB**2+2D0*PTA*PTB*COS(PARU(2)*PYR(0)))
        IF(NTRY.GT.10) PT=PT*0.8D0**(NTRY-10)
      ELSEIF(IABS(MINT(14+ISIDE)).LE.8.OR.MINT(14+ISIDE).EQ.21) THEN
        IF(MSTP(93).LE.0) THEN
          PT=0D0
        ELSEIF(MSTP(93).EQ.1) THEN
          PT=PARP(99)*SQRT(-LOG(PYR(0)))
        ELSEIF(MSTP(93).EQ.2) THEN
          RPT1=PYR(0)
          RPT2=PYR(0)
          PT=-PARP(99)*LOG(RPT1*RPT2)
        ELSEIF(MSTP(93).EQ.3) THEN
          HA=PARP(99)**2
          HB=PARP(100)**2
          PT=SQRT(MAX(0D0,HA*(HA+HB)/(HA+HB-PYR(0)*HB)-HA))
        ELSE
          HA=PARP(99)**2
          HB=PARP(100)**2
          IF(MSTP(93).EQ.5) HB=MIN(VINT(48),PARP(100)**2)
          PT=SQRT(MAX(0D0,HA*((HA+HB)/HA)**PYR(0)-HA))
        ENDIF
        IF(PT.GT.PARP(100)) GOTO 190
      ELSE
        PT=0D0
      ENDIF
      VINT(156+ISIDE)=PT
      PHI=PARU(2)*PYR(0)
      P(IPU3,1)=PT*COS(PHI)
      P(IPU3,2)=PT*SIN(PHI)
      P(IPU3,4)=SQRT(P(IPU3,5)**2+PT**2+P(IPU3,3)**2)
      PMS(3-ISIDE)=P(IPU3,5)**2+P(IPU3,1)**2+P(IPU3,2)**2
      PCP=P(IPU3,4)+ABS(P(IPU3,3))
 
C...Find one or two beam remnants.
      MINT(105)=MINT(102+ISIDE)
      MINT(109)=MINT(106+ISIDE)
      CALL PYSPLI(MINT(10+ISIDE),MINT(12+ISIDE),KFLCH,KFLSP)
      IF(MINT(51).NE.0) THEN
        MINT(51)=0
        GOTO 180
      ENDIF
 
C...Store first remnant parton, with colour info and kinematics.
      I=N+1
      K(I,1)=1
      K(I,2)=KFLSP
      K(I,3)=MINT(83)+ISIDE
      P(I,5)=PYMASS(K(I,2))
      KCOL=KCHG(PYCOMP(KFLSP),2)
      IF(KCOL.NE.0) THEN
        K(I,1)=3
        KFLS=(3-KCOL*ISIGN(1,KFLSP))/2
        K(I,KFLS+3)=MSTU(5)*IPU3
        K(IPU3,6-KFLS)=MSTU(5)*I
        ICOLR=I
      ENDIF
      IF(KFLCH.EQ.0) THEN
        P(I,1)=-P(IPU3,1)
        P(I,2)=-P(IPU3,2)
        PMS(ISIDE)=P(I,5)**2+P(I,1)**2+P(I,2)**2
        P(I,3)=-P(IPU3,3)
        P(I,4)=SQRT(PMS(ISIDE)+P(I,3)**2)
        PRP=P(I,4)+ABS(P(I,3))
 
C...When extra remnant parton or hadron: store extra remnant.
      ELSE
        I=I+1
        K(I,1)=1
        K(I,2)=KFLCH
        K(I,3)=MINT(83)+ISIDE
        P(I,5)=PYMASS(K(I,2))
        KCOL=KCHG(PYCOMP(KFLCH),2)
        IF(KCOL.NE.0) THEN
          K(I,1)=3
          KFLS=(3-KCOL*ISIGN(1,KFLCH))/2
          K(I,KFLS+3)=MSTU(5)*IPU3
          K(IPU3,6-KFLS)=MSTU(5)*I
          ICOLR=I
        ENDIF
 
C...Relative transverse momentum when two remnants.
        LOOP=0
  200   LOOP=LOOP+1
        CALL PYPTDI(1,P(I-1,1),P(I-1,2))
        P(I-1,1)=P(I-1,1)-0.5D0*P(IPU3,1)
        P(I-1,2)=P(I-1,2)-0.5D0*P(IPU3,2)
        PMS(3)=P(I-1,5)**2+P(I-1,1)**2+P(I-1,2)**2
        P(I,1)=-P(IPU3,1)-P(I-1,1)
        P(I,2)=-P(IPU3,2)-P(I-1,2)
        PMS(4)=P(I,5)**2+P(I,1)**2+P(I,2)**2
 
C...Relative distribution of energy for particle into jet plus particle.
        IMB=1
        IF(MOD(MINT(10+ISIDE)/1000,10).NE.0) IMB=2
        IF(MSTP(94).LE.1) THEN
          IF(IMB.EQ.1) CHI=PYR(0)
          IF(IMB.EQ.2) CHI=1D0-SQRT(PYR(0))
          IF(MOD(KFLCH/1000,10).NE.0) CHI=1D0-CHI
        ELSEIF(MSTP(94).EQ.2) THEN
          CHI=1D0-PYR(0)**(1D0/(1D0+PARP(93+2*IMB)))
          IF(MOD(KFLCH/1000,10).NE.0) CHI=1D0-CHI
        ELSEIF(MSTP(94).EQ.3) THEN
          CALL PYZDIS(1,0,PMS(4),ZZ)
          CHI=ZZ
        ELSE
          CALL PYZDIS(1000,0,PMS(4),ZZ)
          CHI=ZZ
        ENDIF
 
C...Construct total transverse mass; reject if too large.
        CHI=MAX(1D-8,MIN(1D0-1D-8,CHI))
        PMS(ISIDE)=PMS(4)/CHI+PMS(3)/(1D0-CHI)
        IF(PMS(ISIDE).GT.P(IPU3,4)**2) THEN
          IF(LOOP.LT.10) GOTO 200
          GOTO 180
        ENDIF
        VINT(158+ISIDE)=CHI
 
C...Subdivide longitudinal momentum according to value selected above.
        PRP=SQRT(PMS(ISIDE)+P(IPU3,3)**2)+ABS(P(IPU3,3))
        PW1=(1D0-CHI)*PRP
        P(I-1,4)=0.5D0*(PW1+PMS(3)/PW1)
        P(I-1,3)=0.5D0*(PW1-PMS(3)/PW1)*SIDESG
        PW2=CHI*PRP
        P(I,4)=0.5D0*(PW2+PMS(4)/PW2)
        P(I,3)=0.5D0*(PW2-PMS(4)/PW2)*SIDESG
      ENDIF
      N=I
 
C...Boost current and remnant systems to correct frame.
      IF(SQRT(PMS(1))+SQRT(PMS(2)).GT.0.99D0*VINT(1)) GOTO 180
      DSQLAM=SQRT(MAX(0D0,(VINT(2)-PMS(1)-PMS(2))**2-4D0*PMS(1)*PMS(2)))
      DRKC=(VINT(2)+PMS(3-ISIDE)-PMS(ISIDE)+DSQLAM)/
     &(2D0*VINT(1)*PCP)
      DRKR=(VINT(2)+PMS(ISIDE)-PMS(3-ISIDE)+DSQLAM)/
     &(2D0*VINT(1)*PRP)
      DBEC=-SIDESG*(DRKC**2-1D0)/(DRKC**2+1D0)
      DBER=SIDESG*(DRKR**2-1D0)/(DRKR**2+1D0)
      CALL PYROBO(IPU3,IPU3,0D0,0D0,0D0,0D0,DBEC)
      CALL PYROBO(IPU3+1,N,0D0,0D0,0D0,0D0,DBER)
 
C...Let current quark shower; recoil but no showering by colour partner.
      QMAX=2D0*SQRT(VINT(309-ISIDE))
      MSTJ48=MSTJ(48)
      MSTJ(48)=1
      PARJ86=PARJ(86)
      PARJ(86)=0D0
      IF(MSTP(71).EQ.1) CALL PYSHOW(IPU3,ICOLR,QMAX)
      MSTJ(48)=MSTJ48
      PARJ(86)=PARJ86
 
      RETURN
      END
