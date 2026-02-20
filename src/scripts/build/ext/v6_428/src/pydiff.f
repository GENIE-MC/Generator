
C*********************************************************************
 
C...PYDIFF
C...Handles diffractive and elastic scattering.
 
      SUBROUTINE PYDIFF
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYJETS/,/PYDAT1/,/PYPARS/,/PYINT1/
 
C...Reset K, P and V vectors. Store incoming particles.
      DO 110 JT=1,MSTP(126)+10
        I=MINT(83)+JT
        DO 100 J=1,5
          K(I,J)=0
          P(I,J)=0D0
          V(I,J)=0D0
  100   CONTINUE
  110 CONTINUE
      N=MINT(84)
      MINT(3)=0
      MINT(21)=0
      MINT(22)=0
      MINT(23)=0
      MINT(24)=0
      MINT(4)=4
      DO 130 JT=1,2
        I=MINT(83)+JT
        K(I,1)=21
        K(I,2)=MINT(10+JT)
        DO 120 J=1,5
          P(I,J)=VINT(285+5*JT+J)
  120   CONTINUE
  130 CONTINUE
      MINT(6)=2
 
C...Subprocess; kinematics.
      SQLAM=(VINT(2)-VINT(63)-VINT(64))**2-4D0*VINT(63)*VINT(64)
      PZ=SQRT(SQLAM)/(2D0*VINT(1))
      DO 200 JT=1,2
        I=MINT(83)+JT
        PE=(VINT(2)+VINT(62+JT)-VINT(65-JT))/(2D0*VINT(1))
        KFH=MINT(102+JT)
 
C...Elastically scattered particle. (Except elastic GVMD states.)
        IF(MINT(16+JT).LE.0.AND.(MINT(10+JT).NE.22.OR.
     &  MINT(106+JT).NE.3)) THEN
          N=N+1
          K(N,1)=1
          K(N,2)=KFH
          K(N,3)=I+2
          P(N,3)=PZ*(-1)**(JT+1)
          P(N,4)=PE
          P(N,5)=SQRT(VINT(62+JT))
 
C...Decay rho from elastic scattering of gamma with sin**2(theta)
C...distribution of decay products (in rho rest frame).
          IF(KFH.EQ.113.AND.MINT(10+JT).EQ.22.AND.MSTP(102).EQ.1) THEN
            NSAV=N
            DBETAZ=P(N,3)/SQRT(P(N,3)**2+P(N,5)**2)
            P(N,3)=0D0
            P(N,4)=P(N,5)
            CALL PYDECY(NSAV)
            IF(N.EQ.NSAV+2.AND.IABS(K(NSAV+1,2)).EQ.211) THEN
              PHI=PYANGL(P(NSAV+1,1),P(NSAV+1,2))
              CALL PYROBO(NSAV+1,NSAV+2,0D0,-PHI,0D0,0D0,0D0)
              THE=PYANGL(P(NSAV+1,3),P(NSAV+1,1))
              CALL PYROBO(NSAV+1,NSAV+2,-THE,0D0,0D0,0D0,0D0)
  140         CTHE=2D0*PYR(0)-1D0
              IF(1D0-CTHE**2.LT.PYR(0)) GOTO 140
              CALL PYROBO(NSAV+1,NSAV+2,ACOS(CTHE),PHI,0D0,0D0,0D0)
            ENDIF
            CALL PYROBO(NSAV,NSAV+2,0D0,0D0,0D0,0D0,DBETAZ)
          ENDIF
 
C...Diffracted particle: low-mass system to two particles.
        ELSEIF(VINT(62+JT).LT.(VINT(66+JT)+PARP(103))**2) THEN
          N=N+2
          K(N-1,1)=1
          K(N,1)=1
          K(N-1,3)=I+2
          K(N,3)=I+2
          PMMAS=SQRT(VINT(62+JT))
          NTRY=0
  150     NTRY=NTRY+1
          IF(NTRY.LT.20) THEN
            MINT(105)=MINT(102+JT)
            MINT(109)=MINT(106+JT)
            CALL PYSPLI(KFH,21,KFL1,KFL2)
            CALL PYKFDI(KFL1,0,KFL3,KF1)
            IF(KF1.EQ.0) GOTO 150
            CALL PYKFDI(KFL2,-KFL3,KFLDUM,KF2)
            IF(KF2.EQ.0) GOTO 150
          ELSE
            KF1=KFH
            KF2=111
          ENDIF
          PM1=PYMASS(KF1)
          PM2=PYMASS(KF2)
          IF(PM1+PM2+PARJ(64).GT.PMMAS) GOTO 150
          K(N-1,2)=KF1
          K(N,2)=KF2
          P(N-1,5)=PM1
          P(N,5)=PM2
          PZP=SQRT(MAX(0D0,(PMMAS**2-PM1**2-PM2**2)**2-
     &    4D0*PM1**2*PM2**2))/(2D0*PMMAS)
          P(N-1,3)=PZP
          P(N,3)=-PZP
          P(N-1,4)=SQRT(PM1**2+PZP**2)
          P(N,4)=SQRT(PM2**2+PZP**2)
          CALL PYROBO(N-1,N,ACOS(2D0*PYR(0)-1D0),PARU(2)*PYR(0),
     &    0D0,0D0,0D0)
          DBETAZ=PZ*(-1)**(JT+1)/SQRT(PZ**2+PMMAS**2)
          CALL PYROBO(N-1,N,0D0,0D0,0D0,0D0,DBETAZ)
 
C...Diffracted particle: valence quark kicked out.
        ELSEIF(MSTP(101).EQ.1.OR.(MSTP(101).EQ.3.AND.PYR(0).LT.
     &    PARP(101))) THEN
          N=N+2
          K(N-1,1)=2
          K(N,1)=1
          K(N-1,3)=I+2
          K(N,3)=I+2
          MINT(105)=MINT(102+JT)
          MINT(109)=MINT(106+JT)
          CALL PYSPLI(KFH,21,K(N,2),K(N-1,2))
          P(N-1,5)=PYMASS(K(N-1,2))
          P(N,5)=PYMASS(K(N,2))
          SQLAM=(VINT(62+JT)-P(N-1,5)**2-P(N,5)**2)**2-
     &    4D0*P(N-1,5)**2*P(N,5)**2
          P(N-1,3)=(PE*SQRT(SQLAM)+PZ*(VINT(62+JT)+P(N-1,5)**2-
     &    P(N,5)**2))/(2D0*VINT(62+JT))*(-1)**(JT+1)
          P(N-1,4)=SQRT(P(N-1,3)**2+P(N-1,5)**2)
          P(N,3)=PZ*(-1)**(JT+1)-P(N-1,3)
          P(N,4)=SQRT(P(N,3)**2+P(N,5)**2)
 
C...Diffracted particle: gluon kicked out.
        ELSE
          N=N+3
          K(N-2,1)=2
          K(N-1,1)=2
          K(N,1)=1
          K(N-2,3)=I+2
          K(N-1,3)=I+2
          K(N,3)=I+2
          MINT(105)=MINT(102+JT)
          MINT(109)=MINT(106+JT)
          CALL PYSPLI(KFH,21,K(N,2),K(N-2,2))
          K(N-1,2)=21
          P(N-2,5)=PYMASS(K(N-2,2))
          P(N-1,5)=0D0
          P(N,5)=PYMASS(K(N,2))
C...Energy distribution for particle into two jets.
  160     IMB=1
          IF(MOD(KFH/1000,10).NE.0) IMB=2
          CHIK=PARP(92+2*IMB)
          IF(MSTP(92).LE.1) THEN
            IF(IMB.EQ.1) CHI=PYR(0)
            IF(IMB.EQ.2) CHI=1D0-SQRT(PYR(0))
          ELSEIF(MSTP(92).EQ.2) THEN
            CHI=1D0-PYR(0)**(1D0/(1D0+CHIK))
          ELSEIF(MSTP(92).EQ.3) THEN
            CUT=2D0*0.3D0/VINT(1)
  170       CHI=PYR(0)**2
            IF((CHI**2/(CHI**2+CUT**2))**0.25D0*(1D0-CHI)**CHIK.LT.
     &      PYR(0)) GOTO 170
          ELSEIF(MSTP(92).EQ.4) THEN
            CUT=2D0*0.3D0/VINT(1)
            CUTR=(1D0+SQRT(1D0+CUT**2))/CUT
  180       CHIR=CUT*CUTR**PYR(0)
            CHI=(CHIR**2-CUT**2)/(2D0*CHIR)
            IF((1D0-CHI)**CHIK.LT.PYR(0)) GOTO 180
          ELSE
            CUT=2D0*0.3D0/VINT(1)
            CUTA=CUT**(1D0-PARP(98))
            CUTB=(1D0+CUT)**(1D0-PARP(98))
  190       CHI=(CUTA+PYR(0)*(CUTB-CUTA))**(1D0/(1D0-PARP(98)))
            IF(((CHI+CUT)**2/(2D0*(CHI**2+CUT**2)))**
     &      (0.5D0*PARP(98))*(1D0-CHI)**CHIK.LT.PYR(0)) GOTO 190
          ENDIF
          IF(CHI.LT.P(N,5)**2/VINT(62+JT).OR.CHI.GT.1D0-P(N-2,5)**2/
     &    VINT(62+JT)) GOTO 160
          SQM=P(N-2,5)**2/(1D0-CHI)+P(N,5)**2/CHI
          PZI=(PE*(VINT(62+JT)-SQM)+PZ*(VINT(62+JT)+SQM))/
     &    (2D0*VINT(62+JT))
          PEI=SQRT(PZI**2+SQM)
          PQQP=(1D0-CHI)*(PEI+PZI)
          P(N-2,3)=0.5D0*(PQQP-P(N-2,5)**2/PQQP)*(-1)**(JT+1)
          P(N-2,4)=SQRT(P(N-2,3)**2+P(N-2,5)**2)
          P(N-1,4)=0.5D0*(VINT(62+JT)-SQM)/(PEI+PZI)
          P(N-1,3)=P(N-1,4)*(-1)**JT
          P(N,3)=PZI*(-1)**(JT+1)-P(N-2,3)
          P(N,4)=SQRT(P(N,3)**2+P(N,5)**2)
        ENDIF
 
C...Documentation lines.
        K(I+2,1)=21
        IF(MINT(16+JT).EQ.0) K(I+2,2)=KFH
        IF(MINT(16+JT).NE.0.OR.(MINT(10+JT).EQ.22.AND.
     &  MINT(106+JT).EQ.3)) K(I+2,2)=ISIGN(9900000,KFH)+10*(KFH/10)
        K(I+2,3)=I
        P(I+2,3)=PZ*(-1)**(JT+1)
        P(I+2,4)=PE
        P(I+2,5)=SQRT(VINT(62+JT))
  200 CONTINUE
 
C...Rotate outgoing partons/particles using cos(theta).
      IF(VINT(23).LT.0.9D0) THEN
        CALL PYROBO(MINT(83)+3,N,ACOS(VINT(23)),VINT(24),0D0,0D0,0D0)
      ELSE
        CALL PYROBO(MINT(83)+3,N,ASIN(VINT(59)),VINT(24),0D0,0D0,0D0)
      ENDIF
 
      RETURN
      END
