 
C*********************************************************************
 
C...PYP
C...Provides various real-valued event related data.
 
      FUNCTION PYP(I,J)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
C...Local array.
      DIMENSION PSUM(4)
 
C...Set default value. For I = 0 sum of momenta or charges,
C...or invariant mass of system.
      PYP=0D0
      IF(I.LT.0.OR.I.GT.MSTU(4).OR.J.LE.0) THEN
      ELSEIF(I.EQ.0.AND.J.LE.4) THEN
        DO 100 I1=1,N
          IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PYP=PYP+P(I1,J)
  100   CONTINUE
      ELSEIF(I.EQ.0.AND.J.EQ.5) THEN
        DO 120 J1=1,4
          PSUM(J1)=0D0
          DO 110 I1=1,N
            IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PSUM(J1)=PSUM(J1)+
     &      P(I1,J1)
  110     CONTINUE
  120   CONTINUE
        PYP=SQRT(MAX(0D0,PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2))
      ELSEIF(I.EQ.0.AND.J.EQ.6) THEN
        DO 130 I1=1,N
          IF(K(I1,1).GT.0.AND.K(I1,1).LE.10) PYP=PYP+PYCHGE(K(I1,2))/3D0
  130   CONTINUE
      ELSEIF(I.EQ.0) THEN
 
C...Direct readout of P matrix.
      ELSEIF(J.LE.5) THEN
        PYP=P(I,J)
 
C...Charge, total momentum, transverse momentum, transverse mass.
      ELSEIF(J.LE.12) THEN
        IF(J.EQ.6) PYP=PYCHGE(K(I,2))/3D0
        IF(J.EQ.7.OR.J.EQ.8) PYP=P(I,1)**2+P(I,2)**2+P(I,3)**2
        IF(J.EQ.9.OR.J.EQ.10) PYP=P(I,1)**2+P(I,2)**2
        IF(J.EQ.11.OR.J.EQ.12) PYP=P(I,5)**2+P(I,1)**2+P(I,2)**2
        IF(J.EQ.8.OR.J.EQ.10.OR.J.EQ.12) PYP=SQRT(PYP)
 
C...Theta and phi angle in radians or degrees.
      ELSEIF(J.LE.16) THEN
        IF(J.LE.14) PYP=PYANGL(P(I,3),SQRT(P(I,1)**2+P(I,2)**2))
        IF(J.GE.15) PYP=PYANGL(P(I,1),P(I,2))
        IF(J.EQ.14.OR.J.EQ.16) PYP=PYP*180D0/PARU(1)
 
C...True rapidity, rapidity with pion mass, pseudorapidity.
      ELSEIF(J.LE.19) THEN
        PMR=0D0
        IF(J.EQ.17) PMR=P(I,5)
        IF(J.EQ.18) PMR=PYMASS(211)
        PR=MAX(1D-20,PMR**2+P(I,1)**2+P(I,2)**2)
        PYP=SIGN(LOG(MIN((SQRT(PR+P(I,3)**2)+ABS(P(I,3)))/SQRT(PR),
     &  1D20)),P(I,3))
 
C...Energy and momentum fractions (only to be used in CM frame).
      ELSEIF(J.LE.25) THEN
        IF(J.EQ.20) PYP=2D0*SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)/PARU(21)
        IF(J.EQ.21) PYP=2D0*P(I,3)/PARU(21)
        IF(J.EQ.22) PYP=2D0*SQRT(P(I,1)**2+P(I,2)**2)/PARU(21)
        IF(J.EQ.23) PYP=2D0*P(I,4)/PARU(21)
        IF(J.EQ.24) PYP=(P(I,4)+P(I,3))/PARU(21)
        IF(J.EQ.25) PYP=(P(I,4)-P(I,3))/PARU(21)
      ENDIF
 
      RETURN
      END
