 
C*********************************************************************
 
C...PYBESQ
C...Calculates the momentum shift in a system of two particles assuming
C...the relative momentum squared should be shifted to Q2NEW. NI is the
C...last position occupied in /PYJETS/.
 
      SUBROUTINE PYBESQ(I1,I2,NI,Q2OLD,Q2NEW)
 
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
      SAVE /PYJETS/,/PYDAT1/
C...Local arrays and data.
      DIMENSION DP(5)
      SAVE HC1
 
      IF(MSTJ(55).EQ.0) THEN
        DQ2=Q2NEW-Q2OLD
        DP2=(P(I1,1)-P(I2,1))**2+(P(I1,2)-P(I2,2))**2+
     &  (P(I1,3)-P(I2,3))**2
        DP12=P(I1,1)**2+P(I1,2)**2+P(I1,3)**2
     &  -P(I2,1)**2-P(I2,2)**2-P(I2,3)**2
        SE=P(I1,4)+P(I2,4)
        DE=P(I1,4)-P(I2,4)
        DQ2SE=DQ2+SE**2
        DA=SE*DE*DP12-DP2*DQ2SE
        DB=DP2*DQ2SE-DP12**2
        HA=(DA+SQRT(MAX(DA**2+DQ2*(DQ2+SE**2-DE**2)*DB,0D0)))/(2D0*DB)
        DO 100 J=1,3
          PD=HA*(P(I1,J)-P(I2,J))
          P(NI+1,J)=PD
          P(NI+2,J)=-PD
  100   CONTINUE
        RETURN
      ENDIF
 
      K(NI+1,1)=1
      K(NI+2,1)=1
      DO 110 J=1,5
        P(NI+1,J)=P(I1,J)
        P(NI+2,J)=P(I2,J)
        DP(J)=P(I1,J)+P(I2,J)
  110 CONTINUE
 
C...Boost to cms and rotate first particle to z-axis
      CALL PYROBO(NI+1,NI+2,0.0D0,0.0D0,
     &-DP(1)/DP(4),-DP(2)/DP(4),-DP(3)/DP(4))
      PHI=PYANGL(P(NI+1,1),P(NI+1,2))
      THE=PYANGL(P(NI+1,3),SQRT(P(NI+1,1)**2+P(NI+1,2)**2))
      S=Q2NEW+(P(I1,5)+P(I2,5))**2
      PZ=0.5D0*SQRT(Q2NEW*(S-(P(I1,5)-P(I2,5))**2)/S)
      P(NI+1,1)=0.0D0
      P(NI+1,2)=0.0D0
      P(NI+1,3)=PZ
      P(NI+1,4)=SQRT(PZ**2+P(I1,5)**2)
      P(NI+2,1)=0.0D0
      P(NI+2,2)=0.0D0
      P(NI+2,3)=-PZ
      P(NI+2,4)=SQRT(PZ**2+P(I2,5)**2)
      DP(4)=SQRT(DP(1)**2+DP(2)**2+DP(3)**2+S)
      CALL PYROBO(NI+1,NI+2,THE,PHI,
     &DP(1)/DP(4),DP(2)/DP(4),DP(3)/DP(4))
 
      DO 120 J=1,3
        P(NI+1,J)=P(NI+1,J)-P(I1,J)
        P(NI+2,J)=P(NI+2,J)-P(I2,J)
  120 CONTINUE
 
      RETURN
      END
