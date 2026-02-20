 
C*********************************************************************
 
C...PY4JTW
C...Auxiliary to PY4JET, to evaluate weight of configuration.
 
      FUNCTION PY4JTW(IA1,IA2,IA3,IA4)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      SAVE /PYJETS/
 
C...First case: when both original partons radiate.
C...IA1 /= 0: N+1 -> IA1 + IA2, N+2 -> IA3 + IA4.
      IF(IA1.NE.0) THEN
        DO 100 J=1,4
          P(N+1,J)=P(IA1,J)+P(IA2,J)
          P(N+2,J)=P(IA3,J)+P(IA4,J)
  100   CONTINUE
        P(N+1,5)=SQRT(MAX(0D0,P(N+1,4)**2-P(N+1,1)**2-P(N+1,2)**2-
     &  P(N+1,3)**2))
        P(N+2,5)=SQRT(MAX(0D0,P(N+2,4)**2-P(N+2,1)**2-P(N+2,2)**2-
     &  P(N+2,3)**2))
        Z1=P(IA1,4)/P(N+1,4)
        WT1=(4D0/3D0)*((1D0+Z1**2)/(1D0-Z1))/(P(N+1,5)**2-P(IA1,5)**2)
        Z2=P(IA3,4)/P(N+2,4)
        WT2=(4D0/3D0)*((1D0+Z2**2)/(1D0-Z2))/(P(N+2,5)**2-P(IA3,5)**2)
 
C...Second case: when one original parton radiates to three.
C...IA1  = 0: N+1 -> IA2 + N+2, N+2 -> IA3 + IA4.
      ELSE
        DO 110 J=1,4
          P(N+2,J)=P(IA3,J)+P(IA4,J)
          P(N+1,J)=P(N+2,J)+P(IA2,J)
  110   CONTINUE
        P(N+1,5)=SQRT(MAX(0D0,P(N+1,4)**2-P(N+1,1)**2-P(N+1,2)**2-
     &  P(N+1,3)**2))
        P(N+2,5)=SQRT(MAX(0D0,P(N+2,4)**2-P(N+2,1)**2-P(N+2,2)**2-
     &  P(N+2,3)**2))
        IF(K(IA2,2).EQ.21) THEN
          Z1=P(N+2,4)/P(N+1,4)
          WT1=(4D0/3D0)*((1D0+Z1**2)/(1D0-Z1))/(P(N+1,5)**2-
     &    P(IA3,5)**2)
        ELSE
          Z1=P(IA2,4)/P(N+1,4)
          WT1=(4D0/3D0)*((1D0+Z1**2)/(1D0-Z1))/(P(N+1,5)**2-
     &    P(IA2,5)**2)
        ENDIF
        Z2=P(IA3,4)/P(N+2,4)
        IF(K(IA2,2).EQ.21) THEN
          WT2=(4D0/3D0)*((1D0+Z2**2)/(1D0-Z2))/(P(N+2,5)**2-
     &    P(IA3,5)**2)
        ELSEIF(K(IA3,2).EQ.21) THEN
          WT2=3D0*((1D0-Z2*(1D0-Z2))**2/(Z2*(1D0-Z2)))/P(N+2,5)**2
        ELSE
          WT2=0.5D0*(Z2**2+(1D0-Z2)**2)
        ENDIF
      ENDIF
 
C...Total weight.
      PY4JTW=WT1*WT2
 
      RETURN
      END
