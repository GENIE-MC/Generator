 
C*********************************************************************
 
C...PY4JTS
C...Auxiliary to PY4JET, to set up chosen configuration.
 
      SUBROUTINE PY4JTS(IA1,IA2,IA3,IA4,IA5,QMAX)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      SAVE /PYJETS/
 
C...Reset info.
      DO 110 I=N+1,N+6
        DO 100 J=1,5
          K(I,J)=0
          V(I,J)=V(IA2,J)
  100   CONTINUE
        K(I,1)=16
  110 CONTINUE
 
C...First case: when both original partons radiate.
C...N+1 -> (IA1=N+3) + (IA2=N+4), N+2 -> (IA3=N+5) + (IA4=N+6).
      IF(IA1.NE.0) THEN
 
C...Set up flavour and history pointers for new partons.
        K(N+1,2)=K(IA1,2)
        K(N+2,2)=K(IA3,2)
        K(N+3,2)=K(IA1,2)
        K(N+4,2)=K(IA2,2)
        K(N+5,2)=K(IA3,2)
        K(N+6,2)=K(IA4,2)
        K(N+1,3)=IA1
        K(N+1,4)=N+3
        K(N+1,5)=N+4
        K(N+2,3)=IA3
        K(N+2,4)=N+5
        K(N+2,5)=N+6
        K(N+3,3)=N+1
        K(N+4,3)=N+1
        K(N+5,3)=N+2
        K(N+6,3)=N+2
 
C...Set up momenta for new partons.
        DO 120 J=1,5
          P(N+1,J)=P(IA1,J)+P(IA2,J)
          P(N+2,J)=P(IA3,J)+P(IA4,J)
          P(N+3,J)=P(IA1,J)
          P(N+4,J)=P(IA2,J)
          P(N+5,J)=P(IA3,J)
          P(N+6,J)=P(IA4,J)
  120   CONTINUE
        P(N+1,5)=SQRT(MAX(0D0,P(N+1,4)**2-P(N+1,1)**2-P(N+1,2)**2-
     &  P(N+1,3)**2))
        P(N+2,5)=SQRT(MAX(0D0,P(N+2,4)**2-P(N+2,1)**2-P(N+2,2)**2-
     &  P(N+2,3)**2))
        QMAX=MIN(P(N+1,5),P(N+2,5))
 
C...Second case: q radiates twice.
C...N+1 -> (IA2=N+4) + N+3, N+3 -> (IA3=N+5) + (IA4=N+6),
C...IA5=N+2 does not radiate.
      ELSEIF(K(IA2,2).EQ.21) THEN
 
C...Set up flavour and history pointers for new partons.
        K(N+1,2)=K(IA3,2)
        K(N+2,2)=K(IA5,2)
        K(N+3,2)=K(IA3,2)
        K(N+4,2)=K(IA2,2)
        K(N+5,2)=K(IA3,2)
        K(N+6,2)=K(IA4,2)
        K(N+1,3)=IA3
        K(N+1,4)=N+3
        K(N+1,5)=N+4
        K(N+2,3)=IA5
        K(N+3,3)=N+1
        K(N+3,4)=N+5
        K(N+3,5)=N+6
        K(N+4,3)=N+1
        K(N+5,3)=N+3
        K(N+6,3)=N+3
 
C...Set up momenta for new partons.
        DO 130 J=1,5
          P(N+1,J)=P(IA2,J)+P(IA3,J)+P(IA4,J)
          P(N+2,J)=P(IA5,J)
          P(N+3,J)=P(IA3,J)+P(IA4,J)
          P(N+4,J)=P(IA2,J)
          P(N+5,J)=P(IA3,J)
          P(N+6,J)=P(IA4,J)
  130   CONTINUE
        P(N+1,5)=SQRT(MAX(0D0,P(N+1,4)**2-P(N+1,1)**2-P(N+1,2)**2-
     &  P(N+1,3)**2))
        P(N+3,5)=SQRT(MAX(0D0,P(N+3,4)**2-P(N+3,1)**2-P(N+3,2)**2-
     &  P(N+3,3)**2))
        QMAX=P(N+3,5)
 
C...Third case: q radiates g, g branches.
C...N+1 -> (IA2=N+3) + N+4, N+4 -> (IA3=N+5) + (IA4=N+6),
C...IA5=N+2 does not radiate.
      ELSE
 
C...Set up flavour and history pointers for new partons.
        K(N+1,2)=K(IA2,2)
        K(N+2,2)=K(IA5,2)
        K(N+3,2)=K(IA2,2)
        K(N+4,2)=21
        K(N+5,2)=K(IA3,2)
        K(N+6,2)=K(IA4,2)
        K(N+1,3)=IA2
        K(N+1,4)=N+3
        K(N+1,5)=N+4
        K(N+2,3)=IA5
        K(N+3,3)=N+1
        K(N+4,3)=N+1
        K(N+4,4)=N+5
        K(N+4,5)=N+6
        K(N+5,3)=N+4
        K(N+6,3)=N+4
 
C...Set up momenta for new partons.
        DO 140 J=1,5
          P(N+1,J)=P(IA2,J)+P(IA3,J)+P(IA4,J)
          P(N+2,J)=P(IA5,J)
          P(N+3,J)=P(IA2,J)
          P(N+4,J)=P(IA3,J)+P(IA4,J)
          P(N+5,J)=P(IA3,J)
          P(N+6,J)=P(IA4,J)
  140   CONTINUE
        P(N+1,5)=SQRT(MAX(0D0,P(N+1,4)**2-P(N+1,1)**2-P(N+1,2)**2-
     &  P(N+1,3)**2))
        P(N+4,5)=SQRT(MAX(0D0,P(N+4,4)**2-P(N+4,1)**2-P(N+4,2)**2-
     &  P(N+4,3)**2))
        QMAX=P(N+4,5)
 
      ENDIF
      N=N+6
 
      RETURN
      END
