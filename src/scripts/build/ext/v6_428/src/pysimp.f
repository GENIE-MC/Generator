 
C*********************************************************************
 
C...PYSIMP
C...Simpson formula for an integral.
 
      FUNCTION PYSIMP(Y,X0,X1,N)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
 
C...Local variables.
      DOUBLE PRECISION Y,X0,X1,H,S
      DIMENSION Y(0:N)
 
      S=0D0
      H=(X1-X0)/N
      DO 100 I=0,N-2,2
        S=S+Y(I)+4D0*Y(I+1)+Y(I+2)
  100 CONTINUE
      PYSIMP=S*H/3D0
 
      RETURN
      END
