C=====================================================================
C     ************* FUNCTION SCALAR PRODUCTS *************************
      DOUBLE PRECISION FUNCTION PYTBHS(A,B)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      DIMENSION A(4),B(4)
      DUM=A(4)*B(4)
      DO 100 ID=1,3
         DUM=DUM-A(ID)*B(ID)
  100 CONTINUE
      PYTBHS=DUM
      RETURN
      END
