 
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
C
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GOTO 110
      R = (DMIN1(DABS(A),DABS(B))/P)**2
  100 CONTINUE
         T = 4.0D0 + R
         IF (T .EQ. 4.0D0) GOTO 110
         S = R/T
         U = 1.0D0 + 2.0D0*S
         P = U*P
         R = (S/U)**2 * R
      GOTO 100
  110 PYTHAG = P
      RETURN
      END
