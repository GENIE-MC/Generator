 
C*********************************************************************
 
C...PYBKSB
C...Auxiliary to PYSIGH, for technicolor corrections to QCD 2 -> 2
C...processes.
 
      SUBROUTINE PYBKSB(A,N,NP,INDX,B)
      IMPLICIT NONE
      INTEGER N,NP,INDX(N)
      COMPLEX*16 A(NP,NP),B(N)
      INTEGER I,II,J,LL
      COMPLEX*16 SUM
 
      II=0
      DO 110 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 100 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
  100     CONTINUE
        ELSE IF (ABS(SUM).NE.0D0) THEN
          II=I
        ENDIF
        B(I)=SUM
  110 CONTINUE
      DO 130 I=N,1,-1
        SUM=B(I)
        DO 120 J=I+1,N
          SUM=SUM-A(I,J)*B(J)
  120   CONTINUE
        B(I)=SUM/A(I,I)
  130 CONTINUE
      RETURN
      END
