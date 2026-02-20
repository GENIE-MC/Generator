 
C*********************************************************************
 
C...PYLDCM
C...Auxiliary to PYSIGH, for technicolor corrections to QCD 2 -> 2
C...processes.
 
      SUBROUTINE PYLDCM(A,N,NP,INDX,D)
      IMPLICIT NONE
      INTEGER N,NP,INDX(N)
      REAL*8 D,TINY
      COMPLEX*16 A(NP,NP)
      PARAMETER (TINY=1.0D-20)
      INTEGER I,IMAX,J,K
      REAL*8 AAMAX,VV(6),DUM
      COMPLEX*16 SUM,DUMC
 
      D=1D0
      DO 110 I=1,N
        AAMAX=0D0
        DO 100 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
  100   CONTINUE
        IF (AAMAX.EQ.0D0) CALL PYERRM(28,'(PYLDCM:) singular matrix')
        VV(I)=1D0/AAMAX
  110 CONTINUE
      DO 180 J=1,N
        DO 130 I=1,J-1
          SUM=A(I,J)
          DO 120 K=1,I-1
            SUM=SUM-A(I,K)*A(K,J)
  120     CONTINUE
          A(I,J)=SUM
  130   CONTINUE
        AAMAX=0D0
        DO 150 I=J,N
          SUM=A(I,J)
          DO 140 K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
  140     CONTINUE
          A(I,J)=SUM
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
  150   CONTINUE
        IF (J.NE.IMAX)THEN
          DO 160 K=1,N
            DUMC=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUMC
  160     CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(ABS(A(J,J)).EQ.0D0) A(J,J)=DCMPLX(TINY,0D0)
        IF(J.NE.N)THEN
          DO 170 I=J+1,N
            A(I,J)=A(I,J)/A(J,J)
  170     CONTINUE
        ENDIF
  180 CONTINUE
 
      RETURN
      END
