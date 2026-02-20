 
C*********************************************************************
 
C...PYCSRT
C...Auxiliary to PYCMQR
C
C     (YR,YI) = COMPLEX DSQRT(XR,XI)
C     BRANCH CHOSEN SO THAT YR .GE. 0.0 AND SIGN(YI) .EQ. SIGN(XI)
C
 
      SUBROUTINE PYCSRT(XR,XI,YR,YI)
 
      DOUBLE PRECISION XR,XI,YR,YI
      DOUBLE PRECISION S,TR,TI,PYTHAG
 
      TR = XR
      TI = XI
      S = DSQRT(0.5D0*(PYTHAG(TR,TI) + DABS(TR)))
      IF (TR .GE. 0.0D0) YR = S
      IF (TI .LT. 0.0D0) S = -S
      IF (TR .LE. 0.0D0) YI = S
      IF (TR .LT. 0.0D0) YR = 0.5D0*(TI/YI)
      IF (TR .GT. 0.0D0) YI = 0.5D0*(TI/YR)
      RETURN
      END
