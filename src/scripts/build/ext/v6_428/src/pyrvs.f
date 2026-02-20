 
C*********************************************************************
 
C...PYRVS
C...Interference function
 
      FUNCTION PYRVS(X,Y,M1,W1,M2,W2)
 
      IMPLICIT NONE
      DOUBLE PRECISION X, Y, PYRVS, PYRVR, M1, M2, W1, W2
      PYRVS = PYRVR(X,M1,W1)*PYRVR(Y,M2,W2)*((X-M1**2)*(Y-M2**2)
     &     +W1*W2*M1*M2)
      RETURN
      END
