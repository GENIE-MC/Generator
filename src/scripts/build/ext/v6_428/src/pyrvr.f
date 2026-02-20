 
C*********************************************************************
 
C...PYRVR
C...Breit-Wigner for resonance contributions
 
      FUNCTION PYRVR(Mab2,RM,RW)
 
      IMPLICIT NONE
      DOUBLE PRECISION Mab2,RM,RW,PYRVR
      PYRVR = 1D0/((Mab2-RM**2)**2+RM**2*RW**2)
      RETURN
      END
