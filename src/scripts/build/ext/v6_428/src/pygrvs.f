 
C*********************************************************************
 
C...PYGRVS
C...Auxiliary for the GRV 94 parton distribution functions
C...for s, c and b sea.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      FUNCTION PYGRVS (X, S, STH, AL, BE, AK, AG, B, D, E, ES)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION (A - Z)
 
C...Evaluation.
      IF(S.LE.STH) THEN
        PYGRVS = 0D0
      ELSE
        DX = SQRT (X)
        LX = LOG (1D0/X)
        PYGRVS = (S - STH)**AL / LX**AK * (1D0+ AG*DX + B*X) *
     &     (1D0- X)**D * EXP (-E + SQRT (ES * S**BE * LX))
      ENDIF
 
      RETURN
      END
