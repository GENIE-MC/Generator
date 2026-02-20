 
C*********************************************************************
 
C...PYGRVV
C...Auxiliary for the GRV 94 parton distribution functions
C...for u and d valence and d-u sea.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      FUNCTION PYGRVV (X, N, AK, BK, A, B, C, D)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION (A - Z)
 
C...Evaluation.
      DX = SQRT (X)
      PYGRVV = N * X**AK * (1D0+ A*X**BK + X * (B + C*DX)) *
     & (1D0- X)**D
 
      RETURN
      END
