 
 
C*********************************************************************
 
C...PYXXGA
C...Calculates chi0_i -> chi0_j + gamma.
 
      FUNCTION PYXXGA(C0,XM1,XM2,XMTR,XMTL)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
 
C...Local variables.
      DOUBLE PRECISION PYXXGA,C0,XM1,XM2,XMTR,XMTL
      DOUBLE PRECISION F1,F2
 
      F1=(1D0+XMTR/(1D0-XMTR)*LOG(XMTR))/(1D0-XMTR)
      F2=(1D0+XMTL/(1D0-XMTL)*LOG(XMTL))/(1D0-XMTL)
      PYXXGA=C0*((XM1**2-XM2**2)/XM1)**3
      PYXXGA=PYXXGA*(2D0/3D0*(F1+F2)-13D0/12D0)**2
 
      RETURN
      END
