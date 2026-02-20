 
C*********************************************************************
 
C...PYGDIR
C...Evaluates the direct contribution, i.e. the C^gamma term,
C...as needed in MSbar parametrizations.
C...Adapted from SaSgam library, authors G.A. Schuler and T. Sjostrand.
 
      SUBROUTINE PYGDIR(X,Q2,P2,Q02,XPGA)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Local array and data.
      DIMENSION XPGA(-6:6)
      DATA PMC/1.3D0/, PMB/4.6D0/, AEM2PI/0.0011614D0/
 
C...Reset output.
      DO 100 KFL=-6,6
        XPGA(KFL)=0D0
  100 CONTINUE
 
C...Evaluate common x-dependent expression.
      XTMP = (X**2+(1D0-X)**2) * (-LOG(X)) - 1D0
      CGAM = 3D0*AEM2PI*X * (XTMP*(1D0+P2/(P2+Q02)) + 6D0*X*(1D0-X))
 
C...d, u, s part by simple charge factor.
      XPGA(1)=(1D0/9D0)*CGAM
      XPGA(2)=(4D0/9D0)*CGAM
      XPGA(3)=(1D0/9D0)*CGAM
 
C...Also fill for antiquarks.
      DO 110 KF=1,5
        XPGA(-KF)=XPGA(KF)
  110 CONTINUE
 
      RETURN
      END
