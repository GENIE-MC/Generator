 
C*********************************************************************
 
C...PYALEM
C...Calculates the running alpha_electromagnetic.
 
      FUNCTION PYALEM(Q2)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
C...Calculate real part of photon vacuum polarization.
C...For leptons simplify by using asymptotic (Q^2 >> m^2) expressions.
C...For hadrons use parametrization of H. Burkhardt et al.
C...See R. Kleiss et al, CERN 89-08, vol. 3, pp. 129-131.
      AEMPI=PARU(101)/(3D0*PARU(1))
      IF(MSTU(101).LE.0.OR.Q2.LT.2D-6) THEN
        RPIGG=0D0
      ELSEIF(MSTU(101).EQ.2.AND.Q2.LT.PARU(104)) THEN
        RPIGG=0D0
      ELSEIF(MSTU(101).EQ.2) THEN
        RPIGG=1D0-PARU(101)/PARU(103)
      ELSEIF(Q2.LT.0.09D0) THEN
        RPIGG=AEMPI*(13.4916D0+LOG(Q2))+0.00835D0*LOG(1D0+Q2)
      ELSEIF(Q2.LT.9D0) THEN
        RPIGG=AEMPI*(16.3200D0+2D0*LOG(Q2))+
     &  0.00238D0*LOG(1D0+3.927D0*Q2)
      ELSEIF(Q2.LT.1D4) THEN
        RPIGG=AEMPI*(13.4955D0+3D0*LOG(Q2))+0.00165D0+
     &  0.00299D0*LOG(1D0+Q2)
      ELSE
        RPIGG=AEMPI*(13.4955D0+3D0*LOG(Q2))+0.00221D0+
     &  0.00293D0*LOG(1D0+Q2)
      ENDIF
 
C...Calculate running alpha_em.
      PYALEM=PARU(101)/(1D0-RPIGG)
      PARU(108)=PYALEM
 
      RETURN
      END
