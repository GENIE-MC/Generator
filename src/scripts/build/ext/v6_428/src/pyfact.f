 
C*********************************************************************
 
C...PYFACT
C...Multiplies histogram contents by factor.
 
      SUBROUTINE PYFACT(ID,F)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Commonblock.
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      SAVE /PYBINS/
 
C...Find initial address in memory. Multiply all contents bins.
      IF(ID.LE.0.OR.ID.GT.IHIST(1)) CALL PYERRM(28,
     &'(PYFACT:) not allowed histogram number')
      IS=INDX(ID)
      IF(IS.EQ.0) CALL PYERRM(28,
     &'(PYFACT:) scaling unbooked histogram')
      DO 100 IX=IS+6,IS+8+NINT(BIN(IS+1))
        BIN(IX)=F*BIN(IX)
  100 CONTINUE
 
      RETURN
      END
