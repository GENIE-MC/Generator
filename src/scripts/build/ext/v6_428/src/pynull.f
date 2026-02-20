 
C*********************************************************************
 
C...PYNULL
C...Resets bin contents of a histogram.
 
      SUBROUTINE PYNULL(ID)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Commonblock.
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      SAVE /PYBINS/
 
      IF(ID.LE.0.OR.ID.GT.IHIST(1)) RETURN
      IS=INDX(ID)
      IF(IS.EQ.0) RETURN
      DO 100 IX=IS+5,IS+8+NINT(BIN(IS+1))
        BIN(IX)=0D0
  100 CONTINUE
 
      RETURN
      END
