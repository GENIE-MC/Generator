 
C*********************************************************************
 
C...PYHIST
C...Prints and resets all histograms.
 
      SUBROUTINE PYHIST
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Commonblock.
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      SAVE /PYBINS/
 
C...Loop over histograms, print and reset used ones.
      DO 100 ID=1,IHIST(1)
        IS=INDX(ID)
        IF(IS.NE.0.AND.NINT(BIN(IS+5)).GT.0) THEN
          CALL PYPLOT(ID)
          CALL PYNULL(ID)
        ENDIF
  100 CONTINUE
 
      RETURN
      END
