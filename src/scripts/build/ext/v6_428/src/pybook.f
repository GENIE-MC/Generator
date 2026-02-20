 
C*********************************************************************
 
C...PYBOOK
C...Books a histogram.
 
      SUBROUTINE PYBOOK(ID,TITLE,NX,XL,XU)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Commonblock.
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      SAVE /PYBINS/
C...Local character variables.
      CHARACTER TITLE*(*), TITFX*60
 
C...Check that input is sensible. Find initial address in memory.
      IF(ID.LE.0.OR.ID.GT.IHIST(1)) CALL PYERRM(28,
     &'(PYBOOK:) not allowed histogram number')
      IF(NX.LE.0.OR.NX.GT.100) CALL PYERRM(28,
     &'(PYBOOK:) not allowed number of bins')
      IF(XL.GE.XU) CALL PYERRM(28,
     &'(PYBOOK:) x limits in wrong order')
      INDX(ID)=IHIST(4)
      IHIST(4)=IHIST(4)+28+NX
      IF(IHIST(4).GT.IHIST(2)) CALL PYERRM(28,
     &'(PYBOOK:) out of histogram space')
      IS=INDX(ID)
 
C...Store histogram size and reset contents.
      BIN(IS+1)=NX
      BIN(IS+2)=XL
      BIN(IS+3)=XU
      BIN(IS+4)=(XU-XL)/NX
      CALL PYNULL(ID)
 
C...Store title by conversion to integer to double precision.
      TITFX=TITLE//' '
      DO 100 IT=1,20
        BIN(IS+8+NX+IT)=256**2*ICHAR(TITFX(3*IT-2:3*IT-2))+
     &  256*ICHAR(TITFX(3*IT-1:3*IT-1))+ICHAR(TITFX(3*IT:3*IT))
  100 CONTINUE
 
      RETURN
      END
