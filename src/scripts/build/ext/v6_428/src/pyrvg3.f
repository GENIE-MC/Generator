 
C*********************************************************************
 
C...PYRVG3
C...Function to do Y integration over true interference contributions
 
      FUNCTION PYRVG3(X)
 
      IMPLICIT NONE
      COMMON/PYRVPM/RM(0:3),A(2),B(2),RESM(2),RESW(2),MFLAG
C...Second Dalitz variable for PYRVG4
      COMMON/PYG2DX/X1
      DOUBLE PRECISION RM, A, B, RESM, RESW, X, X1
      DOUBLE PRECISION E2, E3, C1, SQ1, SR1, SR2, YMIN, YMAX
      DOUBLE PRECISION PYRVG3, PYRVG4, PYGAU2
      LOGICAL MFLAG
      EXTERNAL PYGAU2,PYRVG4
      SAVE/PYRVPM/,/PYG2DX/
      PYRVG3=0D0
      C1=2D0*SQRT(MAX(1D-9,X))
      X1=X
      IF (.NOT.MFLAG) THEN
        E2    = X/C1
        E3    = (RM(0)**2-X)/C1
        YMIN  = 0D0
        YMAX  = 4D0*E2*E3
      ELSE
        E2    = (X-RM(1)**2+RM(2)**2)/C1
        E3    = (RM(0)**2-X-RM(3)**2)/C1
        SQ1   = (E2+E3)**2
        SR1   = SQRT(MAX(0D0,E2**2-RM(2)**2))
        SR2   = SQRT(MAX(0D0,E3**2-RM(3)**2))
        YMIN  = SQ1-(SR1+SR2)**2
        YMAX  = SQ1-(SR1-SR2)**2
      ENDIF
      PYRVG3 = PYGAU2(PYRVG4,YMIN,YMAX,1D-3)
      RETURN
      END
