 
C*********************************************************************
 
C...PYRVG1
C...Integrand for resonance contributions
 
      FUNCTION PYRVG1(X)
 
      IMPLICIT NONE
      COMMON/PYRVPM/RM(0:3),A(2),B(2),RESM(2),RESW(2),MFLAG
      DOUBLE PRECISION X, RM, A, B, RESM, RESW, DELTAY,PYRVR
      DOUBLE PRECISION RVR,PYRVG1,E2,E3,C1,SR1,SR2,A1,A2
      LOGICAL MFLAG
      SAVE/PYRVPM/
      RVR    = PYRVR(X,RESM(1),RESW(1))
      C1     = 2D0*SQRT(MAX(0D0,X))
      IF (.NOT.MFLAG) THEN
        E2     = X/C1
        E3     = (RM(0)**2-X)/C1
        DELTAY = 4D0*E2*E3
        PYRVG1 = DELTAY*RVR*X*(A(1)**2+B(1)**2)*(RM(0)**2-X)
      ELSE
        E2     = (X-RM(1)**2+RM(2)**2)/C1
        E3     = (RM(0)**2-X-RM(3)**2)/C1
        SR1    = SQRT(MAX(0D0,E2**2-RM(2)**2))
        SR2    = SQRT(MAX(0D0,E3**2-RM(3)**2))
        DELTAY = 4D0*SR1*SR2
        A1     = 4.*A(1)*B(1)*RM(3)*RM(0)
        A2     = (A(1)**2+B(1)**2)*(RM(0)**2+RM(3)**2-X)
        PYRVG1 = DELTAY*RVR*(X-RM(1)**2-RM(2)**2)*(A1+A2)
      ENDIF
      RETURN
      END
