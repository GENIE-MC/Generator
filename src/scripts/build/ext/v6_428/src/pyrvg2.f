 
C*********************************************************************
 
C...PYRVG2
C...Integrand for L-R interference contributions
 
      FUNCTION PYRVG2(X)
 
      IMPLICIT NONE
      COMMON/PYRVPM/RM(0:3),A(2),B(2),RESM(2),RESW(2),MFLAG
      DOUBLE PRECISION X, RM, A, B, RESM, RESW, DELTAY, PYRVS
      DOUBLE PRECISION RVS,PYRVG2,E2,E3,C1,SR1,SR2
      LOGICAL MFLAG
      SAVE/PYRVPM/
      C1     = 2D0*SQRT(MAX(0D0,X))
      RVS    = PYRVS(X,X,RESM(1),RESW(1),RESM(2),RESW(2))
      IF (.NOT.MFLAG) THEN
        E2     = X/C1
        E3     = (RM(0)**2-X)/C1
        DELTAY = 4D0*E2*E3
        PYRVG2 = DELTAY*RVS*X*(A(1)*A(2)+B(1)*B(2))*(RM(0)**2-X)
      ELSE
        E2     = (X-RM(1)**2+RM(2)**2)/C1
        E3     = (RM(0)**2-X-RM(3)**2)/C1
        SR1    = SQRT(MAX(0D0,E2**2-RM(2)**2))
        SR2    = SQRT(MAX(0D0,E3**2-RM(3)**2))
        DELTAY = 4D0*SR1*SR2
        PYRVG2 = DELTAY*RVS*(X-RM(1)**2-RM(2)**2)*((A(1)*A(2)
     &       + B(1)*B(2))*(RM(0)**2+RM(3)**2-X)
     &       + 2D0*(A(1)*B(2)+A(2)*B(1))*RM(3)*RM(0))
      ENDIF
      RETURN
      END
