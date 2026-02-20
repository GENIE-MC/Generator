 
C*********************************************************************
 
C...PYRVG4
C...Integrand for true intereference contributions
 
      FUNCTION PYRVG4(Y)
 
      IMPLICIT NONE
      COMMON/PYRVPM/RM(0:3),A(2),B(2),RESM(2),RESW(2),MFLAG
      COMMON/PYG2DX/X
      DOUBLE PRECISION X, Y, PYRVG4, RM, A, B, RESM, RESW, RVS, PYRVS
      LOGICAL MFLAG
      SAVE /PYRVPM/,/PYG2DX/
      PYRVG4=0D0
      RVS=PYRVS(X,Y,RESM(1),RESW(1),RESM(2),RESW(2))
      IF (.NOT.MFLAG) THEN
        PYRVG4 = RVS*B(1)*B(2)*X*Y
      ELSE
        PYRVG4 = RVS*(RM(1)*RM(3)*A(1)*A(2)*(X+Y-RM(1)**2-RM(3)**2)
     &       + RM(1)*RM(0)*B(1)*A(2)*(Y-RM(2)**2-RM(3)**2)
     &       + RM(3)*RM(0)*A(1)*B(2)*(X-RM(1)**2-RM(2)**2)
     &       + B(1)*B(2)*(X*Y-(RM(1)*RM(3))**2-(RM(0)*RM(2))**2))
      ENDIF
      RETURN
      END
