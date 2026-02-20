 
C*********************************************************************
 
C...PYEIG4
C...Finds eigenvalues and eigenvectors to a 4 * 4 matrix.
C...Specific application: mixing in neutralino sector.
 
      SUBROUTINE PYEIG4(A,W,Z)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
 
C...Arrays: in call and local.
      DIMENSION A(4,4),W(4),Z(4,4),X(4),D(4,4),E(4)
 
C...Coefficients of fourth-degree equation from matrix.
C...x**4 + b3 * x**3 + b2 * x**2 + b1 * x + b0 = 0.
      B3=-(A(1,1)+A(2,2)+A(3,3)+A(4,4))
      B2=0D0
      DO 110 I=1,3
        DO 100 J=I+1,4
          B2=B2+A(I,I)*A(J,J)-A(I,J)*A(J,I)
  100   CONTINUE
  110 CONTINUE
      B1=0D0
      B0=0D0
      DO 120 I=1,4
        I1=MOD(I,4)+1
        I2=MOD(I+1,4)+1
        I3=MOD(I+2,4)+1
        B1=B1+A(I,I)*(-A(I1,I1)*A(I2,I2)+A(I1,I2)*A(I2,I1)+
     &  A(I1,I3)*A(I3,I1)+A(I2,I3)*A(I3,I2))-
     &  A(I,I1)*A(I1,I2)*A(I2,I)-A(I,I2)*A(I2,I1)*A(I1,I)
        B0=B0+(-1D0)**(I+1)*A(1,I)*(
     &  A(2,I1)*(A(3,I2)*A(4,I3)-A(3,I3)*A(4,I2))+
     &  A(2,I2)*(A(3,I3)*A(4,I1)-A(3,I1)*A(4,I3))+
     &  A(2,I3)*(A(3,I1)*A(4,I2)-A(3,I2)*A(4,I1)))
  120 CONTINUE
 
C...Coefficients of third-degree equation needed for
C...separation into two second-degree equations.
C...u**3 + c2 * u**2 + c1 * u + c0 = 0.
      C2=-B2
      C1=B1*B3-4D0*B0
      C0=-B1**2-B0*B3**2+4D0*B0*B2
      CQ=C1/3D0-C2**2/9D0
      CR=C1*C2/6D0-C0/2D0-C2**3/27D0
      CQR=CQ**3+CR**2
 
C...Cases with one or three real roots.
      IF(CQR.GE.0D0) THEN
        S1=(CR+SQRT(CQR))**(1D0/3D0)
        S2=(CR-SQRT(CQR))**(1D0/3D0)
        U=S1+S2-C2/3D0
      ELSE
        SABS=SQRT(-CQ)
        THE=ACOS(CR/SABS**3)/3D0
        SRE=SABS*COS(THE)
        U=2D0*SRE-C2/3D0
      ENDIF
 
C...Find and solve two second-degree equations.
      P1=B3/2D0-SQRT(B3**2/4D0+U-B2)
      P2=B3/2D0+SQRT(B3**2/4D0+U-B2)
      Q1=U/2D0+SQRT(U**2/4D0-B0)
      Q2=U/2D0-SQRT(U**2/4D0-B0)
      IF(ABS(P1*Q1+P2*Q2-B1).LT.ABS(P1*Q2+P2*Q1-B1)) THEN
        QSAV=Q1
        Q1=Q2
        Q2=QSAV
      ENDIF
      X(1)=-P1/2D0+SQRT(P1**2/4D0-Q1)
      X(2)=-P1/2D0-SQRT(P1**2/4D0-Q1)
      X(3)=-P2/2D0+SQRT(P2**2/4D0-Q2)
      X(4)=-P2/2D0-SQRT(P2**2/4D0-Q2)
 
C...Order eigenvalues in asceding mass.
      W(1)=X(1)
      DO 150 I1=2,4
        DO 130 I2=I1-1,1,-1
          IF(ABS(X(I1)).GE.ABS(W(I2))) GOTO 140
          W(I2+1)=W(I2)
  130   CONTINUE
  140   W(I2+1)=X(I1)
  150 CONTINUE
 
C...Find equation system for eigenvectors.
      DO 250 I=1,4
        DO 170 J1=1,4
          D(J1,J1)=A(J1,J1)-W(I)
          DO 160 J2=J1+1,4
            D(J1,J2)=A(J1,J2)
            D(J2,J1)=A(J2,J1)
  160     CONTINUE
  170   CONTINUE
 
C...Find largest element in matrix.
        DAMAX=0D0
        DO 190 J1=1,4
          DO 180 J2=1,4
            IF(ABS(D(J1,J2)).LE.DAMAX) GOTO 180
            JA=J1
            JB=J2
            DAMAX=ABS(D(J1,J2))
  180     CONTINUE
  190   CONTINUE
 
C...Subtract others by multiple of row selected above.
        DAMAX=0D0
        DO 210 J3=JA+1,JA+3
          J1=J3-4*((J3-1)/4)
          RL=D(J1,JB)/D(JA,JB)
          DO 200 J2=1,4
            D(J1,J2)=D(J1,J2)-RL*D(JA,J2)
            IF(ABS(D(J1,J2)).LE.DAMAX) GOTO 200
            JC=J1
            JD=J2
            DAMAX=ABS(D(J1,J2))
  200     CONTINUE
  210   CONTINUE
 
C...Do one more subtraction of a row.
        DAMAX=0D0
        DO 230 J3=JC+1,JC+3
          J1=J3-4*((J3-1)/4)
          IF(J1.EQ.JA) GOTO 230
          RL=D(J1,JD)/D(JC,JD)
          DO 220 J2=1,4
            IF(J2.EQ.JB) GOTO 220
            D(J1,J2)=D(J1,J2)-RL*D(JC,J2)
            IF(ABS(D(J1,J2)).LE.DAMAX) GOTO 220
            JE=J1
            DAMAX=ABS(D(J1,J2))
  220     CONTINUE
  230   CONTINUE
 
C...Construct unnormalized eigenvector.
        JF1=JD+1-4*(JD/4)
        JF2=JD+2-4*((JD+1)/4)
        IF(JF1.EQ.JB) JF1=JD+3-4*((JD+2)/4)
        IF(JF2.EQ.JB) JF2=JD+3-4*((JD+2)/4)
        E(JF1)=-D(JE,JF2)
        E(JF2)=D(JE,JF1)
        E(JD)=-(D(JC,JF1)*E(JF1)+D(JC,JF2)*E(JF2))/D(JC,JD)
        E(JB)=-(D(JA,JF1)*E(JF1)+D(JA,JF2)*E(JF2)+D(JA,JD)*E(JD))/
     &  D(JA,JB)
 
C...Normalize and fill in final array.
        EA=SQRT(E(1)**2+E(2)**2+E(3)**2+E(4)**2)
        SGN=(-1D0)**INT(PYR(0)+0.5D0)
        DO 240 J=1,4
          Z(I,J)=SGN*E(J)/EA
  240   CONTINUE
  250 CONTINUE
 
      RETURN
      END
