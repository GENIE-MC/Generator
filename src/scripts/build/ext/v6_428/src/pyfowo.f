 
C*********************************************************************
 
C...PYFOWO
C...Calculates the first few Fox-Wolfram moments.
 
      SUBROUTINE PYFOWO(H10,H20,H30,H40)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
 
C...Copy momenta for particles and calculate H0.
      NP=0
      H0=0D0
      HD=0D0
      DO 110 I=1,N
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 110
        IF(MSTU(41).GE.2) THEN
          KC=PYCOMP(K(I,2))
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.
     &    KC.EQ.18.OR.K(I,2).EQ.KSUSY1+22.OR.K(I,2).EQ.39.OR.
     &    K(I,2).EQ.KSUSY1+39) GOTO 110
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.PYCHGE(K(I,2)).EQ.0)
     &    GOTO 110
        ENDIF
        IF(N+NP.GE.MSTU(4)-MSTU(32)-5) THEN
          CALL PYERRM(11,'(PYFOWO:) no more memory left in PYJETS')
          H10=-1D0
          H20=-1D0
          H30=-1D0
          H40=-1D0
          RETURN
        ENDIF
        NP=NP+1
        DO 100 J=1,3
          P(N+NP,J)=P(I,J)
  100   CONTINUE
        P(N+NP,4)=SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2)
        H0=H0+P(N+NP,4)
        HD=HD+P(N+NP,4)**2
  110 CONTINUE
      H0=H0**2
 
C...Very low multiplicities (0 or 1) not considered.
      IF(NP.LE.1) THEN
        CALL PYERRM(8,'(PYFOWO:) too few particles for analysis')
        H10=-1D0
        H20=-1D0
        H30=-1D0
        H40=-1D0
        RETURN
      ENDIF
 
C...Calculate H1 - H4.
      H10=0D0
      H20=0D0
      H30=0D0
      H40=0D0
      DO 130 I1=N+1,N+NP
        DO 120 I2=I1+1,N+NP
          CTHE=(P(I1,1)*P(I2,1)+P(I1,2)*P(I2,2)+P(I1,3)*P(I2,3))/
     &    (P(I1,4)*P(I2,4))
          H10=H10+P(I1,4)*P(I2,4)*CTHE
          H20=H20+P(I1,4)*P(I2,4)*(1.5D0*CTHE**2-0.5D0)
          H30=H30+P(I1,4)*P(I2,4)*(2.5D0*CTHE**3-1.5D0*CTHE)
          H40=H40+P(I1,4)*P(I2,4)*(4.375D0*CTHE**4-3.75D0*CTHE**2+
     &    0.375D0)
  120   CONTINUE
  130 CONTINUE
 
C...Calculate H1/H0 - H4/H0. Output.
      MSTU(61)=N+1
      MSTU(62)=NP
      H10=(HD+2D0*H10)/H0
      H20=(HD+2D0*H20)/H0
      H30=(HD+2D0*H30)/H0
      H40=(HD+2D0*H40)/H0
 
      RETURN
      END
