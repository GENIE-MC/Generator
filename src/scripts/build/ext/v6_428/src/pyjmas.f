 
C*********************************************************************
 
C...PYJMAS
C...Determines, approximately, the two jet masses that minimize
C...the sum m_H^2 + m_L^2, a la Clavelli and Wyler.
 
      SUBROUTINE PYJMAS(PMH,PML)
 
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
C...Local arrays.
      DIMENSION SM(3,3),SAX(3),PS(3,5)
 
C...Reset.
      NP=0
      DO 120 J1=1,3
        DO 100 J2=J1,3
          SM(J1,J2)=0D0
  100   CONTINUE
        DO 110 J2=1,4
          PS(J1,J2)=0D0
  110   CONTINUE
  120 CONTINUE
      PSS=0D0
      PIMASS=PMAS(PYCOMP(211),1)
 
C...Take copy of particles that are to be considered in mass analysis.
      DO 170 I=1,N
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 170
        IF(MSTU(41).GE.2) THEN
          KC=PYCOMP(K(I,2))
          IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.
     &    KC.EQ.18.OR.K(I,2).EQ.KSUSY1+22.OR.K(I,2).EQ.39.OR.
     &    K(I,2).EQ.KSUSY1+39) GOTO 170
          IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.PYCHGE(K(I,2)).EQ.0)
     &    GOTO 170
        ENDIF
        IF(N+NP+1.GE.MSTU(4)-MSTU(32)-5) THEN
          CALL PYERRM(11,'(PYJMAS:) no more memory left in PYJETS')
          PMH=-2D0
          PML=-2D0
          RETURN
        ENDIF
        NP=NP+1
        DO 130 J=1,5
          P(N+NP,J)=P(I,J)
  130   CONTINUE
        IF(MSTU(42).EQ.0) P(N+NP,5)=0D0
        IF(MSTU(42).EQ.1.AND.K(I,2).NE.22) P(N+NP,5)=PIMASS
        P(N+NP,4)=SQRT(P(N+NP,5)**2+P(I,1)**2+P(I,2)**2+P(I,3)**2)
 
C...Fill information in sphericity tensor and total momentum vector.
        DO 150 J1=1,3
          DO 140 J2=J1,3
            SM(J1,J2)=SM(J1,J2)+P(I,J1)*P(I,J2)
  140     CONTINUE
  150   CONTINUE
        PSS=PSS+(P(I,1)**2+P(I,2)**2+P(I,3)**2)
        DO 160 J=1,4
          PS(3,J)=PS(3,J)+P(N+NP,J)
  160   CONTINUE
  170 CONTINUE
 
C...Very low multiplicities (0 or 1) not considered.
      IF(NP.LE.1) THEN
        CALL PYERRM(8,'(PYJMAS:) too few particles for analysis')
        PMH=-1D0
        PML=-1D0
        RETURN
      ENDIF
      PARU(61)=SQRT(MAX(0D0,PS(3,4)**2-PS(3,1)**2-PS(3,2)**2-
     &PS(3,3)**2))
 
C...Find largest eigenvalue to matrix (third degree equation).
      DO 190 J1=1,3
        DO 180 J2=J1,3
          SM(J1,J2)=SM(J1,J2)/PSS
  180   CONTINUE
  190 CONTINUE
      SQ=(SM(1,1)*SM(2,2)+SM(1,1)*SM(3,3)+SM(2,2)*SM(3,3)-
     &SM(1,2)**2-SM(1,3)**2-SM(2,3)**2)/3D0-1D0/9D0
      SR=-0.5D0*(SQ+1D0/9D0+SM(1,1)*SM(2,3)**2+SM(2,2)*SM(1,3)**2+
     &SM(3,3)*SM(1,2)**2-SM(1,1)*SM(2,2)*SM(3,3))+
     &SM(1,2)*SM(1,3)*SM(2,3)+1D0/27D0
      SP=COS(ACOS(MAX(MIN(SR/SQRT(-SQ**3),1D0),-1D0))/3D0)
      SMA=1D0/3D0+SQRT(-SQ)*MAX(2D0*SP,SQRT(3D0*(1D0-SP**2))-SP)
 
C...Find largest eigenvector by solving equation system.
      DO 210 J1=1,3
        SM(J1,J1)=SM(J1,J1)-SMA
        DO 200 J2=J1+1,3
          SM(J2,J1)=SM(J1,J2)
  200   CONTINUE
  210 CONTINUE
      SMAX=0D0
      DO 230 J1=1,3
        DO 220 J2=1,3
          IF(ABS(SM(J1,J2)).LE.SMAX) GOTO 220
          JA=J1
          JB=J2
          SMAX=ABS(SM(J1,J2))
  220   CONTINUE
  230 CONTINUE
      SMAX=0D0
      DO 250 J3=JA+1,JA+2
        J1=J3-3*((J3-1)/3)
        RL=SM(J1,JB)/SM(JA,JB)
        DO 240 J2=1,3
          SM(J1,J2)=SM(J1,J2)-RL*SM(JA,J2)
          IF(ABS(SM(J1,J2)).LE.SMAX) GOTO 240
          JC=J1
          SMAX=ABS(SM(J1,J2))
  240   CONTINUE
  250 CONTINUE
      JB1=JB+1-3*(JB/3)
      JB2=JB+2-3*((JB+1)/3)
      SAX(JB1)=-SM(JC,JB2)
      SAX(JB2)=SM(JC,JB1)
      SAX(JB)=-(SM(JA,JB1)*SAX(JB1)+SM(JA,JB2)*SAX(JB2))/SM(JA,JB)
 
C...Divide particles into two initial clusters by hemisphere.
      DO 270 I=N+1,N+NP
        PSAX=P(I,1)*SAX(1)+P(I,2)*SAX(2)+P(I,3)*SAX(3)
        IS=1
        IF(PSAX.LT.0D0) IS=2
        K(I,3)=IS
        DO 260 J=1,4
          PS(IS,J)=PS(IS,J)+P(I,J)
  260   CONTINUE
  270 CONTINUE
      PMS=MAX(1D-10,PS(1,4)**2-PS(1,1)**2-PS(1,2)**2-PS(1,3)**2)+
     &MAX(1D-10,PS(2,4)**2-PS(2,1)**2-PS(2,2)**2-PS(2,3)**2)
 
C...Reassign one particle at a time; find maximum decrease of m^2 sum.
  280 PMD=0D0
      IM=0
      DO 290 J=1,4
        PS(3,J)=PS(1,J)-PS(2,J)
  290 CONTINUE
      DO 300 I=N+1,N+NP
        PPS=P(I,4)*PS(3,4)-P(I,1)*PS(3,1)-P(I,2)*PS(3,2)-P(I,3)*PS(3,3)
        IF(K(I,3).EQ.1) PMDI=2D0*(P(I,5)**2-PPS)
        IF(K(I,3).EQ.2) PMDI=2D0*(P(I,5)**2+PPS)
        IF(PMDI.LT.PMD) THEN
          PMD=PMDI
          IM=I
        ENDIF
  300 CONTINUE
 
C...Loop back if significant reduction in sum of m^2.
      IF(PMD.LT.-PARU(48)*PMS) THEN
        PMS=PMS+PMD
        IS=K(IM,3)
        DO 310 J=1,4
          PS(IS,J)=PS(IS,J)-P(IM,J)
          PS(3-IS,J)=PS(3-IS,J)+P(IM,J)
  310   CONTINUE
        K(IM,3)=3-IS
        GOTO 280
      ENDIF
 
C...Final masses and output.
      MSTU(61)=N+1
      MSTU(62)=NP
      PS(1,5)=SQRT(MAX(0D0,PS(1,4)**2-PS(1,1)**2-PS(1,2)**2-PS(1,3)**2))
      PS(2,5)=SQRT(MAX(0D0,PS(2,4)**2-PS(2,1)**2-PS(2,2)**2-PS(2,3)**2))
      PMH=MAX(PS(1,5),PS(2,5))
      PML=MIN(PS(1,5),PS(2,5))
 
      RETURN
      END
