 
C*********************************************************************
 
C...PY4ENT
C...Stores four partons or particles in their CM frame, with
C...the first along the +z axis, the last in the xz plane with x > 0
C...and the second having y < 0 and y > 0 with equal probability.
 
      SUBROUTINE PY4ENT(IP,KF1,KF2,KF3,KF4,PECM,X1,X2,X4,X12,X14)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
 
C...Standard checks.
      MSTU(28)=0
      IF(MSTU(12).NE.12345) CALL PYLIST(0)
      IPA=MAX(1,IABS(IP))
      IF(IPA.GT.MSTU(4)-3) CALL PYERRM(21,
     &'(PY4ENT:) writing outside PYJETS momory')
      KC1=PYCOMP(KF1)
      KC2=PYCOMP(KF2)
      KC3=PYCOMP(KF3)
      KC4=PYCOMP(KF4)
      IF(KC1.EQ.0.OR.KC2.EQ.0.OR.KC3.EQ.0.OR.KC4.EQ.0) CALL PYERRM(12,
     &'(PY4ENT:) unknown flavour code')
 
C...Find masses. Reset K, P and V vectors.
      PM1=0D0
      IF(MSTU(10).EQ.1) PM1=P(IPA,5)
      IF(MSTU(10).GE.2) PM1=PYMASS(KF1)
      PM2=0D0
      IF(MSTU(10).EQ.1) PM2=P(IPA+1,5)
      IF(MSTU(10).GE.2) PM2=PYMASS(KF2)
      PM3=0D0
      IF(MSTU(10).EQ.1) PM3=P(IPA+2,5)
      IF(MSTU(10).GE.2) PM3=PYMASS(KF3)
      PM4=0D0
      IF(MSTU(10).EQ.1) PM4=P(IPA+3,5)
      IF(MSTU(10).GE.2) PM4=PYMASS(KF4)
      DO 110 I=IPA,IPA+3
        DO 100 J=1,5
          K(I,J)=0
          P(I,J)=0D0
          V(I,J)=0D0
  100   CONTINUE
  110 CONTINUE
 
C...Check flavours.
      KQ1=KCHG(KC1,2)*ISIGN(1,KF1)
      KQ2=KCHG(KC2,2)*ISIGN(1,KF2)
      KQ3=KCHG(KC3,2)*ISIGN(1,KF3)
      KQ4=KCHG(KC4,2)*ISIGN(1,KF4)
      IF(MSTU(19).EQ.1) THEN
        MSTU(19)=0
      ELSEIF(KQ1.EQ.0.AND.KQ2.EQ.0.AND.KQ3.EQ.0.AND.KQ4.EQ.0) THEN
      ELSEIF(KQ1.NE.0.AND.KQ2.EQ.2.AND.KQ3.EQ.2.AND.(KQ1+KQ4.EQ.0.OR.
     &  KQ1+KQ4.EQ.4)) THEN
      ELSEIF(KQ1.NE.0.AND.KQ1+KQ2.EQ.0.AND.KQ3.NE.0.AND.KQ3+KQ4.EQ.0D0)
     &  THEN
      ELSE
        CALL PYERRM(2,'(PY4ENT:) unphysical flavour combination')
      ENDIF
      K(IPA,2)=KF1
      K(IPA+1,2)=KF2
      K(IPA+2,2)=KF3
      K(IPA+3,2)=KF4
 
C...Store partons/particles in K vectors for normal case.
      IF(IP.GE.0) THEN
        K(IPA,1)=1
        IF(KQ1.NE.0.AND.(KQ2.NE.0.OR.KQ3.NE.0.OR.KQ4.NE.0)) K(IPA,1)=2
        K(IPA+1,1)=1
        IF(KQ2.NE.0.AND.KQ1+KQ2.NE.0.AND.(KQ3.NE.0.OR.KQ4.NE.0))
     &  K(IPA+1,1)=2
        K(IPA+2,1)=1
        IF(KQ3.NE.0.AND.KQ4.NE.0) K(IPA+2,1)=2
        K(IPA+3,1)=1
 
C...Store partons for parton shower evolution from q-g-g-qbar or
C...g-g-g-g event.
      ELSEIF(KQ1+KQ2.NE.0) THEN
        K(IPA,1)=3
        K(IPA+1,1)=3
        K(IPA+2,1)=3
        K(IPA+3,1)=3
        KCS=4
        IF(KQ1.EQ.-1) KCS=5
        K(IPA,KCS)=MSTU(5)*(IPA+1)
        K(IPA,9-KCS)=MSTU(5)*(IPA+3)
        K(IPA+1,KCS)=MSTU(5)*(IPA+2)
        K(IPA+1,9-KCS)=MSTU(5)*IPA
        K(IPA+2,KCS)=MSTU(5)*(IPA+3)
        K(IPA+2,9-KCS)=MSTU(5)*(IPA+1)
        K(IPA+3,KCS)=MSTU(5)*IPA
        K(IPA+3,9-KCS)=MSTU(5)*(IPA+2)
 
C...Store partons for parton shower evolution from q-qbar-q-qbar event.
      ELSE
        K(IPA,1)=3
        K(IPA+1,1)=3
        K(IPA+2,1)=3
        K(IPA+3,1)=3
        K(IPA,4)=MSTU(5)*(IPA+1)
        K(IPA,5)=K(IPA,4)
        K(IPA+1,4)=MSTU(5)*IPA
        K(IPA+1,5)=K(IPA+1,4)
        K(IPA+2,4)=MSTU(5)*(IPA+3)
        K(IPA+2,5)=K(IPA+2,4)
        K(IPA+3,4)=MSTU(5)*(IPA+2)
        K(IPA+3,5)=K(IPA+3,4)
      ENDIF
 
C...Check kinematics.
      MKERR=0
      IF(0.5D0*X1*PECM.LE.PM1.OR.0.5D0*X2*PECM.LE.PM2.OR.
     &0.5D0*(2D0-X1-X2-X4)*PECM.LE.PM3.OR.0.5D0*X4*PECM.LE.PM4)
     &MKERR=1
      PA1=SQRT(MAX(1D-10,(0.5D0*X1*PECM)**2-PM1**2))
      PA2=SQRT(MAX(1D-10,(0.5D0*X2*PECM)**2-PM2**2))
      PA4=SQRT(MAX(1D-10,(0.5D0*X4*PECM)**2-PM4**2))
      X24=X1+X2+X4-1D0-X12-X14+(PM3**2-PM1**2-PM2**2-PM4**2)/PECM**2
      CTHE4=(X1*X4-2D0*X14)*PECM**2/(4D0*PA1*PA4)
      IF(ABS(CTHE4).GE.1.002D0) MKERR=1
      CTHE4=MAX(-1D0,MIN(1D0,CTHE4))
      STHE4=SQRT(1D0-CTHE4**2)
      CTHE2=(X1*X2-2D0*X12)*PECM**2/(4D0*PA1*PA2)
      IF(ABS(CTHE2).GE.1.002D0) MKERR=1
      CTHE2=MAX(-1D0,MIN(1D0,CTHE2))
      STHE2=SQRT(1D0-CTHE2**2)
      CPHI2=((X2*X4-2D0*X24)*PECM**2-4D0*PA2*CTHE2*PA4*CTHE4)/
     &MAX(1D-8*PECM**2,4D0*PA2*STHE2*PA4*STHE4)
      IF(ABS(CPHI2).GE.1.05D0) MKERR=1
      CPHI2=MAX(-1D0,MIN(1D0,CPHI2))
      IF(MKERR.EQ.1) CALL PYERRM(13,
     &'(PY4ENT:) unphysical kinematical variable setup')
 
C...Store partons/particles in P vectors.
      P(IPA,3)=PA1
      P(IPA,4)=SQRT(PA1**2+PM1**2)
      P(IPA,5)=PM1
      P(IPA+3,1)=PA4*STHE4
      P(IPA+3,3)=PA4*CTHE4
      P(IPA+3,4)=SQRT(PA4**2+PM4**2)
      P(IPA+3,5)=PM4
      P(IPA+1,1)=PA2*STHE2*CPHI2
      P(IPA+1,2)=PA2*STHE2*SQRT(1D0-CPHI2**2)*(-1D0)**INT(PYR(0)+0.5D0)
      P(IPA+1,3)=PA2*CTHE2
      P(IPA+1,4)=SQRT(PA2**2+PM2**2)
      P(IPA+1,5)=PM2
      P(IPA+2,1)=-P(IPA+1,1)-P(IPA+3,1)
      P(IPA+2,2)=-P(IPA+1,2)
      P(IPA+2,3)=-P(IPA,3)-P(IPA+1,3)-P(IPA+3,3)
      P(IPA+2,4)=SQRT(P(IPA+2,1)**2+P(IPA+2,2)**2+P(IPA+2,3)**2+PM3**2)
      P(IPA+2,5)=PM3
 
C...Set N. Optionally fragment/decay.
      N=IPA+3
      IF(IP.EQ.0) CALL PYEXEC
 
      RETURN
      END
