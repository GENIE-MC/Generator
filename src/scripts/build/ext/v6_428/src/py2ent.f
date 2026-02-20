 
C*********************************************************************
 
C...PY2ENT
C...Stores two partons/particles in their CM frame,
C...with the first along the +z axis.
 
      SUBROUTINE PY2ENT(IP,KF1,KF2,PECM)
 
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
      IF(IPA.GT.MSTU(4)-1) CALL PYERRM(21,
     &'(PY2ENT:) writing outside PYJETS memory')
      KC1=PYCOMP(KF1)
      KC2=PYCOMP(KF2)
      IF(KC1.EQ.0.OR.KC2.EQ.0) CALL PYERRM(12,
     &'(PY2ENT:) unknown flavour code')
 
C...Find masses. Reset K, P and V vectors.
      PM1=0D0
      IF(MSTU(10).EQ.1) PM1=P(IPA,5)
      IF(MSTU(10).GE.2) PM1=PYMASS(KF1)
      PM2=0D0
      IF(MSTU(10).EQ.1) PM2=P(IPA+1,5)
      IF(MSTU(10).GE.2) PM2=PYMASS(KF2)
      DO 110 I=IPA,IPA+1
        DO 100 J=1,5
          K(I,J)=0
          P(I,J)=0D0
          V(I,J)=0D0
  100   CONTINUE
  110 CONTINUE
 
C...Check flavours.
      KQ1=KCHG(KC1,2)*ISIGN(1,KF1)
      KQ2=KCHG(KC2,2)*ISIGN(1,KF2)
      IF(MSTU(19).EQ.1) THEN
        MSTU(19)=0
      ELSE
        IF(KQ1+KQ2.NE.0.AND.KQ1+KQ2.NE.4) CALL PYERRM(2,
     &  '(PY2ENT:) unphysical flavour combination')
      ENDIF
      K(IPA,2)=KF1
      K(IPA+1,2)=KF2
 
C...Store partons/particles in K vectors for normal case.
      IF(IP.GE.0) THEN
        K(IPA,1)=1
        IF(KQ1.NE.0.AND.KQ2.NE.0) K(IPA,1)=2
        K(IPA+1,1)=1
 
C...Store partons in K vectors for parton shower evolution.
      ELSE
        K(IPA,1)=3
        K(IPA+1,1)=3
        K(IPA,4)=MSTU(5)*(IPA+1)
        K(IPA,5)=K(IPA,4)
        K(IPA+1,4)=MSTU(5)*IPA
        K(IPA+1,5)=K(IPA+1,4)
      ENDIF
 
C...Check kinematics and store partons/particles in P vectors.
      IF(PECM.LE.PM1+PM2) CALL PYERRM(13,
     &'(PY2ENT:) energy smaller than sum of masses')
      PA=SQRT(MAX(0D0,(PECM**2-PM1**2-PM2**2)**2-(2D0*PM1*PM2)**2))/
     &(2D0*PECM)
      P(IPA,3)=PA
      P(IPA,4)=SQRT(PM1**2+PA**2)
      P(IPA,5)=PM1
      P(IPA+1,3)=-PA
      P(IPA+1,4)=SQRT(PM2**2+PA**2)
      P(IPA+1,5)=PM2
 
C...Set N. Optionally fragment/decay.
      N=IPA+1
      IF(IP.EQ.0) CALL PYEXEC
 
      RETURN
      END
