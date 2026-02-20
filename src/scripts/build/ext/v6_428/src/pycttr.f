C*********************************************************************
 
C...PYCTTR
C...Adapted from PYPREP.
C...Assigns LHA1 colour tags to coloured partons based on
C...K(I,4) and K(I,5) colour connection record.
C...KCS negative signifies that a previous tracing should be continued.
C...(in case the tag to be continued is empty, the routine exits)
C...Starts at I and ends at I or IEND.
C...Special considerations for systems with junctions.
C...Special: if IEND=-1, means trace this parton to its color partner,
C...         then exit. If no partner found, exit with 0. 

      SUBROUTINE PYCTTR(I,KCS,IEND)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYINT1/MINT(400),VINT(400)
C...The common block of colour tags.
      COMMON/PYCTAG/NCT,MCT(4000,2)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYINT1/,/PYCTAG/
      DATA NERRPR/0/
      SAVE NERRPR
 
C...Skip if parton not existing or does not have KCS
      IF (K(I,1).LE.0) GOTO 120
      KC=PYCOMP(K(I,2))
      IF (KC.EQ.0) GOTO 120
      KQ=KCHG(KC,2)
      IF (KQ.EQ.0) GOTO 120
      IF (IABS(KQ).EQ.1.AND.KQ*(9-2*ABS(KCS)).NE.ISIGN(1,K(I,2))) 
     &    GOTO 120
 
      IF (KCS.GT.0) THEN
        NCT=NCT+1
C...Set colour tag of first parton.
        MCT(I,KCS-3)=NCT
        NCS=NCT
      ELSE
        KCS=-KCS
        NCS=MCT(I,KCS-3)
        IF (NCS.EQ.0) GOTO 120
      ENDIF
 
      IA=I
      NSTP=0
  100 NSTP=NSTP+1
      IF(NSTP.GT.4*N) THEN
        CALL PYERRM(14,'(PYCTTR:) caught in infinite loop')
        GOTO 120
      ENDIF
 
C...Finished if reached final-state triplet.
      IF(K(IA,1).EQ.3) THEN
        IF(NSTP.GE.2.AND.KCHG(PYCOMP(K(IA,2)),2).NE.2) GOTO 120
      ENDIF
 
C...Also finished if reached junction.
      IF(K(IA,1).EQ.42) THEN
        GOTO 120
      ENDIF
 
C...GOTO next parton in colour space.
  110 IB=IA
C...If IB's KCS daughter not traced and exists, goto KCS daughter.
      IF(MOD(K(IB,KCS)/MSTU(5)**2,2).EQ.0.AND.MOD(K(IB,KCS),MSTU(5))
     &     .NE.0) THEN
        IA=MOD(K(IB,KCS),MSTU(5))
        K(IB,KCS)=K(IB,KCS)+MSTU(5)**2
        MREV=0
      ELSE
C...If KCS mother traced or KCS mother nonexistent, switch colour.
        IF(K(IB,KCS).GE.2*MSTU(5)**2.OR.MOD(K(IB,KCS)/MSTU(5),
     &       MSTU(5)).EQ.0) THEN
          KCS=9-KCS
          NCT=NCT+1
          NCS=NCT
C...Assign new colour tag on other side of old parton.
          MCT(IB,KCS-3)=NCT
        ENDIF
C...Goto (new) KCS mother, set mother traced tag
        IA=MOD(K(IB,KCS)/MSTU(5),MSTU(5))
        K(IB,KCS)=K(IB,KCS)+2*MSTU(5)**2
        MREV=1
      ENDIF
      IF(IA.LE.0.OR.IA.GT.N) THEN
        IF (IEND.EQ.-1) THEN
          IEND=0
          GOTO 120
        ENDIF
        CALL PYERRM(12,'(PYCTTR:) colour tag tracing failed')
        IF(NERRPR.LT.5) THEN
          write(*,*) 'began at ',I
          write(*,*) 'ended going from', IB, ' to', IA, '  KCS=',KCS,
     &        '  NCS=',NCS,'  MREV=',MREV
          CALL PYLIST(4)
          NERRPR=NERRPR+1
        ENDIF
        MINT(51)=1
        RETURN
      ENDIF
      IF(MOD(K(IA,4)/MSTU(5),MSTU(5)).EQ.IB.OR.MOD(K(IA,5)/MSTU(5),
     &     MSTU(5)).EQ.IB) THEN
        IF(MREV.EQ.1) KCS=9-KCS
        IF(MOD(K(IA,KCS)/MSTU(5),MSTU(5)).NE.IB) KCS=9-KCS
C...Set KSC mother traced tag for IA
        K(IA,KCS)=K(IA,KCS)+2*MSTU(5)**2
      ELSE
        IF(MREV.EQ.0) KCS=9-KCS
        IF(MOD(K(IA,KCS),MSTU(5)).NE.IB) KCS=9-KCS
C...Set KCS daughter traced tag for IA
        K(IA,KCS)=K(IA,KCS)+MSTU(5)**2
      ENDIF
C...Assign new colour tag
      MCT(IA,KCS-3)=NCS
C...Finish if IEND=-1 and found final-state color partner 
      IF (IEND.EQ.-1.AND.K(IA,1).LT.10) THEN
        IEND=IA
        GOTO 120        
      ENDIF
      IF (IA.NE.I.AND.IA.NE.IEND) GOTO 100
 
  120 RETURN
      END
