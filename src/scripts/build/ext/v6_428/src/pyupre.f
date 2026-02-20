 
C*********************************************************************
 
C...PYUPRE
C...Rearranges contents of the HEPEUP commonblock so that
C...mothers precede daughters and daughters of a decay are
C...listed consecutively.
 
      SUBROUTINE PYUPRE
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/
 
C...Local arrays.
      DIMENSION NEWPOS(0:MAXNUP),IDUPT(MAXNUP),ISTUPT(MAXNUP),
     &MOTUPT(2,MAXNUP),ICOUPT(2,MAXNUP),PUPT(5,MAXNUP),
     &VTIUPT(MAXNUP),SPIUPT(MAXNUP)
 
C...Check whether a rearrangement is required.
      NEED=0
      DO 100 IUP=1,NUP
        IF(MOTHUP(1,IUP).GT.IUP) NEED=NEED+1
  100 CONTINUE
      DO 110 IUP=2,NUP
        IF(MOTHUP(1,IUP).LT.MOTHUP(1,IUP-1)) NEED=NEED+1
  110 CONTINUE
 
      IF(NEED.NE.0) THEN
C...Find the new order that particles should have.
        NEWPOS(0)=0
        NNEW=0
        INEW=-1
  120   INEW=INEW+1
        DO 130 IUP=1,NUP
          IF(MOTHUP(1,IUP).EQ.NEWPOS(INEW)) THEN
            NNEW=NNEW+1
            NEWPOS(NNEW)=IUP
          ENDIF
  130   CONTINUE
        IF(INEW.LT.NNEW.AND.INEW.LT.NUP) GOTO 120
        IF(NNEW.NE.NUP) THEN
          CALL PYERRM(2,
     &    '(PYUPRE:) failed to make sense of mother pointers in HEPEUP')
          RETURN
        ENDIF
 
C...Copy old info into temporary storage.
        DO 150 I=1,NUP
          IDUPT(I)=IDUP(I)
          ISTUPT(I)=ISTUP(I)
          MOTUPT(1,I)=MOTHUP(1,I)
          MOTUPT(2,I)=MOTHUP(2,I)
          ICOUPT(1,I)=ICOLUP(1,I)
          ICOUPT(2,I)=ICOLUP(2,I)
          DO 140 J=1,5
            PUPT(J,I)=PUP(J,I)
  140     CONTINUE
          VTIUPT(I)=VTIMUP(I)
          SPIUPT(I)=SPINUP(I)
  150   CONTINUE
 
C...Copy info back into HEPEUP in right order.
        DO 180 I=1,NUP
          IOLD=NEWPOS(I)
          IDUP(I)=IDUPT(IOLD)
          ISTUP(I)=ISTUPT(IOLD)
          MOTHUP(1,I)=0
          MOTHUP(2,I)=0
          DO 160 IMOT=1,I-1
            IF(MOTUPT(1,IOLD).EQ.NEWPOS(IMOT)) MOTHUP(1,I)=IMOT
            IF(MOTUPT(2,IOLD).EQ.NEWPOS(IMOT)) MOTHUP(2,I)=IMOT
  160     CONTINUE
          IF(MOTHUP(2,I).GT.0.AND.MOTHUP(2,I).LT.MOTHUP(1,I)) THEN
            MOTHSW=MOTHUP(1,I)
            MOTHUP(1,I)=MOTHUP(2,I)
            MOTHUP(2,I)=MOTHSW
          ENDIF
          ICOLUP(1,I)=ICOUPT(1,IOLD)
          ICOLUP(2,I)=ICOUPT(2,IOLD)
          DO 170 J=1,5
            PUP(J,I)=PUPT(J,IOLD)
  170     CONTINUE
          VTIMUP(I)=VTIUPT(IOLD)
          SPINUP(I)=SPIUPT(IOLD)
  180   CONTINUE
      ENDIF
 
c...If incoming particles are massive recalculate to put them massless.
      IF(PUP(5,1).NE.0D0.OR.PUP(5,2).NE.0D0) THEN
        PPLUS=(PUP(4,1)+PUP(3,1))+(PUP(4,2)+PUP(3,2))
        PMINUS=(PUP(4,1)-PUP(3,1))+(PUP(4,2)-PUP(3,2))
        PUP(4,1)=0.5D0*PPLUS
        PUP(3,1)=PUP(4,1)
        PUP(5,1)=0D0
        PUP(4,2)=0.5D0*PMINUS
        PUP(3,2)=-PUP(4,2)
        PUP(5,2)=0D0
      ENDIF
 
      RETURN
      END
