 
C*********************************************************************
 
C...PYRSET
C...Reads a state of the random number generator from a file
C...for subsequent generation from this state onwards.
 
      SUBROUTINE PYRSET(LFN,MOVE)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDATR/MRPY(6),RRPY(100)
      SAVE /PYDATR/
C...Local character variable.
      CHARACTER CHERR*8
 
C...Backspace required number of records (or as many as there are).
      IF(MOVE.LT.0) THEN
        NBCK=MIN(MRPY(6),-MOVE)
        DO 100 IBCK=1,NBCK
          BACKSPACE(LFN,ERR=120,IOSTAT=IERR)
  100   CONTINUE
        MRPY(6)=MRPY(6)-NBCK
      ENDIF
 
C...Unformatted read from unit LFN.
      NFOR=1+MAX(0,MOVE)
      DO 110 IFOR=1,NFOR
        READ(LFN,ERR=120,IOSTAT=IERR) (MRPY(I1),I1=1,5),
     &  (RRPY(I2),I2=1,100)
  110 CONTINUE
      MRPY(6)=MRPY(6)+NFOR
      RETURN
 
C...Write error.
  120 WRITE(CHERR,'(I8)') IERR
      CALL PYERRM(18,'(PYRSET:) error when accessing file, IOSTAT ='//
     &CHERR)
 
      RETURN
      END
