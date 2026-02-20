 
C*********************************************************************
 
C...PYSTOP
C...Allows users to handle STOP statemens
 
      SUBROUTINE PYSTOP(MCOD)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/

 
C...Write message, then stop
      WRITE(MSTU(11),5000) MCOD
      STOP

 
C...Formats for output.
 5000 FORMAT(/5X,'PYSTOP called with code: ',I4)
      END
