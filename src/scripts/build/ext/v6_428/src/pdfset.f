 
C*********************************************************************
 
C...PDFSET
C...Dummy routine, to be removed when PDFLIB is to be linked.
 
      SUBROUTINE PDFSET(PARM,VALUE)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
C...Local arrays and character variables.
      CHARACTER*20 PARM(20)
      DOUBLE PRECISION VALUE(20)
 
C...Stop program if this routine is ever called.
      WRITE(MSTU(11),5000)
      CALL PYSTOP(5)
      PARM(20)=PARM(1)
      VALUE(20)=VALUE(1)
 
C...Format for error printout.
 5000 FORMAT(1X,'Error: you did not link PDFLIB correctly.'/
     &1X,'Dummy routine PDFSET in PYTHIA file called instead.'/
     &1X,'Execution stopped!')
 
      RETURN
      END
