 
C*********************************************************************
 
C...VISAJE
C...Dummy function, to be removed when ISAJET (ISASUSY) is to be linked.
 
      FUNCTION VISAJE()
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      CHARACTER*40 VISAJE
 
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
C...Assign default value.
      VISAJE='Undefined'
 
C...Stop program if this routine is ever called.
      WRITE(MSTU(11),5000)
      CALL PYSTOP(110)
 
C...Format for error printout.
 5000 FORMAT(1X,'Error: you did not link ISAJET correctly.'/
     &1X,'Dummy function VISAJE in PYTHIA file called instead.'/
     &1X,'Execution stopped!')
 
      RETURN
      END
