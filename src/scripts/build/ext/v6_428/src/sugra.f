 
C*********************************************************************
 
C...SUGRA
C...Dummy routine, to be removed when ISAJET (ISASUSY) is to be linked.
 
      SUBROUTINE SUGRA(MZERO,MHLF,AZERO,TANB,SGNMU,MTOP,IMODL)
       IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      REAL MZERO,MHLF,AZERO,TANB,SGNMU,MTOP
      INTEGER IMODL
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
C...Stop program if this routine is ever called.
      WRITE(MSTU(11),5000)
      CALL PYSTOP(110)
 
C...Format for error printout.
 5000 FORMAT(1X,'Error: you did not link ISAJET correctly.'/
     &1X,'Dummy routine SUGRA in PYTHIA file called instead.'/
     &1X,'Execution stopped!')
 
      RETURN
      END
