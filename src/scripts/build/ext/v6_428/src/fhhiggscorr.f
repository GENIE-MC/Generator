 
C*********************************************************************
 
C...FHHIGGSCORR
C...Dummy function, to be removed when FEYNHIGGS is to be linked.
 
      SUBROUTINE FHHIGGSCORR(IERR, RMHIGG, SAEFF, UHIGGS)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
C...FeynHiggs variables
      DOUBLE PRECISION RMHIGG(4)
      DOUBLE COMPLEX SAEFF, UHIGGS(3,3)
      DOUBLE COMPLEX DMU,
     &     AE33, AU33, AD33, AE22, AU22, AD22, AE11, AU11, AD11,
     &     DM1, DM2, DM3

C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
C...Stop program if this routine is ever called.
      WRITE(MSTU(11),5000)
      CALL PYSTOP(103)
 
C...Format for error printout.
 5000 FORMAT(1X,'Error: you did not link FEYNHIGGS correctly.'/
     &1X,'Dummy routine FHSETPARA in PYTHIA file called instead.'/
     &1X,'Execution stopped!')
      RETURN
      END
