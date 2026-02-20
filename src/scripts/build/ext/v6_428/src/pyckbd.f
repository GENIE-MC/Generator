 
C*********************************************************************
 
C...PYCKBD
C...Check that BLOCK DATA PYDATA has been loaded.
C...Should not be required, except that some compilers/linkers
C...are pretty buggy in this respect.
 
      SUBROUTINE PYCKBD
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/
 
C...Check a few variables to see they have been sensibly initialized.
      IF(MSTU(4).LT.10.OR.MSTU(4).GT.900000.OR.PMAS(2,1).LT.0.001D0
     &.OR.PMAS(2,1).GT.1D0.OR.CKIN(5).LT.0.01D0.OR.MSTP(1).LT.1.OR.
     &MSTP(1).GT.5) THEN
C...If not, abort the run right away.
        WRITE(*,*) 'Fatal error: BLOCK DATA PYDATA has not been loaded!'
        WRITE(*,*) 'The program execution is stopped now!'
        CALL PYSTOP(8)
      ENDIF
 
      RETURN
      END
