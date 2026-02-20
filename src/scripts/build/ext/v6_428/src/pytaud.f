  
C*********************************************************************
 
C...PYTAUD
C...Dummy routine, to be replaced by user, to handle the decay of a
C...polarized tau lepton.
C...Input:
C...ITAU is the position where the decaying tau is stored in /PYJETS/.
C...IORIG is the position where the mother of the tau is stored;
C...     is 0 when the mother is not stored.
C...KFORIG is the flavour of the mother of the tau;
C...     is 0 when the mother is not known.
C...Note that IORIG=0 does not necessarily imply KFORIG=0;
C...     e.g. in B hadron semileptonic decays the W  propagator
C...     is not explicitly stored but the W code is still unambiguous.
C...Output:
C...NDECAY is the number of decay products in the current tau decay.
C...These decay products should be added to the /PYJETS/ common block,
C...in positions N+1 through N+NDECAY. For each product I you must
C...give the flavour codes K(I,2) and the five-momenta P(I,1), P(I,2),
C...P(I,3), P(I,4) and P(I,5). The rest will be stored automatically.
 
      SUBROUTINE PYTAUD(ITAU,IORIG,KFORIG,NDECAY)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYJETS/,/PYDAT1/
 
C...Stop program if this routine is ever called.
C...You should not copy these lines to your own routine.
      NDECAY=ITAU+IORIG+KFORIG
      WRITE(MSTU(11),5000)
      CALL PYSTOP(10)
 
C...Format for error printout.
 5000 FORMAT(1X,'Error: you did not link your PYTAUD routine ',
     &'correctly.'/1X,'Dummy routine in PYTHIA file called instead.'/
     &1X,'Execution stopped!')
 
      RETURN
      END
