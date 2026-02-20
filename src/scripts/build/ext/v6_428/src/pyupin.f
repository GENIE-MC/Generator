 
C*********************************************************************
 
C...PYUPIN
C...Fills the HEPRUP commonblock with info on incoming beams and allowed
C...processes, and optionally stores that information on file.
 
      SUBROUTINE PYUPIN
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      SAVE /PYJETS/,/PYSUBS/,/PYPARS/,/PYINT5/
 
C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/
 
C...Store info on incoming beams.
      IDBMUP(1)=K(1,2)
      IDBMUP(2)=K(2,2)
      EBMUP(1)=P(1,4)
      EBMUP(2)=P(2,4)
      PDFGUP(1)=0
      PDFGUP(2)=0
      PDFSUP(1)=MSTP(51)
      PDFSUP(2)=MSTP(51)
 
C...Event weighting strategy.
      IDWTUP=3
 
C...Info on individual processes.
      NPRUP=0
      DO 100 ISUB=1,500
        IF(MSUB(ISUB).EQ.1) THEN
          NPRUP=NPRUP+1
          XSECUP(NPRUP)=1D9*XSEC(ISUB,3)
          XERRUP(NPRUP)=XSECUP(NPRUP)/SQRT(MAX(1D0,DBLE(NGEN(ISUB,3))))
          XMAXUP(NPRUP)=1D0
          LPRUP(NPRUP)=ISUB
        ENDIF
  100 CONTINUE
 
C...Write info to file.
      IF(MSTP(161).GT.0) THEN
        WRITE(MSTP(161),5100) IDBMUP(1),IDBMUP(2),EBMUP(1),EBMUP(2),
     &  PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
        DO 110 IPR=1,NPRUP
          WRITE(MSTP(161),5200) XSECUP(IPR),XERRUP(IPR),XMAXUP(IPR),
     &    LPRUP(IPR)
  110   CONTINUE
      ENDIF
 
C...Formats for printout.
 5100 FORMAT(1P,2I8,2E14.6,6I6)
 5200 FORMAT(1P,3E14.6,I6)
 
      RETURN
      END
