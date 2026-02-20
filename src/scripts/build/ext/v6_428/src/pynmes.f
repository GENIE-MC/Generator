 
C*********************************************************************
 
C...PYNMES
C...Generates number of popcorn mesons and stores some relevant
C...parameters.
 
      SUBROUTINE PYNMES(KFDIQ)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 
      MSTU(121)=0
      IF(MSTJ(12).LT.2) RETURN
 
C..Old version: Get 1 or 0 popcorn mesons
      IF(MSTJ(12).LT.5)THEN
         POPWT=PARF(131)
         IF(KFDIQ.NE.0) THEN
            KFDIQA=IABS(KFDIQ)
            KFA=MOD(KFDIQA/1000,10)
            KFB=MOD(KFDIQA/100,10)
            KFS=MOD(KFDIQA,10)
            POPWT=PARF(132)
            IF(KFA.EQ.3) POPWT=PARF(133)
            IF(KFB.EQ.3) POPWT=PARF(134)
            IF(KFS.EQ.1) POPWT=POPWT*SQRT(PARJ(4))
         ENDIF
         MSTU(121)=INT(POPWT/(1D0+POPWT)+PYR(0))
         RETURN
      ENDIF
 
C..New version: Store popcorn- or rank 0 diquark parameters
      MSTU(122)=170
      PARF(193)=PARJ(8)
      PARF(194)=PARF(139)
      IF(KFDIQ.NE.0) THEN
         MSTU(122)=180
         PARF(193)=PARJ(10)
         PARF(194)=PARF(140)
      ENDIF
      IF(PARF(194).LT.1D-5.OR.PARF(194).GT.1D0-1D-5) THEN
         IF(PARF(194).GT.1D0-1D-5) CALL PYERRM(9,
     &        '(PYNMES:) Neglecting too large popcorn possibility')
         RETURN
      ENDIF
 
C..New version: Get number of popcorn mesons
  100 RTST=PYR(0)
      MSTU(121)=-1
  110 MSTU(121)=MSTU(121)+1
      RTST=RTST/PARF(194)
      IF(RTST.LT.1D0) GOTO 110
      IF(KFDIQ.EQ.0.AND.PYR(0)*(2D0+PARF(135)*PARF(161)).GT.
     &     (2D0+PARF(135)*PARF(161)*PARF(138)**MSTU(121))) GOTO 100
      RETURN
      END
