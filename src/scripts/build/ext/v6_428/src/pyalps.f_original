 
C*********************************************************************
 
C...PYALPS
C...Gives the value of alpha_strong.
 
      FUNCTION PYALPS(Q2)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
C...Coefficients for second-order threshold matching.
C...From W.J. Marciano, Phys. Rev. D29 (1984) 580.
      DIMENSION STEPDN(6),STEPUP(6)
c      DATA STEPDN/0D0,0D0,(2D0*107D0/2025D0),(2D0*963D0/14375D0),
c     &(2D0*321D0/3703D0),0D0/
c      DATA STEPUP/0D0,0D0,0D0,(-2D0*107D0/1875D0),
c     &(-2D0*963D0/13225D0),(-2D0*321D0/3381D0)/
      DATA STEPDN/0D0,0D0,0.10568D0,0.13398D0,0.17337D0,0D0/
      DATA STEPUP/0D0,0D0,0D0,-0.11413D0,-0.14563D0,-0.18988D0/
 
C...Constant alpha_strong trivial. Pick artificial Lambda.
      IF(MSTU(111).LE.0) THEN
        PYALPS=PARU(111)
        MSTU(118)=MSTU(112)
        PARU(117)=0.2D0
        IF(Q2.GT.0.04D0) PARU(117)=SQRT(Q2)*EXP(-6D0*PARU(1)/
     &  ((33D0-2D0*MSTU(112))*PARU(111)))
        PARU(118)=PARU(111)
        RETURN
      ENDIF
 
C...Find effective Q2, number of flavours and Lambda.
      Q2EFF=Q2
      IF(MSTU(115).GE.2) Q2EFF=MAX(Q2,PARU(114))
      NF=MSTU(112)
      ALAM2=PARU(112)**2
  100 IF(NF.GT.MAX(3,MSTU(113))) THEN
        Q2THR=PARU(113)*PMAS(NF,1)**2
        IF(Q2EFF.LT.Q2THR) THEN
          NF=NF-1
          Q2RAT=Q2THR/ALAM2
          ALAM2=ALAM2*Q2RAT**(2D0/(33D0-2D0*NF))
          IF(MSTU(111).EQ.2) ALAM2=ALAM2*LOG(Q2RAT)**STEPDN(NF)
          GOTO 100
        ENDIF
      ENDIF
  110 IF(NF.LT.MIN(6,MSTU(114))) THEN
        Q2THR=PARU(113)*PMAS(NF+1,1)**2
        IF(Q2EFF.GT.Q2THR) THEN
          NF=NF+1
          Q2RAT=Q2THR/ALAM2
          ALAM2=ALAM2*Q2RAT**(-2D0/(33D0-2D0*NF))
          IF(MSTU(111).EQ.2) ALAM2=ALAM2*LOG(Q2RAT)**STEPUP(NF)
          GOTO 110
        ENDIF
      ENDIF
      IF(MSTU(115).EQ.1) Q2EFF=Q2EFF+ALAM2
      PARU(117)=SQRT(ALAM2)
 
C...Evaluate first or second order alpha_strong.
      B0=(33D0-2D0*NF)/6D0
      ALGQ=LOG(MAX(1.0001D0,Q2EFF/ALAM2))
      IF(MSTU(111).EQ.1) THEN
        PYALPS=MIN(PARU(115),PARU(2)/(B0*ALGQ))
      ELSE
        B1=(153D0-19D0*NF)/6D0
        PYALPS=MIN(PARU(115),PARU(2)/(B0*ALGQ)*(1D0-B1*LOG(ALGQ)/
     &  (B0**2*ALGQ)))
      ENDIF
      MSTU(118)=NF
      PARU(118)=PYALPS
 
      RETURN
      END
