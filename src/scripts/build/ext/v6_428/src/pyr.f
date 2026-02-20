 
C*********************************************************************
 
C...PYR
C...Generates random numbers uniformly distributed between
C...0 and 1, excluding the endpoints.
 
      FUNCTION PYR(IDUMMY)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDATR/MRPY(6),RRPY(100)
      SAVE /PYDATR/
C...Equivalence between commonblock and local variables.
      EQUIVALENCE (MRPY1,MRPY(1)),(MRPY2,MRPY(2)),(MRPY3,MRPY(3)),
     &(MRPY4,MRPY(4)),(MRPY5,MRPY(5)),(MRPY6,MRPY(6)),
     &(RRPY98,RRPY(98)),(RRPY99,RRPY(99)),(RRPY00,RRPY(100))
 
C...Initialize generation from given seed.
      IF(MRPY2.EQ.0) THEN
        IJ=MOD(MRPY1/30082,31329)
        KL=MOD(MRPY1,30082)
        I=MOD(IJ/177,177)+2
        J=MOD(IJ,177)+2
        K=MOD(KL/169,178)+1
        L=MOD(KL,169)
        DO 110 II=1,97
          S=0D0
          T=0.5D0
          DO 100 JJ=1,48
            M=MOD(MOD(I*J,179)*K,179)
            I=J
            J=K
            K=M
            L=MOD(53*L+1,169)
            IF(MOD(L*M,64).GE.32) S=S+T
            T=0.5D0*T
  100     CONTINUE
          RRPY(II)=S
  110   CONTINUE
        TWOM24=1D0
        DO 120 I24=1,24
          TWOM24=0.5D0*TWOM24
  120   CONTINUE
        RRPY98=362436D0*TWOM24
        RRPY99=7654321D0*TWOM24
        RRPY00=16777213D0*TWOM24
        MRPY2=1
        MRPY3=0
        MRPY4=97
        MRPY5=33
      ENDIF
 
C...Generate next random number.
  130 RUNI=RRPY(MRPY4)-RRPY(MRPY5)
      IF(RUNI.LT.0D0) RUNI=RUNI+1D0
      RRPY(MRPY4)=RUNI
      MRPY4=MRPY4-1
      IF(MRPY4.EQ.0) MRPY4=97
      MRPY5=MRPY5-1
      IF(MRPY5.EQ.0) MRPY5=97
      RRPY98=RRPY98-RRPY99
      IF(RRPY98.LT.0D0) RRPY98=RRPY98+RRPY00
      RUNI=RUNI-RRPY98
      IF(RUNI.LT.0D0) RUNI=RUNI+1D0
      IF(RUNI.LE.0D0.OR.RUNI.GE.1D0) GOTO 130
 
C...Update counters. Random number to output.
      MRPY3=MRPY3+1
      IF(MRPY3.EQ.1000000000) THEN
        MRPY2=MRPY2+1
        MRPY3=0
      ENDIF
      PYR=RUNI
 
      RETURN
      END
