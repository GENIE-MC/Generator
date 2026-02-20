 
C*********************************************************************
 
C...PYOPER
C...Performs operations between histograms.
 
      SUBROUTINE PYOPER(ID1,OPER,ID2,ID3,F1,F2)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Commonblock.
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      SAVE /PYBINS/
C...Character variable.
      CHARACTER OPER*(*)
 
C...Find initial addresses in memory, and histogram size.
      IF(ID1.LE.0.OR.ID1.GT.IHIST(1)) CALL PYERRM(28,
     &'(PYFACT:) not allowed histogram number')
      IS1=INDX(ID1)
      IS2=INDX(MIN(IHIST(1),MAX(1,ID2)))
      IS3=INDX(MIN(IHIST(1),MAX(1,ID3)))
      NX=NINT(BIN(IS3+1))
      IF(OPER.EQ.'M'.AND.ID3.EQ.0) NX=NINT(BIN(IS2+1))
 
C...Update info on number of histogram entries.
      IF(OPER.EQ.'+'.OR.OPER.EQ.'-'.OR.OPER.EQ.'*'.OR.OPER.EQ.'/') THEN
        BIN(IS3+5)=BIN(IS1+5)+BIN(IS2+5)
      ELSEIF(OPER.EQ.'A'.OR.OPER.EQ.'S'.OR.OPER.EQ.'L') THEN
        BIN(IS3+5)=BIN(IS1+5)
      ENDIF
 
C...Operations on pair of histograms: addition, subtraction,
C...multiplication, division.
      IF(OPER.EQ.'+') THEN
        DO 100 IX=6,8+NX
          BIN(IS3+IX)=F1*BIN(IS1+IX)+F2*BIN(IS2+IX)
  100   CONTINUE
      ELSEIF(OPER.EQ.'-') THEN
        DO 110 IX=6,8+NX
          BIN(IS3+IX)=F1*BIN(IS1+IX)-F2*BIN(IS2+IX)
  110   CONTINUE
      ELSEIF(OPER.EQ.'*') THEN
        DO 120 IX=6,8+NX
          BIN(IS3+IX)=F1*BIN(IS1+IX)*F2*BIN(IS2+IX)
  120   CONTINUE
      ELSEIF(OPER.EQ.'/') THEN
        DO 130 IX=6,8+NX
          FA2=F2*BIN(IS2+IX)
          IF(ABS(FA2).LE.1D-20) THEN
            BIN(IS3+IX)=0D0
          ELSE
            BIN(IS3+IX)=F1*BIN(IS1+IX)/FA2
          ENDIF
  130   CONTINUE
 
C...Operations on single histogram: multiplication+addition,
C...square root+addition, logarithm+addition.
      ELSEIF(OPER.EQ.'A') THEN
        DO 140 IX=6,8+NX
          BIN(IS3+IX)=F1*BIN(IS1+IX)+F2
  140   CONTINUE
      ELSEIF(OPER.EQ.'S') THEN
        DO 150 IX=6,8+NX
          BIN(IS3+IX)=F1*SQRT(MAX(0D0,BIN(IS1+IX)))+F2
  150   CONTINUE
      ELSEIF(OPER.EQ.'L') THEN
        ZMIN=1D20
        DO 160 IX=9,8+NX
          IF(BIN(IS1+IX).LT.ZMIN.AND.BIN(IS1+IX).GT.1D-20)
     &    ZMIN=0.8D0*BIN(IS1+IX)
  160   CONTINUE
        DO 170 IX=6,8+NX
          BIN(IS3+IX)=F1*LOG10(MAX(ZMIN,BIN(IS1+IX)))+F2
  170   CONTINUE
 
C...Operation on two or three histograms: average and
C...standard deviation.
      ELSEIF(OPER.EQ.'M') THEN
        DO 180 IX=6,8+NX
          IF(ABS(BIN(IS1+IX)).LE.1D-20) THEN
            BIN(IS2+IX)=0D0
          ELSE
            BIN(IS2+IX)=BIN(IS2+IX)/BIN(IS1+IX)
          ENDIF
          IF(ID3.NE.0) THEN
            IF(ABS(BIN(IS1+IX)).LE.1D-20) THEN
              BIN(IS3+IX)=0D0
            ELSE
              BIN(IS3+IX)=SQRT(MAX(0D0,BIN(IS3+IX)/BIN(IS1+IX)-
     &        BIN(IS2+IX)**2))
            ENDIF
          ENDIF
          BIN(IS1+IX)=F1*BIN(IS1+IX)
  180   CONTINUE
      ENDIF
 
      RETURN
      END
