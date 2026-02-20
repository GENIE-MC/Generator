 
C*********************************************************************
 
C...PYPLOT
C...Prints a histogram (but does not reset it).
 
      SUBROUTINE PYPLOT(ID)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      SAVE /PYDAT1/,/PYBINS/
C...Local arrays and character variables.
      DIMENSION IDATI(6), IROW(100), IFRA(100), DYAC(10)
      CHARACTER TITLE*60, OUT*100, CHA(0:11)*1
 
C...Steps in histogram scale. Character sequence.
      DATA DYAC/.04,.05,.06,.08,.10,.12,.15,.20,.25,.30/
      DATA CHA/'0','1','2','3','4','5','6','7','8','9','X','-'/
 
C...Find initial address in memory; skip if empty histogram.
      IF(ID.LE.0.OR.ID.GT.IHIST(1)) RETURN
      IS=INDX(ID)
      IF(IS.EQ.0) RETURN
      IF(NINT(BIN(IS+5)).LE.0) THEN
        WRITE(MSTU(11),5000) ID
        RETURN
      ENDIF
 
C...Number of histogram lines and x bins.
      LIN=IHIST(3)-18
      NX=NINT(BIN(IS+1))
 
C...Extract title by conversion from double precision via integer.
      DO 100 IT=1,20
        IEQ=NINT(BIN(IS+8+NX+IT))
        TITLE(3*IT-2:3*IT)=CHAR(IEQ/256**2)//CHAR(MOD(IEQ,256**2)/256)
     &  //CHAR(MOD(IEQ,256))
  100 CONTINUE
 
C...Find time; print title.
      CALL PYTIME(IDATI)
      IF(IDATI(1).GT.0) THEN
        WRITE(MSTU(11),5100) ID, TITLE, (IDATI(J),J=1,5)
      ELSE
        WRITE(MSTU(11),5200) ID, TITLE
      ENDIF
 
C...Find minimum and maximum bin content.
      YMIN=BIN(IS+9)
      YMAX=BIN(IS+9)
      DO 110 IX=IS+10,IS+8+NX
        IF(BIN(IX).LT.YMIN) YMIN=BIN(IX)
        IF(BIN(IX).GT.YMAX) YMAX=BIN(IX)
  110 CONTINUE
 
C...Determine scale and step size for y axis.
      IF(YMAX-YMIN.GT.LIN*DYAC(1)*1D-9) THEN
        IF(YMIN.GT.0D0.AND.YMIN.LT.0.1D0*YMAX) YMIN=0D0
        IF(YMAX.LT.0D0.AND.YMAX.GT.0.1D0*YMIN) YMAX=0D0
        IPOT=INT(LOG10(YMAX-YMIN)+10D0)-10
        IF(YMAX-YMIN.LT.LIN*DYAC(1)*10D0**IPOT) IPOT=IPOT-1
        IF(YMAX-YMIN.GT.LIN*DYAC(10)*10D0**IPOT) IPOT=IPOT+1
        DELY=DYAC(1)
        DO 120 IDEL=1,9
          IF(YMAX-YMIN.GE.LIN*DYAC(IDEL)*10D0**IPOT) DELY=DYAC(IDEL+1)
  120   CONTINUE
        DY=DELY*10D0**IPOT
 
C...Convert bin contents to integer form; fractional fill in top row.
        DO 130 IX=1,NX
          CTA=ABS(BIN(IS+8+IX))/DY
          IROW(IX)=SIGN(CTA+0.95D0,BIN(IS+8+IX))
          IFRA(IX)=10D0*(CTA+1.05D0-DBLE(INT(CTA+0.95D0)))
  130   CONTINUE
        IRMI=SIGN(ABS(YMIN)/DY+0.95D0,YMIN)
        IRMA=SIGN(ABS(YMAX)/DY+0.95D0,YMAX)
 
C...Print histogram row by row.
        DO 150 IR=IRMA,IRMI,-1
          IF(IR.EQ.0) GOTO 150
          OUT=' '
          DO 140 IX=1,NX
            IF(IR.EQ.IROW(IX)) OUT(IX:IX)=CHA(IFRA(IX))
            IF(IR*(IROW(IX)-IR).GT.0) OUT(IX:IX)=CHA(10)
  140     CONTINUE
          WRITE(MSTU(11),5300) IR*DELY, IPOT, OUT
  150   CONTINUE
 
C...Print sign and value of bin contents.
        IPOT=INT(LOG10(MAX(YMAX,-YMIN))+10.0001D0)-10
        OUT=' '
        DO 160 IX=1,NX
          IF(BIN(IS+8+IX).LT.-10D0**(IPOT-4)) OUT(IX:IX)=CHA(11)
          IROW(IX)=NINT(10D0**(3-IPOT)*ABS(BIN(IS+8+IX)))
  160   CONTINUE
        WRITE(MSTU(11),5400) OUT
        DO 180 IR=4,1,-1
          DO 170 IX=1,NX
            OUT(IX:IX)=CHA(MOD(IROW(IX),10**IR)/10**(IR-1))
  170     CONTINUE
          WRITE(MSTU(11),5500) IPOT+IR-4, OUT
  180   CONTINUE
 
C...Print sign and value of lower bin edge.
        IPOT=INT(LOG10(MAX(-BIN(IS+2),BIN(IS+3)-BIN(IS+4)))+
     &  10.0001D0)-10
        OUT=' '
        DO 190 IX=1,NX
          IF(BIN(IS+2)+(IX-1)*BIN(IS+4).LT.-10D0**(IPOT-3))
     &    OUT(IX:IX)=CHA(11)
          IROW(IX)=NINT(10D0**(2-IPOT)*ABS(BIN(IS+2)+(IX-1)*BIN(IS+4)))
  190   CONTINUE
        WRITE(MSTU(11),5600) OUT
        DO 210 IR=3,1,-1
          DO 200 IX=1,NX
            OUT(IX:IX)=CHA(MOD(IROW(IX),10**IR)/10**(IR-1))
  200     CONTINUE
          WRITE(MSTU(11),5500) IPOT+IR-3, OUT
  210   CONTINUE
      ENDIF
 
C...Calculate and print statistics.
      CSUM=0D0
      CXSUM=0D0
      CXXSUM=0D0
      DO 220 IX=1,NX
        CTA=ABS(BIN(IS+8+IX))
        X=BIN(IS+2)+(IX-0.5D0)*BIN(IS+4)
        CSUM=CSUM+CTA
        CXSUM=CXSUM+CTA*X
        CXXSUM=CXXSUM+CTA*X**2
  220 CONTINUE
      XMEAN=CXSUM/MAX(CSUM,1D-20)
      XRMS=SQRT(MAX(0D0,CXXSUM/MAX(CSUM,1D-20)-XMEAN**2))
      WRITE(MSTU(11),5700) NINT(BIN(IS+5)),XMEAN,BIN(IS+6),
     &BIN(IS+2),BIN(IS+7),XRMS,BIN(IS+8),BIN(IS+3)
 
C...Formats for output.
 5000 FORMAT(/5X,'Histogram no',I5,' : no entries')
 5100 FORMAT('1'/5X,'Histogram no',I5,6X,A60,5X,I4,'-',I2,'-',I2,1X,
     &I2,':',I2/)
 5200 FORMAT('1'/5X,'Histogram no',I5,6X,A60/)
 5300 FORMAT(2X,F7.2,'*10**',I2,3X,A100)
 5400 FORMAT(/8X,'Contents',3X,A100)
 5500 FORMAT(9X,'*10**',I2,3X,A100)
 5600 FORMAT(/8X,'Low edge',3X,A100)
 5700 FORMAT(/5X,'Entries  =',I12,1P,6X,'Mean =',D12.4,6X,'Underflow ='
     &,D12.4,6X,'Low edge  =',D12.4/5X,'All chan =',D12.4,6X,
     &'Rms  =',D12.4,6X,'Overflow  =',D12.4,6X,'High edge =',D12.4)
 
      RETURN
      END
