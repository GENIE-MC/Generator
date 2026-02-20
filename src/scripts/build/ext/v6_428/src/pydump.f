 
C*********************************************************************
 
C...PYDUMP
C...Dumps histogram contents on file for reading by other program.
C...Can also read back own dump.
 
      SUBROUTINE PYDUMP(MDUMP,LFN,NHI,IHI)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Commonblock.
      COMMON/PYBINS/IHIST(4),INDX(1000),BIN(20000)
      SAVE /PYBINS/
C...Local arrays and character variables.
      DIMENSION IHI(*),ISS(100),VAL(5)
      CHARACTER TITLE*60,FORMAT*13
 
C...Dump all histograms that have been booked,
C...including titles and ranges, one after the other.
      IF(MDUMP.EQ.1) THEN
 
C...Loop over histograms and find which are wanted and booked.
        IF(NHI.LE.0) THEN
          NW=IHIST(1)
        ELSE
          NW=NHI
        ENDIF
        DO 130 IW=1,NW
          IF(NHI.EQ.0) THEN
            ID=IW
          ELSE
            ID=IHI(IW)
          ENDIF
          IS=INDX(ID)
          IF(IS.NE.0) THEN
 
C...Write title, histogram size, filling statistics.
            NX=NINT(BIN(IS+1))
            DO 100 IT=1,20
              IEQ=NINT(BIN(IS+8+NX+IT))
              TITLE(3*IT-2:3*IT)=CHAR(IEQ/256**2)//
     &        CHAR(MOD(IEQ,256**2)/256)//CHAR(MOD(IEQ,256))
  100       CONTINUE
            WRITE(LFN,5100) ID,TITLE
            WRITE(LFN,5200) NX,BIN(IS+2),BIN(IS+3)
            WRITE(LFN,5300) NINT(BIN(IS+5)),BIN(IS+6),BIN(IS+7),
     &      BIN(IS+8)
 
 
C...Write histogram contents, in groups of five.
            DO 120 IXG=1,(NX+4)/5
              DO 110 IXV=1,5
                IX=5*IXG+IXV-5
                IF(IX.LE.NX) THEN
                  VAL(IXV)=BIN(IS+8+IX)
                ELSE
                  VAL(IXV)=0D0
                ENDIF
  110         CONTINUE
              WRITE(LFN,5400) (VAL(IXV),IXV=1,5)
  120       CONTINUE
 
C...Go to next histogram; finish.
          ELSEIF(NHI.GT.0) THEN
            CALL PYERRM(8,'(PYDUMP:) unknown histogram number')
          ENDIF
  130   CONTINUE
 
C...Read back in histograms dumped MDUMP=1.
      ELSEIF(MDUMP.EQ.2) THEN
 
C...Read histogram number, title and range, and book.
  140   READ(LFN,5100,END=170) ID,TITLE
        READ(LFN,5200) NX,XL,XU
        CALL PYBOOK(ID,TITLE,NX,XL,XU)
        IS=INDX(ID)
 
C...Read filling statistics.
        READ(LFN,5300) NENTRY,BIN(IS+6),BIN(IS+7),BIN(IS+8)
        BIN(IS+5)=DBLE(NENTRY)
 
C...Read histogram contents, in groups of five.
        DO 160 IXG=1,(NX+4)/5
          READ(LFN,5400) (VAL(IXV),IXV=1,5)
          DO 150 IXV=1,5
            IX=5*IXG+IXV-5
            IF(IX.LE.NX) BIN(IS+8+IX)=VAL(IXV)
  150     CONTINUE
  160   CONTINUE
 
C...Go to next histogram; finish.
        GOTO 140
  170   CONTINUE
 
C...Write histogram contents in column format,
C...convenient e.g. for GNUPLOT input.
      ELSEIF(MDUMP.EQ.3) THEN
 
C...Find addresses to wanted histograms.
        NSS=0
        IF(NHI.LE.0) THEN
          NW=IHIST(1)
        ELSE
          NW=NHI
        ENDIF
        DO 180 IW=1,NW
          IF(NHI.EQ.0) THEN
            ID=IW
          ELSE
            ID=IHI(IW)
          ENDIF
          IS=INDX(ID)
          IF(IS.NE.0.AND.NSS.LT.100) THEN
            NSS=NSS+1
            ISS(NSS)=IS
          ELSEIF(NSS.GE.100) THEN
            CALL PYERRM(8,'(PYDUMP:) too many histograms requested')
          ELSEIF(NHI.GT.0) THEN
            CALL PYERRM(8,'(PYDUMP:) unknown histogram number')
          ENDIF
  180   CONTINUE
 
C...Check that they have common number of x bins. Fix format.
        NX=NINT(BIN(ISS(1)+1))
        DO 190 IW=2,NSS
          IF(NINT(BIN(ISS(IW)+1)).NE.NX) THEN
            CALL PYERRM(8,'(PYDUMP:) different number of bins')
            RETURN
          ENDIF
  190   CONTINUE
        FORMAT='(1P,000E12.4)'
        WRITE(FORMAT(5:7),'(I3)') NSS+1
 
C...Write histogram contents; first column x values.
        DO 200 IX=1,NX
          X=BIN(ISS(1)+2)+(IX-0.5D0)*BIN(ISS(1)+4)
          WRITE(LFN,FORMAT) X, (BIN(ISS(IW)+8+IX),IW=1,NSS)
  200   CONTINUE
 
      ENDIF
 
C...Formats for output.
 5100 FORMAT(I5,5X,A60)
 5200 FORMAT(I5,1P,2D12.4)
 5300 FORMAT(I12,1P,3D12.4)
 5400 FORMAT(1P,5D12.4)
 
      RETURN
      END
