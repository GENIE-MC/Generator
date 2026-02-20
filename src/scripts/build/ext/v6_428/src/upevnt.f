
C...Old example: handles a simple Pythia 6.4 initialization file.
 
c      SUBROUTINE UPINIT
 
C...Double precision and integer declarations.
c      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
c      IMPLICIT INTEGER(I-N)
 
C...Commonblocks.
c      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
c      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
c      SAVE /PYDAT1/,/PYPARS/
 
C...User process initialization commonblock.
c      INTEGER MAXPUP
c      PARAMETER (MAXPUP=100)
c      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
c      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
c      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
c     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
c     &LPRUP(MAXPUP)
c      SAVE /HEPRUP/
 
C...Read info from file.
c      IF(MSTP(161).GT.0) THEN
c        READ(MSTP(161),*,END=110,ERR=110) IDBMUP(1),IDBMUP(2),EBMUP(1),
c     &  EBMUP(2),PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP
c        DO 100 IPR=1,NPRUP
c          READ(MSTP(161),*,END=110,ERR=110) XSECUP(IPR),XERRUP(IPR),
c     &    XMAXUP(IPR),LPRUP(IPR)
c  100   CONTINUE
c        RETURN
C...Error or prematurely reached end of file.
c  110   WRITE(MSTU(11),5000)
c        STOP
 
C...Else not implemented.
c      ELSE
c        WRITE(MSTU(11),5100)
c        STOP
c      ENDIF
 
C...Format for error printout.
c 5000 FORMAT(1X,'Error: UPINIT routine failed to read information'/
c     &1X,'Execution stopped!')
c 5100 FORMAT(1X,'Error: You have not implemented UPINIT routine'/
c     &1X,'Dummy routine in PYTHIA file called instead.'/
c     &1X,'Execution stopped!')
 
c      RETURN
c      END
 
C*********************************************************************
 
C...UPEVNT
C...Dummy routine, to be replaced by a user implementing external
C...processes. Depending on cross section model chosen, it either has
C...to generate a process of the type IDPRUP requested, or pick a type
C...itself and generate this event. The event is to be stored in the
C...HEPEUP commonblock, including (often) an event weight.

C...New example: handles a standard Les Houches Events File.

      SUBROUTINE UPEVNT
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
C...PYTHIA commonblock: only used to provide read unit MSTP(162).
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /PYPARS/
 
C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/

C...Lines to read in assumed never longer than 200 characters. 
      PARAMETER (MAXLEN=200)
      CHARACTER*(MAXLEN) STRING

C...Format for reading lines.
      CHARACTER*6 STRFMT
      STRFMT='(A000)'
      WRITE(STRFMT(3:5),'(I3)') MAXLEN

C...Loop until finds line beginning with "<event>" or "<event ". 
  100 READ(MSTP(162),STRFMT,END=130,ERR=130) STRING
      IBEG=0
  110 IBEG=IBEG+1
C...Allow indentation.
      IF(STRING(IBEG:IBEG).EQ.' '.AND.IBEG.LT.MAXLEN-6) GOTO 110 
      IF(STRING(IBEG:IBEG+6).NE.'<event>'.AND.
     &STRING(IBEG:IBEG+6).NE.'<event ') GOTO 100

C...Read first line of event info.
      READ(MSTP(162),*,END=130,ERR=130) NUP,IDPRUP,XWGTUP,SCALUP,
     &AQEDUP,AQCDUP

C...Read NUP subsequent lines with information on each particle.
      DO 120 I=1,NUP
        READ(MSTP(162),*,END=130,ERR=130) IDUP(I),ISTUP(I),
     &  MOTHUP(1,I),MOTHUP(2,I),ICOLUP(1,I),ICOLUP(2,I),
     &  (PUP(J,I),J=1,5),VTIMUP(I),SPINUP(I)
  120 CONTINUE
      RETURN

C...Error exit, typically when no more events.
  130 WRITE(*,*) ' Failed to read LHEF event information.'
      WRITE(*,*) ' Will assume end of file has been reached.'
      NUP=0
      MSTI(51)=1
 
      RETURN
      END
