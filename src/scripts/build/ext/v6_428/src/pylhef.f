

C*********************************************************************

C...Combine the two old-style Pythia initialization and event files
C...into a single Les Houches Event File.

      SUBROUTINE PYLHEF
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
C...PYTHIA commonblock: only used to provide read/write units and version.
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /PYPARS/
 
C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &LPRUP(MAXPUP)
      SAVE /HEPRUP/
 
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

C...Rewind initialization and event files. 
      REWIND MSTP(161)
      REWIND MSTP(162)

C...Write header info.
      WRITE(MSTP(163),'(A)') '<LesHouchesEvents version="1.0">'
      WRITE(MSTP(163),'(A)') '<!--'
      WRITE(MSTP(163),'(A,I1,A1,I3)') 'File generated with PYTHIA ',
     &MSTP(181),'.',MSTP(182)
      WRITE(MSTP(163),'(A)') '-->'       

C...Read first line of initialization info and get number of processes.
      READ(MSTP(161),'(A)',END=400,ERR=400) STRING                  
      READ(STRING,*,ERR=400) IDBMUP(1),IDBMUP(2),EBMUP(1),
     &EBMUP(2),PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP

C...Copy initialization lines, omitting trailing blanks. 
C...Embed in <init> ... </init> block.
      WRITE(MSTP(163),'(A)') '<init>' 
      DO 140 IPR=0,NPRUP
        IF(IPR.GT.0) READ(MSTP(161),'(A)',END=400,ERR=400) STRING
        LEN=MAXLEN+1  
  120   LEN=LEN-1
        IF(LEN.GT.1.AND.STRING(LEN:LEN).EQ.' ') GOTO 120
        WRITE(MSTP(163),'(A)',ERR=400) STRING(1:LEN)
  140 CONTINUE
      WRITE(MSTP(163),'(A)') '</init>' 

C...Begin event loop. Read first line of event info or already done.
      READ(MSTP(162),'(A)',END=320,ERR=400) STRING    
  200 CONTINUE

C...Look at first line to know number of particles in event.
      READ(STRING,*,ERR=400) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP

C...Begin an <event> block. Copy event lines, omitting trailing blanks. 
      WRITE(MSTP(163),'(A)') '<event>' 
      DO 240 I=0,NUP
        IF(I.GT.0) READ(MSTP(162),'(A)',END=400,ERR=400) STRING
        LEN=MAXLEN+1  
  220   LEN=LEN-1
        IF(LEN.GT.1.AND.STRING(LEN:LEN).EQ.' ') GOTO 220
        WRITE(MSTP(163),'(A)',ERR=400) STRING(1:LEN)
  240 CONTINUE
              
C...Copy trailing comment lines - with a # in the first column - as is.
  260 READ(MSTP(162),'(A)',END=300,ERR=400) STRING    
      IF(STRING(1:1).EQ.'#') THEN
        LEN=MAXLEN+1  
  280   LEN=LEN-1
        IF(LEN.GT.1.AND.STRING(LEN:LEN).EQ.' ') GOTO 280
        WRITE(MSTP(163),'(A)',ERR=400) STRING(1:LEN)
        GOTO 260
      ENDIF

C..End the <event> block. Loop back to look for next event.
      WRITE(MSTP(163),'(A)') '</event>' 
      GOTO 200

C...Successfully reached end of event loop: write closing tag
C...and remove temporary intermediate files (unless asked not to).
  300 WRITE(MSTP(163),'(A)') '</event>' 
  320 WRITE(MSTP(163),'(A)') '</LesHouchesEvents>' 
      IF(MSTP(164).EQ.1) RETURN
      CLOSE(MSTP(161),ERR=400,STATUS='DELETE')
      CLOSE(MSTP(162),ERR=400,STATUS='DELETE')
      RETURN

C...Error exit.
  400 WRITE(*,*) ' PYLHEF file joining failed!'

      RETURN
      END
