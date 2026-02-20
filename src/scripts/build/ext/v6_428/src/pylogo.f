 
C*********************************************************************
 
C...PYLOGO
C...Writes a logo for the program.
 
      SUBROUTINE PYLOGO
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter for length of information block.
      PARAMETER (IREFER=19)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /PYDAT1/,/PYPARS/
C...Local arrays and character variables.
      INTEGER IDATI(6)
      CHARACTER MONTH(12)*3, LOGO(48)*32, REFER(2*IREFER)*36, LINE*79,
     &VERS*1, SUBV*3, DATE*2, YEAR*4, HOUR*2, MINU*2, SECO*2
 
C...Data on months, logo, titles, and references.
      DATA MONTH/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
     &'Oct','Nov','Dec'/
      DATA (LOGO(J),J=1,19)/
     &'            *......*            ',
     &'       *:::!!:::::::::::*       ',
     &'    *::::::!!::::::::::::::*    ',
     &'  *::::::::!!::::::::::::::::*  ',
     &' *:::::::::!!:::::::::::::::::* ',
     &' *:::::::::!!:::::::::::::::::* ',
     &'  *::::::::!!::::::::::::::::*! ',
     &'    *::::::!!::::::::::::::* !! ',
     &'    !! *:::!!:::::::::::*    !! ',
     &'    !!     !* -><- *         !! ',
     &'    !!     !!                !! ',
     &'    !!     !!                !! ',
     &'    !!                       !! ',
     &'    !!        lh             !! ',
     &'    !!                       !! ',
     &'    !!                 hh    !! ',
     &'    !!    ll                 !! ',
     &'    !!                       !! ',
     &'    !!                          '/
      DATA (LOGO(J),J=20,38)/
     &'Welcome to the Lund Monte Carlo!',
     &'                                ',
     &'PPP  Y   Y TTTTT H   H III   A  ',
     &'P  P  Y Y    T   H   H  I   A A ',
     &'PPP    Y     T   HHHHH  I  AAAAA',
     &'P      Y     T   H   H  I  A   A',
     &'P      Y     T   H   H III A   A',
     &'                                ',
     &'This is PYTHIA version x.xxx    ',
     &'Last date of change: xx xxx 201x',
     &'                                ',
     &'Now is xx xxx 201x at xx:xx:xx  ',
     &'                                ',
     &'Disclaimer: this program comes  ',
     &'without any guarantees. Beware  ',
     &'of errors and use common sense  ',
     &'when interpreting results.      ',
     &'                                ',
     &'Copyright T. Sjostrand (2011)   '/
      DATA (REFER(J),J=1,14)/
     &'An archive of program versions and d',
     &'ocumentation is found on the web:   ',
     &'http://www.thep.lu.se/~torbjorn/Pyth',
     &'ia.html                             ',
     &'                                    ',
     &'                                    ',
     &'When you cite this program, the offi',
     &'cial reference is to the 6.4 manual:',
     &'T. Sjostrand, S. Mrenna and P. Skand',
     &'s, JHEP05 (2006) 026                ',
     &'(LU TP 06-13, FERMILAB-PUB-06-052-CD',
     &'-T) [hep-ph/0603175].               ',
     &'                                    ',
     &'                                    '/
      DATA (REFER(J),J=15,32)/
     &'Also remember that the program, to a',
     &' large extent, represents original  ',
     &'physics research. Other publications',
     &' of special relevance to your       ',
     &'studies may therefore deserve separa',
     &'te mention.                         ',
     &'                                    ',
     &'                                    ',
     &'Main author: Torbjorn Sjostrand; Dep',
     &'artment of Theoretical Physics,     ',
     &'  Lund University, Solvegatan 14A, S',
     &'-223 62 Lund, Sweden;               ',
     &'  phone: + 46 - 46 - 222 48 16; e-ma',
     &'il: torbjorn@thep.lu.se             ',
     &'Author: Stephen Mrenna; Computing Di',
     &'vision, GDS Group,                  ',
     &'  Fermi National Accelerator Laborat',
     &'ory, MS 234, Batavia, IL 60510, USA;'/
      DATA (REFER(J),J=33,2*IREFER)/
     &'  phone: + 1 - 630 - 840 - 2556; e-m',
     &'ail: mrenna@fnal.gov                ',
     &'Author: Peter Skands; CERN/PH-TH, CH',
     &'-1211 Geneva, Switzerland           ',
     &'  phone: + 41 - 22 - 767 24 47; e-ma',
     &'il: peter.skands@cern.ch            '/
 
C...Check that PYDATA linked (check we are in the year 20xx)
      IF(MSTP(183)/100.NE.20) THEN
        WRITE(*,'(1X,A)')
     &  'Error: PYDATA has not been linked.'
        WRITE(*,'(1X,A)') 'Execution stopped!'
        CALL PYSTOP(8)
 
C...Write current version number and current date+time.
      ELSE
        WRITE(VERS,'(I1)') MSTP(181)
        LOGO(28)(24:24)=VERS
        WRITE(SUBV,'(I3)') MSTP(182)
        LOGO(28)(26:28)=SUBV
        IF(MSTP(182).LT.100) LOGO(28)(26:26)='0'
        WRITE(DATE,'(I2)') MSTP(185)
        LOGO(29)(22:23)=DATE
        LOGO(29)(25:27)=MONTH(MSTP(184))
        WRITE(YEAR,'(I4)') MSTP(183)
        LOGO(29)(29:32)=YEAR
        CALL PYTIME(IDATI)
        IF(IDATI(1).LE.0) THEN
          LOGO(31)='                                '
        ELSE
          WRITE(DATE,'(I2)') IDATI(3)
          LOGO(31)(8:9)=DATE
          LOGO(31)(11:13)=MONTH(MAX(1,MIN(12,IDATI(2))))
          WRITE(YEAR,'(I4)') IDATI(1)
          LOGO(31)(15:18)=YEAR
          WRITE(HOUR,'(I2)') IDATI(4)
          LOGO(31)(23:24)=HOUR
          WRITE(MINU,'(I2)') IDATI(5)
          LOGO(31)(26:27)=MINU
          IF(IDATI(5).LT.10) LOGO(31)(26:26)='0'
          WRITE(SECO,'(I2)') IDATI(6)
          LOGO(31)(29:30)=SECO
          IF(IDATI(6).LT.10) LOGO(31)(29:29)='0'
        ENDIF
      ENDIF
 
C...Loop over lines in header. Define page feed and side borders.
      DO 100 ILIN=1,29+IREFER
        LINE=' '
        IF(ILIN.EQ.1) THEN
          LINE(1:1)='1'
        ELSE
          LINE(2:3)='**'
          LINE(78:79)='**'
        ENDIF
 
C...Separator lines and logos.
        IF(ILIN.EQ.2.OR.ILIN.EQ.3.OR.ILIN.GE.28+IREFER) THEN
          LINE(4:77)='***********************************************'//
     &    '***************************'
        ELSEIF(ILIN.GE.6.AND.ILIN.LE.24) THEN
          LINE(6:37)=LOGO(ILIN-5)
          LINE(44:75)=LOGO(ILIN+14)
        ELSEIF(ILIN.GE.26.AND.ILIN.LE.25+IREFER) THEN
          LINE(5:40)=REFER(2*ILIN-51)
          LINE(41:76)=REFER(2*ILIN-50)
        ENDIF
 
C...Write lines to appropriate unit.
        WRITE(MSTU(11),'(A79)') LINE
  100 CONTINUE
 
      RETURN
      END
