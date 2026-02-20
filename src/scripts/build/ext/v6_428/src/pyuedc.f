C*********************************************************************
 
C...PYUEDC
C...Auxiliary to PYXDIN
C...Mass kk states radiative corrections 
C...Radiative corrections are included (hep/ph0204342)

      SUBROUTINE PYUEDC

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP

      PARAMETER(KKPART=25,KKFLA=450)

C...UED Pythia common
      COMMON/PYPUED/IUED(0:99),RUED(0:99)
C...Pythia common: particles properties
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)      
C...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
C...Resonance width and secondary decay treatment.
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

C...Local variables
      DOUBLE PRECISION PI,QUP,QDW
      DOUBLE PRECISION WDTP,WDTE
      DIMENSION WDTP(0:400),WDTE(0:400,0:5)
      DOUBLE PRECISION Q2,ALPHEM,ALPHS,SW2,CW2,RMKK,RMKK2,ZETA3
      DOUBLE PRECISION DSMG2,LOGLAM,DBMG2
      DOUBLE PRECISION DBMQU,DBMQD,DBMQDO,DBMLDO,DBMLE
      DOUBLE PRECISION DSMA2,DSMB2,DBMA2,DBMB2
      DOUBLE PRECISION RFACT,RMW,RMZ,RMZ2,RMW2,A,B,C,SQRDEL,DMB2,DMA2
      DOUBLE PRECISION SWW1,CWW1
      DOUBLE PRECISION RMGST,RMPHST,RMZST,RMWST
      DOUBLE PRECISION RMDQST,RMSQUS,RMSQDS,RMLSLD,RMLSLE
      DOUBLE PRECISION SW21,CW21,SW021,CW021
      COMMON/SW1/SW021,CW021
C...UED related declarations:
C...equivalences between ordered particles (451->475)
C...and UED particle code (5 000 000 + id)
      DIMENSION IUEDEQ(475)
      DATA (IUEDEQ(I),I=451,475)/
C...Singlet quarks      
     & 6100001,6100002,6100003,6100004,6100005,6100006,
C...Doublet quarks
     & 5100001,5100002,5100003,5100004,5100005,5100006, 
C...Singlet leptons
     & 6100011,6100013,6100015,                         
C...Doublet leptons
     & 5100012,5100011,5100014,5100013,5100016,5100015,
C...Gauge boson KK excitations
     & 5100021,5100022,5100023,5100024/                 

C...N.B. rinv=rued(1)
      IF(RUED(1).LE.0.)THEN
         WRITE(MSTU(11),*) 'PYUEDC: RINV < 0 : ',RUED(1)
         WRITE(MSTU(11),*) 'DEFAULT KK STATE MASSES ARE TAKEN '
         RETURN
      ENDIF

      PI=DACOS(-1.D0)
      RMZ  = PMAS(23,1)
      RMZ2 = RMZ**2
      RMW  = PMAS(24,1)
      RMW2 = RMW**2
      ALPHEM = PARU(101)
      QUP = 2./3.
      QDW = -1./3.

c...qt is q-tilde, qs is q-star
c...strong coupling value
      Q2 = RUED(1)**2
      ALPHS=PYALPS(Q2)
      
c...weak mixing angle
      SW2=PARU(102)
      CW2=1D0-PARU(102)
      
c...for the mass corrections
      RMKK = RUED(1)
      RMKK2 = RMKK**2
      ZETA3= 1.2
      
C... Either fix the cutoff scale LAMUED
      IF(IUED(5).EQ.0)THEN
         LOGLAM = DLOG((RUED(3)*(1./RUED(1)))**2)
C... or the ratio LAMUED/RINV (=product Lambda*R)
      ELSEIF(IUED(5).EQ.1)THEN
         LOGLAM = DLOG(RUED(4)**2)
      ELSE
         WRITE(MSTU(11),*) '(PYUEDC:) INVALID VALUE FOR IUED(5)'
         CALL PYSTOP(6000)
      ENDIF

C...Calculate the radiative corrections for the UED KK masses
      IF(IUED(6).EQ.1)THEN
         RFACT=1.D0
C...or induce a minute mass difference
C...keeping the UED KK mass values nearly equal to 1/R
      ELSEIF(IUED(6).EQ.0)THEN
         RFACT=0.01D0
      ELSE
         WRITE(MSTU(11),*) '(PYUEDC:) INVALID VALUE FOR IUED(6)'
         CALL PYSTOP(6001)
      ENDIF

c...Take into account only the strong interactions:

c...The space bulk corrections :
      DSMG2 = RMKK2*(-1.5)*(ALPHS/4./PI)*ZETA3/PI**2
c...The boundary terms:
      DBMG2 = RMKK2*(23./2.)*(ALPHS/4./PI)*LOGLAM

c...Mass corrections for fermions are extracted from 
c...Phys. Rev. D66 036005(2002)9
      DBMQDO=RMKK*(3.*(ALPHS/4./PI)+27./16.*(ALPHEM/4./PI/SW2)
     .     +1./16.*(ALPHEM/4./PI/CW2))*LOGLAM
      DBMQU=RMKK*(3.*(ALPHS/4./PI)
     .     +(ALPHEM/4./PI/CW2))*LOGLAM
      DBMQD=RMKK*(3.*(ALPHS/4./PI)
     .     +0.25*(ALPHEM/4./PI/CW2))*LOGLAM
      
      DBMLDO=RMKK *((27./16.)*(ALPHEM/4./PI/SW2)+9./16.*
     .     (ALPHEM/4./PI/CW2))*LOGLAM
      DBMLE=RMKK *(9./4.*(ALPHEM/4./PI/CW2))*LOGLAM
      
c...Vector boson masss matrix diagonalization
      DBMB2 = RMKK2*(-1./6.)*(ALPHEM/4./PI/CW2)*LOGLAM
      DSMB2 = RMKK2*(-39./2.)*(ALPHEM/4./PI**3/CW2)*ZETA3
      DBMA2 = RMKK2*(15./2.)*(ALPHEM/4./PI/SW2)*LOGLAM
      DSMA2 = RMKK2*(-5./2.)*(ALPHEM/4./PI**3/SW2)*ZETA3
      
c...Elements of the mass matrix
      A = RMZ2*SW2 + DBMB2 + DSMB2
      B = RMZ2*CW2 + DBMA2 + DSMA2
      C = RMZ2*DSQRT(SW2*CW2)
      SQRDEL = DSQRT( (A-B)**2 + 4*C**2 )

c...Eigenvalues: corrections to X1 and Z1 masses
      DMB2 = (A+B-SQRDEL)/2. 
      DMA2 = (A+B+SQRDEL)/2. 
      
c...Rotation angles	
      SWW1 = 2*C
      CWW1 = A-B-SQRDEL
C...Weinberg angle
      SW21= SWW1**2/(SWW1**2 + CWW1**2)
      CW21= 1. - SW21
      
      SW021=SW21
      CW021=CW21
      
c...Masses:
      RMGST = RMKK+RFACT*(DSQRT(RMKK2 + DSMG2 + DBMG2)-RMKK)
      
      RMDQST=RMKK+RFACT*DBMQDO
      RMSQUS=RMKK+RFACT*DBMQU
      RMSQDS=RMKK+RFACT*DBMQD

C...Note: MZ mass is included in ma2
      RMPHST= RMKK+RFACT*(DSQRT(RMKK2 + DMB2)-RMKK)
      RMZST = RMKK+RFACT*(DSQRT(RMKK2 + DMA2)-RMKK)
      RMWST = RMKK+RFACT*(DSQRT(RMKK2 + DBMA2 + DSMA2 + RMW**2)-RMKK)

      RMLSLD=RMKK+RFACT*DBMLDO
      RMLSLE=RMKK+RFACT*DBMLE

      DO 100 IPART=1,5,2
        PMAS(KKFLA+IPART,1)=RMSQDS
 100  CONTINUE
      DO 110 IPART=2,6,2
        PMAS(KKFLA+IPART,1)=RMSQUS
 110  CONTINUE
      DO 120 IPART=7,12
        PMAS(KKFLA+IPART,1)=RMDQST
 120  CONTINUE
      DO 130 IPART=13,15
        PMAS(KKFLA+IPART,1)=RMLSLE
 130  CONTINUE
      DO 140 IPART=16,21
        PMAS(KKFLA+IPART,1)=RMLSLD
 140  CONTINUE
      PMAS(KKFLA+22,1)=RMGST
      PMAS(KKFLA+23,1)=RMPHST
      PMAS(KKFLA+24,1)=RMZST
      PMAS(KKFLA+25,1)=RMWST

      WRITE(MSTU(11),7000) ' PYUEDC: ',
     & 'UED Mass Spectrum (GeV) :'
      WRITE(MSTU(11),7100) '   m(d*_S,s*_S,b*_S) = ',RMSQDS
      WRITE(MSTU(11),7100) '   m(u*_S,c*_S,t*_S) = ',RMSQUS
      WRITE(MSTU(11),7100) '   m(q*_D)           = ',RMDQST
      WRITE(MSTU(11),7100) '   m(l*_S)           = ',RMLSLE
      WRITE(MSTU(11),7100) '   m(l*_D)           = ',RMLSLD
      WRITE(MSTU(11),7100) '   m(g*)             = ',RMGST
      WRITE(MSTU(11),7100) '   m(gamma*)         = ',RMPHST
      WRITE(MSTU(11),7100) '   m(Z*)             = ',RMZST
      WRITE(MSTU(11),7100) '   m(W*)             = ',RMWST
      WRITE(MSTU(11),7000) ' '

C...Initialize widths, branching ratios and life time
      DO 199 IPART=1,25
        KC=KKFLA+IPART
        IF(MWID(KC).EQ.1.AND.MDCY(KC,1).EQ.1)THEN
          CALL PYWIDT(IUEDEQ(KC),PMAS(KC,1)**2,WDTP,WDTE)
          IF(WDTP(0).LE.0)THEN
             WRITE(MSTU(11),*) 
     +             'PYUEDC WARNING: TOTAL WIDTH = 0 --> KC ', KC
             WRITE(MSTU(11),*) 'INITIAL VALUE IS TAKEN',PMAS(KC,2)
             GOTO 199
          ELSE
            DO 180 IDC=1,MDCY(KC,3)
              IC=IDC+MDCY(KC,2)-1
              IF(MDME(IC,1).EQ.1.AND.WDTP(IDC).GT.0.)THEN
C...Life time in cm^{-1}.  paru(3) gev^{-1} -> fm
                PMAS(KC,4)=PARU(3)/WDTP(IDC)*1.D-12
                BRAT(IC)=WDTP(IDC)/WDTP(0)
              ENDIF
 180        CONTINUE
          ENDIF
        ENDIF
 199  CONTINUE

C...Format to use for comments
 7000 FORMAT(' * ',A)
 7100 FORMAT(' * ',A,F12.3)

      END
