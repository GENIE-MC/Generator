C********************************************************************
C...PYXUED
C... Last change: 
C... 13/01/2009 : H. Przysiezniak Frey, P. Skands
C... Original version:
C... M. El Kacimi
C... 05/07/2005
C     Universal Extra Dimensions Subprocess cross sections  
C     The expressions used are from atl-com-phys-2005-003
C     What is coded here is shat**2/pi * dsigma/dt = |M|**2
C     For each UED subprocess, the color flow used is the same 
C     as the equivalent QCD subprocess. Different configuration
C     color flows are considered to have the same probability. 
C
C     The Xsection is calculated following ATL-PHYS-PUB-2005-003
C     by G.Azuelos and P.H.Beauchemin.
C
C     This routine is called from pysigh.

      SUBROUTINE PYXUED(NCHN,SIGS)

C...Double precision and integer declarations
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...
      INTEGER NGRDEC
      COMMON/DECMOD/NGRDEC
C...
      PARAMETER(KKPART=25,KKFLA=450)
C...Commonblocks
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT3/XSFX(2,-40:40),ISIG(1000,3),SIGH(1000)
      COMMON/PYSGCM/ISUB,ISUBSV,MMIN1,MMAX1,MMIN2,MMAX2,MMINA,MMAXA,
     &KFAC(2,-40:40),COMFAC,FACK,FACA,SH,TH,UH,SH2,TH2,UH2,SQM3,SQM4,
     &SHR,SQPTH,TAUP,BE34,CTH,X(2),SQMZ,SQMW,GMMZ,GMMW,
     &AEM,AS,XW,XW1,XWC,XWV,POLL,POLR,POLLL,POLRR
      SAVE /PYDAT2/,/PYINT1/,/PYINT3/,/PYPARS/
C...UED Pythia common
      COMMON/PYPUED/IUED(0:99),RUED(0:99)
C...Local arrays and complex variables
      DOUBLE PRECISION SHAT,SP,THAT,TP,UHAT,UP,ALPHAS
     + ,FAC1,XMNKK,XMUED,SIGS
      INTEGER NCHN

C...Return if UED not switched on
      IF (IUED(1).LE.0) THEN 
        RETURN 
      ENDIF

C...Energy scale of the parton processus
C...taken equal to the mass of the final state kk
c      Q2=XMNKK**2      

C...Default Mandlestam variable (u/t)hatp=(u/t)hatp-xmnkk**2
      XMNKK=PMAS(KKFLA+23,1) 

C...To compare the cross section with phys-pub-2005-03
C...(no radiative corrections), 
C...take xmnkk=rinv  and q2=rinv**2
c++lnk
C...n.b. (rinv=rued(1))
c      IF(NGRDEC.EQ.1)XMNKK=RUED(0)
      IF(NGRDEC.EQ.1)XMNKK=RUED(1)
c--lnk

      SHAT=VINT(44)
      SP=SHAT
      THAT=VINT(45)
      TP=THAT-XMNKK**2
      UHAT=VINT(46)
      UP=UHAT-XMNKK**2
      BETA34=DSQRT(1.D0-4.D0*XMNKK**2/SHAT)
      PI=DACOS(-1.D0)
c++lnk
c      Q2=RUED(0)**2+(TP*UP-RUED(0)**4)/SP
      Q2=RUED(1)**2+(TP*UP-RUED(1)**4)/SP

c      IF(NGRDEC.EQ.1)Q2=RUED(0)**2
      IF(NGRDEC.EQ.1)Q2=RUED(1)**2
c--lnk

C...Strong coupling value
      ALPHAS=PYALPS(Q2)

      IF(ISUB.EQ.311)THEN
C...gg --> g* g*
         FAC1=9./8.*ALPHAS**2/(SP*TP*UP)**2
         XMUED=FAC1*(XMNKK**4*(6.*TP**4+18.*TP**3*UP+
     &        24.*TP**2*UP**2+18.*TP*UP**3+6.*UP**4)
     &        +XMNKK**2*(6.*TP**4*UP+12.*TP**3*UP**2+
     &        12.*TP**2*UP**3+6*TP*UP**4)
     &        +2.*TP**6+6*TP**5*UP+13*TP**4*UP**2+
     &        15.*TP**3*UP**3+13*TP**2*UP**4+
     &        6.*TP*UP**5+2.*UP**6)
         NCHN=NCHN+1
         ISIG(NCHN,1)=21
         ISIG(NCHN,2)=21
C...Three color flow configurations (qcd g+g->g+g)
         XCOL=PYR(0)
         IF(XCOL.LE.1./3.)THEN
            ISIG(NCHN,3)=1
         ELSEIF(XCOL.LE.2./3.)THEN
            ISIG(NCHN,3)=2
         ELSE
            ISIG(NCHN,3)=3
         ENDIF
         SIGH(NCHN)=COMFAC*XMUED
      ELSEIF(ISUB.EQ.312)THEN
C...q + g -> q*_D + g*, q*_S + g*
C...(the two channels have the same cross section)
         FAC1=-1./36.*ALPHAS**2/(SP*TP*UP)**2
         XMUED=FAC1*(12.*SP*UP**5+5.*SP**2*UP**4+22.*SP**3*UP**3+
     &          5.*SP**4*UP**2+12.*SP**5*UP)
         XMUED=COMFAC*2.*XMUED 

          DO 190 I=MMINA,MMAXA
            IF(I.EQ.0.OR.IABS(I).GT.10) GOTO 190
            DO 180 ISDE=1,2

              IF(ISDE.EQ.1.AND.KFAC(1,I)*KFAC(2,21).EQ.0) GOTO 180
              IF(ISDE.EQ.2.AND.KFAC(1,21)*KFAC(2,I).EQ.0) GOTO 180
              NCHN=NCHN+1
              ISIG(NCHN,ISDE)=I
              ISIG(NCHN,3-ISDE)=21
              ISIG(NCHN,3)=1
              SIGH(NCHN)=XMUED
              IF(PYR(0).GT.0.5)ISIG(NCHN,3)=2
  180       CONTINUE
  190     CONTINUE

      ELSEIF(ISUB.EQ.313)THEN
C...qi + qj -> q*_Di + q*_Dj, q*_Si + q*_Sj 
C...(the two channels have the same cross section)
C...qi and qj have the same charge sign 
         DO 100 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.MSTP(58).OR.KFAC(1,I).EQ.0) GOTO 100
            DO 101 J=MMIN2,MMAX2
               JA=IABS(J)
               IF(J.EQ.0.OR.JA.GT.MSTP(58).OR.KFAC(2,J).
     &           EQ.0) GOTO 101
               IF(J*I.LE.0)GOTO 101
               NCHN=NCHN+1
               ISIG(NCHN,1)=I
               ISIG(NCHN,2)=J
               IF(J.EQ.I)THEN
                  FAC1=1./72.*ALPHAS**2/(TP*UP)**2
                  XMUED=FAC1*
     &                  (XMNKK**2*(8*TP**3+4./3.*TP**2*UP+4./3.*TP*UP**2
     &                 +8.*UP**3)+8.*TP**4+56./3.*TP**3*UP+
     &                 20.*TP**2*UP**2+56./3.*
     &                 TP*UP**3+8.*UP**4)
                  SIGH(NCHN)=COMFAC*2.*XMUED
                  ISIG(NCHN,3)=1
                  IF(PYR(0).GT.0.5)ISIG(NCHN,3)=2
               ELSE
                  FAC1=2./9.*ALPHAS**2/TP**2
                  XMUED=FAC1*(-XMNKK**2*SP+SP**2+0.25*TP**2)     
                  SIGH(NCHN)=COMFAC*2.*XMUED
                  ISIG(NCHN,3)=1
               ENDIF
 101       CONTINUE
 100    CONTINUE
      ELSEIF(ISUB.EQ.314)THEN
C...g + g -> q*_D + q*_Dbar, q*_S + q*_Sbar 
C...(the two channels have the same cross section)
         NCHN=NCHN+1
         ISIG(NCHN,1)=21
         ISIG(NCHN,2)=21
         ISIG(NCHN,3)=INT(1.5+PYR(0))

         FAC1=5./6.*ALPHAS**2/(SP*TP*UP)**2
         XMUED=FAC1*(-XMNKK**4*(8.*TP*UP**3+8.*TP**2*UP**2+8.*TP**3*UP
     +          +4.*UP**4+4*TP**4)
     +          -XMNKK**2*(0.5*TP*UP**4+4.*TP**2*UP**3+15./2.*TP**3
     +          *UP**2+ 4.*TP**4*UP)+TP*UP**5-0.25*TP**2*UP**4+
     +          2.*TP**3*UP**3-0.25*TP**4*UP**2+TP**5*UP)
         
         SIGH(NCHN)=COMFAC*XMUED 
C...has been multiplied by 5: all possible quark flavors in final state

      ELSEIF(ISUB.EQ.315)THEN
C...q + qbar -> q*_D + q*_Dbar, q*_S + q*_Sbar
C...(the two channels have the same cross section)
          DO 141 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 141
            DO 142 J=MMIN2,MMAX2
               IF(J.EQ.0.OR.ABS(I).NE.ABS(J).OR.I*J.GE.0) GOTO 142
               FAC1=2./9.*ALPHAS**2*1./(SP*TP)**2
               XMUED=FAC1*(XMNKK**2*SP*(4.*TP**2-SP*TP-SP**2)+
     &              4.*TP**4+3.*SP*TP**3+11./12.*TP**2*SP**2-
     &              2./3.*SP**3*TP+SP**4)                  
               NCHN=NCHN+1
               ISIG(NCHN,1)=I
               ISIG(NCHN,2)=-I
               ISIG(NCHN,3)=1
               SIGH(NCHN)=COMFAC*2.*XMUED
 142        CONTINUE
 141      CONTINUE
      ELSEIF(ISUB.EQ.316)THEN
C...q + qbar' -> q*_D + q*_Sbar' 
         FAC1=2./9.*ALPHAS**2
         DO 300 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.MSTP(58).OR.KFAC(1,I).EQ.0) GOTO 300
            DO 301 J=MMIN2,MMAX2
               JA=IABS(J)
               IF(J.EQ.0.OR.JA.GT.MSTP(58).OR.KFAC(2,J).EQ.0) GOTO 301
               IF(J*I.GE.0.OR.IA.EQ.JA)GOTO 301
               NCHN=NCHN+1
               ISIG(NCHN,1)=I
               ISIG(NCHN,2)=J
               ISIG(NCHN,3)=1
               FAC1=2./9.*ALPHAS**2/TP**2
               XMUED=FAC1*(-XMNKK**2*SP+SP**2+0.25*TP**2)
               SIGH(NCHN)=COMFAC*XMUED 
 301       CONTINUE
 300   CONTINUE
               
      ELSEIF(ISUB.EQ.317)THEN
C...q + qbar' -> q*_D + q*_Dbar' , q*_S + q*_Sbar' 
C...(the two channels have the same cross section)
         DO 400 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.MSTP(58).OR.KFAC(1,I).EQ.0) GOTO 400     
            DO 401 J=MMIN1,MMAX1
               JA=IABS(J)
               IF(J.EQ.0.OR.JA.GT.MSTP(58).OR.KFAC(2,J).EQ.0) GOTO 401
               IF(J*I.GE.0.OR.IA.EQ.JA)GOTO 401
               NCHN=NCHN+1
               ISIG(NCHN,1)=I
               ISIG(NCHN,2)=J
               ISIG(NCHN,3)=1
               FAC1=1./18.*ALPHAS**2/TP**2
               XMUED=FAC1*(4.*XMNKK**2*SP+4.*SP**2+8.*SP*TP+5*TP**2)  
               SIGH(NCHN)=COMFAC*2.*XMUED 
 401       CONTINUE
 400   CONTINUE
      ELSEIF(ISUB.EQ.318)THEN
C...q + q' -> q*_D + q*_S'
         DO 500 I=MMIN1,MMAX1
            IA=IABS(I)
            IF(I.EQ.0.OR.IA.GT.MSTP(58).OR.KFAC(1,I).EQ.0) GOTO 500   
            DO 501 J=MMIN2,MMAX2
               JA=IABS(J)
               IF(J.EQ.0.OR.JA.GT.MSTP(58).OR.KFAC(2,J).EQ.0) GOTO 501 
               IF(J*I.LE.0)GOTO 501
               IF(IA.EQ.JA)THEN
                  NCHN=NCHN+1
                  ISIG(NCHN,1)=I
                  ISIG(NCHN,2)=J
                  ISIG(NCHN,3)=INT(1.5+PYR(0))
                  FAC1=1./36.*ALPHAS**2/(TP*UP)**2
               XMUED=FAC1*(-8.*XMNKK**2*(TP**3+TP**2*UP+TP*UP**2+UP**3)
     &                 +8.*TP**4+4.*TP**2*UP**2+8.*UP**4)
                  SIGH(NCHN)=COMFAC*XMUED              
               ELSE
                  NCHN=NCHN+1
                  ISIG(NCHN,1)=I
                  ISIG(NCHN,2)=J
                  ISIG(NCHN,3)=1
                  FAC1=1./18.*ALPHAS**2/TP**2
                  XMUED=FAC1*(4.*XMNKK**2*SP+4.*SP**2+8.*SP*TP+5*TP**2)
                  SIGH(NCHN)=COMFAC*2.*XMUED
               ENDIF
 501        CONTINUE
 500     CONTINUE
      ELSEIF(ISUB.EQ.319)THEN
C...q + qbar -> q*_D' +q*_Dbar' , q*_S' + q*_Sbar'
C...(the two channels have the same cross section)
          DO 741 I=MMIN1,MMAX1
            IF(I.EQ.0.OR.IABS(I).GT.MSTP(58).OR.
     &      KFAC(1,I)*KFAC(2,-I).EQ.0) GOTO 741
            DO 742 J=MMIN2,MMAX2
               IF(J.EQ.0.OR.IABS(J).NE.IABS(I).OR.J*I.GT.0) GOTO 742
               FAC1=16./9.*ALPHAS**2*1./(SP)**2
               XMUED=FAC1*(2.*XMNKK**2*SP+SP**2+2.*SP*TP+2.*TP**2)
               NCHN=NCHN+1
               ISIG(NCHN,1)=I
               ISIG(NCHN,2)=-I
               ISIG(NCHN,3)=1
               SIGH(NCHN)=COMFAC*2.*XMUED
 742        CONTINUE
 741      CONTINUE   
       
      ENDIF

      RETURN
      END
