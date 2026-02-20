 
C********************************************************************
 
C...PYKFDI
C...Generates a new flavour pair and combines off a hadron
 
      SUBROUTINE PYKFDI(KFL1,KFL2,KFL3,KF)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
C...Local arrays.
      DIMENSION PD(7)
 
      IF(MSTU(123).EQ.0.AND.MSTJ(12).GE.0)  CALL PYKFIN
 
C...Default flavour values. Input consistency checks.
      KF1A=IABS(KFL1)
      KF2A=IABS(KFL2)
      KFL3=0
      KF=0
      IF(KF1A.EQ.0) RETURN
      IF(KF2A.NE.0)THEN
        IF(KF1A.LE.10.AND.KF2A.LE.10.AND.KFL1*KFL2.GT.0) RETURN
        IF(KF1A.GT.10.AND.KF2A.GT.10) RETURN
        IF((KF1A.GT.10.OR.KF2A.GT.10).AND.KFL1*KFL2.LT.0) RETURN
      ENDIF
 
C...Check if tabulated flavour probabilities are to be used.
      IF(MSTJ(15).EQ.1) THEN
        IF(MSTJ(12).GE.5)  CALL PYERRM(29,
     &        '(PYKFDI:) Sorry, option MSTJ(15)=1 not available' //
     &        ' together with MSTJ(12)>=5 modification')
        KTAB1=-1
        IF(KF1A.GE.1.AND.KF1A.LE.6) KTAB1=KF1A
        KFL1A=MOD(KF1A/1000,10)
        KFL1B=MOD(KF1A/100,10)
        KFL1S=MOD(KF1A,10)
        IF(KFL1A.GE.1.AND.KFL1A.LE.4.AND.KFL1B.GE.1.AND.KFL1B.LE.4)
     &  KTAB1=6+KFL1A*(KFL1A-2)+2*KFL1B+(KFL1S-1)/2
        IF(KFL1A.GE.1.AND.KFL1A.LE.4.AND.KFL1A.EQ.KFL1B) KTAB1=KTAB1-1
        IF(KF1A.GE.1.AND.KF1A.LE.6) KFL1A=KF1A
        KTAB2=0
        IF(KF2A.NE.0) THEN
          KTAB2=-1
          IF(KF2A.GE.1.AND.KF2A.LE.6) KTAB2=KF2A
          KFL2A=MOD(KF2A/1000,10)
          KFL2B=MOD(KF2A/100,10)
          KFL2S=MOD(KF2A,10)
          IF(KFL2A.GE.1.AND.KFL2A.LE.4.AND.KFL2B.GE.1.AND.KFL2B.LE.4)
     &    KTAB2=6+KFL2A*(KFL2A-2)+2*KFL2B+(KFL2S-1)/2
          IF(KFL2A.GE.1.AND.KFL2A.LE.4.AND.KFL2A.EQ.KFL2B) KTAB2=KTAB2-1
        ENDIF
        IF(KTAB1.GE.0.AND.KTAB2.GE.0) GOTO 140
      ENDIF
 
C.. Recognize rank 0 diquark case
  100 IRANK=1
      KFDIQ=MAX(KF1A,KF2A)
      IF(KFDIQ.GT.10.AND.KFDIQ.LT.10000) IRANK=0
 
C.. Join two flavours to meson or baryon. Test for popcorn.
      IF(KF2A.GT.0)THEN
        MBARY=0
        IF(KFDIQ.GT.10) THEN
          IF(IRANK.EQ.0.AND.MSTJ(12).LT.5)
     &         CALL PYNMES(KFDIQ)
          IF(MSTU(121).NE.0) THEN
             MSTU(121)=0
             RETURN
          ENDIF
          MBARY=2
        ENDIF
        KFQOLD=KF1A
        KFQVER=KF2A
        GOTO 130
      ENDIF
 
C.. Separate incoming flavours, curtain flavour consistency check
      KFIN=KFL1
      KFQOLD=KF1A
      KFQPOP=KF1A/10000
      IF(KF1A.GT.10)THEN
         KFIN=-KFL1
         KFL1A=MOD(KF1A/1000,10)
         KFL1B=MOD(KF1A/100,10)
         IF(IRANK.EQ.0)THEN
            QAWT=1D0
            IF(KFL1A.GE.3) QAWT=PARF(136+KFL1A/4)
            IF(KFL1B.GE.3) QAWT=QAWT/PARF(136+KFL1B/4)
            KFQPOP=KFL1A+(KFL1B-KFL1A)*INT(1D0/(QAWT+1D0)+PYR(0))
         ENDIF
         IF(KFQPOP.NE.KFL1B.AND.KFQPOP.NE.KFL1A) THEN
             MSTU(121)=0
             RETURN
          ENDIF
         KFQOLD=KFL1A+KFL1B-KFQPOP
      ENDIF
 
C...Meson/baryon choice. Set number of mesons if starting a popcorn
C...system.
  110 MBARY=0
      IF(KF1A.LE.10.AND.MSTJ(12).GT.0)THEN
         IF(MSTU(121).EQ.-1.OR.(1D0+PARJ(1))*PYR(0).GT.1D0)THEN
            MBARY=1
            CALL PYNMES(0)
         ENDIF
      ELSEIF(KF1A.GT.10)THEN
         MBARY=2
         IF(IRANK.EQ.0) CALL PYNMES(KF1A)
         IF(MSTU(121).GT.0) MBARY=-1
      ENDIF
 
C..x->H+q: Choose single vertex quark. Jump to form hadron.
      IF(MBARY.EQ.0.OR.MBARY.EQ.2)THEN
         KFQVER=1+INT((2D0+PARJ(2))*PYR(0))
         KFL3=ISIGN(KFQVER,-KFIN)
         GOTO 130
      ENDIF
 
C..x->H+qq: (IDW=proper PARF position for diquark weights)
      IDW=160
      IF(MBARY.EQ.1)THEN
         IF(MSTU(121).EQ.0) IDW=150
         SQWT=PARF(IDW+1)
         IF(MSTU(121).GT.0) SQWT=SQWT*PARF(135)*PARF(138)**MSTU(121)
         KFQPOP=1+INT((2D0+SQWT)*PYR(0))
C..   Shift to s-curtain parameters if needed
         IF(KFQPOP.GE.3.AND.MSTJ(12).GE.5)THEN
            PARF(194)=PARF(138)*PARF(139)
            PARF(193)=PARJ(8)+PARJ(9)
         ENDIF
      ENDIF
 
C.. x->H+qq: Get vertex quark
      IF(MBARY.EQ.-1.AND.MSTJ(12).GE.5)THEN
         IDW=MSTU(122)
         MSTU(121)=MSTU(121)-1
         IF(IDW.EQ.170) THEN
            IF(MSTU(121).EQ.0)THEN
               IPOS=3*MIN(KFQPOP-1,2)+MIN(KFQOLD-1,2)
            ELSE
               IPOS=3*3+3*MAX(0,MIN(KFQPOP-2,1))+MIN(KFQOLD-1,2)
            ENDIF
         ELSE
            IF(MSTU(121).EQ.0)THEN
               IPOS=3*5+5*MIN(KFQPOP-1,3)+MIN(KFQOLD-1,4)
            ELSE
               IPOS=3*5+5*4+MIN(KFQOLD-1,4)
            ENDIF
         ENDIF
         IPOS=200+30*IPOS+1
 
         IMES=-1
         RMES=PYR(0)*PARF(194)
  120    IMES=IMES+1
         RMES=RMES-PARF(IPOS+IMES)
         IF(IMES.EQ.30) THEN
            MSTU(121)=-1
            KF=-111
            RETURN
         ENDIF
         IF(RMES.GT.0D0) GOTO 120
         KMUL=IMES/5
         KFJ=2*KMUL+1
         IF(KMUL.EQ.2) KFJ=10003
         IF(KMUL.EQ.3) KFJ=10001
         IF(KMUL.EQ.4) KFJ=20003
         IF(KMUL.EQ.5) KFJ=5
         IDIAG=0
         KFQVER=MOD(IMES,5)+1
         IF(KFQVER.GE.KFQOLD) KFQVER=KFQVER+1
         IF(KFQVER.GT.3)THEN
            IDIAG=KFQVER-3
            KFQVER=KFQOLD
         ENDIF
      ELSE
         IF(MBARY.EQ.-1) IDW=170
         SQWT=PARF(IDW+2)
         IF(KFQPOP.EQ.3) SQWT=PARF(IDW+3)
         IF(KFQPOP.GT.3) SQWT=PARF(IDW+3)*(1D0/PARF(IDW+5)+1D0)/2D0
         KFQVER=MIN(3,1+INT((2D0+SQWT)*PYR(0)))
         IF(KFQPOP.LT.3.AND.KFQVER.LT.3)THEN
            KFQVER=KFQPOP
            IF(PYR(0).GT.PARF(IDW+4)) KFQVER=3-KFQPOP
         ENDIF
      ENDIF
 
C..x->H+qq: form outgoing diquark with KFQPOP flag at 10000-pos
      KFLDS=3
      IF(KFQPOP.NE.KFQVER)THEN
         SWT=PARF(IDW+7)
         IF(KFQVER.EQ.3) SWT=PARF(IDW+6)
         IF(KFQPOP.GE.3) SWT=PARF(IDW+5)
         IF((1D0+SWT)*PYR(0).LT.1D0) KFLDS=1
      ENDIF
      KFDIQ=900*MAX(KFQVER,KFQPOP)+100*(KFQVER+KFQPOP)+KFLDS
     &      +10000*KFQPOP
      KFL3=ISIGN(KFDIQ,KFIN)
 
C..x->M+y: flavour for meson.
  130 IF(MBARY.LE.0)THEN
        KFLA=MAX(KFQOLD,KFQVER)
        KFLB=MIN(KFQOLD,KFQVER)
        KFS=ISIGN(1,KFL1)
        IF(KFLA.NE.KFQOLD) KFS=-KFS
C... Form meson, with spin and flavour mixing for diagonal states.
        IF(MBARY.EQ.-1.AND.MSTJ(12).GE.5)THEN
           IF(IDIAG.GT.0) KF=110*IDIAG+KFJ
           IF(IDIAG.EQ.0) KF=(100*KFLA+10*KFLB+KFJ)*KFS*(-1)**KFLA
           RETURN
        ENDIF
        IF(KFLA.LE.2) KMUL=INT(PARJ(11)+PYR(0))
        IF(KFLA.EQ.3) KMUL=INT(PARJ(12)+PYR(0))
        IF(KFLA.GE.4) KMUL=INT(PARJ(13)+PYR(0))
        IF(KMUL.EQ.0.AND.PARJ(14).GT.0D0)THEN
          IF(PYR(0).LT.PARJ(14)) KMUL=2
        ELSEIF(KMUL.EQ.1.AND.PARJ(15)+PARJ(16)+PARJ(17).GT.0D0)THEN
          RMUL=PYR(0)
          IF(RMUL.LT.PARJ(15)) KMUL=3
          IF(KMUL.EQ.1.AND.RMUL.LT.PARJ(15)+PARJ(16)) KMUL=4
          IF(KMUL.EQ.1.AND.RMUL.LT.PARJ(15)+PARJ(16)+PARJ(17)) KMUL=5
        ENDIF
        KFLS=3
        IF(KMUL.EQ.0.OR.KMUL.EQ.3) KFLS=1
        IF(KMUL.EQ.5) KFLS=5
        IF(KFLA.NE.KFLB)THEN
          KF=(100*KFLA+10*KFLB+KFLS)*KFS*(-1)**KFLA
        ELSE
          RMIX=PYR(0)
          IMIX=2*KFLA+10*KMUL
          IF(KFLA.LE.3) KF=110*(1+INT(RMIX+PARF(IMIX-1))+
     &    INT(RMIX+PARF(IMIX)))+KFLS
          IF(KFLA.GE.4) KF=110*KFLA+KFLS
        ENDIF
        IF(KMUL.EQ.2.OR.KMUL.EQ.3) KF=KF+ISIGN(10000,KF)
        IF(KMUL.EQ.4) KF=KF+ISIGN(20000,KF)
 
C..Optional extra suppression of eta and eta'.
C..Allow shift to qq->B+q in old version (set IRANK to 0)
        IF(KF.EQ.221.OR.KF.EQ.331)THEN
           IF(PYR(0).GT.PARJ(25+KF/300))THEN
              IF(KF2A.GT.0) GOTO 130
              IF(MSTJ(12).LT.4) IRANK=0
              GOTO 110
           ENDIF
        ENDIF
        MSTU(121)=0
 
C.. x->B+y: Flavour for baryon
      ELSE
        KFLA=KFQVER
        IF(KF1A.LE.10) KFLA=KFQOLD
        KFLB=MOD(KFDIQ/1000,10)
        KFLC=MOD(KFDIQ/100,10)
        KFLDS=MOD(KFDIQ,10)
        KFLD=MAX(KFLA,KFLB,KFLC)
        KFLF=MIN(KFLA,KFLB,KFLC)
        KFLE=KFLA+KFLB+KFLC-KFLD-KFLF
 
C...  SU(6) factors for formation of baryon.
        KBARY=3
        KDMAX=5
        KFLG=KFLB
        IF(KFLB.NE.KFLC)THEN
           KBARY=2*KFLDS-1
           KDMAX=1+KFLDS/2
           IF(KFLB.GT.2) KDMAX=KDMAX+2
        ENDIF
        IF(KFLA.NE.KFLB.AND.KFLA.NE.KFLC)THEN
           KBARY=KBARY+1
           KFLG=KFLA
        ENDIF
 
        SU6MAX=PARF(140+KDMAX)
        SU6DEC=PARJ(18)
        SU6S  =PARF(146)
        IF(MSTJ(12).GE.5.AND.IRANK.EQ.0) THEN
           SU6MAX=1D0
           SU6DEC=1D0
           SU6S  =1D0
        ENDIF
        SU6OCT=PARF(60+KBARY)
        IF(KFLG.GT.MAX(KFLA+KFLB-KFLG,2))THEN
           SU6OCT=SU6OCT*4*SU6S/(3*SU6S+1)
           IF(KBARY.EQ.2) SU6OCT=PARF(60+KBARY)*4/(3*SU6S+1)
        ELSE
           IF(KBARY.EQ.6) SU6OCT=SU6OCT*(3+SU6S)/(3*SU6S+1)
        ENDIF
        SU6WT=SU6OCT+SU6DEC*PARF(70+KBARY)
 
C..   SU(6) test. Old options enforce new baryon if q->B+qq is rejected.
        IF(SU6WT.LT.PYR(0)*SU6MAX.AND.KF2A.EQ.0)THEN
           MSTU(121)=0
           IF(MSTJ(12).LE.2.AND.MBARY.EQ.1) MSTU(121)=-1
           GOTO 110
        ENDIF
 
C.. Form baryon. Distinguish Lambda- and Sigmalike baryons.
        KSIG=1
        KFLS=2
        IF(SU6WT*PYR(0).GT.SU6OCT) KFLS=4
        IF(KFLS.EQ.2.AND.KFLD.GT.KFLE.AND.KFLE.GT.KFLF)THEN
          KSIG=KFLDS/3
          IF(KFLA.NE.KFLD) KSIG=INT(3*SU6S/(3*SU6S+KFLDS**2)+PYR(0))
        ENDIF
        KF=ISIGN(1000*KFLD+100*KFLE+10*KFLF+KFLS,KFL1)
        IF(KSIG.EQ.0) KF=ISIGN(1000*KFLD+100*KFLF+10*KFLE+KFLS,KFL1)
      ENDIF
      RETURN
 
C...Use tabulated probabilities to select new flavour and hadron.
  140 IF(KTAB2.EQ.0.AND.MSTJ(12).LE.0) THEN
        KT3L=1
        KT3U=6
      ELSEIF(KTAB2.EQ.0.AND.KTAB1.GE.7.AND.MSTJ(12).LE.1) THEN
        KT3L=1
        KT3U=6
      ELSEIF(KTAB2.EQ.0) THEN
        KT3L=1
        KT3U=22
      ELSE
        KT3L=KTAB2
        KT3U=KTAB2
      ENDIF
      RFL=0D0
      DO 160 KTS=0,2
        DO 150 KT3=KT3L,KT3U
          RFL=RFL+PARF(120+80*KTAB1+25*KTS+KT3)
  150   CONTINUE
  160 CONTINUE
      RFL=PYR(0)*RFL
      DO 180 KTS=0,2
        KTABS=KTS
        DO 170 KT3=KT3L,KT3U
          KTAB3=KT3
          RFL=RFL-PARF(120+80*KTAB1+25*KTS+KT3)
          IF(RFL.LE.0D0) GOTO 190
  170   CONTINUE
  180 CONTINUE
  190 CONTINUE
 
C...Reconstruct flavour of produced quark/diquark.
      IF(KTAB3.LE.6) THEN
        KFL3A=KTAB3
        KFL3B=0
        KFL3=ISIGN(KFL3A,KFL1*(2*KTAB1-13))
      ELSE
        KFL3A=1
        IF(KTAB3.GE.8) KFL3A=2
        IF(KTAB3.GE.11) KFL3A=3
        IF(KTAB3.GE.16) KFL3A=4
        KFL3B=(KTAB3-6-KFL3A*(KFL3A-2))/2
        KFL3=1000*KFL3A+100*KFL3B+1
        IF(KFL3A.EQ.KFL3B.OR.KTAB3.NE.6+KFL3A*(KFL3A-2)+2*KFL3B) KFL3=
     &  KFL3+2
        KFL3=ISIGN(KFL3,KFL1*(13-2*KTAB1))
      ENDIF
 
C...Reconstruct meson code.
      IF(KFL3A.EQ.KFL1A.AND.KFL3B.EQ.KFL1B.AND.(KFL3A.LE.3.OR.
     &KFL3B.NE.0)) THEN
        RFL=PYR(0)*(PARF(143+80*KTAB1+25*KTABS)+PARF(144+80*KTAB1+
     &  25*KTABS)+PARF(145+80*KTAB1+25*KTABS))
        KF=110+2*KTABS+1
        IF(RFL.GT.PARF(143+80*KTAB1+25*KTABS)) KF=220+2*KTABS+1
        IF(RFL.GT.PARF(143+80*KTAB1+25*KTABS)+PARF(144+80*KTAB1+
     &  25*KTABS)) KF=330+2*KTABS+1
      ELSEIF(KTAB1.LE.6.AND.KTAB3.LE.6) THEN
        KFLA=MAX(KTAB1,KTAB3)
        KFLB=MIN(KTAB1,KTAB3)
        KFS=ISIGN(1,KFL1)
        IF(KFLA.NE.KF1A) KFS=-KFS
        KF=(100*KFLA+10*KFLB+2*KTABS+1)*KFS*(-1)**KFLA
      ELSEIF(KTAB1.GE.7.AND.KTAB3.GE.7) THEN
        KFS=ISIGN(1,KFL1)
        IF(KFL1A.EQ.KFL3A) THEN
          KFLA=MAX(KFL1B,KFL3B)
          KFLB=MIN(KFL1B,KFL3B)
          IF(KFLA.NE.KFL1B) KFS=-KFS
        ELSEIF(KFL1A.EQ.KFL3B) THEN
          KFLA=KFL3A
          KFLB=KFL1B
          KFS=-KFS
        ELSEIF(KFL1B.EQ.KFL3A) THEN
          KFLA=KFL1A
          KFLB=KFL3B
        ELSEIF(KFL1B.EQ.KFL3B) THEN
          KFLA=MAX(KFL1A,KFL3A)
          KFLB=MIN(KFL1A,KFL3A)
          IF(KFLA.NE.KFL1A) KFS=-KFS
        ELSE
          CALL PYERRM(2,'(PYKFDI:) no matching flavours for qq -> qq')
          GOTO 100
        ENDIF
        KF=(100*KFLA+10*KFLB+2*KTABS+1)*KFS*(-1)**KFLA
 
C...Reconstruct baryon code.
      ELSE
        IF(KTAB1.GE.7) THEN
          KFLA=KFL3A
          KFLB=KFL1A
          KFLC=KFL1B
        ELSE
          KFLA=KFL1A
          KFLB=KFL3A
          KFLC=KFL3B
        ENDIF
        KFLD=MAX(KFLA,KFLB,KFLC)
        KFLF=MIN(KFLA,KFLB,KFLC)
        KFLE=KFLA+KFLB+KFLC-KFLD-KFLF
        IF(KTABS.EQ.0) KF=ISIGN(1000*KFLD+100*KFLF+10*KFLE+2,KFL1)
        IF(KTABS.GE.1) KF=ISIGN(1000*KFLD+100*KFLE+10*KFLF+2*KTABS,KFL1)
      ENDIF
 
C...Check that constructed flavour code is an allowed one.
      IF(KFL2.NE.0) KFL3=0
      KC=PYCOMP(KF)
      IF(KC.EQ.0) THEN
        CALL PYERRM(2,'(PYKFDI:) user-defined flavour probabilities '//
     &  'failed')
        GOTO 100
      ENDIF
 
      RETURN
      END
