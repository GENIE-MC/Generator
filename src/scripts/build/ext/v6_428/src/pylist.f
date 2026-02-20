 
C*********************************************************************
 
C...PYLIST
C...Gives program heading, or lists an event, or particle
C...data, or current parameter values.
 
      SUBROUTINE PYLIST(MLIST)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
 
C...HEPEVT commonblock.
      PARAMETER (NMXHEP=4000)
      COMMON/HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
     &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
      DOUBLE PRECISION PHEP,VHEP
      SAVE /HEPEVT/
 
C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &VTIMUP(MAXNUP),SPINUP(MAXNUP)
      SAVE /HEPEUP/
 
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYCTAG/NCT,MCT(4000,2)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYCTAG/
C...Local arrays, character variables and data.
      CHARACTER CHAP*16,CHAC*16,CHAN*16,CHAD(5)*16,CHDL(7)*4
      DIMENSION PS(6)
      DATA CHDL/'(())',' ','()','!!','<>','==','(==)'/
 
C...Initialization printout: version number and date of last change.
      IF(MLIST.EQ.0.OR.MSTU(12).EQ.1) THEN
        CALL PYLOGO
        MSTU(12)=12345
        IF(MLIST.EQ.0) RETURN
      ENDIF
 
C...List event data, including additional lines after N.
      IF(MLIST.GE.1.AND.MLIST.LE.4) THEN
        IF(MLIST.EQ.1) WRITE(MSTU(11),5100)
        IF(MLIST.EQ.2) WRITE(MSTU(11),5200)
        IF(MLIST.EQ.3) WRITE(MSTU(11),5300)
        IF(MLIST.EQ.4) WRITE(MSTU(11),5400)
        LMX=12
        IF(MLIST.GE.2) LMX=16
        ISTR=0
        IMAX=N
        IF(MSTU(2).GT.0) IMAX=MSTU(2)
        DO 120 I=MAX(1,MSTU(1)),MAX(IMAX,N+MAX(0,MSTU(3)))
          IF(I.GT.IMAX.AND.I.LE.N) GOTO 120
          IF(MSTU(15).EQ.0.AND.K(I,1).LE.0) GOTO 120
          IF(MSTU(15).EQ.1.AND.K(I,1).LT.0) GOTO 120
 
C...Get particle name, pad it and check it is not too long.
          CALL PYNAME(K(I,2),CHAP)
          LEN=0
          DO 100 LEM=1,16
            IF(CHAP(LEM:LEM).NE.' ') LEN=LEM
  100     CONTINUE
          MDL=(K(I,1)+19)/10
          LDL=0
          IF(MDL.EQ.2.OR.MDL.GE.8) THEN
            CHAC=CHAP
            IF(LEN.GT.LMX) CHAC(LMX:LMX)='?'
          ELSE
            LDL=1
            IF(MDL.EQ.1.OR.MDL.EQ.7) LDL=2
            IF(LEN.EQ.0) THEN
              CHAC=CHDL(MDL)(1:2*LDL)//' '
            ELSE
              CHAC=CHDL(MDL)(1:LDL)//CHAP(1:MIN(LEN,LMX-2*LDL))//
     &        CHDL(MDL)(LDL+1:2*LDL)//' '
              IF(LEN+2*LDL.GT.LMX) CHAC(LMX:LMX)='?'
            ENDIF
          ENDIF
 
C...Add information on string connection.
          IF(K(I,1).EQ.1.OR.K(I,1).EQ.2.OR.K(I,1).EQ.11.OR.K(I,1).EQ.12)
     &    THEN
            KC=PYCOMP(K(I,2))
            KCC=0
            IF(KC.NE.0) KCC=KCHG(KC,2)
            IF(IABS(K(I,2)).EQ.39) THEN
              IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='X'
            ELSEIF(KCC.NE.0.AND.ISTR.EQ.0) THEN
              ISTR=1
              IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='A'
            ELSEIF(KCC.NE.0.AND.(K(I,1).EQ.2.OR.K(I,1).EQ.12)) THEN
              IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='I'
            ELSEIF(KCC.NE.0) THEN
              ISTR=0
              IF(LEN+2*LDL+3.LE.LMX) CHAC(LMX-1:LMX-1)='V'
            ENDIF
          ENDIF
          IF((K(I,1).EQ.41.OR.K(I,1).EQ.51).AND.LEN+2*LDL+3.LE.LMX)
     &    CHAC(LMX-1:LMX-1)='I'
 
C...Write data for particle/jet.
          IF(MLIST.EQ.1.AND.ABS(P(I,4)).LT.9999D0) THEN
            WRITE(MSTU(11),5500) I,CHAC(1:12),(K(I,J1),J1=1,3),
     &      (P(I,J2),J2=1,5)
          ELSEIF(MLIST.EQ.1.AND.ABS(P(I,4)).LT.99999D0) THEN
            WRITE(MSTU(11),5600) I,CHAC(1:12),(K(I,J1),J1=1,3),
     &      (P(I,J2),J2=1,5)
          ELSEIF(MLIST.EQ.1) THEN
            WRITE(MSTU(11),5700) I,CHAC(1:12),(K(I,J1),J1=1,3),
     &      (P(I,J2),J2=1,5)
          ELSEIF(MSTU(5).EQ.10000.AND.(K(I,1).EQ.3.OR.K(I,1).EQ.13.OR.
     &      K(I,1).EQ.14.OR.K(I,1).EQ.42.OR.K(I,1).EQ.52)) THEN
            IF(MLIST.NE.4) WRITE(MSTU(11),5800) I,CHAC,(K(I,J1),J1=1,3),
     &      K(I,4)/100000000,MOD(K(I,4)/10000,10000),MOD(K(I,4),10000),
     &      K(I,5)/100000000,MOD(K(I,5)/10000,10000),MOD(K(I,5),10000),
     &      (P(I,J2),J2=1,5)
            IF(MLIST.EQ.4) WRITE(MSTU(11),5900) I,CHAC,(K(I,J1),J1=1,3),
     &      K(I,4)/100000000,MOD(K(I,4)/10000,10000),MOD(K(I,4),10000),
     &           K(I,5)/100000000,MOD(K(I,5)/10000,10000),MOD(K(I,5)
     &           ,10000),MCT(I,1),MCT(I,2)
          ELSE
            IF(MLIST.NE.4) WRITE(MSTU(11),6000) I,CHAC,(K(I,J1),J1=1,5),
     &      (P(I,J2),J2=1,5)
            IF(MLIST.EQ.4) WRITE(MSTU(11),6100) I,CHAC,(K(I,J1),J1=1,5)
     &           ,MCT(I,1),MCT(I,2)
          ENDIF
          IF(MLIST.EQ.3) WRITE(MSTU(11),6200) (V(I,J),J=1,5)
 
C...Insert extra separator lines specified by user.
          IF(MSTU(70).GE.1) THEN
            ISEP=0
            DO 110 J=1,MIN(10,MSTU(70))
              IF(I.EQ.MSTU(70+J)) ISEP=1
  110       CONTINUE
            IF(ISEP.EQ.1) THEN
              IF(MLIST.EQ.1) WRITE(MSTU(11),6300)
              IF(MLIST.EQ.2.OR.MLIST.EQ.3) WRITE(MSTU(11),6400)
              IF(MLIST.EQ.4) WRITE(MSTU(11),6500)
            ENDIF
          ENDIF
  120   CONTINUE
 
C...Sum of charges and momenta.
        DO 130 J=1,6
          PS(J)=PYP(0,J)
  130   CONTINUE
        IF(MLIST.EQ.1.AND.ABS(PS(4)).LT.9999D0) THEN
          WRITE(MSTU(11),6600) PS(6),(PS(J),J=1,5)
        ELSEIF(MLIST.EQ.1.AND.ABS(PS(4)).LT.99999D0) THEN
          WRITE(MSTU(11),6700) PS(6),(PS(J),J=1,5)
        ELSEIF(MLIST.EQ.1) THEN
          WRITE(MSTU(11),6800) PS(6),(PS(J),J=1,5)
        ELSEIF(MLIST.LE.3) THEN
          WRITE(MSTU(11),6900) PS(6),(PS(J),J=1,5)
        ELSE
          WRITE(MSTU(11),7000) PS(6)
        ENDIF
 
C...Simple listing of HEPEVT entries (mainly for test purposes).
      ELSEIF(MLIST.EQ.5) THEN
        WRITE(MSTU(11),7100)
        DO 140 I=1,NHEP
          IF(ISTHEP(I).EQ.0) GOTO 140
          WRITE(MSTU(11),7200) I,ISTHEP(I),IDHEP(I),JMOHEP(1,I),
     &    JMOHEP(2,I),JDAHEP(1,I),JDAHEP(2,I),(PHEP(J,I),J=1,5)
  140   CONTINUE
 
 
C...Simple listing of user-process entries (mainly for test purposes).
      ELSEIF(MLIST.EQ.7) THEN
        WRITE(MSTU(11),7300)
        DO 150 I=1,NUP
          WRITE(MSTU(11),7400) I,ISTUP(I),IDUP(I),MOTHUP(1,I),
     &    MOTHUP(2,I),ICOLUP(1,I),ICOLUP(2,I),(PUP(J,I),J=1,5)
  150   CONTINUE
 
C...Give simple list of KF codes defined in program.
      ELSEIF(MLIST.EQ.11) THEN
        WRITE(MSTU(11),7500)
        DO 160 KF=1,80
          CALL PYNAME(KF,CHAP)
          CALL PYNAME(-KF,CHAN)
          IF(CHAP.NE.' '.AND.CHAN.EQ.' ') WRITE(MSTU(11),7600) KF,CHAP
          IF(CHAN.NE.' ') WRITE(MSTU(11),7600) KF,CHAP,-KF,CHAN
  160   CONTINUE
        DO 190 KFLS=1,3,2
          DO 180 KFLA=1,5
            DO 170 KFLB=1,KFLA-(3-KFLS)/2
              KF=1000*KFLA+100*KFLB+KFLS
              CALL PYNAME(KF,CHAP)
              CALL PYNAME(-KF,CHAN)
              WRITE(MSTU(11),7600) KF,CHAP,-KF,CHAN
  170       CONTINUE
  180     CONTINUE
  190   CONTINUE
        DO 220 KMUL=0,5
          KFLS=3
          IF(KMUL.EQ.0.OR.KMUL.EQ.3) KFLS=1
          IF(KMUL.EQ.5) KFLS=5
          KFLR=0
          IF(KMUL.EQ.2.OR.KMUL.EQ.3) KFLR=1
          IF(KMUL.EQ.4) KFLR=2
          DO 210 KFLB=1,5
            DO 200 KFLC=1,KFLB-1
              KF=10000*KFLR+100*KFLB+10*KFLC+KFLS
              CALL PYNAME(KF,CHAP)
              CALL PYNAME(-KF,CHAN)
              WRITE(MSTU(11),7600) KF,CHAP,-KF,CHAN
              IF(KF.EQ.311) THEN
                KFK=130
                CALL PYNAME(KFK,CHAP)
                WRITE(MSTU(11),7600) KFK,CHAP
                KFK=310
                CALL PYNAME(KFK,CHAP)
                WRITE(MSTU(11),7600) KFK,CHAP
              ENDIF
  200       CONTINUE
            KF=10000*KFLR+110*KFLB+KFLS
            CALL PYNAME(KF,CHAP)
            WRITE(MSTU(11),7600) KF,CHAP
  210     CONTINUE
  220   CONTINUE
        KF=100443
        CALL PYNAME(KF,CHAP)
        WRITE(MSTU(11),7600) KF,CHAP
        KF=100553
        CALL PYNAME(KF,CHAP)
        WRITE(MSTU(11),7600) KF,CHAP
        DO 260 KFLSP=1,3
          KFLS=2+2*(KFLSP/3)
          DO 250 KFLA=1,5
            DO 240 KFLB=1,KFLA
              DO 230 KFLC=1,KFLB
                IF(KFLSP.EQ.1.AND.(KFLA.EQ.KFLB.OR.KFLB.EQ.KFLC))
     &          GOTO 230
                IF(KFLSP.EQ.2.AND.KFLA.EQ.KFLC) GOTO 230
                IF(KFLSP.EQ.1) KF=1000*KFLA+100*KFLC+10*KFLB+KFLS
                IF(KFLSP.GE.2) KF=1000*KFLA+100*KFLB+10*KFLC+KFLS
                CALL PYNAME(KF,CHAP)
                CALL PYNAME(-KF,CHAN)
                WRITE(MSTU(11),7600) KF,CHAP,-KF,CHAN
  230         CONTINUE
  240       CONTINUE
  250     CONTINUE
  260   CONTINUE
        DO 270 KC=1,500
          KF=KCHG(KC,4)
          IF(KF.LT.1000000) GOTO 270
          CALL PYNAME(KF,CHAP)
          CALL PYNAME(-KF,CHAN)
          IF(CHAP.NE.' '.AND.CHAN.EQ.' ') WRITE(MSTU(11),7600) KF,CHAP
          IF(CHAN.NE.' ') WRITE(MSTU(11),7600) KF,CHAP,-KF,CHAN
  270   CONTINUE
 
C...List parton/particle data table. Check whether to be listed.
      ELSEIF(MLIST.EQ.12) THEN
        WRITE(MSTU(11),7700)
        DO 300 KC=1,MSTU(6)
          KF=KCHG(KC,4)
          IF(KF.EQ.0) GOTO 300
          IF(KF.LT.MSTU(1).OR.(MSTU(2).GT.0.AND.KF.GT.MSTU(2)))
     &    GOTO 300
 
C...Find particle name and mass. Print information.
          CALL PYNAME(KF,CHAP)
          IF(KF.LE.100.AND.CHAP.EQ.' '.AND.MDCY(KC,2).EQ.0) GOTO 300
          CALL PYNAME(-KF,CHAN)
          WRITE(MSTU(11),7800) KF,KC,CHAP,CHAN,(KCHG(KC,J1),J1=1,3),
     &    (PMAS(KC,J2),J2=1,4),MDCY(KC,1)
 
C...Particle decay: channel number, branching ratios, matrix element,
C...decay products.
          DO 290 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1
            DO 280 J=1,5
              CALL PYNAME(KFDP(IDC,J),CHAD(J))
  280       CONTINUE
            WRITE(MSTU(11),7900) IDC,MDME(IDC,1),MDME(IDC,2),BRAT(IDC),
     &      (CHAD(J),J=1,5)
  290     CONTINUE
  300   CONTINUE
 
C...List parameter value table.
      ELSEIF(MLIST.EQ.13) THEN
        WRITE(MSTU(11),8000)
        DO 310 I=1,200
          WRITE(MSTU(11),8100) I,MSTU(I),PARU(I),MSTJ(I),PARJ(I),PARF(I)
  310   CONTINUE
      ENDIF
 
C...Format statements for output on unit MSTU(11) (by default 6).
 5100 FORMAT(///28X,'Event listing (summary)'//4X,'I particle/jet KS',
     &5X,'KF  orig    p_x      p_y      p_z       E        m'/)
 5200 FORMAT(///28X,'Event listing (standard)'//4X,'I  particle/jet',
     &'  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)',
     &'       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/)
 5300 FORMAT(///28X,'Event listing (with vertices)'//4X,'I  particle/j',
     &'et  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)       P(I,1)',
     &'       P(I,2)       P(I,3)       P(I,4)       P(I,5)'/73X,
     &'V(I,1)       V(I,2)       V(I,3)       V(I,4)       V(I,5)'/)
 5400 FORMAT(///28X,'Event listing (no momenta)'//4X,'I  particle/jet',
     &     '  K(I,1)   K(I,2) K(I,3)     K(I,4)      K(I,5)',1X
     &     ,'   C tag  AC tag'/)
 5500 FORMAT(1X,I4,1X,A12,1X,I2,I8,1X,I4,5F9.3)
 5600 FORMAT(1X,I4,1X,A12,1X,I2,I8,1X,I4,5F9.2)
 5700 FORMAT(1X,I4,1X,A12,1X,I2,I8,1X,I4,5F9.1)
 5800 FORMAT(1X,I4,2X,A16,1X,I3,1X,I9,1X,I4,2(3X,I1,2I4),5F13.5)
 5900 FORMAT(1X,I4,2X,A16,1X,I3,1X,I9,1X,I4,2(3X,I1,2I4),1X,2I8)
 6000 FORMAT(1X,I4,2X,A16,1X,I3,1X,I9,1X,I4,2(3X,I9),5F13.5)
 6100 FORMAT(1X,I4,2X,A16,1X,I3,1X,I9,1X,I4,2(3X,I9),1X,2I8)
 6200 FORMAT(66X,5(1X,F12.3))
 6300 FORMAT(1X,78('='))
 6400 FORMAT(1X,130('='))
 6500 FORMAT(1X,65('='))
 6600 FORMAT(19X,'sum:',F6.2,5X,5F9.3)
 6700 FORMAT(19X,'sum:',F6.2,5X,5F9.2)
 6800 FORMAT(19X,'sum:',F6.2,5X,5F9.1)
 6900 FORMAT(19X,'sum charge:',F6.2,3X,'sum momentum and inv. mass:',
     &5F13.5)
 7000 FORMAT(19X,'sum charge:',F6.2)
 7100 FORMAT(/10X,'Event listing of HEPEVT common block (simplified)'
     &//'    I IST    ID   Mothers Daughters    p_x      p_y      p_z',
     &'       E        m')
 7200 FORMAT(1X,I4,I2,I8,4I5,5F9.3)
 7300 FORMAT(/10X,'Event listing of user process at input (simplified)'
     &//'   I IST     ID Mothers   Colours    p_x      p_y      p_z',
     &'       E        m')
 7400 FORMAT(1X,I3,I3,I8,2I4,2I5,5F9.3)
 7500 FORMAT(///20X,'List of KF codes in program'/)
 7600 FORMAT(4X,I9,4X,A16,6X,I9,4X,A16)
 7700 FORMAT(///30X,'Particle/parton data table'//8X,'KF',5X,'KC',4X,
     &'particle',8X,'antiparticle',6X,'chg  col  anti',8X,'mass',7X,
     &'width',7X,'w-cut',5X,'lifetime',1X,'decay'/11X,'IDC',1X,'on/off',
     &1X,'ME',3X,'Br.rat.',4X,'decay products')
 7800 FORMAT(/1X,I9,3X,I4,4X,A16,A16,3I5,1X,F12.5,2(1X,F11.5),
     &1X,1P,E13.5,3X,I2)
 7900 FORMAT(10X,I4,2X,I3,2X,I3,2X,F10.6,4X,5A16)
 8000 FORMAT(///20X,'Parameter value table'//4X,'I',3X,'MSTU(I)',
     &8X,'PARU(I)',3X,'MSTJ(I)',8X,'PARJ(I)',8X,'PARF(I)')
 8100 FORMAT(1X,I4,1X,I9,1X,F14.5,1X,I9,1X,F14.5,1X,F14.5)
 
      RETURN
      END
