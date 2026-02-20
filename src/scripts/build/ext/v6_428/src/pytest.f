 
C*********************************************************************
 
C...PYTEST
C...A simple program (disguised as subroutine) to run at installation
C...as a check that the program works as intended.
 
      SUBROUTINE PYTEST(MTEST)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/
C...Local arrays.
      DIMENSION PSUM(5),PINI(6),PFIN(6)
 
C...Save defaults for values that are changed.
      MSTJ1=MSTJ(1)
      MSTJ3=MSTJ(3)
      MSTJ11=MSTJ(11)
      MSTJ42=MSTJ(42)
      MSTJ43=MSTJ(43)
      MSTJ44=MSTJ(44)
      PARJ17=PARJ(17)
      PARJ22=PARJ(22)
      PARJ43=PARJ(43)
      PARJ54=PARJ(54)
      MST101=MSTJ(101)
      MST104=MSTJ(104)
      MST105=MSTJ(105)
      MST107=MSTJ(107)
      MST116=MSTJ(116)
 
C...First part: loop over simple events to be generated.
      IF(MTEST.GE.1) CALL PYTABU(20)
      NERR=0
      DO 180 IEV=1,500
 
C...Reset parameter values. Switch on some nonstandard features.
        MSTJ(1)=1
        MSTJ(3)=0
        MSTJ(11)=1
        MSTJ(42)=2
        MSTJ(43)=4
        MSTJ(44)=2
        PARJ(17)=0.1D0
        PARJ(22)=1.5D0
        PARJ(43)=1D0
        PARJ(54)=-0.05D0
        MSTJ(101)=5
        MSTJ(104)=5
        MSTJ(105)=0
        MSTJ(107)=1
        IF(IEV.EQ.301.OR.IEV.EQ.351.OR.IEV.EQ.401) MSTJ(116)=3
 
C...Ten events each for some single jets configurations.
        IF(IEV.LE.50) THEN
          ITY=(IEV+9)/10
          MSTJ(3)=-1
          IF(ITY.EQ.3.OR.ITY.EQ.4) MSTJ(11)=2
          IF(ITY.EQ.1) CALL PY1ENT(1,1,15D0,0D0,0D0)
          IF(ITY.EQ.2) CALL PY1ENT(1,3101,15D0,0D0,0D0)
          IF(ITY.EQ.3) CALL PY1ENT(1,-2203,15D0,0D0,0D0)
          IF(ITY.EQ.4) CALL PY1ENT(1,-4,30D0,0D0,0D0)
          IF(ITY.EQ.5) CALL PY1ENT(1,21,15D0,0D0,0D0)
 
C...Ten events each for some simple jet systems; string fragmentation.
        ELSEIF(IEV.LE.130) THEN
          ITY=(IEV-41)/10
          IF(ITY.EQ.1) CALL PY2ENT(1,1,-1,40D0)
          IF(ITY.EQ.2) CALL PY2ENT(1,4,-4,30D0)
          IF(ITY.EQ.3) CALL PY2ENT(1,2,2103,100D0)
          IF(ITY.EQ.4) CALL PY2ENT(1,21,21,40D0)
          IF(ITY.EQ.5) CALL PY3ENT(1,2101,21,-3203,30D0,0.6D0,0.8D0)
          IF(ITY.EQ.6) CALL PY3ENT(1,5,21,-5,40D0,0.9D0,0.8D0)
          IF(ITY.EQ.7) CALL PY3ENT(1,21,21,21,60D0,0.7D0,0.5D0)
          IF(ITY.EQ.8) CALL PY4ENT(1,2,21,21,-2,40D0,
     &    0.4D0,0.64D0,0.6D0,0.12D0,0.2D0)
 
C...Seventy events with independent fragmentation and momentum cons.
        ELSEIF(IEV.LE.200) THEN
          ITY=1+(IEV-131)/16
          MSTJ(2)=1+MOD(IEV-131,4)
          MSTJ(3)=1+MOD((IEV-131)/4,4)
          IF(ITY.EQ.1) CALL PY2ENT(1,4,-5,40D0)
          IF(ITY.EQ.2) CALL PY3ENT(1,3,21,-3,40D0,0.9D0,0.4D0)
          IF(ITY.EQ.3) CALL PY4ENT(1,2,21,21,-2,40D0,
     &    0.4D0,0.64D0,0.6D0,0.12D0,0.2D0)
          IF(ITY.GE.4) CALL PY4ENT(1,2,-3,3,-2,40D0,
     &    0.4D0,0.64D0,0.6D0,0.12D0,0.2D0)
 
C...A hundred events with random jets (check invariant mass).
        ELSEIF(IEV.LE.300) THEN
  100     DO 110 J=1,5
            PSUM(J)=0D0
  110     CONTINUE
          NJET=2D0+6D0*PYR(0)
          DO 130 I=1,NJET
            KFL=21
            IF(I.EQ.1) KFL=INT(1D0+4D0*PYR(0))
            IF(I.EQ.NJET) KFL=-INT(1D0+4D0*PYR(0))
            EJET=5D0+20D0*PYR(0)
            THETA=ACOS(2D0*PYR(0)-1D0)
            PHI=6.2832D0*PYR(0)
            IF(I.LT.NJET) CALL PY1ENT(-I,KFL,EJET,THETA,PHI)
            IF(I.EQ.NJET) CALL PY1ENT(I,KFL,EJET,THETA,PHI)
            IF(I.EQ.1.OR.I.EQ.NJET) MSTJ(93)=1
            IF(I.EQ.1.OR.I.EQ.NJET) PSUM(5)=PSUM(5)+PYMASS(KFL)
            DO 120 J=1,4
              PSUM(J)=PSUM(J)+P(I,J)
  120       CONTINUE
  130     CONTINUE
          IF(PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2.LT.
     &    (PSUM(5)+PARJ(32))**2) GOTO 100
 
C...Fifty e+e- continuum events with matrix elements.
        ELSEIF(IEV.LE.350) THEN
          MSTJ(101)=2
          CALL PYEEVT(0,40D0)
 
C...Fifty e+e- continuum event with varying shower options.
        ELSEIF(IEV.LE.400) THEN
          MSTJ(42)=1+MOD(IEV,2)
          MSTJ(43)=1+MOD(IEV/2,4)
          MSTJ(44)=MOD(IEV/8,3)
          CALL PYEEVT(0,90D0)
 
C...Fifty e+e- continuum events with coherent shower.
        ELSEIF(IEV.LE.450) THEN
          CALL PYEEVT(0,500D0)
 
C...Fifty Upsilon decays to ggg or gammagg with coherent shower.
        ELSE
          CALL PYONIA(5,9.46D0)
        ENDIF
 
C...Generate event. Find total momentum, energy and charge.
        DO 140 J=1,4
          PINI(J)=PYP(0,J)
  140   CONTINUE
        PINI(6)=PYP(0,6)
        CALL PYEXEC
        DO 150 J=1,4
          PFIN(J)=PYP(0,J)
  150   CONTINUE
        PFIN(6)=PYP(0,6)
 
C...Check conservation of energy, momentum and charge;
C...usually exact, but only approximate for single jets.
        MERR=0
        IF(IEV.LE.50) THEN
          IF((PFIN(1)-PINI(1))**2+(PFIN(2)-PINI(2))**2.GE.10D0)
     &    MERR=MERR+1
          EPZREM=PINI(4)+PINI(3)-PFIN(4)-PFIN(3)
          IF(EPZREM.LT.0D0.OR.EPZREM.GT.2D0*PARJ(31)) MERR=MERR+1
          IF(ABS(PFIN(6)-PINI(6)).GT.2.1D0) MERR=MERR+1
        ELSE
          DO 160 J=1,4
            IF(ABS(PFIN(J)-PINI(J)).GT.0.0001D0*PINI(4)) MERR=MERR+1
  160     CONTINUE
          IF(ABS(PFIN(6)-PINI(6)).GT.0.1D0) MERR=MERR+1
        ENDIF
        IF(MERR.NE.0) WRITE(MSTU(11),5000) (PINI(J),J=1,4),PINI(6),
     &  (PFIN(J),J=1,4),PFIN(6)
 
C...Check that all KF codes are known ones, and that partons/particles
C...satisfy energy-momentum-mass relation. Store particle statistics.
        DO 170 I=1,N
          IF(K(I,1).GT.20) GOTO 170
          IF(PYCOMP(K(I,2)).EQ.0) THEN
            WRITE(MSTU(11),5100) I
            MERR=MERR+1
          ENDIF
          PD=P(I,4)**2-P(I,1)**2-P(I,2)**2-P(I,3)**2-P(I,5)**2
          IF(ABS(PD).GT.MAX(0.1D0,0.001D0*P(I,4)**2).OR.P(I,4).LT.0D0)
     &    THEN
            WRITE(MSTU(11),5200) I
            MERR=MERR+1
          ENDIF
  170   CONTINUE
        IF(MTEST.GE.1) CALL PYTABU(21)
 
C...List all erroneous events and some normal ones.
        IF(MERR.NE.0.OR.MSTU(24).NE.0.OR.MSTU(28).NE.0) THEN
          IF(MERR.GE.1) WRITE(MSTU(11),6400)
          CALL PYLIST(2)
        ELSEIF(MTEST.GE.1.AND.MOD(IEV-5,100).EQ.0) THEN
          CALL PYLIST(1)
        ENDIF
 
C...Stop execution if too many errors.
        IF(MERR.NE.0) NERR=NERR+1
        IF(NERR.GE.10) THEN
          WRITE(MSTU(11),6300)
          CALL PYLIST(1)
          CALL PYSTOP(9)
        ENDIF
  180 CONTINUE
 
C...Summarize result of run.
      IF(MTEST.GE.1) CALL PYTABU(22)
 
C...Reset commonblock variables changed during run.
      MSTJ(1)=MSTJ1
      MSTJ(3)=MSTJ3
      MSTJ(11)=MSTJ11
      MSTJ(42)=MSTJ42
      MSTJ(43)=MSTJ43
      MSTJ(44)=MSTJ44
      PARJ(17)=PARJ17
      PARJ(22)=PARJ22
      PARJ(43)=PARJ43
      PARJ(54)=PARJ54
      MSTJ(101)=MST101
      MSTJ(104)=MST104
      MSTJ(105)=MST105
      MSTJ(107)=MST107
      MSTJ(116)=MST116
 
C...Second part: complete events of various kinds.
C...Common initial values. Loop over initiating conditions.
      MSTP(122)=MAX(0,MIN(2,MTEST))
      MDCY(PYCOMP(111),1)=0
      DO 230 IPROC=1,8
 
C...Reset process type, kinematics cuts, and the flags used.
        MSEL=0
        DO 190 ISUB=1,500
          MSUB(ISUB)=0
  190   CONTINUE
        CKIN(1)=2D0
        CKIN(3)=0D0
        MSTP(2)=1
        MSTP(11)=0
        MSTP(33)=0
        MSTP(81)=1
        MSTP(82)=1
        MSTP(111)=1
        MSTP(131)=0
        MSTP(133)=0
        PARP(131)=0.01D0
 
C...Prompt photon production at fixed target.
        IF(IPROC.EQ.1) THEN
          PZSUM=300D0
          PESUM=SQRT(PZSUM**2+PYMASS(211)**2)+PYMASS(2212)
          PQSUM=2D0
          MSEL=10
          CKIN(3)=5D0
          CALL PYINIT('FIXT','pi+','p',PZSUM)
 
C...QCD processes at ISR energies.
        ELSEIF(IPROC.EQ.2) THEN
          PESUM=63D0
          PZSUM=0D0
          PQSUM=2D0
          MSEL=1
          CKIN(3)=5D0
          CALL PYINIT('CMS','p','p',PESUM)
 
C...W production + multiple interactions at CERN Collider.
        ELSEIF(IPROC.EQ.3) THEN
          PESUM=630D0
          PZSUM=0D0
          PQSUM=0D0
          MSEL=12
          CKIN(1)=20D0
          MSTP(82)=4
          MSTP(2)=2
          MSTP(33)=3
          CALL PYINIT('CMS','p','pbar',PESUM)
 
C...W/Z gauge boson pairs + pileup events at the Tevatron.
        ELSEIF(IPROC.EQ.4) THEN
          PESUM=1800D0
          PZSUM=0D0
          PQSUM=0D0
          MSUB(22)=1
          MSUB(23)=1
          MSUB(25)=1
          CKIN(1)=200D0
          MSTP(111)=0
          MSTP(131)=1
          MSTP(133)=2
          PARP(131)=0.04D0
          CALL PYINIT('CMS','p','pbar',PESUM)
 
C...Higgs production at LHC.
        ELSEIF(IPROC.EQ.5) THEN
          PESUM=15400D0
          PZSUM=0D0
          PQSUM=2D0
          MSUB(3)=1
          MSUB(102)=1
          MSUB(123)=1
          MSUB(124)=1
          PMAS(25,1)=300D0
          CKIN(1)=200D0
          MSTP(81)=0
          MSTP(111)=0
          CALL PYINIT('CMS','p','p',PESUM)
 
C...Z' production at SSC.
        ELSEIF(IPROC.EQ.6) THEN
          PESUM=40000D0
          PZSUM=0D0
          PQSUM=2D0
          MSEL=21
          PMAS(32,1)=600D0
          CKIN(1)=400D0
          MSTP(81)=0
          MSTP(111)=0
          CALL PYINIT('CMS','p','p',PESUM)
 
C...W pair production at 1 TeV e+e- collider.
        ELSEIF(IPROC.EQ.7) THEN
          PESUM=1000D0
          PZSUM=0D0
          PQSUM=0D0
          MSUB(25)=1
          MSUB(69)=1
          MSTP(11)=1
          CALL PYINIT('CMS','e+','e-',PESUM)
 
C...Deep inelastic scattering at a LEP+LHC ep collider.
        ELSEIF(IPROC.EQ.8) THEN
          P(1,1)=0D0
          P(1,2)=0D0
          P(1,3)=8000D0
          P(2,1)=0D0
          P(2,2)=0D0
          P(2,3)=-80D0
          PESUM=8080D0
          PZSUM=7920D0
          PQSUM=0D0
          MSUB(10)=1
          CKIN(3)=50D0
          MSTP(111)=0
          CALL PYINIT('3MOM','p','e-',PESUM)
        ENDIF
 
C...Generate 20 events of each required type.
        DO 220 IEV=1,20
          CALL PYEVNT
          PESUMM=PESUM
          IF(IPROC.EQ.4) PESUMM=MSTI(41)*PESUM
 
C...Check conservation of energy/momentum/flavour.
          PINI(1)=0D0
          PINI(2)=0D0
          PINI(3)=PZSUM
          PINI(4)=PESUMM
          PINI(6)=PQSUM
          DO 200 J=1,4
            PFIN(J)=PYP(0,J)
  200     CONTINUE
          PFIN(6)=PYP(0,6)
          MERR=0
          DEVE=ABS(PFIN(4)-PINI(4))+ABS(PFIN(3)-PINI(3))
          DEVT=ABS(PFIN(1)-PINI(1))+ABS(PFIN(2)-PINI(2))
          DEVQ=ABS(PFIN(6)-PINI(6))
          IF(DEVE.GT.2D-3*PESUM.OR.DEVT.GT.MAX(0.01D0,1D-4*PESUM).OR.
     &    DEVQ.GT.0.1D0) MERR=1
          IF(MERR.NE.0) WRITE(MSTU(11),5000) (PINI(J),J=1,4),PINI(6),
     &    (PFIN(J),J=1,4),PFIN(6)
 
C...Check that all KF codes are known ones, and that partons/particles
C...satisfy energy-momentum-mass relation.
          DO 210 I=1,N
            IF(K(I,1).GT.20) GOTO 210
            IF(PYCOMP(K(I,2)).EQ.0) THEN
              WRITE(MSTU(11),5100) I
              MERR=MERR+1
            ENDIF
            PD=P(I,4)**2-P(I,1)**2-P(I,2)**2-P(I,3)**2-P(I,5)**2*
     &      SIGN(1D0,P(I,5))
            IF(ABS(PD).GT.MAX(0.1D0,0.002D0*P(I,4)**2,0.002D0*P(I,5)**2)
     &      .OR.(P(I,5).GE.0D0.AND.P(I,4).LT.0D0)) THEN
              WRITE(MSTU(11),5200) I
              MERR=MERR+1
            ENDIF
  210     CONTINUE
 
C...Listing of erroneous events, and first event of each type.
          IF(MERR.GE.1) NERR=NERR+1
          IF(NERR.GE.10) THEN
            WRITE(MSTU(11),6300)
            CALL PYLIST(1)
            CALL PYSTOP(9)
          ENDIF
          IF(MTEST.GE.1.AND.(MERR.GE.1.OR.IEV.EQ.1)) THEN
            IF(MERR.GE.1) WRITE(MSTU(11),6400)
            CALL PYLIST(1)
          ENDIF
  220   CONTINUE
 
C...List statistics for each process type.
        IF(MTEST.GE.1) CALL PYSTAT(1)
  230 CONTINUE
 
C...Summarize result of run.
      IF(NERR.EQ.0) WRITE(MSTU(11),6500)
      IF(NERR.GT.0) WRITE(MSTU(11),6600) NERR
 
C...Format statements for output.
 5000 FORMAT(/' Momentum, energy and/or charge were not conserved ',
     &'in following event'/' sum of',9X,'px',11X,'py',11X,'pz',11X,
     &'E',8X,'charge'/' before',2X,4(1X,F12.5),1X,F8.2/' after',3X,
     &4(1X,F12.5),1X,F8.2)
 5100 FORMAT(/5X,'Entry no.',I4,' in following event not known code')
 5200 FORMAT(/5X,'Entry no.',I4,' in following event has faulty ',
     &'kinematics')
 6300 FORMAT(/5X,'This is the tenth error experienced! Something is ',
     &'wrong.'/5X,'Execution will be stopped after listing of event.')
 6400 FORMAT(5X,'Faulty event follows:')
 6500 FORMAT(//5X,'End result of PYTEST: no errors detected.')
 6600 FORMAT(//5X,'End result of PYTEST:',I2,' errors detected.'/
     &5X,'This should not have happened!')
 
      RETURN
      END
