 
C*********************************************************************
 
C...PYMIRM
C...Picks primordial kT and shares longitudinal momentum among
C...beam remnants.
 
      SUBROUTINE PYMIRM
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...The event record
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Parameters
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
C...The common block of colour tags.
      COMMON/PYCTAG/NCT,MCT(4000,2)
C...The common block of dangling ends
      COMMON/PYINTM/KFIVAL(2,3),NMI(2),IMI(2,800,2),NVC(2,-6:6),
     &     XASSOC(2,-6:6,240),XPSVC(-6:6,-1:240),PVCTOT(2,-1:1),
     &     XMI(2,240),PT2MI(240),IMISEP(0:240)
      SAVE /PYJETS/,/PYDAT1/,/PYPARS/,/PYINT1/,/PYINTM/,/PYCTAG/
C...Local variables
      DIMENSION W(0:2,0:2),VB(3),NNXT(2),IVALQ(2),ICOMQ(2)
C...W(I,J)|  J=0    |   1   |   2   |
C...  I=0 | Wrem**2 |  W+   |  W-   |
C...    1 | W1**2   |  W1+  |  W1-  |
C...    2 | W2**2   |  W2+  |  W2-  |
C...4-product
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)
C...Tentative parametrization of <kT> as a function of Q.
      SIGPT(Q)=MAX(PARJ(21),2.1D0*Q/(7D0+Q))
C      SIGPT(Q)=MAX(0.36D0,4D0*SQRT(Q)/(10D0+SQRT(Q))
C      SIGPT(Q)=MAX(PARJ(21),3D0*SQRT(Q)/(5D0+SQRT(Q))
      GETPT(Q,SIGMA)=MIN(SIGMA*SQRT(-LOG(PYR(0))),PARP(93))
C...Lambda kinematic function.
      FLAM(A,B,C)=A**2+B**2+C**2-2D0*(A*B+B*C+C*A)
 
C...Beginning and end of beam remnant partons
      NOUT=MINT(53)
      ISUB=MINT(1)
 
C...Loopback point if kinematic choices gives impossible configuration.
      NTRY=0
  100 NTRY=NTRY+1
 
C...Assign kT values on each side separately.
      DO 180 JS=1,2
 
C...First zero all kT on this side. Skip if no kT to generate.
        DO 110 IM=1,NMI(JS)
          P(IMI(JS,IM,1),1)=0D0
          P(IMI(JS,IM,1),2)=0D0
  110   CONTINUE
        IF(MSTP(91).LE.0) GOTO 180
 
C...Now assign kT to each (non-collapsed) parton in IMI.
        DO 170 IM=1,NMI(JS)
          I=IMI(JS,IM,1)
C...Select kT according to truncated gaussian or 1/kt6 tails.
C...For first interaction, either use rms width = PARP(91) or fitted.
          IF (IM.EQ.1) THEN
            SIGMA=PARP(91)
            IF (MSTP(91).GE.11.AND.MSTP(91).LE.20) THEN
              Q=SQRT(PT2MI(IM))
              SIGMA=SIGPT(Q)
            ENDIF
          ELSE
C...For subsequent interactions and BR partons use fragmentation width.
            SIGMA=PARJ(21)
          ENDIF
          PHI=PARU(2)*PYR(0)
          PT=0D0
          IF(NTRY.LE.100) THEN
 111        IF (MSTP(91).EQ.1.OR.MSTP(91).EQ.11) THEN
              PT=GETPT(Q,SIGMA)
              PTX=PT*COS(PHI)
              PTY=PT*SIN(PHI)
            ELSEIF (MSTP(91).EQ.2) THEN
              CALL PYERRM(1,'(PYMIRM:) Sorry, MSTP(91)=2 not '//
     &          'available, using MSTP(91)=1.')
              CALL PYGIVE('MSTP(91)=1')
              GOTO 111
            ELSEIF(MSTP(91).EQ.3.OR.MSTP(91).EQ.13) THEN
C...Use distribution with kt**6 tails, rms width = PARP(91).
              EPS=SQRT(3D0/2D0)*SIGMA
C...Generate PTX and PTY separately, each propto 1/KT**6
              DO 119 IXY=1,2
C...Decide which interval to try
 112            P12=1D0/(1D0+27D0/40D0*SIGMA**6/EPS**6)
                IF (PYR(0).LT.P12) THEN
C...Use flat approx with accept/reject up to EPS.
                  PT=PYR(0)*EPS
                  WT=(3D0/2D0*SIGMA**2/(PT**2+3D0/2D0*SIGMA**2))**3
                  IF (PYR(0).GT.WT) GOTO 112
                ELSE
C...Above EPS, use 1/kt**6 approx with accept/reject.
                  PT=EPS/(PYR(0)**(1D0/5D0))
                  WT=PT**6/(PT**2+3D0/2D0*SIGMA**2)**3
                  IF (PYR(0).GT.WT) GOTO 112
                ENDIF
                MSIGN=1
                IF (PYR(0).GT.0.5D0) MSIGN=-1
                IF (IXY.EQ.1) PTX=MSIGN*PT
                IF (IXY.EQ.2) PTY=MSIGN*PT
 119          CONTINUE
            ELSEIF (MSTP(91).EQ.4.OR.MSTP(91).EQ.14) THEN
              PTX=SIGMA*(SQRT(6D0)*PYR(0)-SQRT(3D0/2D0))
              PTY=SIGMA*(SQRT(6D0)*PYR(0)-SQRT(3D0/2D0))
            ENDIF
C...Adjust final PT. Impose upper cutoff, or zero for soft evts.
            PT=SQRT(PTX**2+PTY**2)
            WT=1D0
            IF (PT.GT.PARP(93)) WT=SQRT(PARP(93)/PT)
            IF(ISUB.EQ.95.AND.IM.EQ.1) WT=0D0
            PTX=PTX*WT
            PTY=PTY*WT
            PT=SQRT(PTX**2+PTY**2)
          ENDIF
 
          P(I,1)=P(I,1)+PTX
          P(I,2)=P(I,2)+PTY
 
C...Compensation kicks, with varying degree of local anticorrelations.
          MCORR=MSTP(90)
          IF (MCORR.EQ.0.OR.ISUB.EQ.95) THEN
            PTCX=-PTX/(NMI(JS)-1)
            PTCY=-PTY/(NMI(JS)-1)
            IF(ISUB.EQ.95) THEN
              PTCX=-PTX/(NMI(JS)-2)
              PTCY=-PTY/(NMI(JS)-2)
            ENDIF
            DO 120 IMC=1,NMI(JS)
              IF (IMC.EQ.IM) GOTO 120
              IF(ISUB.EQ.95.AND.IMC.EQ.1) GOTO 120
              P(IMI(JS,IMC,1),1)=P(IMI(JS,IMC,1),1)+PTCX
              P(IMI(JS,IMC,1),2)=P(IMI(JS,IMC,1),2)+PTCY
  120       CONTINUE
          ELSEIF (MCORR.GE.1) THEN
            DO 140 MSID=4,5
              NNXT(MSID-3)=0
C...Count up # of neighbours on either side
              IMO=I
  130         IMO=K(IMO,MSID)/MSTU(5)
              IF (IMO.EQ.0) GOTO 140
              NNXT(MSID-3)=NNXT(MSID-3)+1
C...Stop at quarks and junctions
              IF (MCORR.EQ.1.AND.K(IMO,2).EQ.21) GOTO 130
  140       CONTINUE
C...How should compensation be shared when unequal numbers on the
C...two sides? 50/50 regardless? N1:N2? Assume latter for now.
            NSUM=NNXT(1)+NNXT(2)
            T1=0
            DO 160 MSID=4,5
C...Total momentum to be compensated on this side
              IF (NNXT(MSID-3).EQ.0) GOTO 160
              PTCX=-(NNXT(MSID-3)*PTX)/NSUM
              PTCY=-(NNXT(MSID-3)*PTY)/NSUM
C...RS: compensation supression factor as we go out from parton I.
C...Hardcoded behaviour RS=0.5, i.e. 1/2**n falloff,
C...since (for now) MSTP(90) provides enough variability.
              RS=0.5D0
              FAC=(1D0-RS)/(RS*(1-RS**NNXT(MSID-3)))
              IMO=I
  150         IDA=IMO
              IMO=K(IMO,MSID)/MSTU(5)
              IF (IMO.EQ.0) GOTO 160
              FAC=FAC*RS
              IF (K(IMO,2).NE.88) THEN
                P(IMO,1)=P(IMO,1)+FAC*PTCX
                P(IMO,2)=P(IMO,2)+FAC*PTCY
                IF (MCORR.EQ.1.AND.K(IMO,2).EQ.21) GOTO 150
C...If we reach junction, divide out the kT that would have been
C...assigned to the junction on each of its other legs.
              ELSE
                L1=MOD(K(IMO,4),MSTU(5))
                L2=K(IMO,5)/MSTU(5)
                L3=MOD(K(IMO,5),MSTU(5))
                P(L1,1)=P(L1,1)+0.5D0*FAC*PTCX
                P(L1,2)=P(L1,2)+0.5D0*FAC*PTCY
                P(L2,1)=P(L2,1)+0.5D0*FAC*PTCX
                P(L2,2)=P(L2,2)+0.5D0*FAC*PTCY
                P(L3,1)=P(L3,1)+0.5D0*FAC*PTCX
                P(L3,2)=P(L3,2)+0.5D0*FAC*PTCY
                P(IDA,1)=P(IDA,1)-0.5D0*FAC*PTCX
                P(IDA,2)=P(IDA,2)-0.5D0*FAC*PTCY
              ENDIF
 
  160       CONTINUE
          ENDIF
  170   CONTINUE
C...End assignment of kT values to initiators and remnants.
  180 CONTINUE
 
C...Check kinematics constraints for non-BR partons.
      DO 190 IM=1,MINT(31)
        SHAT=XMI(1,IM)*XMI(2,IM)*VINT(2)
        PT1=SQRT(P(IMI(1,IM,1),1)**2+P(IMI(1,IM,1),2)**2)
        PT2=SQRT(P(IMI(2,IM,1),1)**2+P(IMI(2,IM,1),2)**2)
        PT1PT2=P(IMI(1,IM,1),1)*P(IMI(2,IM,1),1)
     &        +P(IMI(1,IM,1),2)*P(IMI(2,IM,1),2)
        IF (SHAT.LT.2D0*(PT1*PT2-PT1PT2).AND.NTRY.LE.100) THEN
          IF(NTRY.GE.100) THEN
C...Kill this event and start another.
            CALL PYERRM(1,
     &           '(PYMIRM:) No consistent (x,kT) sets found')
            MINT(51)=1
            RETURN
          ENDIF
          GOTO 100
        ENDIF
  190 CONTINUE
 
C...Calculate W+ and W- available for combined remnant system.
      W(0,1)=VINT(1)
      W(0,2)=VINT(1)
      DO 200 IM=1,MINT(31)
        PT2 = (P(IMI(1,IM,1),1)+P(IMI(2,IM,1),1))**2
     &       +(P(IMI(1,IM,1),2)+P(IMI(2,IM,1),2))**2
        ST=XMI(1,IM)*XMI(2,IM)*VINT(2)+PT2
        W(0,1)=W(0,1)-SQRT(XMI(1,IM)/XMI(2,IM)*ST)
        W(0,2)=W(0,2)-SQRT(XMI(2,IM)/XMI(1,IM)*ST)
  200 CONTINUE
C...Also store Wrem**2 = W+ * W-
      W(0,0)=W(0,1)*W(0,2)
 
      IF ((W(0,0).LT.0D0.OR.W(0,1)+W(0,2).LT.0D0).AND.NTRY.LE.100) THEN
          IF(NTRY.GE.100) THEN
C...Kill this event and start another.
            CALL PYERRM(1,
     &    '(PYMIRM:) Negative beam remnant mass squared unavoidable')
            MINT(51)=1
            RETURN
          ENDIF
          GOTO 100
      ENDIF

C...Assign unscaled x values to partons/hadrons in each of the
C...beam remnants and calculate unscaled W+ and W- from them.
      NTRYX=0
  210 NTRYX=NTRYX+1
      DO 280 JS=1,2
        W(JS,1)=0D0
        W(JS,2)=0D0
        DO 270 IM=MINT(31)+1,NMI(JS)
          I=IMI(JS,IM,1)
          KF=K(I,2)
          KFA=IABS(KF)
          ICOMP=IMI(JS,IM,2)
 
C...Skip collapsed gluons and junctions. Reset.
          IF (KFA.EQ.21.AND.K(I,1).EQ.14) GOTO 270
          IF (KFA.EQ.88) GOTO 270
          X=0D0
          IVALQ(1)=0
          IVALQ(2)=0
          ICOMQ(1)=0
          ICOMQ(2)=0
 
C...If gluon then only beam remnant, so takes all.
          IF(KFA.EQ.21) THEN
            X=1D0
C...If valence quark then use parametrized valence distribution.
          ELSEIF(KFA.LE.6.AND.ICOMP.EQ.0) THEN
            IVALQ(1)=KF
C...If companion quark then derive from companion x.
          ELSEIF(KFA.LE.6) THEN
            ICOMQ(1)=ICOMP
C...If valence diquark then use two parametrized valence distributions.
          ELSEIF(KFA.GT.1000.AND.MOD(KFA/10,10).EQ.0.AND.
     &    ICOMP.EQ.0) THEN
            IVALQ(1)=ISIGN(KFA/1000,KF)
            IVALQ(2)=ISIGN(MOD(KFA/100,10),KF)
C...If valence+sea diquark then combine valence + companion choices.
          ELSEIF(KFA.GT.1000.AND.MOD(KFA/10,10).EQ.0.AND.
     &    ICOMP.LT.MSTU(5)) THEN
            IF(KFA/1000.EQ.IABS(K(ICOMP,2))) THEN
              IVALQ(1)=ISIGN(MOD(KFA/100,10),KF)
            ELSE
              IVALQ(1)=ISIGN(KFA/1000,KF)
            ENDIF
            ICOMQ(1)=ICOMP
C...Extra code: workaround for diquark made out of two sea
C...quarks, but where not (yet) ICOMP > MSTU(5).
            DO 220 IM1=1,MINT(31)
              IF(IMI(JS,IM1,2).EQ.I.AND.IMI(JS,IM1,1).NE.ICOMP) THEN
                ICOMQ(2)=IMI(JS,IM1,1)
                IVALQ(1)=0
              ENDIF
  220       CONTINUE
C...If sea diquark then sum of two derived from companion x.
          ELSEIF(KFA.GT.1000.AND.MOD(KFA/10,10).EQ.0) THEN
             ICOMQ(1)=MOD(ICOMP,MSTU(5))
             ICOMQ(2)=ICOMP/MSTU(5)
C...If meson or baryon then use fragmentation function.
C...Somewhat arbitrary split into old and new flavour, but OK normally.
          ELSE
            KFL3=MOD(KFA/10,10)
            IF(MOD(KFA/1000,10).EQ.0) THEN
              KFL1=MOD(KFA/100,10)
            ELSE
              KFL1=MOD(KFA,10000)-10*KFL3-1
              IF(MOD(KFA/1000,10).EQ.MOD(KFA/100,10).AND.
     &        MOD(KFA,10).EQ.2) KFL1=KFL1+2
            ENDIF
            PR=P(I,5)**2+P(I,1)**2+P(I,2)**2
            CALL PYZDIS(KFL1,KFL3,PR,X)
          ENDIF
 
          DO 260 IQ=1,2
C...Calculation of x of valence quark: assume form (1-x)^a/sqrt(x),
C...where a=3.5 for u in proton, =2 for d in proton and =0.8 for meson.
C...In other baryons combine u and d from proton appropriately.
            IF(IVALQ(IQ).NE.0) THEN
              NVAL=0
              IF(KFIVAL(JS,1).EQ.IVALQ(IQ)) NVAL=NVAL+1
              IF(KFIVAL(JS,2).EQ.IVALQ(IQ)) NVAL=NVAL+1
              IF(KFIVAL(JS,3).EQ.IVALQ(IQ)) NVAL=NVAL+1
C...Meson.
              IF(KFIVAL(JS,3).EQ.0) THEN
                MDU=0
C...Baryon with three identical quarks: mix u and d forms.
              ELSEIF(NVAL.EQ.3) THEN
                MDU=INT(PYR(0)+5D0/3D0)
C...Baryon, one of two identical quarks: u form.
              ELSEIF(NVAL.EQ.2) THEN
                MDU=2
C...Baryon with two identical quarks, but not the one picked: d form.
              ELSEIF(KFIVAL(JS,1).EQ.KFIVAL(JS,2).OR.KFIVAL(JS,2).EQ.
     &        KFIVAL(JS,3).OR.KFIVAL(JS,1).EQ.KFIVAL(JS,3)) THEN
                MDU=1
C...Baryon with three nonidentical quarks: mix u and d forms.
              ELSE
                MDU=INT(PYR(0)+5D0/3D0)
              ENDIF
              XPOW=0.8D0
              IF(MDU.EQ.1) XPOW=3.5D0
              IF(MDU.EQ.2) XPOW=2D0
  230         XX=PYR(0)**2
              IF((1D0-XX)**XPOW.LT.PYR(0)) GOTO 230
              X=X+XX
            ENDIF
 
C...Calculation of x of companion quark.
            IF(ICOMQ(IQ).NE.0) THEN
              XCOMP=1D-4
              DO 240 IM1=1,MINT(31)
                IF(IMI(JS,IM1,1).EQ.ICOMQ(IQ)) XCOMP=XMI(JS,IM1)
  240         CONTINUE
              NPOW=MAX(0,MIN(4,MSTP(87)))
  250         XX=XCOMP*(1D0/(1D0-PYR(0)*(1D0-XCOMP))-1D0)
              CORR=((1D0-XCOMP-XX)/(1D0-XCOMP))**NPOW*
     &        (XCOMP**2+XX**2)/(XCOMP+XX)**2
              IF(CORR.LT.PYR(0)) GOTO 250
              X=X+XX
            ENDIF
  260     CONTINUE
 
C...Optionally enchance x of composite systems (e.g. diquarks)
          IF (KFA.GT.100) X=PARP(79)*X
 
C...Store x. Also calculate light cone energies of each system.
          XMI(JS,IM)=X
          W(JS,JS)=W(JS,JS)+X
          W(JS,3-JS)=W(JS,3-JS)+(P(I,5)**2+P(I,1)**2+P(I,2)**2)/X
  270   CONTINUE
        W(JS,JS)=W(JS,JS)*W(0,JS)
        W(JS,3-JS)=W(JS,3-JS)/W(0,JS)
        W(JS,0)=W(JS,1)*W(JS,2)
  280 CONTINUE
 
C...Check W1 W2 < Wrem (can be done before rescaling, since W
C...insensitive to global rescalings of the BR x values).
      IF (SQRT(W(1,0))+SQRT(W(2,0)).GT.SQRT(W(0,0)).AND.NTRYX.LE.100)
     &     THEN
        GOTO 210
      ELSEIF (NTRYX.GT.100.AND.NTRY.LE.100) THEN
        GOTO 100
      ELSEIF (NTRYX.GT.100) THEN
        CALL PYERRM(1,'(PYMIRM:) No consistent (x,kT) sets found')
        MINT(57)=MINT(57)+1
        MINT(51)=1
        RETURN
      ENDIF
 
C...Compute x rescaling factors
      COMTRM=W(0,0)+SQRT(FLAM(W(0,0),W(1,0),W(2,0)))
      R1=(COMTRM+W(1,0)-W(2,0))/(2D0*W(1,1)*W(0,2))
      R2=(COMTRM+W(2,0)-W(1,0))/(2D0*W(2,2)*W(0,1))
 
      IF (R1.LT.0.OR.R2.LT.0) THEN
        CALL PYERRM(19,'(PYMIRM:) negative rescaling factors !')
        MINT(57)=MINT(57)+1
        MINT(51)=1
      ENDIF
 
C...Rescale W(1,*) and W(2,*) (not really necessary, but consistent).
      W(1,1)=W(1,1)*R1
      W(1,2)=W(1,2)/R1
      W(2,1)=W(2,1)/R2
      W(2,2)=W(2,2)*R2
 
C...Rescale BR x values.
      DO 290 IM=MINT(31)+1,MAX(NMI(1),NMI(2))
        XMI(1,IM)=XMI(1,IM)*R1
        XMI(2,IM)=XMI(2,IM)*R2
  290 CONTINUE
 
C...Now we have a consistent set of x and kT values.
C...First set up the initiators and their daughters correctly.
      DO 300 IM=1,MINT(31)
        I1=IMI(1,IM,1)
        I2=IMI(2,IM,1)
        ST=XMI(1,IM)*XMI(2,IM)*VINT(2)+(P(I1,1)+P(I2,1))**2+
     &       (P(I1,2)+P(I2,2))**2
        PT12=P(I1,1)**2+P(I1,2)**2
        PT22=P(I2,1)**2+P(I2,2)**2
C...p_z
        P(I1,3)=SQRT(FLAM(ST,PT12,PT22)/(4D0*ST))
        P(I2,3)=-P(I1,3)
C...Energies (masses should be zero at this stage)
        P(I1,4)=SQRT(PT12+P(I1,3)**2)
        P(I2,4)=SQRT(PT22+P(I2,3)**2)
 
C...Transverse 12 system initiator velocity:
        VB(1)=(P(I1,1)+P(I2,1))/SQRT(ST)
        VB(2)=(P(I1,2)+P(I2,2))/SQRT(ST)
C...Boost to overall initiator system rest frame
        CALL PYROBO(I1,I1,0D0,0D0,-VB(1),-VB(2),0D0)
        CALL PYROBO(I2,I2,0D0,0D0,-VB(1),-VB(2),0D0)

C...Compute phi,theta coordinates of I1 and rotate z axis.
        PHI=PYANGL(P(I1,1),P(I1,2))
        THE=PYANGL(P(I1,3),SQRT(P(I1,1)**2+P(I1,2)**2))
        IMIN=IMISEP(IM-1)+1
C...(include documentation lines if MI = 1)
        IF (IM.EQ.1) IMIN=MINT(83)+5
        IMAX=IMISEP(IM)
C...Rotate entire system in phi
        CALL PYROBO(IMIN,IMAX,0D0,-PHI,0D0,0D0,0D0)
C...Only rotate 12 system in theta
        CALL PYROBO(I1,I1,-THE,0D0,0D0,0D0,0D0)
        CALL PYROBO(I2,I2,-THE,0D0,0D0,0D0,0D0)

C...Now boost entire system back to LAB
        VB(3)=(XMI(1,IM)-XMI(2,IM))/(XMI(1,IM)+XMI(2,IM))
        CALL PYROBO(IMIN,IMAX,THE,PHI,VB(1),VB(2),0D0)
        CALL PYROBO(IMIN,IMAX,0D0,0D0,0D0,0D0,VB(3))

  300 CONTINUE
 
 
C...For the beam remnant partons/hadrons, we only need to set pz and E.
      DO 320 JS=1,2
        DO 310 IM=MINT(31)+1,NMI(JS)
          I=IMI(JS,IM,1)
C...Skip collapsed gluons and junctions.
          IF (K(I,2).EQ.21.AND.K(I,1).EQ.14) GOTO 310
          IF (KFA.EQ.88) GOTO 310
          RMT2=P(I,5)**2+P(I,1)**2+P(I,2)**2
          P(I,4)=0.5D0*(XMI(JS,IM)*W(0,JS)+RMT2/(XMI(JS,IM)*W(0,JS)))
          P(I,3)=0.5D0*(XMI(JS,IM)*W(0,JS)-RMT2/(XMI(JS,IM)*W(0,JS)))
          IF (JS.EQ.2) P(I,3)=-P(I,3)
  310   CONTINUE
  320 CONTINUE
 
 
C...Documentation lines
      DO 340 JS=1,2
        IN=MINT(83)+JS+2
        IO=IMI(JS,1,1)
        K(IN,1)=21
        K(IN,2)=K(IO,2)
        K(IN,3)=MINT(83)+JS
        K(IN,4)=0
        K(IN,5)=0
        DO 330 J=1,5
          P(IN,J)=P(IO,J)
          V(IN,J)=V(IO,J)
  330   CONTINUE
        MCT(IN,1)=MCT(IO,1)
        MCT(IN,2)=MCT(IO,2)
  340 CONTINUE
 
C...Final state colour reconnections.
      IF (MSTP(95).NE.1.OR.MINT(31).LE.1) GOTO 380
 
C...Number of colour tags for which a recoupling will be tried.
      NTOT=NCT
C...Number of recouplings to try
      MINT(34)=0
      NRECP=0
      NITER=0
  350 NRECP=MINT(34)
      NITER=NITER+1
      IITER=0
  360 IITER=IITER+1
      IF (IITER.LE.PARP(78)*NTOT) THEN
C...Select two colour tags at random
C...NB: jj strings do not have colour tags assigned to them,
C...thus they are as yet not affected by anything done here.
        JCT=PYR(0)*NCT+1
        KCT=MOD(INT(JCT+PYR(0)*NCT),NCT)+1
        IJ1=0
        IJ2=0
        IK1=0
        IK2=0
C...Find final state partons with this (anti)colour
        DO 370 I=MINT(84)+1,N
          IF (K(I,1).EQ.3) THEN
            IF (MCT(I,1).EQ.JCT) IJ1=I
            IF (MCT(I,2).EQ.JCT) IJ2=I
            IF (MCT(I,1).EQ.KCT) IK1=I
            IF (MCT(I,2).EQ.KCT) IK2=I
          ENDIF
  370   CONTINUE
C...Only consider recouplings not involving junctions for now.
        IF (IJ1.EQ.0.OR.IJ2.EQ.0.OR.IK1.EQ.0.OR.IK2.EQ.0) GOTO 360
 
        RLO=2D0*FOUR(IJ1,IJ2)*2D0*FOUR(IK1,IK2)
        RLN=2D0*FOUR(IJ1,IK2)*2D0*FOUR(IK1,IJ2)
        IF (RLN.LT.RLO.AND.MCT(IJ2,1).NE.KCT.AND.MCT(IK2,1).NE.JCT) THEN
          MCT(IJ2,2)=KCT
          MCT(IK2,2)=JCT
C...Count up number of reconnections
          MINT(34)=MINT(34)+1
        ENDIF
        IF (MINT(34).LE.1000) THEN
          GOTO 360
        ELSE
          CALL PYERRM(4,'(PYMIRM:) caught in infinite loop')
          GOTO 380
        ENDIF
      ENDIF
      IF (NRECP.LT.MINT(34)) GOTO 350
 
C...Signal PYPREP to use /PYCTAG/ information rather than K(I,KCS).
  380 MINT(33)=1
 
      RETURN
      END
