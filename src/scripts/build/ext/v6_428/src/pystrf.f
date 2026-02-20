 
C*********************************************************************
 
C...PYSTRF
C...Handles the fragmentation of an arbitrary colour singlet
C...jet system according to the Lund string fragmentation model.
 
      SUBROUTINE PYSTRF(IP)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
C...Local arrays. All MOPS variables ends with MO
      DIMENSION DPS(5),KFL(3),PMQ(3),PX(3),PY(3),GAM(3),IE(2),PR(2),
     &IN(9),DHM(4),DHG(4),DP(5,5),IRANK(2),MJU(4),IJU(6),PJU(5,5),
     &TJU(5),KFJH(2),NJS(2),KFJS(2),PJS(4,5),MSTU9T(8),PARU9T(8),
     &INMO(9),PM2QMO(2),XTMO(2),EJSTR(2),IJUORI(2),IBARRK(2),
     &PBST(3,5),TJUOLD(5)
 
C...Function: four-product of two vectors.
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)
      DFOUR(I,J)=DP(I,4)*DP(J,4)-DP(I,1)*DP(J,1)-DP(I,2)*DP(J,2)-
     &DP(I,3)*DP(J,3)
 
C...Reset counters.
      MSTJ(91)=0
      NSAV=N
      MSTU90=MSTU(90)
      NP=0
      KQSUM=0
      DO 100 J=1,5
        DPS(J)=0D0
  100 CONTINUE
      MJU(1)=0
      MJU(2)=0
      NTRYFN=0
      IJUORI(1)=0
      IJUORI(2)=0
 
C...Identify parton system.
      I=IP-1
  110 I=I+1
      IF(I.GT.MIN(N,MSTU(4)-MSTU(32))) THEN
        CALL PYERRM(12,'(PYSTRF:) failed to reconstruct jet system')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      IF(K(I,1).NE.1.AND.K(I,1).NE.2.AND.K(I,1).NE.41) GOTO 110
      KC=PYCOMP(K(I,2))
      IF(KC.EQ.0) GOTO 110
      KQ=KCHG(KC,2)*ISIGN(1,K(I,2))
      IF(KQ.EQ.0.AND.K(I,1).NE.41) GOTO 110
      IF(N+5*NP+11.GT.MSTU(4)-MSTU(32)-5) THEN
        CALL PYERRM(11,'(PYSTRF:) no more memory left in PYJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
 
C...Take copy of partons to be considered. Check flavour sum.
      NP=NP+1
      DO 120 J=1,5
        K(N+NP,J)=K(I,J)
        P(N+NP,J)=P(I,J)
        IF(J.NE.4) DPS(J)=DPS(J)+P(I,J)
  120 CONTINUE
      DPS(4)=DPS(4)+SQRT(P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,5)**2)
      K(N+NP,3)=I
      IF(KQ.NE.2) KQSUM=KQSUM+KQ
      IF(K(I,1).EQ.41) THEN
        IF(MOD(KQSUM,2).EQ.0.AND.MJU(1).EQ.0) THEN
          MJU(1)=N+NP
          IJUORI(1)=I
        ELSE
          MJU(2)=N+NP
          IJUORI(2)=I
        ENDIF
      ENDIF
      IF(K(I,1).EQ.2.OR.K(I,1).EQ.41) GOTO 110
      IF(MOD(KQSUM,3).NE.0) THEN
        CALL PYERRM(12,'(PYSTRF:) unphysical flavour combination')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      IF(MJU(1).GT.0.OR.MJU(2).GT.0) MSTU(29)=1
 
C...Boost copied system to CM frame (for better numerical precision).
      IF(ABS(DPS(3)).LT.0.99D0*DPS(4)) THEN
        MBST=0
        MSTU(33)=1
        CALL PYROBO(N+1,N+NP,0D0,0D0,-DPS(1)/DPS(4),-DPS(2)/DPS(4),
     &  -DPS(3)/DPS(4))
      ELSE
        MBST=1
        HHBZ=SQRT(MAX(1D-6,DPS(4)+DPS(3))/MAX(1D-6,DPS(4)-DPS(3)))
        DO 130 I=N+1,N+NP
          HHPMT=P(I,1)**2+P(I,2)**2+P(I,5)**2
          IF(P(I,3).GT.0D0) THEN
            HHPEZ=MAX(1D-10,(P(I,4)+P(I,3))/HHBZ)
            P(I,3)=0.5D0*(HHPEZ-HHPMT/HHPEZ)
            P(I,4)=0.5D0*(HHPEZ+HHPMT/HHPEZ)
          ELSE
            HHPEZ=MAX(1D-10,(P(I,4)-P(I,3))*HHBZ)
            P(I,3)=-0.5D0*(HHPEZ-HHPMT/HHPEZ)
            P(I,4)=0.5D0*(HHPEZ+HHPMT/HHPEZ)
          ENDIF
  130   CONTINUE
      ENDIF
 
C...Search for very nearby partons that may be recombined.
      NTRYR=0
      NTRYWR=0
      PARU12=PARU(12)
      PARU13=PARU(13)
      MJU(3)=MJU(1)
      MJU(4)=MJU(2)
      NR=NP
      NRMIN=2
      IF(MJU(1).GT.0) NRMIN=NRMIN+2
      IF(MJU(2).GT.0) NRMIN=NRMIN+2
  140 IF(NR.GT.NRMIN) THEN
        PDRMIN=2D0*PARU12
        DO 150 I=N+1,N+NR
          IF(I.EQ.N+NR.AND.IABS(K(N+1,2)).NE.21) GOTO 150
          I1=I+1
          IF(I.EQ.N+NR) I1=N+1
          IF(K(I,1).EQ.41.OR.K(I1,1).EQ.41) GOTO 150
          IF(MJU(1).NE.0.AND.I1.LT.MJU(1).AND.IABS(K(I1,2)).NE.21)
     &    GOTO 150
          IF(MJU(2).NE.0.AND.I.GT.MJU(2).AND.IABS(K(I,2)).NE.21)
     &    GOTO 150
          PAP=SQRT((P(I,1)**2+P(I,2)**2+P(I,3)**2)*(P(I1,1)**2+
     &    P(I1,2)**2+P(I1,3)**2))
          PVP=P(I,1)*P(I1,1)+P(I,2)*P(I1,2)+P(I,3)*P(I1,3)
          PDR=4D0*(PAP-PVP)**2/MAX(1D-6,PARU13**2*PAP+2D0*(PAP-PVP))
          IF(PDR.LT.PDRMIN) THEN
            IR=I
            PDRMIN=PDR
          ENDIF
  150   CONTINUE
 
C...Recombine very nearby partons to avoid machine precision problems.
        IF(PDRMIN.LT.PARU12.AND.IR.EQ.N+NR) THEN
          DO 160 J=1,4
            P(N+1,J)=P(N+1,J)+P(N+NR,J)
  160     CONTINUE
          P(N+1,5)=SQRT(MAX(0D0,P(N+1,4)**2-P(N+1,1)**2-P(N+1,2)**2-
     &    P(N+1,3)**2))
          NR=NR-1
          GOTO 140
        ELSEIF(PDRMIN.LT.PARU12) THEN
          DO 170 J=1,4
            P(IR,J)=P(IR,J)+P(IR+1,J)
  170     CONTINUE
          P(IR,5)=SQRT(MAX(0D0,P(IR,4)**2-P(IR,1)**2-P(IR,2)**2-
     &    P(IR,3)**2))
          IF(MJU(2).NE.0.AND.IR.GT.MJU(2)) K(IR,2)=K(IR+1,2)
          DO 190 I=IR+1,N+NR-1
            K(I,1)=K(I+1,1)
            K(I,2)=K(I+1,2)
            DO 180 J=1,5
              P(I,J)=P(I+1,J)
  180       CONTINUE
  190     CONTINUE
          IF(IR.EQ.N+NR-1) K(IR,2)=K(N+NR,2)
          NR=NR-1
          IF(MJU(1).GT.IR) MJU(1)=MJU(1)-1
          IF(MJU(2).GT.IR) MJU(2)=MJU(2)-1
          GOTO 140
        ENDIF
      ENDIF
      NTRYR=NTRYR+1
 
C...Reset particle counter. Skip ahead if no junctions are present;
C...this is usually the case!
      NRS=MAX(5*NR+11,NP)
      NTRY=0
  200 NTRY=NTRY+1
      IF(NTRY.GT.100.AND.NTRYR.LE.8.AND.NR.GT.NRMIN) THEN
        PARU12=4D0*PARU12
        PARU13=2D0*PARU13
        GOTO 140
      ELSEIF(NTRY.GT.100.OR.NTRYR.GT.100) THEN
        CALL PYERRM(14,'(PYSTRF:) caught in infinite loop')
        IF(MSTU(21).EQ.2) MSTU(90)=0
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      I=N+NRS
      MSTU(90)=MSTU90
      IF(MJU(1).EQ.0.AND.MJU(2).EQ.0) GOTO 650
      IF(MSTJ(12).GE.4) CALL PYERRM(29,'(PYSTRF:) sorry,'//
     &     ' junction strings not handled by MSTJ(12)>3 options')
      DO 640 JT=1,2
        NJS(JT)=0
        IF(MJU(JT).EQ.0) GOTO 640
        JS=3-2*JT
 
C++SKANDS
C...Find and sum up momentum on three sides of junction.
C...Begin with previous boost = zero.
        IJRFIT=0
        DO 210 IX=1,3
          TJUOLD(IX)=0D0
  210   CONTINUE
C...Prevent IJU (specifically IJU(5)) from containing junk below
        DO 215 IU=1,6
          IJU(IU)=0
 215    CONTINUE
        TJUOLD(4)=1D0
  220   IU=0
C...Beginning and end of string system in event record.
        I1BEG=N+1+(JT-1)*(NR-1)
        I1END=N+NR+(JT-1)*(1-NR)
C...Look for junction string piece end points
        DO 230 I1=I1BEG,I1END,JS
          IF(K(I1,2).NE.21.AND.IU.LE.5.AND.IJRFIT.EQ.0) THEN
C...Store junction string piece end points.
C                 1-junction systems        2-junction systems
C           IU :  1     2     3   4     1     2   3     4   5     6
C       IJU(IU):  q-g-g-q-g-g-j-g-q     q-g-g-q-g-j-g-g-j-g-q-g-g-q
            IU=IU+1
            IJU(IU)=I1
          ENDIF
C...Sum over momenta, from junction outwards.
  230   CONTINUE
        DO 280 IU=1,3
          PWT=0D0
C...Initialize junction drag and string piece 4-vectors.
          DO 240 J=1,5
            PBST(IU,J)=0D0
            PJU(IU,J)=0D0
  240     CONTINUE
C...First two branches. Inwards out means opposite direction to JS.
C...(JS is 1 for JT=1, -1 for JT=2)
          IF (IU.LT.3) THEN
            I1A=IJU(IU+1)-JS
            I1B=IJU(IU)
            IDIR=-JS
C...Last branch (gq or gjgqgq). Direction now reversed.
          ELSE
            I1A=IJU(IU)+JS
            I1B=I1END
            IDIR=JS
          ENDIF
          DO 270 I1=I1A,I1B,IDIR
C...Sum up momentum directions with exponential suppression
C...for use in finding junction rest frame below.
            IF (K(I1,2).EQ.88) THEN
C...gjgqgq type system encountered. Use current PWT as start
C...for both strings.
              PWTOLD=PWT
            ELSE
              IF (I1.EQ.IJU(5)+IDIR) PWT=PWTOLD
C...Sum up string piece (boosted) 4-momenta.
              DO 250 J=1,4
                PJU(IU,J)=PJU(IU,J)+P(I1,J)
  250         CONTINUE
C...Compute "junction drag" vectors from (boosted) 4-momenta (initial
C...boost is zero, see above). Skip parton if suppression factor large.
              IF (PWT.GT.10D0) GOTO 270
C...Compute momentum in current frame:
              TDP=TJUOLD(1)*P(I1,1)+TJUOLD(2)*P(I1,2)+TJUOLD(3)*P(I1,3)
              BFC=TDP/(1D0+TJUOLD(4))+P(I1,4)
              DO 260 J=1,3
                PTMP=P(I1,J)+TJUOLD(J)*BFC
                PBST(IU,J)=PBST(IU,J)+PTMP*EXP(-PWT)
  260         CONTINUE
C...Boosted energy
              PTMP=TJUOLD(4)*P(I1,4)+TDP
              PBST(IU,4)=PBST(IU,J)+PTMP*EXP(-PWT)
              PWT=PWT+PTMP/PARJ(48)
            ENDIF
  270     CONTINUE
C...Put |p| rather than m in 5th slot.
          PBST(IU,5)=SQRT(PBST(IU,1)**2+PBST(IU,2)**2+PBST(IU,3)**2)
          PJU(IU,5)=SQRT(PJU(IU,1)**2+PJU(IU,2)**2+PJU(IU,3)**2)
  280   CONTINUE
 
C...Calculate boost from present frame to next JRF candidate.
        IJRFIT=IJRFIT+1
        CALL PYJURF(PBST,TJU)
 
C...After some iterations do not take full step in new direction.
        IF(IJRFIT.GT.5) THEN
          REDUCE=0.8D0**(IJRFIT-5)
          TJU(1)=REDUCE*TJU(1)
          TJU(2)=REDUCE*TJU(2)
          TJU(3)=REDUCE*TJU(3)
          TJU(4)=SQRT(1D0+TJU(1)**2+TJU(2)**2+TJU(3)**2)
        ENDIF
 
C...Combine new boost (TJU) with old boost (TJUOLD)
        TMP=TJU(1)*TJUOLD(1)+TJU(2)*TJUOLD(2)+TJU(3)*TJUOLD(3)
        DO 290 IX=1,3
          TJUOLD(IX)=TJU(IX)+TJUOLD(IX)*(TMP/(1D0+TJUOLD(4))+TJU(4))
  290   CONTINUE
        TJUOLD(4)=SQRT(1D0+TJUOLD(1)**2+TJUOLD(2)**2+TJUOLD(3)**2)
 
C...If last boost small, accept JRF, else iterate.
C...Also prevent possibility of infinite loop.
        IF (ABS((TJU(4)-1D0)/TJUOLD(4)).GT.0.01D0.AND.
     &  IJRFIT.LT.MSTJ(18)) THEN
          GOTO 220
        ELSEIF (IJRFIT.GE.MSTJ(18)) THEN
          CALL PYERRM(1,'(PYSTRF:) failed to converge on JRF')
        ENDIF
 
C...Now store total boost in TJU and change perception.
C...TJUOLD = boost vector from CM of string syst -> JRF. Henceforth,
C...TJU = junction motion vector in string CM, so the sign changes.
        DO 300 J=1,3
          TJU(J)=-TJUOLD(J)
  300   CONTINUE
        TJU(4)=SQRT(1D0+TJU(1)**2+TJU(2)**2+TJU(3)**2)
 
C--SKANDS
 
C...Calculate string piece energies in junction rest frame.
        DO 310 IU=1,3
          PJU(IU,5)=TJU(4)*PJU(IU,4)-TJU(1)*PJU(IU,1)-TJU(2)*PJU(IU,2)-
     &    TJU(3)*PJU(IU,3)
          PBST(IU,5)=TJU(4)*PBST(IU,4)-TJU(1)*PBST(IU,1)-
     &    TJU(2)*PBST(IU,2)-TJU(3)*PBST(IU,3)
  310   CONTINUE
 
C...Start preparing for fragmentation of two strings from junction.
        ISTA=I
        NTRYER=0
  320   NTRYER=NTRYER+1
        MSTU(90)=MSTU90
        I=ISTA
        DO 620 IU=1,2
          NS=IABS(IJU(IU+1)-IJU(IU))
 
C...Junction strings: find longitudinal string directions.
          DO 350 IS=1,NS
            IS1=IJU(IU)+JS*(IS-1)
            IS2=IJU(IU)+JS*IS
            DO 330 J=1,5
              DP(1,J)=0.5D0*P(IS1,J)
              IF(IS.EQ.1) DP(1,J)=P(IS1,J)
              DP(2,J)=0.5D0*P(IS2,J)
              IF(IS.EQ.NS) DP(2,J)=(-PBST(IU,J)+2D0*PBST(IU,5)*TJU(J))*
     &        (PJU(IU,5)/PBST(IU,5))
  330       CONTINUE
            IF(IS.EQ.NS) DP(2,5)=SQRT(MAX(0D0,PJU(IU,4)**2-
     &      PJU(IU,1)**2-PJU(IU,2)**2-PJU(IU,3)**2))
            DP(3,5)=DFOUR(1,1)
            DP(4,5)=DFOUR(2,2)
            DHKC=DFOUR(1,2)
            IF(DP(3,5)+2D0*DHKC+DP(4,5).LE.0D0) THEN
              DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
              DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
              DP(3,5)=0D0
              DP(4,5)=0D0
              DHKC=DFOUR(1,2)
            ENDIF
            DHKS=SQRT(DHKC**2-DP(3,5)*DP(4,5))
            DHK1=0.5D0*((DP(4,5)+DHKC)/DHKS-1D0)
            DHK2=0.5D0*((DP(3,5)+DHKC)/DHKS-1D0)
            IN1=N+NR+4*IS-3
            P(IN1,5)=SQRT(DP(3,5)+2D0*DHKC+DP(4,5))
            DO 340 J=1,4
              P(IN1,J)=(1D0+DHK1)*DP(1,J)-DHK2*DP(2,J)
              P(IN1+1,J)=(1D0+DHK2)*DP(2,J)-DHK1*DP(1,J)
  340       CONTINUE
  350     CONTINUE
 
C...Junction strings: initialize flavour, momentum and starting pos.
          ISAV=I
          MSTU91=MSTU(90)
  360     NTRY=NTRY+1
          IF(NTRY.GT.100.AND.NTRYR.LE.8.AND.NR.GT.NRMIN) THEN
            PARU12=4D0*PARU12
            PARU13=2D0*PARU13
            GOTO 140
          ELSEIF(NTRY.GT.100) THEN
            CALL PYERRM(14,'(PYSTRF:) caught in infinite loop')
            IF(MSTU(21).EQ.2) MSTU(90)=0
            IF(MSTU(21).GE.1) RETURN
          ENDIF
          I=ISAV
          MSTU(90)=MSTU91
          IRANKJ=0
          IE(1)=K(N+1+(JT/2)*(NP-1),3)
          IF (MOD(JT+IU,2).NE.0) THEN
            IE(1)=K(IJU(IU),3)
            IF (NP-NR.NE.0) THEN
C...If gluons have disappeared. Original IJU must be used.
              IT=IP
              NE=1
  370         IT=IT+1
              IF (K(IT,2).NE.21) THEN
                NE=NE+1
              ENDIF
              IF (NE.EQ.IU+4*(JT-1)) THEN
                IE(1)=IT
              ELSEIF (IT.LE.IP+NP) THEN
                GOTO 370
              ELSE
                CALL PYERRM(14,'(PYSTRF:) '//
     &               'Original IJU could not be reconstructed!')
              ENDIF
            ENDIF
          ENDIF
          IN(4)=N+NR+1
          IN(5)=IN(4)+1
          IN(6)=N+NR+4*NS+1
          DO 390 JQ=1,2
            DO 380 IN1=N+NR+2+JQ,N+NR+4*NS-2+JQ,4
              P(IN1,1)=2-JQ
              P(IN1,2)=JQ-1
              P(IN1,3)=1D0
  380       CONTINUE
  390     CONTINUE
          KFL(1)=K(IJU(IU),2)
          PX(1)=0D0
          PY(1)=0D0
          GAM(1)=0D0
          DO 400 J=1,5
            PJU(IU+3,J)=0D0
  400     CONTINUE
 
C...Junction strings: find initial transverse directions.
          DO 410 J=1,4
            DP(1,J)=P(IN(4),J)
            DP(2,J)=P(IN(4)+1,J)
            DP(3,J)=0D0
            DP(4,J)=0D0
  410     CONTINUE
          DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
          DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
          DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4)
          DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4)
          DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4)
          IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1D0
          IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1D0
          IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1D0
          IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1D0
          DHC12=DFOUR(1,2)
          DHCX1=DFOUR(3,1)/DHC12
          DHCX2=DFOUR(3,2)/DHC12
          DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12)
          DHCY1=DFOUR(4,1)/DHC12
          DHCY2=DFOUR(4,2)/DHC12
          DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12
          DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2)
          DO 420 J=1,4
            DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
            P(IN(6),J)=DP(3,J)
            P(IN(6)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     &      DHCYX*DP(3,J))
  420     CONTINUE
 
C...Junction strings: produce new particle, origin.
  430     I=I+1
          IF(2*I-NSAV.GE.MSTU(4)-MSTU(32)-5) THEN
            CALL PYERRM(11,'(PYSTRF:) no more memory left in PYJETS')
            IF(MSTU(21).GE.1) RETURN
          ENDIF
          IRANKJ=IRANKJ+1
          K(I,1)=1
          K(I,3)=IE(1)
          K(I,4)=0
          K(I,5)=0
 
C...Junction strings: generate flavour, hadron, pT, z and Gamma.
  440     CALL PYKFDI(KFL(1),0,KFL(3),K(I,2))
          IF(K(I,2).EQ.0) GOTO 360
          IF(IRANKJ.EQ.1.AND.IABS(KFL(1)).LE.10.AND.
     &    IABS(KFL(3)).GT.10) THEN
            IF(PYR(0).GT.PARJ(19)) GOTO 440
          ENDIF
          P(I,5)=PYMASS(K(I,2))
          CALL PYPTDI(KFL(1),PX(3),PY(3))
          PR(1)=P(I,5)**2+(PX(1)+PX(3))**2+(PY(1)+PY(3))**2
          CALL PYZDIS(KFL(1),KFL(3),PR(1),Z)
          IF(IABS(KFL(1)).GE.4.AND.IABS(KFL(1)).LE.8.AND.
     &    MSTU(90).LT.8) THEN
            MSTU(90)=MSTU(90)+1
            MSTU(90+MSTU(90))=I
            PARU(90+MSTU(90))=Z
          ENDIF
          GAM(3)=(1D0-Z)*(GAM(1)+PR(1)/Z)
          DO 450 J=1,3
            IN(J)=IN(3+J)
  450     CONTINUE
 
C...Junction strings: stepping within 'low' string region.
          IF(IN(1)+1.EQ.IN(2).AND.Z*P(IN(1)+2,3)*P(IN(2)+2,3)*
     &    P(IN(1),5)**2.GE.PR(1)) THEN
            P(IN(1)+2,4)=Z*P(IN(1)+2,3)
            P(IN(2)+2,4)=PR(1)/(P(IN(1)+2,4)*P(IN(1),5)**2)
            DO 460 J=1,4
              P(I,J)=(PX(1)+PX(3))*P(IN(3),J)+(PY(1)+PY(3))*P(IN(3)+1,J)
  460       CONTINUE
            GOTO 560
C...Has used up energy of junction string, i.e. no more hadrons in it.
          ELSEIF(IN(1)+1.EQ.IN(2).AND.IN(1).EQ.N+NR+4*NS-3) THEN
            DO 470 J=1,5
              P(I,J)=0D0
  470       CONTINUE
            GOTO 600
C...Stepping from 'low' string region
          ELSEIF(IN(1)+1.EQ.IN(2)) THEN
            P(IN(2)+2,4)=P(IN(2)+2,3)
            P(IN(2)+2,1)=1D0
            IN(2)=IN(2)+4
            IF(IN(2).GT.N+NR+4*NS) GOTO 360
            IF(FOUR(IN(1),IN(2)).LE.1D-2) THEN
              P(IN(1)+2,4)=P(IN(1)+2,3)
              P(IN(1)+2,1)=0D0
              IN(1)=IN(1)+4
            ENDIF
          ENDIF
 
C...Junction strings: find new transverse directions.
  480     IF(IN(1).GT.N+NR+4*NS.OR.IN(2).GT.N+NR+4*NS.OR.
     &    IN(1).GT.IN(2)) GOTO 360
          IF(IN(1).NE.IN(4).OR.IN(2).NE.IN(5)) THEN
            DO 490 J=1,4
              DP(1,J)=P(IN(1),J)
              DP(2,J)=P(IN(2),J)
              DP(3,J)=0D0
              DP(4,J)=0D0
  490       CONTINUE
            DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
            DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
            DHC12=DFOUR(1,2)
            IF(DHC12.LE.1D-2) THEN
              P(IN(1)+2,4)=P(IN(1)+2,3)
              P(IN(1)+2,1)=0D0
              IN(1)=IN(1)+4
              GOTO 480
            ENDIF
            IN(3)=N+NR+4*NS+5
            DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4)
            DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4)
            DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4)
            IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1D0
            IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1D0
            IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1D0
            IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1D0
            DHCX1=DFOUR(3,1)/DHC12
            DHCX2=DFOUR(3,2)/DHC12
            DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12)
            DHCY1=DFOUR(4,1)/DHC12
            DHCY2=DFOUR(4,2)/DHC12
            DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12
            DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2)
            DO 500 J=1,4
              DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
              P(IN(3),J)=DP(3,J)
              P(IN(3)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     &        DHCYX*DP(3,J))
  500       CONTINUE
C...Express pT with respect to new axes, if sensible.
            PXP=-(PX(3)*FOUR(IN(6),IN(3))+PY(3)*FOUR(IN(6)+1,IN(3)))
            PYP=-(PX(3)*FOUR(IN(6),IN(3)+1)+PY(3)*FOUR(IN(6)+1,IN(3)+1))
            IF(ABS(PXP**2+PYP**2-PX(3)**2-PY(3)**2).LT.0.01D0) THEN
              PX(3)=PXP
              PY(3)=PYP
            ENDIF
          ENDIF
 
C...Junction strings: sum up known four-momentum, coefficients for m2.
          DO 530 J=1,4
            DHG(J)=0D0
            P(I,J)=PX(1)*P(IN(6),J)+PY(1)*P(IN(6)+1,J)+PX(3)*P(IN(3),J)+
     &      PY(3)*P(IN(3)+1,J)
            DO 510 IN1=IN(4),IN(1)-4,4
              P(I,J)=P(I,J)+P(IN1+2,3)*P(IN1,J)
  510       CONTINUE
            DO 520 IN2=IN(5),IN(2)-4,4
              P(I,J)=P(I,J)+P(IN2+2,3)*P(IN2,J)
  520       CONTINUE
  530     CONTINUE
          DHM(1)=FOUR(I,I)
          DHM(2)=2D0*FOUR(I,IN(1))
          DHM(3)=2D0*FOUR(I,IN(2))
          DHM(4)=2D0*FOUR(IN(1),IN(2))
 
C...Junction strings: find coefficients for Gamma expression.
          DO 550 IN2=IN(1)+1,IN(2),4
            DO 540 IN1=IN(1),IN2-1,4
              DHC=2D0*FOUR(IN1,IN2)
              DHG(1)=DHG(1)+P(IN1+2,1)*P(IN2+2,1)*DHC
              IF(IN1.EQ.IN(1)) DHG(2)=DHG(2)-P(IN2+2,1)*DHC
              IF(IN2.EQ.IN(2)) DHG(3)=DHG(3)+P(IN1+2,1)*DHC
              IF(IN1.EQ.IN(1).AND.IN2.EQ.IN(2)) DHG(4)=DHG(4)-DHC
  540       CONTINUE
  550     CONTINUE
 
C...Junction strings: solve (m2, Gamma) equation system for energies.
          DHS1=DHM(3)*DHG(4)-DHM(4)*DHG(3)
          IF(ABS(DHS1).LT.1D-4) GOTO 360
          DHS2=DHM(4)*(GAM(3)-DHG(1))-DHM(2)*DHG(3)-DHG(4)*
     &    (P(I,5)**2-DHM(1))+DHG(2)*DHM(3)
          DHS3=DHM(2)*(GAM(3)-DHG(1))-DHG(2)*(P(I,5)**2-DHM(1))
          P(IN(2)+2,4)=0.5D0*(SQRT(MAX(0D0,DHS2**2-4D0*DHS1*DHS3))/
     &    ABS(DHS1)-DHS2/DHS1)
          IF(DHM(2)+DHM(4)*P(IN(2)+2,4).LE.0D0) GOTO 360
          P(IN(1)+2,4)=(P(I,5)**2-DHM(1)-DHM(3)*P(IN(2)+2,4))/
     &    (DHM(2)+DHM(4)*P(IN(2)+2,4))
 
C...Junction strings: step to new region if necessary.
          IF(P(IN(2)+2,4).GT.P(IN(2)+2,3)) THEN
            P(IN(2)+2,4)=P(IN(2)+2,3)
            P(IN(2)+2,1)=1D0
            IN(2)=IN(2)+4
            IF(IN(2).GT.N+NR+4*NS) GOTO 360
            IF(FOUR(IN(1),IN(2)).LE.1D-2) THEN
              P(IN(1)+2,4)=P(IN(1)+2,3)
              P(IN(1)+2,1)=0D0
              IN(1)=IN(1)+4
            ENDIF
            GOTO 480
          ELSEIF(P(IN(1)+2,4).GT.P(IN(1)+2,3)) THEN
            P(IN(1)+2,4)=P(IN(1)+2,3)
            P(IN(1)+2,1)=0D0
            IN(1)=IN(1)+4
            GOTO 480
          ENDIF
 
C...Junction strings: particle four-momentum, remainder, loop back.
  560     DO 570 J=1,4
            P(I,J)=P(I,J)+P(IN(1)+2,4)*P(IN(1),J)+
     &      P(IN(2)+2,4)*P(IN(2),J)
            PJU(IU+3,J)=PJU(IU+3,J)+P(I,J)
  570     CONTINUE
          IF(P(I,4).LT.P(I,5)) GOTO 360
          PJU(IU+3,5)=TJU(4)*PJU(IU+3,4)-TJU(1)*PJU(IU+3,1)-
     &    TJU(2)*PJU(IU+3,2)-TJU(3)*PJU(IU+3,3)
          IF(PJU(IU+3,5).LT.PJU(IU,5)) THEN
            KFL(1)=-KFL(3)
            PX(1)=-PX(3)
            PY(1)=-PY(3)
            GAM(1)=GAM(3)
            IF(IN(3).NE.IN(6)) THEN
              DO 580 J=1,4
                P(IN(6),J)=P(IN(3),J)
                P(IN(6)+1,J)=P(IN(3)+1,J)
  580         CONTINUE
            ENDIF
            DO 590 JQ=1,2
              IN(3+JQ)=IN(JQ)
              P(IN(JQ)+2,3)=P(IN(JQ)+2,3)-P(IN(JQ)+2,4)
              P(IN(JQ)+2,1)=P(IN(JQ)+2,1)-(3-2*JQ)*P(IN(JQ)+2,4)
  590       CONTINUE
            GOTO 430
          ENDIF
 
C...Junction strings: save quantities left after each string.
          IF(IABS(KFL(1)).GT.10) GOTO 360
  600     I=I-1
          IF(MSTU(90+MSTU(90)).EQ.I+1) MSTU(90)=MSTU(90)-1 
          KFJH(IU)=KFL(1)
          DO 610 J=1,4
            PJU(IU+3,J)=PJU(IU+3,J)-P(I+1,J)
  610     CONTINUE
 
C...Junction strings: loopback if much unused energy in both strings.
          PJU(IU+3,5)=TJU(4)*PJU(IU+3,4)-TJU(1)*PJU(IU+3,1)-
     &    TJU(2)*PJU(IU+3,2)-TJU(3)*PJU(IU+3,3)
          EJSTR(IU)=PJU(IU,5)-PJU(IU+3,5)
  620   CONTINUE
        IF((MIN(EJSTR(1),EJSTR(2)).GT.PARJ(49).OR.
     &  EJSTR(1).GT.PARJ(49)+PYR(0)*PARJ(50).OR.
     &  EJSTR(2).GT.PARJ(49)+PYR(0)*PARJ(50))
     &  .AND.NTRYER.LT.10) GOTO 320
 
C...Junction strings: put together to new effective string endpoint.
        NJS(JT)=I-ISTA
        KFLS=2*INT(PYR(0)+3D0*PARJ(4)/(1D0+3D0*PARJ(4)))+1
        IF(KFJH(1).EQ.KFJH(2)) KFLS=3
        KFJS(JT)=ISIGN(1000*MAX(IABS(KFJH(1)),IABS(KFJH(2)))+
     &  100*MIN(IABS(KFJH(1)),IABS(KFJH(2)))+KFLS,KFJH(1))
        DO 630 J=1,4
          PJS(JT,J)=PJU(1,J)+PJU(2,J)+P(MJU(JT),J)
          PJS(JT+2,J)=PJU(4,J)+PJU(5,J)
  630   CONTINUE
        PJS(JT,5)=SQRT(MAX(0D0,PJS(JT,4)**2-PJS(JT,1)**2-PJS(JT,2)**2-
     &  PJS(JT,3)**2))
        PJS(JT+2,5)=0D0
  640 CONTINUE
 
C...Open versus closed strings. Choose breakup region for latter.
  650 IF(MJU(1).NE.0.AND.MJU(2).NE.0) THEN
        NS=MJU(2)-MJU(1)
        NB=MJU(1)-N
      ELSEIF(MJU(1).NE.0) THEN
        NS=N+NR-MJU(1)
        NB=MJU(1)-N
      ELSEIF(MJU(2).NE.0) THEN
        NS=MJU(2)-N
        NB=1
      ELSEIF(IABS(K(N+1,2)).NE.21) THEN
        NS=NR-1
        NB=1
      ELSE
        NS=NR+1
        W2SUM=0D0
        DO 660 IS=1,NR
          P(N+NR+IS,1)=0.5D0*FOUR(N+IS,N+IS+1-NR*(IS/NR))
          W2SUM=W2SUM+P(N+NR+IS,1)
  660   CONTINUE
        W2RAN=PYR(0)*W2SUM
        NB=0
  670   NB=NB+1
        W2SUM=W2SUM-P(N+NR+NB,1)
        IF(W2SUM.GT.W2RAN.AND.NB.LT.NR) GOTO 670
      ENDIF
 
C...Find longitudinal string directions (i.e. lightlike four-vectors).
      DO 700 IS=1,NS
        IS1=N+IS+NB-1-NR*((IS+NB-2)/NR)
        IS2=N+IS+NB-NR*((IS+NB-1)/NR)
        DO 680 J=1,5
          DP(1,J)=P(IS1,J)
          IF(IABS(K(IS1,2)).EQ.21) DP(1,J)=0.5D0*DP(1,J)
          IF(IS1.EQ.MJU(1)) DP(1,J)=PJS(1,J)-PJS(3,J)
          DP(2,J)=P(IS2,J)
          IF(IABS(K(IS2,2)).EQ.21) DP(2,J)=0.5D0*DP(2,J)
          IF(IS2.EQ.MJU(2)) DP(2,J)=PJS(2,J)-PJS(4,J)
  680   CONTINUE
        IF(IS1.EQ.MJU(1)) DP(1,5)=SQRT(MAX(0D0,DP(1,4)**2-DP(1,1)**2-
     &  DP(1,2)**2-DP(1,3)**2))
        IF(IS2.EQ.MJU(2)) DP(2,5)=SQRT(MAX(0D0,DP(2,4)**2-DP(2,1)**2-
     &  DP(2,2)**2-DP(2,3)**2))
        DP(3,5)=DFOUR(1,1)
        DP(4,5)=DFOUR(2,2)
        DHKC=DFOUR(1,2)
        IF(DP(3,5)+2D0*DHKC+DP(4,5).LE.0D0) GOTO 200
        DHKS=SQRT(DHKC**2-DP(3,5)*DP(4,5))
        DHK1=0.5D0*((DP(4,5)+DHKC)/DHKS-1D0)
        DHK2=0.5D0*((DP(3,5)+DHKC)/DHKS-1D0)
        IN1=N+NR+4*IS-3
        P(IN1,5)=SQRT(DP(3,5)+2D0*DHKC+DP(4,5))
        DO 690 J=1,4
          P(IN1,J)=(1D0+DHK1)*DP(1,J)-DHK2*DP(2,J)
          P(IN1+1,J)=(1D0+DHK2)*DP(2,J)-DHK1*DP(1,J)
  690   CONTINUE
  700 CONTINUE
 
C...Begin initialization: sum up energy, set starting position.
      ISAV=I
      MSTU91=MSTU(90)
  710 NTRY=NTRY+1
      IF(NTRY.GT.100.AND.NTRYR.LE.8.AND.NR.GT.NRMIN) THEN
        PARU12=4D0*PARU12
        PARU13=2D0*PARU13
        GOTO 140
      ELSEIF(NTRY.GT.100) THEN
        CALL PYERRM(14,'(PYSTRF:) caught in infinite loop')
        IF(MSTU(21).EQ.2) MSTU(90)=0
        IF(MSTU(21).GE.1) RETURN
      ENDIF
      I=ISAV
      MSTU(90)=MSTU91
      DO 730 J=1,4
        P(N+NRS,J)=0D0
        DO 720 IS=1,NR
          P(N+NRS,J)=P(N+NRS,J)+P(N+IS,J)
  720   CONTINUE
  730 CONTINUE
      DO 750 JT=1,2
        IRANK(JT)=0
        IF(MJU(JT).NE.0) IRANK(JT)=NJS(JT)
        IF(NS.GT.NR) IRANK(JT)=1
        IBARRK(JT)=0
        IE(JT)=K(N+1+(JT/2)*(NP-1),3)
        IN(3*JT+1)=N+NR+1+4*(JT/2)*(NS-1)
        IN(3*JT+2)=IN(3*JT+1)+1
        IN(3*JT+3)=N+NR+4*NS+2*JT-1
        DO 740 IN1=N+NR+2+JT,N+NR+4*NS-2+JT,4
          P(IN1,1)=2-JT
          P(IN1,2)=JT-1
          P(IN1,3)=1D0
  740   CONTINUE
  750 CONTINUE
 
C.. MOPS variables and switches
      NRVMO=0
      XBMO=1D0
      MSTU(121)=0
      MSTU(122)=0
 
C...Initialize flavour and pT variables for open string.
      IF(NS.LT.NR) THEN
        PX(1)=0D0
        PY(1)=0D0
        IF(NS.EQ.1.AND.MJU(1)+MJU(2).EQ.0) CALL PYPTDI(0,PX(1),PY(1))
        PX(2)=-PX(1)
        PY(2)=-PY(1)
        DO 760 JT=1,2
          KFL(JT)=K(IE(JT),2)
          IF(MJU(JT).NE.0) KFL(JT)=KFJS(JT)
          IF(MJU(JT).NE.0.AND.IABS(KFL(JT)).GT.1000) IBARRK(JT)=1
          MSTJ(93)=1
          PMQ(JT)=PYMASS(KFL(JT))
          GAM(JT)=0D0
  760   CONTINUE
 
C...Closed string: random initial breakup flavour, pT and vertex.
      ELSE
        KFL(3)=INT(1D0+(2D0+PARJ(2))*PYR(0))*(-1)**INT(PYR(0)+0.5D0)
        IBMO=0
  770   CALL PYKFDI(KFL(3),0,KFL(1),KDUMP)
C.. Closed string: first vertex diq attempt => enforced second
C.. vertex diq
        IF(IABS(KFL(1)).GT.10)THEN
           IBMO=1
           MSTU(121)=0
           GOTO 770
        ENDIF
        IF(IBMO.EQ.1) MSTU(121)=-1
        KFL(2)=-KFL(1)
        CALL PYPTDI(KFL(1),PX(1),PY(1))
        PX(2)=-PX(1)
        PY(2)=-PY(1)
        PR3=MIN(25D0,0.1D0*P(N+NR+1,5)**2)
  780   CALL PYZDIS(KFL(1),KFL(2),PR3,Z)
        ZR=PR3/(Z*P(N+NR+1,5)**2)
        IF(ZR.GE.1D0) GOTO 780
        DO 790 JT=1,2
          MSTJ(93)=1
          PMQ(JT)=PYMASS(KFL(JT))
          GAM(JT)=PR3*(1D0-Z)/Z
          IN1=N+NR+3+4*(JT/2)*(NS-1)
          P(IN1,JT)=1D0-Z
          P(IN1,3-JT)=JT-1
          P(IN1,3)=(2-JT)*(1D0-Z)+(JT-1)*Z
          P(IN1+1,JT)=ZR
          P(IN1+1,3-JT)=2-JT
          P(IN1+1,3)=(2-JT)*(1D0-ZR)+(JT-1)*ZR
  790   CONTINUE
      ENDIF
C.. MOPS variables
      DO 800 JT=1,2
         XTMO(JT)=1D0
         PM2QMO(JT)=PMQ(JT)**2
         IF(IABS(KFL(JT)).GT.10) PM2QMO(JT)=0D0
  800 CONTINUE
 
C...Find initial transverse directions (i.e. spacelike four-vectors).
      DO 840 JT=1,2
        IF(JT.EQ.1.OR.NS.EQ.NR-1.OR.MJU(1)+MJU(2).NE.0) THEN
          IN1=IN(3*JT+1)
          IN3=IN(3*JT+3)
          DO 810 J=1,4
            DP(1,J)=P(IN1,J)
            DP(2,J)=P(IN1+1,J)
            DP(3,J)=0D0
            DP(4,J)=0D0
  810     CONTINUE
          DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
          DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
          DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4)
          DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4)
          DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4)
          IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1D0
          IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1D0
          IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1D0
          IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1D0
          DHC12=DFOUR(1,2)
          DHCX1=DFOUR(3,1)/DHC12
          DHCX2=DFOUR(3,2)/DHC12
          DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12)
          DHCY1=DFOUR(4,1)/DHC12
          DHCY2=DFOUR(4,2)/DHC12
          DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12
          DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2)
          DO 820 J=1,4
            DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
            P(IN3,J)=DP(3,J)
            P(IN3+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     &      DHCYX*DP(3,J))
  820     CONTINUE
        ELSE
          DO 830 J=1,4
            P(IN3+2,J)=P(IN3,J)
            P(IN3+3,J)=P(IN3+1,J)
  830     CONTINUE
        ENDIF
  840 CONTINUE
 
C...Remove energy used up in junction string fragmentation.
      IF(MJU(1)+MJU(2).GT.0) THEN
        DO 860 JT=1,2
          IF(NJS(JT).EQ.0) GOTO 860
          DO 850 J=1,4
            P(N+NRS,J)=P(N+NRS,J)-PJS(JT+2,J)
  850     CONTINUE
  860   CONTINUE
        PARJST=PARJ(33)
        IF(MSTJ(11).EQ.2) PARJST=PARJ(34)
        WMIN=PARJST+PMQ(1)+PMQ(2)
        WREM2=FOUR(N+NRS,N+NRS)
        IF(P(N+NRS,4).LT.0D0.OR.WREM2.LT.WMIN**2) THEN
          NTRYWR=NTRYWR+1
          IF(MOD(NTRYWR,20).NE.0) NTRYR=NTRYR-1
          GOTO 140
        ENDIF
      ENDIF
 
C...Produce new particle: side, origin.
  870 I=I+1
      IF(2*I-NSAV.GE.MSTU(4)-MSTU(32)-5) THEN
        CALL PYERRM(11,'(PYSTRF:) no more memory left in PYJETS')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
C.. New side priority for popcorn systems
      IF(MSTU(121).LE.0)THEN
         JT=1.5D0+PYR(0)
         IF(IABS(KFL(3-JT)).GT.10) JT=3-JT
         IF(IABS(KFL(3-JT)).GE.4.AND.IABS(KFL(3-JT)).LE.8) JT=3-JT
      ENDIF
      JR=3-JT
      JS=3-2*JT
      IRANK(JT)=IRANK(JT)+1
      K(I,1)=1
      K(I,4)=0
      K(I,5)=0
 
C...Generate flavour, hadron and pT.
  880 K(I,3)=IE(JT)
      CALL PYKFDI(KFL(JT),0,KFL(3),K(I,2))
      IF(K(I,2).EQ.0) GOTO 710
      MU90MO=MSTU(90)
      IF(MSTU(121).EQ.-1) GOTO 910
      IF(IRANK(JT).EQ.1.AND.IABS(KFL(JT)).LE.10.AND.
     &IABS(KFL(3)).GT.10) THEN
        IF(PYR(0).GT.PARJ(19)) GOTO 880
      ENDIF
      IF(IBARRK(JT).EQ.1.AND.MOD(IABS(K(I,2)),10000).GT.1000)
     &K(I,3)=IJUORI(JT)
      P(I,5)=PYMASS(K(I,2))
      CALL PYPTDI(KFL(JT),PX(3),PY(3))
      PR(JT)=P(I,5)**2+(PX(JT)+PX(3))**2+(PY(JT)+PY(3))**2
 
C...Final hadrons for small invariant mass.
      MSTJ(93)=1
      PMQ(3)=PYMASS(KFL(3))
      PARJST=PARJ(33)
      IF(MSTJ(11).EQ.2) PARJST=PARJ(34)
      WMIN=PARJST+PMQ(1)+PMQ(2)+PARJ(36)*PMQ(3)
      IF(IABS(KFL(JT)).GT.10.AND.IABS(KFL(3)).GT.10) WMIN=
     &WMIN-0.5D0*PARJ(36)*PMQ(3)
      WREM2=FOUR(N+NRS,N+NRS)
      IF(WREM2.LT.0.10D0) GOTO 710
      IF(WREM2.LT.MAX(WMIN*(1D0+(2D0*PYR(0)-1D0)*PARJ(37)),
     &PARJ(32)+PMQ(1)+PMQ(2))**2) GOTO 1080
 
C...Choose z, which gives Gamma. Shift z for heavy flavours.
      CALL PYZDIS(KFL(JT),KFL(3),PR(JT),Z)
      IF(IABS(KFL(JT)).GE.4.AND.IABS(KFL(JT)).LE.8.AND.
     &MSTU(90).LT.8) THEN
        MSTU(90)=MSTU(90)+1
        MSTU(90+MSTU(90))=I
        PARU(90+MSTU(90))=Z
      ENDIF
      KFL1A=IABS(KFL(1))
      KFL2A=IABS(KFL(2))
      IF(MAX(MOD(KFL1A,10),MOD(KFL1A/1000,10),MOD(KFL2A,10),
     &MOD(KFL2A/1000,10)).GE.4) THEN
        PR(JR)=(PMQ(JR)+PMQ(3))**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2
        PW12=SQRT(MAX(0D0,(WREM2-PR(1)-PR(2))**2-4D0*PR(1)*PR(2)))
        Z=(WREM2+PR(JT)-PR(JR)+PW12*(2D0*Z-1D0))/(2D0*WREM2)
        PR(JR)=(PMQ(JR)+PARJST)**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2
        IF((1D0-Z)*(WREM2-PR(JT)/Z).LT.PR(JR)) GOTO 1080
      ENDIF
      GAM(3)=(1D0-Z)*(GAM(JT)+PR(JT)/Z)
 
C.. MOPS baryon model modification
      XTMO3=(1D0-Z)*XTMO(JT)
      IF(IABS(KFL(3)).LE.10) NRVMO=0
      IF(IABS(KFL(3)).GT.10.AND.MSTJ(12).GE.4) THEN
         GTSTMO=1D0
         PTSTMO=1D0
         RTSTMO=PYR(0)
         IF(IABS(KFL(JT)).LE.10)THEN
            XBMO=MIN(XTMO3,1D0-(2D-10))
            GBMO=GAM(3)
            PMMO=0D0
            PGMO=GBMO+LOG(1D0-XBMO)*PM2QMO(JT)
            GTSTMO=1D0-PARF(192)**PGMO
         ELSE
            IF(IRANK(JT).EQ.1) THEN
               GBMO=GAM(JT)
               PMMO=0D0
               XBMO=1D0
            ENDIF
            IF(XBMO.LT.1D0-(1D-10))THEN
               PGNMO=GBMO*XTMO3/XBMO+PM2QMO(JT)*LOG(1D0-XTMO3)
               GTSTMO=(1D0-PARF(192)**PGNMO)/(1D0-PARF(192)**PGMO)
               PGMO=PGNMO
            ENDIF
            IF(MSTJ(12).GE.5)THEN
               PMNMO=SQRT((XBMO-XTMO3)*(GAM(3)/XTMO3-GBMO/XBMO))
               PMMO=PMMO+PMAS(PYCOMP(K(I,2)),1)-PMAS(PYCOMP(K(I,2)),3)
               PTSTMO=EXP((PMMO-PMNMO)*PARF(193))
               PMMO=PMNMO
            ENDIF
         ENDIF
 
C.. MOPS Accepting popcorn system hadron.
         IF(PTSTMO*GTSTMO.GT.RTSTMO) THEN
            IF(IRANK(JT).EQ.1.OR.IABS(KFL(JT)).LE.10) THEN
               NRVMO=I-N-NR
               IF(I+NRVMO.GT.MSTU(4)-MSTU(32)-5) THEN
                  CALL PYERRM(11,
     &                 '(PYSTRF:) no more memory left in PYJETS')
                  IF(MSTU(21).GE.1) RETURN
               ENDIF
               IMO=I
               KFLMO=KFL(JT)
               PMQMO=PMQ(JT)
               PXMO=PX(JT)
               PYMO=PY(JT)
               GAMMO=GAM(JT)
               IRMO=IRANK(JT)
               XMO=XTMO(JT)
               DO 900 J=1,9
                  IF(J.LE.5) THEN
                     DO 890 LINE=1,I-N-NR
                        P(MSTU(4)-MSTU(32)-LINE,J)=P(N+NR+LINE,J)
                        K(MSTU(4)-MSTU(32)-LINE,J)=K(N+NR+LINE,J)
  890                CONTINUE
                  ENDIF
                  INMO(J)=IN(J)
  900          CONTINUE
            ENDIF
         ELSE
C..Reject popcorn system, flag=-1 if enforcing new one
            MSTU(121)=-1
            IF(PTSTMO.GT.RTSTMO) MSTU(121)=-2
         ENDIF
      ENDIF
 
 
C..Lift restoring string outside MOPS block
  910 IF(MSTU(121).LT.0) THEN
         IF(MSTU(121).EQ.-2) MSTU(121)=0
         MSTU(90)=MU90MO
         NRVMO=0
         IF(IRANK(JT).EQ.1.OR.IABS(KFL(JT)).LE.10) GOTO 880
         I=IMO
         KFL(JT)=KFLMO
         PMQ(JT)=PMQMO
         PX(JT)=PXMO
         PY(JT)=PYMO
         GAM(JT)=GAMMO
         IRANK(JT)=IRMO
         XTMO(JT)=XMO
         DO 930 J=1,9
            IF(J.LE.5) THEN
               DO 920 LINE=1,I-N-NR
                  P(N+NR+LINE,J)=P(MSTU(4)-MSTU(32)-LINE,J)
                  K(N+NR+LINE,J)=K(MSTU(4)-MSTU(32)-LINE,J)
  920          CONTINUE
            ENDIF
            IN(J)=INMO(J)
  930    CONTINUE
         GOTO 880
      ENDIF
      XTMO(JT)=XTMO3
C.. MOPS end of modification
 
      DO 940 J=1,3
        IN(J)=IN(3*JT+J)
  940 CONTINUE
 
C...Stepping within or from 'low' string region easy.
      IF(IN(1)+1.EQ.IN(2).AND.Z*P(IN(1)+2,3)*P(IN(2)+2,3)*
     &P(IN(1),5)**2.GE.PR(JT)) THEN
        P(IN(JT)+2,4)=Z*P(IN(JT)+2,3)
        P(IN(JR)+2,4)=PR(JT)/(P(IN(JT)+2,4)*P(IN(1),5)**2)
        DO 950 J=1,4
          P(I,J)=(PX(JT)+PX(3))*P(IN(3),J)+(PY(JT)+PY(3))*P(IN(3)+1,J)
  950   CONTINUE
        GOTO 1040
      ELSEIF(IN(1)+1.EQ.IN(2)) THEN
        P(IN(JR)+2,4)=P(IN(JR)+2,3)
        P(IN(JR)+2,JT)=1D0
        IN(JR)=IN(JR)+4*JS
        IF(JS*IN(JR).GT.JS*IN(4*JR)) GOTO 710
        IF(FOUR(IN(1),IN(2)).LE.1D-2) THEN
          P(IN(JT)+2,4)=P(IN(JT)+2,3)
          P(IN(JT)+2,JT)=0D0
          IN(JT)=IN(JT)+4*JS
        ENDIF
      ENDIF
 
C...Find new transverse directions (i.e. spacelike string vectors).
  960 IF(JS*IN(1).GT.JS*IN(3*JR+1).OR.JS*IN(2).GT.JS*IN(3*JR+2).OR.
     &IN(1).GT.IN(2)) GOTO 710
      IF(IN(1).NE.IN(3*JT+1).OR.IN(2).NE.IN(3*JT+2)) THEN
        DO 970 J=1,4
          DP(1,J)=P(IN(1),J)
          DP(2,J)=P(IN(2),J)
          DP(3,J)=0D0
          DP(4,J)=0D0
  970   CONTINUE
        DP(1,4)=SQRT(DP(1,1)**2+DP(1,2)**2+DP(1,3)**2)
        DP(2,4)=SQRT(DP(2,1)**2+DP(2,2)**2+DP(2,3)**2)
        DHC12=DFOUR(1,2)
        IF(DHC12.LE.1D-2) THEN
          P(IN(JT)+2,4)=P(IN(JT)+2,3)
          P(IN(JT)+2,JT)=0D0
          IN(JT)=IN(JT)+4*JS
          GOTO 960
        ENDIF
        IN(3)=N+NR+4*NS+5
        DP(5,1)=DP(1,1)/DP(1,4)-DP(2,1)/DP(2,4)
        DP(5,2)=DP(1,2)/DP(1,4)-DP(2,2)/DP(2,4)
        DP(5,3)=DP(1,3)/DP(1,4)-DP(2,3)/DP(2,4)
        IF(DP(5,1)**2.LE.DP(5,2)**2+DP(5,3)**2) DP(3,1)=1D0
        IF(DP(5,1)**2.GT.DP(5,2)**2+DP(5,3)**2) DP(3,3)=1D0
        IF(DP(5,2)**2.LE.DP(5,1)**2+DP(5,3)**2) DP(4,2)=1D0
        IF(DP(5,2)**2.GT.DP(5,1)**2+DP(5,3)**2) DP(4,3)=1D0
        DHCX1=DFOUR(3,1)/DHC12
        DHCX2=DFOUR(3,2)/DHC12
        DHCXX=1D0/SQRT(1D0+2D0*DHCX1*DHCX2*DHC12)
        DHCY1=DFOUR(4,1)/DHC12
        DHCY2=DFOUR(4,2)/DHC12
        DHCYX=DHCXX*(DHCX1*DHCY2+DHCX2*DHCY1)*DHC12
        DHCYY=1D0/SQRT(1D0+2D0*DHCY1*DHCY2*DHC12-DHCYX**2)
        DO 980 J=1,4
          DP(3,J)=DHCXX*(DP(3,J)-DHCX2*DP(1,J)-DHCX1*DP(2,J))
          P(IN(3),J)=DP(3,J)
          P(IN(3)+1,J)=DHCYY*(DP(4,J)-DHCY2*DP(1,J)-DHCY1*DP(2,J)-
     &    DHCYX*DP(3,J))
  980   CONTINUE
C...Express pT with respect to new axes, if sensible.
        PXP=-(PX(3)*FOUR(IN(3*JT+3),IN(3))+PY(3)*
     &  FOUR(IN(3*JT+3)+1,IN(3)))
        PYP=-(PX(3)*FOUR(IN(3*JT+3),IN(3)+1)+PY(3)*
     &  FOUR(IN(3*JT+3)+1,IN(3)+1))
        IF(ABS(PXP**2+PYP**2-PX(3)**2-PY(3)**2).LT.0.01D0) THEN
          PX(3)=PXP
          PY(3)=PYP
        ENDIF
      ENDIF
 
C...Sum up known four-momentum. Gives coefficients for m2 expression.
      DO 1010 J=1,4
        DHG(J)=0D0
        P(I,J)=PX(JT)*P(IN(3*JT+3),J)+PY(JT)*P(IN(3*JT+3)+1,J)+
     &  PX(3)*P(IN(3),J)+PY(3)*P(IN(3)+1,J)
        DO 990 IN1=IN(3*JT+1),IN(1)-4*JS,4*JS
          P(I,J)=P(I,J)+P(IN1+2,3)*P(IN1,J)
  990   CONTINUE
        DO 1000 IN2=IN(3*JT+2),IN(2)-4*JS,4*JS
          P(I,J)=P(I,J)+P(IN2+2,3)*P(IN2,J)
 1000   CONTINUE
 1010 CONTINUE
      DHM(1)=FOUR(I,I)
      DHM(2)=2D0*FOUR(I,IN(1))
      DHM(3)=2D0*FOUR(I,IN(2))
      DHM(4)=2D0*FOUR(IN(1),IN(2))
 
C...Find coefficients for Gamma expression.
      DO 1030 IN2=IN(1)+1,IN(2),4
        DO 1020 IN1=IN(1),IN2-1,4
          DHC=2D0*FOUR(IN1,IN2)
          DHG(1)=DHG(1)+P(IN1+2,JT)*P(IN2+2,JT)*DHC
          IF(IN1.EQ.IN(1)) DHG(2)=DHG(2)-JS*P(IN2+2,JT)*DHC
          IF(IN2.EQ.IN(2)) DHG(3)=DHG(3)+JS*P(IN1+2,JT)*DHC
          IF(IN1.EQ.IN(1).AND.IN2.EQ.IN(2)) DHG(4)=DHG(4)-DHC
 1020   CONTINUE
 1030 CONTINUE
 
C...Solve (m2, Gamma) equation system for energies taken.
      DHS1=DHM(JR+1)*DHG(4)-DHM(4)*DHG(JR+1)
      IF(ABS(DHS1).LT.1D-4) GOTO 710
      DHS2=DHM(4)*(GAM(3)-DHG(1))-DHM(JT+1)*DHG(JR+1)-DHG(4)*
     &(P(I,5)**2-DHM(1))+DHG(JT+1)*DHM(JR+1)
      DHS3=DHM(JT+1)*(GAM(3)-DHG(1))-DHG(JT+1)*(P(I,5)**2-DHM(1))
      P(IN(JR)+2,4)=0.5D0*(SQRT(MAX(0D0,DHS2**2-4D0*DHS1*DHS3))/
     &ABS(DHS1)-DHS2/DHS1)
      IF(DHM(JT+1)+DHM(4)*P(IN(JR)+2,4).LE.0D0) GOTO 710
      P(IN(JT)+2,4)=(P(I,5)**2-DHM(1)-DHM(JR+1)*P(IN(JR)+2,4))/
     &(DHM(JT+1)+DHM(4)*P(IN(JR)+2,4))
 
C...Step to new region if necessary.
      IF(P(IN(JR)+2,4).GT.P(IN(JR)+2,3)) THEN
        P(IN(JR)+2,4)=P(IN(JR)+2,3)
        P(IN(JR)+2,JT)=1D0
        IN(JR)=IN(JR)+4*JS
        IF(JS*IN(JR).GT.JS*IN(4*JR)) GOTO 710
        IF(FOUR(IN(1),IN(2)).LE.1D-2) THEN
          P(IN(JT)+2,4)=P(IN(JT)+2,3)
          P(IN(JT)+2,JT)=0D0
          IN(JT)=IN(JT)+4*JS
        ENDIF
        GOTO 960
      ELSEIF(P(IN(JT)+2,4).GT.P(IN(JT)+2,3)) THEN
        P(IN(JT)+2,4)=P(IN(JT)+2,3)
        P(IN(JT)+2,JT)=0D0
        IN(JT)=IN(JT)+4*JS
        GOTO 960
      ENDIF
 
C...Four-momentum of particle. Remaining quantities. Loop back.
 1040 DO 1050 J=1,4
        P(I,J)=P(I,J)+P(IN(1)+2,4)*P(IN(1),J)+P(IN(2)+2,4)*P(IN(2),J)
        P(N+NRS,J)=P(N+NRS,J)-P(I,J)
 1050 CONTINUE
      IF(P(IN(1)+2,4).GT.1D0+PARU(14).OR.P(IN(1)+2,4).LT.-PARU(14).OR.
     &P(IN(2)+2,4).GT.1D0+PARU(14).OR.P(IN(2)+2,4).LT.-PARU(14))
     &GOTO 200
      IF(P(I,4).LT.P(I,5)) GOTO 710
      KFL(JT)=-KFL(3)
      PMQ(JT)=PMQ(3)
      PX(JT)=-PX(3)
      PY(JT)=-PY(3)
      GAM(JT)=GAM(3)
      IF(IN(3).NE.IN(3*JT+3)) THEN
        DO 1060 J=1,4
          P(IN(3*JT+3),J)=P(IN(3),J)
          P(IN(3*JT+3)+1,J)=P(IN(3)+1,J)
 1060   CONTINUE
      ENDIF
      DO 1070 JQ=1,2
        IN(3*JT+JQ)=IN(JQ)
        P(IN(JQ)+2,3)=P(IN(JQ)+2,3)-P(IN(JQ)+2,4)
        P(IN(JQ)+2,JT)=P(IN(JQ)+2,JT)-JS*(3-2*JQ)*P(IN(JQ)+2,4)
 1070 CONTINUE
      IF(IBARRK(JT).EQ.1.AND.MOD(IABS(K(I,2)),10000).GT.1000)
     &IBARRK(JT)=0
      GOTO 870
 
C...Final hadron: side, flavour, hadron, mass.
 1080 I=I+1
      K(I,1)=1
      K(I,3)=IE(JR)
      K(I,4)=0
      K(I,5)=0
      CALL PYKFDI(KFL(JR),-KFL(3),KFLDMP,K(I,2))
      IF(K(I,2).EQ.0) GOTO 710
      IF(IBARRK(JT).EQ.1.AND.MOD(IABS(K(I-1,2)),10000).GT.1000)
     &IBARRK(JT)=0
      IF(IBARRK(JT).EQ.1.AND.MOD(IABS(K(I,2)),10000).GT.1000)
     &K(I,3)=IJUORI(JT)
      IF(IBARRK(JR).EQ.1.AND.MOD(IABS(K(I,2)),10000).GT.1000)
     &K(I,3)=IJUORI(JR)
      P(I,5)=PYMASS(K(I,2))
      PR(JR)=P(I,5)**2+(PX(JR)-PX(3))**2+(PY(JR)-PY(3))**2
 
C...Final two hadrons: find common setup of four-vectors.
      JQ=1
      IF(P(IN(4)+2,3)*P(IN(5)+2,3)*FOUR(IN(4),IN(5)).LT.
     &P(IN(7)+2,3)*P(IN(8)+2,3)*FOUR(IN(7),IN(8))) JQ=2
      DHC12=FOUR(IN(3*JQ+1),IN(3*JQ+2))
      DHR1=FOUR(N+NRS,IN(3*JQ+2))/DHC12
      DHR2=FOUR(N+NRS,IN(3*JQ+1))/DHC12
      IF(IN(4).NE.IN(7).OR.IN(5).NE.IN(8)) THEN
        PX(3-JQ)=-FOUR(N+NRS,IN(3*JQ+3))-PX(JQ)
        PY(3-JQ)=-FOUR(N+NRS,IN(3*JQ+3)+1)-PY(JQ)
        PR(3-JQ)=P(I+(JT+JQ-3)**2-1,5)**2+(PX(3-JQ)+(2*JQ-3)*JS*
     &  PX(3))**2+(PY(3-JQ)+(2*JQ-3)*JS*PY(3))**2
      ENDIF
 
C...Solve kinematics for final two hadrons, if possible.
      WREM2=2D0*DHR1*DHR2*DHC12
      FD=(SQRT(PR(1))+SQRT(PR(2)))/SQRT(WREM2)
      IF(MJU(1)+MJU(2).NE.0.AND.I.EQ.ISAV+2.AND.FD.GE.1D0) GOTO 200
      IF(FD.GE.1D0) GOTO 710
      FA=WREM2+PR(JT)-PR(JR)
      FB=SQRT(MAX(0D0,FA**2-4D0*WREM2*PR(JT)))
      PREVCF=PARJ(42)
      IF(MSTJ(11).EQ.2) PREVCF=PARJ(39)
      PREV=1D0/(1D0+EXP(MIN(50D0,PREVCF*FB*PARJ(40))))
      FB=SIGN(FB,JS*(PYR(0)-PREV))
      KFL1A=IABS(KFL(1))
      KFL2A=IABS(KFL(2))
      IF(MAX(MOD(KFL1A,10),MOD(KFL1A/1000,10),MOD(KFL2A,10),
     &MOD(KFL2A/1000,10)).GE.6) FB=SIGN(SQRT(MAX(0D0,FA**2-
     &4D0*WREM2*PR(JT))),DBLE(JS))
      DO 1090 J=1,4
        P(I-1,J)=(PX(JT)+PX(3))*P(IN(3*JQ+3),J)+(PY(JT)+PY(3))*
     &  P(IN(3*JQ+3)+1,J)+0.5D0*(DHR1*(FA+FB)*P(IN(3*JQ+1),J)+
     &  DHR2*(FA-FB)*P(IN(3*JQ+2),J))/WREM2
        P(I,J)=P(N+NRS,J)-P(I-1,J)
 1090 CONTINUE
      IF(P(I-1,4).LT.P(I-1,5).OR.P(I,4).LT.P(I,5)) GOTO 710
      DM2F1=P(I-1,4)**2-P(I-1,1)**2-P(I-1,2)**2-P(I-1,3)**2-P(I-1,5)**2
      DM2F2=P(I,4)**2-P(I,1)**2-P(I,2)**2-P(I,3)**2-P(I,5)**2
      IF(DM2F1.GT.1D-10*P(I-1,4)**2.OR.DM2F2.GT.1D-10*P(I,4)**2) THEN
        NTRYFN=NTRYFN+1
        IF(NTRYFN.LT.100) GOTO 140
        CALL PYERRM(13,'(PYSTRF:) bad energies for final two hadrons')
      ENDIF
 
C...Mark jets as fragmented and give daughter pointers.
      N=I-NRS+1
      DO 1100 I=NSAV+1,NSAV+NP
        IM=K(I,3)
        K(IM,1)=K(IM,1)+10
        IF(MSTU(16).NE.2) THEN
          K(IM,4)=NSAV+1
          K(IM,5)=NSAV+1
        ELSE
          K(IM,4)=NSAV+2
          K(IM,5)=N
        ENDIF
 1100 CONTINUE
 
C...Document string system. Move up particles.
      NSAV=NSAV+1
      K(NSAV,1)=11
      K(NSAV,2)=92
      K(NSAV,3)=IP
      K(NSAV,4)=NSAV+1
      K(NSAV,5)=N
      DO 1110 J=1,4
        P(NSAV,J)=DPS(J)
        V(NSAV,J)=V(IP,J)
 1110 CONTINUE
      P(NSAV,5)=SQRT(MAX(0D0,DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2))
      V(NSAV,5)=0D0
      DO 1130 I=NSAV+1,N
        DO 1120 J=1,5
          K(I,J)=K(I+NRS-1,J)
          P(I,J)=P(I+NRS-1,J)
          V(I,J)=0D0
 1120   CONTINUE
 1130 CONTINUE
      MSTU91=MSTU(90)
      DO 1140 IZ=MSTU90+1,MSTU91
        MSTU9T(IZ)=MSTU(90+IZ)-NRS+1-NSAV+N
        PARU9T(IZ)=PARU(90+IZ)
 1140 CONTINUE
      MSTU(90)=MSTU90
 
C...Order particles in rank along the chain. Update mother pointer.
      DO 1160 I=NSAV+1,N
        DO 1150 J=1,5
          K(I-NSAV+N,J)=K(I,J)
          P(I-NSAV+N,J)=P(I,J)
 1150   CONTINUE
 1160 CONTINUE
      I1=NSAV
      DO 1190 I=N+1,2*N-NSAV
        IF(K(I,3).NE.IE(1).AND.K(I,3).NE.IJUORI(1)) GOTO 1190
        I1=I1+1
        DO 1170 J=1,5
          K(I1,J)=K(I,J)
          P(I1,J)=P(I,J)
 1170   CONTINUE
        IF(MSTU(16).NE.2) K(I1,3)=NSAV
        DO 1180 IZ=MSTU90+1,MSTU91
          IF(MSTU9T(IZ).EQ.I) THEN
            MSTU(90)=MSTU(90)+1
            MSTU(90+MSTU(90))=I1
            PARU(90+MSTU(90))=PARU9T(IZ)
          ENDIF
 1180   CONTINUE
 1190 CONTINUE
      DO 1220 I=2*N-NSAV,N+1,-1
        IF(K(I,3).EQ.IE(1).OR.K(I,3).EQ.IJUORI(1)) GOTO 1220
        I1=I1+1
        DO 1200 J=1,5
          K(I1,J)=K(I,J)
          P(I1,J)=P(I,J)
 1200   CONTINUE
        IF(MSTU(16).NE.2) K(I1,3)=NSAV
        DO 1210 IZ=MSTU90+1,MSTU91
          IF(MSTU9T(IZ).EQ.I) THEN
            MSTU(90)=MSTU(90)+1
            MSTU(90+MSTU(90))=I1
            PARU(90+MSTU(90))=PARU9T(IZ)
          ENDIF
 1210   CONTINUE
 1220 CONTINUE
 
C...Boost back particle system. Set production vertices.
      IF(MBST.EQ.0) THEN
        MSTU(33)=1
        CALL PYROBO(NSAV+1,N,0D0,0D0,DPS(1)/DPS(4),DPS(2)/DPS(4),
     &  DPS(3)/DPS(4))
      ELSE
        DO 1230 I=NSAV+1,N
          HHPMT=P(I,1)**2+P(I,2)**2+P(I,5)**2
          IF(P(I,3).GT.0D0) THEN
            HHPEZ=(P(I,4)+P(I,3))*HHBZ
            P(I,3)=0.5D0*(HHPEZ-HHPMT/HHPEZ)
            P(I,4)=0.5D0*(HHPEZ+HHPMT/HHPEZ)
          ELSE
            HHPEZ=(P(I,4)-P(I,3))/HHBZ
            P(I,3)=-0.5D0*(HHPEZ-HHPMT/HHPEZ)
            P(I,4)=0.5D0*(HHPEZ+HHPMT/HHPEZ)
          ENDIF
 1230   CONTINUE
      ENDIF
      DO 1250 I=NSAV+1,N
        DO 1240 J=1,4
          V(I,J)=V(IP,J)
 1240   CONTINUE
 1250 CONTINUE
 
      RETURN
      END
