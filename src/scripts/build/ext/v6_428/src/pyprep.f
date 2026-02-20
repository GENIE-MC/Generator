 
C*********************************************************************
 
C...PYPREP
C...Rearranges partons along strings.
C...Special considerations for systems with junctions, with
C...possibility of junction-antijunction annihilation.
C...Allows small systems to collapse into one or two particles.
C...Checks flavours and colour singlet invariant masses.
 
      SUBROUTINE PYPREP(IP)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYINT1/MINT(400),VINT(400)
C...The common block of colour tags.
      COMMON/PYCTAG/NCT,MCT(4000,2)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYINT1/,/PYCTAG/,
     &/PYPARS/
      DATA NERRPR/0/
      SAVE NERRPR
C...Local arrays.
      DIMENSION DPS(5),DPC(5),UE(3),PG(5),E1(3),E2(3),E3(3),E4(3),
     &ECL(3),IJUNC(10,0:4),IPIECE(30,0:4),KFEND(4),KFQ(4),
     &IJUR(4),PJU(4,6),IRNG(4,2),TJJ(2,5),T(5),PUL(3,5),
     &IJCP(0:6),TJUOLD(5)
      CHARACTER CHTMP*6
 
C...Function to give four-product.
      FOUR(I,J)=P(I,4)*P(J,4)-P(I,1)*P(J,1)-P(I,2)*P(J,2)-P(I,3)*P(J,3)
 
C...Rearrange parton shower product listing along strings: begin loop.
      MSTU(24)=0
      NOLD=N
      I1=N
      NJUNC=0
      NPIECE=0
      NJJSTR=0
      MSTU32=MSTU(32)+1
      DO 100 I=MAX(1,IP),N
C...First store junction positions.
        IF(K(I,1).EQ.42) THEN
          NJUNC=NJUNC+1
          IJUNC(NJUNC,0)=I
          IJUNC(NJUNC,4)=0
        ENDIF
  100 CONTINUE
 
      DO 250 MQGST=1,3
        DO 240 I=MAX(1,IP),N
C...Special treatment for junctions
          IF (K(I,1).LE.0) GOTO 240
          IF(K(I,1).EQ.42) THEN
C...MQGST=2: Look for junction-junction strings (not detected in the
C...main search below).
            IF (MQGST.EQ.2.AND.NPIECE.NE.3*NJUNC) THEN
              IF (NJJSTR.EQ.0) THEN
                NJJSTR = (3*NJUNC-NPIECE)/2
              ENDIF
C...Check how many already identified strings end on this junction
              ILC=0
              DO 110 J=1,NPIECE
                IF (IPIECE(J,4).EQ.I) ILC=ILC+1
  110         CONTINUE
C...If less than 3, remaining must be to another junction
              IF (ILC.LT.3) THEN
                IF (ILC.NE.2) THEN
C...Multiple j-j connections not handled yet.
                  CALL PYERRM(2,
     &            '(PYPREP:) Too many junction-junction strings.')
                  MINT(51)=1
                  RETURN
                ENDIF
C...The colour information in the junction is unreadable for the
C...colour space search further down in this routine, so we must
C...start on the colour mother of this junction and then "artificially"
C...prevent the colour mother from connecting here again.
                ITJUNC=MOD(K(I,4)/MSTU(5),MSTU(5))
                KCS=4
                IF (MOD(ITJUNC,2).EQ.0) KCS=5
C...Switch colour if the junction-junction leg is presumably a
C...junction mother leg rather than a junction daughter leg.
                IF (ITJUNC.GE.3) KCS=9-KCS
                IF (MINT(33).EQ.0) THEN
C...Find the unconnected leg and reorder junction daughter pointers so
C...MOD(K(I,4),MSTU(5)) always points to the junction-junction string
C...piece.
                  IA=MOD(K(I,4),MSTU(5))
                  IF (K(IA,KCS)/MSTU(5)**2.GE.2) THEN
                    ITMP=MOD(K(I,5),MSTU(5))
                    IF (K(ITMP,KCS)/MSTU(5)**2.GE.2) THEN
                      ITMP=MOD(K(I,5)/MSTU(5),MSTU(5))
                      K(I,5)=K(I,5)+(IA-ITMP)*MSTU(5)
                    ELSE
                      K(I,5)=K(I,5)+(IA-ITMP)
                    ENDIF
                    K(I,4)=K(I,4)+(ITMP-IA)
                    IA=ITMP
                  ENDIF
                  IF (ITJUNC.LE.2) THEN
C...Beam baryon junction
                    K(IA,KCS)   = K(IA,KCS) + 2*MSTU(5)**2
                    K(I,KCS)    = K(I,KCS) + 1*MSTU(5)**2
C...Else 1 -> 2 decay junction
                  ELSE
                    K(IA,KCS)   = K(IA,KCS) + MSTU(5)**2
                    K(I,KCS)    = K(I,KCS) + 2*MSTU(5)**2
                  ENDIF
                  I1BEG = I1
                  NSTP = 0
                  GOTO 170
C...Alternatively use colour tag information.
                ELSE
C...Find a final state parton with appropriate dangling colour tag.
                  JCT=0
                  IA=0
                  IJUMO=K(I,3)
                  DO 140 J1=MAX(1,IP),N
                    IF (K(J1,1).NE.3) GOTO 140
C...Check for matching final-state colour tag
                    IMATCH=0
                    DO 120 J2=MAX(1,IP),N
                      IF (K(J2,1).NE.3) GOTO 120
                      IF (MCT(J1,KCS-3).EQ.MCT(J2,6-KCS)) IMATCH=1
  120               CONTINUE
                    IF (IMATCH.EQ.1) GOTO 140
C...Check whether this colour tag belongs to the present junction
C...by seeing whether any parton with this colour tag has the same
C...mother as the junction.
                    JCT=MCT(J1,KCS-3)
                    IMATCH=0
                    DO 130 J2=MINT(84)+1,N
                      IMO2=K(J2,3)
C...First scattering partons have IMO1 = 3 and 4.
                      IF (IMO2.EQ.MINT(83)+3.OR.IMO2.EQ.MINT(83)+4)
     &                     IMO2=IMO2-2
                      IF (MCT(J2,KCS-3).EQ.JCT.AND.IMO2.EQ.IJUMO)
     &                     IMATCH=1
  130               CONTINUE
                    IF (IMATCH.EQ.0) GOTO 140
                    IA=J1
  140             CONTINUE
C...Check for junction-junction strings without intermediate final state
C...glue (not detected above).
                  IF (IA.EQ.0) THEN
                    DO 160 MJU=1,NJUNC
                      IJU2=IJUNC(MJU,0)
                      IF (IJU2.EQ.I) GOTO 160
                      ITJU2=MOD(K(IJU2,4)/MSTU(5),MSTU(5))
C...Only opposite types of junctions can connect to each other.
                      IF (MOD(ITJU2,2).EQ.MOD(ITJUNC,2)) GOTO 160
                      IS=0
                      DO 150 J=1,NPIECE
                        IF (IPIECE(J,4).EQ.IJU2) IS=IS+1
  150                 CONTINUE
                      IF (IS.EQ.3) GOTO 160
                      IB=I
                      IA=IJU2
  160               CONTINUE
                  ENDIF
C...Switch to other side of adjacent parton and step from there.
                  KCS=9-KCS
                  I1BEG = I1
                  NSTP = 0
                  GOTO 170
                ENDIF
              ELSE IF (ILC.NE.3) THEN
              ENDIF
            ENDIF
          ENDIF
 
C...Look for coloured string endpoint, or (later) leftover gluon.
          IF(K(I,1).NE.3) GOTO 240
          KC=PYCOMP(K(I,2))
          IF(KC.EQ.0) GOTO 240
          KQ=KCHG(KC,2)
          IF(KQ.EQ.0.OR.(MQGST.LE.2.AND.KQ.EQ.2)) GOTO 240
 
C...Pick up loose string end.
          KCS=4
          IF(KQ*ISIGN(1,K(I,2)).LT.0) KCS=5
          IA=I
          IB=I
          I1BEG=I1
          NSTP=0
  170     NSTP=NSTP+1
          IF(NSTP.GT.4*N) THEN
            CALL PYERRM(14,'(PYPREP:) caught in infinite loop')
            MINT(51)=1
            RETURN
          ENDIF
 
C...Copy undecayed parton. Finished if reached string endpoint.
          IF(K(IA,1).EQ.3) THEN
            IF(I1.GE.MSTU(4)-MSTU32-5) THEN
              CALL PYERRM(11,'(PYPREP:) no more memory left in PYJETS')
              MINT(51)=1
              MSTU(24)=1
              RETURN
            ENDIF
            I1=I1+1
            K(I1,1)=2
            IF(NSTP.GE.2.AND.KCHG(PYCOMP(K(IA,2)),2).NE.2) K(I1,1)=1
            K(I1,2)=K(IA,2)
            K(I1,3)=IA
            K(I1,4)=0
            K(I1,5)=0
            DO 180 J=1,5
              P(I1,J)=P(IA,J)
              V(I1,J)=V(IA,J)
  180       CONTINUE
            K(IA,1)=K(IA,1)+10
            IF(K(I1,1).EQ.1) GOTO 240
          ENDIF
 
C...Also finished (for now) if reached junction; then copy to end.
          IF(K(IA,1).EQ.42) THEN
            NCOPY=I1-I1BEG
            IF(I1.GE.MSTU(4)-MSTU32-NCOPY-5) THEN
              CALL PYERRM(11,'(PYPREP:) no more memory left in PYJETS')
              MINT(51)=1
              MSTU(24)=1
              RETURN
            ENDIF
            IF (MQGST.LE.2.AND.NCOPY.NE.0) THEN
              DO 200 ICOPY=1,NCOPY
                DO 190 J=1,5
                  K(MSTU(4)-MSTU32-ICOPY,J)=K(I1BEG+ICOPY,J)
                  P(MSTU(4)-MSTU32-ICOPY,J)=P(I1BEG+ICOPY,J)
                  V(MSTU(4)-MSTU32-ICOPY,J)=V(I1BEG+ICOPY,J)
  190           CONTINUE
  200         CONTINUE
            ENDIF
C...For junction-junction strings, find end leg and reorder junction
C...daughter pointers so MOD(K(I,4),MSTU(5)) always points to the
C...junction-junction string piece.
            IF (K(I,1).EQ.42.AND.MINT(33).EQ.0) THEN
              ITMP=MOD(K(IA,4),MSTU(5))
              IF (ITMP.NE.IB) THEN
                IF (MOD(K(IA,5),MSTU(5)).EQ.IB) THEN
                  K(IA,5)=K(IA,5)+(ITMP-IB)
                ELSE
                  K(IA,5)=K(IA,5)+(ITMP-IB)*MSTU(5)
                ENDIF
                K(IA,4)=K(IA,4)+(IB-ITMP)
              ENDIF
            ENDIF
            NPIECE=NPIECE+1
C...IPIECE:
C...0: endpoint in original ER
C...1:
C...2:
C...3: Parton immediately next to junction
C...4: Junction
            IPIECE(NPIECE,0)=I
            IPIECE(NPIECE,1)=MSTU32+1
            IPIECE(NPIECE,2)=MSTU32+NCOPY
            IPIECE(NPIECE,3)=IB
            IPIECE(NPIECE,4)=IA
            MSTU32=MSTU32+NCOPY
            I1=I1BEG
            GOTO 240
          ENDIF
 
C...GOTO next parton in colour space.
          IB=IA
          IF (MINT(33).EQ.0) THEN
            IF(MOD(K(IB,KCS)/MSTU(5)**2,2).EQ.0.AND.MOD(K(IB,KCS),MSTU(5
     &           )).NE.0) THEN
              IA=MOD(K(IB,KCS),MSTU(5))
              K(IB,KCS)=K(IB,KCS)+MSTU(5)**2
              MREV=0
            ELSE
              IF(K(IB,KCS).GE.2*MSTU(5)**2.OR.MOD(K(IB,KCS)/MSTU(5),
     &             MSTU(5)).EQ.0) KCS=9-KCS
              IA=MOD(K(IB,KCS)/MSTU(5),MSTU(5))
              K(IB,KCS)=K(IB,KCS)+2*MSTU(5)**2
              MREV=1
            ENDIF
            IF(IA.LE.0.OR.IA.GT.N) THEN
              CALL PYERRM(12,'(PYPREP:) colour rearrangement failed')
              IF(NERRPR.LT.5) THEN
                NERRPR=NERRPR+1
                WRITE(MSTU(11),*) 'started at:', I
                WRITE(MSTU(11),*) 'ended going from',IB,' to',IA
                WRITE(MSTU(11),*) 'MQGST =',MQGST
                CALL PYLIST(4)
              ENDIF
              MINT(51)=1
              RETURN
            ENDIF
            IF(MOD(K(IA,4)/MSTU(5),MSTU(5)).EQ.IB.OR.MOD(K(IA,5)/MSTU(5)
     &           ,MSTU(5)).EQ.IB) THEN
              IF(MREV.EQ.1) KCS=9-KCS
              IF(MOD(K(IA,KCS)/MSTU(5),MSTU(5)).NE.IB) KCS=9-KCS
              K(IA,KCS)=K(IA,KCS)+2*MSTU(5)**2
            ELSE
              IF(MREV.EQ.0) KCS=9-KCS
              IF(MOD(K(IA,KCS),MSTU(5)).NE.IB) KCS=9-KCS
              K(IA,KCS)=K(IA,KCS)+MSTU(5)**2
            ENDIF
            IF(IA.NE.I) GOTO 170
C...Use colour tag information
          ELSE
C...First create colour tags starting on IB if none already present.
            IF (MCT(IB,KCS-3).EQ.0) THEN
              CALL PYCTTR(IB,KCS,IB)
              IF(MINT(51).NE.0) RETURN
            ENDIF
            JCT=MCT(IB,KCS-3)
            IFOUND=0
C...Find final state tag partner
            DO 210 IT=MAX(1,IP),N
              IF (IT.EQ.IB) GOTO 210
              IF (MCT(IT,6-KCS).EQ.JCT.AND.K(IT,1).LT.10.AND.K(IT,1).GT
     &             .0) THEN
                IFOUND=IFOUND+1
                IA=IT
              ENDIF
  210       CONTINUE
C...Just copy and goto next if exactly one partner found.
            IF (IFOUND.EQ.1) THEN
              GOTO 170
C...When no match found, match is presumably junction.
            ELSEIF (IFOUND.EQ.0.AND.MQGST.LE.2) THEN
C...Check whether this colour tag matches a junction
C...by seeing whether any parton with this colour tag has the same
C...mother as a junction.
C...NB: Only type 1 and 2 junctions handled presently.
              DO 230 IJU=1,NJUNC
                IJUMO=K(IJUNC(IJU,0),3)
                ITJUNC=MOD(K(IJUNC(IJU,0),4)/MSTU(5),MSTU(5))
C...Colours only connect to junctions, anti-colours to antijunctions:
                IF (MOD(ITJUNC+1,2)+1.NE.KCS-3) GOTO 230
                IMATCH=0
                DO 220 J1=MAX(1,IP),N
                  IF (K(J1,1).LE.0) GOTO 220
C...First scattering partons have IMO1 = 3 and 4.
                  IMO=K(J1,3)
                  IF (IMO.EQ.MINT(83)+3.OR.IMO.EQ.MINT(83)+4)
     &                 IMO=IMO-2
                  IF (MCT(J1,KCS-3).EQ.JCT.AND.IMO.EQ.IJUMO.AND.MOD(K(J1
     &                 ,3+ITJUNC)/MSTU(5),MSTU(5)).EQ.IJUNC(IJU,0))
     &                 IMATCH=1
C...Attempt at handling type > 3 junctions also. Not tested.
                  IF (ITJUNC.GE.3.AND.MCT(J1,6-KCS).EQ.JCT.AND.IMO.EQ
     &                 .IJUMO) IMATCH=1
  220           CONTINUE
                IF (IMATCH.EQ.0) GOTO 230
                IA=IJUNC(IJU,0)
                IFOUND=IFOUND+1
  230         CONTINUE
 
              IF (IFOUND.EQ.1) THEN
                GOTO 170
              ELSEIF (IFOUND.EQ.0) THEN
                WRITE(CHTMP,'(I6)') JCT
                CALL PYERRM(12,'(PYPREP:) no matching colour tag: '
     &               //CHTMP)
                IF(NERRPR.LT.5) THEN
                  NERRPR=NERRPR+1
                  CALL PYLIST(4)
                ENDIF
                MINT(51)=1
                RETURN
              ENDIF
            ELSEIF (IFOUND.GE.2) THEN
              WRITE(CHTMP,'(I6)') JCT
              CALL PYERRM(12
     &             ,'(PYPREP:) too many occurences of colour line: '//
     &             CHTMP)
              IF(NERRPR.LT.5) THEN
                NERRPR=NERRPR+1
                CALL PYLIST(4)
              ENDIF
              MINT(51)=1
              RETURN
            ENDIF
          ENDIF
          K(I1,1)=1
  240   CONTINUE
  250 CONTINUE
 
C...Junction systems remain.
      IJU=0
      IJUS=0
      IJUCNT=0
      MREV=0
      IJJSTR=0
  260 IJUCNT=IJUCNT+1
      IF (IJUCNT.LE.NJUNC) THEN
C...If we are not processing a j-j string, treat this junction as new.
        IF (IJJSTR.EQ.0) THEN
          IJU=IJUNC(IJUCNT,0)
          MREV=0
C...If junction has already been read, ignore it.
          IF (IJUNC(IJUCNT,4).EQ.1) GOTO 260
C...If we are on a j-j string, goto second j-j junction.
        ELSE
          IJUCNT=IJUCNT-1
          IJU=IJUS
        ENDIF
C...Mark selected junction read.
        DO 270 J=1,NJUNC
          IF (IJUNC(J,0).EQ.IJU) IJUNC(J,4)=1
  270   CONTINUE
C...Determine junction type
        ITJUNC = MOD(K(IJU,4)/MSTU(5),MSTU(5))
C...Type 1 and 2 junctions: ~chi -> q q q, ~chi -> qbar,qbar,qbar
C...Type 3 and 4 junctions: ~qbar -> q q , ~q -> qbar qbar
C...Type 5 and 6 junctions: ~g -> q q q, ~g -> qbar qbar qbar
        IF (ITJUNC.GE.1.AND.ITJUNC.LE.6) THEN
          IHK=0
  280     IHK=IHK+1
C...Find which quarks belong to given junction.
          IHF=0
          DO 290 IPC=1,NPIECE
            IF (IPIECE(IPC,4).EQ.IJU) THEN
              IHF=IHF+1
              IF (IHF.EQ.IHK) IEND=IPIECE(IPC,3)
            ENDIF
            IF (IHK.EQ.3.AND.IPIECE(IPC,0).EQ.IJU) IEND=IPIECE(IPC,3)
  290     CONTINUE
C...IHK = 3 is special. Either normal string piece, or j-j string.
          IF(IHK.EQ.3) THEN
            IF (MREV.NE.1) THEN
              DO 300 IPC=1,NPIECE
C...If there is a j-j string starting on the present junction which has
C...zero length, insert next junction immediately.
                IF (IPIECE(IPC,0).EQ.IJU.AND.K(IPIECE(IPC,4),1)
     &          .EQ.42.AND.IPIECE(IPC,1)-1-IPIECE(IPC,2).EQ.0) THEN
                  IJJSTR = 1
                  GOTO 340
                ENDIF
  300         CONTINUE
              MREV = 1
C...If MREV is 1 and IHK is 3 we are finished with this system.
            ELSE
              MREV=0
              GOTO 260
            ENDIF
          ENDIF
 
C...If we've gotten this far, then either IHK < 3, or
C...an interjunction string exists, or just a third normal string.
          IJUNC(IJUCNT,IHK)=0
          IJJSTR = 0
C..Order pieces belonging to this junction. Also look for j-j.
          DO 310 IPC=1,NPIECE
            IF (IPIECE(IPC,3).EQ.IEND) IJUNC(IJUCNT,IHK)=IPC
            IF (IHK.EQ.3.AND.IPIECE(IPC,0).EQ.IJUNC(IJUCNT,0)
     &      .AND.K(IPIECE(IPC,4),1).EQ.42) THEN
              IJUNC(IJUCNT,IHK)=IPC
              IJJSTR = 1
              MREV = 0
            ENDIF
  310     CONTINUE
C...Copy back chains in proper order. MREV=0/1 : descending/ascending
          IPC=IJUNC(IJUCNT,IHK)
C...Temporary solution to cover for bug.
          IF(IPC.LE.0) THEN
            CALL PYERRM(12,'(PYPREP:) fails to hook up junctions')
            MINT(51)=1
            RETURN
          ENDIF
          DO 330 ICP=IPIECE(IPC,1+MREV),IPIECE(IPC,2-MREV),1-2*MREV
            I1=I1+1
            DO 320 J=1,5
              K(I1,J)=K(MSTU(4)-ICP,J)
              P(I1,J)=P(MSTU(4)-ICP,J)
              V(I1,J)=V(MSTU(4)-ICP,J)
  320       CONTINUE
  330     CONTINUE
          K(I1,1)=2
C...Mark last quark.
          IF (MREV.EQ.1.AND.IHK.GE.2) K(I1,1)=1
C...Do not insert junctions at wrong places.
          IF(IHK.LT.2.OR.MREV.NE.0) GOTO 360
C...Insert junction.
  340     IJUS = IJU
          IF (IHK.EQ.3) THEN
C...Shift to end junction if a j-j string has been processed.
            IF (IJJSTR.NE.0) IJUS = IPIECE(IPC,4)
            MREV= 1
          ENDIF
          I1=I1+1
          DO 350 J=1,5
            K(I1,J)=0
            P(I1,J)=0.
            V(I1,J)=0.
  350     CONTINUE
          K(I1,1)=41
          K(IJUS,1)=K(IJUS,1)+10
          K(I1,2)=K(IJUS,2)
          K(I1,3)=IJUS
  360     IF (IHK.LT.3) GOTO 280
        ELSE
          CALL PYERRM(12,'(PYPREP:) Unknown junction type')
          MINT(51)=1
          RETURN
        ENDIF
        IF (IJUCNT.NE.NJUNC) GOTO 260
      ENDIF
      N=I1
 
C...Rearrange three strings from junction, e.g. in case one has been
C...shortened by shower, so the last is the largest-energy one.
      IF(NJUNC.GE.1) THEN
C...Find systems with exactly one junction.
        MJUN1=0
        NBEG=NOLD+1
        DO 470 I=NOLD+1,N
          IF(K(I,1).NE.1.AND.K(I,1).NE.41) THEN
          ELSEIF(K(I,1).EQ.41) THEN
            MJUN1=MJUN1+1
          ELSEIF(K(I,1).EQ.1.AND.MJUN1.NE.1) THEN
            MJUN1=0
            NBEG=I+1
          ELSE
            NEND=I
C...Sum up energy-momentum in each junction string.
            DO 370 J=1,5
              PJU(1,J)=0D0
              PJU(2,J)=0D0
              PJU(3,J)=0D0
  370       CONTINUE
            NJU=0
            DO 390 I1=NBEG,NEND
              IF(K(I1,2).NE.21) THEN
                NJU=NJU+1
                IJUR(NJU)=I1
              ENDIF
              DO 380 J=1,5
                PJU(MIN(NJU,3),J)=PJU(MIN(NJU,3),J)+P(I1,J)
  380         CONTINUE
  390       CONTINUE
C...Find which of them has highest energy (minus mass) in rest frame.
            DO 400 J=1,5
              PJU(4,J)=PJU(1,J)+PJU(2,J)+PJU(3,J)
  400       CONTINUE
            PMJU=SQRT(MAX(0D0,PJU(4,4)**2-PJU(4,1)**2-PJU(4,2)**2-
     &      PJU(4,3)**2))
            DO 410 I2=1,3
              PJU(I2,6)=(PJU(4,4)*PJU(I2,4)-PJU(4,1)*PJU(I2,1)-
     &        PJU(4,2)*PJU(I2,2)-PJU(4,3)*PJU(I2,3))/PMJU-PJU(I2,5)
  410       CONTINUE
            IF(PJU(3,6).LT.MIN(PJU(1,6),PJU(2,6))) THEN
C...Decide how to rearrange so that new last has highest energy.
              IF(PJU(1,6).LT.PJU(2,6)) THEN
                IRNG(1,1)=IJUR(1)
                IRNG(1,2)=IJUR(2)-1
                IRNG(2,1)=IJUR(4)
                IRNG(2,2)=IJUR(3)+1
                IRNG(4,1)=IJUR(3)-1
                IRNG(4,2)=IJUR(2)
              ELSE
                IRNG(1,1)=IJUR(4)
                IRNG(1,2)=IJUR(3)+1
                IRNG(2,1)=IJUR(2)
                IRNG(2,2)=IJUR(3)-1
                IRNG(4,1)=IJUR(2)-1
                IRNG(4,2)=IJUR(1)
              ENDIF
              IRNG(3,1)=IJUR(3)
              IRNG(3,2)=IJUR(3)
C...Copy in correct order below bottom of current event record.
              I2=N
              DO 440 II=1,4
                DO 430 I1=IRNG(II,1),IRNG(II,2),
     &          ISIGN(1,IRNG(II,2)-IRNG(II,1))
                  I2=I2+1
                  IF(I2.GE.MSTU(4)-MSTU32-5) THEN
                    CALL PYERRM(11,
     &              '(PYPREP:) no more memory left in PYJETS')
                    MINT(51)=1
                    MSTU(24)=1
                    RETURN
                  ENDIF
                  DO 420 J=1,5
                    K(I2,J)=K(I1,J)
                    P(I2,J)=P(I1,J)
                    V(I2,J)=V(I1,J)
  420             CONTINUE
                  IF(K(I2,1).EQ.1) K(I2,1)=2
  430           CONTINUE
  440         CONTINUE
              K(I2,1)=1
C...Copy back up, overwriting but now in correct order.
              DO 460 I1=NBEG,NEND
                I2=I1-NBEG+N+1
                DO 450 J=1,5
                  K(I1,J)=K(I2,J)
                  P(I1,J)=P(I2,J)
                  V(I1,J)=V(I2,J)
  450           CONTINUE
  460         CONTINUE
            ENDIF
            MJUN1=0
            NBEG=I+1
          ENDIF
  470   CONTINUE
 
C...Check whether q-q-j-j-qbar-qbar systems should be collapsed
C...to two q-qbar systems.
C...(MSTJ(19)=1 forces q-q-j-j-qbar-qbar.)
        IF (MSTJ(19).NE.1) THEN
          MJUN1  = 0
          JJGLUE = 0
          NBEG   = NOLD+1
C...Force collapse when MSTJ(19)=2.
          IF (MSTJ(19).EQ.2) THEN
            DELMJJ = 1D9
            DELMQQ = 0D0
          ENDIF
C...Find systems with exactly two junctions.
          DO 700 I=NOLD+1,N
C...Count junctions
            IF (K(I,1).EQ.41) THEN
              MJUN1 = MJUN1+1
C...Check for interjunction gluons
              IF (MJUN1.EQ.2.AND.K(I-1,1).NE.41) THEN
                JJGLUE = 1
              ENDIF
            ELSEIF(K(I,1).EQ.1.AND.(MJUN1.NE.2)) THEN
C...If end of system reached with either zero or one junction, restart
C...with next system.
              MJUN1  = 0
              JJGLUE = 0
              NBEG   = I+1
            ELSEIF(K(I,1).EQ.1) THEN
C...If end of system reached with exactly two junctions, compute string
C...length measure for the (q-q-j-j-qbar-qbar) topology and compare with
C...length measure for the (q-qbar)(q-qbar) topology.
              NEND=I
C...Loop down through chain.
              ISID=0
              DO 480 I1=NBEG,NEND
C...Store string piece division locations in event record
                IF (K(I1,2).NE.21) THEN
                  ISID       = ISID+1
                  IJCP(ISID) = I1
                ENDIF
  480         CONTINUE
C...Randomly choose between (1,3)(2,4) and (1,4)(2,3) topologies.
              ISW=0
              IF (PYR(0).LT.0.5D0) ISW=1
C...Randomly choose which qqbar string gets the jj gluons.
              IGS=1
              IF (PYR(0).GT.0.5D0) IGS=2
C...Only compute string lengths when no topology forced.
              IF (MSTJ(19).EQ.0) THEN
C...Repeat following for each junction
                DO 570 IJU=1,2
C...Initialize iterative procedure for finding JRF
                  IJRFIT=0
                  DO 490 IX=1,3
                    TJUOLD(IX)=0D0
  490             CONTINUE
                  TJUOLD(4)=1D0
C...Start iteration. Sum up momenta in string pieces
  500             DO 540 IJS=1,3
C...JD=-1 for first junction, +1 for second junction.
C...Find out where piece starts and ends and which direction to go.
                    JD=2*IJU-3
                    IF (IJS.LE.2) THEN
                      IA = IJCP((IJU-1)*7 - JD*(IJS+1)) + JD
                      IB = IJCP((IJU-1)*7 - JD*IJS)
                    ELSEIF (IJS.EQ.3) THEN
                      JD =-JD
                      IA = IJCP((IJU-1)*7 + JD*(IJS)) + JD
                      IB = IJCP((IJU-1)*7 + JD*(IJS+3))
                    ENDIF
C...Initialize junction pull 4-vector.
                    DO 510 J=1,5
                      PUL(IJS,J)=0D0
  510               CONTINUE
C...Initialize weight
                    PWT = 0D0
                    PWTOLD = 0D0
C...Sum up (weighted) momenta along each string piece
                    DO 530 ISP=IA,IB,JD
C...If present parton not last in chain
                      IF (ISP.NE.IA.AND.ISP.NE.IB) THEN
C...If last parton was a junction, store present weight
                        IF (K(ISP-JD,2).EQ.88) THEN
                          PWTOLD = PWT
C...If last parton was a quark, reset to stored weight.
                        ELSEIF (K(ISP-JD,2).NE.21) THEN
                          PWT = PWTOLD
                        ENDIF
                      ENDIF
C...Skip next parton if weight already large
                      IF (PWT.GT.10D0) GOTO 530
C...Compute momentum in TJUOLD frame:
                      TDP=TJUOLD(1)*P(ISP,1)+TJUOLD(2)*P(ISP,2)+TJUOLD(3
     &                     )*P(ISP,3)
                      BFC=TDP/(1D0+TJUOLD(4))+P(ISP,4)
                      DO 520 J=1,3
                        TMP=P(ISP,J)+TJUOLD(J)*BFC
                        PUL(IJS,J)=PUL(IJS,J)+TMP*EXP(-PWT)
  520                 CONTINUE
C...Boosted energy
                      TMP=TJUOLD(4)*P(ISP,4)+TDP
                      PUL(IJS,4)=PUL(IJS,J)+TMP*EXP(-PWT)
C...Update weight
                      PWT=PWT+TMP/PARJ(48)
C...Put |p| rather than m in 5th slot
                      PUL(IJS,5)=SQRT(PUL(IJS,1)**2+PUL(IJS,2)**2
     &                     +PUL(IJS,3)**2)
  530               CONTINUE
  540             CONTINUE
C...Compute boost
                  IJRFIT=IJRFIT+1
                  CALL PYJURF(PUL,T)
C...Combine new boost (T) with old boost (TJUOLD)
                  TMP=T(1)*TJUOLD(1)+T(2)*TJUOLD(2)+T(3)*TJUOLD(3)
                  DO 550 IX=1,3
                    TJUOLD(IX)=T(IX)+TJUOLD(IX)*(TMP/(1D0+TJUOLD(4))+T(4
     &                   ))
  550             CONTINUE
                  TJUOLD(4)=SQRT(1D0+TJUOLD(1)**2+TJUOLD(2)**2+TJUOLD(3)
     &                 **2)
C...If last boost small, accept JRF, else iterate.
C...Also prevent possibility of infinite loop.
                  IF (ABS((T(4)-1D0)/TJUOLD(4)).GT.0.01D0.AND.
     &                 IJRFIT.LT.MSTJ(18))THEN
                    GOTO 500
                  ELSEIF (IJRFIT.GE.MSTJ(18)) THEN
                    CALL PYERRM(1,'(PYPREP:) failed to converge on JRF')
                  ENDIF
C...Store final boost, with change of sign since TJJ motion vector.
                  DO 560 IX=1,3
                    TJJ(IJU,IX)=-TJUOLD(IX)
  560             CONTINUE
                  TJJ(IJU,4)=SQRT(1D0+TJJ(IJU,1)**2+TJJ(IJU,2)**2
     &                 +TJJ(IJU,3)**2)
  570           CONTINUE
C...String length measure for (q-qbar)(q-qbar) topology.
C...Note only momenta of nearest partons used (since rest of system
C...identical).
                IF (JJGLUE.EQ.0) THEN
                  DELMQQ=4D0*FOUR(IJCP(2)-1,IJCP(4+ISW)+1)*FOUR(IJCP(3)
     &                 -1,IJCP(5-ISW)+1)
                ELSE
C...Put jj gluons on selected string (IGS selected randomly above).
                  IF (IGS.EQ.1) THEN
                    DELMQQ=8D0*FOUR(IJCP(2)-1,IJCP(4)-1)*FOUR(IJCP(3)+1
     &                   ,IJCP(4+ISW)+1)*FOUR(IJCP(3)-1,IJCP(5-ISW)+1)
                  ELSE
                    DELMQQ=8D0*FOUR(IJCP(2)-1,IJCP(4+ISW)+1)
     &                   *FOUR(IJCP(3)-1,IJCP(4)-1)*FOUR(IJCP(3)+1
     &                   ,IJCP(5-ISW)+1)
                  ENDIF
                ENDIF
C...String length measure for q-q-j-j-q-q topology.
                T1G1=0D0
                T2G2=0D0
                T1T2=0D0
                T1P1=0D0
                T1P2=0D0
                T2P3=0D0
                T2P4=0D0
                ISGN=-1
C...Note only momenta of nearest partons used (since rest of system
C...identical).
                DO 580 IX=1,4
                  IF (IX.EQ.4) ISGN=1
                  T1P1=T1P1+ISGN*TJJ(1,IX)*P(IJCP(2)-1,IX)
                  T1P2=T1P2+ISGN*TJJ(1,IX)*P(IJCP(3)-1,IX)
                  T2P3=T2P3+ISGN*TJJ(2,IX)*P(IJCP(4)+1,IX)
                  T2P4=T2P4+ISGN*TJJ(2,IX)*P(IJCP(5)+1,IX)
                  IF (JJGLUE.EQ.0) THEN
C...Junction motion vector dot product gives length when inter-junction
C...gluons absent.
                    T1T2=T1T2+ISGN*TJJ(1,IX)*TJJ(2,IX)
                  ELSE
C...Junction motion vector dot products with gluon momenta give length
C...when inter-junction gluons present.
                    T1G1=T1G1+ISGN*TJJ(1,IX)*P(IJCP(3)+1,IX)
                    T2G2=T2G2+ISGN*TJJ(2,IX)*P(IJCP(4)-1,IX)
                  ENDIF
  580           CONTINUE
                DELMJJ=16D0*T1P1*T1P2*T2P3*T2P4
                IF (JJGLUE.EQ.0) THEN
                  DELMJJ=DELMJJ*(T1T2+SQRT(T1T2**2-1))
                ELSE
                  DELMJJ=DELMJJ*4D0*T1G1*T2G2
                ENDIF
              ENDIF
C...If delmjj > delmqq collapse string system to q-qbar q-qbar
C...(Always the case for MSTJ(19)=2 due to initialization above)
              IF (DELMJJ.GT.DELMQQ) THEN
C...Put new system at end of event record
                NCOP=N
                DO 650 IST=1,2
                  DO 600 ICOP=IJCP(IST),IJCP(IST+1)-1
                    NCOP=NCOP+1
                    DO 590 IX=1,5
                      P(NCOP,IX)=P(ICOP,IX)
                      K(NCOP,IX)=K(ICOP,IX)
  590               CONTINUE
  600             CONTINUE
                  IF (JJGLUE.NE.0.AND.IST.EQ.IGS) THEN
C...Insert inter-junction gluon string piece (reversed)
                    NJJGL=0
                    DO 620 ICOP=IJCP(4)-1,IJCP(3)+1,-1
                      NJJGL=NJJGL+1
                      NCOP=NCOP+1
                      DO 610 IX=1,5
                        P(NCOP,IX)=P(ICOP,IX)
                        K(NCOP,IX)=K(ICOP,IX)
  610                 CONTINUE
  620               CONTINUE
                    ENDIF
                  IFC=-2*IST+3
                  DO 640 ICOP=IJCP(IST+IFC*ISW+3)+1,IJCP(IST+IFC*ISW+4)
                    NCOP=NCOP+1
                    DO 630 IX=1,5
                      P(NCOP,IX)=P(ICOP,IX)
                      K(NCOP,IX)=K(ICOP,IX)
  630               CONTINUE
  640             CONTINUE
                  K(NCOP,1)=1
  650           CONTINUE
C...Copy system back in right order
                DO 670 ICOP=NBEG,NEND-2
                  DO 660 IX=1,5
                    P(ICOP,IX)=P(N+ICOP-NBEG+1,IX)
                    K(ICOP,IX)=K(N+ICOP-NBEG+1,IX)
  660             CONTINUE
  670           CONTINUE
C...Shift down rest of event record
                DO 690 ICOP=NEND+1,N
                  DO 680 IX=1,5
                    P(ICOP-2,IX)=P(ICOP,IX)
                    K(ICOP-2,IX)=K(ICOP,IX)
  680             CONTINUE
  690             CONTINUE
C...Update length of event record.
                N=N-2
              ENDIF
              MJUN1=0
              NBEG=I+1
            ENDIF
  700     CONTINUE
        ENDIF
      ENDIF
 
C...Done if no checks on small-mass systems.
      IF(MSTJ(14).LT.0) RETURN
      IF(MSTJ(14).EQ.0) GOTO 1140
 
C...Find lowest-mass colour singlet jet system.
      NS=N
  710 NSIN=N-NS
      PDMIN=1D0+PARJ(32)
      IC=0
      DO 770 I=MAX(1,IP),N
        IF(K(I,1).NE.1.AND.K(I,1).NE.2) THEN
        ELSEIF(K(I,1).EQ.2.AND.IC.EQ.0) THEN
          NSIN=NSIN+1
          IC=I
          DO 720 J=1,4
            DPS(J)=P(I,J)
  720     CONTINUE
          MSTJ(93)=1
          DPS(5)=PYMASS(K(I,2))
        ELSEIF(K(I,1).EQ.2.AND.K(I,2).NE.21) THEN
          DO 730 J=1,4
            DPS(J)=DPS(J)+P(I,J)
  730     CONTINUE
          MSTJ(93)=1
          DPS(5)=DPS(5)+PYMASS(K(I,2))
        ELSEIF(K(I,1).EQ.2) THEN
          DO 740 J=1,4
            DPS(J)=DPS(J)+P(I,J)
  740     CONTINUE
        ELSEIF(IC.NE.0.AND.KCHG(PYCOMP(K(I,2)),2).NE.0) THEN
          DO 750 J=1,4
            DPS(J)=DPS(J)+P(I,J)
  750     CONTINUE
          MSTJ(93)=1
          DPS(5)=DPS(5)+PYMASS(K(I,2))
          PD=SQRT(MAX(0D0,DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2))-
     &    DPS(5)
          IF(PD.LT.PDMIN) THEN
            PDMIN=PD
            DO 760 J=1,5
              DPC(J)=DPS(J)
  760       CONTINUE
            IC1=IC
            IC2=I
          ENDIF
          IC=0
        ELSE
          NSIN=NSIN+1
        ENDIF
  770 CONTINUE
 
C...Done if lowest-mass system above threshold for string frag.
      IF(PDMIN.GE.PARJ(32)) GOTO 1140
 
C...Fill small-mass system as cluster.
      NSAV=N
      PECM=SQRT(MAX(0D0,DPC(4)**2-DPC(1)**2-DPC(2)**2-DPC(3)**2))
      K(N+1,1)=11
      K(N+1,2)=91
      K(N+1,3)=IC1
      P(N+1,1)=DPC(1)
      P(N+1,2)=DPC(2)
      P(N+1,3)=DPC(3)
      P(N+1,4)=DPC(4)
      P(N+1,5)=PECM
 
C...Set up history, assuming cluster -> 2 hadrons.
      NBODY=2
      K(N+1,4)=N+2
      K(N+1,5)=N+3
      K(N+2,1)=1
      K(N+3,1)=1
      IF(MSTU(16).NE.2) THEN
        K(N+2,3)=N+1
        K(N+3,3)=N+1
      ELSE
        K(N+2,3)=IC1
        K(N+3,3)=IC2
      ENDIF
      K(N+2,4)=0
      K(N+3,4)=0
      K(N+2,5)=0
      K(N+3,5)=0
      V(N+1,5)=0D0
      V(N+2,5)=0D0
      V(N+3,5)=0D0
 
C...Find total flavour content - complicated by presence of junctions.
      NQ=0
      NDIQ=0
      DO 780 I=IC1,IC2
        IF((K(I,1).EQ.1.OR.K(I,1).EQ.2).AND.K(I,2).NE.21) THEN
          NQ=NQ+1
          KFQ(NQ)=K(I,2)
          IF(IABS(K(I,2)).GT.1000) NDIQ=NDIQ+1
        ENDIF
  780 CONTINUE
 
C...If several diquarks, split up one to give even number of flavours.
      IF(NQ.EQ.3.AND.NDIQ.GE.2) THEN
        I1=3
        IF(IABS(KFQ(3)).LT.1000) I1=1
        KFQ(4)=ISIGN(MOD(IABS(KFQ(I1))/100,10),KFQ(I1))
        KFQ(I1)=KFQ(I1)/1000
        NQ=4
        NDIQ=NDIQ-1
      ENDIF
 
C...If four quark ends, join two to diquark.
      IF(NQ.EQ.4.AND.NDIQ.EQ.0) THEN
        I1=1
        I2=2
        IF(KFQ(I1)*KFQ(I2).LT.0) I2=3
        IF(I2.EQ.3.AND.KFQ(I1)*KFQ(I2).LT.0) I2=4
        KFLS=2*INT(PYR(0)+3D0*PARJ(4)/(1D0+3D0*PARJ(4)))+1
        IF(KFQ(I1).EQ.KFQ(I2)) KFLS=3
        KFQ(I1)=ISIGN(1000*MAX(IABS(KFQ(I1)),IABS(KFQ(I2)))+
     &  100*MIN(IABS(KFQ(I1)),IABS(KFQ(I2)))+KFLS,KFQ(I1))
        KFQ(I2)=KFQ(4)
        NQ=3
        NDIQ=1
      ENDIF
 
C...If two quark ends, plus quark or diquark, join quarks to diquark.
      IF(NQ.EQ.3) THEN
        I1=1
        I2=2
        IF(IABS(KFQ(I1)).GT.1000) I1=3
        IF(IABS(KFQ(I2)).GT.1000) I2=3
        KFLS=2*INT(PYR(0)+3D0*PARJ(4)/(1D0+3D0*PARJ(4)))+1
        IF(KFQ(I1).EQ.KFQ(I2)) KFLS=3
        KFQ(I1)=ISIGN(1000*MAX(IABS(KFQ(I1)),IABS(KFQ(I2)))+
     &  100*MIN(IABS(KFQ(I1)),IABS(KFQ(I2)))+KFLS,KFQ(I1))
        KFQ(I2)=KFQ(3)
        NQ=2
        NDIQ=NDIQ+1
      ENDIF
 
C...Form two particles from flavours of lowest-mass system, if feasible.
      NTRY = 0
  790 NTRY = NTRY + 1
 
C...Open string with two specified endpoint flavours.
      IF(NQ.EQ.2) THEN
        KC1=PYCOMP(KFQ(1))
        KC2=PYCOMP(KFQ(2))
        IF(KC1.EQ.0.OR.KC2.EQ.0) GOTO 1140
        KQ1=KCHG(KC1,2)*ISIGN(1,KFQ(1))
        KQ2=KCHG(KC2,2)*ISIGN(1,KFQ(2))
        IF(KQ1+KQ2.NE.0) GOTO 1140
C...Start with qq, if there is one. Only allow for rank 1 popcorn meson
  800   K1=KFQ(1)
        IF(IABS(KFQ(2)).GT.1000) K1=KFQ(2)
        MSTU(125)=0
        CALL PYDCYK(K1,0,KFLN,K(N+2,2))
        CALL PYDCYK(KFQ(1)+KFQ(2)-K1,-KFLN,KFLDMP,K(N+3,2))
        IF(K(N+2,2).EQ.0.OR.K(N+3,2).EQ.0) GOTO 800
 
C...Open string with four specified flavours.
      ELSEIF(NQ.EQ.4) THEN
        KC1=PYCOMP(KFQ(1))
        KC2=PYCOMP(KFQ(2))
        KC3=PYCOMP(KFQ(3))
        KC4=PYCOMP(KFQ(4))
        IF(KC1.EQ.0.OR.KC2.EQ.0.OR.KC3.EQ.0.OR.KC4.EQ.0) GOTO 1140
        KQ1=KCHG(KC1,2)*ISIGN(1,KFQ(1))
        KQ2=KCHG(KC2,2)*ISIGN(1,KFQ(2))
        KQ3=KCHG(KC3,2)*ISIGN(1,KFQ(3))
        KQ4=KCHG(KC4,2)*ISIGN(1,KFQ(4))
        IF(KQ1+KQ2+KQ3+KQ4.NE.0) GOTO 1140
C...Combine flavours pairwise to form two hadrons.
  810   I1=1
        I2=2
        IF(KQ1*KQ2.GT.0.OR.(IABS(KFQ(1)).GT.1000.AND.
     &  IABS(KFQ(2)).GT.1000)) I2=3
        IF(I2.EQ.3.AND.(KQ1*KQ3.GT.0.OR.(IABS(KFQ(1)).GT.1000.AND.
     &  IABS(KFQ(3)).GT.1000))) I2=4
        I3=3
        IF(I2.EQ.3) I3=2
        I4=10-I1-I2-I3
        CALL PYDCYK(KFQ(I1),KFQ(I2),KFLDMP,K(N+2,2))
        CALL PYDCYK(KFQ(I3),KFQ(I4),KFLDMP,K(N+3,2))
        IF(K(N+2,2).EQ.0.OR.K(N+3,2).EQ.0) GOTO 810
 
C...Closed string.
      ELSE
        IF(IABS(K(IC2,2)).NE.21) GOTO 1140
C...No room for popcorn mesons in closed string -> 2 hadrons.
        MSTU(125)=0
  820   CALL PYDCYK(1+INT((2D0+PARJ(2))*PYR(0)),0,KFLN,KFDMP)
        CALL PYDCYK(KFLN,0,KFLM,K(N+2,2))
        CALL PYDCYK(-KFLN,-KFLM,KFLDMP,K(N+3,2))
        IF(K(N+2,2).EQ.0.OR.K(N+3,2).EQ.0) GOTO 820
      ENDIF
      P(N+2,5)=PYMASS(K(N+2,2))
      P(N+3,5)=PYMASS(K(N+3,2))
 
C...If it does not work: try again (a number of times), give up (if no
C...place to shuffle momentum or too many flavours), or form one hadron.
      IF(P(N+2,5)+P(N+3,5)+PARJ(64).GE.PECM) THEN
        IF(NTRY.LT.MSTJ(17).OR.(NQ.EQ.4.AND.NTRY.LT.5*MSTJ(17))) THEN
          GOTO 790
        ELSEIF(NSIN.EQ.1.OR.NQ.EQ.4) THEN
          GOTO 1140
        ELSE
          GOTO 890
        END IF
      END IF
 
C...Perform two-particle decay of jet system.
C...First step: find reference axis in decaying system rest frame.
C...(Borrow slot N+2 for temporary direction.)
      DO 830 J=1,4
        P(N+2,J)=P(IC1,J)
  830 CONTINUE
      DO 850 I=IC1+1,IC2-1
        IF((K(I,1).EQ.1.OR.K(I,1).EQ.2).AND.
     &  KCHG(PYCOMP(K(I,2)),2).NE.0) THEN
          FRAC1=FOUR(IC2,I)/(FOUR(IC1,I)+FOUR(IC2,I))
          DO 840 J=1,4
            P(N+2,J)=P(N+2,J)+FRAC1*P(I,J)
  840     CONTINUE
        ENDIF
  850 CONTINUE
      CALL PYROBO(N+2,N+2,0D0,0D0,-DPC(1)/DPC(4),-DPC(2)/DPC(4),
     &-DPC(3)/DPC(4))
      THE1=PYANGL(P(N+2,3),SQRT(P(N+2,1)**2+P(N+2,2)**2))
      PHI1=PYANGL(P(N+2,1),P(N+2,2))
 
C...Second step: generate isotropic/anisotropic decay.
      PA=SQRT((PECM**2-(P(N+2,5)+P(N+3,5))**2)*(PECM**2-
     &(P(N+2,5)-P(N+3,5))**2))/(2D0*PECM)
  860 UE(3)=PYR(0)
      IF(PARJ(21).LE.0.01D0) UE(3)=1D0
      PT2=(1D0-UE(3)**2)*PA**2
      IF(MSTJ(16).LE.0) THEN
        PREV=0.5D0
      ELSE
        IF(EXP(-PT2/(2D0*MAX(0.01D0,PARJ(21))**2)).LT.PYR(0)) GOTO 860
        PR1=P(N+2,5)**2+PT2
        PR2=P(N+3,5)**2+PT2
        ALAMBD=SQRT(MAX(0D0,(PECM**2-PR1-PR2)**2-4D0*PR1*PR2))
        PREVCF=PARJ(42)
        IF(MSTJ(11).EQ.2) PREVCF=PARJ(39)
        PREV=1D0/(1D0+EXP(MIN(50D0,PREVCF*ALAMBD*PARJ(40))))
      ENDIF
      IF(PYR(0).LT.PREV) UE(3)=-UE(3)
      PHI=PARU(2)*PYR(0)
      UE(1)=SQRT(1D0-UE(3)**2)*COS(PHI)
      UE(2)=SQRT(1D0-UE(3)**2)*SIN(PHI)
      DO 870 J=1,3
        P(N+2,J)=PA*UE(J)
        P(N+3,J)=-PA*UE(J)
  870 CONTINUE
      P(N+2,4)=SQRT(PA**2+P(N+2,5)**2)
      P(N+3,4)=SQRT(PA**2+P(N+3,5)**2)
 
C...Third step: move back to event frame and set production vertex.
      CALL PYROBO(N+2,N+3,THE1,PHI1,DPC(1)/DPC(4),DPC(2)/DPC(4),
     &DPC(3)/DPC(4))
      DO 880 J=1,4
        V(N+1,J)=V(IC1,J)
        V(N+2,J)=V(IC1,J)
        V(N+3,J)=V(IC2,J)
  880 CONTINUE
      N=N+3
      GOTO 1120
 
C...Else form one particle, if possible.
  890 NBODY=1
      K(N+1,5)=N+2
      DO 900 J=1,4
        V(N+1,J)=V(IC1,J)
        V(N+2,J)=V(IC1,J)
  900 CONTINUE
 
C...Select hadron flavour from available quark flavours.
  910 IF(NQ.EQ.2.AND.IABS(KFQ(1)).GT.100.AND.IABS(KFQ(2)).GT.100) THEN
        GOTO 1140
      ELSEIF(NQ.EQ.2) THEN
        CALL PYKFDI(KFQ(1),KFQ(2),KFLDMP,K(N+2,2))
      ELSE
        KFLN=1+INT((2D0+PARJ(2))*PYR(0))
        CALL PYKFDI(KFLN,-KFLN,KFLDMP,K(N+2,2))
      ENDIF
      IF(K(N+2,2).EQ.0) GOTO 910
      P(N+2,5)=PYMASS(K(N+2,2))
 
C...Use old algorithm for E/p conservation? (EN)
      IF (MSTJ(16).LE.0) GOTO 1080
 
C...Find the string piece closest to the cluster by a loop
C...over the undecayed partons not in present cluster. (EN)
      DGLOMI=1D30
      IBEG=0
      I0=0
      NJUNC=0
      DO 940 I1=MAX(1,IP),N-1
        IF(K(I1,1).EQ.1) NJUNC=0
        IF(K(I1,1).EQ.41) NJUNC=NJUNC+1
        IF(K(I1,1).EQ.41) GOTO 940
        IF(I1.GE.IC1-1.AND.I1.LE.IC2) THEN
          I0=0
        ELSEIF(K(I1,1).EQ.2) THEN
          IF(I0.EQ.0) I0=I1
          I2=I1
  920     I2=I2+1
          IF(K(I2,1).EQ.41) GOTO 940
          IF(K(I2,1).GT.10) GOTO 920
          IF(KCHG(PYCOMP(K(I2,2)),2).EQ.0) GOTO 920
          IF(K(I1,2).EQ.21.AND.K(I2,2).NE.21.AND.K(I2,1).NE.1.AND.
     &    NJUNC.EQ.0) GOTO 940
          IF(K(I1,2).NE.21.AND.K(I2,2).EQ.21.AND.NJUNC.NE.0) GOTO 940
          IF(K(I1,2).NE.21.AND.K(I2,2).NE.21.AND.(I1.GT.I0.OR.
     &    K(I2,1).NE.1)) GOTO 940
 
C...Define velocity vectors e1, e2, ecl and differences e3, e4.
          DO 930 J=1,3
            E1(J)=P(I1,J)/P(I1,4)
            E2(J)=P(I2,J)/P(I2,4)
            ECL(J)=P(N+1,J)/P(N+1,4)
            E3(J)=E2(J)-E1(J)
            E4(J)=ECL(J)-E1(J)
  930     CONTINUE
 
C...Calculate minimal D=(e4-alpha*e3)**2 for 0<alpha<1.
          E3S=E3(1)**2+E3(2)**2+E3(3)**2
          E4S=E4(1)**2+E4(2)**2+E4(3)**2
          E34=E3(1)*E4(1)+E3(2)*E4(2)+E3(3)*E4(3)
          IF(E34.LE.0D0) THEN
            DDMIN=E4S
          ELSEIF(E34.LT.E3S) THEN
            DDMIN=E4S-E34**2/E3S
          ELSE
            DDMIN=E4S-2D0*E34+E3S
          ENDIF
 
C...Is this the smallest so far?
          IF(DDMIN.LT.DGLOMI) THEN
            DGLOMI=DDMIN
            IBEG=I0
            IPCS=I1
          ENDIF
        ELSEIF(K(I1,1).EQ.1.AND.KCHG(PYCOMP(K(I1,2)),2).NE.0) THEN
          I0=0
        ENDIF
  940 CONTINUE
 
C... Check if there are any strings to connect to the new gluon. (EN)
      IF (IBEG.EQ.0) GOTO 1080
 
C...Delta_m = m_clus - m_had > 0: emit a 'gluon' (EN)
      IF (P(N+1,5).GE.P(N+2,5)) THEN
 
C...Construct 'gluon' that is needed to put hadron on the mass shell.
        FRAC=P(N+2,5)/P(N+1,5)
        DO 950 J=1,5
          P(N+2,J)=FRAC*P(N+1,J)
          PG(J)=(1D0-FRAC)*P(N+1,J)
  950   CONTINUE
 
C... Copy string with new gluon put in.
        N=N+2
        I=IBEG-1
  960   I=I+1
        IF(K(I,1).NE.1.AND.K(I,1).NE.2.AND.K(I,1).NE.41) GOTO 960
        IF(KCHG(PYCOMP(K(I,2)),2).EQ.0.AND.K(I,1).NE.41) GOTO 960
        N=N+1
        DO 970 J=1,5
          K(N,J)=K(I,J)
          P(N,J)=P(I,J)
          V(N,J)=V(I,J)
  970   CONTINUE
        K(I,1)=K(I,1)+10
        K(I,4)=N
        K(I,5)=N
        K(N,3)=I
        IF(I.EQ.IPCS) THEN
          N=N+1
          DO 980 J=1,5
            K(N,J)=K(N-1,J)
            P(N,J)=PG(J)
            V(N,J)=V(N-1,J)
  980     CONTINUE
          K(N,2)=21
          K(N,3)=NSAV+1
        ENDIF
        IF(K(I,1).EQ.12.OR.K(I,1).EQ.51) GOTO 960
        GOTO 1120
 
C...Delta_m = m_clus - m_had < 0: have to absorb a 'gluon' instead,
C...from string piece endpoints.
      ELSE
 
C...Begin by copying string that should give energy to cluster.
        N=N+2
        I=IBEG-1
  990   I=I+1
        IF(K(I,1).NE.1.AND.K(I,1).NE.2.AND.K(I,1).NE.41) GOTO 990
        IF(KCHG(PYCOMP(K(I,2)),2).EQ.0.AND.K(I,1).NE.41) GOTO 990
        N=N+1
        DO 1000 J=1,5
          K(N,J)=K(I,J)
          P(N,J)=P(I,J)
          V(N,J)=V(I,J)
 1000   CONTINUE
        K(I,1)=K(I,1)+10
        K(I,4)=N
        K(I,5)=N
        K(N,3)=I
        IF(I.EQ.IPCS) I1=N
        IF(K(I,1).EQ.12.OR.K(I,1).EQ.51) GOTO 990
        I2=I1+1
 
C...Set initial Phad.
        DO 1010 J=1,4
          P(NSAV+2,J)=P(NSAV+1,J)
 1010   CONTINUE
 
C...Calculate Pg, a part of which will be added to Phad later. (EN)
 1020   IF(MSTJ(16).EQ.1) THEN
          ALPHA=1D0
          BETA=1D0
        ELSE
          ALPHA=FOUR(NSAV+1,I2)/FOUR(I1,I2)
          BETA=FOUR(NSAV+1,I1)/FOUR(I1,I2)
        ENDIF
        DO 1030 J=1,4
          PG(J)=ALPHA*P(I1,J)+BETA*P(I2,J)
 1030   CONTINUE
        PG(5)=SQRT(MAX(1D-20,PG(4)**2-PG(1)**2-PG(2)**2-PG(3)**2))
 
C..Solve 2nd order equation, use the best (smallest) solution. (EN)
        PMSCOL=P(NSAV+2,4)**2-P(NSAV+2,1)**2-P(NSAV+2,2)**2-
     &  P(NSAV+2,3)**2
        PCLPG=(P(NSAV+2,4)*PG(4)-P(NSAV+2,1)*PG(1)-
     &  P(NSAV+2,2)*PG(2)-P(NSAV+2,3)*PG(3))/PG(5)**2
        DELTA=SQRT(PCLPG**2+(P(NSAV+2,5)**2-PMSCOL)/PG(5)**2)-PCLPG
 
C...If all gluon energy eaten, zero it and take a step back.
        ITER=0
        IF(DELTA*ALPHA.GT.1D0.AND.I1.GT.NSAV+3.AND.K(I1,2).EQ.21) THEN
          ITER=1
          DO 1040 J=1,4
            P(NSAV+2,J)=P(NSAV+2,J)+P(I1,J)
            P(I1,J)=0D0
 1040     CONTINUE
          P(I1,5)=0D0
          K(I1,1)=K(I1,1)+10
          I1=I1-1
          IF(K(I1,1).EQ.41) ITER=-1
        ENDIF
        IF(DELTA*BETA.GT.1D0.AND.I2.LT.N.AND.K(I2,2).EQ.21) THEN
          ITER=1
          DO 1050 J=1,4
            P(NSAV+2,J)=P(NSAV+2,J)+P(I2,J)
            P(I2,J)=0D0
 1050     CONTINUE
          P(I2,5)=0D0
          K(I2,1)=K(I2,1)+10
          I2=I2+1
          IF(K(I2,1).EQ.41) ITER=-1
        ENDIF
        IF(ITER.EQ.1) GOTO 1020
 
C...If also all endpoint energy eaten, revert to old procedure.
        IF((1D0-DELTA*ALPHA)*P(I1,4).LT.P(I1,5).OR.
     &  (1D0-DELTA*BETA)*P(I2,4).LT.P(I2,5).OR.ITER.EQ.-1) THEN
          DO 1060 I=NSAV+3,N
            IM=K(I,3)
            K(IM,1)=K(IM,1)-10
            K(IM,4)=0
            K(IM,5)=0
 1060     CONTINUE
          N=NSAV
          GOTO 1080
        ENDIF
 
C... Construct the collapsed hadron and modified string partons.
        DO 1070 J=1,4
          P(NSAV+2,J)=P(NSAV+2,J)+DELTA*PG(J)
          P(I1,J)=(1D0-DELTA*ALPHA)*P(I1,J)
          P(I2,J)=(1D0-DELTA*BETA)*P(I2,J)
 1070   CONTINUE
          P(I1,5)=(1D0-DELTA*ALPHA)*P(I1,5)
          P(I2,5)=(1D0-DELTA*BETA)*P(I2,5)
 
C...Finished with string collapse in new scheme.
        GOTO 1120
      ENDIF
 
C... Use old algorithm; by choice or when in trouble.
 1080 CONTINUE
C...Find parton/particle which combines to largest extra mass.
      IR=0
      HA=0D0
      HSM=0D0
      DO 1100 MCOMB=1,3
        IF(IR.NE.0) GOTO 1100
        DO 1090 I=MAX(1,IP),N
          IF(K(I,1).LE.0.OR.K(I,1).GT.10.OR.(I.GE.IC1.AND.I.LE.IC2
     &    .AND.K(I,1).GE.1.AND.K(I,1).LE.2)) GOTO 1090
          IF(MCOMB.EQ.1) KCI=PYCOMP(K(I,2))
          IF(MCOMB.EQ.1.AND.KCI.EQ.0) GOTO 1090
          IF(MCOMB.EQ.1.AND.KCHG(KCI,2).EQ.0.AND.I.LE.NS) GOTO 1090
          IF(MCOMB.EQ.2.AND.IABS(K(I,2)).GT.10.AND.IABS(K(I,2)).LE.100)
     &    GOTO 1090
          HCR=DPC(4)*P(I,4)-DPC(1)*P(I,1)-DPC(2)*P(I,2)-DPC(3)*P(I,3)
          HSR=2D0*HCR+PECM**2-P(N+2,5)**2-2D0*P(N+2,5)*P(I,5)
          IF(HSR.GT.HSM) THEN
            IR=I
            HA=HCR
            HSM=HSR
          ENDIF
 1090   CONTINUE
 1100 CONTINUE
 
C...Shuffle energy and momentum to put new particle on mass shell.
      IF(IR.NE.0) THEN
        HB=PECM**2+HA
        HC=P(N+2,5)**2+HA
        HD=P(IR,5)**2+HA
        HK2=0.5D0*(HB*SQRT(MAX(0D0,((HB+HC)**2-4D0*(HB+HD)*P(N+2,5)**2)/
     &  (HA**2-(PECM*P(IR,5))**2)))-(HB+HC))/(HB+HD)
        HK1=(0.5D0*(P(N+2,5)**2-PECM**2)+HD*HK2)/HB
        DO 1110 J=1,4
          P(N+2,J)=(1D0+HK1)*DPC(J)-HK2*P(IR,J)
          P(IR,J)=(1D0+HK2)*P(IR,J)-HK1*DPC(J)
 1110   CONTINUE
        N=N+2
      ELSE
        CALL PYERRM(3,'(PYPREP:) no match for collapsing cluster')
        RETURN
      ENDIF
 
C...Mark collapsed system and store daughter pointers. Iterate.
 1120 DO 1130 I=IC1,IC2
        IF((K(I,1).EQ.1.OR.K(I,1).EQ.2).AND.
     &  KCHG(PYCOMP(K(I,2)),2).NE.0) THEN
          K(I,1)=K(I,1)+10
          IF(MSTU(16).NE.2) THEN
            K(I,4)=NSAV+1
            K(I,5)=NSAV+1
          ELSE
            K(I,4)=NSAV+2
            K(I,5)=NSAV+1+NBODY
          ENDIF
        ENDIF
        IF(K(I,1).EQ.41) K(I,1)=K(I,1)+10
 1130 CONTINUE
      IF(N.LT.MSTU(4)-MSTU(32)-5) GOTO 710
 
C...Check flavours and invariant masses in parton systems.
 1140 NP=0
      KFN=0
      KQS=0
      NJU=0
      DO 1150 J=1,5
        DPS(J)=0D0
 1150 CONTINUE
      DO 1180 I=MAX(1,IP),N
        IF(K(I,1).EQ.41) NJU=NJU+1
        IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 1180
        KC=PYCOMP(K(I,2))
        IF(KC.EQ.0) GOTO 1180
        KQ=KCHG(KC,2)*ISIGN(1,K(I,2))
        IF(KQ.EQ.0) GOTO 1180
        NP=NP+1
        IF(KQ.NE.2) THEN
          KFN=KFN+1
          KQS=KQS+KQ
          MSTJ(93)=1
          DPS(5)=DPS(5)+PYMASS(K(I,2))
        ENDIF
        DO 1160 J=1,4
          DPS(J)=DPS(J)+P(I,J)
 1160   CONTINUE
        IF(K(I,1).EQ.1) THEN
          NFERR=0
          IF(NJU.EQ.0.AND.NP.NE.1) THEN
            IF(KFN.EQ.1.OR.KFN.GE.3.OR.KQS.NE.0) NFERR=1
          ELSEIF(NJU.EQ.1) THEN
            IF(KFN.NE.3.OR.IABS(KQS).NE.3) NFERR=1
          ELSEIF(NJU.EQ.2) THEN
            IF(KFN.NE.4.OR.KQS.NE.0) NFERR=1
          ELSEIF(NJU.GE.3) THEN
            NFERR=1
          ENDIF
          IF(NFERR.EQ.1) THEN
            CALL PYERRM(2,'(PYPREP:) unphysical flavour combination')
            MINT(51)=1
            RETURN
          ENDIF
          IF(NP.NE.1.AND.DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2.LT.
     &    (0.9D0*PARJ(32)+DPS(5))**2) CALL PYERRM(3,
     &    '(PYPREP:) too small mass in jet system')
          NP=0
          KFN=0
          KQS=0
          NJU=0
          DO 1170 J=1,5
            DPS(J)=0D0
 1170     CONTINUE
        ENDIF
 1180 CONTINUE
 
      RETURN
      END
