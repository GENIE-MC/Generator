 
C***************************************************************
 
C...PYKFIN
C...Precalculates a set of diquark and popcorn weights.
 
      SUBROUTINE PYKFIN
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYDAT1/,/PYDAT2/
 
      DIMENSION SU6(12),SU6M(7),QBB(7),QBM(7),DMB(14)
 
 
      MSTU(123)=1
C..Diquark indices for dimensional variables
      IUD1=1
      IUU1=2
      IUS0=3
      ISU0=4
      IUS1=5
      ISU1=6
      ISS1=7
 
C.. *** SU(6) factors **
C..Modify with decuplet- (and Sigma/Lambda-) suppression.
      PARF(146)=1D0
      IF(MSTJ(12).GE.5) PARF(146)=3D0*PARJ(18)/(2D0*PARJ(18)+1D0)
      IF(PARJ(18).LT.1D0-1D-5.AND.MSTJ(12).LT.5) CALL PYERRM(9,
     &     '(PYKFIN:) PARJ(18)<1 combined with 0<MSTJ(12)<5 option')
      DO 100 I=1,6
         SU6(I)=PARF(60+I)
         SU6(6+I)=SU6(I)*4*PARF(146)/(3*PARF(146)+1)
  100 CONTINUE
      SU6(8)=SU6(2)*4/(3*PARF(146)+1)
      SU6(6)=SU6(6)*(3+PARF(146))/(3*PARF(146)+1)
      DO 110 I=1,6
         SU6(I)=SU6(I)+PARJ(18)*PARF(70+I)
         SU6(6+I)=SU6(6+I)+PARJ(18)*PARF(70+I)
  110 CONTINUE
 
C..SU(6)max            q       q'     s,c,b
      SU6MUD    =MAX(SU6(1) ,       SU6(8) )
      SU6M(IUD1)=MAX(SU6(5) ,       SU6(12))
      SU6M(ISU0)=MAX(SU6(7) ,SU6(2),SU6MUD )
      SU6M(IUU1)=MAX(SU6(3) ,SU6(4),SU6(10))
      SU6M(ISU1)=MAX(SU6(11),SU6(6),SU6M(IUD1))
      SU6M(IUS0)=SU6M(ISU0)
      SU6M(ISS1)=SU6M(IUU1)
      SU6M(IUS1)=SU6M(ISU1)
 
C..Store SU(6)max, in order UD0,UD1,US0,US1,QQ1
      PARF(141)=SU6MUD
      PARF(142)=SU6M(IUD1)
      PARF(143)=SU6M(ISU0)
      PARF(144)=SU6M(ISU1)
      PARF(145)=SU6M(ISS1)
 
C..diquark SU(6) survival =
C..sum over quark (quark tunnel weight)*(SU(6)).
      PUD0=(2D0*SU6(1)+PARJ(2)*SU6(8))
      DMB(ISU0)=(SU6(7)+SU6(2)+PARJ(2)*SU6(1))/PUD0
      DMB(IUS0)=DMB(ISU0)
      DMB(ISS1)=(2D0*SU6(4)+PARJ(2)*SU6(3))/PUD0
      DMB(IUU1)=(SU6(3)+SU6(4)+PARJ(2)*SU6(10))/PUD0
      DMB(ISU1)=(SU6(11)+SU6(6)+PARJ(2)*SU6(5))/PUD0
      DMB(IUS1)=DMB(ISU1)
      DMB(IUD1)=(2D0*SU6(5)+PARJ(2)*SU6(12))/PUD0
 
C.. *** Tunneling factors for Diquark production***
C.. T: half a curtain pair = sqrt(curtain pair factor)
      IF(MSTJ(12).GE.5) THEN
         PMUD0=PYMASS(2101)
         PMUD1=PYMASS(2103)-PMUD0
         PMUS0=PYMASS(3201)-PMUD0
         PMUS1=PYMASS(3203)-PMUS0-PMUD0
         PMSS1=PYMASS(3303)-PMUS0-PMUD0
         QBB(ISU0)=EXP(-(PARJ(9)+PARJ(8))*PMUS0-PARJ(9)*PARF(191))
         QBB(IUS0)=EXP(-PARJ(8)*PMUS0)
         QBB(ISS1)=EXP(-(PARJ(9)+PARJ(8))*PMSS1)*QBB(ISU0)
         QBB(IUU1)=EXP(-PARJ(8)*PMUD1)
         QBB(ISU1)=EXP(-(PARJ(9)+PARJ(8))*PMUS1)*QBB(ISU0)
         QBB(IUS1)=EXP(-PARJ(8)*PMUS1)*QBB(IUS0)
         QBB(IUD1)=QBB(IUU1)
      ELSE
         PAR2M=SQRT(PARJ(2))
         PAR3M=SQRT(PARJ(3))
         PAR4M=SQRT(PARJ(4))
         QBB(ISU0)=PAR2M*PAR3M
         QBB(IUS0)=PAR3M
         QBB(ISS1)=PAR2M*PARJ(3)*PAR4M
         QBB(IUU1)=PAR4M
         QBB(ISU1)=PAR4M*QBB(ISU0)
         QBB(IUS1)=PAR4M*QBB(IUS0)
         QBB(IUD1)=PAR4M
      ENDIF
 
C.. tau: spin*(vertex factor)*(T = half-curtain factor)
      QBM(ISU0)=QBB(ISU0)
      QBM(IUS0)=PARJ(2)*QBB(IUS0)
      QBM(ISS1)=PARJ(2)*6D0*QBB(ISS1)
      QBM(IUU1)=6D0*QBB(IUU1)
      QBM(ISU1)=3D0*QBB(ISU1)
      QBM(IUS1)=PARJ(2)*3D0*QBB(IUS1)
      QBM(IUD1)=3D0*QBB(IUD1)
 
C.. Combine T and tau to diquark weight for q-> B+B+..
      DO 120 I=1,7
         QBB(I)=QBB(I)*QBM(I)
  120 CONTINUE
 
      IF(MSTJ(12).GE.5)THEN
C..New version: tau  for rank 0 diquark.
         DMB(7+ISU0)=EXP(-PARJ(10)*PMUS0)
         DMB(7+IUS0)=PARJ(2)*DMB(7+ISU0)
         DMB(7+ISS1)=6D0*PARJ(2)*EXP(-PARJ(10)*PMSS1)*DMB(7+ISU0)
         DMB(7+IUU1)=6D0*EXP(-PARJ(10)*PMUD1)
         DMB(7+ISU1)=3D0*EXP(-PARJ(10)*PMUS1)*DMB(7+ISU0)
         DMB(7+IUS1)=PARJ(2)*DMB(7+ISU1)
         DMB(7+IUD1)=DMB(7+IUU1)/2D0
 
C..New version: curtain flavour ratios.
C.. s/u for q->B+M+...
C.. s/u for rank 0 diquark: su -> ...M+B+...
C.. Q/q for heavy rank 0 diquark: Qu -> ...M+B+...
         WU=1D0+QBM(IUD1)+QBM(IUS0)+QBM(IUS1)+QBM(IUU1)
         PARF(135)=(2D0*(QBM(ISU0)+QBM(ISU1))+QBM(ISS1))/WU
         WU=1D0+DMB(7+IUD1)+DMB(7+IUS0)+DMB(7+IUS1)+DMB(7+IUU1)
         PARF(136)=(2D0*(DMB(7+ISU0)+DMB(7+ISU1))+DMB(7+ISS1))/WU
         PARF(137)=(DMB(7+ISU0)+DMB(7+ISU1))*
     &        (2D0+DMB(7+ISS1)/(2D0*DMB(7+ISU1)))/WU
      ELSE
C..Old version: reset unused rank 0 diquark weights and
C..             unused diquark SU(6) survival weights
         DO 130 I=1,7
            IF(MSTJ(12).LT.3) DMB(I)=1D0
            DMB(7+I)=1D0
  130    CONTINUE
 
C..Old version: Shuffle PARJ(7) into tau
         QBM(IUS0)=QBM(IUS0)*PARJ(7)
         QBM(ISS1)=QBM(ISS1)*PARJ(7)
         QBM(IUS1)=QBM(IUS1)*PARJ(7)
 
C..Old version: curtain flavour ratios.
C.. s/u for q->B+M+...
C.. s/u for rank 0 diquark: su -> ...M+B+...
C.. Q/q for heavy rank 0 diquark: Qu -> ...M+B+...
         WU=1D0+QBM(IUD1)+QBM(IUS0)+QBM(IUS1)+QBM(IUU1)
         PARF(135)=(2D0*(QBM(ISU0)+QBM(ISU1))+QBM(ISS1))/WU
         PARF(136)=PARF(135)*PARJ(6)*QBM(ISU0)/QBM(IUS0)
         PARF(137)=(1D0+QBM(IUD1))*(2D0+QBM(IUS0))/WU
      ENDIF
 
C..Combine diquark SU(6) survival, SU(6)max, tau and T into factors for:
C..  rank0 D->M+B+..; D->M+B+..; q->B+M+..; q->B+B..
      DO 140 I=1,7
         DMB(7+I)=DMB(7+I)*DMB(I)
         DMB(I)=DMB(I)*QBM(I)
         QBM(I)=QBM(I)*SU6M(I)/SU6MUD
         QBB(I)=QBB(I)*SU6M(I)/SU6MUD
  140 CONTINUE
 
C.. *** Popcorn factors ***
 
      IF(MSTJ(12).LT.5)THEN
C.. Old version: Resulting popcorn weights.
         PARF(138)=PARJ(6)
         WS=PARF(135)*PARF(138)
         WQ=WU*PARJ(5)/3D0
         PARF(132)=WQ*QBM(IUD1)/QBB(IUD1)
         PARF(133)=WQ*
     &        (QBM(IUS1)/QBB(IUS1)+WS*QBM(ISU1)/QBB(ISU1))/2D0
         PARF(134)=WQ*WS*QBM(ISS1)/QBB(ISS1)
         PARF(131)=WQ*(1D0+QBM(IUD1)+QBM(IUU1)+QBM(IUS0)+QBM(IUS1)+
     &                 WS*(QBM(ISU0)+QBM(ISU1)+QBM(ISS1)/2D0))/
     &        (1D0+QBB(IUD1)+QBB(IUU1)+
     &        2D0*(QBB(IUS0)+QBB(IUS1))+QBB(ISS1)/2D0)
      ELSE
C..New version: Store weights for popcorn mesons,
C..get prel. popcorn weights.
         DO 150 IPOS=201,1400
            PARF(IPOS)=0D0
  150    CONTINUE
         DO 160 I=138,140
            PARF(I)=0D0
  160    CONTINUE
         IPOS=200
         PARF(193)=PARJ(8)
         DO 240 MR=0,7,7
           IF(MR.EQ.7) PARF(193)=PARJ(10)
           SQWT=2D0*(DMB(MR+IUS0)+DMB(MR+IUS1))/
     &          (1D0+DMB(MR+IUD1)+DMB(MR+IUU1))
           QQWT=DMB(MR+IUU1)/(1D0+DMB(MR+IUD1)+DMB(MR+IUU1))
           DO 230 NMES=0,1
             IF(NMES.EQ.1) SQWT=PARJ(2)
             DO 220 KFQPOP=1,4
               IF(MR.EQ.0.AND.KFQPOP.GT.3) GOTO 220
               IF(NMES.EQ.0.AND.KFQPOP.GE.3)THEN
                  SQWT=DMB(MR+ISS1)/(DMB(MR+ISU0)+DMB(MR+ISU1))
                  QQWT=0.5D0
                  IF(MR.EQ.0) PARF(193)=PARJ(8)+PARJ(9)
                  IF(KFQPOP.EQ.4) SQWT=SQWT*(1D0/DMB(7+ISU1)+1D0)/2D0
               ENDIF
               DO 210 KFQOLD =1,5
                  IF(MR.EQ.0.AND.KFQOLD.GT.3) GOTO 210
                  IF(NMES.EQ.1) THEN
                     IF(MR.EQ.0.AND.KFQPOP.EQ.1) GOTO 210
                     IF(MR.EQ.7.AND.KFQPOP.NE.1) GOTO 210
                  ENDIF
                  WTTOT=0D0
                  WTFAIL=0D0
      DO 190 KMUL=0,5
         PJWT=PARJ(12+KMUL)
         IF(KMUL.EQ.0) PJWT=1D0-PARJ(14)
         IF(KMUL.EQ.1) PJWT=1D0-PARJ(15)-PARJ(16)-PARJ(17)
         IF(PJWT.LE.0D0) GOTO 190
         IF(PJWT.GT.1D0) PJWT=1D0
         IMES=5*KMUL
         IMIX=2*KFQOLD+10*KMUL
         KFJ=2*KMUL+1
         IF(KMUL.EQ.2) KFJ=10003
         IF(KMUL.EQ.3) KFJ=10001
         IF(KMUL.EQ.4) KFJ=20003
         IF(KMUL.EQ.5) KFJ=5
         DO 180 KFQVER =1,3
            KFLA=MAX(KFQOLD,KFQVER)
            KFLB=MIN(KFQOLD,KFQVER)
            SWT=PARJ(11+KFLA/3+KFLA/4)
            IF(KMUL.EQ.0.OR.KMUL.EQ.2) SWT=1D0-SWT
            SWT=SWT*PJWT
            QWT=SQWT/(2D0+SQWT)
            IF(KFQVER.LT.3)THEN
               IF(KFQVER.EQ.KFQPOP) QWT=(1D0-QWT)*QQWT
               IF(KFQVER.NE.KFQPOP) QWT=(1D0-QWT)*(1D0-QQWT)
            ENDIF
            IF(KFQVER.NE.KFQOLD)THEN
               IMES=IMES+1
               KFM=100*KFLA+10*KFLB+KFJ
               PMM=PMAS(PYCOMP(KFM),1)-PMAS(PYCOMP(KFM),3)
               PARF(IPOS+IMES)=QWT*SWT*EXP(-PARF(193)*PMM)
               WTTOT=WTTOT+PARF(IPOS+IMES)
            ELSE
               DO 170 ID=3,5
                  IF(ID.EQ.3) DWT=1D0-PARF(IMIX-1)
                  IF(ID.EQ.4) DWT=PARF(IMIX-1)-PARF(IMIX)
                  IF(ID.EQ.5) DWT=PARF(IMIX)
                  KFM=110*(ID-2)+KFJ
                  PMM=PMAS(PYCOMP(KFM),1)-PMAS(PYCOMP(KFM),3)
                  PARF(IPOS+5*KMUL+ID)=QWT*SWT*DWT*EXP(-PARF(193)*PMM)
                  IF(KMUL.EQ.0.AND.ID.GT.3) THEN
                     WTFAIL=WTFAIL+QWT*SWT*DWT*(1D0-PARJ(21+ID))
                     PARF(IPOS+5*KMUL+ID)=
     &                    PARF(IPOS+5*KMUL+ID)*PARJ(21+ID)
                  ENDIF
                  WTTOT=WTTOT+PARF(IPOS+5*KMUL+ID)
  170          CONTINUE
            ENDIF
  180    CONTINUE
  190 CONTINUE
                  DO 200 IMES=1,30
                     PARF(IPOS+IMES)=PARF(IPOS+IMES)/(1D0-WTFAIL)
  200             CONTINUE
                  IF(MR.EQ.7) PARF(140)=
     &                 MAX(PARF(140),WTTOT/(1D0-WTFAIL))
                  IF(MR.EQ.0) PARF(139-KFQPOP/3)=
     &                 MAX(PARF(139-KFQPOP/3),WTTOT/(1D0-WTFAIL))
                  IPOS=IPOS+30
  210           CONTINUE
  220         CONTINUE
  230       CONTINUE
  240    CONTINUE
         IF(PARF(139).GT.1D-10) PARF(138)=PARF(138)/PARF(139)
         MSTU(121)=0
 
      ENDIF
 
C..Recombine diquark weights to flavour and spin ratios
      PARF(151)=(2D0*(QBB(ISU0)+QBB(ISU1))+QBB(ISS1))/
     &        (1D0+QBB(IUD1)+QBB(IUU1)+QBB(IUS0)+QBB(IUS1))
      PARF(152)=2D0*(QBB(IUS0)+QBB(IUS1))/(1D0+QBB(IUD1)+QBB(IUU1))
      PARF(153)=QBB(ISS1)/(QBB(ISU0)+QBB(ISU1))
      PARF(154)=QBB(IUU1)/(1D0+QBB(IUD1)+QBB(IUU1))
      PARF(155)=QBB(ISU1)/QBB(ISU0)
      PARF(156)=QBB(IUS1)/QBB(IUS0)
      PARF(157)=QBB(IUD1)
 
      PARF(161)=(2D0*(QBM(ISU0)+QBM(ISU1))+QBM(ISS1))/
     &        (1D0+QBM(IUD1)+QBM(IUU1)+QBM(IUS0)+QBM(IUS1))
      PARF(162)=2D0*(QBM(IUS0)+QBM(IUS1))/(1D0+QBM(IUD1)+QBM(IUU1))
      PARF(163)=QBM(ISS1)/(QBM(ISU0)+QBM(ISU1))
      PARF(164)=QBM(IUU1)/(1D0+QBM(IUD1)+QBM(IUU1))
      PARF(165)=QBM(ISU1)/QBM(ISU0)
      PARF(166)=QBM(IUS1)/QBM(IUS0)
      PARF(167)=QBM(IUD1)
 
      PARF(171)=(2D0*(DMB(ISU0)+DMB(ISU1))+DMB(ISS1))/
     &        (1D0+DMB(IUD1)+DMB(IUU1)+DMB(IUS0)+DMB(IUS1))
      PARF(172)=2D0*(DMB(IUS0)+DMB(IUS1))/(1D0+DMB(IUD1)+DMB(IUU1))
      PARF(173)=DMB(ISS1)/(DMB(ISU0)+DMB(ISU1))
      PARF(174)=DMB(IUU1)/(1D0+DMB(IUD1)+DMB(IUU1))
      PARF(175)=DMB(ISU1)/DMB(ISU0)
      PARF(176)=DMB(IUS1)/DMB(IUS0)
      PARF(177)=DMB(IUD1)
 
      PARF(185)=DMB(7+ISU1)/DMB(7+ISU0)
      PARF(186)=DMB(7+IUS1)/DMB(7+IUS0)
      PARF(187)=DMB(7+IUD1)
 
      RETURN
      END
