 
C*********************************************************************
 
C...PYTHRG
C...Calculates the mass eigenstates of the third generation sfermions.
C...Created:  5-31-96
 
      SUBROUTINE PYTHRG
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/,/PYSSMT/
 
C...Local variables.
      DOUBLE PRECISION BETA
      DOUBLE PRECISION AM2(2,2),RT(2,2),DI(2,2)
      DOUBLE PRECISION XMZ2,XMW2,TANB,XMU,COS2B,XMQL2,XMQR2
      DOUBLE PRECISION XMF,XMF2,DIFF,SAME,XMF12,XMF22,SMALL
      DOUBLE PRECISION ATR,AMQR,AMQL
      INTEGER ID1(3),ID2(3),ID3(3),ID4(3)
      INTEGER IF,I,J,II,JJ,IT,L
      LOGICAL DTERM
      DATA SMALL/1D-3/
      DATA ID1/10,10,13/
      DATA ID2/5,6,15/
      DATA ID3/15,16,17/
      DATA ID4/11,12,14/
      DATA DTERM/.TRUE./
 
      XMZ2=PMAS(23,1)**2
      XMW2=PMAS(24,1)**2
      TANB=RMSS(5)
      XMU=-RMSS(4)
      BETA=ATAN(TANB)
      COS2B=COS(2D0*BETA)
 
C...OPTION TO FIX T1, T2, B1 MASSES AND MIXINGS
 
      IOPT=IMSS(5)
      IF(IOPT.EQ.1) THEN
        CTT=DCOS(RMSS(27))
        CTT2=CTT**2
        STT=DSIN(RMSS(27))
        STT2=STT**2
        XM12=RMSS(10)**2
        XM22=RMSS(12)**2
        XMQL2=CTT2*XM12+STT2*XM22
        XMQR2=STT2*XM12+CTT2*XM22
        XMF2=PYMRUN(6,PMAS(6,1)**2)**2
        ATOP=-XMU/TANB+CTT*STT*(XM12-XM22)/SQRT(XMF2)
        RMSS(16)=ATOP
C......SUBTRACT OUT D-TERM AND FERMION MASS
        XMQL2=XMQL2-XMF2-(4D0*XMW2-XMZ2)*COS2B/6D0
        XMQR2=XMQR2-XMF2+(XMW2-XMZ2)*COS2B*2D0/3D0
        IF(XMQL2.GE.0D0) THEN
          RMSS(10)=SQRT(XMQL2)
        ELSE
          RMSS(10)=-SQRT(-XMQL2)
        ENDIF
        IF(XMQR2.GE.0D0) THEN
          RMSS(12)=SQRT(XMQR2)
        ELSE
          RMSS(12)=-SQRT(-XMQR2)
        ENDIF
 
C SAME FOR BOTTOM SQUARK
        CTT=DCOS(RMSS(26))
        CTT2=CTT**2
        STT=DSIN(RMSS(26))
        STT2=STT**2
        XM22=RMSS(11)**2
        XMF2=PYMRUN(5,PMAS(6,1)**2)**2
        XMQL2=SIGN(RMSS(10)**2,RMSS(10))-(2D0*XMW2+XMZ2)*COS2B/6D0+XMF2
        IF(ABS(CTT).GE..9999D0) THEN
          ABOT=-XMU*TANB
          XMQR2=RMSS(11)**2
        ELSEIF(ABS(CTT).LE.1D-4) THEN
          ABOT=-XMU*TANB
          XMQR2=RMSS(11)**2
        ELSE
          XM12=(XMQL2-STT2*XM22)/CTT2
          XMQR2=STT2*XM12+CTT2*XM22
          ABOT=-XMU*TANB+CTT*STT*(XM12-XM22)/SQRT(XMF2)
        ENDIF
        RMSS(15)=ABOT
C......SUBTRACT OUT D-TERM AND FERMION MASS
        XMQR2=XMQR2-(XMW2-XMZ2)*COS2B/3D0-XMF2
        IF(XMQR2.GE.0D0) THEN
          RMSS(11)=SQRT(XMQR2)
        ELSE
          RMSS(11)=-SQRT(-XMQR2)
        ENDIF
C SAME FOR TAU SLEPTON
        CTT=DCOS(RMSS(28))
        CTT2=CTT**2
        STT=DSIN(RMSS(28))
        STT2=STT**2
        XM12=RMSS(13)**2
        XM22=RMSS(14)**2
        XMQL2=CTT2*XM12+STT2*XM22
        XMQR2=STT2*XM12+CTT2*XM22
        XMFR=PMAS(15,1)
        XMF2=XMFR**2
        ATAU=-XMU*TANB+CTT*STT*(XM12-XM22)/SQRT(XMF2)
        RMSS(17)=ATAU
C......SUBTRACT OUT D-TERM AND FERMION MASS
        XMQL2=XMQL2-XMF2+(-.5D0*XMZ2+XMW2)*COS2B
        XMQR2=XMQR2-XMF2+(XMZ2-XMW2)*COS2B
        IF(XMQL2.GE.0D0) THEN
          RMSS(13)=SQRT(XMQL2)
        ELSE
          RMSS(13)=-SQRT(-XMQL2)
        ENDIF
        IF(XMQR2.GE.0D0) THEN
          RMSS(14)=SQRT(XMQR2)
        ELSE
          RMSS(14)=-SQRT(-XMQR2)
        ENDIF
      ENDIF
      DO 170 L=1,3
        AMQL=RMSS(ID1(L))
        IF(AMQL.LT.0D0) THEN
          XMQL2=-AMQL**2
        ELSE
          XMQL2=AMQL**2
        ENDIF
        ATR=RMSS(ID3(L))
        AMQR=RMSS(ID4(L))
        IF(AMQR.LT.0D0) THEN
          XMQR2=-AMQR**2
        ELSE
          XMQR2=AMQR**2
        ENDIF
        IF=ID2(L)
        XMF=PYMRUN(IF,PMAS(6,1)**2)
        XMF2=XMF**2
        AM2(1,1)=XMQL2+XMF2
        AM2(2,2)=XMQR2+XMF2
        IF(AM2(1,1).EQ.AM2(2,2)) AM2(2,2)=AM2(2,2)*1.00001D0
        IF(DTERM) THEN
          IF(L.EQ.1) THEN
            AM2(1,1)=AM2(1,1)-(2D0*XMW2+XMZ2)*COS2B/6D0
            AM2(2,2)=AM2(2,2)+(XMW2-XMZ2)*COS2B/3D0
            AM2(1,2)=XMF*(ATR+XMU*TANB)
          ELSEIF(L.EQ.2) THEN
            AM2(1,1)=AM2(1,1)+(4D0*XMW2-XMZ2)*COS2B/6D0
            AM2(2,2)=AM2(2,2)-(XMW2-XMZ2)*COS2B*2D0/3D0
            AM2(1,2)=XMF*(ATR+XMU/TANB)
          ELSEIF(L.EQ.3) THEN
            IF(IMSS(8).EQ.1) THEN
              AM2(1,1)=RMSS(6)**2
              AM2(2,2)=RMSS(7)**2
              AM2(1,2)=0D0
              RMSS(13)=RMSS(6)
              RMSS(14)=RMSS(7)
            ELSE
              AM2(1,1)=AM2(1,1)-(-.5D0*XMZ2+XMW2)*COS2B
              AM2(2,2)=AM2(2,2)-(XMZ2-XMW2)*COS2B
              AM2(1,2)=XMF*(ATR+XMU*TANB)
            ENDIF
          ENDIF
        ENDIF
        AM2(2,1)=AM2(1,2)
        DETM=AM2(1,1)*AM2(2,2)-AM2(2,1)**2
        IF(DETM.LT.0D0) THEN
          WRITE(MSTU(11),*) ID2(L),DETM,AM2
          CALL PYERRM(30,' NEGATIVE**2 MASS FOR SFERMION IN PYTHRG ')
        ENDIF
        SAME=0.5D0*(AM2(1,1)+AM2(2,2))
        DIFF=0.5D0*SQRT((AM2(1,1)-AM2(2,2))**2+4D0*AM2(1,2)*AM2(2,1))
        XMF12=SAME-DIFF
        XMF22=SAME+DIFF
        IT=0
        IF(XMF22-XMF12.GT.0D0) THEN
          RT(1,1) = SQRT(MAX(0D0,(XMF22-AM2(1,1))/(XMF22-XMF12)))
          RT(2,2) = RT(1,1)
          RT(1,2) = -SIGN(SQRT(MAX(0D0,1D0-RT(1,1)**2)),
     &    AM2(1,2)/(XMF22-XMF12))
          RT(2,1) = -RT(1,2)
        ELSE
          RT(1,1) = 1D0
          RT(2,2) = RT(1,1)
          RT(1,2) = 0D0
          RT(2,1) = -RT(1,2)
        ENDIF
  100   CONTINUE
        IT=IT+1
 
        DO 140 I=1,2
          DO 130 JJ=1,2
            DI(I,JJ)=0D0
            DO 120 II=1,2
              DO 110 J=1,2
                DI(I,JJ)=DI(I,JJ)+RT(I,J)*AM2(J,II)*RT(JJ,II)
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
  140   CONTINUE
 
        IF(DI(1,1).GT.DI(2,2)) THEN
          WRITE(MSTU(11),*) ' ERROR IN DIAGONALIZATION '
          WRITE(MSTU(11),*) L,SQRT(XMF12),SQRT(XMF22)
          WRITE(MSTU(11),*) AM2
          WRITE(MSTU(11),*) DI
          WRITE(MSTU(11),*) RT
          DI(1,1)=-RT(2,1)
          DI(2,2)=RT(1,2)
          DI(1,2)=-RT(2,2)
          DI(2,1)=RT(1,1)
          DO 160 I=1,2
            DO 150 J=1,2
              RT(I,J)=DI(I,J)
  150       CONTINUE
  160     CONTINUE
          GOTO 100
        ELSEIF(ABS(DI(1,2)*DI(2,1)/DI(1,1)/DI(2,2)).GT.SMALL) THEN
          WRITE(MSTU(11),*) ' ERROR IN DIAGONALIZATION,'//
     &    ' OFF DIAGONAL ELEMENTS '
          WRITE(MSTU(11),*) 'MASSES = ',L,SQRT(XMF12),SQRT(XMF22)
          WRITE(MSTU(11),*) DI
          WRITE(MSTU(11),*) ' ROTATION = ',RT
C...STOP
        ELSEIF(DI(1,1).LT.0D0.OR.DI(2,2).LT.0D0) THEN
          WRITE(MSTU(11),*) ' ERROR IN DIAGONALIZATION,'//
     &    ' NEGATIVE MASSES '
          CALL PYSTOP(111)
        ENDIF
        PMAS(PYCOMP(KSUSY1+IF),1)=SQRT(XMF12)
        PMAS(PYCOMP(KSUSY2+IF),1)=SQRT(XMF22)
        SFMIX(IF,1)=RT(1,1)
        SFMIX(IF,2)=RT(1,2)
        SFMIX(IF,3)=RT(2,1)
        SFMIX(IF,4)=RT(2,2)
  170 CONTINUE
 
C.....TAU SNEUTRINO MASS...L=3
 
      XARG=AM2(1,1)+XMW2*COS2B
      IF(XARG.LT.0D0) THEN
        WRITE(MSTU(11),*) ' PYTHRG:: TAU SNEUTRINO MASS IS NEGATIVE'//
     &  ' FROM THE SUM RULE. '
        WRITE(MSTU(11),*) '  TRY A SMALLER VALUE OF TAN(BETA). '
        RETURN
      ELSE
        PMAS(PYCOMP(KSUSY1+16),1)=SQRT(XARG)
      ENDIF
 
      RETURN
      END
