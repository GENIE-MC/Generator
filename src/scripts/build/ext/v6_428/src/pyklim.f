 
C***********************************************************************
 
C...PYKLIM
C...Checks generated variables against pre-set kinematical limits;
C...also calculates limits on variables used in generation.
 
      SUBROUTINE PYKLIM(ILIM)
 
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
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT2/ISET(500),KFPR(500,2),COEF(500,20),ICOL(40,4,2)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYDAT3/,/PYSUBS/,/PYPARS/,
     &/PYINT1/,/PYINT2/
 
C...Common kinematical expressions.
      MINT(51)=0
      ISUB=MINT(1)
      ISTSB=ISET(ISUB)
      IF(ISUB.EQ.96) GOTO 100
      SQM3=VINT(63)
      SQM4=VINT(64)
      IF(ILIM.NE.0) THEN
        IF(ABS(SQM3).LT.1D-4.AND.ABS(SQM4).LT.1D-4) THEN
          CKIN09=MAX(CKIN(9),CKIN(13))
          CKIN10=MIN(CKIN(10),CKIN(14))
          CKIN11=MAX(CKIN(11),CKIN(15))
          CKIN12=MIN(CKIN(12),CKIN(16))
        ELSE
          CKIN09=MAX(CKIN(9),MIN(0D0,CKIN(13)))
          CKIN10=MIN(CKIN(10),MAX(0D0,CKIN(14)))
          CKIN11=MAX(CKIN(11),MIN(0D0,CKIN(15)))
          CKIN12=MIN(CKIN(12),MAX(0D0,CKIN(16)))
        ENDIF
      ENDIF
      IF(ILIM.NE.1) THEN
        TAU=VINT(21)
        RM3=SQM3/(TAU*VINT(2))
        RM4=SQM4/(TAU*VINT(2))
        BE34=SQRT(MAX(1D-20,(1D0-RM3-RM4)**2-4D0*RM3*RM4))
      ENDIF
      PTHMIN=CKIN(3)
      IF(MIN(SQM3,SQM4).LT.CKIN(6)**2.AND.ISTSB.NE.1.AND.ISTSB.NE.3)
     &PTHMIN=MAX(CKIN(3),CKIN(5))
 
      IF(ILIM.EQ.0) THEN
C...Check generated values of tau, y*, cos(theta-hat), and tau' against
C...pre-set kinematical limits.
        YST=VINT(22)
        CTH=VINT(23)
        TAUP=VINT(26)
        TAUE=TAU
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=TAUP
        X1=SQRT(TAUE)*EXP(YST)
        X2=SQRT(TAUE)*EXP(-YST)
        XF=X1-X2
        IF(MINT(47).NE.1) THEN
          IF(TAU*VINT(2).LT.CKIN(1)**2) MINT(51)=1
          IF(CKIN(2).GE.0D0.AND.TAU*VINT(2).GT.CKIN(2)**2) MINT(51)=1
          IF(YST.LT.CKIN(7).OR.YST.GT.CKIN(8)) MINT(51)=1
          IF(XF.LT.CKIN(25).OR.XF.GT.CKIN(26)) MINT(51)=1
        ENDIF
        IF(MINT(45).NE.1) THEN
          IF(X1.LT.CKIN(21).OR.X1.GT.CKIN(22)) MINT(51)=1
        ENDIF
        IF(MINT(46).NE.1) THEN
          IF(X2.LT.CKIN(23).OR.X2.GT.CKIN(24)) MINT(51)=1
        ENDIF
        IF(MINT(45).EQ.2) THEN
          IF(X1.GT.1D0-2D0*PARP(111)/VINT(1)) MINT(51)=1
        ENDIF
        IF(MINT(46).EQ.2) THEN
          IF(X2.GT.1D0-2D0*PARP(111)/VINT(1)) MINT(51)=1
        ENDIF
        IF(ISTSB.EQ.2.OR.ISTSB.EQ.4) THEN
          PTH=0.5D0*BE34*SQRT(TAU*VINT(2)*MAX(0D0,1D0-CTH**2))
          EXPY3=MAX(1D-20,(1D0+RM3-RM4+BE34*CTH)/
     &    MAX(1D-20,(1D0+RM3-RM4-BE34*CTH)))
          EXPY4=MAX(1D-20,(1D0-RM3+RM4-BE34*CTH)/
     &    MAX(1D-20,(1D0-RM3+RM4+BE34*CTH)))
          Y3=YST+0.5D0*LOG(EXPY3)
          Y4=YST+0.5D0*LOG(EXPY4)
          YLARGE=MAX(Y3,Y4)
          YSMALL=MIN(Y3,Y4)
          ETALAR=20D0
          ETASMA=-20D0
          STH=SQRT(MAX(0D0,1D0-CTH**2))
          EXSQ3=SQRT(MAX(1D-20,((1D0+RM3-RM4)*COSH(YST)+BE34*SINH(YST)*
     &    CTH)**2-4D0*RM3))
          EXSQ4=SQRT(MAX(1D-20,((1D0-RM3+RM4)*COSH(YST)-BE34*SINH(YST)*
     &    CTH)**2-4D0*RM4))
          IF(STH.GE.1D-10) THEN
            EXPET3=((1D0+RM3-RM4)*SINH(YST)+BE34*COSH(YST)*CTH+EXSQ3)/
     &      (BE34*STH)
            EXPET4=((1D0-RM3+RM4)*SINH(YST)-BE34*COSH(YST)*CTH+EXSQ4)/
     &      (BE34*STH)
            ETA3=LOG(MIN(1D10,MAX(1D-10,EXPET3)))
            ETA4=LOG(MIN(1D10,MAX(1D-10,EXPET4)))
            ETALAR=MAX(ETA3,ETA4)
            ETASMA=MIN(ETA3,ETA4)
          ENDIF
          CTS3=((1D0+RM3-RM4)*SINH(YST)+BE34*COSH(YST)*CTH)/EXSQ3
          CTS4=((1D0-RM3+RM4)*SINH(YST)-BE34*COSH(YST)*CTH)/EXSQ4
          CTSLAR=MIN(1D0,MAX(-1D0,CTS3,CTS4))
          CTSSMA=MAX(-1D0,MIN(1D0,CTS3,CTS4))
          SH=TAU*VINT(2)
          RPTS=4D0*VINT(71)**2/SH
          BE34L=SQRT(MAX(0D0,(1D0-RM3-RM4)**2-4D0*RM3*RM4-RPTS))
          RM34=MAX(1D-20,2D0*RM3*RM4)
          IF(2D0*VINT(71)**2/(VINT(21)*VINT(2)).LT.0.0001D0)
     &    RM34=MAX(RM34,2D0*VINT(71)**2/(VINT(21)*VINT(2)))
          RTHM=(4D0*RM3*RM4+RPTS)/(1D0-RM3-RM4+BE34L)
          THA=0.5D0*SH*MAX(RTHM,1D0-RM3-RM4-BE34*CTH)
          UHA=0.5D0*SH*MAX(RTHM,1D0-RM3-RM4+BE34*CTH)
          IF(PTH.LT.PTHMIN) MINT(51)=1
          IF(CKIN(4).GE.0D0.AND.PTH.GT.CKIN(4)) MINT(51)=1
          IF(YLARGE.LT.CKIN(9).OR.YLARGE.GT.CKIN(10)) MINT(51)=1
          IF(YSMALL.LT.CKIN(11).OR.YSMALL.GT.CKIN(12)) MINT(51)=1
          IF(ETALAR.LT.CKIN(13).OR.ETALAR.GT.CKIN(14)) MINT(51)=1
          IF(ETASMA.LT.CKIN(15).OR.ETASMA.GT.CKIN(16)) MINT(51)=1
          IF(CTSLAR.LT.CKIN(17).OR.CTSLAR.GT.CKIN(18)) MINT(51)=1
          IF(CTSSMA.LT.CKIN(19).OR.CTSSMA.GT.CKIN(20)) MINT(51)=1
          IF(CTH.LT.CKIN(27).OR.CTH.GT.CKIN(28)) MINT(51)=1
          IF(THA.LT.CKIN(35)) MINT(51)=1
          IF(CKIN(36).GE.0D0.AND.THA.GT.CKIN(36)) MINT(51)=1
          IF(UHA.LT.CKIN(37)) MINT(51)=1
          IF(CKIN(38).GE.0D0.AND.UHA.GT.CKIN(38)) MINT(51)=1
        ENDIF
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) THEN
          IF(TAUP*VINT(2).LT.CKIN(31)**2) MINT(51)=1
          IF(CKIN(32).GE.0D0.AND.TAUP*VINT(2).GT.CKIN(32)**2) MINT(51)=1
        ENDIF
 
C...Additional cuts on W2 (approximately) in DIS.
        IF(ISUB.EQ.10.AND.MINT(43).GE.2) THEN
          XBJ=X2
          IF(IABS(MINT(12)).LT.20) XBJ=X1
          Q2BJ=THA
          W2BJ=Q2BJ*(1D0-XBJ)/XBJ
          IF(W2BJ.LT.CKIN(39)) MINT(51)=1
          IF(CKIN(40).GT.0D0.AND.W2BJ.GT.CKIN(40)) MINT(51)=1
        ENDIF
 
      ELSEIF(ILIM.EQ.1) THEN
C...Calculate limits on tau
C...0) due to definition
        TAUMN0=0D0
        TAUMX0=1D0
C...1) due to limits on subsystem mass
        TAUMN1=CKIN(1)**2/VINT(2)
        TAUMX1=1D0
        IF(CKIN(2).GE.0D0) TAUMX1=CKIN(2)**2/VINT(2)
C...2) due to limits on pT-hat (and non-overlapping rapidity intervals)
        TM3=SQRT(SQM3+PTHMIN**2)
        TM4=SQRT(SQM4+PTHMIN**2)
        YDCOSH=1D0
        IF(CKIN09.GT.CKIN12) YDCOSH=COSH(CKIN09-CKIN12)
        TAUMN2=(TM3**2+2D0*TM3*TM4*YDCOSH+TM4**2)/VINT(2)
        TAUMX2=1D0
C...3) due to limits on pT-hat and cos(theta-hat)
        CTH2MN=MIN(CKIN(27)**2,CKIN(28)**2)
        CTH2MX=MAX(CKIN(27)**2,CKIN(28)**2)
        TAUMN3=0D0
        IF(CKIN(27)*CKIN(28).GT.0D0) TAUMN3=
     &  (SQRT(SQM3+PTHMIN**2/(1D0-CTH2MN))+
     &  SQRT(SQM4+PTHMIN**2/(1D0-CTH2MN)))**2/VINT(2)
        TAUMX3=1D0
        IF(CKIN(4).GE.0D0.AND.CTH2MX.LT.1D0) TAUMX3=
     &  (SQRT(SQM3+CKIN(4)**2/(1D0-CTH2MX))+
     &  SQRT(SQM4+CKIN(4)**2/(1D0-CTH2MX)))**2/VINT(2)
C...4) due to limits on x1 and x2
        TAUMN4=CKIN(21)*CKIN(23)
        TAUMX4=CKIN(22)*CKIN(24)
C...5) due to limits on xF
        TAUMN5=0D0
        TAUMX5=MAX(1D0-CKIN(25),1D0+CKIN(26))
C...6) due to limits on that and uhat
        TAUMN6=(SQM3+SQM4+CKIN(35)+CKIN(37))/VINT(2)
        TAUMX6=1D0
        IF(CKIN(36).GT.0D0.AND.CKIN(38).GT.0D0) TAUMX6=
     &  (SQM3+SQM4+CKIN(36)+CKIN(38))/VINT(2)
 
C...Net effect of all separate limits.
        VINT(11)=MAX(TAUMN0,TAUMN1,TAUMN2,TAUMN3,TAUMN4,TAUMN5,TAUMN6)
        VINT(31)=MIN(TAUMX0,TAUMX1,TAUMX2,TAUMX3,TAUMX4,TAUMX5,TAUMX6)
        IF(MINT(47).EQ.1.AND.(ISTSB.EQ.1.OR.ISTSB.EQ.2)) THEN
          VINT(11)=1D0-1D-9
          VINT(31)=1D0+1D-9
        ELSEIF(MINT(47).EQ.5) THEN
          VINT(31)=MIN(VINT(31),1D0-2D-10)
        ELSEIF(MINT(47).GE.6) THEN
          VINT(31)=MIN(VINT(31),1D0-1D-10)
        ENDIF
        IF(VINT(31).LE.VINT(11)) MINT(51)=1
 
      ELSEIF(ILIM.EQ.2) THEN
C...Calculate limits on y*
        TAUE=TAU
        IF(ISTSB.GE.3.AND.ISTSB.LE.5) TAUE=VINT(26)
        TAURT=SQRT(TAUE)
C...0) due to kinematics
        YSTMN0=LOG(TAURT)
        YSTMX0=-YSTMN0
C...1) due to explicit limits
        YSTMN1=CKIN(7)
        YSTMX1=CKIN(8)
C...2) due to limits on x1
        YSTMN2=LOG(MAX(TAUE,CKIN(21))/TAURT)
        YSTMX2=LOG(MAX(TAUE,CKIN(22))/TAURT)
C...3) due to limits on x2
        YSTMN3=-LOG(MAX(TAUE,CKIN(24))/TAURT)
        YSTMX3=-LOG(MAX(TAUE,CKIN(23))/TAURT)
C...4) due to limits on xF
        YEPMN4=0.5D0*ABS(CKIN(25))/TAURT
        YSTMN4=SIGN(LOG(MAX(1D-20,SQRT(1D0+YEPMN4**2)+YEPMN4)),CKIN(25))
        YEPMX4=0.5D0*ABS(CKIN(26))/TAURT
        YSTMX4=SIGN(LOG(MAX(1D-20,SQRT(1D0+YEPMX4**2)+YEPMX4)),CKIN(26))
C...5) due to simultaneous limits on y-large and y-small
        YEPSMN=(RM3-RM4)*SINH(CKIN09-CKIN11)
        YEPSMX=(RM3-RM4)*SINH(CKIN10-CKIN12)
        YDIFMN=ABS(LOG(MAX(1D-20,SQRT(1D0+YEPSMN**2)-YEPSMN)))
        YDIFMX=ABS(LOG(MAX(1D-20,SQRT(1D0+YEPSMX**2)-YEPSMX)))
        YSTMN5=0.5D0*(CKIN09+CKIN11-YDIFMN)
        YSTMX5=0.5D0*(CKIN10+CKIN12+YDIFMX)
C...6) due to simultaneous limits on cos(theta-hat) and y-large or
C...   y-small
        CTHLIM=SQRT(MAX(0D0,1D0-4D0*PTHMIN**2/(BE34**2*TAUE*VINT(2))))
        RZMN=BE34*MAX(CKIN(27),-CTHLIM)
        RZMX=BE34*MIN(CKIN(28),CTHLIM)
        YEX3MX=(1D0+RM3-RM4+RZMX)/MAX(1D-10,1D0+RM3-RM4-RZMX)
        YEX4MX=(1D0+RM4-RM3-RZMN)/MAX(1D-10,1D0+RM4-RM3+RZMN)
        YEX3MN=MAX(1D-10,1D0+RM3-RM4+RZMN)/(1D0+RM3-RM4-RZMN)
        YEX4MN=MAX(1D-10,1D0+RM4-RM3-RZMX)/(1D0+RM4-RM3+RZMX)
        YSTMN6=CKIN09-0.5D0*LOG(MAX(YEX3MX,YEX4MX))
        YSTMX6=CKIN12-0.5D0*LOG(MIN(YEX3MN,YEX4MN))
 
C...Net effect of all separate limits.
        VINT(12)=MAX(YSTMN0,YSTMN1,YSTMN2,YSTMN3,YSTMN4,YSTMN5,YSTMN6)
        VINT(32)=MIN(YSTMX0,YSTMX1,YSTMX2,YSTMX3,YSTMX4,YSTMX5,YSTMX6)
        IF(MINT(47).EQ.1) THEN
          VINT(12)=-1D-9
          VINT(32)=1D-9
        ELSEIF(MINT(47).EQ.2.OR.MINT(47).EQ.6) THEN
          VINT(12)=(1D0-1D-9)*YSTMX0
          VINT(32)=(1D0+1D-9)*YSTMX0
        ELSEIF(MINT(47).EQ.3.OR.MINT(47).EQ.7) THEN
          VINT(12)=-(1D0+1D-9)*YSTMX0
          VINT(32)=-(1D0-1D-9)*YSTMX0
        ELSEIF(MINT(47).EQ.5) THEN
          YSTEE=LOG((1D0-1D-10)/TAURT)
          VINT(12)=MAX(VINT(12),-YSTEE)
          VINT(32)=MIN(VINT(32),YSTEE)
        ENDIF
        IF(VINT(32).LE.VINT(12)) MINT(51)=1
 
      ELSEIF(ILIM.EQ.3) THEN
C...Calculate limits on cos(theta-hat)
        YST=VINT(22)
C...0) due to definition
        CTNMN0=-1D0
        CTNMX0=0D0
        CTPMN0=0D0
        CTPMX0=1D0
C...1) due to explicit limits
        CTNMN1=MIN(0D0,CKIN(27))
        CTNMX1=MIN(0D0,CKIN(28))
        CTPMN1=MAX(0D0,CKIN(27))
        CTPMX1=MAX(0D0,CKIN(28))
C...2) due to limits on pT-hat
        CTNMN2=-SQRT(MAX(0D0,1D0-4D0*PTHMIN**2/(BE34**2*TAU*VINT(2))))
        CTPMX2=-CTNMN2
        CTNMX2=0D0
        CTPMN2=0D0
        IF(CKIN(4).GE.0D0) THEN
          CTNMX2=-SQRT(MAX(0D0,1D0-4D0*CKIN(4)**2/
     &    (BE34**2*TAU*VINT(2))))
          CTPMN2=-CTNMX2
        ENDIF
C...3) due to limits on y-large and y-small
        CTNMN3=MIN(0D0,MAX((1D0+RM3-RM4)/BE34*TANH(CKIN11-YST),
     &  -(1D0-RM3+RM4)/BE34*TANH(CKIN10-YST)))
        CTNMX3=MIN(0D0,(1D0+RM3-RM4)/BE34*TANH(CKIN12-YST),
     &  -(1D0-RM3+RM4)/BE34*TANH(CKIN09-YST))
        CTPMN3=MAX(0D0,(1D0+RM3-RM4)/BE34*TANH(CKIN09-YST),
     &  -(1D0-RM3+RM4)/BE34*TANH(CKIN12-YST))
        CTPMX3=MAX(0D0,MIN((1D0+RM3-RM4)/BE34*TANH(CKIN10-YST),
     &  -(1D0-RM3+RM4)/BE34*TANH(CKIN11-YST)))
C...4) due to limits on that
        CTNMN4=-1D0
        CTNMX4=0D0
        CTPMN4=0D0
        CTPMX4=1D0
        SH=TAU*VINT(2)
        IF(CKIN(35).GT.0D0) THEN
          CTLIM=(1D0-RM3-RM4-2D0*CKIN(35)/SH)/BE34
          IF(CTLIM.GT.0D0) THEN
            CTPMX4=CTLIM
          ELSE
            CTPMX4=0D0
            CTNMX4=CTLIM
          ENDIF
        ENDIF
        IF(CKIN(36).GT.0D0) THEN
          CTLIM=(1D0-RM3-RM4-2D0*CKIN(36)/SH)/BE34
          IF(CTLIM.LT.0D0) THEN
            CTNMN4=CTLIM
          ELSE
            CTNMN4=0D0
            CTPMN4=CTLIM
          ENDIF
        ENDIF
C...5) due to limits on uhat
        CTNMN5=-1D0
        CTNMX5=0D0
        CTPMN5=0D0
        CTPMX5=1D0
        IF(CKIN(37).GT.0D0) THEN
          CTLIM=(2D0*CKIN(37)/SH-(1D0-RM3-RM4))/BE34
          IF(CTLIM.LT.0D0) THEN
            CTNMN5=CTLIM
          ELSE
            CTNMN5=0D0
            CTPMN5=CTLIM
          ENDIF
        ENDIF
        IF(CKIN(38).GT.0D0) THEN
          CTLIM=(2D0*CKIN(38)/SH-(1D0-RM3-RM4))/BE34
          IF(CTLIM.GT.0D0) THEN
            CTPMX5=CTLIM
          ELSE
            CTPMX5=0D0
            CTNMX5=CTLIM
          ENDIF
        ENDIF
 
C...Net effect of all separate limits.
        VINT(13)=MAX(CTNMN0,CTNMN1,CTNMN2,CTNMN3,CTNMN4,CTNMN5)
        VINT(33)=MIN(CTNMX0,CTNMX1,CTNMX2,CTNMX3,CTNMX4,CTNMX5)
        VINT(14)=MAX(CTPMN0,CTPMN1,CTPMN2,CTPMN3,CTPMN4,CTPMN5)
        VINT(34)=MIN(CTPMX0,CTPMX1,CTPMX2,CTPMX3,CTPMX4,CTPMX5)
        IF(VINT(33).LE.VINT(13).AND.VINT(34).LE.VINT(14)) MINT(51)=1

        IF(VINT(14).GT.VINT(34)) VINT(34)=VINT(14)
        IF(VINT(13).GT.VINT(33)) VINT(33)=VINT(13)

      ELSEIF(ILIM.EQ.4) THEN
C...Calculate limits on tau'
C...0) due to kinematics
        TAPMN0=TAU
        IF(ISTSB.EQ.5.AND.VINT(201).GT.0D0) THEN
          PQRAT=(VINT(201)+VINT(206))/VINT(1)
          TAPMN0=(SQRT(TAU)+PQRAT)**2
        ENDIF
        TAPMX0=1D0
C...1) due to explicit limits
        TAPMN1=CKIN(31)**2/VINT(2)
        TAPMX1=1D0
        IF(CKIN(32).GE.0D0) TAPMX1=CKIN(32)**2/VINT(2)
 
C...Net effect of all separate limits.
        VINT(16)=MAX(TAPMN0,TAPMN1)
        VINT(36)=MIN(TAPMX0,TAPMX1)
        IF(MINT(47).EQ.1) THEN
          VINT(16)=1D0-1D-9
          VINT(36)=1D0+1D-9
        ELSEIF(MINT(47).EQ.5) THEN
          VINT(36)=MIN(VINT(36),1D0-2D-10)
        ELSEIF(MINT(47).EQ.6.OR.MINT(47).EQ.7) THEN
          VINT(36)=MIN(VINT(36),1D0-1D-10)
        ENDIF
        IF(VINT(36).LE.VINT(16)) MINT(51)=1
 
      ENDIF
      RETURN
 
C...Special case for low-pT and multiple interactions:
C...effective kinematical limits for tau, y*, cos(theta-hat).
  100 IF(ILIM.EQ.0) THEN
      ELSEIF(ILIM.EQ.1) THEN
        IF(MSTP(82).LE.1) THEN
          VINT(11)=4D0*(PARP(81)*(VINT(1)/PARP(89))**PARP(90))**2/
     &    VINT(2)
        ELSE
          VINT(11)=(PARP(82)*(VINT(1)/PARP(89))**PARP(90))**2/VINT(2)
        ENDIF
        VINT(31)=1D0
      ELSEIF(ILIM.EQ.2) THEN
        VINT(12)=0.5D0*LOG(VINT(21))
        VINT(32)=-VINT(12)
      ELSEIF(ILIM.EQ.3) THEN
        IF(MSTP(82).LE.1) THEN
          ST2EFF=4D0*(PARP(81)*(VINT(1)/PARP(89))**PARP(90))**2/
     &    (VINT(21)*VINT(2))
        ELSE
          ST2EFF=0.01D0*(PARP(82)*(VINT(1)/PARP(89))**PARP(90))**2/
     &    (VINT(21)*VINT(2))
        ENDIF
        VINT(13)=-SQRT(MAX(0D0,1D0-ST2EFF))
        VINT(33)=0D0
        VINT(14)=0D0
        VINT(34)=-VINT(13)
      ENDIF
 
      RETURN
      END
