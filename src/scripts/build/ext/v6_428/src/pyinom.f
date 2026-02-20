C*********************************************************************
 
C...PYINOM
C...Finds the mass eigenstates and mixing matrices for neutralinos
C...and charginos.
 
      SUBROUTINE PYINOM
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYCOMP
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
      DOUBLE PRECISION XMW,XMZ,XM(4)
      DOUBLE PRECISION AR(5,5),WR(5),ZR(5,5),ZI(5,5),AI(5,5)
      DOUBLE PRECISION WI(5),FV1(5),FV2(5),FV3(5)
      DOUBLE PRECISION COSW,SINW
      DOUBLE PRECISION XMU
      DOUBLE PRECISION TANB,COSB,SINB
      DOUBLE PRECISION XM1,XM2,XM3,BETA
      DOUBLE PRECISION Q2,AEM,A1,A2,AQ,RM1,RM2
      DOUBLE PRECISION ARG,X0,X1,AX0,AX1,AT,BT
      DOUBLE PRECISION Y0,Y1,AMGX0,AM1X0,AMGX1,AM1X1
      DOUBLE PRECISION ARGX0,AR1X0,ARGX1,AR1X1
      DOUBLE PRECISION PYALPS,PYALEM
      DOUBLE PRECISION PYRNM3
      COMPLEX*16 CAR(4,4),CAI(4,4),CA1,CA2
      INTEGER IERR,INDEX(4),I,J,K,IOPT,ILR,KFNCHI(4)
      DATA KFNCHI/1000022,1000023,1000025,1000035/
 
      IOPT=IMSS(2)
      IF(IMSS(1).EQ.2) THEN
        IOPT=1
      ENDIF
C...M1, M2, AND M3 ARE INDEPENDENT
      IF(IOPT.EQ.0) THEN
        XM1=RMSS(1)
        XM2=RMSS(2)
        XM3=RMSS(3)
      ELSEIF(IOPT.GE.1) THEN
        Q2=PMAS(23,1)**2
        AEM=PYALEM(Q2)
        A2=AEM/PARU(102)
        A1=AEM/(1D0-PARU(102))
        XM1=RMSS(1)
        XM2=RMSS(2)
        IF(IMSS(1).EQ.2) XM1=RMSS(1)/RMSS(20)*A1*5D0/3D0
        IF(IOPT.EQ.1) THEN
          XM2=XM1*A2/A1*3D0/5D0
          RMSS(2)=XM2
        ELSEIF(IOPT.EQ.3) THEN
          XM1=XM2*5D0/3D0*A1/A2
          RMSS(1)=XM1
        ENDIF
        XM3=PYRNM3(XM2/A2)
        RMSS(3)=XM3
        IF(XM3.LE.0D0) THEN
          WRITE(MSTU(11),*) ' ERROR WITH M3 = ',XM3
          CALL PYSTOP(105)
        ENDIF
      ENDIF
 
C...GLUINO MASS
      IF(IMSS(3).EQ.1) THEN
        PMAS(PYCOMP(KSUSY1+21),1)=ABS(XM3)
      ELSE
        AQ=0D0
        DO 110 I=1,4
          DO 100 ILR=1,2
            RM1=PMAS(PYCOMP(ILR*KSUSY1+I),1)**2/XM3**2
            AQ=AQ+0.5D0*((2D0-RM1)*(RM1*LOG(RM1)-1D0)
     &      +(1D0-RM1)**2*LOG(ABS(1D0-RM1)))
  100     CONTINUE
  110   CONTINUE
 
        DO 130 I=5,6
          DO 120 ILR=1,2
            RM1=PMAS(PYCOMP(ILR*KSUSY1+I),1)**2/XM3**2
            RM2=PMAS(I,1)**2/XM3**2
            ARG=(RM1-RM2-1D0)**2-4D0*RM2**2
            IF(ARG.GE.0D0) THEN
              X0=0.5D0*(1D0+RM2-RM1-SQRT(ARG))
              AX0=ABS(X0)
              X1=0.5D0*(1D0+RM2-RM1+SQRT(ARG))
              AX1=ABS(X1)
              IF(X0.EQ.1D0) THEN
                AT=-1D0
                BT=0.25D0
              ELSEIF(X0.EQ.0D0) THEN
                AT=0D0
                BT=-0.25D0
              ELSE
                AT=0.5D0*LOG(ABS(1D0-X0))*(1D0-X0**2)+
     &          0.5D0*X0**2*LOG(AX0)
                BT=(-1D0-2D0*X0)/4D0
              ENDIF
              IF(X1.EQ.1D0) THEN
                AT=-1D0+AT
                BT=0.25D0+BT
              ELSEIF(X1.EQ.0D0) THEN
                AT=0D0+AT
                BT=-0.25D0+BT
              ELSE
                AT=0.5D0*LOG(ABS(1D0-X1))*(1D0-X1**2)+0.5D0*
     &          X1**2*LOG(AX1)+AT
                BT=(-1D0-2D0*X1)/4D0+BT
              ENDIF
              AQ=AQ+AT+BT
            ELSE
              X0=0.5D0*(1D0+RM2-RM1)
              Y0=-0.5D0*SQRT(-ARG)
              AMGX0=SQRT(X0**2+Y0**2)
              AM1X0=SQRT((1D0-X0)**2+Y0**2)
              ARGX0=ATAN2(-X0,-Y0)
              AR1X0=ATAN2(1D0-X0,Y0)
              X1=X0
              Y1=-Y0
              AMGX1=AMGX0
              AM1X1=AM1X0
              ARGX1=ATAN2(-X1,-Y1)
              AR1X1=ATAN2(1D0-X1,Y1)
              AT=0.5D0*LOG(AM1X0)*(1D0-X0**2+3D0*Y0**2)
     &        +0.5D0*(X0**2-Y0**2)*LOG(AMGX0)
              BT=(-1D0-2D0*X0)/4D0+X0*Y0*( AR1X0-ARGX0 )
              AT=AT+0.5D0*LOG(AM1X1)*(1D0-X1**2+3D0*Y1**2)
     &        +0.5D0*(X1**2-Y1**2)*LOG(AMGX1)
              BT=BT+(-1D0-2D0*X1)/4D0+X1*Y1*( AR1X1-ARGX1 )
              AQ=AQ+AT+BT
            ENDIF
  120     CONTINUE
  130   CONTINUE
        PMAS(PYCOMP(KSUSY1+21),1)=ABS(XM3)*(1D0+PYALPS(XM3**2)
     &  /(2D0*PARU(2))*(15D0+AQ))
      ENDIF
 
C...NEUTRALINO MASSES
      DO 150 I=1,4
        DO 140 J=1,4
          AI(I,J)=0D0
  140   CONTINUE
  150 CONTINUE
      XMZ=PMAS(23,1)/100D0
      XMW=PMAS(24,1)/100D0
      XMU=RMSS(4)/100D0
      SINW=SQRT(PARU(102))
      COSW=SQRT(1D0-PARU(102))
      TANB=RMSS(5)
      BETA=ATAN(TANB)
      COSB=COS(BETA)
      SINB=TANB*COSB

      XM2=XM2/100D0
      XM1=XM1/100D0
      
 
C... Definitions:
C...    psi^0 =(-i bino^0, -i wino^0, h_d^0(=H_1^0), h_u^0(=H_2^0))
C... => L_neutralino = -1/2*(psi^0)^T * [AR] * psi^0 + h.c.
      AR(1,1) = XM1*COS(RMSS(30))
      AI(1,1) = XM1*SIN(RMSS(30))
      AR(2,2) = XM2*COS(RMSS(31))
      AI(2,2) = XM2*SIN(RMSS(31))
      AR(3,3) = 0D0
      AR(4,4) = 0D0
      AR(1,2) = 0D0
      AR(2,1) = 0D0
      AR(1,3) = -XMZ*SINW*COSB
      AR(3,1) = AR(1,3)
      AR(1,4) = XMZ*SINW*SINB
      AR(4,1) = AR(1,4)
      AR(2,3) = XMZ*COSW*COSB
      AR(3,2) = AR(2,3)
      AR(2,4) = -XMZ*COSW*SINB
      AR(4,2) = AR(2,4)
      AR(3,4) = -XMU*COS(RMSS(33))
      AI(3,4) = -XMU*SIN(RMSS(33))
      AR(4,3) = -XMU*COS(RMSS(33))
      AI(4,3) = -XMU*SIN(RMSS(33))
C      CALL PYEIG4(AR,WR,ZR)
      CALL PYEICG(5,4,AR,AI,WR,WI,1,ZR,ZI,FV1,FV2,FV3,IERR)
      IF(IERR.NE.0) CALL PYERRM(18,'(PYINOM:) '//
     & 'PROBLEM WITH PYEICG IN PYINOM ')
      DO 160 I=1,4
        INDEX(I)=I
        XM(I)=ABS(WR(I))
  160 CONTINUE
      DO 180 I=2,4
        K=I
        DO 170 J=I-1,1,-1
          IF(XM(K).LT.XM(J)) THEN
            ITMP=INDEX(J)
            XTMP=XM(J)
            INDEX(J)=INDEX(K)
            XM(J)=XM(K)
            INDEX(K)=ITMP
            XM(K)=XTMP
            K=K-1
          ELSE
            GOTO 180
          ENDIF
  170   CONTINUE
  180 CONTINUE
 
 
      DO 210 I=1,4
        K=INDEX(I)
        SMZ(I)=WR(K)*100D0
        PMAS(PYCOMP(KFNCHI(I)),1)=ABS(SMZ(I))
        S=0D0
        DO 190 J=1,4
          S=S+ZR(J,K)**2+ZI(J,K)**2
  190   CONTINUE
        DO 200 J=1,4
          ZMIX(I,J)=ZR(J,K)/SQRT(S)
          ZMIXI(I,J)=ZI(J,K)/SQRT(S)
          IF(ABS(ZMIX(I,J)).LT.1D-6) ZMIX(I,J)=0D0
          IF(ABS(ZMIXI(I,J)).LT.1D-6) ZMIXI(I,J)=0D0
  200   CONTINUE
  210 CONTINUE
 
C...CHARGINO MASSES
C.....Find eigenvectors of X X^*
      DO I=1,4
        DO J=1,4
          AR(I,J)=0D0
          AI(I,J)=0D0
        ENDDO
      ENDDO
      AI(1,1) = 0D0
      AI(2,2) = 0D0
      AR(1,1) = XM2**2+2D0*XMW**2*SINB**2
      AR(2,2) = XMU**2+2D0*XMW**2*COSB**2
      AR(1,2) = SQRT(2D0)*XMW*(XM2*COS(RMSS(31))*COSB+
     &XMU*COS(RMSS(33))*SINB)
      AI(1,2) = SQRT(2D0)*XMW*(XM2*SIN(RMSS(31))*COSB-
     &XMU*SIN(RMSS(33))*SINB)
      AR(2,1) = SQRT(2D0)*XMW*(XM2*COS(RMSS(31))*COSB+
     &XMU*COS(RMSS(33))*SINB)
      AI(2,1) = SQRT(2D0)*XMW*(-XM2*SIN(RMSS(31))*COSB+
     &XMU*SIN(RMSS(33))*SINB)
      CALL PYEICG(5,2,AR,AI,WR,WI,1,ZR,ZI,FV1,FV2,FV3,IERR)
      IF(IERR.NE.0) CALL PYERRM(18,'(PYINOM:) '//
     & 'PROBLEM WITH PYEICG IN PYINOM ')
      INDEX(1)=1
      INDEX(2)=2
      IF(WR(2).LT.WR(1)) THEN
        INDEX(1)=2
        INDEX(2)=1
      ENDIF

 
      DO 240 I=1,2
        K=INDEX(I)
        SMW(I)=SQRT(WR(K))*100D0
        S=0D0
        DO 220 J=1,2
          S=S+ZR(J,K)**2+ZI(J,K)**2
  220   CONTINUE
        DO 230 J=1,2
          UMIX(I,J)=ZR(J,K)/SQRT(S)
          UMIXI(I,J)=-ZI(J,K)/SQRT(S)
          IF(ABS(UMIX(I,J)).LT.1D-6) UMIX(I,J)=0D0
          IF(ABS(UMIXI(I,J)).LT.1D-6) UMIXI(I,J)=0D0
  230   CONTINUE
  240 CONTINUE
C...Force chargino mass > neutralino mass
      IFRC=0
      IF(ABS(SMW(1)).LT.ABS(SMZ(1))+2D0*PMAS(PYCOMP(111),1)) THEN
        CALL PYERRM(8,'(PYINOM:) '//
     &      'forcing m(~chi+_1) > m(~chi0_1) + 2m(pi0)')
        SMW(1)=SIGN(ABS(SMZ(1))+2D0*PMAS(PYCOMP(111),1),SMW(1))
        IFRC=1
      ENDIF
      PMAS(PYCOMP(KSUSY1+24),1)=SMW(1)
      PMAS(PYCOMP(KSUSY1+37),1)=SMW(2)
 
C.....Find eigenvectors of X^* X
      DO I=1,4
        DO J=1,4
          AR(I,J)=0D0
          AI(I,J)=0D0
          ZR(I,J)=0D0
          ZI(I,J)=0D0
        ENDDO
      ENDDO
      AI(1,1) = 0D0
      AI(2,2) = 0D0
      AR(1,1) = XM2**2+2D0*XMW**2*COSB**2
      AR(2,2) = XMU**2+2D0*XMW**2*SINB**2
      AR(1,2) = SQRT(2D0)*XMW*(XM2*COS(RMSS(31))*SINB+
     &XMU*COS(RMSS(33))*COSB)
      AI(1,2) = SQRT(2D0)*XMW*(-XM2*SIN(RMSS(31))*SINB+
     &XMU*SIN(RMSS(33))*COSB)
      AR(2,1) = SQRT(2D0)*XMW*(XM2*COS(RMSS(31))*SINB+
     &XMU*COS(RMSS(33))*COSB)
      AI(2,1) = SQRT(2D0)*XMW*(XM2*SIN(RMSS(31))*SINB-
     &XMU*SIN(RMSS(33))*COSB)
      CALL PYEICG(5,2,AR,AI,WR,WI,1,ZR,ZI,FV1,FV2,FV3,IERR)
      IF(IERR.NE.0) CALL PYERRM(18,'(PYINOM:) '//
     & 'PROBLEM WITH PYEICG IN PYINOM ')
      INDEX(1)=1
      INDEX(2)=2
      IF(WR(2).LT.WR(1)) THEN
        INDEX(1)=2
        INDEX(2)=1
      ENDIF
 
      SIMAG=0D0
      DO 270 I=1,2
        K=INDEX(I)
        S=0D0
        DO 250 J=1,2
          S=S+ZR(J,K)**2+ZI(J,K)**2
          SIMAG=SIMAG+ZI(J,K)**2
  250   CONTINUE
        DO 260 J=1,2
          VMIX(I,J)=ZR(J,K)/SQRT(S)
          VMIXI(I,J)=-ZI(J,K)/SQRT(S)
          IF(ABS(VMIX(I,J)).LT.1D-6) VMIX(I,J)=0D0
          IF(ABS(VMIXI(I,J)).LT.1D-6) VMIXI(I,J)=0D0
  260   CONTINUE
  270 CONTINUE

C.....Simplify if no phases
      IF(SIMAG.LT.1D-6) THEN
        AR(1,1) = XM2*COS(RMSS(31))
        AR(2,2) = XMU*COS(RMSS(33))
        AR(1,2) = SQRT(2D0)*XMW*SINB
        AR(2,1) = SQRT(2D0)*XMW*COSB
        IKNT=0
 300    CONTINUE
        DO I=1,2
          DO J=1,2
            ZR(I,J)=0D0
          ENDDO
        ENDDO

        DO I=1,2
          DO J=1,2
            DO K=1,2
              DO L=1,2
                ZR(I,J)=ZR(I,J)+UMIX(I,K)*AR(K,L)*VMIX(J,L)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        VMIX(1,1)=VMIX(1,1)*SMW(1)/ZR(1,1)/100D0
        VMIX(1,2)=VMIX(1,2)*SMW(1)/ZR(1,1)/100D0
        VMIX(2,1)=VMIX(2,1)*SMW(2)/ZR(2,2)/100D0
        VMIX(2,2)=VMIX(2,2)*SMW(2)/ZR(2,2)/100D0
        IF(IKNT.EQ.2.AND.IFRC.EQ.0) THEN
          CALL PYERRM(18,'(PYINOM:) Problem with Charginos')
        ELSEIF(ZR(1,1).LT.0D0.OR.ZR(2,2).LT.0D0) THEN
          IKNT=IKNT+1
          GOTO 300
        ENDIF
C.....Must deal with phases
      ELSE
        CAR(1,1) = XM2*CMPLX(COS(RMSS(31)),SIN(RMSS(31)))
        CAR(2,2) = XMU*CMPLX(COS(RMSS(33)),SIN(RMSS(33)))
        CAR(1,2) = SQRT(2D0)*XMW*SINB*CMPLX(1D0,0D0)
        CAR(2,1) = SQRT(2D0)*XMW*COSB*CMPLX(1D0,0D0)

        IKNT=0
 310    CONTINUE
        DO I=1,2
          DO J=1,2
            CAI(I,J)=CMPLX(0D0,0D0)
          ENDDO
        ENDDO

        DO I=1,2
          DO J=1,2
            DO K=1,2
              DO L=1,2
                CAI(I,J)=CAI(I,J)+CMPLX(UMIX(I,K),-UMIXI(I,K))*CAR(K,L)*
     &           CMPLX(VMIX(J,L),VMIXI(J,L))
              ENDDO
            ENDDO
          ENDDO
        ENDDO

        CA1=SMW(1)*CAI(1,1)/ABS(CAI(1,1))**2/100D0
        CA2=SMW(2)*CAI(2,2)/ABS(CAI(2,2))**2/100D0
        TEMPR=VMIX(1,1)
        TEMPI=VMIXI(1,1)
        VMIX(1,1)=TEMPR*DBLE(CA1)-TEMPI*DIMAG(CA1)
        VMIXI(1,1)=TEMPI*DBLE(CA1)+TEMPR*DIMAG(CA1)
        TEMPR=VMIX(1,2)
        TEMPI=VMIXI(1,2)
        VMIX(1,2)=TEMPR*DBLE(CA1)-TEMPI*DIMAG(CA1)
        VMIXI(1,2)=TEMPI*DBLE(CA1)+TEMPR*DIMAG(CA1)
        TEMPR=VMIX(2,1)
        TEMPI=VMIXI(2,1)
        VMIX(2,1)=TEMPR*DBLE(CA2)-TEMPI*DIMAG(CA2)
        VMIXI(2,1)=TEMPI*DBLE(CA2)+TEMPR*DIMAG(CA2)
        TEMPR=VMIX(2,2)
        TEMPI=VMIXI(2,2)
        VMIX(2,2)=TEMPR*DBLE(CA2)-TEMPI*DIMAG(CA2)
        VMIXI(2,2)=TEMPI*DBLE(CA2)+TEMPR*DIMAG(CA2)
        IF(IKNT.EQ.2.AND.IFRC.EQ.0) THEN
          CALL PYERRM(18,'(PYINOM:) Problem with Charginos')
        ELSEIF(DBLE(CA1).LT.0D0.OR.DBLE(CA2).LT.0D0.OR.
     &   ABS(IMAG(CA1)).GT.1D-3.OR.ABS(IMAG(CA2)).GT.1D-3) THEN
          IKNT=IKNT+1
          GOTO 310
        ENDIF
      ENDIF 
      RETURN
      END
