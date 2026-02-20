 
C*********************************************************************
 
C...PYPOLE
C...This subroutine computes the CP-even higgs and CP-odd pole
c...Higgs masses and mixing angles.
 
C...Program based on the work by M. Carena, M. Quiros
C...and C.E.M. Wagner, "Effective potential methods and
C...the Higgs mass spectrum in the MSSM", CERN-TH/95-157
 
C...Inputs: IHIGGS(explained below),MCHI,MA,TANB,MQ,MUR,MDR,MTOP,
C...AT,AB,MU
C...where MCHI is the largest chargino mass, MA is the running
C...CP-odd higgs mass, TANB is the value of the ratio of vacuum
C...expectaion values at the scale MTOP, MQ is the third generation
C...left handed squark mass parameter, MUR is the third generation
C...right handed stop mass parameter, MDR is the third generation
C...right handed sbottom mass parameter, MTOP is the pole top quark
C...mass; AT,AB are the soft supersymmetry breaking trilinear
C...couplings of the stop and sbottoms, respectively, and MU is the
C...supersymmetric mass parameter
 
C...The parameter IHIGGS=0,1,2,3 corresponds to the number of
C...Higgses whose pole mass is computed. If IHIGGS=0 only running
C...masses are given, what makes the running of the program
c...much faster and it is quite generally a good approximation
c...(for a theoretical discussion see ref. above). If IHIGGS=1,
C...only the pole mass for H is computed. If IHIGGS=2, then h and H,
c...and if IHIGGS=3, then h,H,A polarizations are computed
 
C...Output: MH and MHP which are the lightest CP-even Higgs running
C...and pole masses, respectively; HM and HMP are the heaviest CP-even
C...Higgs running and pole masses, repectively; SA and CA are the
C...SIN(ALPHA) and COS(ALPHA) where ALPHA is the Higgs mixing angle
C...AMP is the CP-odd Higgs pole mass. STOP1,STOP2,SBOT1 and SBOT2
C...are the stop and sbottom mass eigenvalues. Finally, TANBA is
C...the value of TANB at the CP-odd Higgs mass scale
 
C...This subroutine makes use of CERN library subroutine
C...integration package, which makes the computation of the
C...pole Higgs masses somewhat faster. We thank P. Janot for this
C...improvement. Those who are not able to call the CERN
C...libraries, please use the subroutine SUBHPOLE2.F, which
C...although somewhat slower, gives identical results
 
      SUBROUTINE PYPOLE(IHIGGS,XMC,XMA,TANB,XMQ,XMUR,XMDR,XMT,AT,AB,XMU,
     &XMH,XMHP,HM,HMP,AMP,SA,CA,STOP1,STOP2,SBOT1,SBOT2,TANBA,XMG,DT,DB)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
C...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
      INTEGER PYK,PYCHGE,PYCOMP
 
C...Local variables.
      DIMENSION DELTA(2,2),COUPT(2,2),T(2,2),SSTOP2(2),
     &SSBOT2(2),B(2,2),COUPB(2,2),
     &HCOUPT(2,2),HCOUPB(2,2),
     &ACOUPT(2,2),ACOUPB(2,2),PR(3), POLAR(3)
 
      DELTA(1,1) = 1D0
      DELTA(2,2) = 1D0
      DELTA(1,2) = 0D0
      DELTA(2,1) = 0D0
      V = 174.1D0
      XMZ=91.18D0
      PI=PARU(1)
      RXMT=PYMRUN(6,XMT**2)
      CALL PYRGHM(XMC,XMA,TANB,XMQ,XMUR,XMDR,XMT,AT,AB,
     &XMU,XMH,HM,XMCH,SA,CA,SAB,CAB,TANBA,XMG,DT,DB)
 
      SINB = TANB/(TANB**2+1D0)**0.5D0
      COSB = 1D0/(TANB**2+1D0)**0.5D0
      COS2B = SINB**2 - COSB**2
      SINBPA = SINB*CA + COSB*SA
      COSBPA = COSB*CA - SINB*SA
      RMBOT = PYMRUN(5,XMT**2)
      XMQ2 = XMQ**2
      XMUR2 = XMUR**2
      IF(XMUR.LT.0D0) XMUR2=-XMUR2
      XMDR2 = XMDR**2
      XMST11 = RXMT**2 + XMQ2  - 0.35D0*XMZ**2*COS2B
      XMST22 = RXMT**2 + XMUR2 - 0.15D0*XMZ**2*COS2B
      IF(XMST11.LT.0D0) GOTO 500
      IF(XMST22.LT.0D0) GOTO 500
      XMSB11 = RMBOT**2 + XMQ2  + 0.42D0*XMZ**2*COS2B
      XMSB22 = RMBOT**2 + XMDR2 + 0.08D0*XMZ**2*COS2B
      IF(XMSB11.LT.0D0) GOTO 500
      IF(XMSB22.LT.0D0) GOTO 500
C      WMST11 = RXMT**2 + XMQ2
C      WMST22 = RXMT**2 + XMUR2
      XMST12 = RXMT*(AT - XMU/TANB)
      XMSB12 = RMBOT*(AB - XMU*TANB)
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...STOP EIGENVALUES CALCULATION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      STOP12 = 0.5D0*(XMST11+XMST22) +
     &0.5D0*((XMST11+XMST22)**2 -
     &4D0*(XMST11*XMST22 - XMST12**2))**0.5D0
      STOP22 = 0.5D0*(XMST11+XMST22) -
     &0.5D0*((XMST11+XMST22)**2 - 4D0*(XMST11*XMST22 -
     &XMST12**2))**0.5D0
 
      IF(STOP22.LT.0D0) GOTO 500
      SSTOP2(1) = STOP12
      SSTOP2(2) = STOP22
      STOP1 = STOP12**0.5D0
      STOP2 = STOP22**0.5D0
C      STOP1W = STOP1
C      STOP2W = STOP2
 
      IF(XMST12.EQ.0D0) XST11 = 1D0
      IF(XMST12.EQ.0D0) XST12 = 0D0
      IF(XMST12.EQ.0D0) XST21 = 0D0
      IF(XMST12.EQ.0D0) XST22 = 1D0
 
      IF(XMST12.EQ.0D0) GOTO 110
 
  100 XST11 = XMST12/(XMST12**2+(XMST11-STOP12)**2)**0.5D0
      XST12 = - (XMST11-STOP12)/(XMST12**2+(XMST11-STOP12)**2)**0.5D0
      XST21 = XMST12/(XMST12**2+(XMST11-STOP22)**2)**0.5D0
      XST22 = - (XMST11-STOP22)/(XMST12**2+(XMST11-STOP22)**2)**0.5D0
 
  110 T(1,1) = XST11
      T(2,2) = XST22
      T(1,2) = XST12
      T(2,1) = XST21
 
      SBOT12 = 0.5D0*(XMSB11+XMSB22) +
     &0.5D0*((XMSB11+XMSB22)**2 -
     &4D0*(XMSB11*XMSB22 - XMSB12**2))**0.5D0
      SBOT22 = 0.5D0*(XMSB11+XMSB22) -
     &0.5D0*((XMSB11+XMSB22)**2 - 4D0*(XMSB11*XMSB22 -
     &XMSB12**2))**0.5D0
      IF(SBOT22.LT.0D0) GOTO 500
      SBOT1 = SBOT12**0.5D0
      SBOT2 = SBOT22**0.5D0
 
      SSBOT2(1) = SBOT12
      SSBOT2(2) = SBOT22
 
      IF(XMSB12.EQ.0D0) XSB11 = 1D0
      IF(XMSB12.EQ.0D0) XSB12 = 0D0
      IF(XMSB12.EQ.0D0) XSB21 = 0D0
      IF(XMSB12.EQ.0D0) XSB22 = 1D0
 
      IF(XMSB12.EQ.0D0) GOTO 130
 
  120 XSB11 = XMSB12/(XMSB12**2+(XMSB11-SBOT12)**2)**0.5D0
      XSB12 = - (XMSB11-SBOT12)/(XMSB12**2+(XMSB11-SBOT12)**2)**0.5D0
      XSB21 = XMSB12/(XMSB12**2+(XMSB11-SBOT22)**2)**0.5D0
      XSB22 = - (XMSB11-SBOT22)/(XMSB12**2+(XMSB11-SBOT22)**2)**0.5D0
 
  130 B(1,1) = XSB11
      B(2,2) = XSB22
      B(1,2) = XSB12
      B(2,1) = XSB21
 
 
      SINT = 0.2320D0
      SQR = DSQRT(2D0)
      VP = 174.1D0*SQR
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...STARTING OF LIGHT HIGGS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      IF(IHIGGS.EQ.0) GOTO 490
 
      DO 150 I = 1,2
        DO 140 J = 1,2
          COUPT(I,J) =
     &    SINT*XMZ**2*2D0*SQR/174.1D0/3D0*SINBPA*(DELTA(I,J) +
     &    (3D0 - 8D0*SINT)/4D0/SINT*T(1,I)*T(1,J))
     &    -RXMT**2/174.1D0**2*VP/SINB*CA*DELTA(I,J)
     &    -RXMT/VP/SINB*(AT*CA + XMU*SA)*(T(1,I)*T(2,J) +
     &    T(1,J)*T(2,I))
  140   CONTINUE
  150 CONTINUE
 
 
      DO 170 I = 1,2
        DO 160 J = 1,2
          COUPB(I,J) =
     &    -SINT*XMZ**2*2D0*SQR/174.1D0/6D0*SINBPA*(DELTA(I,J) +
     &    (3D0 - 4D0*SINT)/2D0/SINT*B(1,I)*B(1,J))
     &    +RMBOT**2/174.1D0**2*VP/COSB*SA*DELTA(I,J)
     &    +RMBOT/VP/COSB*(AB*SA + XMU*CA)*(B(1,I)*B(2,J) +
     &    B(1,J)*B(2,I))
  160   CONTINUE
  170 CONTINUE
 
      PRUN = XMH
      EPS = 1D-4*PRUN
      ITER = 0
  180 ITER = ITER + 1
      DO 230  I3 = 1,3
 
        PR(I3)=PRUN+(I3-2)*EPS/2
        P2=PR(I3)**2
        POLT = 0D0
        DO 200 I = 1,2
          DO 190 J = 1,2
            POLT = POLT + COUPT(I,J)**2*3D0*
     &      PYFINT(P2,SSTOP2(I),SSTOP2(J))/16D0/PI**2
  190     CONTINUE
  200   CONTINUE
 
        POLB = 0D0
        DO 220 I = 1,2
          DO 210 J = 1,2
            POLB = POLB + COUPB(I,J)**2*3D0*
     &      PYFINT(P2,SSBOT2(I),SSBOT2(J))/16D0/PI**2
  210     CONTINUE
  220   CONTINUE
C        RXMT2 = RXMT**2
        XMT2=XMT**2
 
        POLTT =
     &  3D0*RXMT**2/8D0/PI**2/  V  **2*
     &  CA**2/SINB**2 *
     &  (-2D0*XMT**2+0.5D0*P2)*
     &  PYFINT(P2,XMT2,XMT2)
 
        POL = POLT + POLB + POLTT
        POLAR(I3) = P2 - XMH**2 - POL
  230 CONTINUE
      DERIV = (POLAR(3)-POLAR(1))/EPS
      DRUN = - POLAR(2)/DERIV
      PRUN = PRUN + DRUN
      P2 = PRUN**2
      IF( ABS(DRUN) .LT. 1D-4 .OR.ITER.GT.500) GOTO 240
      GOTO 180
  240 CONTINUE
 
      XMHP = DSQRT(P2)
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...END OF LIGHT HIGGS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
  250 IF(IHIGGS.EQ.1) GOTO 490
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C... STARTING OF HEAVY HIGGS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      DO 270 I = 1,2
        DO 260 J = 1,2
          HCOUPT(I,J) =
     &    -SINT*XMZ**2*2D0*SQR/174.1D0/3D0*COSBPA*(DELTA(I,J) +
     &    (3D0 - 8D0*SINT)/4D0/SINT*T(1,I)*T(1,J))
     &    -RXMT**2/174.1D0**2*VP/SINB*SA*DELTA(I,J)
     &    -RXMT/VP/SINB*(AT*SA - XMU*CA)*(T(1,I)*T(2,J) +
     &    T(1,J)*T(2,I))
  260   CONTINUE
  270 CONTINUE
 
      DO 290 I = 1,2
        DO 280 J = 1,2
          HCOUPB(I,J) =
     &    SINT*XMZ**2*2D0*SQR/174.1D0/6D0*COSBPA*(DELTA(I,J) +
     &    (3D0 - 4D0*SINT)/2D0/SINT*B(1,I)*B(1,J))
     &    -RMBOT**2/174.1D0**2*VP/COSB*CA*DELTA(I,J)
     &    -RMBOT/VP/COSB*(AB*CA - XMU*SA)*(B(1,I)*B(2,J) +
     &    B(1,J)*B(2,I))
          HCOUPB(I,J)=0D0
  280   CONTINUE
  290 CONTINUE
 
      PRUN = HM
      EPS = 1D-4*PRUN
      ITER = 0
  300 ITER = ITER + 1
      DO 350 I3 = 1,3
        PR(I3)=PRUN+(I3-2)*EPS/2
        HP2=PR(I3)**2
 
        HPOLT = 0D0
        DO 320 I = 1,2
          DO 310 J = 1,2
            HPOLT = HPOLT + HCOUPT(I,J)**2*3D0*
     &      PYFINT(HP2,SSTOP2(I),SSTOP2(J))/16D0/PI**2
  310     CONTINUE
  320   CONTINUE
 
        HPOLB = 0D0
        DO 340 I = 1,2
          DO 330 J = 1,2
            HPOLB = HPOLB + HCOUPB(I,J)**2*3D0*
     &      PYFINT(HP2,SSBOT2(I),SSBOT2(J))/16D0/PI**2
  330     CONTINUE
  340   CONTINUE
 
C        RXMT2 = RXMT**2
        XMT2  = XMT**2
 
        HPOLTT =
     &  3D0*RXMT**2/8D0/PI**2/  V  **2*
     &  SA**2/SINB**2 *
     &  (-2D0*XMT**2+0.5D0*HP2)*
     &  PYFINT(HP2,XMT2,XMT2)
 
        HPOL = HPOLT + HPOLB + HPOLTT
        POLAR(I3) =HP2-HM**2-HPOL
  350 CONTINUE
      DERIV = (POLAR(3)-POLAR(1))/EPS
      DRUN = - POLAR(2)/DERIV
      PRUN = PRUN + DRUN
      HP2 = PRUN**2
      IF( ABS(DRUN) .LT. 1D-4 .OR.ITER.GT.500) GOTO 360
      GOTO 300
  360 CONTINUE
 
 
  370 CONTINUE
      HMP = HP2**0.5D0
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C... END OF HEAVY HIGGS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      IF(IHIGGS.EQ.2) GOTO 490
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...BEGINNING OF PSEUDOSCALAR HIGGS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      DO 390 I = 1,2
        DO 380 J = 1,2
          ACOUPT(I,J) =
     &    -RXMT/VP/SINB*(AT*COSB + XMU*SINB)*
     &    (T(1,I)*T(2,J) -T(1,J)*T(2,I))
  380   CONTINUE
  390 CONTINUE
      DO 410 I = 1,2
        DO 400 J = 1,2
          ACOUPB(I,J) =
     &    RMBOT/VP/COSB*(AB*SINB + XMU*COSB)*
     &    (B(1,I)*B(2,J) -B(1,J)*B(2,I))
  400   CONTINUE
  410 CONTINUE
 
      PRUN = XMA
      EPS = 1D-4*PRUN
      ITER = 0
  420 ITER = ITER + 1
      DO 470 I3 = 1,3
        PR(I3)=PRUN+(I3-2)*EPS/2
        AP2=PR(I3)**2
        APOLT = 0D0
        DO 440 I = 1,2
          DO 430 J = 1,2
            APOLT = APOLT + ACOUPT(I,J)**2*3D0*
     &      PYFINT(AP2,SSTOP2(I),SSTOP2(J))/16D0/PI**2
  430     CONTINUE
  440   CONTINUE
        APOLB = 0D0
        DO 460 I = 1,2
          DO 450 J = 1,2
            APOLB = APOLB + ACOUPB(I,J)**2*3D0*
     &      PYFINT(AP2,SSBOT2(I),SSBOT2(J))/16D0/PI**2
  450     CONTINUE
  460   CONTINUE
C        RXMT2 = RXMT**2
        XMT2=XMT**2
        APOLTT =
     &  3D0*RXMT**2/8D0/PI**2/  V  **2*
     &  COSB**2/SINB**2 *
     &  (-0.5D0*AP2)*
     &  PYFINT(AP2,XMT2,XMT2)
        APOL = APOLT + APOLB + APOLTT
        POLAR(I3) = AP2 - XMA**2 -APOL
  470 CONTINUE
      DERIV = (POLAR(3)-POLAR(1))/EPS
      DRUN = - POLAR(2)/DERIV
      PRUN = PRUN + DRUN
      AP2 = PRUN**2
      IF( ABS(DRUN) .LT. 1D-4 .OR.ITER.GT.500) GOTO 480
      GOTO 420
  480 CONTINUE
 
      AMP = DSQRT(AP2)
 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C...END OF PSEUDOSCALAR HIGGS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 
      IF(IHIGGS.EQ.3) GOTO 490
 
  490 CONTINUE
      RETURN
  500 CONTINUE
      WRITE(MSTU(11),*) ' EXITING IN PYPOLE '
      WRITE(MSTU(11),*) ' XMST11,XMST22 = ',XMST11,XMST22
      WRITE(MSTU(11),*) ' XMSB11,XMSB22 = ',XMSB11,XMSB22
      WRITE(MSTU(11),*) ' STOP22,SBOT22 = ',STOP22,SBOT22
      CALL PYSTOP(107)
      END
