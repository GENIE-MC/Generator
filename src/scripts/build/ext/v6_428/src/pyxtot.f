 
C*********************************************************************
 
C...PYXTOT
C...Parametrizes total, elastic and diffractive cross-sections
C...for different energies and beams. Donnachie-Landshoff for
C...total and Schuler-Sjostrand for elastic and diffractive.
C...Process code IPROC:
C...=  1 : p + p;
C...=  2 : pbar + p;
C...=  3 : pi+ + p;
C...=  4 : pi- + p;
C...=  5 : pi0 + p;
C...=  6 : phi + p;
C...=  7 : J/psi + p;
C...= 11 : rho + rho;
C...= 12 : rho + phi;
C...= 13 : rho + J/psi;
C...= 14 : phi + phi;
C...= 15 : phi + J/psi;
C...= 16 : J/psi + J/psi;
C...= 21 : gamma + p (DL);
C...= 22 : gamma + p (VDM).
C...= 23 : gamma + pi (DL);
C...= 24 : gamma + pi (VDM);
C...= 25 : gamma + gamma (DL);
C...= 26 : gamma + gamma (VDM).
 
      SUBROUTINE PYXTOT
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
      COMMON/PYINT7/SIGT(0:6,0:6,0:5)
      SAVE /PYDAT1/,/PYDAT2/,/PYPARS/,/PYINT1/,/PYINT5/,/PYINT7/
C...Local arrays.
      DIMENSION NPROC(30),XPAR(30),YPAR(30),IHADA(20),IHADB(20),
     &PMHAD(4),BHAD(4),BETP(4),IFITSD(20),IFITDD(20),CEFFS(10,8),
     &CEFFD(10,9),SIGTMP(6,0:5)
 
C...Common constants.
      DATA EPS/0.0808D0/, ETA/-0.4525D0/, ALP/0.25D0/, CRES/2D0/,
     &PMRC/1.062D0/, SMP/0.880D0/, FACEL/0.0511D0/, FACSD/0.0336D0/,
     &FACDD/0.0084D0/
 
C...Number of multiple processes to be evaluated (= 0 : undefined).
      DATA NPROC/7*1,3*0,6*1,4*0,4*3,2*6,4*0/
C...X and Y parameters of sigmatot = X * s**epsilon + Y * s**(-eta).
      DATA XPAR/2*21.70D0,3*13.63D0,10.01D0,0.970D0,3*0D0,
     &8.56D0,6.29D0,0.609D0,4.62D0,0.447D0,0.0434D0,4*0D0,
     &0.0677D0,0.0534D0,0.0425D0,0.0335D0,2.11D-4,1.31D-4,4*0D0/
      DATA YPAR/
     &56.08D0,98.39D0,27.56D0,36.02D0,31.79D0,-1.51D0,-0.146D0,3*0D0,
     &13.08D0,-0.62D0,-0.060D0,0.030D0,-0.0028D0,0.00028D0,4*0D0,
     &0.129D0,0.115D0,0.081D0,0.072D0,2.15D-4,1.70D-4,4*0D0/
 
C...Beam and target hadron class:
C...= 1 : p/n ; = 2 : pi/rho/omega; = 3 : phi; = 4 : J/psi.
      DATA IHADA/2*1,3*2,3,4,3*0,3*2,2*3,4,4*0/
      DATA IHADB/7*1,3*0,2,3,4,3,2*4,4*0/
C...Characteristic class masses, slope parameters, beta = sqrt(X).
      DATA PMHAD/0.938D0,0.770D0,1.020D0,3.097D0/
      DATA BHAD/2.3D0,1.4D0,1.4D0,0.23D0/
      DATA BETP/4.658D0,2.926D0,2.149D0,0.208D0/
 
C...Fitting constants used in parametrizations of diffractive results.
      DATA IFITSD/2*1,3*2,3,4,3*0,5,6,7,8,9,10,4*0/
      DATA IFITDD/2*1,3*2,3,4,3*0,5,6,7,8,9,10,4*0/
      DATA ((CEFFS(J1,J2),J2=1,8),J1=1,10)/
     &0.213D0, 0.0D0, -0.47D0, 150D0, 0.213D0, 0.0D0, -0.47D0, 150D0,
     &0.213D0, 0.0D0, -0.47D0, 150D0, 0.267D0, 0.0D0, -0.47D0, 100D0,
     &0.213D0, 0.0D0, -0.47D0, 150D0, 0.232D0, 0.0D0, -0.47D0, 110D0,
     &0.213D0, 7.0D0, -0.55D0, 800D0, 0.115D0, 0.0D0, -0.47D0, 110D0,
     &0.267D0, 0.0D0, -0.46D0,  75D0, 0.267D0, 0.0D0, -0.46D0,  75D0,
     &0.232D0, 0.0D0, -0.46D0,  85D0, 0.267D0, 0.0D0, -0.48D0, 100D0,
     &0.115D0, 0.0D0, -0.50D0,  90D0, 0.267D0, 6.0D0, -0.56D0, 420D0,
     &0.232D0, 0.0D0, -0.48D0, 110D0, 0.232D0, 0.0D0, -0.48D0, 110D0,
     &0.115D0, 0.0D0, -0.52D0, 120D0, 0.232D0, 6.0D0, -0.56D0, 470D0,
     &0.115D0, 5.5D0, -0.58D0, 570D0, 0.115D0, 5.5D0, -0.58D0, 570D0/
      DATA ((CEFFD(J1,J2),J2=1,9),J1=1,10)/
     &3.11D0, -7.34D0,  9.71D0, 0.068D0, -0.42D0,  1.31D0,
     &-1.37D0,  35.0D0,  118D0,  3.11D0, -7.10D0,  10.6D0,
     &0.073D0, -0.41D0, 1.17D0, -1.41D0,  31.6D0,   95D0,
     &3.12D0, -7.43D0,  9.21D0, 0.067D0, -0.44D0,  1.41D0,
     &-1.35D0,  36.5D0,  132D0,  3.13D0, -8.18D0, -4.20D0,
     &0.056D0, -0.71D0, 3.12D0, -1.12D0,  55.2D0, 1298D0,
     &3.11D0, -6.90D0,  11.4D0, 0.078D0, -0.40D0,  1.05D0,
     &-1.40D0,  28.4D0,   78D0,  3.11D0, -7.13D0,  10.0D0,
     &0.071D0, -0.41D0, 1.23D0, -1.34D0,  33.1D0,  105D0,
     &3.12D0, -7.90D0, -1.49D0, 0.054D0, -0.64D0,  2.72D0,
     &-1.13D0,  53.1D0,  995D0,  3.11D0, -7.39D0,  8.22D0,
     &0.065D0, -0.44D0, 1.45D0, -1.36D0,  38.1D0,  148D0,
     &3.18D0, -8.95D0, -3.37D0, 0.057D0, -0.76D0,  3.32D0,
     &-1.12D0,  55.6D0, 1472D0,  4.18D0, -29.2D0,  56.2D0,
     &0.074D0, -1.36D0, 6.67D0, -1.14D0, 116.2D0, 6532D0/
 
C...Parameters. Combinations of the energy.
      AEM=PARU(101)
      PMTH=PARP(102)
      S=VINT(2)
      SRT=VINT(1)
      SEPS=S**EPS
      SETA=S**ETA
      SLOG=LOG(S)
 
C...Ratio of gamma/pi (for rescaling in parton distributions).
      VINT(281)=(XPAR(22)*SEPS+YPAR(22)*SETA)/
     &(XPAR(5)*SEPS+YPAR(5)*SETA)
      VINT(317)=1D0
      IF(MINT(50).NE.1) RETURN
 
C...Order flavours of incoming particles: KF1 < KF2.
      IF(IABS(MINT(11)).LE.IABS(MINT(12))) THEN
        KF1=IABS(MINT(11))
        KF2=IABS(MINT(12))
        IORD=1
      ELSE
        KF1=IABS(MINT(12))
        KF2=IABS(MINT(11))
        IORD=2
      ENDIF
      ISGN12=ISIGN(1,MINT(11)*MINT(12))
 
C...Find process number (for lookup tables).
      IF(KF1.GT.1000) THEN
        IPROC=1
        IF(ISGN12.LT.0) IPROC=2
      ELSEIF(KF1.GT.100.AND.KF2.GT.1000) THEN
        IPROC=3
        IF(ISGN12.LT.0) IPROC=4
        IF(KF1.EQ.111) IPROC=5
      ELSEIF(KF1.GT.100) THEN
        IPROC=11
      ELSEIF(KF2.GT.1000) THEN
        IPROC=21
        IF(MINT(123).EQ.2.OR.MINT(123).EQ.3) IPROC=22
      ELSEIF(KF2.GT.100) THEN
        IPROC=23
        IF(MINT(123).EQ.2.OR.MINT(123).EQ.3) IPROC=24
      ELSE
        IPROC=25
        IF(MINT(123).EQ.2.OR.MINT(123).EQ.3.OR.MINT(123).EQ.7) IPROC=26
      ENDIF
 
C... Number of multiple processes to be stored; beam/target side.
      NPR=NPROC(IPROC)
      MINT(101)=1
      MINT(102)=1
      IF(NPR.EQ.3) THEN
        MINT(100+IORD)=4
      ELSEIF(NPR.EQ.6) THEN
        MINT(101)=4
        MINT(102)=4
      ENDIF
      N1=0
      IF(MINT(101).EQ.4) N1=4
      N2=0
      IF(MINT(102).EQ.4) N2=4
 
C...Do not do any more for user-set or undefined cross-sections.
      IF(MSTP(31).LE.0) RETURN
      IF(NPR.EQ.0) CALL PYERRM(26,
     &'(PYXTOT:) cross section for this process not yet implemented')
 
C...Parameters. Combinations of the energy.
      AEM=PARU(101)
      PMTH=PARP(102)
      S=VINT(2)
      SRT=VINT(1)
      SEPS=S**EPS
      SETA=S**ETA
      SLOG=LOG(S)
 
C...Loop over multiple processes (for VDM).
      DO 110 I=1,NPR
        IF(NPR.EQ.1) THEN
          IPR=IPROC
        ELSEIF(NPR.EQ.3) THEN
          IPR=I+4
          IF(KF2.LT.1000) IPR=I+10
        ELSEIF(NPR.EQ.6) THEN
          IPR=I+10
        ENDIF
 
C...Evaluate hadron species, mass, slope contribution and fit number.
        IHA=IHADA(IPR)
        IHB=IHADB(IPR)
        PMA=PMHAD(IHA)
        PMB=PMHAD(IHB)
        BHA=BHAD(IHA)
        BHB=BHAD(IHB)
        ISD=IFITSD(IPR)
        IDD=IFITDD(IPR)
 
C...Skip if energy too low relative to masses.
        DO 100 J=0,5
          SIGTMP(I,J)=0D0
  100   CONTINUE
        IF(SRT.LT.PMA+PMB+PARP(104)) GOTO 110
 
C...Total cross-section. Elastic slope parameter and cross-section.
        SIGTMP(I,0)=XPAR(IPR)*SEPS+YPAR(IPR)*SETA
        BEL=2D0*BHA+2D0*BHB+4D0*SEPS-4.2D0
        SIGTMP(I,1)=FACEL*SIGTMP(I,0)**2/BEL
 
C...Diffractive scattering A + B -> X + B.
        BSD=2D0*BHB
        SQML=(PMA+PMTH)**2
        SQMU=S*CEFFS(ISD,1)+CEFFS(ISD,2)
        SUM1=LOG((BSD+2D0*ALP*LOG(S/SQML))/
     &  (BSD+2D0*ALP*LOG(S/SQMU)))/(2D0*ALP)
        BXB=CEFFS(ISD,3)+CEFFS(ISD,4)/S
        SUM2=CRES*LOG(1D0+((PMA+PMRC)/(PMA+PMTH))**2)/
     &  (BSD+2D0*ALP*LOG(S/((PMA+PMTH)*(PMA+PMRC)))+BXB)
        SIGTMP(I,2)=FACSD*XPAR(IPR)*BETP(IHB)*MAX(0D0,SUM1+SUM2)
 
C...Diffractive scattering A + B -> A + X.
        BSD=2D0*BHA
        SQML=(PMB+PMTH)**2
        SQMU=S*CEFFS(ISD,5)+CEFFS(ISD,6)
        SUM1=LOG((BSD+2D0*ALP*LOG(S/SQML))/
     &  (BSD+2D0*ALP*LOG(S/SQMU)))/(2D0*ALP)
        BAX=CEFFS(ISD,7)+CEFFS(ISD,8)/S
        SUM2=CRES*LOG(1D0+((PMB+PMRC)/(PMB+PMTH))**2)/
     &  (BSD+2D0*ALP*LOG(S/((PMB+PMTH)*(PMB+PMRC)))+BAX)
        SIGTMP(I,3)=FACSD*XPAR(IPR)*BETP(IHA)*MAX(0D0,SUM1+SUM2)
 
C...Order single diffractive correctly.
        IF(IORD.EQ.2) THEN
          SIGSAV=SIGTMP(I,2)
          SIGTMP(I,2)=SIGTMP(I,3)
          SIGTMP(I,3)=SIGSAV
        ENDIF
 
C...Double diffractive scattering A + B -> X1 + X2.
        YEFF=LOG(S*SMP/((PMA+PMTH)*(PMB+PMTH))**2)
        DEFF=CEFFD(IDD,1)+CEFFD(IDD,2)/SLOG+CEFFD(IDD,3)/SLOG**2
        SUM1=(DEFF+YEFF*(LOG(MAX(1D-10,YEFF/DEFF))-1D0))/(2D0*ALP)
        IF(YEFF.LE.0) SUM1=0D0
        SQMU=S*(CEFFD(IDD,4)+CEFFD(IDD,5)/SLOG+CEFFD(IDD,6)/SLOG**2)
        SLUP=LOG(MAX(1.1D0,S/(ALP*(PMA+PMTH)**2*(PMB+PMTH)*(PMB+PMRC))))
        SLDN=LOG(MAX(1.1D0,S/(ALP*SQMU*(PMB+PMTH)*(PMB+PMRC))))
        SUM2=CRES*LOG(1D0+((PMB+PMRC)/(PMB+PMTH))**2)*LOG(SLUP/SLDN)/
     &  (2D0*ALP)
        SLUP=LOG(MAX(1.1D0,S/(ALP*(PMB+PMTH)**2*(PMA+PMTH)*(PMA+PMRC))))
        SLDN=LOG(MAX(1.1D0,S/(ALP*SQMU*(PMA+PMTH)*(PMA+PMRC))))
        SUM3=CRES*LOG(1D0+((PMA+PMRC)/(PMA+PMTH))**2)*LOG(SLUP/SLDN)/
     &  (2D0*ALP)
        BXX=CEFFD(IDD,7)+CEFFD(IDD,8)/SRT+CEFFD(IDD,9)/S
        SLRR=LOG(S/(ALP*(PMA+PMTH)*(PMA+PMRC)*(PMB+PMTH)*(PMB+PMRC)))
        SUM4=CRES**2*LOG(1D0+((PMA+PMRC)/(PMA+PMTH))**2)*
     &  LOG(1D0+((PMB+PMRC)/(PMB+PMTH))**2)/MAX(0.1D0,2D0*ALP*SLRR+BXX)
        SIGTMP(I,4)=FACDD*XPAR(IPR)*MAX(0D0,SUM1+SUM2+SUM3+SUM4)
 
C...Non-diffractive by unitarity.
        SIGTMP(I,5)=SIGTMP(I,0)-SIGTMP(I,1)-SIGTMP(I,2)-SIGTMP(I,3)-
     &  SIGTMP(I,4)
  110 CONTINUE
 
C...Put temporary results in output array: only one process.
      IF(MINT(101).EQ.1.AND.MINT(102).EQ.1) THEN
        DO 120 J=0,5
          SIGT(0,0,J)=SIGTMP(1,J)
  120   CONTINUE
 
C...Beam multiple processes.
      ELSEIF(MINT(101).EQ.4.AND.MINT(102).EQ.1) THEN
        IF(MINT(107).EQ.2) THEN
          VINT(317)=(PMHAD(2)**2/(PMHAD(2)**2+VINT(307)))**2
        ELSE
          VINT(317)=16D0*PARP(15)**2*VINT(154)**2/
     &    ((4D0*PARP(15)**2+VINT(307))*(4D0*VINT(154)**2+VINT(307)))
        ENDIF
        IF(MSTP(20).GT.0) THEN
          VINT(317)=VINT(317)*(VINT(2)/(VINT(2)+VINT(307)))**MSTP(20)
        ENDIF
        DO 140 I=1,4
          IF(MINT(107).EQ.2) THEN
            CONV=(AEM/PARP(160+I))*VINT(317)
          ELSEIF(VINT(154).GT.PARP(15)) THEN
            CONV=(AEM/PARU(1))*(KCHG(I,1)/3D0)**2*PARP(18)**2*
     &      (1D0/PARP(15)**2-1D0/VINT(154)**2)*VINT(317)
          ELSE
            CONV=0D0
          ENDIF
          I1=MAX(1,I-1)
          DO 130 J=0,5
            SIGT(I,0,J)=CONV*SIGTMP(I1,J)
  130     CONTINUE
  140   CONTINUE
        DO 150 J=0,5
          SIGT(0,0,J)=SIGT(1,0,J)+SIGT(2,0,J)+SIGT(3,0,J)+SIGT(4,0,J)
  150   CONTINUE
 
C...Target multiple processes.
      ELSEIF(MINT(101).EQ.1.AND.MINT(102).EQ.4) THEN
        IF(MINT(108).EQ.2) THEN
          VINT(317)=(PMHAD(2)**2/(PMHAD(2)**2+VINT(308)))**2
        ELSE
          VINT(317)=16D0*PARP(15)**2*VINT(154)**2/
     &    ((4D0*PARP(15)**2+VINT(308))*(4D0*VINT(154)**2+VINT(308)))
        ENDIF
        IF(MSTP(20).GT.0) THEN
          VINT(317)=VINT(317)*(VINT(2)/(VINT(2)+VINT(308)))**MSTP(20)
        ENDIF
        DO 170 I=1,4
          IF(MINT(108).EQ.2) THEN
            CONV=(AEM/PARP(160+I))*VINT(317)
          ELSEIF(VINT(154).GT.PARP(15)) THEN
            CONV=(AEM/PARU(1))*(KCHG(I,1)/3D0)**2*PARP(18)**2*
     &      (1D0/PARP(15)**2-1D0/VINT(154)**2)*VINT(317)
          ELSE
            CONV=0D0
          ENDIF
          IV=MAX(1,I-1)
          DO 160 J=0,5
            SIGT(0,I,J)=CONV*SIGTMP(IV,J)
  160     CONTINUE
  170   CONTINUE
        DO 180 J=0,5
          SIGT(0,0,J)=SIGT(0,1,J)+SIGT(0,2,J)+SIGT(0,3,J)+SIGT(0,4,J)
  180   CONTINUE
 
C...Both beam and target multiple processes.
      ELSE
        IF(MINT(107).EQ.2) THEN
          VINT(317)=(PMHAD(2)**2/(PMHAD(2)**2+VINT(307)))**2
        ELSE
          VINT(317)=16D0*PARP(15)**2*VINT(154)**2/
     &    ((4D0*PARP(15)**2+VINT(307))*(4D0*VINT(154)**2+VINT(307)))
        ENDIF
        IF(MINT(108).EQ.2) THEN
          VINT(317)=VINT(317)*(PMHAD(2)**2/(PMHAD(2)**2+VINT(308)))**2
        ELSE
          VINT(317)=VINT(317)*16D0*PARP(15)**2*VINT(154)**2/
     &    ((4D0*PARP(15)**2+VINT(308))*(4D0*VINT(154)**2+VINT(308)))
        ENDIF
        IF(MSTP(20).GT.0) THEN
          VINT(317)=VINT(317)*(VINT(2)/(VINT(2)+VINT(307)+
     &    VINT(308)))**MSTP(20)
        ENDIF
        DO 210 I1=1,4
          DO 200 I2=1,4
            IF(MINT(107).EQ.2) THEN
              CONV=(AEM/PARP(160+I1))*VINT(317)
            ELSEIF(VINT(154).GT.PARP(15)) THEN
              CONV=(AEM/PARU(1))*(KCHG(I1,1)/3D0)**2*PARP(18)**2*
     &        (1D0/PARP(15)**2-1D0/VINT(154)**2)*VINT(317)
            ELSE
              CONV=0D0
            ENDIF
            IF(MINT(108).EQ.2) THEN
              CONV=CONV*(AEM/PARP(160+I2))
            ELSEIF(VINT(154).GT.PARP(15)) THEN
              CONV=CONV*(AEM/PARU(1))*(KCHG(I2,1)/3D0)**2*PARP(18)**2*
     &        (1D0/PARP(15)**2-1D0/VINT(154)**2)
            ELSE
              CONV=0D0
            ENDIF
            IF(I1.LE.2) THEN
              IV=MAX(1,I2-1)
            ELSEIF(I2.LE.2) THEN
              IV=MAX(1,I1-1)
            ELSEIF(I1.EQ.I2) THEN
              IV=2*I1-2
            ELSE
              IV=5
            ENDIF
            DO 190 J=0,5
              JV=J
              IF(I2.GT.I1.AND.(J.EQ.2.OR.J.EQ.3)) JV=5-J
              SIGT(I1,I2,J)=CONV*SIGTMP(IV,JV)
  190       CONTINUE
  200     CONTINUE
  210   CONTINUE
        DO 230 J=0,5
          DO 220 I=1,4
            SIGT(I,0,J)=SIGT(I,1,J)+SIGT(I,2,J)+SIGT(I,3,J)+SIGT(I,4,J)
            SIGT(0,I,J)=SIGT(1,I,J)+SIGT(2,I,J)+SIGT(3,I,J)+SIGT(4,I,J)
  220     CONTINUE
          SIGT(0,0,J)=SIGT(1,0,J)+SIGT(2,0,J)+SIGT(3,0,J)+SIGT(4,0,J)
  230   CONTINUE
      ENDIF
 
C...Scale up uniformly for Donnachie-Landshoff parametrization.
      IF(IPROC.EQ.21.OR.IPROC.EQ.23.OR.IPROC.EQ.25) THEN
        RFAC=(XPAR(IPROC)*SEPS+YPAR(IPROC)*SETA)/SIGT(0,0,0)
        DO 260 I1=0,N1
          DO 250 I2=0,N2
            DO 240 J=0,5
              SIGT(I1,I2,J)=RFAC*SIGT(I1,I2,J)
  240       CONTINUE
  250     CONTINUE
  260   CONTINUE
      ENDIF
 
      RETURN
      END
