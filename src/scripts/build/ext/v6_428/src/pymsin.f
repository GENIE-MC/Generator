 
C*********************************************************************
 
C...PYMSIN
C...Initializes supersymmetry: finds sparticle masses and
C...branching ratios and stores this information.
C...AUTHOR: STEPHEN MRENNA
C...Author: P. Skands (SLHA + RPV + ISASUSY Interface, NMSSM)
 
      SUBROUTINE PYMSIN
 
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
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYDAT4/CHAF(500,2)
      CHARACTER CHAF*16
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT4/MWID(500),WIDS(500,5)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYMSRV/RVLAM(3,3,3), RVLAMP(3,3,3), RVLAMB(3,3,3)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
      COMMON/PYHTRI/HHH(7)
      COMMON/PYQNUM/NQNUM,NQDUM,KQNUM(500,0:9)
      SAVE /PYDAT1/,/PYDAT2/,/PYDAT3/,/PYDAT4/,/PYPARS/,/PYINT4/,
     &/PYMSSM/,/PYMSRV/,/PYSSMT/
 
C...Local variables.
      DOUBLE PRECISION ALFA,BETA
      DOUBLE PRECISION TANB,AL,BE,COSA,COSB,SINA,SINB,XW
      INTEGER I,J,J1,I1,K1
      INTEGER KC,LKNT,IDLAM(400,3)
      DOUBLE PRECISION XLAM(0:400)
      DOUBLE PRECISION WDTP(0:400),WDTE(0:400,0:5)
      DOUBLE PRECISION XARG,COS2B,XMW2,XMZ2
      DOUBLE PRECISION DELM,XMDIF
      DOUBLE PRECISION DX,DY,DS,DMU2,DMA2,DQ2,DU2,DD2,DL2,DE2,DHU2,DHD2
      DOUBLE PRECISION ARG,SGNMU,R
      INTEGER IMSSM
      INTEGER IRPRTY
      INTEGER KFSUSY(50),MWIDSU(36),MDCYSU(36)
      SAVE MWIDSU,MDCYSU
      DATA KFSUSY/
     &1000001,2000001,1000002,2000002,1000003,2000003,
     &1000004,2000004,1000005,2000005,1000006,2000006,
     &1000011,2000011,1000012,2000012,1000013,2000013,
     &1000014,2000014,1000015,2000015,1000016,2000016,
     &1000021,1000022,1000023,1000025,1000035,1000024,
     &1000037,1000039,     25,     35,     36,     37,
     &      6,     24,     45,     46,1000045, 9*0/
      DATA INIT/0/
 
C...Automatically read QNUMBERS, MASS, and DECAY tables      
      IF (IMSS(21).NE.0.OR.MSTP(161).NE.0) THEN
        NQNUM=0
        CALL PYSLHA(0,0,IFAIL)
        CALL PYSLHA(5,0,IFAIL)
      ENDIF
      IF (IMSS(22).NE.0.OR.MSTP(161).NE.0) CALL PYSLHA(2,0,IFAIL)

C...Do nothing further if SUSY not requested
      IMSSM=IMSS(1)
      IF(IMSSM.EQ.0) RETURN
      
C...Save copy of MWID(KC) and MDCY(KC,1) values before
C...they are set to zero for the LSP.
      IF(INIT.EQ.0) THEN
        INIT=1
        DO 100 I=1,36
          KF=KFSUSY(I)
          KC=PYCOMP(KF)
          MWIDSU(I)=MWID(KC)
          MDCYSU(I)=MDCY(KC,1)
  100   CONTINUE
      ENDIF
 
C...Restore MWID(KC) and MDCY(KC,1) values previously zeroed for LSP.
      DO 110 I=1,36
        KF=KFSUSY(I)
        KC=PYCOMP(KF)
        IF(MDCY(KC,1).EQ.0.AND.MDCYSU(I).NE.0) THEN
          MWID(KC)=MWIDSU(I)
          MDCY(KC,1)=MDCYSU(I)
        ENDIF
  110 CONTINUE
 
C...First part of routine: set masses and couplings.
 
C...Reset mixing values in sfermion sector to pure left/right.
      DO 120 I=1,16
        SFMIX(I,1)=1D0
        SFMIX(I,4)=1D0
        SFMIX(I,2)=0D0
        SFMIX(I,3)=0D0
  120 CONTINUE
 
C...Add NMSSM states if NMSSM switched on, and change old names.
      IF (IMSS(13).NE.0.AND.PYCOMP(1000045).EQ.0) THEN
C...  Switch on NMSSM
        WRITE(MSTU(11),*) '(PYMSIN:) switching on NMSSM'
 
        KFN=25
        KCN=KFN
        CHAF(KCN,1)='h_10'
        CHAF(KCN,2)=' '
 
        KFN=35
        KCN=KFN
        CHAF(KCN,1)='h_20'
        CHAF(KCN,2)=' '
 
        KFN=45
        KCN=KFN
        CHAF(KCN,1)='h_30'
        CHAF(KCN,2)=' '
 
        KFN=36
        KCN=KFN
        CHAF(KCN,1)='A_10'
        CHAF(KCN,2)=' '
 
        KFN=46
        KCN=KFN
        CHAF(KCN,1)='A_20'
        CHAF(KCN,2)=' '
 
        KFN=1000045
        KCN=PYCOMP(KFN)
        IF (KCN.EQ.0) THEN
          DO 123 KCT=100,MSTU(6)
            IF(KCHG(KCT,4).GT.100) KCN=KCT
 123      CONTINUE
          KCN=KCN+1
          KCHG(KCN,4)=KFN
          MSTU(20)=0
        ENDIF
C...  Set stable for now
        PMAS(KCN,2)=1D-6
        MWID(KCN)=0
        MDCY(KCN,1)=0
        MDCY(KCN,2)=0
        MDCY(KCN,3)=0
        CHAF(KCN,1)='~chi_50'
        CHAF(KCN,2)=' '
      ENDIF
 
C...Read spectrum from SLHA file.
      IF (IMSSM.EQ.11) THEN
        CALL PYSLHA(1,0,IFAIL)
      ENDIF
 
C...Common couplings.
      TANB=RMSS(5)
      BETA=ATAN(TANB)
      COSB=COS(BETA)
      SINB=TANB*COSB
      COS2B=COS(2D0*BETA)
      ALFA=RMSS(18)
      XMW2=PMAS(24,1)**2
      XMZ2=PMAS(23,1)**2
      XW=PARU(102)
 
C...Define sparticle masses for a general MSSM simulation.
      IF(IMSSM.EQ.1) THEN
        IF(IMSS(9).EQ.0) RMSS(22)=RMSS(9)
        DO 130 I=1,5,2
          KC=PYCOMP(KSUSY1+I)
          PMAS(KC,1)=SQRT(RMSS(8)**2-(2D0*XMW2+XMZ2)*COS2B/6D0)
          KC=PYCOMP(KSUSY2+I)
          PMAS(KC,1)=SQRT(RMSS(9)**2+(XMW2-XMZ2)*COS2B/3D0)
          KC=PYCOMP(KSUSY1+I+1)
          PMAS(KC,1)=SQRT(RMSS(8)**2+(4D0*XMW2-XMZ2)*COS2B/6D0)
          KC=PYCOMP(KSUSY2+I+1)
          PMAS(KC,1)=SQRT(RMSS(22)**2-(XMW2-XMZ2)*COS2B*2D0/3D0)
  130   CONTINUE
        XARG=RMSS(6)**2-PMAS(24,1)**2*ABS(COS(2D0*BETA))
        IF(XARG.LT.0D0) THEN
          WRITE(MSTU(11),*) ' SNEUTRINO MASS IS NEGATIVE'//
     &    ' FROM THE SUM RULE. '
          WRITE(MSTU(11),*) '  TRY A SMALLER VALUE OF TAN(BETA). '
          RETURN
        ELSE
          XARG=SQRT(XARG)
        ENDIF
        DO 140 I=11,15,2
          PMAS(PYCOMP(KSUSY1+I),1)=RMSS(6)
          PMAS(PYCOMP(KSUSY2+I),1)=RMSS(7)
          PMAS(PYCOMP(KSUSY1+I+1),1)=XARG
          PMAS(PYCOMP(KSUSY2+I+1),1)=9999D0
  140   CONTINUE
        IF(IMSS(8).EQ.1) THEN
          RMSS(13)=RMSS(6)
          RMSS(14)=RMSS(7)
        ENDIF
 
C...Alternatively derive masses from SUGRA relations.
      ELSEIF(IMSSM.EQ.2) THEN
        RMSS(36)=RMSS(16)
        CALL PYAPPS
C...Or use ISASUSY
      ELSEIF(IMSSM.EQ.12.OR.IMSSM.EQ.13) THEN
        RMSS(36)=RMSS(16)
        CALL PYSUGI
        ALFA=RMSS(18)
        GOTO 170
      ELSE
        GOTO 170
      ENDIF
 
C...Add in extra D-term contributions.
      IF(IMSS(7).EQ.1) THEN
        R=0.43D0
        DX=RMSS(23)
        DY=RMSS(24)
        DS=RMSS(25)
        WRITE(MSTU(11),*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        WRITE(MSTU(11),*) 'C  NEW DTERMS ADDED TO SCALAR MASSES   '
        WRITE(MSTU(11),*) 'C   IN A U(B-L) THEORY                 '
        WRITE(MSTU(11),*) 'C   DX = ',DX
        WRITE(MSTU(11),*) 'C   DY = ',DY
        WRITE(MSTU(11),*) 'C   DS = ',DS
        WRITE(MSTU(11),*) 'C                                      '
        DY=R*DY-4D0/33D0*(1D0-R)*DX+(1D0-R)/33D0*DS
        WRITE(MSTU(11),*) 'C   DY AT THE WEAK SCALE = ',DY
        WRITE(MSTU(11),*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        DQ2=DY/6D0-DX/3D0-DS/3D0
        DU2=-2D0*DY/3D0-DX/3D0-DS/3D0
        DD2=DY/3D0+DX-2D0*DS/3D0
        DL2=-DY/2D0+DX-2D0*DS/3D0
        DE2=DY-DX/3D0-DS/3D0
        DHU2=DY/2D0+2D0*DX/3D0+2D0*DS/3D0
        DHD2=-DY/2D0-2D0*DX/3D0+DS
        DMU2=(-DY/2D0-2D0/3D0*DX+(COSB**2-2D0*SINB**2/3D0)*DS)
     &  /ABS(COS2B)
        DMA2 = 2D0*DMU2+DHU2+DHD2
        DO 150 I=1,5,2
          KC=PYCOMP(KSUSY1+I)
          PMAS(KC,1)=SQRT(PMAS(KC,1)**2+DQ2)
          KC=PYCOMP(KSUSY2+I)
          PMAS(KC,1)=SQRT(PMAS(KC,1)**2+DD2)
          KC=PYCOMP(KSUSY1+I+1)
          PMAS(KC,1)=SQRT(PMAS(KC,1)**2+DQ2)
          KC=PYCOMP(KSUSY2+I+1)
          PMAS(KC,1)=SQRT(PMAS(KC,1)**2+DU2)
  150   CONTINUE
        DO 160 I=11,15,2
          KC=PYCOMP(KSUSY1+I)
          PMAS(KC,1)=SQRT(PMAS(KC,1)**2+DL2)
          KC=PYCOMP(KSUSY2+I)
          PMAS(KC,1)=SQRT(PMAS(KC,1)**2+DE2)
          KC=PYCOMP(KSUSY1+I+1)
          PMAS(KC,1)=SQRT(PMAS(KC,1)**2+DL2)
  160   CONTINUE
        IF(RMSS(4)**2+DMU2.LT.0D0) THEN
          WRITE(MSTU(11),*) ' MU2 DRIVEN NEGATIVE '
          CALL PYSTOP(104)
        ENDIF
        SGNMU=SIGN(1D0,RMSS(4))
        RMSS(4)=SGNMU*SQRT(RMSS(4)**2+DMU2)
        ARG=RMSS(10)**2*SIGN(1D0,RMSS(10))+DQ2
        RMSS(10)=SIGN(SQRT(ABS(ARG)),ARG)
        ARG=RMSS(11)**2*SIGN(1D0,RMSS(11))+DD2
        RMSS(11)=SIGN(SQRT(ABS(ARG)),ARG)
        ARG=RMSS(12)**2*SIGN(1D0,RMSS(12))+DU2
        RMSS(12)=SIGN(SQRT(ABS(ARG)),ARG)
        ARG=RMSS(13)**2*SIGN(1D0,RMSS(13))+DL2
        RMSS(13)=SIGN(SQRT(ABS(ARG)),ARG)
        ARG=RMSS(14)**2*SIGN(1D0,RMSS(14))+DE2
        RMSS(14)=SIGN(SQRT(ABS(ARG)),ARG)
        IF( RMSS(19)**2 + DMA2 .LE. 50D0 ) THEN
          WRITE(MSTU(11),*) ' MA DRIVEN TOO LOW '
          CALL PYSTOP(104)
        ENDIF
        RMSS(19)=SQRT(RMSS(19)**2+DMA2)
        RMSS(6)=SQRT(RMSS(6)**2+DL2)
        RMSS(7)=SQRT(RMSS(7)**2+DE2)
        WRITE(MSTU(11),*) ' MTL = ',RMSS(10)
        WRITE(MSTU(11),*) ' MBR = ',RMSS(11)
        WRITE(MSTU(11),*) ' MTR = ',RMSS(12)
        WRITE(MSTU(11),*) ' SEL = ',RMSS(6),RMSS(13)
        WRITE(MSTU(11),*) ' SER = ',RMSS(7),RMSS(14)
      ENDIF
 
C...Fix the third generation sfermions.
      CALL PYTHRG
 
C...Fix the neutralino--chargino--gluino sector.
      CALL PYINOM
 
C...Fix the Higgs sector.
      CALL PYHGGM(ALFA)
 
C...Choose the Gunion-Haber convention.
      ALFA=-ALFA
      RMSS(18)=ALFA
 
C...Print information on mass parameters.
      IF(IMSSM.EQ.2.AND.MSTP(122).GT.0) THEN
        WRITE(MSTU(11),*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        WRITE(MSTU(11),*) ' USING APPROXIMATE SUGRA RELATIONS '
        WRITE(MSTU(11),*) ' M0 = ',RMSS(8)
        WRITE(MSTU(11),*) ' M1/2=',RMSS(1)
        WRITE(MSTU(11),*) ' TANB=',RMSS(5)
        WRITE(MSTU(11),*) ' MU = ',RMSS(4)
        WRITE(MSTU(11),*) ' AT = ',RMSS(16)
        WRITE(MSTU(11),*) ' MA = ',RMSS(19)
        WRITE(MSTU(11),*) ' MTOP=',PMAS(6,1)
        WRITE(MSTU(11),*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      ENDIF
      IF(IMSS(20).EQ.1) THEN
        WRITE(MSTU(11),*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        WRITE(MSTU(11),*) ' DEBUG MODE '
        WRITE(MSTU(11),*) ' UMIX = ',UMIX(1,1),UMIX(1,2),
     &  UMIX(2,1),UMIX(2,2)
        WRITE(MSTU(11),*) ' UMIXI = ',UMIXI(1,1),UMIXI(1,2),
     &  UMIXI(2,1),UMIXI(2,2)
        WRITE(MSTU(11),*) ' VMIX = ',VMIX(1,1),VMIX(1,2),
     &  VMIX(2,1),VMIX(2,2)
        WRITE(MSTU(11),*) ' VMIXI = ',VMIXI(1,1),VMIXI(1,2),
     &  VMIXI(2,1),VMIXI(2,2)
        WRITE(MSTU(11),*) ' ZMIX = ',(ZMIX(1,I),I=1,4)
        WRITE(MSTU(11),*) ' ZMIXI = ',(ZMIXI(1,I),I=1,4)
        WRITE(MSTU(11),*) ' ZMIX = ',(ZMIX(2,I),I=1,4)
        WRITE(MSTU(11),*) ' ZMIXI = ',(ZMIXI(2,I),I=1,4)
        WRITE(MSTU(11),*) ' ZMIX = ',(ZMIX(3,I),I=1,4)
        WRITE(MSTU(11),*) ' ZMIXI = ',(ZMIXI(3,I),I=1,4)
        WRITE(MSTU(11),*) ' ZMIX = ',(ZMIX(4,I),I=1,4)
        WRITE(MSTU(11),*) ' ZMIXI = ',(ZMIXI(4,I),I=1,4)
        WRITE(MSTU(11),*) ' ALFA = ',ALFA
        WRITE(MSTU(11),*) ' BETA = ',BETA
        WRITE(MSTU(11),*) ' STOP = ',(SFMIX(6,I),I=1,4)
        WRITE(MSTU(11),*) ' SBOT = ',(SFMIX(5,I),I=1,4)
        WRITE(MSTU(11),*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      ENDIF
 
C...Set up the Higgs couplings - needed here since initialization
C...in PYINRE did not yet occur when PYWIDT is called below.
  170 AL=ALFA
      BE=BETA
      SINA=SIN(AL)
      COSA=COS(AL)
      COSB=COS(BE)
      SINB=TANB*COSB
      SBMA=SIN(BE-AL)
      SAPB=SIN(AL+BE)
      CAPB=COS(AL+BE)
      CBMA=COS(BE-AL)
      C2A=COS(2D0*AL)
      C2B=COSB**2-SINB**2
C...tanb (used for H+)
      PARU(141)=TANB
 
C...Firstly: h
C...Coupling to d-type quarks
      PARU(161)=SINA/COSB
C...Coupling to u-type quarks
      PARU(162)=-COSA/SINB
C...Coupling to leptons
      PARU(163)=PARU(161)
C...Coupling to Z
      PARU(164)=SBMA
C...Coupling to W
      PARU(165)=PARU(164)
 
C...Secondly: H
C...Coupling to d-type quarks
      PARU(171)=-COSA/COSB
C...Coupling to u-type quarks
      PARU(172)=-SINA/SINB
C...Coupling to leptons
      PARU(173)=PARU(171)
C...Coupling to Z
      PARU(174)=CBMA
C...Coupling to W
      PARU(175)=PARU(174)
C...Coupling to h
      IF(IMSS(4).GE.2) THEN
        PARU(176)=COS(2D0*AL)*COS(BE+AL)-2D0*SIN(2D0*AL)*SIN(BE+AL)
      ELSE
        HHH(3)=HHH(3)+HHH(4)+HHH(5)
        PARU(176)=-3D0/HHH(1)*(HHH(1)*SINA**2*COSB*COSA+
     1  HHH(2)*COSA**2*SINB*SINA+HHH(3)*(SINA**3*SINB+COSA**3*COSB-
     2  2D0/3D0*CBMA)-HHH(6)*SINA*(COSB*C2A+COSA*CAPB)+
     3  HHH(7)*COSA*(SINB*C2A+SINA*CAPB))
      ENDIF
C...Coupling to H+
C...Define later
      IF(IMSS(4).GE.2) THEN
        PARU(168)=-SBMA-COS(2D0*BE)*SAPB/2D0/(1D0-XW)
      ELSE
        PARU(168)=1D0/HHH(1)*(HHH(1)*SINB**2*COSB*SINA-
     1 HHH(2)*COSB**2*SINB*COSA-HHH(3)*(SINB**3*COSA-COSB**3*SINA)+
     2 2D0*HHH(5)*SBMA-HHH(6)*SINB*(COSB*SAPB+SINA*C2B)-
     3 HHH(7)*COSB*(COSA*C2B-SINB*SAPB)-(HHH(5)-HHH(4))*SBMA)
      ENDIF
C...Coupling to A
      IF(IMSS(4).GE.2) THEN
        PARU(177)=COS(2D0*BE)*COS(BE+AL)
      ELSE
        PARU(177)=-1D0/HHH(1)*(HHH(1)*SINB**2*COSB*COSA+
     1 HHH(2)*COSB**2*SINB*SINA+HHH(3)*(SINB**3*SINA+COSB**3*COSA)-
     2 2D0*HHH(5)*CBMA-HHH(6)*SINB*(COSB*CAPB+COSA*C2B)+
     3 HHH(7)*COSB*(SINB*CAPB+SINA*C2B))
      ENDIF
C...Coupling to H+
      IF(IMSS(4).GE.2) THEN
        PARU(178)=PARU(177)
      ELSE
        PARU(178)=PARU(177)-(HHH(5)-HHH(4))/HHH(1)*CBMA
      ENDIF
C...Thirdly, A
C...Coupling to d-type quarks
      PARU(181)=TANB
C...Coupling to u-type quarks
      PARU(182)=1D0/PARU(181)
C...Coupling to leptons
      PARU(183)=PARU(181)
      PARU(184)=0D0
      PARU(185)=0D0
C...Coupling to Z h
      PARU(186)=COS(BE-AL)
C...Coupling to Z H
      PARU(187)=SIN(BE-AL)
      PARU(188)=0D0
      PARU(189)=0D0
      PARU(190)=0D0
 
C...Finally: H+
C...Coupling to W h
      PARU(195)=COS(BE-AL)
 
C...Tell that all Higgs couplings have been set.
      MSTP(4)=1
 
C...Set R-Violating couplings.
C...Set lambda couplings to common value or "natural values".
      IF ((IMSS(51).NE.3).AND.(IMSS(51).NE.0)) THEN
        VIR3=1D0/(126D0)**3
        DO 200 IRK=1,3
          DO 190 IRI=1,3
            DO 180 IRJ=1,3
              IF (IRI.NE.IRJ) THEN
                IF (IRI.LT.IRJ) THEN
                  RVLAM(IRI,IRJ,IRK)=RMSS(51)
                  IF (IMSS(51).EQ.2) RVLAM(IRI,IRJ,IRK)=RMSS(51)*
     &              SQRT(PMAS(9+2*IRI,1)*PMAS(9+2*IRJ,1)*
     &              PMAS(9+2*IRK,1)*VIR3)
                ELSE
                  RVLAM(IRI,IRJ,IRK)=-RVLAM(IRJ,IRI,IRK)
                ENDIF
              ELSE
                RVLAM(IRI,IRJ,IRK)=0D0
              ENDIF
  180       CONTINUE
  190     CONTINUE
  200   CONTINUE
      ENDIF
C...Set lambda' couplings to common value or "natural values".
      IF ((IMSS(52).NE.3).AND.(IMSS(52).NE.0)) THEN
        VIR3=1D0/(126D0)**3
        DO 230 IRI=1,3
          DO 220 IRJ=1,3
            DO 210 IRK=1,3
              RVLAMP(IRI,IRJ,IRK)=RMSS(52)
              IF (IMSS(52).EQ.2) RVLAMP(IRI,IRJ,IRK)=RMSS(52)*
     &          SQRT(PMAS(9+2*IRI,1)*0.5D0*(PMAS(2*IRJ,1)+
     &          PMAS(2*IRJ-1,1))*PMAS(2*IRK-1,1)*VIR3)
  210       CONTINUE
  220     CONTINUE
  230   CONTINUE
      ENDIF
C...Set lambda'' couplings to common value or "natural values".
      IF ((IMSS(53).NE.3).AND.(IMSS(53).NE.0)) THEN
        VIR3=1D0/(126D0)**3
        DO 260 IRI=1,3
          DO 250 IRJ=1,3
            DO 240 IRK=1,3
              IF (IRJ.NE.IRK) THEN
                IF (IRJ.LT.IRK) THEN
                  RVLAMB(IRI,IRJ,IRK)=RMSS(53)
                  IF (IMSS(53).EQ.2) RVLAMB(IRI,IRJ,IRK)=
     &              RMSS(53)*SQRT(PMAS(2*IRI,1)*PMAS(2*IRJ-1,1)*
     &              PMAS(2*IRK-1,1)*VIR3)
                ELSE
                  RVLAMB(IRI,IRJ,IRK)=-RVLAMB(IRI,IRK,IRJ)
                ENDIF
              ELSE
                RVLAMB(IRI,IRJ,IRK) = 0D0
              ENDIF
  240       CONTINUE
  250     CONTINUE
  260   CONTINUE
      ENDIF
 
C...Antisymmetrize couplings set by user
      IF (IMSS(51).EQ.3.OR.IMSS(53).EQ.3) THEN
        DO 290 IRI=1,3
          DO 280 IRJ=1,3
            DO 270 IRK=1,3
              IF (RVLAM(IRI,IRJ,IRK).NE.-RVLAM(IRJ,IRI,IRK)) THEN
                RVLAM(IRJ,IRI,IRK)=-RVLAM(IRI,IRJ,IRK)
                IF (IRI.EQ.IRJ) RVLAM(IRI,IRJ,IRK)=0D0
              ENDIF
              IF (RVLAMB(IRI,IRJ,IRK).NE.-RVLAMB(IRI,IRK,IRJ)) THEN
                RVLAMB(IRI,IRK,IRJ)=-RVLAMB(IRI,IRJ,IRK)
                IF (IRJ.EQ.IRK) RVLAMB(IRI,IRJ,IRK)=0D0
              ENDIF
  270       CONTINUE
  280     CONTINUE
  290   CONTINUE
      ENDIF
 
C...Write spectrum to SLHA file
      IF (IMSS(23).NE.0) THEN
	IFAIL=0
        CALL PYSLHA(3,0,IFAIL)
      ENDIF
 
C...Second part of routine: set decay modes and branching ratios.
 
C...Allow chi10 -> gravitino + gamma or not.
      KC=PYCOMP(KSUSY1+39)
      IF( IMSS(11) .NE. 0 ) THEN
        PMAS(KC,1)=RMSS(21)/1D9
        PMAS(KC,2)=0D0
        IRPRTY=0
        WRITE(MSTU(11),*) ' ALLOWING DECAYS TO GRAVITINOS '
      ELSE IF (IMSS(51).GE.1.OR.IMSS(52).GE.1.OR.IMSS(53).GE.1) THEN
        IRPRTY=0
        IF (IMSS(51).GE.1) WRITE(MSTU(11),*)
     &       ' ALLOWING SUSY LLE DECAYS'
        IF (IMSS(52).GE.1) WRITE(MSTU(11),*)
     &       ' ALLOWING SUSY LQD DECAYS'
        IF (IMSS(53).GE.1) WRITE(MSTU(11),*)
     &       ' ALLOWING SUSY UDD DECAYS'
        IF (IMSS(53).GE.1.AND.IMSS(52).GE.1) WRITE(MSTU(11),*)
     &   ' --- Warning: R-Violating couplings possibly',
     &       ' incompatible with proton decay'
      ELSE
        PMAS(KC,1)=9999D0
        IRPRTY=1
      ENDIF
 
C...Loop over sparticle and Higgs species.
      PMCHI1=PMAS(PYCOMP(KSUSY1+22),1)
C...Find the LSP or NLSP for a gravitino LSP
      ILSP=0
      PMLSP=1D20
      DO 300 I=1,36
        KF=KFSUSY(I)
        IF(KF.EQ.1000039) GOTO 300
        KC=PYCOMP(KF)
        IF(PMAS(KC,1).LT.PMLSP) THEN
          ILSP=I
          PMLSP=PMAS(KC,1)
        ENDIF
  300 CONTINUE
      DO 370 I=1,50
        IF (I.GT.39.AND.IMSS(13).NE.1) GOTO 370
        KF=KFSUSY(I)
        IF (KF.EQ.0) GOTO 370
        KC=PYCOMP(KF)
        LKNT=0
 
C...Check if there are any decays listed for this sparticle
C...in a file
        IF (IMSS(22).NE.0.OR.MSTP(161).NE.0) THEN
          IFAIL=0
          CALL PYSLHA(2,KF,IFAIL)
          IF (IFAIL.EQ.0.OR.KF.EQ.6.OR.KF.EQ.24) GOTO 370
        ELSEIF (I.GE.37) THEN
          GOTO 370
        ENDIF
 
C...Sfermion decays.
        IF(I.LE.24) THEN
C...First check to see if sneutrino is lighter than chi10.
          IF((I.EQ.15.OR.I.EQ.19.OR.I.EQ.23).AND.
     &    PMAS(KC,1).LT.PMCHI1) THEN
          ELSE
            CALL PYSFDC(KF,XLAM,IDLAM,LKNT)
          ENDIF
 
C...Gluino decays.
        ELSEIF(I.EQ.25) THEN
          CALL PYGLUI(KF,XLAM,IDLAM,LKNT)
          IF(I.EQ.ILSP.AND.IRPRTY.EQ.1) LKNT=0
 
C...Neutralino decays.
        ELSEIF(I.GE.26.AND.I.LE.29) THEN
          CALL PYNJDC(KF,XLAM,IDLAM,LKNT)
C...chi10 stable or chi10 -> gravitino + gamma.
          IF(I.EQ.26.AND.IRPRTY.EQ.1) THEN
            PMAS(KC,2)=1D-6
            MDCY(KC,1)=0
            MWID(KC)=0
          ENDIF
 
C...Chargino decays.
        ELSEIF(I.GE.30.AND.I.LE.31) THEN
          CALL PYCJDC(KF,XLAM,IDLAM,LKNT)
 
C...Gravitino is stable.
        ELSEIF(I.EQ.32) THEN
          MDCY(KC,1)=0
          MWID(KC)=0
 
C...Higgs decays.
        ELSEIF(I.GE.33.AND.I.LE.36) THEN
C...Calculate decays to non-SUSY particles.
          CALL PYWIDT(KF,PMAS(KC,1)**2,WDTP,WDTE)
          LKNT=0
          DO 310 I1=0,100
            XLAM(I1)=0D0
  310     CONTINUE
          DO 330 I1=1,MDCY(KC,3)
            K1=MDCY(KC,2)+I1-1
            IF(IABS(KFDP(K1,1)).GT.KSUSY1.OR.
     &      IABS(KFDP(K1,2)).GT.KSUSY1) GOTO 330
            XLAM(I1)=WDTP(I1)
            XLAM(0)=XLAM(0)+XLAM(I1)
            DO 320 J1=1,3
              IDLAM(I1,J1)=KFDP(K1,J1)
  320       CONTINUE
            LKNT=LKNT+1
  330     CONTINUE
C...Add the decays to SUSY particles.
          CALL PYHEXT(KF,XLAM,IDLAM,LKNT)
        ENDIF
C...Zero the branching ratios for use in loop mode
C...thanks to K. Matchev (FNAL)
        DO 340 IDC=MDCY(KC,2),MDCY(KC,2)+MDCY(KC,3)-1
          BRAT(IDC)=0D0
  340   CONTINUE
 
C...Set stable particles.
        IF(LKNT.EQ.0) THEN
          MDCY(KC,1)=0
          MWID(KC)=0
          PMAS(KC,2)=1D-6
          PMAS(KC,3)=1D-5
          PMAS(KC,4)=0D0
 
C...Store branching ratios in the standard tables.
        ELSE
          IDC=MDCY(KC,2)+MDCY(KC,3)-1
          DELM=1D6
          DO 360 IL=1,LKNT
            IDCSV=IDC
  350       IDC=IDC+1
            BRAT(IDC)=0D0
            IF(IDC.EQ.MDCY(KC,2)+MDCY(KC,3)) IDC=MDCY(KC,2)
            IF(IDLAM(IL,1).EQ.KFDP(IDC,1).AND.IDLAM(IL,2).EQ.
     &      KFDP(IDC,2).AND.IDLAM(IL,3).EQ.KFDP(IDC,3)) THEN
              BRAT(IDC)=XLAM(IL)/XLAM(0)
              XMDIF=PMAS(KC,1)
              IF(MDME(IDC,1).GE.1) THEN
                XMDIF=XMDIF-PMAS(PYCOMP(KFDP(IDC,1)),1)-
     &          PMAS(PYCOMP(KFDP(IDC,2)),1)
                IF(KFDP(IDC,3).NE.0) XMDIF=XMDIF-
     &          PMAS(PYCOMP(KFDP(IDC,3)),1)
              ENDIF
              IF(I.LE.32) THEN
                IF(XMDIF.GE.0D0) THEN
                  DELM=MIN(DELM,XMDIF)
                ELSE
                  WRITE(MSTU(11),*) ' ERROR WITH DELM ',DELM,XMDIF
                  WRITE(MSTU(11),*) ' KF = ',KF
                  WRITE(MSTU(11),*) ' KF(decay) = ',(KFDP(IDC,J),J=1,3)
                ENDIF
              ENDIF
              GOTO 360
            ELSEIF(IDC.EQ.IDCSV) THEN
              WRITE(MSTU(11),*) ' Error in PYMSIN: SUSY decay ',
     &        'channel not recognized:'
              WRITE(MSTU(11),*) KF,' -> ',(IDLAM(IL,J),J=1,3)
              GOTO 360
            ELSE
              GOTO 350
            ENDIF
  360     CONTINUE
 
C...Store width, cutoff and lifetime.
          PMAS(KC,2)=XLAM(0)
          IF(PMAS(KC,2).LT.0.1D0*DELM) THEN
            PMAS(KC,3)=PMAS(KC,2)*10D0
          ELSE
            PMAS(KC,3)=0.95D0*DELM
          ENDIF
          IF(PMAS(KC,2).NE.0D0) THEN
            PMAS(KC,4)=PARU(3)/PMAS(KC,2)*1D-12
          ENDIF
C...Write decays to SLHA file
	  IF (IMSS(24).NE.0) THEN
            IFAIL=0
            CALL PYSLHA(4,KF,IFAIL)
          ENDIF
 
        ENDIF
  370 CONTINUE
 
      RETURN
      END
