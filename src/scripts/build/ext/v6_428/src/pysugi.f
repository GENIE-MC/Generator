 
C*********************************************************************
 
C...PYSUGI
C...Interface to ISASUSY version 7.71.
C...Warning: this interface should not be used with earlier versions
C...of ISASUSY, since common block incompatibilities may then arise.
C...Calls SUGRA (in ISAJET) to perform RGE evolution.
C...Then converts to Gunion-Haber conventions.
 
      SUBROUTINE PYSUGI
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
 
      INTEGER PYK,PYCHGE,PYCOMP
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
 
C...Date of Change
      CHARACTER DOC*11
      PARAMETER (DOC='01 May 2006')
 
C...ISASUGRA Input:
      REAL MZERO,MHLF,AZERO,TANB,SGNMU,MTOP
C...XISAIN contains the MSSMi inputs in natural order.
      COMMON /SUGXIN/ XISAIN(24),XSUGIN(7),XGMIN(14),XNRIN(4),
     $XAMIN(7)
      REAL XISAIN,XSUGIN,XGMIN,XNRIN,XAMIN
      SAVE /SUGXIN/
C...ISASUGRA Output
      CHARACTER*40 ISAVER,VISAJE
      REAL SUPER
      COMMON /SSPAR/ SUPER(72)
      COMMON /SUGMG/ MSS(32),GSS(31),MGUTSS,GGUTSS,AGUTSS,FTGUT,
     $FBGUT,FTAGUT,FNGUT
      REAL MSS,GSS,MGUTSS,GGUTSS,AGUTSS,FTGUT,FBGUT,FTAGUT,FNGUT
      COMMON /SUGPAS/ XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,NOGOOD,IAL3UN,ITACHY,MHPNEG,ASM3,
     $VUMT,VDMT,ASMTP,ASMSS,M3Q
      REAL XTANB,MSUSY,AMT,MGUT,MU,G2,GP,V,VP,XW,
     $A1MZ,A2MZ,ASMZ,FTAMZ,FBMZ,B,SIN2B,FTMT,G3MT,VEV,HIGFRZ,
     $FNMZ,AMNRMJ,ASM3,VUMT,VDMT,ASMTP,ASMSS,M3Q
      INTEGER NOGOOD,IAL3UN,ITACHY,MHPNEG
      INTEGER IALLOW
      SAVE /SUGMG/,/SSPAR/
C SUPER: Filled by ISASUGRA.
C SUPER(1)        = mass of ~g
C SUPER(2:17)     = mass of ~u_L,~u_R,~d_L,~d_R,~s_L,~s_R,~c_L,~c_R,~b_L
C                          ,~b_R,~b_1,~b_2,~t_L,~t_R,~t_1,~t_2
C SUPER(18:25)    = mass of ~e_L,~e_R,~mu_L,~mu_R,~tau_L,~tau_R,~tau_1
C                          ,~tau_2
C SUPER(26:28)    = mass of ~nu_e,~nu_mu,~nu_tau
C SUPER(29)       = Higgsino mass = - mu
C SUPER(30)       = ratio v2/v1 of vev's
C SUPER(31:34)    = Signed neutralino masses
C SUPER(35:50)    = Neutralino mixing matrix
C SUPER(51:52)    = Signed chargino masses
C SUPER(53:54)    = Chargino left, right mixing angles
C SUPER(55:58)    = mass of h0, H0, A0, H+
C SUPER(59)       = Higgs mixing angle alpha
C SUPER(60:65)    = A_t, theta_t, A_b, theta_b, A_tau, theta_tau
C SUPER(66)       = Gravitino mass
C SUPER(67:69)    = Top,Bottom, and Tau masses at MSUSY (not used)
C SUPER(70)       = b-Yukawa at mA scale (not used)
C SUPER(71:72)    = H_u, H_d vev's at MSUSY (not used)
C GSS: Filled by ISASUGRA
C     GSS( 1) = g_1        GSS( 2) = g_2        GSS( 3) = g_3
C     GSS( 4) = y_tau      GSS( 5) = y_b        GSS( 6) = y_t
C     GSS( 7) = M_1        GSS( 8) = M_2        GSS( 9) = M_3
C     GSS(10) = A_tau      GSS(11) = A_b        GSS(12) = A_t
C     GSS(13) = M_h12     GSS(14) = M_h22     GSS(15) = M_er2
C     GSS(16) = M_el2     GSS(17) = M_dnr2    GSS(18) = M_upr2
C     GSS(19) = M_upl2    GSS(20) = M_taur2   GSS(21) = M_taul2
C     GSS(22) = M_btr2    GSS(23) = M_tpr2    GSS(24) = M_tpl2
C     GSS(25) = mu         GSS(26) = B          GSS(27) = Y_N
C     GSS(28) = M_nr       GSS(29) = A_n        GSS(30) = log(vdq)
C     GSS(31) = log(vuq)
C MSS: Filled by ISASUGRA
C     MSS( 1) = glss     MSS( 2) = upl      MSS( 3) = upr
C     MSS( 4) = dnl      MSS( 5) = dnr      MSS( 6) = stl
C     MSS( 7) = str      MSS( 8) = chl      MSS( 9) = chr
C     MSS(10) = b1       MSS(11) = b2       MSS(12) = t1
C     MSS(13) = t2       MSS(14) = nuel     MSS(15) = numl
C     MSS(16) = nutl     MSS(17) = el-      MSS(18) = er-
C     MSS(19) = mul-     MSS(20) = mur-     MSS(21) = tau1
C     MSS(22) = tau2     MSS(23) = z1ss     MSS(24) = z2ss
C     MSS(25) = z3ss     MSS(26) = z4ss     MSS(27) = w1ss
C     MSS(28) = w2ss     MSS(29) = hl0      MSS(30) = hh0
C     MSS(31) = ha0      MSS(32) = h+
C Unification, filled by ISASUGRA if applicable.
C     MGUTSS  = M_GUT    GGUTSS  = g_GUT    AGUTSS  = alpha_GUTC
 
C...SPYTHIA Input/Output
      INTEGER IMSS
      DOUBLE PRECISION RMSS
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
C...SLHA Input/Output
      COMMON/PYLH3P/MODSEL(200),PARMIN(100),PAREXT(200),RMSOFT(0:100),
     &     AU(3,3),AD(3,3),AE(3,3)
C...PYTHIA common blocks
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
 
      SAVE  /PYMSSM/,/PYSSMT/,/PYLH3P/,/PYDAT1/,/PYPARS/,/PYDAT2/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INTEGER IMODEL
      REAL M0,MHF,A0,MT
      CHARACTER*20 CHMOD(5)
      CHARACTER*32 FNAME
 
      COMMON /SUGNU/ XNUSUG(18)
      REAL XNUSUG
      SAVE /SUGNU/
 
      DATA CHMOD/'mSUGRA','mGMSB','non-universal SUGRA',
     &     'truly unified SUGRA', 'non-minimal GMSB'/
 
C...Start by checking for incompatibilities/inconsistencies:
      DO 100 ICHK=2,9
        IF (ICHK.NE.8.AND.ICHK.NE.4.AND.IMSS(ICHK).NE.0) THEN
          WRITE (MSTU(11),*) '(PYSUGI:) IMSS(',ICHK,')=',IMSS(ICHK)
     &         ,' option not used by PYSUGI'
        ENDIF
  100 CONTINUE
C...ISAJET works with REAL numbers.
      MZERO=REAL(RMSS(8))
      MHLF=REAL(RMSS(1))
      AZERO=REAL(RMSS(16))
      TANB=REAL(RMSS(5))
      SGNMU=REAL(RMSS(4))
      MTOP=REAL(PMAS(6,1))
      IMODEL=0
      IF (IMSS(1).EQ.12) THEN
        IMODEL=1
        GOTO 130
      ELSEIF(IMSS(1).EQ.13) THEN
C...Read from isajet par file in IMSS(20)
        LFN=IMSS(20)
C...STOP IF LFN IS ZERO (i.e. if no LFN was given).
        IF (LFN.EQ.0) THEN
          WRITE(MSTU(11),*) '(PYSUGI:) No valid unit given in IMSS(20)'
          GOTO 9999
        ENDIF
        WRITE(MSTU(11),*) 'READING SUSY MODEL FROM FILE...'
CMrenna change to allow any susy model
        WRITE(MSTU(11),*) 'ENTER 1 for mSUGRA:'
        WRITE(MSTU(11),*) 'ENTER 2 for mGMSB:'
        WRITE(MSTU(11),*) 'ENTER 3 for non-universal SUGRA:'
        WRITE(MSTU(11),*) 'ENTER 4 for SUGRA with truly unified'//
     &       ' gauge couplings:'
        WRITE(MSTU(11),*) 'ENTER 5 for non-minimal GMSB:'
        READ(LFN,*) IMODEL
        IF (IMODEL.EQ.4) THEN
          IAL3UN=1
          IMODEL=1
        ENDIF
        IF (IMODEL.EQ.1.OR.IMODEL.EQ.3) THEN
          WRITE(MSTU(11),*) 'ENTER M_0, M_(1/2), A_0, tan(beta),'
     &         //' sgn(mu), M_t:'
          READ(LFN,*) M0,MHF,A0,TANB,SGNMU,MT
          IF (IMODEL.EQ.3) THEN
            IMODEL=1
 110        WRITE(MSTU(11),*) ' ENTER 1,...,5 for NUSUGx keyword;'
     &           //' 0 to continue:'
            WRITE(MSTU(11),*) ' NUSUG1 = GUT scale gaugino masses'
            WRITE(MSTU(11),*) ' NUSUG2 = GUT scale A terms'
            WRITE(MSTU(11),*) ' NUSUG3 = GUT scale Higgs masses'
            WRITE(MSTU(11),*) ' NUSUG4 = GUT scale 1st/2nd'
     &           //' generation masses'
            WRITE(MSTU(11),*)
     &           ' NUSUG5 = GUT scale 3rd generation masses'
            READ(LFN,*) INUSUG
            IF (INUSUG.EQ.0) THEN
              GOTO 120
            ELSEIF (INUSUG.EQ.1) THEN
              WRITE(MSTU(11),*) 'Enter GUT scale M_1, M_2, M_3:'
              READ(LFN,*) XNUSUG(1),XNUSUG(2),XNUSUG(3)
              IF (XNUSUG(3).LE.0.) THEN
                WRITE(MSTU(11),*) ' NEGATIVE M_3 IS NOT ALLOWED'
                CALL PYSTOP(109)
              END IF
            ELSEIF (INUSUG.EQ.2) THEN
              WRITE(MSTU(11),*) 'Enter GUT scale A_t, A_b, A_tau:'
              READ(LFN,*) XNUSUG(6),XNUSUG(5),XNUSUG(4)
            ELSEIF (INUSUG.EQ.3) THEN
              WRITE(MSTU(11),*) 'Enter GUT scale m_Hd, m_Hu:'
              READ(LFN,*) XNUSUG(7),XNUSUG(8)
            ELSEIF (INUSUG.EQ.4) THEN
              WRITE(MSTU(11),*) 'Enter GUT scale M(ul), M(dr),'
     &             //' M(ur), M(el), M(er):'
              READ(LFN,*) XNUSUG(13),XNUSUG(11),XNUSUG(12),
     &             XNUSUG(10),XNUSUG(9)
            ELSEIF (INUSUG.EQ.5) THEN
              WRITE(MSTU(11),*) 'Enter GUT scale M(tl), M(br), M(tr),'
     &              //' M(Ll), M(Lr):'
              READ(LFN,*) XNUSUG(18),XNUSUG(16),XNUSUG(17),
     &             XNUSUG(15),XNUSUG(14)
            ENDIF
            GOTO 110
          ENDIF
        ELSEIF (IMODEL.EQ.2.OR.IMODEL.EQ.5) THEN
          IMSS(11)=1
          WRITE(MSTU(11),*) 'ENTER Lambda, M_mes, N_5, tan(beta),'
     &         ,' sgn(mu), M_t, C_gv:'
          READ(LFN,*) M0,MHF,A0,TANB,SGNMU,MT,XCMGV
          XGMIN(7)=XCMGV
          XGMIN(8)=1.
C...Planck scale: AMPL = 2.4 E18 GeV = {8 pi G_newton}^{1/2}
          AMPL=2.4D18
          AMGVSS=M0*MHF*XCMGV/SQRT(3D0)/AMPL
          IF (IMODEL.EQ.5) THEN
            IMODEL=2
            WRITE(MSTU(11),*) 'Rsl = factor multiplying gaugino'
     &           ,' masses at M_mes'
            WRITE(MSTU(11),*) 'dmH_d2, dmH_u2 = Higgs mass**2'
     &           ,' shifts at M_mes'
            WRITE(MSTU(11),*) 'd_Y = mass**2 shifts proportional to',
     &           ' Y at M_mes'
            WRITE(MSTU(11),*) 'n5_1,n5_2,n5_3 = n5 values for U(1),'
     &           ,'SU(2),SU(3)'
            WRITE(MSTU(11),*) 'ENTER Rsl, dmH_d2, dmH_u2, d_Y, n5_1,'
     &           ,' n5_2, n5_3'
            READ(LFN,*) XGMIN(8),XGMIN(9),XGMIN(10),XGMIN(11),XGMIN(12),
     $           XGMIN(13),XGMIN(14)
          ENDIF
        ELSE
          WRITE(MSTU(11),*) 'Invalid model choice.'
          GOTO 9999
        ENDIF
      ENDIF
 
 120  MZERO=M0
      MHLF=MHF
      AZERO=A0
C     TANB=REAL(RMSS(5))
C     SGNMU=REAL(RMSS(4))
      MTOP=MT
 
C...Initialize MSSM parameter array
 130  DO 140 IPAR=1,72
        SUPER(IPAR)=0.0
 140  CONTINUE
C...Call ISASUGRA
      CALL SUGRA(MZERO,MHLF,AZERO,TANB,SGNMU,MTOP,IMODEL)
C...Check whether ISASUSY thought the model was OK.
      IF (NOGOOD.NE.0) THEN
        IF (NOGOOD.EQ.1) CALL PYERRM(26
     &       ,'(PYSUGI:) SUSY parameters give tachyonic particles.')
        IF (NOGOOD.EQ.2) CALL PYERRM(26
     &       ,'(PYSUGI:) SUSY parameters give no EWSB.')
        IF (NOGOOD.EQ.3) CALL PYERRM(26
     &       ,'(PYSUGI:) SUSY parameters give m(A0) < 0.')
        IF (NOGOOD.EQ.4) CALL PYERRM(26
     &       ,'(PYSUGI:) SUSY parameters give Yukawa > 100.')
        IF (NOGOOD.EQ.7) CALL PYERRM(26
     &       ,'(PYSUGI:) SUSY parameters give x_T EWSB bad.')
        IF (NOGOOD.EQ.8) CALL PYERRM(26
     &       ,'(PYSUGI:) SUSY parameters give m(h0)2 < 0.')
C...Give warning, but don't stop, if LSP not ~chi_10.
        IF (NOGOOD.EQ.5) CALL PYERRM(16
     &       ,'(PYSUGI:) SUSY parameters give ~chi_10 not LSP.')
      ENDIF
C...Warn about possible GUT scale tachyons.
      IF (ITACHY.NE.0) CALL PYERRM(16,
     &       '(PYSUGI:) Tachyonic sleptons at GUT scale.')
C...Finalize spectrum (last iteration)
C...(Thanks to A. Raklev for pointing this out.)
C...NB: SSMSSM also calculates decays, but these are not used by Pythia.
      CALL SSMSSM(XISAIN(1),XISAIN(2),XISAIN(3),
     $ XISAIN(4),XISAIN(5),XISAIN(6),XISAIN(7),XISAIN(8),XISAIN(9),
     $ XISAIN(10),XISAIN(11),XISAIN(12),XISAIN(13),XISAIN(14),
     $ XISAIN(15),XISAIN(16),XISAIN(17),XISAIN(18),XISAIN(19),
     $ XISAIN(20),XISAIN(21),XISAIN(22),XISAIN(23),XISAIN(24),
     $ MTOP,IALLOW,1)
 
C...M1, M2, M3.
      RMSS(1)=dble(GSS(7))
      RMSS(2)=dble(GSS(8))
      RMSS(3)=dble(GSS(9))
      RMSOFT(1)=dble(GSS(7))
      RMSOFT(2)=dble(GSS(8))
      RMSOFT(3)=dble(GSS(9))
C...Mu = - Higgsino mass.
      RMSS(4)=-SUPER(29)
      RMSS(5)=TANB
C...Slepton and squark masses. 2 first generations.
      RMSS(6)=0.5*(SUPER(18)+SUPER(20))
      RMSS(7)=0.5*(SUPER(19)+SUPER(21))
      RMSS(8)=0.25*(SUPER(2)+SUPER(4)+SUPER(6)+SUPER(8))
      RMSS(9)=0.25*(SUPER(3)+SUPER(5)+SUPER(7)+SUPER(9))
C...Third generation.
      RMSS(10)=0.5*(SUPER(14)+SUPER(10))
      RMSS(11)=SUPER(11)
      RMSS(12)=SUPER(15)
      RMSS(13)=SUPER(22)
      RMSS(14)=SUPER(23)
C...SLHA: store exact soft spectrum in RMSOFT
      RMSOFT(31)=SUPER(18)
      RMSOFT(32)=SUPER(20)
      RMSOFT(33)=SUPER(22)
      RMSOFT(34)=SUPER(19)
      RMSOFT(35)=SUPER(21)
      RMSOFT(36)=SUPER(23)
      RMSOFT(41)=0.5D0*(SUPER(2)+SUPER(4))
      RMSOFT(42)=0.5D0*(SUPER(6)+SUPER(8))
      RMSOFT(43)=0.5D0*(SUPER(10)+SUPER(14))
      RMSOFT(44)=SUPER(3)
      RMSOFT(45)=SUPER(9)
      RMSOFT(46)=SUPER(15)
      RMSOFT(47)=SUPER(5)
      RMSOFT(48)=SUPER(7)
      RMSOFT(49)=SUPER(11)
 
C...~b, ~t, and ~tau trilinear couplings and mixing angles.
      RMSS(15)=SUPER(62)
      RMSS(16)=SUPER(60)
      RMSS(17)=SUPER(64)
      RMSS(26)=SUPER(63)
      RMSS(27)=SUPER(61)
      RMSS(28)=SUPER(65)
C...SLHA trilinears
      DO 142 K1=1,3
        DO 141 K2=1,3
          AE(K1,K2)=0D0
          AU(K1,K2)=0D0
          AD(K1,K2)=0D0
 141    CONTINUE
 142  CONTINUE
      AE(3,3)=SUPER(64)
      AU(3,3)=SUPER(60)
      AD(3,3)=SUPER(62)
C...Higgs mixing angle alpha (Gunion-Haber convention).
      RMSS(18)=-SUPER(59)
C...A0 mass.
      RMSS(19)=SUPER(57)
C...GUT scale coupling
      RMSS(20)=AGUTSS
C...Gravitino mass (for future compatibility)
      RMSS(21)=MAX(RMSS(21),DBLE(SUPER(66)))
 
C...Now we're done with RMSS. Time to fill PMAS (m > 0 required).
C...Higgs sector.
      PMAS(PYCOMP(25),1)=ABS(SUPER(55))
      PMAS(PYCOMP(35),1)=ABS(SUPER(56))
      PMAS(PYCOMP(36),1)=ABS(SUPER(57))
      PMAS(PYCOMP(37),1)=ABS(SUPER(58))
C...Gluino.
      PMAS(PYCOMP(KSUSY1+21),1)=ABS(SUPER(1))
C...Squarks and Sleptons.
      DO 150 ILR=1,2
        ILRM=ILR-1
        PMAS(PYCOMP(ILR*KSUSY1+1),1)=ABS(SUPER(4+ILRM))
        PMAS(PYCOMP(ILR*KSUSY1+2),1)=ABS(SUPER(2+ILRM))
        PMAS(PYCOMP(ILR*KSUSY1+3),1)=ABS(SUPER(6+ILRM))
        PMAS(PYCOMP(ILR*KSUSY1+4),1)=ABS(SUPER(8+ILRM))
        PMAS(PYCOMP(ILR*KSUSY1+5),1)=ABS(SUPER(12+ILRM))
        PMAS(PYCOMP(ILR*KSUSY1+6),1)=ABS(SUPER(16+ILRM))
        PMAS(PYCOMP(ILR*KSUSY1+11),1)=ABS(SUPER(18+ILRM))
        PMAS(PYCOMP(ILR*KSUSY1+13),1)=ABS(SUPER(20+ILRM))
        PMAS(PYCOMP(ILR*KSUSY1+15),1)=ABS(SUPER(24+ILRM))
  150 CONTINUE
      PMAS(PYCOMP(KSUSY1+12),1)=ABS(SUPER(26))
      PMAS(PYCOMP(KSUSY1+14),1)=ABS(SUPER(27))
      PMAS(PYCOMP(KSUSY1+16),1)=ABS(SUPER(28))
C...Neutralinos.
      PMAS(PYCOMP(KSUSY1+22),1)=ABS(SUPER(31))
      PMAS(PYCOMP(KSUSY1+23),1)=ABS(SUPER(32))
      PMAS(PYCOMP(KSUSY1+25),1)=ABS(SUPER(33))
      PMAS(PYCOMP(KSUSY1+35),1)=ABS(SUPER(34))
C...Signed masses (extra minus from going to G-H convention).
      SMZ(1)=-SUPER(31)
      SMZ(2)=-SUPER(32)
      SMZ(3)=-SUPER(33)
      SMZ(4)=-SUPER(34)
C...Charginos
      PMAS(PYCOMP(KSUSY1+24),1)=ABS(SUPER(51))
      PMAS(PYCOMP(KSUSY1+37),1)=ABS(SUPER(52))
C...Signed masses (extra minus from going to G-H convention).
      SMW(1)=-SUPER(51)
      SMW(2)=-SUPER(52)
 
C... Neutralino Mixing.
      DO 160 IN=1,4
        ZMIX(IN,1)= SUPER(38+4*(IN-1))
        ZMIX(IN,2)= SUPER(37+4*(IN-1))
        ZMIX(IN,3)=-SUPER(36+4*(IN-1))
        ZMIX(IN,4)=-SUPER(35+4*(IN-1))
  160 CONTINUE
C...Chargino Mixing (PYTHIA same angle as HERWIG).
      THX=1D0
      THY=1D0
      IF (SUPER(53).GT.0) THX=-1D0
      IF (SUPER(54).GT.0) THY=-1D0
      UMIX(1,1) = -SIN(SUPER(53))
      UMIX(1,2) = -COS(SUPER(53))
      UMIX(2,1) = -THX*COS(SUPER(53))
      UMIX(2,2) = THX*SIN(SUPER(53))
      VMIX(1,1) = -SIN(SUPER(54))
      VMIX(1,2) = -COS(SUPER(54))
      VMIX(2,1) = -THY*COS(SUPER(54))
      VMIX(2,2) = THY*SIN(SUPER(54))
C...Sfermion mixing (PYTHIA same angle as ISAJET)
      SFMIX(5,1)=COS(SUPER(63))
      SFMIX(5,2)=SIN(SUPER(63))
      SFMIX(5,3)=-SIN(SUPER(63))
      SFMIX(5,4)=COS(SUPER(63))
      SFMIX(6,1)=COS(SUPER(61))
      SFMIX(6,2)=SIN(SUPER(61))
      SFMIX(6,3)=-SIN(SUPER(61))
      SFMIX(6,4)=COS(SUPER(61))
      SFMIX(15,1)=COS(SUPER(65))
      SFMIX(15,2)=SIN(SUPER(65))
      SFMIX(15,3)=-SIN(SUPER(65))
      SFMIX(15,4)=COS(SUPER(65))
 
      IF (MSTP(122).NE.0) THEN
C...Print a few lines to make the user know what's happening
        ISAVER=VISAJE()
        WRITE(MSTU(11),5000) DOC, ISAVER
        WRITE(MSTU(11),5100)
        IF (IMODEL.EQ.1) THEN
          WRITE(MSTU(11),5200) MZERO, MHLF, AZERO, TANB, NINT(SGNMU),
     &         MTOP
          WRITE(MSTU(11),5300)
        ENDIF
        WRITE(MSTU(11),5500) 'Pole masses'
        WRITE(MSTU(11),5700) (SUPER(IP),IP=2,16,2),(SUPER(IP),IP=3,17,2)
        WRITE(MSTU(11),5800) (SUPER(IP),IP=18,24,2),(SUPER(IP),IP=26,28)
     &       ,(SUPER(IP),IP=19,25,2)
        WRITE(MSTU(11),5900) SUPER(1),(SMZ(IP),IP=1,4), (SMW(IP)
     &       ,IP=1,2)
        WRITE(MSTU(11),5400)
        WRITE(MSTU(11),6000) (SUPER(IP),IP=55,58)
        WRITE(MSTU(11),5400)
        WRITE(MSTU(11),5500) 'EW scale mixing structure'
        WRITE(MSTU(11),6100) ((ZMIX(I,J), J=1,4),I=1,4)
        WRITE(MSTU(11),6200) (UMIX(1,J), J=1,2),(VMIX(1,J),J=1,2)
     &       ,(UMIX(2,J), J=1,2),(VMIX(2,J),J=1,2)
        WRITE(MSTU(11),6300) (SFMIX(5,J), J=1,2),(SFMIX(6,J),J=1,2)
     &       ,(SFMIX(15,J), J=1,2),(SFMIX(5,J),J=3,4),(SFMIX(6,J), J=3,4
     &       ),(SFMIX(15,J),J=3,4)
        WRITE(MSTU(11),5400)
        WRITE(MSTU(11),6450) RMSS(18)
        WRITE(MSTU(11),5400)
        WRITE(MSTU(11),5500) 'Couplings'
        WRITE(MSTU(11),6400) RMSS(15),RMSS(16),RMSS(17),RMSS(20)
        WRITE(MSTU(11),5400)
      ENDIF
 
C...Call FeynHiggs to improve Higgs sector if requested
      IF (IMSS(4).EQ.3) THEN
        IF (MSTP(122).NE.0) WRITE(MSTU(11),'(1x,"*"/1x,"*",A)')
     &       ' (PYSUGI:) Now calling FeynHiggs.'
        CALL PYFEYN(IERR)
        IF (IERR.EQ.0) THEN
          IMSS(4)=2
          IF (MSTP(122).NE.0) THEN
            WRITE(MSTU(11),5400)
            WRITE(MSTU(11),5500)
     &           'Corrected Higgs masses and mixing'
            WRITE(MSTU(11),6000) PMAS(25,1),PMAS(35,1),PMAS(36,1),
     &           PMAS(37,1)
            WRITE(MSTU(11),6450) RMSS(18)
            WRITE(MSTU(11),5400)
          ENDIF
        ENDIF
      ENDIF
 
      IF (MSTP(122).NE.0) WRITE(MSTU(11),6500)
 
C...Fix the higgs sector (in PYMSIN) using the masses and mixing angle
C...output by ISASUSY.
      IMSS(4)=MAX(2,IMSS(4))
 
 5000 FORMAT(1x,19('*'),1x,'PYSUGI v1.52: PYTHIA/ISASUSY '
     &     ,'INTERFACE',1x,19('*')/1x,'*',3x,'PYSUGI: Last Change',1x,A
     &     ,1x,'-',1x,'P. Skands / S. Mrenna'/1x,'*',2x,A/1x,'*')
 5100 FORMAT(1x,'*',1x,'ISASUSY Input:'/1x,'*',1x,'----------------')
 5200 FORMAT(1x,'*',1x,3x,'M_0',6x,'M_1/2',5x,'A_0',3x,'Tan(beta)',
     &     3x,'Sgn(mu)',3x,'M_t'/1x,'*',1x,4(F8.2,1x),I8,2x,F8.2)
 5300 FORMAT(1x,'*'/1x,'*',1x,'ISASUSY Output:'/1x,'*',1x
     &     ,'----------------')
 5400 FORMAT(1x,'*',1x,A)
 5500 FORMAT(1x,'*',1x,A,':')
 5600 FORMAT(1x,'*',2x,2x,'M_GUT',2x,2x,'g_GUT',2x,1x,'alpha_GUT'/
     &       1x,'*',2x,1P,2(1x,E8.2),2x,E8.2)
 5700 FORMAT(1x,'*',4x,4x,'~u',2x,1x,4x,'~d',2x,1x,4x,'~s',2x,1x,
     &     4x,'~c',2x,1x,4x,'~b',2x,1x,2x,'~b(12)',1x,4x,'~t',2x,1x, 2x,
     &     '~t(12)'/1x,'*',2x,'L',1x,8(F8.2,1x)/1x,'*',2x,'R',1x,8(F8.2
     &     ,1x))
 5800 FORMAT(1x,'*'/1x,'*',4x,4x,'~e',2x,1x,3x,'~mu',2x,1x,3x,'~tau',1x
     &     ,1x,'~tau(12)',1x,2x,'~nu_e',1x,1x,1x,'~nu_mu',1x,1x,1x
     &     ,'~nu_tau'/1x,'*',2x,'L',1x,7(F8.2,1x)/1x,'*',2x,'R',1x,4(F8
     &     .2,1x))
 5900 FORMAT(1x,'*'/1x,'*',4x,4x,'~g',2x,1x,1x,'~chi_10',1x,1x,'~chi_20'
     &     ,1x,1x,'~chi_30',1x,1x,'~chi_40',1x,1x,'~chi_1+',1x
     &     ,1x,'~chi_2+'/1x,'*',3x,1x,7(F8.2,1x))
 6000 FORMAT(1x,'*',4x,4x,'h0',2x,1x,4x,'H0',2x,1x,4x,'A0',2x
     &     ,1x,4x,'H+'/1x,'*',3x,1x,5(F8.2,1x))
 6050 FORMAT(1x,'*'/1x,'*',4x,4x,'h0',2x,1x,4x,'H0',2x,1x,4x,'A0',2x
     &     ,1x,4x,'H+'/1x,'*',3x,1x,5(F8.2,1x),3x,'(Before FeynHiggs)')
 6100 FORMAT(1x,'*',11x,'|',3x,'~B',3x,'|',2x,'~W_3',2x,'|',2x
     &     ,'~H_1',2x,'|',2x,'~H_2',2x,'|'/1x,'*',3x,'~chi_10',1x,4('|'
     &     ,1x,F6.3,1x),'|'/1x,'*',3x,'~chi_20',1x,4('|'
     &     ,1x,F6.3,1x),'|'/1x,'*',3x,'~chi_30',1x,4('|'
     &     ,1x,F6.3,1x),'|'/1x,'*',3x,'~chi_40',1x,4('|'
     &     ,1x,F6.3,1x),'|')
 6200 FORMAT(1x,'*'/1x,'*',6x,'L',4x,'|',3x,'~W',3x,'|',3x,'~H',3x,'|'
     &     ,12x,'R',4x,'|',3x,'~W',3x,'|',3x,'~H',3x,'|'/1x,'*',3x
     &     ,'~chi_1+',1x,2('|',1x,F6.3,1x),'|',9x,'~chi_1+',1x,2('|',1x
     &     ,F6.3,1x),'|'/1x,'*',3x,'~chi_2+',1x,2('|',1x,F6.3,1x),'|',9x
     &     ,'~chi_2+',1x,2('|',1x,F6.3,1x),'|')
 6300 FORMAT(1x,'*'/1x,'*',8x,'|',2x,'~b_L',2x,'|',2x,'~b_R',2x,'|',8x
     &     ,'|',2x,'~t_L',2x,'|',2x,'~t_R',2x,'|',10x
     &     ,'|',1x,'~tau_L',1x,'|',1x,'~tau_R',1x,'|'/
     &     1x,'*',3x,'~b_1',1x,2('|',1x,F6.3,1x),'|',3x,'~t_1',1x,2('|'
     &     ,1x,F6.3,1x),'|',3x,'~tau_1',1x,2('|',1x,F6.3,1x),'|'/
     &     1x,'*',3x,'~b_2',1x,2('|',1x,F6.3,1x),'|',3x,'~t_2',1x,2('|'
     &     ,1x,F6.3,1x),'|',3x,'~tau_2',1x,2('|',1x,F6.3,1x),'|')
 6400 FORMAT(1x,'*',3x,'A_b = ',F8.2,4x,'A_t = ',F8.2,4x,'A_tau = ',F8.2
     &     ,4x,'Alpha_GUT = ',F8.2)
 6450 FORMAT(1x,'*',3x,'Alpha_Higgs = ',F8.4)
 6500 FORMAT(1x,32('*'),1x,'END OF PYSUGI',1x,31('*'))
 
 9999 RETURN
      END
