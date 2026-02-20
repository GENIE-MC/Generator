 
C*********************************************************************
 
C...PYSUBH
C...This routine computes the renormalization group improved
C...values of Higgs masses and couplings in the MSSM.
 
C...Program based on the work by M. Carena, J.R. Espinosa,
c...M. Quiros and C.E.M. Wagner, CERN-preprint CERN-TH/95-45
 
C...Input: MA,TANB = TAN(BETA),MQ,MUR,MTOP,AU,AD,MU
C...All masses in GeV units. MA is the CP-odd Higgs mass,
C...MTOP is the physical top mass, MQ and MUR are the soft
C...supersymmetry breaking mass parameters of left handed
C...and right handed stops respectively, AU and AD are the
C...stop and sbottom trilinear soft breaking terms,
C...respectively,  and MU is the supersymmetric
C...Higgs mass parameter. We use the  conventions from
C...the physics report of Haber and Kane: left right
C...stop mixing term proportional to (AU - MU/TANB)
C...We use as input TANB defined at the scale MTOP
 
C...Output: MH,HM,MHCH, SA = SIN(ALPHA), CA= COS(ALPHA), TANBA
C...where MH and HM are the lightest and heaviest CP-even
C...Higgs masses, MHCH is the charged Higgs mass and
C...ALPHA is the Higgs mixing angle
C...TANBA is the angle TANB at the CP-odd Higgs mass scale
 
C...Range of validity:
C...(STOP1**2 - STOP2**2)/(STOP2**2 + STOP1**2) < 0.5
C...(SBOT1**2 - SBOT2**2)/(SBOT2**2 + SBOT2**2) < 0.5
C...where STOP1, STOP2, SBOT1 and SBOT2 are the stop and
C...are the sbottom  mass eigenvalues, respectively. This
C...range automatically excludes the existence of tachyons.
C...For the charged Higgs mass computation, the method is
C...valid if
C...2 * |MB * AD* TANB|  < M_SUSY**2,  2 * |MTOP * AU| < M_SUSY**2
C...2 * |MB * MU * TANB| < M_SUSY**2,  2 * |MTOP * MU| < M_SUSY**2
C...where M_SUSY**2 is the average of the squared stop mass
C...eigenvalues, M_SUSY**2 = (STOP1**2 + STOP2**2)/2. The sbottom
C...masses have been assumed to be of order of the stop ones
C...M_SUSY**2 = (MQ**2 + MUR**2)*0.5 + MTOP**2
 
      SUBROUTINE PYSUBH (XMA,TANB,XMQ,XMUR,XMTOP,AU,AD,XMU,XMH,XHM,
     &XMHCH,SA,CA,TANBA)
 
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
      COMMON/PYHTRI/HHH(7)
      SAVE /PYDAT1/,/PYDAT2/
 
C...Local variables.
      DOUBLE PRECISION PYALEM,PYALPS
      DOUBLE PRECISION TANB,XMQ,XMUR,XMTOP,AU,AD,XMU,XMH,XHM
      DOUBLE PRECISION XMHCH,SA,CA
      DOUBLE PRECISION XMA,AEM,ALP1,ALP2,ALPH3Z,V,PI
      DOUBLE PRECISION Q02
      DOUBLE PRECISION TANBA,TANBT,XMB,ALP3
      DOUBLE PRECISION RMTOP,XMS,T,SINB,COSB
      DOUBLE PRECISION XLAM1,XLAM2,XLAM3,XLAM4,XLAM5,XLAM6
      DOUBLE PRECISION XLAM7,XAU,XAD,G1,G2,G3,HU,HD,HU2
      DOUBLE PRECISION HD2,HU4,HD4,SINBT,COSBT
      DOUBLE PRECISION TRM2,DETM2,XMH2,XHM2,XMHCH2
      DOUBLE PRECISION SINALP,COSALP,AUD,PI2,XMS2,XMS4,AD2
      DOUBLE PRECISION AU2,XMU2,XMZ,XMS3
 
      XMZ = PMAS(23,1)
      Q02=XMZ**2
      AEM=PYALEM(Q02)
      ALP1=AEM/(1D0-PARU(102))
      ALP2=AEM/PARU(102)
      ALPH3Z=PYALPS(Q02)
 
      ALP1 = 0.0101D0
      ALP2 = 0.0337D0
      ALPH3Z = 0.12D0
 
      V = 174.1D0
      PI = PARU(1)
      TANBA = TANB
      TANBT = TANB
 
C...MBOTTOM(MTOP) = 3. GEV
      XMB = PYMRUN(5,XMTOP**2)
      ALP3 = ALPH3Z/(1D0 +(11D0 - 10D0/3D0)/4D0/PI*ALPH3Z*
     &LOG(XMTOP**2/XMZ**2))
 
C...RMTOP= RUNNING TOP QUARK MASS
      RMTOP = XMTOP/(1D0+4D0*ALP3/3D0/PI)
      XMS = ((XMQ**2 + XMUR**2)/2D0 + XMTOP**2)**0.5D0
      T = LOG(XMS**2/XMTOP**2)
      SINB = TANB/((1D0 + TANB**2)**0.5D0)
      COSB = SINB/TANB
C...IF(MA.LE.XMTOP) TANBA = TANBT
      IF(XMA.GT.XMTOP)
     &TANBA = TANBT*(1D0-3D0/32D0/PI**2*
     &(RMTOP**2/V**2/SINB**2-XMB**2/V**2/COSB**2)*
     &LOG(XMA**2/XMTOP**2))
 
      SINBT = TANBT/SQRT(1D0 + TANBT**2)
      COSBT = 1D0/SQRT(1D0 + TANBT**2)
C      COS2BT = (TANBT**2 - 1D0)/(TANBT**2 + 1D0)
      G1 = SQRT(ALP1*4D0*PI)
      G2 = SQRT(ALP2*4D0*PI)
      G3 = SQRT(ALP3*4D0*PI)
      HU = RMTOP/V/SINBT
      HD =  XMB/V/COSBT
      HU2=HU*HU
      HD2=HD*HD
      HU4=HU2*HU2
      HD4=HD2*HD2
      AU2=AU**2
      AD2=AD**2
      XMS2=XMS**2
      XMS3=XMS**3
      XMS4=XMS2*XMS2
      XMU2=XMU*XMU
      PI2=PI*PI
 
      XAU = (2D0*AU2/XMS2)*(1D0 - AU2/12D0/XMS2)
      XAD = (2D0*AD2/XMS2)*(1D0 - AD2/12D0/XMS2)
      AUD = (-6D0*XMU2/XMS2 - ( XMU2- AD*AU)**2/XMS4
     &+ 3D0*(AU + AD)**2/XMS2)/6D0
      XLAM1 = ((G1**2 + G2**2)/4D0)*(1D0-3D0*HD2*T/8D0/PI2)
     &+(3D0*HD4/8D0/PI2) * (T + XAD/2D0 + (3D0*HD2/2D0 + HU2/2D0
     &- 8D0*G3**2) * (XAD*T + T**2)/16D0/PI2)
     &-(3D0*HU4* XMU**4/96D0/PI2/XMS4) * (1+ (9D0*HU2 -5D0* HD2
     &-  16D0*G3**2) *T/16D0/PI2)
      XLAM2 = ((G1**2 + G2**2)/4D0)*(1D0-3D0*HU2*T/8D0/PI2)
     &+(3D0*HU4/8D0/PI2) * (T + XAU/2D0 + (3D0*HU2/2D0 + HD2/2D0
     &- 8D0*G3**2) * (XAU*T + T**2)/16D0/PI2)
     &-(3D0*HD4* XMU**4/96D0/PI2/XMS4) * (1+ (9D0*HD2 -5D0* HU2
     &-  16D0*G3**2) *T/16D0/PI2)
      XLAM3 = ((G2**2 - G1**2)/4D0)*(1D0-3D0*
     &(HU2 + HD2)*T/16D0/PI2)
     &+(6D0*HU2*HD2/16D0/PI2) * (T + AUD/2D0 + (HU2 + HD2
     &- 8D0*G3**2) * (AUD*T + T**2)/16D0/PI2)
     &+(3D0*HU4/96D0/PI2) * (3D0*XMU2/XMS2 - XMU2*AU2/
     &XMS4)* (1D0+ (6D0*HU2 -2D0* HD2/2D0
     &-  16D0*G3**2) *T/16D0/PI2)
     &+(3D0*HD4/96D0/PI2) * (3D0*XMU2/XMS2 - XMU2*AD2/
     &XMS4)*(1D0+ (6D0*HD2 -2D0* HU2
     &-  16D0*G3**2) *T/16D0/PI2)
      XLAM4 = (- G2**2/2D0)*(1D0-3D0*(HU2 + HD2)*T/16D0/PI2)
     &-(6D0*HU2*HD2/16D0/PI2) * (T + AUD/2D0 + (HU2 + HD2
     &- 8D0*G3**2) * (AUD*T + T**2)/16D0/PI2)
     &+(3D0*HU4/96D0/PI2) * (3D0*XMU2/XMS2 - XMU2*AU2/
     &XMS4)*
     &(1+ (6D0*HU2 -2D0* HD2
     &-  16D0*G3**2) *T/16D0/PI2)
     &+(3D0*HD4/96D0/PI2) * (3D0*XMU2/XMS2 - XMU2*AD2/
     &XMS4)*
     &(1+ (6D0*HD2 -2D0* HU2/2D0
     &-  16D0*G3**2) *T/16D0/PI2)
      XLAM5 = -(3D0*HU4* XMU2*AU2/96D0/PI2/XMS4) *
     &(1- (2D0*HD2 -6D0* HU2 + 16D0*G3**2) *T/16D0/PI2)
     &-(3D0*HD4* XMU2*AD2/96D0/PI2/XMS4) *
     &(1- (2D0*HU2 -6D0* HD2 + 16D0*G3**2) *T/16D0/PI2)
      XLAM6 = (3D0*HU4* XMU**3*AU/96D0/PI2/XMS4) *
     &(1- (7D0*HD2/2D0 -15D0* HU2/2D0 + 16D0*G3**2) *T/16D0/PI2)
     &+(3D0*HD4* XMU *(AD**3/XMS3 - 6D0*AD/XMS )/96D0/PI2/XMS) *
     &(1- (HU2/2D0 -9D0* HD2/2D0 + 16D0*G3**2) *T/16D0/PI2)
      XLAM7 = (3D0*HD4* XMU**3*AD/96D0/PI2/XMS4) *
     &(1- (7D0*HU2/2D0 -15D0* HD2/2D0 + 16D0*G3**2) *T/16D0/PI2)
     &+(3D0*HU4* XMU *(AU**3/XMS3 - 6D0*AU/XMS )/96D0/PI2/XMS) *
     &(1- (HD2/2D0 -9D0* HU2/2D0 + 16D0*G3**2) *T/16D0/PI2)
      HHH(1)=XLAM1
      HHH(2)=XLAM2
      HHH(3)=XLAM3
      HHH(4)=XLAM4
      HHH(5)=XLAM5
      HHH(6)=XLAM6
      HHH(7)=XLAM7
      TRM2 = XMA**2 + 2D0*V**2* (XLAM1* COSBT**2 +
     &2D0* XLAM6*SINBT*COSBT
     &+ XLAM5*SINBT**2 + XLAM2* SINBT**2 + 2D0* XLAM7*SINBT*COSBT
     &+ XLAM5*COSBT**2)
      DETM2 = 4D0*V**4*(-(SINBT*COSBT*(XLAM3 + XLAM4) +
     &XLAM6*COSBT**2
     &+ XLAM7* SINBT**2)**2 + (XLAM1* COSBT**2 +
     &2D0* XLAM6* COSBT*SINBT
     &+ XLAM5*SINBT**2)*(XLAM2* SINBT**2 +2D0* XLAM7* COSBT*SINBT
     &+ XLAM5*COSBT**2)) + XMA**2*2D0*V**2 *
     &((XLAM1* COSBT**2 +2D0*
     &XLAM6* COSBT*SINBT + XLAM5*SINBT**2)*COSBT**2 +
     &(XLAM2* SINBT**2 +2D0* XLAM7* COSBT*SINBT + XLAM5*COSBT**2)
     &*SINBT**2
     &+2D0*SINBT*COSBT* (SINBT*COSBT*(XLAM3
     &+ XLAM4) + XLAM6*COSBT**2
     &+ XLAM7* SINBT**2))
 
      XMH2 = (TRM2 - SQRT(TRM2**2 - 4D0* DETM2))/2D0
      XHM2 = (TRM2 + SQRT(TRM2**2 - 4D0* DETM2))/2D0
      XHM = SQRT(XHM2)
      XMH = SQRT(XMH2)
      XMHCH2 = XMA**2 + (XLAM5 - XLAM4)* V**2
      XMHCH = SQRT(XMHCH2)
 
      SINALP = SQRT(((TRM2**2 - 4D0* DETM2)**0.5D0) -
     &((2D0*V**2*(XLAM1* COSBT**2 + 2D0*
     &XLAM6* COSBT*SINBT
     &+ XLAM5*SINBT**2) + XMA**2*SINBT**2)
     &- (2D0*V**2*(XLAM2* SINBT**2 +2D0* XLAM7* COSBT*SINBT
     &+ XLAM5*COSBT**2) + XMA**2*COSBT**2)))/
     &SQRT(((TRM2**2 - 4D0* DETM2)**0.5D0))/2D0**0.5D0
 
      COSALP = (2D0*(2D0*V**2*(SINBT*COSBT*(XLAM3 + XLAM4) +
     &XLAM6*COSBT**2 + XLAM7* SINBT**2) -
     &XMA**2*SINBT*COSBT))/2D0**0.5D0/
     &SQRT(((TRM2**2 - 4D0* DETM2)**0.5D0)*
     &(((TRM2**2 - 4D0* DETM2)**0.5D0) -
     &((2D0*V**2*(XLAM1* COSBT**2 + 2D0*
     &XLAM6* COSBT*SINBT
     &+ XLAM5*SINBT**2) + XMA**2*SINBT**2)
     &- (2D0*V**2*(XLAM2* SINBT**2 +2D0* XLAM7* COSBT*SINBT
     &+ XLAM5*COSBT**2) + XMA**2*COSBT**2))))
 
      SA = -SINALP
      CA = -COSALP
 
  100 CONTINUE
 
      RETURN
      END
