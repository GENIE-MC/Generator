 
C*********************************************************************
 
C...PYGGAM
C...Constructs the F2 and parton distributions of the photon
C...by summing homogeneous (VMD) and inhomogeneous (anomalous) terms.
C...For F2, c and b are included by the Bethe-Heitler formula;
C...in the 'MSbar' scheme additionally a Cgamma term is added.
C...Contains the SaS sets 1D, 1M, 2D and 2M.
C...Adapted from SaSgam library, authors G.A. Schuler and T. Sjostrand.
 
      SUBROUTINE PYGGAM(ISET,X,Q2,P2,IP2,F2GM,XPDFGM)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYINT8/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6),
     &XPDIR(-6:6)
      COMMON/PYINT9/VXPVMD(-6:6),VXPANL(-6:6),VXPANH(-6:6),VXPDGM(-6:6)
      SAVE /PYINT8/,/PYINT9/
C...Local arrays.
      DIMENSION XPDFGM(-6:6),XPGA(-6:6), VXPGA(-6:6)
C...Charm and bottom masses (low to compensate for J/psi etc.).
      DATA PMC/1.3D0/, PMB/4.6D0/
C...alpha_em and alpha_em/(2*pi).
      DATA AEM/0.007297D0/, AEM2PI/0.0011614D0/
C...Lambda value for 4 flavours.
      DATA ALAM/0.20D0/
C...Mixture u/(u+d), = 0.5 for incoherent and = 0.8 for coherent sum.
      DATA FRACU/0.8D0/
C...VMD couplings f_V**2/(4*pi).
      DATA FRHO/2.20D0/, FOMEGA/23.6D0/, FPHI/18.4D0/
C...Masses for rho (=omega) and phi.
      DATA PMRHO/0.770D0/, PMPHI/1.020D0/
C...Number of points in integration for IP2=1.
      DATA NSTEP/100/
 
C...Reset output.
      F2GM=0D0
      DO 100 KFL=-6,6
        XPDFGM(KFL)=0D0
        XPVMD(KFL)=0D0
        XPANL(KFL)=0D0
        XPANH(KFL)=0D0
        XPBEH(KFL)=0D0
        XPDIR(KFL)=0D0
        VXPVMD(KFL)=0D0
        VXPANL(KFL)=0D0
        VXPANH(KFL)=0D0
        VXPDGM(KFL)=0D0
  100 CONTINUE
 
C...Set Q0 cut-off parameter as function of set used.
      IF(ISET.LE.2) THEN
        Q0=0.6D0
      ELSE
        Q0=2D0
      ENDIF
      Q02=Q0**2
 
C...Scale choice for off-shell photon; common factors.
      Q2A=Q2
      FACNOR=1D0
      IF(IP2.EQ.1) THEN
        P2MX=P2+Q02
        Q2A=Q2+P2*Q02/MAX(Q02,Q2)
        FACNOR=LOG(Q2/Q02)/NSTEP
      ELSEIF(IP2.EQ.2) THEN
        P2MX=MAX(P2,Q02)
      ELSEIF(IP2.EQ.3) THEN
        P2MX=P2+Q02
        Q2A=Q2+P2*Q02/MAX(Q02,Q2)
      ELSEIF(IP2.EQ.4) THEN
        P2MX=Q2*(Q02+P2)/(Q2+P2)*EXP(P2*(Q2-Q02)/
     &  ((Q2+P2)*(Q02+P2)))
      ELSEIF(IP2.EQ.5) THEN
        P2MXA=Q2*(Q02+P2)/(Q2+P2)*EXP(P2*(Q2-Q02)/
     &  ((Q2+P2)*(Q02+P2)))
        P2MX=Q0*SQRT(P2MXA)
        FACNOR=LOG(Q2/P2MXA)/LOG(Q2/P2MX)
      ELSEIF(IP2.EQ.6) THEN
        P2MX=Q2*(Q02+P2)/(Q2+P2)*EXP(P2*(Q2-Q02)/
     &  ((Q2+P2)*(Q02+P2)))
        P2MX=MAX(0D0,1D0-P2/Q2)*P2MX+MIN(1D0,P2/Q2)*MAX(P2,Q02)
      ELSE
        P2MXA=Q2*(Q02+P2)/(Q2+P2)*EXP(P2*(Q2-Q02)/
     &  ((Q2+P2)*(Q02+P2)))
        P2MX=Q0*SQRT(P2MXA)
        P2MXB=P2MX
        P2MX=MAX(0D0,1D0-P2/Q2)*P2MX+MIN(1D0,P2/Q2)*MAX(P2,Q02)
        P2MXB=MAX(0D0,1D0-P2/Q2)*P2MXB+MIN(1D0,P2/Q2)*P2MXA
        IF(ABS(Q2-Q02).GT.1D-6) THEN
          FACNOR=LOG(Q2/P2MXA)/LOG(Q2/P2MXB)
        ELSEIF(P2.LT.Q02) THEN
          FACNOR=Q02**3/(Q02+P2)/(Q02**2-P2**2/2D0)
        ELSE
          FACNOR=1D0
        ENDIF
      ENDIF
 
C...Call VMD parametrization for d quark and use to give rho, omega,
C...phi. Note dipole dampening for off-shell photon.
      CALL PYGVMD(ISET,1,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
      XFVAL=VXPGA(1)
      XPGA(1)=XPGA(2)
      XPGA(-1)=XPGA(-2)
      FACUD=AEM*(1D0/FRHO+1D0/FOMEGA)*(PMRHO**2/(PMRHO**2+P2))**2
      FACS=AEM*(1D0/FPHI)*(PMPHI**2/(PMPHI**2+P2))**2
      DO 110 KFL=-5,5
        XPVMD(KFL)=(FACUD+FACS)*XPGA(KFL)
  110 CONTINUE
      XPVMD(1)=XPVMD(1)+(1D0-FRACU)*FACUD*XFVAL
      XPVMD(2)=XPVMD(2)+FRACU*FACUD*XFVAL
      XPVMD(3)=XPVMD(3)+FACS*XFVAL
      XPVMD(-1)=XPVMD(-1)+(1D0-FRACU)*FACUD*XFVAL
      XPVMD(-2)=XPVMD(-2)+FRACU*FACUD*XFVAL
      XPVMD(-3)=XPVMD(-3)+FACS*XFVAL
      VXPVMD(1)=(1D0-FRACU)*FACUD*XFVAL
      VXPVMD(2)=FRACU*FACUD*XFVAL
      VXPVMD(3)=FACS*XFVAL
      VXPVMD(-1)=(1D0-FRACU)*FACUD*XFVAL
      VXPVMD(-2)=FRACU*FACUD*XFVAL
      VXPVMD(-3)=FACS*XFVAL
 
      IF(IP2.NE.1) THEN
C...Anomalous parametrizations for different strategies
C...for off-shell photons; except full integration.
 
C...Call anomalous parametrization for d + u + s.
        CALL PYGANO(-3,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
        DO 120 KFL=-5,5
          XPANL(KFL)=FACNOR*XPGA(KFL)
          VXPANL(KFL)=FACNOR*VXPGA(KFL)
  120   CONTINUE
 
C...Call anomalous parametrization for c and b.
        CALL PYGANO(4,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
        DO 130 KFL=-5,5
          XPANH(KFL)=FACNOR*XPGA(KFL)
          VXPANH(KFL)=FACNOR*VXPGA(KFL)
  130   CONTINUE
        CALL PYGANO(5,X,Q2A,P2MX,ALAM,XPGA,VXPGA)
        DO 140 KFL=-5,5
          XPANH(KFL)=XPANH(KFL)+FACNOR*XPGA(KFL)
          VXPANH(KFL)=VXPANH(KFL)+FACNOR*VXPGA(KFL)
  140   CONTINUE
 
      ELSE
C...Special option: loop over flavours and integrate over k2.
        DO 170 KF=1,5
          DO 160 ISTEP=1,NSTEP
            Q2STEP=Q02*(Q2/Q02)**((ISTEP-0.5D0)/NSTEP)
            IF((KF.EQ.4.AND.Q2STEP.LT.PMC**2).OR.
     &      (KF.EQ.5.AND.Q2STEP.LT.PMB**2)) GOTO 160
            CALL PYGVMD(0,KF,X,Q2,Q2STEP,ALAM,XPGA,VXPGA)
            FACQ=AEM2PI*(Q2STEP/(Q2STEP+P2))**2*FACNOR
            IF(MOD(KF,2).EQ.0) FACQ=FACQ*(8D0/9D0)
            IF(MOD(KF,2).EQ.1) FACQ=FACQ*(2D0/9D0)
            DO 150 KFL=-5,5
              IF(KF.LE.3) XPANL(KFL)=XPANL(KFL)+FACQ*XPGA(KFL)
              IF(KF.GE.4) XPANH(KFL)=XPANH(KFL)+FACQ*XPGA(KFL)
              IF(KF.LE.3) VXPANL(KFL)=VXPANL(KFL)+FACQ*VXPGA(KFL)
              IF(KF.GE.4) VXPANH(KFL)=VXPANH(KFL)+FACQ*VXPGA(KFL)
  150       CONTINUE
  160     CONTINUE
  170   CONTINUE
      ENDIF
 
C...Call Bethe-Heitler term expression for charm and bottom.
      CALL PYGBEH(4,X,Q2,P2,PMC**2,XPBH)
      XPBEH(4)=XPBH
      XPBEH(-4)=XPBH
      CALL PYGBEH(5,X,Q2,P2,PMB**2,XPBH)
      XPBEH(5)=XPBH
      XPBEH(-5)=XPBH
 
C...For MSbar subtraction call C^gamma term expression for d, u, s.
      IF(ISET.EQ.2.OR.ISET.EQ.4) THEN
        CALL PYGDIR(X,Q2,P2,Q02,XPGA)
        DO 180 KFL=-5,5
          XPDIR(KFL)=XPGA(KFL)
  180   CONTINUE
      ENDIF
 
C...Store result in output array.
      DO 190 KFL=-5,5
        CHSQ=1D0/9D0
        IF(IABS(KFL).EQ.2.OR.IABS(KFL).EQ.4) CHSQ=4D0/9D0
        XPF2=XPVMD(KFL)+XPANL(KFL)+XPBEH(KFL)+XPDIR(KFL)
        IF(KFL.NE.0) F2GM=F2GM+CHSQ*XPF2
        XPDFGM(KFL)=XPVMD(KFL)+XPANL(KFL)+XPANH(KFL)
        VXPDGM(KFL)=VXPVMD(KFL)+VXPANL(KFL)+VXPANH(KFL)
  190 CONTINUE
 
      RETURN
      END
