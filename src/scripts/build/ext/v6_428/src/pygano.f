 
C*********************************************************************
 
C...PYGANO
C...Evaluates the parton distributions of the anomalous photon,
C...inhomogeneously evolved from a scale P2 (where it vanishes) to Q2.
C...KF=0 gives the sum over (up to) 5 flavours,
C...KF<0 limits to flavours up to abs(KF),
C...KF>0 is for flavour KF only.
C...ALAM is the 4-flavour Lambda, which is automatically converted
C...to 3- and 5-flavour equivalents as needed.
C...Adapted from SaSgam library, authors G.A. Schuler and T. Sjostrand.
 
      SUBROUTINE PYGANO(KF,X,Q2,P2,ALAM,XPGA,VXPGA)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Local arrays and data.
      DIMENSION XPGA(-6:6), VXPGA(-6:6), ALAMSQ(3:5)
      DATA PMC/1.3D0/, PMB/4.6D0/, AEM/0.007297D0/, AEM2PI/0.0011614D0/
 
C...Reset output.
      DO 100 KFL=-6,6
        XPGA(KFL)=0D0
        VXPGA(KFL)=0D0
  100 CONTINUE
      IF(Q2.LE.P2) RETURN
      KFA=IABS(KF)
 
C...Calculate Lambda; protect against unphysical Q2 and P2 input.
      ALAMSQ(3)=(ALAM*(PMC/ALAM)**(2D0/27D0))**2
      ALAMSQ(4)=ALAM**2
      ALAMSQ(5)=(ALAM*(ALAM/PMB)**(2D0/23D0))**2
      P2EFF=MAX(P2,1.2D0*ALAMSQ(3))
      IF(KF.EQ.4) P2EFF=MAX(P2EFF,PMC**2)
      IF(KF.EQ.5) P2EFF=MAX(P2EFF,PMB**2)
      Q2EFF=MAX(Q2,P2EFF)
      XL=-LOG(X)
 
C...Find number of flavours at lower and upper scale.
      NFP=4
      IF(P2EFF.LT.PMC**2) NFP=3
      IF(P2EFF.GT.PMB**2) NFP=5
      NFQ=4
      IF(Q2EFF.LT.PMC**2) NFQ=3
      IF(Q2EFF.GT.PMB**2) NFQ=5
 
C...Define range of flavour loop.
      IF(KF.EQ.0) THEN
        KFLMN=1
        KFLMX=5
      ELSEIF(KF.LT.0) THEN
        KFLMN=1
        KFLMX=KFA
      ELSE
        KFLMN=KFA
        KFLMX=KFA
      ENDIF
 
C...Loop over flavours the photon can branch into.
      DO 110 KFL=KFLMN,KFLMX
 
C...Light flavours: calculate t range and (approximate) s range.
        IF(KFL.LE.3.AND.(KFL.EQ.1.OR.KFL.EQ.KF)) THEN
          TDIFF=LOG(Q2EFF/P2EFF)
          S=(6D0/(33D0-2D0*NFQ))*LOG(LOG(Q2EFF/ALAMSQ(NFQ))/
     &    LOG(P2EFF/ALAMSQ(NFQ)))
          IF(NFQ.GT.NFP) THEN
            Q2DIV=PMB**2
            IF(NFQ.EQ.4) Q2DIV=PMC**2
            SNFQ=(6D0/(33D0-2D0*NFQ))*LOG(LOG(Q2DIV/ALAMSQ(NFQ))/
     &      LOG(P2EFF/ALAMSQ(NFQ)))
            SNFP=(6D0/(33D0-2D0*(NFQ-1)))*LOG(LOG(Q2DIV/ALAMSQ(NFQ-1))/
     &      LOG(P2EFF/ALAMSQ(NFQ-1)))
            S=S+(LOG(Q2DIV/P2EFF)/LOG(Q2EFF/P2EFF))*(SNFP-SNFQ)
          ENDIF
          IF(NFQ.EQ.5.AND.NFP.EQ.3) THEN
            Q2DIV=PMC**2
            SNF4=(6D0/(33D0-2D0*4))*LOG(LOG(Q2DIV/ALAMSQ(4))/
     &      LOG(P2EFF/ALAMSQ(4)))
            SNF3=(6D0/(33D0-2D0*3))*LOG(LOG(Q2DIV/ALAMSQ(3))/
     &      LOG(P2EFF/ALAMSQ(3)))
            S=S+(LOG(Q2DIV/P2EFF)/LOG(Q2EFF/P2EFF))*(SNF3-SNF4)
          ENDIF
 
C...u and s quark do not need a separate treatment when d has been done.
        ELSEIF(KFL.EQ.2.OR.KFL.EQ.3) THEN
 
C...Charm: as above, but only include range above c threshold.
        ELSEIF(KFL.EQ.4) THEN
          IF(Q2.LE.PMC**2) GOTO 110
          P2EFF=MAX(P2EFF,PMC**2)
          Q2EFF=MAX(Q2EFF,P2EFF)
          TDIFF=LOG(Q2EFF/P2EFF)
          S=(6D0/(33D0-2D0*NFQ))*LOG(LOG(Q2EFF/ALAMSQ(NFQ))/
     &    LOG(P2EFF/ALAMSQ(NFQ)))
          IF(NFQ.EQ.5.AND.NFP.EQ.4) THEN
            Q2DIV=PMB**2
            SNFQ=(6D0/(33D0-2D0*NFQ))*LOG(LOG(Q2DIV/ALAMSQ(NFQ))/
     &      LOG(P2EFF/ALAMSQ(NFQ)))
            SNFP=(6D0/(33D0-2D0*(NFQ-1)))*LOG(LOG(Q2DIV/ALAMSQ(NFQ-1))/
     &      LOG(P2EFF/ALAMSQ(NFQ-1)))
            S=S+(LOG(Q2DIV/P2EFF)/LOG(Q2EFF/P2EFF))*(SNFP-SNFQ)
          ENDIF
 
C...Bottom: as above, but only include range above b threshold.
        ELSEIF(KFL.EQ.5) THEN
          IF(Q2.LE.PMB**2) GOTO 110
          P2EFF=MAX(P2EFF,PMB**2)
          Q2EFF=MAX(Q2,P2EFF)
          TDIFF=LOG(Q2EFF/P2EFF)
          S=(6D0/(33D0-2D0*NFQ))*LOG(LOG(Q2EFF/ALAMSQ(NFQ))/
     &    LOG(P2EFF/ALAMSQ(NFQ)))
        ENDIF
 
C...Evaluate flavour-dependent prefactor (charge^2 etc.).
        CHSQ=1D0/9D0
        IF(KFL.EQ.2.OR.KFL.EQ.4) CHSQ=4D0/9D0
        FAC=AEM2PI*2D0*CHSQ*TDIFF
 
C...Evaluate parton distributions (normalized to unit momentum sum).
        IF(KFL.EQ.1.OR.KFL.EQ.4.OR.KFL.EQ.5.OR.KFL.EQ.KF) THEN
          XVAL= ((1.5D0+2.49D0*S+26.9D0*S**2)/(1D0+32.3D0*S**2)*X**2 +
     &    (1.5D0-0.49D0*S+7.83D0*S**2)/(1D0+7.68D0*S**2)*(1D0-X)**2 +
     &    1.5D0*S/(1D0-3.2D0*S+7D0*S**2)*X*(1D0-X)) *
     &    X**(1D0/(1D0+0.58D0*S)) * (1D0-X**2)**(2.5D0*S/(1D0+10D0*S))
          XGLU= 2D0*S/(1D0+4D0*S+7D0*S**2) *
     &    X**(-1.67D0*S/(1D0+2D0*S)) * (1D0-X**2)**(1.2D0*S) *
     &    ((4D0*X**2+7D0*X+4D0)*(1D0-X)/3D0 - 2D0*X*(1D0+X)*XL)
          XSEA= 0.333D0*S**2/(1D0+4.90D0*S+4.69D0*S**2+21.4D0*S**3) *
     &    X**(-1.18D0*S/(1D0+1.22D0*S)) * (1D0-X)**(1.2D0*S) *
     &    ((8D0-73D0*X+62D0*X**2)*(1D0-X)/9D0 +
     &    (3D0-8D0*X**2/3D0)*X*XL + (2D0*X-1D0)*X*XL**2)
 
C...Threshold factors for c and b sea.
          SLL=LOG(LOG(Q2EFF/ALAM**2)/LOG(P2EFF/ALAM**2))
          XCHM=0D0
          IF(Q2.GT.PMC**2.AND.Q2.GT.1.001D0*P2EFF) THEN
            SCH=MAX(0D0,LOG(LOG(PMC**2/ALAM**2)/LOG(P2EFF/ALAM**2)))
            XCHM=XSEA*(1D0-(SCH/SLL)**3)
          ENDIF
          XBOT=0D0
          IF(Q2.GT.PMB**2.AND.Q2.GT.1.001D0*P2EFF) THEN
            SBT=MAX(0D0,LOG(LOG(PMB**2/ALAM**2)/LOG(P2EFF/ALAM**2)))
            XBOT=XSEA*(1D0-(SBT/SLL)**3)
          ENDIF
        ENDIF
 
C...Add contribution of each valence flavour.
        XPGA(0)=XPGA(0)+FAC*XGLU
        XPGA(1)=XPGA(1)+FAC*XSEA
        XPGA(2)=XPGA(2)+FAC*XSEA
        XPGA(3)=XPGA(3)+FAC*XSEA
        XPGA(4)=XPGA(4)+FAC*XCHM
        XPGA(5)=XPGA(5)+FAC*XBOT
        XPGA(KFL)=XPGA(KFL)+FAC*XVAL
        VXPGA(KFL)=VXPGA(KFL)+FAC*XVAL
  110 CONTINUE
      DO 120 KFL=1,5
        XPGA(-KFL)=XPGA(KFL)
        VXPGA(-KFL)=VXPGA(KFL)
  120 CONTINUE
 
      RETURN
      END
