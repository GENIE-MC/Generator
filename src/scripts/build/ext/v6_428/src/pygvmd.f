 
C*********************************************************************
 
C...PYGVMD
C...Evaluates the VMD parton distributions of a photon,
C...evolved homogeneously from an initial scale P2 to Q2.
C...Does not include dipole suppression factor.
C...ISET is parton distribution set, see above;
C...additionally ISET=0 is used for the evolution of an anomalous photon
C...which branched at a scale P2 and then evolved homogeneously to Q2.
C...ALAM is the 4-flavour Lambda, which is automatically converted
C...to 3- and 5-flavour equivalents as needed.
C...Adapted from SaSgam library, authors G.A. Schuler and T. Sjostrand.
 
      SUBROUTINE PYGVMD(ISET,KF,X,Q2,P2,ALAM,XPGA,VXPGA)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Local arrays and data.
      DIMENSION XPGA(-6:6), VXPGA(-6:6)
      DATA PMC/1.3D0/, PMB/4.6D0/, AEM/0.007297D0/, AEM2PI/0.0011614D0/
 
C...Reset output.
      DO 100 KFL=-6,6
        XPGA(KFL)=0D0
        VXPGA(KFL)=0D0
  100 CONTINUE
      KFA=IABS(KF)
 
C...Calculate Lambda; protect against unphysical Q2 and P2 input.
      ALAM3=ALAM*(PMC/ALAM)**(2D0/27D0)
      ALAM5=ALAM*(ALAM/PMB)**(2D0/23D0)
      P2EFF=MAX(P2,1.2D0*ALAM3**2)
      IF(KFA.EQ.4) P2EFF=MAX(P2EFF,PMC**2)
      IF(KFA.EQ.5) P2EFF=MAX(P2EFF,PMB**2)
      Q2EFF=MAX(Q2,P2EFF)
 
C...Find number of flavours at lower and upper scale.
      NFP=4
      IF(P2EFF.LT.PMC**2) NFP=3
      IF(P2EFF.GT.PMB**2) NFP=5
      NFQ=4
      IF(Q2EFF.LT.PMC**2) NFQ=3
      IF(Q2EFF.GT.PMB**2) NFQ=5
 
C...Find s as sum of 3-, 4- and 5-flavour parts.
      S=0D0
      IF(NFP.EQ.3) THEN
        Q2DIV=PMC**2
        IF(NFQ.EQ.3) Q2DIV=Q2EFF
        S=S+(6D0/27D0)*LOG(LOG(Q2DIV/ALAM3**2)/LOG(P2EFF/ALAM3**2))
      ENDIF
      IF(NFP.LE.4.AND.NFQ.GE.4) THEN
        P2DIV=P2EFF
        IF(NFP.EQ.3) P2DIV=PMC**2
        Q2DIV=Q2EFF
        IF(NFQ.EQ.5) Q2DIV=PMB**2
        S=S+(6D0/25D0)*LOG(LOG(Q2DIV/ALAM**2)/LOG(P2DIV/ALAM**2))
      ENDIF
      IF(NFQ.EQ.5) THEN
        P2DIV=PMB**2
        IF(NFP.EQ.5) P2DIV=P2EFF
        S=S+(6D0/23D0)*LOG(LOG(Q2EFF/ALAM5**2)/LOG(P2DIV/ALAM5**2))
      ENDIF
 
C...Calculate frequent combinations of x and s.
      X1=1D0-X
      XL=-LOG(X)
      S2=S**2
      S3=S**3
      S4=S**4
 
C...Evaluate homogeneous anomalous parton distributions below or
C...above threshold.
      IF(ISET.EQ.0) THEN
        IF(Q2.LE.P2.OR.(KFA.EQ.4.AND.Q2.LT.PMC**2).OR.
     &  (KFA.EQ.5.AND.Q2.LT.PMB**2)) THEN
          XVAL = X * 1.5D0 * (X**2+X1**2)
          XGLU = 0D0
          XSEA = 0D0
        ELSE
          XVAL = (1.5D0/(1D0-0.197D0*S+4.33D0*S2)*X**2 +
     &    (1.5D0+2.10D0*S)/(1D0+3.29D0*S)*X1**2 +
     &    5.23D0*S/(1D0+1.17D0*S+19.9D0*S3)*X*X1) *
     &    X**(1D0/(1D0+1.5D0*S)) * (1D0-X**2)**(2.667D0*S)
          XGLU = 4D0*S/(1D0+4.76D0*S+15.2D0*S2+29.3D0*S4) *
     &    X**(-2.03D0*S/(1D0+2.44D0*S)) * (X1*XL)**(1.333D0*S) *
     &    ((4D0*X**2+7D0*X+4D0)*X1/3D0 - 2D0*X*(1D0+X)*XL)
          XSEA = S2/(1D0+4.54D0*S+8.19D0*S2+8.05D0*S3) *
     &    X**(-1.54D0*S/(1D0+1.29D0*S)) * X1**(2.667D0*S) *
     &    ((8D0-73D0*X+62D0*X**2)*X1/9D0 + (3D0-8D0*X**2/3D0)*X*XL +
     &    (2D0*X-1D0)*X*XL**2)
        ENDIF
 
C...Evaluate set 1D parton distributions below or above threshold.
      ELSEIF(ISET.EQ.1) THEN
        IF(Q2.LE.P2.OR.(KFA.EQ.4.AND.Q2.LT.PMC**2).OR.
     &  (KFA.EQ.5.AND.Q2.LT.PMB**2)) THEN
          XVAL = 1.294D0 * X**0.80D0 * X1**0.76D0
          XGLU = 1.273D0 * X**0.40D0 * X1**1.76D0
          XSEA = 0.100D0 * X1**3.76D0
        ELSE
          XVAL = 1.294D0/(1D0+0.252D0*S+3.079D0*S2) *
     &    X**(0.80D0-0.13D0*S) * X1**(0.76D0+0.667D0*S) * XL**(2D0*S)
          XGLU = 7.90D0*S/(1D0+5.50D0*S) * EXP(-5.16D0*S) *
     &    X**(-1.90D0*S/(1D0+3.60D0*S)) * X1**1.30D0 *
     &    XL**(0.50D0+3D0*S) + 1.273D0 * EXP(-10D0*S) *
     &    X**0.40D0 * X1**(1.76D0+3D0*S)
          XSEA = (0.1D0-0.397D0*S2+1.121D0*S3)/
     &    (1D0+5.61D0*S2+5.26D0*S3) * X**(-7.32D0*S2/(1D0+10.3D0*S2)) *
     &    X1**((3.76D0+15D0*S+12D0*S2)/(1D0+4D0*S))
          XSEA0 = 0.100D0 * X1**3.76D0
        ENDIF
 
C...Evaluate set 1M parton distributions below or above threshold.
      ELSEIF(ISET.EQ.2) THEN
        IF(Q2.LE.P2.OR.(KFA.EQ.4.AND.Q2.LT.PMC**2).OR.
     &  (KFA.EQ.5.AND.Q2.LT.PMB**2)) THEN
          XVAL = 0.8477D0 * X**0.51D0 * X1**1.37D0
          XGLU = 3.42D0 * X**0.255D0 * X1**2.37D0
          XSEA = 0D0
        ELSE
          XVAL = 0.8477D0/(1D0+1.37D0*S+2.18D0*S2+3.73D0*S3) *
     &    X**(0.51D0+0.21D0*S) * X1**1.37D0 * XL**(2.667D0*S)
          XGLU = 24D0*S/(1D0+9.6D0*S+0.92D0*S2+14.34D0*S3) *
     &    EXP(-5.94D0*S) * X**((-0.013D0-1.80D0*S)/(1D0+3.14D0*S)) *
     &    X1**(2.37D0+0.4D0*S) * XL**(0.32D0+3.6D0*S) + 3.42D0 *
     &    EXP(-12D0*S) * X**0.255D0 * X1**(2.37D0+3D0*S)
          XSEA = 0.842D0*S/(1D0+21.3D0*S-33.2D0*S2+229D0*S3) *
     &    X**((0.13D0-2.90D0*S)/(1D0+5.44D0*S)) * X1**(3.45D0+0.5D0*S) *
     &    XL**(2.8D0*S)
          XSEA0 = 0D0
        ENDIF
 
C...Evaluate set 2D parton distributions below or above threshold.
      ELSEIF(ISET.EQ.3) THEN
        IF(Q2.LE.P2.OR.(KFA.EQ.4.AND.Q2.LT.PMC**2).OR.
     &  (KFA.EQ.5.AND.Q2.LT.PMB**2)) THEN
          XVAL = X**0.46D0 * X1**0.64D0 + 0.76D0 * X
          XGLU = 1.925D0 * X1**2
          XSEA = 0.242D0 * X1**4
        ELSE
          XVAL = (1D0+0.186D0*S)/(1D0-0.209D0*S+1.495D0*S2) *
     &    X**(0.46D0+0.25D0*S) *
     &    X1**((0.64D0+0.14D0*S+5D0*S2)/(1D0+S)) * XL**(1.9D0*S) +
     &    (0.76D0+0.4D0*S) * X * X1**(2.667D0*S)
          XGLU = (1.925D0+5.55D0*S+147D0*S2)/(1D0-3.59D0*S+3.32D0*S2) *
     &    EXP(-18.67D0*S) *
     &    X**((-5.81D0*S-5.34D0*S2)/(1D0+29D0*S-4.26D0*S2))
     &    * X1**((2D0-5.9D0*S)/(1D0+1.7D0*S)) *
     &    XL**(9.3D0*S/(1D0+1.7D0*S))
          XSEA = (0.242D0-0.252D0*S+1.19D0*S2)/
     &    (1D0-0.607D0*S+21.95D0*S2) *
     &    X**(-12.1D0*S2/(1D0+2.62D0*S+16.7D0*S2)) * X1**4 * XL**S
          XSEA0 = 0.242D0 * X1**4
        ENDIF
 
C...Evaluate set 2M parton distributions below or above threshold.
      ELSEIF(ISET.EQ.4) THEN
        IF(Q2.LE.P2.OR.(KFA.EQ.4.AND.Q2.LT.PMC**2).OR.
     &  (KFA.EQ.5.AND.Q2.LT.PMB**2)) THEN
          XVAL = 1.168D0 * X**0.50D0 * X1**2.60D0 + 0.965D0 * X
          XGLU = 1.808D0 * X1**2
          XSEA = 0.209D0 * X1**4
        ELSE
          XVAL = (1.168D0+1.771D0*S+29.35D0*S2) * EXP(-5.776D0*S) *
     &    X**((0.5D0+0.208D0*S)/(1D0-0.794D0*S+1.516D0*S2)) *
     &    X1**((2.6D0+7.6D0*S)/(1D0+5D0*S)) *
     &    XL**(5.15D0*S/(1D0+2D0*S)) +
     &    (0.965D0+22.35D0*S)/(1D0+18.4D0*S) * X * X1**(2.667D0*S)
          XGLU = (1.808D0+29.9D0*S)/(1D0+26.4D0*S) * EXP(-5.28D0*S) *
     &    X**((-5.35D0*S-10.11D0*S2)/(1D0+31.71D0*S)) *
     &    X1**((2D0-7.3D0*S+4D0*S2)/(1D0+2.5D0*S)) *
     &    XL**(10.9D0*S/(1D0+2.5D0*S))
          XSEA = (0.209D0+0.644D0*S2)/(1D0+0.319D0*S+17.6D0*S2) *
     &    X**((-0.373D0*S-7.71D0*S2)/(1D0+0.815D0*S+11.0D0*S2)) *
     &    X1**(4D0+S) * XL**(0.45D0*S)
          XSEA0 = 0.209D0 * X1**4
        ENDIF
      ENDIF
 
C...Threshold factors for c and b sea.
      SLL=LOG(LOG(Q2EFF/ALAM**2)/LOG(P2EFF/ALAM**2))
      XCHM=0D0
      IF(Q2.GT.PMC**2.AND.Q2.GT.1.001D0*P2EFF) THEN
        SCH=MAX(0D0,LOG(LOG(PMC**2/ALAM**2)/LOG(P2EFF/ALAM**2)))
        IF(ISET.EQ.0) THEN
          XCHM=XSEA*(1D0-(SCH/SLL)**2)
        ELSE
          XCHM=MAX(0D0,XSEA-XSEA0*X1**(2.667D0*S))*(1D0-SCH/SLL)
        ENDIF
      ENDIF
      XBOT=0D0
      IF(Q2.GT.PMB**2.AND.Q2.GT.1.001D0*P2EFF) THEN
        SBT=MAX(0D0,LOG(LOG(PMB**2/ALAM**2)/LOG(P2EFF/ALAM**2)))
        IF(ISET.EQ.0) THEN
          XBOT=XSEA*(1D0-(SBT/SLL)**2)
        ELSE
          XBOT=MAX(0D0,XSEA-XSEA0*X1**(2.667D0*S))*(1D0-SBT/SLL)
        ENDIF
      ENDIF
 
C...Fill parton distributions.
      XPGA(0)=XGLU
      XPGA(1)=XSEA
      XPGA(2)=XSEA
      XPGA(3)=XSEA
      XPGA(4)=XCHM
      XPGA(5)=XBOT
      XPGA(KFA)=XPGA(KFA)+XVAL
      DO 110 KFL=1,5
        XPGA(-KFL)=XPGA(KFL)
  110 CONTINUE
      VXPGA(KFA)=XVAL
      VXPGA(-KFA)=XVAL
 
      RETURN
      END
