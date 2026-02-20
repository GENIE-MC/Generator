 
 
C*********************************************************************
 
C...PYGBEH
C...Evaluates the Bethe-Heitler cross section for heavy flavour
C...production.
C...Adapted from SaSgam library, authors G.A. Schuler and T. Sjostrand.
 
      SUBROUTINE PYGBEH(KF,X,Q2,P2,PM2,XPBH)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
 
C...Local data.
      DATA AEM2PI/0.0011614D0/
 
C...Reset output.
      XPBH=0D0
      SIGBH=0D0
 
C...Check kinematics limits.
      IF(X.GE.Q2/(4D0*PM2+Q2+P2)) RETURN
      W2=Q2*(1D0-X)/X-P2
      BETA2=1D0-4D0*PM2/W2
      IF(BETA2.LT.1D-10) RETURN
      BETA=SQRT(BETA2)
      RMQ=4D0*PM2/Q2
 
C...Simple case: P2 = 0.
      IF(P2.LT.1D-4) THEN
        IF(BETA.LT.0.99D0) THEN
          XBL=LOG((1D0+BETA)/(1D0-BETA))
        ELSE
          XBL=LOG((1D0+BETA)**2*W2/(4D0*PM2))
        ENDIF
        SIGBH=BETA*(8D0*X*(1D0-X)-1D0-RMQ*X*(1D0-X))+
     &  XBL*(X**2+(1D0-X)**2+RMQ*X*(1D0-3D0*X)-0.5D0*RMQ**2*X**2)
 
C...Complicated case: P2 > 0, based on approximation of
C...C.T. Hill and G.G. Ross, Nucl. Phys. B148 (1979) 373
      ELSE
        RPQ=1D0-4D0*X**2*P2/Q2
        IF(RPQ.GT.1D-10) THEN
          RPBE=SQRT(RPQ*BETA2)
          IF(RPBE.LT.0.99D0) THEN
            XBL=LOG((1D0+RPBE)/(1D0-RPBE))
            XBI=2D0*RPBE/(1D0-RPBE**2)
          ELSE
            RPBESN=4D0*PM2/W2+(4D0*X**2*P2/Q2)*BETA2
            XBL=LOG((1D0+RPBE)**2/RPBESN)
            XBI=2D0*RPBE/RPBESN
          ENDIF
          SIGBH=BETA*(6D0*X*(1D0-X)-1D0)+
     &    XBL*(X**2+(1D0-X)**2+RMQ*X*(1D0-3D0*X)-0.5D0*RMQ**2*X**2)+
     &    XBI*(2D0*X/Q2)*(PM2*X*(2D0-RMQ)-P2*X)
        ENDIF
      ENDIF
 
C...Multiply by charge-squared etc. to get parton distribution.
      CHSQ=1D0/9D0
      IF(IABS(KF).EQ.2.OR.IABS(KF).EQ.4) CHSQ=4D0/9D0
      XPBH=3D0*CHSQ*AEM2PI*X*SIGBH
 
      RETURN
      END
