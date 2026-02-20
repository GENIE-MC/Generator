 
C*********************************************************************
 
C...PYPCMP: Auxiliary to PYPDFU.
C...Giving the momentum integral of a companion quark, with its
C...partner at XS, using an approximate gluon density like (1-x)^NPOW/x.
C...The value corresponds to an unrescaled range between 0 and 1-XS.
 
      FUNCTION PYPCMP(XS,NPOW)
      IMPLICIT NONE
      DOUBLE PRECISION XS, PYPCMP
      INTEGER NPOW
      IF (XS.GE.1D0.OR.XS.LE.0D0) THEN
        PYPCMP=0D0
      ELSEIF (NPOW.LE.0) THEN
        PYPCMP=XS*(5D0+XS*(-9D0-2D0*XS*(-3D0+XS))+3D0*LOG(XS))
        PYPCMP=PYPCMP/((-1D0+XS)*(2D0+XS*(-1D0+2D0*XS)))
      ELSEIF (NPOW.EQ.1) THEN
        PYPCMP=-1D0-3D0*XS+(2D0*(-1D0+XS)**2*(1D0+XS+XS**2))
     &       /(2D0+XS**2*(XS-3D0)+3D0*XS*LOG(XS))
      ELSEIF (NPOW.EQ.2) THEN
        PYPCMP=XS*((1D0-XS)*(19D0+XS*(43D0+4D0*XS))
     &       +6D0*LOG(XS)*(1D0+6D0*XS+4D0*XS**2))
        PYPCMP=PYPCMP/(4D0*((XS-1D0)*(1D0+XS*(4D0+XS))
     &       -3D0*XS*LOG(XS)*(1+XS)))
      ELSEIF (NPOW.EQ.3) THEN
        PYPCMP=3D0*XS*((XS-1)*(7D0+XS*(28D0+13D0*XS))
     &       -2D0*LOG(XS)*(1D0+XS*(9D0+2D0*XS*(6D0+XS))))
        PYPCMP=PYPCMP/(4D0+27D0*XS-31D0*XS**3
     &       +6D0*XS*LOG(XS)*(3D0+2D0*XS*(3D0+XS)))
      ELSE
        PYPCMP=(-9D0*XS*(XS**2-1D0)*(5D0+XS*(24D0+XS))+12D0*XS*LOG(XS)
     &       *(1D0+2D0*XS)*(1D0+2D0*XS*(5D0+2D0*XS)))
        PYPCMP=PYPCMP/(8D0*(1D0+2D0*XS)*((XS-1D0)*(1D0+XS*(10D0+XS))
     &       -6D0*XS*LOG(XS)*(1D0+XS)))
      ENDIF
      RETURN
      END
