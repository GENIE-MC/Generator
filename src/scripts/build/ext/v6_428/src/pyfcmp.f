 
C*********************************************************************
 
C...PYFCMP: Auxiliary to PYPDFU and PYPTIS.
C...Giving the x*f pdf of a companion quark, with its partner at XS,
C...using an approximate gluon density like (1-X)^NPOW/X. The value
C...corresponds to an unrescaled range between 0 and 1-X.
 
      FUNCTION PYFCMP(XC,XS,NPOW)
      IMPLICIT NONE
      DOUBLE PRECISION XC, XS, Y, PYFCMP,FAC
      INTEGER NPOW
 
      PYFCMP=0D0
C...Parent gluon momentum fraction
      Y=XC+XS
      IF (Y.GE.1D0) RETURN
C...Common factor (includes factor XC, since PYFCMP=x*f)
      FAC=3D0*XC*XS*(XC**2+XS**2)/(Y**4)
C...Store normalized companion x*f distribution.
      IF (NPOW.LE.0) THEN
        PYFCMP=FAC/(2D0-XS*(3D0-XS*(3D0-2D0*XS)))
      ELSEIF (NPOW.EQ.1) THEN
        PYFCMP=FAC*(1D0-Y)/(2D0+XS**2*(-3D0+XS)+3D0*XS*LOG(XS))
      ELSEIF (NPOW.EQ.2) THEN
        PYFCMP=FAC*(1D0-Y)**2/(2D0*((1D0-XS)*(1D0+XS*(4D0+XS))
     &       +3D0*XS*(1D0+XS)*LOG(XS)))
      ELSEIF (NPOW.EQ.3) THEN
        PYFCMP=FAC*(1D0-Y)**3*2D0/(4D0+27D0*XS-31D0*XS**3
     &       +6D0*XS*LOG(XS)*(3D0+2D0*XS*(3D0+XS)))
      ELSEIF (NPOW.GE.4) THEN
        PYFCMP=FAC*(1D0-Y)**4/(2D0*(1D0+2D0*XS)*((1D0-XS)*(1D0+
     &       XS*(10D0+XS))+6D0*XS*LOG(XS)*(1D0+XS)))
      ENDIF
      RETURN
      END
