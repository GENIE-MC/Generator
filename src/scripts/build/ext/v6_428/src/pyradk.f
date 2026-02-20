 
C*********************************************************************
 
C...PYRADK
C...Generates initial state photon radiation.
 
      SUBROUTINE PYRADK(ECM,MK,PAK,THEK,PHIK,ALPK)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYDAT1/
 
C...Function: cumulative hard photon spectrum in QFD case.
      FXK(XX)=2D0*LOG(XX)+PARJ(161)*LOG(1D0-XX)+PARJ(162)*XX+
     &PARJ(163)*LOG((XX-SZM)**2+SZW**2)+PARJ(164)*ATAN((XX-SZM)/SZW)
 
C...Determine whether radiative photon or not.
      MK=0
      PAK=0D0
      IF(PARJ(160).LT.PYR(0)) RETURN
      MK=1
 
C...Photon energy range. Find photon momentum in QED case.
      XKL=PARJ(135)
      XKU=MIN(PARJ(136),1D0-(2D0*PARJ(127)/ECM)**2)
      IF(MSTJ(102).LE.1) THEN
  100   XK=1D0/(1D0+(1D0/XKL-1D0)*((1D0/XKU-1D0)/(1D0/XKL-1D0))**PYR(0))
        IF(1D0+(1D0-XK)**2.LT.2D0*PYR(0)) GOTO 100
 
C...Ditto in QFD case, by numerical inversion of integrated spectrum.
      ELSE
        SZM=1D0-(PARJ(123)/ECM)**2
        SZW=PARJ(123)*PARJ(124)/ECM**2
        FXKL=FXK(XKL)
        FXKU=FXK(XKU)
        FXKD=1D-4*(FXKU-FXKL)
        FXKR=FXKL+PYR(0)*(FXKU-FXKL)
        NXK=0
  110   NXK=NXK+1
        XK=0.5D0*(XKL+XKU)
        FXKV=FXK(XK)
        IF(FXKV.GT.FXKR) THEN
          XKU=XK
          FXKU=FXKV
        ELSE
          XKL=XK
          FXKL=FXKV
        ENDIF
        IF(NXK.LT.15.AND.FXKU-FXKL.GT.FXKD) GOTO 110
        XK=XKL+(XKU-XKL)*(FXKR-FXKL)/(FXKU-FXKL)
      ENDIF
      PAK=0.5D0*ECM*XK
 
C...Photon polar and azimuthal angle.
      PME=2D0*(PYMASS(11)/ECM)**2
  120 CTHM=PME*(2D0/PME)**PYR(0)
      IF(1D0-(XK**2*CTHM*(1D0-0.5D0*CTHM)+2D0*(1D0-XK)*PME/MAX(PME,
     &CTHM*(1D0-0.5D0*CTHM)))/(1D0+(1D0-XK)**2).LT.PYR(0)) GOTO 120
      CTHE=1D0-CTHM
      IF(PYR(0).GT.0.5D0) CTHE=-CTHE
      STHE=SQRT(MAX(0D0,(CTHM-PME)*(2D0-CTHM)))
      THEK=PYANGL(CTHE,STHE)
      PHIK=PARU(2)*PYR(0)
 
C...Rotation angle for hadronic system.
      SGN=1D0
      IF(0.5D0*(2D0-XK*(1D0-CTHE))**2/((2D0-XK)**2+(XK*CTHE)**2).GT.
     &PYR(0)) SGN=-1D0
      ALPK=ASIN(SGN*STHE*(XK-SGN*(2D0*SQRT(1D0-XK)-2D0+XK)*CTHE)/
     &(2D0-XK*(1D0-SGN*CTHE)))
 
      RETURN
      END
