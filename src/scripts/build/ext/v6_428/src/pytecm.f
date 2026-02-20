 
 
C*********************************************************************
 
C...PYTECM
C...Finds the s-hat dependent eigenvalues of the inverse propagator
C...matrix for gamma, Z, techni-rho, and techni-omega to optimize the
C...phase space generation.  Extended to include techni-a meson, and
C...to return the width.
 
      SUBROUTINE PYTECM(SMIN,SMOU,WIDO,IOPT)
 
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
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYTCSM/ITCM(0:99),RTCM(0:99)
      SAVE /PYDAT1/,/PYDAT2/,/PYPARS/,/PYTCSM/
 
C...Local variables.
      DOUBLE PRECISION AR(5,5),WR(5),ZR(5,5),ZI(5,5),WORK(12,12),
     &AT(5,5),WI(5),FV1(5),FV2(5),FV3(5),SH,AEM,TANW,CT2W,QUPD,ALPRHT,
     &FAR,FAO,FZR,FZO,SHR,R1,R2,S1,S2,WDTP(0:400),WDTE(0:400,0:5),WX(5)
      INTEGER i,j,ierr

      SH=SMIN
      SHR=SQRT(SH)
      AEM=PYALEM(SH)
 
      SINW=MIN(SQRT(PARU(102)),1D0)
      COSW=SQRT(1D0-SINW**2)
      TANW=SINW/COSW
      CT2W=(1D0-2D0*PARU(102))/(2D0*PARU(102)/TANW)
      QUPD=2D0*RTCM(2)-1D0

      ALPRHT=2.16D0*(3D0/DBLE(ITCM(1)))
      FAR=SQRT(AEM/ALPRHT)
      FAO=FAR*QUPD
      FZR=FAR*CT2W
      FZO=-FAO*TANW
      FZX=-FAR/RTCM(47)/(2D0*SINW*COSW)
      FWR=FAR/(2D0*SINW)
      FWX=-FWR/RTCM(47)

      DO 110 I=1,5
        DO 100 J=1,5
          AT(I,J)=0D0
  100   CONTINUE
  110 CONTINUE

C...NC
      IF(IOPT.EQ.1) THEN
        AR(1,1) = SH
        AR(2,2) = SH-PMAS(23,1)**2
        AR(3,3) = SH-PMAS(PYCOMP(KTECHN+113),1)**2
        AR(4,4) = SH-PMAS(PYCOMP(KTECHN+223),1)**2
        AR(5,5) = SH-PMAS(PYCOMP(KTECHN+115),1)**2
        AR(1,2) = 0D0
        AR(2,1) = 0D0
        AR(1,3) = SH*FAR
        AR(3,1) = AR(1,3)
        AR(1,4) = SH*FAO
        AR(4,1) = AR(1,4)
        AR(2,3) = SH*FZR
        AR(3,2) = AR(2,3)
        AR(2,4) = SH*FZO
        AR(4,2) = AR(2,4)
        AR(3,4) = 0D0
        AR(4,3) = 0D0
        AR(2,5) = SH*FZX
        AR(5,2) = AR(2,5)
        AR(1,5) = 0D0
        AR(5,1) = AR(1,5)
        AR(3,5) = 0D0
        AR(5,3) = AR(3,5)
        AR(4,5) = 0D0
        AR(5,4) = AR(4,5)
        CALL PYWIDT(23,SH,WDTP,WDTE)
        AT(2,2) = WDTP(0)*SHR
        CALL PYWIDT(KTECHN+113,SH,WDTP,WDTE)
        AT(3,3) = WDTP(0)*SHR
        CALL PYWIDT(KTECHN+223,SH,WDTP,WDTE)
        AT(4,4) = WDTP(0)*SHR
        CALL PYWIDT(KTECHN+115,SH,WDTP,WDTE)
        AT(5,5) = WDTP(0)*SHR
        IDIM=5
C...CC
      ELSE
        AR(1,1) = SH-PMAS(24,1)**2
        AR(2,2) = SH-PMAS(PYCOMP(KTECHN+213),1)**2
        AR(3,3) = SH-PMAS(PYCOMP(KTECHN+215),1)**2
        AR(1,2) = SH*FWR
        AR(2,1) = AR(1,2)
        AR(1,3) = SH*FWX
        AR(3,1) = AR(1,3)
        AR(2,3) = 0D0
        AR(3,2) = 0D0
        CALL PYWIDT(24,SH,WDTP,WDTE)
        AT(1,1) = WDTP(0)*SHR
        CALL PYWIDT(KTECHN+213,SH,WDTP,WDTE)
        AT(2,2) = WDTP(0)*SHR
        CALL PYWIDT(KTECHN+215,SH,WDTP,WDTE)
        AT(3,3) = WDTP(0)*SHR
        IDIM=3
      ENDIF
      CALL PYEICG(IDIM,IDIM,AR,AT,WR,WI,0,ZR,ZI,FV1,FV2,FV3,IERR)

      IMIN=1
      SXMN=1D20
      DO 120 I=1,IDIM
        WX(I)=SQRT(ABS(SH-WR(I)))
        WR(I)=ABS(WR(I))
        IF(WR(I).LT.SXMN) THEN
          SXMN=WR(I)
          IMIN=I
        ENDIF
  120 CONTINUE
      SMOU=WX(IMIN)**2
      WIDO=WI(IMIN)/SHR

      RETURN
      END
