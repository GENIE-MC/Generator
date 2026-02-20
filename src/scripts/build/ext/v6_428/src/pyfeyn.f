 
C*********************************************************************
 
C...PYFEYN
C...Interface to FeynHiggs for MSSM Higgs sector.
C...Pythia6.402: Updated to FeynHiggs v.2.3.0+ w/ DOUBLE COMPLEX
C...P. Skands
 
      SUBROUTINE PYFEYN(IERR)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C...SUSY blocks
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
C...FeynHiggs variables
      DOUBLE PRECISION RMHIGG(4)
      DOUBLE COMPLEX SAEFF, UHIGGS(3,3)
      DOUBLE COMPLEX DMU,
     &     AE33, AU33, AD33, AE22, AU22, AD22, AE11, AU11, AD11,
     &     DM1, DM2, DM3
C...SLHA Common Block
      COMMON/PYLH3P/MODSEL(200),PARMIN(100),PAREXT(200),RMSOFT(0:100),
     &     AU(3,3),AD(3,3),AE(3,3)
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/,/PYLH3P/
 
      IERR=0
      CALL FHSETFLAGS(IERR,4,0,0,2,0,2,1,1)
      IF (IERR.NE.0) THEN
        CALL PYERRM(11,'(PYHGGM:) Caught error from FHSETFLAGS.'
     &       //'Will not use FeynHiggs for this run.')
        RETURN
      ENDIF
      Q=RMSOFT(0)
      DMB=PMAS(5,1)
      DMT=PMAS(6,1)
      DMZ=PMAS(23,1)
      DMW=PMAS(24,1)
      DMA=PMAS(36,1)
      DM1=RMSOFT(1)
      DM2=RMSOFT(2)
      DM3=RMSOFT(3)
      DTANB=RMSS(5)
      DMU=RMSS(4)
      DM3SL=RMSOFT(33)
      DM3SE=RMSOFT(36)
      DM3SQ=RMSOFT(43)
      DM3SU=RMSOFT(46)
      DM3SD=RMSOFT(49)
      DM2SL=RMSOFT(32)
      DM2SE=RMSOFT(35)
      DM2SQ=RMSOFT(42)
      DM2SU=RMSOFT(45)
      DM2SD=RMSOFT(48)
      DM1SL=RMSOFT(31)
      DM1SE=RMSOFT(34)
      DM1SQ=RMSOFT(41)
      DM1SU=RMSOFT(44)
      DM1SD=RMSOFT(47)
      AE33=AE(3,3)
      AE22=AE(2,2)
      AE11=AE(1,1)
      AU33=AU(3,3)
      AU22=AU(2,2)
      AU11=AU(1,1)
      AD33=AD(3,3)
      AD22=AD(2,2)
      AD11=AD(1,1)
      CALL FHSETPARA(IERR, 1D0, DMT, DMB, DMW, DMZ, DTANB,
     &     DMA,0D0, DM3SL, DM3SE, DM3SQ, DM3SU, DM3SD,
     &     DM2SL, DM2SE, DM2SQ, DM2SU, DM2SD,
     &     DM1SL, DM1SE, DM1SQ, DM1SU, DM1SD,DMU,
     &     AE33, AU33, AD33, AE22, AU22, AD22, AE11, AU11, AD11,
     &     DM1, DM2, DM3, 0D0, 0D0,Q,Q,Q)
      IF (IERR.NE.0) THEN
        CALL PYERRM(11,'(PYHGGM:) Caught error from FHSETPARA.'
     &       //' Will not use FeynHiggs for this run.')
        RETURN
      ENDIF
C...  Get Higgs masses & alpha_eff. (UHIGGS redundant here, only for CPV)
      SAEFF=0D0
      CALL FHHIGGSCORR(IERR, RMHIGG, SAEFF, UHIGGS)
      IF (IERR.NE.0) THEN
        CALL PYERRM(11,'(PYFEYN:) Caught error from FHHIG'//
     &       'GSCORR. Will not use FeynHiggs for this run.')
        RETURN
      ENDIF
      ALPHA = ASIN(DBLE(SAEFF))
      R=RMSS(18)/ALPHA
      IF (R.LT.0D0.OR.ABS(R).GT.1.2D0.OR.ABS(R).LT.0.8D0) THEN
        CALL PYERRM(1,'(PYFEYN:) Large corrections in Higgs sector.')
        WRITE(MSTU(11),*) '   Old Alpha:', RMSS(18)
        WRITE(MSTU(11),*) '   New Alpha:', ALPHA
      ENDIF
      IF (RMHIGG(1).LT.0.85D0*PMAS(25,1).OR.RMHIGG(1).GT.
     &       1.15D0*PMAS(25,1)) THEN
        CALL PYERRM(1,'(PYFEYN:) Large corrections in Higgs sector.')
        WRITE(MSTU(11),*) '   Old m(h0):', PMAS(25,1)
        WRITE(MSTU(11),*) '   New m(h0):', RMHIGG(1)
      ENDIF
      RMSS(18)=ALPHA
      PMAS(25,1)=RMHIGG(1)
      PMAS(35,1)=RMHIGG(2)
      PMAS(36,1)=RMHIGG(3)
      PMAS(37,1)=RMHIGG(4)
 
      RETURN
      END
