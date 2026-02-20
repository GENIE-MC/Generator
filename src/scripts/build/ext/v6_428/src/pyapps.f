
 
C*********************************************************************
 
C...PYAPPS
C...Uses approximate analytical formulae to determine the full set of
C...MSSM parameters from SUGRA input.
C...See M. Drees and S.P. Martin, hep-ph/9504124
 
      SUBROUTINE PYAPPS
 
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
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      SAVE /PYDAT1/,/PYDAT2/,/PYMSSM/

      WRITE(MSTU(11),*) '(PYAPPS:) approximate mSUGRA relations'//
     &' not intended for serious physics studies'
      IMSS(5)=0
      IMSS(8)=0
      XMT=PMAS(6,1)
      XMZ2=PMAS(23,1)**2
      XMW2=PMAS(24,1)**2
      TANB=RMSS(5)
      BETA=ATAN(TANB)
      XW=PARU(102)
      XMG=RMSS(1)
      XMG2=XMG*XMG
      XM0=RMSS(8)
      XM02=XM0*XM0
C...Temporary sign change for AT. Others unchanged.
      AT=-RMSS(16)
      RMSS(15)=RMSS(16)
      RMSS(17)=RMSS(16)
      SINB=TANB/SQRT(TANB**2+1D0)
      COSB=SINB/TANB
 
      DTERM=XMZ2*COS(2D0*BETA)
      XMER=SQRT(XM02+0.15D0*XMG2-XW*DTERM)
      XMEL=SQRT(XM02+0.52D0*XMG2-(0.5D0-XW)*DTERM)
      RMSS(6)=XMEL
      RMSS(7)=XMER
      XMUR=SQRT(PYRNMQ(2,2D0/3D0*XW*DTERM))
      XMDR=SQRT(PYRNMQ(3,-1D0/3D0*XW*DTERM))
      XMUL=SQRT(PYRNMQ(1,(0.5D0-2D0/3D0*XW)*DTERM))
      XMDL=SQRT(PYRNMQ(1,-(0.5D0-1D0/3D0*XW)*DTERM))
      DO 100 I=1,5,2
        PMAS(PYCOMP(KSUSY1+I),1)=XMDL
        PMAS(PYCOMP(KSUSY2+I),1)=XMDR
        PMAS(PYCOMP(KSUSY1+I+1),1)=XMUL
        PMAS(PYCOMP(KSUSY2+I+1),1)=XMUR
  100 CONTINUE
      XARG=XMEL**2-XMW2*ABS(COS(2D0*BETA))
      IF(XARG.LT.0D0) THEN
        WRITE(MSTU(11),*) ' SNEUTRINO MASS IS NEGATIVE'//
     &  ' FROM THE SUM RULE. '
        WRITE(MSTU(11),*) '  TRY A SMALLER VALUE OF TAN(BETA). '
        RETURN
      ELSE
        XARG=SQRT(XARG)
      ENDIF
      DO 110 I=11,15,2
        PMAS(PYCOMP(KSUSY1+I),1)=XMEL
        PMAS(PYCOMP(KSUSY2+I),1)=XMER
        PMAS(PYCOMP(KSUSY1+I+1),1)=XARG
        PMAS(PYCOMP(KSUSY2+I+1),1)=9999D0
  110 CONTINUE
      RMT=PYMRUN(6,PMAS(6,1)**2)
      XTOP=(RMT/150D0/SINB)**2*(.9D0*XM02+2.1D0*XMG2+
     &(1D0-(RMT/190D0/SINB)**3)*(.24D0*AT**2+AT*XMG))
      RMB=PYMRUN(5,PMAS(6,1)**2)
      XBOT=(RMB/150D0/COSB)**2*(.9D0*XM02+2.1D0*XMG2+
     &(1D0-(RMB/190D0/COSB)**3)*(.24D0*AT**2+AT*XMG))
      XTAU=1D-4/COSB**2*(XM02+0.15D0*XMG2+AT**2/3D0)
      ATP=AT*(1D0-(RMT/190D0/SINB)**2)+XMG*(3.47D0-1.9D0*(RMT/190D0/
     &SINB)**2)
      RMSS(16)=-ATP
      XMU2=-.5D0*XMZ2+(SINB**2*(XM02+.52D0*XMG2-XTOP)-
     &COSB**2*(XM02+.52D0*XMG2-XBOT-XTAU/3D0))/(COSB**2-SINB**2)
      XMA2=2D0*(XM02+.52D0*XMG2+XMU2)-XTOP-XBOT-XTAU/3D0
      XMU=SIGN(SQRT(XMU2),RMSS(4))
      RMSS(4)=XMU
      IF(XMA2.GT.0D0) THEN
        RMSS(19)=SQRT(XMA2)
      ELSE
        WRITE(MSTU(11),*) ' PYAPPS:: PSEUDOSCALAR MASS**2 < 0 '
        CALL PYSTOP(102)
      ENDIF
      ARG=XM02+0.15D0*XMG2-2D0*XTAU/3D0-XW*DTERM
      IF(ARG.GT.0D0) THEN
        RMSS(14)=SQRT(ARG)
      ELSE
        WRITE(MSTU(11),*) ' PYAPPS:: RIGHT STAU MASS**2 < 0 '
        CALL PYSTOP(102)
      ENDIF
      ARG=XM02+0.52D0*XMG2-XTAU/3D0-(0.5D0-XW)*DTERM
      IF(ARG.GT.0D0) THEN
        RMSS(13)=SQRT(ARG)
      ELSE
        WRITE(MSTU(11),*) ' PYAPPS::  LEFT STAU MASS**2 < 0 '
        CALL PYSTOP(102)
      ENDIF
      ARG=PYRNMQ(1,-(XBOT+XTOP)/3D0)
      IF(ARG.GT.0D0) THEN
        RMSS(10)=SQRT(ARG)
      ELSE
        RMSS(10)=-SQRT(-ARG)
      ENDIF
      ARG=PYRNMQ(2,-2D0*XTOP/3D0)
      IF(ARG.GT.0D0) THEN
        RMSS(12)=SQRT(ARG)
      ELSE
        RMSS(12)=-SQRT(-ARG)
      ENDIF
      ARG=PYRNMQ(3,-2D0*XBOT/3D0)
      IF(ARG.GT.0D0) THEN
        RMSS(11)=SQRT(ARG)
      ELSE
        RMSS(11)=-SQRT(-ARG)
      ENDIF
 
      RETURN
      END
