 
C*********************************************************************
 
C...PYJURF
C...From three given input vectors in PJU the boost VJU from
C...the "lab frame" to the junction rest frame is constructed.
 
      SUBROUTINE PYJURF(PJU,VJU)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
C...Input, output and local arrays.
      DIMENSION PJU(3,5),VJU(5),PSUM(5),A(3,3),PENEW(3),PCM(5,5)
      DATA TWOPI/6.283186D0/
 
C...Calculate masses and other invariants.
      DO 100 J=1,4
        PSUM(J)=PJU(1,J)+PJU(2,J)+PJU(3,J)
  100 CONTINUE
      PSUM2=PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-PSUM(3)**2
      PSUM(5)=SQRT(PSUM2)
      DO 120 I=1,3
        DO 110 J=1,3
          A(I,J)=PJU(I,4)*PJU(J,4)-PJU(I,1)*PJU(J,1)-
     &    PJU(I,2)*PJU(J,2)-PJU(I,3)*PJU(J,3)
  110   CONTINUE
  120 CONTINUE
 
C...Pick I to be most massive parton and J to be the one closest to I.
      ITRY=0
      I=1
      IF(A(2,2).GT.A(1,1)) I=2
      IF(A(3,3).GT.MAX(A(1,1),A(2,2))) I=3
  130 ITRY=ITRY+1
      J=1+MOD(I,3)
      K=1+MOD(J,3)
      IF(A(I,K)**2*A(J,J).LT.A(I,J)**2*A(K,K)) THEN
        K=1+MOD(I,3)
        J=1+MOD(K,3)
      ENDIF
      PMI2=A(I,I)
      PMJ2=A(J,J)
      PMK2=A(K,K)
      AIJ=A(I,J)
      AIK=A(I,K)
      AJK=A(J,K)
 
C...Trivial find new parton energies if all three partons are massless.
      IF(PMI2.LT.1D-4) THEN
        PEI=SQRT(2D0*AIK*AIJ/(3D0*AJK))
        PEJ=SQRT(2D0*AJK*AIJ/(3D0*AIK))
        PEK=SQRT(2D0*AIK*AJK/(3D0*AIJ))
 
C...Else find momentum range for parton I and values at extremes.
      ELSE
        PAIMIN=0D0
        PEIMIN=SQRT(PMI2)
        PEJMIN=AIJ/PEIMIN
        PEKMIN=AIK/PEIMIN
        PAJMIN=SQRT(MAX(0D0,PEJMIN**2-PMJ2))
        PAKMIN=SQRT(MAX(0D0,PEKMIN**2-PMK2))
        FMIN=PEJMIN*PEKMIN+0.5D0*PAJMIN*PAKMIN-AJK
        PEIMAX=(AIJ+AIK)/SQRT(PMJ2+PMK2+2D0*AJK)
        IF(PMJ2.GT.1D-4) PEIMAX=AIJ/SQRT(PMJ2)
        PAIMAX=SQRT(MAX(0D0,PEIMAX**2-PMI2))
        HI=PEIMAX**2-0.25D0*PAIMAX**2
        PAJMAX=(PEIMAX*SQRT(MAX(0D0,AIJ**2-PMJ2*HI))-
     &  0.5D0*PAIMAX*AIJ)/HI
        PAKMAX=(PEIMAX*SQRT(MAX(0D0,AIK**2-PMK2*HI))-
     &  0.5D0*PAIMAX*AIK)/HI
        PEJMAX=SQRT(PAJMAX**2+PMJ2)
        PEKMAX=SQRT(PAKMAX**2+PMK2)
        FMAX=PEJMAX*PEKMAX+0.5D0*PAJMAX*PAKMAX-AJK
 
C...If unexpected values at upper endpoint then pick another parton.
        IF(FMAX.GT.0D0.AND.ITRY.LE.2) THEN
          I1=1+MOD(I,3)
          IF(A(I1,I1).GE.1D-4) THEN
            I=I1
            GOTO 130
          ENDIF
          ITRY=ITRY+1
          I1=1+MOD(I,3)
          IF(ITRY.LE.2.AND.A(I1,I1).GE.1D-4) THEN
            I=I1
            GOTO 130
          ENDIF
        ENDIF
 
C..Start binary + linear search to find solution inside range.
        ITER=0
        ITMIN=0
        ITMAX=0
        PAI=0.5D0*(PAIMIN+PAIMAX)
  140   ITER=ITER+1
 
C...Derive momentum of other two partons and distance to root.
        PEI=SQRT(PAI**2+PMI2)
        HI=PEI**2-0.25D0*PAI**2
        PAJ=(PEI*SQRT(MAX(0D0,AIJ**2-PMJ2*HI))-0.5D0*PAI*AIJ)/HI
        PEJ=SQRT(PAJ**2+PMJ2)
        PAK=(PEI*SQRT(MAX(0D0,AIK**2-PMK2*HI))-0.5D0*PAI*AIK)/HI
        PEK=SQRT(PAK**2+PMK2)
        FNOW=PEJ*PEK+0.5D0*PAJ*PAK-AJK
 
C...Pick next I momentum to explore, hopefully closer to root.
        IF(FNOW.GT.0D0) THEN
          PAIMIN=PAI
          FMIN=FNOW
          ITMIN=ITMIN+1
        ELSE
          PAIMAX=PAI
          FMAX=FNOW
          ITMAX=ITMAX+1
        ENDIF
        IF((ITER.LT.10.OR.ITMIN.LE.1.OR.ITMAX.LE.1).AND.ITER.LT.20)
     &  THEN
          PAI=0.5D0*(PAIMIN+PAIMAX)
          GOTO 140
        ELSEIF(ITER.LT.40.AND.FMIN.GT.0D0.AND.FMAX.LT.0D0.AND.
     &  ABS(FNOW).GT.1D-12*PSUM2) THEN
          PAI=PAIMIN+(PAIMAX-PAIMIN)*FMIN/(FMIN-FMAX)
          GOTO 140
        ENDIF
      ENDIF
 
C...Now know energies in junction rest frame.
      PENEW(I)=PEI
      PENEW(J)=PEJ
      PENEW(K)=PEK
 
C...Boost (copy of) partons to their rest frame.
      VXCM=-PSUM(1)/PSUM(5)
      VYCM=-PSUM(2)/PSUM(5)
      VZCM=-PSUM(3)/PSUM(5)
      GAMCM=SQRT(1D0+VXCM**2+VYCM**2+VZCM**2)
      DO 150 I=1,3
        FAC1=PJU(I,1)*VXCM+PJU(I,2)*VYCM+PJU(I,3)*VZCM
        FAC2=FAC1/(1D0+GAMCM)+PJU(I,4)
        PCM(I,1)=PJU(I,1)+FAC2*VXCM
        PCM(I,2)=PJU(I,2)+FAC2*VYCM
        PCM(I,3)=PJU(I,3)+FAC2*VZCM
        PCM(I,4)=PJU(I,4)*GAMCM+FAC1
        PCM(I,5)=SQRT(PCM(I,1)**2+PCM(I,2)**2+PCM(I,3)**2)
  150 CONTINUE
 
C...Construct difference vectors and boost to junction rest frame.
      DO 160 J=1,3
        PCM(4,J)=PCM(1,J)/PCM(1,4)-PCM(2,J)/PCM(2,4)
        PCM(5,J)=PCM(1,J)/PCM(1,4)-PCM(3,J)/PCM(3,4)
  160 CONTINUE
      PCM(4,4)=PENEW(1)/PCM(1,4)-PENEW(2)/PCM(2,4)
      PCM(5,4)=PENEW(1)/PCM(1,4)-PENEW(3)/PCM(3,4)
      PCM4S=PCM(4,1)**2+PCM(4,2)**2+PCM(4,3)**2
      PCM5S=PCM(5,1)**2+PCM(5,2)**2+PCM(5,3)**2
      PCM45=PCM(4,1)*PCM(5,1)+PCM(4,2)*PCM(5,2)+PCM(4,3)*PCM(5,3)
      C4=(PCM5S*PCM(4,4)-PCM45*PCM(5,4))/(PCM4S*PCM5S-PCM45**2)
      C5=(PCM4S*PCM(5,4)-PCM45*PCM(4,4))/(PCM4S*PCM5S-PCM45**2)
      VXJU=C4*PCM(4,1)+C5*PCM(5,1)
      VYJU=C4*PCM(4,2)+C5*PCM(5,2)
      VZJU=C4*PCM(4,3)+C5*PCM(5,3)
      GAMJU=SQRT(1D0+VXJU**2+VYJU**2+VZJU**2)
 
C...Add two boosts, giving final result.
      FCM=(VXJU*VXCM+VYJU*VYCM+VZJU*VZCM)/(1+GAMCM)+GAMJU
      VJU(1)=VXJU+FCM*VXCM
      VJU(2)=VYJU+FCM*VYCM
      VJU(3)=VZJU+FCM*VZCM
      VJU(4)=SQRT(1D0+VJU(1)**2+VJU(2)**2+VJU(3)**2)
      VJU(5)=1D0
 
C...In case of error in reconstruction: revert to CM frame of system.
      CTH12=(PCM(1,1)*PCM(2,1)+PCM(1,2)*PCM(2,2)+PCM(1,3)*PCM(2,3))/
     &(PCM(1,5)*PCM(2,5))
      CTH13=(PCM(1,1)*PCM(3,1)+PCM(1,2)*PCM(3,2)+PCM(1,3)*PCM(3,3))/
     &(PCM(1,5)*PCM(3,5))
      CTH23=(PCM(2,1)*PCM(3,1)+PCM(2,2)*PCM(3,2)+PCM(2,3)*PCM(3,3))/
     &(PCM(2,5)*PCM(3,5))
      ERRCCM=(CTH12+0.5D0)**2+(CTH13+0.5D0)**2+(CTH23+0.5D0)**2
      ERRTCM=TWOPI-ACOS(CTH12)-ACOS(CTH13)-ACOS(CTH23)
      DO 170 I=1,3
        FAC1=PJU(I,1)*VJU(1)+PJU(I,2)*VJU(2)+PJU(I,3)*VJU(3)
        FAC2=FAC1/(1D0+VJU(4))+PJU(I,4)
        PCM(I,1)=PJU(I,1)+FAC2*VJU(1)
        PCM(I,2)=PJU(I,2)+FAC2*VJU(2)
        PCM(I,3)=PJU(I,3)+FAC2*VJU(3)
        PCM(I,4)=PJU(I,4)*VJU(4)+FAC1
        PCM(I,5)=SQRT(PCM(I,1)**2+PCM(I,2)**2+PCM(I,3)**2)
  170 CONTINUE
      CTH12=(PCM(1,1)*PCM(2,1)+PCM(1,2)*PCM(2,2)+PCM(1,3)*PCM(2,3))/
     &(PCM(1,5)*PCM(2,5))
      CTH13=(PCM(1,1)*PCM(3,1)+PCM(1,2)*PCM(3,2)+PCM(1,3)*PCM(3,3))/
     &(PCM(1,5)*PCM(3,5))
      CTH23=(PCM(2,1)*PCM(3,1)+PCM(2,2)*PCM(3,2)+PCM(2,3)*PCM(3,3))/
     &(PCM(2,5)*PCM(3,5))
      ERRCJU=(CTH12+0.5D0)**2+(CTH13+0.5D0)**2+(CTH23+0.5D0)**2
      ERRTJU=TWOPI-ACOS(CTH12)-ACOS(CTH13)-ACOS(CTH23)
      IF(ERRCJU+ERRTJU.GT.ERRCCM+ERRTCM) THEN
        VJU(1)=VXCM
        VJU(2)=VYCM
        VJU(3)=VZCM
        VJU(4)=GAMCM
      ENDIF
 
      RETURN
      END
