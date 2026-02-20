 
C*********************************************************************
 
C...PYEDIT
C...Performs global manipulations on the event record, in particular
C...to exclude unstable or undetectable partons/particles.
 
      SUBROUTINE PYEDIT(MEDIT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYCTAG/NCT,MCT(4000,2)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/,/PYCTAG/
C...Local arrays.
      DIMENSION NS(2),PTS(2),PLS(2)
 
C...Remove unwanted partons/particles.
      IF((MEDIT.GE.0.AND.MEDIT.LE.3).OR.MEDIT.EQ.5) THEN
        IMAX=N
        IF(MSTU(2).GT.0) IMAX=MSTU(2)
        I1=MAX(1,MSTU(1))-1
        DO 110 I=MAX(1,MSTU(1)),IMAX
          IF(K(I,1).EQ.0.OR.(K(I,1).GE.21.AND.K(I,1).LE.40)) GOTO 110
          IF(MEDIT.EQ.1) THEN
            IF(K(I,1).GT.10.AND.K(I,1).NE.41.AND.K(I,1).NE.42) GOTO 110
          ELSEIF(MEDIT.EQ.2) THEN
            IF(K(I,1).GT.10.AND.K(I,1).NE.41.AND.K(I,1).NE.42) GOTO 110
            KC=PYCOMP(K(I,2))
            IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.
     &      KC.EQ.18.OR.K(I,2).EQ.KSUSY1+22.OR.K(I,2).EQ.39.OR.
     &      K(I,2).EQ.KSUSY1+39) GOTO 110
          ELSEIF(MEDIT.EQ.3) THEN
            IF(K(I,1).GT.10.AND.K(I,1).NE.41.AND.K(I,1).NE.42) GOTO 110
            KC=PYCOMP(K(I,2))
            IF(KC.EQ.0) GOTO 110
            IF(KCHG(KC,2).EQ.0.AND.PYCHGE(K(I,2)).EQ.0) GOTO 110
          ELSEIF(MEDIT.EQ.5) THEN
            IF(K(I,1).EQ.13.OR.K(I,1).EQ.14.OR.K(I,1).EQ.52) GOTO 110
            KC=PYCOMP(K(I,2))
            IF(KC.EQ.0) GOTO 110
            IF(K(I,1).GT.10.AND.K(I,1).NE.41.AND.K(I,1).NE.42.AND.
     &      KCHG(KC,2).EQ.0) GOTO 110
          ENDIF
 
C...Pack remaining partons/particles. Origin no longer known.
          I1=I1+1
          DO 100 J=1,5
            K(I1,J)=K(I,J)
            P(I1,J)=P(I,J)
            V(I1,J)=V(I,J)
  100     CONTINUE
          K(I1,3)=0
  110   CONTINUE
        IF(I1.LT.N) MSTU(3)=0
        IF(I1.LT.N) MSTU(70)=0
        N=I1
 
C...Selective removal of class of entries. New position of retained.
      ELSEIF(MEDIT.GE.11.AND.MEDIT.LE.15) THEN
        I1=0
        DO 120 I=1,N
          K(I,3)=MOD(K(I,3),MSTU(5))
          IF(MEDIT.EQ.11.AND.K(I,1).LT.0) GOTO 120
          IF(MEDIT.EQ.12.AND.K(I,1).EQ.0) GOTO 120
          IF(MEDIT.EQ.13.AND.(K(I,1).EQ.11.OR.K(I,1).EQ.12.OR.
     &    K(I,1).EQ.15.OR.K(I,1).EQ.51).AND.K(I,2).NE.94) GOTO 120
          IF(MEDIT.EQ.14.AND.(K(I,1).EQ.13.OR.K(I,1).EQ.14.OR.
     &    K(I,1).EQ.52.OR.K(I,2).EQ.94)) GOTO 120
          IF(MEDIT.EQ.15.AND.K(I,1).GE.21.AND.K(I,1).LE.40) GOTO 120
          I1=I1+1
          K(I,3)=K(I,3)+MSTU(5)*I1
  120   CONTINUE
 
C...Find new event history information and replace old.
        DO 140 I=1,N
          IF(K(I,1).LE.0.OR.(K(I,1).GE.21.AND.K(I,1).LE.40).OR.
     &    K(I,3)/MSTU(5).EQ.0) GOTO 140
          ID=I
  130     IM=MOD(K(ID,3),MSTU(5))
          IF(MEDIT.EQ.13.AND.IM.GT.0.AND.IM.LE.N) THEN
            IF((K(IM,1).EQ.11.OR.K(IM,1).EQ.12.OR.K(IM,1).EQ.15.OR.
     &      K(IM,1).EQ.51).AND.K(IM,2).NE.94) THEN
              ID=IM
              GOTO 130
            ENDIF
          ELSEIF(MEDIT.EQ.14.AND.IM.GT.0.AND.IM.LE.N) THEN
            IF(K(IM,1).EQ.13.OR.K(IM,1).EQ.14.OR.K(IM,1).EQ.52.OR.
     &      K(IM,2).EQ.94) THEN
              ID=IM
              GOTO 130
            ENDIF
          ENDIF
          K(I,3)=MSTU(5)*(K(I,3)/MSTU(5))
          IF(IM.NE.0) K(I,3)=K(I,3)+K(IM,3)/MSTU(5)
          IF(K(I,1).NE.3.AND.K(I,1).NE.13.AND.K(I,1).NE.14.AND.
     &      K(I,1).NE.42.AND.K(I,1).NE.52) THEN
            IF(K(I,4).GT.0.AND.K(I,4).LE.MSTU(4)) K(I,4)=
     &      K(K(I,4),3)/MSTU(5)
            IF(K(I,5).GT.0.AND.K(I,5).LE.MSTU(4)) K(I,5)=
     &      K(K(I,5),3)/MSTU(5)
          ELSE
            KCM=MOD(K(I,4)/MSTU(5),MSTU(5))
            IF(KCM.GT.0.AND.KCM.LE.MSTU(4).AND.K(I,1).NE.42.AND.
     &      K(I,1).NE.52) KCM=K(KCM,3)/MSTU(5)
            KCD=MOD(K(I,4),MSTU(5))
            IF(KCD.GT.0.AND.KCD.LE.MSTU(4)) KCD=K(KCD,3)/MSTU(5)
            K(I,4)=MSTU(5)**2*(K(I,4)/MSTU(5)**2)+MSTU(5)*KCM+KCD
            KCM=MOD(K(I,5)/MSTU(5),MSTU(5))
            IF(KCM.GT.0.AND.KCM.LE.MSTU(4)) KCM=K(KCM,3)/MSTU(5)
            KCD=MOD(K(I,5),MSTU(5))
            IF(KCD.GT.0.AND.KCD.LE.MSTU(4)) KCD=K(KCD,3)/MSTU(5)
            K(I,5)=MSTU(5)**2*(K(I,5)/MSTU(5)**2)+MSTU(5)*KCM+KCD
          ENDIF
  140   CONTINUE
 
C...Pack remaining entries.
        I1=0
        MSTU90=MSTU(90)
        MSTU(90)=0
        DO 170 I=1,N
          IF(K(I,3)/MSTU(5).EQ.0) GOTO 170
          I1=I1+1
          DO 150 J=1,5
            K(I1,J)=K(I,J)
            P(I1,J)=P(I,J)
            V(I1,J)=V(I,J)
  150     CONTINUE
C...Also update LHA1 colour tags
          MCT(I1,1)=MCT(I,1)
          MCT(I1,2)=MCT(I,2)
          K(I1,3)=MOD(K(I1,3),MSTU(5))
          DO 160 IZ=1,MSTU90
            IF(I.EQ.MSTU(90+IZ)) THEN
              MSTU(90)=MSTU(90)+1
              MSTU(90+MSTU(90))=I1
              PARU(90+MSTU(90))=PARU(90+IZ)
            ENDIF
  160     CONTINUE
  170   CONTINUE
        IF(I1.LT.N) MSTU(3)=0
        IF(I1.LT.N) MSTU(70)=0
        N=I1
 
C...Fill in some missing daughter pointers (lost in colour flow).
      ELSEIF(MEDIT.EQ.16) THEN
        DO 220 I=1,N
          IF(K(I,1).LE.10.OR.(K(I,1).GE.21.AND.K(I,1).LE.50)) GOTO 220
          IF(K(I,4).NE.0.OR.K(I,5).NE.0) GOTO 220
C...Find daughters who point to mother.
          DO 180 I1=I+1,N
            IF(K(I1,3).NE.I) THEN
            ELSEIF(K(I,4).EQ.0) THEN
              K(I,4)=I1
            ELSE
              K(I,5)=I1
            ENDIF
  180     CONTINUE
          IF(K(I,5).EQ.0) K(I,5)=K(I,4)
          IF(K(I,4).NE.0) GOTO 220
C...Find daughters who point to documentation version of mother.
          IM=K(I,3)
          IF(IM.LE.0.OR.IM.GE.I) GOTO 220
          IF(K(IM,1).LE.20.OR.K(IM,1).GT.30) GOTO 220
          IF(K(IM,2).NE.K(I,2).OR.ABS(P(IM,5)-P(I,5)).GT.1D-2) GOTO 220
          DO 190 I1=I+1,N
            IF(K(I1,3).NE.IM) THEN
            ELSEIF(K(I,4).EQ.0) THEN
              K(I,4)=I1
            ELSE
              K(I,5)=I1
            ENDIF
  190     CONTINUE
          IF(K(I,5).EQ.0) K(I,5)=K(I,4)
          IF(K(I,4).NE.0) GOTO 220
C...Find daughters who point to documentation daughters who,
C...in their turn, point to documentation mother.
          ID1=IM
          ID2=IM
          DO 200 I1=IM+1,I-1
            IF(K(I1,3).EQ.IM.AND.K(I1,1).GE.21.AND.K(I1,1).LE.30) THEN
              ID2=I1
              IF(ID1.EQ.IM) ID1=I1
            ENDIF
  200     CONTINUE
          DO 210 I1=I+1,N
            IF(K(I1,3).NE.ID1.AND.K(I1,3).NE.ID2) THEN
            ELSEIF(K(I,4).EQ.0) THEN
              K(I,4)=I1
            ELSE
              K(I,5)=I1
            ENDIF
  210     CONTINUE
          IF(K(I,5).EQ.0) K(I,5)=K(I,4)
  220   CONTINUE
 
C...Save top entries at bottom of PYJETS commonblock.
      ELSEIF(MEDIT.EQ.21) THEN
        IF(2*N.GE.MSTU(4)) THEN
          CALL PYERRM(11,'(PYEDIT:) no more memory left in PYJETS')
          RETURN
        ENDIF
        DO 240 I=1,N
          DO 230 J=1,5
            K(MSTU(4)-I,J)=K(I,J)
            P(MSTU(4)-I,J)=P(I,J)
            V(MSTU(4)-I,J)=V(I,J)
  230     CONTINUE
  240   CONTINUE
        MSTU(32)=N
 
C...Restore bottom entries of commonblock PYJETS to top.
      ELSEIF(MEDIT.EQ.22) THEN
        DO 260 I=1,MSTU(32)
          DO 250 J=1,5
            K(I,J)=K(MSTU(4)-I,J)
            P(I,J)=P(MSTU(4)-I,J)
            V(I,J)=V(MSTU(4)-I,J)
  250     CONTINUE
  260   CONTINUE
        N=MSTU(32)
 
C...Mark primary entries at top of commonblock PYJETS as untreated.
      ELSEIF(MEDIT.EQ.23) THEN
        I1=0
        DO 270 I=1,N
          KH=K(I,3)
          IF(KH.GE.1) THEN
            IF(K(KH,1).GE.21.AND.K(KH,1).LE.30) KH=0
          ENDIF
          IF(KH.NE.0) GOTO 280
          I1=I1+1
          IF(K(I,1).GE.11.AND.K(I,1).LE.20) K(I,1)=K(I,1)-10
          IF(K(I,1).GE.51.AND.K(I,1).LE.60) K(I,1)=K(I,1)-10
  270   CONTINUE
  280   N=I1
 
C...Place largest axis along z axis and second largest in xy plane.
      ELSEIF(MEDIT.EQ.31.OR.MEDIT.EQ.32) THEN
        CALL PYROBO(1,N+MSTU(3),0D0,-PYANGL(P(MSTU(61),1),
     &  P(MSTU(61),2)),0D0,0D0,0D0)
        CALL PYROBO(1,N+MSTU(3),-PYANGL(P(MSTU(61),3),
     &  P(MSTU(61),1)),0D0,0D0,0D0,0D0)
        CALL PYROBO(1,N+MSTU(3),0D0,-PYANGL(P(MSTU(61)+1,1),
     &  P(MSTU(61)+1,2)),0D0,0D0,0D0)
        IF(MEDIT.EQ.31) RETURN
 
C...Rotate to put slim jet along +z axis.
        DO 290 IS=1,2
          NS(IS)=0
          PTS(IS)=0D0
          PLS(IS)=0D0
  290   CONTINUE
        DO 300 I=1,N
          IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 300
          IF(MSTU(41).GE.2) THEN
            KC=PYCOMP(K(I,2))
            IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.
     &      KC.EQ.18.OR.K(I,2).EQ.KSUSY1+22.OR.K(I,2).EQ.39.OR.
     &      K(I,2).EQ.KSUSY1+39) GOTO 300
            IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.PYCHGE(K(I,2))
     &      .EQ.0) GOTO 300
          ENDIF
          IS=2D0-SIGN(0.5D0,P(I,3))
          NS(IS)=NS(IS)+1
          PTS(IS)=PTS(IS)+SQRT(P(I,1)**2+P(I,2)**2)
  300   CONTINUE
        IF(NS(1)*PTS(2)**2.LT.NS(2)*PTS(1)**2)
     &  CALL PYROBO(1,N+MSTU(3),PARU(1),0D0,0D0,0D0,0D0)
 
C...Rotate to put second largest jet into -z,+x quadrant.
        DO 310 I=1,N
          IF(P(I,3).GE.0D0) GOTO 310
          IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 310
          IF(MSTU(41).GE.2) THEN
            KC=PYCOMP(K(I,2))
            IF(KC.EQ.0.OR.KC.EQ.12.OR.KC.EQ.14.OR.KC.EQ.16.OR.
     &      KC.EQ.18.OR.K(I,2).EQ.KSUSY1+22.OR.K(I,2).EQ.39.OR.
     &      K(I,2).EQ.KSUSY1+39) GOTO 310
            IF(MSTU(41).GE.3.AND.KCHG(KC,2).EQ.0.AND.PYCHGE(K(I,2))
     &      .EQ.0) GOTO 310
          ENDIF
          IS=2D0-SIGN(0.5D0,P(I,1))
          PLS(IS)=PLS(IS)-P(I,3)
  310   CONTINUE
        IF(PLS(2).GT.PLS(1)) CALL PYROBO(1,N+MSTU(3),0D0,PARU(1),
     &  0D0,0D0,0D0)
      ENDIF
 
      RETURN
      END
