 
C*********************************************************************
 
C...PYADSH
C...Administers the generation of successive final-state showers
C...in external processes.
 
      SUBROUTINE PYADSH(NFIN)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Parameter statement for maximum size of showers.
      PARAMETER (MAXNUR=1000)
C...Commonblocks.
      COMMON/PYPART/NPART,NPARTD,IPART(MAXNUR),PTPART(MAXNUR)
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYCTAG/NCT,MCT(4000,2)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYPART/,/PYJETS/,/PYCTAG/,/PYDAT1/,/PYPARS/,/PYINT1/
C...Local array.
      DIMENSION IBEG(100),KSAV(100,5),PSUM(4),BETA(3)
 
C...Set primary vertex.
      DO 100 J=1,5
        V(MINT(83)+5,J)=0D0
        V(MINT(83)+6,J)=0D0
        V(MINT(84)+1,J)=0D0
        V(MINT(84)+2,J)=0D0
  100 CONTINUE
 
C...Isolate systems of particles with the same mother.
      NSYS=0
      IMS=-1
      DO 140 I=MINT(84)+3,NFIN
        IM=K(I,3)
        IF(IM.GT.0.AND.IM.LE.MINT(84)) IM=K(IM,3)
        IF(IM.NE.IMS) THEN
          NSYS=NSYS+1
          IBEG(NSYS)=I
          IMS=IM
        ENDIF
 
C...Set production vertices.
        IF(IM.LE.MINT(83)+6.OR.(IM.GT.MINT(84).AND.IM.LE.MINT(84)+2))
     &  THEN
          DO 110 J=1,4
            V(I,J)=0D0
  110     CONTINUE
        ELSE
          DO 120 J=1,4
            V(I,J)=V(IM,J)+V(IM,5)*P(IM,J)/P(IM,5)
  120     CONTINUE
        ENDIF
        IF(MSTP(125).GE.1) THEN
          IDOC=I-MSTP(126)+4
          DO 130 J=1,5
            V(IDOC,J)=V(I,J)
  130     CONTINUE
        ENDIF
  140 CONTINUE
 
C...End loop over systems. Return if no showers to be performed.
      IBEG(NSYS+1)=NFIN+1
      IF(MSTP(71).LE.0) RETURN
 
C...Loop through systems of particles; check that sensible size.
      DO 270 ISYS=1,NSYS
        NSIZ=IBEG(ISYS+1)-IBEG(ISYS)
        IF(MINT(35).LE.2) THEN
          IF(NSIZ.EQ.1.AND.ISYS.EQ.1) THEN
            GOTO 270
          ELSEIF(NSIZ.LE.1) THEN
            CALL PYERRM(2,'(PYADSH:) only one particle in system')
            GOTO 270
          ELSEIF(NSIZ.GT.80) THEN
            CALL PYERRM(2,'(PYADSH:) more than 80 particles in system')
            GOTO 270
          ENDIF
        ENDIF
 
C...Save status codes and daughters of showering particles; reset them.
        DO 150 J=1,4
          PSUM(J)=0D0
  150   CONTINUE
        DO 170 II=1,NSIZ
          I=IBEG(ISYS)-1+II
          KSAV(II,1)=K(I,1)
          IF(K(I,1).GT.10) THEN
            K(I,1)=1
            IF(KSAV(II,1).EQ.14) K(I,1)=3
          ENDIF
          IF(KSAV(II,1).LE.10) THEN
          ELSEIF(K(I,1).EQ.1) THEN
            KSAV(II,4)=K(I,4)
            KSAV(II,5)=K(I,5)
            K(I,4)=0
            K(I,5)=0
          ELSE
            KSAV(II,4)=MOD(K(I,4),MSTU(5))
            KSAV(II,5)=MOD(K(I,5),MSTU(5))
            K(I,4)=K(I,4)-KSAV(II,4)
            K(I,5)=K(I,5)-KSAV(II,5)
          ENDIF
          DO 160 J=1,4
            PSUM(J)=PSUM(J)+P(I,J)
  160     CONTINUE
  170   CONTINUE
 
C...Perform shower.
        QMAX=SQRT(MAX(0D0,PSUM(4)**2-PSUM(1)**2-PSUM(2)**2-
     &  PSUM(3)**2))
        IF(ISYS.EQ.1) QMAX=MIN(QMAX,SQRT(PARP(71))*VINT(55))
        NSAV=N
        IF(MINT(35).LE.2) THEN
          IF(NSIZ.EQ.2) THEN
            CALL PYSHOW(IBEG(ISYS),IBEG(ISYS)+1,QMAX)
          ELSE
            CALL PYSHOW(IBEG(ISYS),-NSIZ,QMAX)
          ENDIF
 
C...For external processes, first call, also ISR partons radiate.
C...Can use existing PYPART list, removing partons that radiate later.
        ELSEIF(ISYS.EQ.1) THEN
          NPARTN=0
          DO 175 II=1,NPART
            IF(IPART(II).LT.IBEG(2).OR.IPART(II).GE.IBEG(NSYS+1)) THEN
              NPARTN=NPARTN+1
              IPART(NPARTN)=IPART(II)
              PTPART(NPARTN)=PTPART(II)
            ENDIF
 175      CONTINUE
          NPART=NPARTN
          CALL PYPTFS(1,0.5D0*QMAX,0D0,PTGEN)
        ELSE
C...For subsequent calls use the systems excluded above.
          NPART=NSIZ
          NPARTD=0
          DO 180 II=1,NSIZ
            I=IBEG(ISYS)-1+II
            IPART(II)=I
            PTPART(II)=0.5D0*QMAX
  180     CONTINUE
          CALL PYPTFS(2,0.5D0*QMAX,0D0,PTGEN)
        ENDIF
 
C...Look up showered copies of original showering particles.
        DO 260 II=1,NSIZ
          I=IBEG(ISYS)-1+II
          IMV=I
C...Particles without daughters need not be studied.
          IF(KSAV(II,1).LE.10) GOTO 260
          IF(N.EQ.NSAV.OR.K(I,1).LE.10) THEN
          ELSEIF(K(I,1).EQ.11) THEN
  190       IMV=MOD(K(IMV,4),MSTU(5))
            IF(K(IMV,1).EQ.11) GOTO 190
          ELSE
            KDA1=MOD(K(I,4),MSTU(5))
            IF(KDA1.GT.0) THEN
              IF(K(KDA1,2).EQ.21) KDA1=K(KDA1,5)/MSTU(5)
            ENDIF
            KDA2=MOD(K(I,5),MSTU(5))
            IF(KDA2.GT.0) THEN
              IF(K(KDA2,2).EQ.21) KDA2=K(KDA2,4)/MSTU(5)
            ENDIF
            DO 200 I3=I+1,N
              IF(K(I3,2).EQ.K(I,2).AND.(I3.EQ.KDA1.OR.I3.EQ.KDA2))
     &        THEN
                IMV=I3
                KDA1=MOD(K(I3,4),MSTU(5))
                IF(KDA1.GT.0) THEN
                  IF(K(KDA1,2).EQ.21) KDA1=K(KDA1,5)/MSTU(5)
                ENDIF
                KDA2=MOD(K(I3,5),MSTU(5))
                IF(KDA2.GT.0) THEN
                  IF(K(KDA2,2).EQ.21) KDA2=K(KDA2,4)/MSTU(5)
                ENDIF
              ENDIF
  200       CONTINUE
          ENDIF
 
C...Restore daughter info of original partons to showered copies.
          IF(KSAV(II,1).GT.10) K(IMV,1)=KSAV(II,1)
          IF(KSAV(II,1).LE.10) THEN
          ELSEIF(K(I,1).EQ.1) THEN
            K(IMV,4)=KSAV(II,4)
            K(IMV,5)=KSAV(II,5)
          ELSE
            K(IMV,4)=K(IMV,4)+KSAV(II,4)
            K(IMV,5)=K(IMV,5)+KSAV(II,5)
          ENDIF
 
C...Reset mother info of existing daughters to showered copies.
          DO 210 I3=IBEG(ISYS+1),NFIN
            IF(K(I3,3).EQ.I) K(I3,3)=IMV
            IF(K(I3,1).EQ.3.OR.K(I3,1).EQ.14) THEN
              IF(K(I3,4)/MSTU(5).EQ.I) K(I3,4)=K(I3,4)+MSTU(5)*(IMV-I)
              IF(K(I3,5)/MSTU(5).EQ.I) K(I3,5)=K(I3,5)+MSTU(5)*(IMV-I)
            ENDIF
  210     CONTINUE
 
C...Boost all original daughters to new frame of showered copy.
C...Also update their colour tags.
          IF(IMV.NE.I) THEN
            DO 220 J=1,3
              BETA(J)=(P(IMV,J)-P(I,J))/(P(IMV,4)+P(I,4))
  220       CONTINUE
            FAC=2D0/(1D0+BETA(1)**2+BETA(2)**2+BETA(3)**2)
            DO 230 J=1,3
              BETA(J)=FAC*BETA(J)
  230       CONTINUE
            DO 250 I3=IBEG(ISYS+1),NFIN
              IMO=I3
  240         IMO=K(IMO,3)
              IF(MSTP(128).LE.0) THEN
                IF(IMO.GT.0.AND.IMO.NE.I.AND.IMO.NE.K(I,3)) GOTO 240
                IF(IMO.EQ.I.OR.(K(I,3).LE.MINT(84).AND.IMO.EQ.K(I,3)))
     &          THEN
                  CALL PYROBO(I3,I3,0D0,0D0,BETA(1),BETA(2),BETA(3))
                  IF(MCT(I3,1).EQ.MCT(I,1)) MCT(I3,1)=MCT(IMV,1)
                  IF(MCT(I3,2).EQ.MCT(I,2)) MCT(I3,2)=MCT(IMV,2)
                ENDIF
              ELSE
                IF(IMO.EQ.IMV) THEN
                  CALL PYROBO(I3,I3,0D0,0D0,BETA(1),BETA(2),BETA(3))
                  IF(MCT(I3,1).EQ.MCT(I,1)) MCT(I3,1)=MCT(IMV,1)
                  IF(MCT(I3,2).EQ.MCT(I,2)) MCT(I3,2)=MCT(IMV,2)
                ELSEIF(IMO.GT.0.AND.IMO.NE.I.AND.IMO.NE.K(I,3)) THEN
                  GOTO 240
                ENDIF
              ENDIF
  250       CONTINUE
          ENDIF
  260   CONTINUE
 
C...End of loop over showering systems
  270 CONTINUE
 
      RETURN
      END
