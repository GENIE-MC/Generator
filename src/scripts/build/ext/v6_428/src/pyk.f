 
C*********************************************************************
 
C...PYK
C...Provides various integer-valued event related data.
 
      FUNCTION PYK(I,J)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
 
C...Default value. For I=0 number of entries, number of stable entries
C...or 3 times total charge.
      PYK=0
      IF(I.LT.0.OR.I.GT.MSTU(4).OR.J.LE.0) THEN
      ELSEIF(I.EQ.0.AND.J.EQ.1) THEN
        PYK=N
      ELSEIF(I.EQ.0.AND.(J.EQ.2.OR.J.EQ.6)) THEN
        DO 100 I1=1,N
          IF(J.EQ.2.AND.K(I1,1).GE.1.AND.K(I1,1).LE.10) PYK=PYK+1
          IF(J.EQ.6.AND.K(I1,1).GE.1.AND.K(I1,1).LE.10) PYK=PYK+
     &    PYCHGE(K(I1,2))
  100   CONTINUE
      ELSEIF(I.EQ.0) THEN
 
C...For I > 0 direct readout of K matrix or charge.
      ELSEIF(J.LE.5) THEN
        PYK=K(I,J)
      ELSEIF(J.EQ.6) THEN
        PYK=PYCHGE(K(I,2))
 
C...Status (existing/fragmented/decayed), parton/hadron separation.
      ELSEIF(J.LE.8) THEN
        IF(K(I,1).GE.1.AND.K(I,1).LE.10) PYK=1
        IF(J.EQ.8) PYK=PYK*K(I,2)
      ELSEIF(J.LE.12) THEN
        KFA=IABS(K(I,2))
        KC=PYCOMP(KFA)
        KQ=0
        IF(KC.NE.0) KQ=KCHG(KC,2)
        IF(J.EQ.9.AND.KC.NE.0.AND.KQ.NE.0) PYK=K(I,2)
        IF(J.EQ.10.AND.KC.NE.0.AND.KQ.EQ.0) PYK=K(I,2)
        IF(J.EQ.11) PYK=KC
        IF(J.EQ.12) PYK=KQ*ISIGN(1,K(I,2))
 
C...Heaviest flavour in hadron/diquark.
      ELSEIF(J.EQ.13) THEN
        KFA=IABS(K(I,2))
        PYK=MOD(KFA/100,10)*(-1)**MOD(KFA/100,10)
        IF(KFA.LT.10) PYK=KFA
        IF(MOD(KFA/1000,10).NE.0) PYK=MOD(KFA/1000,10)
        PYK=PYK*ISIGN(1,K(I,2))
 
C...Particle history: generation, ancestor, rank.
      ELSEIF(J.LE.15) THEN
        I2=I
        I1=I
  110   PYK=PYK+1
        I2=I1
        I1=K(I1,3)
        IF(I1.GT.0) THEN
          IF(K(I1,1).GT.0.AND.K(I1,1).LE.20) GOTO 110
        ENDIF
        IF(J.EQ.15) PYK=I2
      ELSEIF(J.EQ.16) THEN
        KFA=IABS(K(I,2))
        IF(K(I,1).LE.20.AND.((KFA.GE.11.AND.KFA.LE.20).OR.KFA.EQ.22.OR.
     &  (KFA.GT.100.AND.MOD(KFA/10,10).NE.0))) THEN
          I1=I
  120     I2=I1
          I1=K(I1,3)
          IF(I1.GT.0) THEN
            KFAM=IABS(K(I1,2))
            ILP=1
            IF(KFAM.NE.0.AND.KFAM.LE.10) ILP=0
            IF(KFAM.EQ.21.OR.KFAM.EQ.91.OR.KFAM.EQ.92.OR.KFAM.EQ.93)
     &      ILP=0
            IF(KFAM.GT.100.AND.MOD(KFAM/10,10).EQ.0) ILP=0
            IF(ILP.EQ.1) GOTO 120
          ENDIF
          IF(K(I1,1).EQ.12) THEN
            DO 130 I3=I1+1,I2
              IF(K(I3,3).EQ.K(I2,3).AND.K(I3,2).NE.91.AND.K(I3,2).NE.92
     &        .AND.K(I3,2).NE.93) PYK=PYK+1
  130       CONTINUE
          ELSE
            I3=I2
  140       PYK=PYK+1
            I3=I3+1
            IF(I3.LT.N.AND.K(I3,3).EQ.K(I2,3)) GOTO 140
          ENDIF
        ENDIF
 
C...Particle coming from collapsing jet system or not.
      ELSEIF(J.EQ.17) THEN
        I1=I
  150   PYK=PYK+1
        I3=I1
        I1=K(I1,3)
        I0=MAX(1,I1)
        KC=PYCOMP(K(I0,2))
        IF(I1.EQ.0.OR.K(I0,1).LE.0.OR.K(I0,1).GT.20.OR.KC.EQ.0) THEN
          IF(PYK.EQ.1) PYK=-1
          IF(PYK.GT.1) PYK=0
          RETURN
        ENDIF
        IF(KCHG(KC,2).EQ.0) GOTO 150
        IF(K(I1,1).NE.12) PYK=0
        IF(K(I1,1).NE.12) RETURN
        I2=I1
  160   I2=I2+1
        IF(I2.LT.N.AND.K(I2,1).NE.11) GOTO 160
        K3M=K(I3-1,3)
        IF(K3M.GE.I1.AND.K3M.LE.I2) PYK=0
        K3P=K(I3+1,3)
        IF(I3.LT.N.AND.K3P.GE.I1.AND.K3P.LE.I2) PYK=0
 
C...Number of decay products. Colour flow.
      ELSEIF(J.EQ.18) THEN
        IF(K(I,1).EQ.11.OR.K(I,1).EQ.12) PYK=MAX(0,K(I,5)-K(I,4)+1)
        IF(K(I,4).EQ.0.OR.K(I,5).EQ.0) PYK=0
      ELSEIF(J.LE.22) THEN
        IF(K(I,1).NE.3.AND.K(I,1).NE.13.AND.K(I,1).NE.14) RETURN
        IF(J.EQ.19) PYK=MOD(K(I,4)/MSTU(5),MSTU(5))
        IF(J.EQ.20) PYK=MOD(K(I,5)/MSTU(5),MSTU(5))
        IF(J.EQ.21) PYK=MOD(K(I,4),MSTU(5))
        IF(J.EQ.22) PYK=MOD(K(I,5),MSTU(5))
      ELSE
      ENDIF
 
      RETURN
      END
