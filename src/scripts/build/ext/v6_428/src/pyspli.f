 
C*********************************************************************
 
C...PYSPLI
C...Splits a hadron remnant into two (partons or hadron + parton)
C...in case it is more complicated than just a quark or a diquark.
 
      SUBROUTINE PYSPLI(KF,KFLIN,KFLCH,KFLSP)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks. PYDAT1 temporary
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYINT1/MINT(400),VINT(400)
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      SAVE /PYPARS/,/PYINT1/,/PYDAT1/
C...Local array.
      DIMENSION KFL(3)
 
C...Preliminaries. Parton composition.
      KFA=IABS(KF)
      KFS=ISIGN(1,KF)
      KFL(1)=MOD(KFA/1000,10)
      KFL(2)=MOD(KFA/100,10)
      KFL(3)=MOD(KFA/10,10)
      IF(KFA.EQ.22.AND.MINT(109).EQ.2) THEN
        KFL(2)=INT(1.5D0+PYR(0))
        IF(MINT(105).EQ.333) KFL(2)=3
        IF(MINT(105).EQ.443) KFL(2)=4
        KFL(3)=KFL(2)
      ELSEIF((KFA.EQ.111.OR.KFA.EQ.113).AND.PYR(0).GT.0.5D0) THEN
        KFL(2)=2
        KFL(3)=2
      ELSEIF(KFA.EQ.223.AND.PYR(0).GT.0.5D0) THEN
        KFL(2)=1
        KFL(3)=1
      ELSEIF((KFA.EQ.130.OR.KFA.EQ.310).AND.PYR(0).GT.0.5D0) THEN
        KFL(2)=MOD(KFA/10,10)
        KFL(3)=MOD(KFA/100,10)
      ENDIF
      IF(KFLIN.NE.21.AND.KFLIN.NE.22.AND.KFLIN.NE.23) THEN
        KFLR=KFLIN*KFS
      ELSE
        KFLR=KFLIN
      ENDIF
      KFLCH=0
 
C...Subdivide lepton.
      IF(KFA.GE.11.AND.KFA.LE.18) THEN
        IF(KFLR.EQ.KFA) THEN
          KFLSP=KFS*22
        ELSEIF(KFLR.EQ.22) THEN
          KFLSP=KFA
        ELSEIF(KFLR.EQ.-24.AND.MOD(KFA,2).EQ.1) THEN
          KFLSP=KFA+1
        ELSEIF(KFLR.EQ.24.AND.MOD(KFA,2).EQ.0) THEN
          KFLSP=KFA-1
        ELSEIF(KFLR.EQ.21) THEN
          KFLSP=KFA
          KFLCH=KFS*21
        ELSE
          KFLSP=KFA
          KFLCH=-KFLR
        ENDIF
 
C...Subdivide photon.
      ELSEIF(KFA.EQ.22.AND.MINT(109).NE.2) THEN
        IF(KFLR.NE.21) THEN
          KFLSP=-KFLR
        ELSE
          RAGR=0.75D0*PYR(0)
          KFLSP=1
          IF(RAGR.GT.0.125D0) KFLSP=2
          IF(RAGR.GT.0.625D0) KFLSP=3
          IF(PYR(0).GT.0.5D0) KFLSP=-KFLSP
          KFLCH=-KFLSP
        ENDIF
 
C...Subdivide Reggeon or Pomeron.
      ELSEIF(KFA.EQ.110.OR.KFA.EQ.990) THEN
        IF(KFLIN.EQ.21) THEN
          KFLSP=KFS*21
        ELSE
          KFLSP=-KFLIN
        ENDIF
 
C...Subdivide meson.
      ELSEIF(KFL(1).EQ.0) THEN
        KFL(2)=KFL(2)*(-1)**KFL(2)
        KFL(3)=-KFL(3)*(-1)**IABS(KFL(2))
        IF(KFLR.EQ.KFL(2)) THEN
          KFLSP=KFL(3)
        ELSEIF(KFLR.EQ.KFL(3)) THEN
          KFLSP=KFL(2)
        ELSEIF(KFLR.EQ.21.AND.PYR(0).GT.0.5D0) THEN
          KFLSP=KFL(2)
          KFLCH=KFL(3)
        ELSEIF(KFLR.EQ.21) THEN
          KFLSP=KFL(3)
          KFLCH=KFL(2)
        ELSEIF(KFLR*KFL(2).GT.0) THEN
          NTRY=0
  100     NTRY=NTRY+1
          CALL PYKFDI(-KFLR,KFL(2),KFDUMP,KFLCH)
          IF(KFLCH.EQ.0.AND.NTRY.LT.100) THEN
            GOTO 100
          ELSEIF(KFLCH.EQ.0) THEN
            CALL PYERRM(14,'(PYSPLI:) caught in infinite loop')
            MINT(51)=1
            RETURN
          ENDIF
          KFLSP=KFL(3)
        ELSE
          NTRY=0
  110     NTRY=NTRY+1
          CALL PYKFDI(-KFLR,KFL(3),KFDUMP,KFLCH)
          IF(KFLCH.EQ.0.AND.NTRY.LT.100) THEN
            GOTO 110
          ELSEIF(KFLCH.EQ.0) THEN
            CALL PYERRM(14,'(PYSPLI:) caught in infinite loop')
            MINT(51)=1
            RETURN
          ENDIF
          KFLSP=KFL(2)
        ENDIF

C...Special case for extracting photon from baryon without splitting
C...the latter. (Currently only used by external programs.)
      ELSEIF(KFLIN.EQ.22.AND.MSTP(98).EQ.1) then
        KFLSP=KFA
        KFLCH=0
 
C...Subdivide baryon.
      ELSE
        NAGR=0
        DO 120 J=1,3
          IF(KFLR.EQ.KFL(J)) NAGR=NAGR+1
  120   CONTINUE
        IF(NAGR.GE.1) THEN
          RAGR=0.00001D0+(NAGR-0.00002D0)*PYR(0)
          IAGR=0
          DO 130 J=1,3
            IF(KFLR.EQ.KFL(J)) RAGR=RAGR-1D0
            IF(IAGR.EQ.0.AND.RAGR.LE.0D0) IAGR=J
  130     CONTINUE
        ELSE
          IAGR=1.00001D0+2.99998D0*PYR(0)
        ENDIF
        ID1=1
        IF(IAGR.EQ.1) ID1=2
        IF(IAGR.EQ.1.AND.KFL(3).GT.KFL(2)) ID1=3
        ID2=6-IAGR-ID1
        KSP=3
        IF(MOD(KFA,10).EQ.2.AND.KFL(1).EQ.KFL(2)) THEN
          IF(IAGR.NE.3.AND.PYR(0).GT.0.25D0) KSP=1
        ELSEIF(MOD(KFA,10).EQ.2.AND.KFL(2).GE.KFL(3)) THEN
          IF(IAGR.NE.1.AND.PYR(0).GT.0.25D0) KSP=1
        ELSEIF(MOD(KFA,10).EQ.2) THEN
          IF(IAGR.EQ.1) KSP=1
          IF(IAGR.NE.1.AND.PYR(0).GT.0.75D0) KSP=1
        ENDIF
        KFLSP=1000*KFL(ID1)+100*KFL(ID2)+KSP
        IF(KFLR.EQ.21) THEN
          KFLCH=KFL(IAGR)
        ELSEIF(NAGR.EQ.0.AND.KFLR.GT.0) THEN
          NTRY=0
  140     NTRY=NTRY+1
          CALL PYKFDI(-KFLR,KFL(IAGR),KFDUMP,KFLCH)
          IF(KFLCH.EQ.0.AND.NTRY.LT.100) THEN
            GOTO 140
          ELSEIF(KFLCH.EQ.0) THEN
            CALL PYERRM(14,'(PYSPLI:) caught in infinite loop')
            MINT(51)=1
            RETURN
          ENDIF
        ELSEIF(NAGR.EQ.0) THEN
          NTRY=0
  150     NTRY=NTRY+1
          CALL PYKFDI(10000*KFL(ID1)+KFLSP,-KFLR,KFDUMP,KFLCH)
          IF(KFLCH.EQ.0.AND.NTRY.LT.100) THEN
            GOTO 150
          ELSEIF(KFLCH.EQ.0) THEN
            CALL PYERRM(14,'(PYSPLI:) caught in infinite loop')
            MINT(51)=1
            RETURN
          ENDIF
          KFLSP=KFL(IAGR)
        ENDIF
      ENDIF
 
C...Add on correct sign for result.
      KFLCH=KFLCH*KFS
      KFLSP=KFLSP*KFS
 
      RETURN
      END
