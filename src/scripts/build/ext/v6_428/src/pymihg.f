 
*********************************************************************
 
C...PYMIHG
C...Collapse JCP1 and connecting tags to JCG1.
C...Collapse JCP2 and connecting tags to JCG2.
 
      SUBROUTINE PYMIHG(JCP1,JCG1,JCP2,JCG2)
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...The event record
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Parameters
      COMMON/PYINT1/MINT(400),VINT(400)
      SAVE /PYJETS/,/PYINT1/
C...Local variables
      COMMON /PYCBLS/MCO(4000,2),NCC,JCCO(4000,2),JCCN(4000,2),MACCPT
      COMMON /PYCTAG/NCT,MCT(4000,2)
      SAVE /PYCBLS/,/PYCTAG/
 
C...Break up JCP1<->JCP2 tag and create JCP1<->JCG1 and JCP2<->JCG2 tags
C...in temporary tag collapse array JCCN. Only break up one connection.
      MACCPT=1
      MCLPS=0
      DO 100 ICC=1,NCC
        JCCN(ICC,1)=JCCO(ICC,1)
        JCCN(ICC,2)=JCCO(ICC,2)
C...If there was a mother, it was previously connected to JCP1.
C...Should be changed to JCP2.
        IF (MCLPS.EQ.0) THEN
          IF (JCCN(ICC,1).EQ.MAX(JCP1,JCP2).AND.JCCN(ICC,2).EQ.MIN(JCP1
     &         ,JCP2)) THEN
            JCCN(ICC,1)=MAX(JCG2,JCP2)
            JCCN(ICC,2)=MIN(JCG2,JCP2)
            MCLPS=1
          ENDIF
        ENDIF
  100 CONTINUE
C...Also collapse colours on JCP1 side of JCG1
      IF (JCP1.NE.0) THEN
        JCCN(NCC+1,1)=MAX(JCP1,JCG1)
        JCCN(NCC+1,2)=MIN(JCP1,JCG1)
      ELSE
        JCCN(NCC+1,1)=MAX(JCP2,JCG2)
        JCCN(NCC+1,2)=MIN(JCP2,JCG2)
      ENDIF
 
C...Initialize event record colour tag array MCT array to MCO.
       DO 110 I=MINT(84)+1,N
        MCT(I,1)=MCO(I,1)
        MCT(I,2)=MCO(I,2)
  110 CONTINUE
 
C...Collapse tags:
C...IS = 1 : All tags connecting to JCG1 on JCG1 side -> JCG1
C...IS = 2 : All tags connecting to JCG2 on JCG2 side -> JCG2
C...IS = 3 : All tags connecting to JCG1 on JCP1 side -> JCG1
C...IS = 4 : All tags connecting to JCG2 on JCP2 side -> JCG2
      DO 160 IS=1,4
C...Skip if junction.
        IF ((IS.EQ.4.AND.JCP2.EQ.0).OR.(IS.EQ.3).AND.JCP1.EQ.0) GOTO 160
C...Define starting point in tag space.
C...JCA = previous tag
C...JCO = present tag
C...JCN = new tag
        IF (MOD(IS,2).EQ.1) THEN
          JCO=JCP1
          JCN=JCG1
          JCALL=JCG1
        ELSEIF (MOD(IS,2).EQ.0) THEN
          JCO=JCP2
          JCN=JCG2
          JCALL=JCG2
        ENDIF
        ITRACE=0
  120   ITRACE=ITRACE+1
        IF (ITRACE.GT.1000) THEN
C...NB: Proper error message should be defined here.
          CALL PYERRM(14
     &         ,'(PYMIHG:) Inf loop when collapsing colours.')
          MINT(57)=MINT(57)+1
          MINT(51)=1
          RETURN
        ENDIF
C...Collapse all JCN tags to JCALL
        DO 130 I=MINT(84)+1,N
          IF (MCO(I,1).EQ.JCN) MCT(I,1)=JCALL
          IF (MCO(I,2).EQ.JCN) MCT(I,2)=JCALL
  130   CONTINUE
C...IS = 1,2: first step forward. IS = 3,4: first step backward.
        IF (IS.GT.2.AND.(JCN.EQ.JCALL)) THEN
          JCA=JCN
          JCN=JCO
        ELSE
          JCA=JCO
          JCO=JCN
        ENDIF
C...If possible, step from JCO to new tag JCN not equal to JCA.
        DO 140 ICC=1,NCC+1
          IF (JCCN(ICC,1).EQ.JCO.AND.JCCN(ICC,2).NE.JCA) JCN=
     &         JCCN(ICC,2)
          IF (JCCN(ICC,2).EQ.JCO.AND.JCCN(ICC,1).NE.JCA) JCN=
     &         JCCN(ICC,1)
  140   CONTINUE
C...Iterate if new colour was arrived at, but don't go in circles.
        IF (JCN.NE.JCO.AND.JCN.NE.JCALL) GOTO 120
C...Change all JCN tags in MCO to JCALL in MCT.
        DO 150 I=MINT(84)+1,N
          IF (MCO(I,1).EQ.JCN) MCT(I,1)=JCALL
          IF (MCO(I,2).EQ.JCN) MCT(I,2)=JCALL
C...If gluon and colour tag = anticolour tag (and not = 0) try again.
          IF (K(I,2).EQ.21.AND.MCT(I,1).EQ.MCT(I,2).AND.MCT(I,1)
     &         .NE.0) MACCPT=0
  150   CONTINUE
  160 CONTINUE
 
      DO 200 JCL=NCT,1,-1
        JCA=0
        JCN=JCL
  170   JCO=JCN
        DO 180 ICC=1,NCC+1
          IF (JCCN(ICC,1).EQ.JCO.AND.JCCN(ICC,2).NE.JCA) JCN
     &         =JCCN(ICC,2)
          IF (JCCN(ICC,2).EQ.JCO.AND.JCCN(ICC,1).NE.JCA) JCN
     &         =JCCN(ICC,1)
  180   CONTINUE
C...Overpaint all JCN with JCL
        IF (JCN.NE.JCO.AND.JCN.NE.JCL) THEN
          DO 190 I=MINT(84)+1,N
            IF (MCT(I,1).EQ.JCN) MCT(I,1)=JCL
            IF (MCT(I,2).EQ.JCN) MCT(I,2)=JCL
C...If gluon and colour tag = anticolour tag (and not = 0) try again.
            IF (K(I,2).EQ.21.AND.MCT(I,1).EQ.MCT(I,2).AND.MCT(I,1)
     &           .NE.0) MACCPT=0
  190     CONTINUE
          JCA=JCO
          GOTO 170
        ENDIF
  200 CONTINUE
 
      RETURN
      END
