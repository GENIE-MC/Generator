 
C*********************************************************************
 
C...PYRVSF
C...Calculates R-violating decays of sfermions.
C...P. Z. Skands
 
      SUBROUTINE PYRVSF(KFIN,XLAM,IDLAM,LKNT)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
C...Parameter statement to help give large particle numbers.
      PARAMETER (KSUSY1=1000000,KSUSY2=2000000,KTECHN=3000000,
     &KEXCIT=4000000,KDIMEN=5000000)
C...Commonblocks.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYMSSM/IMSS(0:99),RMSS(0:99)
      COMMON/PYSSMT/ZMIX(4,4),UMIX(2,2),VMIX(2,2),SMZ(4),SMW(2),
     &SFMIX(16,4),ZMIXI(4,4),UMIXI(2,2),VMIXI(2,2)
      COMMON/PYMSRV/RVLAM(3,3,3), RVLAMP(3,3,3), RVLAMB(3,3,3)
C...Local variables.
      DOUBLE PRECISION XLAM(0:400)
      INTEGER IDLAM(400,3), PYCOMP
      SAVE /PYMSRV/,/PYSSMT/,/PYMSSM/,/PYDAT2/
 
C...IS R-VIOLATION ON ?
      IF ((IMSS(51).GE.1).OR.(IMSS(52).GE.1).OR.(IMSS(53).GE.1)) THEN
C...Mass eigenstate counter
        ICNT=INT(KFIN/KSUSY1)
C...SM KF code of SUSY particle
        KFSM=KFIN-ICNT*KSUSY1
C...Squared Sparticle Mass
        SM=PMAS(PYCOMP(KFIN),1)**2
C... Squared mass of top quark
        SMT=PMAS(PYCOMP(6),1)**2
C...IS L-VIOLATION ON ?
        IF ((IMSS(51).GE.1).OR.(IMSS(52).GE.1)) THEN
C...SLEPTON -> NU(BAR) + LEPTON and UBAR + D
          IF(ICNT.NE.0.AND.(KFSM.EQ.11.OR.KFSM.EQ.13.OR.KFSM.EQ.15))
     &         THEN
            K=INT((KFSM-9)/2)
            DO 110 I=1,3
              DO 100 J=1,3
                IF(I.NE.J) THEN
C...~e,~mu,~tau -> nu_I + lepton-_J
                  LKNT = LKNT+1
                  IDLAM(LKNT,1)= 12 +2*(I-1)
                  IDLAM(LKNT,2)= 11 +2*(J-1)
                  IDLAM(LKNT,3)= 0
                  XLAM(LKNT)=0D0
                  RM2=RVLAM(I,J,K)**2*SFMIX(KFSM,2*ICNT)**2 * SM
                  IF (IMSS(51).NE.0) XLAM(LKNT) =
     &                 PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
C...KINEMATICS CHECK
                  IF (XLAM(LKNT).EQ.0D0) THEN
                    LKNT=LKNT-1
                  ENDIF
                ENDIF
  100         CONTINUE
  110       CONTINUE
C...~e,~mu,~tau -> nu_Ibar + lepton-_K
            J=INT((KFSM-9)/2)
            DO 130 I=1,3
              IF(I.NE.J) THEN
                DO 120 K=1,3
                  LKNT = LKNT+1
                  IDLAM(LKNT,1)=-12 -2*(I-1)
                  IDLAM(LKNT,2)= 11 +2*(K-1)
                  IDLAM(LKNT,3)= 0
                  XLAM(LKNT)=0D0
                  RM2=RVLAM(I,J,K)**2*SFMIX(KFSM,2*ICNT-1)**2 * SM
                  IF (IMSS(51).NE.0) XLAM(LKNT) =
     &                 PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
C...KINEMATICS CHECK
                  IF (XLAM(LKNT).EQ.0D0) THEN
                    LKNT=LKNT-1
                  ENDIF
  120           CONTINUE
              ENDIF
  130       CONTINUE
C...~e,~mu,~tau -> u_Jbar + d_K
            I=INT((KFSM-9)/2)
            DO 150 J=1,3
              DO 140 K=1,3
                LKNT = LKNT+1
                IDLAM(LKNT,1)=-2 -2*(J-1)
                IDLAM(LKNT,2)= 1 +2*(K-1)
                IDLAM(LKNT,3)= 0
                XLAM(LKNT)=0
                IF (IMSS(52).NE.0) THEN
C...Use massive top quark
                  IF (IDLAM(LKNT,1).EQ.-6) THEN
                    RM2=3*RVLAMP(I,J,K)**2*SFMIX(KFSM,2*ICNT-1)**2
     &                   * (SM-SMT)
                    XLAM(LKNT) =
     &                   PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,3)
C...If no top quark, all decay products massless
                  ELSE
                    RM2=3*RVLAMP(I,J,K)**2*SFMIX(KFSM,2*ICNT-1)**2 * SM
                    XLAM(LKNT) =
     &                   PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
                  ENDIF
C...KINEMATICS CHECK
                  IF (XLAM(LKNT).EQ.0D0) THEN
                    LKNT=LKNT-1
                  ENDIF
                ENDIF
  140         CONTINUE
  150       CONTINUE
          ENDIF
C * SNEUTRINO -> LEPTON+ + LEPTON- and DBAR + D
C...No right-handed neutrinos
          IF(ICNT.EQ.1) THEN
            IF(KFSM.EQ.12.OR.KFSM.EQ.14.OR.KFSM.EQ.16) THEN
              J=INT((KFSM-10)/2)
              DO 170 I=1,3
                DO 160 K=1,3
                  IF (I.NE.J) THEN
C...~nu_J -> lepton+_I + lepton-_K
                    LKNT = LKNT+1
                    IDLAM(LKNT,1)=-11 -2*(I-1)
                    IDLAM(LKNT,2)= 11 +2*(K-1)
                    IDLAM(LKNT,3)=  0
                    XLAM(LKNT)=0D0
                    RM2=RVLAM(I,J,K)**2 * SM
                    IF (IMSS(51).NE.0) XLAM(LKNT) =
     &                   PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
C...KINEMATICS CHECK
                    IF (XLAM(LKNT).EQ.0D0) THEN
                      LKNT=LKNT-1
                    ENDIF
                  ENDIF
  160           CONTINUE
  170         CONTINUE
C...~nu_I -> dbar_J + d_K
              I=INT((KFSM-10)/2)
              DO 190 J=1,3
                DO 180 K=1,3
                  LKNT = LKNT+1
                  IDLAM(LKNT,1)=-1 -2*(J-1)
                  IDLAM(LKNT,2)= 1 +2*(K-1)
                  IDLAM(LKNT,3)= 0
                  XLAM(LKNT)=0D0
                  RM2=3*RVLAMP(I,J,K)**2 * SM
                  IF (IMSS(52).NE.0) XLAM(LKNT) =
     &                 PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
C...KINEMATICS CHECK
                  IF (XLAM(LKNT).EQ.0D0) THEN
                    LKNT=LKNT-1
                  ENDIF
  180           CONTINUE
  190         CONTINUE
            ENDIF
          ENDIF
C * SDOWN -> NU(BAR) + D and LEPTON- + U
          IF(ICNT.NE.0.AND.(KFSM.EQ.1.OR.KFSM.EQ.3.OR.KFSM.EQ.5)) THEN
            J=INT((KFSM+1)/2)
            DO 210 I=1,3
              DO 200 K=1,3
C...~d_J -> nu_Ibar + d_K
                LKNT = LKNT+1
                IDLAM(LKNT,1)=-12 -2*(I-1)
                IDLAM(LKNT,2)=  1 +2*(K-1)
                IDLAM(LKNT,3)=  0
                XLAM(LKNT)=0D0
                RM2=RVLAMP(I,J,K)**2*SFMIX(KFSM,2*ICNT-1)**2 * SM
                IF (IMSS(52).NE.0) XLAM(LKNT) =
     &               PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
C...KINEMATICS CHECK
                IF (XLAM(LKNT).EQ.0D0) THEN
                  LKNT=LKNT-1
                ENDIF
  200         CONTINUE
  210       CONTINUE
            K=INT((KFSM+1)/2)
            DO 240 I=1,3
              DO 230 J=1,3
C...~d_K -> nu_I + d_J
                LKNT = LKNT+1
                IDLAM(LKNT,1)= 12 +2*(I-1)
                IDLAM(LKNT,2)=  1 +2*(J-1)
                IDLAM(LKNT,3)=  0
                XLAM(LKNT)=0D0
                RM2=RVLAMP(I,J,K)**2*SFMIX(KFSM,2*ICNT)**2 * SM
                IF (IMSS(52).NE.0) XLAM(LKNT) =
     &               PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
C...KINEMATICS CHECK
                IF (XLAM(LKNT).EQ.0D0) THEN
                  LKNT=LKNT-1
                ENDIF
C...~d_K -> lepton_I- + u_J
  220           LKNT = LKNT+1
                IDLAM(LKNT,1)= 11 +2*(I-1)
                IDLAM(LKNT,2)=  2 +2*(J-1)
                IDLAM(LKNT,3)=  0
                XLAM(LKNT)=0D0
                IF (IMSS(52).NE.0) THEN
C...Use massive top quark
                  IF (IDLAM(LKNT,2).EQ.6) THEN
                    RM2=RVLAMP(I,J,K)**2*SFMIX(KFSM,2*ICNT)**2*(SM-SMT)
                    XLAM(LKNT) =
     &                   PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,2)
C...If no top quark, all decay products massless
                  ELSE
                    RM2=RVLAMP(I,J,K)**2*SFMIX(KFSM,2*ICNT)**2 * SM
                    XLAM(LKNT) =
     &                   PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
                  ENDIF
C...KINEMATICS CHECK
                  IF (XLAM(LKNT).EQ.0D0) THEN
                    LKNT=LKNT-1
                  ENDIF
                ENDIF
  230         CONTINUE
  240       CONTINUE
          ENDIF
C * SUP -> LEPTON+ + D
          IF(ICNT.NE.0.AND.(KFSM.EQ.2.OR.KFSM.EQ.4.OR.KFSM.EQ.6)) THEN
            J=NINT(KFSM/2.)
            DO 260 I=1,3
              DO 250 K=1,3
C...~u_J -> lepton_I+ + d_K
                LKNT = LKNT+1
                IDLAM(LKNT,1)=-11 -2*(I-1)
                IDLAM(LKNT,2)=  1 +2*(K-1)
                IDLAM(LKNT,3)=  0
                XLAM(LKNT)=0D0
                RM2=RVLAMP(I,J,K)**2*SFMIX(KFSM,2*ICNT-1)**2 * SM
                IF (IMSS(52).NE.0) XLAM(LKNT) =
     &               PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
C...KINEMATICS CHECK
                IF (XLAM(LKNT).EQ.0D0) THEN
                  LKNT=LKNT-1
                ENDIF
  250         CONTINUE
  260       CONTINUE
          ENDIF
        ENDIF
C...BARYON NUMBER VIOLATING DECAYS
        IF (IMSS(53).GE.1) THEN
C * SUP -> DBAR + DBAR
          IF(ICNT.NE.0.AND.(KFSM.EQ.2.OR.KFSM.EQ.4.OR.KFSM.EQ.6)) THEN
            I = KFSM/2
            DO 280 J=1,3
              DO 270 K=1,3
C...~u_I -> dbar_J + dbar_K
                IF (J.LT.K) THEN
C...(anti-) symmetry J <-> K.
                  LKNT = LKNT + 1
                  IDLAM(LKNT,1) = -1 -2*(J-1)
                  IDLAM(LKNT,2) = -1 -2*(K-1)
                  IDLAM(LKNT,3) =  0
                  XLAM(LKNT)    =  0D0
                  RM2 = 2.*(RVLAMB(I,J,K)**2)
     &                 * SFMIX(KFSM,2*ICNT)**2 * SM
                  XLAM(LKNT)    =
     &                 PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
C...KINEMATICS CHECK
                  IF (XLAM(LKNT).EQ.0D0) THEN
                    LKNT = LKNT-1
                  ENDIF
                ENDIF
  270         CONTINUE
  280       CONTINUE
          ENDIF
C * SDOWN -> UBAR + DBAR
          IF(ICNT.NE.0.AND.(KFSM.EQ.1.OR.KFSM.EQ.3.OR.KFSM.EQ.5)) THEN
            K=(KFSM+1)/2
            DO 300 I=1,3
              DO 290 J=1,3
C...LAMB coupling antisymmetric in J and K.
                IF (J.NE.K) THEN
C...~d_K -> ubar_I + dbar_K
                  LKNT = LKNT + 1
                  IDLAM(LKNT,1)= -2 -2*(I-1)
                  IDLAM(LKNT,2)= -1 -2*(J-1)
                  IDLAM(LKNT,3)=  0
                  XLAM(LKNT)=0D0
C...Use massive top quark
                  IF (IDLAM(LKNT,1).EQ.-6) THEN
                    RM2=2*RVLAMB(I,J,K)**2*SFMIX(KFSM,2*ICNT)**2*(SM-SMT
     &                   )
                    XLAM(LKNT) =
     &                   PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,3)
C...If no top quark, all decay products massless
                  ELSE
                    RM2=2*RVLAMB(I,J,K)**2*SFMIX(KFSM,2*ICNT)**2 * SM
                    XLAM(LKNT) =
     &                   PYRVSB(KFIN,IDLAM(LKNT,1),IDLAM(LKNT,2),RM2,4)
                  ENDIF
C...KINEMATICS CHECK
                  IF (XLAM(LKNT).EQ.0D0) THEN
                    LKNT=LKNT-1
                  ENDIF
                ENDIF
  290         CONTINUE
  300       CONTINUE
          ENDIF
        ENDIF
      ENDIF
 
      RETURN
      END
