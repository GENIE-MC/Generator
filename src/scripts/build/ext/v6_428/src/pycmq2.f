 
C*********************************************************************
 
C...PYCMQ2
C...Auxiliary to PYEICG.
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR2, NUM. MATH. 16, 181-204(1970) BY PETERS
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A COMPLEX UPPER HESSENBERG MATRIX BY THE QR
C     METHOD.  THE EIGENVECTORS OF A COMPLEX GENERAL MATRIX
C     CAN ALSO BE FOUND IF  CORTH  HAS BEEN USED TO REDUCE
C     THIS GENERAL MATRIX TO HESSENBERG FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        ORTR AND ORTI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  CORTH, IF PERFORMED.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, SET ORTR(J) AND
C          ORTI(J) TO 0.0D0 FOR THESE ELEMENTS.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN FURTHER
C          INFORMATION ABOUT THE TRANSFORMATIONS WHICH WERE USED IN THE
C          REDUCTION BY  CORTH, IF PERFORMED.  IF THE EIGENVECTORS OF
C          THE HESSENBERG MATRIX ARE DESIRED, THESE ELEMENTS MAY BE
C          ARBITRARY.
C
C     ON OUTPUT
C
C        ORTR, ORTI, AND THE UPPER HESSENBERG PORTIONS OF HR AND HI
C          HAVE BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS.  THE EIGENVECTORS
C          ARE UNNORMALIZED.  IF AN ERROR EXIT IS MADE, NONE OF
C          THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS PYCDIV FOR COMPLEX DIVISION.
C     CALLS PYCSRT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED OCTOBER 1989.
C
C  MESHED OVERFLOW CONTROL WITH VECTORS OF ISOLATED ROOTS (10/19/89 BSG)
C  MESHED OVERFLOW CONTROL WITH TRIANGULAR MULTIPLY (10/30/89 BSG)
C
 
      SUBROUTINE PYCMQ2(NM,N,LOW,IGH,ORTR,ORTI,HR,HI,WR,WI,ZR,ZI,IERR)
 
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,NM,NN,IGH,IP1,
     X        ITN,ITS,LOW,LP1,ENM1,IEND,IERR
      DOUBLE PRECISION HR(5,5),HI(5,5),WR(5),WI(5),ZR(5,5),ZI(5,5),
     X       ORTR(5),ORTI(5)
      DOUBLE PRECISION SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2,
     X       PYTHAG
 
      IERR = 0
C     .......... INITIALIZE EIGENVECTOR MATRIX ..........
      DO 110 J = 1, N
C
         DO 100 I = 1, N
            ZR(I,J) = 0.0D0
            ZI(I,J) = 0.0D0
  100    CONTINUE
         ZR(J,J) = 1.0D0
  110 CONTINUE
C     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
C                FROM THE INFORMATION LEFT BY CORTH ..........
      IEND = IGH - LOW - 1
      IF (IEND.LT.0) GOTO 220
      IF (IEND.EQ.0) GOTO 170
C     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 160 II = 1, IEND
         I = IGH - II
         IF (ORTR(I) .EQ. 0.0D0 .AND. ORTI(I) .EQ. 0.0D0) GOTO 160
         IF (HR(I,I-1) .EQ. 0.0D0 .AND. HI(I,I-1) .EQ. 0.0D0) GOTO 160
C     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
         NORM = HR(I,I-1) * ORTR(I) + HI(I,I-1) * ORTI(I)
         IP1 = I + 1
C
         DO 120 K = IP1, IGH
            ORTR(K) = HR(K,I-1)
            ORTI(K) = HI(K,I-1)
  120    CONTINUE
C
         DO 150 J = I, IGH
            SR = 0.0D0
            SI = 0.0D0
C
            DO 130 K = I, IGH
               SR = SR + ORTR(K) * ZR(K,J) + ORTI(K) * ZI(K,J)
               SI = SI + ORTR(K) * ZI(K,J) - ORTI(K) * ZR(K,J)
  130       CONTINUE
C
            SR = SR / NORM
            SI = SI / NORM
C
            DO 140 K = I, IGH
               ZR(K,J) = ZR(K,J) + SR * ORTR(K) - SI * ORTI(K)
               ZI(K,J) = ZI(K,J) + SR * ORTI(K) + SI * ORTR(K)
  140       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  170 L = LOW + 1
C
      DO 210 I = L, IGH
         LL = MIN0(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0D0) GOTO 210
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0D0
C
         DO 180 J = I, N
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  180    CONTINUE
C
         DO 190 J = 1, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  190    CONTINUE
C
         DO 200 J = LOW, IGH
            SI = YR * ZI(J,I) + YI * ZR(J,I)
            ZR(J,I) = YR * ZR(J,I) - YI * ZI(J,I)
            ZI(J,I) = SI
  200    CONTINUE
C
  210 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  220 DO 230 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GOTO 230
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  230 CONTINUE
C
      EN = IGH
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  240 IF (EN .LT. LOW) GOTO 430
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  250 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GOTO 270
         TST1 = DABS(HR(L-1,L-1)) + DABS(HI(L-1,L-1))
     X            + DABS(HR(L,L)) + DABS(HI(L,L))
         TST2 = TST1 + DABS(HR(L,L-1))
         IF (TST2 .EQ. TST1) GOTO 270
  260 CONTINUE
C     .......... FORM SHIFT ..........
  270 IF (L .EQ. EN) GOTO 420
      IF (ITN .EQ. 0) GOTO 550
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GOTO 290
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0D0 .AND. XI .EQ. 0.0D0) GOTO 300
      YR = (HR(ENM1,ENM1) - SR) / 2.0D0
      YI = (HI(ENM1,ENM1) - SI) / 2.0D0
      CALL PYCSRT(YR**2-YI**2+XR,2.0D0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0D0) GOTO 280
      ZZR = -ZZR
      ZZI = -ZZI
  280 CALL PYCDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GOTO 300
C     .......... FORM EXCEPTIONAL SHIFT ..........
  290 SR = DABS(HR(EN,ENM1)) + DABS(HR(ENM1,EN-2))
      SI = 0.0D0
C
  300 DO 310 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  310 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
      DO 330 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0D0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0D0
         HI(I,I-1) = SR / NORM
C
         DO 320 J = I, N
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  320    CONTINUE
C
  330 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0D0) GOTO 350
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0D0
      IF (EN .EQ. N) GOTO 350
      IP1 = EN + 1
C
      DO 340 J = IP1, N
         YR = HR(EN,J)
         YI = HI(EN,J)
         HR(EN,J) = SR * YR + SI * YI
         HI(EN,J) = SR * YI - SI * YR
  340 CONTINUE
C     .......... INVERSE OPERATION (COLUMNS) ..........
  350 DO 390 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 370 I = 1, J
            YR = HR(I,J-1)
            YI = 0.0D0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GOTO 360
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  360       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  370    CONTINUE
C
         DO 380 I = LOW, IGH
            YR = ZR(I,J-1)
            YI = ZI(I,J-1)
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
            ZR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            ZI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
            ZR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            ZI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  380    CONTINUE
C
  390 CONTINUE
C
      IF (SI .EQ. 0.0D0) GOTO 250
C
      DO 400 I = 1, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  400 CONTINUE
C
      DO 410 I = LOW, IGH
         YR = ZR(I,EN)
         YI = ZI(I,EN)
         ZR(I,EN) = SR * YR - SI * YI
         ZI(I,EN) = SR * YI + SI * YR
  410 CONTINUE
C
      GOTO 250
C     .......... A ROOT FOUND ..........
  420 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GOTO 240
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  430 NORM = 0.0D0
C
      DO 440 I = 1, N
C
         DO 440 J = I, N
            TR = DABS(HR(I,J)) + DABS(HI(I,J))
            IF (TR .GT. NORM) NORM = TR
  440 CONTINUE
C
      IF (N .EQ. 1 .OR. NORM .EQ. 0.0D0) GOTO 560
C     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
      DO 500 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         HR(EN,EN) = 1.0D0
         HI(EN,EN) = 0.0D0
         ENM1 = EN - 1
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 490 II = 1, ENM1
            I = EN - II
            ZZR = 0.0D0
            ZZI = 0.0D0
            IP1 = I + 1
C
            DO 450 J = IP1, EN
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  450       CONTINUE
C
            YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR .NE. 0.0D0 .OR. YI .NE. 0.0D0) GOTO 470
               TST1 = NORM
               YR = TST1
  460          YR = 0.01D0 * YR
               TST2 = NORM + YR
               IF (TST2 .GT. TST1) GOTO 460
  470       CONTINUE
            CALL PYCDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN))
C     .......... OVERFLOW CONTROL ..........
            TR = DABS(HR(I,EN)) + DABS(HI(I,EN))
            IF (TR .EQ. 0.0D0) GOTO 490
            TST1 = TR
            TST2 = TST1 + 1.0D0/TST1
            IF (TST2 .GT. TST1) GOTO 490
            DO 480 J = I, EN
               HR(J,EN) = HR(J,EN)/TR
               HI(J,EN) = HI(J,EN)/TR
  480       CONTINUE
C
  490    CONTINUE
C
  500 CONTINUE
C     .......... END BACKSUBSTITUTION ..........
C     .......... VECTORS OF ISOLATED ROOTS ..........
      DO 520 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GOTO 520
C
         DO 510 J = I, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  510    CONTINUE
C
  520 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW DO -- ..........
      DO 540 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 540 I = LOW, IGH
            ZZR = 0.0D0
            ZZI = 0.0D0
C
            DO 530 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  530       CONTINUE
C
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  540 CONTINUE
C
      GOTO 560
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
  550 IERR = EN
  560 RETURN
      END
