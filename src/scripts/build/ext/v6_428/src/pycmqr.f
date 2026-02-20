 
C*********************************************************************
 
C...PYCMQR
C...Auxiliary to PYEICG.
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR, NUM. MATH. 12, 369-376(1968) BY MARTIN
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A COMPLEX
C     UPPER HESSENBERG MATRIX BY THE QR METHOD.
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
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN
C          INFORMATION ABOUT THE UNITARY TRANSFORMATIONS USED IN
C          THE REDUCTION BY  CORTH, IF PERFORMED.
C
C     ON OUTPUT
C
C        THE UPPER HESSENBERG PORTIONS OF HR AND HI HAVE BEEN
C          DESTROYED.  THEREFORE, THEY MUST BE SAVED BEFORE
C          CALLING  COMQR  IF SUBSEQUENT CALCULATION OF
C          EIGENVECTORS IS TO BE PERFORMED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
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
C     THIS VERSION DATED AUGUST 1983.
C
 
      SUBROUTINE PYCMQR(NM,N,LOW,IGH,HR,HI,WR,WI,IERR)
 
      INTEGER I,J,L,N,EN,LL,NM,IGH,ITN,ITS,LOW,LP1,ENM1,IERR
      DOUBLE PRECISION HR(5,5),HI(5,5),WR(5),WI(5)
      DOUBLE PRECISION SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2,
     X       PYTHAG
 
      IERR = 0
      IF (LOW .EQ. IGH) GOTO 130
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
      L = LOW + 1
C
      DO 120 I = L, IGH
         LL = MIN0(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0D0) GOTO 120
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0D0
C
         DO 100 J = I, IGH
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  100    CONTINUE
C
         DO 110 J = LOW, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  110    CONTINUE
C
  120 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  130 DO 140 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GOTO 140
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  140 CONTINUE
C
      EN = IGH
      TR = 0.0D0
      TI = 0.0D0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  150 IF (EN .LT. LOW) GOTO 320
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW D0 -- ..........
  160 DO 170 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GOTO 180
         TST1 = DABS(HR(L-1,L-1)) + DABS(HI(L-1,L-1))
     X            + DABS(HR(L,L)) + DABS(HI(L,L))
         TST2 = TST1 + DABS(HR(L,L-1))
         IF (TST2 .EQ. TST1) GOTO 180
  170 CONTINUE
C     .......... FORM SHIFT ..........
  180 IF (L .EQ. EN) GOTO 300
      IF (ITN .EQ. 0) GOTO 310
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GOTO 200
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0D0 .AND. XI .EQ. 0.0D0) GOTO 210
      YR = (HR(ENM1,ENM1) - SR) / 2.0D0
      YI = (HI(ENM1,ENM1) - SI) / 2.0D0
      CALL PYCSRT(YR**2-YI**2+XR,2.0D0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0D0) GOTO 190
      ZZR = -ZZR
      ZZI = -ZZI
  190 CALL PYCDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GOTO 210
C     .......... FORM EXCEPTIONAL SHIFT ..........
  200 SR = DABS(HR(EN,ENM1)) + DABS(HR(ENM1,EN-2))
      SI = 0.0D0
C
  210 DO 220 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  220 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
      DO 240 I = LP1, EN
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
         DO 230 J = I, EN
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  230    CONTINUE
C
  240 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0D0) GOTO 250
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0D0
C     .......... INVERSE OPERATION (COLUMNS) ..........
  250 DO 280 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 270 I = L, J
            YR = HR(I,J-1)
            YI = 0.0D0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GOTO 260
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  260       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  270    CONTINUE
C
  280 CONTINUE
C
      IF (SI .EQ. 0.0D0) GOTO 160
C
      DO 290 I = L, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  290 CONTINUE
C
      GOTO 160
C     .......... A ROOT FOUND ..........
  300 WR(EN) = HR(EN,EN) + TR
      WI(EN) = HI(EN,EN) + TI
      EN = ENM1
      GOTO 150
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
  310 IERR = EN
  320 RETURN
      END
