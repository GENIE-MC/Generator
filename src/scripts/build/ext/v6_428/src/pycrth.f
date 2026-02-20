 
C*********************************************************************
 
C...PYCRTH
C...Auxiliary to PYEICG.
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE ORTHES, NUM. MATH. 12, 349-368(1968)
C     BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A COMPLEX GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     UNITARY SIMILARITY TRANSFORMATIONS.
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
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX INPUT MATRIX.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.  INFORMATION
C          ABOUT THE UNITARY TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLES UNDER THE
C          HESSENBERG MATRIX.
C
C        ORTR AND ORTI CONTAIN FURTHER INFORMATION ABOUT THE
C          TRANSFORMATIONS.  ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     CALLS PYTHAG FOR  DSQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
 
      SUBROUTINE PYCRTH(NM,N,LOW,IGH,AR,AI,ORTR,ORTI)
 
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      DOUBLE PRECISION AR(5,5),AI(5,5),ORTR(5),ORTI(5)
      DOUBLE PRECISION F,G,H,FI,FR,SCALE,PYTHAG
 
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GOTO 210
C
      DO 200 M = KP1, LA
         H = 0.0D0
         ORTR(M) = 0.0D0
         ORTI(M) = 0.0D0
         SCALE = 0.0D0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
         DO 100 I = M, IGH
  100    SCALE = SCALE + DABS(AR(I,M-1)) + DABS(AI(I,M-1))
C
         IF (SCALE .EQ. 0.0D0) GOTO 200
         MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
         DO 110 II = M, IGH
            I = MP - II
            ORTR(I) = AR(I,M-1) / SCALE
            ORTI(I) = AI(I,M-1) / SCALE
            H = H + ORTR(I) * ORTR(I) + ORTI(I) * ORTI(I)
  110    CONTINUE
C
         G = DSQRT(H)
         F = PYTHAG(ORTR(M),ORTI(M))
         IF (F .EQ. 0.0D0) GOTO 120
         H = H + F * G
         G = G / F
         ORTR(M) = (1.0D0 + G) * ORTR(M)
         ORTI(M) = (1.0D0 + G) * ORTI(M)
         GOTO 130
C
  120    ORTR(M) = G
         AR(M,M-1) = SCALE
C     .......... FORM (I-(U*UT)/H) * A ..........
  130    DO 160 J = M, N
            FR = 0.0D0
            FI = 0.0D0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
            DO 140 II = M, IGH
               I = MP - II
               FR = FR + ORTR(I) * AR(I,J) + ORTI(I) * AI(I,J)
               FI = FI + ORTR(I) * AI(I,J) - ORTI(I) * AR(I,J)
  140       CONTINUE
C
            FR = FR / H
            FI = FI / H
C
            DO 150 I = M, IGH
               AR(I,J) = AR(I,J) - FR * ORTR(I) + FI * ORTI(I)
               AI(I,J) = AI(I,J) - FR * ORTI(I) - FI * ORTR(I)
  150       CONTINUE
C
  160    CONTINUE
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
         DO 190 I = 1, IGH
            FR = 0.0D0
            FI = 0.0D0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
            DO 170 JJ = M, IGH
               J = MP - JJ
               FR = FR + ORTR(J) * AR(I,J) - ORTI(J) * AI(I,J)
               FI = FI + ORTR(J) * AI(I,J) + ORTI(J) * AR(I,J)
  170       CONTINUE
C
            FR = FR / H
            FI = FI / H
C
            DO 180 J = M, IGH
               AR(I,J) = AR(I,J) - FR * ORTR(J) - FI * ORTI(J)
               AI(I,J) = AI(I,J) + FR * ORTI(J) - FI * ORTR(J)
  180       CONTINUE
C
  190    CONTINUE
C
         ORTR(M) = SCALE * ORTR(M)
         ORTI(M) = SCALE * ORTI(M)
         AR(M,M-1) = -G * AR(M,M-1)
         AI(M,M-1) = -G * AI(M,M-1)
  200 CONTINUE
C
  210 RETURN
      END
