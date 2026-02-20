 
C*********************************************************************
 
C...PYGAUS
C...Integration by adaptive Gaussian quadrature.
C...Adapted from the CERNLIB DGAUSS routine by K.S. Kolbig.
 
      FUNCTION PYGAUS(F, A, B, EPS)
 
C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
 
C...Local declarations.
      EXTERNAL F
      DOUBLE PRECISION F,W(12), X(12)
      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/
 
C...The Gaussian quadrature algorithm.
      H = 0D0
      IF(B .EQ. A) GOTO 140
      CONST = 5D-3 / ABS(B-A)
      BB = A
  100 CONTINUE
      AA = BB
      BB = B
  110 CONTINUE
      C1 = 0.5D0*(BB+AA)
      C2 = 0.5D0*(BB-AA)
      S8 = 0D0
      DO 120 I = 1, 4
        U = C2*X(I)
        S8 = S8 + W(I) * (F(C1+U) + F(C1-U))
  120 CONTINUE
      S16 = 0D0
      DO 130 I = 5, 12
        U = C2*X(I)
        S16 = S16 + W(I) * (F(C1+U) + F(C1-U))
  130 CONTINUE
      S16 = C2*S16
      IF(DABS(S16-C2*S8) .LE. EPS*(1D0+DABS(S16))) THEN
        H = H + S16
        IF(BB .NE. B) GOTO 100
      ELSE
        BB = C1
        IF(1D0 + CONST*ABS(C2) .NE. 1D0) GOTO 110
        H = 0D0
        CALL PYERRM(18,'(PYGAUS:) too high accuracy required')
        GOTO 140
      ENDIF
  140 CONTINUE
      PYGAUS = H
 
      RETURN
      END
