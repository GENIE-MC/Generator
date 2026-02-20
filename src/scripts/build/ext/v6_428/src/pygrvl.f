 
C*********************************************************************
 
C...PYGRVL
C...Gives the GRV 94 L (leading order) parton distribution function set
C...in parametrized form.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      SUBROUTINE PYGRVL (X, Q2, UV, DV, DEL, UDB, SB, CHM, BOT, GL)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION (A - Z)
 
C...Common expressions.
      MU2  = 0.23D0
      LAM2 = 0.2322D0 * 0.2322D0
      S  = LOG (LOG(Q2/LAM2) / LOG(MU2/LAM2))
      DS = SQRT (S)
      S2 = S * S
      S3 = S2 * S
 
C...uv :
      NU  =  2.284D0 + 0.802D0 * S + 0.055D0 * S2
      AKU =  0.590D0 - 0.024D0 * S
      BKU =  0.131D0 + 0.063D0 * S
      AU  = -0.449D0 - 0.138D0 * S - 0.076D0 * S2
      BU  =  0.213D0 + 2.669D0 * S - 0.728D0 * S2
      CU  =  8.854D0 - 9.135D0 * S + 1.979D0 * S2
      DU  =  2.997D0 + 0.753D0 * S - 0.076D0 * S2
      UV  = PYGRVV (X, NU, AKU, BKU, AU, BU, CU, DU)
 
C...dv :
      ND  =  0.371D0 + 0.083D0 * S + 0.039D0 * S2
      AKD =  0.376D0
      BKD =  0.486D0 + 0.062D0 * S
      AD  = -0.509D0 + 3.310D0 * S - 1.248D0 * S2
      BD  =  12.41D0 - 10.52D0 * S + 2.267D0 * S2
      CD  =  6.373D0 - 6.208D0 * S + 1.418D0 * S2
      DD  =  3.691D0 + 0.799D0 * S - 0.071D0 * S2
      DV  = PYGRVV (X, ND, AKD, BKD, AD, BD, CD, DD)
 
C...del :
      NE  =  0.082D0 + 0.014D0 * S + 0.008D0 * S2
      AKE =  0.409D0 - 0.005D0 * S
      BKE =  0.799D0 + 0.071D0 * S
      AE  = -38.07D0 + 36.13D0 * S - 0.656D0 * S2
      BE  =  90.31D0 - 74.15D0 * S + 7.645D0 * S2
      CE  =  0.0D0
      DE  =  7.486D0 + 1.217D0 * S - 0.159D0 * S2
      DEL = PYGRVV (X, NE, AKE, BKE, AE, BE, CE, DE)
 
C...udb :
      ALX =  1.451D0
      BEX =  0.271D0
      AKX =  0.410D0 - 0.232D0 * S
      BKX =  0.534D0 - 0.457D0 * S
      AGX =  0.890D0 - 0.140D0 * S
      BGX = -0.981D0
      CX  =  0.320D0 + 0.683D0 * S
      DX  =  4.752D0 + 1.164D0 * S + 0.286D0 * S2
      EX  =  4.119D0 + 1.713D0 * S
      ESX =  0.682D0 + 2.978D0 * S
      UDB = PYGRVW (X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX,
     & DX, EX, ESX)
 
C...sb :
      STS =  0D0
      ALS =  0.914D0
      BES =  0.577D0
      AKS =  1.798D0 - 0.596D0 * S
      AS  = -5.548D0 + 3.669D0 * DS - 0.616D0 * S
      BS  =  18.92D0 - 16.73D0 * DS + 5.168D0 * S
      DST =  6.379D0 - 0.350D0 * S  + 0.142D0 * S2
      EST =  3.981D0 + 1.638D0 * S
      ESS =  6.402D0
      SB  = PYGRVS (X, S, STS, ALS, BES, AKS, AS, BS, DST, EST, ESS)
 
C...cb :
      STC =  0.888D0
      ALC =  1.01D0
      BEC =  0.37D0
      AKC =  0D0
      AC  =  0D0
      BC  =  4.24D0  - 0.804D0 * S
      DCT =  3.46D0  - 1.076D0 * S
      ECT =  4.61D0  + 1.49D0  * S
      ESC =  2.555D0 + 1.961D0 * S
      CHM = PYGRVS (X, S, STC, ALC, BEC, AKC, AC, BC, DCT, ECT, ESC)
 
C...bb :
      STB =  1.351D0
      ALB =  1.00D0
      BEB =  0.51D0
      AKB =  0D0
      AB  =  0D0
      BB  =  1.848D0
      DBT =  2.929D0 + 1.396D0 * S
      EBT =  4.71D0  + 1.514D0 * S
      ESB =  4.02D0  + 1.239D0 * S
      BOT = PYGRVS (X, S, STB, ALB, BEB, AKB, AB, BB, DBT, EBT, ESB)
 
C...gl :
      ALG =  0.524D0
      BEG =  1.088D0
      AKG =  1.742D0 - 0.930D0 * S
      BKG =                         - 0.399D0 * S2
      AG  =  7.486D0 - 2.185D0 * S
      BG  =  16.69D0 - 22.74D0 * S  + 5.779D0 * S2
      CG  = -25.59D0 + 29.71D0 * S  - 7.296D0 * S2
      DG  =  2.792D0 + 2.215D0 * S  + 0.422D0 * S2 - 0.104D0 * S3
      EG  =  0.807D0 + 2.005D0 * S
      ESG =  3.841D0 + 0.316D0 * S
      GL  = PYGRVW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG,
     & DG, EG, ESG)
 
      RETURN
      END
