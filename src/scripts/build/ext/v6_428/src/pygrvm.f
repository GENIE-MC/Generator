 
C*********************************************************************
 
C...PYGRVM
C...Gives the GRV 94 M (MSbar) parton distribution function set
C...in parametrized form.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      SUBROUTINE PYGRVM (X, Q2, UV, DV, DEL, UDB, SB, CHM, BOT, GL)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION (A - Z)
 
C...Common expressions.
      MU2  = 0.34D0
      LAM2 = 0.248D0 * 0.248D0
      S  = LOG (LOG(Q2/LAM2) / LOG(MU2/LAM2))
      DS = SQRT (S)
      S2 = S * S
      S3 = S2 * S
 
C...uv :
      NU  =  1.304D0 + 0.863D0 * S
      AKU =  0.558D0 - 0.020D0 * S
      BKU =          0.183D0 * S
      AU  = -0.113D0 + 0.283D0 * S - 0.321D0 * S2
      BU  =  6.843D0 - 5.089D0 * S + 2.647D0 * S2 - 0.527D0 * S3
      CU  =  7.771D0 - 10.09D0 * S + 2.630D0 * S2
      DU  =  3.315D0 + 1.145D0 * S - 0.583D0 * S2 + 0.154D0 * S3
      UV  = PYGRVV (X, NU, AKU, BKU, AU, BU, CU, DU)
 
C...dv :
      ND  =  0.102D0 - 0.017D0 * S + 0.005D0 * S2
      AKD =  0.270D0 - 0.019D0 * S
      BKD =  0.260D0
      AD  =  2.393D0 + 6.228D0 * S - 0.881D0 * S2
      BD  =  46.06D0 + 4.673D0 * S - 14.98D0 * S2 + 1.331D0 * S3
      CD  =  17.83D0 - 53.47D0 * S + 21.24D0 * S2
      DD  =  4.081D0 + 0.976D0 * S - 0.485D0 * S2 + 0.152D0 * S3
      DV  = PYGRVV (X, ND, AKD, BKD, AD, BD, CD, DD)
 
C...del :
      NE  =  0.070D0 + 0.042D0 * S - 0.011D0 * S2 + 0.004D0 * S3
      AKE =  0.409D0 - 0.007D0 * S
      BKE =  0.782D0 + 0.082D0 * S
      AE  = -29.65D0 + 26.49D0 * S + 5.429D0 * S2
      BE  =  90.20D0 - 74.97D0 * S + 4.526D0 * S2
      CE  =  0.0D0
      DE  =  8.122D0 + 2.120D0 * S - 1.088D0 * S2 + 0.231D0 * S3
      DEL = PYGRVV (X, NE, AKE, BKE, AE, BE, CE, DE)
 
C...udb :
      ALX =  0.877D0
      BEX =  0.561D0
      AKX =  0.275D0
      BKX =  0.0D0
      AGX =  0.997D0
      BGX =  3.210D0 - 1.866D0 * S
      CX  =  7.300D0
      DX  =  9.010D0 + 0.896D0 * DS + 0.222D0 * S2
      EX  =  3.077D0 + 1.446D0 * S
      ESX =  3.173D0 - 2.445D0 * DS + 2.207D0 * S
      UDB = PYGRVW (X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX,
     & DX, EX, ESX)
 
C...sb :
      STS =  0D0
      ALS =  0.756D0
      BES =  0.216D0
      AKS =  1.690D0 + 0.650D0 * DS - 0.922D0 * S
      AS  = -4.329D0 + 1.131D0 * S
      BS  =  9.568D0 - 1.744D0 * S
      DST =  9.377D0 + 1.088D0 * DS - 1.320D0 * S + 0.130D0 * S2
      EST =  3.031D0 + 1.639D0 * S
      ESS =  5.837D0 + 0.815D0 * S
      SB  = PYGRVS (X, S, STS, ALS, BES, AKS, AS, BS, DST, EST, ESS)
 
C...cb :
      STC =  0.820D0
      ALC =  0.98D0
      BEC =  0D0
      AKC = -0.625D0 - 0.523D0 * S
      AC  =  0D0
      BC  =  1.896D0 + 1.616D0 * S
      DCT =  4.12D0  + 0.683D0 * S
      ECT =  4.36D0  + 1.328D0 * S
      ESC =  0.677D0 + 0.679D0 * S
      CHM = PYGRVS (X, S, STC, ALC, BEC, AKC, AC, BC, DCT, ECT, ESC)
 
C...bb :
      STB =  1.297D0
      ALB =  0.99D0
      BEB =  0D0
      AKB =          - 0.193D0 * S
      AB  =  0D0
      BB  =  0D0
      DBT =  3.447D0 + 0.927D0 * S
      EBT =  4.68D0  + 1.259D0 * S
      ESB =  1.892D0 + 2.199D0 * S
      BOT = PYGRVS (X, S, STB, ALB, BEB, AKB, AB, BB, DBT, EBT, ESB)
 
C...gl :
       ALG =  1.014D0
       BEG =  1.738D0
       AKG =  1.724D0 + 0.157D0 * S
       BKG =  0.800D0 + 1.016D0 * S
       AG  =  7.517D0 - 2.547D0 * S
       BG  =  34.09D0 - 52.21D0 * DS + 17.47D0 * S
       CG  =  4.039D0 + 1.491D0 * S
       DG  =  3.404D0 + 0.830D0 * S
       EG  = -1.112D0 + 3.438D0 * S  - 0.302D0 * S2
       ESG =  3.256D0 - 0.436D0 * S
       GL  = PYGRVW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
 
       RETURN
       END
