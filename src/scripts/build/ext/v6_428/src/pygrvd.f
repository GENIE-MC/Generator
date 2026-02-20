 
C*********************************************************************
 
C...PYGRVD
C...Gives the GRV 94 D (DIS) parton distribution function set
C...in parametrized form.
C...Authors: M. Glueck, E. Reya and A. Vogt.
 
      SUBROUTINE PYGRVD (X, Q2, UV, DV, DEL, UDB, SB, CHM, BOT, GL)
 
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
      NU  =  2.484D0 + 0.116D0 * S + 0.093D0 * S2
      AKU =  0.563D0 - 0.025D0 * S
      BKU =  0.054D0 + 0.154D0 * S
      AU  = -0.326D0 - 0.058D0 * S - 0.135D0 * S2
      BU  = -3.322D0 + 8.259D0 * S - 3.119D0 * S2 + 0.291D0 * S3
      CU  =  11.52D0 - 12.99D0 * S + 3.161D0 * S2
      DU  =  2.808D0 + 1.400D0 * S - 0.557D0 * S2 + 0.119D0 * S3
      UV  = PYGRVV (X, NU, AKU, BKU, AU, BU, CU, DU)
 
C...dv :
      ND  =  0.156D0 - 0.017D0 * S
      AKD =  0.299D0 - 0.022D0 * S
      BKD =  0.259D0 - 0.015D0 * S
      AD  =  3.445D0 + 1.278D0 * S + 0.326D0 * S2
      BD  = -6.934D0 + 37.45D0 * S - 18.95D0 * S2 + 1.463D0 * S3
      CD  =  55.45D0 - 69.92D0 * S + 20.78D0 * S2
      DD  =  3.577D0 + 1.441D0 * S - 0.683D0 * S2 + 0.179D0 * S3
      DV  = PYGRVV (X, ND, AKD, BKD, AD, BD, CD, DD)
 
C...del :
      NE  =  0.099D0 + 0.019D0 * S + 0.002D0 * S2
      AKE =  0.419D0 - 0.013D0 * S
      BKE =  1.064D0 - 0.038D0 * S
      AE  = -44.00D0 + 98.70D0 * S - 14.79D0 * S2
      BE  =  28.59D0 - 40.94D0 * S - 13.66D0 * S2 + 2.523D0 * S3
      CE  =  84.57D0 - 108.8D0 * S + 31.52D0 * S2
      DE  =  7.469D0 + 2.480D0 * S - 0.866D0 * S2
      DEL = PYGRVV (X, NE, AKE, BKE, AE, BE, CE, DE)
 
C...udb :
      ALX =  1.215D0
      BEX =  0.466D0
      AKX =  0.326D0 + 0.150D0 * S
      BKX =  0.956D0 + 0.405D0 * S
      AGX =  0.272D0
      BGX =  3.794D0 - 2.359D0 * DS
      CX  =  2.014D0
      DX  =  7.941D0 + 0.534D0 * DS - 0.940D0 * S + 0.410D0 * S2
      EX  =  3.049D0 + 1.597D0 * S
      ESX =  4.396D0 - 4.594D0 * DS + 3.268D0 * S
      UDB = PYGRVW (X, S, ALX, BEX, AKX, BKX, AGX, BGX, CX,
     & DX, EX, ESX)
 
C...sb :
      STS =  0D0
      ALS =  0.175D0
      BES =  0.344D0
      AKS =  1.415D0 - 0.641D0 * DS
      AS  =  0.580D0 - 9.763D0 * DS + 6.795D0 * S  - 0.558D0 * S2
      BS  =  5.617D0 + 5.709D0 * DS - 3.972D0 * S
      DST =  13.78D0 - 9.581D0 * S  + 5.370D0 * S2 - 0.996D0 * S3
      EST =  4.546D0 + 0.372D0 * S2
      ESS =  5.053D0 - 1.070D0 * S  + 0.805D0 * S2
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
      ALG =  1.258D0
      BEG =  1.846D0
      AKG =  2.423D0
      BKG =  2.427D0 + 1.311D0 * S  - 0.153D0 * S2
      AG  =  25.09D0 - 7.935D0 * S
      BG  = -14.84D0 - 124.3D0 * DS + 72.18D0 * S
      CG  =  590.3D0 - 173.8D0 * S
      DG  =  5.196D0 + 1.857D0 * S
      EG  = -1.648D0 + 3.988D0 * S  - 0.432D0 * S2
      ESG =  3.232D0 - 0.542D0 * S
      GL  = PYGRVW (X, S, ALG, BEG, AKG, BKG, AG, BG, CG, DG, EG, ESG)
 
      RETURN
      END
