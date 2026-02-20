 
C*********************************************************************
 
C...PYCT5L
C...Auxiliary function for parametrization of CTEQ5L.
C...Author: J. Pumplin 9/99.
 
C...CTEQ5M1 and CTEQ5L Parton Distribution Functions
C...in Parametrized Form
C...            September 15, 1999
C
C...Ref: "GLOBAL QCD ANALYSIS OF PARTON STRUCTURE OF THE NUCLEON:
C...      CTEQ5 PPARTON DISTRIBUTIONS"
C...hep-ph/9903282
 
C...The CTEQ5M1 set given here is an updated version of the original
C...CTEQ5M set posted, in the table version, on the Web page of CTEQ.
C...The differences between CTEQ5M and CTEQ5M1 are insignificant for
C...almost all applications.
C...The improvement is in the QCD evolution which is now more
C...accurate, and which agrees completely with the benchmark work
C...of the HERA 96/97 Workshop.
C...The differences between the parametrized and the corresponding
C...table versions (on which it is based) are of similar order as
C...between the two version.
 
C...!! Because accurate parametrizations over a wide range of (x,Q)
C...is hard to obtain, only the most widely used sets CTEQ5M and
C...CTEQ5L are available in parametrized form for now.
 
C...These parametrizations were obtained by Jon Pumplin.
 
C  Iset   PDF        Description              Alpha_s(Mz)  Lam4  Lam5
C -------------------------------------------------------------------
C   1    CTEQ5M1  Standard NLO MSbar scheme      0.118     326   226
C   3    CTEQ5L   Leading Order                  0.127     192   146
C -------------------------------------------------------------------
C...Note the Qcd-lambda values given for CTEQ5L is for the leading
C...order form of Alpha_s!!  Alpha_s(Mz) gives the absolute
C...calibration.
 
C...The two Iset value are adopted to agree with the standard table
C...versions.
 
C...Range of validity:
C...The range of (x, Q) covered by this parametrization of the QCD
C...evolved parton distributions is 1E-6 < x < 1 ;
C...1.1 GeV < Q < 10 TeV.  Of course, the PDFs are constrained by
C...data only in a subset of that region; and the assumed DGLAP
C...evolution is unlikely to be valid for all of it either.
 
C...The range of (x, Q) used in the CTEQ5 round of global analysis is
C...approximately 0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2 for
C...fixed target experiments; 0.0001 < x < 0.3 from HERA data; and
C...Q^2 up to 40,000 GeV^2 from Tevatron inclusive Jet data.
 
      FUNCTION PYCT5L(IFL,X,Q)
 
C...Double precision declaration.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
 
      PARAMETER (NEX=8, NLF=2)
      DIMENSION AM(0:NEX,0:NLF,-5:2)
      DIMENSION ALFVEC(-5:2), QMAVEC(-5:2)
      DIMENSION MEXVEC(-5:2), MLFVEC(-5:2)
      DIMENSION UT1VEC(-5:2), UT2VEC(-5:2)
      DIMENSION AF(0:NEX)
 
      DATA MEXVEC( 2) / 8 /
      DATA MLFVEC( 2) / 2 /
      DATA UT1VEC( 2) /  0.4971265E+01 /
      DATA UT2VEC( 2) / -0.1105128E+01 /
      DATA ALFVEC( 2) /  0.2987216E+00 /
      DATA QMAVEC( 2) /  0.0000000E+00 /
      DATA (AM( 0,K, 2),K=0, 2)
     & /  0.5292616E+01, -0.2751910E+01, -0.2488990E+01 /
      DATA (AM( 1,K, 2),K=0, 2)
     & /  0.9714424E+00,  0.1011827E-01, -0.1023660E-01 /
      DATA (AM( 2,K, 2),K=0, 2)
     & / -0.1651006E+02,  0.7959721E+01,  0.8810563E+01 /
      DATA (AM( 3,K, 2),K=0, 2)
     & / -0.1643394E+02,  0.5892854E+01,  0.9348874E+01 /
      DATA (AM( 4,K, 2),K=0, 2)
     & /  0.3067422E+02,  0.4235796E+01, -0.5112136E+00 /
      DATA (AM( 5,K, 2),K=0, 2)
     & /  0.2352526E+02, -0.5305168E+01, -0.1169174E+02 /
      DATA (AM( 6,K, 2),K=0, 2)
     & / -0.1095451E+02,  0.3006577E+01,  0.5638136E+01 /
      DATA (AM( 7,K, 2),K=0, 2)
     & / -0.1172251E+02, -0.2183624E+01,  0.4955794E+01 /
      DATA (AM( 8,K, 2),K=0, 2)
     & /  0.1662533E-01,  0.7622870E-02, -0.4895887E-03 /
 
      DATA MEXVEC( 1) / 8 /
      DATA MLFVEC( 1) / 2 /
      DATA UT1VEC( 1) /  0.2612618E+01 /
      DATA UT2VEC( 1) / -0.1258304E+06 /
      DATA ALFVEC( 1) /  0.3407552E+00 /
      DATA QMAVEC( 1) /  0.0000000E+00 /
      DATA (AM( 0,K, 1),K=0, 2)
     & /  0.9905300E+00, -0.4502235E+00,  0.1624441E+00 /
      DATA (AM( 1,K, 1),K=0, 2)
     & /  0.8867534E+00,  0.1630829E-01, -0.4049085E-01 /
      DATA (AM( 2,K, 1),K=0, 2)
     & /  0.8547974E+00,  0.3336301E+00,  0.1371388E+00 /
      DATA (AM( 3,K, 1),K=0, 2)
     & /  0.2941113E+00, -0.1527905E+01,  0.2331879E+00 /
      DATA (AM( 4,K, 1),K=0, 2)
     & /  0.3384235E+02,  0.3715315E+01,  0.8276930E+00 /
      DATA (AM( 5,K, 1),K=0, 2)
     & /  0.6230115E+01,  0.3134639E+01, -0.1729099E+01 /
      DATA (AM( 6,K, 1),K=0, 2)
     & / -0.1186928E+01, -0.3282460E+00,  0.1052020E+00 /
      DATA (AM( 7,K, 1),K=0, 2)
     & / -0.8545702E+01, -0.6247947E+01,  0.3692561E+01 /
      DATA (AM( 8,K, 1),K=0, 2)
     & /  0.1724598E-01,  0.7120465E-02,  0.4003646E-04 /
 
      DATA MEXVEC( 0) / 8 /
      DATA MLFVEC( 0) / 2 /
      DATA UT1VEC( 0) / -0.4656819E+00 /
      DATA UT2VEC( 0) / -0.2742390E+03 /
      DATA ALFVEC( 0) /  0.4491863E+00 /
      DATA QMAVEC( 0) /  0.0000000E+00 /
      DATA (AM( 0,K, 0),K=0, 2)
     & /  0.1193572E+03, -0.3886845E+01, -0.1133965E+01 /
      DATA (AM( 1,K, 0),K=0, 2)
     & / -0.9421449E+02,  0.3995885E+01,  0.1607363E+01 /
      DATA (AM( 2,K, 0),K=0, 2)
     & /  0.4206383E+01,  0.2485954E+00,  0.2497468E+00 /
      DATA (AM( 3,K, 0),K=0, 2)
     & /  0.1210557E+03, -0.3015765E+01, -0.1423651E+01 /
      DATA (AM( 4,K, 0),K=0, 2)
     & / -0.1013897E+03, -0.7113478E+00,  0.2621865E+00 /
      DATA (AM( 5,K, 0),K=0, 2)
     & / -0.1312404E+01, -0.9297691E+00, -0.1562531E+00 /
      DATA (AM( 6,K, 0),K=0, 2)
     & /  0.1627137E+01,  0.4954111E+00, -0.6387009E+00 /
      DATA (AM( 7,K, 0),K=0, 2)
     & /  0.1537698E+00, -0.2487878E+00,  0.8305947E+00 /
      DATA (AM( 8,K, 0),K=0, 2)
     & /  0.2496448E-01,  0.2457823E-02,  0.8234276E-03 /
 
      DATA MEXVEC(-1) / 8 /
      DATA MLFVEC(-1) / 2 /
      DATA UT1VEC(-1) /  0.3862583E+01 /
      DATA UT2VEC(-1) / -0.1265969E+01 /
      DATA ALFVEC(-1) /  0.2457668E+00 /
      DATA QMAVEC(-1) /  0.0000000E+00 /
      DATA (AM( 0,K,-1),K=0, 2)
     & /  0.2647441E+02,  0.1059277E+02, -0.9176654E+00 /
      DATA (AM( 1,K,-1),K=0, 2)
     & /  0.1990636E+01,  0.8558918E-01,  0.4248667E-01 /
      DATA (AM( 2,K,-1),K=0, 2)
     & / -0.1476095E+02, -0.3276255E+02,  0.1558110E+01 /
      DATA (AM( 3,K,-1),K=0, 2)
     & / -0.2966889E+01, -0.3649037E+02,  0.1195914E+01 /
      DATA (AM( 4,K,-1),K=0, 2)
     & / -0.1000519E+03, -0.2464635E+01,  0.1964849E+00 /
      DATA (AM( 5,K,-1),K=0, 2)
     & /  0.3718331E+02,  0.4700389E+02, -0.2772142E+01 /
      DATA (AM( 6,K,-1),K=0, 2)
     & / -0.1872722E+02, -0.2291189E+02,  0.1089052E+01 /
      DATA (AM( 7,K,-1),K=0, 2)
     & / -0.1628146E+02, -0.1823993E+02,  0.2537369E+01 /
      DATA (AM( 8,K,-1),K=0, 2)
     & / -0.1156300E+01, -0.1280495E+00,  0.5153245E-01 /
 
      DATA MEXVEC(-2) / 7 /
      DATA MLFVEC(-2) / 2 /
      DATA UT1VEC(-2) /  0.1895615E+00 /
      DATA UT2VEC(-2) / -0.3069097E+01 /
      DATA ALFVEC(-2) /  0.5293999E+00 /
      DATA QMAVEC(-2) /  0.0000000E+00 /
      DATA (AM( 0,K,-2),K=0, 2)
     & / -0.6556775E+00,  0.2490190E+00,  0.3966485E-01 /
      DATA (AM( 1,K,-2),K=0, 2)
     & /  0.1305102E+01, -0.1188925E+00, -0.4600870E-02 /
      DATA (AM( 2,K,-2),K=0, 2)
     & / -0.2371436E+01,  0.3566814E+00, -0.2834683E+00 /
      DATA (AM( 3,K,-2),K=0, 2)
     & / -0.6152826E+01,  0.8339877E+00, -0.7233230E+00 /
      DATA (AM( 4,K,-2),K=0, 2)
     & / -0.8346558E+01,  0.2892168E+01,  0.2137099E+00 /
      DATA (AM( 5,K,-2),K=0, 2)
     & /  0.1279530E+02,  0.1021114E+00,  0.5787439E+00 /
      DATA (AM( 6,K,-2),K=0, 2)
     & /  0.5858816E+00, -0.1940375E+01, -0.4029269E+00 /
      DATA (AM( 7,K,-2),K=0, 2)
     & / -0.2795725E+02, -0.5263392E+00,  0.1290229E+01 /
 
      DATA MEXVEC(-3) / 7 /
      DATA MLFVEC(-3) / 2 /
      DATA UT1VEC(-3) /  0.3753257E+01 /
      DATA UT2VEC(-3) / -0.1113085E+01 /
      DATA ALFVEC(-3) /  0.3713141E+00 /
      DATA QMAVEC(-3) /  0.0000000E+00 /
      DATA (AM( 0,K,-3),K=0, 2)
     & /  0.1580931E+01, -0.2273826E+01, -0.1822245E+01 /
      DATA (AM( 1,K,-3),K=0, 2)
     & /  0.2702644E+01,  0.6763243E+00,  0.7231586E-02 /
      DATA (AM( 2,K,-3),K=0, 2)
     & / -0.1857924E+02,  0.3907500E+01,  0.5850109E+01 /
      DATA (AM( 3,K,-3),K=0, 2)
     & / -0.3044793E+02,  0.2639332E+01,  0.5566644E+01 /
      DATA (AM( 4,K,-3),K=0, 2)
     & / -0.4258011E+01, -0.5429244E+01,  0.4418946E+00 /
      DATA (AM( 5,K,-3),K=0, 2)
     & /  0.3465259E+02, -0.5532604E+01, -0.4904153E+01 /
      DATA (AM( 6,K,-3),K=0, 2)
     & / -0.1658858E+02,  0.2923275E+01,  0.2266286E+01 /
      DATA (AM( 7,K,-3),K=0, 2)
     & / -0.1149263E+02,  0.2877475E+01, -0.7999105E+00 /
 
      DATA MEXVEC(-4) / 7 /
      DATA MLFVEC(-4) / 2 /
      DATA UT1VEC(-4) /  0.4400772E+01 /
      DATA UT2VEC(-4) / -0.1356116E+01 /
      DATA ALFVEC(-4) /  0.3712017E-01 /
      DATA QMAVEC(-4) /  0.1300000E+01 /
      DATA (AM( 0,K,-4),K=0, 2)
     & / -0.8293661E+00, -0.3982375E+01, -0.6494283E-01 /
      DATA (AM( 1,K,-4),K=0, 2)
     & /  0.2754618E+01,  0.8338636E+00, -0.6885160E-01 /
      DATA (AM( 2,K,-4),K=0, 2)
     & / -0.1657987E+02,  0.1439143E+02, -0.6887240E+00 /
      DATA (AM( 3,K,-4),K=0, 2)
     & / -0.2800703E+02,  0.1535966E+02, -0.7377693E+00 /
      DATA (AM( 4,K,-4),K=0, 2)
     & / -0.6460216E+01, -0.4783019E+01,  0.4913297E+00 /
      DATA (AM( 5,K,-4),K=0, 2)
     & /  0.3141830E+02, -0.3178031E+02,  0.7136013E+01 /
      DATA (AM( 6,K,-4),K=0, 2)
     & / -0.1802509E+02,  0.1862163E+02, -0.4632843E+01 /
      DATA (AM( 7,K,-4),K=0, 2)
     & / -0.1240412E+02,  0.2565386E+02, -0.1066570E+02 /
 
      DATA MEXVEC(-5) / 6 /
      DATA MLFVEC(-5) / 2 /
      DATA UT1VEC(-5) /  0.5562568E+01 /
      DATA UT2VEC(-5) / -0.1801317E+01 /
      DATA ALFVEC(-5) /  0.4952010E-02 /
      DATA QMAVEC(-5) /  0.4500000E+01 /
      DATA (AM( 0,K,-5),K=0, 2)
     & / -0.6031237E+01,  0.1992727E+01, -0.1076331E+01 /
      DATA (AM( 1,K,-5),K=0, 2)
     & /  0.2933912E+01,  0.5839674E+00,  0.7509435E-01 /
      DATA (AM( 2,K,-5),K=0, 2)
     & / -0.8284919E+01,  0.1488593E+01, -0.8251678E+00 /
      DATA (AM( 3,K,-5),K=0, 2)
     & / -0.1925986E+02,  0.2805753E+01, -0.3015446E+01 /
      DATA (AM( 4,K,-5),K=0, 2)
     & / -0.9480483E+01, -0.9767837E+00, -0.1165544E+01 /
      DATA (AM( 5,K,-5),K=0, 2)
     & /  0.2193195E+02, -0.1788518E+02,  0.9460908E+01 /
      DATA (AM( 6,K,-5),K=0, 2)
     & / -0.1327377E+02,  0.1201754E+02, -0.6277844E+01 /
 
      IF(Q .LE. QMAVEC(IFL)) THEN
         PYCT5L = 0.D0
         RETURN
      ENDIF
 
      IF(X .GE. 1.D0) THEN
         PYCT5L = 0.D0
         RETURN
      ENDIF
 
      TMP = LOG(Q/ALFVEC(IFL))
      IF(TMP .LE. 0.D0) THEN
         PYCT5L = 0.D0
         RETURN
      ENDIF
 
      SB = LOG(TMP)
      SB1 = SB - 1.2D0
      SB2 = SB1*SB1
 
      DO 110 I = 0, NEX
         AF(I) = 0.D0
         SBX = 1.D0
         DO 100 K = 0, MLFVEC(IFL)
            AF(I) = AF(I) + SBX*AM(I,K,IFL)
            SBX = SB1*SBX
  100    CONTINUE
  110 CONTINUE
 
      Y = -LOG(X)
      U = LOG(X/0.00001D0)
 
      PART1 = AF(1)*Y**(1.D0+0.01D0*AF(4))*(1.D0+ AF(8)*U)
      PART2 = AF(0)*(1.D0 - X) + AF(3)*X
      PART3 = X*(1.D0-X)*(AF(5)+AF(6)*(1.D0-X)+AF(7)*X*(1.D0-X))
      PART4 = UT1VEC(IFL)*LOG(1.D0-X) +
     &	      AF(2)*LOG(1.D0+EXP(UT2VEC(IFL))-X)
 
      PYCT5L = EXP(LOG(X) + PART1 + PART2 + PART3 + PART4)
 
C...Include threshold factor.
      PYCT5L = PYCT5L * (1.D0 - QMAVEC(IFL)/Q)
 
      RETURN
      END
