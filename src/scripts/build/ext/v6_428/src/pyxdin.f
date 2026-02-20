C*********************************************************************
 
C...PYXDIN
C...Universal Extra Dimensions Model (UED)
C...Initialize the xd masses and widths
C...M. ELKACIMI 4/03/2006
C...Modified for inclusion in Pythia Apr 2008, H. Przysiezniak, P. Skands

      SUBROUTINE PYXDIN

C...Double precision and integer declarations.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      IMPLICIT INTEGER(I-N)
      INTEGER PYK,PYCHGE,PYCOMP
C...Commonblocks.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
C...UED Pythia common
      COMMON/PYPUED/IUED(0:99),RUED(0:99)

C...SAVE statements
      SAVE /PYDAT1/,/PYDAT3/,/PYSUBS/,/PYPUED/

C...Print out some info about the UED model
      WRITE(MSTU(11),7000) 
     &    ' ',
     &    '********** PYXDIN: initialization of UED ******************',
     &    ' ',
     &    'Universal Extra Dimensions (UED) switched on ',
     &    ' ',
     &    'This implementation is courtesy of',
     &    '       M.Elkacimi, D.Goujdami, H.Przysiezniak,  ', 
     &    '       see [hep-ph/0602198] (Les Houches 2005) ',
     &    ' ',
     &    'The model follows [hep-ph/0012100] (Appelquist, Cheng,   ',
     &    'Dobrescu), with gravity-mediated decay widths calculated in',
     &    '[hep-ph/0001335] (DeRujula, Donini, Gavela, Rigolin) and ',
     &    'radiative corrections to the KK masses from [hep/ph0204342]',
     &    '(Cheng, Matchev, Schmaltz).'
      WRITE(MSTU(11),7000) 
     &    ' ',
     &    'SM particles can propagate into one small extra dimension  ',
     &    'of size 1/R = RUED(1) GeV. For gravity-mediated decays, the',
     &    'graviton is further allowed to propagate into N = IUED(4)', 
     &    'large (eV^-1) extra dimensions.'
      WRITE(MSTU(11),7000) 
     &    ' ',
     &    'The switches and parameters for UED are:',
     &    '    IUED(1): (D=0) main UED ON(=1)/OFF(=0) switch ',
     &    '    IUED(2): (D=0) Grav. med. decays are set ON(=1)/OFF(=0)',
     &    '    IUED(3): (D=5) number of quark flavours',
     &    '    IUED(4): (D=6) number of large extra dimensions into',
     &    '                   which the graviton propagates',
     &    '    IUED(5): (D=0) Lambda (=0) or Lambda*R (=1) is used',
     &    '    IUED(6): (D=1) With/without rad.corrs. (=1/0)',
     &    '                                                 ',
     &    '    RUED(1): (D=1000.) curvature 1/R of the UED (in GeV)',
     &    '    RUED(2): (D=5000.) gravity mediated (GM) scale (in GeV)',
     &    '    RUED(3): (D=20000.) Lambda cutoff scale (in GeV). Used',
     &    '                        when IUED(5)=0',
     &    '    RUED(4): (D=20.) Lambda*R. Used when IUED(5)=1'
      WRITE(MSTU(11),7000) 
     &    ' ',
     &    'N.B.: the Higgs mass is also a free parameter of the UED ',
     &    'model, but is set through pmas(25,1).',
     &    ' '

C...Hardcoded switch, required by current implementation     
      CALL PYGIVE('MSTP(42)=0')

C...Turn the gravity mediated decay (for the KK pphoton) ON or OFF
      IF(IUED(2).EQ.0) CALL PYGIVE('MDCY(C5100022,1)=0')

C...Calculated the radiative corrections to the KK particle masses
      CALL PYUEDC

C...Initialize the graviton mass
C...only if the KK particles decays gravitationally
      IF(IUED(2).EQ.1) CALL PYGRAM(0)

      WRITE(MSTU(11),7000) 
     &    '********** PYXDIN: UED initialization completed  ***********'

C...Format to use for comments
 7000 FORMAT(' * ',A)

      RETURN
      END
