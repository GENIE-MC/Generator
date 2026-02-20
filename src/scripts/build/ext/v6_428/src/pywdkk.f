C*********************************************************************

C...PYWDKK
C...Universal Extra Dimensions Model (UED)
C...
C...Multiplied by the square modulus of a form factor
C...(see GRADEN in function PYGRAW)
C...PYWDKK is the KK boson -> SM boson + graviton
C...gravity mediated partial decay width Gamma(xx, yy)
C...  where xx is exclusive to gravity
C...  yy=m_Graviton/m_bosonKK denotes the Universal extra dimension
C...  and xxa=sqrt(xx**2+yy**2) refers to all of the extra dimensions
C...
C...N.B. The Feynman rules for the couplings of the graviton fields
C...to the UED fields are related to the corresponding couplings of
C...the graviton fields to the SM fields by the form factor.

      DOUBLE PRECISION FUNCTION PYWDKK(X)

C...Double precision and integer declarations
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)

C...Pythia commonblocks
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)

C...Local UED commonblocks and variables
      COMMON/UEDGRA/XMPLNK,XMD,RINV,NDIM
      COMMON/KAPPA/XKAPPA

C...SAVE statements
      SAVE /PYDAT1/,/PYDAT2/,/UEDGRA/,/KAPPA/

      PI=PARU(1)

C...gamma* mass 473
      KCQKK=473
      XMNKK=PMAS(KCQKK,1)

C...Bosons partial width Macesanu hep-ph/0201300
      PYWDKK=XKAPPA**2/(96.*PI)*XMNKK**3/X**4*
     +          ((1.-X**2)**2*(1.+3.*X**2+6.*X**4))

      RETURN
      END
