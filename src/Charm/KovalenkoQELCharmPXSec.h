//____________________________________________________________________________
/*!

\class    genie::KovalenkoQELCharmPXSec
         
\brief    Computes the QEL Charm Production Differential Cross Section
          using \b Kovalenko's duality model approach. 

          The computed differential cross section is the Dxsec = dxsec/dQ^2
          where \n
            \li \c Q2 is the momentum transfer.
                    
          It models the differential cross sections for: \n
             \li v + n \rightarrow mu- + Lambda_{c}^{+} (2285)
             \li v + n \rightarrow mu- + Sigma_{c}^{+}  (2455)
             \li v + p \rightarrow mu- + Sigma_{c}^{++} (2455)

          Is a concrete implementation of the XSecAlgorithmI interface. 
           
\ref      S.G.Kovalenko, Sov.J.Nucl.Phys.52:934 (1990)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 10, 2004

*/
//____________________________________________________________________________

#ifndef _KOVALENKO_QEL_CHARM_PARTIAL_XSEC_H_
#define _KOVALENKO_QEL_CHARM_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "Numerical/IntegratorI.h"

namespace genie {

class KovalenkoQELCharmPXSec : public XSecAlgorithmI {

public:

  KovalenkoQELCharmPXSec();
  KovalenkoQELCharmPXSec(string config);
  virtual ~KovalenkoQELCharmPXSec();

  //-- XSecAlgorithmI interface implementation
  
  double XSec (const Interaction * interaction) const;

private:

  void AssertProcessValidity(const Interaction * interaction) const;
  void Q2Cuts(double & Q2min, double & Q2max) const;

  const IntegratorI * Integrator (void) const;
  
  double ZR       (const Interaction * interaction) const;
  double DR       (const Interaction * interaction, bool norm = false) const;                              
  double MRes     (const Interaction * interaction) const;
  double MRes2    (const Interaction * interaction) const;
  double vR_minus (const Interaction * interaction) const;
  double vR_plus  (const Interaction * interaction) const;
  double SumF2    (const Interaction * interaction) const;
  double ResDM    (const Interaction * interaction) const;
  double xiBar    (const Interaction * interaction, double v) const;  
  double MoScale  (void) const;
};

}       // genie namespace

#endif  // _KOVALENKO_QEL_CHARM_PARTIAL_XSEC_H_
