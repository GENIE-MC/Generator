//____________________________________________________________________________
/*!

\class    genie::AivazisCharmPXSecLO

\brief    Computes, at Leading Order (LO), the differential cross section for
          neutrino charm production using the \b Aivazis,Olness,Tung model.

          The computed cross section is the D2xsec = d^2(xsec) / dy dx \n
          where \n
            \li \c y is the inelasticity, and
            \li \c x is the Bjorken scaling variable \c x

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      M.A.G.Aivazis, F.I.Olness and W.K.Tung
          Phys.Rev.D50, 3085-3101 (1994)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  June 10, 2004

*/
//____________________________________________________________________________

#ifndef _AIVAZIS_CHARM_PARTIAL_XSEC_LO_H_
#define _AIVAZIS_CHARM_PARTIAL_XSEC_LO_H_

#include "Base/XSecAlgorithmI.h"

namespace genie {

class PDFModelI;

class AivazisCharmPXSecLO : public XSecAlgorithmI {

public:

  AivazisCharmPXSecLO();
  AivazisCharmPXSecLO(const char * param_set);
  virtual ~AivazisCharmPXSecLO();

  //-- XSecAlgorithmI interface implementation
  
  double XSec (const Interaction * interaction) const;

private:

  const PDFModelI * PdfModel(void) const;

  //-- methods that allow constant values of charm mass and CKM elements Vcd,
  //   Vcs to be overriden by XML/registry config values if the user wishes to
  double CharmMass (void) const;
  double Vcd       (void) const;
  double Vcs       (void) const;
};

}       // genie namespace

#endif  // _AIVAZIS_CHARM_PARTIAL_XSEC_LO_H_
