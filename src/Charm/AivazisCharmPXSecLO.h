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
  AivazisCharmPXSecLO(string config);
  virtual ~AivazisCharmPXSecLO();

  //-- XSecAlgorithmI interface implementation
  double XSec            (const Interaction * interaction) const;
  bool   ValidProcess    (const Interaction * interaction) const;
  bool   ValidKinematics (const Interaction * interaction) const;

  //-- override the Algorithm::Configure methods to load configuration
  //   data to private data members
  void Configure (const Registry & config);
  void Configure (string param_set);

private:

  void LoadConfigData (void);
  void LoadSubAlg     (void);

  const PDFModelI* fPDFModel;

  bool   fDContributes;
  bool   fSContributes;
  double fMc;
  double fVcd;
  double fVcs;
  double fMc2;
  double fVcd2;
  double fVcs2;
};

}       // genie namespace

#endif  // _AIVAZIS_CHARM_PARTIAL_XSEC_LO_H_
