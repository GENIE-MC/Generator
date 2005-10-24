//____________________________________________________________________________
/*!

\class    genie::BardinDISPXSec

\brief    Computes, the differential cross section for neutrino DIS including
          radiative corrections according to the \b Bardin-Dokuchaeva model.

          The computed xsec is the double differential d^2(xsec) / dy dx \n
          where \n
           \li \c y is the inelasticity, and
           \li \c x is the Bjorken scaling variable \c x

          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      D.Yu.Bardin, V.A.Dokuchaeva, "On the radiative corrections to the
          Neutrino Deep Inelastic Scattering", JINR-E2-86-260, Apr. 1986

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  July 06, 2004
 
*/
//____________________________________________________________________________

#ifndef _BARDIN_DIS_PARTIAL_XSEC_H_
#define _BARDIN_DIS_PARTIAL_XSEC_H_

#include "Base/XSecAlgorithmI.h"
#include "PDF/PDF.h"

namespace genie {

class BardinDISPXSec : public XSecAlgorithmI {

public:

  BardinDISPXSec();
  BardinDISPXSec(string config);
  virtual ~BardinDISPXSec();

  //-- XSecAlgorithmI interface implementation
  
  double XSec (const Interaction * interaction) const;

private:

  //-- Auxiliary methods used for the Bardin-Dokuchaeva xsec calculation.
  
  //   xi is the running variable in the integration for computing
  //   the differential cross section d^2(xsec) / dy dx
  
  double PhiCCi   (double xi, const Interaction * interaction) const;
  double Ii       (double xi, const Interaction * interaction) const;
  double U        (double xi, const Interaction * interaction) const;
  double tau      (double xi, const Interaction * interaction) const;
  double St       (double xi, const Interaction * interaction) const;
  double Su       (double xi, const Interaction * interaction) const;
  
  double S        (const Interaction * interaction) const;
  double DeltaCCi (const Interaction * interaction) const;
  double Sq       (const Interaction * interaction) const;

  double InitQuarkMass (const Interaction * interaction) const;
  double InitQuarkMass (double xi) const;
  double PDFFunc       (const PDF & pdf, int pgdc) const;
};

}       // genie namespace

#endif  // _BARDIN_DIS_PARTIAL_XSEC_H_
