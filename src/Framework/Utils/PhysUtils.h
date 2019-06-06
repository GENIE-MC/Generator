//____________________________________________________________________________
/*!

\namespace  genie::utils::phys

\brief      Various physics formulas & utilities

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    January 22, 2008

\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PHYS_UTILS_H_
#define _PHYS_UTILS_H_

#include <TLorentzVector.h>

namespace genie {
namespace utils {

namespace phys
{
  // Formation zone in fm
  double FormationZone(
     double m, const TLorentzVector & p, const TVector3 & p3hadr, double ct0 /*in fm*/, double K);

  // Longitudinal to transverse cross section ratio (R) parametrizations
  double R99118   (double x, double Q2); ///< PRL 98, 142301, 2007
  double RWhitlow (double x, double Q2);

  // Extract F1, F2, xF3 from a three d^2sigma/dxdy cross section values
  // evaluated at different (E,y) for fixed (x,Q2)
  // See H.Gallagher, Nucl.Phys.Proc.Suppl.159:229-234,2006
/*
  void ExtractStructFunc (
    double x, double Q2, double dxs[3], double& F1, double& F2, double& xF3);
*/

} // phys  namespace
} // utils namespace
} // genie namespace

#endif // _PHYS_UTILS_H_
