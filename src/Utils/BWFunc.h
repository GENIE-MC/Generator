//____________________________________________________________________________
/*!

\namespace genie::utils::bwfunc

\brief     Breit Wigner functions

\author    Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
           STFC, Rutherford Appleton Laboratory

\created   November 22, 2004

\cpright   Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
           All rights reserved.
           For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _BREIT_WIGNER_UTILS_H_
#define _BREIT_WIGNER_UTILS_H_

namespace genie  {
namespace utils  {
namespace bwfunc {

  //-- A realistic Breit-Wigner distribution with L-dependent width.
  double BreitWignerL(
             double W, int L, double mass, double width0, double norm);

  //-- A simple Breit-Wigner distribution.
  double BreitWigner(double W, double mass, double width, double norm);

} // bwfunc namespace
} // utils  namespace
} // genie  namespace

#endif   // _BREIT_WIGNER_UTILS_H_
