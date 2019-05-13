//____________________________________________________________________________
/*!

\namespace genie::utils::bwfunc

\brief     Breit Wigner functions

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           University of Liverpool & STFC Rutherford Appleton Lab

\created   November 22, 2004

\cpright   Copyright (c) 2003-2019, The GENIE Collaboration
           For the full text of the license visit http://copyright.genie-mc.org
           or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _BREIT_WIGNER_UTILS_H_
#define _BREIT_WIGNER_UTILS_H_

namespace genie  {
namespace utils  {
namespace bwfunc {
  //-- A realistic Breit-Wigner distribution with L-dependent and Wlimit
  double BreitWignerLGamma(
             double W, int L, double mass, double width0, double norm);
  


  //-- A realistic Breit-Wigner distribution with L-dependent width.
  double BreitWignerL(
             double W, int L, double mass, double width0, double norm);

  //-- A simple Breit-Wigner distribution.
  double BreitWigner(double W, double mass, double width, double norm);

} // bwfunc namespace
} // utils  namespace
} // genie  namespace

#endif   // _BREIT_WIGNER_UTILS_H_
