//____________________________________________________________________________
/*!

\class    genie::ScatteringParams

\brief    A Registry subclass to hold scattering parameters.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 08, 2004

*/
//____________________________________________________________________________

#ifndef _SCATTERING_PARAMS_H_
#define _SCATTERING_PARAMS_H_

#include <TMath.h>

#include "Registry/Registry.h"

namespace genie {

class ScatteringParams : public Registry {

public:

  ScatteringParams();
  ScatteringParams(const Registry & params);
  virtual ~ScatteringParams();

  double x    (void) const;
  double y    (void) const;
  double Q2   (void) const;
  double q2   (void) const;
  double W    (void) const;
  double lnQ2 (void) const;
  double lnW  (void) const;

};

}       // genie namespace

#endif  // _SCATTERING_PARAMS_H_
