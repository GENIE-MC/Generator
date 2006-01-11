//____________________________________________________________________________
/*!

\class    genie::XSecAlgorithmI

\brief    Cross Section Calculation Interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 03, 2004

*/
//____________________________________________________________________________

#ifndef _XSEC_ALGORITHM_I_H_
#define _XSEC_ALGORITHM_I_H_

#include "Algorithm/Algorithm.h"
#include "Interaction/Interaction.h"

namespace genie {

class XSecAlgorithmI : public Algorithm {

public:

  virtual ~XSecAlgorithmI();

  virtual double XSec            (const Interaction * in) const = 0;
  virtual bool   ValidProcess    (const Interaction * in) const = 0;
  virtual bool   ValidKinematics (const Interaction * in) const = 0;

protected:

  XSecAlgorithmI();
  XSecAlgorithmI(string name);
  XSecAlgorithmI(string name, string config);
};

}       // genie namespace

#endif  // _XSEC_ALGORITHM_I_H_
