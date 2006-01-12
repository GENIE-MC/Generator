//____________________________________________________________________________
/*!

\class    genie::NuclMomentumModelI

\brief    Pure abstract base class.
          Defines the NuclMomentumModelI interface to be implemented by
          any algorithmic class generating momenta for nucleons within nuclei

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 09, 2004

*/
//____________________________________________________________________________

#ifndef _NUCLEAR_MOMENTUM_MODEL_I_H_
#define _NUCLEAR_MOMENTUM_MODEL_I_H_

#include <TH1D.h>

#include "Algorithm/Algorithm.h"
#include "Interaction/Target.h"

namespace genie {

class NuclMomentumModelI : public Algorithm {

public:

  virtual ~NuclMomentumModelI();

  virtual TH1D * ProbabilityDistribution(const Target & target) const = 0;

protected:

  NuclMomentumModelI();
  NuclMomentumModelI(string name);
  NuclMomentumModelI(string name, string config);
};

}         // genie namespace

#endif    // _NUCLEAR_MOMENTUM_DISTRIBUTION_MODEL_I_H_

