//____________________________________________________________________________
/*!

\class    genie::BodekRitchie

\brief    The Bodek-Ritchie model for the probability distribution of nucleon
          momenta within a nucleus.

          Implements the NuclearPDistributionModelI interface.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 09, 2004

*/
//____________________________________________________________________________

#ifndef _BODEK_RITCHIE_H_
#define _BODEK_RITCHIE_H_

#include "Nuclear/NuclearPDistributionModelI.h"

namespace genie {

class BodekRitchie : public NuclearPDistributionModelI {

public:

  BodekRitchie();
  BodekRitchie(const char * param_set);
  virtual ~BodekRitchie();

  TH1D * ProbabilityDistribution(const Target & target) const;
};

}         // genie namespace

#endif    // _BODEK_RITCHIE_H_

