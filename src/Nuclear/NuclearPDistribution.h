//____________________________________________________________________________
/*!

\class    genie::NuclearPDistribution

\brief    Describes a Nucleon Momentum Probability Distribution.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  October 07, 2004

*/
//____________________________________________________________________________

#ifndef _NUCLEON_MOMENTUM_DISTRIBUTION_H_
#define _NUCLEON_MOMENTUM_DISTRIBUTION_H_

#include <TH1D.h>
#include <TVector3.h>

#include "Algorithm/Algorithm.h"
#include "Interaction/Target.h"
#include "Nuclear/NuclearPDistributionModelI.h"

namespace genie {

class NuclearPDistribution {

public:

  NuclearPDistribution();
  NuclearPDistribution(const NuclearPDistribution & nmd);
  virtual ~NuclearPDistribution();

  void     AttachModel           (const NuclearPDistributionModelI * model);
  void     BuildProbDistribution (const Target & target);
  double   Probability           (double p)           const;
  double   Probability           (const TVector3 & p) const;
  double   RandomMomentum        (void) const;
  TVector3 RandomMomentum3       (void) const;

private:

  void   Init (void);

  TH1D * fProbDistribution;

  const NuclearPDistributionModelI * fProbModel;
};

}         // genie namespace

#endif    // _NUCLEON_MOMENTUM_DISTRIBUTION_H_

