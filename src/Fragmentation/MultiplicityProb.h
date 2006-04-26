//____________________________________________________________________________
/*!

\class    genie::MultiplicityProb

\brief    Describes a Multiplicity Probability Distribution.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 21, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _MULTIPLICITY_PROBABILITY_DISTRIBUTION_H_
#define _MULTIPLICITY_PROBABILITY_DISTRIBUTION_H_

#include <TH1D.h>

#include "Algorithm/Algorithm.h"
#include "Fragmentation/MultiplicityProbModelI.h"

namespace genie {

class MultiplicityProb {

public:

  MultiplicityProb();
  MultiplicityProb(const MultiplicityProb & mpd);
  virtual ~MultiplicityProb();

  void         AttachModel           (const MultiplicityProbModelI * model);
  void         BuildProbDistribution (const Interaction * interaction);
  double       Probability           (int n) const;
  unsigned int RandomMultiplicity    (unsigned int min = 2,
                                      unsigned int max = 20) const;

private:

  void   Init (void);

  TH1D *                         fProbDistribution;
  const MultiplicityProbModelI * fMultProbModel;

};

}         // genie namespace

#endif    // _MULTIPLICITY_PROBABILITY_DISTRIBUTION_H_

