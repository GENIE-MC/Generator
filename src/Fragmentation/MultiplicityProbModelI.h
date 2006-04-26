//____________________________________________________________________________
/*!

\class    genie::MultiplicityProbModelI

\brief    Pure abstract base class.
          Defines the MultiplicityProbModelI interface to be implemented by
          any algorithmic class computing multiplicity probability
          distributions for a hadronization model.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 21, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _MULTIPLICITY_PROBABILITY_MODEL_I_H_
#define _MULTIPLICITY_PROBABILITY_MODEL_I_H_

#include <TH1D.h>

#include "Algorithm/Algorithm.h"
#include "Interaction/Interaction.h"

namespace genie {

class MultiplicityProbModelI : public Algorithm {

public:

  virtual ~MultiplicityProbModelI();

  virtual TH1D * ProbabilityDistribution(const Interaction * inter) const = 0;

protected:

  MultiplicityProbModelI();
  MultiplicityProbModelI(string name);
  MultiplicityProbModelI(string name, string config);
};

}         // genie namespace

#endif    // _MULTIPLICITY_PROBABILITY_MODEL_I_H_

