//____________________________________________________________________________
/*!

\class    genie::SchmitzMultiplicityModel

\brief    The 'Schmitz' multiplicity probability model as used in NeuGEN.

          Is a concerete implementation of the MultiplicityProbModelI interface.

\ref      N. Schmitz, Proc. Intl. Symp. on Lepton & Photon Interactions at
          High Energies, Bonn 1981 p.527

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 21, 2004

*/
//____________________________________________________________________________

#ifndef _SCHMITZ_MULTIPLICITY_MODEL_H_
#define _SCHMITZ_MULTIPLICITY_MODEL_H_

#include "Fragmentation/MultiplicityProbModelI.h"

namespace genie {

class SchmitzMultiplicityModel : public MultiplicityProbModelI {

public:

  SchmitzMultiplicityModel();
  SchmitzMultiplicityModel(const char * param_set);
  virtual ~SchmitzMultiplicityModel();

  TH1D * ProbabilityDistribution(const Interaction * interaction) const;

private:

  double SelectOffset(const Interaction * interaction) const;
  
};

}         // genie namespace

#endif    // _SCHMITZ_MULTIPLICITY_MODEL_H_

