//____________________________________________________________________________
/*!

\class    genie::SchmitzMultiplicityModel

\brief    The 'Schmitz' multiplicity probability model as used in NeuGEN.
          Is a concerete implementation of the MultiplicityProbModelI interface.

\ref      N. Schmitz, Proc. Intl. Symp. on Lepton & Photon Interactions at
          High Energies, Bonn 1981 p.527
          The probability scaling factors for low multiplicity states (m=2,3)
          have been taken from NeuGEN (H.Gallagher et al.)

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  August 21, 2004

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _SCHMITZ_MULTIPLICITY_MODEL_H_
#define _SCHMITZ_MULTIPLICITY_MODEL_H_

#include <TH1D.h>

#include "Fragmentation/MultiplicityProbModelI.h"

namespace genie {

class Spline;

class SchmitzMultiplicityModel : public MultiplicityProbModelI {

public:
  SchmitzMultiplicityModel();
  SchmitzMultiplicityModel(string config);
  virtual ~SchmitzMultiplicityModel();

  TH1D * ProbabilityDistribution(const Interaction * interaction) const;

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig   (void);
  double SelectOffset (const Interaction * i) const;
  void   ApplyRijk    (const Interaction * i, TH1D * p, bool norm=true) const;

  //! configuration data

  // KNO distribution
  Spline * fKNO;

  // parameters controling the average multiplicity for given W
  double fAvp;
  double fAvn;
  double fAvbp;
  double fAvbn;
  double fB;

  // NEUGEN's Rijk parameters for controling reductions in the
  // probabilities for low multiplicity states (m=2,3) needed
  // to avoid double-counting with resonance channels
  double fRvpCCm1;
  double fRvpCCm2;
  double fRvpNCm1;
  double fRvpNCm2;
  double fRvnCCm1;
  double fRvnCCm2;
  double fRvnNCm1;
  double fRvnNCm2;
  double fRvbpCCm1;
  double fRvbpCCm2;
  double fRvbpNCm1;
  double fRvbpNCm2;
  double fRvbnCCm1;
  double fRvbnCCm2;
  double fRvbnNCm1;
  double fRvbnNCm2;

  // If true forces a hard upper limit in hadronic multiplicity equal
  // to 10, to be consistent with NEUGEN
  bool fForceNeuGenLimit;
};

}         // genie namespace

#endif    // _SCHMITZ_MULTIPLICITY_MODEL_H_

