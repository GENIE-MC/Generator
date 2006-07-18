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

  const TH1D & ProbabilityDistribution(const Interaction * interaction) const;

  //! overload the Algorithm::Configure() methods to load private data
  //! members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  void   LoadConfig     (void);
  void   CreateProbHist (double maxmult) const;
  double AverageChMult  (int nu_pdg, int nuc_pdg, double W) const;
  double KNO            (int nu_pdg, int nuc_pdg, double z) const;
  double SelectOffset   (const Interaction * i) const;
  void   ApplyRijk      (const Interaction * i, bool norm=true) const;

  // computed multiplicity probability distribution
  mutable TH1D * fMultProb;

  //! configuration data

  // KNO distribution
  Spline * fKNO;

  //! flags
  bool fForceNeuGenLimit;   ///< force upper hadronic multiplicity to NeuGEN limit
  bool fApplyRijk;          ///< apply multiplicity probability scaling factors
  bool fRenormalize;        ///< re-normalize after applying scaling factors
  bool fUseLegacyKNOSpline; ///< use legacy spline instead of Levy

  //! parameters controling the average multiplicity for given W
  double fAvp;
  double fAvn;
  double fAvbp;
  double fAvbn;
  double fB;

  //! Levy function parameter c
  double fCvp;
  double fCvn;
  double fCvbp;
  double fCvbn;

  //! under the DIS/RES joining scheme, multiplicity probability scaling 
  //! factors would be applied for W<Wcut
  double fWcut;

  //! NEUGEN's low-multiplicity probability scaling parameters
  double fRvpCCm2;   ///< vp,  CC, multiplicity = 2
  double fRvpCCm3;   ///< vp,  CC, multiplicity = 3
  double fRvpNCm2;   ///< vp,  NC, multiplicity = 2
  double fRvpNCm3;   ///< vp,  NC, multiplicity = 3
  double fRvnCCm2;   ///< vn,  CC, multiplicity = 2
  double fRvnCCm3;   ///< vn,  CC, multiplicity = 3
  double fRvnNCm2;   ///< vn,  NC, multiplicity = 2
  double fRvnNCm3;   ///< vn,  NC, multiplicity = 3
  double fRvbpCCm2;  ///< vbp, CC, multiplicity = 2
  double fRvbpCCm3;  ///< vbp, CC, multiplicity = 3
  double fRvbpNCm2;  ///< vbp, NC, multiplicity = 2
  double fRvbpNCm3;  ///< vbp, NC, multiplicity = 3
  double fRvbnCCm2;  ///< vbn, CC, multiplicity = 2
  double fRvbnCCm3;  ///< vbn, CC, multiplicity = 3
  double fRvbnNCm2;  ///< vbn, NC, multiplicity = 2
  double fRvbnNCm3;  ///< vbn, NC, multiplicity = 3
};

}         // genie namespace
#endif    // _SCHMITZ_MULTIPLICITY_MODEL_H_

