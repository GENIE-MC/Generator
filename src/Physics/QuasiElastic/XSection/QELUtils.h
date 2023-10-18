//____________________________________________________________________________
/*!

\brief    QE utilities

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermilab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
*/
//____________________________________________________________________________

#ifndef _QEL_UTILS_H_
#define _QEL_UTILS_H_

#include "Framework/Interaction/Interaction.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/EventGen/XSecAlgorithmI.h"

#include "TLorentzVector.h"
#include "TVector3.h"

#include <string>

namespace genie {

  // Enumerated type used to specify the method for determining the off-shell energy
  // of the hit nucleon for quasielastic events
  typedef enum EQELEvGenBindingMode {

    // Use removal energy from the nuclear model
    kUseNuclearModel,

    // Calculate binding energy assuming that the remnant nucleus is left in its
    // ground state
    kUseGroundStateRemnant,

    // Leave the struck nucleon on shell, effectively ignoring its binding energy
    kOnShell
  } QELEvGen_BindingMode_t;

  namespace utils {

    double EnergyDeltaFunctionSolutionQEL(const Interaction& inter);

    QELEvGen_BindingMode_t StringToQELBindingMode( const std::string& mode_str );

    double ComputeFullQELPXSec(Interaction* interaction,
      const NuclearModelI* nucl_model, const XSecAlgorithmI* xsec_model,
      double cos_theta_0, double phi_0, double& Eb,
      QELEvGen_BindingMode_t hitNucleonBindingMode, double fMinAngleEM,bool bind_nucleon = true);

    double CosTheta0Max(const genie::Interaction& interaction);

    void BindHitNucleon(Interaction& interaction, const NuclearModelI& nucl_model,
      double& Eb, QELEvGen_BindingMode_t hitNucleonBindingMode);

    double ComputeTestQELPXSec(Interaction* interaction,
      const XSecAlgorithmI* xsec_model, double omega, double ctl,
      double phiLep, double phiNi, double pNi, double Eremove, KinePhaseSpace_t kps = genie::kPSQELEvGen);
  
    void Rotate_qvec_alongZ(TLorentzVector &probe_leptonP4, 
      TLorentzVector &out_leptonP4, std::vector<TLorentzVector> &otherP4);
  }
}

#endif
