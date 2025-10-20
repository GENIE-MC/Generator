//____________________________________________________________________________
/*!

\brief    QE utilities

\author   Steven Gardiner <gardiner \at fnal.gov>
          Fermilab

\created  May 04, 2004

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
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

    // Leave the struck nucleon on shell, effectively ignoring its binding
    // energy
    kOnShell,

    // Use a prescription equivalent to that of the Valencia model (see
    // Eq. (43) in https://arxiv.org/abs/nucl-th/0408005). In this case,
    // the effective energy transfer implies an off-shell hit nucleon
    // total energy in the initial state.
    kValenciaStyleQValue,
  } QELEvGen_BindingMode_t;

  namespace utils {

    double EnergyDeltaFunctionSolutionQEL(const Interaction& inter);

    QELEvGen_BindingMode_t StringToQELBindingMode( const std::string& mode_str );

    double ComputeFullQELPXSec(Interaction* interaction,
      const NuclearModelI* nucl_model, const XSecAlgorithmI* xsec_model,
      double cos_theta_0, double phi_0, double& Eb,
      QELEvGen_BindingMode_t hitNucleonBindingMode, double min_angle_EM = 0.,
      bool bind_nucleon = true);

    double CosTheta0Max(const genie::Interaction& interaction);

    void BindHitNucleon(Interaction& interaction, const NuclearModelI& nucl_model,
      double& Eb, QELEvGen_BindingMode_t hitNucleonBindingMode);
  }
}

#endif
