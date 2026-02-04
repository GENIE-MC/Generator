//____________________________________________________________________________
/*!

\class    genie::LeptoHadPythia6

\brief    Provides access to the PYTHIA6 hadronization. \n

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC (Valencia)

\created  December 12, 2024

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _LEPTO_HAD_PYTHIA6__H_
#define _LEPTO_HAD_PYTHIA6__H_

#include "Framework/EventGen/EventRecordVisitorI.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/EventGen/EVGThreadException.h"

#ifdef __GENIE_PYTHIA6_ENABLED__
#include <TPythia6.h>
#endif

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::math;

namespace genie {

class GHepParticle;

class LeptoHadPythia6 : public EventRecordVisitorI {

public:
  LeptoHadPythia6();
  LeptoHadPythia6(string config);
  virtual ~LeptoHadPythia6();

  //-- implement the HadronizationModelI interface
  void ProcessEventRecord(GHepRecord * event) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  bool           Hadronize        (GHepRecord * event) const;
  void           Initialize       (void)               const;
  void           LoadConfig       (void);

  //-- configuration parameters
  int    fMaxIterHad;         // Maxmium number of iterations to look for a combination of hadrons
  double fPrimordialKT;       // Width of Gaussian distribution for the primordial transverse momentum kT of partons in the nucleon.
  double fRemnantPT;          // Width of Gaussian distribution in transverse momentum when a non-trivial target remnant is split into two particles
  double fMinESinglet;        // It is, with quark masses added, used to define the minimum allowable energy of a colour-singlet parton system.
  double fWmin;               // Minimum value of W

  // PYTHIA physics configuration parameters used
  double fSSBarSuppression;       ///< ssbar suppression
  double fGaussianPt2;            ///< gaussian pt2 distribution width
  double fNonGaussianPt2Tail;     ///< non gaussian pt2 tail parameterization
  double fRemainingECutoff;       ///< remaining E cutoff stopping fragmentation
  double fDiQuarkSuppression;     ///< di-quark suppression parameter
  double fLightVMesonSuppression; ///< light vector meson suppression
  double fSVMesonSuppression;     ///< strange vector meson suppression
  double fLunda;                  ///< Lund a parameter
  double fLundb;                  ///< Lund b parameter
  double fLundaDiq;               ///< adjustment of Lund a for di-quark

#ifdef __GENIE_PYTHIA6_ENABLED__
  mutable TPythia6 * fPythia;   ///< PYTHIA6 wrapper class
#endif

};

}         // genie namespace

#endif    // _LEPTO_HAD_PYTHIA6__H_

