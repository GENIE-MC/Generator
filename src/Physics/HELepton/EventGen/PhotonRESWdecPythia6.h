//____________________________________________________________________________
/*!

\class    genie::PhotonRESWdecPythia6

\brief    PhotonRES W decay with pythia6.

\author   Alfonso Garcia <aagarciasoto \at km3net.de>
          IFIC (Valencia)

\created  Dec 12, 2024

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PHOTON_RES_WDEC_PYTHIA6_H_
#define _PHOTON_RES_WDEC_PYTHIA6_H_

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
#include "Framework/Utils/StringUtils.h"
#include "Physics/HELepton/XSection/Born.h"

#ifdef __GENIE_PYTHIA6_ENABLED__
#include <TPythia6.h>
#endif

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::math;

namespace genie {

class PhotonRESWdecPythia6 : public EventRecordVisitorI {

public :
  PhotonRESWdecPythia6();
  PhotonRESWdecPythia6(string config);
 ~PhotonRESWdecPythia6();

  //-- implement the HadronizationModelI interface
  void ProcessEventRecord(GHepRecord * event) const;

  //-- overload the Algorithm::Configure() methods to load private data
  //   members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:

  bool           Wdecay           (GHepRecord * event) const;

  void           Initialize       (void)               const;
  void           LoadConfig       (void);

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

  double fQ2PDFmin;

  Born * born;

#ifdef __GENIE_PYTHIA6_ENABLED__
  mutable TPythia6 * fPythia;   ///< PYTHIA6 wrapper class
#endif

};

}      // genie namespace
#endif // _PHOTON_RES_WDEC_PYTHIA6_H_
