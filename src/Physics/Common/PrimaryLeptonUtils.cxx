//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TMatrixD.h"
 #include "TDecompSVD.h"

#include "Framework/Conventions/Constants.h"
#include "Physics/Common/PrimaryLeptonUtils.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/ParticleData/PDGLibrary.h"

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;

//___________________________________________________________________________
void genie::utils::SetPrimaryLeptonPolarization( GHepRecord * ev )
{
// Moved out of the PrimaryLeptonGenerator class to make the same treatment
// accessible for generators that use a more unified approach (e.g.,
// QELEventGenerator and MECGenerator). -- S. Gardiner

  // get the final state primary lepton
  GHepParticle * fsl = ev->FinalStatePrimaryLepton();
  if ( !fsl ) {
    LOG("LeptonicVertex", pERROR)
      << "Final state lepton not set yet! \n" << *ev;
    return;
  }
  //-- Get the interaction
  Interaction * interaction = ev->Summary();
  //-- Access cross section algorithm for running thread
  RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
  const EventGeneratorI * evg = rtinfo->RunningThread();
  const XSecAlgorithmI * xsec_alg = evg->CrossSectionAlg();
  fsl->SetPolarization(xsec_alg->FinalLeptonPolarization(interaction));
    
  LOG("LeptonicVertex", pINFO)
    << "Setting polarization angles for particle: " << fsl->Name();

  if ( fsl->PolzIsSet() ) {
    LOG("LeptonicVertex", pINFO)
      << "Polarization (rad): Polar = "  << fsl->PolzPolarAngle()
      << ", Azimuthal = " << fsl->PolzAzimuthAngle();
  }
  
}
