//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Steven Gardiner <gardiner \at fnal.gov>
 Fermi National Accelerator Laboratory
*/
//____________________________________________________________________________

#include "TVector3.h"

#include "Physics/Common/PrimaryLeptonUtils.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/EventGen/EventGeneratorI.h"
#include "Framework/EventGen/RunningThreadInfo.h"
#include "Framework/EventGen/XSecAlgorithmI.h"

using namespace genie;
using namespace genie::utils;

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
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  // Polarization of final charged lepton according to Kuzmin, Lyubushkin and Naumov, hep-ph/0312107
  if ( proc_info.IsWeakCC() )
  {
     //-- Access cross section algorithm for running thread
     RunningThreadInfo * rtinfo = RunningThreadInfo::Instance();
     const EventGeneratorI * evg = rtinfo->RunningThread();
     const XSecAlgorithmI * xsec_alg = evg->CrossSectionAlg();
     interaction->SetBit(kISkipProcessChk);
     Kinematics * kine = interaction->KinePtr();
     kine->UseSelectedKinematics();
     if ( proc_info.IsQuasiElastic() )
     {
        std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
        double xsec = xsec_alg->XSec(interaction, kPSyfEx);
        std::cout << "xsec = " << xsec << "\n";
        xsec = xsec_alg->XSec(interaction, kPSQELEvGen);
        std::cout << "xsec = " << xsec << "\n";
        std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
     }
     return;
  }
  
  // Get (px,py,pz) @ LAB
  TVector3 plab( fsl->Px(), fsl->Py(), fsl->Pz() );

  // In the limit m/E->0: leptons are left-handed and their anti-particles
  // are right-handed
  int pdgc = fsl->Pdg();
  if ( pdg::IsNeutrino(pdgc) || pdg::IsElectron(pdgc) || 
       pdg::IsMuon(pdgc)     || pdg::IsTau(pdgc) )
  {
    plab *= -1; // left-handed
  }

  LOG("LeptonicVertex", pINFO)
    << "Setting polarization angles for particle: " << fsl->Name();

  fsl->SetPolarization( plab );

  if ( fsl->PolzIsSet() ) {
    LOG("LeptonicVertex", pINFO)
      << "Polarization (rad): Polar = "  << fsl->PolzPolarAngle()
      << ", Azimuthal = " << fsl->PolzAzimuthAngle();
  }
  
  // reset trust bits
  interaction->ResetBit(kISkipProcessChk);
  interaction->ResetBit(kISkipKinematicChk);
  
}
