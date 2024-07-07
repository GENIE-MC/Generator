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
     const AlgId & xsec_alg_id = xsec_alg->Id();
     std::string xsec_alg_name = xsec_alg_id.Name();
     std::string xsec_alg_config = xsec_alg_id.Config();
     std::string xsec_alg_key = xsec_alg_id.Key();
     interaction->SetBit(kISkipProcessChk);
     Kinematics * kine = interaction->KinePtr();
     kine->UseSelectedKinematics();
     std::cout << "\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
     std::cout << "Name: " << xsec_alg_name << ", Config: " << xsec_alg_config << ", Key: " << xsec_alg_key << "\n";
     double xsec;
     if ( proc_info.IsQuasiElastic() )
     {
        if (xsec_alg_name != "genie::SmithMonizQELCCPXSec" && xsec_alg_name != "genie::SuSAv2QELPXSec" && xsec_alg_key != "genie::HybridXSecAlgorithm/SuSAv2-QEL")
        {
            xsec = xsec_alg->XSec(interaction, kPSyphi0fEx);
        }
     }
     xsec = xsec_alg->XSec(interaction, kPSxyfE);
     //if ( proc_info.IsMEC() )
     //{
         //if (xsec_alg_name == "genie::NievesSimoVacasMECPXSec2016" || xsec_alg_name == "genie::SuSAv2MECPXSec" || xsec_alg_name == "genie::EmpiricalMECPXSec2015")
         
     //}
     std::cout << "xsec = " << std::scientific << xsec << "\n";
     std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n" << std::endl;
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
