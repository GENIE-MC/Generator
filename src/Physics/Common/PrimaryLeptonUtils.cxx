//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
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

using namespace genie;
using namespace genie::utils;

//___________________________________________________________________________
void genie::utils::SetPrimaryLeptonPolarization( GHepRecord * ev )
{
// Moved out of the PrimaryLeptonGenerator class to make the same treatment
// accessible for generators that use a more unified approach (e.g.,
// QELEventGenerator and MECGenerator). -- S. Gardiner

// Set the final state lepton polarization. A mass-less lepton would be fully
// polarized. This would be exact for neutrinos and a very good approximation
// for electrons for the energies this generator is going to be used. This is
// not the case for muons and, mainly, for taus. I need to refine this later.
// How? See Kuzmin, Lyubushkin and Naumov, hep-ph/0312107

  // get the final state primary lepton
  GHepParticle * fsl = ev->FinalStatePrimaryLepton();
  if ( !fsl ) {
    LOG("LeptonicVertex", pERROR)
      << "Final state lepton not set yet! \n" << *ev;
    return;
  }

  // Get (px,py,pz) @ LAB
  TVector3 plab( fsl->Px(), fsl->Py(), fsl->Pz() );

  // In the limit m/E->0: leptons are left-handed and their anti-particles
  // are right-handed
  int pdgc = fsl->Pdg();
  if ( pdg::IsNeutrino(pdgc) || pdg::IsElectron(pdgc) ||
    pdg::IsMuon(pdgc) || pdg::IsTau(pdgc) )
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
}
