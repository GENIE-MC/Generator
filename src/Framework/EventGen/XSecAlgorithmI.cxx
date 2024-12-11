//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include "Framework/EventGen/XSecAlgorithmI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGUtils.h"

using namespace genie;

//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI() :
Algorithm()
{
    fFinalLeptonPolarization = TVector3(0, 0, 0);
    fIsPreciseLeptonPolarization = false;
}
//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI(string name) :
Algorithm(name)
{
    fFinalLeptonPolarization = TVector3(0, 0, 0);
    fIsPreciseLeptonPolarization = false;
}
//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI(string name, string config) :
Algorithm(name, config)
{
    fFinalLeptonPolarization = TVector3(0, 0, 0);
    fIsPreciseLeptonPolarization = false;
}
//___________________________________________________________________________
XSecAlgorithmI::~XSecAlgorithmI()
{

}
//___________________________________________________________________________
bool XSecAlgorithmI::ValidKinematics(const Interaction* interaction) const
{
// can offer common implementation for all concrete x-section models because
// the input interaction is aware of its kinematic limits

  if ( interaction->TestBit(kISkipKinematicChk) ) return true;

  const KPhaseSpace& kps = interaction->PhaseSpace();

  if ( ! kps.IsAboveThreshold() ) {
     LOG("XSecBase", pINFO)  << "*** Below energy threshold";
     return false;
  }
  if ( ! kps.IsAllowed() ) {
     LOG("XSecBase", pINFO)  << "*** Not in allowed kinematical space";
     return false;
  }
  return true;
}
//___________________________________________________________________________
const TVector3 & XSecAlgorithmI::FinalLeptonPolarization (const Interaction* i) const
{
    int pdg = i->FSPrimLeptonPdg();
    if ( pdg::IsNeutrino(pdg) || pdg::IsElectron(pdg) || pdg::IsMuon(pdg) || pdg::IsTau(pdg) )
    {
        fFinalLeptonPolarization = TVector3(0, 0, -1);
    }
    else
    {
        fFinalLeptonPolarization = TVector3(0, 0, 1);
    }
    return fFinalLeptonPolarization;
}
//___________________________________________________________________________
