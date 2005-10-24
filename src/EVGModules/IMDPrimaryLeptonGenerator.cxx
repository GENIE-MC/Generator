//____________________________________________________________________________
/*!

\class   genie::IMDPrimaryLeptonGenerator

\brief   Generates the final state primary lepton in inverse muon decay (IMD)
         neutrino interaction events.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created July 13, 2005

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "EVGModules/IMDPrimaryLeptonGenerator.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
IMDPrimaryLeptonGenerator::IMDPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::IMDPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
IMDPrimaryLeptonGenerator::IMDPrimaryLeptonGenerator(string config):
PrimaryLeptonGenerator("genie::IMDPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
IMDPrimaryLeptonGenerator::~IMDPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void IMDPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton

  //-- Get the interaction & initial state objects & get the neutrino 4P and
  //   the kinematical variable
  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  TLorentzVector * p4nu = init_state.GetProbeP4(kRfLab); //ignore e velocity
  double y = interaction->GetScatteringParams().y();

  assert(1);
  assert(y>=0 && y<=1);

  // Final state primary lepton PDG code
  int pdgc = kPdgMuon;

  //-- Kinematics: Compute the lepton energy and the scattering angle with
  //   respect to the incoming neutrino

  double Ev    = p4nu->Energy();
  double Emu   = y*Ev;
  double Emu2  = TMath::Power(Emu,2);
  double me    = kElectronMass;
  double mmu2  = kMuonMass_2;
  double pmu   = TMath::Sqrt(Emu2-mmu2);      assert(Emu2>=mmu2);
  double Q2    = 2*(Ev-Emu)*me;
  double cThSc = (Emu-0.5*(Q2+mmu2)/Ev)/pmu;

  //warn about overflow in costheta and ignore it if it is small or abort
  if( TMath::Abs(cThSc)>1 ) {
     LOG("LeptonicVertex", pWARN)
       << "El = " << Emu << ", Ev = " << Ev << ", cos(theta) = " << cThSc;
     if(TMath::Abs(cThSc)-1<0.01) cThSc = 1.0;
  }
  assert(TMath::Abs(cThSc)<=1);

  //-- The computed scattering angle is with respect to the incoming neutrino
  //   direction. Rotate the lepton 4-momentum to the LAB according to
  //   unit' = R(Theta0,Phi0) * R(ThetaSc,PhiSc) * R^-1(Theta0,Phi0) * unit
  TLorentzVector * p4l = this->Rotate4P(p4nu, pdgc, cThSc, Emu);

  //-- Create a GHepParticle and add it to the event record
  //   (use the insertion method at the base PrimaryLeptonGenerator visitor)
  this->AddToEventRecord(evrec, pdgc, p4l);

  delete p4nu;
  delete p4l;
}
//___________________________________________________________________________
