//____________________________________________________________________________
/*!

\class   genie::COHPrimaryLeptonGenerator

\brief   Generates the final state primary lepton in v COH NC interactions.

         Is a concrete implementation of the EventRecordVisitorI interface.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 26, 2005

*/
//____________________________________________________________________________

#include "EVGModules/COHPrimaryLeptonGenerator.h"
#include "GHEP/GHepOrder.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Interaction/Interaction.h"

using namespace genie;

//___________________________________________________________________________
COHPrimaryLeptonGenerator::COHPrimaryLeptonGenerator() :
PrimaryLeptonGenerator("genie::COHPrimaryLeptonGenerator")
{

}
//___________________________________________________________________________
COHPrimaryLeptonGenerator::COHPrimaryLeptonGenerator(string config) :
PrimaryLeptonGenerator("genie::COHPrimaryLeptonGenerator", config)
{

}
//___________________________________________________________________________
COHPrimaryLeptonGenerator::~COHPrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
void COHPrimaryLeptonGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// This method generates the final state primary lepton

  //-- Get the attached interaction & get initial state & kinematics objects
  Interaction * interaction = evrec->GetInteraction();

  const InitialState & init_state = interaction->GetInitialState();
  const Kinematics &   kinematics = interaction->GetKinematics();

  //-- Coherent Scattering Kinematics: Compute the lepton energy and the
  //   scattering angle with respect to the incoming neutrino

  //auxiliary vars
  TParticlePDG * fsl = interaction->GetFSPrimaryLepton();
  int    pdgc = fsl->PdgCode();
  double Ev   = init_state.GetProbeE(kRfLab);
  double x    = kinematics.x();
  double y    = kinematics.y();
  double M    = init_state.GetTarget().Mass();
  double ml   = fsl->Mass();
  double ml2  = ml*ml;
  double Q2   = 2*x*y*M*Ev;

  //Compute outgoing lepton energy & momentum
  double El = (1-y)*Ev;
  double pl = TMath::Sqrt( TMath::Max(0., El*El-ml2) );

  //Compute outgoing lepton scat. angle with respect to the incoming v
  double cThSc = (El - 0.5*(Q2+ml2)/Ev) / pl; // cos(theta-scat) [-1,1]
  assert( TMath::Abs(cThSc) <= 1 );

  //-- Get the neutrino 4-p in LAB
  int nupos = GHepOrder::ProbePosition();
  TLorentzVector * p4nu = evrec->GetParticle(nupos)->GetP4();

  //-- Rotate its 4-momentum to the LAB
  //   unit' = R(Theta0,Phi0) * R(ThetaSc,PhiSc) * R^-1(Theta0,Phi0) * unit
  TLorentzVector * p4l = this->Rotate4P(p4nu, pdgc, cThSc, El);

  //-- Create a GHepParticle and add it to the event record
  //   (use the insertion method at the base PrimaryLeptonGenerator visitor)
  this->AddToEventRecord(evrec, pdgc, p4l);

  delete p4l;
  delete p4nu;
}
//___________________________________________________________________________
