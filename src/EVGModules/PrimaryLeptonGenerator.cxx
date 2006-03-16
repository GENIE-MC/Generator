//____________________________________________________________________________
/*!

\class   genie::PrimaryLeptonGenerator

\brief   Abstract class. Is used to pass common implementation to concrete
         implementations of the EventRecordVisitorI interface generating the
         primary lepton for a specific processes (QEL,DIS,RES,IMD,...)

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created October 03, 2004

*/
//____________________________________________________________________________

#include <TVector3.h>
#include <TLorentzVector.h>

#include "Conventions/Constants.h"
#include "EVGModules/PrimaryLeptonGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepOrder.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGUtils.h"
#include "Numerical/RandomGen.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
PrimaryLeptonGenerator::PrimaryLeptonGenerator() :
EventRecordVisitorI()
{

}
//___________________________________________________________________________
PrimaryLeptonGenerator::PrimaryLeptonGenerator(string name) :
EventRecordVisitorI(name)
{

}
//___________________________________________________________________________
PrimaryLeptonGenerator::PrimaryLeptonGenerator(string name, string config) :
EventRecordVisitorI(name, config)
{

}
//___________________________________________________________________________
PrimaryLeptonGenerator::~PrimaryLeptonGenerator()
{

}
//___________________________________________________________________________
TVector3 * PrimaryLeptonGenerator::NucRestFrame2Lab(GHepRecord * evrec) const
{
// Velocity for an active Lorentz transform taking the final state primary
// lepton from the [nucleon rest frame] --> [LAB]

  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  TLorentzVector * pnuc4 = init_state.GetTarget().StruckNucleonP4(); //[LAB]

  double bx = pnuc4->Px() / pnuc4->Energy();
  double by = pnuc4->Py() / pnuc4->Energy();
  double bz = pnuc4->Pz() / pnuc4->Energy();

  TVector3 * b = new TVector3(bx,by,bz);
  return b;
}
//___________________________________________________________________________
TLorentzVector * PrimaryLeptonGenerator::P4InNucRestFrame(
                           GHepRecord * evrec, double cThSc, double El) const
{
// Takes the final state primary lepton scattering angle (with respect to the
// incoming neutrino direction) and its generated energy and computes its 4-P
// nucleon rest frame.
// Inputs:
//   -- evrec: current event record
//   -- cThSc: generated cosine(theta-angle) relative to the v direction
//   -- El:    generated lepton energy
//
// cThSc and El depend on the actual kinematics and are calculated by more
// specialized, higher level EventRecordVisitors that subclass this ABC.

  Interaction * interaction = evrec->GetInteraction();
  const InitialState & init_state = interaction->GetInitialState();

  // Figure out the final state primary lepton
  int pdgc = interaction->GetFSPrimaryLepton()->PdgCode();

  // Take the incoming neutrino 4-momentum in the nucleon rest frame
  TLorentzVector * p4nu = init_state.GetProbeP4(kRfStruckNucAtRest);

  // Rotate the final state lepton 4-p to the reference frame of the
  // input neutrino direction
  TLorentzVector * p4l = this->Rotate4P(p4nu, pdgc, cThSc, El);

  delete p4nu;
  return p4l;
}
//___________________________________________________________________________
TLorentzVector * PrimaryLeptonGenerator::Rotate4P(
              TLorentzVector * p4nu, int pdgc, double cThSc, double El) const
{
// Rotate the final state 4-P from the reference frame where the scattering
// angle is measured with respect to the incoming neutrino direction, to the
// reference frame in which the input neutrino 4-P is given (eg nucleon rest
// frame, LAB,...)
// It needs the:
//   -- p4nu:  neutrino 4-P
//   -- pdgc:  final state primary lepton pdg code
//   -- cThSc: generated cosine(theta-angle) relative to the v direction
//   -- El:    generated lepton energy
//
// The direction of the final state lepton is given by rotating the unit
// vector along the input neutrino direction as:
//
//   u' = R(Theta0,Phi0) * R(ThetaSc,PhiSc) * R^-1(Theta0,Phi0) * u
//
// where
//   Theta0, Phi0 are the v zenith and azimuth angle in the frame where p4nu
//   is given, ThetaSc, PhiSc are the angles of the emerging final state
//   lepton with respect to the incoming v, R is a rotation matrix and R^-1
//   its inverse.
// For simplicity the rotation matrix multiplications have already been
// carried out.


  RandomGen * rnd = RandomGen::Instance();

  // Compute azimuthal angle [uniform over [0, 2*pi]:
  double PhiSc = 2 * kPi * (rnd->Random1().Rndm());

  // Compute the lepton |momentum|
  double ml  = PDGLibrary::Instance()->Find(pdgc)->Mass();
  double ml2 = TMath::Power(ml,2);
  double pl  = TMath::Sqrt( TMath::Max(0., El*El-ml2) );

  LOG("LeptonicVertex", pINFO)
              << "f/s prim. lepton: E = " << El << ", |p| = " << pl;
  LOG("LeptonicVertex", pINFO)
      << "f/s prim. lepton: "
           << "cos(theta-sc) = " << cThSc << ", phi-sc = " << PhiSc;

  //-- Compute the remaining needed trigonometric numbers
  double sThSc = TMath::Sqrt(1. - cThSc*cThSc); // sin(theta-scattering)
  double sPhSc = TMath::Sin(PhiSc);             // sin(phi-scattering)
  double cPhSc = TMath::Cos(PhiSc);             // cos(phi-scattering)

  //-- Go from the (theta,fi)-coordinates relative to the v direction to
  //   (theta,fi)-coordinates

  // (theta,fi) coordinates of the neutrino
  double Theta0 = p4nu->Theta();
  double Phi0   = p4nu->Phi();

  // trigonometric numbers involved in rotation of reference frames
  double sTh0 = TMath::Sin( Theta0 );
  double cPh0 = TMath::Cos( Phi0   );
  double sPh0 = TMath::Sin( Phi0   );

  // unit vector along the direction of the neutrino
  TVector3 unit = p4nu->BoostVector().Unit();

  // unit' = R(Theta0,Phi0) * R(ThetaSc,PhiSc) * R^-1(Theta0,Phi0) * unit
  //
  double plx = unit.x() * cThSc + sThSc * (unit.z()*cPhSc*cPh0 - sPhSc*sPh0);
  double ply = unit.y() * cThSc + sThSc * (unit.z()*cPhSc*sPh0 + sPhSc*cPh0);
  double plz = unit.z() * cThSc - sThSc * sTh0 * cPhSc;

  LOG("LeptonicVertex", pINFO)
         << "i/s lepton direction [NRF]: u = (" << unit.x()
                            << ", " << unit.y()  << ", " << unit.z() << ")";
  LOG("LeptonicVertex", pINFO)
         << "f/s lepton direction [NRF]: u = ("
                               << plx << ", " << ply  << ", " << plz << ")";

  TVector3 plv(plx, ply, plz);
  plv.SetMag(pl);

  //-- Output 4-momentum
  TLorentzVector * p4l = new TLorentzVector(plv, El);
  return p4l;
}
//___________________________________________________________________________
void PrimaryLeptonGenerator::AddToEventRecord(
              GHepRecord * evrec, int pdgc, const TLorentzVector * p4) const
{
// Adds the final state primary lepton GHepParticle to the event record.
// To be called by all concrete PrimaryLeptonGenerators before exiting.

  int mom = GHepOrder::ProbePosition();
  TLorentzVector vdummy(0,0,0,0); // position 4-vector

  evrec->AddParticle(pdgc, kIStStableFinalState, mom,-1,-1,-1, *p4, vdummy);
}
//___________________________________________________________________________
void PrimaryLeptonGenerator::SetPolarization(GHepRecord * ev) const
{
// Set the final state lepton polarization. A mass-less lepton would be fully
// polarized. This would be exact for neutrinos and a very good approximation
// for electrons for the energies this generator is going to be used. This is
// not the case for muons and, mainly, for taus. I need to refine this later.
// How? See Kuzmin, Lyubushkin and Naumov, hep-ph/0312107

  // get the final state primary lepton
  GHepParticle * prb = ev->GetParticle( GHepOrder::ProbePosition() );
  assert(prb);
  GHepParticle * fsl = ev->GetParticle( prb->FirstDaughter() );

  if(!fsl) {
     LOG("LeptonicVertex", pERROR)
                    << "Final state lepton not set yet! \n" << *ev;
     return;
  }

  // get (px,py,pz) @ LAB
  TVector3 plab(fsl->Px(), fsl->Py(), fsl->Pz());

  // in the limit m/E->0: leptons are left-handed and their anti-particles
  // are right-handed
  int pdgc = fsl->PdgCode();
  if(pdg::IsNeutrino(pdgc) || pdg::IsElectron(pdgc) ||
                    pdg::IsMuon(pdgc) || pdg::IsTau(pdgc) ) {
    plab *= -1; // left-handed
  }

  LOG("LeptonicVertex", pINFO) 
            << "Setting polarization angles for particle: " << fsl->Name();

  fsl->SetPolarization(plab);

  if(fsl->PolzIsSet()) {
     LOG("LeptonicVertex", pINFO) 
          << "Polarization (rad): Polar = "  << fsl->PolzPolarAngle() 
                           << ", Azimuthal = " << fsl->PolzAzimuthAngle();
  }
}
//___________________________________________________________________________

