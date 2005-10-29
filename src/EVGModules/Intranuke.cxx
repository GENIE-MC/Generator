//____________________________________________________________________________
/*!

\class   genie::Intranuke

\brief   The INTRANUKE cascading MC for intranuclear rescattering.

         add description here
         this is a EventRecordVisitorI template

         Is a concrete implementation of the EventRecordVisitorI interface.

\author

\created Month xx, yyyy

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "EVGModules/Intranuke.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
Intranuke::Intranuke() :
EventRecordVisitorI("genie::Intranuke")
{

}
//___________________________________________________________________________
Intranuke::Intranuke(string config) :
EventRecordVisitorI("genie::Intranuke", config)
{

}
//___________________________________________________________________________
Intranuke::~Intranuke()
{

}
//___________________________________________________________________________
void Intranuke::ProcessEventRecord(GHepRecord * event_rec) const
{
  LOG("Intranuke", pDEBUG) << "Running INTRANUKE";

  // get the Interaction attached to the event record & get its InitialState
  // and Target objects
  Interaction * interaction = event_rec->GetInteraction();

  const InitialState & init_state = interaction->GetInitialState();
  const Target &       target     = init_state.GetTarget();

  // return if the neutrino was not scatterred off a nuclear target
  if (! target.IsNucleus()) {
    LOG("Intranuke", pINFO) << "No nuclear target found - INTRANUKE exits";
    return;
  }

  // Intranuke configuration (uncomment when intranuke is reallt coded)
  //double t0  = fConfig->GetDoubleDef("t0",    kInukeFormationTime);
  //double K   = fConfig->GetDoubleDef("Kpt2",  kInukeKpt2);
  //double Ro  = fConfig->GetDoubleDef("Ro",    kInukeNuclRadius);
  double t0  = 1.; // tmp
  double K   = 1.; // tmp
  double Ro  = 1.; // tmp

  // Get the random number generator

  RandomGen * rnd = RandomGen::Instance();

  // Get hadronic system's momentum vector

  TLorentzVector * p4hadronic = new TLorentzVector;
  // .. .. .. fill p4hadronic from GHEP
  TVector3 p3hadronic = p4hadronic->Vect(); // get px,py,pz

  // loop over the event record entries and look for final state pi+,pi-,pi0
  TObjArrayIter piter(event_rec);
  GHepParticle * p = 0;

  while( (p = (GHepParticle *) piter.Next()) ) {

    if( p->Status() == kIStStableFinalState ) {

      int  pdgc = p->PdgCode();
      bool ispi = (pdgc==kPdgPiPlus || pdgc==kPdgPiMinus || pdgc==kPdgPi0);

      if (ispi) {

         LOG("Intranuke", pINFO) 
                  << "Intranuke will attempt rescattering a " << p->Name();

         TLorentzVector * p4 = p->GetP4();
         TVector3         p3 = p4->Vect();

         double m   = p->Mass();          // hadron's m
         double m2  = m*m;                //          m^2
         double P   = p4->P();            //          |p|
         double Pt  = p3.Pt(p3hadronic);  //          pT
         double Pt2 = Pt*Pt;              //          pT^2

         double fz  = P*t0*m/(m2+K*Pt2);  // formation zone

         LOG("Intranuke", pINFO) 
            << "|P| = " << P << ", Pt = " << Pt << ", FormationZone = " << fz;

         // generate a random vertex within the nuclear radius

         double R        = Ro * rnd->Random2().Rndm();
         double costheta = -1. + 2. * rnd->Random2().Rndm();
         double sintheta = TMath::Sqrt(1.-costheta*costheta);
         double fi       = 2 * kPi * rnd->Random2().Rndm();
         double cosfi    = TMath::Cos(fi);
         double sinfi    = TMath::Sin(fi);

         double vtxx = R*sintheta*cosfi;
         double vtxy = R*sintheta*sinfi;
         double vtxz = R*costheta;

         LOG("Intranuke", pINFO) 
            << "(vtxx = " << vtxx << ", vtxy = " << vtxy 
                                               << ", vtxz = " << vtxz << ")";

         // ...

      } // is pi+,pi-.pi0
    } // is stable-final-state
  }// stdhep entries

  //...
}
//___________________________________________________________________________
