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
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::utils;
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
  LOG("Intranuke", pINFO) << "Running INTRANUKE";

  // remove this line to continue
  //LOG("Intranuke", pWARN) << "(not ready yet)"; return; 

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

  // Get Intranuke configuration (or set defaults)
  // (uncomment when intranuke is really coded)

  bool   opaque = fConfig->GetDoubleDef("opaque", false);
  //double t0     = fConfig->GetDoubleDef("t0",    kInukeFormationTime);
  //double K      = fConfig->GetDoubleDef("Kpt2",  kInukeKpt2);
  //double Ro     = fConfig->GetDoubleDef("Ro",    kInukeNuclRadius);
  double t0  = 1.; // tmp
  double K   = 1.; // tmp
  double Ro  = 1.; // tmp

  // Get hadronic system's momentum vector

  TLorentzVector p4hadronic(0,0,0,0);
  // (.. .. .. fill p4hadronic from GHEP)
  TVector3 p3hadronic = p4hadronic.Vect(); // get px,py,pz

  // Get the random number generator
  RandomGen * rnd = RandomGen::Instance();

  // Generate a random vertex within the nuclear radius
  double R        = Ro * rnd->Random2().Rndm();
  double costheta = -1. + 2. * rnd->Random2().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->Random2().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  TVector3 vtx(R*sintheta*cosfi, R*sintheta*sinfi, R*costheta);

  LOG("Intranuke", pINFO) << "Vtx = " << print::Vec3AsString(&vtx);

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

         // check hadron's position after it is advanced by the 'formation
         // length' in 3-D

         TVector3 dr = p3.Unit(); // take a unit vector along momentum
         dr.SetMag(fz);           // set its length to be te formation length
         TVector3 r = vtx + dr;   // current position = vtx + dr

         LOG("Intranuke", pDEBUG) 
              << "The " << p->Name() << " stepped by " << print::Vec3AsString(&dr);
         LOG("Intranuke", pDEBUG) 
              << "The " << p->Name() << " is now at  " << print::Vec3AsString(&r);

         bool is_in = (r.Mag() < Ro);
         if(!is_in) {
            LOG("Intranuke", pINFO) 
                     << "Hadron is out of the nuclear radius. Done with it.";
            delete p4;
            continue;
         }

         // Start Intranuclear rescattering for current hadron

         // Check if the "opaque" config var is true in which case absorb 
         // the hadron anyway

         if(opaque) {
            // ?? The particle has a kIstStableFinalState status
            // ?? Now, mark it as kIstHadronInTheNucleus

            p->SetStatus(kIstHadronInTheNucleus);

            delete p4;
            continue;
         }

         // ... ... ... ... ... 

      } // is pi+,pi-.pi0
    } // is stable-final-state
  }// stdhep entries

  //...

  LOG("Intranuke", pINFO) << "Done with this event";
}
//___________________________________________________________________________
