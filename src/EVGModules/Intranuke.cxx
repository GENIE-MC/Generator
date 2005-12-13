//____________________________________________________________________________
/*!

\class   genie::Intranuke

\brief   The INTRANUKE cascading MC for intranuclear rescattering.

         Is a concrete implementation of the EventRecordVisitorI interface.

\ref     R.Merenyi et al., Phys.Rev.D45 (1992)
         R.D.Ransome, Nucl.Phys.B 139 (2005)

         The original INTRANUKE cascade MC was developed (in fortran) for the
         NeuGEN MC by G.F.Pearce, R.Edgecock, W.A.Mann and H.Gallagher.

\author  Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts University
         Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk> CCLRC, Rutherford Lab

\created September 20, 2005

*/
//____________________________________________________________________________

#include "Conventions/Constants.h"
#include "EVGModules/Intranuke.h"
#include "EVGModules/IntranukeConstants.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepOrder.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/NuclearUtils.h"

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

  // comment out this line to continue
  LOG("Intranuke", pWARN) << "(not ready yet)"; return;

  // Get the Interaction attached to the event record & get its InitialState
  // and Target objects
  Interaction * interaction = event_rec->GetInteraction();

  const InitialState & init_state = interaction->GetInitialState();
  const Target &       target     = init_state.GetTarget();

  // Return if the neutrino was not scatterred off a nuclear target
  if (! target.IsNucleus()) {
    LOG("Intranuke", pINFO) << "No nuclear target found - INTRANUKE exits";
    return;
  }
  if(target.A() != 56) {
    LOG("Intranuke", pWARN)
                    << "Can only select pion interactions for iron target";
    return;
  }

  // Get Intranuke configuration (or set defaults)

  double RoDef  = nuclear::Radius(target); // = R*A^(1/3), GeV

  bool   opaque = fConfig->GetDoubleDef("opaque", false);
  double ct0    = fConfig->GetDoubleDef("ct0",    kInukeFormationL); //GeV^-1
  double K      = fConfig->GetDoubleDef("Kpt2",   kInukeKpt2);
  double R0     = fConfig->GetDoubleDef("R0",     RoDef); // GeV

  R0  /= units::m; // GeV -> m
  ct0 /= units::m; // GeV -> m

  LOG("Intranuke", pDEBUG) << "ct0     = " << ct0 << " m";
  LOG("Intranuke", pDEBUG) << "K(pt^2) = " << K;
  LOG("Intranuke", pDEBUG) << "R0      = " << R0  << " m";

  // Get hadronic system's momentum vector
  TVector3 p3hadronic = this->Hadronic3P(event_rec);
  LOG("Intranuke", pINFO) << "P3Hadronic = " << print::P3AsString(&p3hadronic);

  // Get the random number generator
  RandomGen * rnd = RandomGen::Instance();

  // Generate a random vertex within the nuclear radius
  double R        = R0 * rnd->Random2().Rndm();
  double costheta = -1. + 2. * rnd->Random2().Rndm();
  double sintheta = TMath::Sqrt(1.-costheta*costheta);
  double fi       = 2 * kPi * rnd->Random2().Rndm();
  double cosfi    = TMath::Cos(fi);
  double sinfi    = TMath::Sin(fi);

  TVector3 vtx(R*sintheta*cosfi, R*sintheta*sinfi, R*costheta);

  LOG("Intranuke", pINFO) << "Vtx (in m) = " << print::Vec3AsString(&vtx);

  // Loop over the event record entries and run intranuclear rescattering
  // code for selected entries
  TObjArrayIter piter(event_rec);
  GHepParticle * p = 0;

  while( (p = (GHepParticle *) piter.Next()) ) {

    GHepStatus_t ist  = p->Status();
    int          pdgc = p->PdgCode();

    // Look for final state pi+,pi0,pi-
    bool infs   = (ist == kIStStableFinalState);
    bool ispi   = (pdgc==kPdgPiPlus || pdgc==kPdgPiMinus || pdgc==kPdgPi0);
    bool handle = (infs && ispi);

    if(!handle) continue; // <-- skip to next GHEP entry

    LOG("Intranuke", pINFO)
              << "*** Intranuke will attempt to rescatter a " << p->Name();

    TLorentzVector p4(*p->P4()); // current particle's 4-momentum (in GeV)
    TLorentzVector x4(*p->V4()); // current particle's 4-position (in m,sec)

    TVector3 p3  = p4.Vect();          // hadron's: p (px,py,pz)
    double   m   = p->Mass();          //           m
    double   m2  = m*m;                //           m^2
    double   P   = p4.P();             //           |p|
    double   Pt  = p3.Pt(p3hadronic);  //           pT
    double   Pt2 = Pt*Pt;              //           pT^2
    double   fz  = P*ct0*m/(m2+K*Pt2); //           formation zone, in m

    LOG("Intranuke", pINFO)
          << "|P| = " << P << " GeV, Pt = " << Pt
                              << " GeV, Formation Zone = " << fz << " m";

    // Advanced the intranuked hadron by the 'formation length' in 4-D
    LOG("Intranuke", pDEBUG) << "Advancing " << p->Name()
                             << " to account for its formation length";
    x4 = this->StepParticle(x4, p4, fz);

    // Check whether the hadron's position is still within the nucleus
    bool is_in = (x4.Vect().Mag() < R0);
    if(!is_in) {
       LOG("Intranuke", pINFO)
             << "Hadron is out of the nuclear radius. Done with it.";

       LOG("Intranuke", pDEBUG) << "Updating 4-position";
       p->SetVertex(x4);

       continue;
    }

    // Check if the "opaque" config var is true in which case absorb
    // the hadron anyway
    if(opaque) {
       // ?? Change status kIstStableFinalState -> kIstHadronInTheNucleus
       LOG("Intranuke", pINFO)
             << "In *OPAQUE* mode. Absorbing hadron & done with it.";
       p->SetStatus(kIstHadronInTheNucleus);
       p->SetVertex(x4);

       continue;
    }

    // ** Start stepping through nucleus & generating hadron interactions **

    bool go_on = true;
    while (go_on) {

      // Compute mean free path L and generate an 'interaction'
      // distance d from a exp(-d/L) distribution
      double L = this->MeanFreePath(target, p4, pdgc);
      double d = -1.*L * TMath::Log(rnd->Random2().Rndm());

      LOG("Intranuke", pDEBUG)
                << "Mean free path = " << L << " m, "
                              << "Generated path length = " << d << " m";

      // Check whether the hadron interacts (is it still inside the
      // nucleus after stepping it by d)
      x4 = this->StepParticle(x4, p4, d);

      bool interacts = (x4.Vect().Mag() < R0);

      // If it does not interact, advance, check whether it exits
      // the nucleus & continue...
      if(!interacts) {
        LOG("Intranuke", pINFO)
                 << "The intranuked " << p->Name()
                                << " has escaped the nucleus. Stopping!";
        go_on = false;
        continue;
      }

      // If it interacts, decide the hadron interaction type
      INukeProc_t inukp = this->PionFate(target, p4, pdgc);
      LOG("Intranuke", pINFO)
           << "Selected intranuke interaction for " << p->Name()
                                  << ": " << INukeProc::AsString(inukp);

      // Compute 4-p of the interacted hadron
      switch (inukp) {
        case kINukAbsorption:
             LOG("Intranuke", pINFO) << "Absorbing hadron & done with it.";
             p->SetStatus(kIstHadronInTheNucleus);
             p->SetVertex(x4);
             go_on = false;
             break;

        case kINukChargeExchange:
             LOG("Intranuke", pWARN) << "Can not do charge exchange yet.";
             break;
        case kINukInelastic:
             LOG("Intranuke", pWARN) << "Can not do inelastic scatering yet.";
             break;
        case kINukElastic:
             LOG("Intranuke", pWARN) << "Can not do elastic scatering yet.";
             break;
        default:
             LOG("Intranuke", pWARN) << "Oops! Unknown interaction type.";
             break;
      }

      // ... ... ... ... ...

    } // stepping
  }// stdhep entries

  //...

  LOG("Intranuke", pINFO) << "Done with this event";
}
//___________________________________________________________________________
double Intranuke::MeanFreePath(
            const Target & target, const TLorentzVector & p4, int pdgc) const
{
// Mean free path for the input GHEP record particle p in the input nuclear
// target.
// Adapted from NeuGEN's intranuke_piscat
// Original documentation:
//   Determine the scattering length LAMDA from yellowbook data and scale
//   it to the nuclear density.

  // compute pion kinetic energy K
  double E = p4.E();
  double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
  double K = 1000*(E-m); // in MeV

  // collision length in fermis
  double rlmbda;
  if      (K < 25 ) rlmbda = 200.;
  else if (K > 550) rlmbda = 25.;
  else              rlmbda = 1.0 / (-2.444/K + .1072 - .00011810*K);

  // additional correction for iron
  if(target.A() == 56) rlmbda = rlmbda/2.297;

  rlmbda *= (units::fermi/units::m); // fermi -> m

  // scale to nuclear density.
  //double density = kNucDensity;
  double density = 3.4; // units?

  double L  = (rlmbda/density); // units?
  return L;
}
//___________________________________________________________________________
INukeProc_t   Intranuke::PionFate(
            const Target & target, const TLorentzVector & p4, int pdgc) const
{
// Selects interaction type for the particle that is currently rescaterred.
// Adapted from NeuGEN's intranuke_pifate

  // compute pion kinetic energy K
  double E = p4.E();
  double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
  double K = 1000*(E-m); // in MeV

  int Kbin = TMath::Min (int(K/50.), intranuke::kPNDataPoints-1);

  RandomGen * rnd = RandomGen::Instance();
  double t = rnd->Random2().Rndm();

  if      ( t < intranuke::kPElastic[Kbin]    ) return kINukElastic;
  else if ( t < intranuke::kPInelastic[Kbin]  ) return kINukInelastic;
  else if ( t < intranuke::kPAbsorption[Kbin] ) return kINukAbsorption;
  else                                          return kINukChargeExchange;

  return kINukUndefined;
}
//___________________________________________________________________________
TLorentzVector Intranuke::StepParticle(
     const TLorentzVector & x4, const TLorentzVector & p4, double step) const
{
// Steps a particle starting at position x4 (in m, sec) moving along the
// direction of the input momentum by the input step (in m)

  LOG("Intranuke", pDEBUG) << "Stepping particle in nucleus:";

  TVector3 dr = p4.Vect().Unit();  // unit vector along momentum direction

  LOG("Intranuke", pDEBUG)
                << "dr = " << step
                         << " m, direction = " << print::Vec3AsString(&dr);
  LOG("Intranuke", pDEBUG)
                     << "X4[init] (in m,sec) = " << print::X4AsString(&x4);

  // spatial step size:
  dr.SetMag(step);

  // temporal step:
  double c  = kLightSpeed / (units::m/units::s); // c in m/sec
  double dt = step/c;

  TLorentzVector dx4(dr,dt);       // 4-vector step
  TLorentzVector x4new = x4 + dx4; // new position

  LOG("Intranuke", pDEBUG)
                  << "X4[new] (in m,sec) = " << print::X4AsString(&x4new);
  return x4new;
}
//___________________________________________________________________________
TVector3 Intranuke::Hadronic3P(GHepRecord * event_rec) const
{
  Interaction * interaction = event_rec->GetInteraction();

  int ipos = GHepOrder::StruckNucleonPosition(interaction);
  assert(ipos>0);

  GHepParticle * struck_nucleon = event_rec->GetParticle(ipos);

  int daughter1 = struck_nucleon->FirstDaughter();
  int daughter2 = struck_nucleon->LastDaughter();

  TLorentzVector p4had(0.,0.,0.,0.);
  for(int i=daughter1; i<=daughter2; i++) {
    TLorentzVector & p4 = *(event_rec->GetParticle(i)->P4());
    p4had = p4had + p4;
  }

  return p4had.Vect();  // (px,py,pz)
}
//___________________________________________________________________________

