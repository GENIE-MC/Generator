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
  LOG("Intranuke", pINFO) << "************ Running INTRANUKE ************";

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
  RoDef  /= units::m; // GeV -> m

  double R0 = (fR0>0) ? fR0 : RoDef;

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
  this->SetVtxPosition(event_rec, vtx);

  LOG("Intranuke", pINFO) << "Vtx (in m) = " << print::Vec3AsString(&vtx);

  // Loop over the event record entries and run intranuclear rescattering
  // code for selected entries
  TObjArrayIter piter(event_rec);
  GHepParticle * p = 0;
  int icurr = -1;

  while( (p = (GHepParticle *) piter.Next()) )
  {
    icurr++;

    // check whether the particle can be rescattered or continue
    if(!this->CanRescatter(p)) continue; // <-- skip to next GHEP entry

    // check whether the particle is already outside the nucleus
    // (eg, if it is a f/s pion added to GHEP at previous intranuke steps)
    if(!this->IsInNucleus(p, R0)) continue; // <-- skip to next GHEP entry

    LOG("Intranuke", pINFO)
           << "*** Intranuke will attempt to rescatter a " << p->Name();

    p->SetStatus(kIstHadronInTheNucleus);     // <-- mark it & done with it
    GHepParticle * sp = new GHepParticle(*p); // <-- rescater a clone

    sp->SetFirstMother(icurr); // set its unscattered clone as its mom

    TVector3 p3  = sp->P4()->Vect();     // hadron's: p (px,py,pz)
    double   m   = sp->Mass();           //           m
    double   m2  = m*m;                  //           m^2
    double   P   = sp->P4()->P();        //           |p|
    double   Pt  = p3.Pt(p3hadronic);    //           pT
    double   Pt2 = Pt*Pt;                //           pT^2
    double   fz  = P*fct0*m/(m2+fK*Pt2); //           formation zone, in m

    LOG("Intranuke", pINFO)
          << "|P| = " << P << " GeV, Pt = " << Pt
                              << " GeV, Formation Zone = " << fz << " m";

    // Advanced the intranuked hadron by the 'formation length' in 4-D
    this->StepParticle(sp, fz);

    // Check whether the hadron's position is still within the nucleus
    if(!this->IsInNucleus(sp, R0)) {
       LOG("Intranuke", pINFO)
                  << "Hadron went out of the nucleus! Done with it.";
       sp->SetStatus(kIStStableFinalState);
       event_rec->AddParticle(*sp);
    }

    // Check if the "opaque" config var is true in which case absorb
    // the hadron anyway
    if(fIsOpaque) {
       // ?? Change status kIstStableFinalState -> kIstHadronInTheNucleus
       LOG("Intranuke", pINFO)
             << "In *OPAQUE* mode. Absorbing hadron & done with it.";
       sp->SetStatus(kIstHadronInTheNucleus);
       event_rec->AddParticle(*sp);
       continue;
    }

    // ** Start stepping through nucleus & generating hadron interactions **

    bool go_on = true;
    while (go_on) {
      // Compute mean free path L and generate an 'interaction'
      // distance d from a exp(-d/L) distribution
      double K = sp->KinE();
      double L = this->MeanFreePath(K);
      double d = -1.*L * TMath::Log(rnd->Random2().Rndm());

      LOG("Intranuke", pDEBUG)
                << "Mean free path = " << L << " m, "
                              << "Generated path length = " << d << " m";

      // Check whether the hadron interacts (is it still inside the
      // nucleus after stepping it by d)
      this->StepParticle(sp, d);

      bool interacts = this->IsInNucleus(sp, R0);

      // If it does not interact, it must exit the nucleus.
      // Done with it! Add it to GHEP and continue.
      if(!interacts) {
        LOG("Intranuke", pINFO)
                 << "The intranuked " << sp->Name()
                              << " has escaped the nucleus. Stopping!";
        event_rec->AddParticle(*sp);
        go_on = false;
        continue;
      }

      // If it interacts, decide the hadron interaction type
      INukeProc_t inukp = this->ParticleFate(sp);

      // Compute 4-p of the interacted hadron
      switch (inukp) {
        case kINukAbsorption:
             LOG("Intranuke", pINFO) << "Absorbing hadron & done with it.";
             sp->SetStatus(kIstHadronInTheNucleus);
             event_rec->AddParticle(*sp);
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
bool Intranuke::CanRescatter(const GHepParticle * p) const
{
// checks whether the particle can be rescattered by this cascade MC

  GHepStatus_t ist  = p->Status();
  int          pdgc = p->PdgCode();

  // is it a final state particle
  bool in_fs = (ist == kIStStableFinalState);

  // is it a pi+,pi-,pi0?
  bool is_pi = (pdgc==kPdgPiPlus || pdgc==kPdgPiMinus || pdgc==kPdgPi0);

  // rescatter only f/s pions
  bool ican = (in_fs && is_pi);

  return ican;
}
//___________________________________________________________________________
double Intranuke::MeanFreePath(double KinE) const
{
// Mean free path for the pions with the input kinetic energy.
// Adapted from NeuGEN's intranuke_piscat
// Original documentation:
//   Determine the scattering length LAMDA from yellowbook data and scale
//   it to the nuclear density.

  // pion kinetic energy in MeV
  double K = (KinE/units::MeV);

  // collision length in fermis
  double rlmbda;
  if      (K < 25 ) rlmbda = 200.;
  else if (K > 550) rlmbda = 25.;
  else              rlmbda = 1.0 / (-2.444/K + .1072 - .00011810*K);

  // additional correction for iron
  rlmbda = rlmbda/2.297;
  rlmbda *= (units::fermi/units::m); // fermi -> m

  // scale to nuclear density.
  //double density = kNucDensity;
  double density = 3.4; // units?

  double L  = (rlmbda/density); // units?
  return L;
}
//___________________________________________________________________________
INukeProc_t Intranuke::ParticleFate(const GHepParticle * p) const
{
// Selects interaction type for the particle that is currently rescaterred.
// Adapted from NeuGEN's intranuke_pifate

  // compute pion kinetic energy K in MeV
  double K = (p->KinE() / units::MeV);

  // find the kinetic energy bin in the cummulative interaction prob array
  int Kbin = TMath::Min (int(K/50.), intranuke::kPNDataPoints-1);

  // select rescattering type
  RandomGen * rnd = RandomGen::Instance();
  double t = rnd->Random2().Rndm();

  INukeProc_t inukp = kINukUndefined;

  if      ( t < intranuke::kPElastic[Kbin]    ) inukp = kINukElastic;
  else if ( t < intranuke::kPInelastic[Kbin]  ) inukp = kINukInelastic;
  else if ( t < intranuke::kPAbsorption[Kbin] ) inukp = kINukAbsorption;
  else                                          inukp = kINukChargeExchange;

  LOG("Intranuke", pINFO)
           << "Selected intranuke interaction for " << p->Name()
                                     << ": " << INukeProc::AsString(inukp);
  return inukp;
}
//___________________________________________________________________________
void Intranuke::StepParticleMF(GHepParticle * p, double step, int /*Z*/) const
{
// Steps a particle starting from its current position (in m, sec) and moving
// along the direction of its current momentum by the input step (in m).
// If the particle is neutral it is stepped in a straight line.
// If the particle is charged it is stepped in a curved trajectory depending
// on its q/m and the mean field of the nucleus.

  LOG("Intranuke", pDEBUG)
      << "Stepping particle [" << p->Name() << "] by dr = " << step << " m";

  //..
  //..
}
//___________________________________________________________________________
void Intranuke::StepParticle(GHepParticle * p, double step) const
{
// Steps a particle starting from its current position (in m, sec) and moving
// along the direction of its current momentum by the input step (in m).
// The particle is stepped in a straight line.

  LOG("Intranuke", pDEBUG)
      << "Stepping particle [" << p->Name() << "] by dr = " << step << " m";

  TVector3 dr = p->P4()->Vect().Unit();  // unit vector along its direction

  LOG("Intranuke", pDEBUG)
               << "Init direction = " << print::Vec3AsString(&dr);
  LOG("Intranuke", pDEBUG)
       << "Init position (in m,sec) = " << print::X4AsString(p->X4());

  // spatial step size:
  dr.SetMag(step);

  // temporal step:
  double c  = kLightSpeed / (units::m/units::s); // c in m/sec
  double dt = step/c;

  TLorentzVector dx4(dr,dt);       // 4-vector step
  TLorentzVector x4new = *(p->X4()) + dx4; // new position

  LOG("Intranuke", pDEBUG)
                  << "X4[new] (in m,sec) = " << print::X4AsString(&x4new);
  p->SetPosition(x4new);
}
//___________________________________________________________________________
bool Intranuke::IsInNucleus(const GHepParticle * p, double R0) const
{
  return (p->X4()->Vect().Mag() < R0);
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
void Intranuke::SetVtxPosition(GHepRecord * event_rec, TVector3 & v) const
{
// set the interaction vtx position in the target nucleus coordinate system

  Interaction * interaction = event_rec->GetInteraction();

  int ip = GHepOrder::ProbePosition();
  int in = GHepOrder::StruckNucleonPosition(interaction);

  GHepParticle * probe = event_rec->GetParticle(ip);
  assert(probe);
  probe->SetPosition(v.x(), v.y(), v.z(), 0.);
  GHepParticle * nucleon = event_rec->GetParticle(in);
  if(nucleon) {
    nucleon->SetPosition(v.x(), v.y(), v.z(), 0.);
  }
}
//___________________________________________________________________________
void Intranuke::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
}
//___________________________________________________________________________
void Intranuke::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfigData();
}
//___________________________________________________________________________
void Intranuke::LoadConfigData(void)
{
  fIsOpaque = fConfig->GetBoolDef   ("opaque", false);
  fct0      = fConfig->GetDoubleDef ("ct0",    kInukeFormationL); //GeV^-1
  fK        = fConfig->GetDoubleDef ("Kpt2",   kInukeKpt2);
  fR0       = fConfig->GetDoubleDef ("R0",     -1.); // GeV

  fR0  /= units::m; // GeV -> m
  fct0 /= units::m; // GeV -> m

  LOG("Intranuke", pDEBUG) << "IsOpaque = " << fIsOpaque;
  LOG("Intranuke", pDEBUG) << "ct0      = " << fct0 << " m";
  LOG("Intranuke", pDEBUG) << "K(pt^2)  = " << fK;
  LOG("Intranuke", pDEBUG) << "R0       = " << fR0  << " m";
}
//___________________________________________________________________________


