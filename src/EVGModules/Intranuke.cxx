//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts University
         Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, CCLRC, Rutherford Lab
         September 20, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "EVGModules/Intranuke.h"
#include "EVGModules/IntranukeConstants.h"
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
void Intranuke::ProcessEventRecord(GHepRecord * evrec) const
{
  LOG("Intranuke", pNOTICE) << "************ Running INTRANUKE ************";

  // Return if the neutrino was not scatterred off a nuclear target
  GHepParticle * nucltgt = evrec->TargetNucleus();
  if (!nucltgt) {
    LOG("Intranuke", pINFO) << "No nuclear target found - INTRANUKE exits";
    return;
  }

  // Generate and set a vertex in the nucleus coordinate system
  this->GenerateVertex(evrec);

  // If in transparent mode, take all particle outside the nucleus and exit
  if(fIsTransparent) {
    this->TransportInTransparentNuc(evrec);
    return;
  }

  this->TransportInPhysicalNuc(evrec);

  LOG("Intranuke", pINFO) << "Done with this event";
}
//___________________________________________________________________________
void Intranuke::GenerateVertex(GHepRecord * evrec) const
{
// generate a vtx and set it to all GHEP physical particles

  GHepParticle * nucltgt = evrec->TargetNucleus();
  assert(nucltgt);
  this->SetNuclearRadius(nucltgt);

  RandomGen * rnd = RandomGen::Instance();

  double R     = fNuclRadius * rnd->Random1().Rndm();
  double cos9  = -1. + 2. * rnd->Random1().Rndm();    
  double sin9  = TMath::Sqrt(1.-cos9*cos9);   
  double fi    = 2 * kPi * rnd->Random1().Rndm();
  double cosfi = TMath::Cos(fi);
  double sinfi = TMath::Sin(fi);

  TVector3 vtx(R*sin9*cosfi,R*sin9*sinfi,R*cos9);

  LOG("Intranuke", pINFO) << "Vtx (in fm) = " << print::Vec3AsString(&vtx);

  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) )
  {
    if(p->IsFake()) continue;
    p->SetPosition(vtx.x(), vtx.y(), vtx.z(), 0.);
  }
}
//___________________________________________________________________________
void Intranuke::SetNuclearRadius(const GHepParticle * p) const
{
  assert(p && p->IsNucleus());

  int A = p->A();
  if(fR0>0) {
     fNuclRadius = fR0 * TMath::Power(A, 1./3.);
  } else {
     fNuclRadius  = nuclear::Radius(A); 
  }
}
//___________________________________________________________________________
bool Intranuke::NeedsRescattering(const GHepParticle * p) const
{
// checks whether the particle should be rescattered

  assert(p);
  return (p->Status() == kIStHadronInTheNucleus);
}
//___________________________________________________________________________
bool Intranuke::CanRescatter(const GHepParticle * p) const
{
// checks whether a particle that needs to be rescattered, can in fact be 
// rescattered by this cascade MC

  assert(p);
  return  ( p->PdgCode() == kPdgPiPlus  || 
            p->PdgCode() == kPdgPiMinus || 
            p->PdgCode() == kPdgPi0
          );
}
//___________________________________________________________________________
bool Intranuke::IsInNucleus(const GHepParticle * p) const
{
  return (p->X4()->Vect().Mag() < fNuclRadius);
}
//___________________________________________________________________________
void Intranuke::TransportInTransparentNuc(GHepRecord * evrec) const
{
// transport all hadrons assuming a transparent nucleus

  LOG("Intranuke", pNOTICE) << "INTRANUKE is set in **transparent** mode";

  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  int icurr=-1;

  while( (p = (GHepParticle *) piter.Next()) )
  {
    icurr++;

    // Check whether the particle needs rescattering, otherwise skip it
    if( ! this->NeedsRescattering(p) ) continue;

    LOG("Intranuke", pINFO)
       << "Transporting " << p->Name() << " out of the nuclear target";

    // move it outside the nucleus
    GHepParticle * cp = new GHepParticle(*p); // create a clone

    cp->SetFirstMother(icurr);                // clone's mother
    cp->SetStatus(kIStStableFinalState);      // mark it & done with it

    evrec->AddParticle(*cp); // add it at the event record

    //LOG("Intranuke", pDEBUG) << *evrec;
  }
}
//___________________________________________________________________________
void Intranuke::TransportInPhysicalNuc(GHepRecord * evrec) const
{
// transport all hadrons 

  // Loop over GHEP and run intranuclear rescattering on handled particles
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  int icurr = -1;

  while( (p = (GHepParticle *) piter.Next()) )
  {
    icurr++;

    // Check whether the particle needs rescattering, otherwise skip it
    if( ! this->NeedsRescattering(p) ) continue;

    LOG("Intranuke", pINFO)
             << "*** Intranuke will attempt to rescatter a " << p->Name();

    // Rescatter a clone, not the original particle
    GHepParticle * sp = new GHepParticle(*p); 

    // Check whether the particle can be rescattered 
    if(!this->CanRescatter(sp)) {

       // if I can't rescatter it, I will just take it out of the nucleus
       LOG("Intranuke", pINFO)
                   << "Current version can't rescatter a " << sp->Name();
       sp->SetFirstMother(icurr); 
       sp->SetStatus(kIStStableFinalState);
       evrec->AddParticle(*sp);
       continue; // <-- skip to next GHEP entry
    }

    // Check whether it is a 'fresh' hadron 
    bool fresh = this->IsFreshHadron(evrec,sp);

    // Set clone's mom to be the hadron that was cloned
    sp->SetFirstMother(icurr); 

    // Advance a 'fresh' hadron by the formation zone, store and continue 
    if (fresh) {
       this->AdvanceFreshHadron(evrec,sp);

       // Check whether the hadron's position is still within the nucleus
       if(!this->IsInNucleus(sp)) {
           LOG("Intranuke", pINFO)
                  << "Hadron went out of the nucleus! Done with it.";
           sp->SetStatus(kIStStableFinalState);
       }
       evrec->AddParticle(*sp);
       continue; // <-- skip to next GHEP entry
    }

    // Step in the nucleus
    double d = this->GenerateStep(evrec,sp);
    this->StepParticle(sp, d);

    // Check whether it interacts
    bool interacts = this->IsInNucleus(sp);

    // If it does not interact, it must exit the nucleus. Done with it! 
    if(!interacts) {
       LOG("Intranuke", pINFO) << "*** Hadron escaped the nucleus!";
       sp->SetStatus(kIStStableFinalState);
       evrec->AddParticle(*sp);
       continue; 
    }

    // The stepped particle interacts. Simulate the hadronic interaction.
    this->SimHadronicInteraction(evrec,sp);

  }// GHEP entries
}
//___________________________________________________________________________
bool Intranuke::IsFreshHadron(GHepRecord* evrec, GHepParticle* p) const
{
// Decide whether the particle p is a 'fresh' hadron (direct descendant of 
// the 'HadronicSystem' GHEP entry)

  int imom = p->FirstMother();
  assert(imom>0);

  if(evrec->Particle(imom)->PdgCode() == kPdgHadronicSyst) return true;

  return false;
}
//___________________________________________________________________________
double Intranuke::FormationZone(GHepRecord* evrec, GHepParticle* p) const
{
// Compute the formation zone for the particle p of the input event

  // Get hadronic system's 3-momentum
  GHepParticle * hadronic_system = evrec->FinalStateHadronicSystem();
  TVector3 p3hadr = hadronic_system->P4()->Vect(); // (px,py,pz)

  // Compute formation zone
  TVector3 p3  = p->P4()->Vect();      // hadron's: p (px,py,pz)
  double   m   = p->Mass();            //           m
  double   m2  = m*m;                  //           m^2
  double   P   = p->P4()->P();         //           |p|
  double   Pt  = p3.Pt(p3hadr);        //           pT
  double   Pt2 = Pt*Pt;                //           pT^2
  double   fz  = P*fct0*m/(m2+fK*Pt2); //           formation zone, in m

  LOG("Intranuke", pINFO)
          << "|P| = " << P << " GeV, Pt = " << Pt
                              << " GeV, Formation Zone = " << fz << " m";
  return fz;
}
//___________________________________________________________________________
void Intranuke::AdvanceFreshHadron(GHepRecord* evrec, GHepParticle* cp) const
{
// Advance 'fresh' hadrons by a formation zone. Note that the input particle
// may be modified but the event record not (ie the particle has been cloned 
// from an event entry to undergo rescattering but it should not have been
// added at the event just yet)

  double fzone = this->FormationZone(evrec, cp);
  this->StepParticle(cp, fzone);
}
//___________________________________________________________________________
double Intranuke::GenerateStep(GHepRecord* evrec, GHepParticle* p) const
{
// Generate a step (in fermis) for particle p in the input event.
// Computes the mean free path L and generate an 'interaction' distance d from 
// an exp(-d/L) distribution

  RandomGen * rnd = RandomGen::Instance();

  double L = this->MeanFreePath(evrec, p);
  double d = -1.*L * TMath::Log(rnd->Random1().Rndm());

  LOG("Intranuke", pDEBUG)
            << "Mean free path = " << L << " fm, "
                              << "Generated path length = " << d << " fm";
  return d;
}
//___________________________________________________________________________
double Intranuke::MeanFreePath(GHepRecord* evrec, GHepParticle* p) const
{
// Mean free path for the pions with the input kinetic energy.
// Adapted from NeuGEN's intranuke_piscat
// Original documentation:
//   Determine the scattering length LAMDA from yellowbook data and scale
//   it to the nuclear density.

  // compute pion kinetic energy K in MeV
  double K = (p->KinE() / units::MeV);

  // collision length in fermis
  double rlmbda;
  if      (K < 25 ) rlmbda = 200.;
  else if (K > 550) rlmbda = 25.;
  else              rlmbda = 1.0 / (-2.444/K + .1072 - .00011810*K);

  // additional correction for iron
  GHepParticle * target = evrec->TargetNucleus();
  if(target->Z()==26) {
      rlmbda = rlmbda/2.297;
  }

  // scale to nuclear density.
  //double density = kNucDensity;
  double density = 3.4; // units?

  double L  = (rlmbda/density); // units?
  return L;
}
//___________________________________________________________________________
void Intranuke::SimHadronicInteraction(
                                    GHepRecord* evrec, GHepParticle* p) const
{
  HadroProc_t proc = this->HadronFate(p);

  switch (proc) {
    case kHPcAbsorption:
        this->SimAbsorption(evrec, p);
        break;
    case kHPcChargeExchange:
        this->SimChargeExchange(evrec, p);
        break;
    case kHPcInelastic:
        this->SimInelasticScattering(evrec, p);
        break;
    case kHPcElastic:
        this->SimElasticScattering(evrec, p);
        break;
    default:
        LOG("Intranuke", pFATAL) << "Unknown INTRANUKE interaction type.";
        exit(1);
  }
}
//___________________________________________________________________________
void Intranuke::SimAbsorption(GHepRecord* evrec, GHepParticle* p) const
{
   LOG("Intranuke", pWARN) << "Can not do hadron absorption yet.";
}
//___________________________________________________________________________
void Intranuke::SimChargeExchange(
			   GHepRecord* /*evrec*/, GHepParticle* /*p*/) const
{
   LOG("Intranuke", pWARN) << "Can not do charge exchange yet.";
}
//___________________________________________________________________________
void Intranuke::SimInelasticScattering(
			   GHepRecord* /*evrec*/, GHepParticle* /*p*/) const
{
   LOG("Intranuke", pWARN) << "Can not do inelastic scatering yet.";
}
//___________________________________________________________________________
void Intranuke::SimElasticScattering(
			   GHepRecord* /*evrec*/, GHepParticle* /*p*/) const
{
   LOG("Intranuke", pWARN) << "Can not do elastic scatering yet.";
}
//___________________________________________________________________________
HadroProc_t Intranuke::HadronFate(const GHepParticle * p) const
{
// Selects interaction type for the particle that is currently rescaterred.
// Adapted from NeuGEN's intranuke_pifate

  // compute pion kinetic energy K in MeV
  double K = (p->KinE() / units::MeV);

  // find the kinetic energy bin in the cummulative interaction prob array
  int Kbin = TMath::Min (int(K/50.), intranuke::kPNDataPoints-1);

  // select rescattering type
  RandomGen * rnd = RandomGen::Instance();
  double t = rnd->Random1().Rndm();

  HadroProc_t proc = kHPcUndefined;

  if      ( t < intranuke::kPElastic[Kbin]    ) proc = kHPcElastic;
  else if ( t < intranuke::kPInelastic[Kbin]  ) proc = kHPcInelastic;
  else if ( t < intranuke::kPAbsorption[Kbin] ) proc = kHPcAbsorption;
  else                                          proc = kHPcChargeExchange;

  LOG("Intranuke", pINFO)
        << "Selected hadronic process for " 
                        << p->Name() << ": " << HadroProc::AsString(proc);
  return proc;
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
void Intranuke::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void Intranuke::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//___________________________________________________________________________
void Intranuke::LoadConfig(void)
{
  AlgConfigPool * confp = AlgConfigPool::Instance();
  const Registry * gc = confp->GlobalParameterList();

  fIsTransparent = fConfig->GetBoolDef("transparent-mode", false);

  fct0 = fConfig->GetDoubleDef ("ct0",  gc->GetDouble("INUKE-FormationZone")); // fermi
  fK   = fConfig->GetDoubleDef ("Kpt2", gc->GetDouble("INUKE-KPt2"));
  fR0  = fConfig->GetDoubleDef ("R0",   gc->GetDouble("INUKE-Ro")); // fermi

  LOG("Intranuke", pDEBUG) << "IsTransparent = " << fIsTransparent;
  LOG("Intranuke", pDEBUG) << "ct0           = " << fct0 << " fermi";
  LOG("Intranuke", pDEBUG) << "K(pt^2)       = " << fK;
  LOG("Intranuke", pDEBUG) << "R0            = " << fR0  << " fermi";
}
//___________________________________________________________________________


