//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Steve Dytman <dytman+@pitt.edu>, Pittsburgh Univ.
         Hugh Gallagher <gallag@minos.phy.tufts.edu>, Tufts Univ.
         Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>, Rutherford Lab.
         September 20, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <cstdlib>
#include <sstream>

#include <TMath.h>

#include "Algorithm/AlgConfigPool.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepParticle.h"
#include "HadronTransport/Intranuke.h"
#include "HadronTransport/INukeHadroData.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/NuclearUtils.h"

using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;

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

  // Transport all particles outside the nucleus and exit
  this->TransportHadrons(evrec);

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
  TVector3 vtx(0.,0.,0.);

  // *** For h+A events (test mode): 
  // Generate vtx on the nuclear radius. Put the vtx on the
  // 'correct' half of the sphere so that the hadron does not
  // point away from the nucleus

  if(fInTestMode) {
    // generate a point in spherical coordinates
    double R     = fNuclRadius;
    double cos9  = -1. + 2. * rnd->RndFsi().Rndm();    
    double sin9  = TMath::Sqrt(1.-cos9*cos9);   
    double fi    = 2 * kPi * rnd->RndFsi().Rndm();
    double cosfi = TMath::Cos(fi);
    double sinfi = TMath::Sin(fi);

    // polar -> rectangular coordinates assuming the unit vector along 
    // the hadron direction is (0,0,1)
    vtx.SetXYZ(R*sin9*cosfi,R*sin9*sinfi, -1. * TMath::Abs(R*cos9));

    // get the actual unit vector along the incoming hadron direction
    TVector3 direction = evrec->Probe()->P4()->Vect().Unit();

    // rotate the vtx position
    vtx.RotateUz(direction);
  } 

  // *** For v+A events:
  // Use const probability per unit volume within the sphere
  // (don't do the silly mistake to generate vertices uniformly in R)
  //

  else {
    while(vtx.Mag() > fNuclRadius) {
      vtx.SetX(-fNuclRadius + 2*fNuclRadius * rnd->RndFsi().Rndm());
      vtx.SetY(-fNuclRadius + 2*fNuclRadius * rnd->RndFsi().Rndm());
      vtx.SetZ(-fNuclRadius + 2*fNuclRadius * rnd->RndFsi().Rndm());
    }
  }

  LOG("Intranuke", pNOTICE) 
     << "Generated vtx @ R = " << vtx.Mag() << " fm / " 
     << print::Vec3AsString(&vtx);

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

  if(fInTestMode) {
    // in test mode (initial state: hadron+nucleus) rescatter the initial
    // state hadron
    return (p->Status() == kIStInitialState && !pdg::IsIon(p->Pdg()));
  }
  else {
   // attempt to rescatter anything marked as 'hadron in the nucleus'
   return (p->Status() == kIStHadronInTheNucleus);
  }
}
//___________________________________________________________________________
bool Intranuke::CanRescatter(const GHepParticle * p) const
{
// checks whether a particle that needs to be rescattered, can in fact be 
// rescattered by this cascade MC

  assert(p);
  return  ( p->Pdg() == kPdgPiP     || 
            p->Pdg() == kPdgPiM     || 
            p->Pdg() == kPdgPi0     ||
            p->Pdg() == kPdgProton  ||
            p->Pdg() == kPdgNeutron
          );
}
//___________________________________________________________________________
bool Intranuke::IsInNucleus(const GHepParticle * p) const
{
// check whether the input particle is still within the nucleus
//
  return (p->X4()->Vect().Mag() < fNuclRadius);
}
//___________________________________________________________________________
void Intranuke::TransportHadrons(GHepRecord * evrec) const
{
// transport all hadrons outside the nucleus

  // Loop over GHEP and run intranuclear rescattering on handled particles
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  int icurr = -1;

  while( (p = (GHepParticle *) piter.Next()) )
  {
    icurr++;

    // Check whether the particle needs rescattering, otherwise skip it
    if( ! this->NeedsRescattering(p) ) continue;

    LOG("Intranuke", pNOTICE)
      << " >> Stepping a " << p->Name() 
                        << " with kinetic E = " << p->KinE() << " GeV";

    // Rescatter a clone, not the original particle
    GHepParticle * sp = new GHepParticle(*p); 

    // Set clone's mom to be the hadron that was cloned
    sp->SetFirstMother(icurr); 

    // Check whether the particle can be rescattered 
    if(!this->CanRescatter(sp)) {

       // if I can't rescatter it, I will just take it out of the nucleus
       LOG("Intranuke", pNOTICE)
              << "... Current version can't rescatter a " << sp->Name();
       sp->SetFirstMother(icurr); 
       sp->SetStatus(kIStStableFinalState);
       evrec->AddParticle(*sp);
       continue; // <-- skip to next GHEP entry
    }

    // ** Special treatment for hadrons handled for the very first time
    // ** 

    // All 'fresh' hadrons would be advanced by the formation zone. Then
    // their corresponding (cloned) GHEP entry would be stored and would
    // move on to the next particle...

    bool fresh = this->IsFreshHadron(evrec,sp);
    if (fresh) {
       this->AdvanceFreshHadron(evrec,sp);

       // Check whether the hadron's position is still within the nucleus
       if(!this->IsInNucleus(sp)) {
           LOG("Intranuke", pNOTICE)
                  << "*** Hadron escaped the nucleus! Done with it.";
           sp->SetStatus(kIStStableFinalState);
       }

       evrec->AddParticle(*sp); // <-- add it to the GHEP record
       continue;                // <-- skip to next GHEP entry
    }

    // ** Standard intranuclear rescattering step
    // ** 

    // Generate a step and advance the hadron by it
    double d = this->GenerateStep(evrec,sp);
    this->StepParticle(sp, d);

    // Check whether it interacts
    bool interacts = this->IsInNucleus(sp);

    // If it does not interact, it must exit the nucleus. Done with it! 
    if(!interacts) {
       LOG("Intranuke", pNOTICE) 
                    << "*** Hadron escaped the nucleus! Done with it.";
       sp->SetStatus(kIStStableFinalState);
       evrec->AddParticle(*sp);
       continue; 
    }

    // The stepped particle interacts. Simulate the hadronic interaction.
    this->SimHadroProc(evrec,sp);

    // Current snapshot
    //LOG("Intranuke", pINFO) << "Current event record snapshot: " << *evrec;

  }// GHEP entries
}
//___________________________________________________________________________
bool Intranuke::IsFreshHadron(GHepRecord* evrec, GHepParticle* p) const
{
// Decide whether the particle p is a 'fresh' hadron (direct descendant of 
// the 'HadronicSystem' GHEP entry -in case od the KNO model- or descedant of 
// the JETSET special particles -cluster,string,indep-)

  // in test mode the initial state is hadron+nucleus and that initial
  // hadron does not have a formation zone
  if(fInTestMode) return false;

  // mom
  int imom = p->FirstMother();
  assert(imom>0);

  // grand-mom
  int igmom = evrec->Particle(imom)->FirstMother();
  assert(igmom>0);

  // grand-mom pdgc
  int gmom_pdg = evrec->Particle(igmom)->Pdg();

  if (gmom_pdg == kPdgHadronicSyst ||
      gmom_pdg == kPdgCluster      ||
      gmom_pdg == kPdgString       ||
      gmom_pdg == kPdgIndep) return true;

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

  LOG("Intranuke", pNOTICE)
    << "Hadron handled for first time is advanced by the formation zone = "
    << fzone << " m";
}
//___________________________________________________________________________
double Intranuke::GenerateStep(GHepRecord* evrec, GHepParticle* p) const
{
// Generate a step (in fermis) for particle p in the input event.
// Computes the mean free path L and generate an 'interaction' distance d 
// from an exp(-d/L) distribution

  RandomGen * rnd = RandomGen::Instance();

  double L = this->MeanFreePath(evrec, p);
  double d = -1.*L * TMath::Log(rnd->RndFsi().Rndm());

  LOG("Intranuke", pNOTICE)
            << "Mean free path = " << L << " fm / "
                              << "Generated path length = " << d << " fm";
  return d;
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
double Intranuke::MeanFreePath(GHepRecord* evrec, GHepParticle* p) const
{
// [adapted from neugen3 mean_free_path.F]
//
// Mean free path for the pions with the input kinetic energy.
// Determines the scattering length lamda from the SAID pi+N PWA fit ands
// scales it to the nuclear density.
//
// output mean free path units: fm

  // compute pion kinetic energy K in MeV
  double K = (p->KinE() / units::MeV);

  // get the nuclear target mass number
  GHepParticle * nucltgt = evrec->TargetNucleus();
  double A = nucltgt->A();

  // get current radial position within the nucleus
  double rnow = p->X4()->Vect().Mag();

  // get the nuclear density at the current position
  double rrms = (A>=20) ? 3.76 /* heavy:Fe */ : 2.44 /* light:C */;
  double req  = rrms*1.26;
  double req3 = TMath::Power(req,3);
  double rho  = A/(4*kPi*req3/3.);

  if(A>20) rho = A * this->DensityWoodsSaxon(rnow);
  else     rho = A * this->DensityGaus(rnow);

  // get total xsection for the incident hadron at its current 
  // kinetic energy
  int    pdgc   = p->Pdg();
  double sigtot = 0;

  if (pdgc == kPdgPiP) 
      sigtot = fHadroData -> XSecPipN_Tot() -> Evaluate(K);
  else if (pdgc == kPdgPi0) 
      sigtot = fHadroData -> XSecPi0N_Tot() -> Evaluate(K);
  else if (pdgc == kPdgPiM) 
      sigtot = fHadroData -> XSecPimN_Tot() -> Evaluate(K);
  else if (pdgc == kPdgProton) 
      sigtot = fHadroData -> XSecPN_Tot()   -> Evaluate(K);
  else if (pdgc == kPdgNeutron) 
      sigtot = fHadroData -> XSecNN_Tot()   -> Evaluate(K);
  else {
     LOG("Intranuke", pWARN)
         << "Can't compute mean free path for particle code = " << pdgc;
     return 0;
  }

  // the xsection splines in INukeHadroData return the hadron x-section in
  // mb -> convert to fm^2
  sigtot *= (units::mb / units::fm2);
  LOG("Intranuke", pDEBUG) 
       << "SigmaTotal(KineticE = " << K << " MeV) = " << sigtot << " fm2";

  // compute the mean free path
  double lamda = 1. / (rho * sigtot);
  LOG("Intranuke", pDEBUG) << "Mean free path = " << lamda << " fm";

  return lamda;
}
//___________________________________________________________________________
double Intranuke::DensityGaus(double r) const
{
// [adapted from neugen3 density_gaus.F]
//
// Modified harmonic osc density distribution. Norm gives normalization to 1
// Only C12 here - can be expalned to others later
//
// input  : radial distance in nucleus [units: fm]
// output : nuclear density            [units: fm^-3]

  double norm = 0.0132;
  double a    = 1.635;
  double alf  = 1.403;
  double b    = TMath::Power(r/a,2.);
  double dens = norm * (1. + alf*b) * TMath::Exp(-b);
  return dens;
}
//___________________________________________________________________________
double Intranuke::DensityWoodsSaxon(double r) const
{
// [adapted from neugen3 density_ws.F]
//
// Woods-Saxon desity distribution. Norn gives normalization to 1
// Onlt Fe here - can be expalned to others later
//
// input  : radial distance in nucleus [units: fm]
// output : nuclear density            [units: fm^-3]

  double norm = 0.002925;
  double c    = 4.10;
  double z    = 0.56;
  double dens = norm / (1 + TMath::Exp((r-c)/z));
  return dens;
}
//___________________________________________________________________________
void Intranuke::SimHadroProc(GHepRecord* ev, GHepParticle* p) const
{
// Simulate a hadron interaction for the input particle p
//
  switch(fMode) {
  case (kIMdUndefined) :
          LOG("Intranuke", pERROR) << "Undefined INTRANUKE mode";
          break;
  case (kIMdHN) :
          // simulate the hadron interaction in INTRANUKE's hN mode
          this->SimHadroProcHN(ev,p);
          break;
  case (kIMdHA) :
          // simulate the hadron interaction in INTRANUKE's hA mode
          this->SimHadroProcHA(ev,p);
          break;
  default :
          LOG("Intranuke", pERROR) << "Undefined INTRANUKE mode";
          break;
  }
}
//___________________________________________________________________________
//___________________________________________________________________________
// Methods specific to INTRANUKE's HA-mode
//___________________________________________________________________________
//___________________________________________________________________________
void Intranuke::SimHadroProcHA(GHepRecord* ev, GHepParticle* p) const
{
// Simulate a hadron interaction for the input particle p in HA mode
//
  // check inputs
  if(!p || !ev) {
     LOG("Intranuke", pERROR) << "** Null input!";
     return;
  }

  // check particle id
  int  pdgc = p->Pdg();
  bool is_pion    = (pdgc==kPdgPiP || pdgc==kPdgPiM || pdgc==kPdgPi0);
  bool is_baryon  = (pdgc==kPdgProton || pdgc==kPdgNeutron);
  bool is_handled = (is_baryon || is_pion);
  if(!is_handled) {
     LOG("Intranuke", pERROR) << "** Can not handle particle: " << p->Name();
     return;     
  }

  // select a fate for the input particle
  INukeFateHA_t fate = this->HadronFateHA(p);
  if(fate == kIHAFtUndefined) {
     LOG("Intranuke", pERROR) << "** Couldn't select a fate";
     p->SetStatus(kIStStableFinalState);
     ev->AddParticle(*p);
     return;     
  }
  LOG("Intranuke", pNOTICE)
    << "Selected " << p->Name() << " fate: " << INukeHadroFates::AsString(fate);

  // get the reaction products for the selected fate
  if (fate == kIHAFtCEx || fate == kIHAFtElas  || fate == kIHAFtInelas) {
   if (is_pion)   this->PiSlam(ev,p,fate);
   if (is_baryon) this->PnSlam(ev,p,fate);
  }
  else if (fate == kIHAFtAbsNP   ||
           fate == kIHAFtAbsPP   ||
           fate == kIHAFtAbsNPP  ||
           fate == kIHAFtAbsNNP  ||
           fate == kIHAFtAbs2N2P ||
           fate == kIHAFtAbs2N3P ||
           fate == kIHAFtNPip    || 
           fate == kIHAFtNPipPi0) {
    this->Inelastic(ev,p,fate);
  }
}
//___________________________________________________________________________
INukeFateHA_t Intranuke::HadronFateHA(const GHepParticle * p) const
{
// Select a hadron fate in HA mode
//
  RandomGen * rnd = RandomGen::Instance();

  // get pdgc code & kinetic energy in MeV
  int    pdgc = p->Pdg();
  double ke   = p->KinE() / units::MeV;
 
  // try to generate a hadron fate
  unsigned int iter = 0;
  while(iter++ < kRjMaxIterations) {

    // handle pions
    if (pdgc==kPdgPiP || pdgc==kPdgPiM || pdgc==kPdgPi0) {
       double tot_frac = 0;
       double r = rnd->RndFsi().Rndm();
       LOG("Intranuke", pDEBUG) << "r = " << r;

       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtCEx,     ke))) return kIHAFtCEx;     // cex
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtElas,    ke))) return kIHAFtElas;    // elas
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtInelas,  ke))) return kIHAFtInelas;  // inelas
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtAbsNP,   ke))) return kIHAFtAbsNP;   // abs np
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtAbsPP,   ke))) return kIHAFtAbsPP;   // abs pp
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtAbsNPP,  ke))) return kIHAFtAbsNPP;  // abs npp
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtAbsNNP,  ke))) return kIHAFtAbsNNP;  // abs nnp
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtAbs2N2P, ke))) return kIHAFtAbs2N2P; // abs 2n2p 
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtNPipPi0, ke))) return kIHAFtNPipPi0; // pi production : n pi+ pi0

       LOG("Intranuke", pWARN) 
         << "No selection after going through all fates! "
         << "Total fraction = " << tot_frac << " (r = " << r << ")";
    }

    // handle nucleons
    else if (pdgc==kPdgProton || pdgc==kPdgNeutron) {
       double tot_frac = 0;
       double r = rnd->RndFsi().Rndm();
       LOG("Intranuke", pDEBUG) << "r = " << r;

       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtCEx,     ke))) return kIHAFtCEx;     // cex
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtElas,    ke))) return kIHAFtElas;    // elas
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtInelas,  ke))) return kIHAFtInelas;  // inelas
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtAbsNP,   ke))) return kIHAFtAbsNP;   // abs np
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtAbsPP,   ke))) return kIHAFtAbsPP;   // abs pp
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtAbsNPP,  ke))) return kIHAFtAbsNPP;  // abs npp
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtAbsNNP,  ke))) return kIHAFtAbsNNP;  // abs nnp
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtAbs2N3P, ke))) return kIHAFtAbs2N3P; // abs 2n3p
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtNPip,    ke))) return kIHAFtNPip;    // pi production : n pi+
       if(r < (tot_frac += fHadroData->Frac(pdgc, kIHAFtNPipPi0, ke))) return kIHAFtNPipPi0; // pi production : n pi+ pi0

       LOG("Intranuke", pWARN) 
         << "No selection after going through all fates! "
         << "Total fraction = " << tot_frac << " (r = " << r << ")";
    }
  }//iterations

  return kIHAFtUndefined; 
}
//___________________________________________________________________________
void Intranuke::PiSlam(
         GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const
{
// [adapted from neugen3 intranuke_pislam.F]
//
// Scatters charged pions on iron nuclei.
// Note: At pi momenta below a few GeV, beta for the Lab to CM tranformation
// is small due to the huge mass of an iron nucleus. Here, the Lab and CM 
// are, for our purposes, identical.
//
  LOG("Intranuke", pDEBUG) 
      << "PiSlam() is invoked for a : " << p->Name() 
      << " whose fate is : " << INukeHadroFates::AsString(fate);

  if (fate!=kIHAFtCEx && fate!=kIHAFtElas && fate!=kIHAFtInelas) {
     LOG("Intranuke", pWARN) 
         << "Can not handle fate: " << INukeHadroFates::AsString(fate);
     return;
  }
  RandomGen * rnd = RandomGen::Instance();

  // compute the scattering angle
  //
  double costheta = 0.; 
  if(fate == kIHAFtElas) {  
    // Elastic scattering
    double theta = this->PiBounce();
    costheta = TMath::Cos(theta);
  } else {   
    // CEx & inelastic scattering
    costheta = 2.0 * rnd->RndFsi().Rndm() - 1.0;
  }
  double sintheta = TMath::Sqrt(1.0 - TMath::Power(costheta,2));

  // degrade the momentum for CEx and Inelastic scattering
  //
  double ptot = 0;
  if (fate == kIHAFtCEx || fate == kIHAFtInelas) {
     double ke     = p->KinE();
     double ke_fin = 0.;
     if(costheta < 0.) { ke_fin = this->Degrade(ke); }
     else              { ke_fin = ke * (0.5 * rnd->RndFsi().Rndm() + 0.5); }

     // handle pdg code change for CEx
     if (fate == kIHAFtCEx) {
       int pdgc = p->Pdg();
       if(pdgc != kPdgPi0) pdgc = kPdgPi0; // pi0
       else {
         double r = rnd->RndFsi().Rndm();
         pdgc = (r > 0.5) ? kPdgPiP : kPdgPiM; // 50% pi+ - 50% pi-
       }
       p->SetPdgCode(pdgc);
     }
     ptot = TMath::Sqrt(ke_fin*ke_fin + 2.0 * p->Mass() * ke_fin);
  }
  else { ptot = p->P4()->Mag(); }

  // update pi momentum & energy
  //
  double pz  = ptot * costheta;
  double pt  = ptot * sintheta;
  double phi = 2 * kPi * rnd->RndFsi().Rndm();

  TVector3 plab(0,pt,pz);
  TVector3 dir = p->P4()->Vect().Unit(); // unit vector along incoming pi

  plab.RotateUz(dir);  // align pz with the direction of the incoming pi
  plab.RotateZ(phi);   // randomize tranverse components

  double ptot2  = TMath::Power(ptot,      2);
  double pmass2 = TMath::Power(p->Mass(), 2);
  double energy = TMath::Sqrt(ptot2 + pmass2);

  TLorentzVector p4lab(plab, energy);
  p->SetMomentum(p4lab);

  // now add the modified particle to the GHEP record
  //
  p->SetStatus(kIStStableFinalState); // done with it in hA?
  ev->AddParticle(*p);
}
//___________________________________________________________________________
void Intranuke::PnSlam(
         GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const
{
// [adapted from neugen3 intranuke_pnslam.F]
//
// Scatters charged protons off nuclei.
//
  LOG("Intranuke", pDEBUG) 
      << "PnSlam() is invoked for a : " << p->Name() 
      << " whose fate is : " << INukeHadroFates::AsString(fate);

  if (fate!=kIHAFtCEx && fate!=kIHAFtElas && fate!=kIHAFtInelas) {
     LOG("Intranuke", pWARN) 
         << "Can not handle fate: " << INukeHadroFates::AsString(fate);
     return;
  }

  RandomGen * rnd = RandomGen::Instance();

  // compute the scattering angle
  //
  double theta    = this->PnBounce();
  double costheta = TMath::Cos(theta);
  double sintheta = TMath::Sqrt(1.0 - TMath::Power(costheta,2));

  // handle pdg code change for CEx
  //
  if (fate == kIHAFtCEx) {
     int pdgc = p->Pdg();
     pdgc = pdg::SwitchProtonNeutron(pdgc);
     p->SetPdgCode(pdgc);
  }
     
  // update nucleon momentum & energy
  //
  double ptot = p->P4()->Mag();  
  double pz   = ptot * costheta;
  double pt   = ptot * sintheta;
  double phi  = 2 * kPi * rnd->RndFsi().Rndm();

  TVector3 plab(0,pt,pz);
  TVector3 dir = p->P4()->Vect().Unit(); // unit vector along incoming nucleon

  plab.RotateUz(dir);  // align pz with the direction of the incoming nucleon
  plab.RotateZ(phi);   // randomize tranverse components

  double ptot2  = TMath::Power(ptot,      2);
  double pmass2 = TMath::Power(p->Mass(), 2);
  double energy = TMath::Sqrt(ptot2 + pmass2);

  TLorentzVector p4lab(plab, energy);
  p->SetMomentum(p4lab);

  // now add the modified particle to the GHEP record
  //
  p->SetStatus(kIStStableFinalState); // done with it in hA?
  ev->AddParticle(*p);
}
//___________________________________________________________________________
double Intranuke::PiBounce(void) const
{
// [adapted from neugen3 intranuke_bounce.F]
// [is a fortran stub / difficult to understand - needs to be improved]
//
// Generates theta in radians for elastic pion-nucleus scattering/
// Lookup table is based on Fig 17 of Freeman, Miller and Henley, Nucl.Phys.
// A389, 457 (1982)
//
  const int nprob = 25;
  double dintor = 0.0174533;
  double denom  = 47979.453;
  double rprob[nprob] = {
    5000., 4200., 3000., 2600., 2100., 1800., 1200., 750., 500., 230., 120., 
    35., 9., 3., 11., 18., 29., 27., 20., 14., 10., 6., 2., 0.14, 0.19 };

  double angles[nprob];
  for(int i=0; i<nprob; i++) angles[i] = 2.5*i;

  RandomGen * rnd = RandomGen::Instance();
  double r = rnd->RndFsi().Rndm();
  
  double xsum  = 0.;
  double theta = 0.;
  double binl  = 0.;
  double binh  = 0.;
  int tj = 0;
  for(int i=0; i<60; i++) {
   theta = i+0.5;
   for(int j=0; j < nprob-1; j++) {
     binl = angles[j];
     binh = angles[j+1];
     tj=j;
     if(binl<=theta && binh>=theta) break;
     tj=0;
   }//j
   int itj = tj;
   double tfract = (theta-binl)/2.5;
   double delp   = rprob[itj+1] - rprob[itj];
   xsum += (rprob[itj] + tfract*delp)/denom;
   if(xsum>r) break;
   theta = 0.;
  }//i

  theta *= dintor;

  LOG("Intranuke", pINFO) 
     << "Generated pi+A elastic scattering angle = " << theta << " radians";

  return theta;
}
//___________________________________________________________________________
double Intranuke::PnBounce(void) const
{
// [adapted from neugen3 intranuke_pnbounce.F]
// [is a fortran stub / difficult to understand - needs to be improved]
//
// Generates theta in radians for elastic nucleon-nucleus scattering.
// Use 800 MeV p+O16 as template in same (highly simplified) spirit as pi+A
// from table in Adams et al., PRL 1979. Guess value at 0-2 deg based on Ni
// data.
//
  const int nprob = 20;
  double dintor = 0.0174533;
  double denom  = 11967.0;
  double rprob[nprob] = {
    2400., 2350., 2200., 2000., 1728., 1261., 713., 312., 106., 35., 
    6., 5., 10., 12., 11., 9., 6., 1., 1., 1. };

  double angles[nprob];
  for(int i=0; i<nprob; i++) angles[i] = 1.0*i;

  RandomGen * rnd = RandomGen::Instance();
  double r = rnd->RndFsi().Rndm();

  double xsum  = 0.;
  double theta = 0.;
  double binl  = 0.;
  double binh  = 0.;
  int tj = 0;
  for(int i=0; i < nprob; i++) {
    theta = i + 0.5;
    for(int j=0; j < nprob-1; j++) {
      binl = angles[j];
      binh = angles[j+1];
      tj = j;
      if(binl<theta && binh>theta) break;
      tj = 0;
    } //j 
    double tfract = (theta-binl)/1.;
    int itj = tj;
    double delp = rprob[itj+1] - rprob[itj];
    xsum += (rprob[itj] + tfract*delp)/denom;
    if(xsum > r) break;
    theta = 0.0;
  } //i

  theta *= dintor;

  LOG("Intranuke", pINFO) 
     << "Generated N+A elastic scattering angle = " << theta << " radians";

  return theta;
}
//___________________________________________________________________________
double Intranuke::Degrade(double ke) const
{
// [adapted from neugen3 intranuke_degrade.F]
// [is a fortran stub / difficult to understand - needs to be improved]
//
// The method receives the incident pion kinetic energy for inelastic scatt.
// process and returns the degraded final-state kinetic energy.
// Is based on inelastic scattering distribution Fig.18, of C.H.Q.Ingram,
// Nucl.Phys.A374, 319 (1982)
//
  const int np = 15;
  double efract[np] = { 
   0.13, 0.19, 0.25, 0.31, 0.38, 0.44, 0.50, 0.56, 0.63, 0.69, 0.75, 0.81,
   0.88, 0.94, 1.00 };
  double rprob[np] = {
   0.00, 0.33, 0.54, 0.67, 0.72, 0.76, 0.71, 0.59, 0.41, 0.32, 0.25, 0.19, 
   0.11, 0.08, 0.02 };
  double ekin[np] = {
    20.,  30. , 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 
    140., 150., 160. };
  double distn = 56.95;

  RandomGen * rnd = RandomGen::Instance();
  double r = rnd->RndFsi().Rndm();

  double xsum   = 0.;
  double fractn = 0.;
  for(int i=0; i<70; i++) {
    double x = 1. + i;
    double tmev = 20. + 2.*x - 1.;
    int ibinl = -1;
    for(int j=0; j<np; j++) {
       ibinl = j;
       if(ekin[j]>tmev) break;
    }//j
    if(ibinl>0) {
       int ibinh = ibinl+1;
       double delt   = (tmev-ekin[ibinl])/10.0;
       double dprob  = rprob[ibinh]  - rprob[ibinl];
       double dfract = efract[ibinh] - efract[ibinl];
       double probi  = rprob[ibinl] + delt * dprob;
       xsum += (2.*probi)/distn;
       if(xsum>r) {
         fractn = efract[ibinl] + delt * dfract;
         break;
       }
    } else {
       fractn = 0.;
       break;
    }    
  }//i

  double ke_fin = fractn * ke;

  LOG("Intranuke", pINFO) 
     << "Degrading pion kinetic energy for pi+N inelastic scattering : " 
     << ke << " -> " << ke_fin;

  return ke_fin;
}
//___________________________________________________________________________
void Intranuke::Inelastic(
          GHepRecord* ev, GHepParticle* p, INukeFateHA_t fate) const
{
// [adapted from neugen3 intranuke_pi_inelastic.F, intranuke_pn_inelastic.F] 
//

  LOG("Intranuke", pDEBUG) 
      << "Inelastic() is invoked for a : " << p->Name() 
      << " whose fate is : " << INukeHadroFates::AsString(fate);

  bool allow_dup = true;
  PDGCodeList list(allow_dup); // list of final state particles

  // figure out the final state according to the fate
  //
  switch (fate) {
   case (kIHAFtAbsNP)   : // -> n+p
         list.push_back(kPdgNeutron);
         list.push_back(kPdgProton);
         break;
   case (kIHAFtAbsPP)   : // -> p+p
         list.push_back(kPdgProton);
         list.push_back(kPdgProton);
         break;
   case (kIHAFtAbsNPP)  : // -> n+p+p
         list.push_back(kPdgNeutron);
         list.push_back(kPdgProton);
         list.push_back(kPdgProton);
         break;
   case (kIHAFtAbsNNP)  : // -> n+n+p
         list.push_back(kPdgNeutron);
         list.push_back(kPdgNeutron);
         list.push_back(kPdgProton);
         break;
   case (kIHAFtAbs2N2P) : // -> n+n+p+p
         list.push_back(kPdgNeutron);
         list.push_back(kPdgNeutron);
         list.push_back(kPdgProton);
         list.push_back(kPdgProton);
         break;
   case (kIHAFtAbs2N3P) : // -> n+n+p+p+p
         list.push_back(kPdgNeutron);
         list.push_back(kPdgNeutron);
         list.push_back(kPdgProton);
         list.push_back(kPdgProton);
         list.push_back(kPdgProton);
         break;
   case (kIHAFtNPip)  : // -> n + pi+ 
         list.push_back(kPdgNeutron);
         list.push_back(kPdgPiP);
         break;
   case (kIHAFtNPipPi0)  : // -> n + pi+ + pi0 
         list.push_back(kPdgNeutron);
         list.push_back(kPdgPiP);
         list.push_back(kPdgPi0);
         break;
   default : 
         LOG("Intranuke", pWARN) 
             << "Can not handle fate: " << INukeHadroFates::AsString(fate);
         return;
  }

  // do the phase space decay & save all f/s particles to the event record
  //
  this->PhaseSpaceDecay(ev,p,list);
}
//___________________________________________________________________________
//___________________________________________________________________________
// Methods specific to INTRANUKE's HN-mode
//___________________________________________________________________________
//___________________________________________________________________________
void Intranuke::SimHadroProcHN(GHepRecord* /*ev*/, GHepParticle* /*p*/) const
{
// Simulate a hadron interaction for the input particle p in HN mode
//
  LOG("Intranuke", pERROR) << "** HN-mode not implemented yet";
}
//___________________________________________________________________________
// Generic Phase Space Decay method
//___________________________________________________________________________
bool Intranuke::PhaseSpaceDecay(
             GHepRecord* ev, GHepParticle* p, const PDGCodeList & pdgv) const
{
// General method decaying the input particle system 'pdgv' with available 4-p
// given by 'pd'. The decayed system is used to populate the input TMCParticle
// array starting from the slot 'offset'.
//
  LOG("Intranuke", pINFO) << "*** Performing a Phase Space Decay";
  assert(pdgv.size() > 1);

  // Get the decay product masses & names

  ostringstream state_sstream;
  state_sstream << "( ";
  vector<int>::const_iterator pdg_iter;
  int i = 0;
  double * mass = new double[pdgv.size()];
  double   mass_sum = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    double m  = PDGLibrary::Instance()->Find(pdgc)->Mass();
    string nm = PDGLibrary::Instance()->Find(pdgc)->GetName();
    mass[i++] = m;
    mass_sum += m;
    state_sstream << nm << " ";
  }
  state_sstream << ")";

  TLorentzVector * pd = p->GetP4(); // incident particle 4p

  // update available energy -> init (mass + kinetic) + sum of f/s masses
  double availE = pd->Energy() + mass_sum; 
  availE += 0.01; // steve's recipe
  pd->SetE(availE);

  LOG("Intranuke", pINFO)
    << "Final state = " << state_sstream.str() << " has N = " << pdgv.size() 
    << " particles / total mass = " << mass_sum;
  LOG("Intranuke", pINFO)
    << "Decaying system p4 = " << utils::print::P4AsString(pd);

  // Set the decay
  bool permitted = fGenPhaseSpace.SetDecay(*pd, pdgv.size(), mass);
  if(!permitted) {
     LOG("Intranuke", pERROR)
       << " *** Phase space decay is not permitted \n"
       << " Total particle mass = " << mass_sum << "\n"
       << " Decaying system p4 = " << utils::print::P4AsString(pd);

     // clean-up and return
     delete [] mass;
     delete pd;
     return false;
  }

  // The decay is permitted - add the incident particle at the event record
  // and mark is as 'decayed state'
  p->SetStatus(kIStDecayedState);
  ev->AddParticle(*p);

  // Get the maximum weight
  double wmax = -1;
  for(int i=0; i<200; i++) {
     double w = fGenPhaseSpace.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);

  LOG("Intranuke", pNOTICE)
   << "Max phase space gen. weight @ current hadronic interaction: " << wmax;

  // Generate an unweighted decay

  RandomGen * rnd = RandomGen::Instance();
  wmax *= 1.2;

  bool accept_decay=false;
  unsigned int itry=0;

  while(!accept_decay)
  {
    itry++;

    if(itry>kMaxUnweightDecayIterations) {
       // report, clean-up and return
       LOG("Intranuke", pWARN)
             << "Couldn't generate an unweighted phase space decay after "
             << itry << " attempts";
       delete [] mass;
       delete pd;
       return false;
    }

    double w  = fGenPhaseSpace.Generate();
    double gw = wmax * rnd->RndFsi().Rndm();

    LOG("Intranuke", pINFO) << "Decay weight = " << w << " / R = " << gw;
    accept_decay = (gw<=w);
  }

  // Insert final state products into the event record
  // - the particles are added as daughters of the decayed state
  // - the particles are marked as final stable state (in hA mode)
  i=0;
  int mom = ev->ParticlePosition(p);
  GHepStatus_t ist = kIStStableFinalState;

  TLorentzVector * v4 = p->GetX4();

  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {

     //-- current PDG code
     int pdgc = *pdg_iter;

     //-- get the 4-momentum of the i-th final state particle
     TLorentzVector * p4fin = fGenPhaseSpace.GetDecay(i++);

     //-- add the particle at the event record
     //   (if it is a cascade nucleon set the binding energy to be 
     //    removed later on)
    
     GHepParticle new_particle(pdgc, ist, mom,-1,-1,-1, *p4fin, *v4);
     if(pdg::IsNeutronOrProton(pdgc)) 
             new_particle.SetRemovalEnergy(fNucRmvE);

     ev->AddParticle(new_particle);
  }
  // Clean-up
  delete [] mass;
  delete pd;
  delete v4;

  return true;
}
//___________________________________________________________________________
//___________________________________________________________________________
// Methods for configuring the INTRANUKE module:
//___________________________________________________________________________
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

  //-- load hadronic cross sections
  fHadroData = INukeHadroData::Instance();

  //-- intranuke mode (h+N or h+A)
  string mode = fConfig->GetStringDef ("mode", gc->GetString("INUKE-Mode"));
  if      (mode == "hA" || mode == "HA") fMode = kIMdHA;
  else if (mode == "hN" || mode == "HN") fMode = kIMdHN;
  else                                   fMode = kIMdUndefined;

  //-- in test mode? (def:no)
  fInTestMode = fConfig->GetBoolDef ("test-mode", false);

  //-- other intranuke config params
  fct0     = fConfig->GetDoubleDef ("ct0",  gc->GetDouble("INUKE-FormationZone")); // fm
  fK       = fConfig->GetDoubleDef ("Kpt2", gc->GetDouble("INUKE-KPt2"));
  fR0      = fConfig->GetDoubleDef ("R0",   gc->GetDouble("INUKE-Ro")); // fm
  fNucRmvE = fConfig->GetDoubleDef ("nucleon-removal-energy", gc->GetDouble("INUKE-NucRemovalE")); // GeV

  //-- report
  LOG("Intranuke", pDEBUG) << "mode    = " << INukeMode::AsString(fMode);
  LOG("Intranuke", pDEBUG) << "ct0     = " << fct0 << " fermi";
  LOG("Intranuke", pDEBUG) << "K(pt^2) = " << fK;
  LOG("Intranuke", pDEBUG) << "R0      = " << fR0  << " fermi";
}
//___________________________________________________________________________


