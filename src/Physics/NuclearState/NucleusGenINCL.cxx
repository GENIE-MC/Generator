//____________________________________________________________________________
/*!

\class    genie::NucleusGenINCL

\brief    It visits the event record & computes a Fermi motion momentum for
          initial state nucleons bound in nuclei.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Liang Liu <liangliu \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  October 17, 2024

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__
#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Physics/NuclearState/NucleusGenINCL.h"

#include "Physics/NuclearState/NuclearModel.h"
#include "Physics/NuclearState/NuclearModelI.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/NuclearState/FermiMomentumTablePool.h"
#include "Physics/NuclearState/FermiMomentumTable.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NuclearState/INCLNucleus.h"

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
NucleusGenINCL::NucleusGenINCL() :
EventRecordVisitorI("genie::NucleusGenINCL")
{

}
//___________________________________________________________________________
NucleusGenINCL::NucleusGenINCL(string config) :
EventRecordVisitorI("genie::NucleusGenINCL", config)
{

}
//___________________________________________________________________________
NucleusGenINCL::~NucleusGenINCL()
{

}

//___________________________________________________________________________
void NucleusGenINCL::ProcessEventRecord(GHepRecord * evrec) const
{
  // skip if not a nuclear target
  if(! evrec->Summary()->InitState().Tgt().IsNucleus()) return;

  // skip if no hit nucleon is set
  if(! evrec->HitNucleon()) return;

  LOG("NucleusGenINCL", pINFO) << "Adding final state nucleus";
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  LOG("NucleusGenINCL", pINFO) << "Adding final state nucleus";
  incl_nucleus->initialize(evrec);
  incl_nucleus->reset(evrec);
  incl_nucleus->initialize(evrec);
  LOG("NucleusGenINCL", pINFO) << "Adding final state nucleus";

  // give hit nucleon a vertex
  this->setInitialStateVertex(evrec);
  LOG("NucleusGenINCL", pINFO) << "Adding final state nucleus";
  // give hit nucleon a Fermi momentum
  this->setInitialStateMomentum(evrec);
  LOG("NucleusGenINCL", pINFO) << "Adding final state nucleus";

  // handle the addition of the recoil nucleon
  // TODO:  INCL has it own SRC model
//  if ( fSecondEmitter ) fSecondEmitter -> ProcessEventRecord( evrec ) ;

  // add a recoiled nucleus remnant
  this->setTargetNucleusRemnant(evrec);
  LOG("NucleusGenINCL", pINFO) << "Adding final state nucleus";
}

//___________________________________________________________________________
//  using INCL model to get the position and momentum of 
//  Hit  nucleon
void NucleusGenINCL::setInitialStateVertex(GHepRecord * evrec) const{

  INCLNucleus *incl_nucleus = INCLNucleus::Instance();

// generate a vtx and set it to all GHEP physical particles
  Interaction * interaction = evrec->Summary();
  GHepParticle * nucltgt = evrec->TargetNucleus();
  TVector3 vtx(9999999.,999999.,999999.);
  if(!nucltgt){
    vtx.SetXYZ(0.,0.,0.);
  }else{
    double A = nucltgt->A();
 //   vtx = GenerateVertex(interaction,A);
 	
  const ProcessInfo & proc_info = interaction->ProcInfo();
  bool is_coh = proc_info.IsCoherentProduction() || proc_info.IsCoherentElastic();
  bool is_ve  = proc_info.IsInverseMuDecay() ||
    proc_info.IsIMDAnnihilation() ||
    proc_info.IsNuElectronElastic() ||
    proc_info.IsGlashowResonance() ||
    proc_info.IsPhotonResonance() ||
    proc_info.IsPhotonCoherent();


  if(is_coh||is_ve) {
    // ** For COH or ve- set a vertex positon on the nuclear boundary
    //
  //  LOG("Vtx", pINFO)  << "Setting vertex on the nuclear boundary";
  //  double phi      = 2*kPi * rnd->RndFsi().Rndm();
  //  double cosphi   = TMath::Cos(phi);
  //  double sinphi   = TMath::Sin(phi);
  //  double costheta = -1 + 2 * rnd->RndFsi().Rndm();
  //  double sintheta = TMath::Sqrt(1-costheta*costheta);
  //  vtx.SetX(R*sintheta*cosphi);
  //  vtx.SetY(R*sintheta*sinphi);
  //  vtx.SetZ(R*costheta);
  }
  else {
    vtx = incl_nucleus->getHitNucleonPosition();
  }

  }

  LOG("NucleusGenINCL", pINFO) << "Position";
  vtx.Print();
  // Copy the vertex info to the particles already in the event  record
  //
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) )
  {
    if(pdg::IsPseudoParticle(p->Pdg())) continue;
    if(pdg::IsIon           (p->Pdg())) continue;

    LOG("NucleusGenINCL", pINFO) << "Setting vertex position for: " << p->Name();
    p->SetPosition(vtx.x(), vtx.y(), vtx.z(), 0.);
  }



}

void NucleusGenINCL::setInitialStateMomentum(GHepRecord * evrec) const{
  Interaction *  interaction = evrec       -> Summary();
  InitialState * init_state  = interaction -> InitStatePtr();
  Target *       tgt         = init_state  -> TgtPtr();

  // do nothing for non-nuclear targets
  if(!tgt->IsNucleus()) return;

  TLorentzVector * p4 = tgt->HitNucP4Ptr();

  // do nothing if the struct nucleon 4-momentum was set (eg as part of the
  // initial state selection)
  if(p4->Px()>0 || p4->Py()>0 || p4->Pz()>0) return;

  // access the hit nucleon and target nucleus at the GHEP record
  GHepParticle * nucleon = evrec->HitNucleon();
  GHepParticle * nucleus = evrec->TargetNucleus();
  assert(nucleon);
  assert(nucleus);

  // initialize INCL nucleus model
  // INCL nucleus model sample all nucleons with r-p correlation
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  G4INCL::Nucleus *incl_nuc = incl_nucleus->getNuclues();
  TLorentzVector p4tgt;
  p4tgt.SetPx(incl_nuc->getMomentum().getX() / 1000.);
  p4tgt.SetPy(incl_nuc->getMomentum().getY() / 1000. );
  p4tgt.SetPz(incl_nuc->getMomentum().getZ() / 1000. );
  p4tgt.SetE(incl_nuc->getEnergy() / 1000.);
  init_state->SetTgtP4(p4tgt);
  nucleus->SetMomentum(p4tgt);


  // get a random nucleon with respect to the isospin of evrec->HitNucleon();
  // the removal energy maybe not necessary
  TVector3 p3 = incl_nucleus->getHitNucleonMomentum();
  double   hit_nucleon_energy = incl_nucleus->getHitNucleonEnergy();
  double   w  = incl_nucleus->getRemovalEnergy();
  //-- update the struck nucleon 4p at the interaction summary and at
  // the GHEP record
  p4->SetPx(p3.Px()/1000.);
  p4->SetPy(p3.Py()/1000.);
  p4->SetPz(p3.Pz()/1000.);
  p4->SetE (hit_nucleon_energy/1000.);
  LOG("NucleusGenINCL", pINFO) << "Momentum";
  p4->Print();

  nucleon->SetMomentum(*p4);  // update GHEP value
  nucleon->SetRemovalEnergy(w);  // FIXME this may be not necessary


  // Sometimes, for interactions near threshold, Fermi momentum might bring
  // the neutrino energy in the nucleon rest frame below threshold (for the
  // selected interaction). In this case mark the event as unphysical and
  // abort the current thread.
  const KPhaseSpace & kps = interaction->PhaseSpace();
  if(!kps.IsAboveThreshold()) {
     LOG("NucleusGenINCL", pNOTICE)
                  << "Event below threshold after generating Fermi momentum";

     double Ethr = kps.Threshold();
     double Ev   = init_state->ProbeE(kRfHitNucRest);
     LOG("NucleusGenINCL", pNOTICE)
              << "Ev (@ nucleon rest frame) = " << Ev << ", Ethr = " << Ethr;

     evrec->EventFlags()->SetBitNumber(kBelowThrNRF, true);
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("E < Ethr after generating nucleon Fermi momentum");
     exception.SwitchOnFastForward();
     throw exception;
  }

}

void NucleusGenINCL::setTargetNucleusRemnant(GHepRecord * evrec)const{
// add the remnant nuclear target at the GHEP record

  LOG("NucleusGenINCL", pINFO) << "Adding final state nucleus";

  double Px = 0;
  double Py = 0;
  double Pz = 0;
  double E  = 0;

  GHepParticle * nucleus = evrec->TargetNucleus();
  int A = nucleus->A();
  int Z = nucleus->Z();

  int fd = nucleus->FirstDaughter();
  int ld = nucleus->LastDaughter();

  INCLNucleus *incl_nucleus = INCLNucleus::Instance();

  G4INCL::Nucleus *incl_nuc = incl_nucleus->getNuclues();


  for(int id = fd; id <= ld; id++) {

    // compute A,Z for final state nucleus & get its PDG code and its mass
    GHepParticle * particle = evrec->Particle(id);
    assert(particle);
    int  pdgc = particle->Pdg();
    bool is_p  = pdg::IsProton (pdgc);
    bool is_n  = pdg::IsNeutron(pdgc);

    if (is_p) Z--;
    if (is_p || is_n) A--;

    Px += particle->Px();
    Py += particle->Py();
    Pz += particle->Pz();
    E  += particle->E();

  }//daughters

  TParticlePDG * remn = 0;
  int ipdgc = pdg::IonPdgCode(A, Z);
  remn = PDGLibrary::Instance()->Find(ipdgc);
  if(!remn) {
    LOG("HadronicVtx", pFATAL)
          << "No particle with [A = " << A << ", Z = " << Z
                            << ", pdgc = " << ipdgc << "] in PDGLibrary!";
    assert(remn);
  }

  double Mi = nucleus->Mass();
  Px = incl_nuc->getMomentum().getX()/1000. - Px;
  Py = incl_nuc->getMomentum().getY()/1000. - Py;
  Pz = incl_nuc->getMomentum().getZ()/1000. - Pz;
  E = incl_nuc->getEnergy()/1000. - E;

  // Add the nucleus to the event record
  LOG("FermiMover", pINFO)
     << "Adding nucleus [A = " << A << ", Z = " << Z
     << ", pdgc = " << ipdgc << "]";

  int imom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
       ipdgc,kIStStableFinalState, imom,-1,-1,-1, Px,Py,Pz,E, 0,0,0,0);
}


//___________________________________________________________________________
void NucleusGenINCL::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGenINCL::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGenINCL::LoadConfig(void)
{

}
//____________________________________________________________________________


#endif // end  __GENIE_INCL_ENABLED__
