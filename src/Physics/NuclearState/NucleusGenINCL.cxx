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
#include <cstdlib>  // For getenv
#include <string>
#include <regex>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>
#include "G4INCLKinematicsUtils.hh"

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
NucleusGenI("genie::NucleusGenINCL")
{

}
//___________________________________________________________________________
NucleusGenINCL::NucleusGenINCL(string config) :
NucleusGenI("genie::NucleusGenINCL", config)
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


  LOG("NucleusGenINCL", pINFO) << "Initialize a nucleus for a new event!";

  // give hit nucleon a vertex
  this->setInitialStateVertex(evrec);
  // give hit nucleon a Fermi momentum
  this->setInitialStateMomentum(evrec);
}

//___________________________________________________________________________

void NucleusGenINCL::GenerateCluster(GHepRecord * evrec) const{

  LOG("NucleusGenINCL", pINFO) << "Initialize a nucleus for a new event!";

  GHepParticle * target_nucleus = evrec->TargetNucleus();
  assert(target_nucleus);
  GHepParticle * nucleon_cluster = evrec->HitNucleon();
  assert(nucleon_cluster);
  Target tgt(target_nucleus->Pdg());
  tgt.SetHitNucPdg(nucleon_cluster->Pdg());

  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  incl_nucleus->reset(&tgt);

  std::shared_ptr<G4INCL::Cluster> incl_cluster =  incl_nucleus->getHitNNCluster();
  G4INCL::Nucleus *nucleus =  incl_nucleus->getNuclues();
  LOG("NucleusGenINCL", pINFO) << incl_cluster->print();
  LOG("NucleusGenINCL", pINFO) << incl_cluster->getMomentum().print();
  LOG("NucleusGenINCL", pINFO) << incl_cluster->getPosition().print();


  G4INCL::ThreeVector cluster_posi = incl_cluster->getPosition();
  TVector3 vtx(cluster_posi.getX(), cluster_posi.getY(), cluster_posi.getZ());
  // Copy the vertex info to the particles already in the event  record
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) )
  {
    if(pdg::IsPseudoParticle(p->Pdg())) continue;
    if(pdg::IsIon           (p->Pdg())) continue;

    LOG("NucleusGenINCL", pDEBUG) << "Setting vertex position for: " << p->Name();
    p->SetPosition(vtx.x(), vtx.y(), vtx.z(), 0.);
  }

  evrec->HitNucleon()->SetPosition(vtx.x(), vtx.y(), vtx.z(), 0.);


  G4INCL::ThreeVector cluster_mom(0, 0, 0);
  double cluster_energy = 0;
  G4INCL::ParticleList particles = incl_cluster->getParticleList();
  for(G4INCL::ParticleIter i=particles.begin(), e=particles.end(); i!=e; ++i) {
    double localEnergy = G4INCL::KinematicsUtils::getLocalEnergy(nucleus, (*i));
    double oldEnergy = (*i)->getEnergy();
    (*i)->setEnergy(oldEnergy - localEnergy);
    (*i)->adjustMomentumFromEnergy();
    cluster_mom+=(*i)->getMomentum();
    cluster_energy += (*i)->getEnergy();
    (*i)->setEnergy(oldEnergy);
    (*i)->adjustMomentumFromEnergy();
  }

  LOG("NucleusGenINCL", pINFO) << cluster_mom.print();
  LOG("NucleusGenINCL", pINFO) << cluster_energy;

  //G4INCL::ThreeVector cluster_mom = incl_cluster->getMomentum();
  TLorentzVector p4nclust   (   cluster_mom.getX() / 1000.,    
      cluster_mom.getY() / 1000.,
      cluster_mom.getZ() / 1000.,
      cluster_energy / 1000.   );

  nucleon_cluster->SetMomentum(p4nclust);
  cluster_bind.SetPxPyPzE(0.,0.,0.,0.);
  LOG("NucleusGenINCL", pINFO) << "success!";
//  abort();

}

//___________________________________________________________________________
//  using INCL model to get the position and momentum of 
//  Hit  nucleon
void NucleusGenINCL::setInitialStateVertex(GHepRecord * evrec) const{

  // get the target and use it to initialize the incl nucleus
  Target* tgt = evrec->Summary()->InitState().TgtPtr();
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  incl_nucleus->reset(tgt);

// generate a vtx and set it to all GHEP physical particles
  Interaction * interaction = evrec->Summary();
  GHepParticle * nucltgt = evrec->TargetNucleus();
  TVector3 vtx(9999999.,999999.,999999.);
  if(!nucltgt){
    vtx.SetXYZ(0.,0.,0.);
  }else{
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

  // Copy the vertex info to the particles already in the event  record
  //
  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) )
  {
    if(pdg::IsPseudoParticle(p->Pdg())) continue;
    if(pdg::IsIon           (p->Pdg())) continue;

    LOG("NucleusGenINCL", pDEBUG) << "Setting vertex position for: " << p->Name();
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

  // get a random nucleon with respect to the isospin of evrec->HitNucleon();
  // the removal energy maybe not necessary
  TVector3 p3 = incl_nucleus->getHitNucleonMomentum();
  double   hit_nucleon_energy = incl_nucleus->getHitNucleonEnergy();
  // double   w  = incl_nucleus->getRemovalEnergy();
  //-- update the struck nucleon 4p at the interaction summary and at
  // the GHEP record
  p4->SetPx(p3.Px()/1000.);
  p4->SetPy(p3.Py()/1000.);
  p4->SetPz(p3.Pz()/1000.);
  p4->SetE ((hit_nucleon_energy)/1000.);

  nucleon->SetMomentum(*p4);  // update GHEP value
  nucleon->SetRemovalEnergy(0);  // FIXME this may be not necessary


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

  this->setTargetNucleusRemnant(evrec);
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

  Px *= -1;
  Py *= -1;
  Pz *= -1;
  E = Mi-E;

  // Add the nucleus to the event record
  LOG("FermiMover", pINFO)
     << "Adding nucleus [A = " << A << ", Z = " << Z
     << ", pdgc = " << ipdgc << "]";

  int imom = evrec->TargetNucleusPosition();
  evrec->AddParticle(
       ipdgc, kIStStableFinalState, imom,-1,-1,-1, Px,Py,Pz,E, 0,0,0,0);
}

//___________________________________________________________________________
void NucleusGenINCL::BindHitNucleon(Interaction& interaction, double& Eb, QELEvGen_BindingMode_t hitNucleonBindingMode) const {
  // BindHitNucleon assign the value of four momentum to class Interaction
  // Eb and hitNucleonBindingMode will not be used in INCL nuclear model
  (void)Eb; (void)hitNucleonBindingMode;
  Target* tgt = interaction.InitState().TgtPtr();
  TLorentzVector* p4Ni = tgt->HitNucP4Ptr();

  if(!flag_isRadius){
    p4Ni->SetVect(TVector3(0.,0.,0.));
    // Look up the (on-shell) mass of the initial nucleon
    TDatabasePDG* tb = TDatabasePDG::Instance();
    double mNi = tb->GetParticle( tgt->HitNucPdg() )->Mass();
    p4Ni->SetE(mNi);
    return;
  }


  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  // get a random nucleon with respect to the isospin of evrec->HitNucleon();
  // the removal energy maybe not necessary
  TVector3 p3 = incl_nucleus->getHitNucleonMomentum();
  double   hit_nucleon_energy = incl_nucleus->getHitNucleonEnergy();
  // Update the initial nucleon lab-frame 4-momentum in the interaction with
  // its current components
  p4Ni->SetVect(TVector3(p3.X()/1000., p3.Y()/1000., p3.Z()/1000.));
  p4Ni->SetE(hit_nucleon_energy/1000.);

}
//___________________________________________________________________________


void NucleusGenINCL::GenerateNucleon(Interaction* interaction, ResamplingHitNucleon_t resampling_mode) const{

  // isRadius:
  // isRadius == 0 (false) put the nucleon in origin
  // isRadius == 1 (true) sampling the radius for nucleon
  //
  // Call the GenerateNucleon will reset the INCLNucleus and generate a new nucleus
  if(! interaction->InitState().Tgt().IsNucleus()) return;
  Target* tgt = interaction->InitState().TgtPtr();
  flag_isRadius = true; // initialize true 
  if(resampling_mode == isOrigin){
    tgt->SetHitNucPosition(0.);
    flag_isRadius = false;
    return;
  }
  else if(resampling_mode == fixRadius){
    INCLNucleus *incl_nucleus = INCLNucleus::Instance();
    incl_nucleus->ResamplingHitNucleon();
    //TVector3 vertex_pos = incl_nucleus->getHitNucleonPosition();
    //double radius = vertex_pos.Mag();
    //tgt->SetHitNucPosition( radius );
  }
  else if(resampling_mode == BothRPResamping){
    INCLNucleus *incl_nucleus = INCLNucleus::Instance();
    incl_nucleus->reset(tgt);
    TVector3 vertex_pos = incl_nucleus->getHitNucleonPosition();
    double radius = vertex_pos.Mag();
    tgt->SetHitNucPosition( radius );
  }
}

//___________________________________________________________________________

// is function is a util function for QEL
bool NucleusGenINCL::isRPValid(double r, double p, const Target & tgt) const {
  (void) tgt;
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  return incl_nucleus->isRPValid(r, p);
}

void NucleusGenINCL::SetHitNucleonOnShellMom(TVector3 p3) const {
  (void) p3;

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

  std::string inclxxpath;
  GetParamDef( "inclxx-data-dir", inclxxpath, std::string(""));
  LOG("NucleusGenINCL", pINFO) << this->expandEnvironmentPath(inclxxpath);

  std::string ablaxxpath;
  GetParamDef( "inclxx-ablaxx-dir", ablaxxpath, std::string(""));
  LOG("NucleusGenINCL", pINFO) << this->expandEnvironmentPath(ablaxxpath);

  std::string abla07path;
  GetParamDef( "inclxx-abla07-dir", abla07path, std::string(""));
  LOG("NucleusGenINCL", pINFO) << this->expandEnvironmentPath(abla07path);

  std::string geminixxpath;
  GetParamDef( "inclxx-geminixx-dir", geminixxpath, std::string(""));
  LOG("NucleusGenINCL", pINFO) << this->expandEnvironmentPath(geminixxpath);

  std::string deExType;
  G4INCL::DeExcitationType deExcitationType = G4INCL::DeExcitationABLA07;
  GetParamDef( "inclxx-de-excitation", deExType, std::string(""));
  LOG("NucleusGenINCL", pINFO) << "inclxx-de-excitation : " << deExType;
  if(!deExType.compare("ABLA07")){
    deExcitationType = G4INCL::DeExcitationABLA07;
  } else if(!deExType.compare("ABLAXX")) {
    deExcitationType = G4INCL::DeExcitationABLAXX;
  } else if(!deExType.compare("GEMINIXX")) {
    deExcitationType = G4INCL::DeExcitationGEMINIXX;
  } else {
    std::stringstream ss;
    ss<< "########################################################\n"
      << "###              WARNING WARNING WARNING             ###\n"
      << "###                                                  ###\n"
      << "### You are running the code without any coupling to ###\n"
      << "###              a de-excitation model!              ###\n"
      << "###    Results will be INCOMPLETE and UNPHYSICAL!    ###\n"
      << "###    Are you sure this is what you want to do?     ###\n"
      << "########################################################\n";
    LOG("NucleusGenINCL", pWARN) << '\n' << ss.str();
  }
  G4INCL::PotentialType potentialType = G4INCL::IsospinEnergyPotential;
  std::string potentialType_str;
  GetParamDef( "inclxx-potential", potentialType_str, std::string(""));
  LOG("NucleusGenINCL", pINFO) << "inclxx-potential : " << potentialType_str;
  if(!potentialType_str.compare("IsospinEnergyPotential")){
    potentialType = G4INCL::IsospinEnergyPotential;
  LOG("NucleusGenINCL", pINFO) << "inclxx-potential : " << potentialType;
  }
  else if(!potentialType_str.compare("IsospinEnergySmoothPotential")){
    potentialType = G4INCL::IsospinEnergySmoothPotential;
  LOG("NucleusGenINCL", pINFO) << "inclxx-potential : " << potentialType;
  }
  else if(!potentialType_str.compare("IsospinPotential")){
    potentialType = G4INCL::IsospinPotential;
  LOG("NucleusGenINCL", pINFO) << "inclxx-potential : " << potentialType;
  }
  else if(!potentialType_str.compare("ConstantPotential")){
    potentialType = G4INCL::ConstantPotential;
  LOG("NucleusGenINCL", pINFO) << "inclxx-potential : " << potentialType;
  }
  LOG("NucleusGenINCL", pINFO) << "inclxx-potential : " << potentialType;

  G4INCL::PauliType pauliType = G4INCL::StrictStatisticalPauli;
  std::string pauliString;
  GetParamDef( "inclxx-pauli-block", pauliString, std::string(""));
  LOG("NucleusGenINCL", pINFO) << "inclxx-pauli-block : " << pauliString;
  if(!pauliString.compare("StrictStatisticalPauli")){
    pauliType = G4INCL::StrictStatisticalPauli;
  }
  else if(!pauliString.compare("StatisticalPauli")){
    pauliType = G4INCL::StatisticalPauli;
  }
  else if(!pauliString.compare("StrictPauli")){
    pauliType = G4INCL::StrictPauli;
  }
  else if(!pauliString.compare("GlobalPauli")){
    pauliType = G4INCL::GlobalPauli;
  }
  else if(!pauliString.compare("NoPauli")){
    pauliType = G4INCL::NoPauli;
  }
  LOG("NucleusGenINCL", pINFO) << "inclxx-pauli-block : " << pauliType;

  // Local Energy type BB
  G4INCL::LocalEnergyType localEnergyTypeBB = G4INCL::FirstCollisionLocalEnergy;
  std::string localEnergyStringBB;
  GetParamDef( "local-energy-BB", localEnergyStringBB, std::string(""));
  LOG("NucleusGenINCL", pINFO) << "local-energy-BB: " << localEnergyStringBB;
  if(!localEnergyStringBB.compare("always")){
    localEnergyTypeBB = G4INCL::AlwaysLocalEnergy;
  }
  else if(!localEnergyStringBB.compare("first-collision")){
    localEnergyTypeBB = G4INCL::FirstCollisionLocalEnergy;
  }
  else if(!localEnergyStringBB.compare("never")){
    localEnergyTypeBB = G4INCL::FirstCollisionLocalEnergy;
  }

  // Local Energy type pi
  G4INCL::LocalEnergyType localEnergyTypepi = G4INCL::FirstCollisionLocalEnergy;
  std::string localEnergyStringpi;
  GetParamDef( "local-energy-pi", localEnergyStringpi, std::string(""));
  LOG("NucleusGenINCL", pINFO) << "local-energy-pi: " << localEnergyStringpi;
  if(!localEnergyStringpi.compare("always")){
    localEnergyTypepi = G4INCL::AlwaysLocalEnergy;
  }
  else if(!localEnergyStringpi.compare("first-collision")){
    localEnergyTypepi = G4INCL::FirstCollisionLocalEnergy;
  }
  else if(!localEnergyStringpi.compare("never")){
    localEnergyTypepi = G4INCL::FirstCollisionLocalEnergy;
  }

  double hadronizationTime = 0.0;
  GetParamDef( "hadronizationTime", hadronizationTime, 0.0);


  // set the Cluster Algorithm
  std::string clusterAlgorithmString = "intercomparison";
  G4INCL::ClusterAlgorithmType clusterAlgorithmType = G4INCL::IntercomparisonClusterAlgorithm;
  GetParamDef( "cluster-algorithm", clusterAlgorithmString, std::string("intercomparison"));
  if(!clusterAlgorithmString.compare("intercomparison")){
    clusterAlgorithmType = G4INCL::IntercomparisonClusterAlgorithm;
  }
  else if(!clusterAlgorithmString.compare("none")){
    clusterAlgorithmType = G4INCL::NoClusterAlgorithm;
  }

  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  incl_nucleus->setINCLXXDataFilePath(this->expandEnvironmentPath(inclxxpath));
  incl_nucleus->setABLAXXDataFilePath(this->expandEnvironmentPath(ablaxxpath));
  incl_nucleus->setABLA07DataFilePath(this->expandEnvironmentPath(abla07path));
  incl_nucleus->setGEMINIXXDataFilePath(this->expandEnvironmentPath(geminixxpath));
  incl_nucleus->setDeExcitationType(deExcitationType);
  incl_nucleus->setPotentialType(potentialType);
  incl_nucleus->setPauliType(pauliType);
  incl_nucleus->setPauliString(pauliString);

  incl_nucleus->setLocalEnergyBBType(localEnergyTypeBB);
  incl_nucleus->setLocalEnergyPiType(localEnergyTypepi);
  incl_nucleus->setHadronizationTime(hadronizationTime);

  incl_nucleus->setClusterAlgorithmType(clusterAlgorithmType);
  incl_nucleus->setClusterAlgorithmString(clusterAlgorithmString);

  incl_nucleus->configure();
}
//____________________________________________________________________________

std::string NucleusGenINCL::expandEnvironmentPath(const std::string& path){
  std::regex env_var_pattern(R"(\$\{([^}]+)\})"); // Regex to match $VAR_NAME
  std::smatch match;

  std::string expanded_path = path;

  while (std::regex_search(expanded_path, match, env_var_pattern)) {
    std::string env_var = match[1].str(); // Extract variable name
    const char* env_value = std::getenv(env_var.c_str());

    if (env_value) {
      expanded_path.replace(match.position(0), match.length(0), env_value);
    } else {
      LOG("NucleusGenINCL", pERROR) << "Warning: Environment variable " << env_var << " is not set!";
      expanded_path.replace(match.position(0), match.length(0), ""); // Remove unresolved variables
    }
  }

  return expanded_path;
}



#endif // end  __GENIE_INCL_ENABLED__
