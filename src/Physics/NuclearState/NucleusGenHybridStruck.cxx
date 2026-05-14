//____________________________________________________________________________
/*!

\class    genie::NucleusGenHybridStruck

\brief    It visits the event record & combines the FermoMover and VertexGenerator
          computes a Fermi motion momentum and position for initial state nucleons 
	  bound in nuclei. It is configured in NucleusGenerator. This new interface
	  is compatible with previous implments.
          Is a concrete implementation of the EventRecordVisitorI interface.

\author   Liang Liu <liangliu \at fnal.gov>
          Fermi National Accelerator Laboratory

\created  Jan 17, 2025

\cpright  Copyright (c) 2003-2024, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          
*/
//____________________________________________________________________________

#include "Framework/Conventions/GBuild.h"
#include <cstdlib>

#include <TLorentzVector.h>
#include <TVector3.h>
#include <TParticlePDG.h>
#include <TMath.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Physics/NuclearState/NucleusGenHybridStruck.h"

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
#include "Physics/Common/VertexGenerator.h"
#ifdef __GENIE_INCL_ENABLED__
#include "Physics/NuclearState/INCLNucleus.h"
#endif

using namespace genie;
using namespace genie::constants;

//___________________________________________________________________________
NucleusGenHybridStruck::NucleusGenHybridStruck() :
NucleusGenI("genie::NucleusGenHybridStruck")
{

}
//___________________________________________________________________________
NucleusGenHybridStruck::NucleusGenHybridStruck(string config) :
NucleusGenI("genie::NucleusGenHybridStruck", config)
{

}
//___________________________________________________________________________
NucleusGenHybridStruck::~NucleusGenHybridStruck()
{

}

//___________________________________________________________________________
void NucleusGenHybridStruck::ProcessEventRecord(GHepRecord * evrec) const
{
  // skip if not a nuclear target
  if(! evrec->Summary()->InitState().Tgt().IsNucleus()) return;
  this->setInitialStateVertex(evrec);
  this->setInitialStateMomentum(evrec);
}

void NucleusGenHybridStruck::setInitialStateVertex(GHepRecord * evrec) const{
#ifdef __GENIE_INCL_ENABLED__
  if(fINCLVertex){
    // randomly pick up a nucleon from INCL nucleus as the struck nucleon
    fNucleusGen->setInitialStateVertex(evrec);
  }
  else{
#endif
    // using the GENIE vertex model
    fVertexGenerator->ProcessEventRecord(evrec);
#ifdef __GENIE_INCL_ENABLED__
    if(fINCLFSI){
      Target* tgt = evrec->Summary()->InitState().TgtPtr();
      INCLNucleus *incl_nucleus = INCLNucleus::Instance();
      // TODO: simplify this initialization
      // initialization will randomly pick a nucleon as a hit nucleon
      // it will be replaced by GENIE simulation later
      // in the hitParticle->setPosition(ip->X3()); in G4INCLGENIEQELChannel.cxx
      incl_nucleus->reset(tgt);
    }
  }
#endif
}


void NucleusGenHybridStruck::setINCLVertex(GHepRecord * evrec) const {
#ifdef __GENIE_INCL_ENABLED__
  GHepParticle * nucleon = evrec->HitNucleon();
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  TVector3 vtx = incl_nucleus->ResamplingVertex(nucleon->Pdg());

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
#endif
}

void NucleusGenHybridStruck::setInitialStateMomentum(GHepRecord * evrec) const{
  fFermiMover->ProcessEventRecord(evrec);
}

void NucleusGenHybridStruck::setClusterVertex(GHepRecord * evrec) const{
  #ifdef __GENIE_INCL_ENABLED__
  if(fINCLVertex){
    // randomly pick up a nucleon from INCL nucleus as the struck nucleon
    GHepParticle * target_nucleus = evrec->TargetNucleus();
    assert(target_nucleus);
    GHepParticle * nucleon_cluster = evrec->HitNucleon();
    assert(nucleon_cluster);
    Target tgt(target_nucleus->Pdg());
    tgt.SetHitNucPdg(nucleon_cluster->Pdg());

    INCLNucleus *incl_nucleus = INCLNucleus::Instance();
    incl_nucleus->reset(&tgt);

    std::shared_ptr<G4INCL::Cluster> incl_cluster =  incl_nucleus->getHitNNCluster();
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
    
  }
  else{
#endif
    // using the GENIE vertex model
    fVertexGenerator->ProcessEventRecord(evrec);
    evrec->HitNucleon()->SetPosition(*evrec->Probe()->X4());
#ifdef __GENIE_INCL_ENABLED__
    if(fINCLFSI){
      // find the closet nucleon in INCL as the struck nucleon
      // get the target and use it to initialize the incl nucleus
      Target* tgt = evrec->Summary()->InitState().TgtPtr();
      INCLNucleus *incl_nucleus = INCLNucleus::Instance();
      incl_nucleus->reset(tgt);
      // Get the position of hit nucelon
      GHepParticle * nucleon = evrec->HitNucleon();
      TVector3 posi = nucleon->X4()->Vect();
      int nucleon_pdg = tgt->HitNucPdg();
      if(nucleon_pdg == kPdgClusterNN){
        incl_nucleus->setHitNNCluster(kPdgNeutron, kPdgNeutron, posi);
      }
      else if(nucleon_pdg == kPdgClusterNP){
        incl_nucleus->setHitNNCluster(kPdgNeutron, kPdgProton, posi);
      }
      else if(nucleon_pdg == kPdgClusterPP){
        incl_nucleus->setHitNNCluster(kPdgProton, kPdgProton, posi);
      }
      else{
        LOG("INCLNucleus", pFATAL) << "Can't get";
        exit(1);
      }
    }
  }
#endif
}

void NucleusGenHybridStruck::GenerateCluster(GHepRecord * evrec) const{
  // This is function will be used in MEC channela
  // first, generate vertex for hit cluster
  if(! evrec->Summary()->InitState().Tgt().IsNucleus()) return;
  this->setClusterVertex(evrec);

  GHepParticle * target_nucleus = evrec->TargetNucleus();
  assert(target_nucleus);
  GHepParticle * nucleon_cluster = evrec->HitNucleon();
  assert(nucleon_cluster);
  GHepParticle * remnant_nucleus = evrec->RemnantNucleus();
  assert(remnant_nucleus);

  // generate a Fermi momentum for each nucleon

  Target tgt(target_nucleus->Pdg());

  bool allowdup = true;
  PDGCodeList pdgv(allowdup);

  if(nucleon_cluster->Pdg() == kPdgClusterNN) {
    pdgv.push_back(kPdgNeutron);
    pdgv.push_back(kPdgNeutron);
  }
  else
    if(nucleon_cluster->Pdg() == kPdgClusterNP) {
      pdgv.push_back(kPdgNeutron);
      pdgv.push_back(kPdgProton);
    }
    else
      if(nucleon_cluster->Pdg() == kPdgClusterPP) {
        pdgv.push_back(kPdgProton);
        pdgv.push_back(kPdgProton);
      }
      else
      {
        LOG("NucleusGenHybridStruck", pERROR)
          << "Unknown di-nucleon cluster PDG code (" << nucleon_cluster->Pdg() << ")";
        exit(1);
      }

  assert(pdgv.size()==2);

  TVector3 p3a, p3b;
  double removalenergy1, removalenergy2;
  fNuclModel->GenerateCluster(tgt, pdgv, &p3a, &p3b, &removalenergy1, &removalenergy2);

  LOG("NucleusGenHybridStruck", pINFO)
    << "1st nucleon (code = " << pdgv[0] << ") generated momentum: ("
    << p3a.Px() << ", " << p3a.Py() << ", " << p3a.Pz() << "), "
    << "|p| = " << p3a.Mag();
  LOG("NucleusGenHybridStruck", pINFO)
    << "2nd nucleon (code = " << pdgv[1] << ") generated momentum: ("
    << p3b.Px() << ", " << p3b.Py() << ", " << p3b.Pz() << "), "
    << "|p| = " << p3b.Mag();

  TVector3 p3 = p3a + p3b;
  double M2n = PDGLibrary::Instance()->Find(nucleon_cluster->Pdg())-> Mass(); // nucleon cluster mass
  // nucleon cluster energy
  double EN = TMath::Sqrt(p3.Mag2() + M2n*M2n);

  TLorentzVector p4nclust   (   p3.Px(),    p3.Py(),    p3.Pz(),  EN   );
  nucleon_cluster->SetMomentum(p4nclust);
  cluster_bind.SetPxPyPzE(0., 0., 0., -1.0 * (removalenergy1 + removalenergy2));

}
//___________________________________________________________________________
void NucleusGenHybridStruck::BindHitNucleon(Interaction& interaction, double& Eb, QELEvGen_BindingMode_t hitNucleonBindingMode) const {

  Target* tgt = interaction.InitState().TgtPtr();
  TLorentzVector* p4Ni = tgt->HitNucP4Ptr();

  // Initial nucleon 3-momentum (lab frame)
  TVector3 p3Ni = fNuclModel->Momentum3();

  // Look up the (on-shell) mass of the initial nucleon
  TDatabasePDG* tb = TDatabasePDG::Instance();
  double mNi = tb->GetParticle( tgt->HitNucPdg() )->Mass();

  // Set the (possibly off-shell) initial nucleon energy based on
  // the selected binding energy mode. Always put the initial nucleon
  // on shell if it is not part of a composite nucleus
  double ENi = 0.;
  if ( tgt->IsNucleus() && hitNucleonBindingMode != genie::kOnShell ) {

    // For a nuclear target with a bound initial struck nucleon, take binding
    // energy effects and Pauli blocking into account when computing QE
    // differential cross sections
    interaction.ResetBit( kIAssumeFreeNucleon );

    // Initial nucleus mass
    double Mi = tgt->Mass();

    // Final nucleus mass
    double Mf = 0.;

    // If we use the removal energy reported by the nuclear
    // model, then it implies a certain value for the final
    // nucleus mass
    if ( hitNucleonBindingMode == genie::kUseNuclearModel ) {
      Eb = fNuclModel->RemovalEnergy();
      // This equation is the definition that we assume
      // here for the "removal energy" (Eb) returned by the
      // nuclear model. It matches GENIE's convention for
      // the Bodek/Ritchie Fermi gas model.
      Mf = Mi + Eb - mNi;
    }
    // We can also assume that the final nucleus is in its
    // ground state. In this case, we can just look up its
    // mass from the standard table. This implies a particular
    // binding energy for the hit nucleon.
    else if ( hitNucleonBindingMode == genie::kUseGroundStateRemnant ) {
      // Determine the mass and proton numbers for the remnant nucleus
      int Af = tgt->A() - 1;
      int Zf = tgt->Z();
      if ( genie::pdg::IsProton( tgt->HitNucPdg()) ) --Zf;
      Mf = genie::PDGLibrary::Instance()->Find( genie::pdg::IonPdgCode(Af, Zf) )->Mass();

      // Deduce the binding energy from the final nucleus mass
      Eb = Mf - Mi + mNi;
    }
    // A third option is to assign the binding energy and final nuclear
    // mass in a way that is equivalent to the original Valencia model
    // treatment
    else if ( hitNucleonBindingMode == genie::kValenciaStyleQValue ) {
      // Compute the Q-value needed for the ground-state-to-ground-state
      // transition without nucleon removal (see Sec. III B of
      // https://arxiv.org/abs/nucl-th/0408005). Note that this will be
      // zero for NC and EM scattering since the proton and nucleon numbers
      // are unchanged. Also get Q_LFG, the difference in Fermi energies
      // between the final and initial struck nucleon species. This will
      // likewise be zero for NC and EM interactions.
      const genie::ProcessInfo& info = interaction.ProcInfo();
      double Qvalue = 0.;
      double Q_LFG = 0.;
      if ( info.IsWeakCC() ) {
        // Without nucleon removal, the final nucleon number remains the same
        int Af = tgt->A();
        // For CC interactions, the proton number will change by one. Choose
        // the right change by checking whether we're working with a neutrino
        // or an antineutrino.
        int Zf = tgt->Z();
        int probe_pdg = interaction.InitState().ProbePdg();
        if ( genie::pdg::IsAntiNeutrino(probe_pdg) ) {
          --Zf;
        }
        else if ( genie::pdg::IsNeutrino(probe_pdg) ) {
          ++Zf;
        }
        else {
          LOG( "QELEvent", pFATAL ) << "Unhandled probe PDG code " << probe_pdg
            << " encountered for a CC interaction in"
            << " genie::utils::BindHitNucleon()";
          gAbortingInErr = true;
          std::exit( 1 );
        }

        // We have what we need to get the Q-value. Get the final nuclear
        // mass (without nucleon removal)
        double mf_keep_nucleon = genie::PDGLibrary::Instance()
          ->Find( genie::pdg::IonPdgCode(Af, Zf) )->Mass();

        Qvalue = mf_keep_nucleon - Mi;

        // Get the Fermi energies for the initial and final nucleons. Include
        // the radial dependence if using the LFG.
        double hit_nucleon_radius = tgt->HitNucPosition();

        // Average of the proton and neutron masses. It may actually be better
        // to use the exact on-shell masses here. However, the original paper
        // uses the same value for protons and neutrons when computing the
        // difference in Fermi energies. For consistency with the Valencia
        // model publication, I'll do the same.
        const double mN = genie::constants::kNucleonMass;

        double kF_Ni = fNuclModel->LocalFermiMomentum( *tgt,
            tgt->HitNucPdg(), hit_nucleon_radius );
        double EFermi_Ni = std::sqrt( std::max(0., mN*mN + kF_Ni*kF_Ni) );

        double kF_Nf = fNuclModel->LocalFermiMomentum( *tgt,
            interaction.RecoilNucleonPdg(), hit_nucleon_radius );
        // (On-shell) final nucleon mass
        //double mNf = interaction.RecoilNucleon()->Mass();
        double EFermi_Nf = std::sqrt( std::max(0., mN*mN + kF_Nf*kF_Nf) );

        // The difference in Fermi energies is Q_LFG
        Q_LFG = EFermi_Nf - EFermi_Ni;
      }
      // Fail here for interactions that aren't one of EM/NC/CC. This is
      // intended to help avoid silently doing the wrong thing in the future.
      else if ( !info.IsEM() && !info.IsWeakNC() ) {
        LOG( "QELEvent", pFATAL ) << "Unhandled process type \"" << info
          << "\" encountered in genie::utils::BindHitNucleon()";
        gAbortingInErr = true;
        std::exit( 1 );
      }

      // On-shell total energy of the initial struck nucleon
      double ENi_OnShell = std::sqrt( std::max(0., mNi*mNi + p3Ni.Mag2()) );

      // Total energy of the remnant nucleus (with the hit nucleon removed)
      double Ef = Mi - ENi_OnShell + Qvalue - Q_LFG;

      // Mass and kinetic energy of the remnant nucleus (with the hit nucleon
      // removed)
      Mf = std::sqrt( std::max(0., Ef*Ef - p3Ni.Mag2()) );

      // Deduce the binding energy from the final nucleus mass
      Eb = Mf - Mi + mNi;

      LOG( "QELEvent", pDEBUG ) << "Qvalue = " << Qvalue
        << ", Q_LFG = " << Q_LFG << " at radius = " << tgt->HitNucPosition();
    }

    // The (lab-frame) off-shell initial nucleon energy is the difference
    // between the lab frame total energies of the initial and remnant nuclei
    ENi = Mi - std::sqrt( Mf*Mf + p3Ni.Mag2() );
  }
  else {
    // Keep the struck nucleon on shell either because
    // hitNucleonBindingMode == kOnShell or because
    // the target is a single nucleon
    ENi = std::sqrt( p3Ni.Mag2() + std::pow(mNi, 2) );
    Eb = 0.;

    // If we're dealing with a nuclear target but using the on-shell
    // binding mode, set the "assume free nucleon" flag. This turns off
    // Pauli blocking and the requirement that q0 > 0 in the QE cross section
    // models (an on-shell nucleon *is* a free nucleon)
    if ( tgt->IsNucleus() ) interaction.SetBit( kIAssumeFreeNucleon );
  }

  LOG( "QELEvent", pDEBUG ) << "Eb = " << Eb << ", pNi = " << p3Ni.Mag()
    << ", ENi = " << ENi;

  // Update the initial nucleon lab-frame 4-momentum in the interaction with
  // its current components
  p4Ni->SetVect( p3Ni );
  p4Ni->SetE( ENi );



} 

//___________________________________________________________________________

void NucleusGenHybridStruck::GenerateNucleon(Interaction* interaction, ResamplingHitNucleon_t resampling_mode) const {
  // isRadius:
  // isRadius == 0 (false) put the nucleon in origin
  // isRadius == 1 (true) sampling the radius for nucleon
  Target* tgt = interaction->InitState().TgtPtr();
  if(resampling_mode == BothRPResamping){
    //const VertexGenerator* vtx_gen = dynamic_cast<const VertexGenerator*>(fVertexGenerator);
    //TVector3 vertex_pos = vtx_gen->GenerateVertex( interaction, tgt->A() );
    TVector3 vertex_pos = this->GetVertex(interaction);
    double radius = vertex_pos.Mag();

    LOG("NucleusGenHybridStruck", pDEBUG) << "radius : " << radius;

    tgt->SetHitNucPosition( radius );
    fNuclModel->GenerateNucleon(*tgt, radius);
  }
  else if(resampling_mode == isOrigin){
    tgt->SetHitNucPosition(0.);
    if ( tgt->IsNucleus() ){
      fNuclModel->GenerateNucleon(*tgt, 0.);
    }
    else{
      fNuclModel->SetRemovalEnergy(0.);
      interaction->SetBit( kIAssumeFreeNucleon );
    }
    fNuclModel->SetMomentum3( TVector3(0., 0., 0.) );
  }
  else if(resampling_mode == fixRadius){
    LOG("NucleusGenHybridStruck", pNOTICE) << "radius : " << tgt->HitNucPosition();
    fNuclModel->GenerateNucleon(*tgt, tgt->HitNucPosition());
  }
}

TVector3 NucleusGenHybridStruck::GetVertex(Interaction* interaction) const{
  Target* tgt = interaction->InitState().TgtPtr();
#ifdef __GENIE_INCL_ENABLED__
  if(fINCLVertex){
    // randomly pick up a nucleon from INCL nucleus as the struck nucleon
    INCLNucleus *incl_nucleus = INCLNucleus::Instance();
    incl_nucleus->reset(tgt);
    TVector3 vertex_pos = incl_nucleus->getHitNucleonPosition();
    return vertex_pos;

  }
  else{
#endif
    // using the GENIE vertex model
    LOG("NucleusGenHybridStruck", pDEBUG) << "Get new vertex";
    const VertexGenerator* vtx_gen = dynamic_cast<const VertexGenerator*>(fVertexGenerator);
    TVector3 vertex_pos = vtx_gen->GenerateVertex( interaction, tgt->A() );
    return vertex_pos;
#ifdef __GENIE_INCL_ENABLED__
  }
#endif
}


//___________________________________________________________________________

// is function is a util function for QEL
bool NucleusGenHybridStruck::isRPValid(double r, double p, const Target & tgt) const {
  // TODO: document this, won't work for spectral functions
  double dummy_w = -1.;
  double prob = fNuclModel->Prob(p, dummy_w, tgt, r);
  return (prob > 0.);
}

//___________________________________________________________________________

void NucleusGenHybridStruck::SetHitNucleonOnShellMom(TVector3 p3) const {
  // Set the nucleon we're using to be upstream at max energy and unbound
  fNuclModel->SetMomentum3( p3 );
  fNuclModel->SetRemovalEnergy( 0. );
}

//___________________________________________________________________________
void NucleusGenHybridStruck::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGenHybridStruck::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void NucleusGenHybridStruck::LoadConfig(void)
{
  bool is_configed = true;
  fFermiMover = nullptr;
  fVertexGenerator = nullptr;
  fINCLVertex = true;

  this->GetParam("INCL-vertex", fINCLVertex);
  std::cout << "DEBUG : " << fINCLVertex << std::endl;

  fFermiMover = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("FermiMover"));
  if(!fFermiMover){
    is_configed = false;
    LOG("NucleusGenHybridStruck", pERROR) << "The SubAlg FermiMover is not configed";
  }


  fVertexGenerator = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("VertexGenerator"));
  if(!fVertexGenerator){
    is_configed = false;
    LOG("NucleusGenHybridStruck", pERROR) << "The SubAlg VertexGenerator is not configed";
  }

  if(!is_configed){
    LOG("NucleusGenHybridStruck", pERROR) << "Configuration has failed.";
    exit(78);
  }

  RgKey nuclkey = "NuclearModel";
  fNuclModel = nullptr;
  fNuclModel = dynamic_cast<const NuclearModelI *>(fFermiMover->SubAlg(nuclkey));
  assert(fNuclModel);

  fINCLFSI = false;
  Registry * algos = AlgConfigPool::Instance() -> GlobalParameterList() ;
  algos->Print(std::cout);
  bool fHadronTranspEnable = algos -> GetBool("HadronTransp-Enable");
  string fsi_model_name = algos -> GetAlg("HadronTransp-Model").name;
  std::cout << "DEBUG : " << fsi_model_name << std::endl;
  if(fHadronTranspEnable && fsi_model_name == string("genie::INCLCascadeIntranuke")){
    std::cout << "DEBUG : " << fsi_model_name << "  OK!!" << std::endl;
    fINCLFSI = true;
  }
  std::cout << "DEBUG : " << fINCLFSI << std::endl;
  std::cout << "DEBUG : " << fINCLVertex << std::endl;

  // Using NucleusGenINCL to setup the configuration of INCL nuclear model.
  // Using INCLNucleus directly instead of NucleusGenINCL 
#ifdef __GENIE_INCL_ENABLED__
  if(fINCLVertex || fINCLFSI){
    std::cout << "DEBUG : " << fINCLFSI << std::endl;
    std::cout << "DEBUG : " << fINCLVertex << std::endl;
    RgKey nuclgenkey = "NucleusGenerator";
    fNucleusGen = nullptr;
    fNucleusGen = dynamic_cast<const NucleusGenI *> (this->SubAlg(nuclgenkey));
    assert(fNucleusGen);
    INCLNucleus *incl_nucleus = INCLNucleus::Instance();
    incl_nucleus->setHybridModel(fNuclModel->ModelType(Target()));
  }
#endif
}
//____________________________________________________________________________

