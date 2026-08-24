#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_INCL_ENABLED__

//---------------------
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

// GENIE
#include "Physics/HadronTransport/INCLCascadeIntranuke.h"
#include "Physics/HadronTransport/G4INCLGENIECascadeAction.h"
#include "Framework/ParticleData/BaryonResUtils.h"
#include "Framework/Algorithm/AlgConfigPool.h"

// INCL++
#include "G4INCLConfig.hh"
#include "G4INCLCascade.hh"
#include "G4INCLConfigEnums.hh"
#include "G4INCLParticle.hh"
// signal handler (for Linux and GCC)
#include "G4INCLSignalHandling.hh"


#include "G4INCLPauliBlocking.hh"
#include "G4INCLCrossSections.hh"
#include "G4INCLDecayAvatar.hh"
#include "G4INCLPhaseSpaceGenerator.hh"
#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"
#include "G4INCLNuclearDensityFactory.hh"
#include "G4INCLINuclearPotential.hh"
#include "G4INCLCoulombDistortion.hh"
#include "G4INCLClustering.hh"
#include "G4INCLIntersection.hh"
#include "G4INCLBinaryCollisionAvatar.hh"
#include "G4INCLCascadeAction.hh"
#include "G4INCLAvatarDumpAction.hh"
#include "G4INCLClusterDecay.hh"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/NuclearState/NucleusGenI.h"

// Generic de-excitation interface
#include "G4INCLIDeExcitation.hh"

// ABLA v3p de-excitation
#ifdef INCL_DEEXCITATION_ABLAXX
#include "G4INCLAblaInterface.hh"
#endif

// ABLA07 de-excitation
#ifdef INCL_DEEXCITATION_ABLA07
#include "G4INCLAbla07Interface.hh"
#endif

// SMM de-excitation
#ifdef INCL_DEEXCITATION_SMM
#include "G4INCLSMMInterface.hh"
#endif

// GEMINIXX de-excitation
#ifdef INCL_DEEXCITATION_GEMINIXX
#include "G4INCLGEMINIXXInterface.hh"
#endif

// --------------------------------------Include for GENIE---------------------
// GENIE

#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Conventions/Units.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Physics/NuclearState/INCLNucleus.h"
#include "Physics/HadronTransport/G4INCLGENIEAvatar.h"
#include "Physics/HadronTransport/G4INCLGENIEParticleRecord.h"

using namespace genie;
using namespace genie::utils;
using namespace G4INCL;
using namespace std;

INCLCascadeIntranuke::INCLCascadeIntranuke() :
  EventRecordVisitorI("genie::INCLCascadeIntranuke"),
  theINCLConfig(0), theINCLModel(0), minRemnantSize(4), fRemnantFullyDecayed(false),
  cascadeAction(std::make_unique<G4INCL::GENIECascadeAction>())
{
  LOG("INCLCascadeIntranuke", pDEBUG)
    << "default ctor";
}

//______________________________________________________________________________
INCLCascadeIntranuke::INCLCascadeIntranuke(string config) :
  EventRecordVisitorI("genie::INCLCascadeIntranuke", config),
  theINCLConfig(0), theINCLModel(0), minRemnantSize(4), fRemnantFullyDecayed(false),
  cascadeAction(std::make_unique<G4INCL::GENIECascadeAction>())
{
  LOG("INCLCascadeIntranuke", pDEBUG)
    << "ctor from config " << config;
}

//______________________________________________________________________________
INCLCascadeIntranuke::~INCLCascadeIntranuke()
{

  // Config is owned by model once handed over
  if ( theINCLConfig   ) { theINCLConfig=0;   }
  if ( theINCLModel    ) { delete theINCLModel;    theINCLModel=0;    }

}

//______________________________________________________________________________
void INCLCascadeIntranuke::LoadConfig(void)
{
  fResonanceDecayer = nullptr;
  fResonanceDecayer = dynamic_cast<const EventRecordVisitorI *> (this->SubAlg("Decayer"));
  assert(fResonanceDecayer);

  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  theConfig = incl_nucleus->getConfig();
  // if the inclxx data file path is not setup, then initialize the genie::NucleusGenINCL
  if(theConfig->getINCLXXDataFilePath().size() == 0){ 
    AlgFactory * algf = AlgFactory::Instance();
    // initialize INCL configuration for hadron model, will never use it.
    auto incl = dynamic_cast<const NucleusGenI*> (algf->GetAlgorithm("genie::NucleusGenINCL", "Default")); 
    assert(incl);
  }
}

//______________________________________________________________________________
// a helper function to add final state from INCL hadron nucleus scattering to GENIE GHEP Event Record.
void INCLCascadeIntranuke::AddINCLParticle(int i, G4INCL::EventInfo &result, GHepRecord * evrec, int first_mother, int second_mother) const {
  double EKin  = result.EKin[i];
  double px  = result.px[i];
  double py  = result.py[i];
  double pz  = result.pz[i];
  int pdg = result.PDGCode[i];
  //short origin = result.origin[i];
  double mass = ((px*px + py*py + pz*pz) - EKin * EKin) / (2.0 * EKin);
  double E = EKin + mass;

  TLorentzVector p4(px* units::MeV / units::GeV,
      py * units::MeV / units::GeV, 
      pz * units::MeV / units::GeV, 
      E * units::MeV / units::GeV);
  TLorentzVector x4(0,0,0,0);
  if(pdg > 1000){
    if(pdg != 2212 && pdg != 2112 && pdg != 3122){
      int S = pdg / 1000000;
      int pdg_no_s = pdg % 1000000;
      int A = pdg_no_s % 1000;
      int Z = pdg_no_s / 1000;
      pdg = genie::pdg::IonPdgCode( A , Z, S, 0);
      pdg = this->INCLPDG_to_GHEPPDG(pdg, A, Z, S);
    }
  }
  evrec->AddParticle(pdg, kIStStableFinalState, first_mother, second_mother, -1, -1, p4, x4);
}
//______________________________________________________________________________
int INCLCascadeIntranuke::doCascade(GHepRecord * evrec) const {

  // do hadron nucleus cascade 
  // re-implement the interface from inclxx/main/src/INCLCascade.cc
  int tpos = evrec->TargetNucleusPosition();
  GHepParticle * target = evrec->Particle(tpos);
  GHepParticle * pprobe = evrec->Probe();
  G4INCL::ParticleType theType = this->PDG_to_INCLType(pprobe->Pdg());
  G4INCL::ParticleSpecies  theSpecies(theType);
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  // setup INCL config 
  theConfig = incl_nucleus->getConfig();
  theConfig->setProjectileSpecies(theSpecies);
  theConfig->setProjectileKineticEnergy((pprobe->E() - pprobe->Mass()) * units::MeV / units::GeV);
  theConfig->setTargetA(target->A());
  theConfig->setTargetZ(target->Z());
  theConfig->setTargetS(0);


  // initialize INCL model
  theINCLModel = new G4INCL::INCL(theConfig);
  G4INCL::EventInfo result;
  // process INCL event
  result = theINCLModel->processEvent();
  if( result.transparent ) {
    evrec->AddParticle(pprobe->Pdg(), kIStStableFinalState, 0, -1, -1, -1, *pprobe->P4(),*pprobe->X4());
    evrec->AddParticle(target->Pdg(), kIStFinalStateNuclearRemnant, 1, -1, -1, -1, *target->P4(),*target->X4());
  }
  else{

    int n_outgoing = 0;
    // add INCL final states from cascade to GHEP Event Record
    for(int i=0; i<result.nParticles; ++i) {
      INCLCascadeIntranuke::AddINCLParticle(i, result, evrec, 0, 1);
      n_outgoing++;
    }

    // add nuclear remnant from cascade to GHEP Event Record
    // do we always have only one remnant?
    double Rem_px = result.pxRem[0]* units::MeV / units::GeV;
    double Rem_py = result.pyRem[0]* units::MeV / units::GeV;
    double Rem_pz = result.pzRem[0]* units::MeV / units::GeV;
    double Rem_Kin = result.EKinRem[0]* units::MeV / units::GeV;
    double Rem_EStar = result.EStarRem[0]* units::MeV / units::GeV;
    int A = result.ARem[0];
    int Z = result.ZRem[0];
    int S = result.SRem[0];
    int pdg = genie::pdg::IonPdgCode( A , Z, S, 0 );
    TParticlePDG * prem = PDGLibrary::Instance()->Find(pdg);
    int PreDeExPDG = pdg;
    double Rem_E = 0;
    if(!prem){ 
      PreDeExPDG = kPdgHadronicBlob; 
      double Rem_p2 = std::sqrt(Rem_px*Rem_px + Rem_py*Rem_py + Rem_pz*Rem_pz);
      double Rem_mass = (Rem_p2*Rem_p2 - Rem_Kin*Rem_Kin) / ( 2.0 * Rem_Kin);
      Rem_E = Rem_Kin + Rem_mass;
    }
    else{
      double Rem_MStar = Rem_EStar + prem->Mass();
      Rem_E  = std::sqrt(Rem_px*Rem_px + Rem_py*Rem_py + Rem_pz*Rem_pz + Rem_MStar*Rem_MStar);
    }
    TLorentzVector p4mom(Rem_px, Rem_py, Rem_pz, Rem_E);
    TLorentzVector p4posi(0,0,0,0);
    evrec->AddParticle(PreDeExPDG, kIStPreDeExNuclearRemnant, 1, -1, -1, -1, p4mom, p4posi);

    // running de-excitation model

    switch(theConfig->getDeExcitationType()){
      case G4INCL::DeExcitationABLAXX:
        {
          std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLAblaInterface>(theConfig);
          theDeExcitation->deExcite(&result);
          break;
        }
      case G4INCL::DeExcitationABLA07:
        {
          std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<ABLA07CXX::Abla07Interface>(theConfig);
          theDeExcitation->deExcite(&result);
          break;
        }
      case G4INCL::DeExcitationGEMINIXX:
        {
          std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLGEMINIXXInterface>(theConfig);
          theDeExcitation->deExcite(&result);
          break;
        }
      default:
        {
          // FIXME: error message
          exit(1);
        }
    }

    // add final states from de-excitation to GHEP Event Record
    for(int i=n_outgoing; i<result.nParticles; ++i) {
      INCLCascadeIntranuke::AddINCLParticle(i, result, evrec, n_outgoing + 2);
      if(i == (result.nParticles - 1) && ( pdg > 1000000000 && pdg%10000/10 > 4)){
        evrec->Particle(i)->SetStatus(kIStFinalStateNuclearRemnant);
      }
    }
  }
  return 0;
}


void INCLCascadeIntranuke::ProcessEventRecord(GHepRecord * evrec)  const {
  LOG("INCLCascadeIntranuke", pINFO) << "Start with this event";

  fGMode = evrec->EventGenerationMode();
  if(fGMode == kGMdHadronNucleus || fGMode == kGMdPhotonNucleus){
    this->doCascade(evrec);
    return;
  }

  this->PreparePrimaryVertex(evrec);

  INCLNucleus *incl_nucleus = INCLNucleus::Instance();
  theConfig = incl_nucleus->getConfig();
  incl_target = incl_nucleus->getNuclues();
  propagationModel =  incl_nucleus->getPropagationModel();
  incl_target->setParticleNucleusCollision();
  fRemnantFullyDecayed = false;
  std::unique_ptr<G4INCL::FinalState> finalState(new FinalState);

  // passing the GHepRecord to cascadeAction
  cascadeAction->setGHepRecord(evrec);
  cascadeAction->beforeRunUserAction(theConfig);
  prob = evrec->Probe();

  // Set the minimum remnant size
  int theA = incl_target->getA();
  minRemnantSize = std::min(theA - 1, 4);

  const bool canRunCascade = this->preCascade();
  if(!canRunCascade){
    LOG("INCLCascadeIntranuke", pNOTICE) << "Never happened!";
    exit(0);
  }


  // TODO: stopping time: 
  // INCL don't have the stopping time for neutrino.
  // We can calculate the longest stopping time for all the daughters from primary interaction.
  // We need to tune this maximumTime.
  double maximumTime = 29.8 * std::pow(incl_target->getA(), 0.16);
  propagationModel->setStoppingTime(maximumTime);
  propagationModel->setNucleus(incl_target);
  propagationModel->generateAllAvatars();

  // before neutrino cascade action
  cascadeAction->beforeCascadeUserAction(propagationModel);
  std::shared_ptr<G4INCL::IAvatar> npv_avatar = this->fillFinalState(evrec, finalState.get(), cascadeAction->getEventRecord());
  cascadeAction->afterNPVAvatarUserAction(npv_avatar.get(), incl_target, finalState.get());

  incl_target->applyFinalState(finalState.get());

  //    incl_target->getStore()->getBook().incrementCascading();   // FIXME
  incl_target->getStore()->getBook().incrementAcceptedCollisions();
  unsigned long loopCounter = 0;
  const unsigned long maxLoopCounter = 10000000;
  while(loopCounter < maxLoopCounter && continueCascade()){
    // Run book keeping actions that should take place before propagation:
    cascadeAction->beforePropagationUserAction(propagationModel);

    // Get the avatar with the smallest time and propagate particles
    // to that point in time.
    IAvatar *avatar = propagationModel->propagate(finalState.get());

    finalState->reset();

    // Run book keeping actions that should take place after propagation:
    cascadeAction->afterPropagationUserAction(propagationModel, avatar);

    if(avatar == 0) break; // No more avatars in the avatar list.

    cascadeAction->beforeAvatarUserAction(avatar, incl_target);

    avatar->fillFinalState(finalState.get());
    // Must fill event record before incl nucleus applyFinalState
    // applyFinalState will delete destroyed particles.

    cascadeAction->afterAvatarUserAction(avatar, incl_target, finalState.get());

    incl_target->applyFinalState(finalState.get());

    LOG("INCLCascadeIntranuke", pNOTICE) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
    delete avatar;
    loopCounter++;
  }
  LOG("INCLCascadeIntranuke", pDEBUG) << "cascade loops: " << loopCounter;
  cascadeAction->afterCascadeUserAction(incl_target);


  primarylepton = evrec->FinalStatePrimaryLepton();

  // put the nuclear remnant in the event record
  LOG("INCLCascadeIntranuke", pDEBUG) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  this->postCascade(finalState.get());

  TObjArrayIter piter(evrec);
  GHepParticle * fsp = nullptr;
  int stable_finalstate = 0;
  while ( (fsp = (GHepParticle *) piter.Next() ) ) {
    if(fsp->Status() == kIStStableFinalState){
      stable_finalstate++;
    }
  }
  const int n_outgoing = theEventInfo.nParticles;
  LOG("INCLCascadeIntranuke", pDEBUG) << "number of stable final state from cascade in INCL and GENIE record: stable_finalstate and n_outgoing: " << stable_finalstate << "  " << n_outgoing;
  

  if(fRemnantFullyDecayed){
    // decayMe() phase-space decayed the target remnant (Z==0 or N==0): the pre-de-excitation
    // entry and every nucleon are already in the record; nothing is left to de-excite.
    LOG("INCLCascadeIntranuke", pINFO) << "Target remnant fully decayed by decayMe(); skipping the remnant/de-excitation block";
  }
  else{
  // Get the 4-momentum of excited remnant
  double Rem_p2 = theEventInfo.pxRem[0]*theEventInfo.pxRem[0]
    + theEventInfo.pyRem[0]*theEventInfo.pyRem[0]
    + theEventInfo.pzRem[0]*theEventInfo.pzRem[0];
  double Rem_Kin = theEventInfo.EKinRem[0];
  double Rem_M = (Rem_p2 - Rem_Kin*Rem_Kin) / ( 2.0 * Rem_Kin);
  double Rem_px = theEventInfo.pxRem[0] * units::MeV / units::GeV;
  double Rem_py = theEventInfo.pyRem[0] * units::MeV / units::GeV;
  double Rem_pz = theEventInfo.pzRem[0] * units::MeV / units::GeV;
  double Rem_E =  std::sqrt(Rem_p2 + Rem_M*Rem_M) * units::MeV / units::GeV;


  // 4 momentum of excited remnant 
  TLorentzVector p4ex(Rem_px, Rem_py, Rem_pz, Rem_E);


  int pdg = genie::pdg::IonPdgCode( incl_target->getA() ,
      incl_target->getZ(),
      incl_target->getS(), 0);
  LOG("INCLCascadeIntranuke", pDEBUG) << "remnant pdg: " <<  pdg;

  int PreDeExPDG = pdg;
  TParticlePDG * prem = PDGLibrary::Instance()->Find(pdg);
  if(!prem){ PreDeExPDG = kPdgHadronicBlob; }
  int rem_mom_id = (evrec->Particle(3)->FirstDaughter() == -1) ? 3 : evrec->Particle(3)->FirstDaughter();

  // position of nucleus at origin
  TLorentzVector p4posi(0,0,0,0);
  evrec->AddParticle(PreDeExPDG, kIStPreDeExNuclearRemnant, rem_mom_id, -1, -1, -1, p4ex, p4posi);


  // processing the deexcitation
  switch(theConfig->getDeExcitationType()){
    case G4INCL::DeExcitationABLAXX:
      {
        std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLAblaInterface>(theConfig);
        theDeExcitation->deExcite(&theEventInfo);
        break;
      }
    case G4INCL::DeExcitationABLA07:
      {
        std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<ABLA07CXX::Abla07Interface>(theConfig);
        theDeExcitation->deExcite(&theEventInfo);
        break;
      }
    case G4INCL::DeExcitationGEMINIXX:
      {
        std::unique_ptr<G4INCL::IDeExcitation> theDeExcitation = std::make_unique<G4INCLGEMINIXXInterface>(theConfig);
        theDeExcitation->deExcite(&theEventInfo);
        break;
      }
    default:
      {
        LOG("INCLCascadeIntranuke", pERROR) << "No set up for de-excitation after INCL FSI!!!";
        exit(1);
      }
  }

  // particle from de-excitation
  const int remnant_id = evrec->GetEntries() - 1;
  if(theEventInfo.nParticles > n_outgoing){
    for(int i = n_outgoing; i < theEventInfo.nParticles; i++){
      LOG("INCLCascadeIntranuke", pDEBUG) << "Final State Particles INCL PDG: " << theEventInfo.PDGCode[i];
      int depdg = this->INCLPDG_to_GHEPPDG(theEventInfo.PDGCode[i], theEventInfo.A[i], theEventInfo.Z[i], theEventInfo.S[i]);
      LOG("INCLCascadeIntranuke", pDEBUG) << "Final State Particles GENIE PDG: " << depdg;
      TParticlePDG * p = PDGLibrary::Instance()->Find(depdg);

      double p2 = theEventInfo.px[i]*theEventInfo.px[i]
        + theEventInfo.py[i]*theEventInfo.py[i]
        + theEventInfo.pz[i]*theEventInfo.pz[i];
      double EKin = theEventInfo.EKin[i];
      double mass = (p2 - EKin*EKin) / (2.0 * EKin);

      double M = p->Mass();
      if(std::fabs(mass * units::MeV / units::GeV - M) > 0.05){
        LOG("INCLCascadeIntranuke", pERROR) << " particle from de-excitation is unphysical for (" << depdg << "), mass = " << mass/1000. << "(" << M <<")"; 
      }
      double E = (EKin + mass) * units::MeV / units::GeV;
      TLorentzVector p4de(theEventInfo.px[i] * units::MeV / units::GeV,
          theEventInfo.py[i] * units::MeV / units::GeV,
          theEventInfo.pz[i] * units::MeV / units::GeV, 
          E);

      EGHepStatus ptype;
      if(theEventInfo.A[i] > minRemnantSize){
        ptype = kIStFinalStateNuclearRemnant;
      } else {
        ptype = kIStStableFinalState;
      }
      evrec->AddParticle(depdg, ptype, remnant_id, -1, -1, -1, p4de, TLorentzVector(0,0,0,0));
    }
  }
  else{
    int rem_pdg = genie::pdg::IonPdgCode( incl_target->getA(),
        incl_target->getZ(),
        std::abs(incl_target->getS()),
        0 );
    TParticlePDG * p = PDGLibrary::Instance()->Find(rem_pdg);
    double M = p->Mass();
    Rem_E = sqrt(Rem_p2/1000000. + M*M);
    TLorentzVector p4rem(Rem_px, Rem_py, Rem_pz, Rem_E);
    evrec->AddParticle(rem_pdg, kIStFinalStateNuclearRemnant, remnant_id, -1, -1, -1, p4rem, TLorentzVector(0,0,0,0));
  }
  } // end if(!fRemnantFullyDecayed)

  // check the baryon number conservation
  // FIXME: we don't consider any beyond SM in INCL FSI
  bool is_conservation = this->BaryonNumberConservation(evrec);
  if(!is_conservation){
    exit(1);
  }
  LOG("INCLCascadeIntranuke", pINFO) << "Done with this event";
}


//______________________________________________________________________________
void INCLCascadeIntranuke::Configure(const Registry & config) {

  LOG("INCLCascadeIntranuke", pDEBUG)
    << "Configure from Registry: '" << config.Name() << "'\n"
    << config;

  Algorithm::Configure(config);
  this->LoadConfig();

}

//___________________________________________________________________________
void INCLCascadeIntranuke::Configure(string param_set) {

  LOG("INCLCascadeIntranuke", pDEBUG)
    << "Configure from param_set name: " << param_set;

  Algorithm::Configure(param_set);
  this->LoadConfig();

}

bool INCLCascadeIntranuke::continueCascade() const{

  bool continueCascade_ = true;

  if(propagationModel->getCurrentTime() > propagationModel->getStoppingTime()){
    continueCascade_ = false;
    LOG("INCLCascadeIntranuke", pWARN) << "stop time : " << propagationModel->getCurrentTime() << " : " << propagationModel->getStoppingTime();

  }
  if(incl_target->getStore()->getBook().getCascading()==0 &&
      incl_target->getStore()->getIncomingParticles().empty()){
    continueCascade_ = false;
    LOG("INCLCascadeIntranuke", pWARN) << "stop cascading ";
  }
  if(incl_target->getA() <= minRemnantSize) {
    continueCascade_ = false;
    LOG("INCLCascadeIntranuke", pWARN) << "stop min size ";
  }

  if(incl_target->getTryCompoundNucleus()) {
    continueCascade_ = false;
    LOG("INCLCascadeIntranuke", pWARN) << "stop Compound ";
  }

  return continueCascade_;

}

void INCLCascadeIntranuke::postCascade(G4INCL::FinalState * finalState) const {
  // Fill in the event information
  theEventInfo.stoppingTime = propagationModel->getCurrentTime();

  // The event bias
  theEventInfo.eventBias = (Double_t) Particle::getTotalBias();

  // Check if the nucleus contains strange particles
  theEventInfo.sigmasInside = incl_target->containsSigma();
  theEventInfo.antikaonsInside = incl_target->containsAntiKaon();
  theEventInfo.lambdasInside = incl_target->containsLambda();
  theEventInfo.kaonsInside = incl_target->containsKaon();

  LOG("INCLCascadeIntranuke", pNOTICE) << "n Sigma and anti-Kaon " << incl_target->containsSigma() << "  " << incl_target->containsAntiKaon() << "  " << incl_target->containsLambda() << "  " << incl_target->containsKaon();

  // FIXME: It seems the INCL will transfer Sigma and anti-Kaon
  // into Lambda via interacting with the first nucleon in Store
  // list
  // - However, de-excitation
  // Capture antiKaons and Sigmas and produce Lambda instead
  theEventInfo.absorbedStrangeParticle = this->decayInsideStrangeParticles(finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();

  // Emit strange particles still inside the nucleus
  this->emitInsideStrangeParticles(finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  theEventInfo.emitKaon = this->emitInsideKaon(finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  theEventInfo.emitLambda = this->emitInsideLambda(finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();

  // Check if the nucleus contains deltas
  theEventInfo.deltasInside = incl_target->containsDeltas();

  // Take care of any remaining deltas

  //if(theEventInfo.deltasInside){
  //  G4INCL::ParticleList const &inside = incl_target->getStore()->getOutgoingParticles();
  //}
  //theEventInfo.forcedDeltasOutside = incl_target->decayOutgoingDeltas();   //FIXME: leave resonances to pythia
  theEventInfo.forcedDeltasInside = this->decayInsideDeltas(finalState);
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();

  // Take care of any remaining etas, omegas, neutral Sigmas and/or neutral kaons
  double timeThreshold=theConfig->getDecayTimeThreshold();
  theEventInfo.forcedPionResonancesOutside = this->decayOutgoingPionResonances(timeThreshold, finalState);  //FIXME: leave resonances to pythia
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  this->decayOutgoingSigmaZero(timeThreshold, finalState); //FIXME: leave resonances to pythia
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  this->decayOutgoingNeutralKaon(finalState); //FIXME: leave resonances to pythia
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();

  //this->emitInsidePions(evrec, finalState); // FIXME: is it valid?

  // Apply Coulomb distortion, if appropriate
  // Note that this will apply Coulomb distortion also on pions emitted by
  // unphysical remnants (see decayInsideDeltas). This is at variance with
  // what INCL4.6 does, but these events are (should be!) so rare that
  // whatever we do doesn't (shouldn't!) make any noticeable difference.
  // CoulombDistortion::distortOut(incl_target->getStore()->getOutgoingParticles(), incl_target);

  // If the normal cascade predicted complete fusion, use the tabulated
  // masses to compute the excitation energy, the recoil, etc.
  //if(incl_target->getStore()->getOutgoingParticles().size()==0
  //    && (!incl_target->getProjectileRemnant()
  //	|| incl_target->getProjectileRemnant()->getParticles().size()==0)) {
  if( false && (!incl_target->getProjectileRemnant()
        || incl_target->getProjectileRemnant()->getParticles().size()==0)) {

    LOG("INCLCascadeIntranuke", pINFO) << "Cascade resulted in complete fusion, using realistic fusion kinematics";

    incl_target->useFusionKinematics();

    if(incl_target->getExcitationEnergy()<0.) {
      // Complete fusion is energetically impossible, return a transparent
      LOG("INCLCascadeIntranuke", pINFO) << "Complete-fusion kinematics yields negative excitation energy, returning a transparent!";
      theEventInfo.transparent = true;
      return;
    }

  }
  else {

    // Set the excitation energy
    incl_target->setExcitationEnergy(incl_target->computeExcitationEnergy());

    // Make a projectile pre-fragment out of the geometrical and dynamical
    // spectators
    //theEventInfo.nUnmergedSpectators = makeProjectileRemnant();

    // Compute recoil momentum, energy and spin of the nucleus
    if(incl_target->getA()==1 && minRemnantSize>1) {
      LOG("INCLCascadeIntranuke", pINFO) << "Computing one-nucleon recoil kinematics. We should never be here nowadays, cascade should stop earlier than this.";
    }
    incl_target->computeRecoilKinematics();

    // Make room for the remnant recoil by rescaling the energies of the
    // outgoing particles. FIXME
    if(incl_target->hasRemnant()) this->rescaleOutgoingForRecoil();
  }

  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  theEventInfo.clusterDecay = this->decayOutgoingClusters(finalState) || this->decayMe(finalState); 
  LOG("INCLCascadeIntranuke", pWARN) << "A and Z: " << incl_target->getA() << "  " << incl_target->getZ();
  incl_target->fillEventInfo(&theEventInfo);

}

bool INCLCascadeIntranuke::preCascade() const {
  // Reset theEventInfo
  theEventInfo.reset();
  EventInfo::eventNumber++;

  // Fill in the event information
  // neutrino projectile
  // theEventInfo.projectileType
  // Projectile/target bookkeeping (no annihilation logic)
  // FIXME: don't have the type of projectile yet
  theEventInfo.Ep = prob->P4()->E() * units::GeV / units::MeV;
  theEventInfo.Ap = 0;
  theEventInfo.Zp = 0;
  theEventInfo.Sp = 0;
  theEventInfo.At = incl_target->getA();
  theEventInfo.Zt = incl_target->getZ();
  theEventInfo.St = incl_target->getS();

  // have no impact parameter for neutrino interaction
  // randomly sampling a hit nucleon
  theEventInfo.impactParameter = 0.;
  theEventInfo.effectiveImpactParameter = 0.;

  // transparent always false since the primary vertex success
  // return true for cascade
  theEventInfo.transparent = false;
  return true;
}

void INCLCascadeIntranuke::rescaleOutgoingForRecoil() const {
  TLorentzVector *p4prob = prob->P4();
  TLorentzVector *p4lep  = primarylepton->P4();
  G4INCL::ThreeVector transferQ((p4prob->Px() - p4lep->Px()) * units::GeV / units::MeV,
      (p4prob->Py() - p4lep->Py()) * units::GeV / units::MeV,
      (p4prob->Pz() - p4lep->Pz()) * units::GeV / units::MeV);

  RecoilCMFunctor theRecoilFunctor(incl_target, theEventInfo, transferQ);

  // Apply the root-finding algorithm
  const G4INCL::RootFinder::Solution theSolution = G4INCL::RootFinder::solve(&theRecoilFunctor, 1.0);
  if(theSolution.success) {
    theRecoilFunctor(theSolution.x); // Apply the solution
  } else {
    LOG("INCLCascadeIntranuke", pINFO) << "Couldn't accommodate remnant recoil while satisfying energy conservation, root-finding algorithm failed.";
  }
}

int INCLCascadeIntranuke::INCLPDG_to_GHEPPDG(int pdg, int A, int Z, int S) const{
  // Convert the INCL pdg id to GHep pdg id
  TDatabasePDG * fDatabasePDG = TDatabasePDG::Instance();
  TParticlePDG * p = fDatabasePDG->GetParticle(pdg);
  // pdg code is not in the data base
  if(!p){     
    // It is a nucleus, using A, Z, S to construct the nucleus ID.
    if(A != 0){  
      int ion_pdg =  genie::pdg::IonPdgCode( A , Z, std::abs(S), 0 );
      TParticlePDG *ion = fDatabasePDG->GetParticle(ion_pdg);
      // the hyper-nuclei is not in the database, add a temp one 
      if(S != 0 && !ion){ 
        PDGLibrary *pdg_library = PDGLibrary::Instance();
        pdg_library->AddHypernucleus(ion_pdg);
      }
      return ion_pdg;
    }
    else{
      LOG("INCLCascadeIntranuke", pERROR) << "Particle is not identified: pdg = " << pdg << "; A, Z, S => " 
        << A << " " << Z << " " << S;
      exit(1);
    }
  }
  else{ // find the pdg code in database, return it.
    return pdg;
  }
}



std::shared_ptr<G4INCL::IAvatar> INCLCascadeIntranuke::fillFinalState(GHepRecord * evrec, G4INCL::FinalState * finalState, std::vector<G4INCL::GENIEParticleRecord> *eventRecord) const{

  // neutrino interaction info
  const ProcessInfo & proc_info = evrec->Summary()->ProcInfo();
  INCLNucleus *incl_nucleus = INCLNucleus::Instance();

  // get the nuclear model for neutrino primary vertex
  // GENIE model (if we use GENIE nuclear for the n.p.v., it is a hybrid model)
  //   - kNucmFermiGas
  //   - kNucmLocalFermiGas
  //   - kNucmSpectralFunc
  //   - kNucmEffSpectralFunc
  // INCL model
  //   - kNucmINCL
  NuclearModel_t nucl_model = incl_nucleus->getHybridModel();
  LOG("INCLCascadeIntranuke", pWARN) << "Nuclear model " << NuclearModel::AsString(nucl_model);
  bool isHybridModel = true;
  switch(nucl_model){
    case kNucmINCL: isHybridModel = false; break;
    case kNucmUndefined: 
        LOG("INCLCascadeIntranuke", pERROR) << "Nuclear model is not setup correctly!";
        exit(1);
    default: break;
  }


  std::shared_ptr<G4INCL::IAvatar> avatar;
  if(proc_info.IsMEC()){ // MEC 2p2h channel
    avatar = std::make_shared<G4INCL::GENIEAvatar>(0, (incl_nucleus->getHitNNCluster()).get(), incl_nucleus->getNuclues(), eventRecord, isHybridModel);
    avatar->fillFinalState(finalState);
  }
  else{ // QE, RES, DIS channel, one hit nucleon
    avatar = std::make_shared<G4INCL::GENIEAvatar>(0, incl_nucleus->getHitParticle(), incl_nucleus->getNuclues(), eventRecord, isHybridModel);
    avatar->fillFinalState(finalState);
  }

  if(finalState->getValidity() != G4INCL::ValidFS){
    LOG("INCLCascadeIntranuke", pWARN)
      << "Enforcing energy conservation: failed! ";
    evrec->EventFlags()->SetBitNumber(kKineGenErr, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Couldn't select kinematics");
    exception.SwitchOnFastForward();
    throw exception;
  }
  return avatar;
}

G4INCL::ParticleType INCLCascadeIntranuke::PDG_to_INCLType(int pdg) const {
  switch(pdg){
    case 2212: return G4INCL::Proton;
    case 2112: return G4INCL::Neutron;
    case 211: return G4INCL::PiPlus;
    case -211: return G4INCL::PiMinus;
    case 111: return G4INCL::PiZero;
    case 2224: return G4INCL::DeltaPlusPlus;
    case 2214: return G4INCL::DeltaPlus;
    case 2114: return G4INCL::DeltaZero;
    case 1114: return G4INCL::DeltaMinus;
    case 221: return G4INCL::Eta;
    case 223: return G4INCL::Omega;
    case 331: return G4INCL::EtaPrime;
    case 22: return G4INCL::Photon;
    case 3122: return G4INCL::Lambda;
    case 3222: return G4INCL::SigmaPlus;
    case 3212: return G4INCL::SigmaZero;
    case 3112: return G4INCL::SigmaMinus;
    case -2212: return G4INCL::antiProton;
    case 3312: return G4INCL::XiMinus;
    case 3322: return G4INCL::XiZero;
    case -2112: return G4INCL::antiNeutron;
    case -3122: return G4INCL::antiLambda;
    case -3222: return G4INCL::antiSigmaPlus;
    case -3212: return G4INCL::antiSigmaZero;
    case -3112: return G4INCL::antiSigmaMinus;
    case -3312: return G4INCL::antiXiMinus;
    case -3322: return G4INCL::antiXiZero;
    case 321: return G4INCL::KPlus;
    case 311: return G4INCL::KZero;
    case -311: return G4INCL::KZeroBar;
    case -321: return G4INCL::KMinus;
    case 310: return G4INCL::KShort;
    case 130: return G4INCL::KLong;
    default:
              return G4INCL::UnknownParticle;
  }

}

void INCLCascadeIntranuke::PreparePrimaryVertex(GHepRecord * event_rec) const{
  // 1. for CCQE
  int inucl = -1;
  inucl = event_rec->RemnantNucleusPosition();
  GHepParticle * nucl = event_rec->Particle(inucl);
  nucl->SetStatus(kIStIntermediateState);

  // INCL has INCL::prepareReaction
  // Prepare the target nucleus and reaction parameters 
  // before simulating a projectile (proton, neutron, 
  // antiproton, antineutron, antideuteron, etc.) interacting 
  // with a nucleus of mass number A, atomic number Z, and 
  // strangeness S. 
  //
  // For genie simulation, we only simulate neutrino with transparent falese
  //
  // don't have impact parameter
  //
  // minRemnantSize = std::min(theA-1, 4);
}

void INCLCascadeIntranuke::DecayResonance(GHepRecord *evrec) const{
  if(evrec->Summary()->ProcInfo().IsResonant()){
    bool decay_flag = true;
    while(decay_flag){
      TObjArrayIter piter(evrec);
      GHepParticle * p = nullptr;
      decay_flag = false;
      while ( (p = (GHepParticle *) piter.Next() ) ) {
        if(p->Status() == kIStPreDecayResonantState && p->FirstDaughter() == -1 && PDG_to_INCLType(p->Pdg()) == G4INCL::UnknownParticle){
          decay_flag = true;
        }
      }
      if(decay_flag){
        fResonanceDecayer->ProcessEventRecord(evrec);
      }
    }
  }
}

bool INCLCascadeIntranuke::BaryonNumberConservation(GHepRecord *evrec) const {
  int final_C = 0, final_A = 0;
  Target *tgt = evrec->Summary()->InitStatePtr()->TgtPtr();
  // GHepParticle::Charge() is in |e|/3 (TParticlePDG convention), Target::Charge() in +e
  int inital_C = evrec->Probe()->Charge()/3. + tgt->Charge();
  int target_A = tgt->A();
  LOG("INCLCascadeIntranuke", pINFO) << "Inital states charge and A: " << inital_C << " " << target_A;
  TObjArrayIter piter(evrec);
  GHepParticle * p = nullptr;
  while ( (p = (GHepParticle *) piter.Next() ) ) {
    // the code of the particles in primary neutrino interaction
    std::string particle_class = std::string(PDGLibrary::Instance()->Find(p->Pdg())->ParticleClass());
    if(p->Status() == kIStStableFinalState || p->Status() == kIStFinalStateNuclearRemnant){
      if(genie::pdg::IsNeutronOrProton(p->Pdg())){
        final_A++;
      }
      else if(genie::pdg::IsBaryonResonance(p->Pdg())){
        final_A++;
      }
      else if(particle_class.find("Baryon") != std::string::npos){
        if(p->Pdg() < 0){
          final_A--;
        }
        else{
          final_A++;
        }
      }
      else if(genie::pdg::IsIon(p->Pdg())){
        final_A+=p->A();
      }
      final_C += p->Charge()/3;
      LOG("INCLCascadeIntranuke", pINFO) << "PDG and Charge: " << p->Pdg() << "  " << p->Charge() << "  " << p->A() << "  " << final_A << "  " << final_C;
    }
  }
  LOG("INCLCascadeIntranuke", pINFO) << "Final states charge and A: " << final_C << " " << final_A;
  if(inital_C != final_C || target_A != final_A){
    evrec->Print(std::cout);
    LOG("INCLCascadeIntranuke", pFATAL) << "Violating charge number and baryon number!";
    LOG("INCLCascadeIntranuke", pFATAL) << "Inital states charge and A: " << inital_C << " " << target_A;
    LOG("INCLCascadeIntranuke", pFATAL) << "Final states charge and A: " << final_C << " " << final_A;
    exit(1);
    return false;
  }
  return true;
}

#include "INCLPostCascade.icc"

#endif // __GENIE_INCL_ENABLED__
