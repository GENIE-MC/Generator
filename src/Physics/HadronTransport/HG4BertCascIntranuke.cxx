#include "Framework/Conventions/GBuild.h"
#ifdef __GENIE_GEANT4_INTERFACE_ENABLED__
//____________________________________________________________________________
/*

 Author: Dennis Wright <dwright@slac.stanford.edu>
         31 January 2017

 See header file for class documentation.

*/
//____________________________________________________________________________

#include <cstdlib>
#include <sstream>
#include <exception>

#include <TMath.h>

//rwh//#include "Framework/Algorithm/AlgConfigPool.h"
//rwh//#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Physics/HadronTransport/HG4BertCascIntranuke.h"
//rwh//#include "Physics/HadronTransport/INukeHadroData.h"
#include "Physics/HadronTransport/INukeUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"

// Geant4 headers
#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4LeptonConstructor.hh"
#include "Geant4/G4MesonConstructor.hh"
#include "Geant4/G4GenericIon.hh"
#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Fancy3DNucleus.hh"
#include "Geant4/G4InuclParticle.hh"
#include "Geant4/G4InuclCollider.hh"
#include "Geant4/G4InuclElementaryParticle.hh"
#include "Geant4/G4InuclNuclei.hh"
#include "Geant4/G4KineticTrackVector.hh"


using std::ostringstream;

using namespace genie;
using namespace genie::utils;
using namespace genie::utils::intranuke;
using namespace genie::constants;
using namespace genie::controls;

//___________________________________________________________________________
HG4BertCascIntranuke::HG4BertCascIntranuke()
: EventRecordVisitorI("genie::HG4BertCascIntranuke")
{
  InitG4Particles();
}

//___________________________________________________________________________
HG4BertCascIntranuke::HG4BertCascIntranuke(string /* config */)
: EventRecordVisitorI("genie::HG4BertCascIntranuke")
{
  InitG4Particles();
}

//___________________________________________________________________________
HG4BertCascIntranuke::~HG4BertCascIntranuke()
{
}
//___________________________________________________________________________
int HG4BertCascIntranuke::G4BertCascade(GHepRecord * evrec) const{
  GHepParticle* probe = evrec->Probe();                      // incoming particle
  GHepParticle* tgtNucl = evrec->TargetNucleus();             // target

  //
  G4ParticleDefinition* incidentDef = PDGtoG4Particle(probe->Pdg() );

  int Zinit = tgtNucl->Z();
  int Ainit = tgtNucl->A();

  G4InuclNuclei *theNucleus = new G4InuclNuclei(Ainit,Zinit);


  //rwh: unused //const  G4ParticleDefinition* tgtDef = theNucleus->getDefinition();

  const G4ThreeVector tgtmomentum(0.,0.,0.);
  const G4double tgtEnergy =tgtNucl->Energy()*1000;
  G4LorentzVector tgt4momentum(tgtmomentum,tgtEnergy);

  //
  TLorentzVector  *p4Genie = probe->P4();
  //

  G4ThreeVector incidentDir(p4Genie->Vect().Unit().Px(),p4Genie->Vect().Unit().Py(),p4Genie->Vect().Unit().Pz());

  double dynamicMass = std::sqrt(p4Genie->M2() );
  double incidentKE = p4Genie->E() - dynamicMass;
  const G4DynamicParticle p(incidentDef, incidentDir, incidentKE/units::MeV, dynamicMass/units::MeV);
  G4InuclElementaryParticle* incident = new G4InuclElementaryParticle(p,G4InuclParticle::bullet);

    // Set up output and start the cascade
  G4CollisionOutput cascadeOutput;
  G4InuclCollider bertCollider;
  //collide
  bertCollider.collide(incident,theNucleus,cascadeOutput);
  bertCollider.useCascadeDeexcitation();
  delete incident;
  delete theNucleus;

//
// Add Geant4 generated particles to the event record
//
  TLorentzVector remX(0.,0.,0.,0.);
  int Nfrag = cascadeOutput.numberOfOutgoingNuclei();
  const std::vector<G4InuclNuclei>& outgoingFragments = cascadeOutput.getOutgoingNuclei();
  int npdg = 0;
  fRemnZ = 0;
  fRemnA = 0;

  if (Nfrag > 0) {
    // Get largest nuclear fragment in output and call it the remnant
    int maxA = 0;
    int rem_index = 0;
    for (int j = 0; j < Nfrag; j++) {
      if (outgoingFragments[j].getA() > maxA) {
        maxA = outgoingFragments[j].getA();
        rem_index = j;
      }
    }

    fRemnZ = outgoingFragments[rem_index].getZ();
    fRemnA = outgoingFragments[rem_index].getA();

    // Get remnant momentum from cascade
    G4LorentzVector nucP = outgoingFragments[rem_index].getMomentum();
    TLorentzVector remP(nucP.px(), nucP.py(), nucP.pz(), nucP.e() );
    npdg = outgoingFragments[rem_index].getDefinition()->GetPDGEncoding();
    GHepParticle largest_Fragment(npdg, kIStStableFinalState,
     1,0,-1,-1, remP, remX);
    evrec->AddParticle(largest_Fragment);

    // If any nuclear fragments left, add them to the event
    for (G4int k = 0; k < Nfrag; k++) {
      if (k != rem_index) {
        npdg = outgoingFragments[k].getDefinition()->GetPDGEncoding();
        nucP = outgoingFragments[k].getMomentum();  // need to boost by fRemnP4
        TLorentzVector tempP(nucP.px(), nucP.py(), nucP.pz(), nucP.e() );
        GHepParticle nuclear_Fragment(npdg, kIStStableFinalState, 0, 1,-1,-1,
          tempP, remX);
        evrec->AddParticle(nuclear_Fragment);
      }
    }
  } // Nfrag > 0

  // Now the single hadrons
  int Nhad = cascadeOutput.numberOfOutgoingParticles();
  const std::vector<G4InuclElementaryParticle>& outgoingHadrons = cascadeOutput.getOutgoingParticles();
  for (int l = 0; l < Nhad; l++) {
    npdg = outgoingHadrons[l].getDefinition()->GetPDGEncoding();
    G4LorentzVector hadP = outgoingHadrons[l].getMomentum();
    TLorentzVector tempP(hadP.px(), hadP.py(), hadP.pz(), hadP.e() );
    GHepParticle new_particle(npdg, kIStStableFinalState, 0, 1,-1,-1,
      tempP, remX);
    evrec->AddParticle(new_particle);
  }

  return 0; // rwh: not sure what return code should imply
            // usage seems to discard any value
}
//___________________________________________________________________________
void HG4BertCascIntranuke::ProcessEventRecord(GHepRecord* evrec) const
{
  // Check the event generation mode: should be lepton-nucleus
  fGMode = evrec->EventGenerationMode();
  if ( fGMode== kGMdLeptonNucleus) {
    LOG("HG4BertCascIntranuke", pINFO)
      << " Lepton-nucleus event generation mode ";
    GHepParticle* nucltgt = evrec->TargetNucleus();
    if (nucltgt) {
    // Decide tracking radius for the current nucleus (few * R0 * A^1/3)
      SetTrackingRadius(nucltgt);
    } else {
      LOG("HG4BertCascIntranuke", pINFO) << "No nuclear target found - exit ";
      return;
    }

    // Collect hadrons from initial interaction and track them through the
    // nucleus using Bertini cascade
    TransportHadrons(evrec);
  } else if ( fGMode == kGMdHadronNucleus ||
              fGMode == kGMdPhotonNucleus    ) {
    G4BertCascade(evrec);
  } else{
    LOG("HG4BertCascIntranuke", pINFO) << " Inappropriate event generation mode - exit ";
    return;
  }

  LOG("HG4BertCascIntranuke", pINFO) << "Done with this event";
}

//___________________________________________________________________________
void HG4BertCascIntranuke::InitG4Particles() const
{
  G4LeptonConstructor::ConstructParticle();
  G4MesonConstructor::ConstructParticle();
  G4He3::He3();
  G4Gamma::Gamma();
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  pTable->SetReadiness(true);
  G4GenericIon* gIon = G4GenericIon::Definition();
  gIon->SetProcessManager(new G4ProcessManager(gIon) );
}

//___________________________________________________________________________
void HG4BertCascIntranuke::GenerateVertex(GHepRecord * /* evrec */) const
{
}

//___________________________________________________________________________
void HG4BertCascIntranuke::SetTrackingRadius(const GHepParticle * p) const
{
  assert(p && pdg::IsIon(p->Pdg()));
  double A = p->A();
  fTrackingRadius = fR0*TMath::Power(A, 1./3.);

  // multiply that by some input factor so that hadrons are tracked
  // beyond the nuclear 'boundary' since the nuclear density distribution
  // is not zero there
  fTrackingRadius *= fNR;

  LOG("HG4BertCascIntranuke", pNOTICE)
  << "Setting tracking radius to R = " << fTrackingRadius;
}

//___________________________________________________________________________
bool HG4BertCascIntranuke::IsInNucleus(const GHepParticle * /* p */) const
{
  // check whether the input particle is still within the nucleus
  //
  return false;
  //  return (p->X4()->Vect().Mag() < fTrackingRadius + fHadStep);
}
//___________________________________________________________________________
void HG4BertCascIntranuke::TransportHadrons(GHepRecord * evrec) const
{
  // Set up nucleus through which hadrons will be propagated
  // In frozen nucleus approximation it is the remnant nucleus corrected for
  // final state lepton charge

  GHepParticle* probe = evrec->Probe();                      // incoming neutrino
  GHepParticle* tgtNucl = evrec->TargetNucleus();
  GHepParticle* remNucl = evrec->RemnantNucleus();
  GHepParticle* outLept = evrec->FinalStatePrimaryLepton();  // outgoing lepton
  GHepParticle* struckNucleon = evrec->HitNucleon();

//  double sepEnergy = tgtNucl->Mass() - remNucl->Mass() - struckNucleon->Mass();
//  LOG("HG4BertCascIntranuke", pINFO)
//     << " sepEnergy = " << sepEnergy;

  // Get outgoing nucleon from initial vertex
  // It is the incident particle for rescattering
  GHepParticle* p = 0;
  G4ParticleDefinition* incidentDef = 0;
  GHepParticle* incidentBaryon = 0;
  TObjArrayIter piter(evrec);
  while( (p = (GHepParticle*) piter.Next()) ) {
    if (p->Status() == kIStHadronInTheNucleus) {
      if (p->Pdg() == kPdgProton || p->Pdg() == kPdgNeutron || p->Pdg() == kPdgLambda ||
        p->Pdg() == kPdgSigmaP || p->Pdg() == kPdgSigma0  || p->Pdg() == kPdgSigmaM) {
        incidentBaryon = p;
      incidentDef = PDGtoG4Particle(p->Pdg() );
      break;
      }
  }
}
if (!incidentBaryon) LOG("HG4BertCascIntranuke::TransportHadrons", pWARN)
 << "Unrecognized baryon in nucleus";

int Zinit = remNucl->Z() - outLept->Charge()/3;

//  if (incidentBaryon->Pdg() != struckNucleon->Pdg() ) Zinit--;
Zinit += (struckNucleon->Charge() - incidentBaryon->Charge() )/3;
int Ainit = remNucl->A();
 LOG("HG4BertCascIntranuke", pINFO) << " Zinit = " << Zinit
                                    << " Ainit = " << Ainit;

G4Fancy3DNucleus* g4Nucleus = new G4Fancy3DNucleus();
g4Nucleus->Init(Ainit, Zinit);

double EE = struckNucleon->E() - tgtNucl->Mass() + g4Nucleus->GetMass()*units::MeV;
TLorentzVector struckMomentum(struckNucleon->Px(), struckNucleon->Py(), struckNucleon->Pz(), EE);

TLorentzVector pIncident = *(tgtNucl->P4()) - *(remNucl->P4()) + *(probe->P4()) - *(outLept->P4()) - struckMomentum;
G4ThreeVector incidentDir(pIncident.Vect().Unit().Px(),
  pIncident.Vect().Unit().Py(),
  pIncident.Vect().Unit().Pz() );

double dynamicMass = std::sqrt(pIncident.M2() );
double incidentKE = pIncident.E() - dynamicMass;

  // Create pseudo-particle to supply Bertini collider with bullet
G4DynamicParticle dp(incidentDef, incidentDir, incidentKE/units::MeV, dynamicMass/units::MeV);
G4InuclElementaryParticle* incident = new G4InuclElementaryParticle(dp,G4InuclParticle::bullet);

  // Get hadronic secondaries and convert them to G4KineticTracks
G4KineticTrackVector* g4secondaries = ConvertGenieSecondariesToG4(evrec);
int Nsec = g4secondaries->size();

  // Set up output and start the cascade
G4CollisionOutput cascadeOutput;
G4InuclCollider bertCollider;
bertCollider.useCascadeDeexcitation();
  // bertCollider.setVerboseLevel(3);
bertCollider.rescatter(incident, g4secondaries, g4Nucleus, cascadeOutput);
delete incident;
delete g4Nucleus;
for (int n = 0; n < Nsec; n++) delete (*g4secondaries)[n];
  delete g4secondaries;

  //
  // Add Geant4 generated particles to the event record
  //
TLorentzVector remX(tgtNucl->Vx(), tgtNucl->Vy(), tgtNucl->Vz(), tgtNucl->Vt() );

  int rem_nucl = evrec->RemnantNucleusPosition();  // position in array
  int Nfrag = cascadeOutput.numberOfOutgoingNuclei();
  const std::vector<G4InuclNuclei>& outgoingFragments = cascadeOutput.getOutgoingNuclei();
  int npdg = 0;
  fRemnZ = 0;
  fRemnA = 0;

  if (Nfrag > 0) {
    // Get largest nuclear fragment in output and call it the remnant
    int maxA = 0;
    int rem_index = 0;
    for (int j = 0; j < Nfrag; j++) {
      if (outgoingFragments[j].getA() > maxA) {
        maxA = outgoingFragments[j].getA();
        rem_index = j;
      }
    }

    fRemnZ = outgoingFragments[rem_index].getZ();
    fRemnA = outgoingFragments[rem_index].getA();

    // Get remnant momentum from cascade
    G4LorentzVector nucP = outgoingFragments[rem_index].getMomentum();
    TLorentzVector remP(nucP.px(), nucP.py(), nucP.pz(), nucP.e() );
    npdg = outgoingFragments[rem_index].getDefinition()->GetPDGEncoding();
    GHepParticle largest_Fragment(npdg, kIStStableFinalState,
     rem_nucl,-1,-1,-1, remP, remX);
    evrec->AddParticle(largest_Fragment);

    // If any nuclear fragments left, add them to the event
    for (G4int k = 0; k < Nfrag; k++) {
      if (k != rem_index) {
        npdg = outgoingFragments[k].getDefinition()->GetPDGEncoding();
        nucP = outgoingFragments[k].getMomentum();  // need to boost by fRemnP4
        TLorentzVector tempP(nucP.px(), nucP.py(), nucP.pz(), nucP.e() );
        GHepParticle nuclear_Fragment(npdg, kIStStableFinalState, rem_nucl, -1,-1,-1,
          tempP, remX);
        evrec->AddParticle(nuclear_Fragment);
      }
    }
  } // Nfrag > 0

  // Now the single hadrons
  int Nhad = cascadeOutput.numberOfOutgoingParticles();
  const std::vector<G4InuclElementaryParticle>& outgoingHadrons = cascadeOutput.getOutgoingParticles();
  for (int l = 0; l < Nhad; l++) {
    npdg = outgoingHadrons[l].getDefinition()->GetPDGEncoding();
    G4LorentzVector hadP = outgoingHadrons[l].getMomentum();
    TLorentzVector tempP(hadP.px(), hadP.py(), hadP.pz(), hadP.e() );
    GHepParticle new_particle(npdg, kIStStableFinalState, rem_nucl, -1,-1,-1,
      tempP, remX);
    evrec->AddParticle(new_particle);
  }

  // Geant4 conservation test
  if (Conserve4Momentum(evrec) ) {
    LOG("HG4BertCascIntranuke", pINFO)
      << " Energy conservation test ";
  }
}

//___________________________________________________________________________
G4ParticleDefinition* HG4BertCascIntranuke::PDGtoG4Particle(int pdg) const
{
  G4ParticleDefinition* pDef = 0;

  switch(pdg) {
    case kPdgNuE: pDef = G4NeutrinoE::Definition(); break;
    case kPdgAntiNuE: pDef = G4AntiNeutrinoE::Definition(); break;
    case kPdgNuMu: pDef = G4NeutrinoMu::Definition(); break;
    case kPdgAntiNuMu: pDef = G4AntiNeutrinoMu::Definition(); break;
    case kPdgNuTau: pDef = G4NeutrinoTau::Definition(); break;
    case kPdgAntiNuTau: pDef = G4AntiNeutrinoTau::Definition(); break;

    case kPdgMuon: pDef = G4MuonMinus::Definition(); break;

    case kPdgProton: pDef = G4Proton::Definition(); break;
    case kPdgNeutron: pDef = G4Neutron::Definition(); break;
    case kPdgPiP: pDef = G4PionPlus::Definition(); break;
    case kPdgPiM: pDef = G4PionMinus::Definition(); break;
    case kPdgPi0: pDef = G4PionZero::Definition(); break;
    case kPdgGamma: pDef = G4Gamma::Definition(); break;
    case kPdgKP: pDef = G4KaonPlus::Definition(); break;
    case kPdgKM: pDef = G4KaonMinus::Definition(); break;
    case kPdgK0: pDef = G4KaonZero::Definition(); break;
    case kPdgAntiK0: pDef = G4AntiKaonZero::Definition(); break;
    case kPdgK0L:
    case kPdgK0S: {
      pDef = G4KaonZero::Definition();
      if (RandomGen::Instance()->RndFsi().Rndm() > 0.5)
        pDef = G4AntiKaonZero::Definition();
      break;
    }
    case kPdgLambda: pDef = G4Lambda::Definition(); break;
    case kPdgSigmaP: pDef = G4SigmaPlus::Definition(); break;
    case kPdgSigma0: pDef = G4SigmaZero::Definition(); break;
    case kPdgSigmaM: pDef = G4SigmaMinus::Definition(); break;
    case kPdgXi0: pDef = G4XiZero::Definition(); break;
    case kPdgXiM: pDef = G4XiMinus::Definition(); break;
    case kPdgOmegaM: pDef = G4OmegaMinus::Definition(); break;
    default:
    LOG("HG4BertCascIntranuke::PDGtoG4Particle", pWARN)
    << "Unrecognized Bertini particle type: " << pdg;
  }

  return pDef;
}

//___________________________________________________________________________
G4KineticTrackVector* HG4BertCascIntranuke::ConvertGenieSecondariesToG4(GHepRecord* evrec) const
{
  // Get hadronic secondaries from event record and convert them to
  // G4KineticTracks for use in cascade

  // distance units: Geant4 mm = 1.0, GENIE fermi = 1.0 => conversion 1e-12
  // nuclear radius parameter R0: Geant4 2.81967*1.2, GENIE = 1.4 => conversion 2.41686
  static double GeToG4length = 2.81967*1.0e-12*1.2/1.4;

  TObjArrayIter piter(evrec);
  GHepParticle* p = 0;
  G4ParticleDefinition* pDef = 0;
  G4KineticTrackVector* g4secondaries = new G4KineticTrackVector;
  G4KineticTrack* kt = 0;

  while( (p = (GHepParticle*) piter.Next()) ) {
    if (p->Status() == kIStHadronInTheNucleus) {
      pDef = PDGtoG4Particle(p->Pdg() );
      double formationTime = p->Vt();
      G4ThreeVector formationPosition(p->Vx()*GeToG4length,
        p->Vy()*GeToG4length,
        p->Vz()*GeToG4length);
      G4LorentzVector formationMomentum(p->Px()/units::MeV, p->Py()/units::MeV,
        p->Pz()/units::MeV, p->E()/units::MeV );

      kt = new G4KineticTrack(pDef, formationTime,
        formationPosition, formationMomentum);
      kt->SetDefinition(pDef);    // defeat reset to K0L/K0S in ctor
      kt->SetState(G4KineticTrack::inside);
      g4secondaries->push_back(kt);
    }
  }

  return g4secondaries;
}

//___________________________________________________________________________
bool HG4BertCascIntranuke::Conserve4Momentum(GHepRecord* evrec) const
{
  GHepParticle* probe = evrec->Probe(); // incoming neutrino
  GHepParticle* targetNucleus = evrec->TargetNucleus();
  GHepParticle* finalLepton = evrec->FinalStatePrimaryLepton();
  // rwh: unused //GHepParticle* remnantNucleus = evrec->RemnantNucleus();

  int Nentries = evrec->GetEntries();

  // Collect 4-momenta of final state particles
  GHepParticle* p = 0;
  TLorentzVector sum4mom(0.0, 0.0, 0.0, 0.0);
  for (int j = 0; j < Nentries; j++) {
    p = (GHepParticle*) (*evrec)[j];
    if (p->FirstMother() == evrec->RemnantNucleusPosition() ) {
      sum4mom += *(p->P4() );
      LOG("HG4BertCascIntranuke", pINFO)
        << " Final state 4-momentum = ("
        << p->P4()->Px() << ", " << p->P4()->Py() << ", "
        << p->P4()->Pz() << ", " << p->P4()->E() << ") ";
    }
  }
  sum4mom += *(finalLepton->P4() );
  LOG("HG4BertCascIntranuke", pINFO)
    << " Total Final state 4-momentum = (" << sum4mom.Px()
    << ", " << sum4mom.Py()
    << ", " << sum4mom.Pz()
    << ", " << sum4mom.E() << ") ";

  TLorentzVector initial4mom = *(targetNucleus->P4() ) + *(probe->P4() );
  LOG("HG4BertCascIntranuke", pINFO)
    << " Total Initial state 4-momentum = (" << initial4mom.Px()
    << ", " << initial4mom.Py()
    << ", " << initial4mom.Pz()
    << ", " << initial4mom.E() << ") ";

  // Charge conservation
  double Qinit = targetNucleus->Charge();
  double Qfinal = finalLepton->Charge();

  for (int j = 0; j < Nentries; j++) {
    p = (GHepParticle*) (*evrec)[j];
    if (p->FirstMother() == evrec->RemnantNucleusPosition() ) {
      Qfinal += p->Charge();
    }
  }

  if (Qinit != Qfinal) {
    LOG("HG4BertCascIntranuke", pINFO)
      << " Charge not conserved: \n"
      << " Qfinal = " << Qfinal << '\n'
      << " Qinit = " << Qinit;
  }

  return false;
}

//___________________________________________________________________________
void HG4BertCascIntranuke::LoadConfig(void)
{
  /* was for R-2
  AlgConfigPool* confp = AlgConfigPool::Instance();
  const Registry* gc = confp->GlobalParameterList();

  fR0 = fConfig->GetDoubleDef("R0", gc->GetDouble("NUCL-R0"));  // fm
  fNR = fConfig->GetDoubleDef("NR", gc->GetDouble("NUCL-NR"));
  */

  GetParam("NUCL-R0", fR0); // fm
  GetParam("NUCL-NR", fNR);

  // report
  LOG("HG4BertCascIntranuke", pINFO) << "Settings for INTRANUKE mode: " << INukeMode::AsString(kIMdHA);
  LOG("HG4BertCascIntranuke", pINFO) << "R0          = " << fR0 << " fermi";
  LOG("HG4BertCascIntranuke", pINFO) << "NR          = " << fNR;
}

//___________________________________________________________________________

void HG4BertCascIntranuke::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  LoadConfig();
}

//___________________________________________________________________________
void HG4BertCascIntranuke::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  LoadConfig();
}

//___________________________________________________________________________
#endif // __GENIE_GEANT4_INTERFACE_ENABLED__
