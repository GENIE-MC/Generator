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

#include "TMath.h"

#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/HadronTransport/HG4BertCascIntranuke.h"
#include "Physics/HadronTransport/INukeUtils.h"
#include "Physics/HadronTransport/INukeUtils2018.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"


// Geant4 headers
#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4IonTable.hh"
#include "Geant4/G4LeptonConstructor.hh"
#include "Geant4/G4MesonConstructor.hh"
#include "Geant4/G4BaryonConstructor.hh"
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
#include "Geant4/G4Diproton.hh"
#include "Geant4/G4Dineutron.hh"
#include "Geant4/G4UnboundPN.hh"


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
HG4BertCascIntranuke::HG4BertCascIntranuke(string config)
: EventRecordVisitorI("genie::HG4BertCascIntranuke", config)
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


  const G4ParticleDefinition* incidentDef = PDGtoG4Particle(probe->Pdg() );

  int Zinit = tgtNucl->Z();
  int Ainit = tgtNucl->A();

  G4InuclNuclei *theNucleus = new G4InuclNuclei(Ainit,Zinit);


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


  this->SetTrackingRadius(tgtNucl);
  this->GenerateVertex(evrec);
  fRemnA=Ainit;
  fRemnZ=Zinit;

  GHepParticle * sp = new GHepParticle(*probe);
  bool has_interacted = false;
  while ( this-> IsInNucleus(sp) ) {
    // advance the hadron by a step
    utils::intranuke2018::StepParticle(sp, fHadStep);

    // check whether it interacts
    double d = this->GenerateStep(evrec,sp);
    has_interacted = (d<fHadStep);
    if(has_interacted) break;
  } //stepping

  if ( has_interacted ) {

    // Set up output and start the cascade
    G4CollisionOutput cascadeOutput;
    G4InuclCollider bertCollider;
    bertCollider.useCascadeDeexcitation();
    //collide
    bertCollider.collide(incident,theNucleus,cascadeOutput);

    delete incident;
    delete theNucleus;

    //
    // Add Geant4 generated particles to the event record
    //
    TLorentzVector remX(0.,0.,0.,0.);
    int Nfrag = cascadeOutput.numberOfOutgoingNuclei();
    const std::vector<G4InuclNuclei>& outgoingFragments =
      cascadeOutput.getOutgoingNuclei();
    int npdg = 0;
    fRemnZ = 0;
    fRemnA = 0;


    // Now the single hadrons
    int Nhad = cascadeOutput.numberOfOutgoingParticles();
    const std::vector<G4InuclElementaryParticle>& outgoingHadrons =
      cascadeOutput.getOutgoingParticles();
    for (int l = 0; l < Nhad; l++) {
      npdg = outgoingHadrons[l].getDefinition()->GetPDGEncoding();
      G4LorentzVector hadP = outgoingHadrons[l].getMomentum();
      TLorentzVector tempP(hadP.px(), hadP.py(), hadP.pz(), hadP.e() );
      GHepParticle new_particle(npdg, kIStStableFinalState, 0, 1,-1,-1,tempP, remX);
      evrec->AddParticle(new_particle);
    }

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
      //Checks if the remnant is present i PDGLibrary
      TParticlePDG * pdgRemn=PDGLibrary::Instance()->Find(npdg,false);
      
      if(!pdgRemn)
      {
        LOG("HG4BertCascIntranuke", pINFO)
        << "NO Particle with pdg = " << npdg << " in PDGLibrary!";
              // Add the particle with status id=15 and change it to HadroBlob
        GHepParticle largest_Fragment(kPdgHadronicBlob, kIStFinalStateNuclearRemnant,
          1,0,-1,-1, remP, remX);
        evrec->AddParticle(largest_Fragment);        
      }
      else
      {
        GHepParticle largest_Fragment(npdg, kIStStableFinalState,1,-1,-1,-1, remP, remX);
        evrec->AddParticle(largest_Fragment);
      }
      // If any nuclear fragments left, add them to the event
      for (G4int k = 0; k < Nfrag; k++) {
        if (k != rem_index) {
          npdg = outgoingFragments[k].getDefinition()->GetPDGEncoding();
          nucP = outgoingFragments[k].getMomentum();  // need to boost by fRemnP4
          TLorentzVector tempP(nucP.px(), nucP.py(), nucP.pz(), nucP.e() );
          GHepParticle nuclear_Fragment(npdg, kIStStableFinalState, 0, 1,-1,-1, tempP, remX);
          evrec->AddParticle(nuclear_Fragment);
        }
      }
    } // Nfrag > 0

  } else { // transparent event

    TLorentzVector p4h   (0.,0.,probe->Pz(),probe->E());
    TLorentzVector x4null(0.,0.,0.,0.);
    TLorentzVector p4tgt (0.,0.,0.,tgtNucl->Mass());
    evrec->AddParticle(probe->Pdg(),   kIStStableFinalState, 0,-1,-1,-1, p4h,x4null);
    evrec->AddParticle(tgtNucl->Pdg(), kIStStableFinalState, 1,-1,-1,-1,p4tgt,x4null);
  }
  delete sp;
  return 0;
}
//___________________________________________________________________________
double HG4BertCascIntranuke::GenerateStep(GHepRecord*  /*evrec*/, GHepParticle* p) const {
  // Generate a step (in fermis) for particle p in the input event.
  // Computes the mean free path L and generate an 'interaction' distance d
  // from an exp(-d/L) distribution

  RandomGen * rnd = RandomGen::Instance();

  double LL = utils::intranuke2018::MeanFreePath(p->Pdg(), *p->X4(), *p->P4(),
                                                 fRemnA, fRemnZ,
                                                 fDelRPion, fDelRNucleon,
                                                 fUseOset, fAltOset,
                                                 fXsecNNCorr);

  double d = -1.*LL * TMath::Log(rnd->RndFsi().Rndm());

  return d;
}
//___________________________________________________________________________
void HG4BertCascIntranuke::ProcessEventRecord(GHepRecord* evrec) const
{
  // Check the event generation mode: should be lepton-nucleus
  fGMode = evrec->EventGenerationMode();
  if ( fGMode== kGMdLeptonNucleus) {
    LOG("HG4BertCascIntranuke", pINFO) << " Lepton-nucleus event generation mode ";
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
  } else if(fGMode == kGMdHadronNucleus || fGMode == kGMdPhotonNucleus){
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
  G4BaryonConstructor::ConstructParticle();
  G4He3::He3();
  G4Gamma::Gamma();
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  pTable->SetReadiness(true);
  G4GenericIon* gIon = G4GenericIon::Definition();
  gIon->SetProcessManager(new G4ProcessManager(gIon) );

  // testing
  const G4ParticleDefinition* electron = PDGtoG4Particle(11);
  const G4ParticleDefinition* proton   = PDGtoG4Particle(2212);
  const G4ParticleDefinition* neutron  = PDGtoG4Particle(2112);
  const G4ParticleDefinition* piplus   = PDGtoG4Particle(211);
  LOG("HG4BertCascIntranuke", pINFO)
    << "testing initialization of G4 particles \n"
    << " e   0x" << electron << "\n"
    << " p   0x" << proton << "\n"
    << " n   0x" << neutron << "\n"
    << " pi+ 0x" << piplus << "\n"
    << "...InitG4Particles complete";
  if ( electron == 0 || proton == 0 || neutron == 0 || piplus == 0 ) {
    LOG("HG4BertCascIntranuke", pFATAL)
      << "something is seriously wrong with g4 particle lookup";
    exit(42);
  }
}

//___________________________________________________________________________
void HG4BertCascIntranuke::GenerateVertex(GHepRecord * evrec) const
{
  // Sets a vertex in the nucleus periphery
  // Called onlt in hadron/photon-nucleus interactions.

  GHepParticle * nucltgt = evrec->TargetNucleus();
  assert(nucltgt);

  RandomGen * rnd = RandomGen::Instance();
  TVector3 vtx(999999.,999999.,999999.);

  // *** For h+A events (test mode):
  // Assume a hadron beam with uniform intensity across an area,
  // so we need to choose events uniformly within that area.
  double x=999999., y=999999., epsilon = 0.001;
  double R2  = TMath::Power(fTrackingRadius,2.);
  double rp2 = TMath::Power(x,2.) + TMath::Power(y,2.);
  while(rp2 > R2-epsilon) {
    x = (fTrackingRadius-epsilon) * rnd->RndFsi().Rndm();
    y = -fTrackingRadius + 2*fTrackingRadius * rnd->RndFsi().Rndm();
    y -= ((y>0) ? epsilon : -epsilon);
    rp2 = TMath::Power(x,2.) + TMath::Power(y,2.);
  }
  vtx.SetXYZ(x,y, -1.*TMath::Sqrt(TMath::Max(0.,R2-rp2)) + epsilon);

  // get the actual unit vector along the incoming hadron direction
  TVector3 direction = evrec->Probe()->P4()->Vect().Unit();

  // rotate the vtx position
  vtx.RotateUz(direction);

  LOG("HG4BertCascIntranuke", pNOTICE)
    << "Generated vtx @ R = " << vtx.Mag() << " fm / "
    << print::Vec3AsString(&vtx);

  TObjArrayIter piter(evrec);
  GHepParticle * p = 0;
  while( (p = (GHepParticle *) piter.Next()) ) {
    if ( pdg::IsPseudoParticle(p->Pdg()) ) continue;
    if ( pdg::IsIon(p->Pdg())            ) continue;

    p->SetPosition(vtx.x(), vtx.y(), vtx.z(), 0.);
  }
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
bool HG4BertCascIntranuke::IsInNucleus(const GHepParticle * p) const
{
  // check whether the input particle is still within the nucleus
  return (p->X4()->Vect().Mag() < fTrackingRadius + fHadStep);
}
//___________________________________________________________________________
void HG4BertCascIntranuke::TransportHadrons(GHepRecord * evrec) const
{
  // Set up nucleus through which hadrons will be propagated
  // In frozen nucleus approximation it is the remnant nucleus corrected for
  // final state lepton charge

  GHepParticle* tgtNucl = evrec->TargetNucleus();
  GHepParticle* remNucl = evrec->RemnantNucleus();
  GHepParticle* struckNucleon = evrec->HitNucleon();

  int inucl = evrec->RemnantNucleusPosition();

  //  double sepEnergy = tgtNucl->Mass() - remNucl->Mass() - struckNucleon->Mass();
  //  std::cout << " sepEnergy = " << sepEnergy << std::endl;

  GHepParticle* p = 0;
  const G4ParticleDefinition* incidentDef = 0;
  GHepParticle* incidentBaryon = 0;
  TObjArrayIter piter(evrec);
  TObjArrayIter pitter(evrec);
  TObjArrayIter pittter(evrec);
  int icurr =-1;
  bool has_incidentBaryon(false), has_secondaries(false);
  bool has_remnant(false), has_incidentparticle(false);

  fRemnA=remNucl->A();
  fRemnZ=remNucl->Z();

  while( (p = (GHepParticle*) piter.Next()) ) {
    icurr++;
    bool has_interacted(false);
    if ( ! this->NeedsRescattering(p) ) continue;

    GHepParticle * sp = new GHepParticle(*p);
    sp->SetFirstMother(icurr);

    if ( ! this->CanRescatter(sp) ) {
      sp->SetStatus(kIStStableFinalState);
      evrec->AddParticle(*sp);
      evrec->Particle(sp->FirstMother())->SetRescatterCode(1);
      delete sp;
      continue;
    }

    while ( this-> IsInNucleus(sp) ) {
      // advance the hadron by a step
      utils::intranuke2018::StepParticle(sp, fHadStep);
      double d = this->GenerateStep(evrec,sp);
      has_interacted = (d<fHadStep);
      if ( has_interacted ) {
        has_secondaries=true;
        break;
      }
    } //stepping
    if ( ! has_interacted ) {
      sp->SetStatus(kIStStableFinalState);
      evrec->AddParticle(*sp);
      evrec->Particle(sp->FirstMother())->SetRescatterCode(1);
      //rwh-unused//transparents=true;
    }
    if ( ! has_incidentBaryon && sp->Status() == kIStHadronInTheNucleus ) {
      if ( IsBaryon(sp) ) {
        incidentBaryon = new GHepParticle(*sp);
        incidentDef = PDGtoG4Particle(sp->Pdg() );
        has_incidentBaryon=true;
      } else {
        if (sp->Pdg() == kPdgProton ||
            sp->Pdg() == kPdgNeutron||
            sp->Pdg() == kPdgLambda ||
            sp->Pdg() == kPdgSigmaP ||
            sp->Pdg() == kPdgSigma0 ||
            sp->Pdg() == kPdgSigmaM ) {
          LOG("G4BertCascInterface::TransportHadrons", pWARN)
            << "IsBaryon failed to tag " << *sp;
        }
      }
    }
    delete sp;
  }

  int Nsec(0);
  std::vector<int> Postion_evrec; 
  if ( has_secondaries ) {
    if ( ! incidentBaryon ) LOG("G4BertCascInterface::TransportHadrons", pINFO)
                              << "Unrecognized baryon in nucleus";

    delete incidentBaryon;
    G4Fancy3DNucleus* g4Nucleus = new G4Fancy3DNucleus();

    TLorentzVector pIncident;

    g4Nucleus->Init(remNucl->A(),remNucl->Z());
    double EE = struckNucleon->E() - tgtNucl->Mass() + g4Nucleus->GetMass()*units::MeV;
    TLorentzVector struckMomentum(struckNucleon->Px(), struckNucleon->Py(), struckNucleon->Pz(), EE);
    Double_t PxI(0),PyI(0),PzI(0),EEI(0), Q(0), P(0), N(0);
    int icccur=-1;
    int pos_in_evrec(0);
    while( (p = (GHepParticle*) pitter.Next()) ) {
      icccur++;
      if (p->Status() == kIStHadronInTheNucleus && this->CanRescatter(p) && p->RescatterCode()!=1) {
        PxI+=p->P4()->Px();
        PyI+=p->P4()->Py();
        PzI+=p->P4()->Pz();
        EEI+=p->P4()->E();
        Q+=p->Charge()/3;
        if ( pos_in_evrec==0 ) pos_in_evrec = icccur;
        Postion_evrec.push_back(icccur);  
        if (genie::pdg::IsProton(p->Pdg())) P++;
        if (genie::pdg::IsNeutron(p->Pdg())) N++;
        if ( pos_in_evrec==0 ) pos_in_evrec = icccur;
        if ( ! has_incidentparticle ) { // take the baryon as incident particle
          if ( IsBaryon(p) ) {
            incidentDef = PDGtoG4Particle(p->Pdg() );
            has_incidentparticle=true;
          }
        }
      }
    }
    if ( ! has_incidentparticle) {
      GHepParticle * pinN = evrec->Particle(pos_in_evrec);
      incidentDef=PDGtoG4Particle(pinN->Pdg()); // if no baryon among the secondaries
    }
    if (P == 2) incidentDef = PDGtoG4Particle(kPdgClusterPP);
    else if (N == 2) incidentDef = PDGtoG4Particle(kPdgClusterNN);
    else if (P == 1 && N == 1) incidentDef = PDGtoG4Particle(kPdgClusterNP);
    
    
    pIncident.SetPxPyPzE(PxI,PyI,PzI,EEI);


    G4ThreeVector incidentDir(pIncident.Vect().Unit().Px(),
                              pIncident.Vect().Unit().Py(),
                              pIncident.Vect().Unit().Pz());

    double dynamicMass = std::sqrt(pIncident.M2() );
    double incidentKE = pIncident.E() - dynamicMass;
    // Create pseudo-particle to supply Bertini collider with bullet

    G4DynamicParticle dp(incidentDef, incidentDir, incidentKE/units::MeV, dynamicMass/units::MeV);
    dp.SetCharge(Q);

    G4InuclElementaryParticle* incident = new G4InuclElementaryParticle(dp,G4InuclParticle::bullet);

    // Get hadronic secondaries and convert them to G4KineticTracks

    G4KineticTrackVector* g4secondaries = ConvertGenieSecondariesToG4(evrec);

    Nsec = g4secondaries->size();

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


    // Now the single hadrons
    int Nhad = cascadeOutput.numberOfOutgoingParticles();
    const std::vector<G4InuclElementaryParticle>& outgoingHadrons = cascadeOutput.getOutgoingParticles();

    int mother1=Postion_evrec.at(0);
    int mother2(0);
    if(Nsec==1)mother2=-1;
    else if(Nsec>1)mother2=Postion_evrec.at(Nsec-1);
    for (int l = 0; l < Nhad; l++) {
      npdg = outgoingHadrons[l].getDefinition()->GetPDGEncoding();

      G4LorentzVector hadP = outgoingHadrons[l].getMomentum();
      TLorentzVector tempP(hadP.px(), hadP.py(), hadP.pz(), hadP.e() );

      GHepParticle new_particle(npdg, kIStStableFinalState,mother1, mother2,-1,-1,tempP, remX);
      evrec->AddParticle(new_particle);
    }

    if (Nfrag > 0) {
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

      // If any nuclear fragments left, add them to the event
      for (G4int k = 0; k < Nfrag; k++) {
        if (k != rem_index) {
          npdg = outgoingFragments[k].getDefinition()->GetPDGEncoding();
          nucP = outgoingFragments[k].getMomentum();  // need to boost by fRemnP4
          TLorentzVector tempP(nucP.px(), nucP.py(), nucP.pz(), nucP.e() );

          GHepParticle nuclear_Fragment(npdg, kIStStableFinalState, rem_nucl,-1,-1,-1,tempP, remX);
          evrec->AddParticle(nuclear_Fragment);
        }
      }

      // Get largest nuclear fragment in output and call it the remnant
      npdg = outgoingFragments[rem_index].getDefinition()->GetPDGEncoding();
      remP.SetPx(remP.Px()+remNucl->P4()->Px());
      remP.SetPy(remP.Py()+remNucl->P4()->Py());
      remP.SetPz(remP.Pz()+remNucl->P4()->Pz());

      //Checks if the remnant is present i PDGLibrary
      TParticlePDG * pdgRemn=PDGLibrary::Instance()->Find(npdg,false);
      if(!pdgRemn)
      {
        LOG("HG4BertCascIntranuke", pINFO)
        << "NO Particle with pdg = " << npdg << " in PDGLibrary!";
        // Add the particle with status id=15 and change it to HadroBlob
        GHepParticle largest_Fragment(kPdgHadronicBlob, kIStFinalStateNuclearRemnant,
          rem_nucl,-1,-1,-1, remP, remX);
        evrec->AddParticle(largest_Fragment);
      }
      else
      {
        GHepParticle largest_Fragment(npdg, kIStStableFinalState,rem_nucl,-1,-1,-1, remP, remX);
        evrec->AddParticle(largest_Fragment);
      }
    } // Nfrag > 0
    has_remnant=true;
  }
  // Mark the initial remnant nucleus as an intermediate state
  if(!has_remnant){
    GHepParticle * sp = new GHepParticle(*evrec->Particle(inucl));
    sp->SetFirstMother(inucl);
    sp->SetStatus(kIStStableFinalState);
    evrec->AddParticle(*sp);
    delete sp;
  }
  evrec->Particle(inucl)->SetStatus(kIStIntermediateState);
  // Tests
  int dau1(0), dau2(0);
  if(Nsec>1){
      GHepParticle * pinN = evrec->Particle(Postion_evrec.at(0));
      dau1=pinN->FirstDaughter();
      dau2=pinN->LastDaughter(); 
    for(int ii=1;ii<Nsec;ii++){
      evrec->Particle(Postion_evrec.at(ii))->SetFirstDaughter(dau1);
      evrec->Particle(Postion_evrec.at(ii))->SetLastDaughter(dau2);
    }
  }


  // Geant4 conservation test
  // this probably isn't configured right ... skip it for now
  /*
  if ( Conserve4Momentum(evrec) ) {
    std::cout << " Energy conservation test " << std::endl;
  }
  */

}

//____________________________________________________________________________
bool HG4BertCascIntranuke::NeedsRescattering(const GHepParticle * p) const {
  // checks whether the particle should be rescattered
  assert(p);
  // attempt to rescatter anything marked as 'hadron in the nucleus'  (14)
  return ( p->Status() == kIStHadronInTheNucleus );
}

//___________________________________________________________________________
bool HG4BertCascIntranuke::CanRescatter(const GHepParticle * p) const
{
  assert(p);
  return  ( p->Pdg() == kPdgPiP           ||
            p->Pdg() == kPdgPiM           ||
            p->Pdg() == kPdgPi0           ||
            p->Pdg() == kPdgProton        ||
            p->Pdg() == kPdgAntiProton    ||
            p->Pdg() == kPdgNeutron       ||
            p->Pdg() == kPdgKP            ||
            p->Pdg() == kPdgKM            ||
            p->Pdg() == kPdgK0            ||
            p->Pdg() == kPdgK0L           ||
            p->Pdg() == kPdgSigma0        ||
            p->Pdg() == kPdgSigmaM        ||
            p->Pdg() == kPdgSigmaP        ||
          //p->Pdg() == kPdgSigmaPPc      ||
            p->Pdg() == kPdgXiM           ||
            p->Pdg() == kPdgXi0           ||
            p->Pdg() == kPdgLambda
            );
}

//___________________________________________________________________________
bool HG4BertCascIntranuke::IsBaryon(const GHepParticle * p) const
{
  // use PDGLibrary to determine if the particle is baryon
  assert(p);

  TParticlePDG * ppdg = PDGLibrary::Instance()->Find(p->Pdg());
  if ( ! ppdg ) {
    LOG("G4BertCascInterface", pWARN)
      << "no entry for PDG " << p->Pdg() << " in PDGLibrary";
  } else {
    if ( std::string(ppdg->ParticleClass()) == std::string("Baryon") ) {
      return true;
    }
  }
  return false;
}
//___________________________________________________________________________
const G4ParticleDefinition* HG4BertCascIntranuke::PDGtoG4Particle(int pdg) const
{
  const G4ParticleDefinition* pDef = 0;

  if (pdg == kPdgClusterPP) return G4Diproton::Diproton();
  if (pdg == kPdgClusterNN) return G4Dineutron::Dineutron();
  if (pdg == kPdgClusterNP) return G4UnboundPN::UnboundPN();
 
  if ( abs(pdg) < 1000000000 ) {
    pDef = G4ParticleTable::GetParticleTable()->FindParticle(pdg);
  } else if ( pdg < 2000000000 ) {
    pDef = G4IonTable::GetIonTable()->GetIon(pdg);
  }
  
  if ( ! pDef ) {
    LOG("HG4BertCascIntranuke", pWARN)
      << "Unrecognized Bertini particle type: " << pdg;
  }

  return pDef;
}

//___________________________________________________________________________
G4KineticTrackVector* HG4BertCascIntranuke::ConvertGenieSecondariesToG4(std::vector<GHepParticle> partList) const
{
  static double GeToG4length = 2.81967*1.0e-12*1.2/1.4;

  GHepParticle* p = 0;
  const G4ParticleDefinition* pDef = 0;
  G4KineticTrackVector* g4secondaries = new G4KineticTrackVector;
  G4KineticTrack* kt = 0;

  for (size_t it=0 ; it<partList.size(); it++) {
    p= &partList.at(it);
    pDef = PDGtoG4Particle(p->Pdg() );
    double formationTime = p->Vt();
    G4ThreeVector formationPosition(p->Vx()*GeToG4length,
                                    p->Vy()*GeToG4length,
                                    p->Vz()*GeToG4length);
    G4LorentzVector formationMomentum(p->Px()/units::MeV, p->Py()/units::MeV,
                                      p->Pz()/units::MeV, p->E()/units::MeV );
    kt = new G4KineticTrack(pDef, formationTime, formationPosition, formationMomentum);
    kt->SetDefinition(pDef);    // defeat reset to K0L/K0S in ctor
    kt->SetState(G4KineticTrack::inside);
    g4secondaries->push_back(kt);
  }

  return g4secondaries;
}

//___________________________________________________________________________
G4KineticTrackVector* HG4BertCascIntranuke::ConvertGenieSecondariesToG4(GHepRecord* evrec) const {
  // Get hadronic secondaries from event record and convert them to
  // G4KineticTracks for use in cascade

  // distance units: Geant4 mm = 1.0, GENIE fermi = 1.0 => conversion 1e-12
  // nuclear radius parameter R0: Geant4 2.81967*1.2, GENIE = 1.4 => conversion 2.41686
  static double GeToG4length = 2.81967*1.0e-12*1.2/1.4;

  TObjArrayIter piter(evrec);
  GHepParticle* p = 0;
  const G4ParticleDefinition* pDef = 0;
  G4KineticTrackVector* g4secondaries = new G4KineticTrackVector;
  G4KineticTrack* kt = 0;

  while( (p = (GHepParticle*) piter.Next()) ) {
    if ( p->Status() == kIStHadronInTheNucleus &&
         this->CanRescatter(p) && ( p->RescatterCode() != 1) ) {
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
  bool failed = false;

  GHepParticle* probe = evrec->Probe(); // incoming neutrino
  GHepParticle* targetNucleus = evrec->TargetNucleus();
  GHepParticle* finalLepton = evrec->FinalStatePrimaryLepton();

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

  TLorentzVector initial4mom = *(targetNucleus->P4() ) + *(probe->P4() );

  TLorentzVector diff = initial4mom - sum4mom;
  // rwh need some threshold for acceptable differences
  const double maxdiff  = 1.0e-9;  // set crazy small for now
  double   diffE    = diff.E();
  TVector3 diffp3   = diff.Vect();
  double   diffpmag = diffp3.Mag();
  if ( TMath::Abs(diffE) > maxdiff || diffpmag > maxdiff ) {
    failed = true;
    LOG("HG4BertCascIntranuke", pWARN)
      << "final state - initial state differs by > " << maxdiff << "\n"
      << " dE = " << diffE << ", d|p3| = " << diffpmag;

    LOG("HG4BertCascIntranuke", pWARN)
      << " Total Final state 4-momentum = (" << sum4mom.Px()
      << ", " << sum4mom.Py()
      << ", " << sum4mom.Pz()
      << ", " << sum4mom.E() << ") ";
    LOG("HG4BertCascIntranuke", pWARN)
      << " Total Initial state 4-momentum = (" << initial4mom.Px()
      << ", " << initial4mom.Py()
      << ", " << initial4mom.Pz()
      << ", " << initial4mom.E() << ") ";
  }

  // Charge conservation
  double Qinit = targetNucleus->Charge();
  double Qfinal = finalLepton->Charge();

  for (int j = 0; j < Nentries; j++) {
    p = (GHepParticle*) (*evrec)[j];
    if ( p->FirstMother() == evrec->RemnantNucleusPosition() ) {
      Qfinal += p->Charge();
    }
  }

  if (Qinit != Qfinal) {
    failed = true;
    LOG("HG4BertCascIntranuke", pWARN)
      << " Charge not conserved: \n"
      << " Qfinal = " << Qfinal
      << " Qinit = " << Qinit;
  }

  return ( ! failed );
}


//___________________________________________________________________________
void HG4BertCascIntranuke::LoadConfig(void)
{
  // fermi momentum setup
  // this is specifically set in Intranuke2018::Configure(string)
  fNuclmodel = dynamic_cast<const NuclearModelI *>( this -> SubAlg("NuclearModel") ) ;

  // other intranuke config params
  GetParam( "NUCL-R0",             fR0 );             // fm
  GetParam( "NUCL-NR",             fNR );

  GetParam( "INUKE-NucRemovalE",   fNucRmvE );        // GeV
  GetParam( "INUKE-HadStep",       fHadStep ) ;
  GetParam( "INUKE-NucAbsFac",     fNucAbsFac ) ;
  GetParam( "INUKE-NucCEXFac",     fNucCEXFac ) ;
  GetParam( "INUKE-Energy_Pre_Eq", fEPreEq ) ;
  GetParam( "INUKE-FermiFac",      fFermiFac ) ;
  GetParam( "INUKE-FermiMomentum", fFermiMomentum ) ;
  GetParam( "INUKE-DoFermi",       fDoFermi ) ;
  GetParam( "INUKE-XsecNNCorr",    fXsecNNCorr ) ;
  GetParamDef( "UseOset",          fUseOset, false ) ;
  GetParamDef( "AltOset",          fAltOset, false ) ;
  GetParam( "HNINUKE-DelRPion",    fDelRPion ) ;
  GetParam( "HNINUKE-DelRNucleon", fDelRNucleon ) ;

  // report
  LOG("HG4BertCascIntranuke", pNOTICE)
    << "Settings for INTRANUKE mode: " << INukeMode::AsString(kIMdHA) << "\n"
    << "R0          = " << fR0 << " fermi" << "\n"
    << "NR          = " << fNR << "\n"
    << "DelRPion    = " << fDelRPion << "\n"
    << "DelRNucleon = " << fDelRNucleon << "\n"
    << "HadStep     = " << fHadStep << " fermi" << "\n"
    << "EPreEq      = " << fHadStep << " fermi" << "\n"
    << "NucAbsFac   = " << fNucAbsFac << "\n"
    << "NucCEXFac   = " << fNucCEXFac << "\n"
    << "FermiFac    = " << fFermiFac << "\n"
    << "FermiMomtm  = " << fFermiMomentum << "\n"
    << "DoFermi?    = " << ((fDoFermi)?(true):(false)) << "\n"
    << "XsecNNCorr? = " << ((fXsecNNCorr)?(true):(false));

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
