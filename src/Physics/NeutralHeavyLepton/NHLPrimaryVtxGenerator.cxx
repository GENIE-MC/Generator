//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
 John Plows <komninos-john.plows \at physics.ox.ac.uk>
 University of Oxford
*/
//____________________________________________________________________________

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/GMode.h"
#include "Framework/Interaction/Target.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"
#include "Physics/NuclearState/NuclearUtils.h"
#include "Physics/NeutralHeavyLepton/NHLPrimaryVtxGenerator.h"
#include "Physics/NeutralHeavyLepton/NHLDecayUtils.h"
#include "Physics/NeutralHeavyLepton/NHLDecayMode.h"

using namespace genie;
using namespace genie::NHL;

//____________________________________________________________________________
NHLPrimaryVtxGenerator::NHLPrimaryVtxGenerator() :
EventRecordVisitorI("genie::NHLPrimaryVtxGenerator")
{

}
//____________________________________________________________________________
NHLPrimaryVtxGenerator::NHLPrimaryVtxGenerator(string config) :
EventRecordVisitorI("genie::NHLPrimaryVtxGenerator",config)
{

}
//____________________________________________________________________________
NHLPrimaryVtxGenerator::~NHLPrimaryVtxGenerator()
{

}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::ProcessEventRecord(GHepRecord * event) const
{
  Interaction * interaction = event->Summary();

  fCurrInitStatePdg = interaction->InitState().ProbePdg();
  fCurrDecayMode = (NHLDecayMode_t) interaction->ExclTag().DecayMode();

  LOG("NHL", pNOTICE)
    << "Simulating NHL decay " << utils::nhl::AsString(fCurrDecayMode)
    << " for an initial state with PDG code " << fCurrInitStatePdg;

  this->AddInitialState(event);
  //this->GenerateDecayPosition(event);
  this->GenerateDecayProducts(event);
  this->UpdateEventRecord(event);

}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::AddInitialState(GHepRecord * event) const
{
  if( fUe42 == -1.0 && fUm42 == -1.0 && fUt42 == -1.0 ){
    LOG( "NHL", pINFO )
      << "Setting couplings to (1,1,0).";
    fUe42 = 1.0;
    fUm42 = 1.0;
    fUt42 = 0.0;
  }

  std::vector< double > * prodVtx = this->GenerateDecayPosition( event );
  std::vector< double > * p3NHL = this->GenerateMomentum( event );

  // RETHERE don't sample production vtx if user isn't asking for geom! It's pointless.
  TLorentzVector v4( prodVtx->at(0), prodVtx->at(1), prodVtx->at(2), 0.0 );

  Interaction * interaction = event->Summary();

  double px = p3NHL->at(0);
  double py = p3NHL->at(1);
  double pz = p3NHL->at(2);
  double E = interaction->InitState().ProbeE(kRfLab);
  
  TLorentzVector p4( px, py, pz, E );

  LOG( "NHL", pDEBUG )
    << "Probe p4 = ( " << px << ", " << py << ", " << pz << ", " << E << " )";

  InitialState * init_state = interaction->InitStatePtr();
  init_state->SetProbeP4( p4 );
  
  //init_state->SetTgtP4( v4 );
  SetProdVtxPosition( v4 );

  int hpdg = interaction->InitState().ProbePdg();
  event->AddParticle(hpdg, kIStInitialState, 0,-1,-1,-1, p4, v4);
}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::GenerateDecayProducts(GHepRecord * event) const
{
  LOG("NHL", pINFO) << "Generating decay...";

  // do we have nubar?
  PDGCodeList pdgv0 = utils::nhl::DecayProductList(fCurrDecayMode);
  int typeMod = ( fCurrInitStatePdg >= 0 ) ? 1 : -1; 
  PDGCodeList pdgv(true);
  for( std::vector<int>::iterator it = pdgv0.begin(); it != pdgv0.end(); ++it ){
    int pdgc = *it; 
    int newpdgc = ( pdgc == genie::kPdgPi0 ) ? pdgc : typeMod * pdgc; // pi-0 is its own antiparticle
    LOG("NHL", pDEBUG) << "Adding " << pdgc << " --> " << newpdgc;
    pdgv.push_back( newpdgc );
  }

  LOG("NHL", pINFO) << "Decay product IDs: " << pdgv;
  assert ( pdgv.size() > 1);

  // RETHERE may not want a phase space decay!
  LOG("NHL", pINFO) << "Performing a phase space decay...";

  // Get the decay product masses

  vector<int>::const_iterator pdg_iter;
  int i = 0;
  double * mass = new double[pdgv.size()];
  double   sum  = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
    mass[i++] = m;
    sum += m;
  }

  LOG("NHL", pINFO)
    << "Decaying N = " << pdgv.size() << " particles / total mass = " << sum;

  int nhl_id = 0;
  GHepParticle * nhl = event->Particle(nhl_id);
  assert(nhl);
  TLorentzVector * p4d = nhl->GetP4();
  TLorentzVector * v4d = nhl->GetX4();

  LOG("NHL", pINFO)
    << "Decaying system p4 = " << utils::print::P4AsString(p4d);

  // Set the decay
  bool permitted = fPhaseSpaceGenerator.SetDecay(*p4d, pdgv.size(), mass);
  if(!permitted) {
     LOG("NHL", pERROR)
       << " *** Phase space decay is not permitted \n"
       << " Total particle mass = " << sum << "\n"
       << " Decaying system p4 = " << utils::print::P4AsString(p4d);
     // clean-up
     delete [] mass;
     delete p4d;
     delete v4d;
     // throw exception
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Decay not permitted kinematically");
     exception.SwitchOnFastForward();
     throw exception;
  }

  // Get the maximum weight
  //double wmax = fPhaseSpaceGenerator.GetWtMax();
  double wmax = -1;
  for(int idec=0; idec<200; idec++) {
     double w = fPhaseSpaceGenerator.Generate();
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  wmax *= 2;

  LOG("NHL", pNOTICE)
     << "Max phase space gen. weight @ current hadronic system: " << wmax;

  // Generate an unweighted decay
  RandomGen * rnd = RandomGen::Instance();

  bool accept_decay=false;
  unsigned int itry=0;
  while(!accept_decay)
  {
     itry++;

     if(itry > controls::kMaxUnweightDecayIterations) {
       // report, clean-up and return
       LOG("NHL", pWARN)
           << "Couldn't generate an unweighted phase space decay after "
           << itry << " attempts";
       // clean up
       delete [] mass;
       delete p4d;
       delete v4d;
       // throw exception
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select decay after N attempts");
       exception.SwitchOnFastForward();
       throw exception;
     }
     double w  = fPhaseSpaceGenerator.Generate();
     if(w > wmax) {
        LOG("NHL", pWARN)
           << "Decay weight = " << w << " > max decay weight = " << wmax;
     }
     double gw = wmax * rnd->RndHadro().Rndm();
     accept_decay = (gw<=w);

     LOG("NHL", pINFO)
        << "Decay weight = " << w << " / R = " << gw
        << " - accepted: " << accept_decay;

  } //!accept_decay

  // Insert final state products into a TClonesArray of GHepParticle's
  TLorentzVector v4(*v4d);
  int idp = 0; int npip = 0, npi0 = 0, npim = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
     int pdgc = *pdg_iter;
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(idp);
     GHepStatus_t ist = kIStStableFinalState;
     event->AddParticle(pdgc, ist, nhl_id,-1,-1,-1, *p4fin, v4);
     switch( pdgc ){
     case kPdgPiP: npip++; break;
     case kPdgPi0: npi0++; break;
     case kPdgPiM: npim++; break;
     }
     idp++;
  }
  Interaction * interaction = event->Summary();
  interaction->ExclTagPtr()->SetNPions( npip, npi0, npim );

  // Manually set up some mother-daughter links
  // find primary (=leading) lepton
  double Elead = -1.0; int ilead = -1;
  std::vector< int > lpdgv = { 11, 12, 13, 14, 15, 16 };
  for( std::vector<int>::iterator lit = lpdgv.begin(); lit != lpdgv.end(); ++lit ) {
    GHepParticle * tmpPart = event->FindParticle( (*lit), kIStStableFinalState, 1 );
    if( !tmpPart ) tmpPart = event->FindParticle( -1 * (*lit), kIStStableFinalState, 1 ); //antiparticle?
    if( tmpPart ){
      double tmpE = tmpPart->E();
      if( tmpE > Elead ){ Elead = tmpE; ilead = event->ParticlePosition( tmpPart ); }
    }
  }
  event->Particle( 0 )->SetFirstDaughter( ilead );
  event->Particle( 0 )->SetFirstMother(-1); // why do I need to do this explicitly?
  event->Particle( 0 )->SetLastMother(-1);

  assert( event->Probe() );
  assert( event->FinalStatePrimaryLepton() );
  
  // loop over all FS particles and set their mother to NHL
  int itmp = 1, ilast = 1;
  while( event->Particle( itmp ) ){
    if( event->Particle(itmp)->Status() != kIStStableFinalState ){ itmp++; continue; }
    event->Particle(itmp)->SetFirstMother(0);
    event->Particle(itmp)->SetLastMother(-1);
    if( itmp != ilead ) ilast = itmp;
    itmp++;
  }
  event->Particle( 0 )->SetLastDaughter( ilast );
  // "last daughter" of NHL means last non-primary-FS-lepton, so can be less than "first" daughter

  LOG("NHL", pNOTICE)
    << "Finished with decay products. Clean up and exit!";

  // Clean-up
  delete [] mass;
  delete p4d;
  delete v4d;
}
//____________________________________________________________________________
std::vector< double > * NHLPrimaryVtxGenerator::GenerateDecayPosition( GHepRecord * event ) const
{
  // let's query *where* the NHL decayed from.
  // RETHERE - perhaps should return to GCylindTH1Flux-like implementation?
  if( !fProdVtxHist || fProdVtxHist == 0 ){
    std::string pvPath = "/GENIEv2/Generator/data/flux/HNL/HNL_vertex_positions.root"; // RETHERE - need to fix this!
    std::string pdName = "";
    std::string pvName = "hHNLVtxPos";
    fProdVtxHist = NHLFluxReader::getFluxHist3D( pvPath, pdName, pvName );
  }
  assert( fProdVtxHist );
  LOG( "NHL", pDEBUG )
    << "Found production vertex histo with " << fProdVtxHist->GetEntries() << " entries. Good!";
  
  std::vector< double > * prodVtx = NHLFluxReader::generateVtx3X( fProdVtxHist );
  LOG( "NHL", pDEBUG )
    << "Production vertex at: ( " << prodVtx->at(0) << ", " << prodVtx->at(1) << ", " << prodVtx->at(2) << ") [cm]";

  return prodVtx;
}
//____________________________________________________________________________
std::vector< double > * NHLPrimaryVtxGenerator::GenerateMomentum( GHepRecord * event ) const
{
  Interaction * interaction = event->Summary();
  double E = interaction->InitState().ProbeE(kRfLab);
  double M = PDGLibrary::Instance()->Find(kPdgNHL)->Mass();
  double p = TMath::Sqrt(E*E-M*M);

  // set some initial deviation from beam axis due to collimation effect
  double thetaDev = fAngularDeviation; // deg
  thetaDev *= genie::constants::kPi / 180.0; // rad
  RandomGen * Rng = RandomGen::Instance();
  double theta = Rng->RndGen().Gaus(0.0, thetaDev);
  if( theta < 0.0 ) theta *= -1.0;
  double phi = Rng->RndGen().Uniform(0.0, 2.0 * genie::constants::kPi);

  double px = p * std::sin(theta) * std::cos(phi);
  double py = p * std::sin(theta) * std::sin(phi);
  double pz = p * std::cos(theta);

  std::vector< double > * p3NHL = new std::vector< double >();
  p3NHL->emplace_back(px);
  p3NHL->emplace_back(py);
  p3NHL->emplace_back(pz);

  LOG( "NHL", pDEBUG )
    << "Generated momentum: ( " << px << ", " << py << ", " << pz << " )"; 

  return p3NHL;
}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::UpdateEventRecord(GHepRecord * event) const
{
  Interaction * interaction = event->Summary();

  interaction->KinePtr()->Sett( 0.0, true );
  interaction->KinePtr()->SetW( interaction->InitState().Probe()->Mass(), true );
  TLorentzVector * p4NHL = interaction->InitState().GetProbeP4( genie::kRfLab ); assert( p4NHL );
  // primary lepton is FirstDaughter() of Probe()
  // need Probe() as a GHepParticle(), not a TParticlePDG()!
  // get from event record position 0
  LOG( "NHL", pDEBUG ) << "Particle(0) has PDG code " << event->Particle(0)->Pdg();
  int iFSL = event->Particle(0)->FirstDaughter();
  LOG( "NHL", pDEBUG ) << "First daughter = " << iFSL << " with status " 
		       << (int) (event->Particle( iFSL ))->Status();
  assert( event->Particle( iFSL ) );
  TLorentzVector * p4FSL = ( event->Particle( iFSL ) )->GetP4(); 
  assert( p4FSL );
  TLorentzVector p4DIF( p4NHL->Px() - p4FSL->Px(),
			p4NHL->Py() - p4FSL->Py(),
			p4NHL->Pz() - p4FSL->Pz(),
			p4NHL->E() - p4FSL->E() );
  interaction->KinePtr()->SetQ2( p4DIF.M2(), true );

  LOG( "NHL", pDEBUG )
    << "\nNHL p4 = ( " << p4NHL->E() << ", " << p4NHL->Px() << ", " << p4NHL->Py() << ", " << p4NHL->Pz() << " )"
    << "\nFSL p4 = ( " << p4FSL->E() << ", " << p4FSL->Px() << ", " << p4FSL->Py() << ", " << p4FSL->Pz() << " )"
    << "\nDIF p4 = ( " << p4DIF.E() << ", " << p4DIF.Px() << ", " << p4DIF.Py() << ", " << p4DIF.Pz() << " )";

  // Set probe
  interaction->InitStatePtr()->SetProbePdg( event->Particle(0)->Pdg() );
  interaction->InitStatePtr()->SetProbeP4( *(event->Particle(0)->P4()) );
  
  // Set target: always Particle(1)
  // This is the charged pion in channels that have it, the pi0 in N --> pi0 pi0 v,
  // and the SM neutrino in 3-lepton channels (for N --> v v v it is the one marked "nu_e")
  interaction->InitStatePtr()->SetTgtPdg( event->Particle(1)->Pdg() );
  interaction->InitStatePtr()->SetTgtP4( *(event->Particle(1)->P4()) );
  
  LOG( "NHL", pDEBUG )
    << "Target info: " << Target();
  
  // clean up
  delete p4NHL;
  delete p4FSL;
}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void NHLPrimaryVtxGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void NHLPrimaryVtxGenerator::LoadConfig(void)
{
  LOG("NHL", pDEBUG)
    << "Loading configuration from file...";

  this->GetParam( "NHL-Mass", fMass );
  std::vector< double > U4l2s;
  this->GetParamVect( "NHL-LeptonMixing", U4l2s );
  SetNHLCouplings( U4l2s.at(0), U4l2s.at(1), U4l2s.at(2) );
  this->GetParam( "NHL-Majorana", fIsMajorana );
  this->GetParam( "NHL-Type", fType );

  this->GetParam( "NHL-angular_deviation", fAngularDeviation );

  this->GetParamVect( "Beam2User_T", fB2UTranslation );
  this->GetParamVect( "Beam2User_R", fB2URotation );
  SetBeam2User( fB2UTranslation, fB2URotation );

  fIntChannels = {}; bool itChan = false;
  // RETHERE for now I parse the channels manually... need to add automatic recognition?
  this->GetParam( "NHL-2B_mu_pi",  itChan ); if( itChan ) fIntChannels.push_back( kNHLDcyPiMu );
  this->GetParam( "NHL-2B_e_pi",   itChan ); if( itChan ) fIntChannels.push_back( kNHLDcyPiE );
  this->GetParam( "NHL-2B_nu_pi0", itChan ); if( itChan ) fIntChannels.push_back( kNHLDcyPi0Nu );
  this->GetParam( "NHL-3B_nu_nu_nu",   itChan ); if( itChan ) fIntChannels.push_back( kNHLDcyNuNuNu );
  this->GetParam( "NHL-3B_nu_mu_mu",   itChan ); if( itChan ) fIntChannels.push_back( kNHLDcyNuMuMu );
  this->GetParam( "NHL-3B_nu_e_e",     itChan ); if( itChan ) fIntChannels.push_back( kNHLDcyNuEE );
  this->GetParam( "NHL-3B_nu_mu_e",    itChan ); if( itChan ) fIntChannels.push_back( kNHLDcyNuMuE );
  this->GetParam( "NHL-3B_e_pi_pi0",   itChan ); if( itChan ) fIntChannels.push_back( kNHLDcyPiPi0E );
  this->GetParam( "NHL-3B_mu_pi_pi0",  itChan ); if( itChan ) fIntChannels.push_back( kNHLDcyPiPi0Mu );
  this->GetParam( "NHL-3B_nu_pi0_pi0", itChan ); if( itChan ) fIntChannels.push_back( kNHLDcyPi0Pi0Nu );

  LOG("NHL", pDEBUG)
    << "Read the following params:"
    << "\nMass = " << fMass << " GeV"
    << "\nECoup = " << fUe42
    << "\nMCoup = " << fUm42
    << "\nTCoup = " << fUt42
    << "\nIsMajorana = " << fIsMajorana
    << "\nType = " << fType
    << "\nAngular deviation = " << fAngularDeviation << " deg";

  fIsConfigLoaded = true;
}
//___________________________________________________________________________
void NHLPrimaryVtxGenerator::SetNHLCouplings( double Ue42, double Um42, double Ut42 ) const
{
  fUe42 = Ue42;
  fUm42 = Um42;
  fUt42 = Ut42;
}
//___________________________________________________________________________
void NHLPrimaryVtxGenerator::SetBeam2User( std::vector< double > translation, std::vector< double > rotation ) const
{
  fTx = translation.at(0);
  fTy = translation.at(1);
  fTz = translation.at(2);

  fR1 = rotation.at(0);
  fR2 = rotation.at(1);
  fR3 = rotation.at(2);

  fRM11 = std::cos( fR2 );
  fRM12 = -std::cos( fR3 ) * std::sin( fR2 );
  fRM13 = std::sin( fR2 ) * std::sin( fR3 );
  fRM21 = std::cos( fR1 ) * std::sin( fR2 );
  fRM22 = std::cos( fR1 ) * std::cos( fR2 ) * std::cos( fR3 ) - std::sin( fR1 ) * std::sin( fR3 );
  fRM23 = -std::cos( fR3 ) * std::sin( fR1 ) - std::cos( fR1 ) * std::cos( fR2 ) * std::sin( fR3 );
  fRM31 = std::sin( fR1 ) * std::sin( fR2 );
  fRM32 = std::cos( fR1 ) * std::sin( fR3 ) + std::cos( fR2 ) * std::cos( fR3 ) * std::sin( fR1 );
  fRM33 = std::cos( fR1 ) * std::cos( fR3 ) - std::cos( fR2 ) * std::sin( fR1 ) * std::sin( fR3 );

  fRTx = fTx * fRM11 + fTy * fRM12 + fTz * fRM13;
  fRTy = fTx * fRM21 + fTy * fRM22 + fTz * fRM23;
  fRTz = fTx * fRM31 + fTy * fRM32 + fTz * fRM33;

  LOG( "NHL", pDEBUG )
    << "Set BEAM origin = (0,0,0) to UNROTATED USER coordinates = ( " << fTx << ", " << fTy << ", " << fTz << " ) [m]"
    << "\nSet Euler (extrinsic x-z-x) angles to ( " << fR1 << ", " << fR2 << ", " << fR3 << " ) [rad]"
    << "\nThe rotation matrix is as follows:"
    << "\nROW 1 = ( " << fRM11 << ", " << fRM12 << ", " << fRM13 << " )"
    << "\nROW 2 = ( " << fRM21 << ", " << fRM22 << ", " << fRM23 << " )"
    << "\nROW 3 = ( " << fRM31 << ", " << fRM32 << ", " << fRM33 << " )"
    << "\nROTATED USER corrdinates = ( " << fRTx << ", " << fRTy << ", " << fRTz << " )";  
}
//___________________________________________________________________________
double NHLPrimaryVtxGenerator::GetNHLMass(string config)
{
  if( !fIsConfigLoaded ) this->Configure(config);
  return fMass;
}
//___________________________________________________________________________
std::vector< double > NHLPrimaryVtxGenerator::GetNHLCouplings(string config)
{
  if( !fIsConfigLoaded ) this->Configure(config);
  std::vector< double > coupVec = { fUe42, fUm42, fUt42 };
  return coupVec;
}
//___________________________________________________________________________
SimpleNHL NHLPrimaryVtxGenerator::GetNHLInstance(string config)
{
  if( !fIsConfigLoaded ) this->Configure(config);
  SimpleNHL sh = SimpleNHL( "NHLInstance", 0, genie::kPdgNHL, genie::kPdgKP,
			    fMass, fUe42, fUm42, fUt42, fIsMajorana );
  sh.SetType( fType );
  sh.SetInterestingChannelsVec( fIntChannels );
  sh.SetAngularDeviation( fAngularDeviation );
  sh.SetBeam2UserTranslation( fTx, fTy, fTz );
  sh.SetBeam2UserRotation( fR1, fR2, fR3 );
  return sh;
}
//____________________________________________________________________________
void NHLPrimaryVtxGenerator::SetProdVtxPosition(const TLorentzVector & v4) const
{
  TLorentzVector * pv4 = new TLorentzVector();
  pv4->SetXYZT( v4.X(), v4.Y(), v4.Z(), v4.T() );
  fProdVtx = pv4;
}
//____________________________________________________________________________
TLorentzVector * NHLPrimaryVtxGenerator::GetProdVtxPosition(void)
{
  return fProdVtx;
}
