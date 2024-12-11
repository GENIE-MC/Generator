//____________________________________________________________________________
/*
 Copyright (c) 2003-2024, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
*/
//____________________________________________________________________________

#include <cstring>

#include <RVersion.h>
#ifdef __GENIE_PYTHIA6_ENABLED__
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#endif
#include <TClonesArray.h>
#include <TMath.h>

#include "Physics/HELepton/EventGen/PhotonRESGenerator.h"
#include "Physics/HEDIS/EventGen/HEDISGenerator.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/StringUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::math;

//___________________________________________________________________________
PhotonRESGenerator::PhotonRESGenerator() :
EventRecordVisitorI("genie::PhotonRESGenerator")
{
  this->Initialize();
  born = new Born();
}
//___________________________________________________________________________
PhotonRESGenerator::PhotonRESGenerator(string config) :
EventRecordVisitorI("genie::PhotonRESGenerator", config)
{
  this->Initialize();
}
//___________________________________________________________________________
PhotonRESGenerator::~PhotonRESGenerator()
{

}
//____________________________________________________________________________
void PhotonRESGenerator::Initialize(void) const
{
#ifdef __GENIE_PYTHIA6_ENABLED__
  fPythia = TPythia6::Instance();

  // sync GENIE/PYTHIA6 seed number
  RandomGen::Instance();
#endif

#ifdef __GENIE_PYTHIA8_ENABLED__
  //One cool feature of having independent Pythia8 instances
  //We can define our intial state with the proper setting (i.e. disabling decays)
  //without affecting other classes of GENIE that also use Pythia8

  fPythiaP = new Pythia8::Pythia();
  fPythiaP->readString("Print:quiet = on");
  fPythiaP->readString("Random:setSeed = on");
  fPythiaP->readString("WeakSingleBoson:ffbar2ffbar(s:W) = on");
  fPythiaP->readString("PDF:lepton = off");
  fPythiaP->readString("24:onMode = off");
  fPythiaP->readString("24:onIfAny = 1 2 3 4 5"); //enable W->hadron only 
  fPythiaP->readString("Beams:idA = 12");
  fPythiaP->readString("Beams:idB = -11");

  fPythiaN = new Pythia8::Pythia();
  fPythiaN->readString("Print:quiet = on");
  fPythiaN->readString("Random:setSeed = on");
  fPythiaN->readString("WeakSingleBoson:ffbar2ffbar(s:W) = on");
  fPythiaN->readString("PDF:lepton = off");
  fPythiaN->readString("24:onMode = off");
  fPythiaN->readString("24:onIfAny = 1 2 3 4 5"); //enable W->hadron only 
  fPythiaN->readString("Beams:idA = -12");
  fPythiaN->readString("Beams:idB = 11");
#endif

}
//___________________________________________________________________________
void PhotonRESGenerator::ProcessEventRecord(GHepRecord * evrec) const
{

  Interaction * interaction = evrec->Summary();
  const InitialState & init_state = interaction->InitState();

  //incoming v & struck particle & remnant nucleus
  GHepParticle * nu = evrec->Probe();
  GHepParticle * el = evrec->HitNucleon();

  GHepParticle * target = evrec -> TargetNucleus();
  if(target) evrec->AddParticle(target->Pdg(), kIStFinalStateNuclearRemnant, 1,-1,-1,-1, *(target->P4()), *(target->X4()) );

  TVector3 unit_nu = nu->P4()->Vect().Unit();

  int probepdg = init_state.ProbePdg();

  long double Mtarget = init_state.Tgt().HitNucMass();
  long double mlout = interaction->FSPrimLepton()->Mass(); //mass of charged lepton

  long double Enuin = init_state.ProbeE(kRfLab);
  long double s = born->GetS(Mtarget,Enuin);

  long double n1 = interaction->Kine().GetKV(kKVn1);
  long double n2 = interaction->Kine().GetKV(kKVn2);

  long double costhCM = n1;
  long double sinthCM = sqrtl(1-costhCM*costhCM);

  long double xmin = fQ2PDFmin/2./Enuin/Mtarget;
  long double x = expl( logl(xmin) + (logl(1.0)-logl(xmin))*n2 );
  long double s_r = s*x;

  // Boost velocity CM -> LAB
  long double EnuinCM = sqrtl(s_r)/2.;
  long double beta = (powl(Enuin,2)-powl(EnuinCM,2))/(powl(Enuin,2)+powl(EnuinCM,2));

  // Final state primary lepton PDG code
  int pdgl = interaction->FSPrimLeptonPdg();
  assert(pdgl!=0);

  if ( pdg::IsElectron(TMath::Abs(pdgl)) || pdg::IsMuon(TMath::Abs(pdgl)) || pdg::IsTau(TMath::Abs(pdgl)) ) {

    long double ElpoutCM = (s_r+mlout*mlout)/sqrtl(s_r)/2.;
    long double EnuoutCM = (s_r-mlout*mlout)/sqrtl(s_r)/2.;
    LongLorentzVector p4_lpout( 0.,  EnuoutCM*sinthCM,  EnuoutCM*costhCM, ElpoutCM );
    LongLorentzVector p4_nuout( 0., -EnuoutCM*sinthCM, -EnuoutCM*costhCM, EnuoutCM );

    p4_lpout.BoostZ(beta);
    p4_nuout.BoostZ(beta);

    TLorentzVector p4lp_o( (double)p4_lpout.Px(), (double)p4_lpout.Py(), (double)p4_lpout.Pz(), (double)p4_lpout.E() );
    TLorentzVector p4nu_o( (double)p4_nuout.Px(), (double)p4_nuout.Py(), (double)p4_nuout.Pz(), (double)p4_nuout.E() );

    // Randomize transverse components
    RandomGen * rnd = RandomGen::Instance();
    double phi  = 2* kPi * rnd->RndLep().Rndm();
    p4lp_o.RotateZ(phi);
    p4nu_o.RotateZ(phi);

    //rotate from LAB=[0,0,Ev,Ev]->[px,py,pz,E]
    p4lp_o.RotateUz(unit_nu);
    p4nu_o.RotateUz(unit_nu);

    int pdgvout = 0;
    if      ( pdg::IsElectron(pdgl) ) pdgvout = kPdgAntiNuE;
    else if ( pdg::IsPositron(pdgl) ) pdgvout = kPdgNuE;
    else if ( pdg::IsMuon(pdgl)     ) pdgvout = kPdgAntiNuMu;
    else if ( pdg::IsAntiMuon(pdgl) ) pdgvout = kPdgNuMu;
    else if ( pdg::IsTau(pdgl)      ) pdgvout = kPdgAntiNuTau;
    else if ( pdg::IsAntiTau(pdgl)  ) pdgvout = kPdgNuTau;

    int pdgboson = pdg::IsNeutrino(probepdg) ? kPdgWP : kPdgWM;

    // Create a GHepParticle and add it to the event record
    evrec->AddParticle( pdgboson, kIStDecayedState,     0, -1,  5,  6, *nu->P4()+*el->P4(), *(nu->X4()) ); //W [mothers: nuebar_in,e_in][daugthers: nulbar_out,lp_out]
    evrec->AddParticle( pdgl,     kIStStableFinalState, 4, -1, -1, -1, p4lp_o,              *(nu->X4()) );
    evrec->AddParticle( pdgvout,  kIStStableFinalState, 4, -1, -1, -1, p4nu_o,              *(nu->X4()) );
    evrec->Summary()->KinePtr()->SetFSLeptonP4(p4lp_o);

  }
  else {

    char p6frame[10];
    strcpy(p6frame, "CMS"    );

    char p6nu[10], p6tgt[10];
    if      (pdg::IsAntiNeutrino(nu->Pdg())) { strcpy(p6nu, "nu_ebar");   strcpy(p6tgt, "e-");   }
    else if (pdg::IsNeutrino    (nu->Pdg())) { strcpy(p6nu, "nu_e");      strcpy(p6tgt, "e+");   }

#ifdef __GENIE_PYTHIA6_ENABLED__
    int def61  = fPythia->GetMSTP(61);
    int def71  = fPythia->GetMSTP(71);
    int def206 = fPythia->GetMDME(206,1);
    int def207 = fPythia->GetMDME(207,1);
    int def208 = fPythia->GetMDME(208,1);
    fPythia->SetMSTP(61,0);    // (Default=2) master switch for initial-state QCD and QED radiation.
    fPythia->SetMSTP(71,0);    // (Default=2) master switch for initial-state QCD and QED radiation.
    fPythia->SetMDME(206,1,0); //swicht off W decay leptonic modes
    fPythia->SetMDME(207,1,0);
    fPythia->SetMDME(208,1,0);

    fPythia->Pyinit(p6frame, p6nu, p6tgt, sqrtl(s_r));
    fPythia->Pyevnt();

    fPythia->SetMSTP(61,def61);
    fPythia->SetMSTP(71,def71);
    fPythia->SetMDME(206,1,def206);
    fPythia->SetMDME(207,1,def207);
    fPythia->SetMDME(208,1,def208);

    // get LUJETS record
    fPythia->GetPrimaries();
    TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");
    int np = pythia_particles->GetEntries();
    assert(np>0);

    TMCParticle * particle = 0;
    TIter piter(pythia_particles);
    while( (particle = (TMCParticle *) piter.Next()) ) {

      int pdgc = particle->GetKF();
      int ks = particle->GetKS();

      if ( ks==21 ) { continue; } //we dont want to save first particles from pythia (init states)

      LongLorentzVector p4longo(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetEnergy());
      p4longo.BoostZ(beta);

      TLorentzVector p4o( (double)p4longo.Px(), (double)p4longo.Py(), (double)p4longo.Pz(), (double)p4longo.E() );
      p4o.RotateUz(unit_nu);

      TParticlePDG * part = PDGLibrary::Instance()->Find(pdgc);
      if ( (ks==1 || ks==4) && p4o.E() < part->Mass() ) {
        LOG("PhotonRESGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
        LOG("PhotonRESGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
        LOG("PhotonRESGenerator", pWARN) << "E = " << p4o.E() << " // |p| = " << TMath::Sqrt(p4o.P());
        LOG("PhotonRESGenerator", pWARN) << "p = [ " << p4o.Px() << " , "  << p4o.Py() << " , "  << p4o.Pz() << " ]";
        LOG("PhotonRESGenerator", pWARN) << "m    = " << p4o.M() << " // mpdg = " << part->Mass();
        p4o.SetXYZT(0,0,0,part->Mass());
      }

      // copy final state particles to the event record
      GHepStatus_t ist = (ks==1 || ks==4) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

      // fix numbering scheme used for mother/daughter assignments
      int firstmother = -1;
      int lastmother  = -1;
      int firstchild  = -1;
      int lastchild   = -1;

      if ( particle->GetParent() < 10 ) {
        if      ( TMath::Abs(pdgc)<7 ) {   //outgoing quarks: mother will be the boson (saved in position 4)
          firstmother = 4;
          firstchild  = particle->GetFirstChild() - 6;
          lastchild   = particle->GetLastChild()  - 6;
        }
        else if ( TMath::Abs(pdgc)==24 ) { //produced W boson: mother will be the incoming neutrino
          firstmother = 0;
          firstchild  = particle->GetFirstChild() - 6;
          lastchild   = particle->GetLastChild()  - 6;
        }
        else if ( pdgc==22 ) {             //radiative photons: mother will be the incoming electron
          firstmother = 2;
        }
      }
      else { //rest
        firstmother = particle->GetParent()     - 6; //shift to match boson position
        firstchild  = (particle->GetFirstChild()==0) ? particle->GetFirstChild() - 1 : particle->GetFirstChild() - 6;
        lastchild   = (particle->GetLastChild()==0)  ? particle->GetLastChild()  - 1 : particle->GetLastChild()  - 6;
      }

      double lightspeed = 299792458e3; //c in mm/s. Used for time in PYTHIA t[s]=t_pythia[mm]/c[mm/s]
      double vx = nu->X4()->X() + particle->GetVx()*1e12; //pythia gives position in [mm] while genie uses [fm]
      double vy = nu->X4()->Y() + particle->GetVy()*1e12;
      double vz = nu->X4()->Z() + particle->GetVz()*1e12;
      double vt = nu->X4()->T() + particle->GetTime()/lightspeed;
      TLorentzVector pos( vx, vy, vz, vt );

      evrec->AddParticle(pdgc, ist, firstmother, lastmother, firstchild, lastchild, p4o, pos );

    }

    delete particle;
    pythia_particles->Clear("C");

#endif // __GENIE_PYTHIA6_ENABLED__

#ifdef __GENIE_PYTHIA8_ENABLED__

    RandomGen * rnd = RandomGen::Instance();

    Pythia8::Event fEvent;

    //We have to initialize Pythia8 everytime we generate a new event
    //And we have to change the seed because it would regenerate the same event
    //if we dont do it.
    if ( pdg::IsNeutrino(nu->Pdg()) ) {
      fPythiaP->settings.mode("Random:seed", rnd->RndLep().Integer(100000000));
      fPythiaP->settings.parm("Beams:eCM",sqrtl(s_r));
      fPythiaP->init(); 
      fPythiaP->next();
      // fPythiaP->event.list();
      // fPythiaP->stat();
      fEvent = fPythiaP->event;
    }
    else if ( pdg::IsAntiNeutrino(nu->Pdg()) ) {
      fPythiaN->settings.mode("Random:seed", rnd->RndLep().Integer(100000000));
      fPythiaN->settings.parm("Beams:eCM",sqrtl(s_r));
      fPythiaN->init(); 
      fPythiaN->next();
      // fPythiaN->event.list();
      // fPythiaN->stat();
      fEvent = fPythiaN->event;      
    }

    int np = fEvent.size();
    assert(np>0);

    for (int i = 5; i < np; i++) { // ignore firt entry -> system + input particles

      int pdgc = fEvent[i].id();
      if (!PDGLibrary::Instance()->Find(pdgc)) continue; // some intermidatie particles not part of genie tables

      int ks = fEvent[i].status();

      LongLorentzVector p4longo(fEvent[i].px(), fEvent[i].py(), fEvent[i].pz(), fEvent[i].e());
      p4longo.BoostZ(beta); 

      TLorentzVector p4o( (double)p4longo.Px(), (double)p4longo.Py(), (double)p4longo.Pz(), (double)p4longo.E() );
      p4o.RotateUz(unit_nu); 
     
      double massPDG = PDGLibrary::Instance()->Find(pdgc)->Mass();
      if ( ks>0 && p4o.E()<massPDG ) {
        LOG("GLRESGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
        LOG("GLRESGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
        LOG("GLRESGenerator", pWARN) << "E = " << p4o.E() << " // |p| = " << TMath::Sqrt(p4o.P()); 
        LOG("GLRESGenerator", pWARN) << "p = [ " << p4o.Px() << " , "  << p4o.Py() << " , "  << p4o.Pz() << " ]";
        LOG("GLRESGenerator", pWARN) << "m    = " << p4o.M() << " // mpdg = " << massPDG;
        p4o.SetXYZT(0,0,0,massPDG);
      }

      // copy final state particles to the event record
      GHepStatus_t ist = (ks>0) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

      // fix numbering scheme used for mother/daughter assignments
      int firstmother = -1;
      int lastmother  = -1;
      int firstchild  = -1;
      int lastchild   = -1;

      if ( fEvent[i].mother1() == 3 ) { //produced W boson: mother will be the incoming neutrino
          firstmother = 0;
          firstchild  = fEvent[i].daughter1() - 3;
          lastchild   = fEvent[i].daughter2() - 3;        
      }
      else { //rest
        firstmother = fEvent[i].mother1() - 3; //shift to match boson position
        firstchild  = fEvent[i].daughter1() - 3;
        lastchild   = fEvent[i].daughter2() - 3;
      }

      double vx = nu->X4()->X() + fEvent[i].xProd()*1e12; //pythia gives position in [mm] while genie uses [fm]
      double vy = nu->X4()->Y() + fEvent[i].yProd()*1e12;
      double vz = nu->X4()->Z() + fEvent[i].zProd()*1e12;
      double vt = nu->X4()->T() + fEvent[i].tProd()*(units::millimeter/units::second);
      TLorentzVector pos( vx, vy, vz, vt );

      evrec->AddParticle(pdgc, ist, firstmother, lastmother, firstchild, lastchild, p4o, pos );
    }

#endif

  }

}
//___________________________________________________________________________
void PhotonRESGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonRESGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonRESGenerator::LoadConfig(void)
{

  // PYTHIA parameters only valid for HEDIS
  double wmin;        GetParam( "Xsec-Wmin",     wmin ) ;
  int warnings;       GetParam( "Warnings",      warnings ) ;
  int errors;         GetParam( "Errors",        errors ) ;
  int qrk_mass;       GetParam( "QuarkMass",     qrk_mass ) ;
#ifdef __GENIE_PYTHIA6_ENABLED__
  fPythia->SetPARP(2,  wmin);         //(D = 10. GeV) lowest c.m. energy for the event as a whole that the program will accept to simulate. (bellow 2GeV pythia crashes)
  fPythia->SetMSTU(26, warnings);     // (Default=10) maximum number of warnings that are printed
  fPythia->SetMSTU(22, errors);       // (Default=10) maximum number of errors that are printed
  fPythia->SetMSTJ(93, qrk_mass);     // light (d, u, s, c, b) quark masses are taken from PARF(101) - PARF(105) rather than PMAS(1,1) - PMAS(5,1). Diquark masses are given as sum of quark masses, without spin splitting term.

  fPythia->SetPMAS(24,1,kMw);  //mass of the W boson (pythia=80.450 // genie=80.385)
  fPythia->SetPMAS(24,2,0.);   //set to 0 the width of the W boson to avoid problems with energy conservation
  fPythia->SetPMAS(6,2,0.);    //set to 0 the width of the top to avoid problems with energy conservation
#endif // __GENIE_PYTHIA6_ENABLED__

#ifdef __GENIE_PYTHIA8_ENABLED__
  // Same default mass of the W boson in pythia8 and genie, so no need to change in pythia8
  // No problem with energy conservation W and top decays, so no need to set the width to 0

  // Important for not decaying too early!!! 
  // Copy of the method from Decayer class to have a consistent treatment.
  // The idea is to not decay in this stage what we dont want to decay in Decayer.
  // Hadronization  ->   Decayer
  //------------------------------
  // Decay          ->   Decay
  // Undecay        ->   Undecay
  // Undecay        ->   Decay
  // Decay          ->   Undecay -> This is the problematic one
  // So we loop over all the particles for which we inhibit decay in the config file
  // for the Decayer stage and we inhibit it also at hadronization.
  // This is not need in LeptonHadronization and PhotonCOH because they use the 
  // Pythia8Singleton class which already defines the decays consistently
  //
  RgKeyList klist = GetConfig().FindKeys("DecayParticleWithCode=");
  RgKeyList::const_iterator kiter = klist.begin();
  for( ; kiter != klist.end(); ++kiter) {
    RgKey key = *kiter;
    bool decay = GetConfig().GetBool(key);
    vector<string> kv = utils::str::Split(key,"=");
    int pdg_code = atoi(kv[1].c_str());
    if(!decay) {
      LOG("GLRESGenerator", pDEBUG) << "Configured to inhibit decays for  " <<  pdg_code;
      
      // Pythia8 class for nue+e+>W+
      auto pdentryP = fPythiaP->particleData.particleDataEntryPtr(pdg_code);
      for (int ichan=0; ichan<=pdentryP->sizeChannels()-1; ++ichan) {
        int onMode = pdentryP->channel(ichan).onMode();
        bool is_particle = (pdg_code>0);
        switch ( onMode ) {
        case 0: // off for particles & antiparticles
          break; // already off
        case 1: // on for both particle & antiparticle
          onMode = (is_particle ? 3 : 2);
          break;
        case 2: // on for particle only
          onMode = (is_particle ? 0 : 2);
          break;
        case 3: // on for antiparticle only
          onMode = (is_particle ? 3 : 0);
          break;
        }
        pdentryP->channel(ichan).onMode(onMode);
      }

      // Pythia8 class for anue+e->W-
      auto pdentryN = fPythiaN->particleData.particleDataEntryPtr(pdg_code);
      for (int ichan=0; ichan<=pdentryN->sizeChannels()-1; ++ichan) {
        int onMode = pdentryN->channel(ichan).onMode();
        bool is_particle = (pdg_code>0);
        switch ( onMode ) {
        case 0: // off for particles & antiparticles
          break; // already off
        case 1: // on for both particle & antiparticle
          onMode = (is_particle ? 3 : 2);
          break;
        case 2: // on for particle only
          onMode = (is_particle ? 0 : 2);
          break;
        case 3: // on for antiparticle only
          onMode = (is_particle ? 3 : 0);
          break;
        }
        pdentryN->channel(ichan).onMode(onMode);
      }

    }
  }

  //Important: We dont initialize Pythia8 here because it must be done event by event. See above
#endif

  GetParam( "Q2Grid-Min", fQ2PDFmin );

}
//____________________________________________________________________________
