//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC (Valencia)
*/
//____________________________________________________________________________

#include "Physics/HELepton/EventGen/PhotonRESWdecPythia8.h"

//___________________________________________________________________________
PhotonRESWdecPythia8::PhotonRESWdecPythia8() :
EventRecordVisitorI("genie::PhotonRESWdecPythia8")
{
  this->Initialize();
  born = new Born();
}
//___________________________________________________________________________
PhotonRESWdecPythia8::PhotonRESWdecPythia8(string config) :
EventRecordVisitorI("genie::PhotonRESWdecPythia8", config)
{
  this->Initialize();
  born = new Born();
}
//___________________________________________________________________________
PhotonRESWdecPythia8::~PhotonRESWdecPythia8()
{

}
//____________________________________________________________________________
void PhotonRESWdecPythia8::ProcessEventRecord(GHepRecord * event) const
{

  if(!this->Wdecay(event)) {
    LOG("PhotonRESGenerator", pWARN) << "W decayer failed!";
    event->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Could not simulate the W decay system");
    exception.SwitchOnFastForward();
    throw exception;
    return;
  }

}
//___________________________________________________________________________
bool PhotonRESWdecPythia8::Wdecay(GHepRecord * 
#ifdef __GENIE_PYTHIA8_ENABLED__
  event // avoid unused variable warning if PYTHIA8 is not enabled
#endif
) const
{

#ifdef __GENIE_PYTHIA8_ENABLED__
  Interaction * interaction = event->Summary();
  const InitialState & init_state = interaction->InitState();

  //incoming v & struck particle & remnant nucleus
  GHepParticle * nu = event->Probe();
  GHepParticle * el = event->HitNucleon();

  GHepParticle * target = event -> TargetNucleus();
  if(target) event->AddParticle(target->Pdg(), kIStFinalStateNuclearRemnant, 1,-1,-1,-1, *(target->P4()), *(target->X4()) );

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
    event->AddParticle( pdgboson, kIStDecayedState,     0, -1,  5,  6, *nu->P4()+*el->P4(), *(nu->X4()) ); //W [mothers: nuebar_in,e_in][daugthers: nulbar_out,lp_out]
    event->AddParticle( pdgl,     kIStStableFinalState, 4, -1, -1, -1, p4lp_o,              *(nu->X4()) );
    event->AddParticle( pdgvout,  kIStStableFinalState, 4, -1, -1, -1, p4nu_o,              *(nu->X4()) );
    event->Summary()->KinePtr()->SetFSLeptonP4(p4lp_o);

  }
  else {


    //We have to initialize Pythia8 everytime we generate a new event
    //And we have to change the seed because it would regenerate the same event
    //if we dont do it.
    RandomGen * rnd = RandomGen::Instance();
    Pythia8::Event fEvent;
    if (pdg::IsNeutrino(probepdg)) {
      fPythiaP->settings.mode("Random:seed", rnd->RndLep().Integer(100000000));
      fPythiaP->settings.parm("Beams:eCM",sqrtl(s_r));
      fPythiaP->init(); 
      fPythiaP->next();
      // fPythiaP->event.list();
      // fPythiaP->stat();
      fEvent = fPythiaP->event;
    }
    else {
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

    for (int i = 5; i < np; i++) { // ignore first entry -> system + input particles

      int pdgc = fEvent[i].id();
      if (!PDGLibrary::Instance()->Find(pdgc)) continue; // some intermediate particles not part of genie tables

      int ks = fEvent[i].status();

      LongLorentzVector p4longo(fEvent[i].px(), fEvent[i].py(), fEvent[i].pz(), fEvent[i].e());
      p4longo.BoostZ(beta); 

      TLorentzVector p4o( (double)p4longo.Px(), (double)p4longo.Py(), (double)p4longo.Pz(), (double)p4longo.E() );
      p4o.RotateUz(unit_nu); 
     
      double massPDG = PDGLibrary::Instance()->Find(pdgc)->Mass();
      if ( ks>0 && p4o.E()<massPDG ) {
        LOG("PhotonRESGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
        LOG("PhotonRESGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
        LOG("PhotonRESGenerator", pWARN) << "E = " << p4o.E() << " // |p| = " << TMath::Sqrt(p4o.P()); 
        LOG("PhotonRESGenerator", pWARN) << "p = [ " << p4o.Px() << " , "  << p4o.Py() << " , "  << p4o.Pz() << " ]";
        LOG("PhotonRESGenerator", pWARN) << "m    = " << p4o.M() << " // mpdg = " << massPDG;
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

      event->AddParticle(pdgc, ist, firstmother, lastmother, firstchild, lastchild, p4o, pos );
    }

  }

  return true;
#else
  return false;
#endif

}
//___________________________________________________________________________
void PhotonRESWdecPythia8::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonRESWdecPythia8::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonRESWdecPythia8::LoadConfig(void)
{

#ifdef __GENIE_PYTHIA8_ENABLED__
  GetParam( "SSBarSuppression",       fSSBarSuppression       );
  GetParam( "GaussianPt2",            fGaussianPt2            );
  GetParam( "NonGaussianPt2Tail",     fNonGaussianPt2Tail     );
  GetParam( "RemainingEnergyCutoff",  fRemainingECutoff       );
  GetParam( "DiQuarkSuppression",     fDiQuarkSuppression     );
  GetParam( "LightVMesonSuppression", fLightVMesonSuppression );
  GetParam( "SVMesonSuppression",     fSVMesonSuppression     );
  GetParam( "Lunda",                  fLunda                  );
  GetParam( "Lundb",                  fLundb                  );
  GetParam( "LundaDiq",               fLundaDiq               );

  GetParam( "Q2Grid-Min", fQ2PDFmin );

  fPythiaP->settings.parm("StringFlav:probStoUD",         fSSBarSuppression);
  fPythiaP->settings.parm("Diffraction:primKTwidth",      fGaussianPt2);
  fPythiaP->settings.parm("StringPT:enhancedFraction",    fNonGaussianPt2Tail);
  fPythiaP->settings.parm("StringFragmentation:stopMass", fRemainingECutoff);
  fPythiaP->settings.parm("StringFlav:probQQtoQ",         fDiQuarkSuppression);
  fPythiaP->settings.parm("StringFlav:mesonUDvector",     fLightVMesonSuppression);
  fPythiaP->settings.parm("StringFlav:mesonSvector",      fSVMesonSuppression);
  fPythiaP->settings.parm("StringZ:aLund",                fLunda);
  fPythiaP->settings.parm("StringZ:bLund",                fLundb);
  fPythiaP->settings.parm("StringZ:aExtraDiquark",        fLundaDiq);

  fPythiaN->settings.parm("StringFlav:probStoUD",         fSSBarSuppression);
  fPythiaN->settings.parm("Diffraction:primKTwidth",      fGaussianPt2);
  fPythiaN->settings.parm("StringPT:enhancedFraction",    fNonGaussianPt2Tail);
  fPythiaN->settings.parm("StringFragmentation:stopMass", fRemainingECutoff);
  fPythiaN->settings.parm("StringFlav:probQQtoQ",         fDiQuarkSuppression);
  fPythiaN->settings.parm("StringFlav:mesonUDvector",     fLightVMesonSuppression);
  fPythiaN->settings.parm("StringFlav:mesonSvector",      fSVMesonSuppression);
  fPythiaN->settings.parm("StringZ:aLund",                fLunda);
  fPythiaN->settings.parm("StringZ:bLund",                fLundb);
  fPythiaN->settings.parm("StringZ:aExtraDiquark",        fLundaDiq);

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
      LOG("PhotonRESGenerator", pDEBUG) << "Configured to inhibit decays for  " <<  pdg_code;
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

}
//____________________________________________________________________________
void PhotonRESWdecPythia8::Initialize(void) const
{

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
