//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC (Valencia)
*/
//____________________________________________________________________________

#include "Physics/HELepton/EventGen/GLRESWdecPythia8.h"

//___________________________________________________________________________
GLRESWdecPythia8::GLRESWdecPythia8() :
EventRecordVisitorI("genie::GLRESWdecPythia8")
{
  this->Initialize();
  born = new Born();
}
//___________________________________________________________________________
GLRESWdecPythia8::GLRESWdecPythia8(string config) :
EventRecordVisitorI("genie::GLRESWdecPythia8", config)
{
  this->Initialize();
  born = new Born();
}
//___________________________________________________________________________
GLRESWdecPythia8::~GLRESWdecPythia8()
{

}
//____________________________________________________________________________
void GLRESWdecPythia8::ProcessEventRecord(GHepRecord * event) const
{

  if(!this->Wdecay(event)) {
    LOG("GLRESGenerator", pWARN) << "W decayer failed!";
    event->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Could not simulate the W decay system");
    exception.SwitchOnFastForward();
    throw exception;
    return;
  }

}
//___________________________________________________________________________
bool GLRESWdecPythia8::Wdecay(GHepRecord * 
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
  GHepParticle * el = event->HitElectron();

  GHepParticle * target = event -> TargetNucleus();
  if(target) event->AddParticle(target->Pdg(), kIStFinalStateNuclearRemnant, 1,-1,-1,-1, *(target->P4()), *(target->X4()) );

  TVector3 unit_nu = nu->P4()->Vect().Unit();

  long double mlout = interaction->FSPrimLepton()->Mass(); //mass of charged lepton
  long double mlin  = kElectronMass;                       //mass of incoming charged lepton

  long double Enuin = init_state.ProbeE(kRfLab);
  long double s = born->GetS(mlin,Enuin);

  long double n1 = interaction->Kine().GetKV(kKVn1);
  long double n2 = interaction->Kine().GetKV(kKVn2);

  long double costhCM = n1;
  long double sinthCM = sqrtl(1-costhCM*costhCM);
  
  long double t = born->GetT(mlin,mlout,s,n1);
  long double zeta  = born->GetReAlpha()/kPi*(2.0*logl(sqrtl(-t)/kElectronMass)-1.0);
  long double omx   = powl(n2, 1.0/zeta );
  long double s_r = s*( 1.-omx );

  // Boost velocity CM -> LAB
  long double EnuinCM = (s_r-mlin*mlin)/sqrtl(s_r)/2.;
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

    int pdgboson = pdg::IsNeutrino(init_state.ProbePdg()) ? kPdgWP : kPdgWM;

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
    fPythia->settings.mode("Random:seed", rnd->RndLep().Integer(100000000));
    fPythia->settings.parm("Beams:eCM",sqrtl(s_r));
    fPythia->init(); 
    fPythia->next();

    // fPythia->event.list();
    // fPythia->stat();

    Pythia8::Event &fEvent = fPythia->event;
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
          firstchild  = fEvent[i].daughter1() - 1;
          lastchild   = fEvent[i].daughter2() - 1;        
      }
      else { //rest
        firstmother = fEvent[i].mother1() - 1; //shift to match boson position
        firstchild  = fEvent[i].daughter1() - 1;
        lastchild   = fEvent[i].daughter2() - 1;
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
void GLRESWdecPythia8::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESWdecPythia8::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESWdecPythia8::LoadConfig(void)
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

  fPythia->settings.parm("StringFlav:probStoUD",         fSSBarSuppression);
  fPythia->settings.parm("Diffraction:primKTwidth",      fGaussianPt2);
  fPythia->settings.parm("StringPT:enhancedFraction",    fNonGaussianPt2Tail);
  fPythia->settings.parm("StringFragmentation:stopMass", fRemainingECutoff);
  fPythia->settings.parm("StringFlav:probQQtoQ",         fDiQuarkSuppression);
  fPythia->settings.parm("StringFlav:mesonUDvector",     fLightVMesonSuppression);
  fPythia->settings.parm("StringFlav:mesonSvector",      fSVMesonSuppression);
  fPythia->settings.parm("StringZ:aLund",                fLunda);
  fPythia->settings.parm("StringZ:bLund",                fLundb);
  fPythia->settings.parm("StringZ:aExtraDiquark",        fLundaDiq);

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
      auto pdentry = fPythia->particleData.particleDataEntryPtr(pdg_code);
      for (int ichan=0; ichan<=pdentry->sizeChannels()-1; ++ichan) {
        int onMode = pdentry->channel(ichan).onMode();
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
        pdentry->channel(ichan).onMode(onMode);
      }
    }
  }

  //Important: We dont initialize Pythia8 here because it must be done event by event. See above
#endif

}
//____________________________________________________________________________
void GLRESWdecPythia8::Initialize(void) const
{

#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia = new Pythia8::Pythia();

  fPythia->readString("Print:quiet = on");
  fPythia->readString("Random:setSeed = on");

  //One cool feature of having independent Pythia8 instances
  //We can define our intial state with the proper setting (i.e. disabling decays)
  //without affecting other classes of GENIE that also use Pythia8
  fPythia->readString("WeakSingleBoson:ffbar2ffbar(s:W) = on");
  fPythia->readString("PDF:lepton = off");
  fPythia->readString("24:onMode = off");
  fPythia->readString("24:onIfAny = 1 2 3 4 5"); //enable W->hadron only 
  fPythia->readString("Beams:idA = -12");
  fPythia->readString("Beams:idB = 11");
#endif

}
