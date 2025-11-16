//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC (Valencia)
*/
//____________________________________________________________________________

#include "Physics/HELepton/EventGen/PhotonCOHWdecPythia8.h"

//___________________________________________________________________________
PhotonCOHWdecPythia8::PhotonCOHWdecPythia8() :
EventRecordVisitorI("genie::PhotonCOHWdecPythia8")
{
  this->Initialize();
}
//___________________________________________________________________________
PhotonCOHWdecPythia8::PhotonCOHWdecPythia8(string config) :
EventRecordVisitorI("genie::PhotonCOHWdecPythia8", config)
{
  this->Initialize();
}
//___________________________________________________________________________
PhotonCOHWdecPythia8::~PhotonCOHWdecPythia8()
{

}
//____________________________________________________________________________
void PhotonCOHWdecPythia8::ProcessEventRecord(GHepRecord * event) const
{

  if(!this->Wdecay(event)) {
    LOG("PhotonCOHGenerator", pWARN) << "W decayer failed!";
    event->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Could not simulate the W decay system");
    exception.SwitchOnFastForward();
    throw exception;
    return;
  }

}
//___________________________________________________________________________
bool PhotonCOHWdecPythia8::Wdecay(GHepRecord * 
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

  GHepParticle * target = event -> TargetNucleus();
  if(target) event->AddParticle(target->Pdg(), kIStFinalStateNuclearRemnant, 1,-1,-1,-1, *(target->P4()), *(target->X4()) );

  TVector3 unit_nu = nu->P4()->Vect().Unit();

  long double Ev = init_state.ProbeE(kRfLab);

  long double Mtgt = init_state.Tgt().Z()*kProtonMass + init_state.Tgt().N()*kNeutronMass;

  long double n1 = interaction->Kine().GetKV(kKVn1);
  long double n2 = interaction->Kine().GetKV(kKVn2);
  long double n3 = interaction->Kine().GetKV(kKVn3);

  long double costh = n1;
  long double sinth = sqrtl(1-costh*costh);

  long double mlout = 0;
  if      (pdg::IsNuE  (TMath::Abs(nu->Pdg()))) mlout = kElectronMass;
  else if (pdg::IsNuMu (TMath::Abs(nu->Pdg()))) mlout = kMuonMass;
  else if (pdg::IsNuTau(TMath::Abs(nu->Pdg()))) mlout = kTauMass;
  long double mlout2 = mlout*mlout;

  long double mL      = mlout+kMw;
  long double Delta   = sqrtl( powl(2.*Ev*Mtgt-mL*mL,2)-4.*powl(Mtgt*mL,2) );
  long double s12_min = Ev/(2.*Ev+Mtgt)*(mL*mL+2.*Ev*Mtgt-Delta);
  long double s12_max = Ev/(2.*Ev+Mtgt)*(mL*mL+2.*Ev*Mtgt+Delta);
  long double s12     = expl( logl(s12_min)+(logl(s12_max)-logl(s12_min))*n2);
  long double Q2_min  = powl(s12,2)*Mtgt/(2.*Ev*(2.*Ev*Mtgt-s12));
  long double Q2_max  = s12 - mL*mL;
  double Q2  = expl( logl(Q2_min) + (logl(Q2_max)-logl(Q2_min))*n3 );
  double s_r = s12 - Q2;

  long double EW = (s_r+kMw2-mlout2)/sqrtl(s_r)/2.;
  long double El = (s_r-kMw2+mlout2)/sqrtl(s_r)/2.;
  long double p  = sqrtl( EW*EW - kMw2 );
  LongLorentzVector p4_lout( 0., -p*sinth, -p*costh, El );

  long double bz = 4.*Ev*Mtgt*Q2/(Q2+s_r)/(2.*Ev*Mtgt-Q2) - (2.*Ev*Mtgt+Q2)/(2.*Ev*Mtgt-Q2);
  long double by = sqrtl(Mtgt*powl(Q2+s_r,2)/(2.*Ev*Q2*(s_r+Q2-2.*Ev*Mtgt))+1.);

  p4_lout.BoostZ(-bz);
  p4_lout.BoostY(-by);

  TLorentzVector p4l_o(p4_lout.Px(),p4_lout.Py(),p4_lout.Pz(),p4_lout.E());
  p4l_o.RotateX((double)acosl(by)-kPi/2.);


  // Randomize transverse components
  RandomGen * rnd = RandomGen::Instance();
  double phi  = 2 * kPi * rnd->RndLep().Rndm();

  p4l_o.RotateZ(phi);

  //rotate from LAB=[0,0,Ev,Ev]->[px,py,pz,E]
  p4l_o.RotateUz(unit_nu);

  int pdglout = 0;
  if      (pdg::IsAntiNuE  (nu->Pdg())) pdglout = kPdgPositron;
  else if (pdg::IsNuE      (nu->Pdg())) pdglout = kPdgElectron;
  else if (pdg::IsAntiNuMu (nu->Pdg())) pdglout = kPdgAntiMuon;
  else if (pdg::IsNuMu     (nu->Pdg())) pdglout = kPdgMuon;
  else if (pdg::IsAntiNuTau(nu->Pdg())) pdglout = kPdgAntiTau;
  else if (pdg::IsNuTau    (nu->Pdg())) pdglout = kPdgTau;

  // Create a GHepParticle and add it to the event record
  event->AddParticle( pdglout,  kIStStableFinalState, 0, -1, -1, -1, p4l_o, *(nu->X4()) );

  int pdgboson = pdg::IsNeutrino(init_state.ProbePdg()) ? kPdgWP : kPdgWM;

  fPythia->event.reset();
  fPythia->event.append(pdgboson, 23, 0, 0, 0., p*sinth, p*costh, EW, kMw);

  fPythia->next();

  // fPythia->event.list();
  // fPythia->stat();

  Pythia8::Event &fEvent = fPythia->event;
  int np = fEvent.size();
  assert(np>0);

  for (int i = 1; i < np; i++) {

    int pdgc = fEvent[i].id();
    if (!PDGLibrary::Instance()->Find(pdgc)) continue; // some intermediate particles not part of genie tables

    int ks = fEvent[i].status();

    LongLorentzVector p4longo(fEvent[i].px(), fEvent[i].py(), fEvent[i].pz(), fEvent[i].e());
    p4longo.BoostZ(-bz);
    p4longo.BoostY(-by);

    TLorentzVector p4o(p4longo.Px(),p4longo.Py(),p4longo.Pz(),p4longo.E());
    p4o.RotateX((double)acosl(by)-kPi/2.);
    p4o.RotateZ(phi);
    p4o.RotateUz(unit_nu);

    double massPDG = PDGLibrary::Instance()->Find(pdgc)->Mass();
    if ( ks>0 && p4o.E()<massPDG ) {
      LOG("PhotonCOHGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
      LOG("PhotonCOHGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
      LOG("PhotonCOHGenerator", pWARN) << "E = " << p4o.E() << " // |p| = " << TMath::Sqrt(p4o.P());
      LOG("PhotonCOHGenerator", pWARN) << "p = [ " << p4o.Px() << " , "  << p4o.Py() << " , "  << p4o.Pz() << " ]";
      LOG("PhotonCOHGenerator", pWARN) << "m    = " << p4o.M() << " // mpdg = " << massPDG;
      p4o.SetXYZT(0,0,0,massPDG);
    }

    // copy final state particles to the event record
    GHepStatus_t ist = (ks>0) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

    // fix numbering scheme used for mother/daughter assignments
    int firstmother = -1;
    int lastmother  = -1;
    int firstchild  = -1;
    int lastchild   = -1;

    if (fEvent[i].mother1()==0) {
      firstmother = 0;
    }
    else {
      firstmother = fEvent[i].mother1() + 3;
      if (fEvent[i].daughter1()!=0) firstchild  = fEvent[i].daughter1() + 3;
      if (fEvent[i].daughter2()!=0) lastchild   = fEvent[i].daughter2() + 3;

    }

    double vx = nu->X4()->X() + fEvent[i].xProd()*1e12; //pythia gives position in [mm] while genie uses [fm]
    double vy = nu->X4()->Y() + fEvent[i].yProd()*1e12;
    double vz = nu->X4()->Z() + fEvent[i].zProd()*1e12;
    double vt = nu->X4()->T() + fEvent[i].tProd()*(units::millimeter/units::second);
    TLorentzVector pos( vx, vy, vz, vt );

    event->AddParticle(pdgc, ist, firstmother, lastmother, firstchild, lastchild, p4o, pos );

  }

  return true;
#else
  return false;
#endif

}
//___________________________________________________________________________
void PhotonCOHWdecPythia8::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonCOHWdecPythia8::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonCOHWdecPythia8::LoadConfig(void)
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

  LOG("PhotonCOHGenerator", pINFO) << "Initialising PYTHIA..." ;
  fPythia->init(); 
#endif

}
//____________________________________________________________________________
void PhotonCOHWdecPythia8::Initialize(void) const
{

#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia = Pythia8Singleton::Instance()->Pythia8();

  fPythia->readString("Print:quiet = on");

  // sync GENIE and PYTHIA8 seeds
  RandomGen * rnd = RandomGen::Instance();
  long int seed = rnd->GetSeed();
  fPythia->readString("Random:setSeed = on");
  fPythia->settings.mode("Random:seed", seed);
  LOG("LeptoHad", pINFO) << "PYTHIA8  seed = " << fPythia->settings.mode("Random:seed");

  //needed to only do hadronization
  fPythia->readString("ProcessLevel:all = off");
#endif

}
