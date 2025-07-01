//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC (Valencia)
*/
//____________________________________________________________________________

#include "Physics/HELepton/EventGen/PhotonCOHWdecPythia6.h"

#ifdef __GENIE_PYTHIA6_ENABLED__
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#endif

//___________________________________________________________________________
PhotonCOHWdecPythia6::PhotonCOHWdecPythia6() :
EventRecordVisitorI("genie::PhotonCOHWdecPythia6")
{
  this->Initialize();
}
//___________________________________________________________________________
PhotonCOHWdecPythia6::PhotonCOHWdecPythia6(string config) :
EventRecordVisitorI("genie::PhotonCOHWdecPythia6", config)
{
  this->Initialize();
}
//___________________________________________________________________________
PhotonCOHWdecPythia6::~PhotonCOHWdecPythia6()
{

}
//____________________________________________________________________________
void PhotonCOHWdecPythia6::ProcessEventRecord(GHepRecord * event) const
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
bool PhotonCOHWdecPythia6::Wdecay(GHepRecord * 
#ifdef __GENIE_PYTHIA6_ENABLED__
  event // avoid unused variable warning if PYTHIA6 is not enabled
#endif
) const
{

#ifdef __GENIE_PYTHIA6_ENABLED__
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

  int def61  = fPythia->GetMSTP(61);
  int def71  = fPythia->GetMSTP(71);
  fPythia->SetMSTP(61,0);    // (Default=2) master switch for initial-state QCD and QED radiation.
  fPythia->SetMSTP(71,0);    // (Default=2) master switch for initial-state QCD and QED radiation.

  fPythia->Py1ent( -1, pdgboson, EW, acosl(costh), kPi/2. ); //k(1,2) = 2
  fPythia->Pyexec();

  fPythia->SetMSTP(61,def61);
  fPythia->SetMSTP(71,def71);

  fPythia->GetPrimaries();
  TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");
  int np = pythia_particles->GetEntries();
  assert(np>0);

  TMCParticle * particle = 0;
  TIter piter(pythia_particles);
  while( (particle = (TMCParticle *) piter.Next()) ) {

    int pdgc = particle->GetKF();
    int ks = particle->GetKS();

    LongLorentzVector p4longo(particle->GetPx(), particle->GetPy(), particle->GetPz(), particle->GetEnergy());
    p4longo.BoostZ(-bz);
    p4longo.BoostY(-by);

    TLorentzVector p4o(p4longo.Px(),p4longo.Py(),p4longo.Pz(),p4longo.E());
    p4o.RotateX((double)acosl(by)-kPi/2.);
    p4o.RotateZ(phi);
    p4o.RotateUz(unit_nu);

    double massPDG = PDGLibrary::Instance()->Find(pdgc)->Mass();
    if ( (ks==1 || ks==4) && p4o.E()<massPDG ) {
      LOG("PhotonCOHGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
      LOG("PhotonCOHGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
      LOG("PhotonCOHGenerator", pWARN) << "E = " << p4o.E() << " // |p| = " << TMath::Sqrt(p4o.P());
      LOG("PhotonCOHGenerator", pWARN) << "p = [ " << p4o.Px() << " , "  << p4o.Py() << " , "  << p4o.Pz() << " ]";
      LOG("PhotonCOHGenerator", pWARN) << "m    = " << p4o.M() << " // mpdg = " << massPDG;
      p4o.SetXYZT(0,0,0,massPDG);
    }

    // copy final state particles to the event record
    GHepStatus_t ist = (ks==1 || ks==4) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

    // fix numbering scheme used for mother/daughter assignments
    int firstmother = -1;
    int lastmother  = -1;
    int firstchild  = -1;
    int lastchild   = -1;

    if (particle->GetParent()==0) {
      firstmother = 0;
    }
    else {
      firstmother = particle->GetParent() + 3;
      if (particle->GetFirstChild()!=0) firstchild  = particle->GetFirstChild() + 3;
      if (particle->GetLastChild() !=0) lastchild   = particle->GetLastChild() + 3;

    }

    double vx = nu->X4()->X() + particle->GetVx()*1e12; //pythia gives position in [mm] while genie uses [fm]
    double vy = nu->X4()->Y() + particle->GetVy()*1e12;
    double vz = nu->X4()->Z() + particle->GetVz()*1e12;
    double vt = nu->X4()->T() + particle->GetTime()*(units::millimeter/units::second);
    TLorentzVector pos( vx, vy, vz, vt );

    event->AddParticle(pdgc, ist, firstmother, lastmother, firstchild, lastchild, p4o, pos );

  }

  delete particle;
  pythia_particles->Clear("C");

  return true;
#else
  return false;
#endif

}
//___________________________________________________________________________
void PhotonCOHWdecPythia6::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonCOHWdecPythia6::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void PhotonCOHWdecPythia6::LoadConfig(void)
{

#ifdef __GENIE_PYTHIA6_ENABLED__
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

  // PYTHIA parameters only valid for HEDIS
  double wmin;        GetParam( "Xsec-Wmin",     wmin ) ;
  int warnings;       GetParam( "Warnings",      warnings ) ;
  int errors;         GetParam( "Errors",        errors ) ;
  int qrk_mass;       GetParam( "QuarkMass",     qrk_mass ) ;
  fPythia->SetPARP(2,  wmin);         //(D = 10. GeV) lowest c.m. energy for the event as a whole that the program will accept to simulate. (bellow 2GeV pythia crashes)
  fPythia->SetMSTU(26, warnings);     // (Default=10) maximum number of warnings that are printed
  fPythia->SetMSTU(22, errors);       // (Default=10) maximum number of errors that are printed
  fPythia->SetMSTJ(93, qrk_mass);     // light (d, u, s, c, b) quark masses are taken from PARF(101) - PARF(105) rather than PMAS(1,1) - PMAS(5,1). Diquark masses are given as sum of quark masses, without spin splitting term.
  fPythia->SetPMAS(24,1,kMw);  //mass of the W boson (pythia=80.450 // genie=80.385)
  fPythia->SetPMAS(24,2,0.);   //set to 0 the width of the W boson to avoid problems with energy conservation
  fPythia->SetPMAS(6,2,0.);    //set to 0 the width of the top to avoid problems with energy conservation
  fPythia->SetMDME(192,1,0);   // W->dbar+t decay off 
  fPythia->SetMDME(196,1,0);   // W->cbar+t decay off 
  fPythia->SetMDME(200,1,0);   // W->cbar+t decay off 

  // PYTHIA tuned parameters
  fPythia->SetPARJ(2,  fSSBarSuppression       );
  fPythia->SetPARJ(21, fGaussianPt2            );
  fPythia->SetPARJ(23, fNonGaussianPt2Tail     );
  fPythia->SetPARJ(33, fRemainingECutoff       );
  fPythia->SetPARJ(1,  fDiQuarkSuppression     );
  fPythia->SetPARJ(11, fLightVMesonSuppression );
  fPythia->SetPARJ(12, fSVMesonSuppression     );
  fPythia->SetPARJ(41, fLunda                  );
  fPythia->SetPARJ(42, fLundb                  );
  fPythia->SetPARJ(45, fLundaDiq               );
#endif

}
//____________________________________________________________________________
void PhotonCOHWdecPythia6::Initialize(void) const
{

#ifdef __GENIE_PYTHIA6_ENABLED__
  fPythia = TPythia6::Instance();

  // sync GENIE/PYTHIA6 seed number
  RandomGen::Instance();
#endif

}
