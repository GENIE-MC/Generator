//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC (Valencia)
*/
//____________________________________________________________________________

#include "Physics/HELepton/EventGen/GLRESWdecPythia6.h"

#ifdef __GENIE_PYTHIA6_ENABLED__
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#endif // __GENIE_PYTHIA6_ENABLED__

//___________________________________________________________________________
GLRESWdecPythia6::GLRESWdecPythia6() :
EventRecordVisitorI("genie::GLRESWdecPythia6")
{
  this->Initialize();
  born = new Born();
}
//___________________________________________________________________________
GLRESWdecPythia6::GLRESWdecPythia6(string config) :
EventRecordVisitorI("genie::GLRESWdecPythia6", config)
{
  this->Initialize();
  born = new Born();
}
//___________________________________________________________________________
GLRESWdecPythia6::~GLRESWdecPythia6()
{

}
//____________________________________________________________________________
void GLRESWdecPythia6::ProcessEventRecord(GHepRecord * event) const
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
bool GLRESWdecPythia6::Wdecay(GHepRecord * 
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

    char p6frame[10],p6nu[10],p6tgt[10];
    strcpy(p6frame, "CMS"     );
    strcpy(p6nu,    "nu_ebar" );   
    strcpy(p6tgt,   "e-"      );
    
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
     
      double massPDG = PDGLibrary::Instance()->Find(pdgc)->Mass();
      if ( (ks==1 || ks==4) && p4o.E()<massPDG ) {
        LOG("GLRESGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
        LOG("GLRESGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
        LOG("GLRESGenerator", pWARN) << "E = " << p4o.E() << " // |p| = " << TMath::Sqrt(p4o.P()); 
        LOG("GLRESGenerator", pWARN) << "p = [ " << p4o.Px() << " , "  << p4o.Py() << " , "  << p4o.Pz() << " ]";
        LOG("GLRESGenerator", pWARN) << "m    = " << p4o.M() << " // mpdg = " << massPDG;
        p4o.SetXYZT(0,0,0,massPDG);
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

      double vx = nu->X4()->X() + particle->GetVx()*1e12; //pythia gives position in [mm] while genie uses [fm]
      double vy = nu->X4()->Y() + particle->GetVy()*1e12;
      double vz = nu->X4()->Z() + particle->GetVz()*1e12;
      double vt = nu->X4()->T() + particle->GetTime()*(units::millimeter/units::second);
      TLorentzVector pos( vx, vy, vz, vt );

      event->AddParticle(pdgc, ist, firstmother, lastmother, firstchild, lastchild, p4o, pos );

    }

    delete particle;
    pythia_particles->Clear("C");

  }

  return true;
#else
  return false;
#endif


}
//___________________________________________________________________________
void GLRESWdecPythia6::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESWdecPythia6::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESWdecPythia6::LoadConfig(void)
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
#endif // __GENIE_PYTHIA6_ENABLED__

}
//____________________________________________________________________________
void GLRESWdecPythia6::Initialize(void) const
{

#ifdef __GENIE_PYTHIA6_ENABLED__
  fPythia = TPythia6::Instance();

  // sync GENIE/PYTHIA6 seed number
  RandomGen::Instance();
#endif

}
