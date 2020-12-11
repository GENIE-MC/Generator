//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
*/
//____________________________________________________________________________

#include <cstring>

#include <RVersion.h>
#include <TClonesArray.h>
#include <TMath.h>

#include "Framework/Conventions/Constants.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h"
#include "Framework/EventGen/EVGThreadException.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Physics/GlashowResonance/EventGen/GLRESGenerator.h"

#ifdef __GENIE_PYTHIA6_ENABLED__
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#endif // __GENIE_PYTHIA6_ENABLED__

using namespace genie;
using namespace genie::constants;
using namespace genie::utils::math;

//___________________________________________________________________________
GLRESGenerator::GLRESGenerator() :
EventRecordVisitorI("genie::GLRESGenerator")
{

}
//___________________________________________________________________________
GLRESGenerator::GLRESGenerator(string config) :
EventRecordVisitorI("genie::GLRESGenerator", config)
{

}
//___________________________________________________________________________
GLRESGenerator::~GLRESGenerator()
{

}
//___________________________________________________________________________
void GLRESGenerator::ProcessEventRecord(GHepRecord *
  #ifdef __GENIE_PYTHIA6_ENABLED__
    event // avoid unused variable warning if PYTHIA6 is not enabled
  #endif
) const
{

#ifdef __GENIE_PYTHIA6_ENABLED__

  Interaction * interaction = event->Summary();
  const InitialState & init_state = interaction->InitState();

  //incoming v & struck particle & remnant nucleus
  GHepParticle * nu     = event -> Probe();
  GHepParticle * el     = event -> HitElectron();
  GHepParticle * target = event -> TargetNucleus();
  assert(nu);
  assert(el);
  assert(target);

  if(target) event->AddParticle(target->Pdg(), kIStStableFinalState, 1,-1,-1,-1, *(target->P4()), *(target->X4()) );

  LongLorentzVector p4_nu( * event->Probe()->P4()       );
  LongLorentzVector p4_el( * event->HitElectron()->P4() );
  LongLorentzVector p4_Wlong( p4_nu.Px()+p4_el.Px(), p4_nu.Py()+p4_el.Py(), p4_nu.Pz()+p4_el.Pz(), p4_nu.E()+p4_el.E() );

  long double Wmass = p4_Wlong.M();
  LOG("GLRESGenerator", pINFO) << "Wmass : " << Wmass;

  if(Wmass < fWmin) {
    LOG("GLRESGenerator", pWARN) << "Low invariant mass, W = " << Wmass << " GeV!!";
    event->EventFlags()->SetBitNumber(kHadroSysGenErr, true);
    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Could not simulate the hadronic system");
    exception.SwitchOnFastForward();
    throw exception;
    return;
  }

  // Final state primary lepton PDG code
  int pdgl = interaction->FSPrimLeptonPdg();
  assert(pdgl!=0);

  LOG("GLRESGenerator", pINFO) << "Channel : " << interaction->FSPrimLeptonPdg();

  if ( pdg::IsElectron(pdgl) || pdg::IsMuon(pdgl) || pdg::IsTau(pdgl) ) {

    // Get selected kinematics
    long double y = interaction->Kine().y(true);
    assert(y>0 && y<1);

    // Compute the neutrino and muon energy
    long double Ev  = init_state.ProbeE(kRfLab);
    long double El  = y*Ev;

    LOG("GLRESGenerator", pINFO) << "Ev = " << Ev << ", y = " << y << ", -> El = " << El;

    // Compute the momentum transfer and scattering angle
    long double El2   = powl(El,2);
    long double ml    = interaction->FSPrimLepton()->Mass();
    long double ml2   = powl(ml,2);
    long double pl    = sqrtl(El2-ml2);

    long double Q2    = 2*(Ev-El)*kElectronMass + kElectronMass*kElectronMass;
    long double costh = (El-0.5*(Q2+ml2)/Ev)/pl;
    long double sinth = sqrtl( 1-powl(costh,2.) );
    // Randomize transverse components
    RandomGen * rnd = RandomGen::Instance();
    long double phi  = 2* constants::kPi * rnd->RndLep().Rndm();

    LOG("GLRESGenerator", pINFO) << "Q2 = " << Q2 << ", cos(theta) = " << costh << ", phi = " << phi;

    long double plx = pl*sinth*cosl(phi);
    long double ply = pl*sinth*sinl(phi);
    long double plz = pl*costh;

    // Lepton 4-momentum in the LAB frame
    LongLorentzVector p4lplong( plx, ply, plz, El );
    p4lplong.Rotate(p4_nu);
    LongLorentzVector p4nulong( p4_Wlong.Px()-p4lplong.Px(), p4_Wlong.Py()-p4lplong.Py(), p4_Wlong.Pz()-p4lplong.Pz(), p4_Wlong.E()-p4lplong.E() );

    TLorentzVector p4_W    ( (double)p4_Wlong.Px(), (double)p4_Wlong.Py(), (double)p4_Wlong.Pz(), (double)p4_Wlong.E() );
    TLorentzVector p4_lpout( (double)p4lplong.Px(), (double)p4lplong.Py(), (double)p4lplong.Pz(), (double)p4lplong.E() );
    TLorentzVector p4_nuout( (double)p4nulong.Px(), (double)p4nulong.Py(), (double)p4nulong.Pz(), (double)p4nulong.E() );

    int pdgvout = 0;
    if      ( pdg::IsElectron(pdgl) ) pdgvout = kPdgAntiNuE;
    else if ( pdg::IsMuon(pdgl)     ) pdgvout = kPdgAntiNuMu;
    else if ( pdg::IsTau(pdgl)      ) pdgvout = kPdgAntiNuTau;

    // Create a GHepParticle and add it to the event record
    event->AddParticle( kPdgWM,  kIStDecayedState,     0, -1,  5,  6, p4_W,     *(nu->X4()) ); //W [mothers: nuebar_in,e_in][daugthers: nulbar_out,lp_out]
    event->AddParticle( pdgvout, kIStStableFinalState, 4, -1, -1, -1, p4_nuout, *(nu->X4()) );
    event->AddParticle( pdgl,    kIStStableFinalState, 4, -1, -1, -1, p4_lpout, *(nu->X4()) );
    event->Summary()->KinePtr()->SetFSLeptonP4(p4_lpout);

  }
  else {

    char p6frame[10], p6nu[10], p6tgt[10];
    strcpy(p6frame, "CMS"    );
    strcpy(p6nu,    "nu_ebar");
    strcpy(p6tgt,   "e-"     );


    int mstp61  = fPythia->GetMSTP(61);
    int mstp71  = fPythia->GetMSTP(71);
    int mdme206 = fPythia->GetMDME(206,1);
    int mdme207 = fPythia->GetMDME(207,1);
    int mdme208 = fPythia->GetMDME(208,1);
    double pmas1_W = fPythia->GetPMAS(24,1);
    double pmas2_W = fPythia->GetPMAS(24,2);
    double pmas2_t = fPythia->GetPMAS(6,2);
    fPythia->SetMSTP(61,0);    // (Default=2) master switch for initial-state QCD and QED radiation.
    fPythia->SetMSTP(71,0);    // (Default=2) master switch for initial-state QCD and QED radiation.
    fPythia->SetMDME(206,1,0); //swicht off W decay leptonic modes
    fPythia->SetMDME(207,1,0);
    fPythia->SetMDME(208,1,0);
    fPythia->SetPMAS(24,1,kMw); //mass of the W boson (pythia=80.450 // genie=80.385)
    fPythia->SetPMAS(24,2,0.);  //set to 0 the width of the W boson to avoid problems with energy conservation
    fPythia->SetPMAS(6,2,0.);  //set to 0 the width of the top to avoid problems with energy conservation

    fPythia->Pyinit(p6frame, p6nu, p6tgt, (double)Wmass);

    fPythia->Pyevnt();

    fPythia->SetMSTP(61,mstp61);
    fPythia->SetMSTP(71,mstp71);
    fPythia->SetMDME(206,1,mdme206);
    fPythia->SetMDME(207,1,mdme207);
    fPythia->SetMDME(208,1,mdme208);
    fPythia->SetPMAS(24,1,pmas1_W);
    fPythia->SetPMAS(24,2,pmas2_W);
    fPythia->SetPMAS(6,2,pmas2_t);


    // Use for debugging purposes
    //fPythia->Pylist(3);

    // get LUJETS record
    fPythia->GetPrimaries();
    TClonesArray * pythia_particles = (TClonesArray *) fPythia->ImportParticles("All");
    // copy PYTHIA container to a new TClonesArray so as to transfer ownership
    // of the container and of its elements to the calling method
    int np = pythia_particles->GetEntries();
    assert(np>0);
    TClonesArray * particle_list = new TClonesArray("genie::GHepParticle", np);
    particle_list->SetOwner(true);

    // Boost velocity LAB' -> Resonance rest frame
    long double beta = p4_Wlong.P()/p4_Wlong.E();

    TMCParticle * p = 0;
    TIter particle_iter(pythia_particles);
    while( (p = (TMCParticle *) particle_iter.Next()) ) {

      int pdgc   = p->GetKF();
      int ks     = p->GetKS();
      int parent = p->GetParent();

      if ( ks==21 ) { continue; } //we dont want to save first particles from pythia (init states)

      LongLorentzVector p4long( p->GetPx(), p->GetPy(), p->GetPz(), p->GetEnergy()  );
      p4long.Boost(beta);
      p4long.Rotate(p4_Wlong);

      TLorentzVector p4( (double)p4long.Px(), (double)p4long.Py(), (double)p4long.Pz(), (double)p4long.E() );

      double massPDG = PDGLibrary::Instance()->Find(pdgc)->Mass();
      if ( (ks==1 || ks==4) && p4.E() < massPDG ) {
        LOG("GLRESGenerator", pWARN) << "Putting at rest one stable particle generated by PYTHIA because E < m";
        LOG("GLRESGenerator", pWARN) << "PDG = " << pdgc << " // State = " << ks;
        LOG("GLRESGenerator", pWARN) << "E = " << p4.E() << " // |p| = " << TMath::Sqrt(p4.P());
        LOG("GLRESGenerator", pWARN) << "p = [ " << p4.Px() << " , "  << p4.Py() << " , "  << p4.Pz() << " ]";
        LOG("GLRESGenerator", pWARN) << "m    = " << p4.M() << " // mpdg = " << massPDG;
        p4.SetXYZT(0,0,0,massPDG);
      }

      // copy final state particles to the event record
      GHepStatus_t ist = (ks==1 || ks==4) ? kIStStableFinalState : kIStDISPreFragmHadronicState;

      // fix numbering scheme used for mother/daughter assignments
      int firstmother = -1;
      int lastmother  = -1;
      int firstchild  = -1;
      int lastchild   = -1;

      if ( parent < 10 ) {  // I=10 is the position were W boson is stored
        if      ( TMath::Abs(pdgc)<7 ) {   //outgoing quarks: mother will be the boson (saved in position 4)
          firstmother = 4;
          firstchild  = p->GetFirstChild() - 6;  //boson is stored in I=10 with pythia and I=4 for GENIE => shift 6
          lastchild   = p->GetLastChild()  - 6;
        }
        else if ( TMath::Abs(pdgc)==24 ) { //produced W boson: mother will be the incoming neutrino
          firstmother = 0;
          firstchild  = p->GetFirstChild() - 6;
          lastchild   = p->GetLastChild()  - 6;
        }
        else if ( pdgc==22 ) {             //radiative photons: mother will be the incoming electron
          firstmother = 2;
        }
      }
      else { //rest
        firstmother = parent - 6; //shift to match boson position
        firstchild  = (p->GetFirstChild()==0) ? p->GetFirstChild() - 1 : p->GetFirstChild() - 6;
        lastchild   = (p->GetLastChild()==0)  ? p->GetLastChild()  - 1 : p->GetLastChild()  - 6;
      }

      double vx = nu->X4()->X() + p->GetVx()*1e12; //pythia gives position in [mm] while genie uses [fm]
      double vy = nu->X4()->Y() + p->GetVy()*1e12;
      double vz = nu->X4()->Z() + p->GetVz()*1e12;
      double vt = nu->X4()->T() + p->GetTime()*(units::millimeter/units::second);
      TLorentzVector pos( vx, vy, vz, vt );

      event->AddParticle(pdgc, ist, firstmother, lastmother, firstchild, lastchild, p4, pos );

    }

  }

#else
  LOG("GLRESGenerator", pFATAL)
    << "Calling GENIE/PYTHIA6 without enabling PYTHIA6";
  gAbortingInErr = true;
  std::exit(1);
#endif

}
//___________________________________________________________________________
void GLRESGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void GLRESGenerator::LoadConfig(void)
{
#ifdef __GENIE_PYTHIA6_ENABLED__
  fPythia = TPythia6::Instance();

  // sync GENIE/PYTHIA6 seed number
  RandomGen::Instance();

  // PYTHIA parameters only valid for HEDIS
  GetParam( "Xsec-Wmin", fWmin ) ;
  fPythia->SetPARP(2,  fWmin); //(D = 10. GeV) lowest c.m. energy for the event as a whole that the program will accept to simulate. (bellow 2GeV pythia crashes)

  int warnings;       GetParam( "Warnings",      warnings ) ;
  int errors;         GetParam( "Errors",        errors ) ;
  int qrk_mass;       GetParam( "QuarkMass",     qrk_mass ) ;
  fPythia->SetMSTU(26, warnings);     // (Default=10) maximum number of warnings that are printed
  fPythia->SetMSTU(22, errors);       // (Default=10) maximum number of errors that are printed
  fPythia->SetMSTJ(93, qrk_mass);     // light (d, u, s, c, b) quark masses are taken from PARF(101) - PARF(105) rather than PMAS(1,1) - PMAS(5,1). Diquark masses are given as sum of quark masses, without spin splitting term.

#endif // __GENIE_PYTHIA6_ENABLED__

}
//____________________________________________________________________________
