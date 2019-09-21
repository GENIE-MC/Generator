//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)

         Shivesh Mandalia <s.p.mandalia@qmul.ac.uk>
         Queen Mary University of London
*/
//____________________________________________________________________________

#include <RVersion.h>
#include <TVector3.h>
#include <TF1.h>
#include <TROOT.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/GHEP/GHepStatus.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepRecord.h"
#include "Framework/GHEP/GHepFlags.h" 
#include "Framework/EventGen/EVGThreadException.h"
#include "Physics/Hadronization/CharmP8Hadronization.h"
#include "Physics/Hadronization/FragmentationFunctionI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/KineUtils.h"
#include "Framework/Utils/StringUtils.h"
#include "Physics/Hadronization/FragmRecUtils.h"
#include "Framework/Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

//____________________________________________________________________________
CharmP8Hadronization::CharmP8Hadronization() :
EventRecordVisitorI("genie::CharmP8Hadronization")
{
  this->Initialize();
}
//____________________________________________________________________________
CharmP8Hadronization::CharmP8Hadronization(string config) :
EventRecordVisitorI("genie::CharmP8Hadronization", config)
{
  this->Initialize();
}
//____________________________________________________________________________
CharmP8Hadronization::~CharmP8Hadronization()
{
  delete fCharmPT2pdf;
  fCharmPT2pdf = 0;

  delete fD0FracSpl;
  fD0FracSpl = 0;

  delete fDpFracSpl;
  fDpFracSpl = 0;

  delete fDsFracSpl;
  fDsFracSpl = 0;
}
//____________________________________________________________________________
void CharmP8Hadronization::Initialize(void) const
{
#ifdef __GENIE_PYTHIA8_ENABLED__
  fPythia = Pythia8Singleton::Instance();

  fPythia->Pythia8()->readString("ProcessLevel:all = off");
  fPythia->Pythia8()->readString("Print:quiet      = on");
  fPythia->Pythia8()->init();
#endif
}
//____________________________________________________________________________
void CharmP8Hadronization::ProcessEventRecord(GHepRecord * event) const
{
  Interaction * interaction = event->Summary();
  TClonesArray * particle_list = this->Hadronize(interaction);

  if(! particle_list ) {
    LOG("CharmP8Hadronization", pWARN) << "Got an empty particle list. Hadronizer failed!";
    LOG("CharmP8Hadronization", pWARN) << "Quitting the current event generation thread";
    
    event->EventFlags()->SetBitNumber(kHadroSysGenErr, true);

    genie::exceptions::EVGThreadException exception;
    exception.SetReason("Could not simulate the hadronic system");
    exception.SwitchOnFastForward();
    throw exception;
    
    return;
  }
  
  int mom = event->FinalStateHadronicSystemPosition();
  assert(mom!=-1);

  // find the proper status for the particles we are going to put in event record
  bool is_nucleus = interaction->InitState().Tgt().IsNucleus();
  GHepStatus_t istfin = (is_nucleus) ?
    kIStHadronInTheNucleus : kIStStableFinalState ;

  // retrieve the hadronic blob lorentz boost
  // Because Hadronize() returned particles not in the LAB reference frame
  const TLorentzVector * had_syst = event -> Particle(mom) -> P4() ;
  TVector3 beta( 0., 0., had_syst -> P()/ had_syst -> Energy() ) ; 

  // Vector defining rotation from LAB to LAB' (z:= \vec{phad})
  TVector3 unitvq = had_syst -> Vect().Unit();
  
  GHepParticle * neutrino  = event->Probe();                                                                                                                                                 
  const TLorentzVector & vtx = *(neutrino->X4());                                                                                                                                            

  GHepParticle * particle = 0;
  TIter particle_iter(particle_list);
  while ((particle = (GHepParticle *) particle_iter.Next()))  {

    int pdgc = particle -> Pdg() ;

    //  bring the particle in the LAB reference frame
    particle -> P4() -> Boost( beta ) ;
    particle -> P4() -> RotateUz( unitvq ) ; 

    // set the proper status according to a number of things:
    // interaction on a nucleaus or nucleon, particle type
    GHepStatus_t ist = ( particle -> Status() ==1 ) ? istfin : kIStDISPreFragmHadronicState;

    // handle gammas, and leptons that might come from internal pythia decays
    // mark them as final state particles
    bool not_hadron = ( pdgc == kPdgGamma ||
			pdg::IsNeutralLepton(pdgc) ||
			pdg::IsChargedLepton(pdgc) ) ;

    if(not_hadron)  { ist = kIStStableFinalState; }
    particle -> SetStatus( ist ) ;

    int im  = mom + 1 + particle -> FirstMother() ;
    int ifc = ( particle -> FirstDaughter() == -1) ? -1 : mom + 1 + particle -> FirstDaughter();
    int ilc = ( particle -> LastDaughter()  == -1) ? -1 : mom + 1 + particle -> LastDaughter();

    particle -> SetFirstMother( im ) ;
    if ( ifc > -1 ) {
      particle -> SetFirstDaughter( ifc ) ; 
      particle -> SetLastDaughter( ilc ) ; 
    }
    
    // the Pythia particle position is overridden    
    particle -> SetPosition( vtx ) ;         
    
    event->AddParticle(*particle);
  }

  particle_list -> Delete() ;
  delete particle_list ; 
  
  // update the weight of the event
  event -> SetWeight ( Weight() * event->Weight() );

}
//____________________________________________________________________________
TClonesArray * CharmP8Hadronization::Hadronize(
                                        const Interaction * interaction) const
{
  LOG("CharmHad", pNOTICE) << "** Running CHARM hadronizer";

  PDGLibrary * pdglib = PDGLibrary::Instance();
  RandomGen *  rnd    = RandomGen::Instance();
  
  // ....................................................................
  // Get information on the input event
  //
  const InitialState & init_state = interaction -> InitState();
  const Kinematics &   kinematics = interaction -> Kine();
  const Target &       target     = init_state.Tgt();

  const TLorentzVector & p4Had = kinematics.HadSystP4();

  double Ev = init_state.ProbeE(kRfLab);
  double W  = kinematics.W(true);

  TVector3 beta = -1 * p4Had.BoostVector(); // boost vector for LAB' -> HCM'
  TLorentzVector p4H(0,0,0,W);              // hadronic system 4p @ HCM'

  double Eh = p4Had.Energy();

  LOG("CharmHad", pNOTICE) << "Ehad (LAB) = " << Eh << ", W = " << W;

  int  nu_pdg  = init_state.ProbePdg();
  int  nuc_pdg = target.HitNucPdg();
//int  qpdg    = target.HitQrkPdg();
//bool sea     = target.HitSeaQrk();
  bool isp     = pdg::IsProton (nuc_pdg);
  bool isn     = pdg::IsNeutron(nuc_pdg);
  bool isnu    = pdg::IsNeutrino(nu_pdg);
  bool isnub   = pdg::IsAntiNeutrino(nu_pdg);
  bool isdm    = pdg::IsDarkMatter(nu_pdg);

  // ....................................................................
  // Attempt to generate a charmed hadron & its 4-momentum
  //

  TLorentzVector p4C(0,0,0,0);
  int ch_pdg = -1;

  bool got_charmed_hadron = false;
  unsigned int itry=0;

  while(itry++ < kRjMaxIterations && !got_charmed_hadron) {

     // Generate a charmed hadron PDG code
     int    pdg = this->GenerateCharmHadron(nu_pdg,Ev); // generate hadron
     double mc  = pdglib->Find(pdg)->Mass();           // lookup mass
     
     LOG("CharmHad", pNOTICE) 
         << "Trying charm hadron = " << pdg << "(m = " << mc << ")";

     if(mc>=W) continue; // dont' accept

     // Generate the charmed hadron energy based on the input
     // fragmentation function
     double z   = fFragmFunc->GenerateZ();  // generate z(=Eh/Ev)
     double Ec  = z*Eh;                     // @ LAB'
     double mc2 = TMath::Power(mc,2);
     double Ec2 = TMath::Power(Ec,2);
     double pc2 = Ec2-mc2;

     LOG("CharmHad", pINFO) 
         << "Trying charm hadron z = " << z << ", E = " << Ec;
     
     if(pc2<=0) continue; 

     // Generate the charm hadron pT^2 and pL^2 (with respect to the
     // hadronic system direction @ the LAB)
     double ptc2 = fCharmPT2pdf->GetRandom();   
     double plc2 = Ec2 - ptc2 - mc2;
     LOG("CharmHad", pINFO) 
           << "Trying charm hadron pT^2 (tranv to pHad) = " << ptc2;
     if(plc2<0) continue;

     // Generate the charm hadron momentum components (@ LAB', z:\vec{pHad})
     double ptc = TMath::Sqrt(ptc2);
     double plc = TMath::Sqrt(plc2);   
     double phi = (2*kPi) * rnd->RndHadro().Rndm();
     double pxc = ptc * TMath::Cos(phi);
     double pyc = ptc * TMath::Sin(phi);
     double pzc = plc; 

     p4C.SetPxPyPzE(pxc,pyc,pzc,Ec); // @ LAB'

     // Boost charm hadron 4-momentum from the LAB' to the HCM' frame
     //
     LOG("CharmHad", pDEBUG) 
       << "Charm hadron p4 (@LAB') = " << utils::print::P4AsString(&p4C);

     p4C.Boost(beta);

     LOG("CharmHad", pDEBUG) 
        << "Charm hadron p4 (@HCM') = " << utils::print::P4AsString(&p4C);

     // Hadronic non-charm remnant 4p at HCM'
     TLorentzVector p4 = p4H - p4C;
     double wr = p4.M();
     LOG("CharmHad", pINFO) 
             << "Invariant mass of remnant hadronic system= " << wr;
     if(wr < kNucleonMass + kPionMass + kASmallNum) {
        LOG("CharmHad", pINFO) << "Too small hadronic remnant mass!";
        continue;
     }
    
     ch_pdg = pdg;
     got_charmed_hadron = true;

     LOG("CharmHad", pNOTICE) 
           << "Generated charm hadron = " << pdg << "(m = " <<  mc << ")";
     LOG("CharmHad", pNOTICE) 
           << "Generated charm hadron z = " << z << ", E = " << Ec;
  }

  // ....................................................................
  // Check whether the code above had difficulty generating the charmed
  // hadron near the W - threshold.
  // If yes, attempt a phase space decay of a low mass charm hadron + nucleon
  // pair that maintains the charge.
  // That's a desperate solution but don't want to quit too early as that 
  // would distort the generated dsigma/dW distribution near threshold.
  //
  bool used_lowW_strategy = false;
  int  fs_nucleon_pdg = -1;
  if(ch_pdg==-1 && W < 3.){
     LOG("CharmHad", pNOTICE) 
         << "Had difficulty generating charm hadronic system near W threshold";
     LOG("CharmHad", pNOTICE) 
         << "Trying an alternative strategy";

     double qfsl  = interaction->FSPrimLepton()->Charge() / 3.;
     double qinit = pdglib->Find(nuc_pdg)->Charge() / 3.;
     int qhad  = (int) (qinit - qfsl);

     int remn_pdg = -1; 
     int chrm_pdg = -1; 

     //cc-only: qhad(nu) = +1,+2, qhad(nubar)= -1,0
     //
     if(qhad == 2) { 
         chrm_pdg = kPdgDP; remn_pdg = kPdgProton; 
     } else if(qhad ==  1) { 
         if(rnd->RndHadro().Rndm() > 0.5) { 
            chrm_pdg = kPdgD0; remn_pdg = kPdgProton; 
         } else { 
            chrm_pdg = kPdgDP; remn_pdg = kPdgNeutron; 
         }
     } else if(qhad ==  0) { 
         chrm_pdg = kPdgAntiD0; remn_pdg = kPdgNeutron; 
     } else if(qhad == -1) { 
         chrm_pdg = kPdgDM; remn_pdg = kPdgNeutron; 
     } 

     double mc  = pdglib->Find(chrm_pdg)->Mass();           
     double mn  = pdglib->Find(remn_pdg)->Mass();          

     if(mc+mn < W) {
        // Set decay
        double mass[2] = {mc, mn};
        bool permitted = fPhaseSpaceGenerator.SetDecay(p4H, 2, mass);
        assert(permitted);
 
        // Get the maximum weight
        double wmax = -1;
        for(int i=0; i<200; i++) {
           double w = fPhaseSpaceGenerator.Generate();
           wmax = TMath::Max(wmax,w);
        }

        if(wmax>0) {
           wmax *= 2;

           // Generate unweighted decay
           bool accept_decay=false;
           unsigned int idecay_try=0;
           while(!accept_decay)
           {
             idecay_try++;

             if(idecay_try>kMaxUnweightDecayIterations) {
                LOG("CharmHad", pWARN)
                  << "Couldn't generate an unweighted phase space decay after "
                  << idecay_try << " attempts";
             }
             double w  = fPhaseSpaceGenerator.Generate();
             if(w > wmax) {
                LOG("CharmHad", pWARN)
                 << "Decay weight = " << w << " > max decay weight = " << wmax;
             }
             double gw = wmax * rnd->RndHadro().Rndm();
             accept_decay = (gw<=w);

             if(accept_decay) {
                 used_lowW_strategy = true;
                 TLorentzVector * p4 = fPhaseSpaceGenerator.GetDecay(0);
                 p4C            = *p4;
                 ch_pdg         = chrm_pdg;
                 fs_nucleon_pdg = remn_pdg;
             } 
           } // decay loop
        }//wmax>0

     }// allowed decay
  } // alt low-W strategy

  // ....................................................................
  // Check success in generating the charm hadron & compute 4p for
  // remnant system
  //
  if(ch_pdg==-1){
     LOG("CharmHad", pWARN) 
         << "Couldn't generate charm hadron for: " << *interaction;
     return 0;
  }

  TLorentzVector p4R = p4H - p4C;
  double WR = p4R.M();
  //double MC = pdglib->Find(ch_pdg)->Mass();

  LOG("CharmHad", pNOTICE) << "Remnant hadronic system mass = " << WR;

  // ....................................................................
  // Handle case where the user doesn't want to remnant system to be
  // hadronized (add as 'hadronic blob') 
  //
  if(fCharmOnly) {
    // Create particle list (fragmentation record)
    TClonesArray * particle_list = new TClonesArray("genie::GHepParticle", 2);
    particle_list->SetOwner(true);

    // insert the generated particles
    new ((*particle_list)[0]) GHepParticle (ch_pdg,kIStStableFinalState,
	     -1,-1,-1,-1, p4C.Px(),p4C.Py(),p4C.Pz(),p4C.E(), 0,0,0,0);
    new ((*particle_list)[1]) GHepParticle (kPdgHadronicBlob,kIStStableFinalState,
             -1,-1,-1,-1, p4R.Px(),p4R.Py(),p4R.Pz(),p4R.E(), 0,0,0,0);

    return particle_list;
  }

  // ....................................................................
  // Handle case where the remnant system is already known and doesn't
  // have to be hadronized. That happens when (close to the W threshold)
  // the hadronic system was generated by a simple 2-body decay
  //
  if(used_lowW_strategy) {
    // Create particle list (fragmentation record)
    TClonesArray * particle_list = new TClonesArray("genie::GHepParticle", 3);
    particle_list->SetOwner(true);

    // insert the generated particles
    new ((*particle_list)[0]) GHepParticle (ch_pdg,kIStStableFinalState,
	    -1,-1,-1,-1, p4C.Px(),p4C.Py(),p4C.Pz(),p4C.E(), 0,0,0,0);
    new ((*particle_list)[1]) GHepParticle (kPdgHadronicBlob,kIStNucleonTarget,
            -1,-1,2,2, p4R.Px(),p4R.Py(),p4R.Pz(),p4R.E(), 0,0,0,0);
    new ((*particle_list)[2]) GHepParticle (fs_nucleon_pdg,kIStStableFinalState,
            1,1,-1,-1, p4R.Px(),p4R.Py(),p4R.Pz(),p4R.E(), 0,0,0,0);

    return particle_list;
  }
  
  // ....................................................................
  // --------------------------------------------------------------------
  // Hadronize non-charm hadronic blob using PYTHIA/JETSET
  // --------------------------------------------------------------------
  // ....................................................................

  // Create output event record
  // Insert the generated charm hadron & the hadronic (non-charm) blob.
  // In this case the hadronic blob is entered as a pre-fragm. state.

  TClonesArray * particle_list = new TClonesArray("genie::GHepParticle");
  particle_list->SetOwner(true);

  new ((*particle_list)[0]) GHepParticle (ch_pdg,kIStStableFinalState,
          -1,-1,-1,-1,  p4C.Px(),p4C.Py(),p4C.Pz(),p4C.E(), 0,0,0,0);
  new ((*particle_list)[1]) GHepParticle (kPdgHadronicBlob,kIStNucleonTarget,
          -1,-1,2,3, p4R.Px(),p4R.Py(),p4R.Pz(),p4R.E(), 0,0,0,0);

#ifdef __GENIE_PYTHIA8_ENABLED__
  unsigned int rpos =2; // offset in event record

  bool use_pythia = (WR>1.5);

/*
  // Determining quark systems to input to PYTHIA based on simple quark model 
  // arguments
  //
  // Neutrinos
  // ------------------------------------------------------------------
  // Scattering off valence q
  // ..................................................................
  // p: [uu]+d 
  //         |--> c --> D0       <c+\bar(u)>  : [u]
  //                --> D+       <c+\bar(d)>  : [d]
  //                --> Ds+      <c+\bar(s)>  : [s]
  //                --> Lamda_c+ <c+ud     >  : [\bar(ud)]
  //
  // (for n: [uu] -> 50%[ud]_{0} + 50%[ud]_{1})
  //
  // Scattering off sea q
  // ..................................................................
  // p: [uud] + [\bar(d)]d (or)
  //            [\bar(s)]s 
  //                     |--> c --> D0       <c+\bar(u)>  : [u]
  //                            --> D+       <c+\bar(d)>  : [d]
  //                            --> Ds+      <c+\bar(s)>  : [s]
  //                            --> Lamda_c+ <c+ud     >  : [\bar(ud)]
  // Anti-Neutrinos
  // ------------------------------------------------------------------
  // Scattering off sea q
  // ..................................................................
  // p: [uud] + [d] \bar(d) (or)
  //            [s] \bar(s) 
  //                   |----> \bar(c) --> \bar(D0) <\bar(c)+u>  : [\bar(u)]
  //                                  --> D-       <\bar(c)+d>  : [\bar(d)]
  //                                  --> Ds-      <\bar(c)+s>  : [\bar(s)]
  // [Summary]
  //                                                                                    Qq
  // | v        + p  [val/d]        -->    D0       +  {           u          uu       }(+2) / u,uu
  // | v        + p  [val/d]        -->    D+       +  {           d          uu       }(+1) / d,uu
  // | v        + p  [val/d]        -->    Ds+      +  {           s          uu       }(+1) / s,uu
  // | v        + p  [val/d]        -->    Lc+      +  {           \bar(ud)   uu       }(+1) / \bar(d),u
  // | v        + n  [val/d]        -->    D0       +  {           u          ud       }(+1) / u,ud
  // | v        + n  [val/d]        -->    D+       +  {           d          ud       }( 0) / d,ud
  // | v        + n  [val/d]        -->    Ds+      +  {           s          ud       }( 0) / s,ud
  // | v        + n  [val/d]        -->    Lc+      +  {           \bar(ud)   ud       }( 0) / \bar(d),d
  // | v        + p  [sea/d]        -->    D0       +  {  uud      \bar(d)    u        }(+2) / u,uu
  // | v        + p  [sea/d]        -->    D+       +  {  uud      \bar(d)    d        }(+1) / d,uu
  // | v        + p  [sea/d]        -->    Ds+      +  {  uud      \bar(d)    s        }(+1) / s,uu
  // | v        + p  [sea/d]        -->    Lc+      +  {  uud      \bar(d)    \bar(ud) }(+1) / \bar(d),u
  // | v        + n  [sea/d]        -->    D0       +  {  udd      \bar(d)    u        }(+1) / u,ud
  // | v        + n  [sea/d]        -->    D+       +  {  udd      \bar(d)    d        }( 0) / d,ud
  // | v        + n  [sea/d]        -->    Ds+      +  {  udd      \bar(d)    s        }( 0) / s,ud
  // | v        + n  [sea/d]        -->    Lc+      +  {  udd      \bar(d)    \bar(ud) }( 0) / \bar(d),d
  // | v        + p  [sea/s]        -->    D0       +  {  uud      \bar(s)    u        }(+2) / u,uu
  // | v        + p  [sea/s]        -->    D+       +  {  uud      \bar(s)    d        }(+1) / d,uu
  // | v        + p  [sea/s]        -->    Ds+      +  {  uud      \bar(s)    s        }(+1) / s,uu
  // | v        + p  [sea/s]        -->    Lc+      +  {  uud      \bar(s)    \bar(ud) }(+1) / \bar(d),u
  // | v        + n  [sea/s]        -->    D0       +  {  udd      \bar(s)    u        }(+1) / u,ud
  // | v        + n  [sea/s]        -->    D+       +  {  udd      \bar(s)    d        }( 0) / d,ud
  // | v        + n  [sea/s]        -->    Ds+      +  {  udd      \bar(s)    s        }( 0) / s,ud
  // | v        + n  [sea/s]        -->    Lc+      +  {  udd      \bar(s)    \bar(ud) }( 0) / \bar(d),d

  // | \bar(v)  + p  [sea/\bar(d)]  -->    \bar(D0) +  {  uud      d          \bar(u)  }( 0) / d,ud
  // | \bar(v)  + p  [sea/\bar(d)]  -->    D-       +  {  uud      d          \bar(d)  }(+1) / d,uu
  // | \bar(v)  + p  [sea/\bar(d)]  -->    Ds-      +  {  uud      d          \bar(s)  }(+1) / d,uu
  // | \bar(v)  + n  [sea/\bar(d)]  -->    \bar(D0) +  {  udd      d          \bar(u)  }(-1) / d,dd
  // | \bar(v)  + n  [sea/\bar(d)]  -->    D-       +  {  udd      d          \bar(d)  }( 0) / d,ud
  // | \bar(v)  + n  [sea/\bar(d)]  -->    Ds-      +  {  udd      d          \bar(s)  }( 0) / d,ud
  // | \bar(v)  + p  [sea/\bar(s)]  -->    \bar(D0) +  {  uud      s          \bar(u)  }( 0) / d,ud
  // | \bar(v)  + p  [sea/\bar(s)]  -->    D-       +  {  uud      s          \bar(d)  }(+1) / d,uu
  // | \bar(v)  + p  [sea/\bar(s)]  -->    Ds-      +  {  uud      s          \bar(s)  }(+1) / d,uu
  // | \bar(v)  + n  [sea/\bar(s)]  -->    \bar(D0) +  {  udd      s          \bar(u)  }(-1) / d,dd
  // | \bar(v)  + n  [sea/\bar(s)]  -->    D-       +  {  udd      s          \bar(d)  }( 0) / d,ud
  // | \bar(v)  + n  [sea/\bar(s)]  -->    Ds-      +  {  udd      s          \bar(s)  }( 0) / d,ud
  //
  //
  // Taking some short-cuts below :
  // In the process of obtaining 2 q systems (while conserving the charge) I might tread 
  // d,s or \bar(d),\bar(s) as the same
  // In the future I should perform the first steps of the multi-q hadronization manualy
  // (to reduce the number of q's input to PYTHIA) or use py3ent_(), py4ent_() ...
  //
*/

  if(use_pythia) {
    int  qrkSyst1 = 0;
    int  qrkSyst2 = 0;
    if(isnu||isdm) { // neutrinos
       if(ch_pdg==kPdgLambdaPc) {
          if(isp) { qrkSyst1 = kPdgAntiDQuark; qrkSyst2 = kPdgUQuark; }
          if(isn) { qrkSyst1 = kPdgAntiDQuark; qrkSyst2 = kPdgDQuark; }
       } else {
          if(isp) { qrkSyst2 = kPdgUUDiquarkS1; }
          if(isn) { qrkSyst2 = (rnd->RndHadro().Rndm()<0.5) ?  kPdgUDDiquarkS0 :  kPdgUDDiquarkS1; }
          if (ch_pdg==kPdgD0      ) { qrkSyst1 = kPdgUQuark; }
          if (ch_pdg==kPdgDP      ) { qrkSyst1 = kPdgDQuark; }  
          if (ch_pdg==kPdgDPs     ) { qrkSyst1 = kPdgSQuark; }  
       }
     }
     if(isnub) { // antineutrinos
       qrkSyst1 = kPdgDQuark;    
       if (isp && ch_pdg==kPdgAntiD0) { qrkSyst2 = (rnd->RndHadro().Rndm()<0.5) ? kPdgUDDiquarkS0 : kPdgUDDiquarkS1; }
       if (isp && ch_pdg==kPdgDM    ) { qrkSyst2 = kPdgUUDiquarkS1; }
       if (isp && ch_pdg==kPdgDMs   ) { qrkSyst2 = kPdgUUDiquarkS1; }
       if (isn && ch_pdg==kPdgAntiD0) { qrkSyst2 = kPdgDDDiquarkS1; }
       if (isn && ch_pdg==kPdgDM    ) { qrkSyst2 = (rnd->RndHadro().Rndm()<0.5) ? kPdgUDDiquarkS0 : kPdgUDDiquarkS1; }
       if (isn && ch_pdg==kPdgDMs   ) { qrkSyst2 = (rnd->RndHadro().Rndm()<0.5) ? kPdgUDDiquarkS0 : kPdgUDDiquarkS1; }
     }
     if(qrkSyst1 == 0 && qrkSyst2 == 0) {
         LOG("CharmHad", pWARN) 
             << "Couldn't generate quark systems for PYTHIA in: " << *interaction;
         return 0;
     }

     //
     // Run PYTHIA for the hadronization of remnant system
     //
     fPythia->Pythia8()->particleData.mayDecay(kPdgPi0,              false); // don't decay pi0
     fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_DeltaM,  true); // decay Delta+
     fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_Delta0,  true); // decay Delta++
     fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_DeltaP,  true); // decay Delta++
     fPythia->Pythia8()->particleData.mayDecay(kPdgP33m1232_DeltaPP, true); // decay Delta++
//   fPythia->Pythia8()->particleData.mayDecay(kPdgDeltaP,  true); // decay Delta+
//   fPythia->Pythia8()->particleData.mayDecay(kPdgDeltaPP, true); // decay Delta++

     // -- hadronize --

     double mA    = fPythia->Pythia8()->particleData.m0(qrkSyst1);
     double mB    = fPythia->Pythia8()->particleData.m0(qrkSyst2);
     double pzAcm = 0.5 * Pythia8::sqrtpos( (WR + mA + mB) * (WR - mA - mB) * (WR - mA + mB) * (WR + mA - mB) ) / WR;
     double pzBcm = -pzAcm;
     double eA    = sqrt(mA*mA + pzAcm*pzAcm);
     double eB    = sqrt(mB*mB + pzBcm*pzBcm);

     fPythia->Pythia8()->event.reset();

      // Pythia8 status code for outgoing particles of the hardest subprocesses is 23
      // anti/colour tags for these 2 particles must complement each other
      // antiparticles must have positive anticolour to avoid PYTHIA errors
     if (qrkSyst1 > 0) {
         fPythia->Pythia8()->event.append(qrkSyst1, 23, 101, 0, 0., 0., pzAcm, eA, mA);
         fPythia->Pythia8()->event.append(qrkSyst2, 23, 0, 101, 0., 0., pzBcm, eB, mB);
     } else {
         fPythia->Pythia8()->event.append(qrkSyst1, 23, 0, 101, 0., 0., pzAcm, eA, mA);
         fPythia->Pythia8()->event.append(qrkSyst2, 23, 101, 0, 0., 0., pzBcm, eB, mB);
     }
     fPythia->Pythia8()->next();

      // List the event information
     fPythia->Pythia8()->event.list();
     fPythia->Pythia8()->stat();

     fPythia->Pythia8()->particleData.mayDecay(kPdgPi0, true); // restore

     //-- Get PYTHIA's event record
     Pythia8::Event &fEvent = fPythia->Pythia8()->event;
     int numpart = fEvent.size();
     assert(numpart>0);

     // Offset the initial (system) particle
     int ioff = 0;
     if (fEvent[0].id() == 90) ioff = -1;

     // Hadronic 4vec
     TLorentzVector p4Had = kinematics.HadSystP4();
 
     // Vector defining rotation from LAB to LAB' (z:= \vec{phad})
     TVector3 unitvq = p4Had.Vect().Unit();
 
     // Boost velocity LAB' -> HCM
     TVector3 beta(0,0,p4Had.P()/p4Had.Energy());

     for (int i = 1; i < numpart; ++i) {
        /*
         * Convert Pythia8 status code to Pythia6
         * Initial quark has a pythia6 status code of 12
         * The initial diquark and the fragmented particles have a pythia6 code
         * of 11 (kIStNucleonTarget)
         * Final state particles have a positive pythia8 code and a pythia6 code of
         * 1 (kIStStableFinalState)
         * The fragmentation products are generated in the hadronic CM frame
         * where the z>0 axis is the \vec{phad} direction. For each particle 
         * returned by the hadronizer:
         * - boost it back to LAB' frame {z:=\vec{phad}} / doesn't affect pT
         * - rotate its 3-momentum from LAB' to LAB
         */
        GHepStatus_t gStatus;
        if (i == 1) gStatus = kIStDISPreFragmHadronicState;
        else gStatus = (fEvent[i].status()>0) ? kIStStableFinalState : kIStNucleonTarget;

        LOG("CharmP8Had", pDEBUG)
            << "Adding final state particle pdgc = " << fEvent[i].id()
            << " with status = " << gStatus;

        if (fEvent[i].status() > 0){
          if( pdg::IsQuark  (fEvent[i].id()) || 
                  pdg::IsDiQuark(fEvent[i].id()) ) {
            LOG("CharmP8Had", pERROR)
                << "Hadronization failed! Bare quark/di-quarks appear in final state!";
            particle_list->Delete();
            delete particle_list;
            return 0;            
          }
        }

        TLorentzVector p4o(fEvent[i].px(), fEvent[i].py(), fEvent[i].pz(), fEvent[i].e());
        p4o.Boost(beta); 
        TVector3 p3 = p4o.Vect();
        p3.RotateUz(unitvq); 
        TLorentzVector p4(p3,p4o.Energy());

        new((*particle_list)[i++]) GHepParticle(
                fEvent[i].id(),               // pdg
                gStatus,                      // status
                fEvent[i].mother1()   + ioff, // first parent
                fEvent[i].mother2()   + ioff, // second parent
                fEvent[i].daughter1() + ioff, // first daughter
                fEvent[i].daughter2() + ioff, // second daughter
                p4.Px(),                      // px [GeV/c]
                p4.Py(),                      // py [GeV/c]
                p4.Pz(),                      // pz [GeV/c]
                p4.Energy(),                  // e  [GeV]
                fEvent[i].xProd(),            // x  [mm]
                fEvent[i].yProd(),            // y  [mm]
                fEvent[i].zProd(),            // z  [mm]
                fEvent[i].tProd()             // t  [mm/c]
        );
     }
  } // use_pythia

  // ....................................................................
  // Hadronizing low-W non-charm hadronic blob using a phase space decay
  // ....................................................................

  else {
     // Just a small fraction of events (low-W remnant syste) causing trouble
     // to PYTHIA/JETSET
     // Set a 2-body N+pi system that matches the remnant system charge and
     // do a simple phase space decay
     //
     // q(remn)  remn/syst
     // +2    :  (p pi+)
     // +1    :  50%(p pi0) + 50% n pi+
     //  0    :  50%(p pi-) + 50% n pi0
     // -1    :  (n pi-)
     //
     double qfsl  = interaction->FSPrimLepton()->Charge() / 3.;
     double qinit = pdglib->Find(nuc_pdg)->Charge() / 3.;
     double qch   = pdglib->Find(ch_pdg)->Charge() / 3.;
     int Q = (int) (qinit - qfsl - qch); // remnant hadronic system charge

     bool allowdup=true;
     PDGCodeList pd(allowdup);
     if(Q==2) {
            pd.push_back(kPdgProton);  pd.push_back(kPdgPiP);  } 
     else if (Q==1) {
       if(rnd->RndHadro().Rndm()<0.5) {
            pd.push_back(kPdgProton);  pd.push_back(kPdgPi0);  }
       else {
            pd.push_back(kPdgNeutron); pd.push_back(kPdgPiP);  }
     } 
     else if (Q==0) {
        if(rnd->RndHadro().Rndm()<0.5) {
            pd.push_back(kPdgProton);  pd.push_back(kPdgPiM);  }
        else {
           pd.push_back(kPdgNeutron);  pd.push_back(kPdgPi0);  }
     } 
     else if (Q==-1) {
           pd.push_back(kPdgNeutron);  pd.push_back(kPdgPiM);  }

     double mass[2] = {
       pdglib->Find(pd[0])->Mass(), pdglib->Find(pd[1])->Mass()
     };

     // Set the decay
     bool permitted = fPhaseSpaceGenerator.SetDecay(p4R, 2, mass);
     if(!permitted) {
       LOG("CharmHad", pERROR) << " *** Phase space decay is not permitted";
       return 0;
     }
     // Get the maximum weight
     double wmax = -1;
     for(int i=0; i<200; i++) {
       double w = fPhaseSpaceGenerator.Generate();
       wmax = TMath::Max(wmax,w);
     }
     if(wmax<=0) {
       LOG("CharmHad", pERROR) << " *** Non-positive maximum weight";
       LOG("CharmHad", pERROR) << " *** Can not generate an unweighted phase space decay";
       return 0;
     }

     LOG("CharmHad", pINFO)
        << "Max phase space gen. weight @ current hadronic system: " << wmax;

      // *** generating an un-weighted decay ***
      wmax *= 1.3;
      bool accept_decay=false;
      unsigned int idectry=0;
      while(!accept_decay)
      {
       idectry++;
       if(idectry>kMaxUnweightDecayIterations) {
         // report, clean-up and return
         LOG("Char,Had", pWARN)
             << "Couldn't generate an unweighted phase space decay after "
             << itry << " attempts";
         return 0;
      }
      double w = fPhaseSpaceGenerator.Generate();
      if(w > wmax) {
          LOG("CharmHad", pWARN)
             << "Decay weight = " << w << " > max decay weight = " << wmax;
       }
       double gw = wmax * rnd->RndHadro().Rndm();
       accept_decay = (gw<=w);
       LOG("CharmHad", pDEBUG)
          << "Decay weight = " << w << " / R = " << gw << " - accepted: " << accept_decay;
     }
     for(unsigned int i=0; i<2; i++) {
        int pdgc = pd[i];
        TLorentzVector * p4d = fPhaseSpaceGenerator.GetDecay(i);
        new ( (*particle_list)[rpos+i] ) GHepParticle(
           pdgc,kIStStableFinalState,1,1,-1,-1,p4d->Px(),p4d->Py(),p4d->Pz(),p4d->Energy(),
           0,0,0,0);
     }
  } 

#endif
  //-- Print & return the fragmentation record
  //utils::fragmrec::Print(particle_list);
  return particle_list;
}
//____________________________________________________________________________
int CharmP8Hadronization::GenerateCharmHadron(int nu_pdg, double EvLab) const
{
  // generate a charmed hadron pdg code using a charm fraction table

  RandomGen * rnd = RandomGen::Instance();
  double r  = rnd->RndHadro().Rndm();

  if(pdg::IsNeutrino(nu_pdg)) {
     double tf = 0;
     if      (r < (tf+=fD0FracSpl->Evaluate(EvLab)))  return kPdgD0;       // D^0
     else if (r < (tf+=fDpFracSpl->Evaluate(EvLab)))  return kPdgDP;       // D^+
     else if (r < (tf+=fDsFracSpl->Evaluate(EvLab)))  return kPdgDPs;      // Ds^+
     else                                             return kPdgLambdaPc; // Lamda_c^+

  } else if(pdg::IsAntiNeutrino(nu_pdg)) {
     if      (r < fD0BarFrac)          return kPdgAntiD0;
     else if (r < fD0BarFrac+fDmFrac)  return kPdgDM;
     else                              return kPdgDMs;
  }

  LOG("CharmHad", pERROR) << "Could not generate a charm hadron!";
  return 0;
}
//____________________________________________________________________________
double CharmP8Hadronization::Weight(void) const 
{
  return 1. ;
}

//____________________________________________________________________________
void CharmP8Hadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CharmP8Hadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CharmP8Hadronization::LoadConfig(void)
{

  bool hadronize_remnants ; 
  GetParamDef( "HadronizeRemnants", hadronize_remnants, true ) ;

  fCharmOnly = ! hadronize_remnants ;

  //-- Get a fragmentation function
  fFragmFunc = dynamic_cast<const FragmentationFunctionI *> (
    this->SubAlg("FragmentationFunc"));
  assert(fFragmFunc);

  string pt_function ;
  this -> GetParam( "PTFunction", pt_function ) ;

  fCharmPT2pdf = new TF1("fCharmPT2pdf", pt_function.c_str(),0,0.6);

  // stop ROOT from deleting this object of its own volition
  gROOT->GetListOfFunctions()->Remove(fCharmPT2pdf);

  // neutrino charm fractions: D^0, D^+, Ds^+ (remainder: Lamda_c^+)
  std::vector<double> ec, d0frac, dpfrac, dsfrac ;

  std::string raw ;
  std::vector<std::string> bits ;

  bool invalid_configuration = false ;

  // load energy points
  this -> GetParam( "CharmFrac-E", raw ) ;
  bits = utils::str::Split( raw, ";" ) ;

  if ( ! utils::str::Convert(bits, ec) ) {
    LOG("CharmP8Hadronization", pFATAL) <<
    		"Failed to decode CharmFrac-E string: ";
    LOG("CharmP8Hadronization", pFATAL) << "string: "<< raw ;
    invalid_configuration = true ;
  }

  // load D0 fractions
  this -> GetParam( "CharmFrac-D0", raw ) ;
  bits = utils::str::Split( raw, ";" ) ;

  if ( ! utils::str::Convert(bits, d0frac) ) {
    LOG("CharmP8Hadronization", pFATAL) <<
    		"Failed to decode CharmFrac-D0 string: ";
    LOG("CharmP8Hadronization", pFATAL) << "string: "<< raw ;
    invalid_configuration = true ;
  }

  // check the size
  if ( d0frac.size() != ec.size() ) {
	  LOG("CharmP8Hadronization", pFATAL) << "E entries don't match D0 fraction entries";
	  LOG("CharmP8Hadronization", pFATAL) << "E:  " << ec.size() ;
	  LOG("CharmP8Hadronization", pFATAL) << "D0: " << d0frac.size() ;
	  invalid_configuration = true ;
  }

  // load D+ fractions
    this -> GetParam( "CharmFrac-D+", raw ) ;
    bits = utils::str::Split( raw, ";" ) ;

    if ( ! utils::str::Convert(bits, dpfrac) ) {
      LOG("CharmP8Hadronization", pFATAL) <<
      		"Failed to decode CharmFrac-D+ string: ";
      LOG("CharmP8Hadronization", pFATAL) << "string: "<< raw ;
      invalid_configuration = true ;
    }

    // check the size
    if ( dpfrac.size() != ec.size() ) {
  	  LOG("CharmP8Hadronization", pFATAL) << "E entries don't match D+ fraction entries";
  	  LOG("CharmP8Hadronization", pFATAL) << "E:  " << ec.size() ;
  	  LOG("CharmP8Hadronization", pFATAL) << "D+: " << dpfrac.size() ;
  	  invalid_configuration = true ;
    }

    // load D_s fractions
    this -> GetParam( "CharmFrac-Ds", raw ) ;
    bits = utils::str::Split( raw, ";" ) ;

    if ( ! utils::str::Convert(bits, dsfrac) ) {
    	LOG("CharmP8Hadronization", pFATAL) <<
    			"Failed to decode CharmFrac-Ds string: ";
    	LOG("CharmP8Hadronization", pFATAL) << "string: "<< raw ;
    	invalid_configuration = true ;
    }

    // check the size
    if ( dsfrac.size() != ec.size() ) {
    	LOG("CharmP8Hadronization", pFATAL) << "E entries don't match Ds fraction entries";
    	LOG("CharmP8Hadronization", pFATAL) << "E:  " << ec.size() ;
    	LOG("CharmP8Hadronization", pFATAL) << "Ds: " << dsfrac.size() ;
    	invalid_configuration = true ;
    }

  fD0FracSpl = new Spline( ec.size(), & ec[0], & d0frac[0] );
  fDpFracSpl = new Spline( ec.size(), & ec[0], & dpfrac[0] );
  fDsFracSpl = new Spline( ec.size(), & ec[0], & dsfrac[0] );

  // anti-neutrino charm fractions: bar(D^0), D^-, (remainder: Ds^-)

  this -> GetParam( "CharmFrac-D0bar", fD0BarFrac ) ;
  this -> GetParam( "CharmFrac-D-",    fDmFrac ) ;

  if ( invalid_configuration ) {

	    LOG("CharmP8Hadronization", pFATAL)
	      << "Invalid configuration: Exiting" ;

	    // From the FreeBSD Library Functions Manual
	    //
	    // EX_CONFIG (78)   Something was found in an unconfigured or miscon-
	    //                  figured state.

	    exit( 78 ) ;
  }
}
//____________________________________________________________________________
