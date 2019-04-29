//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)

*/
//____________________________________________________________________________

#include <RVersion.h>
#if ROOT_VERSION_CODE >= ROOT_VERSION(5,15,6)
#include <TMCParticle.h>
#else
#include <TMCParticle6.h>
#endif
#include <TPythia6.h>
#include <TVector3.h>
#include <TF1.h>
#include <TROOT.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Controls.h"
#include "Physics/Hadronization/CharmHadronization.h"
#include "Physics/Hadronization/FragmentationFunctionI.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/Utils/KineUtils.h"
#include "Physics/Hadronization/FragmRecUtils.h"
#include "Framework/Utils/PrintUtils.h"

using namespace genie;
using namespace genie::constants;
using namespace genie::controls;

extern "C" void py1ent_(int *,  int *, double *, double *, double *);
extern "C" void py2ent_(int *,  int *, int *, double *);

//____________________________________________________________________________
CharmHadronization::CharmHadronization() :
HadronizationModelI("genie::CharmHadronization")
{
  this->Initialize();
}
//____________________________________________________________________________
CharmHadronization::CharmHadronization(string config) :
HadronizationModelI("genie::CharmHadronization", config)
{
  this->Initialize();
}
//____________________________________________________________________________
CharmHadronization::~CharmHadronization()
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
void CharmHadronization::Initialize(void) const
{
  fPythia = TPythia6::Instance();
}
//____________________________________________________________________________
TClonesArray * CharmHadronization::Hadronize(
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
  double MC = pdglib->Find(ch_pdg)->Mass();

  LOG("CharmHad", pNOTICE) << "Remnant hadronic system mass = " << WR;

  // ....................................................................
  // Handle case where the user doesn't want to remnant system to be
  // hadronized (add as 'hadronic blob') 
  //
  if(fCharmOnly) {
    // Create particle list (fragmentation record)
    TClonesArray * particle_list = new TClonesArray("TMCParticle", 2);
    particle_list->SetOwner(true);

    // insert the generated particles
    new ((*particle_list)[0]) TMCParticle (1,ch_pdg,
	     -1,-1,-1, p4C.Px(),p4C.Py(),p4C.Pz(),p4C.E(),MC, 0,0,0,0,0);
    new ((*particle_list)[1]) TMCParticle (1,kPdgHadronicBlob,
             -1,-1,-1, p4R.Px(),p4R.Py(),p4R.Pz(),p4R.E(),WR, 0,0,0,0,0);

    return particle_list;
  }

  // ....................................................................
  // Handle case where the remnant system is already known and doesn't
  // have to be hadronized. That happens when (close to the W threshold)
  // the hadronic system was generated by a simple 2-body decay
  //
  if(used_lowW_strategy) {
    // Create particle list (fragmentation record)
    TClonesArray * particle_list = new TClonesArray("TMCParticle", 3);
    particle_list->SetOwner(true);

    // insert the generated particles
    new ((*particle_list)[0]) TMCParticle (1,ch_pdg,
	     -1,-1,-1, p4C.Px(),p4C.Py(),p4C.Pz(),p4C.E(),MC, 0,0,0,0,0);
    new ((*particle_list)[1]) TMCParticle (11,kPdgHadronicBlob,
               -1,2,2, p4R.Px(),p4R.Py(),p4R.Pz(),p4R.E(),WR, 0,0,0,0,0);
    new ((*particle_list)[2]) TMCParticle (1,fs_nucleon_pdg,
              1,-1,-1, p4R.Px(),p4R.Py(),p4R.Pz(),p4R.E(),WR, 0,0,0,0,0);

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

  TClonesArray * particle_list = new TClonesArray("TMCParticle");
  particle_list->SetOwner(true);

  new ((*particle_list)[0]) TMCParticle (1,ch_pdg,
                  -1,-1,-1,  p4C.Px(),p4C.Py(),p4C.Pz(),p4C.E(),MC, 0,0,0,0,0);
  new ((*particle_list)[1]) TMCParticle (11,kPdgHadronicBlob,
                     -1,2,3, p4R.Px(),p4R.Py(),p4R.Pz(),p4R.E(),WR, 0,0,0,0,0);

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
     fPythia->SetMDCY(fPythia->Pycomp(kPdgPi0),              1,0); // don't decay pi0
     fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaM),  1,1); // decay Delta+
     fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_Delta0),  1,1); // decay Delta++
     fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaP),  1,1); // decay Delta++
     fPythia->SetMDCY(fPythia->Pycomp(kPdgP33m1232_DeltaPP), 1,1); // decay Delta++
//   fPythia->SetMDCY(fPythia->Pycomp(kPdgDeltaP),  1,1); // decay Delta+
//   fPythia->SetMDCY(fPythia->Pycomp(kPdgDeltaPP), 1,1); // decay Delta++

     int ip = 0;
     py2ent_(&ip, &qrkSyst1, &qrkSyst2, &WR); // hadronize

     fPythia->SetMDCY(fPythia->Pycomp(kPdgPi0),1,1); // restore

     //-- Get PYTHIA's LUJETS event record
     TClonesArray * remnants = 0;
     fPythia->GetPrimaries();
     remnants = dynamic_cast<TClonesArray *>(fPythia->ImportParticles("All"));
     if(!remnants) {
         LOG("CharmHad", pWARN) << "Couldn't hadronize (non-charm) remnants!";
         return 0;
      }

      // PYTHIA performs the hadronization at the *remnant hadrons* centre of mass 
      // frame  (not the hadronic centre of mass frame). 
      // Boost all hadronic blob fragments to the HCM', fix their mother/daughter 
      // assignments and add them to the fragmentation record.

      TVector3 rmnbeta = +1 * p4R.BoostVector(); // boost velocity

      TMCParticle * remn  = 0; // remnant
      TMCParticle * bremn = 0; // boosted remnant
      TIter remn_iter(remnants);
      while( (remn = (TMCParticle *) remn_iter.Next()) ) {

         // insert and get a pointer to inserted object for mods
         bremn = new ((*particle_list)[rpos++]) TMCParticle (*remn);

         // boost 
         TLorentzVector p4(remn->GetPx(),remn->GetPy(),remn->GetPz(),remn->GetEnergy());
         p4.Boost(rmnbeta);
         bremn -> SetPx     (p4.Px());
         bremn -> SetPy     (p4.Py());
         bremn -> SetPz     (p4.Pz());
         bremn -> SetEnergy (p4.E() );

         // handle insertion of charmed hadron 
         int jp  = bremn->GetParent();
         int ifc = bremn->GetFirstChild();
         int ilc = bremn->GetLastChild();
         bremn -> SetParent     ( (jp  == 0 ?  1 : jp +1) );
         bremn -> SetFirstChild ( (ifc == 0 ? -1 : ifc+1) );
         bremn -> SetLastChild  ( (ilc == 0 ? -1 : ilc+1) );
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
        new ( (*particle_list)[rpos+i] ) TMCParticle(
           1,pdgc,1,-1,-1,p4d->Px(),p4d->Py(),p4d->Pz(),p4d->Energy(),
           mass[i],0,0,0,0,0);
     }
  } 

  //-- Print & return the fragmentation record
  //utils::fragmrec::Print(particle_list);
  return particle_list;
}
//____________________________________________________________________________
double CharmHadronization::Weight(void) const
{
  return 1.; // does not generate weighted events
}
//____________________________________________________________________________
PDGCodeList * CharmHadronization::SelectParticles(
                               const Interaction * /*interaction*/) const
{
  return 0;
}
//____________________________________________________________________________
TH1D * CharmHadronization::MultiplicityProb(
           const Interaction * /*interaction*/, Option_t * /*opt*/)  const
{
  return 0;
}
//____________________________________________________________________________
int CharmHadronization::GenerateCharmHadron(int nu_pdg, double EvLab) const
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
void CharmHadronization::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CharmHadronization::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void CharmHadronization::LoadConfig(void)
{

  bool hadronize_remnants = true ;
  GetParam( "HadronizeRemnants", hadronize_remnants, false ) ;

  fCharmOnly = ! hadronize_remnants ;

  //-- Get a fragmentation function
  fFragmFunc = dynamic_cast<const FragmentationFunctionI *> (
    this->SubAlg("FragmentationFunc"));
  assert(fFragmFunc);

  fCharmPT2pdf = new TF1("fCharmPT2pdf", "exp(-0.213362-6.62464*x)",0,0.6);
  // stop ROOT from deleting this object of its own volition
  gROOT->GetListOfFunctions()->Remove(fCharmPT2pdf);

  // neutrino charm fractions: D^0, D^+, Ds^+ (remainder: Lamda_c^+)
  //
  const int nc = 15;
  double ec[nc] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,50.,60.,70.,80.,90.,100.};

  double d0frac[nc] = { .000, .320, .460,  .500,  .520,  .530, .540, .540,
                        .540, .550,  .550,  .560,  .570, .580, .600  };
  double dpfrac[nc] = { .000, .120, .180,  .200,  .200,  .210, .210, .210,
                        .210, .210,  .220,  .220,  .220, .230, .230  };
  double dsfrac[nc] = { .000, .054, .078,  .130,  .130,  .140, .140, .140,
                        .140, .140,  .140,  .140,  .140, .150, .150  };
 
  fD0FracSpl = new Spline(nc, ec, d0frac);
  fDpFracSpl = new Spline(nc, ec, dpfrac);
  fDsFracSpl = new Spline(nc, ec, dsfrac);

  // anti-neutrino charm fractions: bar(D^0), D^-, (remainder: Ds^-)
  //
  fD0BarFrac = 0.667;
  fDmFrac    = 0.222;
}
//____________________________________________________________________________
