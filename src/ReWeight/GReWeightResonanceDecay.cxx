//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 10, 2010 - CA
   First included in v2.7.1.

*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtupleD.h>
#include <TH1D.h>
#include <TParticlePDG.h>
#include <TDecayChannel.h>

#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGLibrary.h"
#include "ReWeight/GReWeightResonanceDecay.h"
#include "ReWeight/GSystUncertainty.h"

//#define _G_REWEIGHT_RESDEC_DEBUG_

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightResonanceDecay::GReWeightResonanceDecay() :
GReWeightI()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightResonanceDecay::~GReWeightResonanceDecay()
{
#ifdef _G_REWEIGHT_RESDEC_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;  
#endif
}
//_______________________________________________________________________________________
bool GReWeightResonanceDecay::IsHandled(GSyst_t syst)
{
  switch(syst) {
     case ( kRDcyTwkDial_BR1gamma         ) :
     case ( kRDcyTwkDial_BR1eta           ) :
     case ( kRDcyTwkDial_Theta_Delta2Npi  ) :
        return true;
        break;
     default:
        return false;
        break;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightResonanceDecay::SetSystematic(GSyst_t syst, double twk_dial)
{
  switch(syst) {
     case ( kRDcyTwkDial_BR1gamma ) :
        fBR1gammaTwkDial = twk_dial;
        break;
     case ( kRDcyTwkDial_BR1eta ) :
        fBR1etaTwkDial = twk_dial;
        break;
     case ( kRDcyTwkDial_Theta_Delta2Npi  ) :
        fThetaDelta2NpiTwkDial = twk_dial;
        break;
     default:
        return;
        break;
  }
}
//_______________________________________________________________________________________
void GReWeightResonanceDecay::Reset(void)
{
  fBR1gammaTwkDial       = 0.;
  fBR1etaTwkDial         = 0.;
  fThetaDelta2NpiTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightResonanceDecay::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightResonanceDecay::CalcWeight(const EventRecord & event) 
{  
  Interaction * interaction = event.Summary();

  bool is_res = interaction->ProcInfo().IsResonant();
  if(!is_res) return 1.;

  bool is_cc  = interaction->ProcInfo().IsWeakCC();
  bool is_nc  = interaction->ProcInfo().IsWeakNC();
  if(is_cc && !fRewCC) return 1.;
  if(is_nc && !fRewNC) return 1.;

  int nupdg = interaction->InitState().ProbePdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  double wght = 
    this->RewBR(event) * 
    this->RewThetaDelta2Npi(event);

  return wght;
}
//_______________________________________________________________________________________
double GReWeightResonanceDecay::CalcChisq(void)
{
  double chisq =
    TMath::Power(fBR1gammaTwkDial,       2.) +
    TMath::Power(fBR1etaTwkDial,         2.) +
    TMath::Power(fThetaDelta2NpiTwkDial, 2.);

  return chisq;
}
//_______________________________________________________________________________________
double GReWeightResonanceDecay::RewBR(const EventRecord & event)
{
  bool tweaked_br1gamma = (TMath::Abs(fBR1gammaTwkDial) > controls::kASmallNum);
  bool tweaked_br1eta   = (TMath::Abs(fBR1etaTwkDial)   > controls::kASmallNum);

  bool tweaked = (tweaked_br1gamma || tweaked_br1eta);
  if(!tweaked) return 1.;

  double wght = 1.;

  GSystUncertainty * uncertainty = GSystUncertainty::Instance();

  LOG("ReW", pDEBUG) << "Checking resonance decay mode.";

  GHepParticle * p = 0;
  TIter iter(&event);
  while((p=(GHepParticle*)iter.Next())) {
    if(utils::res::IsBaryonResonance(p->Pdg())) {
      int fd = p->FirstDaughter();
      int ld = p->LastDaughter();
      int ngamma = 0;
      int neta   = 0;
      for(int i=fd; i<=ld; i++) {
        int dpdg = event.Particle(i)->Pdg();
        if(dpdg==kPdgGamma) ngamma++; 
        if(dpdg==kPdgEta  ) neta++; 
      }//i
      bool is_1gamma = (ngamma==1);
      bool is_1eta   = (neta  ==1);

      if(is_1gamma) {
         LOG("ReW", pDEBUG) << "A resonance -> X + 1gamma event";
      }
      if(is_1eta) {
         LOG("ReW", pDEBUG) << "A resonance -> X + 1eta event";
      }

      //
      // For Resonance -> X + gamma, tweak the branching ratio as:
      // BR_{tweaked} = BR_{default} * (1 + dial * fractional_error)
      // so that, basically, weight = 1 + dial * fractional_error.
      // For all other decay modes adjust the branching ratio so that Sum{BR} = 1.
      //
      if(tweaked_br1gamma) {
        double frerr = uncertainty->OneSigmaErr(kRDcyTwkDial_BR1gamma);
        double dial  = fBR1gammaTwkDial;
        double w = (1. + dial*frerr);
        w = TMath::Max(0.,w);
        TH1D * brfw  = fMpBR1gammaDef[p->Pdg()];
        double mass = p->P4()->M();
        mass = TMath::Min(mass, brfw->GetXaxis()->GetXmax());
        int ibin = brfw->FindBin(mass);
        double brdef = brfw->GetBinContent(ibin);
        double brtwk = brdef*w;
        if(brtwk>1) {
         brtwk = 1.;
         w = brtwk/brdef;
        }
        if(is_1gamma) {
         wght *= w;
        } else {         
         wght *= ((1-brtwk)/(1-brdef));
        }
      }
      // Similarly for Resonance -> X + eta
      //
      if(tweaked_br1eta) {
        double frerr = uncertainty->OneSigmaErr(kRDcyTwkDial_BR1eta);
        double dial  = fBR1etaTwkDial;
        double w = (1. + dial*frerr);
        w = TMath::Max(0.,w);
        TH1D * brfw  = fMpBR1etaDef[p->Pdg()];
        double mass = p->P4()->M();
        mass = TMath::Min(mass, brfw->GetXaxis()->GetXmax());
        int ibin = brfw->FindBin(mass);
        double brdef = brfw->GetBinContent(ibin);
        double brtwk = brdef*w;
        if(brtwk>1) {
         brtwk = 1.;
         w = brtwk/brdef;
        }
        if(is_1eta) {
         wght *= w;
        } else {         
         wght *= ((1-brtwk)/(1-brdef));
        }
      }
      // Similarly for other modes with tweaked BRs
      //
      //...
      //...

    }//res?
  }//p

  return wght;
}
//_______________________________________________________________________________________
double GReWeightResonanceDecay::RewThetaDelta2Npi(const EventRecord & event)
{
  bool tweaked = (TMath::Abs(fThetaDelta2NpiTwkDial) > controls::kASmallNum);
  if(!tweaked) return 1.;

  bool is_Delta_1pi = false;
  int ir = -1; // resonance position
  int ip = -1; // pion position
  int i  =  0;
  GHepParticle * p = 0;
  TIter iter(&event);
  while((p=(GHepParticle*)iter.Next())) {
    bool is_Deltapp = (p->Pdg()==kPdgP33m1232_DeltaPP);
    if(is_Deltapp) {
      ir = i;
      int fd = p->FirstDaughter();
      int ld = p->LastDaughter();
      int nd = 1 + ld - fd;
      if(nd==2) {
        int fpdg = event.Particle(fd)->Pdg();
        int lpdg = event.Particle(ld)->Pdg();
        if(fpdg==kPdgProton && lpdg==kPdgPiP) { is_Delta_1pi = true; ip = ld; }
        if(lpdg==kPdgProton && fpdg==kPdgPiP) { is_Delta_1pi = true; ip = fd; }
      }
    }

    bool is_Deltap = (p->Pdg()==kPdgP33m1232_DeltaP);
    if(is_Deltap) {
      ir = i;
      int fd = p->FirstDaughter();
      int ld = p->LastDaughter();
      int nd = 1 + ld - fd;
      if(nd==2) {
        int fpdg = event.Particle(fd)->Pdg();
        int lpdg = event.Particle(ld)->Pdg();
        if(fpdg==kPdgProton && lpdg==kPdgPi0) { is_Delta_1pi = true; ip = ld; }
        if(lpdg==kPdgProton && fpdg==kPdgPi0) { is_Delta_1pi = true; ip = fd; }
        if(fpdg==kPdgNeutron && lpdg==kPdgPiP) { is_Delta_1pi = true; ip = ld; }
        if(lpdg==kPdgNeutron && fpdg==kPdgPiP) { is_Delta_1pi = true; ip = fd; }
      }
    }

    bool is_Delta0 = (p->Pdg()==kPdgP33m1232_Delta0);
    if(is_Delta0) {
      ir = i;
      int fd = p->FirstDaughter();
      int ld = p->LastDaughter();
      int nd = 1 + ld - fd;
      if(nd==2) {
        int fpdg = event.Particle(fd)->Pdg();
        int lpdg = event.Particle(ld)->Pdg();
        if(fpdg==kPdgProton && lpdg==kPdgPiM) { is_Delta_1pi = true; ip = ld; }
        if(lpdg==kPdgProton && fpdg==kPdgPiM) { is_Delta_1pi = true; ip = fd; }
        if(fpdg==kPdgNeutron && lpdg==kPdgPi0) { is_Delta_1pi = true; ip = ld; }
        if(lpdg==kPdgNeutron && fpdg==kPdgPi0) { is_Delta_1pi = true; ip = fd; }
      }
    }

    bool is_DeltaM = (p->Pdg()==kPdgP33m1232_DeltaM);
    if(is_DeltaM) {
      ir = i;
      int fd = p->FirstDaughter();
      int ld = p->LastDaughter();
      int nd = 1 + ld - fd;
      if(nd==2) {
        int fpdg = event.Particle(fd)->Pdg();
        int lpdg = event.Particle(ld)->Pdg();
        if(fpdg==kPdgNeutron && lpdg==kPdgPiM) { is_Delta_1pi = true; ip = ld; }
        if(lpdg==kPdgNeutron && fpdg==kPdgPiM) { is_Delta_1pi = true; ip = fd; }
      }
    }

    if(is_Delta_1pi) break;
    i++;
  }

  if(!is_Delta_1pi) return 1.;

  LOG("ReW", pDEBUG) << "A Delta++ -> p pi+ event:";
  LOG("ReW", pDEBUG) << "Resonance is at position: " << ir;
  LOG("ReW", pDEBUG) << "Pion is at position: " << ip;
//LOG("ReW", pDEBUG) << event;

  // Get Delta and pi+ 4-momentum vectors
  TLorentzVector p4res(*event.Particle(ir)->P4());
  TLorentzVector p4pip(*event.Particle(ip)->P4());

  // Boost pi+ to the Delta CM
  TVector3 bv = -1*p4res.BoostVector();
  p4pip.Boost(bv);

  //
  // Calculate weight.
  // Angular distribution for pi+, given 2 possible sub-states mi=1/2,3/2
  // is: Wpi(theta) = 1 - p(3/2)*P2(costheta) + p(1/2)*P2(costheta).
  // mi is the projection of Delta angular momentum along quantization axis.
  // P2(costheta) is the 2nd order Legendre polynomial.
  // Isotropy implies equal populations: p(3/2)=p(1/2)=0.5
  // R/S predict: p(3/2)=0.75, p(1/2)=0.25
  // The above follow the notation used in Sam Zeller's invited talk in a 
  // MINOS Physics Simulation mtg on Sep 28th, 2007.
  // Sam's talk cited a paper by G.Garvey and O.Lalakulich which I haven't
  // been able to locate.
  //

  double p32iso   = 0.50;
  double p12iso   = 0.50;
  double p32rs    = 0.75;
  double p12rs    = 0.25;
  double costheta = p4pip.Vect().CosTheta();
  double P2       = 0.5 * (3.*costheta*costheta - 1.);
  double Wiso     = 1 - p32iso * P2 + p12iso * P2; // = 1.0
  double Wrs      = 1 - p32rs  * P2 + p12rs  * P2;
  double dial     = fThetaDelta2NpiTwkDial;
  double Wdef     = Wiso;
  double Wtwk     = dial*Wrs + (1-dial)*Wiso;

  double wght = 1.;
  if(Wdef>0. && Wtwk>0.) {
     wght = Wtwk/Wdef;
  }

  LOG("ReW", pDEBUG) 
       << "Pion Cos(ThetaCM) = " << costheta << ", weight = " << wght;

#ifdef _G_REWEIGHT_RESDEC_DEBUG_
  fTestNtp->Fill(costheta,wght);
#endif

  return wght;
}
//_______________________________________________________________________________________
void GReWeightResonanceDecay::Init(void)
{
  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);
  this->RewCC     (true);
  this->RewNC     (true);

  fBR1gammaTwkDial       = 0.;
  fBR1etaTwkDial         = 0.;
  fThetaDelta2NpiTwkDial = 0.;

  const int respdgarray[] = {
    kPdgP33m1232_DeltaM, kPdgP33m1232_Delta0, kPdgP33m1232_DeltaP, kPdgP33m1232_DeltaPP,
    kPdgS31m1620_DeltaM, kPdgS31m1620_Delta0, kPdgS31m1620_DeltaP, kPdgS31m1620_DeltaPP,
    kPdgD33m1700_DeltaM, kPdgD33m1700_Delta0, kPdgD33m1700_DeltaP, kPdgD33m1700_DeltaPP,
    kPdgF35m1905_DeltaM, kPdgF35m1905_Delta0, kPdgF35m1905_DeltaP, kPdgF35m1905_DeltaPP,
    kPdgP31m1910_DeltaM, kPdgP31m1910_Delta0, kPdgP31m1910_DeltaP, kPdgP31m1910_DeltaPP,
    kPdgP33m1920_DeltaM, kPdgP33m1920_Delta0, kPdgP33m1920_DeltaP, kPdgP33m1920_DeltaPP,
    kPdgF37m1950_DeltaM, kPdgF37m1950_Delta0, kPdgF37m1950_DeltaP, kPdgF37m1950_DeltaPP,
    kPdgP11m1440_N0, kPdgP11m1440_NP,    
    kPdgD13m1520_N0, kPdgD13m1520_NP,   
    kPdgS11m1535_N0, kPdgS11m1535_NP,  
    kPdgS11m1650_N0, kPdgS11m1650_NP,  
    kPdgD15m1675_N0, kPdgD15m1675_NP,    
    kPdgF15m1680_N0, kPdgF15m1680_NP,     
    kPdgD13m1700_N0, kPdgD13m1700_NP,  
    kPdgP11m1710_N0, kPdgP11m1710_NP,      
    kPdgP13m1720_N0, kPdgP13m1720_NP,
    0
  };

  // init
  const int    kNW   = 30;
  const double kWmin = 0.9;
  const double kWmax = 5.0;
  unsigned int ires=0;
  int respdg = 0;
  while((respdg = respdgarray[ires++])) {
    fMpBR1gammaDef [respdg] = new TH1D("","",kNW,kWmin,kWmax);
    fMpBR1etaDef   [respdg] = new TH1D("","",kNW,kWmin,kWmax);
    fMpBR1gammaDef [respdg] -> SetDirectory(0);
    fMpBR1etaDef   [respdg] -> SetDirectory(0);
  }

  // find corresponding decay channels and store default BR
  PDGLibrary * pdglib = PDGLibrary::Instance();
  ires=0;
  while((respdg = respdgarray[ires++])) {
    TParticlePDG * res = pdglib->Find(respdg);
    if(!res) continue;
    for(int j=0; j<res->NDecayChannels(); j++) {
        TDecayChannel * dch = res->DecayChannel(j);
        int ngamma = 0;
        int neta   = 0;
        double br = dch->BranchingRatio();
        double mt = 0.;
        for(int k=0; k<dch->NDaughters(); k++) {
          int dpdg = dch->DaughterPdgCode(k);
          if(dpdg==kPdgGamma) ngamma++; 
          if(dpdg==kPdgEta  ) neta++; 
          mt += 0;
        }//decay channel f/s particles
        bool is_1gamma = (ngamma==1);
        bool is_1eta   = (neta  ==1);
        for(int ibin = 1; ibin <= fMpBR1gammaDef[respdg]->GetNbinsX(); ibin++) {
          double W = fMpBR1gammaDef[respdg]->GetBinLowEdge(ibin) + 
                     fMpBR1gammaDef[respdg]->GetBinWidth(ibin);
          bool is_allowed = (W>mt);
          if(is_allowed && is_1gamma) { fMpBR1gammaDef[respdg]->Fill(W, br); }
          if(is_allowed && is_1eta  ) { fMpBR1etaDef  [respdg]->Fill(W, br); }
        }//W bins
    }//decay channels
  }//resonances

#ifdef _G_REWEIGHT_RESDEC_DEBUG_
  fTestFile = new TFile("./resdec_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","costheta:wght");
#endif
}
//_______________________________________________________________________________________
