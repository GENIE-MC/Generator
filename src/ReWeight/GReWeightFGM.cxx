//____________________________________________________________________________
/*
 Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Apr 27, 2010 - CA
   First included in v2.7.1.

*/
//____________________________________________________________________________

#include <TFile.h>
#include <TNtupleD.h>

#include "Algorithm/AlgFactory.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "Nuclear/FermiMomentumTablePool.h"
#include "Nuclear/FermiMomentumTable.h"
#include "Nuclear/NuclearModelI.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightFGM.h"
#include "ReWeight/GSystUncertainty.h"
#include "Utils/NuclearUtils.h"

//#define _G_REWEIGHT_FGM_DEBUG_

using namespace genie;
using namespace genie::rew;
using namespace genie::utils;

//_______________________________________________________________________________________
GReWeightFGM::GReWeightFGM() :
GReWeightI()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightFGM::~GReWeightFGM()
{
#ifdef _G_REWEIGHT_FGM_DEBUG_
  fTestFile->cd();
  fTestNtp ->Write();
  fTestFile->Close();
  delete fTestFile;
#endif
}
//_______________________________________________________________________________________
bool GReWeightFGM::IsHandled(GSyst_t syst)
{
  switch(syst) {
    case ( kSystNucl_CCQEPauliSupViaKF   ) : 
    case ( kSystNucl_CCQEMomDistroFGtoSF ) : 
        return true;
        break;
    default:
        return false;
        break;
  }
  return false;
}
//_______________________________________________________________________________________
void GReWeightFGM::SetSystematic(GSyst_t syst, double val)
{
  switch(syst) {
    case ( kSystNucl_CCQEPauliSupViaKF   ) : 
        fKFTwkDial = val;
        break;
    case ( kSystNucl_CCQEMomDistroFGtoSF ) : 
        fMomDistroTwkDial = val;
        break;
    default:
        return;
        break;
  }
}
//_______________________________________________________________________________________
void GReWeightFGM::Reset(void)
{
  fKFTwkDial        = 0.;
  fMomDistroTwkDial = 0.;
}
//_______________________________________________________________________________________
void GReWeightFGM::Reconfigure(void)
{
  
}
//_______________________________________________________________________________________
double GReWeightFGM::CalcWeight(const EventRecord & event) 
{ 
  double wght =
    this->RewCCQEPauliSupViaKF   (event) *
    this->RewCCQEMomDistroFGtoSF (event);
 
  return wght;
}
//_______________________________________________________________________________________
double GReWeightFGM::CalcChisq(void)
{
  return 0.;
}
//_______________________________________________________________________________________
double GReWeightFGM::RewCCQEPauliSupViaKF(const EventRecord & event) 
{
  bool kF_tweaked = (TMath::Abs(fKFTwkDial) > controls::kASmallNum);
  if(!kF_tweaked) return 1.;

  bool is_qe = event.Summary()->ProcInfo().IsQuasiElastic();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if(!is_qe || !is_cc) return 1.;

  const Target & target = event.Summary()->InitState().Tgt();
  if (!target.IsNucleus()) {
     return 1.; 
  }

  double A = target.A();
  if(A <= 4) {
     return 1.; 
  }

  int target_pdgc         = target.Pdg();
  int struck_nucleon_pdgc = target.HitNucPdg();
  int final_nucleon_pdgc  = pdg::SwitchProtonNeutron(struck_nucleon_pdgc);
             
  FermiMomentumTablePool * kftp = FermiMomentumTablePool::Instance();
  const FermiMomentumTable * kft  = kftp->GetTable("Default");

  // Fermi momentum for initial, final nucleons
  double kFi = kft->FindClosestKF(target_pdgc, struck_nucleon_pdgc);
  double kFf = kft->FindClosestKF(target_pdgc, final_nucleon_pdgc );
 
  // hit nucleon mass
  double Mn = target.HitNucP4Ptr()->M(); // can be off m/shell

  // momentum transfer
  const Kinematics & kine = event.Summary()->Kine();
  bool selected = true;
  double q2 = kine.q2(selected);
  
  const double pmax = 0.5; // GeV

  // default nuclear suppression 
  double R = utils::nuclear::RQEFG_generic(q2, Mn, kFi, kFf, pmax);

  // tweak kF
  GSystUncertainty * uncertainty = GSystUncertainty::Instance();
  double kF_fracerr = uncertainty->OneSigmaErr(kSystNucl_CCQEPauliSupViaKF);

  double kFi_twk = kFi * (1 + fKFTwkDial * kF_fracerr);
  double kFf_twk = kFf * (1 + fKFTwkDial * kF_fracerr);

  // calculate tweaked nuclear suppression factor
  double Rtwk = nuclear::RQEFG_generic(q2, Mn, kFi_twk, kFf_twk, pmax);
  
  // calculate weight (ratio of suppression factors)

  double wght = 1.0;

  if(R>0 && Rtwk>0) {
    wght = Rtwk/R;
  }

#ifdef _G_REWEIGHT_FGM_DEBUG_
  fTestNtp->Fill(-q2,wght);
#endif

  return wght;
}
//_______________________________________________________________________________________
double GReWeightFGM::RewCCQEMomDistroFGtoSF(const EventRecord & event) 
{
  bool momdistro_tweaked = (TMath::Abs(fMomDistroTwkDial) > controls::kASmallNum);
  if(!momdistro_tweaked) return 1.;

  bool is_qe = event.Summary()->ProcInfo().IsQuasiElastic();
  bool is_cc = event.Summary()->ProcInfo().IsWeakCC();
  if(!is_qe || !is_cc) return 1.;

  GHepParticle * tgt = event.TargetNucleus();
  if(!tgt) return 1.; // scattering off free-nucleon 

  GHepParticle * hitnuc = event.HitNucleon();
  if(!hitnuc) return 1.;

  const double kPmax = 0.5;
  double p = hitnuc->P4()->Vect().Mag();
  if(p > kPmax) return 1.;

  TH1D * hfg = 0;
  TH1D * hsf = 0;

  int tgtpdg = tgt    -> Pdg();
  int nucpdg = hitnuc -> Pdg();

  map<int, TH1D*> & mapfg = pdg::IsNeutron(nucpdg) ? fMapFGn : fMapFGp;
  map<int, TH1D*> & mapsf = pdg::IsNeutron(nucpdg) ? fMapSFn : fMapSFp;

  map<int, TH1D*>::iterator it;
  it = mapfg.find(tgtpdg);
  if(it != mapfg.end()) { hfg = it->second; }
  it = mapsf.find(tgtpdg);
  if(it != mapsf.end()) { hsf = it->second; }

  bool have_weight_func = (hfg!=0) && (hsf!=0);
  if(!have_weight_func) {
     const int kNEv  = 20000;
     const int kNP   = 500;
     hfg = new TH1D("","",kNP,0.,kPmax);
     hsf = new TH1D("","",kNP,0.,kPmax);
     hfg -> SetDirectory(0);
     hsf -> SetDirectory(0);
     const Target & tgt = event.Summary()->InitState().Tgt();
     bool ok = true;
     for(int iev=0; iev<kNEv; iev++) {
       ok = fFG->GenerateNucleon(tgt);
       if(!ok) {
         delete hfg;
         delete hsf;
         return 1.;
       }
       hfg->Fill(fFG->Momentum());
     }//fg
     for(int iev=0; iev<kNEv; iev++) {
       ok = fSF->GenerateNucleon(tgt);
       if(!ok) {
         delete hfg;
         delete hsf;
         return 1.;
       }
       hsf->Fill(fSF->Momentum());
     }//sf
     hfg->Scale(1. / hfg->Integral("width"));
     hsf->Scale(1. / hsf->Integral("width"));
     mapfg.insert(map<int,TH1D*>::value_type(tgtpdg,hfg));
     mapsf.insert(map<int,TH1D*>::value_type(tgtpdg,hsf));
  }//create & store momentum distributions


  double f_fg = hfg->GetBinContent( hfg->FindBin(p) );
  double f_sf = hsf->GetBinContent( hsf->FindBin(p) );
  double dial = fMomDistroTwkDial;
  double wght = (f_sf * dial + f_fg * (1-dial)) / f_fg;

  return wght;
}
//_______________________________________________________________________________________
void GReWeightFGM::Init(void)
{
  fKFTwkDial        = 0.;
  fMomDistroTwkDial = 0.;

  AlgFactory * algf = AlgFactory::Instance();

  fFG = dynamic_cast<const NuclearModelI*> (
    algf->GetAlgorithm("genie::FGMBodekRitchie","Default"));
  fSF = dynamic_cast<const NuclearModelI*> (
    algf->GetAlgorithm("genie::SpectralFunc","Default"));

#ifdef _G_REWEIGHT_FGM_DEBUG_
  fTestFile = new TFile("./fgm_reweight_test.root","recreate");
  fTestNtp  = new TNtupleD("testntp","","Q2:wght");
#endif
}
//_______________________________________________________________________________________

