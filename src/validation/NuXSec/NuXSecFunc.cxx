//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ 01/05/2012 - CA
   Added in v2.7.1
*/
//____________________________________________________________________________

#include <sstream>

#include <TH1D.h>
#include <TMath.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "PDG/PDGLibrary.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightI.h"
#include "ReWeight/GSystSet.h"
#include "ReWeight/GReWeightNuXSecNCEL.h"
#include "ReWeight/GReWeightNuXSecCCQE.h"
#include "ReWeight/GReWeightNuXSecCCRES.h"
#include "ReWeight/GReWeightNuXSecCOH.h"
#include "ReWeight/GReWeightNonResonanceBkg.h"
#include "ReWeight/GReWeightFGM.h"
#include "ReWeight/GReWeightDISNuclMod.h"
#include "ReWeight/GReWeightResonanceDecay.h"
#include "ReWeight/GReWeightFZone.h"
#include "ReWeight/GReWeightINuke.h"
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GReWeightNuXSecCCQEvec.h"
#include "ReWeight/GReWeightNuXSecNCRES.h"  
#include "ReWeight/GReWeightNuXSecDIS.h"    
#include "validation/NuXSec/NuXSecFunc.h"

using std::ostringstream;

using namespace genie;
using namespace genie::rew;
using namespace genie::mc_vs_data;

//____________________________________________________________________________
string NuXSecFunc::BuildXSecDirectoryName(int nupdg, int tgtpdg)
{
// build the name of TDirectory of the input cross-section file 
// that corresponds to the specified initial state
  ostringstream xsec_dir_name;
  PDGLibrary * pdglib = PDGLibrary::Instance();
  string probe_name = pdglib->Find(nupdg)->GetName();
  string tgt_name = (tgtpdg==1000000010) ? "n" : pdglib->Find(tgtpdg)->GetName();
  xsec_dir_name << probe_name << "_" << tgt_name;
  return xsec_dir_name.str();
}
//____________________________________________________________________________
XSecForModeX::XSecForModeX() :
NuXSecFunc()
{
  fNuisanceParams.clear();
  fXSecDirectoryName = "";
  fModeXSplineName   = "";
  fXSecScaleFactor   = 1.;
}
//.............................................................................
XSecForModeX::~XSecForModeX()
{

}
//.............................................................................
TGraphAsymmErrors * XSecForModeX::ExtractFromEventSample(
     int imodel, double Emin, double Emax,
     int n, bool inlogE, bool scale_with_E, bool incl_err_band)
{
  if(!fGenieInputs) {
     LOG("gvldtest", pERROR) << "No GENIE MC inputs";
     return 0;
  }
  TFile * genie_xsec_file = fGenieInputs->XSecFile(imodel);
  if(!genie_xsec_file) {
     LOG("gvldtest", pERROR) << "No input GENIE cross-section file";
     return 0;
  }
  TChain * genie_event_tree = fGenieInputs->EvtChain(imodel);
  if(!genie_event_tree) {
     LOG("gvldtest", pERROR) << "No input GENIE event tree";
     return 0;
  }

  // Get xsec directory and retrieve inclusive CC xsec
  TDirectory * xsec_dir = 
      (TDirectory *) genie_xsec_file->Get(fXSecDirectoryName.c_str());
  if(!xsec_dir) {
     LOG("gvldtest", pERROR) 
        << "Can't find cross-section directory: " << fXSecDirectoryName;
     return 0;
  }
  TGraph * incl_xsec_spline = (TGraph*) xsec_dir->Get("tot_cc");
  if(!incl_xsec_spline) {
     LOG("gvldtest", pERROR) 
        << "Can't find results of inclusive CC cross-section calculation";
     return 0;
  }
  // Check whether the cross-section for mode X is pre-computed and is
  // available in the input cross-section file. This is possible for
  // some of the simpler moder such as QE or coherent, while for most
  // exclusive inelastic channels it can only be calculated via MC.
  TGraph * modex_xsec_spline = 0;
  if(fModeXSplineName.size()>0) {
    modex_xsec_spline = (TGraph*) xsec_dir->Get(fModeXSplineName.c_str());
    if(!modex_xsec_spline) {
       LOG("gvldtest", pERROR) 
	 << "Can't find cross-section spline: " << fModeXSplineName;
       return 0;
    }
  }

  // Set event tree branch address
  genie_event_tree->SetBranchStatus("gmcrec", 1);   
  Long64_t nmax = genie_event_tree->GetEntries();
  if (nmax<0) {
     LOG("gvldtest", pERROR) << "Number of events = 0";
     return 0;
  }
  LOG("gvldtest", pNOTICE) 
    << "Found " << nmax << " entries in the event tree";
  NtpMCEventRecord * mcrec = 0;
  genie_event_tree->SetBranchAddress("gmcrec", &mcrec);   
  if (!mcrec) {
     LOG("gvldtest", pERROR) << "Null MC record";
     return 0;
  }

  // Book histograms for nominal and systematically varied event fractions
  const unsigned int knsyst_max = 10;
  const unsigned int kntwk = 2;
  const double twk[kntwk] = { +1., -1. };
  double xmin = (inlogE) ? TMath::Log10(Emin) : Emin;
  double xmax = (inlogE) ? TMath::Log10(Emax) : Emax;
  TH1D * hmx   = new TH1D("hmx",   "", n, xmin, xmax);
  TH1D * hincl = new TH1D("hincl", "", n, xmin, xmax);
  TH1D * hmxtwk[knsyst_max][kntwk];
  for(unsigned int isyst = 0; isyst < knsyst_max; isyst++) {    
    for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
       hmxtwk[isyst][itwk] = 0;
       if(isyst < fNuisanceParams.size()) {
          hmxtwk[isyst][itwk] = new TH1D( Form("hmx_%d_%d",isyst,itwk), "", n, xmin, xmax);
       }       
    }
  }

  // Fill in nominal event fractions
  // If I am going to show pre-tabulated cross-section calculations and not 
  // cross-section computed from the input sample, then this is only needed
  // if I am going to show error envelopes.
  if(!modex_xsec_spline || incl_err_band) {
    int curr_tree_num = -1;
    for(Long64_t iev = 0; iev < nmax; iev++) {
      genie_event_tree->GetEntry(iev); 
      EventRecord & event = *(mcrec->event);
      // Generated event files used as inputs in the validation programs
      // correspond to given neutrino+target initial states.
      // If the current file in the chain doesn't correspond to a desired 
      // initial state (as determined by looking at the first event in the file)
      // then skip this file althogether to save processing time.
      int tree_num = genie_event_tree->GetTreeNumber();
      if(curr_tree_num != tree_num) {
         curr_tree_num = tree_num;
         bool skip = !this->UseTree(event);
         if(skip) {
           LOG("gvldtest", pNOTICE) << "Skipping tree # : " << tree_num;
           iev += (genie_event_tree->GetTree()->GetEntries() - 1);
           continue;
         }//skip?
      }//new tree in chain
      GHepParticle * neutrino = event.Probe();
      assert(neutrino);
      double E = neutrino->P4()->E();
      double x = (inlogE) ? TMath::Log10(E) : E;
      if(this->IsModeX (event)) { hmx  ->Fill(x); }
      if(this->IsCC    (event)) { hincl->Fill(x); }
    }
    hmx->Divide(hincl);
    hmx->Smooth(2);
  }

  // Include uncertainties if requested
  if(incl_err_band) {
     GSystSet & systlist = fRew.Systematics();
     for(unsigned int isyst = 0; isyst < fNuisanceParams.size(); isyst++) {
       systlist.Init(fNuisanceParams[isyst]);
       for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
          systlist.Set(fNuisanceParams[isyst], twk[itwk]);
          fRew.Reconfigure();
          int curr_tree_num = -1;
          for(Long64_t iev = 0; iev < nmax; iev++) {
              genie_event_tree->GetEntry(iev);
              EventRecord & event = *(mcrec->event);
              int tree_num = genie_event_tree->GetTreeNumber();
              if(curr_tree_num != tree_num) {
                curr_tree_num = tree_num;
                 bool skip = !this->UseTree(event);
                 if(skip) {           
                   LOG("gvldtest", pNOTICE) << "Skipping tree # : " << tree_num;
                   iev += (genie_event_tree->GetTree()->GetEntries() - 1);
                   continue;
                 }//skip?
              }//new tree in chain
              if(this->IsModeX(event)) {
                 GHepParticle * neutrino = event.Probe();
                 assert(neutrino);
                 double E = neutrino->P4()->E();
                 double x = (inlogE) ? TMath::Log10(E) : E;
                 double wght = fRew.CalcWeight(event);
                 hmxtwk[isyst][itwk]->Fill(x, wght);
              }//X?
          }//iev
          hmxtwk[isyst][itwk]->Divide(hincl);
          hmxtwk[isyst][itwk]->Smooth(2);
       }//itwk
       systlist.Set(fNuisanceParams[isyst], 0.);//reset
       fRew.Reconfigure();
       systlist.Remove(fNuisanceParams[isyst]);
     }//isyst
  }//incl_err_band?


  // Calculate cross-section = f(E) and corresponding uncertainty.
  // Sources of uncertainty are taken to be uncorellated.
  double * energy_array    = new double [n];
  double * xsec_array      = new double [n];
  double * xsec_array_errp = new double [n];
  double * xsec_array_errm = new double [n];
  for(int i = 0; i < n; i++) {
    int ibin = i+1;
    double energy = (inlogE) ? 
       TMath::Power(10., hmx->GetBinCenter(ibin)) : hmx->GetBinCenter(ibin);
    assert(energy>0.);
    double incl_xsec = incl_xsec_spline->Eval(energy);
    double evt_frac  = hmx->GetBinContent(ibin);
    double xsec      = evt_frac * incl_xsec;
    double xsec_errp = 0;
    double xsec_errm = 0;
    if(incl_err_band) {
       double xsec_errp_s[knsyst_max];
       double xsec_errm_s[knsyst_max];
       for(unsigned int isyst = 0; isyst < fNuisanceParams.size(); isyst++) {
          xsec_errp_s[isyst] = 0;
          xsec_errm_s[isyst] = 0;
          double max_evt_frac = -999999;
          double min_evt_frac =  999999;
          for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
            max_evt_frac = TMath::Max(max_evt_frac, hmxtwk[isyst][itwk]->GetBinContent(ibin));
            min_evt_frac = TMath::Min(min_evt_frac, hmxtwk[isyst][itwk]->GetBinContent(ibin));
          }
          double xsec_hi = max_evt_frac * incl_xsec;
          double xsec_lo = min_evt_frac * incl_xsec;
          xsec_errp_s[isyst] = TMath::Max(0., xsec_hi-xsec);
          xsec_errm_s[isyst] = TMath::Max(0., xsec-xsec_lo);
          xsec_errp += TMath::Power(xsec_errp_s[isyst],2.);
          xsec_errm += TMath::Power(xsec_errm_s[isyst],2.);
       }//isyst
       xsec_errp = TMath::Sqrt(xsec_errp);
       xsec_errm = TMath::Sqrt(xsec_errm);       
    }//incl_err_band?
    // 
    if(modex_xsec_spline) {
       xsec = modex_xsec_spline->Eval(energy);
    }
    // apply scaling factor, typically to normalize nuclear cross-section
    // to number of nucleons as expt cross-sections are typically quoted as such.
    xsec      *= fXSecScaleFactor;        
    xsec_errp *= fXSecScaleFactor;        
    xsec_errm *= fXSecScaleFactor;        
    LOG("gvldtest", pNOTICE)  
        << "xsec (E = " << energy << " GeV) = " << xsec << " 1E-38 cm^2 "
        << "+ " << ((xsec>0) ? 100*xsec_errp/xsec : 0) << "% "
        << "- " << ((xsec>0) ? 100*xsec_errm/xsec : 0) << "% ";    
    energy_array[i]    = energy;
    xsec_array[i]      = (scale_with_E) ? xsec/energy       : xsec;
    xsec_array_errp[i] = (scale_with_E) ? xsec_errp/energy  : xsec_errp;
    xsec_array_errm[i] = (scale_with_E) ? xsec_errm/energy  : xsec_errm;
  }

  // Build cross-section graph
  TGraphAsymmErrors * model = new TGraphAsymmErrors(
     n,energy_array,xsec_array,0,0,xsec_array_errm, xsec_array_errp);
  assert(model);
 
  // Clean-up
  delete hmx;
  delete hincl;
  for(unsigned int isyst = 0; isyst < knsyst_max; isyst++) {
    for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
       if(hmxtwk[isyst][itwk]) { delete hmxtwk[isyst][itwk]; }
    }
  }
  delete [] energy_array;
  delete [] xsec_array;
  delete [] xsec_array_errp;
  delete [] xsec_array_errm;

  LOG("gvldtest", pNOTICE) << "Returning model prediction";

  return model;
}
//____________________________________________________________________________
CCQEXSec::CCQEXSec(int nupdg, int tgtpdg, int hitnucpdg, double xsec_scale) :
XSecForModeX(),
fNuPdg(nupdg),
fTgtPdg(tgtpdg),
fHitNucPdg(hitnucpdg)
{
  fXSecScaleFactor   = xsec_scale; 
  fXSecDirectoryName = this->BuildXSecDirectoryName(nupdg, tgtpdg);

  fModeXSplineName = "";
  if (hitnucpdg == kPdgNeutron) { fModeXSplineName = "qel_cc_n"; }
  if (hitnucpdg == kPdgProton ) { fModeXSplineName = "qel_cc_p"; }

  // List of nuisance params considered for CCQE cross-section
  fNuisanceParams.push_back(kXSecTwkDial_MaCCQE);
  if(pdg::IonPdgCodeToA(tgtpdg) > 2) {
     fNuisanceParams.push_back(kSystNucl_CCQEPauliSupViaKF);
  }
  // Add weight calculation engines
  fRew.AdoptWghtCalc( "xsec_ccqe",  new GReWeightNuXSecCCQE );
  fRew.AdoptWghtCalc( "nuclear_qe", new GReWeightFGM        );
  // Fine-tuning reweighting engines
  GReWeightNuXSecCCQE * rwccqe =      
    dynamic_cast<GReWeightNuXSecCCQE *> (fRew.WghtCalc("xsec_ccqe"));  
  rwccqe->SetMode(GReWeightNuXSecCCQE::kModeMa);

}
//.............................................................................
CCQEXSec::~CCQEXSec()
{

}
//.............................................................................
string CCQEXSec::Name(void) const 
{ 
  return Form("CCQE;nu=%d;tgt=%d;hitnuc=%d;scale=%.4f",
               fNuPdg,fTgtPdg,fHitNucPdg,fXSecScaleFactor);
}
//.............................................................................
bool CCQEXSec::IsCC(EventRecord & event)
{
  Interaction * in = event.Summary();
  if(!in->ProcInfo().IsWeakCC()) return false;
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;
  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  // For free-nucleon targets, the usual nucleon pdg codes are used,
  // not the ion codes 1000010010 and 1000000010.
  // In this case, compare with the specified hit nucleon code.
  if(tgtpdg == kPdgNeutron || tgtpdg == kPdgProton) {
    if(tgtpdg != fHitNucPdg) return false;
  }
  else {
    if(tgtpdg != fTgtPdg) return false;
  }
  return true;
}
//.............................................................................
bool CCQEXSec::IsModeX(EventRecord & event)
{
  Interaction * in = event.Summary();
  if(!in->ProcInfo().IsWeakCC()) return false;
  if(!in->ProcInfo().IsQuasiElastic()) return false;
  if(in->ExclTag().IsCharmEvent()) return false; // skip charm QE
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;
  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  // For free-nucleon targets, the usual nucleon pdg codes are used,
  // not the ion codes 1000010010 and 1000000010.
  // In this case, compare with the specified hit nucleon code.
  if(tgtpdg == kPdgNeutron || tgtpdg == kPdgProton) {
    if(tgtpdg != fHitNucPdg) return false;
  }
  else {
    if(tgtpdg != fTgtPdg) return false;
  }
  return true;
}
//.............................................................................
bool CCQEXSec::UseTree(EventRecord & event)
{
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;

  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  if(pdg::IonPdgCodeToA(fTgtPdg) > 1) {
    if(tgtpdg == fTgtPdg) return true;
  }
  else
  if(pdg::IonPdgCodeToA(fTgtPdg) == 1) {
    bool match =
       (tgtpdg == kPdgProton  && fTgtPdg == kPdgTgtFreeP) ||
       (tgtpdg == kPdgNeutron && fTgtPdg == kPdgTgtFreeN);
    if(match) return true;
  }
  return false;
}
//____________________________________________________________________________
CCPionXSec::CCPionXSec(
 int nupdg, int tgtpdg, int hitnucpdg, 
 int npip, int npi0, int npim, int np, int nn, double xsec_scale) :
XSecForModeX(),
fNuPdg(nupdg),
fTgtPdg(tgtpdg),
fHitNucPdg(hitnucpdg),
fNpip(npip),
fNpi0(npi0),
fNpim(npim),
fNp(np),
fNn(nn)
{
  fXSecScaleFactor   = xsec_scale; 
  fXSecDirectoryName = this->BuildXSecDirectoryName(nupdg, tgtpdg);

  // List of nuisance params considered for CCQE cross-section

  fNuisanceParams.push_back(kXSecTwkDial_MaCCRES);
//fNuisanceParams.push_back(kXSecTwkDial_AhtBY);
//fNuisanceParams.push_back(kXSecTwkDial_BhtBY);
//fNuisanceParams.push_back(kXSecTwkDial_CV1uBY);
//fNuisanceParams.push_back(kXSecTwkDial_CV2uBY);

  int npi = npip + npi0 + npim;
  if( pdg::IsNeutrino(nupdg) && pdg::IsProton(hitnucpdg) && npi==1) {
     fNuisanceParams.push_back(kXSecTwkDial_RvpCC1pi);
  }
  if( pdg::IsNeutrino(nupdg) && pdg::IsProton(hitnucpdg) && npi==2) {
     fNuisanceParams.push_back(kXSecTwkDial_RvpCC2pi);
  }
  if( pdg::IsNeutrino(nupdg) && pdg::IsNeutron(hitnucpdg) && npi==1) {
     fNuisanceParams.push_back(kXSecTwkDial_RvnCC1pi);
  }
  if( pdg::IsNeutrino(nupdg) && pdg::IsNeutron(hitnucpdg) && npi==2) {
     fNuisanceParams.push_back(kXSecTwkDial_RvnCC2pi);
  }
  if( pdg::IsAntiNeutrino(nupdg) && pdg::IsProton(hitnucpdg) && npi==1) {
     fNuisanceParams.push_back(kXSecTwkDial_RvbarpCC1pi);
  }
  if( pdg::IsAntiNeutrino(nupdg) && pdg::IsProton(hitnucpdg) && npi==2) {
     fNuisanceParams.push_back(kXSecTwkDial_RvbarpCC2pi);
  }
  if( pdg::IsAntiNeutrino(nupdg) && pdg::IsNeutron(hitnucpdg) && npi==1) {
     fNuisanceParams.push_back(kXSecTwkDial_RvbarnCC1pi);
  }
  if( pdg::IsAntiNeutrino(nupdg) && pdg::IsNeutron(hitnucpdg) && npi==2) {
     fNuisanceParams.push_back(kXSecTwkDial_RvbarnCC2pi);
  }

  // Add weight calculation engines
  fRew.AdoptWghtCalc( "xsec_ccres",      new GReWeightNuXSecCCRES     );
  fRew.AdoptWghtCalc( "xsec_nonresbkg",  new GReWeightNonResonanceBkg );
//fRew.AdoptWghtCalc( "xsec_dis",        new GReWeightNuXSecDIS       );

  // Fine-tune weight calculation engines
  GReWeightNuXSecCCRES * rwccres = 
      dynamic_cast<GReWeightNuXSecCCRES *> (fRew.WghtCalc("xsec_ccres"));  
  rwccres->SetMode(GReWeightNuXSecCCRES::kModeMaMv);
}
//.............................................................................
CCPionXSec::~CCPionXSec()
{

}
//.............................................................................
string CCPionXSec::Name(void) const 
{ 
  return Form(
   "CCpi;nu=%d;tgt=%d;hitnuc=%d;fspi=%d,%d,%d;fsnuc:%d,%d;scale=%.4f",
    fNuPdg,fTgtPdg,fHitNucPdg,fNpip,fNpi0,fNpim,fNp,fNn,fXSecScaleFactor);
}
//.............................................................................
bool CCPionXSec::IsCC(EventRecord & event)
{
  Interaction * in = event.Summary();
  if(!in->ProcInfo().IsWeakCC()) return false;
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;
  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  if(tgtpdg == kPdgNeutron || tgtpdg == kPdgProton) {
    if(tgtpdg != fHitNucPdg) return false;
  }
  else {
    if(tgtpdg != fTgtPdg) return false;
  }
  return true;
}
//.............................................................................
bool CCPionXSec::IsModeX(EventRecord & event)
{
  Interaction * in = event.Summary();
  if(!in->ProcInfo().IsWeakCC()) return false;
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;
  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  if(tgtpdg == kPdgNeutron || tgtpdg == kPdgProton) {
    if(tgtpdg != fHitNucPdg) return false;
  }
  else {
    if(tgtpdg != fTgtPdg) return false;
    GHepParticle * hitnuc = event.Particle(2); 
    if(hitnuc->Pdg() != fHitNucPdg) return false;
  }
  int npip = 0;
  int npi0 = 0;
  int npim = 0;
  int np   = 0;
  int nn   = 0;
  TObjArrayIter piter(&event);
  GHepParticle * p = 0;
  while ((p = (GHepParticle *) piter.Next())) {
    if(p->Status() != kIStStableFinalState) continue;
    if(p->Pdg() == kPdgPiP    ) npip++;
    if(p->Pdg() == kPdgPi0    ) npi0++;
    if(p->Pdg() == kPdgPiM    ) npim++;
    if(p->Pdg() == kPdgProton ) np++;
    if(p->Pdg() == kPdgNeutron) nn++;
  }
  bool accept = (npip==fNpip && npi0==fNpi0 && npim==fNpim && np==fNp && nn==fNn);
  return accept;
}
//.............................................................................
bool CCPionXSec::UseTree(EventRecord & event)
{
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;

  GHepParticle * target = event.Particle(1);
  int tgtpdg = target->Pdg();
  if(pdg::IonPdgCodeToA(fTgtPdg) > 1) {
    if(tgtpdg == fTgtPdg) return true; 
  }
  else
  if(pdg::IonPdgCodeToA(fTgtPdg) == 1) {
    bool match =
       (tgtpdg == kPdgProton  && fTgtPdg == kPdgTgtFreeP) ||
       (tgtpdg == kPdgNeutron && fTgtPdg == kPdgTgtFreeN);
    if(match) return true;
  }
  return false;
}
//____________________________________________________________________________
CohPionXSec::CohPionXSec(int nupdg, int tgtpdg, int pipdg, double xsec_scale) :
XSecForModeX(),
fNuPdg(nupdg),
fTgtPdg(tgtpdg),
fPiPdg(pipdg)
{
  fXSecScaleFactor   = xsec_scale; 
  fXSecDirectoryName = this->BuildXSecDirectoryName(nupdg, tgtpdg);

  fModeXSplineName = "";
  if (pipdg == kPdgPi0) { 
     fModeXSplineName = "coh_nc"; 
  } else
  if (pipdg == kPdgPiP || pipdg == kPdgPiM) { 
     fModeXSplineName = "coh_cc"; 
  }

  // List of nuisance params considered for coherent pion cross-section
// fNuisanceParams.push_back(kXSecTwkDial_MaCOHpi);
// fNuisanceParams.push_back(kXSecTwkDial_R0COHpi);
  // Add weight calculation engines
// fRew.AdoptWghtCalc( "xsec_coh",  new GReWeightNuXSecCOH);

}
//.............................................................................
CohPionXSec::~CohPionXSec()
{

}
//.............................................................................
string CohPionXSec::Name(void) const 
{
  return Form("COHpi;nu=%d;tgt=%d;pi=%d;scale=%.4f",
               fNuPdg,fTgtPdg,fPiPdg,fXSecScaleFactor);
}
//.............................................................................
bool CohPionXSec::IsCC(EventRecord & event)
{
  Interaction * in = event.Summary();
  if(!in->ProcInfo().IsWeakCC()) return false;
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;
  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  if(tgtpdg != fTgtPdg) return false;

  return true;
}
//.............................................................................
bool CohPionXSec::IsModeX(EventRecord & event)
{
  Interaction * in = event.Summary();
  if(!in->ProcInfo().IsCoherent()) return false;
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;
  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  if(tgtpdg != fTgtPdg) return false;
  bool allowed_mode =
    (pdg::IsNeutrino(fNuPdg)     && in->ProcInfo().IsWeakCC() && fPiPdg==kPdgPiP) ||
    (pdg::IsNeutrino(fNuPdg)     && in->ProcInfo().IsWeakNC() && fPiPdg==kPdgPi0) ||
    (pdg::IsAntiNeutrino(fNuPdg) && in->ProcInfo().IsWeakCC() && fPiPdg==kPdgPiM) ||
    (pdg::IsAntiNeutrino(fNuPdg) && in->ProcInfo().IsWeakNC() && fPiPdg==kPdgPi0); 
  if(!allowed_mode) return false;

  return true;
}
//.............................................................................
bool CohPionXSec::UseTree(EventRecord & event)
{
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;

  GHepParticle * target = event.Particle(1);
  int tgtpdg = target->Pdg();
  if(tgtpdg == kPdgNeutron || tgtpdg == kPdgProton) return false;
  
  if(tgtpdg == fTgtPdg) return true; 

  return false;
}
//____________________________________________________________________________
XSecRatioForModesXY::XSecRatioForModesXY() :
NuXSecFunc()
{
  fNuisanceParams.clear();
}
//____________________________________________________________________________
XSecRatioForModesXY::~XSecRatioForModesXY()
{

}
//.............................................................................
TGraphAsymmErrors * XSecRatioForModesXY::ExtractFromEventSample(
     int imodel, double Emin, double Emax,
     int n, bool inlogE, bool scale_with_E, bool incl_err_band)
{
  if(!fGenieInputs) return 0;
  TChain * genie_event_tree = fGenieInputs->EvtChain(imodel);
  if(!genie_event_tree) return 0;

  // Set event tree branch address
  genie_event_tree->SetBranchStatus("gmcrec", 1);   
  Long64_t nmax = genie_event_tree->GetEntries();
  if (nmax<0) {
     LOG("gvldtest", pERROR) << "Number of events = 0";
     return 0;
  }
  LOG("gvldtest", pNOTICE) 
    << "Found " << nmax << " entries in the event tree";
  NtpMCEventRecord * mcrec = 0;
  genie_event_tree->SetBranchAddress("gmcrec", &mcrec);   
  if (!mcrec) {
     LOG("gvldtest", pERROR) << "Null MC record";
     return 0;
  }

  // Book histograms for nominal and systematically varied event fractions
  const unsigned int knsyst_max = 10;
  const unsigned int kntwk = 2;
  const double twk[kntwk] = { +1., -1. };
  double xmin = (inlogE) ? TMath::Log10(Emin) : Emin;
  double xmax = (inlogE) ? TMath::Log10(Emax) : Emax;
  TH1D * hmx   = new TH1D("hmx",   "", n, xmin, xmax);
  TH1D * hmy   = new TH1D("hmy",   "", n, xmin, xmax);
  TH1D * hmxtwk[knsyst_max][kntwk];
  TH1D * hmytwk[knsyst_max][kntwk];
  for(unsigned int isyst = 0; isyst < knsyst_max; isyst++) {    
    for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
       hmxtwk[isyst][itwk] = 0;
       hmytwk[isyst][itwk] = 0;
       if(isyst < fNuisanceParams.size()) {
          hmxtwk[isyst][itwk] = new TH1D( Form("hmx_%d_%d",isyst,itwk), "", n, xmin, xmax);
          hmytwk[isyst][itwk] = new TH1D( Form("hmy_%d_%d",isyst,itwk), "", n, xmin, xmax);
       }       
    }
  }

  // Fill in nominal event fractions
  int curr_tree_num = -1;
  for(Long64_t iev = 0; iev < nmax; iev++) {
      genie_event_tree->GetEntry(iev); 
      EventRecord & event = *(mcrec->event);
      // Generated event files used as inputs in the validation programs
      // correspond to given neutrino+target initial states.
      // If the current file in the chain doesn't correspond to a desired 
      // initial state (as determined by looking at the first event in the file)
      // then skip this file althogether to save processing time.
      int tree_num = genie_event_tree->GetTreeNumber();
      if(curr_tree_num != tree_num) {
         curr_tree_num = tree_num;
         bool skip = !this->UseTree(event);
         if(skip) {
           LOG("gvldtest", pNOTICE) << "Skipping tree # : " << tree_num;
           iev += (genie_event_tree->GetTree()->GetEntries() - 1);
           continue;
         }//skip?
      }//new tree in chain
      GHepParticle * neutrino = event.Probe();
      assert(neutrino);
      double E = neutrino->P4()->E();
      double x = (inlogE) ? TMath::Log10(E) : E;
      if(this->IsModeX (event)) { hmx  ->Fill(x); }
      if(this->IsModeY (event)) { hmy  ->Fill(x); }
  }
  hmx->Divide(hmy);
  hmx->Smooth(2);

  // Include uncertainties if requested
  if(incl_err_band) {
     GSystSet & systlist = fRew.Systematics();
     for(unsigned int isyst = 0; isyst < fNuisanceParams.size(); isyst++) {
       systlist.Init(fNuisanceParams[isyst]);
       for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
          systlist.Set(fNuisanceParams[isyst], twk[itwk]);
          fRew.Reconfigure();
          int curr_tree_num = -1;
          for(Long64_t iev = 0; iev < nmax; iev++) {
              genie_event_tree->GetEntry(iev);
              EventRecord & event = *(mcrec->event);
              int tree_num = genie_event_tree->GetTreeNumber();
              if(curr_tree_num != tree_num) {
                curr_tree_num = tree_num;
                 bool skip = !this->UseTree(event);
                 if(skip) {           
                   LOG("gvldtest", pNOTICE) << "Skipping tree # : " << tree_num;
                   iev += (genie_event_tree->GetTree()->GetEntries() - 1);
                   continue;
                 }//skip?
              }//new tree in chain
              if(this->IsModeX(event)) {
                 GHepParticle * neutrino = event.Probe();
                 assert(neutrino);
                 double E = neutrino->P4()->E();
                 double x = (inlogE) ? TMath::Log10(E) : E;
                 double wght = fRew.CalcWeight(event);
                 hmxtwk[isyst][itwk]->Fill(x, wght);
              }//X?
              if(this->IsModeY(event)) {
                 GHepParticle * neutrino = event.Probe();
                 assert(neutrino);
                 double E = neutrino->P4()->E();
                 double x = (inlogE) ? TMath::Log10(E) : E;
                 double wght = fRew.CalcWeight(event);
                 hmytwk[isyst][itwk]->Fill(x, wght);
              }//X?
          }//iev
          hmxtwk[isyst][itwk]->Divide(hmytwk[isyst][itwk]);
          hmxtwk[isyst][itwk]->Smooth(2);
          systlist.Set(fNuisanceParams[isyst], 0.);//reset
          fRew.Reconfigure();
       }//itwk
     }//isyst
  }//incl_err_band?

  // Calculate cross-section ratio = f(E) and corresponding uncertainty.
  // Sources of uncertainty are taken to be uncorellated.
  double * energy_array = new double [n];
  double * R_array      = new double [n];
  double * R_array_errp = new double [n];
  double * R_array_errm = new double [n];
  for(int i = 0; i < n; i++) {
    int ibin = i+1;
    double energy = (inlogE) ? 
       TMath::Power(10., hmx->GetBinCenter(ibin)) : hmx->GetBinCenter(ibin);
    assert(energy>0.);
    double R = hmx->GetBinContent(ibin);
    double R_errp = 0;
    double R_errm = 0;
    if(incl_err_band) {
       double R_errp_s[knsyst_max];
       double R_errm_s[knsyst_max];
       for(unsigned int isyst = 0; isyst < fNuisanceParams.size(); isyst++) {
          R_errp_s[isyst] = 0;
          R_errm_s[isyst] = 0;
          double Rmax = -999999;
          double Rmin =  999999;
          for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
            Rmax = TMath::Max(Rmax, hmxtwk[isyst][itwk]->GetBinContent(ibin));
            Rmin = TMath::Min(Rmin, hmxtwk[isyst][itwk]->GetBinContent(ibin));
          }
          R_errp_s[isyst] = TMath::Max(0., Rmax-R);
          R_errm_s[isyst] = TMath::Max(0., R-Rmin);
          R_errp += TMath::Power(R_errp_s[isyst],2.);
          R_errm += TMath::Power(R_errm_s[isyst],2.);
       }//isyst
       R_errp = TMath::Sqrt(R_errp);
       R_errm = TMath::Sqrt(R_errm);       
    }//incl_err_band?
    // apply scaling factor, typically to normalize nuclear cross-section
    // to number of nucleons as expt cross-sections are typically quoted as such.
    LOG("gvldtest", pNOTICE)  
        << "R (E = " << energy << " GeV) = " << R 
        << " +" << ((R>0) ? 100*R_errp/R : 0) << "% "
        << " -" << ((R>0) ? 100*R_errm/R : 0) << "% ";    
    energy_array[i] = energy;
    R_array[i]      = (scale_with_E) ? R/energy       : R;
    R_array_errp[i] = (scale_with_E) ? R_errp/energy  : R_errp;
    R_array_errm[i] = (scale_with_E) ? R_errm/energy  : R_errm;
  }

  // Build cross-section graph
  TGraphAsymmErrors * model = new TGraphAsymmErrors(
     n,energy_array,R_array,0,0,R_array_errm,R_array_errp);
 
  // Clean-up
  delete hmx;
  delete hmy;
  for(unsigned int isyst = 0; isyst < knsyst_max; isyst++) {
    for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
       if(hmxtwk[isyst][itwk]) { delete hmxtwk[isyst][itwk]; }
       if(hmytwk[isyst][itwk]) { delete hmytwk[isyst][itwk]; }
    }
  }
  delete [] energy_array;
  delete [] R_array;
  delete [] R_array_errp;
  delete [] R_array_errm;

  return model;
}
//____________________________________________________________________________
CCpi0_CCQE::CCpi0_CCQE(int nupdg, int tgtpdg) :
XSecRatioForModesXY(),
fNuPdg  (nupdg),
fTgtPdg (tgtpdg)
{
  // List of nuisance params considered for CCpi0/CCQE cross-section

  fNuisanceParams.push_back(kXSecTwkDial_MaCCQE);
  if(pdg::IonPdgCodeToA(tgtpdg) > 2) {
     fNuisanceParams.push_back(kSystNucl_CCQEPauliSupViaKF);
  }
  fNuisanceParams.push_back(kXSecTwkDial_MaCCRES);
//fNuisanceParams.push_back(kXSecTwkDial_AhtBY);
//fNuisanceParams.push_back(kXSecTwkDial_BhtBY);
//fNuisanceParams.push_back(kXSecTwkDial_CV1uBY);
//fNuisanceParams.push_back(kXSecTwkDial_CV2uBY);
  fNuisanceParams.push_back(kXSecTwkDial_RvpCC1pi);
  fNuisanceParams.push_back(kXSecTwkDial_RvpCC2pi);
  fNuisanceParams.push_back(kXSecTwkDial_RvnCC1pi);
  fNuisanceParams.push_back(kXSecTwkDial_RvnCC2pi);

  // Add weight calculation engines
  fRew.AdoptWghtCalc( "xsec_ccqe",       new GReWeightNuXSecCCQE      );
  fRew.AdoptWghtCalc( "nuclear_qe",      new GReWeightFGM             );
  fRew.AdoptWghtCalc( "xsec_ccres",      new GReWeightNuXSecCCRES     );
  fRew.AdoptWghtCalc( "xsec_nonresbkg",  new GReWeightNonResonanceBkg );
//fRew.AdoptWghtCalc( "xsec_dis",        new GReWeightNuXSecDIS       );

  // Fine-tuning reweighting engines
  GReWeightNuXSecCCQE * rwccqe =      
    dynamic_cast<GReWeightNuXSecCCQE *> (fRew.WghtCalc("xsec_ccqe"));  
  rwccqe->SetMode(GReWeightNuXSecCCQE::kModeMa);
  GReWeightNuXSecCCRES * rwccres = 
      dynamic_cast<GReWeightNuXSecCCRES *> (fRew.WghtCalc("xsec_ccres"));  
  rwccres->SetMode(GReWeightNuXSecCCRES::kModeMaMv);
}
//.............................................................................
CCpi0_CCQE::~CCpi0_CCQE()
{
}
//.............................................................................
bool CCpi0_CCQE::IsModeX (EventRecord & event)
{
  Interaction * in = event.Summary();
  if(!in->ProcInfo().IsWeakCC()) return false;
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;
  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  if(tgtpdg == kPdgNeutron) {
    if(fTgtPdg != 1000000010) return false;
  }
  else if(tgtpdg == kPdgProton) {
    if(fTgtPdg != 1000010010) return false;
  }
  else {
    if(tgtpdg != fTgtPdg) return false;
  }

  bool has_pi0 = false;
  TObjArrayIter piter(&event);
  GHepParticle * p = 0;
  while ((p = (GHepParticle *) piter.Next())) {
    if(p->Status() != kIStStableFinalState) continue;
    if(p->Pdg() == kPdgPi0) {
       has_pi0 = true;
       break;
    }
  }
  if(!has_pi0) return false;

  return true;
}
//.............................................................................
bool CCpi0_CCQE::IsModeY (EventRecord & event)
{
  Interaction * in = event.Summary();
  if(!in->ProcInfo().IsWeakCC()) return false;
  if(!in->ProcInfo().IsQuasiElastic()) return false;
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;
  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  if(tgtpdg == kPdgNeutron) {
    if(fTgtPdg != 1000000010) return false;
  }
  else if(tgtpdg == kPdgProton) {
    if(fTgtPdg != 1000010010) return false;
  }
  else {
    if(tgtpdg != fTgtPdg) return false;
  }

  return true;
}
//.............................................................................
bool CCpi0_CCQE::UseTree(EventRecord & event)
{
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;

  GHepParticle * target = event.Particle(1);
  int tgtpdg = target->Pdg();
  if(pdg::IonPdgCodeToA(fTgtPdg) > 1) {
    if(tgtpdg == fTgtPdg) return true; 
  }
  else
  if(pdg::IonPdgCodeToA(fTgtPdg) == 1) {
    bool match =
       (tgtpdg == kPdgProton  && fTgtPdg == kPdgTgtFreeP) ||
       (tgtpdg == kPdgNeutron && fTgtPdg == kPdgTgtFreeN);
    if(match) return true;
  }
  return false;
}
//.............................................................................
string CCpi0_CCQE::Name (void) const
{
  return Form("CCpi0/CCQE;nu=%d;tgt=%d",fNuPdg,fTgtPdg);
}
//____________________________________________________________________________
CCIsoInclXSec::CCIsoInclXSec(int nupdg) :
NuXSecFunc(),
fNuPdg(nupdg)
{
  // List of nuisance params considered for CCQE cross-section
  fNuisanceParams.push_back(kXSecTwkDial_MaCCQE);
  fNuisanceParams.push_back(kXSecTwkDial_MaCCRES);
//fNuisanceParams.push_back(kXSecTwkDial_AhtBY);
//fNuisanceParams.push_back(kXSecTwkDial_BhtBY);
//fNuisanceParams.push_back(kXSecTwkDial_CV1uBY);
//fNuisanceParams.push_back(kXSecTwkDial_CV2uBY);
  if( pdg::IsNeutrino(nupdg) ) {
     fNuisanceParams.push_back(kXSecTwkDial_RvpCC1pi);
     fNuisanceParams.push_back(kXSecTwkDial_RvpCC2pi);
     fNuisanceParams.push_back(kXSecTwkDial_RvnCC1pi);
     fNuisanceParams.push_back(kXSecTwkDial_RvnCC2pi);
  }
  if( pdg::IsAntiNeutrino(nupdg) ) {
     fNuisanceParams.push_back(kXSecTwkDial_RvbarpCC1pi);
     fNuisanceParams.push_back(kXSecTwkDial_RvbarpCC2pi);
     fNuisanceParams.push_back(kXSecTwkDial_RvbarnCC1pi);
     fNuisanceParams.push_back(kXSecTwkDial_RvbarnCC2pi);
  }

  // Add weight calculation engines
  fRew.AdoptWghtCalc( "xsec_ccqe",       new GReWeightNuXSecCCQE      );
  fRew.AdoptWghtCalc( "xsec_ccres",      new GReWeightNuXSecCCRES     );
  fRew.AdoptWghtCalc( "xsec_nonresbkg",  new GReWeightNonResonanceBkg );
//fRew.AdoptWghtCalc( "xsec_dis",        new GReWeightNuXSecDIS       );

  // Fine-tune weight calculation engines
  GReWeightNuXSecCCQE * rwccqe = 
      dynamic_cast<GReWeightNuXSecCCQE *> (fRew.WghtCalc("xsec_ccqe"));  
  rwccqe->SetMode(GReWeightNuXSecCCQE::kModeMa);
  GReWeightNuXSecCCRES * rwccres = 
      dynamic_cast<GReWeightNuXSecCCRES *> (fRew.WghtCalc("xsec_ccres"));  
  rwccres->SetMode(GReWeightNuXSecCCRES::kModeMaMv);

}
//.............................................................................
CCIsoInclXSec::~CCIsoInclXSec()
{

}
//.............................................................................
string CCIsoInclXSec::Name(void) const 
{
  return Form("CCIsoIncl;nu=%d",fNuPdg);
}
//.............................................................................
TGraphAsymmErrors * CCIsoInclXSec::ExtractFromEventSample(
     int imodel, double Emin, double Emax,
     int n, bool inlogE, bool scale_with_E, bool incl_err_band)
{
  if(!fGenieInputs) {
     LOG("gvldtest", pERROR) << "No GENIE MC inputs";
     return 0;
  }
  TFile * genie_xsec_file = fGenieInputs->XSecFile(imodel);
  if(!genie_xsec_file) {
     LOG("gvldtest", pERROR) << "No input GENIE cross-section file";
     return 0;
  }
  TChain * genie_event_tree = fGenieInputs->EvtChain(imodel);
  if(!genie_event_tree) {
     LOG("gvldtest", pERROR) << "No input GENIE event tree";
     return 0;
  }

  // Get xsec directory and retrieve inclusive CC xsec
  string xsec_dir_name_vp = BuildXSecDirectoryName(fNuPdg, kPdgTgtFreeP);
  string xsec_dir_name_vn = BuildXSecDirectoryName(fNuPdg, kPdgTgtFreeN);
  TDirectory * xsec_dir_vp = 
      (TDirectory *) genie_xsec_file->Get(xsec_dir_name_vp.c_str());
  if(!xsec_dir_vp) {
     LOG("gvldtest", pERROR) 
        << "Can't find cross-section directory: " << xsec_dir_name_vp;
     return 0;
  }
  TGraph * cc_incl_xsec_vp = (TGraph*) xsec_dir_vp->Get("tot_cc");
  if(!cc_incl_xsec_vp) {
     LOG("gvldtest", pERROR) 
        << "Can't find inclusive CC cross-section calculation";
     return 0;
  }
  TDirectory * xsec_dir_vn = 
      (TDirectory *) genie_xsec_file->Get(xsec_dir_name_vn.c_str());
  if(!xsec_dir_vn) {
     LOG("gvldtest", pERROR) 
        << "Can't find cross-section directory: " << xsec_dir_name_vn;
     return 0;
  }
  TGraph * cc_incl_xsec_vn = (TGraph*) xsec_dir_vn->Get("tot_cc");
  if(!cc_incl_xsec_vn) {
     LOG("gvldtest", pERROR) 
        << "Can't find inclusive CC cross-section calculation";
     return 0;
  }

  // Set event tree branch address
  genie_event_tree->SetBranchStatus("gmcrec", 1);   
  Long64_t nmax = genie_event_tree->GetEntries();
  if (nmax<0) {
     LOG("gvldtest", pERROR) << "Number of events = 0";
     return 0;
  }
  LOG("gvldtest", pNOTICE) 
    << "Found " << nmax << " entries in the event tree";
  NtpMCEventRecord * mcrec = 0;
  genie_event_tree->SetBranchAddress("gmcrec", &mcrec);   
  if (!mcrec) {
     LOG("gvldtest", pERROR) << "Null MC record";
     return 0;
  }

  // Book histograms
  const unsigned int knsyst_max = 10;
  const unsigned int kntwk = 2;
  const double twk[kntwk] = { +1., -1. };
  double xmin = (inlogE) ? TMath::Log10(Emin) : Emin;
  double xmax = (inlogE) ? TMath::Log10(Emax) : Emax;
  TH1D * hnom_ccvp =  new TH1D("hnom_ccvp", "", n, xmin, xmax);
  TH1D * hnom_ccvn =  new TH1D("hnom_ccvn", "", n, xmin, xmax);
  TH1D * htwk_ccvp[knsyst_max][kntwk];
  TH1D * htwk_ccvn[knsyst_max][kntwk];
  for(unsigned int isyst = 0; isyst < knsyst_max; isyst++) {    
    for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
       htwk_ccvp[isyst][itwk] = 0;
       htwk_ccvn[isyst][itwk] = 0;
       if(isyst < fNuisanceParams.size()) {
          htwk_ccvp[isyst][itwk] = new TH1D( Form("htwk_ccvp_%d_%d",isyst,itwk), "", n, xmin, xmax);
          htwk_ccvn[isyst][itwk] = new TH1D( Form("htwk_ccvn_%d_%d",isyst,itwk), "", n, xmin, xmax);
       }       
    }
  }

  if(incl_err_band) {
     int curr_tree_num = -1;
     for(Long64_t iev = 0; iev < nmax; iev++) {
       genie_event_tree->GetEntry(iev);
       EventRecord & event = *(mcrec->event);
       int tree_num = genie_event_tree->GetTreeNumber();
       if(curr_tree_num != tree_num) {
          curr_tree_num = tree_num;
          bool skip = !this->UseTree(event);
          if(skip) {
            LOG("gvldtest", pNOTICE) << "Skipping tree # : " << tree_num;
            iev += (genie_event_tree->GetTree()->GetEntries() - 1);
            continue;
          }//skip?
       }//new tree in chain
       Interaction * in = event.Summary();
       if(!in->ProcInfo().IsWeakCC()) continue;
       GHepParticle * neutrino = event.Probe();
       if(neutrino->Pdg() != fNuPdg) continue;
       double E = neutrino->P4()->E();
       double x = (inlogE) ? TMath::Log10(E) : E;
       GHepParticle * target = event.Particle(1); 
       int tgtpdg = target->Pdg();
       if(tgtpdg == kPdgProton ) { hnom_ccvp->Fill(x); }
       if(tgtpdg == kPdgNeutron) { hnom_ccvn->Fill(x); }
     }//iev
     GSystSet & systlist = fRew.Systematics();
     for(unsigned int isyst = 0; isyst < fNuisanceParams.size(); isyst++) {
       systlist.Init(fNuisanceParams[isyst]);
       for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
          systlist.Set(fNuisanceParams[isyst], twk[itwk]);
          fRew.Reconfigure();
          curr_tree_num = -1;
          for(Long64_t iev = 0; iev < nmax; iev++) {
              genie_event_tree->GetEntry(iev);
              EventRecord & event = *(mcrec->event);
              int tree_num = genie_event_tree->GetTreeNumber();
              if(curr_tree_num != tree_num) {
                 curr_tree_num = tree_num;
                 bool skip = !this->UseTree(event);
                 if(skip) {
                   LOG("gvldtest", pNOTICE) << "Skipping tree # : " << tree_num;
                   iev += (genie_event_tree->GetTree()->GetEntries() - 1);
                   continue;
                 }//skip?
              }//new tree in chain
              Interaction * in = event.Summary();
              if(!in->ProcInfo().IsWeakCC()) continue;
              GHepParticle * neutrino = event.Probe();
              if(neutrino->Pdg() != fNuPdg) continue;
              GHepParticle * target = event.Particle(1); 
              int tgtpdg = target->Pdg();
              if(tgtpdg==kPdgProton || tgtpdg==kPdgNeutron) {
                 double E = neutrino->P4()->E();
                 double x = (inlogE) ? TMath::Log10(E) : E;
                 double wght = fRew.CalcWeight(event);              
                 if(tgtpdg == kPdgProton ) { htwk_ccvp[isyst][itwk]->Fill(x,wght); }
                 if(tgtpdg == kPdgNeutron) { htwk_ccvn[isyst][itwk]->Fill(x,wght); }
              }
          }//iev
          htwk_ccvp[isyst][itwk]->Divide(hnom_ccvp);
          htwk_ccvn[isyst][itwk]->Divide(hnom_ccvn);
          systlist.Set(fNuisanceParams[isyst], 0.);//reset
          fRew.Reconfigure();
       }//itwk
     }//isyst
  }//incl_err_band?

  // Calculate cross-section = f(E) and corresponding uncertainty.
  // Sources of uncertainty are taken to be uncorellated but they are 100% correlated
  // between the nu+p and nu+n components of the isoscalar neutrino cross-section.
  double * energy_array    = new double [n];
  double * xsec_array      = new double [n];
  double * xsec_array_errp = new double [n];
  double * xsec_array_errm = new double [n];
  for(int i = 0; i < n; i++) {
    int ibin = i+1;
    double energy = (inlogE) ? 
       TMath::Power(10., hnom_ccvp->GetBinCenter(ibin)) : 
       hnom_ccvp->GetBinCenter(ibin);
    assert(energy>0.);
    double xsec = 0.5 * (
        cc_incl_xsec_vp->Eval(energy) + cc_incl_xsec_vn->Eval(energy));
    double xsec_errp = 0;
    double xsec_errm = 0;
    if(incl_err_band) {
       double xsec_errp_s[knsyst_max];
       double xsec_errm_s[knsyst_max];
       for(unsigned int isyst = 0; isyst < fNuisanceParams.size(); isyst++) {
          double xsec_max = -99999;
          double xsec_min = +99999;
          for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
            double xsec_curr = 0.5 * (
                    cc_incl_xsec_vp->Eval(energy) * htwk_ccvp[isyst][itwk]->GetBinContent(ibin) + 
                    cc_incl_xsec_vn->Eval(energy) * htwk_ccvn[isyst][itwk]->GetBinContent(ibin)
            );
            xsec_max = TMath::Max(xsec_max, xsec_curr);
            xsec_min = TMath::Min(xsec_min, xsec_curr);
          }
          xsec_errp_s[isyst] = TMath::Max(0., xsec_max-xsec);
          xsec_errm_s[isyst] = TMath::Max(0., xsec-xsec_min);
          xsec_errp += TMath::Power(xsec_errp_s[isyst],2.);
          xsec_errm += TMath::Power(xsec_errm_s[isyst],2.);
       }//isyst
       xsec_errp = TMath::Sqrt(xsec_errp);
       xsec_errm = TMath::Sqrt(xsec_errm);       
    }//incl_err_band?
    LOG("gvldtest", pNOTICE)  
        << "xsec (E = " << energy << " GeV) = " << xsec << " 1E-38 cm^2 "
        << "+ " << ((xsec>0) ? 100*xsec_errp/xsec : 0) << "% "
        << "- " << ((xsec>0) ? 100*xsec_errm/xsec : 0) << "% ";
    energy_array[i]    = energy;
    xsec_array[i]      = (scale_with_E) ? xsec/energy       : xsec;
    xsec_array_errp[i] = (scale_with_E) ? xsec_errp/energy  : xsec_errp;
    xsec_array_errm[i] = (scale_with_E) ? xsec_errm/energy  : xsec_errm;
  }

  // Build cross-section graph
  TGraphAsymmErrors * model = new TGraphAsymmErrors(
     n,energy_array,xsec_array,0,0,xsec_array_errm, xsec_array_errp);

  // Clean-up
  delete hnom_ccvp;
  delete hnom_ccvn;
  for(unsigned int isyst = 0; isyst < knsyst_max; isyst++) {
    for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
      if (htwk_ccvp[isyst][itwk]) delete htwk_ccvp[isyst][itwk];
      if (htwk_ccvn[isyst][itwk]) delete htwk_ccvn[isyst][itwk];
    }
  }
  delete [] energy_array;
  delete [] xsec_array;
  delete [] xsec_array_errp;
  delete [] xsec_array_errm;

  LOG("gvldtest", pNOTICE) << "Returning CC inclusive model prediction";

  return model;
}
//.............................................................................
bool CCIsoInclXSec::UseTree(EventRecord & event)
{
  GHepParticle * neutrino = event.Probe();
  if(neutrino->Pdg() != fNuPdg) return false;
  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  if(!pdg::IsNucleon(tgtpdg)) return false;
  return true;
}
//____________________________________________________________________________
r::r() :
NuXSecFunc()
{
  // List of nuisance params considered for CCQE cross-section
  fNuisanceParams.push_back(kXSecTwkDial_MaCCQE);
  fNuisanceParams.push_back(kXSecTwkDial_MaCCRES);
//fNuisanceParams.push_back(kXSecTwkDial_AhtBY);
//fNuisanceParams.push_back(kXSecTwkDial_BhtBY);
//fNuisanceParams.push_back(kXSecTwkDial_CV1uBY);
//fNuisanceParams.push_back(kXSecTwkDial_CV2uBY);
  fNuisanceParams.push_back(kXSecTwkDial_RvpCC1pi);
  fNuisanceParams.push_back(kXSecTwkDial_RvpCC2pi);
  fNuisanceParams.push_back(kXSecTwkDial_RvnCC1pi);
  fNuisanceParams.push_back(kXSecTwkDial_RvnCC2pi);
  fNuisanceParams.push_back(kXSecTwkDial_RvbarpCC1pi);
  fNuisanceParams.push_back(kXSecTwkDial_RvbarpCC2pi);
  fNuisanceParams.push_back(kXSecTwkDial_RvbarnCC1pi);
  fNuisanceParams.push_back(kXSecTwkDial_RvbarnCC2pi);

  // Add weight calculation engines
  fRew.AdoptWghtCalc( "xsec_ccqe",       new GReWeightNuXSecCCQE      );
  fRew.AdoptWghtCalc( "xsec_ccres",      new GReWeightNuXSecCCRES     );
  fRew.AdoptWghtCalc( "xsec_nonresbkg",  new GReWeightNonResonanceBkg );
//fRew.AdoptWghtCalc( "xsec_dis",        new GReWeightNuXSecDIS       );

  // Fine-tune weight calculation engines
  GReWeightNuXSecCCQE * rwccqe = 
      dynamic_cast<GReWeightNuXSecCCQE *> (fRew.WghtCalc("xsec_ccqe"));  
  rwccqe->SetMode(GReWeightNuXSecCCQE::kModeMa);
  GReWeightNuXSecCCRES * rwccres = 
      dynamic_cast<GReWeightNuXSecCCRES *> (fRew.WghtCalc("xsec_ccres"));  
  rwccres->SetMode(GReWeightNuXSecCCRES::kModeMaMv);

}
//.............................................................................
r::~r()
{

}
//.............................................................................
string r::Name(void) const 
{
  return "r";
}
//.............................................................................
TGraphAsymmErrors * r::ExtractFromEventSample(
     int imodel, double Emin, double Emax,
     int n, bool inlogE, bool /*scale_with_E*/, bool incl_err_band)
{
  if(!fGenieInputs) return 0;
  TFile * genie_xsec_file = fGenieInputs->XSecFile(imodel);
  if(!genie_xsec_file) return 0;
  TChain * genie_event_tree = fGenieInputs->EvtChain(imodel);
  if(!genie_event_tree) return 0;

  // Get xsec directory and retrieve inclusive CC xsec
  string xsec_dir_name_vp    = BuildXSecDirectoryName(kPdgNuMu,     kPdgTgtFreeP);
  string xsec_dir_name_vn    = BuildXSecDirectoryName(kPdgNuMu,     kPdgTgtFreeN);
  string xsec_dir_name_vbarp = BuildXSecDirectoryName(kPdgAntiNuMu, kPdgTgtFreeP);
  string xsec_dir_name_vbarn = BuildXSecDirectoryName(kPdgAntiNuMu, kPdgTgtFreeN);

  TDirectory * xsec_dir_vp = 
      (TDirectory *) genie_xsec_file->Get(xsec_dir_name_vp.c_str());
  if(!xsec_dir_vp) {
     LOG("gvldtest", pERROR) 
        << "Can't find cross-section directory: " << xsec_dir_name_vp;
     return 0;
  }
  TGraph * cc_incl_xsec_vp = (TGraph*) xsec_dir_vp->Get("tot_cc");
  if(!cc_incl_xsec_vp) {
     LOG("gvldtest", pERROR) 
        << "Can't find inclusive CC cross-section calculation";
     return 0;
  }
  TDirectory * xsec_dir_vn = 
      (TDirectory *) genie_xsec_file->Get(xsec_dir_name_vn.c_str());
  if(!xsec_dir_vn) {
     LOG("gvldtest", pERROR) 
        << "Can't find cross-section directory: " << xsec_dir_name_vn;
     return 0;
  }
  TGraph * cc_incl_xsec_vn = (TGraph*) xsec_dir_vn->Get("tot_cc");
  if(!cc_incl_xsec_vn) {
     LOG("gvldtest", pERROR) 
        << "Can't find inclusive CC cross-section calculation";
     return 0;
  }


  TDirectory * xsec_dir_vbarp = 
      (TDirectory *) genie_xsec_file->Get(xsec_dir_name_vbarp.c_str());
  if(!xsec_dir_vbarp) {
     LOG("gvldtest", pERROR) 
        << "Can't find cross-section directory: " << xsec_dir_name_vbarp;
     return 0;
  }
  TGraph * cc_incl_xsec_vbarp = (TGraph*) xsec_dir_vbarp->Get("tot_cc");
  if(!cc_incl_xsec_vbarp) {
     LOG("gvldtest", pERROR) 
        << "Can't find inclusive CC cross-section calculation";
     return 0;
  }
  TDirectory * xsec_dir_vbarn = 
      (TDirectory *) genie_xsec_file->Get(xsec_dir_name_vbarn.c_str());
  if(!xsec_dir_vbarn) {
     LOG("gvldtest", pERROR) 
        << "Can't find cross-section directory: " << xsec_dir_name_vbarn;
     return 0;
  }
  TGraph * cc_incl_xsec_vbarn = (TGraph*) xsec_dir_vbarn->Get("tot_cc");
  if(!cc_incl_xsec_vbarn) {
     LOG("gvldtest", pERROR) 
        << "Can't find inclusive CC cross-section calculation";
     return 0;
  }

  // Set event tree branch address
  genie_event_tree->SetBranchStatus("gmcrec", 1);   
  Long64_t nmax = genie_event_tree->GetEntries();
  if (nmax<0) {
     LOG("gvldtest", pERROR) << "Number of events = 0";
     return 0;
  }
  LOG("gvldtest", pNOTICE) 
    << "Found " << nmax << " entries in the event tree";
  NtpMCEventRecord * mcrec = 0;
  genie_event_tree->SetBranchAddress("gmcrec", &mcrec);   
  if (!mcrec) {
     LOG("gvldtest", pERROR) << "Null MC record";
     return 0;
  }

  // Book histograms
  const unsigned int knsyst_max = 10;
  const unsigned int kntwk = 2;
  const double twk[kntwk] = { +1., -1. };
  double xmin = (inlogE) ? TMath::Log10(Emin) : Emin;
  double xmax = (inlogE) ? TMath::Log10(Emax) : Emax;
  TH1D * hnom_ccvp    =  new TH1D("hnom_ccvp",    "", n, xmin, xmax);
  TH1D * hnom_ccvn    =  new TH1D("hnom_ccvn",    "", n, xmin, xmax);
  TH1D * hnom_ccvbarp =  new TH1D("hnom_ccvbarp", "", n, xmin, xmax);
  TH1D * hnom_ccvbarn =  new TH1D("hnom_ccvbarn", "", n, xmin, xmax);
  TH1D * htwk_ccvp   [knsyst_max][kntwk];
  TH1D * htwk_ccvn   [knsyst_max][kntwk];
  TH1D * htwk_ccvbarp[knsyst_max][kntwk];
  TH1D * htwk_ccvbarn[knsyst_max][kntwk];
  for(unsigned int isyst = 0; isyst < knsyst_max; isyst++) {    
    for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
       htwk_ccvp   [isyst][itwk] = 0;
       htwk_ccvn   [isyst][itwk] = 0;
       htwk_ccvbarp[isyst][itwk] = 0;
       htwk_ccvbarn[isyst][itwk] = 0;
       if(isyst < fNuisanceParams.size()) {
          htwk_ccvp   [isyst][itwk] = new TH1D( Form("htwk_ccvp_%d_%d",isyst,itwk),    "", n, xmin, xmax);
          htwk_ccvn   [isyst][itwk] = new TH1D( Form("htwk_ccvn_%d_%d",isyst,itwk),    "", n, xmin, xmax);
          htwk_ccvbarp[isyst][itwk] = new TH1D( Form("htwk_ccvbarp_%d_%d",isyst,itwk), "", n, xmin, xmax);
          htwk_ccvbarn[isyst][itwk] = new TH1D( Form("htwk_ccvbarn_%d_%d",isyst,itwk), "", n, xmin, xmax);
       }       
    }
  }

  if(incl_err_band) {
     int curr_tree_num = -1;
     for(Long64_t iev = 0; iev < nmax; iev++) {
       genie_event_tree->GetEntry(iev);
       EventRecord & event = *(mcrec->event);
       int tree_num = genie_event_tree->GetTreeNumber();
       if(curr_tree_num != tree_num) {
          curr_tree_num = tree_num;
          bool skip = !this->UseTree(event);
          if(skip) {
            LOG("gvldtest", pNOTICE) << "Skipping tree # : " << tree_num;
            iev += (genie_event_tree->GetTree()->GetEntries() - 1);
            continue;
          }//skip?
       }//new tree in chain
       Interaction * in = event.Summary();
       if(!in->ProcInfo().IsWeakCC()) continue;
       GHepParticle * neutrino = event.Probe();
       int nupdg = neutrino->Pdg();
       if(nupdg != kPdgNuMu && nupdg != kPdgAntiNuMu) continue;
       double E = neutrino->P4()->E();
       double x = (inlogE) ? TMath::Log10(E) : E;
       GHepParticle * target = event.Particle(1); 
       int tgtpdg = target->Pdg();
       if(tgtpdg == kPdgProton  && nupdg == kPdgNuMu    ) { hnom_ccvp    -> Fill(x); }
       if(tgtpdg == kPdgNeutron && nupdg == kPdgNuMu    ) { hnom_ccvn    -> Fill(x); }
       if(tgtpdg == kPdgProton  && nupdg == kPdgAntiNuMu) { hnom_ccvbarp -> Fill(x); }
       if(tgtpdg == kPdgNeutron && nupdg == kPdgAntiNuMu) { hnom_ccvbarn -> Fill(x); }
     }//iev
     GSystSet & systlist = fRew.Systematics();
     for(unsigned int isyst = 0; isyst < fNuisanceParams.size(); isyst++) {
       systlist.Init(fNuisanceParams[isyst]);
       for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
          systlist.Set(fNuisanceParams[isyst], twk[itwk]);
          fRew.Reconfigure();
          curr_tree_num = -1;
          for(Long64_t iev = 0; iev < nmax; iev++) {
              genie_event_tree->GetEntry(iev);
              EventRecord & event = *(mcrec->event);
              int tree_num = genie_event_tree->GetTreeNumber();
              if(curr_tree_num != tree_num) {
                 curr_tree_num = tree_num;
                 bool skip = !this->UseTree(event);
                 if(skip) {
                   LOG("gvldtest", pNOTICE) << "Skipping tree # : " << tree_num;
                   iev += (genie_event_tree->GetTree()->GetEntries() - 1);
                   continue;
                 }//skip?
              }//new tree in chain
              Interaction * in = event.Summary();
              if(!in->ProcInfo().IsWeakCC()) continue;
              GHepParticle * neutrino = event.Probe();
              int nupdg = neutrino->Pdg();
              if(nupdg != kPdgNuMu && nupdg != kPdgAntiNuMu) continue;
              GHepParticle * target = event.Particle(1); 
              int tgtpdg = target->Pdg();
              if(tgtpdg==kPdgProton || tgtpdg==kPdgNeutron) {
                 double E = neutrino->P4()->E();
                 double x = (inlogE) ? TMath::Log10(E) : E;
                 double wght = fRew.CalcWeight(event);              
                 if(tgtpdg == kPdgProton  && nupdg == kPdgNuMu    ) { htwk_ccvp   [isyst][itwk]->Fill(x,wght); }
                 if(tgtpdg == kPdgNeutron && nupdg == kPdgNuMu    ) { htwk_ccvn   [isyst][itwk]->Fill(x,wght); }
                 if(tgtpdg == kPdgProton  && nupdg == kPdgAntiNuMu) { htwk_ccvbarp[isyst][itwk]->Fill(x,wght); }
                 if(tgtpdg == kPdgNeutron && nupdg == kPdgAntiNuMu) { htwk_ccvbarn[isyst][itwk]->Fill(x,wght); }
              }
          }//iev
          htwk_ccvp   [isyst][itwk]->Divide(hnom_ccvp);
          htwk_ccvn   [isyst][itwk]->Divide(hnom_ccvn);
          htwk_ccvbarp[isyst][itwk]->Divide(hnom_ccvbarp);
          htwk_ccvbarn[isyst][itwk]->Divide(hnom_ccvbarn);
          systlist.Set(fNuisanceParams[isyst], 0.);//reset
          fRew.Reconfigure();
       }//itwk
     }//isyst
  }//incl_err_band?

  // Calculate cross-section = f(E) and corresponding uncertainty.
  // Sources of uncertainty are taken to be uncorellated (even Rijk params for neutrinos and anti-neutrinos) 
  // but they are 100% correlated between the nu+p, nu+n, nubar+p, nubar+n components of r
  double * energy_array = new double [n];
  double * r_array      = new double [n];
  double * r_array_errp = new double [n];
  double * r_array_errm = new double [n];
  for(int i = 0; i < n; i++) {
    int ibin = i+1;
    double energy = (inlogE) ? 
       TMath::Power(10., hnom_ccvp->GetBinCenter(ibin)) : 
       hnom_ccvp->GetBinCenter(ibin);
    assert(energy>0.);
    double xsec_nu    = 0.5 * (cc_incl_xsec_vp   ->Eval(energy) + cc_incl_xsec_vn   ->Eval(energy));
    double xsec_nubar = 0.5 * (cc_incl_xsec_vbarp->Eval(energy) + cc_incl_xsec_vbarn->Eval(energy));
    double r = (xsec_nu>0) ? xsec_nubar/xsec_nu : 0;
    double r_errp = 0;
    double r_errm = 0;
    if(incl_err_band) {
       double r_errp_s[knsyst_max];
       double r_errm_s[knsyst_max];
       for(unsigned int isyst = 0; isyst < fNuisanceParams.size(); isyst++) {
          double r_max = -99999;
          double r_min = +99999;
          for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
            double xsec_nu_curr = 0.5 * (
                    cc_incl_xsec_vp->Eval(energy) * htwk_ccvp[isyst][itwk]->GetBinContent(ibin) + 
                    cc_incl_xsec_vn->Eval(energy) * htwk_ccvn[isyst][itwk]->GetBinContent(ibin)
            );
            double xsec_nubar_curr = 0.5 * (
                    cc_incl_xsec_vbarp->Eval(energy) * htwk_ccvbarp[isyst][itwk]->GetBinContent(ibin) + 
                    cc_incl_xsec_vbarn->Eval(energy) * htwk_ccvbarn[isyst][itwk]->GetBinContent(ibin)
            );
            double r_curr = (xsec_nu_curr>0) ? xsec_nubar_curr/xsec_nu_curr : 0;
            r_max = TMath::Max(r_max, r_curr);
            r_min = TMath::Min(r_min, r_curr);
          }
          r_errp_s[isyst] = TMath::Max(0., r_max-r);
          r_errm_s[isyst] = TMath::Max(0., r-r_min);
          r_errp += TMath::Power(r_errp_s[isyst],2.);
          r_errm += TMath::Power(r_errm_s[isyst],2.);
       }//isyst
       r_errp = TMath::Sqrt(r_errp);
       r_errm = TMath::Sqrt(r_errm);       
    }//incl_err_band?
    LOG("gvldtest", pNOTICE)  
        << "r (E = " << energy << " GeV) = " << r 
        << "+ " << ((r>0) ? 100*r_errp/r : 0) << "% "
        << "- " << ((r>0) ? 100*r_errm/r : 0) << "% ";
    energy_array[i] = energy;
    r_array[i]      = r;
    r_array_errp[i] = r_errp;
    r_array_errm[i] = r_errm;
  }

  // Build cross-section graph
  TGraphAsymmErrors * model = new TGraphAsymmErrors(
     n,energy_array,r_array,0,0,r_array_errm, r_array_errp);

  // Clean-up
  delete hnom_ccvp;
  delete hnom_ccvn;
  delete hnom_ccvbarp;
  delete hnom_ccvbarn;
  for(unsigned int isyst = 0; isyst < knsyst_max; isyst++) {
    for(unsigned int itwk = 0; itwk < kntwk; itwk++) {
      if (htwk_ccvp   [isyst][itwk]) delete htwk_ccvp   [isyst][itwk];
      if (htwk_ccvn   [isyst][itwk]) delete htwk_ccvn   [isyst][itwk];
      if (htwk_ccvbarp[isyst][itwk]) delete htwk_ccvbarp[isyst][itwk];
      if (htwk_ccvbarn[isyst][itwk]) delete htwk_ccvbarn[isyst][itwk];
    }
  }
  delete [] energy_array;
  delete [] r_array;
  delete [] r_array_errp;
  delete [] r_array_errm;

  return model;
}
//.............................................................................
bool r::UseTree(EventRecord & event)
{
  GHepParticle * neutrino = event.Probe();
  int nupdg = neutrino->Pdg();
  if( nupdg != kPdgNuMu && nupdg != kPdgAntiNuMu) return false;
  GHepParticle * target = event.Particle(1); 
  int tgtpdg = target->Pdg();
  if(!pdg::IsNucleon(tgtpdg)) return false;
  return true;
}
//____________________________________________________________________________
