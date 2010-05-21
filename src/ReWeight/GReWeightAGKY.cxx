//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Authors: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
          STFC, Rutherford Appleton Laboratory

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Sep 10, 2009 - CA
   Skeleton first included in v2.5.1.
 @ May 21, 2010 - CA
   Added option to reweiht xF and pT2 for nucleon+pion states
*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TF1.h>
#include <TMath.h>

#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "ReWeight/GReWeightAGKY.h"
#include "ReWeight/GReWeightUtils.h"
#include "ReWeight/GSystUncertainty.h"

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GReWeightAGKY::GReWeightAGKY() :
GReWeightI()
{
  this->Init();
}
//_______________________________________________________________________________________
GReWeightAGKY::~GReWeightAGKY()
{
  delete fBaryonXFpdf;   
  delete fBaryonPT2pdf; 
  delete fBaryonXFpdfTwk;
  delete fBaryonPT2pdfTwk;
}
//_______________________________________________________________________________________
bool GReWeightAGKY::IsHandled(GSyst_t syst)
{
  switch(syst) {
     case ( kHadrAGKYTwkDial_xF1pi ) :
     case ( kHadrAGKYTwkDial_pT1pi ) :
       return true;
       break;
     default:
       return false;
       break;
  }

  return false;
}
//_______________________________________________________________________________________
void GReWeightAGKY::SetSystematic(GSyst_t syst, double twk_dial)
{
  switch(syst) {
     case ( kHadrAGKYTwkDial_xF1pi ) :
       fPeakBaryonXFTwkDial = twk_dial;
       break;
     case ( kHadrAGKYTwkDial_pT1pi ) :
       fAvgPT2TwkDial = twk_dial;
       break;
     default:
       return;
       break;
  }
}
//_______________________________________________________________________________________
void GReWeightAGKY::Reset(void)
{
   fPeakBaryonXFTwkDial = 0.;
   fAvgPT2TwkDial       = 0.;
}
//_______________________________________________________________________________________
void GReWeightAGKY::Reconfigure(void)
{

}
//_______________________________________________________________________________________
double GReWeightAGKY::CalcWeight(const EventRecord & event) 
{  
  Interaction * interaction = event.Summary();

  bool is_cc  = interaction->ProcInfo().IsWeakCC();
  bool is_nc  = interaction->ProcInfo().IsWeakNC();
  if(is_cc && !fRewCC) return 1.;   
  if(is_nc && !fRewNC) return 1.;
       
  int nupdg = interaction->InitState().ProbePdg();
  if(nupdg==kPdgNuMu     && !fRewNumu   ) return 1.;
  if(nupdg==kPdgAntiNuMu && !fRewNumubar) return 1.;
  if(nupdg==kPdgNuE      && !fRewNue    ) return 1.;
  if(nupdg==kPdgAntiNuE  && !fRewNuebar ) return 1.;

  // Skip events not handled by the AGKY hadronization model
  if(! utils::rew::HadronizedByAGKY(event)) return 1.0;

  // Skip high-W events hadronized by AGKY/PYTHIA
  if(utils::rew::HadronizedByAGKYPythia(event)) return 1.0;

  //
  // Reweight events hadronized by AGKY/KNO
  //

  double wght =      
    this->RewxFpT1pi(event);  /* reweight pT2 and xF for `nucleon+pion' states */
    
  return wght;
}
//_______________________________________________________________________________________
double GReWeightAGKY::RewxFpT1pi(const EventRecord & event) 
{
  bool PT2tweaked = (TMath::Abs(fAvgPT2TwkDial)       > controls::kASmallNum);
  bool XFtweaked  = (TMath::Abs(fPeakBaryonXFTwkDial) > controls::kASmallNum);

  bool tweaked = (PT2tweaked || XFtweaked);
  if(!tweaked) return 1.;

  //
  // Did the hadronization code produced a `nucleon+pion' system?
  // If yes, keep the nucleon xF and pT (in HCM) for event reweighting.
  //
  bool is_Npi = false;
  TLorentzVector p4N(0.,0.,0.,0.);

  //
  // ...
  // ...
  // ...
  //

/*
  GHepParticle * p = 0;
  TIter event_iter(&event);
  while ( (p = dynamic_cast<GHepParticle *>(event_iter.Next())) ) {
	int imom = p->FirstMother();
        if(imom==-1) continue;
        if(event.Particle(imom)->Pdg() != kPdgHadronicSyst) continue;

  }
*/

  if(!is_Npi) return 1.;

  // Get the hadronic system 4-momentum and boost velocity (LAB -> HCM)
  TLorentzVector p4hadsyst = utils::rew::Hadronic4pLAB(event);
  TVector3 bv = -1 * p4hadsyst.BoostVector();

  // Get baryon xF and pT2 in HCM
  double XF  = 0;
  double PT2 = 0;
  p4N.Boost(bv);
  
  // ...
  // ...
  // ...

  // Calculate event weight
  
  GSystUncertainty * uncertainty = GSystUncertainty::Instance();

  double XFwght = 1.;
  bool XFinrange = (XF > fXFmin && XF < fXFmax);
  if(XFinrange && XFtweaked) {
    double frcerr = uncertainty->OneSigmaErr(kHadrAGKYTwkDial_xF1pi);
    double XFpeak = fDefPeakBaryonXF * (1. + fPeakBaryonXFTwkDial * frcerr);
    fBaryonXFpdfTwk->SetParameter(1, XFpeak);
    double I = fBaryonXFpdfTwk->Integral(fXFmin,fXFmax);
    if(I>0) {
       double norm = 1/I;
       fBaryonXFpdfTwk->SetParameter(0, norm);
       double def = fBaryonXFpdf    -> Eval(XF);
       double twk = fBaryonXFpdfTwk -> Eval(XF);
       if(def>0 && twk>=0) {
         XFwght = twk/def;
       }	    
    }
  }
  
  double PT2wght = 1.;
  bool PT2inrange = (PT2 > fPT2min && PT2 < fPT2max);
  if(PT2inrange && PT2tweaked) {
    double frcerr = uncertainty->OneSigmaErr(kHadrAGKYTwkDial_pT1pi);
    double PT2avg = fDefAvgPT2* (1. + fAvgPT2TwkDial * frcerr);
    fBaryonPT2pdfTwk->SetParameter(1, PT2avg);
    double I = fBaryonPT2pdfTwk->Integral(fPT2min,fPT2max);
    if(I>0.) {
       double norm = 1/I;
       fBaryonPT2pdfTwk->SetParameter(0, norm);
       double def = fBaryonPT2pdf    -> Eval(PT2);
       double twk = fBaryonPT2pdfTwk -> Eval(PT2);
       if(def>0 && twk>=0) {
         PT2wght = twk/def;
       }	    
    }
  }

  double wght = PT2wght * XFwght;
  return wght;
}
//_______________________________________________________________________________________
double GReWeightAGKY::CalcChisq(void)
{
  double chisq =
    TMath::Power(fPeakBaryonXFTwkDial, 2.) +
    TMath::Power(fAvgPT2TwkDial,       2.);

  return chisq;
}
//_______________________________________________________________________________________
void GReWeightAGKY::Init(void)
{
  this->RewNue    (true);
  this->RewNuebar (true);
  this->RewNumu   (true);
  this->RewNumubar(true);
  this->RewCC     (true);
  this->RewNC     (true);

  // Baryon pT^2 and xF parameterizations used as PDFs in AGKY.
  // Code from KNOHadronization.cxx:
  //

  fXFmin = -1.0;
  fXFmax =  0.5;
  fPT2min = 0;
  fPT2max = 0.6;

  fBaryonXFpdf  = new TF1("fBaryonXFpdf",
                   "0.083*exp(-0.5*pow(x+0.385,2.)/0.131)", fXFmin, fXFmax);  
  fBaryonPT2pdf = new TF1("fBaryonPT2pdf", 
                   "exp(-0.214-6.625*x)", fPT2min, fPT2max);  
  //
  // Now, define same function as above, but insert tweaking dials.
  // - For the xF  PDF add a parameter to shift the peak of the distribution.
  // - For the pT2 PDF add a parameter to shift the <pT2>
  // Include parameter to allow re-normalizing tweaked PDF to 1.

  fDefPeakBaryonXF = -0.385;
  fBaryonXFpdfTwk  = new TF1("fBaryonXFpdf",
                   "[0]*0.083*exp(-0.5*pow(x-[1],2.)/0.131)", fXFmin, fXFmax);
  fBaryonXFpdfTwk->SetParameter(0, 1.); // norm
  fBaryonXFpdfTwk->SetParameter(1, fDefPeakBaryonXF);

  fDefAvgPT2 = 1./6.625;
  fBaryonPT2pdfTwk = new TF1("fBaryonPT2pdf", 
                   "[0]*exp(-0.214-x/[1])", fPT2min, fPT2max);  
  fBaryonPT2pdfTwk->SetParameter(0, 1.); // norm
  fBaryonPT2pdfTwk->SetParameter(1,fDefAvgPT2);

  // init tweaking dials
  fPeakBaryonXFTwkDial = 0.;
  fAvgPT2TwkDial       = 0.;
}
//_______________________________________________________________________________________

