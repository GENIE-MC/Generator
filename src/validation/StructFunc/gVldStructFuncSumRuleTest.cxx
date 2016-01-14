//____________________________________________________________________________
/*!

\program gvld_sf_sum_rule_test

\brief   Structure function sum rule test

         Syntax:
           gvld_sf_sum_rule_test
                [--resonance-xsec-model genie::model_name/model_config]
                [--dis-xsec-model       genie::model_name/model_config]
                [--dis-charm-xsec-model genie::model_name/model_config]

         Options:

           --resonance-xsec-model
              Specify GENIE resonance cross-section model.

           --dis-xsec-model
              Specify GENIE DIS cross-section model.

           --dis-charm-xsec-model
              Specify GENIE DIS charm production cross-section model.

         Example:

            % gvld_sf_sum_rule_test \
                --resonance-xsec-model genie::ReinSehgalRESPXSec/Default \
                --dis-xsec-model genie::QPMDISPXSec/Default \
                --dis-charm-xsec-model genie::AivazisCharmPXSecLO/CC-Default

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created May 30, 2009

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>

#include <TSystem.h>
#include <TFile.h>
#include <TGraph.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPavesText.h>
#include <TText.h>
#include <TLegend.h>

#include "Algorithm/AlgFactory.h"
#include "Base/XSecAlgorithmI.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/StringUtils.h"
#include "Utils/SystemUtils.h"
#include "Utils/Style.h"
#include "validation/StructFunc/StructFunc.h"

using std::ostringstream;
using std::string;

using namespace genie;
using namespace genie::mc_vs_data;

// globals
string gOptRESModelName        = ""; // --resonance-xsec-model command-line argument
string gOptDISModelName        = ""; // --dis-xsec-model command-line argument
string gOptDISCharmModelName   = ""; // --dis-charm-xsec-model command-line argument

const XSecAlgorithmI * gRESXSecModel      = 0; // resonance cross-section model
const XSecAlgorithmI * gDISXSecModel      = 0; // DIS cross-section model
const XSecAlgorithmI * gDISCharmXSecModel = 0; // DIS cross-section model (charm production)

StructFunc gStructFuncCalc;   // utility class extracting structure functions from the cross-section model

TPostScript *  gPS = 0;
TCanvas *      gC  = 0;

// constants

const int    kNQ2      = 50;
const double kQ2min    = 1E-1;
const double kQ2max    = 1E+2;
const double klogQ2min = TMath::Log(kQ2min);
const double klogQ2max = TMath::Log(kQ2max);
const double kdlogQ2   = (klogQ2max-klogQ2min)/(kNQ2-1);
const int    kNx       = 1000;
const double kxmin     = 1E-6;
const double kxmax     = 0.999999;
const double klogxmin  = TMath::Log10(kxmin);
const double klogxmax  = TMath::Log10(kxmax);
const double kdlogx    = (klogxmax-klogxmin)/(kNx-1);

// function prototypes

void GetCommandLineArgs             (int argc, char ** argv);
void PrintSyntax                    (void);
void Init                           (void);
void End                            (void);
void TestAdlerSumRule               (void);
void TestGrossLlewellynSmithSumRule (void);
void TestGottfriedSumRule           (void);
void TestBjorkenSumRule             (void);

//_________________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs (argc,argv);
  Init();
  TestAdlerSumRule();
  TestGrossLlewellynSmithSumRule();
  TestGottfriedSumRule();
  TestBjorkenSumRule();
  End();
  LOG("gvldtest", pINFO)  << "Done!";

  return 0;
}
//_________________________________________________________________________________
void Init(void)
{
  LOG("gvldtest", pNOTICE) << "Initializing...";

  // Set GENIE style
  utils::style::SetDefaultStyle();

  //
  // Get cross-section models
  //
  AlgFactory * algf = AlgFactory::Instance();
  gRESXSecModel = 0;
  if(gOptRESModelName != "none"){
      vector<string> modelv = utils::str::Split(gOptRESModelName,"/");
      assert(modelv.size()==2);
      string model_name = modelv[0];
      string model_conf = modelv[1];
      gRESXSecModel =dynamic_cast<const XSecAlgorithmI *> (
         algf->GetAlgorithm(model_name, model_conf));
  }
  gDISXSecModel = 0;
  if(gOptDISModelName != "none"){
       vector<string> modelv = utils::str::Split(gOptDISModelName,"/");
       assert(modelv.size()==2);
       string model_name = modelv[0];
       string model_conf = modelv[1];
       gDISXSecModel = dynamic_cast<const XSecAlgorithmI *> (
          algf->GetAlgorithm(model_name, model_conf));
  }
  gDISCharmXSecModel = 0;
  if(gOptDISCharmModelName != "none"){
       vector<string> modelv = utils::str::Split(gOptDISCharmModelName,"/");
       assert(modelv.size()==2);
       string model_name = modelv[0];
       string model_conf = modelv[1];
       gDISCharmXSecModel = dynamic_cast<const XSecAlgorithmI *> (
          algf->GetAlgorithm(model_name, model_conf));
  }

  gStructFuncCalc.SetResonanceXSecModel (gRESXSecModel);
  gStructFuncCalc.SetDISXSecModel       (gDISXSecModel);
  gStructFuncCalc.SetDISCharmXSecModel  (gDISCharmXSecModel);

  // Create plot canvas
  gC = new TCanvas("c","",20,20,500,650);
  gC->SetBorderMode(0);
  gC->SetFillColor(0);
  gC->SetGridx();
  gC->SetGridy();

  // Get local time to tag outputs
  string lt_for_filename   = utils::system::LocalTimeAsString("%02d.%02d.%02d_%02d.%02d.%02d");
  string lt_for_cover_page = utils::system::LocalTimeAsString("%02d/%02d/%02d %02d:%02d:%02d");

  // Create output postscript file
  string filename  = Form("genie-structfunc_sumrules-%s.ps",lt_for_filename.c_str());
  gPS = new TPostScript(filename.c_str(), 111);

  // Add cover page
  gPS->NewPage();
  gC->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText(" ");
  hdr.AddText("GENIE structure function sum rule tests");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(lt_for_cover_page.c_str());
  hdr.AddText(" ");
  hdr.Draw();
  gC->Update();

}
//_________________________________________________________________________________
void End(void)
{
  LOG("gvldtest", pNOTICE) << "Cleaning up...";

  gPS->Close();

  delete gC;
  delete gPS;
}
//_________________________________________________________________________________
void TestAdlerSumRule(void)
{
//
// ***********************************************
// * Adler sum rule:                             *
// *                                             *
// *          /1                                 *
// * S_{A} =  |dx  (F2{vn} - F2{vp}) / 2x  = 1.  *
// *          /0                                 *
// ***********************************************
//

  double S  [kNQ2];
  double Q2 [kNQ2];
  for(int i=0; i<kNQ2; i++) {
    Q2[i] = TMath::Exp(klogQ2min + i * kdlogQ2);
    double x [kNx];
    double I [kNx];
    for(int j=0; j<kNx; j++) {
       x[j] = TMath::Exp(klogxmin + j * kdlogx);
       double temp1, temp2;
       double F2vn = 0.;
       double F2vp = 0.;
       gStructFuncCalc.ExtractF1F2xF3(x[j], Q2[i], kPdgNuMu, kPdgNeutron, temp1, F2vn, temp2);
       gStructFuncCalc.ExtractF1F2xF3(x[j], Q2[i], kPdgNuMu, kPdgProton,  temp1, F2vp, temp2);
       I[j] = (F2vn - F2vp)/2.;
    }//j
    S[i] = 0.5*(I[0]+I[kNx-1]);
    for(int j=1; j<kNx-1; j++) {
      S[i] += (I[j] * (j%2+1));
    }
    S[i] *= (2.*kdlogx/3.);
  }//i
  TGraph * gr = new  TGraph(kNQ2,Q2,S);
  gPS-> NewPage();
  gC -> SetFillColor(0);
  gC -> SetBorderMode(0);
  gC -> SetLogx(1);
  TH1F * hframe = (TH1F*) gC->DrawFrame(kQ2min,0,kQ2max,1.5);
  hframe->GetXaxis()->SetTitle("Q^{2} (GeV^{2}/c^{4})");
  hframe->GetYaxis()->SetTitle("S_{A}");
  gr -> Draw("l");
  TLatex * t = new TLatex();
  t->SetTextColor(kBlack);
  t->SetTextFont(42);
  t->SetTextSize(0.03);
  t->DrawLatex(kQ2min,1.6,
    "Adler sum rule: S_{A} = #int_{0}^{1}dx (F_{2}(#nun)-F_{2}(#nup))/2x = 1");
  gC -> Update();
}
//_________________________________________________________________________________
void TestGrossLlewellynSmithSumRule(void)
{
//
// *******************************************
// * Gross - Llewellyn Smith sum rule:       *
// *                                         *
// *            /1                           *
// * S_{GLS} =  | dx xF3{vN} / 2x = 3.       * 
// *           /0                            *
// *******************************************
//
  double S  [kNQ2];
  double Q2 [kNQ2];
  for(int i=0; i<kNQ2; i++) {
    Q2[i] = TMath::Exp(klogQ2min + i * kdlogQ2);
    double x [kNx];
    double I [kNx];
    for(int j=0; j<kNx; j++) {
       x[j] = TMath::Exp(klogxmin + j * kdlogx);
       double temp1, temp2;
       double xF3vn = 0.;
       double xF3vp = 0.;
       gStructFuncCalc.ExtractF1F2xF3(x[j], Q2[i], kPdgNuMu, kPdgNeutron, temp1, temp2, xF3vn);
       gStructFuncCalc.ExtractF1F2xF3(x[j], Q2[i], kPdgNuMu, kPdgProton,  temp1, temp2, xF3vp);
       double xF3 = (xF3vn+xF3vp)/2.;
       I[j] = xF3/2.;
    }//j
    S[i] = 0.5*(I[0]+I[kNx-1]);
    for(int j=1; j<kNx-1; j++) {
      S[i] += (I[j] * (j%2+1));
    }
    S[i] *= (2.*kdlogx/3.);
  }//i
  TGraph * gr = new  TGraph(kNQ2,Q2,S);
  gPS-> NewPage();
  gC -> SetFillColor(0);
  gC -> SetBorderMode(0);
  gC -> SetLogx(1);
  TH1F * hframe = (TH1F*) gC->DrawFrame(kQ2min,0,kQ2max,4);
  hframe->GetXaxis()->SetTitle("Q^{2} (GeV^{2}/c^{4})");
  hframe->GetYaxis()->SetTitle("S_{GLS}");
  gr -> Draw("l");
  TLatex * t = new TLatex();
  t->SetTextColor(kBlack);
  t->SetTextFont(42);
  t->SetTextSize(0.03);
  t->DrawLatex(kQ2min,4.1,
    "Gross - Llewellyn Smith sum rule: S_{GLS} = #int_{0}^{1}dx xF_{3}(#nuN)/2x = 3");
  gC -> Update();
}
//_________________________________________________________________________________
void TestGottfriedSumRule(void)
{
//
// ************************************************
// * Gottfried sum rule:                          *
// *                                              *
// *         /1                                   *
// * S_{G} = | dx (F2{mup} - F2{mun}) / x = 1./3. * 
// *        /0                                    *
// ************************************************
//


}
//_________________________________________________________________________________
void TestBjorkenSumRule(void)
{

}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gvldtest", pNOTICE) << "*** Parsing command line arguments";

  CmdLnArgParser parser(argc,argv);

  // Get GENIE model names to be used
  if(parser.OptionExists("resonance-xsec-model")){
     gOptRESModelName = parser.Arg("resonance-xsec-model");
  }
  if(parser.OptionExists("dis-xsec-model")){
     gOptDISModelName = parser.Arg("dis-xsec-model");
  }
  if(parser.OptionExists("dis-charm-xsec-model")){
     gOptDISCharmModelName = parser.Arg("dis-charm-xsec-model");
  }

}
//_________________________________________________________________________________
void PrintSyntax(void)
{

}
//_________________________________________________________________________________
