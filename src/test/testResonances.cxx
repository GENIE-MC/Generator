//____________________________________________________________________________
/*!

\program testResonances

\brief   test program used for testing/debugging the Breit-Wigner distributions 
         for Baryon Resonances

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created June 20, 2004

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TPad.h>
#include <TPaveLabel.h>

#include "Algorithm/Algorithm.h"
#include "Algorithm/AlgFactory.h"
#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "BaryonResonance/BreitWignerI.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"

using std::string;

using namespace genie;
using namespace genie::constants;
using namespace genie::units;

int main(int /*argc*/, char ** /*argv*/)
{
 AlgFactory * algf = AlgFactory::Instance();

 const int kNResonances = 18;
 const Resonance_t kResonances[kNResonances] = {

    kP33_1232, kS11_1535, kD13_1520,
    kS11_1650, kD13_1700, kD15_1675,
    kS31_1620, kD33_1700, kP11_1440,
    kP33_1600, kP13_1720, kF15_1680,
    kP31_1910, kP33_1920, kF35_1905,
    kF37_1950, kP11_1710, kF17_1970
 };

 const int color[kNResonances] = {1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19};

 //-- Graphs to hold baryon resonances Breit Wigner distributions
 TGraph * gr_simplistic[kNResonances];
 TGraph * gr_realistic [kNResonances];

 //-- Define an output ntuple
 TNtuple resnt("resnt","baryon resonances","res:bwtype:W:BW");

 //-- Loop over baryon resonances
 for(int ires = 0; ires < kNResonances; ires++) {

    //-- Get the current baryon resonance
    Resonance_t resonance = kResonances[ires];

    //-- Get Breit Wigner algorithms
    const BreitWignerI * bw_simplistic =
              dynamic_cast<const BreitWignerI *>
                     (algf->GetAlgorithm("genie::BreitWignerRes","Default"));
    const BreitWignerI * bw_realistic =
              dynamic_cast<const BreitWignerI *>
                    (algf->GetAlgorithm("genie::BreitWignerLRes","Default"));

    double data_bw_simplistic[50], data_bw_realistic[50], W[50];

    for(int iW = 0; iW<50; iW++) {

       W[iW] = 0.01 + iW*(3.5/50.);

       data_bw_simplistic [iW] = bw_simplistic -> Eval (resonance,W[iW]);
       data_bw_realistic  [iW] = bw_realistic  -> Eval (resonance,W[iW]);

       resnt.Fill((float)resonance,0.,W[iW],data_bw_simplistic[iW]);
       resnt.Fill((float)resonance,1.,W[iW],data_bw_realistic[iW]);

       LOG("VldRes", pINFO) << "W = " << W[iW]
                            << " BWs = " << data_bw_simplistic [iW]
                            << " BWr = " << data_bw_realistic  [iW];
    }

    gr_simplistic[ires] = new TGraph(50, W, data_bw_simplistic);
    gr_realistic [ires] = new TGraph(50, W, data_bw_realistic );
 }

 //-- Create & format output canvas

 TCanvas * c = new TCanvas("c","",20,20,500,500);

 c->SetBorderMode(0);
 c->SetFillColor(0);
 c->Draw();

 //-- Define & format TPads

 TPad * upper_data_pad = new TPad("upper_data_pad","",0.05,0.51,0.95,0.99);
 TPad * lower_data_pad = new TPad("lower_data_pad","",0.05,0.01,0.95,0.49);

 upper_data_pad -> Draw();
 lower_data_pad -> Draw();

 upper_data_pad -> SetFillColor(0);
 lower_data_pad -> SetFillColor(0);
 upper_data_pad -> SetBorderMode(0);
 lower_data_pad -> SetBorderMode(0);

 //-- UPPER PAD = simplistic Breit Wigner distributions

 upper_data_pad->cd();
 TH1F * hframe1 = (TH1F*) upper_data_pad->DrawFrame(1.0,0.5,3.5,6.0);
 hframe1->Draw();
 TLegend * legend1 = new TLegend(0.6,0.2,0.8,0.8);

 for(int ires = 0; ires < kNResonances; ires++) {
    gr_simplistic[ires]->SetLineWidth(2);
    gr_simplistic[ires]->SetLineColor( color[ires] );
    gr_simplistic[ires]->Draw("C");
    legend1->AddEntry(gr_simplistic[ires],utils::res::AsString(kResonances[ires]),"L");
 }

 legend1->SetFillColor(0);
 legend1->Draw();

 //-- Format Axis

 hframe1->GetXaxis()->SetTitle("W (GeV)");
 hframe1->GetYaxis()->SetTitle("Breit Wigner Amplitude");
 hframe1->GetXaxis()->SetTitleOffset(1.3);
 hframe1->GetYaxis()->SetTitleOffset(1.3);

 //-- LOWER PAD = realistic Breit Wigner distributions

 lower_data_pad->cd();
 TH1F * hframe2 = (TH1F*) lower_data_pad->DrawFrame(1.0,0.5,3.5,6.0);
 hframe2->Draw();
 TLegend * legend2 = new TLegend(0.6,0.2,0.8,0.8);

 for(int ires = 0; ires < kNResonances; ires++) {
    gr_realistic[ires]->SetLineWidth(2);
    gr_realistic[ires]->SetLineColor( color[ires] );
    gr_realistic[ires]->Draw("C");
    legend2->AddEntry(gr_realistic[ires],utils::res::AsString(kResonances[ires]),"L");
 }

 legend2->SetFillColor(0);
 legend2->Draw();

 //-- Format Axis
 hframe2->GetXaxis()->SetTitle("W (GeV)");
 hframe2->GetYaxis()->SetTitle("Breit Wigner Amplitude");
 hframe2->GetXaxis()->SetTitleOffset(1.3);
 hframe2->GetYaxis()->SetTitleOffset(1.3);

 c->Update();

 TFile f("./genie-resonances.root","recreate");
 c->Write("plot");
 resnt.Write();
 f.Close();

 for(int ires = 0; ires < kNResonances; ires++) {
  delete gr_simplistic[ires];
  delete gr_realistic[ires];
 }
 delete legend1;
 delete legend2;
 delete upper_data_pad;
 delete lower_data_pad;
 delete c;
}
//____________________________________________________________________________

