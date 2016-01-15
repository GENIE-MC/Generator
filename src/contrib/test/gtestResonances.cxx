//____________________________________________________________________________
/*!

\program gtestResonances

\brief   Program used for testing/debugging the Breit-Wigner distributions 
         for Baryon Resonances

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created June 20, 2004

\cpright Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
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

#include "BaryonResonance/BaryonResonance.h"
#include "BaryonResonance/BaryonResUtils.h"
#include "Conventions/Constants.h"
#include "Conventions/Units.h"
#include "Messenger/Messenger.h"
#include "Utils/BWFunc.h"

using std::string;

using namespace genie;
using namespace genie::constants;
using namespace genie::utils;
using namespace genie::units;

int main(int /*argc*/, char ** /*argv*/)
{
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
 TGraph * grbw[kNResonances];

 //-- Define an output ntuple
 TNtuple resnt("resnt","baryon resonances","res:W:BW");

 //-- Loop over baryon resonances
 for(int ires = 0; ires < kNResonances; ires++) {

    //-- Get the current baryon resonance
    Resonance_t resonance = kResonances[ires];

    LOG("test", pINFO) << "@ Resonance = " << utils::res::AsString(resonance);

    int    LR  = utils::res::OrbitalAngularMom (resonance);
    double MR  = utils::res::Mass              (resonance);
    double WR  = utils::res::Width             (resonance);
    double NR  = utils::res::BWNorm            (resonance);

    LOG("test", pINFO) << "- mass    = " << MR << " GeV";
    LOG("test", pINFO) << "- width   = " << WR << " GeV";
    LOG("test", pINFO) << "- BW norm = " << NR;
    LOG("test", pINFO) << "- L       = " << LR;

    double bw[50], W[50];

    for(int iW = 0; iW<50; iW++) {
       W[iW]  = 0.01 + iW*(3.5/50.);
       bw[iW] = utils::bwfunc::BreitWignerL(W[iW],LR,MR,WR,NR);
       resnt.Fill((float)resonance, W[iW], bw[iW]);
       LOG("test", pINFO) 
         << "  BreitWigner(W = " << W[iW] << " GeV) = " << bw[iW];
    }
    grbw[ires] = new TGraph(50, W, bw);
 }

 //-- Create & format output canvas

 TCanvas * c = new TCanvas("c","",20,20,500,500);

 c->SetBorderMode(0);
 c->SetFillColor(0);
 c->Draw();

 TH1F * hframe = (TH1F*) c->DrawFrame(1.0,0.5,3.5,6.0);
 hframe->Draw();
 TLegend * legend = new TLegend(0.6,0.2,0.8,0.8);

 for(int ires = 0; ires < kNResonances; ires++) {
    grbw[ires]->SetLineWidth(2);
    grbw[ires]->SetLineColor( color[ires] );
    grbw[ires]->Draw("C");
    legend->AddEntry(grbw[ires],utils::res::AsString(kResonances[ires]),"L");
 }

 legend->SetFillColor(0);
 legend->Draw();

 hframe->GetXaxis()->SetTitle("W (GeV)");
 hframe->GetYaxis()->SetTitle("Breit Wigner Amplitude");
 hframe->GetXaxis()->SetTitleOffset(1.3);
 hframe->GetYaxis()->SetTitleOffset(1.3);

 c->Update();

 TFile f("./genie-resonances.root","recreate");
 c->Write("plot");
 resnt.Write();
 f.Close();

 for(int ires = 0; ires < kNResonances; ires++) {
  delete grbw[ires];
 }
 delete legend;
 delete c;
}
//____________________________________________________________________________

