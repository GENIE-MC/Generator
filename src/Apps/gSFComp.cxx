//____________________________________________________________________________
/*!

\program gsfcomp

\brief   Structure function comparison tool

\syntax  gsfcomp --structure-func sf_set [-o output]

         --structure-func :
          Specifies a comma separated list of GENIE structure function models.
          The full algorithm name and configuration should be provided for
          each GENIE model as in `genie::model_name/model_config'.
          (unique color and line styles are defined for up to 4 sets)

         -o :
          Specifies a name to be used in the output files.
          Default: sf_comp

\example gsfcomp --structure-func genie::Blah/Default,genie::Blah/Tweaked

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created Feb 10, 2016

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <vector>
#include <string>

#include <TNtuple.h>
#include <TFile.h>
#include <TMath.h>
#include <TPostScript.h>
#include <TPavesText.h>
#include <TText.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TPaletteAxis.h>

#include "Framework/Algorithm/AlgFactory.h"
#include "Framework/Messenger/Messenger.h"
#include "Physics/DeepInelastic/XSection/DISStructureFunc.h"
#include "Physics/DeepInelastic/XSection/DISStructureFuncModelI.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/Style.h"

using namespace std;
using namespace genie;
using namespace genie::utils;

// globals
string        gOptSF      = "";         // --structure-func argument
string        gOptOutFile = "sf_comp";  // -o argument
vector<const DISStructureFuncModelI *> gSFAlgList;

// function prototypes
void GetCommandLineArgs (int argc, char ** argv);
void GetAlgorithms      (void);
void MakePlots          (void);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  utils::style::SetDefaultStyle();

  GetCommandLineArgs (argc,argv);   // Get command line arguments
  GetAlgorithms();                  // Get requested SF algorithms
  MakePlots();   // Produce all output plots and fill output n-tuple

  return 0;
}
//_________________________________________________________________________________
void MakePlots (void)
{
  const unsigned int nm = 4; // number of models
  int    col [nm] = { kBlack, kRed+1,  kBlue-3, kGreen+2 };
  int    sty [nm] = { kSolid, kDashed, kDashed, kDashed  };
  int    mrk [nm] = { 20,     20,      20,      20       };
  double msz [nm] = { 0.7,    0.7,     0.7,     0.7      };
  const char * opt   [nm] = { "ap", "l", "l", "l" };
  const char * lgopt [nm] = { "P",  "L", "L", "L" };

  // Q2 range and values for 1-D plots
  const unsigned int nQ2 = 20;
  const double Q2min = 1E-1; // GeV^2
  const double Q2max = 1E+3; // GeV^2
  const double log10Q2min = TMath::Log10(Q2min);
  const double log10Q2max = TMath::Log10(Q2max);
  const double dlog10Q2 = (log10Q2max-log10Q2min)/(nQ2-1);
  double Q2_arr [nQ2];
  for(unsigned int iq2 = 0; iq2 < nQ2; iq2++) {
     double Q2 = TMath::Power(10, log10Q2min + iq2*dlog10Q2);
     Q2_arr[iq2] = Q2;
  }

  // x values for 1-D plots
  const unsigned int nx = 22;
  double x_arr [nx] = {
    0.0001, 0.0010, 0.0100, 0.0250, 0.0500, 
    0.0750, 0.1000, 0.1500, 0.2000, 0.2500, 
    0.3500, 0.4000, 0.4500, 0.5000, 0.5500, 
    0.6000, 0.7000, 0.7500, 0.8000, 0.8500, 
    0.9000, 0.9500
  };

  // Q2 range and values for 2-D plots
  const unsigned int nQ2_2d = 50;
  const double Q2min_2d = 1E-1; // GeV^2
  const double Q2max_2d = 1E+3; // GeV^2
  const double log10Q2min_2d = TMath::Log10(Q2min_2d);
  const double log10Q2max_2d = TMath::Log10(Q2max_2d);
  const double dlog10Q2_2d = (log10Q2max_2d-log10Q2min_2d)/(nQ2_2d-1);
  double Q2_bin_edges_2d [nQ2_2d];
  for(unsigned int iq2 = 0; iq2 < nQ2_2d; iq2++) {
     double Q2 = TMath::Power(10, log10Q2min_2d + iq2*dlog10Q2_2d);
     Q2_bin_edges_2d[iq2] = Q2;
  }

  // x range and values for 2-D plots
  const unsigned int nx_2d = 50;
  const double xmin_2d = 1E-4; 
  const double xmax_2d = 0.95; 
  const double log10xmin_2d = TMath::Log10(xmin_2d);
  const double log10xmax_2d = TMath::Log10(xmax_2d);
  const double dlog10x_2d = (log10xmax_2d-log10xmin_2d)/(nx_2d-1);
  double x_bin_edges_2d [nx_2d];
  for(unsigned int ix = 0; ix < nx_2d; ix++) {
     double x = TMath::Power(10, log10xmin_2d + ix*dlog10x_2d);
     x_bin_edges_2d[ix] = x;
  }

  // Output ntuple
  TNtuple * ntpl = new TNtuple("nt","structure functions","i:F1:F2:F3:F4:F5:F6:x:Q2");

  // Canvas for output plots
  TCanvas * cnv = new TCanvas("c","",20,20,500,650);
  cnv->SetBorderMode(0);
  cnv->SetFillColor(0);

  // Create output ps file
  string ps_filename = gOptOutFile + ".ps";
  TPostScript * ps = new TPostScript(ps_filename.c_str(), 111);

  // Front page
  ps->NewPage();
  cnv->Range(0,0,100,100);
  TPavesText hdr(10,40,90,70,3,"tr");
  hdr.AddText("GENIE structure function comparisons");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText("Models used:");
  hdr.AddText(" ");
  for(unsigned int im=0; im < gSFAlgList.size(); im++) {
      const char * label = gSFAlgList[im]->Id().Key().c_str();
      hdr.AddText(label);
  }
  hdr.AddText(" ");
  hdr.Draw();
  cnv->Update();

  ps->NewPage();

  TLatex * tex = new TLatex;
  TLegend * lgnd = 0;

  //
  // Plots
  //

  // structure functions = f(Q^2) for selected vales of x (1-D plots)
  TGraph * gr_F1_Q2  [nm] = { 0 };
  TGraph * gr_F2_Q2  [nm] = { 0 }; 
  TGraph * gr_F3_Q2  [nm] = { 0 }; 
  TGraph * gr_F4_Q2  [nm] = { 0 }; 
  TGraph * gr_F5_Q2  [nm] = { 0 }; 
  TGraph * gr_F6_Q2  [nm] = { 0 }; 

  // structure functions = f(x,Q^2) with fine x,Q2 binning (2-D plots)
  TH2D * h2_F1 [nm] = { 0 };
  TH2D * h2_F2 [nm] = { 0 };
  TH2D * h2_F3 [nm] = { 0 };
  TH2D * h2_F4 [nm] = { 0 };
  TH2D * h2_F5 [nm] = { 0 };
  TH2D * h2_F6 [nm] = { 0 };

  // structure functions ratios relative to first one specified = f(x,Q^2) with fine x,Q2 binning (2-D plots)
  TH2D * h2_F1_r [nm] = { 0 };
  TH2D * h2_F2_r [nm] = { 0 };
  TH2D * h2_F3_r [nm] = { 0 };
  TH2D * h2_F4_r [nm] = { 0 };
  TH2D * h2_F5_r [nm] = { 0 };
  TH2D * h2_F6_r [nm] = { 0 };

  //
  // Generate tructure function plots as f(Q^2) for selected vales of x (1-D plots)
  //

  for(unsigned int ix=0; ix < nx; ix++) {
    double x = x_arr[ix];
    double F1_arr [nm][nQ2];
    double F2_arr [nm][nQ2];
    double F3_arr [nm][nQ2];
    double F4_arr [nm][nQ2];
    double F5_arr [nm][nQ2];
    double F6_arr [nm][nQ2];
    for(unsigned int im=0; im < gSFAlgList.size(); im++) {
      DISStructureFunc sf;
      sf.SetModel(gSFAlgList[im]);
      for(unsigned int iq2 = 0; iq2 < nQ2; iq2++) {
        double Q2 = Q2_arr[iq2];
        Interaction * interaction = Interaction::DISCC(kPdgTgtFreeP,kPdgProton,kPdgNuMu);
      //LOG("gsfcomp", pNOTICE) << "Setting x = " << x << ", Q2 = " << Q2;
        interaction->KinePtr()->Setx(x);
        interaction->KinePtr()->SetQ2(Q2);
        sf.Calculate(interaction);
        double F1 = sf.F1();
        double F2 = sf.F2();
        double F3 = sf.F3();
        double F4 = sf.F4();
        double F5 = sf.F5();
        double F6 = sf.F6();
        F1_arr [im][iq2] = F1;
        F2_arr [im][iq2] = F2;
        F3_arr [im][iq2] = F3;
        F4_arr [im][iq2] = F4;
        F5_arr [im][iq2] = F5;
        F6_arr [im][iq2] = F6;
        ntpl->Fill(im,F1,F2,F3,F4,F5,F6,x,Q2);
        delete interaction;
      }//iq2

      gr_F1_Q2 [im] = new TGraph (nQ2, Q2_arr, F1_arr [im]);
      gr_F2_Q2 [im] = new TGraph (nQ2, Q2_arr, F2_arr [im]);
      gr_F3_Q2 [im] = new TGraph (nQ2, Q2_arr, F3_arr [im]);
      gr_F4_Q2 [im] = new TGraph (nQ2, Q2_arr, F4_arr [im]);
      gr_F5_Q2 [im] = new TGraph (nQ2, Q2_arr, F5_arr [im]);
      gr_F6_Q2 [im] = new TGraph (nQ2, Q2_arr, F6_arr [im]);

      genie::utils::style::Format( gr_F1_Q2 [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);
      genie::utils::style::Format( gr_F2_Q2 [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);
      genie::utils::style::Format( gr_F3_Q2 [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);
      genie::utils::style::Format( gr_F4_Q2 [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);
      genie::utils::style::Format( gr_F5_Q2 [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);
      genie::utils::style::Format( gr_F6_Q2 [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);

      gr_F1_Q2 [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");
      gr_F2_Q2 [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");
      gr_F3_Q2 [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");
      gr_F4_Q2 [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");
      gr_F5_Q2 [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");
      gr_F6_Q2 [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");

      gr_F1_Q2 [im] -> GetYaxis() -> SetTitle("F1(x,Q^{2})");
      gr_F2_Q2 [im] -> GetYaxis() -> SetTitle("F2(x,Q^{2})");
      gr_F3_Q2 [im] -> GetYaxis() -> SetTitle("F3(x,Q^{2})");
      gr_F4_Q2 [im] -> GetYaxis() -> SetTitle("F4(x,Q^{2})");
      gr_F5_Q2 [im] -> GetYaxis() -> SetTitle("F5(x,Q^{2})");
      gr_F6_Q2 [im] -> GetYaxis() -> SetTitle("F6(x,Q^{2})");

    }//im

    ps->NewPage();

    if(x<0.3) {
       lgnd = new TLegend(0.20, 0.55, 0.50, 0.85);
    } else {
       lgnd = new TLegend(0.60, 0.55, 0.80, 0.85);
    }
    lgnd -> SetFillColor(0);
    lgnd -> SetFillStyle(0);
    lgnd -> SetBorderSize(0);

    lgnd->Clear();
    for(unsigned int im=0; im < gSFAlgList.size(); im++) {
      const char * label = gSFAlgList[im]->Id().Key().c_str();
      lgnd->AddEntry(gr_F1_Q2 [im], Form("%s",label), lgopt[im]);
    }

    cnv->Divide(2,3);

    cnv->cd(1); gPad->SetLogx();
    cnv->cd(2); gPad->SetLogx();
    cnv->cd(3); gPad->SetLogx();
    cnv->cd(4); gPad->SetLogx();
    cnv->cd(5); gPad->SetLogx();
    cnv->cd(6); gPad->SetLogx();

    for(unsigned int im=0; im < gSFAlgList.size(); im++) {
      cnv->cd(1); gr_F1_Q2 [im]->Draw(opt[im]);
      cnv->cd(2); gr_F2_Q2 [im]->Draw(opt[im]);
      cnv->cd(3); gr_F3_Q2 [im]->Draw(opt[im]);
      cnv->cd(4); gr_F4_Q2 [im]->Draw(opt[im]);
      cnv->cd(5); gr_F5_Q2 [im]->Draw(opt[im]);
      cnv->cd(6); gr_F6_Q2 [im]->Draw(opt[im]);
    }

    cnv->cd(1); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, Form("x = %.3e",x));
    cnv->cd(2); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, Form("x = %.3e",x));
    cnv->cd(3); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, Form("x = %.3e",x));
    cnv->cd(4); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, Form("x = %.3e",x));
    cnv->cd(5); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, Form("x = %.3e",x));
    cnv->cd(6); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, Form("x = %.3e",x));

    cnv->Update();
  }

  //
  // Plot structure functions = f(x,Q2)
  //

  for(unsigned int im=0; im < gSFAlgList.size(); im++) {
    h2_F1 [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_F2 [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_F3 [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_F4 [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_F5 [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_F6 [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    DISStructureFunc sf;
    sf.SetModel(gSFAlgList[im]);
    for(int ibinx = 1; 
            ibinx <= h2_F1[im]->GetXaxis()->GetNbins(); ibinx++) {
      double x = h2_F1[im]->GetXaxis()->GetBinCenter(ibinx);
      for(int ibinq2 = 1; 
              ibinq2 <= h2_F1[im]->GetYaxis()->GetNbins(); ibinq2++) {
         double Q2 = h2_F1[im]->GetYaxis()->GetBinCenter(ibinq2);
         Interaction * interaction = Interaction::DISCC(kPdgTgtFreeP,kPdgProton,kPdgNuMu);
       //LOG("gsfcomp", pNOTICE) << "Setting x (bin " << ibinx << ") = " << x << ", Q2 (bin " << ibinq2 << ") = " << Q2;
         interaction->KinePtr()->Setx(x);
         interaction->KinePtr()->SetQ2(Q2);
         sf.Calculate(interaction);
         double F1 = sf.F1();
         double F2 = sf.F2();
         double F3 = sf.F3();
         double F4 = sf.F4();
         double F5 = sf.F5();
         double F6 = sf.F6();
         h2_F1 [im] -> SetBinContent(ibinx, ibinq2, F1);
         h2_F2 [im] -> SetBinContent(ibinx, ibinq2, F2); 
         h2_F3 [im] -> SetBinContent(ibinx, ibinq2, F3); 
         h2_F4 [im] -> SetBinContent(ibinx, ibinq2, F4); 
         h2_F5 [im] -> SetBinContent(ibinx, ibinq2, F5); 
         h2_F6 [im] -> SetBinContent(ibinx, ibinq2, F6); 
         delete interaction;
      }
    }

    ps->NewPage();

    h2_F1 [im] -> GetXaxis() -> SetTitle("x");
    h2_F2 [im] -> GetXaxis() -> SetTitle("x");
    h2_F3 [im] -> GetXaxis() -> SetTitle("x");
    h2_F4 [im] -> GetXaxis() -> SetTitle("x");
    h2_F5 [im] -> GetXaxis() -> SetTitle("x");
    h2_F6 [im] -> GetXaxis() -> SetTitle("x");

    h2_F1 [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_F2 [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_F3 [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_F4 [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_F5 [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_F6 [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");

    cnv->Divide(2,3);

    TPaletteAxis * palette = 0;
    tex->SetTextSize(0.03);

    cnv->cd(1); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F1 [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F1[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, Form("F1: %s", gSFAlgList[im]->Id().Key().c_str()) );

    cnv->cd(2); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F2 [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F2[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, Form("F2: %s", gSFAlgList[im]->Id().Key().c_str()) );

    cnv->cd(3); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F3 [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F3[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, Form("F3: %s", gSFAlgList[im]->Id().Key().c_str()) );

    cnv->cd(4); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F4 [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F4[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, Form("F4: %s", gSFAlgList[im]->Id().Key().c_str()) );

    cnv->cd(5); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F5[im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F5[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, Form("F5: %s", gSFAlgList[im]->Id().Key().c_str()) );

    cnv->cd(6); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F6[im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F6[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, Form("F6: %s", gSFAlgList[im]->Id().Key().c_str()) );

    cnv->Update();
  }

  //
  // In the case of multiple input structure functions, 
  // plot the ratio of each structure function = f(x,Q2) to the first one
  // 

  for(unsigned int im=1; im < gSFAlgList.size(); im++) {
    h2_F1_r [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_F2_r [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_F3_r [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_F4_r [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_F5_r [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_F6_r [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    for(int ibinx = 1; 
            ibinx <= h2_F1[im]->GetXaxis()->GetNbins(); ibinx++) {
      for(int ibinq2 = 1; 
              ibinq2 <= h2_F1[im]->GetYaxis()->GetNbins(); ibinq2++) {

         double F1  = h2_F1 [im] -> GetBinContent(ibinx, ibinq2);
         double F2  = h2_F2 [im] -> GetBinContent(ibinx, ibinq2);
         double F3  = h2_F3 [im] -> GetBinContent(ibinx, ibinq2);
         double F4  = h2_F4 [im] -> GetBinContent(ibinx, ibinq2);
         double F5  = h2_F5 [im] -> GetBinContent(ibinx, ibinq2);
         double F6  = h2_F6 [im] -> GetBinContent(ibinx, ibinq2);

         double F10 = h2_F1 [0]  -> GetBinContent(ibinx, ibinq2);
         double F20 = h2_F2 [0]  -> GetBinContent(ibinx, ibinq2);
         double F30 = h2_F3 [0]  -> GetBinContent(ibinx, ibinq2);
         double F40 = h2_F4 [0]  -> GetBinContent(ibinx, ibinq2);
         double F50 = h2_F5 [0]  -> GetBinContent(ibinx, ibinq2);
         double F60 = h2_F6 [0]  -> GetBinContent(ibinx, ibinq2);

         double F1r = ((F10 > 0.) ? (F1-F10)/F10 : 0.);
         double F2r = ((F20 > 0.) ? (F2-F20)/F20 : 0.);
         double F3r = ((F30 > 0.) ? (F3-F30)/F30 : 0.);
         double F4r = ((F40 > 0.) ? (F4-F40)/F40 : 0.);
         double F5r = ((F50 > 0.) ? (F5-F50)/F50 : 0.);
         double F6r = ((F60 > 0.) ? (F6-F60)/F60 : 0.);

         h2_F1_r [im] -> SetBinContent(ibinx, ibinq2, F1r);
         h2_F2_r [im] -> SetBinContent(ibinx, ibinq2, F2r); 
         h2_F3_r [im] -> SetBinContent(ibinx, ibinq2, F3r); 
         h2_F4_r [im] -> SetBinContent(ibinx, ibinq2, F4r); 
         h2_F5_r [im] -> SetBinContent(ibinx, ibinq2, F5r); 
         h2_F6_r [im] -> SetBinContent(ibinx, ibinq2, F6r); 
      }
    }

    ps->NewPage();

    h2_F1_r [im] -> GetXaxis() -> SetTitle("x");
    h2_F2_r [im] -> GetXaxis() -> SetTitle("x");
    h2_F3_r [im] -> GetXaxis() -> SetTitle("x");
    h2_F4_r [im] -> GetXaxis() -> SetTitle("x");
    h2_F5_r [im] -> GetXaxis() -> SetTitle("x");
    h2_F6_r [im] -> GetXaxis() -> SetTitle("x");

    h2_F1_r [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_F2_r [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_F3_r [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_F4_r [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_F5_r [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_F6_r [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");

    cnv->Divide(2,3);

    TPaletteAxis * palette = 0;

    cnv->cd(1); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F1_r [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F1_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, Form("F1: (%s-%s)/(%s)", 
      gSFAlgList[im]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str()) );

    cnv->cd(2); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F2_r [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F2_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, Form("F2: (%s-%s)/(%s)", 
      gSFAlgList[im]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str()) );

    cnv->cd(3); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F3_r [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F3_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, Form("F3: (%s-%s)/(%s)", 
      gSFAlgList[im]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str()) );

    cnv->cd(4); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F4_r [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F4_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, Form("F4: (%s-%s)/(%s)", 
      gSFAlgList[im]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str()) );

    cnv->cd(5); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F5_r[im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F5_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, Form("F5: (%s-%s)/(%s)", 
      gSFAlgList[im]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str()) );

    cnv->cd(6); gPad->SetLogx(); gPad->SetLogy(); 
    h2_F6_r[im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_F6_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.8);
       palette->SetX2NDC(0.85);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, Form("F6: (%s-%s)/(%s)", 
      gSFAlgList[im]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str(), gSFAlgList[0]->Id().Key().c_str()) );

    cnv->Update();
  }

  ps->Close();

  delete cnv;
  delete ps;
  delete lgnd;

  string root_filename = gOptOutFile + ".root";
  TFile f(root_filename.c_str(),"recreate");
  ntpl->Write("sfntp");
  f.Close();
  delete ntpl;
}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  CmdLnArgParser parser(argc,argv);

  if(parser.OptionExists("structure-function")){
    gOptSF = parser.Arg("structure-function");
    LOG("gsfcomp", pNOTICE) << "Input structure functions: " << gOptSF; 
  } else {
    LOG("gsfcomp", pFATAL) 
       << "Please specify structure function models sets "
       << "using the --structure-function argument";
    gAbortingInErr = true;
    exit(1);
  }

  if(parser.OptionExists('o')){
    gOptOutFile = parser.Arg('o');
  }

}
//_________________________________________________________________________________
void GetAlgorithms(void)
{
  vector<string> vsfset = str::Split(gOptSF, ",");
  LOG("gsfcomp", pNOTICE) 
      << "Number of input structure functions: " << vsfset.size();
  if(vsfset.size() == 0) {
     LOG("gsfcomp", pFATAL) << "Need at least 1 structure function!";
     gAbortingInErr = true;
     exit(1);
  }
  vector<string>::iterator vsfset_iter = vsfset.begin();
  for( ; vsfset_iter != vsfset.end(); ++vsfset_iter) {
     vector<string> vsf = str::Split(*vsfset_iter, "/"); 
     if(vsf.size() != 2) {
        LOG("gsfcomp", pFATAL) 
           << "Need to specify both a structire function algorithm name and configuration "
           << "as in genie::BYStrucFunc/Default";
        gAbortingInErr = true;
        exit(1);
     } 
     string sf_alg_name = vsf[0];
     string sf_alg_conf = vsf[1];

     AlgFactory * algf = AlgFactory::Instance();
     const DISStructureFuncModelI * sf_alg = 
         dynamic_cast<const DISStructureFuncModelI *> (
                algf->GetAlgorithm(sf_alg_name, sf_alg_conf));

     if(!sf_alg) {
        LOG("gsfcomp", pFATAL) 
          << "Couldn't instantiate " << sf_alg_name << "/" << sf_alg_conf;
        gAbortingInErr = true;
        exit(1);
     } 
     LOG("gsfcomp", pNOTICE) 
       << "\n Instantiated: " << sf_alg->Id()
       << " with the following configuration: " 
       << sf_alg->GetConfig();

     gSFAlgList.push_back(sf_alg);
 }
}
//_________________________________________________________________________________

