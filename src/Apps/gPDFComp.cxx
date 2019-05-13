//____________________________________________________________________________
/*!

\program gpdfcomp

\brief   PDF comparison tool

\syntax  gpdfcomp --pdf-set pdf_set [-o output]

         --pdf-set :
          Specifies a comma separated list of GENIE PDFs.
          The full algorithm name and configuration should be provided for
          each GENIE PDF as in `genie::model_name/model_config'.
          (unique color and line styles are defined for up to 4 sets)

         -o :
          Specifies a name to be used in the output files.
          Default: pdf_comp

\example gpdfcomp --pdf-set genie::GRV98LO/Default,genie::BYPDF/Default 

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
#include "Physics/PartonDistributions/PDFModelI.h"
//#include "Physics/PartonDistributions/LHAPDF5.h"
#include "Physics/PartonDistributions/PDF.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/Style.h"

using namespace std;
using namespace genie;
using namespace genie::utils;

// globals
string        gOptPDFSet  = "";         // --pdf-set argument
string        gOptOutFile = "pdf_comp"; // -o argument
vector<const PDFModelI *> gPDFAlgList;

// function prototypes
void GetCommandLineArgs (int argc, char ** argv);
void GetAlgorithms (void);
void MakePlots     (void);

//___________________________________________________________________
int main(int argc, char ** argv)
{
  utils::style::SetDefaultStyle();

  GetCommandLineArgs (argc,argv);   // Get command line arguments
  GetAlgorithms();                  // Get requested PDF algorithms
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
  TNtuple * ntpl = new TNtuple("nt","pdfs","i:uv:dv:us:ds:s:g:x:Q2");

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
  hdr.AddText("GENIE parton density function (pdf) comparisons");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText(" ");
  hdr.AddText("Models used:");
  hdr.AddText(" ");
  for(unsigned int im=0; im < gPDFAlgList.size(); im++) {
      const char * label = gPDFAlgList[im]->Id().Key().c_str();
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

  // PDFs = f(Q^2) for selected vales of x (1-D plots)
  TGraph * gr_xuv_Q2  [nm] = { 0 };
  TGraph * gr_xdv_Q2  [nm] = { 0 }; 
  TGraph * gr_xus_Q2  [nm] = { 0 }; 
  TGraph * gr_xds_Q2  [nm] = { 0 }; 
  TGraph * gr_xstr_Q2 [nm] = { 0 }; 
  TGraph * gr_xglu_Q2 [nm] = { 0 }; 

  // PDFs = f(x,Q^2) with fine x,Q2 binning (2-D plots)
  TH2D * h2_xuv [nm] = { 0 };
  TH2D * h2_xdv [nm] = { 0 };
  TH2D * h2_xus [nm] = { 0 };
  TH2D * h2_xds [nm] = { 0 };
  TH2D * h2_xstr[nm] = { 0 };
  TH2D * h2_xglu[nm] = { 0 };

  // PDFs ratios relative to first one specified = f(x,Q^2) with fine x,Q2 binning (2-D plots)
  TH2D * h2_xuv_r [nm] = { 0 };
  TH2D * h2_xdv_r [nm] = { 0 };
  TH2D * h2_xus_r [nm] = { 0 };
  TH2D * h2_xds_r [nm] = { 0 };
  TH2D * h2_xstr_r[nm] = { 0 };
  TH2D * h2_xglu_r[nm] = { 0 };

  //
  // Generate PDF plots as f(Q^2) for selected vales of x (1-D plots)
  //

  for(unsigned int ix=0; ix < nx; ix++) {
    double x = x_arr[ix];
    double xuv_arr  [nm][nQ2];
    double xdv_arr  [nm][nQ2];
    double xus_arr  [nm][nQ2];
    double xds_arr  [nm][nQ2];
    double xstr_arr [nm][nQ2];
    double xglu_arr [nm][nQ2];
    
    double max_gr_xuv_Q2  = -9E9;
    double max_gr_xdv_Q2  = -9E9;
    double max_gr_xus_Q2  = -9E9;
    double max_gr_xds_Q2  = -9E9;
    double max_gr_xstr_Q2 = -9E9;
    double max_gr_xglu_Q2 = -9E9;
    
    for(unsigned int im=0; im < gPDFAlgList.size(); im++) {
      PDF pdf;
      pdf.SetModel(gPDFAlgList[im]);
      for(unsigned int iq2 = 0; iq2 < nQ2; iq2++) {
        double Q2 = Q2_arr[iq2];
        pdf.Calculate(x, Q2);
      //LOG("gpdfcomp", pINFO) << "PDFs:\n" << pdf;
        double xuv  = pdf.UpValence();
        double xdv  = pdf.DownValence();
        double xus  = pdf.UpSea();
        double xds  = pdf.DownSea();
        double xstr = pdf.Strange();
        double xglu = pdf.Gluon();
        xuv_arr  [im][iq2] = x * xuv;
        xdv_arr  [im][iq2] = x * xdv;
        xus_arr  [im][iq2] = x * xus;
        xds_arr  [im][iq2] = x * xds;
        xstr_arr [im][iq2] = x * xstr;
        xglu_arr [im][iq2] = x * xglu;
        ntpl->Fill(im,xuv,xdv,xus,xds,xstr,xglu,x,Q2);
      }//iq2

      gr_xuv_Q2  [im] = new TGraph (nQ2, Q2_arr, xuv_arr  [im]);
      gr_xdv_Q2  [im] = new TGraph (nQ2, Q2_arr, xdv_arr  [im]);
      gr_xus_Q2  [im] = new TGraph (nQ2, Q2_arr, xus_arr  [im]);
      gr_xds_Q2  [im] = new TGraph (nQ2, Q2_arr, xds_arr  [im]);
      gr_xstr_Q2 [im] = new TGraph (nQ2, Q2_arr, xstr_arr [im]);
      gr_xglu_Q2 [im] = new TGraph (nQ2, Q2_arr, xglu_arr [im]);

      genie::utils::style::Format( gr_xuv_Q2  [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);
      genie::utils::style::Format( gr_xdv_Q2  [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);
      genie::utils::style::Format( gr_xus_Q2  [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);
      genie::utils::style::Format( gr_xds_Q2  [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);
      genie::utils::style::Format( gr_xstr_Q2 [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);
      genie::utils::style::Format( gr_xglu_Q2 [im], col[im], sty[im], 2, col[im], mrk[im], msz[im]);

      gr_xuv_Q2  [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");
      gr_xdv_Q2  [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");
      gr_xus_Q2  [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");
      gr_xds_Q2  [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");
      gr_xstr_Q2 [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");
      gr_xglu_Q2 [im] -> GetXaxis() -> SetTitle("Q^{2} (GeV^{2}/c^{2})");

      gr_xuv_Q2  [im] -> GetYaxis() -> SetTitle("x*u_{val}(x,Q^{2})");
      gr_xdv_Q2  [im] -> GetYaxis() -> SetTitle("x*d_{val}(x,Q^{2})");
      gr_xus_Q2  [im] -> GetYaxis() -> SetTitle("x*u_{sea}(x,Q^{2})");
      gr_xds_Q2  [im] -> GetYaxis() -> SetTitle("x*d_{sea}(x,Q^{2})");
      gr_xstr_Q2 [im] -> GetYaxis() -> SetTitle("x*s(x,Q^{2})");
      gr_xglu_Q2 [im] -> GetYaxis() -> SetTitle("x*g(x,Q^{2})");
      
      double this_max_gr_xuv_Q2  = TMath::MaxElement(gr_xuv_Q2 [im]->GetN(),gr_xuv_Q2 [im]->GetY());
      double this_max_gr_xdv_Q2  = TMath::MaxElement(gr_xdv_Q2 [im]->GetN(),gr_xdv_Q2 [im]->GetY());
      double this_max_gr_xus_Q2  = TMath::MaxElement(gr_xus_Q2 [im]->GetN(),gr_xus_Q2 [im]->GetY());
      double this_max_gr_xds_Q2  = TMath::MaxElement(gr_xds_Q2 [im]->GetN(),gr_xds_Q2 [im]->GetY());
      double this_max_gr_xstr_Q2 = TMath::MaxElement(gr_xstr_Q2[im]->GetN(),gr_xstr_Q2[im]->GetY());
      double this_max_gr_xglu_Q2 = TMath::MaxElement(gr_xglu_Q2[im]->GetN(),gr_xglu_Q2[im]->GetY());
      max_gr_xuv_Q2  = std::max(max_gr_xuv_Q2 ,this_max_gr_xuv_Q2 );
      max_gr_xdv_Q2  = std::max(max_gr_xdv_Q2 ,this_max_gr_xdv_Q2 );
      max_gr_xus_Q2  = std::max(max_gr_xus_Q2 ,this_max_gr_xus_Q2 );
      max_gr_xds_Q2  = std::max(max_gr_xds_Q2 ,this_max_gr_xds_Q2 );
      max_gr_xstr_Q2 = std::max(max_gr_xstr_Q2,this_max_gr_xstr_Q2);
      max_gr_xglu_Q2 = std::max(max_gr_xglu_Q2,this_max_gr_xglu_Q2);

    }//im
    
    // Now loop to set sensible limits
    for(unsigned int im=0; im < gPDFAlgList.size(); im++) {
      gr_xuv_Q2  [im] -> SetMinimum(0.);
      gr_xdv_Q2  [im] -> SetMinimum(0.);
      gr_xus_Q2  [im] -> SetMinimum(0.);
      gr_xds_Q2  [im] -> SetMinimum(0.);
      gr_xstr_Q2 [im] -> SetMinimum(0.);
      gr_xglu_Q2 [im] -> SetMinimum(0.);
      gr_xuv_Q2  [im] -> SetMaximum(1.1*max_gr_xuv_Q2  );
      gr_xdv_Q2  [im] -> SetMaximum(1.1*max_gr_xdv_Q2  );
      gr_xus_Q2  [im] -> SetMaximum(1.1*max_gr_xus_Q2  );
      gr_xds_Q2  [im] -> SetMaximum(1.1*max_gr_xds_Q2  );
      gr_xstr_Q2 [im] -> SetMaximum(1.1*max_gr_xstr_Q2 );
      gr_xglu_Q2 [im] -> SetMaximum(1.1*max_gr_xglu_Q2 );
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
    for(unsigned int im=0; im < gPDFAlgList.size(); im++) {
      std::string label(gPDFAlgList[im]->Id().Key());
      lgnd->AddEntry(gr_xuv_Q2 [im], label.c_str(), lgopt[im]);
    }

    cnv->Clear();
    cnv->Divide(2,3);

    cnv->cd(1); gPad->SetLogx();
    cnv->cd(2); gPad->SetLogx();
    cnv->cd(3); gPad->SetLogx();
    cnv->cd(4); gPad->SetLogx();
    cnv->cd(5); gPad->SetLogx();
    cnv->cd(6); gPad->SetLogx();

    for(unsigned int im=0; im < gPDFAlgList.size(); im++) {
      cnv->cd(1); gr_xuv_Q2 [im]->Draw(opt[im]);
      cnv->cd(2); gr_xdv_Q2 [im]->Draw(opt[im]);
      cnv->cd(3); gr_xus_Q2 [im]->Draw(opt[im]);
      cnv->cd(4); gr_xds_Q2 [im]->Draw(opt[im]);
      cnv->cd(5); gr_xstr_Q2[im]->Draw(opt[im]);
      cnv->cd(6); gr_xglu_Q2[im]->Draw(opt[im]);
    }

    cnv->cd(1); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, TString::Format("x = %.3e",x).Data());
    cnv->cd(2); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, TString::Format("x = %.3e",x).Data());
    cnv->cd(3); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, TString::Format("x = %.3e",x).Data());
    cnv->cd(4); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, TString::Format("x = %.3e",x).Data());
    cnv->cd(5); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, TString::Format("x = %.3e",x).Data());
    cnv->cd(6); lgnd->Draw(); tex->DrawTextNDC(0.4, 0.95, TString::Format("x = %.3e",x).Data());

    cnv->Update();
  }

  //
  // Plot PDFs = f(x,Q2)
  //

  for(unsigned int im=0; im < gPDFAlgList.size(); im++) {
    h2_xuv [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_xdv [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_xus [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_xds [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_xstr[im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_xglu[im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    PDF pdf;
    pdf.SetModel(gPDFAlgList[im]);
    for(int ibinx = 1; 
            ibinx <= h2_xuv[im]->GetXaxis()->GetNbins(); ibinx++) {
      double x = h2_xuv[im]->GetXaxis()->GetBinCenter(ibinx);
      for(int ibinq2 = 1; 
              ibinq2 <= h2_xuv[im]->GetYaxis()->GetNbins(); ibinq2++) {
         double Q2 = h2_xuv[im]->GetYaxis()->GetBinCenter(ibinq2);
         pdf.Calculate(x, Q2);
       //LOG("gpdfcomp", pINFO) << "PDFs:\n" << pdf;
         double xuv  = x * pdf.UpValence();
         double xdv  = x * pdf.DownValence();
         double xus  = x * pdf.UpSea();
         double xds  = x * pdf.DownSea();
         double xstr = x * pdf.Strange();
         double xglu = x * pdf.Gluon();
         h2_xuv [im] -> SetBinContent(ibinx, ibinq2, xuv );
         h2_xdv [im] -> SetBinContent(ibinx, ibinq2, xdv ); 
         h2_xus [im] -> SetBinContent(ibinx, ibinq2, xus ); 
         h2_xds [im] -> SetBinContent(ibinx, ibinq2, xds ); 
         h2_xstr[im] -> SetBinContent(ibinx, ibinq2, xstr); 
         h2_xglu[im] -> SetBinContent(ibinx, ibinq2, xglu); 
      }
    }

    ps->NewPage();

    h2_xuv [im] -> GetXaxis() -> SetTitle("x");
    h2_xdv [im] -> GetXaxis() -> SetTitle("x");
    h2_xus [im] -> GetXaxis() -> SetTitle("x");
    h2_xds [im] -> GetXaxis() -> SetTitle("x");
    h2_xstr[im] -> GetXaxis() -> SetTitle("x");
    h2_xglu[im] -> GetXaxis() -> SetTitle("x");

    h2_xuv [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_xdv [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_xus [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_xds [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_xstr[im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_xglu[im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    
    h2_xuv [im] -> SetMinimum(0.);
    h2_xdv [im] -> SetMinimum(0.);
    h2_xus [im] -> SetMinimum(0.);
    h2_xds [im] -> SetMinimum(0.);
    h2_xstr[im] -> SetMinimum(0.);
    h2_xglu[im] -> SetMinimum(0.);
    
    

    cnv->Divide(2,3);

    TPaletteAxis * palette = 0;
    tex->SetTextSize(0.03);

    cnv->cd(1); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xuv [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xuv[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, TString::Format("x*uv: %s", gPDFAlgList[im]->Id().Key().c_str()).Data() );

    cnv->cd(2); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xdv [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xdv[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, TString::Format("x*dv: %s", gPDFAlgList[im]->Id().Key().c_str()).Data() );

    cnv->cd(3); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xus [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xus[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, TString::Format("x*us: %s", gPDFAlgList[im]->Id().Key().c_str()).Data() );

    cnv->cd(4); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xds [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xds[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, TString::Format("x*ds: %s", gPDFAlgList[im]->Id().Key().c_str()).Data() );

    cnv->cd(5); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xstr[im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xstr[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, TString::Format("x*str: %s", gPDFAlgList[im]->Id().Key().c_str()).Data() );

    cnv->cd(6); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xglu[im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xglu[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.30, 0.95, TString::Format("x*glu: %s", gPDFAlgList[im]->Id().Key().c_str()).Data() );

    cnv->Update();
  }

  //
  // For multiple PDFs sets, plot the ratio of each PDF = f(x,Q2) to the first one
  // 

  for(unsigned int im=1; im < gPDFAlgList.size(); im++) {
    h2_xuv_r [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_xdv_r [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_xus_r [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_xds_r [im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_xstr_r[im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    h2_xglu_r[im] = new TH2D("","", nx_2d-1, x_bin_edges_2d, nQ2_2d-1, Q2_bin_edges_2d);
    for(int ibinx = 1; 
            ibinx <= h2_xuv[im]->GetXaxis()->GetNbins(); ibinx++) {
      for(int ibinq2 = 1; 
              ibinq2 <= h2_xuv[im]->GetYaxis()->GetNbins(); ibinq2++) {

         double xuv   = h2_xuv [im] -> GetBinContent(ibinx, ibinq2);
         double xdv   = h2_xdv [im] -> GetBinContent(ibinx, ibinq2);
         double xus   = h2_xus [im] -> GetBinContent(ibinx, ibinq2);
         double xds   = h2_xds [im] -> GetBinContent(ibinx, ibinq2);
         double xstr  = h2_xstr[im] -> GetBinContent(ibinx, ibinq2);
         double xglu  = h2_xglu[im] -> GetBinContent(ibinx, ibinq2);

         double xuv0  = h2_xuv [0]  -> GetBinContent(ibinx, ibinq2);
         double xdv0  = h2_xdv [0]  -> GetBinContent(ibinx, ibinq2);
         double xus0  = h2_xus [0]  -> GetBinContent(ibinx, ibinq2);
         double xds0  = h2_xds [0]  -> GetBinContent(ibinx, ibinq2);
         double xstr0 = h2_xstr[0]  -> GetBinContent(ibinx, ibinq2);
         double xglu0 = h2_xglu[0]  -> GetBinContent(ibinx, ibinq2);

         double xuv_r  = ((xuv0  > 0.) ? (xuv -xuv0 )/xuv0  : 0.);
         double xdv_r  = ((xdv0  > 0.) ? (xdv -xdv0 )/xdv0  : 0.);
         double xus_r  = ((xus0  > 0.) ? (xus -xus0 )/xus0  : 0.);
         double xds_r  = ((xds0  > 0.) ? (xds -xds0 )/xds0  : 0.);
         double xstr_r = ((xstr0 > 0.) ? (xstr-xstr0)/xstr0 : 0.);
         double xglu_r = ((xglu0 > 0.) ? (xglu-xglu0)/xglu0 : 0.);

         h2_xuv_r [im] -> SetBinContent(ibinx, ibinq2, xuv_r );
         h2_xdv_r [im] -> SetBinContent(ibinx, ibinq2, xdv_r ); 
         h2_xus_r [im] -> SetBinContent(ibinx, ibinq2, xus_r ); 
         h2_xds_r [im] -> SetBinContent(ibinx, ibinq2, xds_r ); 
         h2_xstr_r[im] -> SetBinContent(ibinx, ibinq2, xstr_r); 
         h2_xglu_r[im] -> SetBinContent(ibinx, ibinq2, xglu_r); 
      }
    }

    ps->NewPage();

    h2_xuv_r [im] -> GetXaxis() -> SetTitle("x");
    h2_xdv_r [im] -> GetXaxis() -> SetTitle("x");
    h2_xus_r [im] -> GetXaxis() -> SetTitle("x");
    h2_xds_r [im] -> GetXaxis() -> SetTitle("x");
    h2_xstr_r[im] -> GetXaxis() -> SetTitle("x");
    h2_xglu_r[im] -> GetXaxis() -> SetTitle("x");

    h2_xuv_r [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_xdv_r [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_xus_r [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_xds_r [im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_xstr_r[im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");
    h2_xglu_r[im] -> GetXaxis() -> SetTitle("Q^2 (GeV^2/c^2)");

    cnv->Divide(2,3);

    TPaletteAxis * palette = 0;

    cnv->cd(1); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xuv_r [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xuv_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, TString::Format("uv: (%s-%s)/(%s)", 
      gPDFAlgList[im]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str()).Data() );

    cnv->cd(2); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xdv_r [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xdv_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, TString::Format("dv: (%s-%s)/(%s)", 
      gPDFAlgList[im]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str()).Data() );

    cnv->cd(3); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xus_r [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xus_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, TString::Format("us: (%s-%s)/(%s)", 
      gPDFAlgList[im]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str()).Data() );

    cnv->cd(4); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xds_r [im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xds_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, TString::Format("ds: (%s-%s)/(%s)", 
      gPDFAlgList[im]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str()).Data() );

    cnv->cd(5); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xstr_r[im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xstr_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, TString::Format("str: (%s-%s)/(%s)", 
      gPDFAlgList[im]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str()).Data() );

    cnv->cd(6); gPad->SetLogx(); gPad->SetLogy(); 
    h2_xglu_r[im] -> Draw("colz");
    gPad->Update();
    palette = (TPaletteAxis*)h2_xglu_r[im]->GetListOfFunctions()->FindObject("palette");
    if(palette) {
       palette->SetX1NDC(0.2);
       palette->SetX2NDC(0.25);
       palette->SetY1NDC(0.4);
       palette->SetY2NDC(0.8);
    }
    tex->DrawTextNDC(0.1, 0.95, TString::Format("glu: (%s-%s)/(%s)", 
      gPDFAlgList[im]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str(), gPDFAlgList[0]->Id().Key().c_str()).Data() );

    cnv->Update();
  }

  ps->Close();

  delete cnv;
  delete ps;
  delete lgnd;

  string root_filename = gOptOutFile + ".root";
  TFile f(root_filename.c_str(),"recreate");
  ntpl->Write("pdflib");

  f.Close();
  delete ntpl;
  
  std::cout<<"Done."<<std::endl;
}
//_________________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  CmdLnArgParser parser(argc,argv);

  if(parser.OptionExists("pdf-set")){
    gOptPDFSet = parser.Arg("pdf-set");
    LOG("gpdfcomp", pNOTICE) << "Input PDF set: " << gOptPDFSet; 
  } else {
    LOG("gpdfcomp", pFATAL) 
       << "Please specify PDF sets using the --pdf-set argument";
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
  vector<string> vpdfset = str::Split(gOptPDFSet, ",");
  LOG("gpdfcomp", pNOTICE) 
      << "Number of input PDF sets: " << vpdfset.size();
  if(vpdfset.size() == 0) {
     LOG("gpdfcomp", pFATAL) << "Need at least 1 PDF set!";
     gAbortingInErr = true;
     exit(1);
  }
  vector<string>::iterator vpdfset_iter = vpdfset.begin();
  for( ; vpdfset_iter != vpdfset.end(); ++vpdfset_iter) {
     vector<string> vpdf = str::Split(*vpdfset_iter, "/"); 
     if(vpdf.size() != 2) {
        LOG("gpdfcomp", pFATAL) 
           << "Need to specify both a PDF algorithm name and configuration "
           << "as in genie::GRV98LO/Default";
        gAbortingInErr = true;
        exit(1);
     } 
     string pdf_alg_name = vpdf[0];
     string pdf_alg_conf = vpdf[1];

     AlgFactory * algf = AlgFactory::Instance();
     const PDFModelI * pdf_alg = 
         dynamic_cast<const PDFModelI *> (
                algf->GetAlgorithm(pdf_alg_name, pdf_alg_conf));

     if(!pdf_alg) {
        LOG("gpdfcomp", pFATAL) 
          << "Couldn't instantiate " << pdf_alg_name << "/" << pdf_alg_conf;
        gAbortingInErr = true;
        exit(1);
     } 
     LOG("gpdfcomp", pNOTICE) 
       << "\n Instantiated: " << pdf_alg->Id()
       << " with the following configuration: " 
       << pdf_alg->GetConfig();

     gPDFAlgList.push_back(pdf_alg);
 }
}
//_________________________________________________________________________________

