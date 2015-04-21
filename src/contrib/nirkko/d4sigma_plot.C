// *********************************************************************
//  Plot differential cross-sections - Martti Nirkko (28th Nov 2014)
//  Compile and run in terminal:     root -l -b -q d4sigma_plot.C+g
// *********************************************************************

#include <TMath.h>
#include <TH3D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraph.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>

// Compile using:   root -l -b -q d4sigma_plot.C+g
void d4sigma_plot() {
  
  // Should be the same as when input rootfile was generated
  const int COMP = 10;
  
  // Initialisation of variables
  const double pi = TMath::Pi();
  const int nsteps1 = 10*COMP, nsteps2 = 10*COMP, nsteps3 = 10*COMP, nsteps4 = 2*COMP;
  double varmin1, varmin2, varmin3, varmin4, varmax1, varmax2, varmax3, varmax4;
  int i,j,k,l;
  
  double Enu;                           // neutrino energy [GeV]
  printf("Please enter neutrino energy: ");
  scanf("%lf", &Enu);
  printf("Trying to find input file for E_nu = %3.1lf GeV...\n", Enu);
  
  std::string fname = Form("data/d4sigma_hist_%3.1lfGeV.root", Enu);
  TFile *file = new TFile(fname.c_str());
  
  if (!file) {
    printf("ERROR: File not found!");
    return;
  }
  
  TH3D *hist[nsteps4];
  for (l=0; l<nsteps4; l++) {
    hist[l] = (TH3D*)file->Get(Form("d4sigma_hist_%d",l));
  }
  
  // Get minimum/maximum of variables
  varmin1 = hist[0]->GetXaxis()->GetXmin(); varmax1 = hist[0]->GetXaxis()->GetXmax();
  varmin2 = hist[0]->GetYaxis()->GetXmin(); varmax2 = hist[0]->GetYaxis()->GetXmax();
  varmin3 = hist[0]->GetZaxis()->GetXmin(); varmax3 = hist[0]->GetZaxis()->GetXmax();
  varmin4 = -pi; varmax4 = pi;
  
  // Get bin edges for each variable
  double binedge1[nsteps1+1], binedge2[nsteps2+1], binedge3[nsteps3+1];
  binedge1[0] = varmin1; binedge2[0] = varmin2; binedge3[0] = varmin3;
  for (i=1; i<nsteps1; i++) binedge1[i] = hist[0]->GetXaxis()->GetBinLowEdge(i+1);
  for (j=1; j<nsteps2; j++) binedge2[j] = hist[0]->GetYaxis()->GetBinLowEdge(j+1);
  for (k=1; k<nsteps3; k++) binedge3[k] = hist[0]->GetZaxis()->GetBinLowEdge(k+1);
  binedge1[nsteps1] = varmax1; binedge2[nsteps2] = varmax2; binedge3[nsteps3] = varmax3;
  
  // Get bin widths for each variable
  double binwidth1[nsteps1], binwidth2[nsteps2], binwidth3[nsteps3], binwidth4[nsteps4];
  for (i=0; i<nsteps1; i++) binwidth1[i] = hist[0]->GetXaxis()->GetBinWidth(i+1);
  for (j=0; j<nsteps2; j++) binwidth2[j] = hist[0]->GetYaxis()->GetBinWidth(j+1);
  for (k=0; k<nsteps3; k++) binwidth3[k] = hist[0]->GetZaxis()->GetBinWidth(k+1);
  for (l=0; l<nsteps4; l++) binwidth4[l] = (varmax4-varmin4)/nsteps4;
  
  // Initialise histograms with supposedly clever bins
  TH1D* dTk = new TH1D("dTk", "d#sigma/dT_{kaon}",           nsteps1, binedge1);
  TH1D* dTl = new TH1D("dTl", "d#sigma/dT_{lepton}",         nsteps2, binedge2);
  TH1D* dct = new TH1D("dct", "d#sigma/dcos(#theta_{#nul})", nsteps3, binedge3);
  TH1D* dph = new TH1D("dph", "d#sigma/d#phi_{kq}",          nsteps4, varmin4, varmax4);
  
  // Get differential cross-section over T_kaon
  double diff1Tk, diff2Tk, diff3Tk;
  for (i=0; i<nsteps1; i++) {
    diff1Tk = 0.;
    for (j=0; j<nsteps2; j++) {
      diff2Tk = 0.;
      for (k=0; k<nsteps3; k++) {
        diff3Tk = 0.;
        for (l=0; l<nsteps4; l++) {
          diff3Tk += hist[l]->GetBinContent(i+1,j+1,k+1)*binwidth4[l];
        }
        diff2Tk += diff3Tk*binwidth3[k];
      }
      diff1Tk += diff2Tk*binwidth2[j];
    }
    dTk->SetBinContent(i+1, diff1Tk*binwidth1[i]);
  }
  
  // Get differential cross-section over T_lepton
  double diff1Tl, diff2Tl, diff3Tl;
  for (j=0; j<nsteps2; j++) {
    diff1Tl = 0.;
    for (i=0; i<nsteps1; i++) {
      diff2Tl = 0.;
      for (k=0; k<nsteps3; k++) {
        diff3Tl = 0.;
        for (l=0; l<nsteps4; l++) {
          diff3Tl += hist[l]->GetBinContent(i+1,j+1,k+1)*binwidth4[l];
        }
        diff2Tl += diff3Tl*binwidth3[k];
      }
      diff1Tl += diff2Tl*binwidth1[i];
    }
    dTl->SetBinContent(j+1, diff1Tl*binwidth2[j]);
  }
  
  // Get differential cross-section over cos(theta)
  double diff1ct, diff2ct, diff3ct;
  for (k=0; k<nsteps3; k++) {
    diff1ct = 0.;
    for (i=0; i<nsteps1; i++) {
      diff2ct = 0.;
      for (j=0; j<nsteps2; j++) {
        diff3ct = 0.;
        for (l=0; l<nsteps4; l++) {
          diff3ct += hist[l]->GetBinContent(i+1,j+1,k+1)*binwidth4[l];
        }
        diff2ct += diff3ct*binwidth2[j];
      }
      diff1ct += diff2ct*binwidth1[i];
    }
    dct->SetBinContent(k+1, diff1ct*binwidth3[k]);
  }
  
  // Get differential cross-section over phi_kq
  double diff1ph, diff2ph, diff3ph;
  for (l=0; l<nsteps4; l++) {
    diff1ph = 0.;
    for (i=0; i<nsteps1; i++) {
      diff2ph = 0.;
      for (j=0; j<nsteps2; j++) {
        diff3ph = 0.;
        for (k=0; k<nsteps3; k++) {
          diff3ph += hist[l]->GetBinContent(i+1,j+1,k+1)*binwidth3[k];
        }
        diff2ph += diff3ph*binwidth2[j];
      }
      diff1ph += diff2ph*binwidth1[i];
    }
    dph->SetBinContent(l+1, diff1ph*binwidth4[l]);
  }
  
  double sum1 = dTk->Integral();
  double sum2 = dTl->Integral();
  double sum3 = dct->Integral();
  double sum4 = dph->Integral();
  
  printf("Integrals are:\n %.6e\n %.6e\n %.6e\n %.6e\n",sum1,sum2,sum3,sum4);
  
  dTk->SetMinimum(0);
  dTl->SetMinimum(0);
  dct->SetMinimum(0);
  dph->SetMinimum(0);
  
  gStyle->SetOptStat(0);
  
  std::string outname = Form("d4sigma_plot_%3.1lfGeV", Enu);
  
  dTl->SetXTitle("T_{lepton} [GeV]");
  dTk->SetXTitle("T_{kaon} [GeV]");
  dct->SetXTitle("cos(#theta_{#nul}) [ ]");
  dph->SetXTitle("#phi_{Kq} [rad]");
  
  TCanvas *c1 = new TCanvas("c1", "Differential cross-sections", 800, 600);
  c1->Divide(2,2);
  c1->cd(1); dTl->Draw();
  c1->cd(2); dTk->Draw();
  c1->cd(3); dct->Draw();
  c1->cd(4); dph->Draw();
  c1->Print(("images/"+outname+".png").c_str()); c1->Close();
  
  TFile* outfile = new TFile(("data/"+outname+".root").c_str(), "RECREATE");
  dTl->Write("dsigma_dTlepton");
  dTk->Write("dsigma_dTkaon");
  dct->Write("dsigma_dcostheta");
  dph->Write("dsigma_dphikq");
  outfile->Close();
  std::cout << std::endl << "Output written to file: " << outname << ".root" << std::endl << std::endl;
  
}

