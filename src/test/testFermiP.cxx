//____________________________________________________________________________
/*!

\program testFermiP

\brief   program used for testing / debugging the Fermi momentum distribution
         models

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created June 20, 2004

\cpright Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1F.h>
#include <TPad.h>
#include <TPaveLabel.h>
#include <TMath.h>
#include <TVector3.h>
#include <TPostScript.h>
#include <TStyle.h>

#include "Algorithm/AlgFactory.h"
#include "Interaction/Target.h"
#include "Messenger/Messenger.h"
#include "Nuclear/NuclearModelI.h"
#include "PDG/PDGCodes.h"
#include "Utils/PrintUtils.h"

using namespace genie;

int main(int /*argc*/, char ** /*argv*/)
{
  AlgFactory * algf = AlgFactory::Instance();
  const NuclearModelI * nuclmodel =
       dynamic_cast<const NuclearModelI *> (algf->GetAlgorithm(
			     "genie::FGMBodekRitchie","Default"));
  //-- Declare targets
  const unsigned int kNTargets = 4;

  Target * nucltgt[kNTargets];

  nucltgt[0] = new Target ( 6,  12); // C12
  nucltgt[1] = new Target ( 8,  16); // O16
  nucltgt[2] = new Target (26,  56); // Fe56
  nucltgt[3] = new Target (82, 208); // Pb208

  //-- Declare histogram arrays
  TH1F * hnucpx[kNTargets];
  TH1F * hnucpy[kNTargets];
  TH1F * hnucpz[kNTargets];
  TH1F * hpprob[kNTargets];

  //-- Set style
  gStyle->SetTitle(0);
  gStyle->SetOptStat(0);

  //-- Create & format output canvas
  TCanvas * c = new TCanvas("c","",20,20,500,500);
  c->SetBorderMode(0);
  c->SetFillColor(0);
  c->Draw();

  TPad * plot_pad = new TPad("plot_pad","",0.05,0.05,0.95,0.95);
  plot_pad -> Draw();
  plot_pad -> SetFillColor(0);
  plot_pad -> SetBorderMode(0);

  for(unsigned int i = 0; i < kNTargets; i++) {

     nucltgt[i]->SetHitNucPdg(kPdgProton);
     const Target & target = *nucltgt[i];

     LOG("Main", pINFO) << "Current nuclear target: " << target;

     //-- create histograms
     hnucpx[i] = new TH1F(Form("px%d",i),"", 40, -1.,1.);
     hnucpy[i] = new TH1F(Form("py%d",i),"", 40, -1.,1.);
     hnucpz[i] = new TH1F(Form("pz%d",i),"", 40, -1.,1.);
     hpprob[i] = new TH1F(Form("P%d",i) ,"", 40,  0.,1.);

     //-- fill in px,py,pz plots
     LOG("Main", pINFO) << "Generating nucleon 4-momenta...";

     for(unsigned int iev = 0; iev < 2000; iev++) {

         nuclmodel->GenerateNucleon(target);
         TVector3 p3 = nuclmodel->Momentum3();
         LOG("Main", pDEBUG)
                      << "Nucleon 4-P = " << utils::print::Vec3AsString(&p3);
         hnucpx[i]->Fill(p3.Px());
         hnucpy[i]->Fill(p3.Py());
         hnucpz[i]->Fill(p3.Pz());
     }

     //-- fill in the probability vs momentum plot
     LOG("Main", pINFO) << "Computing momentrum probability distributions...";

     for(int ibin = 1; ibin < hpprob[i]->GetNbinsX(); ibin++) {

         double p    = double (hpprob[i]->GetBinCenter(ibin));
         double Prob = nuclmodel->Prob(p,-1,target);

         LOG("Main", pDEBUG)
                    << "Nucleon |P| = " << p << " -> Probability = " << Prob;

         hpprob[i]->Fill(p,Prob);
     }
  }//targets


  //-- Format graphs
  LOG("Main", pDEBUG) << "Formatting plots";

  for(unsigned int i = 0; i < kNTargets; i++) {
     hnucpx[i]->SetLineWidth(2);
     hnucpy[i]->SetLineWidth(2);
     hnucpz[i]->SetLineWidth(2);
     hpprob[i]->SetLineWidth(2);

     hnucpx[i]->SetLineColor(i+1);
     hnucpy[i]->SetLineColor(i+1);
     hnucpz[i]->SetLineColor(i+1);
     hpprob[i]->SetLineColor(i+1);
  }

  //-- Divide the plot_pad
  LOG("Main", pDEBUG) << "Dividing plot pad";
  plot_pad->cd();
  plot_pad->Divide(2,2);

  //--- Plot Probabilities

  //plot
  plot_pad->cd(1);
  for(unsigned int i = 0; i < kNTargets; i++) {
    if(i==0) hpprob[i]->Draw();
    else     hpprob[i]->Draw("SAME");
  }
  // add legend
  TLegend * legendP = new TLegend(0.6,0.81,0.9,0.9);
  for(unsigned int i = 0; i < kNTargets; i++) {
    legendP->AddEntry(hpprob[i], nucltgt[i]->AsString().c_str(), "L");
  }
  legendP->SetFillColor(0);
  legendP->Draw();
  // format axis
  hpprob[0]->GetXaxis()->SetTitle("|p| (GeV)");
  hpprob[0]->GetYaxis()->SetTitle("Probability");
  hpprob[0]->GetXaxis()->SetTitleOffset(1.3);
  hpprob[0]->GetYaxis()->SetTitleOffset(1.6);

  //--- Plot Px

  //plot
  plot_pad->cd(2);
  for(unsigned int i = 0; i < kNTargets; i++) {
    if(i==0) hnucpx[i]->Draw();
    else     hnucpx[i]->Draw("SAME");
  }
  // add legend
  TLegend * legendPx = new TLegend(0.6,0.81,0.9,0.9);
  for(unsigned int i = 0; i < kNTargets; i++) {
    legendPx->AddEntry(hnucpx[i], nucltgt[i]->AsString().c_str(), "L");
  }
  legendPx->Draw();
  legendPx->SetFillColor(0);
  // format axis
  hnucpx[0]->GetXaxis()->SetTitle("px (GeV)");
  hnucpx[0]->GetXaxis()->SetTitleOffset(1.3);

  //--- Plot Py

  //plot
  plot_pad->cd(3);
  for(unsigned int i = 0; i < kNTargets; i++) {
    if(i==0) hnucpy[i]->Draw();
    else     hnucpy[i]->Draw("SAME");
  }
  // add legend
  TLegend * legendPy = new TLegend(0.6,0.81,0.9,0.9);
  for(unsigned int i = 0; i < kNTargets; i++) {
    legendPy->AddEntry(hnucpy[i], nucltgt[i]->AsString().c_str(), "L");
  }
  legendPy->SetFillColor(0);
  legendPy->Draw();
  // format axis
  hnucpy[0]->GetXaxis()->SetTitle("py (GeV)");
  hnucpy[0]->GetXaxis()->SetTitleOffset(1.3);

  //--- Plot Pz

  //plot
  plot_pad->cd(4);
  for(unsigned int i = 0; i < kNTargets; i++) {
    if(i==0) hnucpz[i]->Draw();
    else     hnucpz[i]->Draw("SAME");
  }
  // add legend
  TLegend * legendPz = new TLegend(0.6,0.81,0.9,0.9);
  for(unsigned int i = 0; i < kNTargets; i++) {
    legendPz->AddEntry(hnucpz[i], nucltgt[i]->AsString().c_str(), "L");
  }
  legendPz->SetFillColor(0);
  legendPz->Draw();
  // format axis
  hnucpz[0]->GetXaxis()->SetTitle("pz (GeV)");
  hnucpz[0]->GetXaxis()->SetTitleOffset(1.3);

  c->Update();

  TFile f("./genie-fermi-motion.root","recreate");
  c->Write("plot");
  hpprob[0]->Write("C12p");
  hpprob[1]->Write("O16p");
  hpprob[2]->Write("Fe56p");
  hpprob[3]->Write("Pb208p");
  f.Close();

  for(unsigned int i=0; i<kNTargets; i++) {
    delete nucltgt[i];
    delete hnucpx[i];
    delete hnucpy[i];
    delete hnucpz[i];
    delete hpprob[i];
  }
  delete legendP;
  delete legendPx;
  delete legendPy;
  delete legendPz;
  delete plot_pad;
  delete c;
}
//____________________________________________________________________________
