//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Tingjun Yang (Stanford Univ.)

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Mar 05, 2009 - TY
   Was first added in v2.5.1.

*/
//____________________________________________________________________________
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TProfile.h"

#include "ValidationTools/Hadronization/HadPlotter.h"
#include "Riostream.h"
#include <vector>
#include <iostream>

using namespace genie;
using namespace genie::vld;
using namespace std;

const int mccolors[7] = {2,4,3,1,6,7,5};

TH2D* DrawFrame(TCanvas *c1, int dir, double fr_xmin, double fr_xmax, double fr_ymin, double fr_ymax, string frtit, string frxtit, string frytit, bool logx, bool logy);

TGraphErrors* MakeGraph(string file);
void DrawGraph(TGraph* gr, int markerstyle, int markercolor, double markersize=0.8, string opt = "def");

void SetLegend(TLegend *leg){
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
}

//____________________________________________________________________________
void HadPlotter::AddPlots(HadPlots hp)
{
  hadPlots.push_back(hp);
}
//____________________________________________________________________________
void HadPlotter::ShowPlots()
{
  
  TH2D *htemp;
  
  //................charged hadron multiplicity...........................
  TGraphErrors *prd27_47_f4_p = MakeGraph("prd27_47_fig4_vp.dat");
  TGraphErrors *prd27_47_f4_n = MakeGraph("prd27_47_fig4_vn.dat");
  TGraphErrors *npb181_385_fig1_nch_W = MakeGraph("npb181_385_fig1_nch_W.dat");
  TGraphErrors *zpc24_119_vn = MakeGraph("zpc24_119_vn.dat");
  TCanvas *cMulCh = new TCanvas("cMulCh","cMulCh",800,350);
  cMulCh->Divide(2,1);
  DrawFrame(cMulCh,1,1,300,0,10,"","W^{2}(GeV^{2}/c^{4})","<n_{ch}>",true,false);
  DrawGraph(prd27_47_f4_p,20,1,0.8);
  DrawGraph(npb181_385_fig1_nch_W,24,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].nch_w[0],24,mccolors[i],0.8,"c");
  }
  TLegend *lMulCh1 = new TLegend(0.2,0.6,0.5,0.85);
  SetLegend(lMulCh1);
  lMulCh1->AddEntry(prd27_47_f4_p,"15' #nuD_{2} (1983)","p");
  lMulCh1->AddEntry(npb181_385_fig1_nch_W,"BEBC #nuH_{2} (1981)","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lMulCh1->AddEntry(hadPlots[i].nch_w[0],name.c_str(),"l");
  }
  lMulCh1->Draw();
  TLatex *tMulCh1 = new TLatex(0.7,0.3,"#nup#rightarrow#mu^{-}X^{++}");
  tMulCh1->SetNDC();
  tMulCh1->SetTextSize(0.06);
  tMulCh1->Draw();

  DrawFrame(cMulCh,2,1,300,0,10,"","W^{2}(GeV^{2}/c^{4})","<n_{ch}>",true,false);
  DrawGraph(prd27_47_f4_n,20,1,0.8);
  DrawGraph(zpc24_119_vn,24,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].nch_w[1],24,mccolors[i],0.8,"c");
  }
  TLegend *lMulCh2 = new TLegend(0.2,0.6,0.5,0.85);
  SetLegend(lMulCh2);
  lMulCh2->AddEntry(prd27_47_f4_n,"15' #nuD_{2} (1983)","p");
  lMulCh2->AddEntry(zpc24_119_vn,"BEBC #nuD_{2} (1984)","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lMulCh2->AddEntry(hadPlots[i].nch_w[1],name.c_str(),"l");
  }
  lMulCh2->Draw();
  TLatex *tMulCh2 = new TLatex(0.7,0.3,"#nun#rightarrow#mu^{-}X^{+}");
  tMulCh2->SetNDC();
  tMulCh2->Draw();
  tMulCh2->SetTextSize(0.06);

  //..................dispersion (D_)...............................
  TGraphErrors *prd27_47_fig5_D_n_vp = MakeGraph("prd27_47_fig5_D_n_vp.dat");
  TGraphErrors *prd27_47_fig5_D_n_vn = MakeGraph("prd27_47_fig5_D_n_vn.dat");
  TCanvas *cDispCh = new TCanvas("cDispCh","cDispCh",800,350);
  cDispCh->Divide(2,1);
  DrawFrame(cDispCh,1,0,4,0,2,"","<n_{-}>","D_{-}",false,false);
  DrawGraph(prd27_47_fig5_D_n_vp,24,1,0.8);
  DrawGraph(prd27_47_fig5_D_n_vn,20,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].D_nneg[0],24,mccolors[i],0.8,"c");
    hadPlots[i].D_nneg[1]->SetLineStyle(2);
    DrawGraph(hadPlots[i].D_nneg[1],20,mccolors[i],0.8,"c");
  }
  TLegend *lDispCh1 = new TLegend(0.5,0.25,0.9,0.45);
  SetLegend(lDispCh1);
  lDispCh1->AddEntry(prd27_47_fig5_D_n_vp,"#nup 15' #nuD_{2} (1983)","p");
  lDispCh1->AddEntry(prd27_47_fig5_D_n_vn,"#nun 15' #nuD_{2} (1983)","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = "#nup "+hadPlots[i].modelName;
    lDispCh1->AddEntry(hadPlots[i].D_nneg[0],name.c_str(),"l");
    name = "#nun "+hadPlots[i].modelName;
    lDispCh1->AddEntry(hadPlots[i].D_nneg[1],name.c_str(),"l");
  }
  lDispCh1->Draw();
  //.....................D/<n_{ch}>..................................
  TGraphErrors *prd27_47_fig5_D_W_vp = MakeGraph("prd27_47_fig5_D_W_vp.dat");
  TGraphErrors *prd27_47_fig5_D_W_vn = MakeGraph("prd27_47_fig5_D_W_vn.dat");
  DrawFrame(cDispCh,2,1,1000,0,0.8,"","W^{2}(GeV^{2}/c^{4})","D/<n_{ch}>",true,false);
  DrawGraph(prd27_47_fig5_D_W_vp,24,1,0.8);
  DrawGraph(prd27_47_fig5_D_W_vn,20,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].D_W2[0],24,mccolors[i],0.8,"c");
    hadPlots[i].D_W2[1]->SetLineStyle(2);
    DrawGraph(hadPlots[i].D_W2[1],20,mccolors[i],0.8,"c");
  }
  TLegend *lDispCh2 = new TLegend(0.5,0.25,0.9,0.45);
  SetLegend(lDispCh2);
  lDispCh2->AddEntry(prd27_47_fig5_D_W_vp,"#nup 15' #nuD_{2} (1983","p");
  lDispCh2->AddEntry(prd27_47_fig5_D_W_vn,"#nun 15' #nuD_{2} (1983)","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = "#nup "+hadPlots[i].modelName;
    lDispCh2->AddEntry(hadPlots[i].D_W2[0],name.c_str(),"l");
    name = "#nun "+hadPlots[i].modelName;
    lDispCh2->AddEntry(hadPlots[i].D_W2[1],name.c_str(),"l");
  }
  lDispCh2->Draw();

  //.....................KNO distributions..........................
  for (unsigned i = 0; i<hadPlots.size(); i++){
    for (int j = 0; j<2; j++){
      for (int k = 0; k<5; k++){
	hadPlots[i].kno[k][j]->SetMarkerColor(2);
	hadPlots[i].kno[k][j]->SetMarkerSize(0.8);
	hadPlots[i].kno[k][j]->SetLineColor(2);
      }
      hadPlots[i].kno[0][j]->SetMarkerStyle(20);
      hadPlots[i].kno[1][j]->SetMarkerStyle(21);
      hadPlots[i].kno[2][j]->SetMarkerStyle(22);
      hadPlots[i].kno[3][j]->SetMarkerStyle(23);
      hadPlots[i].kno[4][j]->SetMarkerStyle(29);
    }
  }
  TGraphErrors *vn_1_3 = MakeGraph("prd27_47_fig6_vn_1_3.dat");
  TGraphErrors *vn_3_5 = MakeGraph("prd27_47_fig6_vn_3_5.dat");
  TGraphErrors *vn_5_7 = MakeGraph("prd27_47_fig6_vn_5_7.dat");
  TGraphErrors *vn_7_10 = MakeGraph("prd27_47_fig6_vn_7_10.dat");
  TGraphErrors *vn_10_15 = MakeGraph("prd27_47_fig6_vn_10_15.dat");
  TGraphErrors *vn_kno = MakeGraph("prd27_47_fig6_vn.dat");
  
  vn_1_3->SetMarkerStyle(24);
  vn_3_5->SetMarkerStyle(25);
  vn_5_7->SetMarkerStyle(26);
  vn_7_10->SetMarkerStyle(27);
  vn_10_15->SetMarkerStyle(30);

  vn_1_3->SetMarkerSize(0.8);
  vn_3_5->SetMarkerSize(0.8);
  vn_5_7->SetMarkerSize(0.8);
  vn_7_10->SetMarkerSize(0.8);
  vn_10_15->SetMarkerSize(0.8);

  TGraphErrors *vp_1_3 = MakeGraph("prd27_47_fig6_vp_1_3.dat");
  TGraphErrors *vp_3_5 = MakeGraph("prd27_47_fig6_vp_3_5.dat");
  TGraphErrors *vp_5_7 = MakeGraph("prd27_47_fig6_vp_5_7.dat");
  TGraphErrors *vp_7_10 = MakeGraph("prd27_47_fig6_vp_7_10.dat");
  TGraphErrors *vp_10_15 = MakeGraph("prd27_47_fig6_vp_10_15.dat");
  TGraphErrors *vp_kno = MakeGraph("prd27_47_fig6_vp.dat");

  vp_1_3->SetMarkerStyle(24);
  vp_3_5->SetMarkerStyle(25);
  vp_5_7->SetMarkerStyle(26);
  vp_7_10->SetMarkerStyle(27);
  vp_10_15->SetMarkerStyle(30);

  vp_1_3->SetMarkerSize(0.8);
  vp_3_5->SetMarkerSize(0.8);
  vp_5_7->SetMarkerSize(0.8);
  vp_7_10->SetMarkerSize(0.8);
  vp_10_15->SetMarkerSize(0.8);

  TF1 *levy = new TF1("fLevy","2*TMath::Exp(-[0])*pow([0],[0]*x+1)/TMath::Gamma([0]*x+1)",0,3);

  TCanvas *cKNO = new TCanvas("cKNO","cKNO",800,400);
  cKNO->Divide(2,1);
  DrawFrame(cKNO,1,0,3,0.01,20,"","n_{ch}/<n_{ch}>","<n_{ch}>P(n_{ch})",false,true);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    for (int j = 0; j<5; j++){
      hadPlots[i].kno[j][0]->Draw("e same");
    }
  }
  vp_1_3->Draw("p");
  vp_3_5->Draw("p");
  vp_5_7->Draw("p");
  vp_7_10->Draw("p");
  vp_10_15->Draw("p");
  //KNO->Draw("p");
  levy->SetParameter(0,8.67);
  vp_kno->Fit("fLevy","Q");
  levy->SetLineWidth(1);
  levy->SetLineColor(1);
  levy->DrawCopy("same");
  std::cout<<"#nup c="<<levy->GetParameter(0)<<"+/-"<<levy->GetParError(0)<<endl;

  TLegend *lKNO1 = new TLegend(0.65,0.55,0.95,0.9);
  SetLegend(lKNO1);
  lKNO1->SetHeader("Data 15' #nuD_{2}");
  lKNO1->AddEntry(vn_1_3,"1<W<3GeV","p");
  lKNO1->AddEntry(vn_3_5,"3<W<5GeV","p");
  lKNO1->AddEntry(vn_5_7,"5<W<7GeV","p");
  lKNO1->AddEntry(vn_7_10,"7<W<10GeV","p");
  lKNO1->AddEntry(vn_10_15,"10<W<15GeV","p");
  lKNO1->SetBorderSize(0);
  lKNO1->Draw("same");
  TLatex *tKNO1 = new TLatex(0.8,6,"#nup");
  tKNO1->Draw();

  DrawFrame(cKNO,2,0,3,0.01,20,"","n_{ch}/<n_{ch}>","<n_{ch}>P(n_{ch})",false,true);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    for (int j = 0; j<5; j++){
      hadPlots[i].kno[j][1]->Draw("e same");
    }
  }
  vn_1_3->Draw("p");
  vn_3_5->Draw("p");
  vn_5_7->Draw("p");
  vn_7_10->Draw("p");
  vn_10_15->Draw("p");
  vn_kno->Fit("fLevy","Q");
  levy->DrawCopy("same");
  std::cout<<"#nun c="<<levy->GetParameter(0)<<"+/-"<<levy->GetParError(0)<<endl;
  TLegend *lKNO2 = new TLegend(0.65,0.55,0.95,0.9);
  SetLegend(lKNO2);
  lKNO2->SetHeader(hadPlots[0].modelName.c_str());
  lKNO2->SetTextColor(2);
  lKNO2->SetLineColor(2);
  lKNO2->AddEntry(hadPlots[0].kno[0][0],"W<2GeV","p");
  lKNO2->AddEntry(hadPlots[0].kno[1][0],"2<=W<3GeV","p");
  lKNO2->AddEntry(hadPlots[0].kno[2][0],"3<=W<4GeV","p");
  lKNO2->AddEntry(hadPlots[0].kno[3][0],"4<=W<5GeV","p");
  lKNO2->AddEntry(hadPlots[0].kno[4][0],"W>=5GeV","p");
  lKNO2->SetBorderSize(0);
  lKNO2->Draw("same");

  TLatex *tKNO2 = new TLatex(0.8,6,"#nun");
  tKNO2->Draw();

  //.......................pi0 multiplicity..............................
  TGraphErrors *sn41_963_fig2 = MakeGraph("sn41_963_fig2_vA.dat");
  TGraphErrors *npb223_269_fig8 = MakeGraph("npb223_269_fig8_vp.dat");
  TGraphErrors *zpc40_231_fig13_pi0 = MakeGraph("zpc40_231_fig13_pi0.dat");
  TCanvas *cMulPi0 = new TCanvas("cMulPi0","cMulPi0",800,350);
  cMulPi0->Divide(2,1);
  DrawFrame(cMulPi0,1,1,1000,0,4,"","W^{2}(GeV^{2}/c^{4})","<n_{#pi^{0}}>",true,false);
  DrawGraph(sn41_963_fig2,20,1,0.8);
  DrawGraph(npb223_269_fig8,26,1,0.8);
  DrawGraph(zpc40_231_fig13_pi0,27,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].npi0_w[0],24,mccolors[i],0.8,"c");
    hadPlots[i].npi0_w[1]->SetLineStyle(2);
    DrawGraph(hadPlots[i].npi0_w[1],20,mccolors[i],0.8,"c");
  }
  TLegend *lMulPi0_1 = new TLegend(0.6,0.25,0.95,0.55);
  SetLegend(lMulPi0_1);
  lMulPi0_1->AddEntry(sn41_963_fig2,"#nuA SKAT #nuFreon","p");
  lMulPi0_1->AddEntry(npb223_269_fig8,"#nup BEBC #nuH_{2}","p");
  lMulPi0_1->AddEntry(zpc40_231_fig13_pi0,"#nuneon BEBC","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = "#nup "+hadPlots[i].modelName;
    lMulPi0_1->AddEntry(hadPlots[i].npi0_w[0],name.c_str(),"l");
    name = "#nun "+hadPlots[i].modelName;
    lMulPi0_1->AddEntry(hadPlots[i].npi0_w[1],name.c_str(),"l");
  }
  lMulPi0_1->Draw();
  //..................dispersion (D_pi0)...............................
  TGraphErrors *sn41_963_fig4_vA = MakeGraph("sn41_963_fig4_vA.dat");
  DrawFrame(cMulPi0,2,0,3,0,2,"","<n_{#pi^{0}}>","D_{#pi^{0}}",false,false);
  DrawGraph(sn41_963_fig4_vA,20,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].D_npi0[0],24,mccolors[i],0.8,"l");
    hadPlots[i].D_npi0[1]->SetLineStyle(2);
    DrawGraph(hadPlots[i].D_npi0[1],20,mccolors[i],0.8,"l");
  }
  TLegend *lMulPi0_2 = new TLegend(0.6,0.25,0.95,0.6);
  SetLegend(lMulPi0_2);
  lMulPi0_2->AddEntry(sn41_963_fig4_vA,"#nuA SKAT #nuFreon","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = "#nup "+hadPlots[i].modelName;
    lMulPi0_2->AddEntry(hadPlots[i].D_npi0[0],name.c_str(),"l");
    name = "#nun "+hadPlots[i].modelName;
    lMulPi0_2->AddEntry(hadPlots[i].D_npi0[1],name.c_str(),"l");
  }
  lMulPi0_2->Draw();

  //...................<npi0> <nch-> correlation.....................
  TGraphErrors *npb223_269_fig9_vp_3_4 = MakeGraph("npb223_269_fig9_vp_3_4.dat");
  TGraphErrors *npb223_269_fig9_vp_4_5 = MakeGraph("npb223_269_fig9_vp_4_5.dat");
  TGraphErrors *npb223_269_fig9_vp_5_7 = MakeGraph("npb223_269_fig9_vp_5_7.dat");
  TGraphErrors *npb223_269_fig9_vp_7_10 = MakeGraph("npb223_269_fig9_vp_7_10.dat");
  TGraphErrors *npb223_269_fig9_vp_10_14 = MakeGraph("npb223_269_fig9_vp_10_14.dat");

  npb223_269_fig9_vp_3_4->SetMarkerStyle(24);
  npb223_269_fig9_vp_4_5->SetMarkerStyle(24);
  npb223_269_fig9_vp_5_7->SetMarkerStyle(24);
  npb223_269_fig9_vp_7_10->SetMarkerStyle(24);
  npb223_269_fig9_vp_10_14->SetMarkerStyle(24);

  TCanvas *cCorr_Pi0_Ch = new TCanvas("cCorr_Pi0_Ch","cCorr_Pi0_Ch",800,500);
  cCorr_Pi0_Ch->Divide(2,2,0,0);
  cCorr_Pi0_Ch->cd(1);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0);
  //gPad->SetTickx(2);
  htemp = DrawFrame(cCorr_Pi0_Ch,1,-0.2,6.5,-0.2,4.8,"","n_{-}","<n_{#pi^{0}}>",false,false);
  htemp->GetYaxis()->SetTitleSize(0.07);
  htemp->GetYaxis()->SetTitleOffset(0.8);
  DrawGraph(npb223_269_fig9_vp_3_4,24,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    TH1D *h0 = hadPlots[i].npi0_nneg[0][0]->ProjectionX(Form("h0_%d",i));
    TGraph *gr0 = new TGraph(4);
    for (int j = 0; j<4; j++){
      gr0->SetPoint(j,h0->GetBinCenter(j+2),h0->GetBinContent(j+2));
    }
    gr0->SetLineColor(mccolors[i]);
    gr0->SetLineWidth(2);
    gr0->Draw("c");
  }
  TLegend* lCorr_Pi0_Ch = new TLegend(0.55,0.7,0.85,0.9);
  SetLegend(lCorr_Pi0_Ch);
  lCorr_Pi0_Ch->AddEntry(npb223_269_fig9_vp_3_4,"#nup BEBC #nuH_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = "#nup "+hadPlots[i].modelName;
    hadPlots[i].npi0_nneg[0][0]->SetLineColor(mccolors[i]);
    lCorr_Pi0_Ch->AddEntry(hadPlots[i].npi0_nneg[0][0],name.c_str(),"l");
    lCorr_Pi0_Ch->Draw();
  }
  TLatex *tCorr_Pi0_Ch1 = new TLatex(4,1,"(a) 3<W<4GeV/c^{2}");
  tCorr_Pi0_Ch1->SetTextSize(0.08);
  tCorr_Pi0_Ch1->Draw();

  cCorr_Pi0_Ch->cd(2);
  gPad->SetBottomMargin(0);
  gPad->SetLeftMargin(0);
  //gPad->SetTickx(2);
  gPad->SetTicky(2);
  htemp = DrawFrame(cCorr_Pi0_Ch,2,-0.2,6.5,-0.2,4.8,"","n_{-}","",false,false);
  htemp->GetYaxis()->SetLabelSize(0);
  DrawGraph(npb223_269_fig9_vp_4_5,24,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    TH1D *h1 = hadPlots[i].npi0_nneg[1][0]->ProjectionX(Form("h1_%d",i));
    TGraph *gr1 = new TGraph(4);
    for (int j = 0; j<4; j++){
      gr1->SetPoint(j,h1->GetBinCenter(j+2),h1->GetBinContent(j+2));
    }
    gr1->SetLineColor(mccolors[i]);
    gr1->SetLineWidth(2);
    gr1->Draw("c");
  }
  lCorr_Pi0_Ch->Draw();
  TLatex *tCorr_Pi0_Ch2 = new TLatex(4,1,"(b) 4<W<5GeV/c^{2}");
  tCorr_Pi0_Ch2->SetTextSize(0.08);
  tCorr_Pi0_Ch2->Draw();
  
  cCorr_Pi0_Ch->cd(3);
  gPad->SetTopMargin(0);
  gPad->SetRightMargin(0);
  DrawFrame(cCorr_Pi0_Ch,3,-0.2,6.5,-0.2,4.8,"","n_{-}","<n_{#pi^{0}}>",false,false);
  DrawGraph(npb223_269_fig9_vp_5_7,24,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    TH1D *h2 = hadPlots[i].npi0_nneg[2][0]->ProjectionX(Form("h2_%d",i));
    TGraph *gr2 = new TGraph(4);
    for (int j = 0; j<5; j++){
      gr2->SetPoint(j,h2->GetBinCenter(j+2),h2->GetBinContent(j+2));
    }
    gr2->SetLineColor(mccolors[i]);
    gr2->SetLineWidth(2);
    gr2->Draw("c");
  }
  lCorr_Pi0_Ch->Draw();
  TLatex *tCorr_Pi0_Ch3 = new TLatex(2,0.5,"(c) 5<W<7GeV/c^{2}");
  tCorr_Pi0_Ch3->SetTextSize(0.06);
  tCorr_Pi0_Ch3->Draw();

  cCorr_Pi0_Ch->cd(4);
  gPad->SetTopMargin(0);
  gPad->SetLeftMargin(0);
  htemp = DrawFrame(cCorr_Pi0_Ch,4,-0.2,6.5,-0.2,4.8,"","n_{-}","",false,false);
  htemp->GetYaxis()->SetLabelSize(0);
  DrawGraph(npb223_269_fig9_vp_7_10,24,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    TH1D *h3 = hadPlots[i].npi0_nneg[3][0]->ProjectionX(Form("h3_%d",i));
    TGraph *gr3 = new TGraph(4);
    for (int j = 0; j<6; j++){
      gr3->SetPoint(j,h3->GetBinCenter(j+2),h3->GetBinContent(j+2));
    }
    gr3->SetLineColor(mccolors[i]);
    gr3->SetLineWidth(2);
    gr3->Draw("c");
  }
  lCorr_Pi0_Ch->Draw();
  TLatex *tCorr_Pi0_Ch4 = new TLatex(2,0.5,"(d) 7<W<10GeV");
  tCorr_Pi0_Ch4->SetTextSize(0.06);
  tCorr_Pi0_Ch4->Draw();

//...........normalised topological cross sections..................
  TGraphErrors *prd27_47_fig3_vp_2 = MakeGraph("prd27_47_fig3_vp_2.dat");
  TGraphErrors *prd27_47_fig3_vp_4 = MakeGraph("prd27_47_fig3_vp_4.dat");
  TGraphErrors *prd27_47_fig3_vp_6 = MakeGraph("prd27_47_fig3_vp_6.dat");
  TGraphErrors *prd27_47_fig3_vp_8 = MakeGraph("prd27_47_fig3_vp_8.dat");
  TGraphErrors *prd27_47_fig3_vp_10 = MakeGraph("prd27_47_fig3_vp_10.dat");
  TGraphErrors *prd27_47_fig3_vp_12 = MakeGraph("prd27_47_fig3_vp_12.dat");

  TGraphErrors *prd27_47_fig3_vn_1 = MakeGraph("prd27_47_fig3_vn_1.dat");
  TGraphErrors *prd27_47_fig3_vn_3 = MakeGraph("prd27_47_fig3_vn_3.dat");
  TGraphErrors *prd27_47_fig3_vn_5 = MakeGraph("prd27_47_fig3_vn_5.dat");
  TGraphErrors *prd27_47_fig3_vn_7 = MakeGraph("prd27_47_fig3_vn_7.dat");
  TGraphErrors *prd27_47_fig3_vn_9 = MakeGraph("prd27_47_fig3_vn_9.dat");
  TGraphErrors *prd27_47_fig3_vn_11 = MakeGraph("prd27_47_fig3_vn_11.dat");

  TCanvas *cTopo = new TCanvas("cTopo","cTopo",900,400);
  cTopo->Divide(2,1);
  DrawFrame(cTopo,1,1,1000,1e-4,10,"#nup","W^{2}(GeV^{2}/c^{4})","P(n_{ch})",true,true);
  DrawGraph(prd27_47_fig3_vp_2,20,1,0.7);
  DrawGraph(prd27_47_fig3_vp_4,20,2,0.7);
  DrawGraph(prd27_47_fig3_vp_6,20,3,0.7);
  DrawGraph(prd27_47_fig3_vp_8,20,4,0.7);
  DrawGraph(prd27_47_fig3_vp_10,20,7,0.7);
  DrawGraph(prd27_47_fig3_vp_12,20,6,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pnch_w2[2][0]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[4][0]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[6][0]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[8][0]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[10][0]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[12][0]->SetLineStyle(i+1);
    DrawGraph(hadPlots[i].pnch_w2[2][0],20,1,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[4][0],20,2,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[6][0],20,3,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[8][0],20,4,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[10][0],20,7,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[12][0],20,6,0.8,"c");
  }
  TLegend *lTopo1 = new TLegend(0.78,0.45,1,0.85);
  SetLegend(lTopo1);
  lTopo1->SetHeader("#nup 15' #nuD_{2}");
  lTopo1->AddEntry(prd27_47_fig3_vp_2,"n = 2","p");
  lTopo1->AddEntry(prd27_47_fig3_vp_4,"n = 4","p");
  lTopo1->AddEntry(prd27_47_fig3_vp_6,"n = 6","p");
  lTopo1->AddEntry(prd27_47_fig3_vp_8,"n = 8","p");
  lTopo1->AddEntry(prd27_47_fig3_vp_10,"n = 10","p");
  lTopo1->AddEntry(prd27_47_fig3_vp_12,"n = 12","p");
  lTopo1->Draw();

  DrawFrame(cTopo,2,1,1000,1e-4,10,"#nup","W^{2}(GeV^{2}/c^{4})","P(n_{ch})",true,true);
  DrawGraph(prd27_47_fig3_vn_1,20,1,0.7);
  DrawGraph(prd27_47_fig3_vn_3,20,2,0.7);
  DrawGraph(prd27_47_fig3_vn_5,20,3,0.7);
  DrawGraph(prd27_47_fig3_vn_7,20,4,0.7);
  DrawGraph(prd27_47_fig3_vn_9,20,7,0.7);
  DrawGraph(prd27_47_fig3_vn_11,20,6,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pnch_w2[1][1]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[3][1]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[5][1]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[7][1]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[9][1]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[11][1]->SetLineStyle(i+1);
    DrawGraph(hadPlots[i].pnch_w2[1][1],20,1,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[3][1],20,2,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[5][1],20,3,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[7][1],20,4,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[9][1],20,7,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[11][1],20,6,0.8,"c");
  }
  TLegend *lTopo2 = new TLegend(0.78,0.45,1,0.85);
  SetLegend(lTopo2);
  lTopo2->SetHeader("#nun 15' #nuD_{2}");
  lTopo2->AddEntry(prd27_47_fig3_vn_1,"n = 1","p");
  lTopo2->AddEntry(prd27_47_fig3_vn_3,"n = 3","p");
  lTopo2->AddEntry(prd27_47_fig3_vn_5,"n = 5","p");
  lTopo2->AddEntry(prd27_47_fig3_vn_7,"n = 7","p");
  lTopo2->AddEntry(prd27_47_fig3_vn_9,"n = 9","p");
  lTopo2->AddEntry(prd27_47_fig3_vn_11,"n = 11","p");
  lTopo2->Draw();

  //.................F/B multiplicity.........................
  TGraphErrors *prd27_47_fig7_vp_tar = MakeGraph("prd27_47_fig7_vp_tar.dat");
  TGraphErrors *prd27_47_fig7_vn_tar = MakeGraph("prd27_47_fig7_vn_tar.dat");
  TGraphErrors *prd27_47_fig7_vp_cur = MakeGraph("prd27_47_fig7_vp_cur.dat");
  TGraphErrors *prd27_47_fig7_vn_cur = MakeGraph("prd27_47_fig7_vn_cur.dat");
  TGraphErrors *npb223_269_nF_W = MakeGraph("npb223_269_nF_W.dat");
  TGraphErrors *npb223_269_nB_W = MakeGraph("npb223_269_nB_W.dat");
  TGraphErrors *zpc24_119_fig4_vp_F = MakeGraph("zpc24_119_fig4_vp_F.dat");
  TGraphErrors *zpc24_119_fig4_vp_B = MakeGraph("zpc24_119_fig4_vp_B.dat");
  TGraphErrors *zpc24_119_fig4_vn_F = MakeGraph("zpc24_119_fig4_vn_F.dat");
  TGraphErrors *zpc24_119_fig4_vn_B = MakeGraph("zpc24_119_fig4_vn_B.dat");

  TCanvas *cMulFB = new TCanvas("cMulFB","cMulFB",800,500);
  cMulFB->Divide(2,2,0,0);
  cMulFB->cd(1);
  gPad->SetBottomMargin(0);
  gPad->SetRightMargin(0);
  //gPad->SetTickx(2);
  htemp = DrawFrame(cMulFB, 1, 1, 300, -0.2, 5.5, "","W^{2}(GeV^{2}/c^{4})", "<n_{ch}>", true, false);
  htemp->GetYaxis()->SetTitleSize(0.07);
  htemp->GetYaxis()->SetTitleOffset(0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].nch_w_f[0],24,mccolors[i],0.8,"c");
  }
  DrawGraph(prd27_47_fig7_vp_cur,27,1,0.8);
  DrawGraph(npb223_269_nF_W,24,1,0.8);
  DrawGraph(zpc24_119_fig4_vp_F,20,1,0.8);
  TLegend *lMulFB_1 = new TLegend(0.65,0.05,0.85,0.35);
  SetLegend(lMulFB_1);
  lMulFB_1->AddEntry(npb223_269_nF_W,"BEBC #nuH_{2}","p");
  lMulFB_1->AddEntry(zpc24_119_fig4_vp_F,"BEBC #nuD_{2}","p");
  lMulFB_1->AddEntry(prd27_47_fig7_vp_cur,"15' #nuD_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lMulFB_1->AddEntry(hadPlots[i].nch_w_f[0],name.c_str(),"l");
    lMulFB_1->Draw();
  }
  TLatex *tMulFB_1 = new TLatex(3,4,"(a) #nup forward");
  tMulFB_1->SetTextSize(0.085);
  tMulFB_1->Draw();
  cMulFB->cd(2);
  //gPad->SetTopMargin(0);
  gPad->SetBottomMargin(0);
  gPad->SetLeftMargin(0);
  //gPad->SetRightMargin(2);
  //gPad->SetTickx(2);
  gPad->SetTicky(2);
  htemp = DrawFrame(cMulFB, 2, 1, 300, -0.2, 5.5, "","W^{2}(GeV^{2}/c^{4})", "", true, false);
  htemp->GetYaxis()->SetLabelSize(0);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].nch_w_b[0],24,mccolors[i],0.8,"c");
  }
  DrawGraph(prd27_47_fig7_vp_tar,27,1,0.8);
  DrawGraph(npb223_269_nB_W,24,1,0.8);
  DrawGraph(zpc24_119_fig4_vp_B,20,1,0.8);
  lMulFB_1->Draw();
  TLatex *tMulFB_2 = new TLatex(3,4,"(b) #nup backward");
  tMulFB_2->SetTextSize(0.085);
  tMulFB_2->Draw();
  cMulFB->cd(3);
  gPad->SetTopMargin(0);
  gPad->SetRightMargin(0);
  DrawFrame(cMulFB, 3, 1, 300, -0.2, 5.5, "","W^{2}(GeV^{2}/c^{4})", "<n_{ch}>", true, false);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].nch_w_f[1],24,mccolors[i],0.8,"c");
  }
  DrawGraph(prd27_47_fig7_vn_cur,27,1,0.8);
  DrawGraph(zpc24_119_fig4_vn_F,20,1,0.8);
  TLegend *lMulFB_3 = new TLegend(0.65,0.25,0.85,0.45);
  SetLegend(lMulFB_3);
  lMulFB_3->AddEntry(zpc24_119_fig4_vn_F,"BEBC #nuD_{2}","p");
  lMulFB_3->AddEntry(prd27_47_fig7_vn_cur,"15' #nuD_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lMulFB_3->AddEntry(hadPlots[i].nch_w_f[1],name.c_str(),"l");
  }
  lMulFB_3->Draw();
  TLatex *tMulFB_3 = new TLatex(3,4,"(c) #nun forward");
  tMulFB_3->SetTextSize(0.07);
  tMulFB_3->Draw();
  cMulFB->cd(4);
  gPad->SetTopMargin(0);
  gPad->SetLeftMargin(0);
  //gPad->SetRightMargin(2);
  gPad->SetTicky(2);
  htemp = DrawFrame(cMulFB, 4, 1, 300, -0.2, 5.5, "","W^{2}(GeV^{2}/c^{4})", "", true, false);
  htemp->GetYaxis()->SetLabelSize(0);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].nch_w_b[1],20,mccolors[i],0.8,"c");
  }
  DrawGraph(prd27_47_fig7_vn_tar,27,1,0.8);
  DrawGraph(zpc24_119_fig4_vn_B,20,1,0.8);
  lMulFB_3->Draw();
  TLatex *tMulFB_4 = new TLatex(3,4,"(d) #nun backward");
  tMulFB_4->SetTextSize(0.07);
  tMulFB_4->Draw();
  
   //..................xF distribution.................................
  TGraphErrors *npb214_369_fig7_pip = MakeGraph("npb214_369_fig7_pip.dat");
  TGraphErrors *npb214_369_fig7_pim = MakeGraph("npb214_369_fig7_pim.dat");
  TCanvas *cXF = new TCanvas("cXF","cXF",800,400);
  cXF->Divide(2,1);
  DrawFrame(cXF, 1, -1, 1, 0.01, 10, "","X_{F}", "(1/N_{ev})dN^{#pi}/dX_{F}", false, true);
  DrawGraph(npb214_369_fig7_pip,20,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].xf_pip[0]->SetLineColor(mccolors[i]);
    hadPlots[i].xf_pip[0]->SetMarkerStyle(1);
    hadPlots[i].xf_pip[0]->Draw("hist c same");
  }
  TLegend *lXF = new TLegend(0.45,0.25,0.75,0.45);
  SetLegend(lXF);
  lXF->SetHeader("W>3GeV");
  lXF->AddEntry(npb214_369_fig7_pip, "#nup, BEBC #nuH_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lXF->AddEntry(hadPlots[i].xf_pip[0],name.c_str(),"l");
  }
  lXF->Draw();
  TLatex *txf1 = new TLatex(-0.9,5,"#pi^{+} from #nup");
  txf1->SetTextSize(0.06);
  txf1->Draw();
  DrawFrame(cXF, 2, -1, 1, 0.01, 10, "","X_{F}", "(1/N_{ev})dN^{#pi}/dX_{F}", false, true);
  DrawGraph(npb214_369_fig7_pim,20,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].xf_pim[0]->SetLineColor(mccolors[i]);
    hadPlots[i].xf_pim[0]->SetMarkerStyle(1);
    hadPlots[i].xf_pim[0]->Draw("hist c same");
  }
  TLatex *txf2 = new TLatex(-0.9,5,"#pi^{-} from #nup");
  txf2->SetTextSize(0.06);
  txf2->Draw();
  lXF->Draw();

  //...........................z distribution.............................
  TGraphErrors *zpc24_119_fig8_vd_pos = MakeGraph("zpc24_119_fig8_vd_pos.dat");
  TGraphErrors *zpc24_119_fig8_vd_neg = MakeGraph("zpc24_119_fig8_vd_neg.dat");
  TCanvas *cZ = new TCanvas("cZ","cZ",800,400);
  cZ->Divide(2,1);
  DrawFrame(cZ, 1, 0, 1, 0.001, 40, "","z", "D^{+}(z)", false, true);
  DrawGraph(zpc24_119_fig8_vd_pos,20,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].z_pos[0]->SetLineColor(mccolors[i]);
    hadPlots[i].z_pos[0]->SetMarkerStyle(1);
    hadPlots[i].z_pos[0]->Draw("hist c same");
    hadPlots[i].z_pos[1]->SetLineColor(mccolors[i]);
    hadPlots[i].z_pos[1]->SetLineStyle(2);
    hadPlots[i].z_pos[1]->SetMarkerStyle(1);
    hadPlots[i].z_pos[1]->Draw("hist c same");
  }
  TLegend *leg_z = new TLegend(0.2,0.25,0.55,0.45);
  SetLegend(leg_z);
  leg_z->AddEntry(zpc24_119_fig8_vd_pos,"#nuD_{2} BEBC","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = "#nup "+hadPlots[i].modelName;
    leg_z->AddEntry(hadPlots[i].z_pos[0],name.c_str(),"l");
    name = "#nun "+hadPlots[i].modelName;
    leg_z->AddEntry(hadPlots[i].z_pos[1],name.c_str(),"l");
  }
  leg_z->Draw();
  TLatex *tz1 = new TLatex(0.6,10,"h^{+}");
  tz1->SetTextSize(0.06);
  tz1->Draw();
  DrawFrame(cZ, 2, 0, 1, 0.001, 40, "","z", "D^{-}(z)", false, true);
  DrawGraph(zpc24_119_fig8_vd_neg,20,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].z_neg[0]->SetLineColor(mccolors[i]);
    hadPlots[i].z_neg[0]->SetMarkerStyle(1);
    hadPlots[i].z_neg[0]->Draw("hist c same");
    hadPlots[i].z_neg[1]->SetLineColor(mccolors[i]);
    hadPlots[i].z_neg[1]->SetLineStyle(2);
    hadPlots[i].z_neg[1]->SetMarkerStyle(1);
    hadPlots[i].z_neg[1]->Draw("hist c same");
  }
  leg_z->Draw();
  TLatex *tz2 = new TLatex(0.6,10,"h^{-}");
  tz2->SetTextSize(0.06);
  tz2->Draw();

  /*
  //..................pT distribution......................
  TGraphErrors *prd19_1_fig11_xf1 = MakeGraph("prd19_1_fig11_xf1.dat");
  TGraphErrors *prd19_1_fig11_xf2 = MakeGraph("prd19_1_fig11_xf2.dat");
  TGraphErrors *prd19_1_fig11_xf3 = MakeGraph("prd19_1_fig11_xf3.dat");
  TGraphErrors *prd19_1_fig11_xf4 = MakeGraph("prd19_1_fig11_xf4.dat");
  TCanvas *cPt = new TCanvas("cPt","cPt",800,800);
  cPt->Divide(2,2);
  DrawFrame(cPt,1,0,1.2,0.01,30,"0<X_{F}<0.1","P_{T}^{2}[(GeV/c)^{2}]","h^{#pm}Tracks/Event/(GeV/c)^2",false,true);
  DrawGraph(prd19_1_fig11_xf1,20,1,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt2_xf1[0]->SetLineColor(mccolors[i]);
    hadPlots[i].pt2_xf1[0]->SetMarkerStyle(1);
    hadPlots[i].pt2_xf1[0]->Draw("hist c same");
  }
  TLegend *leg_pt = new TLegend(0.45,0.5,0.85,0.7);
  SetLegend(leg_pt);
  leg_pt->AddEntry(prd19_1_fig11_xf1,"#nup 15'#nuH_{2}, W>4GeV","l");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = "#nup "+hadPlots[i].modelName;
    leg_pt->AddEntry(hadPlots[i].pt2_xf1[0],name.c_str(),"l");
  }
  leg_pt->Draw();

  DrawFrame(cPt,2,0,1.2,0.01,30,"0.1<X_{F}<0.3","P_{T}^{2}[(GeV/c)^{2}]","h^{#pm}Tracks/Event/(GeV/c)^2",false,true);
  DrawGraph(prd19_1_fig11_xf2,20,1,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt2_xf2[0]->SetLineColor(mccolors[i]);
    hadPlots[i].pt2_xf2[0]->SetMarkerStyle(1);
    hadPlots[i].pt2_xf2[0]->Draw("hist c same");
  }

  DrawFrame(cPt,3,0,1.2,0.01,30,"0.3<X_{F}<0.6","P_{T}^{2}[(GeV/c)^{2}]","h^{#pm}Tracks/Event/(GeV/c)^2",false,true);
  DrawGraph(prd19_1_fig11_xf3,20,1,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt2_xf3[0]->SetLineColor(mccolors[i]);
    hadPlots[i].pt2_xf3[0]->SetMarkerStyle(1);
    hadPlots[i].pt2_xf3[0]->Draw("hist c same");
  }

  DrawFrame(cPt,4,0,1.2,0.001,3,"0.6<X_{F}<1.0","P_{T}^{2}[(GeV/c)^{2}]","h^{#pm}Tracks/Event/(GeV/c)^2",false,true);
  DrawGraph(prd19_1_fig11_xf4,20,1,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt2_xf4[0]->SetLineColor(mccolors[i]);
    hadPlots[i].pt2_xf4[0]->SetMarkerStyle(1);
    hadPlots[i].pt2_xf4[0]->Draw("hist c same");
  }
  */

  //......................pT2 vs W2....................................
  TGraphErrors *zpc27_239_fig3_vD_F = MakeGraph("zpc27_239_fig3_vD_F.dat");
  TGraphErrors *zpc27_239_fig3_vD_B = MakeGraph("zpc27_239_fig3_vD_B.dat");
  TCanvas *cPt_W2 = new TCanvas("cPt_W2","cPt_W2",800,400);
  cPt_W2->Divide(2,1);
  DrawFrame(cPt_W2,1,1,300,0,1,"","W^{2}(GeV^{2}/c^{4})","<P_{T}^{2}>(GeV/c)^{2}",true,false);
  DrawGraph(zpc27_239_fig3_vD_F,20,1);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt2_W2_F[0]->SetLineColor(mccolors[i]);
    hadPlots[i].pt2_W2_F[0]->Draw("hist c same");
    hadPlots[i].pt2_W2_F[1]->SetLineColor(mccolors[i]);
    hadPlots[i].pt2_W2_F[1]->SetLineStyle(2);
    hadPlots[i].pt2_W2_F[1]->Draw("hist c same");
  }
  TLegend *leg_pt2_w2_F = new TLegend(0.2,0.6,0.6,0.85);
  SetLegend(leg_pt2_w2_F);
  leg_pt2_w2_F->SetHeader("X_{F}>0.3");
  leg_pt2_w2_F->AddEntry(zpc27_239_fig3_vD_F,"#nuD_{2}BEBC","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = "#nup "+hadPlots[i].modelName;
    leg_pt2_w2_F->AddEntry(hadPlots[i].pt2_W2_F[0],name.c_str(),"l");
    name = "#nun "+hadPlots[i].modelName;
    leg_pt2_w2_F->AddEntry(hadPlots[i].pt2_W2_F[1],name.c_str(),"l");
  }
  leg_pt2_w2_F->Draw();
  
  DrawFrame(cPt_W2,2,1,300,0,1,"","W^{2}(GeV^{2}/c^{4})","<P_{T}^{2}>(GeV/c)^{2}",true,false);  
  DrawGraph(zpc27_239_fig3_vD_B,20,1);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt2_W2_B[0]->SetLineColor(mccolors[i]);
    hadPlots[i].pt2_W2_B[0]->Draw("hist c same");
    hadPlots[i].pt2_W2_B[1]->SetLineColor(mccolors[i]);
    hadPlots[i].pt2_W2_B[1]->SetLineStyle(2);
    hadPlots[i].pt2_W2_B[1]->Draw("hist c same");
  }
  TLegend *leg_pt2_w2_B = new TLegend(0.2,0.6,0.6,0.85);
  SetLegend(leg_pt2_w2_B);
  leg_pt2_w2_B->SetHeader("X_{F}<0.3");
  leg_pt2_w2_B->AddEntry(zpc27_239_fig3_vD_B,"#nuD_{2} BEBC","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = "#nup "+hadPlots[i].modelName;
    leg_pt2_w2_B->AddEntry(hadPlots[i].pt2_W2_B[0],name.c_str(),"l");
    name = "#nun "+hadPlots[i].modelName;
    leg_pt2_w2_B->AddEntry(hadPlots[i].pt2_W2_B[1],name.c_str(),"l");
  }
  leg_pt2_w2_B->Draw();

  //..............................pt2 vs xf..............................
  TGraphErrors *zpc27_239_fig6_vp_loW = MakeGraph("zpc27_239_fig6_vp_loW.dat");
  TGraphErrors *zpc27_239_fig6_vn_loW = MakeGraph("zpc27_239_fig6_vn_loW.dat");
  TCanvas *cPt2_xf = new TCanvas("cPt2_xf","cPt2_xf",800,400);
  cPt2_xf->Divide(2,1);
  DrawFrame(cPt2_xf,1,-1,1,0,0.8,"","X_{F}","<P_{T}^{2}>(GeV/c)^{2}",false,false);
  DrawGraph(zpc27_239_fig6_vp_loW,20,1);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt2_xf_loW[0]->SetLineColor(mccolors[i]);
    hadPlots[i].pt2_xf_loW[0]->Draw("hist c same");
  }
  TLegend *leg_pt2_xf_p = new TLegend(0.3,0.6,0.8,0.85);
  SetLegend(leg_pt2_xf_p);
  leg_pt2_xf_p->SetHeader("#nup, 9<W^{2}<25GeV^{2}");
  leg_pt2_xf_p->AddEntry(zpc27_239_fig6_vp_loW,"BEBC #nuD_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    leg_pt2_xf_p->AddEntry(hadPlots[i].pt2_xf_loW[0],name.c_str(),"l");
  }
  leg_pt2_xf_p->Draw();

  DrawFrame(cPt2_xf,2,-1,1,0,0.8,"","X_{F}","<P_{T}^{2}>(GeV/c)^{2}",false,false);
  DrawGraph(zpc27_239_fig6_vn_loW,20,1);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt2_xf_loW[1]->SetLineColor(mccolors[i]);
    hadPlots[i].pt2_xf_loW[1]->Draw("hist c same");
  }
  TLegend *leg_pt2_xf_n = new TLegend(0.3,0.6,0.8,0.85);
  SetLegend(leg_pt2_xf_n);
  leg_pt2_xf_n->SetHeader("#nun, 9<W^{2}<25GeV^{2}");
  leg_pt2_xf_n->AddEntry(zpc27_239_fig6_vn_loW,"BEBC #nuD_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    leg_pt2_xf_n->AddEntry(hadPlots[i].pt2_xf_loW[1],name.c_str(),"l");
  }
  leg_pt2_xf_n->Draw();

  cMulCh->Print("cMulCh.gif");
  cDispCh->Print("cDispCh.gif");
  cKNO->Print("cKNO.gif");
  cMulPi0->Print("cMulPi0.gif");
  cCorr_Pi0_Ch->Print("cCorr_Pi0_Ch.gif");
  cTopo->Print("cTopo.gif");
  cMulFB->Print("cMulFB.gif");
  cXF->Print("cXF.gif");
  cZ->Print("cZ.gif");
  cPt_W2->Print("cPt_W2.gif");
  cPt2_xf->Print("cPt2_xf.gif");

  //______________________anti-neutrino plots____________________________

  //................charged hadron multiplicity...........................
  TGraphErrors *prd25_624_fig5 = MakeGraph("prd25_624_fig5.dat");
  TCanvas *cMulCh_nubar = new TCanvas("cMulCh_nubar","cMulCh_nubar",600,400);
  DrawFrame(cMulCh_nubar,0,1,300,0,10,"","W^{2}(GeV^{2}/c^{4})","<n_{ch}>",true,false);
  DrawGraph(prd25_624_fig5,24,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].nch_w[2],24,mccolors[i],0.8,"c");
  }
  TLegend *lMulCh_nubar = new TLegend(0.65,0.25,0.8,0.45);
  SetLegend(lMulCh_nubar);
  lMulCh_nubar->AddEntry(prd25_624_fig5,"15' #bar{#nu}H_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lMulCh_nubar->AddEntry(hadPlots[i].nch_w[2],name.c_str(),"l");
  }
  lMulCh_nubar->Draw();
  TLatex *tMulCh_nubar = new TLatex(3,7,"#bar{#nu}p#rightarrow#mu^{+}X^{0}");
  tMulCh_nubar->SetTextSize(0.06);
  tMulCh_nubar->Draw();

  //..................dispersion (D_)...............................
  TGraphErrors *prd25_624_fig8 = MakeGraph("prd25_624_fig8.dat");
  TCanvas *cDispCh_nubar = new TCanvas("cDispCh_nubar","cDispCh_nubar",600,400);
  DrawFrame(cDispCh_nubar,0,0,4,0,2,"","<n_{-}>","D_{-}",false,false);
  DrawGraph(prd25_624_fig8,24,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].D_nneg[2],24,mccolors[i],0.8,"c");
  }
  TLegend *lDispCh_nubar = new TLegend(0.75,0.25,0.95,0.4);
  SetLegend(lDispCh_nubar);
  lDispCh_nubar->AddEntry(prd25_624_fig8,"15' #bar{#nu}H_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lDispCh_nubar->AddEntry(hadPlots[i].D_nneg[2],name.c_str(),"l");
  }
  lDispCh_nubar->Draw();
  TLatex *tDispCh_nubar = new TLatex(0.5,1.4,"#bar{#nu}p#rightarrow#mu^{+}X^{0}");
  tDispCh_nubar->SetTextSize(0.06);
  tDispCh_nubar->Draw();

  //...........normalised topological cross sections..................

  TGraphErrors *prd25_624_0 = MakeGraph("prd25_624_0.dat");
  TGraphErrors *prd25_624_2 = MakeGraph("prd25_624_2.dat");
  TGraphErrors *prd25_624_4 = MakeGraph("prd25_624_4.dat");
  TGraphErrors *prd25_624_6 = MakeGraph("prd25_624_6.dat");
  TGraphErrors *prd25_624_8 = MakeGraph("prd25_624_8.dat");
  TGraphErrors *prd25_624_10 = MakeGraph("prd25_624_10.dat");
  TCanvas *cTopo_nubar = new TCanvas("cTopo_nubar","cTopo_nubar",600,370);
  DrawFrame(cTopo_nubar,0,1,1000,1e-4,10,"","W^{2}(GeV^{2}/c^{4})","P(n_{ch})",true,true);
  DrawGraph(prd25_624_0,20,1,0.7);
  DrawGraph(prd25_624_2,20,2,0.7);
  DrawGraph(prd25_624_4,20,3,0.7);
  DrawGraph(prd25_624_6,20,4,0.7);
  DrawGraph(prd25_624_8,20,7,0.7);
  DrawGraph(prd25_624_10,20,6,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pnch_w2[0][2]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[2][2]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[4][2]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[6][2]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[8][2]->SetLineStyle(i+1);
    hadPlots[i].pnch_w2[10][2]->SetLineStyle(i+1);
    DrawGraph(hadPlots[i].pnch_w2[0][2],20,1,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[2][2],20,2,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[4][2],20,3,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[6][2],20,4,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[8][2],20,7,0.8,"c");
    DrawGraph(hadPlots[i].pnch_w2[10][2],20,6,0.8,"c");
  }
  TLegend *lTopo_nubar = new TLegend(0.7,0.45,0.95,0.85);
  SetLegend(lTopo_nubar);
  lTopo_nubar->AddEntry(prd25_624_0,"n = 0, 15' #bar{#nu}H_{2}(MC)","p");
  lTopo_nubar->AddEntry(prd25_624_2,"n = 2, 15' #bar{#nu}H_{2}","p");
  lTopo_nubar->AddEntry(prd25_624_4,"n = 4, 15' #bar{#nu}H_{2}","p");
  lTopo_nubar->AddEntry(prd25_624_6,"n = 6, 15' #bar{#nu}H_{2}","p");
  lTopo_nubar->AddEntry(prd25_624_8,"n = 8, 15' #bar{#nu}H_{2}","p");
  lTopo_nubar->AddEntry(prd25_624_10,"n = 10, 15' #bar{#nu}H_{2}","p");
  lTopo_nubar->Draw();

  //.......................pi0 multiplicity..............................
  TGraphErrors *nc51a_539_fig4 = MakeGraph("nc51a_539_fig4.dat");
  TGraphErrors *npb223_269_fig8b = MakeGraph("npb223_269_fig8b.dat");
  TGraphErrors *zpc40_231_fig13_nubar = MakeGraph("zpc40_231_fig13_nubar.dat");
  TCanvas *cMulPi0_nubar = new TCanvas("cMulPi0_nubar","cMulPi0_nubar");
  DrawFrame(cMulPi0_nubar,0,1,1000,0,4,"","W^{2}(GeV^{2}/c^{4})","<n_{#pi^{0}}>",true,false);
  //DrawGraph(sn41_963_fig2,20,1,0.8);
  DrawGraph(nc51a_539_fig4,24,1,0.8);
  DrawGraph(npb223_269_fig8b,26,1,0.8);
  DrawGraph(zpc40_231_fig13_nubar,20,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].npi0_w[2],24,mccolors[i],0.8,"c");
    hadPlots[i].npi0_w[3]->SetLineStyle(2);
    DrawGraph(hadPlots[i].npi0_w[3],24,mccolors[i],0.8,"c");
  }
  TLegend *lMulPi0_nubar = new TLegend(0.7,0.25,0.95,0.55);
  SetLegend(lMulPi0_nubar);
  lMulPi0_nubar->AddEntry(nc51a_539_fig4,"FNAL 15' #bar{#nu}neon+H_{2}","p");
  lMulPi0_nubar->AddEntry(zpc40_231_fig13_nubar,"BEBC #bar{#nu}neon+H_{2}","p");
  lMulPi0_nubar->AddEntry(npb223_269_fig8b,"BEBC #bar{#nu}H_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = "#bar{#nu}p "+hadPlots[i].modelName;
    lMulPi0_nubar->AddEntry(hadPlots[i].npi0_w[2],name.c_str(),"l");
    name = "#bar{#nu}n "+hadPlots[i].modelName;
    lMulPi0_nubar->AddEntry(hadPlots[i].npi0_w[3],name.c_str(),"l");
  }
  lMulPi0_nubar->Draw();
  
  //.................F/B multiplicity.........................
  TGraphErrors *prd24_1071_fig21a = MakeGraph("prd24_1071_fig21a.dat");
  TGraphErrors *prd24_1071_fig21b = MakeGraph("prd24_1071_fig21b.dat");
  TCanvas *cMulFB_nubar = new TCanvas("cMulFB_nubar","cMulFB_nubar",800,400);
  cMulFB_nubar->Divide(2,1);
  DrawFrame(cMulFB_nubar, 1, 1, 300, 0, 5, "","W^{2}(GeV^{2}/c^{4})", "<n_{ch}>", true, false);
  DrawGraph(prd24_1071_fig21a,24,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].nch_w_f[2],24,mccolors[i],0.8,"c");
  }
  TLegend *lMulFB1_nubar = new TLegend(0.2,0.7,0.5,0.85);
  SetLegend(lMulFB1_nubar);
  lMulFB1_nubar->AddEntry(prd24_1071_fig21a,"15' #nuH_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lMulFB1_nubar->AddEntry(hadPlots[i].nch_w_f[2],name.c_str(),"l");
  }
  lMulFB1_nubar->Draw();
  TLatex *tMulFB1_nubar = new TLatex(0.6,0.3,"#bar{#nu}p forward");
  tMulFB1_nubar->SetNDC();
  tMulFB1_nubar->SetTextSize(0.06);
  tMulFB1_nubar->Draw();

  DrawFrame(cMulFB_nubar, 2, 1, 300, 0, 5, "","W^{2}(GeV^{2}/c^{4})", "<n_{ch}>", true, false);
  DrawGraph(prd24_1071_fig21b,24,1,0.8);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    DrawGraph(hadPlots[i].nch_w_b[2],24,mccolors[i],0.8,"c");
  }
  TLegend *lMulFB2_nubar = new TLegend(0.2,0.7,0.5,0.85);
  SetLegend(lMulFB2_nubar);
  lMulFB2_nubar->AddEntry(prd24_1071_fig21a,"15' #bar{#nu}H_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lMulFB2_nubar->AddEntry(hadPlots[i].nch_w_b[2],name.c_str(),"l");
  }
  lMulFB2_nubar->Draw();
  TLatex *tMulFB2_nubar = new TLatex(0.6,0.3,"#bar{#nu}p backward");
  tMulFB2_nubar->SetNDC();
  tMulFB2_nubar->SetTextSize(0.06);
  tMulFB2_nubar->Draw();

  //..................xf distribution......................
  TGraphErrors *prd24_1071_fig8a = MakeGraph("prd24_1071_fig8a.dat");
  //TGraphErrors *prd24_1071_fig8b = MakeGraph("prd24_1071_fig8b.dat");
  TGraphErrors *prd24_1071_fig9a = MakeGraph("prd24_1071_fig9a.dat");
  //TGraphErrors *prd24_1071_fig9b = MakeGraph("prd24_1071_fig9b.dat");
  TGraphErrors *prd24_1071_fig8a_hi = MakeGraph("prd24_1071_fig8a_hi.dat");
  //TGraphErrors *prd24_1071_fig8b_hi = MakeGraph("prd24_1071_fig8b_hi.dat");
  TGraphErrors *prd24_1071_fig9a_hi = MakeGraph("prd24_1071_fig9a_hi.dat");
  //TGraphErrors *prd24_1071_fig9b_hi = MakeGraph("prd24_1071_fig9b_hi.dat");
  TCanvas *cXf_nubar = new TCanvas("cXf_nubar","cXf_nubar",800,800);
  cXf_nubar->Divide(2,2);
  DrawFrame(cXf_nubar,1,-1,1,0,.5,"","X_{F}=2P_{L}/W","F(X_{F})",false,false);
  DrawGraph(prd24_1071_fig8a,20,1,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].Fxf_pos1[2]->SetLineColor(mccolors[i]);
    hadPlots[i].Fxf_pos1[2]->Draw("hist c same");
  }
  TLatex *tXf1_nubar = new TLatex(-0.8,0.45,"2<W<4GeV,Q^{2}<45(GeV/c)^{2}");
  tXf1_nubar->SetTextSize(0.05);
  tXf1_nubar->Draw();
  TLatex *tXf2_nubar = new TLatex(-0.8,0.4,"Positive Hadrons");
  tXf2_nubar->SetTextSize(0.05);
  tXf2_nubar->Draw();
  TLegend *lXf1_nubar = new TLegend(0.18,0.6,0.48,0.7);
  SetLegend(lXf1_nubar);
  lXf1_nubar->AddEntry(prd24_1071_fig8a,"15' #bar{#nu}H_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lXf1_nubar->AddEntry(hadPlots[i].Fxf_pos1[2],name.c_str(),"l");
  }
  lXf1_nubar->Draw();

  DrawFrame(cXf_nubar,2,-1,1,0,.5,"","X_{F}=2P_{L}/W","F(X_{F})",false,false);
  DrawGraph(prd24_1071_fig9a,20,1,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].Fxf_neg1[2]->SetLineColor(mccolors[i]);
    hadPlots[i].Fxf_neg1[2]->Draw("hist c same");
  }
  tXf1_nubar->Draw();
  TLatex *tXf3_nubar = new TLatex(-0.8,0.4,"Negative Hadrons");
  tXf3_nubar->SetTextSize(0.05);
  tXf3_nubar->Draw();
  lXf1_nubar->Draw();

  DrawFrame(cXf_nubar,3,-1,1,0,0.5,"","X_{F}=2P_{L}/W","F(X_{F})",false,false);
  DrawGraph(prd24_1071_fig8a_hi,20,1,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].Fxf_pos1_hi[2]->SetLineColor(mccolors[i]);
    hadPlots[i].Fxf_pos1_hi[2]->Draw("hist c same");
  }
  TLatex *tXf4_nubar = new TLatex(-0.8,0.45,"4<W<10GeV,Q^{2}<45(GeV/c)^{2}");
  tXf4_nubar->SetTextSize(0.05);
  tXf4_nubar->Draw();
  tXf2_nubar->Draw();
  lXf1_nubar->Draw();

  DrawFrame(cXf_nubar,4,-1,1,0,0.5,"","X_{F}=2P_{L}/W","F(X_{F})",false,false);
  DrawGraph(prd24_1071_fig9a_hi,20,1,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].Fxf_neg1_hi[2]->SetLineColor(mccolors[i]);
    hadPlots[i].Fxf_neg1_hi[2]->Draw("hist c same");
  }
  tXf4_nubar->Draw();
  tXf3_nubar->Draw();
  lXf1_nubar->Draw();

  //Nucl.Phys.B214:369,1983
  TGraphErrors *npb214_369_fig5c = MakeGraph("npb214_369_fig5c.dat");
  TGraphErrors *npb214_369_fig5d = MakeGraph("npb214_369_fig5d.dat");
  TGraphErrors *zpc24_119_fig6g = MakeGraph("zpc24_119_fig6g.dat");
  TGraphErrors *zpc24_119_fig6c = MakeGraph("zpc24_119_fig6c.dat");
  TCanvas *cXf2_nubar = new TCanvas("cXf2_nubar","cXf2_nubar",800,400);
  cXf2_nubar->Divide(2,1);
  DrawFrame(cXf2_nubar,1,-1,1,0,.3,"","X_{F}=2P_{L}/W","F(X_{F})",false,false);
  DrawGraph(npb214_369_fig5c,20,1,0.7);
  DrawGraph(zpc24_119_fig6g,24,1,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].fxf_pip[2]->SetLineColor(mccolors[i]);
    hadPlots[i].fxf_pip[2]->Draw("hist c same");
  }
  TLegend *lXf2_nubar = new TLegend(0.18,0.7,0.47,0.9);
  SetLegend(lXf2_nubar);
  lXf2_nubar->AddEntry(npb214_369_fig5c,"WA21,BEBC #bar{#nu}H_{2}","p");
  lXf2_nubar->AddEntry(zpc24_119_fig6g,"WA25,BEBC #bar{#nu}D_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    lXf2_nubar->AddEntry(hadPlots[0].fxf_pip[2],name.c_str(),"l");
  }
  lXf2_nubar->Draw();
  TLatex *piplus = new TLatex(0.4,0.25,"#pi^{+}");
  piplus->SetTextSize(0.06);
  piplus->Draw();

  DrawFrame(cXf2_nubar,2,-1,1,0,.3,"","X_{F}=2P_{L}/W","F(X_{F})",false,false);
  DrawGraph(npb214_369_fig5d,20,1,0.7);
  DrawGraph(zpc24_119_fig6c,24,1,0.7);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].fxf_pim[2]->SetLineColor(mccolors[i]);
    hadPlots[i].fxf_pim[2]->Draw("hist c same");
  }
  lXf2_nubar->Draw();
  TLatex *piminus = new TLatex(0.4,0.25,"#pi^{-}");
  piminus->SetTextSize(0.06);
  piminus->Draw();

  //........................z from anti-nu..................
  TGraphErrors *prd17_fig14_E1 = MakeGraph("prd17_fig14_E1.dat");
  TGraphErrors *prd17_fig14_E2 = MakeGraph("prd17_fig14_E2.dat");
  TGraphErrors *prd17_fig14_E3 = MakeGraph("prd17_fig14_E3.dat");
  TCanvas *cZ_nubar = new TCanvas("cZ_nubar","cZ_nubar",600,600);
  DrawFrame(cZ_nubar,0,0,1,0.1,100,"","z","(1/N_{E})dN_{E}/dZ",false,true);
  DrawGraph(prd17_fig14_E1,20,1);
  DrawGraph(prd17_fig14_E2,24,2);
  DrawGraph(prd17_fig14_E3,25,4);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].z_E1[2]->SetLineStyle(i+1);
    hadPlots[i].z_E2[2]->SetLineStyle(i+1);
    hadPlots[i].z_E3[2]->SetLineStyle(i+1);
    hadPlots[i].z_E1[2]->SetLineColor(1);
    hadPlots[i].z_E2[2]->SetLineColor(2);
    hadPlots[i].z_E3[2]->SetLineColor(4);
    hadPlots[i].z_E1[2]->Draw("hist c same");
    hadPlots[i].z_E2[2]->Draw("hist c same");
    hadPlots[i].z_E3[2]->Draw("hist c same");
  }
  TLegend *leg_z_nubar = new TLegend(0.4,0.6,0.85,0.85);
  SetLegend(leg_z_nubar);
  leg_z_nubar->AddEntry(prd17_fig14_E1,"#bar{#nu}p, 15'#bar{#nu}H_{2} E_{#nu}<15GeV","p");
  leg_z_nubar->AddEntry(prd17_fig14_E2,"#bar{#nu}p, 15'#bar{#nu}H_{2} 15<E_{#nu}<30GeV","p");
  leg_z_nubar->AddEntry(prd17_fig14_E3,"#bar{#nu}p, 15'#bar{#nu}H_{2} E_{#nu}>30GeV","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName+", E_{#nu}<15GeV";
    leg_z_nubar->AddEntry(hadPlots[i].z_E1[2],name.c_str(),"l");
    name = hadPlots[i].modelName+", 15<E_{#nu}<30GeV";
    leg_z_nubar->AddEntry(hadPlots[i].z_E2[2],name.c_str(),"l");
    name = hadPlots[i].modelName+", E_{#nu}>30GeV";
    leg_z_nubar->AddEntry(hadPlots[i].z_E3[2],name.c_str(),"l");
  }
  leg_z_nubar->Draw();

  //......................pT vs W....................................  
  TGraphErrors *prd24_1071_fig37a = MakeGraph("prd24_1071_fig37a.dat");
  TGraphErrors *prd24_1071_fig37b = MakeGraph("prd24_1071_fig37b.dat");
  TGraphErrors *prd24_1071_fig37c = MakeGraph("prd24_1071_fig37c.dat");
  TCanvas *cPt_W_nubar = new TCanvas("cPt_W_nubar","cPt_W_nubar",900,300);
  cPt_W_nubar->Divide(3,1);
  DrawFrame(cPt_W_nubar,1,1,10,0,0.6,"","W(GeV/c^{2})","<P_{T}>(GeV/c)",false,false);
  DrawGraph(prd24_1071_fig37a,20,1);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt_W_F[2]->SetLineColor(mccolors[i]);
    hadPlots[i].pt_W_F[2]->Draw("hist c same");
  }
  TLegend *leg_pt_w = new TLegend(0.5,0.3,0.8,0.5);
  SetLegend(leg_pt_w);
  leg_pt_w->AddEntry(prd24_1071_fig37a,"#bar{#nu}p 15' #bar{#nu}D_{2}","p");
  for (unsigned i = 0; i<hadPlots.size(); i++){
    string name = hadPlots[i].modelName;
    leg_pt_w->AddEntry(hadPlots[i].pt_W_F[2],name.c_str(),"l");
  }
  leg_pt_w->Draw();
  TLatex *tptw1 = new TLatex(3,0.5,"x_{F}>0");
  tptw1->SetTextSize(0.08);
  tptw1->Draw();

  DrawFrame(cPt_W_nubar,2,1,10,0,0.6,"","W(GeV/c^{2})","<P_{T}>(GeV/c)",false,false);
  DrawGraph(prd24_1071_fig37b,20,1);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt_W_B[2]->SetLineColor(mccolors[i]);
    hadPlots[i].pt_W_B[2]->Draw("hist c same");
  }
  leg_pt_w->Draw();
  TLatex *tptw2 = new TLatex(3,0.5,"x_{F}<0");
  tptw2->SetTextSize(0.08);
  tptw2->Draw();

  DrawFrame(cPt_W_nubar,3,1,10,0,0.6,"","W(GeV/c^{2})","<P_{T}>(GeV/c)",false,false);
  DrawGraph(prd24_1071_fig37c,20,1);
  for (unsigned i = 0; i<hadPlots.size(); i++){
    hadPlots[i].pt_W[2]->SetLineColor(mccolors[i]);
    hadPlots[i].pt_W[2]->Draw("hist c same");
  }
  leg_pt_w->Draw();
  TLatex *tptw3 = new TLatex(3,0.5,"All x_{F}");
  tptw3->SetTextSize(0.08);
  tptw3->Draw();

  cMulCh_nubar->Print("cMulCh_nubar.gif");
  cDispCh_nubar->Print("cDispCh_nubar.gif");
  cTopo_nubar->Print("cTopo_nubar.gif");
  cMulPi0_nubar->Print("cMulPi0_nubar.gif");
  cMulFB_nubar->Print("cMulFB_nubar.gif");
  cXf_nubar->Print("cXf_nubar.gif");
  cXf2_nubar->Print("cXf2_nubar.gif");
  cZ_nubar->Print("cZ_nubar.gif");
  cPt_W_nubar->Print("cPt_W_nubar.gif");

}

TGraphErrors* MakeGraph(string file)
{
  string filename = "data/"+file;
  ifstream in;
  in.open(filename.c_str());
  double x, y, ex, ey;
  vector<double> vx;
  vector<double> vy;
  vector<double> vex;
  vector<double> vey;
  while (1) {
    in >>x >> y >>ey;
    ex = 0;
    if (!in.good()) break;
    vx.push_back(x);
    vy.push_back(y);
    vex.push_back(ex);
    vey.push_back(ey);
  }

  TGraphErrors *gr = new TGraphErrors(vx.size(),&vx[0],&vy[0],&vex[0],&vey[0]);
  return gr;
}

TH2D* DrawFrame(TCanvas *c1, int dir, double fr_xmin, double fr_xmax, double fr_ymin, double fr_ymax, string frtit, string frxtit, string frytit, bool logx, bool logy)
{
  if (dir==0){
    c1->cd();
  }
  else c1->cd(dir);
  static int nfr = -1;
  nfr ++;
  char frname[100];
  sprintf(frname,"fr%d",nfr);
  TH2D *fr1 = new TH2D(frname,frname,100,fr_xmin,fr_xmax,100,fr_ymin,fr_ymax);
  fr1->SetTitle(frtit.c_str());
  fr1->SetStats(0);
  fr1->SetXTitle(frxtit.c_str());
  fr1->SetYTitle(frytit.c_str());
  if (logx) gPad->SetLogx();
  if (logy) gPad->SetLogy();
  fr1->GetXaxis()->CenterTitle();
  fr1->Draw();
  return fr1;
}
void DrawGraph(TGraph* gr, int markerstyle, int markercolor, double markersize, string opt)
{
  gr->SetMarkerStyle(markerstyle);
  gr->SetMarkerColor(markercolor);
  gr->SetLineColor(markercolor);
  gr->SetLineWidth(2);
  gr->SetMarkerSize(markersize);
  if (opt=="def"){
    gr->Draw("p");
  }
  else {
    gr->Draw(opt.c_str());
  }
}
