TFile* inf;
TNtuple* nt;
TCanvas* c1, * c1l, * c2;
TGraph* xsg, * psg, * rsg, * nsg;
TH2* h, * rh;

void drawXsec(const Char_t* infn="VLExsecNT.root") {
   
   gStyle->SetOptStat(0);
   
   inf = TFile::Open(infn);
   if (inf==0 || inf->IsZombie()) {
      Error("drawXsec","Could not open [%s]",infn);
      return;
   }
   
   nt = dynamic_cast<TNtuple*>(inf->Get("nt"));
   if (nt==0) {
      Error("drawXsec","Could not get nt from [%s]",infn);
      return;
   }
   
   // draw the cross sections
   c1 = new TCanvas("c1", "c1: xsec", 700, 500);
   c1->cd();
   
   h = new TH2F("h","Inverse Beta Decay Cross Section"
      ";E_{#nu} (MeV);#sigma #times 10^{-41} cm^{2}",
      110, 0, 220,
      100, 1e-3, 350);
   h->GetXaxis()->SetNdivisions(505);
   
   const Long64_t xents = nt->Draw("xsec:Ev","","goff");
   xsg = new TGraph(xents,nt->GetV2(),nt->GetV1());
   xsg->SetName("xsg");
   xsg->SetMarkerStyle(24);
   xsg->SetMarkerSize(0.8);
   
   const Long64_t xnents = nt->Draw("xsnun:Ev","","goff");
   nsg = new TGraph(xnents,nt->GetV2(),nt->GetV1());
   nsg->SetName("nsg");
   nsg->SetMarkerStyle(27);
   nsg->SetMarkerColor(kRed+1);
   
   const Long64_t pents = nt->Draw("xspaper:Ev","","goff");
   psg = new TGraph(pents,nt->GetV2(),nt->GetV1());
   psg->SetName("psg");
   psg->SetLineColor(kAzure-8);

   h->Draw("axis");
   psg->Draw("c");
   xsg->Draw("p");
   nsg->Draw("p");

   TLegend* xl = new TLegend(.188,.561,.560,.845);
   xl->SetBorderSize(0);
   xl->SetFillColor(0);
   xl->AddEntry(psg,"Strumia (#bar{#nu_{e}}p#rightarrowne^{+})","l");
   xl->AddEntry(xsg,"VLE GENIE (#bar{#nu_{e}}p#rightarrowne^{+})","p");
   xl->AddEntry(nsg,"VLE GENIE (#nu_{e}n#rightarrowpe^{-})","p");
   xl->Draw();
   

   
   // draw the cross sections on log scale
   c1l = new TCanvas("c1l", "c1l: xsec (logy)", 700, 500);
   c1l->cd();
   c1l->SetLogy();
   h->Draw("axis");
   psg->Draw("c");
   xsg->Draw("p");
   nsg->Draw("p");
   
   TLegend* xll = new TLegend(.453,.360,.825,.644);
   xll->SetBorderSize(0);
   xll->SetFillColor(0);
   xll->AddEntry(psg,"Strumia (#bar{#nu_{e}}p#rightarrowne^{+})","l");
   xll->AddEntry(xsg,"VLE GENIE (#bar{#nu_{e}}p#rightarrowne^{+})","p");
   xll->AddEntry(nsg,"VLE GENIE (#nu_{e}n#rightarrowpe^{-})","p");
   xll->Draw();
   
   // draw ratio: xsec / strumia for nu_e_bar + p --> n + e+
   c2 = new TCanvas("c2", "c2: xsec ratio", 700, 500);;
   c2->cd();
   c2->SetGridy();
   
   rsg = new TGraph(*xsg);
   rsg->SetName("rsg");
   rsg->SetMarkerStyle(21);
   rsg->SetMarkerSize(1);
   rsg->SetMarkerColor(kAzure-8);
   const Double_t* py = psg->GetY();
   Double_t*       ry = rsg->GetY();
   const Int_t     np = rsg->GetN();
   for (Int_t i=0; i<np; i++, ry++, py++) {
      if (*py > 1e-10) *ry /= *py;
   }
   
   rh = new TH2F("rh","IBD X-sec Ratio: GENIE / Strumia Paper"
      ";E_{#nu} (MeV);#sigma(GENIE) / #sigma(Strumia Paper)",
      110, 0, 220,
      100, 0.85, 1.05);
   rh->GetXaxis()->SetNdivisions(505);
   rh->GetYaxis()->SetTitleOffset(1.35);
   
   rh->Draw("axis");
   rsg->Draw("p");
   
   TLatex* rl = new TLatex(126,1.023,"#bar{#nu_{e}}p#rightarrowne^{+}");
   rl->SetTextColor(kAzure-8);
   rl->Draw();
   
}
