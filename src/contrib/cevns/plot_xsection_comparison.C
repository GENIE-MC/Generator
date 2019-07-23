{
  TFile * f = new TFile("./cevns.root","read");

  TCanvas * c = new TCanvas("c","",20,20,700,500);
  c->SetLogy();
  TH1F * hf = (TH1F*) c->DrawFrame(0.003,1E-44,0.06,1E-36);
  hf->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hf->GetYaxis()->SetTitle("#sigma_{CEvNS} (GeV)");
  published_xsec_nuAr40->Draw("p");
  genie_xsec_nuAr40->Draw("lsame");

  TLegend * leg = new TLegend(0.5,0.1,0.9,0.4);
  leg->SetBorderSize(0);
  leg->AddEntry(published_xsec_nuAr40,"#nu+Ar40, digitised from Fig.1 in arXiv:1803.09183","P");
  leg->AddEntry(genie_xsec_nuAr40,"#nu+Ar40, GENIE","L");
  leg->Draw();
}
