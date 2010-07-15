{
// Relative normalization of SuperK event samples
//

  // output of $GENIE/src/scripts/production/misc/generate_sk_flux_histograms.C
  const char * skfluxfile = "/opt/ppd/t2k/GENIE/data/job_inputs/t2k_flux/10/sk/sk_flux_histograms.root";

  // output of $GENIE/src/support/t2k/SKNorm/gSKXSecTable.cxx (gSKxsect executable)
  const char xsecfile = "./genie_sk_xsec_table.dat";

  // binning in xsection and flux files
  double Emin =  0.0;
  double Emax = 15.0;
  double dE   =  0.050;
  int    nE   = (Emax-Emin)/dE;

  // load genie cross sections for water
  TTree xsec_water;
  xsec_water.ReadFile("genie_sk_xsec_table.dat", "E/D:xsec_numu/D:xsec_numubar/D:xsec_nue/D:xsec_nuebar/D");
  TH1D * xnumu    = new TH1D("xnumu",   "", nE, Emin, Emax);
  TH1D * xnumubar = new TH1D("xnumubar","", nE, Emin, Emax);
  TH1D * xnue     = new TH1D("xnue",    "", nE, Emin, Emax);
  xsec_water.Draw("E>>xnumu",    "xsec_numu",    "goff");
  xsec_water.Draw("E>>xnumubar", "xsec_numubar", "goff");
  xsec_water.Draw("E>>xnue",     "xsec_nue",     "goff");

  // load SuperK flux histograms
  TFile * flux = new TFile(skfluxfile, "read");
  TH1D * fnumu    = (TH1D*) flux->Get("numu_flux");
  TH1D * fnumubar = (TH1D*) flux->Get("numubar_flux");
  TH1D * fnue     = (TH1D*) flux->Get("nue_flux");

  // fill weight histograms
  double Emaxw = 8;
  int    nEw   = (Emaxw-Emin)/dE;
  TH1D * wnumubar = new TH1D("wnumubar","", nEw, Emin, Emaxw);
  TH1D * wnue     = new TH1D("wnue",    "", nEw, Emin, Emaxw);
  TH1D * wnuesig  = new TH1D("wnuesig", "", nEw, Emin, Emaxw);

  for(int i=1; i<=wnue->GetNbinsX(); i++) {
    double E = wnue->GetBinCenter(i);
    
    double xsec_numu    = xnumu    -> GetBinContent (xnumu    -> FindBin(E));
    double flux_numu    = fnumu    -> GetBinContent (fnumu    -> FindBin(E));
    double xsec_numubar = xnumubar -> GetBinContent (xnumubar -> FindBin(E));
    double flux_numubar = fnumubar -> GetBinContent (fnumubar -> FindBin(E));
    double xsec_nue     = xnue     -> GetBinContent (xnue     -> FindBin(E));
    double flux_nue     = fnue     -> GetBinContent (fnue     -> FindBin(E));
    double xsec_nuesig  = xnue     -> GetBinContent (xnue     -> FindBin(E));
    double flux_nuesig  = fnumu    -> GetBinContent (fnumu    -> FindBin(E)); //numu->nue

    double wght_numubar = (xsec_numubar * flux_numubar) / (xsec_numu * flux_numu);
    double wght_nue     = (xsec_nue     * flux_nue    ) / (xsec_numu * flux_numu);
    double wght_nuesig  = (xsec_nuesig  * flux_nuesig ) / (xsec_numu * flux_numu);

    cout << "** E = " << E;
    cout << "numu:    " << xsec_numu    << ", " << flux_numu    << endl;
    cout << "numubar: " << xsec_numubar << ", " << flux_numubar << endl;
    cout << "nue:     " << xsec_nue     << ", " << flux_nue     << endl;
    cout << "nuesig:  " << xsec_nuesig  << ", " << flux_nuesig  << endl;
    cout << "wght(numubar) = " << wght_numubar << ", wght(nue) = " << wght_nue << ", wght(nuesig) = " << wght_nuesig << endl;

    wnumubar-> SetBinContent(i, wght_numubar);
    wnue    -> SetBinContent(i, wght_nue);
    wnuesig -> SetBinContent(i, wght_nuesig);
  }

  wnuesig  -> SetLineStyle(kSolid);
  wnue     -> SetLineStyle(kDashed);
  wnumubar -> SetLineStyle(kSolid);
  wnumubar -> SetLineColor(kRed);

  TCanvas * c = new TCanvas("c","",20,20,500,500);
  TH1F * hframe = (TH1F*)c->DrawFrame(0,1E-3,Emaxw,10);
  wnuesig  -> Draw("same");
  wnue     -> Draw("same");
  wnumubar -> Draw("same");
  TLegend * legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->SetFillColor(0);
  legend->AddEntry(wnuesig,  "signal #nu_{e} / #nu_{#mu}",  "L");
  legend->AddEntry(wnue,     "#nu_{e} / #nu_{#mu}",         "L");
  legend->AddEntry(wnumubar, "#bar{#nu_{#mu}} / #nu_{#mu}", "L");
  legend->Draw();
  c->SetLogy();

  TFile fout("./skwght.root","recreate");
  c        -> Write();
  wnuesig  -> Write();
  wnue     -> Write();
  wnumubar -> Write();
  fout.Close(); 
}
