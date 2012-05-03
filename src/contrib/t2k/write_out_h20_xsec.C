//
// read the output of make_sk_xsec_table.C used by SuperK for MC job normalization, 
// and convert to ROOT format
//

{
  double Emin =  0.0;
  double Emax = 15.0;
  double dE   =  0.050;
  int    nE   = (Emax-Emin)/dE;

  TFile f("./xsec_H20.root","recreate");

  TH1D * xsec_numuH20    = new TH1D("xsec_numuH20",    "", nE, Emin, Emax);
  TH1D * xsec_numubarH20 = new TH1D("xsec_numubarH20", "", nE, Emin, Emax);
  TH1D * xsec_nueH20     = new TH1D("xsec_nueH20",     "", nE, Emin, Emax);
  TH1D * xsec_nuebarH20  = new TH1D("xsec_nuebarH20",  "", nE, Emin, Emax);

  TTree xsec_water;
  const char * xsec_water_file = "./genie_sk_xsec_table.dat";
  xsec_water.ReadFile(xsec_water_file, "E/D:xsec_numu/D:xsec_numubar/D:xsec_nue/D:xsec_nuebar/D");

  xsec_water.Draw("E>>xsec_numuH20",    "xsec_numu",    "goff");
  xsec_water.Draw("E>>xsec_numubarH20", "xsec_numubar", "goff");
  xsec_water.Draw("E>>xsec_nueH20",     "xsec_nue",     "goff");
  xsec_water.Draw("E>>xsec_nuebarH20",  "xsec_nuebar",  "goff");

  f.Write();
  f.Close();
}
