{
// Absolute normalization of SuperK numu event sample
//
// N = Integral{
//      (d3Flux / dE dS dI) * sig(E) * (Na/A) * rho*L * dE*dS*dI}
//   = Integral{
//      (d3Flux / dE dS dI) * sig(E) * (Na/A) * dE*dM*dI}
// where
//   (d3Flux / dE dS dI): numu flux per energy bin, per unit area, per POT 
//   sigma(E): total numu cross section on water
//   Na: Avogadro's number
//   A: mass number for water 
//   rho: water density 
//   L: path length
//   M: mass
//   I: beam intensity (POT)
//
// SK flux is given in #neutrinos per 0.05 GeV per cm2 per 1E+21 POTs
// Input water cross sections are given in 1E-28 cm2
//
// N = 6.023E+23 x 1E-38 x NF x (Mfv/A) x Ebinsz x Sum_{i} { F_{i} * sig_{i} }
//
// where 
//  Mfv: fiducial volume mass
//  NF : number of 1E+21 POT worth of flux files chained together to produce the flux histogram
//  Ebinsz: energy bin size (0.05 GeV)
//  F_{i): flux in bin i
//  sig_{i}: cross section evaluated at centre of bin i
//
// notes:
//  Ptot = P(H1) + P(O16) = sigma(H1)*w(H1)*rho/A(H1) + sigma(O16)*w(O16)*rho/A(O16)
//  rho: water density and w: mass contribution, w(H1)=2/18, w(O16)=16/18
//  So: Ptot ~ sigma(H1) * (2/18) + sigma(O16) * (16/18) / 16 =>
//      Ptot ~ (2*sigma(H1)+sigma(O16))/18 =>
//      Ptot ~ sigma(H20)/A(H20)
//

  // output of $GENIE/src/scripts/production/misc/generate_sk_flux_histograms.C
  const char * skfluxfile = "/opt/ppd/t2k/GENIE/data/job_inputs/t2k_flux/10/sk/sk_flux_histograms.root";

  // output of $GENIE/src/support/t2k/SKNorm/gSKXSecTable.cxx (gSKxsect executable)
  const char xsecfile = "./genie_sk_xsec_table.dat";

  // consts
  double Na  = 6.023E+23;
  double Mfv = 1E+10; // gr
  double A   = 18; // gr
  double Io  = 1E+23; // pot (all flux files)
  double If  = 1E+21; // pot (per flux file)

  // binning in xsection and flux files
  double Emin =  0.0;
  double Emax = 15.0;
  double dE   =  0.050;
  int    nE   = (Emax-Emin)/dE;

  // load genie cross sections for water
  TTree xsec_water;
  xsec_water.ReadFile("genie_sk_xsec_table.dat", "E/D:xsec_numu/D:xsec_numubar/D:xsec_nue/D:xsec_nuebar/D");
  TH1D * xnumu = new TH1D("xnumu", "", nE, Emin, Emax);
  xsec_water.Draw("E>>xnumu", "xsec_numu", "goff");

  // load SuperK flux histograms
  TFile * flux = new TFile(skfluxfile, "read");
  TH1D * fnumu = (TH1D*) flux->Get("numu_flux");

  double N = 0;

  for(int i=1; i<=fnumu->GetNbinsX(); i++) {
    double E = fnumu->GetBinCenter(i);
    
    double xsec_numu = xnumu -> GetBinContent (xnumu -> FindBin(E));
    double flux_numu = fnumu -> GetBinContent (fnumu -> FindBin(E));

    cout << "E = " << E << ", sigma(numu) = " << xsec_numu << ", flux(numu) = " << flux_numu << endl;

    N += (flux_numu * xsec_numu);
    cout << N << endl;
  }

  N *= ( Na * 1E-38 * (Io/If) * (Mfv/A) * dE );

  cout << "n = " << N << " SK numu events "
       << "per " << Io << " POT per " << Mfv/1E+9 << " kton fiducial"<< endl;
}
