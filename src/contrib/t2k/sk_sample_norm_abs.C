//________________________________________________________________________________________
/*!

\macro   sk_sample_norm_abs.C

\brief   Calculate absolute normalization of SuperK event samples.
         
\details The program calculates the number of events of each flavour per 1E+21 POT
         per 22.5 kton water fiducial.

         N = Integral{
              (d3Flux / dE dS dI) * sig(E) * (Na/A) * rho*L * dE*dS*dI}
           = Integral{
              (d3Flux / dE dS dI) * sig(E) * (Na/A) * dE*dM*dI} (1)

         where
          - (d3Flux / dE dS dI): numu flux per energy bin, per unit area, per POT 
          - sigma(E): total numu cross section on water
          - Na: Avogadro's number
          - A: mass number for water 
          - rho: water density 
          - L: path length
          - M: mass
          - I: beam intensity (POT)
  
         The SuperK flux is given in #neutrinos per 0.05 GeV per cm2 per 1E+21 POTs.
         Input water cross sections are given in 1E-38 cm2.
         Equation (1) becomes

         N = 6.023E+23 x 1E-38 x (Mfv/A) x Ebinsz x Sum_{i} { F_{i} * sig_{i} }
  
         where 
         - N  : expected number of events for NF x 1E+21 POT
         - Mfv: fiducial volume mass
         - NF : number of 1E+21 POT worth of flux files chained together to produce 
                the flux histograms
         - Ebinsz: energy bin size (0.05 GeV)
         - F_{i): flux in bin i (in number of flux neutrinos / Ebinsz / (NF x 1E+21 POT))
         - sig_{i}: cross section evaluated at centre of bin i (1E-38 cm^2)
  
         Notes:
         Ptot = P(H1) + P(O16) = sigma(H1)*w(H1)*rho/A(H1) + sigma(O16)*w(O16)*rho/A(O16)
         rho: water density and w: mass contribution, w(H1)=2/18, w(O16)=16/18
         So: Ptot ~ sigma(H1) * (2/18) + sigma(O16) * (16/18) / 16 =>
             Ptot ~ (2*sigma(H1)+sigma(O16))/18 =>
             Ptot ~ sigma(H20)/A(H20)

\inputs
         - xsecfile:
              Neutrino - water cross section file.
              The output of  $GENIE/src/contrib/t2k/make_sk_xsec_table.C
         - skfluxfile :
              Input neutrino flux file.
              The output of $GENIE/src/scripts/production/misc/generate_sk_flux_histograms.C
         - IF :
              Number of POTs per flux simulation file used for filling-in the flux histograms 
              (typically 1E+21 POT)
         - NF :
              Number of flux files used filling-in the flux histograms 


\author  Costas Andreopoulos <costas.andreopoulos@stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created Nov 24, 2008

\cpright Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//_________________________________________________________________________________________


void sk_sample_norm_abs
(
  const char * xsecfile,      // Neutrino - water cross section file
  const char * skfluxfile,    // Input neutrino flux file.
  double       IF  = 1E+21,   // Number of POTs per flux file
  int          NF  = 500      // Number of flux files used
)

{
  // consts
  double Mfv = 2.25E+10; // want final event numbers shown for this ficucial mass (gr)
  double I0  = 1E+21;    // want final event numbers shown for this POT exposure
  double Na  = 6.023E+23;
  double A   = 18;      // gr

  // binning in xsection and flux files

  double Emin =  0.0;
  double Emax = 15.0;
  double dE   =  0.050;
  int    nE   = (Emax-Emin)/dE;

  // load genie cross sections for water

  TTree xsec_water;
  xsec_water.ReadFile(xsecfile, "E/D:xsec_numu/D:xsec_numubar/D:xsec_nue/D:xsec_nuebar/D");
  TH1D * xnumu    = new TH1D("xnumu",    "", nE, Emin, Emax);
  TH1D * xnumubar = new TH1D("xnumubar", "", nE, Emin, Emax);
  TH1D * xnue     = new TH1D("xnue",     "", nE, Emin, Emax);
  TH1D * xnuebar  = new TH1D("xnuebar",  "", nE, Emin, Emax);
  xsec_water.Draw("E>>xnumu",    "xsec_numu",    "goff");
  xsec_water.Draw("E>>xnumubar", "xsec_numubar", "goff");
  xsec_water.Draw("E>>xnue",     "xsec_nue",     "goff");
  xsec_water.Draw("E>>xnuebar",  "xsec_nuebar",  "goff");

  // load SuperK flux histograms

  TFile * flux = new TFile(skfluxfile, "read");
  TH1D * fnumu    = (TH1D*) flux->Get("numu_flux");
  TH1D * fnumubar = (TH1D*) flux->Get("numubar_flux");
  TH1D * fnue     = (TH1D*) flux->Get("nue_flux");     // instrinsic beam nue
  TH1D * fnuebar  = (TH1D*) flux->Get("nuebar_flux");  // ...

  // integrate flux x cross section and apply exposure / fiducial mass and dimensional factors

  double Nnumu    = 0;
  double Nnumubar = 0;
  double Nnue     = 0;
  double Nnuebar  = 0;
  double Nnuesig  = 0; // numu->nue

  for(int i=1; i<=fnumu->GetNbinsX(); i++) {
    double E = fnumu->GetBinCenter(i);
    
    double xsec_numu    = xnumu    -> GetBinContent (xnumu    -> FindBin(E));
    double flux_numu    = fnumu    -> GetBinContent (fnumu    -> FindBin(E));
    double xsec_numubar = xnumubar -> GetBinContent (xnumubar -> FindBin(E));
    double flux_numubar = fnumubar -> GetBinContent (fnumubar -> FindBin(E));
    double xsec_nue     = xnue     -> GetBinContent (xnue     -> FindBin(E));
    double flux_nue     = fnue     -> GetBinContent (fnue     -> FindBin(E));
    double xsec_nuebar  = xnuebar  -> GetBinContent (xnuebar  -> FindBin(E));
    double flux_nuebar  = fnuebar  -> GetBinContent (fnuebar  -> FindBin(E));

    cout << "E = " << E << " GeV";
    cout << " - numu   : sigma(H20) = " << xsec_numu    << " x1E-38 cm2, flux(@SK) = " << flux_numu    
         << " /" << dE << " GeV /" << (NF*IF) << " POT /cm2" << endl;
    cout << " - numubar: sigma(H20) = " << xsec_numubar << " x1E-38 cm2, flux(@SK) = " << flux_numubar 
         << " /" << dE << " GeV /" << (NF*IF) << " POT /cm2" << endl;
    cout << " - nue    : sigma(H20) = " << xsec_nue     << " x1E-38 cm2, flux(@SK) = "  << flux_nue    
         << " /" << dE << " GeV /" << (NF*IF) << " POT /cm2" << endl;
    cout << " - nuebar : sigma(H20) = " << xsec_nuebar  << " x1E-38 cm2, flux(@SK) = "  << flux_nuebar
         << " /" << dE << " GeV /" << (NF*IF) << " POT /cm2" << endl;

    Nnumu    += ( flux_numu    * xsec_numu    );
    Nnumubar += ( flux_numubar * xsec_numubar );
    Nnue     += ( flux_nue     * xsec_nue     );
    Nnuebar  += ( flux_nuebar  * xsec_nuebar  );
    Nnuesig  += ( flux_numu    * xsec_nue     ); // 100% numu->nue 
  }

  double f =  Na * 1E-38 * (Mfv/A) * I0 / (NF*IF);

  Nnumu    *= f;
  Nnumubar *= f;
  Nnue     *= f;
  Nnuebar  *= f;
  Nnuesig  *= f;

  // print-out results

  cout << endl;
  cout << endl;
  cout << "---------------------------------------------------------------------------" << endl;
  cout << "species                 | # of SK events per " 
       << I0  << " POT per " << Mfv/1E+9 << " kton fiducial"<< endl;
  cout << "---------------------------------------------------------------------------" << endl;
  cout << "numu                    | " << Nnumu    << endl; 
  cout << "numubar                 | " << Nnumubar << endl;
  cout << "nue(bkg)                | " << Nnue     << endl;
  cout << "nuebar(bkg)             | " << Nnuebar  << endl;
  cout << "nue(sig,100% numu->nue) | " << Nnuesig  << endl;
  cout << "---------------------------------------------------------------------------" << endl;
  cout << endl;

}
