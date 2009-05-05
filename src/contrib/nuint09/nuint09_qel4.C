//
// NuINT09 Conference, Benchmark Calculations (GENIE contribution)
//
// QEL.4:
// dSigma / dOmega_p dKE_p at E_nu = 0.5 and 1.0 GeV
//
// Inputs:
// - sample id: 0 (nu_mu+C12), 1 (nu_mu+O16), 2 (nu_mu+Fe56)
//
// Costas Andreopoulos, STFC / Rutherford Appleton Laboratory
//
#include <iomanip>

//
// consts
//
const int kNSamples       = 3;
const int kNWCur          = 1;
const int kNEnergies      = 2;
const int kNRunsPerCase   = 5;
const int kNEvtPerRun     = 100000;

const char * kLabel[kNSamples] = 
{  
 /* 0 */ "nu_mu_C12",
 /* 1 */ "nu_mu_O16",
 /* 2 */ "nu_mu_Fe56"
};

const int kRunNu[kNSamples][kNWCur][kNEnergies][kNRunsPerCase] = 
{  
 /* indices : sample ; cc/nc ; energy */
 {
  /* 0,0,0 (nu_mu C12,  CC, 0.5 GeV) */ { {200100, 200101, 200102, 200103, 200104},
  /* 0,0,1 (nu_mu C12,  CC, 1.0 GeV) */   {200200, 200201, 200202, 200203, 200204} }
 },
 {
  /* 1,0,0 (nu_mu O16,  CC, 0.5 GeV) */ { {210100, 210101, 210102, 210103, 210104},
  /* 1,0,1 (nu_mu O16,  CC, 1.0 GeV) */   {210200, 210201, 210202, 210203, 210204} }
 },
 {
  /* 2,0,0 (nu_mu Fe56, CC, 0.5 GeV) */ { {220100, 220101, 220102, 220103, 220104},
  /* 2,0,1 (nu_mu Fe56, CC, 1.0 GeV) */   {220200, 220201, 220202, 220203, 220204} }
 }
};
int kZ[3] =
{
 // Z for nuclear target at each sample
 /* 0 */  6,
 /* 1 */  8,
 /* 2 */ 26 
};
int kA[3] =
{
 // A for nuclear target at each sample
 /* 0 */  12,
 /* 1 */  16,
 /* 2 */  56 
};

void nuint09_qel4(int isample)
{
  cout << " ***** running: QEL.4" << endl;

  if(isample<0 || isample >= kNSamples) return;

  const char * label = kLabel[isample];

  int A = kA[isample];
  int Z = kZ[isample];
  int N = A-Z;

  // get cross section graphs

  TFile fsig("../sig/splines.root","read");
  TDirectory * sig_dir = (TDirectory *) fsig.Get(label);  

  TGraph * sig_graph_qelcc = (TGraph*) sig_dir->Get("qel_cc_n");

  // range & spacing

  const int    nKEp    = 60;
  const double KEpmin  =  0.01;
  const double KEpmax  =  1.50;

  const int    ncosth   = 30;
  const double costhmin = -1;
  const double costhmax = +1;

  // create output stream

  ostringstream out_filename;
  out_filename << label << ".qel_4.d2sig_dKEpdOmega.data";
  ofstream out_stream(out_filename.str().c_str(), ios::out);

  // write out txt file

  out_stream << "# [" << label << "]" << endl;
  out_stream << "#  " << endl;
  out_stream << "# [QEL.4]:" << endl;
  out_stream << "#  dSigma / dOmega_p dKE_p at E_nu = 0.5 and 1.0 GeV" << endl;
  out_stream << "#  " << endl;
  out_stream << "#  Note:" << endl;
  out_stream << "#   - proton energies are _kinetic_ energies " << endl;
  out_stream << "#   - proton kinetic energy KE in GeV, between KEmin = " << KEpmin << " GeV, KEmax = " << KEpmax << " GeV "  << endl;
  out_stream << "#   - cross sections in 1E-38 cm^2 /GeV /sterad" << endl;
  out_stream << "#   - quoted cross section is nuclear cross section per nucleon contributing in the scattering (eg only neutrons for nu_mu QELCC)" << endl;
  out_stream << "#  Columns:" << endl;
  out_stream << "#  |  KE(p)  |  cos(theta_p)  | sig(QELCC; Ev=0.5GeV)  | sig(QELCC; Ev=1.0GeV) | " << endl;

  out_stream << setiosflags(ios::fixed) << setprecision(6);

  //
  // load event data
  //

  TChain * chain = new TChain("gst");

  // loop over CC/NC cases
  for(int iwkcur=0; iwkcur<kNWCur; iwkcur++) {
    // loop over energies
    for(int ie=0; ie<kNEnergies; ie++) {
       // loop over runs for current case
       for(int ir=0; ir<kNRunsPerCase; ir++) {
	  // build filename
	  ostringstream filename;
          int run_number = kRunNu[isample][iwkcur][ie][ir];
          filename << "../gst/gntp." << run_number << ".gst.root";
          // add to chain
          cout << "Adding " << filename.str() << " to event chain" << endl;
          chain->Add(filename.str().c_str());
       }
    }
  }

  //
  // get QELCC cross sections at given energies for normalization purposes
  //
  double sig_qelcc_0500MeV = sig_graph_qelcc->Eval(0.5) / N; 
  double sig_qelcc_1000MeV = sig_graph_qelcc->Eval(1.0) / N; 

  //
  // book histograms
  //
  TH2D * hst_d2sig_dKEpdOmg_qelcc_0500MeV =  new TH2D("hst_d2sig_dKEpdOmg_qelcc_0500MeV",
           "dsig/d2KEpdOmega, nu_mu QEL CC, Enu=0.5 GeV", nKEp, KEpmin, KEpmax, ncosth, costhmin, costhmax);
  TH2D * hst_d2sig_dKEpdOmg_qelcc_1000MeV = new TH2D("hst_d2sig_dKEpdOmg_qelcc_1000MeV",
           "dsig/d2KEpdOmega, nu_mu QEL CC, Enu=1.0 GeV", nKEp, KEpmin, KEpmax, ncosth, costhmin, costhmax);

  //
  // fill histograms
  //
  chain->Draw("pzi/sqrt(pzi*pzi+pyi*pyi+pxi*pxi):(Ei-0.938272)>>hst_d2sig_dKEpdOmg_qelcc_0500MeV","qel&&cc&&Ev>0.49&&Ev<0.51&&pdgi==2212","GOFF");
  chain->Draw("pzi/sqrt(pzi*pzi+pyi*pyi+pxi*pxi):(Ei-0.938272)>>hst_d2sig_dKEpdOmg_qelcc_1000MeV","qel&&cc&&Ev>0.99&&Ev<1.01&&pdgi==2212","GOFF");

  //
  // normalize
  //
  double norm_qelcc_0500MeV = hst_d2sig_dKEpdOmg_qelcc_0500MeV -> Integral("width") * 2*TMath::Pi() / sig_qelcc_0500MeV;
  double norm_qelcc_1000MeV = hst_d2sig_dKEpdOmg_qelcc_1000MeV -> Integral("width") * 2*TMath::Pi() / sig_qelcc_1000MeV;

  if (norm_qelcc_0500MeV > 0) hst_d2sig_dKEpdOmg_qelcc_0500MeV -> Scale(1./norm_qelcc_0500MeV);
  if (norm_qelcc_1000MeV > 0) hst_d2sig_dKEpdOmg_qelcc_1000MeV -> Scale(1./norm_qelcc_1000MeV);

  for(int i = 1; i <= hst_d2sig_dKEpdOmg_qelcc_1000MeV->GetNbinsX(); i++) {
     for(int j = 1; j <= hst_d2sig_dKEpdOmg_qelcc_1000MeV->GetNbinsY(); j++) {

        double KEp    = hst_d2sig_dKEpdOmg_qelcc_1000MeV -> GetXaxis() -> GetBinCenter(i);
        double costhp = hst_d2sig_dKEpdOmg_qelcc_1000MeV -> GetYaxis() -> GetBinCenter(j);

        double d2sig_dKEpdOmg_qelcc_0500MeV = hst_d2sig_dKEpdOmg_qelcc_0500MeV -> GetBinContent(i,j);
        double d2sig_dKEpdOmg_qelcc_1000MeV = hst_d2sig_dKEpdOmg_qelcc_1000MeV -> GetBinContent(i,j);

        d2sig_dKEpdOmg_qelcc_0500MeV = TMath::Max(0., d2sig_dKEpdOmg_qelcc_0500MeV);
        d2sig_dKEpdOmg_qelcc_1000MeV = TMath::Max(0., d2sig_dKEpdOmg_qelcc_1000MeV);

        out_stream << setw(15) << KEp 
                   << setw(15) << costhp
                   << setw(15) << d2sig_dKEpdOmg_qelcc_0500MeV
                   << setw(15) << d2sig_dKEpdOmg_qelcc_1000MeV
                   << endl;
    }//costh
  }//E

  out_stream.close();

  delete chain;
}
