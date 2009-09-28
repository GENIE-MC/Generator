//
// NuINT09 Conference, Benchmark Calculations (GENIE contribution)
//
// QEL.6:
// CCQE-like cross section (1 nucleon + 0 pions after FSI) at E_nu= 0.5, 1.0 and 1.5 GeV as a function of the nucleon kinetic energy.
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
const int kNEnergies      = 3;
const int kNRunsPerCase   = 5;
const int kNEvtPerRun     = 100000;

const char * kLabel[kNSamples] = 
{  
 // available samples
 /* 0 */ "nu_mu_C12",
 /* 1 */ "nu_mu_O16",
 /* 2 */ "nu_mu_Fe56"
};

const int kRunNuQEL6[kNSamples][kNWCur][kNEnergies][kNRunsPerCase] =
{
 /* indices : sample ; cc/nc ; energy */
 {
  /* 0,0,0 (nu_mu C12,  CC, 0.5 GeV) */ { {900100, 900101, 900102, 900103, 900104},
  /* 0,0,1 (nu_mu C12,  CC, 1.0 GeV) */   {900200, 900201, 900202, 900203, 900204},
  /* 0,0,2 (nu_mu C12,  CC, 1.5 GeV) */   {900300, 900301, 900302, 900303, 900304} }
 },
 {
  /* 1,0,0 (nu_mu O16,  CC, 0.5 GeV) */ { {910100, 910101, 910102, 910103, 910104},
  /* 1,0,1 (nu_mu O16,  CC, 1.0 GeV) */   {910200, 910201, 910202, 910203, 910204},
  /* 1,0,2 (nu_mu O16,  CC, 1.5 GeV) */   {910300, 910301, 910302, 910303, 910304} }
 },
 {
  /* 2,0,0 (nu_mu Fe56, CC, 0.5 GeV) */ { {920100, 920101, 920102, 920103, 920104},
  /* 2,0,1 (nu_mu Fe56, CC, 1.0 GeV) */   {920200, 920201, 920202, 920203, 920204},
  /* 2,0,2 (nu_mu Fe56, CC, 1.5 GeV) */   {920300, 920301, 920302, 920303, 920304} }
 }
};

int kA[kNSamples] = 
{  
 // A for nuclear target at each sample
 /* 0 */  12,
 /* 1 */  16,
 /* 2 */  56
};

void nuint09_qel6(int isample)
{
  cout << " ***** running: QEL.6" << endl;

  if(isample<0 || isample >= kNSamples) return;

  const char * label = kLabel[isample];
  int A = kA[isample];

  // get cross section graph  

  TFile fsig("../sig/splines.root","read");
  TDirectory * sig_dir = (TDirectory *) fsig.Get(label);  

  TGraph * sig_totcc = (TGraph*) sig_dir->Get("tot_cc");

  // range & spacing

  const int    nKEnuc   =  100;
  const double KEnucmin =  0.00;
  const double KEnucmax =  1.30;

  // create output stream

  ostringstream out_filename;
  out_filename << label << ".qel_6.dsigCCQElike_dKEp.data";

  ofstream out_stream(out_filename.str().c_str(), ios::out);

  // write out txt file

  out_stream << "# [" << label << "]" << endl;
  out_stream << "#  " << endl;
  out_stream << "# [QEL.6]:" << endl;
  out_stream << "#  CCQE-like cross section (1 nucleon + 0 pions after FSI) at E_nu = 0.5, 1.0 and 1.5 GeV as a function of the nucleon kinetic energy" << endl;
  out_stream << "#  " << endl;
  out_stream << "#  Note:" << endl;
  out_stream << "#   - nucleon kinetic energy KE in GeV, linear spacing between KEmin = " << KEnucmin << " GeV, KEmax = " << KEnucmax << " GeV " << endl;
  out_stream << "#   - cross sections in 1E-38 cm^2 / GeV" << endl;
  out_stream << "#   - quoted cross section is nuclear cross section divided with number of nucleons A" << endl;
  out_stream << "#  Columns:" << endl;
  out_stream << "#  |  KE(proton)  | dsig(numu CCQE-like; Enu = 0.5 GeV) | dsig(numu CCQE-like; Enu = 1.0 GeV) | dsig(numu CCQE-like; Enu = 1.5 GeV)  | "  << endl;
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
          int run_number = kRunNuQEL6[isample][iwkcur][ie][ir];
          filename << "../gst/gntp." << run_number << ".gst.root";
          // add to chain
          cout << "Adding " << filename.str() << " to event chain" << endl;
          chain->Add(filename.str().c_str());
       }
    }
  } 

  // 
  // get CC cross sections at given energies for normalization purposes
  //
  double sig_totcc_0500MeV = sig_totcc -> Eval (0.5) / A;
  double sig_totcc_1000MeV = sig_totcc -> Eval (1.0) / A;
  double sig_totcc_1500MeV = sig_totcc -> Eval (1.5) / A;
  
  //
  // book histograms
  //
  TH1D * hst_dsig_dKEnuc_0500MeV = new TH1D("hst_dsig_dKEnuc_0500MeV","dsig/dKEp, numu CCQE-like after FSI, Enu=0.5 GeV", nKEnuc, KEnucmin, KEnucmax);
  TH1D * hst_dsig_dKEnuc_1000MeV = new TH1D("hst_dsig_dKEnuc_1000MeV","dsig/dKEp, numu CCQE-like after FSI, Enu=1.0 GeV", nKEnuc, KEnucmin, KEnucmax);
  TH1D * hst_dsig_dKEnuc_1500MeV = new TH1D("hst_dsig_dKEnuc_1500MeV","dsig/dKEp, numu CCQE-like after FSI, Enu=1.5 GeV", nKEnuc, KEnucmin, KEnucmax);

  //
  // fill histograms
  //
  chain->Draw("(Ef-0.938)>>hst_dsig_dKEnuc_0500MeV","cc&&Ev>0.49&&Ev<0.51&&nfp==1&&nfn==0&&nfpip==1&&nfpim==0&&nfpi0==0&&nfother==0","GOFF");
  chain->Draw("(Ef-0.938)>>hst_dsig_dKEnuc_1000MeV","cc&&Ev>0.99&&Ev<1.01&&nfp==1&&nfn==0&&nfpip==1&&nfpim==0&&nfpi0==0&&nfother==0","GOFF");
  chain->Draw("(Ef-0.938)>>hst_dsig_dKEnuc_1500MeV","cc&&Ev>1.49&&Ev<1.51&&nfp==1&&nfn==0&&nfpip==1&&nfpim==0&&nfpi0==0&&nfother==0","GOFF");
                
  //
  // normalize
  //
  double norm_0500MeV = hst_dsig_dKEnuc_0500MeV -> Integral("width") / sig_totcc_0500MeV;
  double norm_1000MeV = hst_dsig_dKEnuc_1000MeV -> Integral("width") / sig_totcc_1000MeV;
  double norm_1500MeV = hst_dsig_dKEnuc_1500MeV -> Integral("width") / sig_totcc_1500MeV;
  
  if (norm_0500MeV > 0) hst_dsig_dKEnuc_0500MeV -> Scale(1./norm_0500MeV);
  if (norm_1000MeV > 0) hst_dsig_dKEnuc_1000MeV -> Scale(1./norm_1000MeV);
  if (norm_1500MeV > 0) hst_dsig_dKEnuc_1500MeV -> Scale(1./norm_1500MeV);

  //
  // write-out
  //
  for(int i=1; i <= hst_dsig_dKEnuc_1000MeV->GetNbinsX(); i++) {

     double KEnuc = hst_dsig_dKEnuc_1000MeV -> GetBinCenter(i);

     double dsig_dKEnuc_0500MeV = hst_dsig_dKEnuc_0500MeV -> GetBinContent(i);
     double dsig_dKEnuc_1000MeV = hst_dsig_dKEnuc_1000MeV -> GetBinContent(i);
     double dsig_dKEnuc_1500MeV = hst_dsig_dKEnuc_1500MeV -> GetBinContent(i);

     dsig_dKEnuc_0500MeV = TMath::Max(0., dsig_dKEnuc_0500MeV);
     dsig_dKEnuc_1000MeV = TMath::Max(0., dsig_dKEnuc_1000MeV);
     dsig_dKEnuc_1500MeV = TMath::Max(0., dsig_dKEnuc_1500MeV);

     out_stream << setw(15) << KEnuc 
                << setw(15) << dsig_dKEnuc_0500MeV
                << setw(15) << dsig_dKEnuc_1000MeV
                << setw(15) << dsig_dKEnuc_1500MeV
                << endl;
  }

  out_stream.close();

  delete chain;
  fsig.Close();
}
