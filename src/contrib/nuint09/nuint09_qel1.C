//
// NuINT09 Conference, Benchmark Calculations (GENIE contribution)
//
// QEL.1:
// QEL total cross section as a function of energy
//
// Inputs:
// - sample id: 0 (nu_mu+C12), 1 (nu_mu+O16), 2 (nu_mu+Fe56)
//
// Costas Andreopoulos, STFC / Rutherford Appleton Laboratory
//
#include <iomanip>

const int kNSamples = 3;

const char * kLabel[kNSamples] = 
{  
 // available samples
 /* 0 */ "nu_mu_C12",
 /* 1 */ "nu_mu_O16",
 /* 2 */ "nu_mu_Fe56"
};
int kZ[kNSamples] = 
{  
 // Z for nuclear target at each sample
 /* 0 */  6,
 /* 1 */  8,
 /* 2 */ 26
};
int kA[kNSamples] = 
{  
 // A for nuclear target at each sample
 /* 0 */  12,
 /* 1 */  16,
 /* 2 */  56
};

void nuint09_qel1(int isample)
{
  cout << " ***** running: QEL.1" << endl;

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

  const int    ne   = 60;
  const double emin =  0.05;
  const double emax = 30.00;

  const double dloge = (TMath::Log10(emax) - TMath::Log10(emin)) / (ne-1);

  // create output stream

  ostringstream out_filename;
  out_filename << label << ".qel_1.sig_vs_Enu.data";
  ofstream out_stream(out_filename.str().c_str(), ios::out);

  // write out txt file

  out_stream << "# [" << label << "]" << endl;
  out_stream << "#  " << endl;
  out_stream << "# [QEL.1]:" << endl;
  out_stream << "#  QEL total cross section as a function of neutrino energy" << endl;
  out_stream << "#  " << endl;
  out_stream << "#  Note:" << endl;
  out_stream << "#   - neutrino energy E in GeV, log spacing between Emin = " << emin << " GeV, Emax = " << emax << " GeV "  << endl;
  out_stream << "#   - cross sections in 1E-38 cm^2 " << endl;
  out_stream << "#   - quoted cross section is nuclear cross section per nucleon contributing in the scattering (eg only neutrons for nu_mu QELCC)" << endl;
  out_stream << "#  Columns:" << endl;
  out_stream << "#  |  Energy     |   sig(nu_mu + n[bound] -> mu- + p)   |  "  << endl;

  out_stream << setiosflags(ios::fixed) << setprecision(6);

  for(int i=0; i < ne; i++) {

     double e = TMath::Power(10., TMath::Log10(emin) + i * dloge);

     double sig_qelcc = sig_graph_qelcc->Eval(e) / N; 
     sig_qelcc = TMath::Max(0., sig_qelcc);

     out_stream << setw(15) << e << setw(15) << sig_qelcc << endl;
  }

  out_stream.close();

  // visual inspection
  TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  sig_graph_qelcc->Draw("alp");
  c1->Update();
}
