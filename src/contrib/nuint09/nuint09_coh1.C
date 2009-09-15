//
// NuINT09 Conference, Benchmark Calculations (GENIE contribution)
//
// COH.1:
// Total coherent cross section (NC and CC) as a function of neutrino energy 
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
 /* 0 */ "nu_mu_C12",
 /* 1 */ "nu_mu_O16",
 /* 2 */ "nu_mu_Fe56"
};

void nuint09_coh1(int isample)
{
  cout << " ***** running: COH.1" << endl;

  if(isample<0 || isample >= kNSamples) return;

  const char * label = kLabel[isample];

  // get cross section graphs

  TFile fsig("../sig/splines.root","read");
  TDirectory * sig_dir = (TDirectory *) fsig.Get(label);  

  TGraph * sig_graph_cohcc = (TGraph*) sig_dir->Get("coh_cc");
  TGraph * sig_graph_cohnc = (TGraph*) sig_dir->Get("coh_nc");

  // range & spacing

  const int    ne   = 60;
  const double emin =  0.05;
  const double emax = 30.00;

  const double dloge = (TMath::Log10(emax) - TMath::Log10(emin)) / (ne-1);

  // create output stream

  ostringstream out_filename;
  out_filename << label << ".coh_1.sig_vs_Enu.data";
  ofstream out_stream(out_filename.str().c_str(), ios::out);

  // write out txt file

  out_stream << "# [" << label << "]" << endl;
  out_stream << "#  " << endl;
  out_stream << "# [COH.1]:" << endl;
  out_stream << "#  Total coherent cross section (NC and CC) as a function of neutrino energy" << endl;
  out_stream << "#  " << endl;
  out_stream << "#  Note:" << endl;
  out_stream << "#   - neutrino energy E in GeV, log spacing between Emin = " << emin << " GeV, Emax = " << emax << " GeV "  << endl;
  out_stream << "#   - cross sections in 1E-38 cm^2 " << endl;
  out_stream << "#   - for coherent scattering we quote _nuclear_ cross section " << endl;
  out_stream << "#  Columns:" << endl;
  out_stream << "#  |  Energy     |   sig(coherent; CC)   |   sig(coherent; NC)  |"  << endl;

  out_stream << setiosflags(ios::fixed) << setprecision(6);

  for(int i=0; i < ne; i++) {

     double e = TMath::Power(10., TMath::Log10(emin) + i * dloge);

     double sig_cohcc = sig_graph_cohcc->Eval(e); 
     double sig_cohnc = sig_graph_cohnc->Eval(e); 

     sig_cohcc = TMath::Max(0., sig_cohcc);
     sig_cohnc = TMath::Max(0., sig_cohnc);

     out_stream << setw(15) << e << setw(15) << sig_cohcc << setw(15) << sig_cohnc << endl;
  }

  out_stream.close();

  // visual inspection

  TCanvas * c1 = new TCanvas("c1","",20,20,500,500);
  sig_graph_cohcc->Draw("alp");
  sig_graph_cohnc->Draw("lp");
  c1->Update();
}
