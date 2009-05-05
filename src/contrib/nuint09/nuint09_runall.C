//
// NuINT09 Conference, Benchmark Calculations (GENIE contribution)
//
// Run all macros
//
// Costas Andreopoulos, STFC / Rutherford Appleton Laboratory
//


//
// sample IDs:
//
// 0 <--> nu_mu + C12
// 1 <--> nu_mu + O16
// 2 <--> nu_mu + Fe56
//


{

  //
  // process nu_mu+C12 samples
  //

  gROOT->ProcessLine(".x nuint09_coh1.C(0)");
  gROOT->ProcessLine(".x nuint09_coh2.C(0)");
  gROOT->ProcessLine(".x nuint09_coh3.C(0)");
  gROOT->ProcessLine(".x nuint09_qel1.C(0)");
  gROOT->ProcessLine(".x nuint09_qel2.C(0)");
  gROOT->ProcessLine(".x nuint09_qel3.C(0)");
  gROOT->ProcessLine(".x nuint09_qel4.C(0)");
  gROOT->ProcessLine(".x nuint09_qel5.C(0)");
  gROOT->ProcessLine(".x nuint09_1pi1.C(0,0)");
  gROOT->ProcessLine(".x nuint09_1pi2.C(0,0)");
  gROOT->ProcessLine(".x nuint09_1pi3.C(0,0)");
  gROOT->ProcessLine(".x nuint09_1pi4.C(0,0)");
  gROOT->ProcessLine(".x nuint09_1pi1.C(0,1)");
  gROOT->ProcessLine(".x nuint09_1pi2.C(0,1)");
  gROOT->ProcessLine(".x nuint09_1pi3.C(0,1)");
  gROOT->ProcessLine(".x nuint09_1pi4.C(0,1)");
  gROOT->ProcessLine(".x nuint09_1pi1.C(0,2)");
  gROOT->ProcessLine(".x nuint09_1pi2.C(0,2)");
  gROOT->ProcessLine(".x nuint09_1pi3.C(0,2)");
  gROOT->ProcessLine(".x nuint09_1pi4.C(0,2)");

}
