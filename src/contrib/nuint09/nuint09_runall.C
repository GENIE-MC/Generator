//
// NuINT09 Conference, Benchmark Calculations (GENIE contribution)
// This macro runs all other macros
//
// Costas Andreopoulos, STFC / Rutherford Appleton Laboratory
// ................................................................
//
// sample:     0 --> nu_mu + C12
//             1 --> nu_mu + O16
//             2 --> nu_mu + Fe56
// 1pi source: 0 --> all
//             1 --> P33(1232) only
//             2 --> all resonances only
// stage:      0 --> primary hadronic state (no FSI)
//             1 --> final state

{
  //
  // process nu_mu+C12 samples
  //
  cout << "Processing the nu_mu+C12 samples" << endl;

  // coherent
  gROOT->ProcessLine(".x nuint09_coh1.C(0)");
  gROOT->ProcessLine(".x nuint09_coh2.C(0)");
  gROOT->ProcessLine(".x nuint09_coh3.C(0)");

  // quasi-elastic
  gROOT->ProcessLine(".x nuint09_qel1.C(0)");
  gROOT->ProcessLine(".x nuint09_qel2.C(0,0)"); // -FSI
  gROOT->ProcessLine(".x nuint09_qel2.C(0,1)"); // +FSI
  gROOT->ProcessLine(".x nuint09_qel3.C(0)");
  gROOT->ProcessLine(".x nuint09_qel4.C(0)");
  gROOT->ProcessLine(".x nuint09_qel5.C(0)");
  gROOT->ProcessLine(".x nuint09_qel6.C(0)");

  // final state {X pi+}
  gROOT->ProcessLine(".x nuint09_1pi1.C(0,0,1)"); // all sources
  gROOT->ProcessLine(".x nuint09_1pi2.C(0,0,1)");
  gROOT->ProcessLine(".x nuint09_1pi3.C(0,0,1)");
  gROOT->ProcessLine(".x nuint09_1pi4.C(0,0,1)");
  gROOT->ProcessLine(".x nuint09_1pi1.C(0,1,1)"); // P33(1232) only
  gROOT->ProcessLine(".x nuint09_1pi2.C(0,1,1)");
  gROOT->ProcessLine(".x nuint09_1pi3.C(0,1,1)");
  gROOT->ProcessLine(".x nuint09_1pi4.C(0,1,1)");
  gROOT->ProcessLine(".x nuint09_1pi1.C(0,2,1)"); // all resonances only
  gROOT->ProcessLine(".x nuint09_1pi2.C(0,2,1)");
  gROOT->ProcessLine(".x nuint09_1pi3.C(0,2,1)");
  gROOT->ProcessLine(".x nuint09_1pi4.C(0,2,1)");

  // primary {X pi+}
  gROOT->ProcessLine(".x nuint09_1pi1.C(0,0,0)"); // all sources
  gROOT->ProcessLine(".x nuint09_1pi2.C(0,0,0)");
  gROOT->ProcessLine(".x nuint09_1pi3.C(0,0,0)");
  gROOT->ProcessLine(".x nuint09_1pi4.C(0,0,0)");
  gROOT->ProcessLine(".x nuint09_1pi1.C(0,1,0)"); // P33(1232) only
  gROOT->ProcessLine(".x nuint09_1pi2.C(0,1,0)");
  gROOT->ProcessLine(".x nuint09_1pi3.C(0,1,0)");
  gROOT->ProcessLine(".x nuint09_1pi4.C(0,1,0)");
  gROOT->ProcessLine(".x nuint09_1pi1.C(0,2,0)"); // all resonances only
  gROOT->ProcessLine(".x nuint09_1pi2.C(0,2,0)");
  gROOT->ProcessLine(".x nuint09_1pi3.C(0,2,0)");
  gROOT->ProcessLine(".x nuint09_1pi4.C(0,2,0)");

  //
  // process nu_mu+O16 samples
  //
  cout << "Processing the nu_mu+O16 samples" << endl;

  // quasi-elastic
  gROOT->ProcessLine(".x nuint09_qel1.C(1)");
  gROOT->ProcessLine(".x nuint09_qel2.C(1,0)"); // -FSI
  gROOT->ProcessLine(".x nuint09_qel2.C(1,1");  // +FSI
  gROOT->ProcessLine(".x nuint09_qel3.C(1)");
  gROOT->ProcessLine(".x nuint09_qel4.C(1)");
  gROOT->ProcessLine(".x nuint09_qel5.C(1)");

  //
  // process nu_mu+Fe56 samples
  //
  cout << "Processing the nu_mu+Fe56 samples" << endl;

  // quasi-elastic
  gROOT->ProcessLine(".x nuint09_qel1.C(2)");
  gROOT->ProcessLine(".x nuint09_qel2.C(2,0)"); // -FSI
  gROOT->ProcessLine(".x nuint09_qel2.C(2,1)"); // +FSI
  gROOT->ProcessLine(".x nuint09_qel3.C(2)");
  gROOT->ProcessLine(".x nuint09_qel4.C(2)");
  gROOT->ProcessLine(".x nuint09_qel5.C(2)");
}
