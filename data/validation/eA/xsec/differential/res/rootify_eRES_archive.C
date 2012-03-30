//
// ROOT-ify data files downloaded from the 
// Resonance Electron-Nucleus Scattering Archive maintained by JLAB
//
// C.Andreopoulos 
//

#include <fstream>
#include <string>
#include <multimap>
#include <iomanip>

//#define _show_debug_mesg_

using namespace std;

int rootify_eRES_archive(void)
{
  const int nfiles = 25;

  string filename[nfiles] = 
  {
    "eRES_1H_JLab_e00_002.data",
    "eRES_1H_JLab_e94_110.data",
    "eRES_1H_JLab_ioanna.data",
    "eRES_1H_SLAC_e133.data",
    "eRES_1H_SLAC_e140.data",
    "eRES_1H_SLAC_e140x.data",
    "eRES_1H_SLAC_e49a10.data",
    "eRES_1H_SLAC_e49a6.data",
    "eRES_1H_SLAC_e49b.data",
    "eRES_1H_SLAC_e61.data",
    "eRES_1H_SLAC_e87.data",
    "eRES_1H_SLAC_e891.data",
    "eRES_1H_SLAC_e8920.data",
    "eRES_1H_SLAC_ne11.data",
    "eRES_1H_SLAC_onen1haf.data",
    "eRES_2H_JLab_ioanna.data",
    "eRES_2H_SLAC_e133.data",
    "eRES_2H_SLAC_e140x.data",
    "eRES_2H_SLAC_e49a10.data",
    "eRES_2H_SLAC_e49a6.data",
    "eRES_2H_SLAC_e49b.data",
    "eRES_2H_SLAC_e61.data",
    "eRES_2H_SLAC_e891.data",
    "eRES_2H_SLAC_e8920.data",
    "eRES_2H_SLAC_ne11.data"
  };
  // Skip a particular file if something seems odd with the data till
  // the issus is clarified.
  bool skip[nfiles] = 
  {
    false, // eRES_1H_JLab_e00_002.data
    false, // eRES_1H_JLab_e94_110.data
    false, // eRES_1H_JLab_ioanna.data
    false, // eRES_1H_SLAC_e133.data
    false, // eRES_1H_SLAC_e140.data
    false, // eRES_1H_SLAC_e140x.data
    false, // eRES_1H_SLAC_e49a10.data
    false, // eRES_1H_SLAC_e49a6.data
    false, // eRES_1H_SLAC_e49b.data
    false, // eRES_1H_SLAC_e61.data
    false, // eRES_1H_SLAC_e87.data
    false, // eRES_1H_SLAC_e891.data
    false, // eRES_1H_SLAC_e8920.data
    false, // eRES_1H_SLAC_ne11.data
    false, // eRES_1H_SLAC_onen1haf.data
    false, // eRES_2H_JLab_ioanna.data
    true,  // eRES_2H_SLAC_e133.data - has entries with negative W^2
    false, // eRES_2H_SLAC_e140x.data
    false, // eRES_2H_SLAC_e49a10.data
    false, // eRES_2H_SLAC_e49a6.data
    false, // eRES_2H_SLAC_e49b.data
    false, // eRES_2H_SLAC_e61.data
    false, // eRES_2H_SLAC_e891.data
    false, // eRES_2H_SLAC_e8920.data
    false  // eRES_2H_SLAC_ne11.data
  };
  // For some data files, the quoted xsec error includes both statistical and systematic uncertainties.
  // But for some data files, the quoted xsec error includes only statistical uncertainties and the
  // systematic uncertainty is quoted separately. This % systematic uncertainty is listed below and is
  // added to the statistical uncertainty.
  double pct_xsec_err_to_add[nfiles] = 
  {
    0.000, // eRES_1H_JLab_e00_002.data
    0.000, // eRES_1H_JLab_e94_110.data
    0.000, // eRES_1H_JLab_ioanna.data
    5.000, // eRES_1H_SLAC_e133.data
    0.000, // eRES_1H_SLAC_e140.data
    5.000, // eRES_1H_SLAC_e140x.data
    5.000, // eRES_1H_SLAC_e49a10.data
    5.000, // eRES_1H_SLAC_e49a6.data
    5.000, // eRES_1H_SLAC_e49b.data
    5.000, // eRES_1H_SLAC_e61.data
    5.000, // eRES_1H_SLAC_e87.data
    5.000, // eRES_1H_SLAC_e891.data
    5.000, // eRES_1H_SLAC_e8920.data
    0.000, // eRES_1H_SLAC_ne11.data
    5.000, // eRES_1H_SLAC_onen1haf.data
    0.000, // eRES_2H_JLab_ioanna.data
    5.000, // eRES_2H_SLAC_e133.data
    5.000, // eRES_2H_SLAC_e140x.data
    5.000, // eRES_2H_SLAC_e49a10.data
    5.000, // eRES_2H_SLAC_e49a6.data
    5.000, // eRES_2H_SLAC_e49b.data
    5.000, // eRES_2H_SLAC_e61.data
    5.000, // eRES_2H_SLAC_e891.data
    5.000, // eRES_2H_SLAC_e8920.data
    0.000  // eRES_2H_SLAC_ne11.data
  };
  double whitlow_norm[nfiles] = 
  {
    1.000, // eRES_1H_JLab_e00_002.data
    1.000, // eRES_1H_JLab_e94_110.data
    1.000, // eRES_1H_JLab_ioanna.data
    1.000, // eRES_1H_SLAC_e133.data
    1.000, // eRES_1H_SLAC_e140.data
    1.000, // eRES_1H_SLAC_e140x.data
    1.000, // eRES_1H_SLAC_e49a10.data
    1.012, // eRES_1H_SLAC_e49a6.data
    0.981, // eRES_1H_SLAC_e49b.data
    1.011, // eRES_1H_SLAC_e61.data
    0.982, // eRES_1H_SLAC_e87.data
    0.989, // eRES_1H_SLAC_e891.data
    0.953, // eRES_1H_SLAC_e8920.data
    1.000, // eRES_1H_SLAC_ne11.data
    1.000, // eRES_1H_SLAC_onen1haf.data
    1.000, // eRES_2H_JLab_ioanna.data
    1.000, // eRES_2H_SLAC_e133.data
    1.000, // eRES_2H_SLAC_e140x.data
    1.000, // eRES_2H_SLAC_e49a10.data
    1.001, // eRES_2H_SLAC_e49a6.data
    0.981, // eRES_2H_SLAC_e49b.data
    1.033, // eRES_2H_SLAC_e61.data
    0.985, // eRES_2H_SLAC_e891.data
    0.949, // eRES_2H_SLAC_e8920.data
    1.000  // eRES_2H_SLAC_ne11.data
  };
  char* expt_name[nfiles] = 
  {
    "jlab_e00_002",  // eRES_1H_JLab_e00_002.data
    "jlab_e94_110",  // eRES_1H_JLab_e94_110.data
    "jlab_ioanna",   // eRES_1H_JLab_ioanna.data
    "jlab_e133",     // eRES_1H_SLAC_e133.data
    "slac_e140",     // eRES_1H_SLAC_e140.data
    "slac_e140x",    // eRES_1H_SLAC_e140x.data
    "slac_e49a10",   // eRES_1H_SLAC_e49a10.data
    "slac_e49a6",    // eRES_1H_SLAC_e49a6.data
    "slac_e49b",     // eRES_1H_SLAC_e49b.data
    "slac_e61",      // eRES_1H_SLAC_e61.data
    "slac_e87",      // eRES_1H_SLAC_e87.data
    "slac_e891",     // eRES_1H_SLAC_e891.data
    "slac_e8920",    // eRES_1H_SLAC_e8920.data
    "slac_ne11",     // eRES_1H_SLAC_ne11.data
    "slac_onen1haf", // eRES_1H_SLAC_onen1haf.data
    "jlab_ioanna",   // eRES_2H_JLab_ioanna.data
    "slac_e133",     // eRES_2H_SLAC_e133.data
    "slac_e140x",    // eRES_2H_SLAC_e140x.data
    "slac_e49a10",   // eRES_2H_SLAC_e49a10.data
    "slac_e49a6",    // eRES_2H_SLAC_e49a6.data
    "slac_e49b",     // eRES_2H_SLAC_e49b.data
    "slac_e61",      // eRES_2H_SLAC_e61.data
    "slac_e891",     // eRES_2H_SLAC_e891.data
    "slac_e8920",    // eRES_2H_SLAC_e8920.data
    "slac_ne11"      // eRES_2H_SLAC_ne11.data
  };

  // Store all unique (E,theta) experimental configurations for each file.
  // This facilitates grouping the data later on as data points are typically
  // plotted for fixed E and theta.
  multimap<double,double> data_summary[nfiles];

  //
  // Read all data and save to ROOT tree
  //

  // define output tree
  int     Z        = 0;   // atomic mass number
  int     A        = 0;   // target mass number
  double  E        = 0.0; // incoming  electron energy, GeV
  double  Ep       = 0.0; // scattered electron energy, GeV
  double  theta    = 0.0; // scattering angle, degrees 
  double  Q2       = 0.0; // momentum transfer, GeV^2
  double  W2       = 0.0; // hadronic invariant mass, GeV^2
  double  v        = 0.0; // energy loss v, GeV
  double  epsilon  = 0.0; // 
  double  gamma    = 0.0; // 
  double  x        = 0.0; // bjorken x
  double  wnorm    = 0.0; // Whitlow norm
  double  xsec     = 0.0; // cross section, nb/sr/GeV
  double  xsec_err = 0.0; // cross section uncertainty, nb/sr/GeV
  char *  expt     = "";  // experiment

  TFile outfile("eRES.root", "RECREATE");
  TTree restree("resnt", "Resonance Electron Nucleus Scattering Archive");
  restree.Branch ("Z",        &Z,        "Z/I"        );
  restree.Branch ("A",        &A,        "A/I"        );
  restree.Branch ("E",        &E,        "E/D"        );
  restree.Branch ("Ep",       &Ep,       "Ep/D"       );
  restree.Branch ("theta",    &theta,    "theta/D"    );
  restree.Branch ("Q2",       &Q2,       "Q2/D"       );
  restree.Branch ("W2",       &W2,       "W2/D"       );
  restree.Branch ("v",        &v,        "v/D"        );
  restree.Branch ("epsilon",  &epsilon,  "epsilon/D"  );
  restree.Branch ("gamma",    &gamma,    "gamma/D"    );
  restree.Branch ("x",        &x,        "x/D"        );
  restree.Branch ("wnorm",    &wnorm,    "wnorm/D"    );
  restree.Branch ("xsec",     &xsec,     "xsec/D"     );
  restree.Branch ("xsec_err", &xsec_err, "xsec_err/D" );
  restree.Branch ("expt", (void*)expt,   "expt/C", 128);

  // temp vars
  double xsec_err_to_add = 0;
  double xsec_err_quoted = 0;

  // loop over files
  for(int i = 0; i < nfiles; i++) {
     // use this file?
     if(skip[i]) continue;
     // open the input data file
     ifstream infile(filename[i].c_str());
     if ( !infile.good() ) {
          cerr << "Can't open file: " << filename[i] << endl;
          return 1;
     }
     cout << "** ROOTify data from: " << filename[i] << endl;
     // get name, whitlow norm and Z,A which are common for all entries
     // of the current data file
     strcpy(expt, expt_name[i]);
     wnorm = whitlow_norm[i];
     Z = -99;
     A = -99;
     if(filename[i].find("_1H_") != string::npos) { Z = 1; A = 1; }
     if(filename[i].find("_2H_") != string::npos) { Z = 1; A = 2; }
     // loop over rows
     while(1) {
          // skip header lines staring with #
          if(infile.peek() == '#') {
            infile.ignore(2048,'\n');
#ifdef _show_debug_mesg_
	    cout << "Skipping header line..." << endl;
#endif
          } else {
   	    infile >> E >> Ep >> theta >> Q2 >> W2 >> v >> epsilon >> gamma >> x >> xsec >> xsec_err_quoted;
            if(infile.eof()) break;            
            // add systematic error by hand if not included in the number quoted in the file
            if(pct_xsec_err_to_add[i] > 0) {
              xsec_err_to_add = xsec * pct_xsec_err_to_add[i]/100.;
              xsec_err = TMath::Sqrt(xsec_err_quoted*xsec_err_quoted + xsec_err_to_add*xsec_err_to_add);
            } else {
              xsec_err = xsec_err_quoted;
            }
#ifdef _show_debug_mesg_
	    cout << "Sigma(eRES; " << expt << "; "
                 << "Z = "         << Z 
                 << ", A = "       << A 
                 << ", E = "       << E     << " GeV"
                 << ", Ep = "      << Ep    << " GeV"
                 << ", theta = "   << theta << " deg"
                 << ", Q2 = "      << Q2    << " GeV^2" 
                 << ", W2 = "      << W2    << " GeV^2" 
                 << ", v = "       << v     << " GeV" 
                 << ", epsilon = " << epsilon
                 << ", gamma = "   << gamma
                 << ", x = "       << x
                 << ", wnorm = "   << wnorm << ") = "
                 << xsec << " +/- " << xsec_err << " nb/sr/GeV"
                 << endl;
#endif
            // add current entry to tree
            restree.Fill();
            // add entry to summary 
            bool found = false;
            multimap<double,double>::const_iterator it = data_summary[i].begin();
            for( ; it != data_summary[i].end(); ++it) {       
               if( TMath::Abs(E - it->first) < 1E-3 && TMath::Abs(theta - it->second) < 5E-2) {
                  found = true;
                  break;
               }
            }
            if(!found) {
              pair<double,double> p(E,theta);
              data_summary[i].insert(p);
            }
          }
     }//!eof
  }//nfiles

  cout << "Wrote " << restree.GetEntries() << " entries in the RES tree"<< endl;  
  outfile.Write();

  // Write-out a summary of the experimental data included in the ROOT tree
  ofstream summary_file;
  summary_file.open("summary.txt");
//  summary_file << fixed;
  summary_file << "# target    experiment   E (GeV)   theta (deg)" << endl;;
  for(int i = 0; i < nfiles; i++) {
    cout << "Expt: " << expt_name[i] << endl;
     int tgtpdg=0;
     if(filename[i].find("_1H_") != string::npos) { tgtpdg = 1000010010; }
     if(filename[i].find("_2H_") != string::npos) { tgtpdg = 1000010020; }

    multimap<double,double>::const_iterator it = data_summary[i].begin();
    for( ; it != data_summary[i].end(); ++it) {       
       double E     = it->first;
       double theta = it->second;
       summary_file << tgtpdg << " "
                    << setw(20) << setfill(' ') << expt_name[i] 
                    << setw(10) << setfill(' ') << E 
                    << setw(10) << setfill(' ') << theta 
                    << endl;
    }//it   
  }//i
  summary_file.close();

  return 0;
}

