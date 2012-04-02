//
// ROOT-ify multiplicity ratio data files 
//
// C.Andreopoulos 
//

#include <fstream>
#include <string>
#include <iomanip>

//#define _show_debug_mesg_

using namespace std;

int rootify(void)
{
  const int nfiles = 12;

  string filename[nfiles] = 
  {
      "hermes/NPB780_2007_fig2-RhA_vsQ2-K+_He.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-K+_Kr.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-K+_Ne.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-K+_Xe.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-p_He.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-p_Kr.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-p_Ne.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-p_Xe.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-pi+_He.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-pi+_Kr.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-pi+_Ne.dat",
      "hermes/NPB780_2007_fig2-RhA_vsQ2-pi+_Xe.dat"
  };


  // define output tree
  int     Z        = 0;   // atomic mass number
  int     A        = 0;   // target mass number
  int     hadron   = 0;   // hadon PDG code
  double  Q2       = 0.0; // momentum transfer, GeV^2
  double  R        = 0.0; // multiplicity ratio
  double  dRp      = 0.0; // uncertainty on multiplicity ratio (+)
  double  dRm      = 0.0; // uncertainty on multiplicity ratio (-)
  char *  citation = "";  // in SPIRES notation

  TFile outfile("Rh.root", "RECREATE");
  TTree hmrnt("hmr", "Hadron Multiplicity Ratios");
  hmrnt.Branch ("Z",         &Z,             "Z/I"        );
  hmrnt.Branch ("A",         &A,             "A/I"        );
  hmrnt.Branch ("hadron",    &hadron,        "hadron/I"   );
  hmrnt.Branch ("Q2",        &Q2,            "Q2/D"       );
  hmrnt.Branch ("Rh",        &R,             "R/D"        );
  hmrnt.Branch ("dRp",       &dRp,           "dRp/D"      );
  hmrnt.Branch ("dRm",       &dRm,           "dRm/D"      );
  hmrnt.Branch ("citation", (void*)citation, "citation/C", 128);

  // loop over files
  for(int i = 0; i < nfiles; i++) {
     ifstream infile(filename[i].c_str());
     if ( !infile.good() ) {
          cerr << "Can't open file: " << filename[i] << endl;
          return 1;
     }
     cout << "** ROOTify data from: " << filename[i] << endl;
     while(1) {
          // skip header lines staring with #
          if(infile.peek() == '#') {
            infile.ignore(2048,'\n');
#ifdef _show_debug_mesg_
	    cout << "Skipping header line..." << endl;
#endif
          } else {
            double Rp = 0, Rm = 0;
   	    infile >> Q2 >> R >> Rp >> Rm >> Z >> A >> hadron >> citation;
            if(infile.eof()) break;            
            dRm = TMath::Abs(R-Rm);
            dRp = TMath::Abs(R-Rp);
#ifdef _show_debug_mesg_
	    cout << "R(h = " << hadron << "; " 
                 << "Z = "         << Z 
                 << ", A = "       << A 
                 << ", Q2 = "      << Q2    << " GeV^2) = " 
                 << R << " +" << dRp << " -dRp" << " [" << citation << "]" 
                 << endl;
#endif
            // add current entry to tree
            hmrnt.Fill();
          }
     }//!eof
  }//nfiles

  cout << "Wrote " << hmrnt.GetEntries() << " entries in the hadron multiplicity ratios tree"<< endl;  
  outfile.Write();

  return 0;
}

