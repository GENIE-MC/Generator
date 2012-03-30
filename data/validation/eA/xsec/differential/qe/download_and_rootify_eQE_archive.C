//
// Download e- QE Scattering Archive maintained by Donal Day 
// and (optionally) save data in a ROOT tree.
//
// The script uses wget to download the data.
// If you are behind a proxy, remember to set the http_proxy environental variable.
//
// C.Andreopoulos 
//

#include <string>
#include <sstream>
#include <fstream>

//#define _show_debug_mesg_

using namespace std;

string kDonaldDayBaseURL = "http://faculty.virginia.edu/qes-archive/";

int download_and_rootify_eQE_archive(bool rootify=true)
{
  const int ntargets = 13;

  const char * remote_filename[ntargets] = 
  {
    "2H/2H.dat",     // deuterium
    "3H/3H.dat",     // tritium
    "3He/3He.dat",   // 3He
    "4He/4He.dat",   // 4He
    "C12/12C.dat",   // 12C
    "O16/O16.dat",   // 16O
    "Al/27Al.dat",   // 27Al
    "Ca40/40Ca.dat", // 40Ca
    "Ca48/Ca48.dat", // 48Ca
    "Fe/56Fe.dat",   // 56Fe
    "Au/197Au.dat",  // 197Au
    "Pb/208Pb.dat",  // 208Pb
    "U238/238U.dat"  // 238U
  };
  const char * target[ntargets] = 
  {
    "2H",   
    "3H",   
    "3He",  
    "4He",  
    "12C",  
    "16O",  
    "27Al", 
    "40Ca", 
    "48Ca", 
    "56Fe", 
    "197Au",
    "208Pb",
    "238U"
  };

  //
  // Download data from D.Day's dbase
  //
  for(int i=0; i < ntargets; i++) {

    ostringstream local_filename;
    local_filename << "eQE_" << target[i] << ".dat";

    // create header
    ofstream header;
    header.open("tmp.header");
    header << "# " << endl;;
    header << "# Quasielastic Electron Nucleus Scattering Archive" << endl;
    header << "# O. Benhar, D. Day and I. Sick, Rev. Mod. Phys. 80, 189-224, 2008" << endl;
    header << "# Source: http://faculty.virginia.edu/qes-archive/" << endl;
    header << "# Target: " << target[i] << endl;
    header << "# " << endl;
    header << "# Z  A   E      Theta     v      xsec        error    citation" << endl;
    header << "#        (GeV)  (deg)    (GeV)  (nb/sr/GeV)                   " << endl;
    header << "# " << endl;;
    header.close();

    // download data
    ostringstream cmd_get;
    cmd_get << "wget --output-document=tmp.dat " << kDonaldDayBaseURL << "/" << remote_filename[i];
    gSystem->Exec(cmd_get.str().c_str());

    // merge
    ostringstream cmd_merge;
    cmd_merge << "cat tmp.header tmp.dat > " << local_filename.str();
    gSystem->Exec(cmd_merge.str().c_str());

    // clean-up
    gSystem->Exec("rm tmp.header");
    gSystem->Exec("rm tmp.dat");
  }

  //
  // Download citations
  //
  gSystem->Exec(Form("wget --output-document=citations.bib %s/data-archive.bib",kDonaldDayBaseURL.c_str()));

  //
  // Read all data and save to ROOT tree
  //
  if(rootify) {
     int     Z        = 0;
     int     A        = 0;
     double  E        = 0.0; // GeV
     double  theta    = 0.0; // degrees 
     double  v        = 0.0; // energy loss v, GeV
     double  xsec     = 0.0; // cross section, nb/sr/GeV
     double  xsec_err = 0.0; // random;
     char *  citation = "";  // in SPIRES notation

     TFile outfile("eQE.root", "RECREATE");
     TTree qetree("qent", "Quasielastic Electron Nucleus Scattering Archive");
     qetree.Branch ("Z",        &Z,        "Z/I"       );
     qetree.Branch ("A",        &A,        "A/I"       );
     qetree.Branch ("E",        &E,        "E/D"       );
     qetree.Branch ("theta",    &theta,    "theta/D"   );
     qetree.Branch ("v",        &v,        "v/D"       );
     qetree.Branch ("xsec",     &xsec,     "xsec/D"    );
     qetree.Branch ("xsec_err", &xsec_err, "xsec_err/D");
     qetree.Branch ("citation", (void*)citation, "citation/C", 128);
      
     for(int i=0; i < ntargets; i++) {
       // open the input data file
       ostringstream local_filename;
       local_filename << "eQE_" << target[i] << ".dat";
       ifstream infile(local_filename.str().c_str());
       if ( !infile.good() ) {
          cerr << "Can't open file: " << local_filename.str() << endl;
          return EXIT_FAILURE;
       }
       cout << "** ROOTify data from: " << local_filename.str() << endl;
       // loop over rows
       while(1) {
          // skip header lines staring with #
          if(infile.peek() == '#') {
            infile.ignore(2048,'\n');
#ifdef _show_debug_mesg_
	    cout << "Skipping header line..." << endl;
#endif
          } else {
   	    infile >> Z >> A >> E >> theta >> v >> xsec >> xsec_err >> citation;
            if(infile.eof()) break;            
            qetree.Fill();
#ifdef _show_debug_mesg_
	    cout << "Sigma(eQE; Z = " << Z << ", A = " << A 
                 << ", E = " << E << " GeV, theta = " << theta << " deg, v = " << v << " GeV) = " 
                 << xsec << " nb/sr/GeV [" << citation << "]" 
                 << endl;
#endif
          }
       }//eof
     }//i
     cout << "Wrote " << qetree.GetEntries() << " entries in the QE tree"<< endl;  
     outfile.Write();
  }//rootify?

  return 0;
}

