//
// ROOT-ify neutrino-nucleus cross-section data files 
//
// root[0] .L rootify.C++
// root[1] rootify()
//
// C.Andreopoulos 
//

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <iomanip>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

//#define _show_debug_mesg_

using namespace std;

const char * trim_spaces(const char * in);

int rootify(void)
{
  const int buffer_size = 100;

  // open input data summary file
  ifstream summary_file("./summary.txt");
  if (!summary_file.good() ) {
     cerr << "Can't open data summary file" << endl;
     return 1;
  }

  // open output ROOT data file
  TFile outfile("nuXSec.root", "RECREATE");

  // define output ROOT n-tuple
  char   dataset         [buffer_size];
  char   citation_spires [buffer_size];
  char   citation        [buffer_size];
  char   measurement     [buffer_size];
  char   comments        [buffer_size];
  double E;
  double Emin;
  double Emax;
  double xsec_datum; // absolute cross-section or cross-section ratio
  double xsec_datum_err_p;
  double xsec_datum_err_m;

  TTree outtree("nuxsnt", "World Neutrino-Nucleus Cross-Section Data");

  outtree.Branch ("dataset",           (void*)dataset,         "dataset/C",          buffer_size);
  outtree.Branch ("citation_spires",   (void*)citation_spires, "citation_spires/C",  buffer_size);
  outtree.Branch ("citation",          (void*)citation,        "citation/C",         buffer_size);
  outtree.Branch ("measurement",       (void*)measurement,     "measurement/C",      buffer_size);
  outtree.Branch ("comments",          (void*)comments,        "comments/C",         buffer_size);
  outtree.Branch ("E",                 &E,                     "E/D"                );
  outtree.Branch ("Emin",              &Emin,                  "Emin/D"             );
  outtree.Branch ("Emax",              &Emax,                  "Emax/D"             );
  outtree.Branch ("xsec_datum",        &xsec_datum,            "xsec_datum/D"       );
  outtree.Branch ("xsec_datum_err_p",  &xsec_datum_err_p,      "xsec_datum_err_p/D" );
  outtree.Branch ("xsec_datum_err_m",  &xsec_datum_err_m,      "xsec_datum_err_m/D" );
  
  // loop over summary file entries
  while(1) {
      // skip header lines staring with #
      if(summary_file.peek() == '#') {
         summary_file.ignore(1000, '\n');
         //cerr << "ignore comment line..." << endl;
       } else {
         // read a line
         char   temp[buffer_size];
         string expt;
         int    id;
         int    np;
         double syst_err;
  	 summary_file.getline(temp, buffer_size, '|'); expt = trim_spaces(temp);
   	 summary_file.getline(temp, buffer_size, '|'); id   = atoi(temp);
	 summary_file.getline(temp, buffer_size, '|'); strcpy(citation_spires, trim_spaces(temp));
	 summary_file.getline(temp, buffer_size, '|'); strcpy(citation,        trim_spaces(temp));
   	 summary_file.getline(temp, buffer_size, '|'); np       = atoi(temp);
   	 summary_file.getline(temp, buffer_size, '|'); syst_err = atof(temp);
    	 summary_file.getline(temp, buffer_size, '|'); strcpy(measurement, trim_spaces(temp));
    	 summary_file.getline(temp, buffer_size,'\n'); strcpy(comments,    trim_spaces(temp));
         if(summary_file.eof()) break;            

         // form dataset key
	 strcpy(dataset, Form("%s,%d", expt.c_str(), id));

         // work-out data file directory & form filename
         string directory = "";
         if      ( string(measurement).find("numu CC pi0 / numu CC QE")                 != string::npos ) directory = "ccpi_ccqe";
         else if ( string(measurement).find("numubar CC inclusive / numu CC inclusive") != string::npos ) directory = "ccr";
         else if ( string(measurement).find("numu CC mu-mu+ / numu CC inclusive")       != string::npos ) directory = "ccchm";
         else if ( string(measurement).find("numu CC mu-e+ / numu CC inclusive")        != string::npos ) directory = "ccchm";
         else if ( string(measurement).find("numu CC mu-l+ / numu CC inclusive")        != string::npos ) directory = "ccchm";
         else if ( string(measurement).find("numubar CC mu+mu- / numubar CC inclusive") != string::npos ) directory = "ccchm";
         else if ( string(measurement).find("numubar CC mu+e- / numubar CC inclusive")  != string::npos ) directory = "ccchm";
         else if ( string(measurement).find("numubar CC mu+l- / numubar CC inclusive")  != string::npos ) directory = "ccchm";
         else if ( string(measurement).find("numu CC charm / numu CC inclusive")        != string::npos ) directory = "ccchm";
         else if ( string(measurement).find("CC inclusive")                             != string::npos ) directory = "ccincl";
         else if ( string(measurement).find("CC QE")                                    != string::npos ) directory = "ccqe";
         else if ( string(measurement).find("CC 1pi+")                                  != string::npos ) directory = "ccpi";
         else if ( string(measurement).find("CC 1pi-")                                  != string::npos ) directory = "ccpi";
         else if ( string(measurement).find("CC 1pi0")                                  != string::npos ) directory = "ccpi";
         else if ( string(measurement).find("CC pi+pi-")                                != string::npos ) directory = "ccpi";
         else if ( string(measurement).find("CC pi+pi0")                                != string::npos ) directory = "ccpi";
         else if ( string(measurement).find("CC pi+pi+")                                != string::npos ) directory = "ccpi";
         else if ( string(measurement).find("CC coherent pi+")                          != string::npos ) directory = "cccoh";
         else if ( string(measurement).find("CC coherent pi-")                          != string::npos ) directory = "cccoh";
         else if ( string(measurement).find("NC coherent pi0")                          != string::npos ) directory = "nccoh";
         else {
            cerr << "Don't know where to find data file for dataset: " << dataset << endl;
            return 1;
         }

         string data_filename = Form("%s/%s-%d.data", directory.c_str(), expt.c_str(), id);

         // print-out summary info for current dataset
         cout << "\n* Dataset: " << dataset << " [" << citation_spires << " : " << citation << "]"
              << "\n Measurement: " << measurement 
              << "\n Comments: " << comments
              << "\n Reading " << np << " data-points from: " << data_filename 
              << " (and adding " << syst_err << "% systematic error by hand...)"
              << endl;

         // open the data file for the current dataset and read data points
         ifstream data_file(data_filename.c_str());
         if (!data_file.good() ) {
	   cerr << "Can't open data file: " << data_filename << endl;
           return 1;
         }

         while(1) {
           // skip header lines staring with #
          if(data_file.peek() == '#') {
            data_file.ignore(1000, '\n');
            //cerr << "ignore comment line..." << endl;
          } else {
            double xsec_datum_err_quoted_p;
            double xsec_datum_err_quoted_m;
  	    data_file >> E >> Emin >> Emax >> xsec_datum >> xsec_datum_err_quoted_p >> xsec_datum_err_quoted_m;
            if(data_file.eof()) break;            

            // add gross systematic error by hand if not included in the number quoted in the file
            if(syst_err > 0) { // this is the % err assigned to the current dataset
              double xsec_datum_err_to_add = xsec_datum * syst_err/100.;
              xsec_datum_err_p = TMath::Sqrt(xsec_datum_err_quoted_p*xsec_datum_err_quoted_p + xsec_datum_err_to_add*xsec_datum_err_to_add);
              xsec_datum_err_m = TMath::Sqrt(xsec_datum_err_quoted_m*xsec_datum_err_quoted_m + xsec_datum_err_to_add*xsec_datum_err_to_add);
            } else {
              xsec_datum_err_p = xsec_datum_err_quoted_p;
              xsec_datum_err_m = xsec_datum_err_quoted_m;
            }
#ifdef _show_debug_mesg_
            // print-out data point
            bool isRatio   = ( string(measurement).find("/")        != string::npos );
            bool isCohXSec = ( string(measurement).find("coherent") != string::npos );
            if(isRatio) { cout << "  - R";    }
            else        { cout << "  - xsec"; }
  	    cout << "(E = " << E 
                 << "; [" << Emin << ", " << Emax << "] GeV) = " 
                 << xsec_datum 
                 << " +" << xsec_datum_err_p 
                 << " -" << xsec_datum_err_m;
            if(!isRatio) {
               if(isCohXSec) { cout << " 1E-38*cm2/GeV/nucleus"; }
               else          { cout << " 1E-38*cm2/GeV/nucleon"; }
            } 
            cout << endl;
#endif

            // add data point in output tree
            outtree.Fill();

          }//data file comment line?
       }//while reading lines from data file

       data_file.close();

    }//summary file comment line?
  }//while reading lines from summary file

  cout << "\nWrote " << outtree.GetEntries() << " entries in the vA cross-section data tree"<< endl;  

  outtree.Write();
  outfile.Close();

  summary_file.close();

  cout << "Done"<< endl;  

  return 0;
}

const char * trim_spaces(const char * in)
{
// trim leading and trailing spaces
  string out = in;

  if( out.find_first_not_of(" \n") != 0)
       out.erase( 0, out.find_first_not_of(" \n")  );

  if( out.find_last_not_of(" \n") != out.length() )
       out.erase( out.find_last_not_of(" \n")+1, out.length() );

  return out.c_str();
}
