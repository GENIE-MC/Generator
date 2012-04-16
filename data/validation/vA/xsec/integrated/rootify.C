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
  char   process         [buffer_size];
  char   target          [buffer_size];
  char   comments        [buffer_size];
  int    nu_pdg;
  double E;
  double Emin;
  double Emax;
  double xsec;
  double xsec_err_p;
  double xsec_err_m;

  TTree outtree("nuxsnt", "World Neutrino-Nucleus Cross-Section Data");

  outtree.Branch ("dataset",         (void*)dataset,         "dataset/C",          buffer_size);
  outtree.Branch ("citation_spires", (void*)citation_spires, "citation_spires/C",  buffer_size);
  outtree.Branch ("citation",        (void*)citation,        "citation/C",         buffer_size);
  outtree.Branch ("process",         (void*)process,         "process/C",          buffer_size);
  outtree.Branch ("target",          (void*)target,          "target/C",           buffer_size);
  outtree.Branch ("comments",        (void*)comments,        "comments/C",         buffer_size);
  outtree.Branch ("nu_pdg",          &nu_pdg,                "nu_pdg/I"    );
  outtree.Branch ("E",               &E,                     "E/D"         );
  outtree.Branch ("Emin",            &Emin,                  "Emin/D"      );
  outtree.Branch ("Emax",            &Emax,                  "Emax/D"      );
  outtree.Branch ("xsec",            &xsec,                  "xsec/D"      );
  outtree.Branch ("xsec_err_p",      &xsec_err_p,            "xsec_err_p/D");
  outtree.Branch ("xsec_err_m",      &xsec_err_m,            "xsec_err_m/D");
  
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
   	 summary_file.getline(temp, buffer_size, '|'); nu_pdg   = atoi(temp);
   	 summary_file.getline(temp, buffer_size, '|'); strcpy(target,   trim_spaces(temp));
    	 summary_file.getline(temp, buffer_size, '|'); strcpy(process,  trim_spaces(temp));
    	 summary_file.getline(temp, buffer_size,'\n'); strcpy(comments, trim_spaces(temp));
         if(summary_file.eof()) break;            

         // form dataset key
	 strcpy(dataset, Form("%s,%d", expt.c_str(), id));

         // work-out data file directory & form filename
         bool isCC  = (string(process).find("CC")        != string::npos);
         bool isNC  = (string(process).find("NC")        != string::npos);
         bool isQE  = (string(process).find("QE")        != string::npos);
         bool isCoh = (string(process).find("coherent")  != string::npos);
         bool isPi  = (string(process).find("pi")        != string::npos);
         bool isInc = (string(process).find("inclusive") != string::npos);
         string directory = "";
         if      (isCC && isCoh) directory = "cccoh";
         else if (isNC && isCoh) directory = "nccoh";
         else if (isCC && isPi ) directory = "ccpi";
         else if (isNC && isPi ) directory = "ncpi";
         else if (isCC && isQE ) directory = "ccqe";
         else if (isNC && isQE ) directory = "ncel";
         else if (isCC && isInc) directory = "ccincl";
         else if (isNC && isInc) directory = "ncincl";
         else {
           cerr << "Don't know where to find data file for dataset: " << dataset << endl;
           return 1;
         }
         string data_filename = Form("%s/%s-%d.data", directory.c_str(), expt.c_str(), id);

         // print-out summary info for current dataset
         cout << "\n* Dataset: " << dataset << " [" << citation_spires << " : " << citation << "]"
              << " - Neutrino: " << nu_pdg << ", target: " << target 
              << ", process: " << process << ", comments: " << comments
              << "\n  Reading " << np << " data-points from: " << data_filename 
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
            double xsec_err_quoted_p;
            double xsec_err_quoted_m;
  	    data_file >> E >> Emin >> Emax >> xsec >> xsec_err_quoted_p >> xsec_err_quoted_m;
            if(data_file.eof()) break;            

            // add gross systematic error by hand if not included in the number quoted in the file
            if(syst_err > 0) { // this is the % err assigned to the current dataset
              double xsec_err_to_add = xsec * syst_err/100.;
              xsec_err_p = TMath::Sqrt(xsec_err_quoted_p*xsec_err_quoted_p + xsec_err_to_add*xsec_err_to_add);
              xsec_err_m = TMath::Sqrt(xsec_err_quoted_m*xsec_err_quoted_m + xsec_err_to_add*xsec_err_to_add);
            } else {
              xsec_err_p = xsec_err_quoted_p;
              xsec_err_m = xsec_err_quoted_m;
            }
#ifdef _show_debug_mesg_
            // print-out data point
  	    cout << "   - xsec (E = " << E << ", E bin = [" << Emin << ", " << Emax << "] GeV) = " 
                 << xsec << "+ " << xsec_err_p << "- " << xsec_err_m
                 << ( (isCoh) ? " 1E-38*cm2/GeV/nucleus" : " 1E-38*cm2/GeV/nucleon" ) 
                 << endl;
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
