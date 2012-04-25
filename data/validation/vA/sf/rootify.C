//
// ROOT-ify structure function data 
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
	TFile outfile("structFunc.root", "RECREATE");

	// define output ROOT n-tuple
	char   dataset         [buffer_size];
	char   citation_spires [buffer_size];
	char   citation        [buffer_size];
	char   process         [buffer_size];
	char   strfunc         [buffer_size];
	char   comments        [buffer_size];
	int    target_pdg;
	int    nu_pdg;
	int    id;
	double x;
	double q2;
	double f;
	double f_err_p;
	double f_err_m;

	TTree outtree("structFunc", "Structure Function Data");

	outtree.Branch ("dataset",         (void*)dataset,         "dataset/C",          buffer_size);
	outtree.Branch ("citation_spires", (void*)citation_spires, "citation_spires/C",  buffer_size);
	outtree.Branch ("citation",        (void*)citation,        "citation/C",         buffer_size);
	outtree.Branch ("process",         (void*)process,         "process/C",          buffer_size);
	outtree.Branch ("comments",        (void*)comments,        "comments/C",         buffer_size);
	outtree.Branch ("strfunc",         (void*)strfunc,         "strfunc/C",          buffer_size);
	outtree.Branch ("nu_pdg",          	&nu_pdg,               "nu_pdg/I"    );
	outtree.Branch ("target_pdg",       &target_pdg,           "target_pdg/I");
	outtree.Branch ("x",               	&x,                    "x/D"      );
	outtree.Branch ("q2",              	&q2,                   "q2/D"     );
	outtree.Branch ("f", 		        &f,                    "f/D"      );
	outtree.Branch ("f_err_p",         	&f_err_p,              "f_err_p/D");
	outtree.Branch ("f_err_m",         	&f_err_m,              "f_err_p/D");

	// loop over summary file entries
	while(1) {
		// skip header lines staring with #
		if(summary_file.peek() == '#') {
			summary_file.ignore(1000, '\n');
			cerr << "ignore comment line..." << endl;
		} else {
			// read a line
			char   temp[buffer_size];
			string expt;
			int    np;//number of points
			double syst_err;
			summary_file.getline(temp, buffer_size, '|'); expt = trim_spaces(temp);
			summary_file.getline(temp, buffer_size, '|'); id   = atoi(temp);
			summary_file.getline(temp, buffer_size, '|'); strcpy(citation_spires, trim_spaces(temp));
			summary_file.getline(temp, buffer_size, '|'); strcpy(citation,        trim_spaces(temp));
			summary_file.getline(temp, buffer_size, '|'); np           = atoi(temp);
			summary_file.getline(temp, buffer_size, '|'); syst_err     = atof(temp);
			summary_file.getline(temp, buffer_size, '|'); nu_pdg       = atoi(temp);
			summary_file.getline(temp, buffer_size, '|'); target_pdg   = atoi(temp);
			summary_file.getline(temp, buffer_size, '|'); strcpy(process,  trim_spaces(temp));
			summary_file.getline(temp, buffer_size, '|'); strcpy(strfunc,  trim_spaces(temp));
			summary_file.getline(temp, buffer_size,'\n'); strcpy(comments, trim_spaces(temp));
			if(summary_file.eof()) break;            

			// form dataset key

			strcpy(dataset, Form("%s,%d", expt.c_str(), id));
			
			// work-out data file directory & form filename
			string directory = "";
			string bob = strfunc;//do not know better way of geting a string from a char........
			if (bob=="F2") directory = "f2";
			else if (bob=="xF3") directory = "xf3";
			else {
				cerr << "Don't know where to find data file for dataset: " << dataset << endl;
				return 1;
			}

			string data_filename = Form("%s/%s_%s.data", directory.c_str(),expt.c_str(), strfunc);
			cout<<data_filename<<endl;
			// print-out summary info for current dataset
			cout << "\n* Dataset: " << dataset << " [" << citation_spires << " : " << citation << "]"
				<< " - Neutrino: " << nu_pdg << ", target: " << target_pdg 
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
					double f_err_quoted_p;
					double f_err_quoted_m;
					data_file >> x  >> q2 >> f >> f_err_quoted_p;
					f_err_quoted_m = f_err_quoted_p;
					if(data_file.eof()) break;            

					// add gross systematic error by hand if not included in the number quoted in the file
					if(syst_err > 0) { // this is the % err assigned to the current dataset
						double f_err_to_add = f * syst_err/100.;
						f_err_p = TMath::Sqrt(f_err_quoted_p*f_err_quoted_p + f_err_to_add*f_err_to_add);
						f_err_m = TMath::Sqrt(f_err_quoted_m*f_err_quoted_m + f_err_to_add*f_err_to_add);
					} else {
						f_err_p = f_err_quoted_p;
						f_err_m = f_err_quoted_m;
					}
#ifdef _show_debug_mesg_
					// print-out data point
					cout << "   x = " << x << ", 
						<< f << "+ " << f_err_p << "- " << f__err_m
						<< endl;
#endif

					// add data point in output tree
					outtree.Fill();

				}//data file comment line?
			}//while reading lines from data file

			data_file.close();

		}//summary file comment line?
	}//while reading lines from summary file

	cout << "\nWrote " << outtree.GetEntries() << " entries in the sf data tree"<< endl;  

	outtree.Write();
	outfile.Close();

	summary_file.close();

	cout << "Done"<< endl;  

	return 0;
}
//===============================
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
