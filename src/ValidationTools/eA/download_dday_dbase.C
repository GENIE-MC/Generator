//
// Download e- QE Scattering Archive maintained by Donald Day
// and save data in a ROOT tree
//
// root[0] .x download_dday_dbase.C(target_pdgc)
//
// where target_pdgc is one of
// - 1000010020  (deuterium)
// - 1000010030  (tritium)
// - 1000020030  (3He)
// - 1000020040  (4He)
// - 1000060120  (12C)
// - 1000080160  (16O)
// - 1000130270  (27Al)
// - 1000200400  (40Ca)
// - 1000260560  (56Fe)
// - 1000791970  (197Au)
// - 1000822080) (208Pb)
//
// The script uses wget to download the data.
// If you are behind a proxy, remember to http_proxy environental variable.
//

string kDonaldDayBaseURL = "http://faculty.virginia.edu/qes-archive/";

int download_dday_dbase(int pdgc)
{
  string data_file_path;
  string data_file_name;

  switch(pdgc) {

   // deuterium
   case(1000010020) :
   {
     data_file_path = "2H";
     data_file_name = "2H.dat";
     break;
   }
   // tritium
   case(1000010030) :
   {
     data_file_path = "3H";
     data_file_name = "3H.dat";
     break;
   }
   // 3He
   case(1000020030) :
   {
     data_file_path = "3He";
     data_file_name = "3He.dat";
     break;
   }
   // 4He
   case(1000020040) :
   {
     data_file_path = "4He";
     data_file_name = "4He.dat";
     break;
   }
   // 12C
   case(1000060120) :
   {
     data_file_path = "C12";
     data_file_name = "C12.dat";
     break;
   }
   // 16O
   case(1000080160) :
   {
     data_file_path = "O16";
     data_file_name = "O16.dat";
     break;
   }
   // 27Al
   case(1000130270) :
   {
     data_file_path = "Al";
     data_file_name = "27Al.dat";
     break;
   }
   // 40Ca
   case(1000200400) :
   {
     data_file_path = "Ca40";
     data_file_name = "40Ca.dat";
     break;
   }
   // 56Fe
   case(1000260560) :
   {
     data_file_path = "Fe";
     data_file_name = "56Fe.dat";
     break;
   }
   // 197Au
   case(1000791970) :
   {
     data_file_path = "Au";
     data_file_name = "197Au.dat";
     break;
   }
   // 208Pb
   case(1000822080) :
   {
     data_file_path = "Pb";
     data_file_name = "208Pb.dat";
     break;
   }
   default :
   {
     cout << "Don't know the path for data on the following target: " << pdgc;
     exit(1);
   }
  }

  ostringstream cmd;
  cmd << "wget " << kDonaldDayBaseURL << "/" << data_file_path << "/" << data_file_name;;

  gSystem->Exec(cmd.str().c_str());

  //
  // convert to a ROOT tree
  //
  int    Z        = 0;
  int    A        = 0;
  double E        = 0.0; // GeV
  double theta    = 0.0; // degrees 
  double v        = 0.0; // energy loss v, GeV
  double xsec     = 0.0; // nb/sr/GeV
  double dxsec    = 0.0; // random;
  string citation = "";  // spires notation

  ostringstream outf;
  outf << "eQE_" << pdgc << ".root";

  TFile outfile(outf.str().c_str(), "RECREATE");
  TTree qetree("qetree", "QE Archive Ntuple");
  qetree.Branch ("Z",     &Z,     "Z/I"    );
  qetree.Branch ("A",     &A,     "A/I"    );
  qetree.Branch ("E",     &E,     "E/D"    );
  qetree.Branch ("theta", &theta, "theta/D");
  qetree.Branch ("v",     &v,     "v/D"    );
  qetree.Branch ("xsec",  &xsec,  "xsec/D" );
  qetree.Branch ("dxsec", &dxsec, "dxsec/D");

  // open the input data file
  std::ifstream infile(data_file_name.c_str());
  if ( !infile.good() ) {
    cerr<<"Can't open file: " << data_file << endl;
    return EXIT_FAILURE;
  }

  while( !infile.eof() ) {
	infile >> Z >> A >> E >> theta >> v >> xsec >> dxsec >> citation;
	qetree.Fill();
	cout << "cite = " << citation << endl;
  }

  cout << "In the new tree there are: " << qetree.GetEntries() << " entries"<< endl;  

  outfile.Write();

  return 0;
}

