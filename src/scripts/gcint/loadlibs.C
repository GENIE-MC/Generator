
Long_t load_single_lib(const char * libname)
{
  // returns
  //  0 on success
  //  1 if lib already loaded
  // -1 error or non-existent
  // -2 if version mismatch
  return gSystem->Load(libname);
}
Long_t add_single_lpath(const char * lpath, bool verbose = false)
{
   if (verbose) std::cout << "AddDynamicPath(" << lpath << ")" << std::endl;
   gSystem->AddDynamicPath(lpath);
   // AddDynamicPath doesn't return a status, so just hope for the best
   return 0;
}

Long_t has_feature(const char * feature_name, bool verbose = true)
{
  // take a feature name, e.g. "pythia6" and
  // Returns 1 if enabled, 0 if not
  Long_t status = 0;
  TString command = TString("genie-config --has-") + feature_name;
  FILE * f = gSystem->OpenPipe(command,"r");
  TString line;
  line.Gets(f);
  if (verbose) {
    std::cout << "command: " << command.Data() << " result: " << line.Data() << std::endl;
  }
  if (line == "yes") status = 1;
  if (verbose) {
    std::cout << "return " << status << std::endl;
  }
  gSystem->ClosePipe(f);
  return status;
}

Long_t load_libs_from_command(const char * list_libs_command, bool verbose = false)
{ // Takes a linker-style list of libs (eg "-lmygreatlib -lmylessgoodlib" )
  // and loads each using gSystem::Load()
  // Returns 0 on success, 1 on failure
  Long_t status = 0;
  FILE * f = gSystem->OpenPipe(list_libs_command,"r");

  TPRegexp relib("^-l([\\d\\w]*)");
  TPRegexp reLpath("^-L([\\d\\w/]*)");

  while (true) {
    TString line;
    if (!line.Gets(f)) {break;}
    TObjArray * tokens = line.Tokenize(" ");

    for (int i = 0 ; i < tokens->GetEntries() ; i++) {
      TObjString * token_os = static_cast<TObjString*>(tokens->At(i));
      if (verbose) {
        std::cout << std::endl << "token_os: '" << token_os->GetString()
                  << "'" << std::endl;
      }
      if (!token_os) {continue;}

      TObjArray * matcheslib = relib.MatchS(token_os->GetString());
      TObjArray * matchesLpath = reLpath.MatchS(token_os->GetString());

      if (matcheslib->GetEntries()==2) {
        TObjString * libname_os = static_cast<TObjString*>(matcheslib->At(1));
        if (verbose && libname_os ) {
          std::cout << "libname_os: '" << libname_os->GetString() << "'"
                    << std::endl;
        } else if (verbose) std::cout << "libname_os: was NULL" << std::endl;

        if (!libname_os) {continue;}
        TString full_libname = "lib"+libname_os->GetString();
        Long_t this_status = load_single_lib(full_libname.Data());
        if (verbose) {
          std::cout << "load_single_lib " << full_libname
                    << " status:\t" << this_status << std::endl;
        }
        status = (status || this_status);
      } else {

        if (verbose) {
          std::cout << "matchesLpath->GetEntries() "
                    << matchesLpath->GetEntries() << std::endl;
        }
        if (matchesLpath->GetEntries()==2) {
          TObjString * libpath_os = static_cast<TObjString*>(matchesLpath->At(1));
          if (verbose && libpath_os ) {
            std::cout << "libpath_os: '" << libpath_os->GetString() << "'"
                      << std::endl;
          } else if (verbose) std::cout << "libpath_os: was NULL" << std::endl;
          if (!libpath_os) {continue;}
          TString full_lpath = libpath_os->GetString();
          Long_t this_status = add_single_lpath(full_lpath.Data());
          if (verbose) {
            std::cout << "add_single_lpath " << full_lpath
                      <<" status:\t" << this_status << std::endl;
          }
          status = (status || this_status);
        }
        delete matcheslib;
        delete matchesLpath;
      }
    }
    delete tokens;
  }
  gSystem->ClosePipe(f);
  return status;
}

int loadlibs()
{
  TString libs0 = gSystem->GetDynamicPath();
  TString libs  = libs0 + ":/usr/local/lib64:/usr/lib64:/lib64:"+
    "/usr/local/lib/:/usr/lib/:/lib:/opt/lib:/opt/local/lib";
  gSystem->SetDynamicPath(libs.Data());

  bool has_pythia6 = has_feature("pythia6");

  // PYTHIA 6 lib
  if (has_pythia6) gSystem->Load("libPythia6");

  // extra ROOT libs
  gSystem->Load("libPhysics");
  gSystem->Load("libEG");
  if (has_pythia6) gSystem->Load("libEGPythia6");
  gSystem->Load("libGeom");
  gSystem->Load("libTree");
  gSystem->Load("libMathMore");
  gSystem->Load("libMinuit");
  gSystem->Load("libGenVector");
  
  // libxml2 and log4cpp libs
  gSystem->Load("libxml2");
  gSystem->Load("liblog4cpp");
  
  // GSL
  gSystem->Load("libgslcblas");
  gSystem->Load("libgsl");
  gSystem->Load("libm");
  //~ load_libs_from_command("gsl-config --libs"); // not guaranteed to be in the right order

  //
  // GENIE libs
  //

  Long_t status = load_libs_from_command("genie-config --libs");

  if (status) {
    std::cout<<"Failed to load GENIE libraries"<<std::endl;
    exit(1);
  }
  return 0;
}

