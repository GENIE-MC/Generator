Long_t load_single_lib(const char * libname)
{
  return gSystem->Load(libname);
}

Long_t load_libs_from_command(const char * list_libs_command, bool verbose = false)
{ // Takes a linker-style list of libs (eg "-lmygreatlib -lmylessgoodlib" )
  // and loads each using gSystem::Load()
  // Returns 0 on success, 1 on failure
  Long_t status = 0;
  FILE * f = gSystem->OpenPipe(list_libs_command,"r");

  TPRegexp re("-l([\\d\\w]*)");
  while (true) {
    TString line;
    if (!line.Gets(f)) {break;}
    TObjArray * tokens = line.Tokenize(" ");
    for (int i = 0 ; i < tokens->GetEntries() ; i++) {
      TObjString * token_os = static_cast<TObjString*>(tokens->At(i));
      if (!token_os) {continue;}
      TObjArray * matches = re.MatchS(token_os->GetString());
      if (matches->GetEntries()!=2) { continue; }
      TObjString * libname_os = static_cast<TObjString*>(matches->At(1));
      if (!libname_os) {continue;}
      TString full_libname = "lib"+libname_os->GetString();
      Long_t this_status = load_single_lib(full_libname.Data());
      if (verbose) std::cout<<full_libname<<":\t"<<this_status<<std::endl;
      status = (status || this_status);
      delete matches;
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

  // PYTHIA 6 lib
  gSystem->Load("libPythia6");

  // extra ROOT libs
  gSystem->Load("libPhysics");
  gSystem->Load("libEG");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libGeom");
  gSystem->Load("libTree");

  // libxml2 and log4cpp libs
  gSystem->Load("libxml2");
  gSystem->Load("liblog4cpp");

  // GSL
  gSystem->Load("libgslcblas");
  gSystem->Load("libgsl");
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

