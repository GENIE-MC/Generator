{
  TString libs0 = gSystem->GetDynamicPath();
  TString libs  = libs0 + ":/usr/lib:/usr/local/lib:/opt/lib:/opt/local/lib";
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

  //
  // GENIE libs
  //
  
  // Read them from genie-config
  TString command = TString::Format("genie-config --libs");
  FILE * f = gSystem->OpenPipe(command.Data(),"r");
  
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
      //~ cerr<<full_libname<<endl;
      gSystem->Load(full_libname.Data());
      delete matches;
    }
    delete tokens;
  }
  gSystem->ClosePipe(f);
}

