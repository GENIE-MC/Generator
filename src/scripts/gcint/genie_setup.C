void change_prompt()
{
    TRint* rint = dynamic_cast<TRint*>(gApplication);
    if (!rint) return;

    rint->SetPrompt("genie [%d] ");
}

void genie_setup()
{
    TString script_dir = gSystem->Getenv("GENIE");
    script_dir += "/src/scripts/gcint/";
    
    int genie_incs_status = -999;
    int genie_libs_status = -999;

    TString loadincs = "#include \""+script_dir+"/loadincs.C\"";
    TString execute_loadincs = TString::Format("*(int*)%p = loadincs();",&genie_incs_status);
    
    TString loadlibs = "#include \""+script_dir+"/loadlibs.C\"";
    TString execute_loadlibs = TString::Format("*(int*)%p = loadlibs();",&genie_libs_status);
    
    gROOT->ProcessLine(loadincs.Data());
    gROOT->ProcessLine(execute_loadincs.Data());
    gROOT->ProcessLine(loadlibs.Data());
    gROOT->ProcessLine(execute_loadlibs.Data());
    
    if (genie_incs_status || genie_libs_status) {
      std::cerr<<"Failed to load GENIE libraries."<<std::endl;
      exit(-1);
    }

    change_prompt();
    
    gROOT->ProcessLine("genie::RunOpt::Instance()->ReadFromCommandLine( 0, 0 );" ) ;
    gROOT->ProcessLine("genie::RunOpt::Instance()->BuildTune();");

}
