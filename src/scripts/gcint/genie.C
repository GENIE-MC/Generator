void change_prompt()
{
    TRint* rint = dynamic_cast<TRint*>(gApplication);
    if (!rint) return;

    rint->SetPrompt("genie [%d] ");
}

void genie()
{   
    TString script_dir = gSystem->Getenv("GENIE");
    script_dir += "/src/scripts/gcint/";

    TString curr_dir = gSystem->pwd();

    gSystem->cd(script_dir.Data());

    gROOT->ProcessLine(".x loadincs.C");    
    gROOT->ProcessLine(".x loadlibs.C");    

    gSystem->cd(curr_dir.Data());

    change_prompt();
}
