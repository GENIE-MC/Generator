{
  TString genie_topdir = gSystem->Getenv("GENIE");

  TString mp = gROOT->GetMacroPath();
  TString ip;

  TString ginc = genie_topdir + "/src";
  const char* p = ginc.Data();
  if (p) {
    mp += ":";
    mp += p;
    ip += " -I";
    ip += p;
  }

  mp += ":/usr/local/include/";
  ip += "  -I/usr/local/include/";

  gROOT->SetMacroPath(mp.Data());
  gSystem->SetIncludePath(ip);

  TString dip = ".include ";
  dip += genie_topdir.Data();
  dip += " /usr/local/include/";
  dip += " /usr/local/include/log4cpp";
  gROOT->ProcessLine(dip.Data());
}

