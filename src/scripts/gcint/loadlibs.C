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

  // framework and utilities
  gSystem->Load("libGMessenger");
  gSystem->Load("libGRegistry");
  gSystem->Load("libGAlgorithm");
  gSystem->Load("libGInteraction");
  gSystem->Load("libGHEP");
  gSystem->Load("libGBase");
  gSystem->Load("libGNumerical");
  gSystem->Load("libGUtils");
  gSystem->Load("libGPDG");
  gSystem->Load("libGBaryonResonance");
  gSystem->Load("libGEVGCore");
  gSystem->Load("libGEVGDrivers");
  gSystem->Load("libGNtuple");
  gSystem->Load("libGGeo");
  gSystem->Load("libGFluxDrivers");

  // physics
  gSystem->Load("libGPDF");
  gSystem->Load("libGElFF");
  gSystem->Load("libGDecay");
  gSystem->Load("libGFragmentation");
  gSystem->Load("libGNuclear");
  gSystem->Load("libGLlewellynSmith");
  gSystem->Load("libGCrossSections");	
  gSystem->Load("libGCharm");
  gSystem->Load("libGElas");
  gSystem->Load("libGGiBUU");
  gSystem->Load("libGReinSehgal");
  gSystem->Load("libGQPM");
  gSystem->Load("libGBodekYang");
  gSystem->Load("libGEVGModules");
  gSystem->Load("libGQEL");
  gSystem->Load("libGRES");
  gSystem->Load("libGDIS");
  gSystem->Load("libGCoh");
  gSystem->Load("libGDfrc");
  gSystem->Load("libGMEC");
  gSystem->Load("libGNuE");
  gSystem->Load("libGNuGamma");
//gSystem->Load("libGVLE");
//gSystem->Load("libGVHE");
  gSystem->Load("libGHadronTransp");
}

