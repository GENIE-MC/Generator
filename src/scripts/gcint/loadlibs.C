{
 TString libs0 = gSystem->GetDynamicPath();
 TString libs  = libs0 + ":/usr/lib:/usr/local/lib:/opt/lib:/opt/local/lib";
 gSystem->SetDynamicPath(libs.Data());

 // PYTHIA 6 lib
 gSystem->Load("libPythia6.so");

 // extra ROOT libs
 gSystem->Load("libPhysics");
 gSystem->Load("libEG");
 gSystem->Load("libEGPythia6");
 gSystem->Load("libGeom");
 gSystem->Load("libTree");

 // libxml2 and log4cpp libs
 gSystem->Load("libxml2.so");
 gSystem->Load("liblog4cpp.so");

 // std GENIE libs
 gSystem->Load("libGMessenger.so");
 gSystem->Load("libGRegistry.so");
 gSystem->Load("libGAlgorithm.so");
 gSystem->Load("libGNumerical.so");
 gSystem->Load("libGInteraction.so");
 gSystem->Load("libGHEP.so");
 gSystem->Load("libGUtils.so");
 gSystem->Load("libGBase.so");
 gSystem->Load("libGPDG.so");
 gSystem->Load("libGPDF.so");
 gSystem->Load("libGElFF.so");
 gSystem->Load("libGDecay.so");
 gSystem->Load("libGFragmentation.so");
 gSystem->Load("libGBaryonResonance.so");
 gSystem->Load("libGNuclear.so");
 gSystem->Load("libGLlewellynSmith.so");
 gSystem->Load("libGCrossSections.so");	
 gSystem->Load("libGCharm.so");
 gSystem->Load("libGElas.so");
 gSystem->Load("libGGiBUU.so");
 gSystem->Load("libGReinSeghal.so");
 gSystem->Load("libGPaschos.so");
 gSystem->Load("libGQPM.so");
 gSystem->Load("libGBodekYang.so");
 gSystem->Load("libGEVGCore.so");
 gSystem->Load("libGEVGModules.so");
 gSystem->Load("libGEVGDrivers.so");
 gSystem->Load("libGNtuple.so");
 gSystem->Load("libGQEL.so");
 gSystem->Load("libGRES.so");
 gSystem->Load("libGDIS.so");
 gSystem->Load("libGCoh.so");
 gSystem->Load("libGDfrc.so");
 gSystem->Load("libGMEC.so");
 gSystem->Load("libGNuE.so");
 gSystem->Load("libGNuGamma.so");
 gSystem->Load("libGVHE.so");
 gSystem->Load("libGHadronTransp.so");
 gSystem->Load("libGGeo.so");
 gSystem->Load("libGFluxDrivers.so");
}

