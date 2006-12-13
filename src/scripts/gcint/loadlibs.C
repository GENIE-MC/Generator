{
 // PYTHIA 6 lib
 gSystem->Load("libPythia6.so");

 // extra ROOT libs
 gSystem->Load("libPhysics");
 gSystem->Load("libEG");
 gSystem->Load("libEGPythia6");
 gSystem->Load("libGeom");

 // libxml2 and log4cpp libs
 gSystem->Load("libxml2.so");
 gSystem->Load("liblog4cpp.so");

 // std GENIE libs
 gSystem->Load("libGMessenger.so");
 gSystem->Load("libGRegistry.so");
 gSystem->Load("libGAlgorithm.so");
 gSystem->Load("libGNumerical.so");
 gSystem->Load("libGInteraction.so");
 gSystem->Load("libGUtils.so");
 gSystem->Load("libGBase.so");
 gSystem->Load("libGPDG.so");
 gSystem->Load("libGPDF.so");
 gSystem->Load("libGDecay.so");
 gSystem->Load("libGFragmentation.so");
 gSystem->Load("libGBaryonResonance.so");
 gSystem->Load("libGNuclear.so");
 gSystem->Load("libGLlewellynSmith.so");
 gSystem->Load("libGCrossSections.so");
 gSystem->Load("libGInvMuDecay.so");
 gSystem->Load("libGCharm.so");
 gSystem->Load("libGElas.so");
 gSystem->Load("libGReinSeghal.so");
 gSystem->Load("libGRadiative.so");
 gSystem->Load("libGPaschos.so");
 gSystem->Load("libGVHE.so");
 gSystem->Load("libGPartonModel.so");
 gSystem->Load("libGBodekYang.so");
 gSystem->Load("libGHEP.so");
 gSystem->Load("libGEVGCore.so");
 gSystem->Load("libGEVGModules.so");
 gSystem->Load("libGEVGDrivers.so");
 gSystem->Load("libGNtuple.so");

 gSystem->Load("libGGeo.so");
 gSystem->Load("libGFluxDrivers.so");
}

