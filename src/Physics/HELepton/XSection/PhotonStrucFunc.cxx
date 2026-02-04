//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Alfonso Garcia <aagarciasoto \at km3net.de>
 IFIC & Harvard University
*/
//____________________________________________________________________________

#include "Physics/HELepton/XSection/PhotonStrucFunc.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"

#include "TSystem.h"

using namespace genie;

//_________________________________________________________________________
PhotonStrucFunc * PhotonStrucFunc::fgInstance = 0;
//_________________________________________________________________________
PhotonStrucFunc::PhotonStrucFunc()
{

  string basedir = "";
  if ( gSystem->Getenv("PHOTON_SF_DATA_PATH")==NULL ) basedir = string(gSystem->Getenv("GENIE")) + "/data/evgen/photon-sf";
  else                                                basedir = string(gSystem->Getenv("PHOTON_SF_DATA_PATH"));
  LOG("PhotonStrucFunc", pWARN) << "Base diretory: " << basedir;

  int nucs[2] = { kPdgProton, kPdgNeutron };
  int pdgs[6] = { kPdgNuE, kPdgAntiNuE, kPdgNuMu, kPdgAntiNuMu, kPdgNuTau, kPdgAntiNuTau };
  
  for (int k=0; k<2; k++) {
    for(int j=0; j<6; j++) {
      string SFname = basedir + "/PhotonSF_hitnuc"+std::to_string(nucs[k])+"_hitlep"+std::to_string(pdgs[j])+".dat";
      if ( gSystem->AccessPathName( SFname.c_str(), kReadPermission ) ) {
        LOG("PhotonStrucFunc", pFATAL) << "File doesnt exist or you dont have read permission.";        
        LOG("PhotonStrucFunc", pFATAL) << "Remember!!!";
        LOG("PhotonStrucFunc", pFATAL) << "Path to base directory is defined with the enviroment variable PHOTON_SF_DATA_PATH.";
        LOG("PhotonStrucFunc", pFATAL) << "If not defined, default location is $GENIE/data/evgen/photon-sf";
        LOG("PhotonStrucFunc", pFATAL) << "Photon SF tables must be computed with gmkglressf.";        
        assert(0);
      }
      fSFTables[nucs[k]].Table[pdgs[j]] = new genie::Spline();
      fSFTables[nucs[k]].Table[pdgs[j]]->LoadFromAsciiFile(SFname);
    }        
  }

  fgInstance = 0;

}
//_________________________________________________________________________
PhotonStrucFunc::~PhotonStrucFunc()
{

}
//_________________________________________________________________________
PhotonStrucFunc * PhotonStrucFunc::Instance()
{
  if(fgInstance == 0) {
    LOG("PhotonStrucFunc", pINFO) << "Late initialization";
    static PhotonStrucFunc::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fgInstance = new PhotonStrucFunc();
  }  
  return fgInstance;
}