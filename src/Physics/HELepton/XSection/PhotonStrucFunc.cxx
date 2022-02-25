//____________________________________________________________________________
/*
 Copyright (c) 2003-2022, The GENIE Collaboration
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

  int nucs[2] = { kPdgProton, kPdgNeutron };
  int pdgs[6] = { kPdgNuE, kPdgAntiNuE, kPdgNuMu, kPdgAntiNuMu, kPdgNuTau, kPdgAntiNuTau };
  
  for (int k=0; k<2; k++) {
    for(int j=0; j<6; j++) {
      string SFname = string(gSystem->Getenv("GENIE")) + "/data/evgen/photon-sf/PhotonSF_hitnuc"+std::to_string(nucs[k])+"_hitlep"+std::to_string(pdgs[j])+".dat";
      if ( gSystem->AccessPathName( SFname.c_str()) ) {
        LOG("PhotonStrucFunc", pWARN) << "File doesnt exist. SF table must be compute with gmkglressf.";        
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