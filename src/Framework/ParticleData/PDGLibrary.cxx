//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)

         Changes required to implement the GENIE Dark Neutrino module
         were installed by Iker de Icaza (Univ. of Sussex)
	 
	 Changes required to implement the GENIE BeamHNL module
	 were installed by John Plows (Univ. of Oxford)
*/
//____________________________________________________________________________

#include <iostream>
#include <string>

#include <TSystem.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"


using std::string;

using namespace genie;

//____________________________________________________________________________
PDGLibrary * PDGLibrary::fInstance = 0;
//____________________________________________________________________________
PDGLibrary::PDGLibrary()
{
  if( ! LoadDBase() ) LOG("PDG", pERROR) << "Could not load PDG data";

#ifdef __GENIE_DARK_NEUTRINO_ENABLED__
  LOG("PDG", pINFO) << "Loading Dark sector Info";
  if ( ! AddDarkSector() ) { 
    LOG("PDG", pFATAL) << "Could not load Dark Neutrino data";
    exit(78);
  }
#endif // __GENIE_DARK_NEUTRINO_ENABLED__

#ifdef __GENIE_HEAVY_NEUTRAL_LEPTON_ENABLED__
  LOG("PDG", pINFO) << "Loading Heavy Neutral Lepton data";
  if( ! AddHNL() ){
    LOG("PDG", pFATAL) << "Could not load Heavy Neutral Lepton data";
    exit(78);
  }
#endif // #ifdef __GENIE_HEAVY_NEUTRAL_LEPTON_ENABLED__
  
  fInstance =  0;
}
//____________________________________________________________________________
PDGLibrary::~PDGLibrary()
{
  fInstance = 0;
}
//____________________________________________________________________________
PDGLibrary * PDGLibrary::Instance()
{
  if(fInstance == 0) {
    LOG("PDG", pINFO) << "PDGLibrary late initialization";

    static PDGLibrary::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new PDGLibrary;
  }
  return fInstance;
}
//____________________________________________________________________________
TDatabasePDG * PDGLibrary::DBase(void)
{
  return fDatabasePDG;
}
//____________________________________________________________________________
TParticlePDG * PDGLibrary::Find(int pdgc, bool must_exist )
{

  TParticlePDG * ret = fDatabasePDG->GetParticle(pdgc);
  if(ret) return ret;

  if ( must_exist ) {
    LOG("PDG", pERROR) << "Requested missing particle with PDG: " << pdgc ;
  }

  return ret ;
}

//____________________________________________________________________________
bool PDGLibrary::LoadDBase(void)
{
  fDatabasePDG = TDatabasePDG::Instance();

  // loading PDG data from $GENIE/config/
  const char* altpdgtable = gSystem->Getenv("GENIE_PDG_TABLE");
  if ( altpdgtable ) {
    if ( ! (gSystem->AccessPathName(altpdgtable) ) ) {
        LOG("PDG", pINFO) << "Load PDG data from $GENIE_PDG_TABLE: "
                          << altpdgtable;
        fDatabasePDG->ReadPDGTable( altpdgtable );
        return true;
    }
  }

  if ( gSystem->Getenv("GENIE") ) {
    string base_dir = string( gSystem->Getenv("GENIE") );
    base_dir += string("/data/evgen/catalogues/pdg/") ; 

    string file_name = "genie_pdg_table.txt" ; 
    const Registry * reg = AlgConfigPool::Instance()->CommonList("Param", "PDG");
    if( reg ) {
      file_name = reg -> GetString("PDG-TableName") ;
      LOG("PDG", pINFO) << "Found file name specification: " << file_name ;

    }
    
    string path = base_dir + file_name ;

    if ( ! (gSystem->AccessPathName(path.c_str()) ) ) {
        LOG("PDG", pINFO) << "Load PDG data from: " << path;
        fDatabasePDG->ReadPDGTable( path.c_str() );
        return true;
    }
  }

  // no PDG data in $GENIE/config/ - Try $ROOTSYS/etc/

  if(gSystem->Getenv("ROOTSYS")) {
    string base_dir  = string( gSystem->Getenv("ROOTSYS") );
    string path = base_dir  + string("/etc/pdg_table.txt");

    if ( !(gSystem->AccessPathName(path.c_str())) ) {
        LOG("PDG", pINFO) << "Load PDG data from: " << path;
        fDatabasePDG->ReadPDGTable( path.c_str() );
        return true;
     }
  }

  LOG("PDG", pERROR) << " *** The PDG extensions will not be loaded!! ***";
  return false;
};
//____________________________________________________________________________
void PDGLibrary::AddDarkMatter(double mass, double med_ratio)
{
// Add dark matter particle to PDG database

  double med_mass = mass*med_ratio;
  TParticlePDG * dm_particle = fDatabasePDG->GetParticle(kPdgDarkMatter);
  TParticlePDG * med_particle = fDatabasePDG->GetParticle(kPdgMediator);
  if (!dm_particle) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("chi_dm","chi_dm",mass,true,0.,0,"DarkMatter",kPdgDarkMatter);
  }
  else {
    assert(dm_particle->Mass() == mass);
  }
  if (!med_particle) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("Z_prime","Z_prime",med_mass,true,0.,0,"DarkMatter",kPdgMediator);
  }
  else {
    assert(med_particle->Mass() == med_mass);
  }
}
//____________________________________________________________________________
bool PDGLibrary::AddHNL()
{
  // Add HNL to PDG database
  const Registry * reg = AlgConfigPool::Instance()->CommonList("HNL", "ParameterSpace");
  if (!reg) {
    LOG("PDG", pERROR) << "Cannot find HNL ParameterSpace param_set";
    return false;
  }
  TParticlePDG * hnl = fDatabasePDG->GetParticle(kPdgHNL);
  if (!hnl) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("HNL","HNL",reg->GetDouble("HNL-Mass"),true,0.,0,"HNL",kPdgHNL);
    fDatabasePDG->AddParticle("HNLBar","HNLBar",reg->GetDouble("HNL-Mass"),true,0.,0,"HNL",-1*kPdgHNL);
  }
  return true;
}
//____________________________________________________________________________
bool PDGLibrary::AddDarkSector()
{
  // Add dark neutrino particles to PDG database

  const Registry * reg = AlgConfigPool::Instance()->CommonList("Dark", "Masses");
  if(!reg) {
    LOG("PDG", pERROR) << "The Dark Sector masses not available.";
    return false;
  }
  TParticlePDG * dnu_particle = fDatabasePDG->GetParticle(kPdgDarkNeutrino);
  TParticlePDG * anti_dnu_particle = fDatabasePDG->GetParticle(kPdgAntiDarkNeutrino);
  TParticlePDG * med_particle = fDatabasePDG->GetParticle(kPdgDNuMediator);
  if (!dnu_particle) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("nu_D","#nu_{D}",reg->GetDouble("Dark-NeutrinoMass"),
                              true,0.,0,"DarkNeutrino",kPdgDarkNeutrino);
  }
  if (!anti_dnu_particle) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("nu_D_bar","#bar{#nu}_{D}",reg->GetDouble("Dark-NeutrinoMass"),
                              true,0.,0,"DarkNeutrino",kPdgAntiDarkNeutrino);
  }
  if (!med_particle) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("Z_D","Z_{D}",reg->GetDouble("Dark-MediatorMass"),
                              true,0.,0,"DarkNeutrino",kPdgDNuMediator);
  }
  return true;
}
//____________________________________________________________________________
// EDIT: need a way to clear and then reload the PDG database
void PDGLibrary::ReloadDBase(void)
{
  if(fDatabasePDG) {
    delete fDatabasePDG;
  }

  if( ! LoadDBase() ) LOG("PDG", pERROR) << "Could not load PDG data";
}
//____________________________________________________________________________
