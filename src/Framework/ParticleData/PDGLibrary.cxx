//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)

         Changes required to implement the GENIE Dark Neutrino module
         were installed by Iker de Icaza (Univ. of Sussex)
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
void PDGLibrary::AddNHL(double mass)
{
// Add NHL to PDG database

  TParticlePDG * nhl = fDatabasePDG->GetParticle(kPdgNHL);
  if (!nhl) {
    // Name Title Mass Stable Width Charge Class PDG
    fDatabasePDG->AddParticle("NHL","NHL",mass,true,0.,0,"NHL",kPdgNHL);
  }
  else {
    assert(nhl->Mass() == mass);
  }
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
