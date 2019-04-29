//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab 

         Changes required to implement the GENIE Boosted Dark Matter module
         were installed by Josh Berger (Univ. of Wisconsin)
*/
//____________________________________________________________________________

#include <iostream>
#include <string>

#include <TSystem.h>

#include "Framework/Messenger/Messenger.h"
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
TParticlePDG * PDGLibrary::Find(int pdgc)
{
// save some typing in the most frequently typed TDatabasePDG method

  return fDatabasePDG->GetParticle(pdgc);
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
    string path = base_dir + 
      string("/data/evgen/catalogues/pdg/genie_pdg_table.txt");

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
// EDIT: need a way to clear and then reload the PDG database
void PDGLibrary::ReloadDBase(void)
{
  if(fDatabasePDG) {
    delete fDatabasePDG;
  }

  if( ! LoadDBase() ) LOG("PDG", pERROR) << "Could not load PDG data";
}
