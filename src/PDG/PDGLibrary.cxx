//____________________________________________________________________________
/*!

\class    genie::PDGLibrary

\brief    Singleton class to load & serve a TDatabasePDG.

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 06, 2004

*/
//____________________________________________________________________________

#include <iostream>
#include <string>

#include <TSystem.h>

#include "Messenger/Messenger.h"
#include "PDG/PDGLibrary.h"

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

  if(gSystem->Getenv("GENIE")) {
    string base_dir = string( gSystem->Getenv("GENIE") );
    string path = base_dir + string("/data/pdg/genie_pdg_table.txt");

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

