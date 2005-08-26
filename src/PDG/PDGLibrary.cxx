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
  bool pdg_data_loaded = false;

  fDatabasePDG = TDatabasePDG::Instance();

  //-- get base GENIE & ROOT base directory from the environment
  string genie_base_dir = string( gSystem->Getenv("GENIE") );
  string root_base_dir  = string( gSystem->Getenv("ROOTSYS") );

  //-- build the full pathnames for possible pdg data file locations
  string path_1 = genie_base_dir + string("/data/pdg/genie_pdg_table.txt");
  string path_2 = root_base_dir  + string("/etc/pdg_table.txt");

  if ( ! (gSystem->AccessPathName( path_1.c_str() ) ) ) {

     // loading PDG data from $GENIE/config/
     LOG("PDG", pINFO)
              << "\n ***** Loading PDG data and extensions from " << path_1;
     fDatabasePDG->ReadPDGTable( path_1.c_str() );
     pdg_data_loaded = true;

  }  else {
      // no PDG data in $GENIE/config/ - Try $ROOTSYS/etc/
     if ( ! (gSystem->AccessPathName( path_2.c_str() ) ) ) {

        // no PDG data in $ROOTSYS/etc either !!
        LOG("PDG", pWARN)
           << "\n ***** Loading standard PDG data from " << path_2
                      << " - The PDG extensions will not be loaded!! [<-bad]";

        fDatabasePDG->ReadPDGTable( path_2.c_str() );
        pdg_data_loaded = true;
     }
  }
  return pdg_data_loaded;
};
//____________________________________________________________________________

