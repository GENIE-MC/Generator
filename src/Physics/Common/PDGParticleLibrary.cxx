//____________________________________________________________________________
/*
 Copyright (c) 2003-2020, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory 

 Marco Roda <mroda \at liverpool.ac.uk>                                                                    
 University of Liverpool                                                                                   
 
*/
//____________________________________________________________________________

#include "Physics/Common/PDGParticleLibrary.h"

#include <string>

using namespace genie;
using std::string ;

PDGParticleLibrary::PDGParticleLibrary() :
  ParticleLibrary( "genie::PDGParticleLibrary", "Default" ),
  fDatabase( nullptr ) {

  TDatabasePDG::Instance() ;

  // The TDatabasePDG has a quasi-singletone implementation which is not great
  // After dedicated experimentats did with ROOT 6.22/06 we figured out that 
  // 1) if the Database singletone access is called first, 
  //    any subsequent new allocation of the TDatbasePDG will create another 
  //    different databaser
  // 2) if a normal constructor is called before the first call of 
  //    TDatabasePDG::Instance(), the singletone will call the first
  //    the first constructed database 
  //    while any new constructor will create a separate database. 
  // Since we want the this class to be independed from the rest of the world
  // as we might want to change it freely 
  // in the constructor we call the singletone so that later we will
  // be able to instantiate separate objects wihtout interfering with 
  // the rest 
}

//___________________________________________________________________________

PDGParticleLibrary::PDGParticleLibrary( const string & config )  :
  ParticleLibrary( "genie::PDGParticleLibrary", config ),
  fDatabase( nullptr ) {

  TPDGDatabase::Instance() ;
  // see default constructor for the explanation of the previous call
} 
//___________________________________________________________________________
PDGParticleLibrary::~PDGParticleLibrary() {
  ;
}
//___________________________________________________________________________
TParticlePDG * PDGParticleLibrary::Find(int pdgc) const {
  
  return fDatabase -> GetParticle(pdgc); 
}
//___________________________________________________________________________
void PDGParticleLibrary::LoadConfig( void ) {

  bool good_configuration = true ;

  string file_name = PDGParticleLibrary::PDGTableFile() ;

  if ( gSystem->AccessPathName( file_name.c_str() ) ) {
    // failure to load as file not existent
    LOG("PDGParticleLibrary", pERROR) << "Invalid file: " << file_name ;
    good_configuration = false ;
  }

  if ( ! good_configuration ) {

    LOG("PDGParticleLibrary", pFATAL) << "Invalid configuration: Exiting" ;
  }

  // actual allocation 
  fDatabase.reset( new TDatabasePDG() ) ;
  fDataset -> ReadPDGTable( file_name.c_str() );

}
//___________________________________________________________________________
std::string PDGParticleLibrary::PDGTableFile() const {

  // loading PDG data from a number of possibilities
  // in order of priority
  // 1) GENIE_PDG_TABLE env variable set
  // 2) in the dedicated $GENIE data directory
  //    using the tune defined file
  // 3) using the bare ROOT one


  // loading PDG data from $GENIE/config/
  const char* altpdgtable = gSystem->Getenv("GENIE_PDG_TABLE");
  if ( altpdgtable ) {
    
    LOG("PDGParticleLibrary", pINFO) << "Trying to load PDG data from $GENIE_PDG_TABLE: "
				     << altpdgtable;

    if ( ! (gSystem->AccessPathName(altpdgtable) ) ) {
      return altpdgtable ;
    }
  }

  if ( gSystem->Getenv("GENIE") ) {
    string base_dir = string( gSystem->Getenv("GENIE") );
    string path = base_dir + string("/data/evgen/catalogues/pdg/");

    string file_name ;
    GetParamDef( "PDG-TableName", file_name, "genie_pdg_table.txt" ) ;

    string file = path + file_name ;
    LOG("PDGParticleLibrary", pINFO) << "Trying to load PDG data from: " << file;

    if ( ! (gSystem->AccessPathName( file.c_str()) ) ) {
      return file ;
    }
  }

  // no PDG data in $GENIE/config/ - Try $ROOTSYS/etc/

  if(gSystem->Getenv("ROOTSYS")) {
    string base_dir  = string( gSystem->Getenv("ROOTSYS") );
    string path = base_dir  + string("/etc/pdg_table.txt");
    
    LOG("PDGParticleLibrary", pINFO) << "Trying to load PDG data from: " << path;

    if ( !(gSystem->AccessPathName(path.c_str())) ) {

      return path;
    }
  }

  LOG("PDGParticleLibrary", pERROR) << " PDG loading failure " ; 
  return "";


}
