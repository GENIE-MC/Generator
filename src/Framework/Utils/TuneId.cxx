//____________________________________________________________________________
/*
 Copyright (c) 2003-2018, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Marco Roda <Marco.Roda \at liverpool.ac.uk>
         University of Liverpool

         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk
         University of Liverpool & STFC Rutherford Appleton Lab
*/
//____________________________________________________________________________

//#include <sstream>

#include "Framework/Utils/TuneId.h"

#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/XmlParserUtils.h"
#include "Framework/Utils/SystemUtils.h"

#include "Framework/Messenger/Messenger.h"

//using std::ostringstream;

using namespace genie;

//____________________________________________________________________________
namespace genie
{
  ostream & operator << (ostream & stream, const TuneId & id)
  {
    id.Print(stream);
    return stream;
  }
  //..........................................................................
  bool operator == (const TuneId & id1, const TuneId & id2)
  {
    return id1.Compare(id2);
  }
  //..........................................................................
  bool operator != (const TuneId & id1, const TuneId & id2)
  {
    return !id1.Compare(id2);
  }
}
//____________________________________________________________________________
TuneId::TuneId(const string & id_str) :
  fName(id_str)
{
  this -> Decode(id_str);

  if ( this -> CheckDirectory() ) {
    LOG("TuneId", pINFO) << Name() <<" Tune configured " ;
  }
  else {
    LOG("TuneId", pFATAL) << "No valid tune directory associated with " << Name() ;
    exit(0) ;
  }
}
//____________________________________________________________________________
TuneId::TuneId(const TuneId & id) 
{
  this->Copy(id);

  if ( ! CheckDirectory() ) {
    LOG("TuneId", pWARN) << "No valid subdirectory associated with " << Name() ;
  }
}
//____________________________________________________________________________
TuneId::~TuneId()
{

}
//____________________________________________________________________________
string TuneId::CGC(void) const {

  string cgc = fPrefix ;
  cgc += fYear ;
  cgc += "_" ;
  cgc += ModelId() ;

  return cgc ;
}
//____________________________________________________________________________
string TuneId::Tail(void) const {

  string tail = fTunedParamSetId ;
  tail += "_" + fFitDataSetId    ;
  return tail ;
}
//____________________________________________________________________________
string TuneId::CGCDirectory(void) const {

  string dir = fBaseDirectory ;
  dir += "/" + CGC() ;

  return dir ;

}
//____________________________________________________________________________
string TuneId::TuneDirectory   (void) const {

  string dir = CGCDirectory() ;
  if ( ! OnlyConfiguration() ) dir += "/" + Name() ;

  return dir ;
}
//____________________________________________________________________________
void TuneId::Decode(string id_str)
{

  std::vector<string> parts = utils::str::Split( id_str, "_" ) ;

  this -> fPrefix        = parts[0].substr( parts[0].size()-2 ) ;
  this -> fYear          = parts[0].substr( 0, parts[0].size()-2 ) ;
  this -> fMajorModelId  = parts[1].substr( 0, 2 ) ;
  this -> fMinorModelId  = parts[1].substr( 2 ) ;
  this->fTunedParamSetId = parts[2] ;
  this ->fFitDataSetId   = parts[3] ;

}
//____________________________________________________________________________
void TuneId::Copy(const TuneId & id)
{
  this->fName            = id.Name();
  this->fPrefix          = id.Prefix();
  this->fYear            = id.Year();
  this->fMajorModelId    = id.MajorModelId();
  this->fMinorModelId    = id.MinorModelId();
  this->fTunedParamSetId = id.TunedParamSetId();
  this->fFitDataSetId    = id.FitDataSetId();
}
//____________________________________________________________________________
bool TuneId::Compare(const TuneId & id) const
{
  return (this->Name() == id.Name());
}
//____________________________________________________________________________
void TuneId::Print(ostream & stream) const
{
  stream << " GENIE tune: " << std::endl;
  stream << " - Prefix ............... : " << this->Prefix()          << std::endl;
  stream << " - Year ................. : " << this->Year()            << std::endl;
  stream << " - Major model ID ....... : " << this->MajorModelId()    << std::endl;
  stream << " - Minor model ID ....... : " << this->MinorModelId()    << std::endl;
  stream << " - Tuned param set ID ... : " << this->TunedParamSetId() << std::endl;
  stream << " - Fit dataset ID ....... : " << this->FitDataSetId()    << std::endl;
}
//____________________________________________________________________________
bool TuneId::CheckDirectory() {

  std::string pathlist = genie::utils::xml::GetXMLPathList(false) ;
  std::vector<std::string> paths = genie::utils::str::Split(pathlist,":;,");

  fBaseDirectory = "" ;

  for ( size_t i=0; i< paths.size(); ++i ) {
     const char* tmppath = paths[i].c_str();
     std::string onepath = gSystem->ExpandPathName(tmppath);
     string test = onepath + "/" + CGC() ;
     LOG("TuneId", pDEBUG) << " Testing  " << test << " directory" ;
     if ( utils::system::DirectoryExists( test.c_str() ) ) {
       fBaseDirectory = onepath ;
       break ;
     }
  }

  if ( fBaseDirectory.size() == 0 ) {
    LOG("TuneId", pWARN) << " No " << CGC() << " subdirectory found in pathlist";
    return false ;
  }

  if ( ! OnlyConfiguration() ) {
    if ( ! utils::system::DirectoryExists( TuneDirectory().c_str() ) ) {
      LOG("TuneId", pWARN) << "No " << Name() << " subdirectory found in " << CGC() ;
      return false ;
    }
  }

  return true ;
}

