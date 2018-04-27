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

#include "TPRegexp.h"
#include "TObjArray.h"
#include "TObjString.h"

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
TuneId::TuneId(const TuneId & id) 
{
  this->Copy(id);

  if ( ! CheckDirectory() ) {
    LOG("TuneId", pWARN) << "No valid subdirectory associated with " << Name() ;
  }
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
void TuneId::Build(const string & name ) {

  if ( name.size() > 0 ) fName = name ;

  this -> Decode( fName );

  if ( this -> CheckDirectory() ) {
    LOG("TuneId", pINFO) << Name() <<" Tune configured " ;
  }
  else {
    LOG("TuneId", pFATAL) << "No valid tune directory associated with " << Name() ;
    exit(0) ;
  }
}
//____________________________________________________________________________
void TuneId::Decode(string id_str)
{
  static TPRegexp pattern("([A-Za-z]+)(\\d{2})_(\\d{2})([a-z])_(\\d{2})_(\\d{3})");
  TString tstr(id_str.c_str());
  TObjArray * matches = pattern.MatchS(tstr);
  if ( matches -> GetEntries() != 7) {
    LOG("TuneId", pFATAL) << "Bad tune pattern "<<id_str<<" - form is eg G18_01a_00_000";
    exit(-1);
  }
  
  this -> fPrefix          = ((TObjString*)matches->At(1))->String().Data();
  this -> fYear            = ((TObjString*)matches->At(2))->String().Data();
  this -> fMajorModelId    = ((TObjString*)matches->At(3))->String().Data();
  this -> fMinorModelId    = ((TObjString*)matches->At(4))->String().Data();
  this -> fTunedParamSetId = ((TObjString*)matches->At(5))->String().Data();
  this -> fFitDataSetId    = ((TObjString*)matches->At(6))->String().Data();
  
  delete matches;
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
  stream << (IsCustom() ? "Custom" : "Standard") << " GENIE tune: " << this -> Name()  << std::endl;
  stream << " - Prefix ............... : " << this->Prefix()          << std::endl;
  stream << " - Year ................. : " << this->Year()            << std::endl;
  stream << " - Major model ID ....... : " << this->MajorModelId()    << std::endl;
  stream << " - Minor model ID ....... : " << this->MinorModelId()    << std::endl;
  stream << " - Tuned param set ID ... : " << this->TunedParamSetId() << std::endl;
  stream << " - Fit dataset ID ....... : " << this->FitDataSetId()    << std::endl;
  stream << " - Tune directory ....... : " << this->fBaseDirectory    << std::endl;
  if ( IsCustom() )
    stream << " - Custom directory ..... : " << this -> fCustomSource   << std::endl;
}
//____________________________________________________________________________
bool TuneId::CheckDirectory() {

  std::string pathlist = utils::xml::GetXMLPathList(false) ;
  std::vector<std::string> paths = utils::str::Split(pathlist,":;,");

  string top_path = gSystem->ExpandPathName( paths[0].c_str() ) ;
  string def_path = gSystem->ExpandPathName( utils::xml::GetXMLDefaultPath().c_str() ) ;

  if ( top_path != def_path ) {
    fCustomSource = top_path ;
  }

  fBaseDirectory = "" ;
  LOG("TuneId",pDEBUG) << "Base dir validation " ;
  
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

  LOG("TuneId",pDEBUG) << fBaseDirectory ;

  return true ;
}

