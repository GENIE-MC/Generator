//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
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
TuneId::TuneId(const string & id_str, bool failOnInvalid)
  : fName(genie::utils::str::TrimSpaces(id_str)) // remove any lead/trailing
  , fIsConfigured(false)
  , fIsValidated(false)
{
  Build(fName);
  if ( failOnInvalid && ! fIsValidated ) {
    //    status & 0377 is returned to parent on exit() call e.g. [0:255]
    // SYSEXITS(3)       FreeBSD Library Functions Manual
    // According to style(9),  it is not a good practice to call exit(3) with
    // arbitrary values to indicate a failure condition when ending a program.
    // Instead, the pre-defined exit codes from sysexits should be used, so the
    // caller of the process can get a rough estimation about the failure class
    // without looking up the source code.
    // EX_USAGE (64)        The command was used incorrectly, e.g., with the
    //                      wrong number of arguments, a bad flag, a bad syntax
    //                      in a parameter, or whatever.
    // EX_UNAVAILABLE (69)  A service is unavailable.  This can occur if a sup­
    //                      port program or file does not exist.  This can also
    //                      be used as a catchall message when something you
    //                      wanted to do doesn't work, but you don't know why.

    //   use 64 when failed Decode name (i.e. ! fIsConfigured )
    //   use 69 when failed to find directory (i.e. ! fIsValidated )
    if ( fIsConfigured ) exit(69);
    else                 exit(64);
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
string TuneId::CMC(void) const {

  string cmc = fPrefix ;
  cmc += fYear ;
  cmc += "_" ;
  cmc += ModelId() ;

  return cmc ;
}
//____________________________________________________________________________
string TuneId::Tail(void) const {

  string tail = fTunedParamSetId ;
  tail += "_" + fFitDataSetId    ;
  return tail ;
}
//____________________________________________________________________________
string TuneId::CMCDirectory(void) const {

  string dir = fBaseDirectory ;
  dir += "/" + CMC() ;

  return dir ;

}
//____________________________________________________________________________
string TuneId::TuneDirectory   (void) const {

  string dir = CMCDirectory() ;
  if ( ! OnlyConfiguration() ) dir += "/" + Name() ;

  return dir ;
}
//____________________________________________________________________________
void TuneId::Build(const string & name ) {
  LOG("TuneId",pDEBUG)<<"Building tune "<<name;
  if ( name.size() > 0 ) fName = name ;

  this -> Decode( fName );
  if ( ! fIsConfigured ) return; // no point going on

  if ( this -> CheckDirectory() ) {
    LOG("TuneId", pINFO) << Name() <<" Tune configured " ;
    fIsValidated = true;
  } else {
    LOG("TuneId", pFATAL) << "No valid tune directory associated with " << Name() ;
    fIsValidated = false;
  }
}
//____________________________________________________________________________
void TuneId::Decode(string id_str)
{
  static TPRegexp pattern("^([A-Za-z]+)(\\d{2})_(\\d{2})([a-z])_([a-z0-9]{2})_([a-z0-9]{3})$");
  // TPRegexp pattern("([A-Za-z]+)(\\d{2})_(\\d{2})([a-z])_(\\d{2})_(\\d{3})");

  TString tstr(id_str.c_str());
  TObjArray * matches = pattern.MatchS(tstr);
  if ( matches -> GetEntries() != 7) {
    LOG("TuneId", pFATAL) << "Bad tune pattern "<<id_str<<" - form is eg G18_01a_00_000";
    fIsConfigured = false;
    return;
  } else {
    fIsConfigured = true;
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

  this->fIsConfigured    = id.IsConfigured();
  this->fIsValidated     = id.IsValidated();
}
//____________________________________________________________________________
bool TuneId::Compare(const TuneId & id) const
{
  return (this->Name() == id.Name());
}
//____________________________________________________________________________
void TuneId::Print(ostream & stream) const
{
  std::string             status = "Standard";
  if ( IsCustom()       ) status = "Custom";
  if ( ! IsValidated()  ) status = "BadDirectory";
  if ( ! IsConfigured() ) status = "BadConfig";
  stream << status << " GENIE tune: " << this -> Name()  << std::endl;
  stream << " - Prefix ............... : " << this->Prefix()          << std::endl;
  stream << " - Year ................. : " << this->Year()            << std::endl;
  stream << " - Major model ID ....... : " << this->MajorModelId()    << std::endl;
  stream << " - Minor model ID ....... : " << this->MinorModelId()    << std::endl;
  stream << " - Tuned param set ID ... : " << this->TunedParamSetId() << std::endl;
  stream << " - Fit dataset ID ....... : " << this->FitDataSetId()    << std::endl;
  stream << " - Tune directory ....... : " << this->TuneDirectory()   << std::endl;
  stream << " - Base directory ....... : " << this->fBaseDirectory    << std::endl;
  if ( IsCustom() )
    stream << " - Custom directory ..... : " << this -> fCustomSource   << std::endl;

  stream << std::flush;
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
    string test = onepath + "/" + CMC() ;
    LOG("TuneId", pDEBUG) << " Testing  " << test << " directory" ;
    if ( utils::system::DirectoryExists( test.c_str() ) ) {
      fBaseDirectory = onepath ;
      break ;
    }
  }

  if ( fBaseDirectory.size() == 0 ) {
    LOG("TuneId", pWARN) << " No " << CMC() << " subdirectory found in pathlist";
    return false ;
  }

  if ( ! OnlyConfiguration() ) {
    if ( ! utils::system::DirectoryExists( TuneDirectory().c_str() ) ) {
      LOG("TuneId", pWARN) << "No " << Name() << " subdirectory found in " << CMC() ;
      return false ;
    }
  }

  LOG("TuneId",pDEBUG) << fBaseDirectory ;

  return true ;
}
