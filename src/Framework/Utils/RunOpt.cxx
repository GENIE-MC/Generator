//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <iostream>
#include <cstdlib>

#include <TMath.h>
#include <TBits.h>

#include "Framework/GHEP/GHepFlags.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Messenger/Messenger.h"

using std::cout;
using std::endl;

namespace genie {

  static const string gDefaultTune = "G18_02a_00_000";

//____________________________________________________________________________
ostream & operator << (ostream & stream, const RunOpt & opt)
{
  opt.Print(stream);
  return stream;
}
//____________________________________________________________________________
RunOpt * RunOpt::fInstance = 0;
//____________________________________________________________________________
RunOpt::RunOpt() : fTune(0)
{
  fInstance = 0;

  this->Init();
}
//____________________________________________________________________________
RunOpt::~RunOpt()
{
  if ( fTune )                    delete fTune ;
  if ( fUnphysEventMask )         delete fUnphysEventMask ;
  fInstance = 0;
}
//____________________________________________________________________________
RunOpt * RunOpt::Instance()
{
  if(fInstance == 0) {
    static RunOpt::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new RunOpt;
  }
  return fInstance;
}
//____________________________________________________________________________
void RunOpt::Init(void)
{
  fTune = 0 ;
  fEnableBareXSecPreCalc = true;
  fCacheFile = "";
  fMesgThresholds = "";
  fUnphysEventMask = new TBits(GHepFlags::NFlags());
//fUnphysEventMask->ResetAllBits(true);
  for(unsigned int i = 0; i < GHepFlags::NFlags(); i++) {
   fUnphysEventMask->SetBitNumber(i, true);
  }
  fMCJobStatusRefreshRate = 50;
  fEventRecordPrintLevel  = 3;
  fEventGeneratorList     = "Default";
  fXMLPath = "";
}
//____________________________________________________________________________
void RunOpt::SetTuneName(string tuneName)
{
  if ( tuneName == "Default" || tuneName == "" ) tuneName = gDefaultTune;
  if ( fTune ) {
    LOG("RunOpt",pNOTICE) << "RunOpt::SetTune() already had " << fTune->Name()
              << ", now being re-set to " << tuneName;
    delete fTune;
  }
  fTune = new TuneId( tuneName ) ;
}
//____________________________________________________________________________
void RunOpt::BuildTune()
{
  LOG("RunOpt",pINFO) << "Building tune "<<Tune()->Name();
  Tune()->Build() ;
  XSecSplineList::Instance()->SetCurrentTune( Tune()->Name() ) ;
}
//____________________________________________________________________________
void RunOpt::ReadFromCommandLine(int argc, char ** argv)
{
  LOG("RunOpt",pDEBUG) << "Reading "<<argc-1<<" command line arguments.";
  CmdLnArgParser parser(argc,argv);

  if( parser.OptionExists("enable-bare-xsec-pre-calc") ) {
    fEnableBareXSecPreCalc = true;
  } else
  if( parser.OptionExists("disable-bare-xsec-pre-calc") ) {
    fEnableBareXSecPreCalc = false;
  }

  if( parser.OptionExists("cache-file") ) {
    fCacheFile = parser.ArgAsString("cache-file");
  }

  if( parser.OptionExists("message-thresholds") ) {
    fMesgThresholds = parser.ArgAsString("message-thresholds");
  }

  if( parser.OptionExists("event-record-print-level") ) {
    fEventRecordPrintLevel = parser.ArgAsInt("event-record-print-level");
  }

  if( parser.OptionExists("mc-job-status-refresh-rate") ) {
    fMCJobStatusRefreshRate = TMath::Max(
        1, parser.ArgAsInt("mc-job-status-refresh-rate"));
  }

  if( parser.OptionExists("event-generator-list") ) {
    SetEventGeneratorList(parser.ArgAsString("event-generator-list"));
  }

  if (parser.OptionExists("xml-path")) {
    fXMLPath = parser.ArgAsString("xml-path");
  }

  if( parser.OptionExists("tune") ) {
    SetTuneName( parser.ArgAsString("tune") ) ;
  }
  else {
    SetTuneName( "Default" );
  }// else ( parser.OptionExists("tune") )

  if( parser.OptionExists("unphysical-event-mask") ) {
    const char * bitfield =
       parser.ArgAsString("unphysical-event-mask").c_str();
    unsigned int n = GHepFlags::NFlags();
    unsigned int i = 0;
    while (i < n) {
        bool flag = (bitfield[i]=='1');
        fUnphysEventMask->SetBitNumber(n-1-i,flag);
        i++;
     } //i
  }

}
//____________________________________________________________________________
std::string RunOpt::RunOptSyntaxString(bool include_generator_specific)
{
  std::ostringstream s;
  s << "\n"
    << "\n // command line args handled by RunOpt:"
    // for v3, all tunes should have a Default event-generator-list
    << "\n         [--event-generator-list list_name] // default \"Default\" "
    // G18_02a_00_000 is currently the default tune
    << "\n         [--tune tune_name]  // default \"" << gDefaultTune << "\" "
    << "\n         [--xml-path path]"
    << "\n         [--message-thresholds xml_file]";

  if (include_generator_specific) {
    // these options are only for generator applications
    //    << "\n // not all options used by all applications "
    s << "\n"
      << "\n         [--event-record-print-level level]"
      << "\n         [--mc-job-status-refresh-rate rate]"
      << "\n         [--cache-file root_file]"
      << "\n         [--enable-bare-xsec-pre-calc]"
      << "\n         [--disable-bare-xsec-pre-calc]"
      << "\n         [--unphysical-event-mask mask]"
      << "\n";
  }

  return s.str();
}
//____________________________________________________________________________
void RunOpt::Print(ostream & stream) const
{
  stream << "Global running options:";
  if ( fTune ) stream << "\n GENIE tune: " << *fTune;
  stream << "\n Event generator list: " << fEventGeneratorList;
  stream << "\n User-specified message thresholds : " << fMesgThresholds;
  stream << "\n Cache file : " << fCacheFile;
  stream << "\n Unphysical event mask (bits: "
         << GHepFlags::NFlags()-1 << " -> 0) : " << *fUnphysEventMask;
  stream << "\n Event record print level : " << fEventRecordPrintLevel;
  stream << "\n MC job status file refresh rate: " << fMCJobStatusRefreshRate;
  stream << "\n Pre-calculate all free-nucleon cross-sections? : "
         << ((fEnableBareXSecPreCalc) ? "Yes" : "No");

  if (fXMLPath.size()) {
    stream << "\n XMLPath over-ride : "<<fXMLPath;
  }

  stream << "\n";
}
//___________________________________________________________________________

} // genie namespace
