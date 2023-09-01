//____________________________________________________________________________
/*
 Copyright (c) 2003-2023, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
 University of Liverpool & STFC Rutherford Appleton Laboratory
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

  desc = std::shared_ptr<boost::program_options::options_description>(new boost::program_options::options_description("RunOpt"));


  desc->add_options()
    ("enable-bare-xsec-pre-calc", po::bool_switch(&fEnableBareXSecPreCalc)->default_value(false), "")
    ("disable-bare-xsec-pre-calc", po::bool_switch()->default_value(false), "")
    ("cache-file", po::value<std::string>(&fCacheFile)->default_value(""), "")
    ("message-thresholds", po::value<std::string>(&fMesgThresholds)->default_value(""), "")
    ("event-record-print-level", po::value<int>(&fEventRecordPrintLevel)->default_value(0), "")
    ("mc-job-status-refresh=rate", po::value<int>(&fMCJobStatusRefreshRate)->default_value(1), "")
    ("event-generator-list", po::value<std::string>()->default_value(""), "")
    ("xml-path", po::value<std::string>(&fXMLPath)->default_value(""), "");
    ("tune", po::value<std::string>()->default_value("Default"), "");
}
//____________________________________________________________________________
std::string RunOpt::HelpString() {
  std::stringstream ss;
  // create some dummy args
  char* argv[] = {"RunOpt",  NULL};
  po::variables_map vm;
  po::store(po::parse_command_line(1, argv, (*desc) ), vm);
  ss << (*desc);
  return ss.str();
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

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(*desc).allow_unregistered().run(), vm);
  po::notify(vm);

  fMCJobStatusRefreshRate = TMath::Max(1, fMCJobStatusRefreshRate);
  if (vm.count("enable-bare-xsec-pre-calc") &&
      !vm["enable-bare-xsec-pre-calc"].defaulted() &&
      vm.count("disable-bare-xsec-pre-calc") &&
      !vm["disable-bare-xsec-pre-calc"].defaulted()) {
    // LOG("RunOpt", pERROR) << "Inconsistent use of bare-xsec-pre-calc!\n";
    exit(-1);
  }
  else if (vm.count("disable-bare-xsec-pre-calc") &&
      vm["disable-bare-xsec-pre-calc"].as<bool>())
    fEnableBareXSecPreCalc = false;

  if( vm.count("event-generator-list") ) {
    SetEventGeneratorList(vm["event-generator-list"].as<std::string>());
  }

  if( vm.count("tune") ) {
    SetTuneName( vm["tune"].as<std::string>() ) ;
  }

  if( vm.count("unphysical-event-mask") ) {
    const char * bitfield =
       vm["unphysical-event-mask"].as<std::string>().c_str();
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
