//____________________________________________________________________________
/*!

\program gevgen

\brief   A simple 'generic' GENIE v+A event generation driver (gevgen).

       It handles:
         a) event generation for a fixed init state (v+A) at fixed energy, or
         b) event generation for simple fluxes (specified either via some
            functional form, tabular text file or a ROOT histogram) and for
            simple 'geometries' (a target mix with its corresponding weights)

         See the GENIE manual for other apps handling experiment-specific
         event generation cases using the outputs of detailed neutrino flux
         simulations and realistic detector geometry descriptions.

         Syntax :
           gevgen [-h]
                  [-r run#]
                   -n nev
                   -e energy (or energy range)
                   -p neutrino_pdg
                   -t target_pdg
                  [-f flux_description]
                  [-o outfile_name]
                  [-w]
                  [--seed random_number_seed]
                  [--cross-sections xml_file]
                  [--event-generator-list list_name]
                  [--tune genie_tune]
                  [--message-thresholds xml_file]
                  [--unphysical-event-mask mask]
                  [--event-record-print-level level]
                  [--mc-job-status-refresh-rate  rate]
                  [--cache-file root_file]
                  [--xml-path config_xml_dir]

         Options :
           [] Denotes an optional argument.
           -h
              Prints-out help on using gevgen and exits.
           -n
              Specifies the number of events to generate.
           -r
              Specifies the MC run number.
           -e
              Specifies the neutrino energy.
              If what follows the -e option is a comma separated pair of values
              it will be interpreted as an energy range for the flux specified
              via the -f option (see below).
           -p
              Specifies the neutrino PDG code.
           -t
              Specifies the target PDG code (pdg format: 10LZZZAAAI) _or_ a target
              mix (pdg codes with corresponding weights) typed as a comma-separated
              list of pdg codes with the corresponding weight fractions in brackets,
              eg code1[fraction1],code2[fraction2],...
              For example, to use a target mix of 95% O16 and 5% H type:
              `-t 1000080160[0.95],1000010010[0.05]'.
           -f
              Specifies the neutrino flux spectrum.
              It can be any of:
              -- A function:
                 eg ` -f x*x+4*exp(-x)'
              -- A vector file:
                 The vector file should contain 2 columns corresponding to
                 energy,flux (see $GENIE/data/flux/ for few examples).
              -- A 1-D ROOT histogram (TH1D):
                 The general syntax is `-f /full/path/file.root,object_name'
           -o
              Specifies the name of the output file events will be saved in.
           -w
              Forces generation of weighted events.
              This option is relevant only if a neutrino flux is specified.
              Note that 'weighted' refers to the selection of the primary
              flux neutrino + target that were forced to interact. A weighting
              scheme for the generated kinematics of individual processes can
              still be in effect if enabled..
              ** Only use that option if you understand what it means **
           --seed
              Random number seed.
           --cross-sections
              Name (incl. full path) of an XML file with pre-computed
              cross-section values used for constructing splines.
           --event-generator-list
              List of event generators to load in event generation drivers.
              [default: "Default"].
           --tune
              Specifies a GENIE comprehensive neutrino interaction model tune.
              [default: "Default"].
           --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.
           --unphysical-event-mask
              Allows users to specify a 16-bit mask to allow certain types of
              unphysical events to be written in the output file.
              [default: all unphysical events are rejected]
           --event-record-print-level
              Allows users to set the level of information shown when the event
              record is printed in the screen. See GHepRecord::Print().
           --mc-job-status-refresh-rate
              Allows users to customize the refresh rate of the status file.
           --cache-file
              Allows users to specify a cache file so that the cache can be
              re-used in subsequent MC jobs.
           --xml-path
              A directory to load XML files from - overrides $GXMLPATH, and $GENIE/config

        ***  See the User Manual for more details and examples. ***

\author  Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created October 05, 2004

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cstdlib>
#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
#include <fenv.h> // for `feenableexcept`
#endif

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TH1.h>
#include <TF1.h>

#include "Framework/Conventions/XmlParserStatus.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/EventGen/GFluxI.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/EventGen/GMCJDriver.h"
#include "Framework/EventGen/GMCJMonitor.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpWriter.h"
#include "Framework/Ntuple/NtpMCFormat.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/Spline.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/CmdLnArgParser.h"

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#define __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
#include "Tools/Flux/GCylindTH1Flux.h"
#include "Tools/Flux/GMonoEnergeticFlux.h"
#include "Tools/Geometry/PointGeomAnalyzer.h"
#endif
#endif

using std::string;
using std::vector;
using std::map;
using std::ostringstream;

using namespace genie;
using namespace genie::controls;

void GetCommandLineArgs (int argc, char ** argv);
void Initialize         (void);
void PrintSyntax        (void);

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
void            GenerateEventsUsingFluxOrTgtMix();
GeomAnalyzerI * GeomDriver              (void);
GFluxI *        FluxDriver              (void);
GFluxI *        MonoEnergeticFluxDriver (void);
GFluxI *        TH1FluxDriver           (void);
#endif

void GenerateEventsAtFixedInitState (void);

//Default options (override them using the command line arguments):
int           kDefOptNevents   = 0;       // n-events to generate
NtpMCFormat_t kDefOptNtpFormat = kNFGHEP; // ntuple format
Long_t        kDefOptRunNu     = 0;       // default run number

//User-specified options:
int             gOptNevents;      // n-events to generate
double          gOptNuEnergy;     // neutrino E, or min neutrino energy in spectrum
double          gOptNuEnergyRange;// energy range in input spectrum
int             gOptNuPdgCode;    // neutrino PDG code
map<int,double> gOptTgtMix;       // target mix (each with its relative weight)
Long_t          gOptRunNu;        // run number
string          gOptFlux;         //
bool            gOptWeighted;     //
bool            gOptUsingFluxOrTgtMix = false;
long int        gOptRanSeed;      // random number seed
string          gOptInpXSecFile;  // cross-section splines
string          gOptOutFileName;  // Optional outfile name
string          gOptStatFileName; // Status file name, set if gOptOutFileName was set.

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc,argv);
  Initialize();

  // throw on NaNs and Infs...
#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
  //
  // Generate neutrino events
  //

  if(gOptUsingFluxOrTgtMix) {
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
        GenerateEventsUsingFluxOrTgtMix();
#else
  LOG("gevgen", pERROR)
    << "\n   To be able to generate neutrino events from a flux and/or a target mix"
    << "\n   you need to add the following config options at your GENIE installation:"
    << "\n   --enable-flux-drivers  --enable-geom-drivers \n" ;
#endif
  } else {
     GenerateEventsAtFixedInitState();
  }
  return 0;
}
//____________________________________________________________________________
void Initialize()
{

  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gevgen", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  // Initialization of random number generators, cross-section table,
  // messenger thresholds, cache file
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
  utils::app_init::CacheFile(RunOpt::Instance()->CacheFile());
  utils::app_init::RandGen(gOptRanSeed);
  utils::app_init::XSecTable(gOptInpXSecFile, false);

  // Set GHEP print level
  GHepRecord::SetPrintLevel(RunOpt::Instance()->EventRecordPrintLevel());
}
//____________________________________________________________________________
void GenerateEventsAtFixedInitState(void)
{
  int neutrino = gOptNuPdgCode;
  int target   = gOptTgtMix.begin()->first;
  double Ev    = gOptNuEnergy;
  TLorentzVector nu_p4(0.,0.,Ev,Ev); // px,py,pz,E (GeV)

  // Create init state
  InitialState init_state(target, neutrino);

  // Create/config event generation driver
  GEVGDriver evg_driver;
  evg_driver.SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
  evg_driver.SetUnphysEventMask(*RunOpt::Instance()->UnphysEventMask());
  evg_driver.Configure(init_state);

  // Initialize an Ntuple Writer
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);

  // If an output file name has been specified... use it
  if (!gOptOutFileName.empty()){
    ntpw.CustomizeFilename(gOptOutFileName);
  }
  ntpw.Initialize();


  // Create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);
  mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());

  // If a status file name has been given... use it
  if (!gOptStatFileName.empty()){
    mcjmonitor.CustomizeFilename(gOptStatFileName);
  }


  LOG("gevgen", pNOTICE)
    << "\n ** Will generate " << gOptNevents << " events for \n"
    << init_state << " at Ev = " << Ev << " GeV";

  // Generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while (ievent < gOptNevents) {
     LOG("gevgen", pNOTICE)
        << " *** Generating event............ " << ievent;

     // generate a single event
     EventRecord * event = evg_driver.GenerateEvent(nu_p4);

     if(!event) {
        LOG("gevgen", pNOTICE)
          << "Last attempt failed. Re-trying....";
        continue;
     }

     LOG("gevgen", pNOTICE)
        << "Generated Event GHEP Record: " << *event;

     // add event at the output ntuple, refresh the mc job monitor & clean up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     ievent++;
     delete event;
  }

  // Save the generated MC events
  ntpw.Save();
}
//____________________________________________________________________________

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
//............................................................................
void GenerateEventsUsingFluxOrTgtMix(void)
{
  // Get flux and geom drivers
  GFluxI *        flux_driver = FluxDriver();
  GeomAnalyzerI * geom_driver = GeomDriver();

  // Create the monte carlo job driver
  GMCJDriver * mcj_driver = new GMCJDriver;
  mcj_driver->SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
  mcj_driver->SetUnphysEventMask(*RunOpt::Instance()->UnphysEventMask());
  mcj_driver->UseFluxDriver(flux_driver);
  mcj_driver->UseGeomAnalyzer(geom_driver);
  mcj_driver->Configure();
  mcj_driver->UseSplines();
  if(!gOptWeighted)
        mcj_driver->ForceSingleProbScale();

  // Initialize an Ntuple Writer to save GHEP records into a TTree
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);

  // If an output file name has been specified... use it
  if (!gOptOutFileName.empty()){
    ntpw.CustomizeFilename(gOptOutFileName);
  }
  ntpw.Initialize();

  // Create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);
  mcjmonitor.SetRefreshRate(RunOpt::Instance()->MCJobStatusRefreshRate());

  // If a status file name has been given... use it
  if (!gOptStatFileName.empty()){
    mcjmonitor.CustomizeFilename(gOptStatFileName);
  }


  // Generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

     LOG("gevgen", pNOTICE) << " *** Generating event............ " << ievent;

     // generate a single event for neutrinos coming from the specified flux
     EventRecord * event = mcj_driver->GenerateEvent();

     LOG("gevgen", pNOTICE) << "Generated Event GHEP Record: " << *event;

     // add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     ievent++;
     delete event;
  }

  // Save the generated MC events
  ntpw.Save();

  delete flux_driver;
  delete geom_driver;
  delete mcj_driver;;
}
//____________________________________________________________________________
GeomAnalyzerI * GeomDriver(void)
{
// create a trivial point geometry with the specified target or target mix

  GeomAnalyzerI * geom_driver = new geometry::PointGeomAnalyzer(gOptTgtMix);
  return geom_driver;
}
//____________________________________________________________________________
GFluxI * FluxDriver(void)
{
// create & configure one of the generic flux drivers
//
  GFluxI * flux_driver = 0;

  if(gOptNuEnergyRange<0) flux_driver = MonoEnergeticFluxDriver();
  else flux_driver = TH1FluxDriver();

  return flux_driver;
}
//____________________________________________________________________________
GFluxI * MonoEnergeticFluxDriver(void)
{
//
//
  flux::GMonoEnergeticFlux * flux =
              new flux::GMonoEnergeticFlux(gOptNuEnergy, gOptNuPdgCode);
  GFluxI * flux_driver = dynamic_cast<GFluxI *>(flux);
  return flux_driver;
}
//____________________________________________________________________________
GFluxI * TH1FluxDriver(void)
{
//
//
  flux::GCylindTH1Flux * flux = new flux::GCylindTH1Flux;
  TH1D * spectrum = 0;

  int flux_entries = 100000;

  double emin = gOptNuEnergy;
  double emax = gOptNuEnergy+gOptNuEnergyRange;
  double de   = gOptNuEnergyRange;

  // check whether the input flux is a file or a functional form
  //
  bool input_is_text_file = ! gSystem->AccessPathName(gOptFlux.c_str());
  bool input_is_root_file = gOptFlux.find(".root") != string::npos &&
                            gOptFlux.find(",") != string::npos;
  if(input_is_text_file) {
    //
    // ** generate the flux histogram from the x,y pairs in the input text file
    //
    Spline * input_flux = new Spline(gOptFlux.c_str());
    int  n = 100;
    double estep = (emax-emin)/(n-1);
    double ymax  = -1, ry = -1, gy = -1, e = -1;
    for(int i=0; i<n; i++) {
      e = emin + i*estep;
      ymax = TMath::Max(ymax, input_flux->Evaluate(e));
    }
    ymax *= 1.3;

    RandomGen * r = RandomGen::Instance();
    spectrum  = new TH1D("spectrum","neutrino flux", 300, emin, emax);
    spectrum->SetDirectory(0);

    for(int ientry=0; ientry<flux_entries; ientry++) {
      bool accept = false;
      unsigned int iter=0;
      while(!accept) {
        iter++;
        if(iter > kRjMaxIterations) {
           LOG("gevgen", pFATAL) << "Couldn't generate a flux histogram";
           exit(1);
        }
        e = emin + de * r->RndGen().Rndm();
        gy = ymax * r->RndGen().Rndm();
        ry = input_flux->Evaluate(e);
        accept = gy < ry;
        if(accept) spectrum->Fill(e);
      }
    }
    delete input_flux;
  }
  else if(input_is_root_file) {
    //
    // ** extract specified flux histogram from the input root file
    //
    vector<string> fv = utils::str::Split(gOptFlux,",");
    assert(fv.size()==2);
    assert( !gSystem->AccessPathName(fv[0].c_str()) );

    LOG("gevgen", pNOTICE) << "Getting input flux from root file: " << fv[0];
    TFile * flux_file = new TFile(fv[0].c_str(),"read");

    LOG("gevgen", pNOTICE) << "Flux name: " << fv[1];
    TH1D * hst = (TH1D *)flux_file->Get(fv[1].c_str());
    assert(hst);

    LOG("gevgen", pNOTICE) << hst->GetEntries();

    // Copy in the flux histogram from the root file and remove bins outside the emin,emax range
    spectrum = (TH1D*)hst->Clone();
    spectrum->SetNameTitle("spectrum","neutrino_flux");
    spectrum->SetDirectory(0);
    for(int ibin = 1; ibin <= hst->GetNbinsX(); ibin++) {
      if(hst->GetBinLowEdge(ibin) + hst->GetBinWidth(ibin) > emax ||
         hst->GetBinLowEdge(ibin) < emin) {
        spectrum->SetBinContent(ibin, 0);
      }
    }

    LOG("gevgen", pNOTICE) << spectrum->GetEntries();

    flux_file->Close();
    delete flux_file;

    LOG("gevgen", pNOTICE) << spectrum->GetEntries();

  } else {
    //
    // ** generate the flux histogram from the input functional form
    //
    TF1 *  input_func = new TF1("input_func", gOptFlux.c_str(), emin, emax);
    spectrum  = new TH1D("spectrum","neutrino flux", 300, emin, emax);
    spectrum->SetDirectory(0);
    spectrum->FillRandom("input_func", flux_entries);
    delete input_func;
  }
  // save input flux

  TFile f("./input-flux.root","recreate");
  spectrum->Write();
  f.Close();

  TVector3 bdir (0,0,1);
  TVector3 bspot(0,0,0);

  flux->SetNuDirection      (bdir);
  flux->SetBeamSpot         (bspot);
  flux->SetTransverseRadius (-1);
  flux->AddEnergySpectrum   (gOptNuPdgCode, spectrum);

  GFluxI * flux_driver = dynamic_cast<GFluxI *>(flux);
  return flux_driver;
}
//............................................................................
#endif

//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen", pINFO) << "Parsing command line arguments";

  // Common run options. Set defaults and read.
  RunOpt::Instance()->EnableBareXSecPreCalc(true);
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // help?
  bool help = parser.OptionExists('h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen", pINFO) << "Reading number of events to generate";
    gOptNevents = parser.ArgAsInt('n');
  } else {
    LOG("gevgen", pINFO)
       << "Unspecified number of events to generate - Using default";
    gOptNevents = kDefOptNevents;
  }

  // run number
  if( parser.OptionExists('r') ) {
    LOG("gevgen", pINFO) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen", pINFO) << "Unspecified run number - Using default";
    gOptRunNu = kDefOptRunNu;
  }

  // Output file name
  if( parser.OptionExists('o') ) {
    LOG("gevgen", pINFO) << "Reading output file name";
    gOptOutFileName = parser.ArgAsString('o');

    gOptStatFileName = gOptOutFileName;
    // strip the output file format and replace with .status
    if (gOptOutFileName.find_last_of(".") != string::npos)
      gOptStatFileName =
        gOptStatFileName.substr(0, gOptOutFileName.find_last_of("."));
    gOptStatFileName .append(".status");
  }

  // flux functional form
  bool using_flux = false;
  if( parser.OptionExists('f') ) {
    LOG("gevgen", pINFO) << "Reading flux function";
    gOptFlux = parser.ArgAsString('f');
    using_flux = true;
  }

  if(parser.OptionExists('s')) {
    LOG("gevgen", pWARN)
      << "-s option no longer available. Please read the revised code documentation";
    gAbortingInErr = true;
    exit(1);
  }


  // generate weighted events option (only relevant if using a flux)
  gOptWeighted = parser.OptionExists('w');

  // neutrino energy
  if( parser.OptionExists('e') ) {
    LOG("gevgen", pINFO) << "Reading neutrino energy";
    string nue = parser.ArgAsString('e');

    // is it just a value or a range (comma separated set of values)
    if(nue.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(nue, ",");
       assert(nurange.size() == 2);
       double emin = atof(nurange[0].c_str());
       double emax = atof(nurange[1].c_str());
       assert(emax>emin && emin>=0);
       gOptNuEnergy      = emin;
       gOptNuEnergyRange = emax-emin;
       if(!using_flux) {
          LOG("gevgen", pWARN)
             << "No flux was specified but an energy range was input!";
          LOG("gevgen", pWARN)
             << "Events will be generated at fixed E = " << gOptNuEnergy << " GeV";
          gOptNuEnergyRange = -1;
       }
    } else {
       gOptNuEnergy       = atof(nue.c_str());
       gOptNuEnergyRange = -1;
    }
  } else {
    LOG("gevgen", pFATAL) << "Unspecified neutrino energy - Exiting";
    PrintSyntax();
    exit(1);
  }

  // neutrino PDG code
  if( parser.OptionExists('p') ) {
    LOG("gevgen", pINFO) << "Reading neutrino PDG code";
    gOptNuPdgCode = parser.ArgAsInt('p');
  } else {
    LOG("gevgen", pFATAL) << "Unspecified neutrino PDG code - Exiting";
    PrintSyntax();
    exit(1);
  }

  // target mix (their PDG codes with their corresponding weights)
  bool using_tgtmix = false;
  if( parser.OptionExists('t') ) {
    LOG("gevgen", pINFO) << "Reading target mix";
    string stgtmix = parser.ArgAsString('t');
    gOptTgtMix.clear();
    vector<string> tgtmix = utils::str::Split(stgtmix,",");
    if(tgtmix.size()==1) {
         int    pdg = atoi(tgtmix[0].c_str());
         double wgt = 1.0;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
    } else {
      using_tgtmix = true;
      vector<string>::const_iterator tgtmix_iter = tgtmix.begin();
      for( ; tgtmix_iter != tgtmix.end(); ++tgtmix_iter) {
         string tgt_with_wgt = *tgtmix_iter;
         string::size_type open_bracket  = tgt_with_wgt.find("[");
         string::size_type close_bracket = tgt_with_wgt.find("]");
         string::size_type ibeg = 0;
         string::size_type iend = open_bracket;
         string::size_type jbeg = open_bracket+1;
         string::size_type jend = close_bracket-1;
         int    pdg = atoi(tgt_with_wgt.substr(ibeg,iend).c_str());
         double wgt = atof(tgt_with_wgt.substr(jbeg,jend).c_str());
         LOG("Main", pNOTICE)
            << "Adding to target mix: pdg = " << pdg << ", wgt = " << wgt;
         gOptTgtMix.insert(map<int, double>::value_type(pdg, wgt));
      }//tgtmix_iter
    }//>1

  } else {
    LOG("gevgen", pFATAL) << "Unspecified target PDG code - Exiting";
    PrintSyntax();
    exit(1);
  }

  gOptUsingFluxOrTgtMix = using_flux || using_tgtmix;

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevgen", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  // input cross-section file
  if( parser.OptionExists("cross-sections") ) {
    LOG("gevgen", pINFO) << "Reading cross-section file";
    gOptInpXSecFile = parser.ArgAsString("cross-sections");
  } else {
    LOG("gevgen", pINFO) << "Unspecified cross-section file";
    gOptInpXSecFile = "";
  }

  //
  // print-out the command line options
  //
  LOG("gevgen", pNOTICE)
     << "\n"
     << utils::print::PrintFramedMesg("gevgen job configuration");
  LOG("gevgen", pNOTICE)
     << "MC Run Number: " << gOptRunNu;
  if(gOptRanSeed != -1) {
     LOG("gevgen", pNOTICE)
       << "Random number seed: " << gOptRanSeed;
  } else {
     LOG("gevgen", pNOTICE)
       << "Random number seed was not set, using default";
  }
  LOG("gevgen", pNOTICE)
       << "Number of events requested: " << gOptNevents;
  if(gOptInpXSecFile.size() > 0) {
     LOG("gevgen", pNOTICE)
       << "Using cross-section splines read from: " << gOptInpXSecFile;
  } else {
     LOG("gevgen", pNOTICE)
       << "No input cross-section spline file";
  }
  LOG("gevgen", pNOTICE)
       << "Flux: " << gOptFlux;
  LOG("gevgen", pNOTICE)
       << "Generate weighted events? " << gOptWeighted;
  if(gOptNuEnergyRange>0) {
     LOG("gevgen", pNOTICE)
        << "Neutrino energy: ["
        << gOptNuEnergy << ", " << gOptNuEnergy+gOptNuEnergyRange << "]";
  } else {
     LOG("gevgen", pNOTICE)
        << "Neutrino energy: " << gOptNuEnergy;
  }
  LOG("gevgen", pNOTICE)
      << "Neutrino code (PDG): " << gOptNuPdgCode;
  LOG("gevgen", pNOTICE)
      << "Target code (PDG) & weight fraction (in case of multiple targets): ";
  map<int,double>::const_iterator iter;
  for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
      int    tgtpdgc = iter->first;
      double wgt     = iter->second;
      LOG("gevgen", pNOTICE)
          << " >> " <<  tgtpdgc << " (weight fraction = " << wgt << ")";
  }
  LOG("gevgen", pNOTICE) << "\n";

  LOG("gevgen", pNOTICE) << *RunOpt::Instance();

}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "\n      gevgen [-h]"
    << "\n              [-r run#]"
    << "\n               -n nev"
    << "\n               -e energy (or energy range) "
    << "\n               -p neutrino_pdg"
    << "\n               -t target_pdg "
    << "\n              [-f flux_description]"
    << "\n              [-o outfile_name]"
    << "\n              [-w]"
    << "\n              [--seed random_number_seed]"
    << "\n              [--cross-sections xml_file]"
    << "\n              [--event-generator-list list_name]"
    << "\n              [--message-thresholds xml_file]"
    << "\n              [--unphysical-event-mask mask]"
    << "\n              [--event-record-print-level level]"
    << "\n              [--mc-job-status-refresh-rate  rate]"
    << "\n              [--cache-file root_file]"
    << "\n              [--xml-path config_xml_dir]"
    << "\n";
}
//____________________________________________________________________________
