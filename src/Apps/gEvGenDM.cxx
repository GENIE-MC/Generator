//____________________________________________________________________________
/*!

\program gevgen_dm

\brief   A simple 'generic' GENIE DM+A event generation driver (gevgen_dm).

         It handles:
         a) event generation for a fixed init state (DM+A) at fixed energy, or
         b) event generation for simple fluxes (specified either via some
            functional form, tabular text file or a ROOT histogram) and for
            simple 'geometries' (a target mix with its corresponding weights)

         See the GENIE manual for other apps handling experiment-specific
         event generation cases using the outputs of detailed dark matter flux
         simulations and realistic detector geometry descriptions.

         Syntax :
           gevgen_dm [-h]
                     [-r run#]
                     -n nev
                     -e energy (or energy range)
                     -m mass
                     -t target_pdg
                     [-g zp_coupling]
                     [-z med_ratio]
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

         Options :
           [] Denotes an optional argument.
           -h
              Prints-out help on using gevgen_dm and exits.
           -n
              Specifies the number of events to generate.
           -r
              Specifies the MC run number.
           -e
              Specifies the dark matter energy.
              If what follows the -e option is a comma separated pair of values
              it will be interpreted as an energy range for the flux specified
              via the -f option (see below).
           -m
              Specifies the dark matter mass.
           -t
              Specifies the target PDG code (pdg format: 10LZZZAAAI) _or_ a target
              mix (pdg codes with corresponding weights) typed as a comma-separated
              list of pdg codes with the corresponding weight fractions in brackets,
              eg code1[fraction1],code2[fraction2],...
              For example, to use a target mix of 95% O16 and 5% H type:
              `-t 1000080160[0.95],1000010010[0.05]'.
           -z
              Specifies the ratio of the mediator mass to dark matter mass.
              Default: 0.5
           -g
              Specifies the Z' coupling constant
              Default: Value in UserPhysicsOptions.xml
           -f
              Specifies the dark matter flux spectrum.
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
              This option is relevant only if a dark matter flux is specified.
              Note that 'weighted' refers to the selection of the primary
              flux dark matter + target that were forced to interact. A weighting
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
              Specifies a GENIE comprehensive dark matter interaction model tune.
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

        ***  See the User Manual for more details and examples. ***

\author  Joshua Berger <jberger \at physics.wisc.edu>
         University of Wisconsin-Madison
         Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created September 1, 2017

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


#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/XmlParserStatus.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Controls.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
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
#include "Framework/ParticleData/PDGLibrary.h"
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
using namespace genie::constants;
using namespace genie::units;

void GetCommandLineArgs (int argc, char ** argv);
void Initialize         (void);
void PrintSyntax        (void);
bool CheckUnitarityLimit(InitialState init_state);

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
double          gOptDMEnergy;     // dark matter E, or min dark matter energy in spectrum
double          gOptDMEnergyRange;// energy range in input spectrum
double          gOptDMMass;       // dark matter mass
double          gOptZpCoupling;   // mediator coupling
map<int,double> gOptTgtMix;       // target mix (each with its relative weight)
double          gOptMedRatio;     // ratio of mediator to DM mass
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
  PDGLibrary::Instance()->AddDarkMatter(gOptDMMass,gOptMedRatio);
  if (gOptZpCoupling > 0.) {
      Registry * r = AlgConfigPool::Instance()->CommonList("Param", "BoostedDarkMatter");
      r->UnLock();
      r->Set("ZpCoupling", gOptZpCoupling);
      r->Lock();
  }
  Initialize();


  // throw on NaNs and Infs...
#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
  //
  // Generate dark matter events
  //

  if(gOptUsingFluxOrTgtMix) {
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
        GenerateEventsUsingFluxOrTgtMix();
#else
  LOG("gevgen_dm", pERROR)
    << "\n   To be able to generate dark matter events from a flux and/or a target mix"
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
    LOG("gmkspl", pFATAL) << " No TuneId in RunOption";
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
  int dark_matter = kPdgDarkMatter;
  int target   = gOptTgtMix.begin()->first;
  double Ed    = gOptDMEnergy;
  double Md    = gOptDMMass;
  double pd    = TMath::Sqrt(Ed*Ed - Md*Md);
  assert(pd>=0.);
  TLorentzVector dm_p4(0.,0.,pd,Ed); // px,py,pz,E (GeV)

  // Create init state
  InitialState init_state(target, dark_matter);

  bool unitary = CheckUnitarityLimit(init_state);
  if (!unitary) {
    LOG("gevgen_dm", pFATAL)
      << "Cross-section risks exceeding unitarity limit - Exiting";
    exit(1);
  }


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


  LOG("gevgen_dm", pNOTICE)
    << "\n ** Will generate " << gOptNevents << " events for \n"
    << init_state << " at Ev = " << Ed << " GeV";

  // Generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while (ievent < gOptNevents) {
    LOG("gevgen_dm", pNOTICE)
      << " *** Generating event............ " << ievent;

    // generate a single event
    EventRecord * event = evg_driver.GenerateEvent(dm_p4);

    if(!event) {
      LOG("gevgen_dm", pNOTICE)
        << "Last attempt failed. Re-trying....";
      continue;
    }

    LOG("gevgen_dm", pNOTICE)
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

     LOG("gevgen_dm", pNOTICE) << " *** Generating event............ " << ievent;

     // generate a single event for dark matter particles coming from the specified flux
     EventRecord * event = mcj_driver->GenerateEvent();

     LOG("gevgen_dm", pNOTICE) << "Generated Event GHEP Record: " << *event;

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

  if(gOptDMEnergyRange<0) flux_driver = MonoEnergeticFluxDriver();
  else flux_driver = TH1FluxDriver();

  return flux_driver;
}
//____________________________________________________________________________
GFluxI * MonoEnergeticFluxDriver(void)
{
//
//
  flux::GMonoEnergeticFlux * flux =
              new flux::GMonoEnergeticFlux(gOptDMEnergy, kPdgDarkMatter);
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

  double emin = gOptDMEnergy;
  double emax = gOptDMEnergy+gOptDMEnergyRange;
  double de   = gOptDMEnergyRange;

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
    spectrum  = new TH1D("spectrum","dark matter flux", 300, emin, emax);
    spectrum->SetDirectory(0);

    for(int ientry=0; ientry<flux_entries; ientry++) {
      bool accept = false;
      unsigned int iter=0;
      while(!accept) {
        iter++;
        if(iter > kRjMaxIterations) {
           LOG("gevgen_dm", pFATAL) << "Couldn't generate a flux histogram";
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

    LOG("gevgen_dm", pNOTICE) << "Getting input flux from root file: " << fv[0];
    TFile * flux_file = new TFile(fv[0].c_str(),"read");

    LOG("gevgen_dm", pNOTICE) << "Flux name: " << fv[1];
    TH1D * hst = (TH1D *)flux_file->Get(fv[1].c_str());
    assert(hst);

    LOG("gevgen_dm", pNOTICE) << hst->GetEntries();

    // Copy in the flux histogram from the root file and remove bins outside the emin,emax range
    spectrum = (TH1D*)hst->Clone();
    spectrum->SetNameTitle("spectrum","dark_matter_flux");
    spectrum->SetDirectory(0);
    for(int ibin = 1; ibin <= hst->GetNbinsX(); ibin++) {
      if(hst->GetBinLowEdge(ibin) + hst->GetBinWidth(ibin) > emax ||
         hst->GetBinLowEdge(ibin) < emin) {
        spectrum->SetBinContent(ibin, 0);
      }
    }

    LOG("gevgen_dm", pNOTICE) << spectrum->GetEntries();

    flux_file->Close();
    delete flux_file;

    LOG("gevgen_dm", pNOTICE) << spectrum->GetEntries();

  } else {
    //
    // ** generate the flux histogram from the input functional form
    //
    TF1 *  input_func = new TF1("input_func", gOptFlux.c_str(), emin, emax);
    spectrum  = new TH1D("spectrum","dark matter flux", 300, emin, emax);
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
  flux->AddEnergySpectrum   (kPdgDarkMatter, spectrum);

  GFluxI * flux_driver = dynamic_cast<GFluxI *>(flux);
  return flux_driver;
}
//............................................................................
#endif

//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen_dm", pINFO) << "Parsing command line arguments";

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

  if ( ! parser.OptionExists("tune") ) {
    LOG("gevgen_dm", pFATAL) << "No Dark Matter tune selected, please select one ";
    LOG("gevgen_dm", pFATAL) << "Exiting ";
    exit( 0 ) ;
  }

  // number of events
  if( parser.OptionExists('n') ) {
    LOG("gevgen_dm", pINFO) << "Reading number of events to generate";
    gOptNevents = parser.ArgAsInt('n');
  } else {
    LOG("gevgen_dm", pINFO)
       << "Unspecified number of events to generate - Using default";
    gOptNevents = kDefOptNevents;
  }

  // run number
  if( parser.OptionExists('r') ) {
    LOG("gevgen_dm", pINFO) << "Reading MC run number";
    gOptRunNu = parser.ArgAsLong('r');
  } else {
    LOG("gevgen_dm", pINFO) << "Unspecified run number - Using default";
    gOptRunNu = kDefOptRunNu;
  }

  // Output file name
  if( parser.OptionExists('o') ) {
    LOG("gevgen_dm", pINFO) << "Reading output file name";
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
    LOG("gevgen_dm", pINFO) << "Reading flux function";
    gOptFlux = parser.ArgAsString('f');
    using_flux = true;
  }

  if(parser.OptionExists('s')) {
    LOG("gevgen_dm", pWARN)
      << "-s option no longer available. Please read the revised code documentation";
    gAbortingInErr = true;
    exit(1);
  }


  // generate weighted events option (only relevant if using a flux)
  gOptWeighted = parser.OptionExists('w');

  // dark matter energy
  if( parser.OptionExists('e') ) {
    LOG("gevgen_dm", pINFO) << "Reading dark matter energy";
    string dme = parser.ArgAsString('e');

    // is it just a value or a range (comma separated set of values)
    if(dme.find(",") != string::npos) {
       // split the comma separated list
       vector<string> nurange = utils::str::Split(dme, ",");
       assert(nurange.size() == 2);
       double emin = atof(nurange[0].c_str());
       double emax = atof(nurange[1].c_str());
       assert(emax>emin && emin>=0);
       gOptDMEnergy      = emin;
       gOptDMEnergyRange = emax-emin;
       if(!using_flux) {
          LOG("gevgen_dm", pWARN)
             << "No flux was specified but an energy range was input!";
          LOG("gevgen_dm", pWARN)
             << "Events will be generated at fixed E = " << gOptDMEnergy << " GeV";
          gOptDMEnergyRange = -1;
       }
    } else {
       gOptDMEnergy       = atof(dme.c_str());
       gOptDMEnergyRange = -1;
    }
  } else {
    LOG("gevgen_dm", pFATAL) << "Unspecified dark matter energy - Exiting";
    PrintSyntax();
    exit(1);
  }

  // dark matter mass
  if( parser.OptionExists('m') ) {
    LOG("gevgen_dm", pINFO) << "Reading dark matter mass";
    gOptDMMass = parser.ArgAsDouble('m');
  } else {
    LOG("gevgen_dm", pFATAL) << "Unspecified dark matter mass - Exiting";
    PrintSyntax();
    exit(1);
  }

  // mediator coupling
  if( parser.OptionExists('g') ) {
    LOG("gevgen_dm", pINFO) << "Reading mediator coupling";
    gOptZpCoupling = parser.ArgAsDouble('g');
  } else {
    LOG("gevgen_dm", pINFO) << "Unspecified mediator coupling - Using value from config file";
    gOptZpCoupling = -1.;
  }

  // target mix (their PDG codes with their corresponding weights)
  bool using_tgtmix = false;
  if( parser.OptionExists('t') ) {
    LOG("gevgen_dm", pINFO) << "Reading target mix";
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
    LOG("gevgen_dm", pFATAL) << "Unspecified target PDG code - Exiting";
    PrintSyntax();
    exit(1);
  }

  // mediator mass ratio
  if( parser.OptionExists('z') ) {
    LOG("gevgen_dm", pINFO) << "Reading mediator mass ratio";
    gOptMedRatio = parser.ArgAsDouble('z');
  } else {
    LOG("gevgen_dm", pINFO) << "Unspecified mediator mass ratio - Using default";
    gOptMedRatio = 0.5;
  }

  gOptUsingFluxOrTgtMix = using_flux || using_tgtmix;

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gevgen_dm", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gevgen_dm", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  // input cross-section file
  if( parser.OptionExists("cross-sections") ) {
    LOG("gevgen_dm", pINFO) << "Reading cross-section file";
    gOptInpXSecFile = parser.ArgAsString("cross-sections");
  } else {
    LOG("gevgen_dm", pINFO) << "Unspecified cross-section file";
    gOptInpXSecFile = "";
  }

  //
  // print-out the command line options
  //
  LOG("gevgen_dm", pNOTICE)
     << "\n"
     << utils::print::PrintFramedMesg("gevgen_dm job configuration");
  LOG("gevgen_dm", pNOTICE)
     << "MC Run Number: " << gOptRunNu;
  if(gOptRanSeed != -1) {
     LOG("gevgen_dm", pNOTICE)
       << "Random number seed: " << gOptRanSeed;
  } else {
     LOG("gevgen_dm", pNOTICE)
       << "Random number seed was not set, using default";
  }
  LOG("gevgen_dm", pNOTICE)
       << "Number of events requested: " << gOptNevents;
  if(gOptInpXSecFile.size() > 0) {
     LOG("gevgen_dm", pNOTICE)
       << "Using cross-section splines read from: " << gOptInpXSecFile;
  } else {
     LOG("gevgen_dm", pNOTICE)
       << "No input cross-section spline file";
  }
  LOG("gevgen_dm", pNOTICE)
       << "Flux: " << gOptFlux;
  LOG("gevgen_dm", pNOTICE)
       << "Generate weighted events? " << gOptWeighted;
  if(gOptDMEnergyRange>0) {
     LOG("gevgen_dm", pNOTICE)
        << "Dark matter energy: ["
        << gOptDMEnergy << ", " << gOptDMEnergy+gOptDMEnergyRange << "]";
  } else {
     LOG("gevgen_dm", pNOTICE)
        << "Dark matter energy: " << gOptDMEnergy;
  }
  LOG("gevgen_dm", pNOTICE)
      << "Dark matter mass: " << gOptDMMass;
  LOG("gevgen_dm", pNOTICE)
      << "Target code (PDG) & weight fraction (in case of multiple targets): ";
  LOG("gevgen_dm", pNOTICE)
      << "Mediator mass ratio: " << gOptMedRatio;
  map<int,double>::const_iterator iter;
  for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
      int    tgtpdgc = iter->first;
      double wgt     = iter->second;
      LOG("gevgen_dm", pNOTICE)
          << " >> " <<  tgtpdgc << " (weight fraction = " << wgt << ")";
  }
  LOG("gevgen_dm", pNOTICE) << "\n";

  LOG("gevgen_dm", pNOTICE) << *RunOpt::Instance();

}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen_dm", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "\n      gevgen_dm [-h]"
    << "\n                [-r run#]"
    << "\n                -n nev"
    << "\n                -e energy (or energy range) "
    << "\n                -m mass"
    << "\n                -t target_pdg "
    << "\n                [-g zp_coupling]"
    << "\n                [-z med_ratio]"
    << "\n                [-f flux_description]"
    << "\n                [-o outfile_name]"
    << "\n                [-w]"
    << "\n                [--seed random_number_seed]"
    << "\n                [--cross-sections xml_file]"
    << "\n                [--event-generator-list list_name]"
    << "\n                [--message-thresholds xml_file]"
    << "\n                [--unphysical-event-mask mask]"
    << "\n                [--event-record-print-level level]"
    << "\n                [--mc-job-status-refresh-rate  rate]"
    << "\n                [--cache-file root_file]"
    << "\n";
}
//____________________________________________________________________________
bool CheckUnitarityLimit(InitialState init_state)
{
  // Before generating the events, perform a simple sanity check
  // We estimate the leading divergent piece of the cross-section
  // We make sure it does not exceed the unitarity limit
  double gzp;
  Registry * r = AlgConfigPool::Instance()->CommonList("Param", "BoostedDarkMatter");
  r->Get("ZpCoupling", gzp);
  double gzp4 = TMath::Power(gzp,4);
  double Mzp  = gOptMedRatio * gOptDMMass;
  double Mzp2 = Mzp*Mzp;
  // The leading, forward-dominated piece is the same for both DM models
  double xsec_est = gzp4 / (4. * kPi * Mzp2);
  double ml   = gOptDMMass;
  double ml2  = ml*ml;
  double M    = kNucleonMass;
  double M2   = M*M;
  double Ed   = gOptDMEnergy;
  double Ed2  = Ed*Ed;
  double pcm2 = M2 * (Ed2 - ml2) / (ml2 + M2 + 2.*M*Ed);
  double xsec_lim = kPi / pcm2;
  bool unitary = xsec_lim > xsec_est;
  if (!unitary) {
    LOG("gevgen_dm", pWARN)
      << "Estimated a cross-section " << xsec_est/cm2 << " cm^2";
    LOG("gevgen_dm", pWARN)
      << "Unitarity limit set to " << xsec_lim/cm2 << " cm^2";
  }
  return unitary;
}
//____________________________________________________________________________
