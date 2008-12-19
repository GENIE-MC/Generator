//____________________________________________________________________________
/*!

\program gevgen

\brief   A 'generic' GENIE v+A event generation driver (gevgen)

	 This generic event generation driver can handle:
 	 a) event generation for a fixed init state (v+A) at fixed energy, or
         b) event generation for simple fluxes (specified either via some
            functional form, tabular text file or a ROOT histogram) and for 
            simple 'geometries' (a target mix with its corresponding weights)

         For more complex event generation cases using realistic / detailed 
         neutrino flux simulations and detector geometries see for example 
         the GENIE event generation driver for T2K (gT2Kevgen) at
         $GENIE/src/support/t2k/EvGen/gT2Kevgen.cxx
	 Use that as a template for your own experiment-specific case.

         Syntax :
           gevgen [-h] -n nev [-s] -e E -p nupdg -t tgtpdg [-r run#] [-f flux] [-w]

         Options :
           [] Denotes an optional argument
           -h Prints help and exits
           -n Specifies the number of events to generate
           -r Specifies the MC run number
           -s Turns on cross section spline generation at job initialization.
              ** Always use that option **
              Not using the cross section splines will slow down the event
              generation considerably. Since building the splines is a 
              considerable overhead we recommend building them in advance
              (see gmkspl utility) and loading via the $GSPLOAD env. variable.
           -e Specifies the neutrino energy
	      If what follows the -e option is a comma separated pair of values
              it will be interpreted as an energy range for the flux specified
              via the -f option (see below).
           -p Specifies the neutrino PDG code
           -t Specifies the target PDG code (pdg format: 10LZZZAAAI) _or_ a target
              mix (pdg codes with corresponding weights) typed as a comma-separated 
              list of pdg codes with the corresponding weight fractions in brackets, 
              eg code1[fraction1],code2[fraction2],... For example, to use a target
              mix of 95% O16 and 5% H type: '-t 1000080160[0.95],1000010010[0.05]'.
              See the examples below.
           -f Specifies the neutrino flux spectrum - it can be either:
	      -- A function, eg 'x*x+4*exp(-x)' 
              -- A text file containing 2 columns corresponding to energy,flux
                 (see $GENIE/data/flux/ for few examples). 
              -- A 1-D histogram taken from a ROOT file. It is specified as
                 /full/path/file.root,object_name
              If this option is enabled, the GMCJDriver is partially utilized 
              to the use the specified flux and a trivial 'point geometry'. 
              See $GENIE/src/test/testMCJobDriver.cxx to see how you can generate 
              events for arbitrarily complex fluxes and detector geometries. 
              Note that in order to use the -f option you must have enabled the 
              flux and geometry drivers during the GENIE installation (see 
              installation instructions).
           -w Forces generation of weighted events.
              This option is relevant only if a neutrino flux is specified.
              Note that 'weighted' refers to the selection of the primary
              flux neutrino + target that were forced to interact. A weighting
              scheme for the generated kinematics of individual processes can
              still be in effect if enabled..
              ** Only use that option if you understand what it means **
 
         <Examples>

          (1) gevgen -n 3000 -s -e 6.5 -p 14 -t 1000260560 

           will generate 3000 events of muon neutrinos (pdg=14) on Iron (A=56,Z=26) 
           at E = 6.5 GeV. The -s option will force it to generate cross section
           splines during the job initialization. Only the cross section splines
           not loaded via $GSPLOAD (see below) will be generated. Note that building
           cross section splines for all necessary processes is a time consuming step
           and we do recommend building them in advance (using the gmkspl utility)
           and loading them via $GSPLOAD.

          (2) gevgen -n 30000 -s -e 1,4 -p 14 -t 1000080160 -f 'x*x*exp((-x*x+3)/4)' 

           will generate 30000 _unweighted_ events of muon neutrinos (pdg=14) on O16
           (A=16,Z=8) using a neutrino flux of the specified functional form at the 
           energy range from 1 to 4 GeV. 
           Above comments on cross section splines apply here too.

          (3) gevgen -n 30000 -s -e 1,4 -p 14 -t 1000080160 -f /a/path/flux.data 

           like above except that the neutrino flux is taken from the input vector file 
           rather than the input functional form.

          (4) gevgen -n 30000 -s -e 1,4 -p 14 -t 1000080160 -f /a/path/file.root,flux

           like above except that the flux is taken from a TH1D object called 'flux'
           stored in /a/path/file.root. 
           Note that the event generation driver will use only the input histogram 
           bins that fall within the specified (with -e) energy range. In the above 
           example, all the neutrino flux bins that do not fall in the 1 GeV - 4 GeV 
           range will be neglected (Note that the bins including 1 GeV and 4 GeV will 
           be taken into account - So the actual energy range used is: from the lower 
           edge of the bin containing 1 GeV to the upper edge of the bin containing 
           4 GeV).
           
          (5) gevgen -n 30000 -s -e 1,4 -p 14 -t 1000080160[0.95],1000010010[0.05] 
                      -f /a/path/file.root,flux 

           like above but in this case the target is a mix containing 95% O16 and
           5% H.

         Other control options:

         You can further control the program behaviour by setting the GEVGL,
         GSPLOAD, GSPSAVE, GHEPPRINTLEVEL, GSEED environmental variables.

           - Set the GEVGL environmental variable to contol the list of event
             generator objects that get loaded (at job initialization, the
             EventGeneratorListAssember assembles the list of all requested
             event generator objects -meaning: all concrete implementations of
             the EventGeneratorI interface that know how to generate various
             classes of events-). You can see the possible GEVGL options by
             looking up the names of the EventGeneratorListAssembler configuration
             sets at $GENIE/config/event_generator_list_assembler.xml
             If GEVGL is not set, the 'Default' generator list is used.

           - If you specify the -s option you can set the GSPSAVE env.variable
             to specify an XML filename for saving the computed xsec splines.
             If GSPSAVE is not set the computed splines will not be saved.

           - If you specify the -s option you can set the GSPLOAD env.variable
             to specify an XML filename from which to load xsec splines at
             initialization. GENIE will load the splines and then go on and
             compute only the xsec splines that it needs but were not loaded.
             If GSPLOAD is not set, GENIE will not attempt to load any splines
             and will compute them all.

           - You can set the GMSGCONF env.variable to point to a messenger
             XML configuration file (following the syntax of the default one
             that can be found in $GENIE/config/messenger.xml) and modify the
             verbosity of GENIE output. Both the default and your messenger
             configuration will be read but yours will take precedence in
             case of clashing priority for the same stream.

           - You can set the GHEPPRINTLEVEL to control the GHEP print options

           - You can set the GSEED env variable to define the random number
             seed number (default seed number in (src/Conventions/Controls.h)

             Examples:
             By setting the following env.vars you ask GENIE to generate QEL
             events only, and load the cross section splines from splines.xml
             and also to add the priority levels in mymsg.xml on top of the
             default ones.
             shell% export GEVGL=QEL
             shell% export GSPLOAD=/home/me/GENIE/mydata/splines.xml
             shell% export GMSGCONF=/home/me/GENIE/mydata/mymsg.xml

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory

\created October 05, 2004

\cpright Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TVector3.h>
#include <TH1.h>
#include <TF1.h>

#include "Conventions/XmlParserStatus.h"
#include "Conventions/GBuild.h"
#include "Conventions/Controls.h"
#include "EVGCore/EventRecord.h"
#include "EVGDrivers/GFluxI.h"
#include "EVGDrivers/GEVGDriver.h"
#include "EVGDrivers/GMCJDriver.h"
#include "EVGDrivers/GMCJMonitor.h"
#include "Interaction/Interaction.h"
#include "Messenger/Messenger.h"
#include "Ntuple/NtpWriter.h"
#include "Ntuple/NtpMCFormat.h"
#include "Numerical/RandomGen.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "Utils/XSecSplineList.h"
#include "Utils/StringUtils.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"

#ifdef __GENIE_FLUX_DRIVERS_ENABLED__
#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#define __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
#include "FluxDrivers/GCylindTH1Flux.h"
#include "FluxDrivers/GMonoEnergeticFlux.h"
#include "Geo/PointGeomAnalyzer.h"
#endif
#endif

using std::string;
using std::vector;
using std::map;
using std::ostringstream;

using namespace genie;
using namespace genie::controls;

void GetCommandLineArgs (int argc, char ** argv);
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
bool            gOptBuildSplines; // spline building option
double          gOptNuEnergy;     // neutrino E, or min neutrino energy in spectrum
double          gOptNuEnergyRange;// energy range in input spectrum
int             gOptNuPdgCode;    // neutrino PDG code
map<int,double> gOptTgtMix;       // target mix (each with its relative weight)
Long_t          gOptRunNu;        // run number
string          gOptFlux;         //
bool            gOptWeighted;     // 
bool            gOptUsingFluxOrTgtMix = false;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);
  
  //-- Autoload splines (from the XML file pointed at the $GSPLOAD env. var.,
  //   if the env. var. has been set)
  XSecSplineList * xspl = XSecSplineList::Instance();
  xspl->AutoLoad();

  //-- Generate neutrino events
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
  } else 
	GenerateEventsAtFixedInitState();

  return 0;
}
//____________________________________________________________________________
void GenerateEventsAtFixedInitState(void)
{
  int neutrino = gOptNuPdgCode;
  int target   = gOptTgtMix.begin()->first;
  double Ev    = gOptNuEnergy;
  TLorentzVector nu_p4(0.,0.,Ev,Ev); // px,py,pz,E (GeV)

  //-- Create init state
  InitialState init_state(target, neutrino);

  //-- Create/config event generation driver 
  GEVGDriver evg_driver;
  evg_driver.Configure(init_state);

  //-- If splines are used, then create any spline that is needed but is not
  //   already loaded (can built them all here if no spline at all is loaded)
  if(gOptBuildSplines) 
	evg_driver.CreateSplines();

  //-- initialize an Ntuple Writer
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);
  ntpw.Initialize();

  //-- create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  LOG("gevgen", pNOTICE) 
    << "\n ** Will generate " << gOptNevents << " events for \n" 
    << init_state << " at Ev = " << Ev << " GeV";

  //-- generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while (ievent < gOptNevents) {
     LOG("gevgen", pNOTICE) 
        << " *** Generating event............ " << ievent;

     // generate a single event
     EventRecord * event = evg_driver.GenerateEvent(nu_p4);
     LOG("gevgen", pINFO) 
	<< "Generated Event GHEP Record: " << *event;

     // add event at the output ntuple, refresh the mc job monitor & clean up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     ievent++;
     delete event;
  }

  //-- save the generated MC events
  ntpw.Save();
}
//____________________________________________________________________________

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX_OR_TGTMIX__
//............................................................................
void GenerateEventsUsingFluxOrTgtMix(void)
{
  //-- get flux and geom drivers
  GFluxI *        flux_driver = FluxDriver();
  GeomAnalyzerI * geom_driver = GeomDriver();

  //-- create the monte carlo job driver
  GMCJDriver * mcj_driver = new GMCJDriver;
  mcj_driver->UseFluxDriver(flux_driver);
  mcj_driver->UseGeomAnalyzer(geom_driver);
  mcj_driver->Configure();
  mcj_driver->UseSplines();
  if(!gOptWeighted) 
	mcj_driver->ForceSingleProbScale();

  //-- initialize an Ntuple Writer to save GHEP records into a TTree
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);
  ntpw.Initialize();

  //-- create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //-- generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

     LOG("gevgen", pNOTICE) << " *** Generating event............ " << ievent;

     // generate a single event for neutrinos coming from the specified flux
     EventRecord * event = mcj_driver->GenerateEvent();

     LOG("gevgen", pINFO) << "Generated Event GHEP Record: " << *event;

     // add event at the output ntuple, refresh the mc job monitor & clean-up
     ntpw.AddEventRecord(ievent, event);
     mcjmonitor.Update(ievent,event);
     ievent++;
     delete event;
  }

  //-- save the generated MC events
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

    // copy all bins between emin,emax
    spectrum  = new TH1D("spectrum","neutrino flux",
        hst->GetNbinsX(),  hst->GetXaxis()->GetXmin(), hst->GetXaxis()->GetXmax());
    spectrum->SetDirectory(0);
    for(int ibin = 1; ibin <= hst->GetNbinsX(); ibin++) {
       if(hst->GetBinLowEdge(ibin) < emax &&
         hst->GetBinLowEdge(ibin) + hst->GetBinWidth(ibin) > emin) {
           spectrum->SetBinContent(ibin, hst->GetBinContent(ibin));
	    LOG("gevgen", pNOTICE) << "adding => " << ibin << ": " << hst->GetBinContent(ibin);
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
  LOG("gevgen", pNOTICE) << "Parsing command line arguments";

  // help?
  bool help = genie::utils::clap::CmdLineArgAsBool(argc,argv,'h');
  if(help) {
      PrintSyntax();
      exit(0);
  }

  // number of events:
  try {
    LOG("gevgen", pINFO) << "Reading number of events to generate";
    gOptNevents = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pINFO)
            << "Unspecified number of events to generate - Using default";
      gOptNevents = kDefOptNevents;
    }
  }

  // run number:
  try {
    LOG("gevgen", pINFO) << "Reading MC run number";
    gOptRunNu = genie::utils::clap::CmdLineArgAsInt(argc,argv,'r');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pINFO) << "Unspecified run number - Using default";
      gOptRunNu = kDefOptRunNu;
    }
  }

  // flux functional form
  bool using_flux = false;
  try {
    LOG("gevgen", pINFO) << "Reading flux function";
    gOptFlux = genie::utils::clap::CmdLineArgAsString(argc,argv,'f');
    using_flux = true;
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
    }
  }

  // spline building option
  gOptBuildSplines = genie::utils::clap::CmdLineArgAsBool(argc,argv,'s');

  // generate weighted events option (only relevant if using a flux)
  gOptWeighted = genie::utils::clap::CmdLineArgAsBool(argc,argv,'w');

  //-- Required arguments

  // neutrino energy:
  try {
    LOG("gevgen", pINFO) << "Reading neutrino energy";
    string nue = genie::utils::clap::CmdLineArgAsString(argc,argv,'e');

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
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pFATAL) << "Unspecified neutrino energy - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // neutrino PDG code:
  try {
    LOG("gevgen", pINFO) << "Reading neutrino PDG code";
    gOptNuPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'p');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pFATAL) << "Unspecified neutrino PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // target mix (their PDG codes with their corresponding weights):
  bool using_tgtmix = false;
  try {
    LOG("gevgen", pINFO) << "Reading target mix";
    string stgtmix = genie::utils::clap::CmdLineArgAsString(argc,argv,'t');
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

  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pFATAL) << "Unspecified target PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  gOptUsingFluxOrTgtMix = using_flux || using_tgtmix;

  // print-out the command line options
  //
  LOG("gevgen", pINFO) << "MC Run Number              = " << gOptRunNu;
  LOG("gevgen", pINFO) << "Number of events requested = " << gOptNevents;
  LOG("gevgen", pINFO) << "Building splines at init.  = " << gOptBuildSplines; 
  LOG("gevgen", pINFO) << "Flux                       = " << gOptFlux;
  LOG("gevgen", pINFO) << "Generate weighted events?  = " << gOptWeighted;
  LOG("gevgen", pINFO) << "Neutrino PDG code          = " << gOptNuPdgCode;
  if(gOptNuEnergyRange>0) {
     LOG("gevgen", pINFO) << "Neutrino energy            = [" 
        << gOptNuEnergy << ", " << gOptNuEnergy+gOptNuEnergyRange << "]";
  } else {
     LOG("gevgen", pINFO) << "Neutrino energy            = " << gOptNuEnergy;
  }
  map<int,double>::const_iterator iter;
  for(iter = gOptTgtMix.begin(); iter != gOptTgtMix.end(); ++iter) {
      int    tgtpdgc = iter->first;
      double wgt     = iter->second;
      LOG("gevgen", pINFO) 
          << "Target mix - element = " <<  tgtpdgc << ", wgt = " << wgt;
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gevgen [-h] -n nev [-s] -e energy -p nupdg -t tgtmix [-f flux] [-w]\n";
}
//____________________________________________________________________________
