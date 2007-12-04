//____________________________________________________________________________
/*!

\program gevgen

\brief   GENIE v+A event generation driver 

         Syntax :
           gevgen [-n nev] [-s] -e E -p nupdg -t tgtpdg [-r run#] [-f flux] [-u]

         Options :
           [] Denotes an optional argument
           -n Specifies the number of events to generate
           -s Turns on cross section spline generation at job initialization
              Not using the cross section splines will slow down the event
              generation considerably. Since building the splines is a 
              considerable overhead we recommend building them in advance
              (see gmkspl utility) and loading via the $GSPLOAD env. variable.
           -e Specifies the neutrino energy
	      If what follows the -e option is a comma separated pair of values
              it will be interpreted as an energy range and neutrino energies 
              will be generated uniformly in this range. Note that this is the 
              energy of an "interacting" neutrino, not the energy of a "flux" 
              neutrino (the neutrino is forced to interact anyway)
              However, if a flux has been determined using the -f option then 
              the energy range will be interpreted as the energy range of flux
              neutrinos
           -p Specifies the neutrino PDG code
           -t Specifies the target PDG code (pdg format: 10LZZZAAAI)
           -r Specifies the MC run number
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
           -u Forces generation of unweighted events.
              This option is relevant only if a neutrino flux is specified.
              Note that 'unweighted' refers to the selection of the primary
              flux neutrino + target that were forced to interact. A weighting
              scheme for the generated kinematics of individual processes can
              still be in effect if enabled..
 
         Examples:
         (1)
           gevgen -n 300 -s -e 6.5 -p 14 -t 1000260560 

           will generate 300 events of muon neutrinos (pdg = 14) on Iron (A=56,Z=26) 
           at E = 6.5 GeV. The -s option will force it to generate cross section
           splines during the job initialization. Only the cross section splines
           not loaded via $GSPLOAD (see below) will be generated. Note that building
           cross section splines for all necessary processes is a time consuming step
           and we do recommend building them in advance (using the gmkspl utility)
           and loading them via $GSPLOAD.

         (2)
           gevgen -n 300 -s -e 1,5 -p 14 -t 1000260560 

           like above except that interactions will be generated uniformly at the
           energy range from 1 to 5 GeV (note: uniform spectrum of interacting 
           neutrinos _not_ flux neutrinos).

          (3)
           gevgen -n 30000 -s -e 1,4 -p 14 -t 1000080160 -f 'x*x*exp((-x*x+3)/4)'

           will generate 30000 events of muon neutrinos (pdg = 14) on Oxygen
           (A=16,Z=8) using a neutrino flux of the specified functional form at the 
           energy range from 1 to 4 GeV. Above comments on cross section splines apply
           here too.

          (4)
           gevgen -n 30000 -s -e 1,4 -p 14 -t 1000080160 -f /some/path/flux.data

           like above except that the neutrino flux is taken from the input vector file 
           rather than the input functional form.

          (4)
           gevgen -n 30000 -s -e 1,4 -p 14 -t 1000080160 -f /some/path/my_file.root,my_flux

           like above except that the flux is taken from a TH1D object called my_flux
           stored in /some/path/my_file.root. Note that the event generation driver
           will use only the input histogram bins that fall within the specified
           (with -e) energy range. In the above example, all the neutrino flux bins 
           that do not fall in the 1 GeV - 4 GeV range will be neglected (Note that
           the bins including 1 GeV and 4 GeV will be taken into account - So the 
           actual energy range used is: from the lower edge of the bin containing 
           1 GeV to the upper edge of the bin containing 4 GeV).
           
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

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <sstream>
#include <string>
#include <vector>

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
#define __CAN_GENERATE_EVENTS_USING_A_FLUX__
#endif
#endif

#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__
#include "FluxDrivers/GCylindTH1Flux.h"
#include "Geo/PointGeomAnalyzer.h"
#endif

using std::string;
using std::vector;
using std::ostringstream;

using namespace genie;
using namespace genie::controls;

void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

void GenerateEventsUsingANeutrinoFlux(void);
void GenerateEventsAtFixedEnergies   (void);

//Default options (override them using the command line arguments):
int           kDefOptNevents   = 100;     // n-events to generate
NtpMCFormat_t kDefOptNtpFormat = kNFGHEP; // ntuple format
Long_t        kDefOptRunNu     = 0;       // default run number

//User-specified options:
int           gOptNevents;      // n-events to generate
bool          gOptBuildSplines; // spline building option
double        gOptMinNuEnergy;  // min neutrino energy 
double        gOptNuEnergyRange;// max-min neutrino energy 
int           gOptNuPdgCode;    // neutrino PDG code
int           gOptTgtPdgCode;   // target PDG code
Long_t        gOptRunNu;        // run number
string        gOptFlux;         //
bool          gOptUnweighted;   // 
bool          gOptUsingFlux = false;
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
  if(gOptUsingFlux) GenerateEventsUsingANeutrinoFlux();
  else              GenerateEventsAtFixedEnergies();

  return 0;
}
//____________________________________________________________________________
void GenerateEventsUsingANeutrinoFlux(void)
{
#ifdef __CAN_GENERATE_EVENTS_USING_A_FLUX__

  int flux_entries = 100000;

  //-- create the flux driver
  double emin = gOptMinNuEnergy;
  double emax = gOptMinNuEnergy+gOptNuEnergyRange;
  double de   = gOptNuEnergyRange;

  flux::GCylindTH1Flux * flux = new flux::GCylindTH1Flux;
  TH1D * spectrum = 0;

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

  //-- create a trivial point geometry
  GeomAnalyzerI * geom_driver = new geometry::PointGeomAnalyzer(gOptTgtPdgCode);

  //-- create the monte carlo job driver
  GMCJDriver * mcj_driver = new GMCJDriver;
  mcj_driver->UseFluxDriver(flux_driver);
  mcj_driver->UseGeomAnalyzer(geom_driver);
  mcj_driver->Configure();
  mcj_driver->UseSplines();

  if(gOptUnweighted) mcj_driver->ForceSingleProbScale();

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

#else

  LOG("gevgen", pERROR) 
    << "\n   To be able to generate neutrino events from a flux you need" 
    << "\n   to add the following configuration options at your GENIE installation:" 
    << "\n   --enable-flux-drivers  --enable-geom-drivers \n" ;

#endif
}
//____________________________________________________________________________
void GenerateEventsAtFixedEnergies(void)
{
  //-- create the GENIE high level event generation object for the given 
  //   initial state
  InitialState init_state(gOptTgtPdgCode, gOptNuPdgCode);

  GEVGDriver evg_driver;
  evg_driver.Configure(init_state);

  //-- If splines are used, then create any spline that is needed but is not
  //   already loaded (can built them all here if no spline at all is loaded)
  if(gOptBuildSplines) evg_driver.CreateSplines();

  //-- initialize an Ntuple Writer
  NtpWriter ntpw(kDefOptNtpFormat, gOptRunNu);
  ntpw.Initialize();

  //-- create an MC Job Monitor
  GMCJMonitor mcjmonitor(gOptRunNu);

  //-- get a random number generator
  RandomGen * r = RandomGen::Instance();

  //-- generate events / print the GHEP record / add it to the ntuple
  int ievent = 0;
  while ( ievent < gOptNevents) {

     LOG("gevgen", pNOTICE) << " *** Generating event............ " << ievent;

     // generate neutrino energy (if an energy range was defined)
     double Ev = 0;
     if(gOptNuEnergyRange>0) {
       Ev = gOptMinNuEnergy + gOptNuEnergyRange * r->RndGen().Rndm();
     } else {
       Ev = gOptMinNuEnergy;
     }
     TLorentzVector nu_p4(0.,0.,Ev,Ev); // px,py,pz,E (GeV)

     // generate a single event
     EventRecord * event = evg_driver.GenerateEvent(nu_p4);

     LOG("gevgen", pINFO) << "Generated Event GHEP Record: " << *event;

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
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gevgen", pNOTICE) << "Parsing command line arguments";

  //-- Optional arguments

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
  try {
    LOG("gevgen", pINFO) << "Reading flux function";
    gOptFlux = genie::utils::clap::CmdLineArgAsString(argc,argv,'f');
    gOptUsingFlux = true;
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
    }
  }

  //spline building option
  gOptBuildSplines = genie::utils::clap::CmdLineArgAsBool(argc,argv,'s');

  //generate unweighted option
  gOptUnweighted = genie::utils::clap::CmdLineArgAsBool(argc,argv,'u');

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
       gOptMinNuEnergy   = emin;
       gOptNuEnergyRange = emax-emin;
    } else {
       gOptMinNuEnergy   = atof(nue.c_str());
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

  // target PDG code:
  try {
    LOG("gevgen", pINFO) << "Reading target PDG code";
    gOptTgtPdgCode = genie::utils::clap::CmdLineArgAsInt(argc,argv,'t');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gevgen", pFATAL) << "Unspecified target PDG code - Exiting";
      PrintSyntax();
      exit(1);
    }
  }

  // print the command line options
  LOG("gevgen", pINFO) << "Number of events requested = " << gOptNevents;
  LOG("gevgen", pINFO) << "Building splines at init.  = " << gOptBuildSplines; 
  LOG("gevgen", pINFO) << "Neutrino PDG code          = " << gOptNuPdgCode;
  LOG("gevgen", pINFO) << "Target PDG code            = " << gOptTgtPdgCode;
  LOG("gevgen", pINFO) << "MC Run Number              = " << gOptRunNu;
  if(gOptNuEnergyRange>0) {
    LOG("gevgen", pINFO) << "Neutrino energy            = [" 
        << gOptMinNuEnergy << ", " << gOptMinNuEnergy+gOptNuEnergyRange << "]";
  } else {
    LOG("gevgen", pINFO) << "Neutrino energy            = " << gOptMinNuEnergy;
  }
  if(gOptUsingFlux) {
    LOG("gevgen", pINFO) << "Flux                       = " << gOptFlux;
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gevgen", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gevgen [-n nev] [-s] -e energy -p nupdg -t tgtpdg [-f flux] [-u]\n";
}
//____________________________________________________________________________
