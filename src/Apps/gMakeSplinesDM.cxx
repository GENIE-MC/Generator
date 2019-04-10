//____________________________________________________________________________
/*!

\program gmkspl_dm

\brief   GENIE utility program building XML cross section splines that can
         be loaded into GENIE to speed-up event generation.
         The list of neutrino PDG codes is passed from the command line.
         The list of nuclear target PDG codes is either passed from the
         command line or extracted from the input ROOT/GEANT geometry.

         Syntax :
           gmkspl_dm -m masses
                 <-t target_pdg_codes,
                  -f geometry_file>
                  <-o | --output-cross-sections> output_xml_xsec_file
                  [-g zp_couplings]
                  [-z med_ratios]
                  [-n nknots]
                  [-e max_energy]
                  [--no-copy]
                  [--seed random_number_seed]
                  [--input-cross-sections xml_file]
                  [--event-generator-list list_name]
                  [--tune genie_tune]
                  [--message-thresholds xml_file]

         Note :
           [] marks optional arguments.
           <> marks a list of arguments out of which only one can be
              selected at any given time.

         Options :
           -m
               A comma separated list of DM masses.
           -t
               A comma separated list of tgt PDG codes.
               PDG code format: 10LZZZAAAI
           -f
               A ROOT file containing a ROOT/GEANT geometry description.
           -o, --output-cross-sections
               Name of output XML file containing computed cross-section data.
               Default: `xsec_splines.xml'.
           -g
               A comma separated list of Z' coupling constants
               Default: Value in UserPhysicsOptions.xml
           -z
               A comma separated list of ratios of mediator mass to dark
               matter mass
               Default: 0.5
           -n
               Number of knots per spline.
               Default: 15 knots per decade of energy range with a minimum
               of 30 knots totally.
           -e
               Maximum energy in spline.
               Default: The max energy in the validity range of the spline
               generating thread.
           --no-copy
               Does not write out the input cross-sections in the output file
           --seed
              Random number seed.
           --input-cross-sections
              Name (incl. full path) of an XML file with pre-computed
              free-nucleon cross-section values. If loaded, it can speed-up
              cross-section calculation for nuclear targets.
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

#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>

#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
#include <fenv.h> // for `feenableexcept`
#endif

#include <TSystem.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/EventGen/GEVGDriver.h"
#include "Framework/Interaction/Interaction.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/RunOpt.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/StringUtils.h"
//#include "Framework/Utils/SystemUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/XSecSplineList.h"
#include "Framework/Utils/CmdLnArgParser.h"

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
#include "Tools/Geometry/ROOTGeomAnalyzer.h"
#endif

using std::string;
using std::vector;

using namespace genie;

#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
using namespace genie::geometry;
#endif

// Prototypes:
void          GetCommandLineArgs (int argc, char ** argv);
void          PrintSyntax        (void);
PDGCodeList * GetTargetCodes     (void);

// User-specified options:
string   gOptTgtPdgCodeList = "";
string   gOptGeomFilename   = "";
int      gOptNKnots         = -1;
double   gOptMaxE           = -1.;
vector<double>   gOptDMMasses;
vector<double>   gOptMedRatios;
vector<double>   gOptZpCouplings;
bool     gOptNoCopy         = false;
long int gOptRanSeed        = -1;   // random number seed
string   gOptInpXSecFile    = "";   // input cross-section file
string   gOptOutXSecFile    = "";   // output cross-section file

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  // Parse command line arguments
  GetCommandLineArgs(argc,argv);

  if ( ! RunOpt::Instance()->Tune() ) {
    LOG("gmkspl", pFATAL) << " No TuneId in RunOption";
    exit(-1);
  }
  RunOpt::Instance()->BuildTune();

  for (vector<double>::iterator mass = gOptDMMasses.begin(); mass != gOptDMMasses.end(); ++mass) {
    for (vector<double>::iterator ratio = gOptMedRatios.begin(); ratio != gOptMedRatios.end(); ++ratio) {
      for (vector<double>::iterator coup = gOptZpCouplings.begin(); coup != gOptZpCouplings.end(); ++coup) {
        // Add dark matter to the table
        PDGLibrary::Instance()->ReloadDBase();
        PDGLibrary::Instance()->AddDarkMatter(*mass,*ratio);
        if (*coup > 0.) {
            Registry * r = AlgConfigPool::Instance()->CommonList("Param", "BoostedDarkMatter");
            r->UnLock();
            r->Set("ZpCoupling", *coup);
            r->Lock();
        }


        // throw on NaNs and Infs...
#if defined(HAVE_FENV_H) && defined(HAVE_FEENABLEEXCEPT)
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

        // Init
        utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());
        utils::app_init::RandGen(gOptRanSeed);
        utils::app_init::XSecTable(gOptInpXSecFile, false);

        // Get list of neutrinos and nuclear targets

        //  PDGCodeList * neutrinos = GetNeutrinoCodes();
        PDGCodeList * targets   = GetTargetCodes();

        if(!targets || targets->size() == 0 ) {
          LOG("gmkspl_dm", pFATAL) << "Empty target PDG code list";
          PrintSyntax();
          exit(3);
        }

        LOG("gmkspl_dm", pINFO) << "Targets: "   << *targets;

        // Loop over all possible input init states and ask the GEVGDriver
        // to build splines for all the interactions that its loaded list
        // of event generators can generate.

        //  PDGCodeList::const_iterator nuiter;
        PDGCodeList::const_iterator tgtiter;
        //  for(nuiter = neutrinos->begin(); nuiter != neutrinos->end(); ++nuiter) {
        for(tgtiter = targets->begin(); tgtiter != targets->end(); ++tgtiter) {
          int dmpdgc  = kPdgDarkMatter;
          int tgtpdgc = *tgtiter;
          InitialState init_state(tgtpdgc, dmpdgc);
          GEVGDriver driver;
          driver.SetEventGeneratorList(RunOpt::Instance()->EventGeneratorList());
          driver.Configure(init_state);
          driver.CreateSplines(gOptNKnots, gOptMaxE);
        }
        delete targets;
      }
    }
  }
  //  }

  // Save the splines at the requested XML file
  XSecSplineList * xspl = XSecSplineList::Instance();
  bool save_init = !gOptNoCopy;
  xspl->SaveAsXml(gOptOutXSecFile, save_init);


  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gmkspl_dm", pINFO) << "Parsing command line arguments";

  // Common run options. Set defaults and read.
  RunOpt::Instance()->EnableBareXSecPreCalc(true);
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // output XML file name
  if( parser.OptionExists('o') ||
      parser.OptionExists("output-cross-sections") )
  {
    LOG("gmkspl_dm", pINFO) << "Reading output filename";
    if( parser.OptionExists('o') ) {
      gOptOutXSecFile = parser.ArgAsString('o');
    }
    else {
      gOptOutXSecFile = parser.ArgAsString("output-cross-sections");
    }
  } else {
    LOG("gmkspl_dm", pINFO) << "Unspecified filename - Using default";
    gOptOutXSecFile = "xsec_splines.xml";
  }

  // number of knots
  if( parser.OptionExists('n') ) {
    LOG("gmkspl_dm", pINFO) << "Reading number of knots/spline";
    gOptNKnots = parser.ArgAsInt('n');
  } else {
    LOG("gmkspl_dm", pINFO)
      << "Unspecified number of knots - Using default";
    gOptNKnots = -1;
  }

  // max spline energy (if < max of validity range)
  if( parser.OptionExists('e') ) {
    LOG("gmkspl_dm", pINFO) << "Reading maximum spline energy";
    gOptMaxE = parser.ArgAsDouble('e');
  } else {
    LOG("gmkspl_dm", pINFO)
       << "Unspecified maximum spline energy - Using default";
    gOptMaxE = -1;
  }

  // get the dark matter mass
  if( parser.OptionExists('m') ) {
    LOG("gmkspl_dm", pINFO) << "Reading dark matter mass";
    string dmm = parser.ArgAsString('m');

    if(dmm.find(",") != string::npos) {
      vector<string> massrange = utils::str::Split(dmm, ",");
      assert(massrange.size() > 1);
      for (vector<string>::iterator m = massrange.begin(); m != massrange.end(); ++m) {
        gOptDMMasses.push_back(atof((*m).c_str()));
      }
    }
    else {
      gOptDMMasses.push_back(atof(dmm.c_str()));
    }
  } else {
    LOG("gmkspl_dm", pFATAL)
       << "Unspecified dark matter masses - Exiting";
    PrintSyntax();
    exit(1);
  }

  // get the mediator mass ratio
  if( parser.OptionExists('z') ) {
    LOG("gmkspl_dm", pINFO) << "Reading mediator mass ratios";
    string mr = parser.ArgAsString('z');

    if(mr.find(",") != string::npos) {
      vector<string> ratiorange = utils::str::Split(mr, ",");
      assert(ratiorange.size() > 1);
      for (vector<string>::iterator r = ratiorange.begin(); r != ratiorange.end(); ++r) {
        gOptMedRatios.push_back(atof((*r).c_str()));
      }
    }
    else {
      gOptMedRatios.push_back(atof(mr.c_str()));
    }
  } else {
    LOG("gmkspl_dm", pINFO)
       << "Unspecified mediator mass ratios - Using default";
    gOptMedRatios.push_back(0.5);
  }

  // write out input splines?
  if( parser.OptionExists("no-copy") ) {
    LOG("gmkspl_dm", pINFO) << "Not copying input splines to output";
    gOptNoCopy = true;
  }

  // get the mediator coupling
  if( parser.OptionExists('g') ) {
    LOG("gmkspl_dm", pINFO) << "Reading mediator couplings";
    string gzp = parser.ArgAsString('g');

    if(gzp.find(",") != string::npos) {
      vector<string> couprange = utils::str::Split(gzp, ",");
      assert(couprange.size() > 1);
      for (vector<string>::iterator g = couprange.begin(); g != couprange.end(); ++g) {
        gOptZpCouplings.push_back(atof((*g).c_str()));
      }
    }
    else {
      gOptZpCouplings.push_back(atof(gzp.c_str()));
    }
  } else {
    LOG("gmkspl_dm", pNOTICE)
       << "Unspecified mediator couplings - Using value from config";
    gOptZpCouplings.push_back(-1.);
  }

  // write out input splines?
  if( parser.OptionExists("no-copy") ) {
    LOG("gmkspl_dm", pINFO) << "Not copying input splines to output";
    gOptNoCopy = true;
  }


  // comma-separated target PDG code list or input geometry file
  bool tgt_cmd = true;
  if( parser.OptionExists('t') ) {
    LOG("gmkspl_dm", pINFO) << "Reading target nuclei PDG codes";
    gOptTgtPdgCodeList = parser.ArgAsString('t');
  } else {
    LOG("gmkspl_dm", pINFO) << "No code list specified from the command line";
    tgt_cmd = false;
  }

  bool tgt_geom = true;
  if( parser.OptionExists('f') ) {
    LOG("gmkspl_dm", pINFO) << "Reading ROOT geometry filename";
    gOptGeomFilename = parser.ArgAsString('f');
  } else {
    LOG("gmkspl_dm", pINFO) << "No geometry file was specified";
    tgt_cmd = false;
  }

  bool both =  tgt_geom &&  tgt_cmd;
  bool none = !tgt_geom && !tgt_cmd;
  if(none) {
    LOG("gmkspl_dm", pFATAL)
          << "No geom file or cmd line target list was specified - Exiting";
    PrintSyntax();
    exit(1);
  }
  if(both) {
    LOG("gmkspl_dm", pFATAL)
       << "You specified both a geom file and a cmd line target list "
         << "- Exiting confused";
    PrintSyntax();
    exit(1);
  }

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gmkspl_dm", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gmkspl_dm", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  // input cross-section file
  if( parser.OptionExists("input-cross-sections") ) {
    LOG("gmkspl_dm", pINFO) << "Reading cross-section file";
    gOptInpXSecFile = parser.ArgAsString("input-cross-sections");
  } else {
    LOG("gmkspl_dm", pINFO) << "Unspecified input cross-section file";
    gOptInpXSecFile = "";
  }

  //
  // print the command-line options
  //
  LOG("gmkspl_dm", pNOTICE)
     << "\n"
     << utils::print::PrintFramedMesg("gmkspl_dm job configuration")
     << "\n Target PDG codes : " << gOptTgtPdgCodeList
     << "\n Input ROOT geometry : " << gOptGeomFilename
     << "\n Output cross-section file : " << gOptOutXSecFile
     << "\n Input cross-section file : " << gOptInpXSecFile
     << "\n Random number seed : " << gOptRanSeed
     << "\n";

  // print list of DM masses
  for (vector<double>::iterator m = gOptDMMasses.begin(); m != gOptDMMasses.end(); ++m) {
    LOG("gmkspl_dm",pNOTICE)
      << "Considering DM mass : " << *m;
  }
  for (vector<double>::iterator r = gOptMedRatios.begin(); r != gOptMedRatios.end(); ++r) {
    LOG("gmkspl_dm",pNOTICE)
      << "Considering mediator mass ratio : " << *r;
  }
  for (vector<double>::iterator g = gOptZpCouplings.begin(); g != gOptZpCouplings.end(); ++g) {
    LOG("gmkspl_dm",pNOTICE)
      << "Considering Z' couplings : " << *g;
  }
  LOG("gmkspl_dm", pNOTICE) << *RunOpt::Instance();
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gmkspl_dm", pNOTICE)
    << "\n\n" << "Syntax:" << "\n"
    << "   gmkspl_dm -m dm_masses "
    << " <-t tgtpdg, -f geomfile> "
    << " <-o | --output-cross-section> xsec_xml_file_name"
    << " [-g zp_couplings] "
    << " [-z med_ratios] "
    << " [-n nknots] [-e max_energy] "
    << " [--seed seed_number]"
    << " [--input-cross-section xml_file]"
    << " [--event-generator-list list_name]"
    << " [--message-thresholds xml_file]\n\n";
}
//____________________________________________________________________________
PDGCodeList * GetTargetCodes(void)
{
  bool from_geom_file = ( gOptGeomFilename.size()   > 0 );
  bool from_cmd_line  = ( gOptTgtPdgCodeList.size() > 0 );

  if (from_cmd_line) {
     // split the comma separated list
     vector<string> tgtvec = utils::str::Split(gOptTgtPdgCodeList, ",");

     // fill in the PDG code list
     PDGCodeList * list = new PDGCodeList;
     vector<string>::const_iterator iter;
     for(iter = tgtvec.begin(); iter != tgtvec.end(); ++iter) {
        list->push_back( atoi(iter->c_str()) );
     }
     return list;
  }

  if (from_geom_file) {
#ifdef __GENIE_GEOM_DRIVERS_ENABLED__
     // create/configure a geometry driver
     LOG("gmkspl_dm", pINFO) << "Creating/configuring a ROOT geom. driver";
     ROOTGeomAnalyzer * geom = new ROOTGeomAnalyzer(gOptGeomFilename);

     PDGCodeList * list = new PDGCodeList(geom->ListOfTargetNuclei());

     delete geom;
     return list;
#else
     LOG("gmkspl_dm", pFATAL)
      << "To read-in a ROOT geometry you need to enable the geometry drivers!";
     gAbortingInErr = true;
     exit(1);
     return 0;
#endif

  }
  return 0;
}
//____________________________________________________________________________
