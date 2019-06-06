//____________________________________________________________________________
/*!

\program gmxpl

\brief   GENIE utility program computing the maximum path lengths for any
         given ROOT/GEANT geometry and saving them in an output XML file.

         The maximum path lengths XML file can then be input to the GENIE
         event generation drivers to speed up the job initialization.

         Note that this program actually computes the 'density weighted' path
         lengths required for computing interaction probabilities in the input
         geometry volumes.
         For pure materials, this program computes:
              -> [path length]*[material density]
         whereas,  for the ith element of a mixture, it computes:
              -> [path length]*[mixture density]*[element weight fraction]

         Syntax :
           gmxpl -f geom_file [-L length_units] [-D density_units] 
                 [-t top_vol_name] [-o output_xml_file] [-n np] [-r nr]
                 [-seed random_number_seed]
                 [--message-thresholds xml_file]

         Options :
           -f  
              A ROOT file containing a ROOT/GEANT geometry description
           -L  
               Geometry length units [ default: mm ]
           -D  
               Geometry density units [ default: gr/cm3 ]
           -t  
               Top volume name [ default: "" ]
           -n  
               Number of  scanning points / surface [ default: see geom driver's defaults ]
           -r  
               Number of scanning rays / point [ default: see geom driver's defaults ]
           -o  
               Name of output XML file [ default: maxpl.xml ]
           --seed 
               Random number seed.
          --message-thresholds
              Allows users to customize the message stream thresholds.
              The thresholds are specified using an XML file.
              See $GENIE/config/Messenger.xml for the XML schema.

         Example:

           gmxpl -f mygeometry.root -L m -D kg_m3 -o out.xml -n 1000 -r 1000

           will compute the maximum density weighted path lengths for all the 
           materials of the ROOT geometry at the mygeometry.root input file. 
           The program will use 'm' and 'kg/m^3' as the length and density 
           units of the input geometry. 
           The input geometry will be scanned with 1E+3 points / surface and
           1E+3 rays / surface.
           Results will be saved in the out.xml XML file at SI units.
           See $GENIE/src/Conventions/Units.h for GENIE unit definitions.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab

\created September 27, 2005

\cpright Copyright (c) 2003-2019, The GENIE Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <string>

#include <TMath.h>

#include "Framework/EventGen/PathLengthList.h"
#include "Tools/Geometry/ROOTGeomAnalyzer.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Utils/CmdLnArgParser.h"
#include "Framework/Utils/UnitUtils.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Utils/AppInit.h"
#include "Framework/Utils/RunOpt.h"

using std::string;

using namespace genie;
using namespace genie::geometry;

// Prototypes:
void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

// Defaults for optional options:
string kDefOptXMLFilename  = "maxpl.xml"; // default output xml filename
string kDefOptGeomLUnits   = "mm";        // default geometry length units
string kDefOptGeomDUnits   = "g_cm3";     // default geometry density units

// User-specified options:
string    gOptGeomFilename    = "";          // input geometry file
string    gOptXMLFilename     = "";          // input xml filename
string    gOptRootGeomTopVol  = "";          // input root geometry top vol name
double    gOptGeomLUnits      = 0;           // input geometry length units
double    gOptGeomDUnits      = 0;           // input geometry density units
int       gOptNPoints         = -1;          // input number of points / surf
int       gOptNRays           = -1;          // input number of rays / point
long int  gOptRanSeed         = -1;          // random number seed

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  GetCommandLineArgs(argc,argv);

  utils::app_init::RandGen(gOptRanSeed);
  utils::app_init::MesgThresholds(RunOpt::Instance()->MesgThresholdFiles());

  // Create the geometry driver
  LOG("gmxpl", pINFO)
     << "Creating/configuring a ROOT geom. driver";

  ROOTGeomAnalyzer * geom = new ROOTGeomAnalyzer(gOptGeomFilename);
  geom -> SetLengthUnits       (gOptGeomLUnits);
  geom -> SetDensityUnits      (gOptGeomDUnits);
  geom -> SetWeightWithDensity (true);

  // Set the top volume name
  geom -> SetTopVolName        (gOptRootGeomTopVol);
  geom -> SetWeightWithDensity (true);

  if(gOptNPoints > 0) geom->SetScannerNPoints(gOptNPoints);
  if(gOptNRays   > 0) geom->SetScannerNRays  (gOptNRays);

  // Compute the maximum path lengths
  LOG("gmxpl", pINFO)
      << "Asking input GeomAnalyzerI for the max path-lengths";
  const PathLengthList & plmax = geom->ComputeMaxPathLengths();

  // Print & save the maximum path lengths in XML format
  LOG("gmxpl", pINFO)
      << "Maximum path lengths: " << plmax;
  plmax.SaveAsXml(gOptXMLFilename);

  delete geom;

  return 0;
}
//____________________________________________________________________________
void GetCommandLineArgs(int argc, char ** argv)
{
  LOG("gmxpl", pINFO) << "Parsing command line arguments";

  // Common run options. 
  RunOpt::Instance()->ReadFromCommandLine(argc,argv);

  // Parse run options for this app

  CmdLnArgParser parser(argc,argv);

  // output XML file name
  if( parser.OptionExists('o') ) {
    LOG("gmxpl", pDEBUG) << "Reading output filename";
    gOptXMLFilename = parser.ArgAsString('o');
  } else {
    LOG("gmxpl", pDEBUG) 
       << "Unspecified output filename - Using default";
    gOptXMLFilename = kDefOptXMLFilename;
  } // -o

  // legth & density units
  string lunits, dunits;
  if( parser.OptionExists('L') ) {
    LOG("gmxpl", pDEBUG) << "Checking for input geometry length units";
    lunits = parser.ArgAsString('L');
  } else {
    LOG("gmxpl", pDEBUG) << "Using default geometry length units";
    lunits = kDefOptGeomLUnits;
  } // -L
  if( parser.OptionExists('D') ) {
    LOG("gmxpl", pDEBUG) << "Checking for input geometry density units";
    dunits = parser.ArgAsString('D');
  } else {
    LOG("gmxpl", pDEBUG) << "Using default geometry density units";
    dunits = kDefOptGeomDUnits;
  } // -D
  gOptGeomLUnits = genie::utils::units::UnitFromString(lunits);
  gOptGeomDUnits = genie::utils::units::UnitFromString(dunits);

  // root geometry top volume name
  if( parser.OptionExists('t') ) {
    LOG("gmxpl", pDEBUG) 
       << "Reading root geometry top volume name";
    gOptRootGeomTopVol = parser.ArgAsString('t');
  } else {
    LOG("gmxpl", pDEBUG) 
       << "Unspecified geometry top volume - Using default";
    gOptRootGeomTopVol = "";
  } // -o
  
  // number of scanning points / surface
  if( parser.OptionExists('n') ) {
    LOG("gmxpl", pDEBUG) 
       << "Reading input number of scanning points/surface";
    gOptNPoints = parser.ArgAsInt('n');
  } else {
    LOG("gmxpl", pDEBUG)
      << "Unspecified number of points - Using driver's default";
  } //-n

  // number of scanning rays / point
  if( parser.OptionExists('r') ) {
    LOG("gmxpl", pDEBUG) 
       << "Reading input number of scanning rays/point";
    gOptNRays = parser.ArgAsInt('r');
  } else {
    LOG("gmxpl", pDEBUG)
      << "Unspecified number of rays - Using driver's default";
  } //-r

  // input geometry file
  if( parser.OptionExists('f') ) {
    LOG("gmxpl", pDEBUG) 
       << "Reading ROOT/GEANT geometry filename";
    gOptGeomFilename = parser.ArgAsString('f');
  } else {
    LOG("gmxpl", pFATAL) 
       << "No geometry file was specified - Exiting";
    PrintSyntax();
    exit(1);
  } //-f

  // random number seed
  if( parser.OptionExists("seed") ) {
    LOG("gmxpl", pINFO) << "Reading random number seed";
    gOptRanSeed = parser.ArgAsLong("seed");
  } else {
    LOG("gmxpl", pINFO) << "Unspecified random number seed - Using default";
    gOptRanSeed = -1;
  }

  // print the command line arguments
  LOG("gmxpl", pNOTICE)
     << "\n"
     << utils::print::PrintFramedMesg("gmxpl job inputs");
  LOG("gmxpl", pNOTICE) << "Command line arguments";
  LOG("gmxpl", pNOTICE) << "Input ROOT geometry     : " << gOptGeomFilename;
  LOG("gmxpl", pNOTICE) << "Output XML file         : " << gOptXMLFilename;
  LOG("gmxpl", pNOTICE) << "Geometry length units   : " << gOptGeomLUnits;
  LOG("gmxpl", pNOTICE) << "Geometry density units  : " << gOptGeomDUnits;
  LOG("gmxpl", pNOTICE) << "Scanner points/surface  : " << gOptNPoints;
  LOG("gmxpl", pNOTICE) << "Scanner rays/point      : " << gOptNRays;
  LOG("gmxpl", pNOTICE) << "Random number seed      : " << gOptRanSeed;

  LOG("gmxpl", pNOTICE) << "\n";
  LOG("gmxpl", pNOTICE) << *RunOpt::Instance();
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gmxpl", pNOTICE)
      << "\n\n" << "Syntax:" << "\n"
      << "   gmxpl"
      << " -f geom_file"
      << " [-L length_units]"
      << " [-D density_units]" 
      << " [-t top_volume_name]"
      << " [-o output_xml_file]"
      << " [-seed random_number_seed]"
      << " [--message-thresholds xml_file]\n";

}
//____________________________________________________________________________
