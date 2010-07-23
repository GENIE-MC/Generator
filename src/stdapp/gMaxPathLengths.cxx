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

         Options :
           -f  a ROOT file containing a ROOT/GEANT geometry description
           -L  geometry length units       [ default: mm ]
           -D  geometry density units      [ default: gr/cm3 ]
           -t  top volume name             [ default: "" ]
           -n  n scanning points / surface [ default: see geom driver's defaults ]
           -r  n scanning rays / point     [ default: see geom driver's defaults ]
           -o  output XML filename         [ default: maxpl.xml ]

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
         STFC, Rutherford Appleton Laboratory

\created September 27, 2005

\cpright Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
         For the full text of the license visit http://copyright.genie-mc.org
         or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#include <cassert>
#include <string>

#include <TMath.h>

#include "EVGDrivers/PathLengthList.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "Messenger/Messenger.h"
#include "Utils/CmdLnArgParser.h"
#include "Utils/UnitUtils.h"

using std::string;

using namespace genie;
using namespace genie::geometry;

//Prototypes:
void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

//Defaults for optional options:
string kDefOptXMLFilename  = "maxpl.xml"; // default output xml filename
string kDefOptGeomLUnits   = "mm";        // default geometry length units
string kDefOptGeomDUnits   = "g_cm3";     // default geometry density units

//User-specified options:
string gOptGeomFilename    = "";          // input geometry file
string gOptXMLFilename     = "";          // input xml filename
string gOptRootGeomTopVol  = "";          // input root geometry top vol name
double gOptGeomLUnits      = 0;           // input geometry length units
double gOptGeomDUnits      = 0;           // input geometry density units
int    gOptNPoints         = -1;          // input number of points / surf
int    gOptNRays           = -1;          // input number of rays / point

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);

  //-- create the geometry driver
  LOG("gmxpl", pINFO)
     << "Creating/configuring a ROOT geom. driver";

  ROOTGeomAnalyzer * geom = new ROOTGeomAnalyzer(gOptGeomFilename);
  geom -> SetLengthUnits       (gOptGeomLUnits);
  geom -> SetDensityUnits      (gOptGeomDUnits);
  geom -> SetWeightWithDensity (true);

  //-- set the top volume name
  geom -> SetTopVolName        (gOptRootGeomTopVol);
  geom -> SetWeightWithDensity (true);

  if(gOptNPoints > 0) geom->SetScannerNPoints(gOptNPoints);
  if(gOptNRays   > 0) geom->SetScannerNRays  (gOptNRays);

  //-- compute the maximum path lengths
  LOG("gmxpl", pINFO)
            << "Asking input GeomAnalyzerI for the max path-lengths";
  const PathLengthList & plmax = geom->ComputeMaxPathLengths();

  //-- print & save the maximum path lengths in XML format
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

  CmdLnArgParser parser(argc,argv);

  // output XML file name:
  if( parser.OptionExists('o') ) {
    LOG("gmxpl", pDEBUG) << "Reading output filename";
    gOptXMLFilename = parser.ArgAsString('o');
  } else {
    LOG("gmxpl", pDEBUG) 
       << "Unspecified output filename - Using default";
    gOptXMLFilename = kDefOptXMLFilename;
  } // -o

  // legth & density units:
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

  // root geometry top volume name:
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

  // input geometry file:
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

  // print the command line arguments
  LOG("gmxpl", pINFO) << "Command line arguments:";
  LOG("gmxpl", pINFO) << "Input ROOT geometry     = " << gOptGeomFilename;
  LOG("gmxpl", pINFO) << "Output XML file         = " << gOptXMLFilename;
  LOG("gmxpl", pINFO) << "Geometry length units   = " << gOptGeomLUnits;
  LOG("gmxpl", pINFO) << "Geometry density units  = " << gOptGeomDUnits;
  LOG("gmxpl", pINFO) << "Scanner points/surface  = " << gOptNPoints;
  LOG("gmxpl", pINFO) << "Scanner rays/point      = " << gOptNRays;
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
      << " [-o output_xml_file]";
}
//____________________________________________________________________________
