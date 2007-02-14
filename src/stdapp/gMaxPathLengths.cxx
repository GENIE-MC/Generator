//____________________________________________________________________________
/*!

\program gmxpl

\brief   GENIE utility program computing the maximum path lengths for any
         given ROOT/GEANT geometry and saving them in an output XML file.

         The maximum path lengths XML file can then be input to the GENIE
         MC Job driver to speed up the job initialization.
         In a MC job, the max, path lengths are used in computing the max.
         interaction probabilities that, at event generation, scale the
         computed interaction probabilities. The MC job driver can compute
         these maximum path lengths at initialization but since this is a
         time consuming operation for complex geometries (involving the
         generation of random points at the surface enclosing the detector,
         the generation of random rays crossing the detector and geometry
         navigation along these rays) and you only really need to do it
         just once for any particular detector setup you would be much
         better off feeding these numbers to the MC job driver.

         Note that what [path length], really means is
           for pure materials
              -> [path length]*[material density]
           for the ith element of a mixture
              -> [path length]*[mixture density]*[element weight fraction]

         Syntax :
           gmxpl -f geom_file [-u units] [-o output_xml_file] [-n np] [-r nr]

         Options :
           -f  a ROOT file containing a ROOT/GEANT geometry description
           -u  geometry length units       [ default: meter ]
           -n  n scanning points / surface [ default: see geom driver's defaults ]
           -r  n scanning rays / point     [ default: see geom driver's defaults ]
           -o  output XML filename         [ default: maxpl.xml ]

         Examples:

           gmxpl -f ~/data/mygeometry.root -u cm

           will compute the maximum path lengths for all the materials found
           in the input ROOT geometry (at mygeometry.root), in which the used
           length unit is cm and will save them in the maxpl.xml XML file.
           See $GENIE/src/Conventions/Units.h for GENIE unit definitions.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created September 27, 2005

\cpright Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
         All rights reserved.
         For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#include <cassert>
#include <string>

#include "EVGDrivers/PathLengthList.h"
#include "Geo/ROOTGeomAnalyzer.h"
#include "Messenger/Messenger.h"
#include "Utils/CmdLineArgParserUtils.h"
#include "Utils/CmdLineArgParserException.h"
#include "Utils/UnitUtils.h"

using std::string;

using namespace genie;
using namespace genie::geometry;

//Prototypes:
void GetCommandLineArgs (int argc, char ** argv);
void PrintSyntax        (void);

//Defaults for optional options:
string kDefOptXMLFilename = "maxpl.xml";
string kDefOptGeomUnits   = "meter";

//User-specified options:
string gOptGeomFilename = "";
string gOptXMLFilename  = "";
string gOptGeomUnits    = "";
int    gOptNPoints      = -1;
int    gOptNRays        = -1;

//____________________________________________________________________________
int main(int argc, char ** argv)
{
  //-- parse command line arguments
  GetCommandLineArgs(argc,argv);

  //-- print the options you got from command line arguments
  LOG("gmxpl", pINFO) << "Command line arguments:";
  LOG("gmxpl", pINFO) << "Input ROOT geometry    = " << gOptGeomFilename;
  LOG("gmxpl", pINFO) << "Output XML file        = " << gOptXMLFilename;
  LOG("gmxpl", pINFO) << "Geometry units         = " << gOptGeomUnits;
  LOG("gmxpl", pINFO) << "Scanner points/surface = " << gOptNPoints;
  LOG("gmxpl", pINFO) << "Scanner rays/point     = " << gOptNRays;

  //-- create the geometry driver
  LOG("gmxpl", pINFO)
            << "Creating/configuring a ROOT geom. driver";

  ROOTGeomAnalyzer * geom = new ROOTGeomAnalyzer(gOptGeomFilename);

  geom->SetUnits(genie::utils::units::UnitFromString(gOptGeomUnits));
  geom->SetWeightWithDensity(true);

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
  LOG("gmxpl", pINFO) << "Parsing commad line arguments";

  //-- Optional arguments

  //output XML file name:
  try {
    LOG("gmxpl", pINFO) << "Reading output filename";
    gOptXMLFilename = genie::utils::clap::CmdLineArgAsString(argc,argv,'o');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmxpl", pINFO) << "Unspecified output filename - Using default";
      gOptXMLFilename = kDefOptXMLFilename;
    }
  }
  //geometry units:
  try {
    LOG("gmxpl", pINFO) << "Reading input geometry units";
    gOptGeomUnits = genie::utils::clap::CmdLineArgAsString(argc,argv,'u');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmxpl", pINFO) << "Unspecified geometry units - Using default";
      gOptGeomUnits = kDefOptGeomUnits;
    }
  }

  //number of scanning points / surface
  try {
    LOG("gmxpl", pINFO) << "Reading input number of scanning points/surface";
    gOptNPoints = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmxpl", pINFO)
         << "Unspecified number of points - Using driver's default";
    }
  }
  //number of scanning rays / point
  try {
    LOG("gmxpl", pINFO) << "Reading input number of scanning rays/point";
    gOptNRays = genie::utils::clap::CmdLineArgAsInt(argc,argv,'n');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmxpl", pINFO)
         << "Unspecified number of rays - Using driver's default";
    }
  }

  //-- Required arguments

  //output geometry file:
  try {
    LOG("gmxpl", pINFO) << "Reading ROOT/GEANT geometry filename";
    gOptGeomFilename = genie::utils::clap::CmdLineArgAsString(argc,argv,'f');
  } catch(exceptions::CmdLineArgParserException e) {
    if(!e.ArgumentFound()) {
      LOG("gmxpl", pFATAL) << "No geometry file was specified - Exiting";
      PrintSyntax();
      exit(1);
    }
  }
}
//____________________________________________________________________________
void PrintSyntax(void)
{
  LOG("gmxpl", pNOTICE)
      << "\n\n" << "Syntax:" << "\n"
          << "   gmxpl -f geom_file [-u geom_units] [-o output_xml_file]";
}
//____________________________________________________________________________
