//____________________________________________________________________________
/*!

\program testROOTGeometry

\brief   Tests the ROOT Geometry Analyzer.

\syntax  testROOTGeometry -f /path/to/root/geometry/file.root

         If no file is specified, $GENIE/src/test/TestGeometry.root is used
         as default.

\author  Anselmo Meregaglia <anselmo.meregaglia@cern.ch>
         ETH Zurich

\created August 11, 2005
*/
//____________________________________________________________________________

#include <string>

#include <TLorentzVector.h>

#include "Geo/ROOTGeomAnalyzer.h"
#include "Geo/PathLengthList.h"
#include "Messenger/Messenger.h"

using std::string;

using namespace genie;

int main(int argc, char ** argv)
{
  // get filename from the command line argument (following -f)

  string filename = "$GENIE/src/test/TestGeometry.root"; // default
  for(int iarg = 0; iarg < argc-1; iarg++) {
     string argument(argv[iarg]);
     if( argument.compare("-f") == 0 ) filename = string(argv[++iarg]);
  }

  LOG("Test",pINFO) << "Starting ROOTGeomAnalyzer with geometry from: " << filename;

  ROOTGeomAnalyzer* root_analyzer = new ROOTGeomAnalyzer(filename);
 
  LOG("Test",pINFO) << "Computing path lengths";

  TLorentzVector* x= new TLorentzVector(0,0,0,0);
  TLorentzVector* p= new TLorentzVector(1,0,0,1);
  
  const PathLengthList & pl = root_analyzer->ComputePathLengths(*x,*p);

  LOG("Test",pINFO) << "Printing computed path lengths:";
  LOG("Test",pINFO) << pl;
  
  double matLength = root_analyzer->SetVtxMaterial(1039018000);

  LOG("Test",pINFO) << "Path in selected material ...";
  LOG("Test",pINFO) << matLength;
  return 0;
}

