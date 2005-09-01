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
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <TLorentzVector.h>

#include "Geo/ROOTGeomAnalyzer.h"
#include "Geo/PathLengthList.h"
#include "Messenger/Messenger.h"

using std::string;

using namespace genie;

int main(int argc, char ** argv)
{
  //-- Default geometry
  string base_dir = string( gSystem->Getenv("GENIE") );
  string filename = base_dir+ string("/src/test/TestGeometry.root");
  //-- Scan for filename from the command line argument (following -f)
  for(int iarg = 0; iarg < argc-1; iarg++) {
     string argument(argv[iarg]);
     if( argument.compare("-f") == 0 ) filename = string(argv[++iarg]);
  }

  LOG("Test",pINFO) << "Starting ROOTGeomAnalyzer with geometry from: " << filename;
  ROOTGeomAnalyzer* root_analyzer = new ROOTGeomAnalyzer(filename);


  LOG("Test",pINFO) << "Computing Max path lengths";
  const PathLengthList & maxpl = root_analyzer->ComputeMaxPathLengths();
  LOG("Test",pINFO) << "Printing computed Max path lengths:";
  LOG("Test",pINFO) << maxpl;

  LOG("Test",pINFO) << "Computing path lengths";
  TLorentzVector* x= new TLorentzVector(0,0,0,0);
  TLorentzVector* p= new TLorentzVector(1,0,0,1);
  const PathLengthList & pl = root_analyzer->ComputePathLengths(*x,*p);
  LOG("Test",pINFO) << "Printing computed path lengths:";
  LOG("Test",pINFO) << pl;

  //material selected for the vertex generation
  int pdg(1039018000);
  //number of vertices to be generated
  int numVtx(100);
  ofstream outfileVTX("VtxCoord.txt",std::ios_base::app);

  for(int i=0;i<numVtx;i++)
    {
      const TVector3 & vtx = root_analyzer->GenerateVertex(*x,*p,pdg);
      LOG("Test",pINFO) << "Vertex selected ...";
      LOG("Test",pINFO) << " x "<<vtx.X()<<" y "<<vtx.Y()<<" z "<<vtx.Z();
      outfileVTX<<vtx.X()<<"\t"<<vtx.Y()<<"\t"<<vtx.Z()<<std::endl;
    }

  return 0;
}

