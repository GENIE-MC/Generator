//____________________________________________________________________________
/*!

\program testROOTGeometry

\author  Anselmo Meregaglia <anselmo.meregaglia@cern.ch>
         ETH Zurich

\created August 11, 2005
*/
//____________________________________________________________________________

#include "Geo/ROOTGeomAnalyzer.h"
#include <TLorentzVector.h>

using namespace genie;

int main(int argc, char ** argv)
{

  ROOTGeomAnalyzer* root_analyzer = new ROOTGeomAnalyzer();
  
  //root_analyzer->Load("$GENIE/src/test/TestGeometry.root");
  //root_analyzer->Load("/eth/store_6/home_6/amerega/TEST/vgm.2.03/examples/E01/N03/bin/Linux-g++/Geometry.root");

  //root_analyzer->SetVtxMaterial("Lead");

   
  //TLorentzVector* x= new TLorentzVector(0,0,0,0);
  //TLorentzVector* p= new TLorentzVector(1,0,0,1);

  //root_analyzer->ComputePathLengths(*x,*p);
  root_analyzer->test();
  
  return 0;
}
