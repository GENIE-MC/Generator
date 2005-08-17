//____________________________________________________________________________
/*!

\class   genie::ROOTGeomAnalyzer

\brief   A ROOT/GEANT Geometry Analyzer

\author  Anselmo Meregaglia <anselmo.meregaglia@cern.ch>
         ETH Zurich

\created May 24, 2005

*/
//____________________________________________________________________________

#include <TGeoVolume.h>
#include <TLorentzVector.h>
#include <TVector3.h>

#include "Geo/ROOTGeomAnalyzer.h"
#include "Geo/PathLengthList.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"

using namespace genie;

//___________________________________________________________________________
ROOTGeomAnalyzer::ROOTGeomAnalyzer(TGeoVolume * geom) :
GeomAnalyzerI(),
fGeometry(geom)
{
  this->Initialize();
}
//___________________________________________________________________________
ROOTGeomAnalyzer::~ROOTGeomAnalyzer()
{
  if( fCurrPathLengthList ) delete fCurrPathLengthList;
  if( fCurrPDGCodeList    ) delete fCurrPDGCodeList;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::Initialize(void)
{
  fCurrPathLengthList = 0;
  fCurrPDGCodeList    = 0;

  this->BuildListOfTargetNuclei();

  const PDGCodeList & pdglist = this->ListOfTargetNuclei();

  fCurrPathLengthList = new PathLengthList(pdglist);
  fCurrVertex = new TVector3(0.,0.,0.);
}
//___________________________________________________________________________
const PDGCodeList & ROOTGeomAnalyzer::ListOfTargetNuclei(void)
{
  return *fCurrPDGCodeList;
}
//___________________________________________________________________________
const PathLengthList & ROOTGeomAnalyzer::ComputePathLengths(
                          const TLorentzVector & x, const TLorentzVector & p)
{
// Computes the path-length within each detector material for a neutrino
// starting from point x and travelling along the direction of p

  // reset current list of path-lengths
  fCurrPathLengthList->SetAllToZero();

  LOG("GROOTGeom", pINFO)
       << "\nComputing path-lengths for neutrino: "
       << "\n  with 4-momentum : " << print_utils::P4AsString(&p)
       << "\n  starting from   : " << print_utils::X4AsString(&x);

  //-- fill in list using the input TGeoVolume
  //
  // ... ... ...

  return *fCurrPathLengthList;
}
//___________________________________________________________________________
const TVector3 & ROOTGeomAnalyzer::GenerateVertex(
              const TLorentzVector & x, const TLorentzVector & p, int tgtpdg)
{
// Generates a random vertex, within the detector material with the input
// PDG code, for a neutrino starting from point x and travelling along the
// direction of p

  // reset current interaction vertex
  fCurrVertex->SetXYZ(0.,0.,0.);

  LOG("GROOTGeom", pINFO)
       << "\nGenerating an vertex at material with PDG code = " << tgtpdg
       << "\nfor a neutrino: "
       << "\n  with 4-momentum : " << print_utils::P4AsString(&p)
       << "\n  starting from   : " << print_utils::X4AsString(&x);

  //-- generate the vertex
  //
  // ... ... ...

  return *fCurrVertex;
}
//___________________________________________________________________________
void ROOTGeomAnalyzer::BuildListOfTargetNuclei(void)
{
  fCurrPDGCodeList = new PDGCodeList;

  //-- fill in list using the input TGeoVolume
  //
  //... ...
}
//___________________________________________________________________________
