//____________________________________________________________________________
/*
 Copyright (c) 2003-2007, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - July 14, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TVector3.h>

#include "Geo/PointGeomAnalyzer.h"
#include "EVGDrivers/PathLengthList.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"

using namespace genie;
using namespace genie::geometry;

//___________________________________________________________________________
PointGeomAnalyzer::PointGeomAnalyzer(int tgtpdgc) :
GeomAnalyzerI()
{
  fCurrVertex = new TVector3(0,0,0);

  fCurrPDGCodeList = new PDGCodeList;
  fCurrPDGCodeList->clear();
  fCurrPDGCodeList->push_back(tgtpdgc);

  fCurrPathLengthList = new PathLengthList(*fCurrPDGCodeList);
  fCurrPathLengthList->SetPathLength(tgtpdgc,1.);

  LOG("PointGeom", pNOTICE) << *fCurrPDGCodeList;
  LOG("PointGeom", pNOTICE) << *fCurrPathLengthList;
}
//___________________________________________________________________________
PointGeomAnalyzer::~PointGeomAnalyzer()
{
  if( fCurrPathLengthList ) delete fCurrPathLengthList;
  if( fCurrPDGCodeList    ) delete fCurrPDGCodeList;
}
//___________________________________________________________________________
const PDGCodeList & PointGeomAnalyzer::ListOfTargetNuclei(void)
{
// pdg code list contains a single code corresponding to the material passed
// at the geom analyser ctor

  return *fCurrPDGCodeList;
}
//___________________________________________________________________________
const PathLengthList & PointGeomAnalyzer::ComputeMaxPathLengths(void)
{
// this is irrelevant for the 'point' geometry - return a path length of 1.
// for the only defined material

  return *fCurrPathLengthList;
}
//___________________________________________________________________________
const PathLengthList & PointGeomAnalyzer::ComputePathLengths(
                  const TLorentzVector & /*x*/, const TLorentzVector & /*p*/)
{
// this is irrelevant for the 'point' geometry - return a path length of 1.
// for the only defined material

  return *fCurrPathLengthList;
}
//___________________________________________________________________________
const TVector3 & PointGeomAnalyzer::GenerateVertex(
  const TLorentzVector & /*x*/, const TLorentzVector & /*p*/, int /*tgtpdg*/)
{
// this is irrelevant for the 'point' geometry - return a vtx at (0,0,0)

  return *fCurrVertex;
}
//___________________________________________________________________________
