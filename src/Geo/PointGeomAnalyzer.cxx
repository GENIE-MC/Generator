//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         STFC, Rutherford Appleton Laboratory - July 14, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 05, 2008 - CA
   Support a mix of targets (with their corresponding weights) in the same 
   'point geometry'.
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
PointGeomAnalyzer::PointGeomAnalyzer(int pdg) :
GeomAnalyzerI()
{
/*
  fCurrVertex = new TVector3(0,0,0);

  fCurrPDGCodeList = new PDGCodeList;
  fCurrPDGCodeList->clear();
  fCurrPDGCodeList->push_back(tgtpdgc);

  fCurrPathLengthList = new PathLengthList(*fCurrPDGCodeList);
  fCurrPathLengthList->SetPathLength(tgtpdgc,1.);

  LOG("PointGeom", pNOTICE) << *fCurrPDGCodeList;
  LOG("PointGeom", pNOTICE) << *fCurrPathLengthList;
*/
  const int    tgtpdgc[1] = { pdg };
  const double weight[1]  = { 1.0 };

  this->Initialize(1,tgtpdgc,weight);
}
//___________________________________________________________________________
PointGeomAnalyzer::PointGeomAnalyzer(
           unsigned int n, const int tgtpdgc[], const double weight[]) :
GeomAnalyzerI()
{
  this->Initialize(n,tgtpdgc,weight);
}
//___________________________________________________________________________
PointGeomAnalyzer::~PointGeomAnalyzer()
{
  this->CleanUp();
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
void PointGeomAnalyzer::Initialize(
          unsigned int n, const int tgtpdgc[], const double weight[]) 
{
  fCurrVertex = new TVector3(0,0,0);

  fCurrPDGCodeList = new PDGCodeList;
  fCurrPDGCodeList->clear();
  for(unsigned int i=0; i<n; i++) {
	  fCurrPDGCodeList->push_back(tgtpdgc[i]);
  }
  fCurrPathLengthList = new PathLengthList(*fCurrPDGCodeList);
  for(unsigned int i=0; i<n; i++) {
	  fCurrPathLengthList->SetPathLength(tgtpdgc[i], weight[i]);
  }
  LOG("PointGeom", pNOTICE) << *fCurrPDGCodeList;
  LOG("PointGeom", pNOTICE) << *fCurrPathLengthList;
}
//___________________________________________________________________________
void PointGeomAnalyzer::CleanUp(void)
{
  if( fCurrPathLengthList ) delete fCurrPathLengthList;
  if( fCurrPDGCodeList    ) delete fCurrPDGCodeList;
}
//___________________________________________________________________________

