//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - July 14, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Feb 05, 2008 - CA
   Support a mix of targets (with their corresponding weights) in the same 
   'point geometry'.
*/
//____________________________________________________________________________

#include <TLorentzVector.h>
#include <TVector3.h>

#include "Tools/Geometry/PointGeomAnalyzer.h"
#include "Framework/EventGen/PathLengthList.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"

using namespace genie;
using namespace genie::geometry;

//___________________________________________________________________________
PointGeomAnalyzer::PointGeomAnalyzer(int pdg) :
GeomAnalyzerI()
{
  map<int,double> tgtmap;
  tgtmap.insert( map<int, double>::value_type(pdg, 1.) );

  this->Initialize(tgtmap);
}
//___________________________________________________________________________
PointGeomAnalyzer::PointGeomAnalyzer(
           unsigned int n, const int tgtpdgc[], const double weight[]) :
GeomAnalyzerI()
{
  map<int,double> tgtmap;
  for(unsigned int i=0; i<n; i++) 
     tgtmap.insert( map<int, double>::value_type(tgtpdgc[i], weight[i]) );

  this->Initialize(tgtmap);
}
//___________________________________________________________________________
PointGeomAnalyzer::PointGeomAnalyzer(const map<int,double> & tgtmap) :
GeomAnalyzerI()
{
  this->Initialize(tgtmap);
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
void PointGeomAnalyzer::Initialize(const map<int,double> & tgtmap)
{
  fCurrVertex = new TVector3(0,0,0);

  fCurrPDGCodeList = new PDGCodeList;
  fCurrPDGCodeList->clear();

  map<int,double>::const_iterator iter;
  for(iter = tgtmap.begin(); iter != tgtmap.end(); ++iter) {
	int tgtpdgc = iter->first;
	fCurrPDGCodeList->push_back(tgtpdgc);
  }

  fCurrPathLengthList = new PathLengthList(tgtmap);

  LOG("PointGeom", pNOTICE) << *fCurrPDGCodeList;
  LOG("PointGeom", pNOTICE) << *fCurrPathLengthList;
}
//___________________________________________________________________________
void PointGeomAnalyzer::CleanUp(void)
{
  if( fCurrVertex )         delete fCurrVertex;
  if( fCurrPathLengthList ) delete fCurrPathLengthList;
  if( fCurrPDGCodeList    ) delete fCurrPDGCodeList;
}
//___________________________________________________________________________

