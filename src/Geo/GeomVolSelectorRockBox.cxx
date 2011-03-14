//____________________________________________________________________________
/*
 Copyright (c) 2003-2011, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "Geo/GeomVolSelectorRockBox.h"
#include "Utils/StringUtils.h"

using namespace genie;
using namespace genie::geometry;

#include <TGeoVolume.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>

//____________________________________________________________________________
GeomVolSelectorRockBox::GeomVolSelectorRockBox() 
  : GeomVolSelectorFiducial(), fMinimumWall(0.), fDeDx(1.)
  , fRockBoxShape(0), fROOTGeom(0)
{
  fName = "RockBox";
  // base class' fiducial volume always treated as reverse (if even exists)
  this->SetReverseFiducial(true);
  
  for (int i=0; i<3; ++i) {
    fMinimalXYZMin[i] = 0;
    fMinimalXYZMax[i] = 0;
  }
}

//___________________________________________________________________________
GeomVolSelectorRockBox::~GeomVolSelectorRockBox()
{
  if ( fRockBoxShape ) delete fRockBoxShape;
  fRockBoxShape = 0;
  fROOTGeom = 0;  // was reference only
}

//___________________________________________________________________________
void GeomVolSelectorRockBox::TrimSegment(PathSegment& ps) const
{
  // First trim the segment based on the ray vs. cylinder or box 
  // Then trim futher according to the Basic parameters
  
  if ( ! fInterceptRock.fIsHit ) {
      // want in rock box, ray misses => reject all segments
      ps.fStepRangeSet.clear();  // 
  } else {
    // ray hit rock box volume, some segments steps need rejection, some need splitting...
    // check the steps in this segment
    Double_t dist = ps.fRayDist;
    StepRangeSet::iterator srs_itr = ps.fStepRangeSet.begin();
    StepRangeSet::iterator srs_end = ps.fStepRangeSet.end();
    StepRangeSet modifiedStepRangeSet;
    Bool_t ismod = false;
    
    // loop over steps within this segement 
    for ( ; srs_itr != srs_end; ++srs_itr ) {
      Double_t slo = srs_itr->first;
      Double_t shi = srs_itr->second;
      Bool_t split = false;
      StepRange step1, step2;
      // determine new trimmed or split steps
      ismod |= NewStepPairs(false,dist,slo,shi,
                            fInterceptRock,split,step1,step2);
      // build up new step list
      bool nonzerostep = ( step1.first != step1.second );
      if ( nonzerostep || ! fRemoveEntries ) {
        modifiedStepRangeSet.push_back(step1);
        if (split) {
          modifiedStepRangeSet.push_back(step2);
        }
      }
    } // loop over step range set elements
    if ( ismod ) ps.fStepRangeSet = modifiedStepRangeSet;

  } // fIsHit

  GeomVolSelectorFiducial::TrimSegment(ps);
}

//___________________________________________________________________________
void GeomVolSelectorRockBox::BeginPSList(const PathSegmentList* untrimmed) const
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.

  GeomVolSelectorFiducial::BeginPSList(untrimmed); 

  fCurrPathSegmentList = untrimmed;

  MakeRockBox();

  if ( ! fRockBoxShape ) {
    LOG("GeomVolSel", pFATAL) << "no shape defined";
    fInterceptRock = RayIntercept();
  } else {
    fInterceptRock = 
      fRockBoxShape->Intercept(fCurrPathSegmentList->GetStartPos(),
                               fCurrPathSegmentList->GetDirection());
  }

  //cout << "BeginPSList: " << endl
  //     << " fid:  " << fIntercept << endl
  //     << " rock: " << fInterceptRock << endl;

}

//___________________________________________________________________________
void GeomVolSelectorRockBox::EndPSList() const
{
  // Completed current path segment list processsing
}

//___________________________________________________________________________
void GeomVolSelectorRockBox::ConvertShapeMaster2Top(const ROOTGeomAnalyzer* rgeom)
{
  GeomVolSelectorFiducial::ConvertShapeMaster2Top(rgeom);
  fROOTGeom = rgeom;   // stash away a copy
}
//___________________________________________________________________________
void GeomVolSelectorRockBox::SetRockBoxMinimal(Double_t* xyzmin,
                                               Double_t* xyzmax)
{
  // This sets parameters for a minimal box

  for ( int j = 0; j < 3; ++j ) {
    fMinimalXYZMin[j] = TMath::Min(xyzmin[j],xyzmax[j]);
    fMinimalXYZMax[j] = TMath::Max(xyzmax[j],xyzmax[j]);
  }

}
//___________________________________________________________________________
void GeomVolSelectorRockBox::MakeRockBox() const
{
  // This sets parameters for a box

  // expanded box
  double energy = fP4.Energy();
  double boxXYZMin[3], boxXYZMax[3];
  for ( int j = 0; j < 3; ++j ) {
    double dmin = 0, dmax = 0;
    double dircos = fCurrPathSegmentList->GetDirection()[j];
    if ( dircos > 0 ) dmin =  dircos*energy/fDeDx;  // pad upstream
    else              dmax = -dircos*energy/fDeDx;
    // cout << "MakeRockBox wall " << fMinimumWall << " dm " << dmin << " " << dmax  << " dedx " << fDeDx << endl;
    boxXYZMin[j] = fMinimalXYZMin[j] - TMath::Max(fMinimumWall,dmin);
    boxXYZMax[j] = fMinimalXYZMax[j] + TMath::Max(fMinimumWall,dmax);
  }

  FidPolyhedron* poly = new FidPolyhedron();
  // careful about sign of "d" vs. direction normal
  PlaneParam pln0(-1,0,0, boxXYZMin[0]);  poly->push_back(pln0);
  PlaneParam pln1(0,-1,0, boxXYZMin[1]);  poly->push_back(pln1);
  PlaneParam pln2(0,0,-1, boxXYZMin[2]);  poly->push_back(pln2);
  PlaneParam pln3(+1,0,0,-boxXYZMax[0]);  poly->push_back(pln3);
  PlaneParam pln4(0,+1,0,-boxXYZMax[1]);  poly->push_back(pln4);
  PlaneParam pln5(0,0,+1,-boxXYZMax[2]);  poly->push_back(pln5);

  if ( fRockBoxShape ) delete fRockBoxShape;
  fRockBoxShape = poly;

  //cout << "MakeRockBox min [" 
  //     << boxXYZMin[0] << ","
  //     << boxXYZMin[1] << ","
  //     << boxXYZMin[2] << "]" << endl;
  //cout << "MakeRockBox max [" 
  //     << boxXYZMax[0] << ","
  //     << boxXYZMax[1] << ","
  //     << boxXYZMax[2] << "]" << endl;
  //cout << "rock before:" << *fRockBoxShape << endl;

  if ( fROOTGeom ) fRockBoxShape->ConvertMaster2Top(fROOTGeom);

  //cout << "rock after: " << *fRockBoxShape << endl;
  //cout << "fid  after: " << *fShape << endl;

}

//___________________________________________________________________________


