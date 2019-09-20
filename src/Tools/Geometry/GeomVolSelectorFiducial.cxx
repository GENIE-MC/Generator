//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Tools/Geometry/GeomVolSelectorFiducial.h"
#include "Framework/Utils/StringUtils.h"

using namespace genie;
using namespace genie::geometry;

#include <TGeoVolume.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>

//____________________________________________________________________________
GeomVolSelectorFiducial::GeomVolSelectorFiducial() 
  : GeomVolSelectorBasic(), fSelectReverse(false)
  , fShape(0), fCurrPathSegmentList(0)
{
  fName = "Fiducial";
}

//___________________________________________________________________________
GeomVolSelectorFiducial::~GeomVolSelectorFiducial()
{
  if ( fShape ) delete fShape;
  fShape = 0;
  fCurrPathSegmentList = 0;  // was reference only
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::TrimSegment(PathSegment& ps) const
{
  // First trim the segment based on the ray vs. cylinder or box 
  // Then trim futher according to the Basic parameters
  
  if ( ! fIntercept.fIsHit ) {
    // simple case when ray doesn't intersect the fiducial volume at all
    if ( fSelectReverse ) {
      // fiducial is reversed (ie. only want regions outside) 
      /// so a miss means blindly accept all segments
    } else {
      // want in fiducial, ray misses => reject all segments
      ps.fStepRangeSet.clear();  // 
    }
  } else {
    // ray hit fiducial volume, some segments steps need rejection, some need splitting...
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
      ismod |= NewStepPairs(fSelectReverse,dist,slo,shi,
                            fIntercept,split,step1,step2);
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

  GeomVolSelectorBasic::TrimSegment(ps);
}

//___________________________________________________________________________
Bool_t GeomVolSelectorFiducial::NewStepPairs(Bool_t selectReverse, Double_t raydist, 
                                             Double_t slo, Double_t shi,
                                             const RayIntercept& intercept, 
                                             Bool_t& split,
                                             StepRange& step1, StepRange& step2)
{
  // modifying a step based on a range
  // there seem to be six possible cases:    
  //     step: |===| {slo,shi}
  //
  //   range:   in       out
  //            #        #            normal     reverse
  // 0   |===|  #        #            {0,0}      {slo,shi} 
  // 1       |==#==|     #            {in,shi}   {slo,in}
  // 2       |==#========#==|         {in,out}   {slo,in}+{out,shi}
  // 3          #  |==|  #            {slo,shi}  {0,0}
  // 4          #     |==#==|         {slo,out}  {out,shi} 
  // 5          #        #  |====|    {0,0}      {slo,shi}
  //            #        #
  //
  bool ismodified = true;
  step1 = StepRange(slo,shi);
  split = false;
  // dist in/out are relative to the ray origin
  Double_t sdistin  = intercept.fDistIn  - raydist;
  Double_t sdistout = intercept.fDistOut - raydist;
  // do remaining calculations relative to steps within this segment
  if ( slo < sdistin ) {
    if ( shi < sdistin ) {
      // case 0
      if ( selectReverse ) { ismodified = false; }
      else                 { step1 = StepRange(0,0); }
    } else if ( shi < sdistout ) {
      // case 1
      if ( selectReverse ) { step1.second = sdistin; }
      else                 { step1.first  = sdistin; }
    } else {
      // case 2
      if ( selectReverse ) { split = true;
        step1 = StepRange(slo,sdistin); 
        step2 = StepRange(sdistout,shi); }
      else                 { step1 = StepRange(sdistin,sdistout); }
    }
  } else if ( slo < sdistout ) {
    if ( shi < sdistout ) {
      // case 3
      if ( selectReverse ) { step1 = StepRange(0,0); }
      else                 { ismodified = false; }
    } else {
      // case 4
      if ( selectReverse ) { step1.first  = sdistout; }
      else                 { step1.second = sdistout; }
    }
  } else {
    // case 5
    if ( selectReverse ) { ismodified = false; }
    else                 { step1 = StepRange(0,0); }
  }

  return ismodified;
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::BeginPSList(const PathSegmentList* untrimmed) const
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.

  GeomVolSelectorBasic::BeginPSList(untrimmed);  // initialize base class in case it is needed

  fCurrPathSegmentList = untrimmed;

  if ( ! fShape ) {
    LOG("GeomVolSel", pFATAL) << "no shape defined";
    fIntercept = RayIntercept();
  } else {
    fIntercept = fShape->Intercept(fCurrPathSegmentList->GetStartPos(),
                                   fCurrPathSegmentList->GetDirection());
  }

}

//___________________________________________________________________________
void GeomVolSelectorFiducial::EndPSList() const
{
  // Completed current path segment list processsing
  GeomVolSelectorBasic::EndPSList();
  fCurrPathSegmentList = 0;
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::AdoptFidShape(FidShape* shape)
{
  if ( fShape ) delete fShape;
  fShape = shape;
}
//___________________________________________________________________________
void GeomVolSelectorFiducial::ConvertShapeMaster2Top(const ROOTGeomAnalyzer* rgeom)
{
  if ( fShape ) fShape->ConvertMaster2Top(rgeom);
}
//___________________________________________________________________________
void GeomVolSelectorFiducial::MakeSphere(Double_t x0, Double_t y0, Double_t z0, Double_t radius)
{
  AdoptFidShape(new FidSphere(TVector3(x0,y0,z0),radius));
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::MakeXCylinder(Double_t y0, Double_t z0, Double_t radius,
                                            Double_t xmin, Double_t xmax)
{
  // This sets parameters for a cylinder parallel to the x-axis
  Double_t base[] = {  0, y0, z0 };
  Double_t axis[] = {  1,  0,  0 };
  Double_t cap1[] = { -1,  0,  0,  xmin };
  Double_t cap2[] = { +1,  0,  0, -xmax };  // note sign change
  this->MakeCylinder(base,axis,radius,cap1,cap2);
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::MakeYCylinder(Double_t x0, Double_t z0, Double_t radius,
                                            Double_t ymin, Double_t ymax)
{
  // This sets parameters for a cylinder parallel to the y-axis
  Double_t base[] = { x0,  0, z0 };
  Double_t axis[] = {  0,  1,  0 };
  Double_t cap1[] = {  0, -1,  0,  ymin };
  Double_t cap2[] = {  0, +1,  0, -ymax };  // note sign change
  this->MakeCylinder(base,axis,radius,cap1,cap2);
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::MakeZCylinder(Double_t x0, Double_t y0, Double_t radius,
                                            Double_t zmin, Double_t zmax)
{
  // This sets parameters for a cylinder parallel to the z-axis
  Double_t base[] = { x0, y0,  0 };
  Double_t axis[] = {  0,  0,  1 };
  Double_t cap1[] = {  0,  0, -1,  zmin };
  Double_t cap2[] = {  0,  0, +1, -zmax };  // note sign change
  this->MakeCylinder(base,axis,radius,cap1,cap2);
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::MakeCylinder(Double_t* base, Double_t* axis, Double_t radius,
                                           Double_t* cap1, Double_t* cap2)
{
  // This sets parameters for a cylinder
  AdoptFidShape(new FidCylinder(TVector3(base),TVector3(axis),radius,
                                PlaneParam(cap1),PlaneParam(cap2)));
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::MakeBox(Double_t* xyzmin, Double_t* xyzmax)
{
  // This sets parameters for a box

  double boxXYZmin[3], boxXYZmax[3];
  for ( int j = 0; j < 3; ++j ) {
    boxXYZmin[j] = TMath::Min(xyzmin[j],xyzmax[j]);
    boxXYZmax[j] = TMath::Max(xyzmin[j],xyzmax[j]);
  }

  FidPolyhedron* poly = new FidPolyhedron();
  // careful about sign of "d" vs. direction normal
  PlaneParam pln0(-1,0,0, boxXYZmin[0]);  poly->push_back(pln0);
  PlaneParam pln1(0,-1,0, boxXYZmin[1]);  poly->push_back(pln1);
  PlaneParam pln2(0,0,-1, boxXYZmin[2]);  poly->push_back(pln2);
  PlaneParam pln3(+1,0,0,-boxXYZmax[0]);  poly->push_back(pln3);
  PlaneParam pln4(0,+1,0,-boxXYZmax[1]);  poly->push_back(pln4);
  PlaneParam pln5(0,0,+1,-boxXYZmax[2]);  poly->push_back(pln5);

  AdoptFidShape(poly);
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::MakeZPolygon(Int_t n, Double_t x0, Double_t y0, Double_t inradius, 
                                           Double_t phi0deg, Double_t zmin, Double_t zmax)
{
  // phi=0 will put flat face towards +x
  // inscribed radius

  FidPolyhedron* poly = new FidPolyhedron();

  Double_t dphir = TMath::TwoPi()/n;
  Double_t phi0r = phi0deg * TMath::DegToRad();
  for ( int iface = 0; iface < n; ++iface ) {
    Double_t dx = TMath::Cos(dphir*iface + phi0r);
    Double_t dy = TMath::Sin(dphir*iface + phi0r);
    PlaneParam face(dx,dy,0,-inradius-dx*x0-dy*y0); poly->push_back(face);
  }
  PlaneParam plnz1(0,0,-1, zmin);  poly->push_back(plnz1);
  PlaneParam plnz2(0,0,+1,-zmax);  poly->push_back(plnz2);

  AdoptFidShape(poly);
}

//___________________________________________________________________________


