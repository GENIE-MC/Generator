//____________________________________________________________________________
/*
 Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include "Messenger/Messenger.h"
#include "Geo/GeomVolSelectorFiducial.h"
#include "Geo/PathSegmentList.h"
#include "Utils/StringUtils.h"

using namespace genie;
using namespace genie::geometry;

#include <TGeoVolume.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>

const Double_t gLocalInfPos = 1.0e30;
const Double_t gLocalInfNeg = -gLocalInfPos;

//____________________________________________________________________________
GeomVolSelectorFiducial::GeomVolSelectorFiducial() 
  : GeomVolSelectorBasic(), fSelectReverse(false), fFidType(kUnknown),
    fCurrPathSegmentList(0)
{
  fName = "Fiducial";
}

//___________________________________________________________________________
GeomVolSelectorFiducial::~GeomVolSelectorFiducial()
{
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::TrimSegment(PathSegment& ps) const
{
  // First trim the segment based on the ray vs. cylinder or box 
  // Then trim futher according to the Basic parameters
  
  if ( ! fIsHit ) {
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
      Double_t slo_mod1, shi_mod1, slo_mod2, shi_mod2;
      // determine new trimmed or split steps
      ismod |= NewStepPairs(dist,slo,shi,split,
                            slo_mod1,shi_mod1,slo_mod2,shi_mod2);
      // build up new step list
      bool nonzerostep = ( slo_mod1 != shi_mod1 );
      if ( nonzerostep || ! fRemoveEntries ) {
        StepRange rng1(slo_mod1,shi_mod1);
        modifiedStepRangeSet.push_back(rng1);
        if (split) {
          StepRange rng2(slo_mod2,shi_mod2);
          modifiedStepRangeSet.push_back(rng2);
        }
      }
    } // loop over step range set elements
    if ( ismod ) ps.fStepRangeSet = modifiedStepRangeSet;

  } // fIsHit



  GeomVolSelectorBasic::TrimSegment(ps);
}

//___________________________________________________________________________
Bool_t GeomVolSelectorFiducial::NewStepPairs(Double_t raydist, 
                                             Double_t slo, Double_t shi, 
                                             Bool_t& split,
                                             Double_t& slo_mod1, Double_t& shi_mod1,
                                             Double_t& slo_mod2, Double_t& shi_mod2) const
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
  slo_mod1 = slo;
  shi_mod1 = shi;
  split = false;
  // dist in/out are relative to the ray origin
  Double_t sdistin  = fDistIn  - raydist;
  Double_t sdistout = fDistOut - raydist;
  // do remaining calculations relative to steps within this segment
  if ( slo < sdistin ) {
    if ( shi < sdistin ) {
      // case 0
      if ( fSelectReverse ) { ismodified = false; }
      else                  { slo_mod1 = 0; shi_mod1 = 0; }
    } else if ( shi < sdistout ) {
      // case 1
      if ( fSelectReverse ) { shi_mod1 = sdistin; }
      else                  { slo_mod1 = sdistin; }
    } else {
      // case 2
      if ( fSelectReverse ) { split = true;
        slo_mod1 = slo;      shi_mod1 = sdistin;
        slo_mod2 = sdistout; shi_mod2 = shi; }
      else                  { slo_mod1 = sdistin; shi_mod1 = sdistout; }
    }
  } else if ( slo < sdistout ) {
    if ( shi < sdistout ) {
      // case 3
      if ( fSelectReverse ) { slo_mod1 = 0; shi_mod1 = 0; }
      else                  { ismodified = false; }
    } else {
      // case 4
      if ( fSelectReverse ) { slo_mod1 = sdistout; }
      else                  { shi_mod1 = sdistout; }
    }
  } else {
    // case 5
    if ( fSelectReverse ) { ismodified = false; }
    else                  { slo_mod1 = 0; shi_mod1 = 0; }
  }

  return ismodified;
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::BeginPSList(const PathSegmentList* untrimmed) const
{
  LOG("GeomVolSel", pNOTICE)
    << "GeomVolSelectorFiducial::BeginPSList getting started";

  // A new neutrino ray has been set, calculate the entrance/exit distances.
  GeomVolSelectorBasic::BeginPSList(untrimmed);  // initialize base class in case it is needed

  fCurrPathSegmentList = untrimmed;

  fIsHit   = false;
  fDistIn  = gLocalInfNeg;
  fDistOut = gLocalInfPos;
  fSurfIn  = -1;
  fSurfOut = -1;

  switch ( fFidType ) {
  case kSphere:      CalcRay2Sphere(); break;
  case kZCylinder:   CalcRay2Cyl();    break;
  case kCylinder:    CalcRay2Cyl();    break;
  case kBox:         CalcRay2Poly();   break;
  case kConvexPoly:  CalcRay2Poly();   break;
  default:
    // do nothing
    LOG("GeomVolSel", pFATAL) << "no support for FidType=" << fFidType;
    break;
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
void GeomVolSelectorFiducial::SetSphere(Double_t x0, Double_t y0, Double_t z0, Double_t radius)
{
  fCenter  = TVector3(x0,y0,z0);
  fSRadius = radius;
  fFidType = kSphere;
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::SetZCylinder(Double_t x0, Double_t y0, Double_t radius, 
                                           Double_t zmin, Double_t zmax)
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.
  // This sets parameters for a cylinder parallel to the z-axis
  Double_t base[] = { x0, y0,  0 };
  Double_t axis[] = {  0,  0,  1 };
  Double_t cap1[] = {  0,  0, -1,  zmin };
  Double_t cap2[] = {  0,  0, +1, -zmax };  // note sign change
  this->SetCylinder(base,axis,radius,cap1,cap2);

  fFidType = kZCylinder;
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::SetCylinder(Double_t* base, Double_t* axis, Double_t radius,
                                          Double_t* cap1, Double_t* cap2)
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.
  // This sets parameters for a cylinder
  fFidType = kCylinder;
  fCylBase = TVector3(base);
  fCylAxis = TVector3(axis);
  fCylRadius = radius;
  fCylCap1 = PlaneParam(cap1);
  fCylCap2 = PlaneParam(cap2);
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::SetBox(Double_t* xyzmin, Double_t* xyzmax)
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.
  // This sets parameters for a box

  for ( int j = 0; j < 3; ++j ) {
    fBoxXYZMin[j] = TMath::Min(xyzmin[j],xyzmax[j]);
    fBoxXYZMax[j] = TMath::Max(xyzmax[j],xyzmax[j]);
  }

  fPolyFaces.clear();
  // careful about sign of "d" vs. direction normal
  PlaneParam pln0(-1,0,0, xyzmin[0]);  fPolyFaces.push_back(pln0);
  PlaneParam pln1(0,-1,0, xyzmin[1]);  fPolyFaces.push_back(pln1);
  PlaneParam pln2(0,0,-1, xyzmin[2]);  fPolyFaces.push_back(pln2);
  PlaneParam pln3(+1,0,0,-xyzmax[0]);  fPolyFaces.push_back(pln3);
  PlaneParam pln4(0,+1,0,-xyzmax[1]);  fPolyFaces.push_back(pln4);
  PlaneParam pln5(0,0,+1,-xyzmax[2]);  fPolyFaces.push_back(pln5);

  fFidType = kBox;
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::SetZPolygon(Int_t n, Double_t x0, Double_t y0, Double_t inradius, 
                                          Double_t phi0deg, Double_t zmin, Double_t zmax)
{
  // phi=0 will put flat face towards +x
  // inscribed radius
  fPolyFaces.clear();

  Double_t dphir = TMath::TwoPi()/n;
  Double_t phi0r = phi0deg * TMath::DegToRad();
  for ( int iface = 0; iface < n; ++iface ) {
    Double_t dx = TMath::Cos(dphir*iface + phi0r);
    Double_t dy = TMath::Sin(dphir*iface + phi0r);
    PlaneParam face(dx,dy,0,-inradius-dx*x0-dy*y0); fPolyFaces.push_back(face);
  }
  PlaneParam plnz1(0,0,-1, zmin);  fPolyFaces.push_back(plnz1);
  PlaneParam plnz2(0,0,+1,-zmax);  fPolyFaces.push_back(plnz2);

  fFidType = kConvexPoly;
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::CalcRay2Sphere() const
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.
  // This sets fDistIn/fDistOut for an sphere
  TVector3 oc = fCenter - fCurrPathSegmentList->GetStartPos();
  Double_t loc2 = oc.Mag2();
  Double_t r2 = fSRadius*fSRadius;
  //LOG("GeomVolSel", pNOTICE) << " loc2 = " << loc2 << " r2 " << r2;
  // if ( loc2 > r2 ) ray originates outside the sphere
  const TVector3& d = fCurrPathSegmentList->GetDirection();
  Double_t d2 = d.Mag2();
  Double_t tca = oc.Dot(d)/d2;
  //if ( tca < 0.0 ) sphere _center_ behind the ray orgin
  Double_t lhc2 = ( r2 -loc2 )/d2 + tca*tca;
  if ( lhc2 < 0.0 ) return; // ray misses the sphere
  fIsHit = true;
  Double_t lhc = TMath::Sqrt(lhc2);
  fDistIn  = tca - lhc;
  if ( fDistIn >= 0.0 ) fSurfIn  = 1;    // we hit
  else                  fDistIn  = 0.0;  // started inside
  fDistOut = tca + lhc;
  fSurfOut = 1;
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::CalcRay2CylUncapped() const
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.
  // This sets fDistIn/fDistOut for an infinite cylinder
  // Take as "hit" if the ray is parallel to the axis but inside the radius

  TVector3 rc   = fCurrPathSegmentList->GetStartPos() - fCylBase;
  TVector3 n    = fCurrPathSegmentList->GetDirection().Cross(fCylAxis);
  Double_t len  = n.Mag();
  Double_t dist = 0;
  if ( len == 0.0 ) {
    // ray is parallel to axis
    Double_t dist = rc.Dot(fCylAxis);
    TVector3 d = rc - dist*fCylAxis;
    dist = d.Mag();
    //LOG("GeomVolSel", pNOTICE) << " len = " << len << " is parallel, dist " << dist << ", " << fCylRadius;
    if ( dist <= fCylRadius ) {
      fIsHit   = true;  // inside is considered a hit
      fSurfIn  = 0;
      fSurfOut = 0;
      return;
    }
  }
  // ray is not parallel
  if ( len != 1.0 ) n.SetMag(1.);  // normalize if it isn't already
  dist = TMath::Abs(rc.Dot(n));    // closest approach distance
  //LOG("GeomVolSel", pNOTICE) << " len = " << len << " not parallel, dist " << dist << ", " << fCylRadius;
  if ( dist <= fCylRadius ) {
    fIsHit = true; // yes, it hits
    fSurfIn  = 0;
    fSurfOut = 0;
    TVector3 o = rc.Cross(fCylAxis);
    Double_t t = - o.Dot(n)/len;
    o = n.Cross(fCylAxis);
    o.SetMag(1.);
    Double_t s = TMath::Abs( TMath::Sqrt(fCylRadius*fCylRadius-dist*dist) /
                             fCurrPathSegmentList->GetDirection().Dot(o)    );
    fDistIn  = t - s;
    fDistOut = t + s;
    //LOG("GeomVolSel", pNOTICE) << " hits, t = " << t << " s = " << s;
  }
}

//___________________________________________________________________________
void GeomVolSelectorFiducial::CalcRay2Cyl() const
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.

  CalcRay2CylUncapped();  // find where ray hits the infinite cylinder
  // trim this down by applying the end caps
  if ( fIsHit ) {
    if ( fCylCap1.IsValid() ) {
      Double_t vd1 = fCylCap1.Vd(fCurrPathSegmentList->GetDirection());
      Double_t vn1 = fCylCap1.Vn(fCurrPathSegmentList->GetStartPos());
      if ( vd1 != 0.0 ) {
        Double_t t1 = -vn1 / vd1;
        if ( vd1 < 0.0 && t1 > fDistIn  ) { fDistIn  = t1;  fSurfIn  = 1; }
        if ( vd1 > 0.0 && t1 < fDistOut ) { fDistOut = t1;  fSurfOut = 1; }
      }
    }
    if ( fCylCap2.IsValid() ) {
      Double_t vd2 = fCylCap2.Vd(fCurrPathSegmentList->GetDirection());
      Double_t vn2 = fCylCap2.Vn(fCurrPathSegmentList->GetStartPos());
      if ( vd2 != 0.0 ) {
        Double_t t2 = -vn2 / vd2;
        if ( vd2 < 0.0 && t2 > fDistIn  ) { fDistIn  = t2;  fSurfIn  = 2; }
        if ( vd2 > 0.0 && t2 < fDistOut ) { fDistOut = t2;  fSurfOut = 2; }
      }
    }
  }

}

//___________________________________________________________________________
void GeomVolSelectorFiducial::CalcRay2Poly() const
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.
  // Calculate the point of intersection of a ray (directed line) and a
  // *convex* polyhedron constructed from the intersection of a list
  // of planar equations (no check on convex condition).

  Double_t tnear = gLocalInfNeg;
  Double_t tfar  = gLocalInfPos;
  Int_t surfNear = -1;
  Int_t surfFar  = -1;
  Bool_t parallel = false;

  // test each plane in the polyhedron
  for ( size_t iface=0; iface < fPolyFaces.size(); ++iface ) {
    const PlaneParam& pln = fPolyFaces[iface];
    if ( ! pln.IsValid() ) continue;

    // calculate numerator, denominator to "t" = distance along ray to intersection w/ pln
    Double_t vd = pln.Vd(fCurrPathSegmentList->GetDirection());
    Double_t vn = pln.Vn(fCurrPathSegmentList->GetStartPos());
    
    //LOG("GeomVolSel", pNOTICE)
    //  << " face " << iface << " [" << pln.a << "," << pln.b << "," << pln.c << "," << pln.d
    //  << "] vd=" << vd << " vn=" << vn;

    if ( vd == 0.0 ) {
      // ray is parallel to plane - check if ray origin is inside plane's half-space
      //LOG("GeomVolSel", pNOTICE)
      //  << " vd=0, " << " possibly parallel ";
      if ( vn > 0.0 ) { parallel = true; break; }  // wrong side ... complete miss
    } else {
      // ray is not parallel to plane -- get the distance to the plane
      Double_t t = -vn / vd; // notice negative sign!
      //LOG("GeomVolSel", pNOTICE) << "   t=" << t << " tnear=" << tnear << " tfar=" << tfar;
      if ( vd < 0.0 ) {
        // front face: t is a near point
        if ( t > tnear ) {
          surfNear = iface;
          tnear    = t;
        }
      } else {
        // back face: t is a far point
        if ( t < tfar ) {
          surfFar = iface;
          tfar    = t;
        }
      }
      //LOG("GeomVolSel", pNOTICE) << "     new surf " <<  surfNear << "," << surfFar 
      //                           << " tnear=" << tnear << " tfar=" << tfar;
    }
  }
  if ( ! parallel ) {
    // survived all the tests
    if ( tnear > 0.0 ) {
      if ( tnear < tfar ) {
        //LOG("GeomVolSel", pNOTICE) << "is hit case1 ";
        fIsHit   = true;
        fSurfIn  = surfNear;
        fSurfOut = surfFar;
      }
    } else {
      if ( tfar > 0.0 ) {
        //LOG("GeomVolSel", pNOTICE) << "is hit case2 ";
        fIsHit   = true;
        fSurfIn  = -1;
        fSurfOut = surfFar;
      }
    }
  }
  fDistIn  = tnear;
  fDistOut = tfar;
  //LOG("GeomVolSel", pNOTICE) << " hit? " << (fIsHit?"true":"false") 
  //                           << " dist in " << fDistIn << " out " << fDistOut;
}

//___________________________________________________________________________


