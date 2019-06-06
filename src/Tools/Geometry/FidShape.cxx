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
#include <iomanip>
using namespace std;

#include "Tools/Geometry/FidShape.h"
#include "Tools/Geometry/ROOTGeomAnalyzer.h"

using namespace genie;
using namespace genie::geometry;

//___________________________________________________________________________
std::ostream& 
genie::geometry::operator<< (std::ostream& stream, 
                             const genie::geometry::PlaneParam& pparam)
{
  pparam.Print(stream);
  return stream;
}

std::ostream& 
genie::geometry::operator<< (std::ostream& stream, 
                             const genie::geometry::RayIntercept& ri)
{
  stream << "RayIntercept: dist in/out " << ri.fDistIn << "/" << ri.fDistOut
         << " hit=" << ((ri.fIsHit)?"true":"false")
         << " surf " << ri.fSurfIn << "/" << ri.fSurfOut;
  return stream;
}

std::ostream& 
genie::geometry::operator<< (std::ostream& stream, 
                             const genie::geometry::FidShape& shape)
{
  shape.Print(stream);
  return stream;
}

//___________________________________________________________________________
void PlaneParam::ConvertMaster2Top(const ROOTGeomAnalyzer* rgeom)
{
  // convert a plane equation from master coordinates to "top vol" coordinates
  TVector3 dir(a,b,c);
  rgeom->Master2TopDir(dir);  // new a,b,c
  TVector3 zero(0,0,0);
  rgeom->Master2Top(zero);
  a = dir.X();
  b = dir.Y();
  c = dir.Z();
  d = d - ( a*zero.X() + b*zero.Y() + c*zero.Z() );

}

//___________________________________________________________________________
void PlaneParam::Print(std::ostream& stream) const
{
  stream << "PlaneParam=[" << a << "," << b << "," << c << "," << d << "]"; 
}

//___________________________________________________________________________
RayIntercept FidSphere::Intercept(const TVector3& start, const TVector3& dir) const
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.
  // This sets fDistIn/fDistOut for an sphere
  RayIntercept intercept;

  TVector3 oc = fCenter - start;
  Double_t loc2 = oc.Mag2();
  Double_t r2 = fSRadius*fSRadius;
  //LOG("GeomVolSel", pNOTICE) << " loc2 = " << loc2 << " r2 " << r2;
  // if ( loc2 > r2 ) ray originates outside the sphere
  const TVector3& d = dir;
  Double_t d2 = d.Mag2();
  Double_t tca = oc.Dot(d)/d2;
  //if ( tca < 0.0 ) sphere _center_ behind the ray orgin
  Double_t lhc2 = ( r2 -loc2 )/d2 + tca*tca;
  if ( lhc2 < 0.0 ) return intercept; // ray misses the sphere
  intercept.fIsHit = true;
  Double_t lhc = TMath::Sqrt(lhc2);

  intercept.fDistIn  = tca - lhc;
  intercept.fSurfIn  = 1;
  intercept.fDistOut = tca + lhc;
  intercept.fSurfOut = 1;

  return intercept;
}

//___________________________________________________________________________
void FidSphere::ConvertMaster2Top(const ROOTGeomAnalyzer* rgeom)
{
  rgeom->Master2Top(fCenter);
}

//___________________________________________________________________________
void FidSphere::Print(std::ostream& stream) const
{
  stream << "FidSphere @ ["
         << fCenter.X() << ","
         << fCenter.Y() << ","
         << fCenter.Z() << "]  r = " << fSRadius;
}

//___________________________________________________________________________
RayIntercept FidCylinder::InterceptUncapped(const TVector3& start, const TVector3& dir) const
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.
  // This sets fDistIn/fDistOut for an infinite cylinder
  // Take as "hit" if the ray is parallel to the axis but inside the radius
  RayIntercept intercept;

  TVector3 rc   = start - fCylBase;
  TVector3 n    = dir.Cross(fCylAxis);
  Double_t len  = n.Mag();
  Double_t dist = 0;
  if ( len == 0.0 ) {
    // ray is parallel to axis
    dist = rc.Dot(fCylAxis);
    TVector3 d = rc - dist*fCylAxis;
    dist = d.Mag();
    //LOG("GeomVolSel", pNOTICE) << " len = " << len << " is parallel, dist " << dist << ", " << fCylRadius;
    if ( dist <= fCylRadius ) {
      intercept.fIsHit   = true;  // inside is considered a hit
      intercept.fSurfIn  = 0;
      intercept.fSurfOut = 0;
    }
    return intercept;
  }
  // ray is not parallel
  if ( len != 1.0 ) n.SetMag(1.);  // normalize if it isn't already
  dist = TMath::Abs(rc.Dot(n));    // closest approach distance
  //LOG("GeomVolSel", pNOTICE) << " len = " << len << " not parallel, dist " << dist << ", " << fCylRadius;
  if ( dist <= fCylRadius ) {
    intercept.fIsHit = true; // yes, it hits
    intercept.fSurfIn  = 0;
    intercept.fSurfOut = 0;
    TVector3 o = rc.Cross(fCylAxis);
    Double_t t = - o.Dot(n)/len;
    o = n.Cross(fCylAxis);
    o.SetMag(1.);
    Double_t s = TMath::Abs( TMath::Sqrt(fCylRadius*fCylRadius-dist*dist) /
                             dir.Dot(o)    );
    intercept.fDistIn  = t - s;
    intercept.fDistOut = t + s;
    //LOG("GeomVolSel", pNOTICE) << " hits, t = " << t << " s = " << s;
  }
  return intercept;
}

//___________________________________________________________________________
RayIntercept FidCylinder::Intercept(const TVector3& start, const TVector3& dir) const
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.

  RayIntercept intercept = InterceptUncapped(start,dir);  // find where ray hits the infinite cylinder
  // trim this down by applying the end caps
  if ( ! intercept.fIsHit ) return intercept;
  for ( int icap=1; icap <= 2; ++icap ) {
    const PlaneParam& cap = (icap==1) ? fCylCap1 : fCylCap2;
    if ( ! cap.IsValid() ) continue;
    Double_t vd = cap.Vd(dir);
    Double_t vn = cap.Vn(start);
    //std::cout << "FidCyl::Intercept cap " << icap 
    //          << " vd " << vd << " vn " << vn;
    if ( vd == 0.0 ) { // parallel to surface, is it on the right side?
      //std::cout << " vd=0, vn " << ((vn>0)?"wrong":"right") << "side " << std::endl;
      if ( vn > 0 ) { intercept.fIsHit = false; break; } // wrong side
    } else {
      Double_t t = -vn / vd;
      //std::cout << " t " << t << " in/out "
      //          << intercept.fDistIn << "/" << intercept.fDistOut << std::endl;
      if ( vd < 0.0 ) { // t is the entering point
        if ( t > intercept.fDistIn  ) 
          { intercept.fDistIn  = t;  intercept.fSurfIn  = 1; }
      } else { // t is the exiting point
        if ( t < intercept.fDistOut ) 
          { intercept.fDistOut = t;  intercept.fSurfOut = 1; }
      }
    }
  }
  if ( intercept.fDistIn > intercept.fDistOut ) intercept.fIsHit = false;
  return intercept;
}

//___________________________________________________________________________
void FidCylinder::ConvertMaster2Top(const ROOTGeomAnalyzer* rgeom)
{
  rgeom->Master2Top(fCylBase);
  rgeom->Master2TopDir(fCylAxis);
  fCylCap1.ConvertMaster2Top(rgeom);
  fCylCap2.ConvertMaster2Top(rgeom);
}

//___________________________________________________________________________
void FidCylinder::Print(std::ostream& stream) const
{
  stream << "FidCylinder @ ["
         << fCylBase.X() << ","
         << fCylBase.Y() << ","
         << fCylBase.Z() << "]  dir ["
         << fCylAxis.X() << ","
         << fCylAxis.Y() << ","
         << fCylAxis.Z() << "]  r = " << fCylRadius;
  stream << " cap1=" << fCylCap1 << " cap2=" << fCylCap2;
}

//___________________________________________________________________________
RayIntercept FidPolyhedron::Intercept(const TVector3& start, const TVector3& dir) const
{
  // A new neutrino ray has been set, calculate the entrance/exit distances.
  // Calculate the point of intersection of a ray (directed line) and a
  // *convex* polyhedron constructed from the intersection of a list
  // of planar equations (no check on convex condition).

  RayIntercept intercept;

  Double_t tnear = -DBL_MAX;
  Double_t tfar  =  DBL_MAX;
  Int_t surfNear = -1;
  Int_t surfFar  = -1;
  Bool_t parallel = false;

  // test each plane in the polyhedron
  for ( size_t iface=0; iface < fPolyFaces.size(); ++iface ) {
    const PlaneParam& pln = fPolyFaces[iface];
    if ( ! pln.IsValid() ) continue;

    // calculate numerator, denominator to "t" = distance along ray to intersection w/ pln
    Double_t vd = pln.Vd(dir);
    Double_t vn = pln.Vn(start);
    
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
        intercept.fIsHit   = true;
        intercept.fSurfIn  = surfNear;
        intercept.fSurfOut = surfFar;
      }
    } else {
      if ( tfar > 0.0 ) {
        //LOG("GeomVolSel", pNOTICE) << "is hit case2 ";
        intercept.fIsHit   = true;
        intercept.fSurfIn  = -1;
        intercept.fSurfOut = surfFar;
      }
    }
  }
  intercept.fDistIn  = tnear;
  intercept.fDistOut = tfar;
  //LOG("GeomVolSel", pNOTICE) << " hit? " << (fIsHit?"true":"false") 
  //                           << " dist in " << fDistIn << " out " << fDistOut;
  return intercept;
}

//___________________________________________________________________________
void FidPolyhedron::ConvertMaster2Top(const ROOTGeomAnalyzer* rgeom)
{
  for (unsigned int i = 0; i < fPolyFaces.size(); ++i ) {
    PlaneParam& aplane = fPolyFaces[i];
    aplane.ConvertMaster2Top(rgeom);
  }
}
void FidPolyhedron::Print(std::ostream& stream) const
{
  size_t nfaces = fPolyFaces.size();
  stream << "FidPolyhedron n=" << nfaces;
  for ( size_t i=0; i<nfaces; ++i) {
    const PlaneParam& aface = fPolyFaces[i];
    stream << std::endl << "[" << setw(2) << i << "]" << aface;
  }
}

//___________________________________________________________________________
