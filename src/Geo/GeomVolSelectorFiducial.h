//____________________________________________________________________________
/*!

\class    genie::geometry::GeomVolSelectorFiducial

\brief    GENIE Interface for user-defined volume selector functors
          Trim path segments based on the intersection with a cylinder, box
          or sphere as well as everything the Basic selector can do.

          Assumes that the fiducial volume is defined in the same coords
          and units as the PathSegmentList ("top vol") and that the ray
          always starts outside the defined volume.  (If not user should
          cap the fid volume just down from the flux window or use the
          SetUpstreamZ() in the flux driver to push the ray back to make it so).

\author   Robert Hatcher <rhatcher@fnal.gov>
          FNAL

\created  July 14, 2010

\cpright  Copyright (c) 2003-2010, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

          Some of the algorithms here are (loosely) based on those found in:
          Graphics Gems II, ISBN 0-12-064480-0
          pg. 247 "Fast Ray-Convex Polyhedron Intersection"
          Graphics Gems IV, ed. Paul Heckbert, ISBN 0-12-336156-7  T385.G6974 (1994)
          pg. 356 "Intersecting a Ray with a Cylinder"
*/
//____________________________________________________________________________

#ifndef _GEOM_VOL_SELECTOR_FIDUCIAL_H_
#define _GEOM_VOL_SELECTOR_FIDUCIAL_H_

#include <string>
#include <vector>
#include "TMath.h"
#include "TLorentzVector.h"
#include "Geo/GeomVolSelectorBasic.h"

using namespace std;

namespace genie {
namespace geometry {

class PathSegmentList;

class PlaneParam {
  // A plane is described by the equation a*x +b*y + c*z + d = 0
  // n = [a,b,c] are the plane normal components  (one must be non-zero)
  // d is the distance to the origin
  // for a point "p" on the plane:  d = - p.n    (note the "-")
 public:
  PlaneParam(Double_t ain=0, Double_t bin=0, Double_t cin=0, Double_t din=0)
    { a = ain; b = bin; c = cin; d = din; Normalize(); }
  PlaneParam(Double_t* abcd)
    { a = abcd[0]; b = abcd[1]; c = abcd[2]; d = abcd[3]; Normalize(); }

  void     Normalize()  // make the a,b,c parameters a unit normal
    { Double_t mag = TMath::Sqrt(a*a+b*b+c*c); 
      if (mag>0) { a /= mag; b /= mag; c /= mag; d /= mag; } }
  Double_t Vn(const TVector3& raybase) const
    { return raybase.X()*a + raybase.Y()*b + raybase.Z()*c + d; }
  Double_t Vd(const TVector3& raycos) const
    { return raycos.Px()*a + raycos.Py()*b + raycos.Pz()*c; }
  Bool_t   IsValid() const { return (a != 0 || b != 0 || c != 0 ); }

  Double_t a, b, c, d; // the parameters
};


class GeomVolSelectorFiducial : public GeomVolSelectorBasic {

public :
           GeomVolSelectorFiducial();
  virtual ~GeomVolSelectorFiducial();

  //
  // define the missing part of the GeomVolSelectorI interface:
  //
  void TrimSegment(PathSegment& segment) const;
  void BeginPSList(const PathSegmentList* untrimmed) const;
  void EndPSList() const;

  //
  // set fiducial volume parameter (call only one)
  // in "top vol" coordinates and units
  //
  void SetSphere(Double_t x0, Double_t y0, Double_t z0, Double_t radius);
  void SetZCylinder(Double_t x0, Double_t y0, Double_t radius, Double_t zmin, Double_t zmax);
  void SetCylinder(Double_t* base, Double_t* axis, Double_t radius, Double_t* cap1, Double_t* cap2);
  void SetBox(Double_t* xyzmin, Double_t* xyzmax);
  void SetZPolygon(Int_t n, Double_t x0, Double_t y0, Double_t inradius, Double_t phi0deg, Double_t zmin, Double_t zmax);

  void SetReverseFiducial(Bool_t reverse=true) { fSelectReverse = reverse; }

protected:

  // const (but modify mutable values)
  void CalcRay2Sphere() const;
  void CalcRay2CylUncapped() const;
  void CalcRay2Cyl() const;
  void CalcRay2Poly() const;

  Bool_t NewStepPairs(Double_t raydist, Double_t slo, Double_t shi, Bool_t& split,
                      Double_t& slo_mod1, Double_t& shi_mod1,
                      Double_t& slo_mod2, Double_t& shi_mod2) const;

  enum EFidType {
    kUnknown    = 0,
    kSphere     = 1,
    kZCylinder  = 2,
    kCylinder   = 3,
    kBox        = 4,
    kConvexPoly = 5
  };

  Bool_t fSelectReverse; /// select for "outside" fiducial?
  EFidType  fFidType;    /// what type of fiducial volume
  
  // Sphere parameters
  TVector3   fCenter;
  Double_t   fSRadius;

  // Cylinder parameters  (surface normals point out, especially caps)
  TVector3   fCylBase;   /// base point on cylinder axis
  TVector3   fCylAxis;   /// direction cos of cylinder axis
  Double_t   fCylRadius; /// radius of cylinder
  PlaneParam fCylCap1;   /// define a plane for 1st cylinder cap
  PlaneParam fCylCap2;   /// define a plane for 2nd cylinder cap

  // Box parameters
  Double_t fBoxXYZMin[3];  /// box minimum xyz positions
  Double_t fBoxXYZMax[3];  /// box maximum xyz positions

  // Polyhedron parameters  (surface normals point out)
  std::vector<PlaneParam> fPolyFaces;  /// the collection of planar equations for the faces

  // values calculated during BeginPSList():
  mutable const PathSegmentList* fCurrPathSegmentList;  // reference only, for ray info
  mutable Double_t  fDistIn;   /// distance along ray to enter fid volume
  mutable Double_t  fDistOut;  /// distance along ray to exit fid volume
  mutable Bool_t    fIsHit;    /// was the volume hit
  mutable Int_t     fSurfIn;   /// what surface was hit on way in
  mutable Int_t     fSurfOut;  /// what surface was hit on way out

};

}      // geometry namespace
}      // genie    namespace

#endif // _GEOM_VOL_SELECTOR_FIDUCIAL_H_
