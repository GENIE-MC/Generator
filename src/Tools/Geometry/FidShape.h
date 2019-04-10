//____________________________________________________________________________
/*!

\class    genie::geometry::FidShape

\brief    Some simple volumes that know how to calculate where a ray 
          intercepts them.

\author   Robert Hatcher <rhatcher@fnal.gov>
          FNAL

\created  August 3, 2010

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

          Some of the algorithms here are (loosely) based on those found in:
          Graphics Gems II, ISBN 0-12-064480-0
          pg. 247 "Fast Ray-Convex Polyhedron Intersection"
          Graphics Gems IV, ed. Paul Heckbert, ISBN 0-12-336156-7  T385.G6974 (1994)
          pg. 356 "Intersecting a Ray with a Cylinder"
*/
//____________________________________________________________________________

#ifndef _FID_SHAPE_H_
#define _FID_SHAPE_H_

#include <vector>
#include <cfloat> // for DBL_MAX

#include "TMath.h"
#include "TLorentzVector.h"

namespace genie {
namespace geometry {

class ROOTGeomAnalyzer;

class PlaneParam;
std::ostream& operator<< (std::ostream& stream, 
                           const genie::geometry::PlaneParam& pparam);

class RayIntercept {
  /// A class to hold information about where a ray intercepts a
  /// convex shape.  
  public:
  RayIntercept() : fDistIn(-DBL_MAX), fDistOut(DBL_MAX), 
      fIsHit(false), fSurfIn(-1), fSurfOut(-1) { ; }
  ~RayIntercept() { ; }
  Double_t  fDistIn;   /// distance along ray to enter fid volume
  Double_t  fDistOut;  /// distance along ray to exit fid volume
  Bool_t    fIsHit;    /// was the volume hit
  Int_t     fSurfIn;   /// what surface was hit on way in
  Int_t     fSurfOut;  /// what surface was hit on way out
};
std::ostream& operator<< (std::ostream& stream, 
                          const genie::geometry::RayIntercept& ri);

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
  void     ConvertMaster2Top(const ROOTGeomAnalyzer* rgeom);
  void     Print(std::ostream& stream) const;
  friend std::ostream& operator<< (std::ostream& stream, 
                                   const genie::geometry::PlaneParam& pparam);

  Double_t a, b, c, d; // the parameters
};

class FidShape;
std::ostream& operator<< (std::ostream& stream, 
                          const genie::geometry::FidShape& shape);

class FidShape {
  // generic fiducial shape
  public:
  FidShape() { ; } 
  virtual ~FidShape() { ; }
  /// derived classes must implement the Intercept() method
  /// which calculates the entry/exit point of a ray w/ the shape
  virtual RayIntercept Intercept(const TVector3& start, const TVector3& dir) const = 0; 
  /// derived classes must implement the ConvertMaster2Top() method
  /// which transforms the shape specification from master coordinates to "top vol"
  virtual void ConvertMaster2Top(const ROOTGeomAnalyzer* rgeom) = 0;
  virtual void Print(std::ostream& stream) const = 0;
  friend std::ostream& operator<< (std::ostream& stream, 
                              const genie::geometry::FidShape& shape);

};

class FidSphere : public FidShape {
 public:
 FidSphere(const TVector3& center, Double_t radius) : fCenter(center), fSRadius(radius) { ; }
 RayIntercept Intercept(const TVector3& start, const TVector3& dir) const;
 void         ConvertMaster2Top(const ROOTGeomAnalyzer* rgeom);
 void         Print(std::ostream& stream) const;
 protected:
 TVector3    fCenter;   /// center of the sphere
 Double_t    fSRadius;  /// radius of the sphere
};

class FidCylinder : public FidShape {
 public:
 FidCylinder(const TVector3& base, const TVector3& axis, Double_t radius, 
             const PlaneParam& cap1, const PlaneParam& cap2) 
   : fCylBase(base), fCylAxis(axis), fCylRadius(radius), fCylCap1(cap1), fCylCap2(cap2) { ; }
 RayIntercept Intercept(const TVector3& start, const TVector3& dir) const;
 RayIntercept InterceptUncapped(const TVector3& start, const TVector3& dir) const;
 void         ConvertMaster2Top(const ROOTGeomAnalyzer* rgeom);
 void         Print(std::ostream& stream) const;
 protected:

 TVector3    fCylBase;   /// base point on cylinder axis
 TVector3    fCylAxis;   /// direction cosines of cylinder axis
 Double_t    fCylRadius; /// radius of cylinder
 PlaneParam  fCylCap1;   /// define a plane for 1st cylinder cap
 PlaneParam  fCylCap2;   /// define a plane for 2nd cylinder cap
};

class FidPolyhedron : public FidShape {
  /// convex polyhedron is made of multiple planar equations
 public:
 FidPolyhedron() { ; }
 void push_back(const PlaneParam& pln) { fPolyFaces.push_back(pln); }
 void clear() { fPolyFaces.clear(); }
 RayIntercept Intercept(const TVector3& start, const TVector3& dir) const;
 void         ConvertMaster2Top(const ROOTGeomAnalyzer* rgeom);
 void         Print(std::ostream& stream) const;
 protected:
 std::vector<PlaneParam> fPolyFaces;  /// the collection of planar equations for the faces
};

}      // geometry namespace
}      // genie    namespace

#endif // _FID_SHAPE_H_
