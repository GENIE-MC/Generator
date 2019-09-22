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

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

*/
//____________________________________________________________________________

#ifndef _GEOM_VOL_SELECTOR_FIDUCIAL_H_
#define _GEOM_VOL_SELECTOR_FIDUCIAL_H_

#include <string>
#include <vector>

#include "TMath.h"
#include "TLorentzVector.h"
#include "Tools/Geometry/FidShape.h"
#include "Tools/Geometry/PathSegmentList.h"
#include "Tools/Geometry/GeomVolSelectorBasic.h"

using namespace std;

namespace genie {
namespace geometry {

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

  // allow the selection to be reversed (i.e. exclude "fid" region)
  void SetReverseFiducial(Bool_t reverse=true) { fSelectReverse = reverse; }

  //
  // set fiducial volume parameter (call only one)
  // in "top vol" coordinates and units
  //
  void AdoptFidShape(FidShape* shape);
  void MakeSphere(Double_t x0, Double_t y0, Double_t z0, Double_t radius);
  void MakeXCylinder(Double_t y0, Double_t z0, Double_t radius, Double_t xmin, Double_t xmax);
  void MakeYCylinder(Double_t x0, Double_t z0, Double_t radius, Double_t ymin, Double_t ymax);
  void MakeZCylinder(Double_t x0, Double_t y0, Double_t radius, Double_t zmin, Double_t zmax);
  void MakeCylinder(Double_t* base, Double_t* axis, Double_t radius, Double_t* cap1, Double_t* cap2);
  void MakeBox(Double_t* xyzmin, Double_t* xyzmax);
  void MakeZPolygon(Int_t n, Double_t x0, Double_t y0, Double_t inradius, Double_t phi0deg, Double_t zmin, Double_t zmax);

  // by default shapes are assumed to be in "top vol" coordinates
  // in the case where they are entered in master coordinates
  // ask the configured shape to convert itself  
  // (do this only once for any shape definition)
  virtual void ConvertShapeMaster2Top(const ROOTGeomAnalyzer* rgeom);

protected:

  static Bool_t NewStepPairs(Bool_t selectReverse,
                             Double_t raydist, Double_t slo, Double_t shi, 
                             const RayIntercept& intercept, Bool_t& split,
                             StepRange& step1, StepRange& step2);

  Bool_t fSelectReverse; /// select for "outside" fiducial?

  FidShape* fShape;   /// shape

  // values calculated during BeginPSList():
  mutable const PathSegmentList* fCurrPathSegmentList;  // reference only, for ray info
  mutable RayIntercept fIntercept;  // current intercept parameters

};

}      // geometry namespace
}      // genie    namespace

#endif // _GEOM_VOL_SELECTOR_FIDUCIAL_H_
