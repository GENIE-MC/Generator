//____________________________________________________________________________
/*!

\class    genie::geometry::GeomVolSelectorRockBox

\brief    GENIE Interface for limiting vertex selection in the rock
          to a volume that depends (in part) on the neutrino p4.
          Uses GeomVolSelectorFiducial to possibly exclude an inner region.

\author   Robert Hatcher <rhatcher@fnal.gov>
          FNAL

\created  August 5, 2010

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE

*/
//____________________________________________________________________________

#ifndef _GEOM_VOL_SELECTOR_ROCKBOX_H_
#define _GEOM_VOL_SELECTOR_ROCKBOX_H_


#include <string>
#include <vector>

#include "TMath.h"
#include "TLorentzVector.h"
#include "Tools/Geometry/FidShape.h"
#include "Tools/Geometry/PathSegmentList.h"
#include "Tools/Geometry/GeomVolSelectorFiducial.h"

using namespace std;

namespace genie {
namespace geometry {

class GeomVolSelectorRockBox : public GeomVolSelectorFiducial {

public :
           GeomVolSelectorRockBox();
  virtual ~GeomVolSelectorRockBox();

  //
  // define the missing part of the GeomVolSelectorI interface:
  //
  void TrimSegment(PathSegment& segment) const;
  void BeginPSList(const PathSegmentList* untrimmed) const;
  void EndPSList() const;

  //
  // set fiducial volume parameter (call only once)
  //   in "top vol" coordinates and units
  // set minimal (inner (optionally exclusion)) extent 
  //   before either wall or inclusion extent (which define region 
  //   where events are always accepted)
  void SetRockBoxMinimal(Double_t* xyzmin, Double_t* xyzmax);
  void SetRockBoxInclusion(Double_t* xyzmin, Double_t* xyzmax);
  void SetMinimumWall(Double_t w);
  void SetDeDx(Double_t dedx) { fDeDx = dedx; }
  void SetExpandFromInclusion(bool how=false) { fExpandInclusion = how; }

  // by default shapes are assumed to be in "top vol" coordinates
  // in the case where they are entered in master coordinates
  // ask the configured shape to convert itself  
  // (do this only once for any shape definition)
  virtual void ConvertShapeMaster2Top(const ROOTGeomAnalyzer* rgeom);

protected:

  void MakeRockBox() const;

  Double_t  fMinimalXYZMin[3];   /// interior box lower corner
  Double_t  fMinimalXYZMax[3];   /// interior box upper corner
  Double_t  fMinimumWall;        /// minimum distance around (XYZmin,XYZmax)
  Double_t  fInclusionXYZMin[3]; /// box within which events are always 
  Double_t  fInclusionXYZMax[3]; ///   accepted
  Double_t  fDeDx;               /// how to scale from energy to distance
  Bool_t    fExpandInclusion;    /// expand from minimal or inclusion box?

  mutable FidShape* fRockBoxShape;   /// shape changes for every nu ray
  
  const ROOTGeomAnalyzer* fROOTGeom;  // ref! only (for coordinate transforms, units)

  // values calculated during BeginPSList():
  mutable RayIntercept fInterceptRock;  // current intercept parameters

};

}      // geometry namespace
}      // genie    namespace

#endif // _GEOM_VOL_SELECTOR_ROCKBOX_H_
