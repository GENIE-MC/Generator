//____________________________________________________________________________
/*!

\class    genie::geometry::GeomVolSelectorI

\brief    GENIE Interface for user-defined volume selector functors

\author   Robert Hatcher <rhatcher@fnal.gov>
          FNAL

\created  August 25, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GEOM_VOL_SELECTOR_I_H_
#define _GEOM_VOL_SELECTOR_I_H_

#include <string>
#include "TLorentzVector.h"

namespace genie {
namespace geometry {

class PathSegment;
class PathSegmentList;

class GeomVolSelectorI {

public :
  virtual ~GeomVolSelectorI();
  //
  // define the GeomVolSelectorI interface:
  //

  /// create and return a new PathSegmentList from the old list
  /// relinquishes ownership of returned object
  virtual PathSegmentList* GenerateTrimmedList(const PathSegmentList* untrimmed) const;

  /// This is the method every derived version must implement
  /// To reject a segment outright:  segment.fStepRangeSet.clear()
  virtual void TrimSegment(PathSegment& segment) const = 0;

  /// Every derived version must also respond to a signal that starts
  /// a new path segment list processing and ends it.
  /// In general they can simply ignore the signal.
  /// If the derived class needs to cache something, make it mutable
  virtual void BeginPSList(const PathSegmentList* untrimmed) const = 0;
  virtual void EndPSList() const = 0;

  /// configure for individual neutrino ray
  void SetCurrentRay(const TLorentzVector& x4, const TLorentzVector& p4)
  { fX4 = x4; fP4 = p4; }

  /// set scale factor for SI to "raydist" units of PathSegmentList
  void SetSI2Local(double scale) { fScale = scale; }

  void        SetRemoveEntries(bool rmset) { fRemoveEntries = rmset; }
  bool        GetRemoveEntries()           { return fRemoveEntries; }

  void        SetNeedPath()                { fNeedPath = true; }  /// allow toggle *on* only
  bool        GetNeedPath() const          { return fNeedPath; }

  std::string GetName() const              { return fName; }

protected:

 GeomVolSelectorI();
 GeomVolSelectorI(std::string name);

 TLorentzVector     fX4;             ///< current neutrino ray's start position (global)
 TLorentzVector     fP4;             ///< current neutrino ray's momentum (global)
 double             fScale;          ///< SI->raydist scale factor
 bool               fRemoveEntries;  ///< whether selector should remove entries or set hi=lo
 bool               fNeedPath;       ///< selector needs PathSegment "path" string
 std::string        fName;           ///< volume selector name    

};

}      // geometry namespace
}      // genie    namespace

#endif // _GEOM_VOL_SELECTOR_I_H_
