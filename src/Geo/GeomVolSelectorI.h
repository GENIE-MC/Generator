//____________________________________________________________________________
/*!

\class    genie::geometry::GeomVolSelectorI

\brief    GENIE Interface for user-defined volume selector functors

\author   Robert Hatcher <rhatcher@fnal.gov>
          FNAL

\created  December 3, 2008

\cpright  Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
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

class PathSegmentList;

class GeomVolSelectorI {

public :
  virtual ~GeomVolSelectorI();
  //
  // define the GeomVolSelectorI interface:
  //

  /// create and return a new PathSegmentList from the old list
  /// relinquishes ownership of returned object
  virtual PathSegmentList* GenerateTrimmedList(const PathSegmentList* untrimmed) const = 0;

  /// configure for individual neutrino ray
  void SetCurrentRay(const TLorentzVector& x4, const TLorentzVector& p4)
  { fX4 = x4; fP4 = p4; }

  /// set scale factor for SI to "raydist" units of PathSegmentList
  void SetSI2Local(double scale) { fScale = scale; }

 std::string GetName() const { return fName; }

protected:

 GeomVolSelectorI();
 GeomVolSelectorI(std::string name);

 TLorentzVector     fX4;    ///< current neutrino ray's start position (global)
 TLorentzVector     fP4;    ///< current neutrino ray's momentum (global)
 double             fScale; ///< SI->raydist scale factor
 std::string        fName;  ///< volume selector name                                    

};

}      // geometry namespace
}      // genie    namespace

#endif // _GEOM_VOL_SELECTOR_I_H_
