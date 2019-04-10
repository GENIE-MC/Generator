//____________________________________________________________________________
/*!

\class    genie::geometry::GeomVolSelectorBasic

\brief    GENIE Interface for user-defined volume selector functors
          This basic version allows configurations that depend on PathSegment elements'
          material/media/volume and/or "path"

\author   Robert Hatcher <rhatcher@fnal.gov>
          FNAL

\created  December 3, 2008

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _GEOM_VOL_SELECTOR_BASIC_H_
#define _GEOM_VOL_SELECTOR_BASIC_H_

#include <string>
#include <vector>
#include "TLorentzVector.h"
#include "Tools/Geometry/GeomVolSelectorI.h"

using namespace std;

namespace genie {
namespace geometry {

class PathSegmentList;

class GeomVolSelectorBasic : public GeomVolSelectorI {

public :
           GeomVolSelectorBasic();
  virtual ~GeomVolSelectorBasic();

  ///
  ///  Selections are string based, elements are specified as a list of items separated by
  ///  comma, semicolon or colons.  Elements that start with "-" are rejections; elements
  ///  that start with "+" (or nothing) are required, e.g.
  ///     "+N276B,-air0"
  ///
  void     SetVolumeSelection(string volstr);
  void     SetMediumSelection(string medstr);
  void     SetMaterialSelection(string matstr);
  void     SetPathSelection(string pathstr);

  //
  // define the missing parts of the GeomVolSelectorI interface:
  //
  void TrimSegment(PathSegment& segment) const;
  void BeginPSList(const PathSegmentList* untrimmed) const;
  void EndPSList() const;

protected:

  void ParseSelection(const string& str, vector<string>& required, vector<string>& forbidden);
  bool RejectString(const string& str, const vector<string>& required, const vector<string>& forbidden) const;

  // PathSegment must contain one of the things in these lists (if there are any)
  vector<string> fRequiredVol;
  vector<string> fRequiredMed;
  vector<string> fRequiredMat;
  vector<string> fRequiredPath;

  // PathSegment must not contain any of the things in these lists
  vector<string> fForbiddenVol;
  vector<string> fForbiddenMed;
  vector<string> fForbiddenMat;
  vector<string> fForbiddenPath;

};

}      // geometry namespace
}      // genie    namespace

#endif // _GEOM_VOL_SELECTOR_BASIC_H_
