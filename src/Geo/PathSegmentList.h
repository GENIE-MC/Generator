//____________________________________________________________________________
/*!

\class   genie::PathSegmentList

\brief   Object to be filled with the neutrino path-segments representing
         geometry volume steps (generally boundary-to-boundary) along with
         geometry materials.  Good for a single starting position and 
         travelling along the direction of the neutrino 4-momentum.

\author  Robert Hatcher <rhatcher@fnal.gov>
         FNAL

\created May 26, 2009

\cpright  Copyright (c) 2009, GENIE Neutrino MC Generator Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PATH_SEGMENT_LIST_H_
#define _PATH_SEGMENT_LIST_H_

#include <vector>
#include <ostream>
#include <string>
#include <map>

#include <TVector3.h>
class TGeoVolume;
class TGeoMedium;
class TGeoMaterial;

using std::vector;
using std::ostream;
using std::string;

namespace genie {

class PDGCodeList;

class PathSegment {

 public:
  PathSegment();
 ~PathSegment() { ; }

  void SetEnter(const TVector3 & p3enter, double raydist) 
     { fEnter = p3enter; fRayDist = raydist; }
  void SetEnter(const Double_t * p3enter, double raydist) 
     { fEnter.SetXYZ(p3enter[0],p3enter[1],p3enter[2]); fRayDist = raydist; }
  void SetExit(const TVector3 & p3exit) { fExit = p3exit;  }
  void SetExit(const Double_t * p3exit) 
     { fExit.SetXYZ(p3exit[0],p3exit[1],p3exit[2]); }
  void SetStep(Double_t step ) { fStepLength = step; }
  void SetGeo(const TGeoVolume * gvol, const TGeoMedium * gmed, 
              const TGeoMaterial * gmat)
  { fVolume = gvol; fMedium = gmed; fMaterial = gmat; }

  //void Copy (const PathSegment & ps);
  //PathSegment& operator = (const PathSegment & ps);

  void Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const PathSegment & list);

  TVector3               fEnter;      ///< top vol coordinates and units
  TVector3               fExit;       ///< top vol coordinates and units
  Double_t               fStepLength;
  Double_t               fRayDist;    ///< distance from start of ray
  const TGeoVolume *     fVolume;     ///< ref only ptr to TGeoVolume
  const TGeoMedium *     fMedium;     ///< ref only ptr to TGeoMedium
  const TGeoMaterial *   fMaterial;   ///< ref only ptr to TGeoMaterial
};


class PathSegmentList {

public :
  PathSegmentList();
  PathSegmentList(const PathSegmentList & plist);
 ~PathSegmentList();

  void    SetAllToZero    (void);
  void    SetStartInfo    (const TVector3& pos = TVector3(0,0,1e37), 
                          const TVector3& dir = TVector3(0,0,0)     );
  bool    IsSameStart     (const TVector3& pos, const TVector3& dir) const;
  void    AddSegment      (PathSegment& ps) { fSegmentList.push_back(ps); }

  typedef std::vector<PathSegment> PathSegmentV_t;
  typedef PathSegmentV_t::const_iterator PathSegVCItr_t;

  const   PathSegmentV_t&   GetPathSegmentV (void) { return fSegmentList; }  
  size_t                    size(void) const { return fSegmentList.size(); }

  typedef std::map<const TGeoMaterial*,Double_t> MaterialMap_t;
  typedef MaterialMap_t::const_iterator MaterialMapCItr_t;

  void                      FillMatStepSum   (void);
  const   MaterialMap_t&    GetMatStepSumMap (void) { return fMatStepSum; };

#ifdef UNNEEDED_SEGFUNCS
  void   AddPathLength   (int pdgc, double pl); // path-legth(pdgc) += pl
  void   SetPathLength   (int pdgc, double pl); // path-legth(pdgc)  = pl
  bool   AreAllZero      (void) const;
  void   ScalePathLength (int pdgc, double scale);
  double PathLength      (int pdgc) const;

  //  XmlParserStatus_t LoadFromXml (string filename);
  //  void              SaveAsXml   (string filename) const;

#endif

  void Copy  (const PathSegmentList & plist);
  PathSegmentList & operator =  (const PathSegmentList & list);

  void Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const PathSegmentList & list);

 protected:

  /// Record, for future comparison, the path taken
  TVector3         fStartPos;  ///< starting position (in top vol coords)
  TVector3         fDirection; ///< direction (in top vol coords)

  /// Actual list of segments
  PathSegmentV_t   fSegmentList;

  /// Segment list re-evaluated by material for fast lookup of path lengths
  MaterialMap_t    fMatStepSum;

};

}      // genie namespace

#endif // _PATH_SEGMENT_LIST_H_
