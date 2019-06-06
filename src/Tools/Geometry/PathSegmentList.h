//____________________________________________________________________________
/*!

\class   genie::geometry::PathSegmentList

\brief   Object to be filled with the neutrino path-segments representing
         geometry volume steps (generally boundary-to-boundary) along with
         geometry materials.  Good for a single starting position and 
         travelling along the direction of the neutrino 4-momentum.

\author  Robert Hatcher <rhatcher@fnal.gov>
         FNAL

\created May 26, 2009

\cpright  Copyright (c) 2003-2019, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org
          or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _PATH_SEGMENT_LIST_H_
#define _PATH_SEGMENT_LIST_H_

/// --- for test purposes allow compilation of class without string member
/// fetching/keeping the geometry path seems to add a significant (2x)
/// overhead to swimming through the geometry.
#define PATHSEG_KEEP_PATH
//#undef  PATHSEG_KEEP_PATH

#include <utility>  // for pair<>
#include <vector>
#include <list>
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
using std::pair;

namespace genie {
namespace geometry {

class PathSegment;
ostream & operator << (ostream & stream, const PathSegment & list);

typedef std::pair<Double_t,Double_t> StepRange;
typedef std::vector<StepRange>       StepRangeSet;

class PathSegment {

 public:
  PathSegment();
 ~PathSegment() { ; }

  /// point of entry to geometry element
  void SetEnter(const TVector3 & p3enter, double raydist) 
     { fEnter = p3enter; fRayDist = raydist; }
  void SetEnter(const Double_t * p3enter, double raydist) 
     { fEnter.SetXYZ(p3enter[0],p3enter[1],p3enter[2]); fRayDist = raydist; }

  /// point of exit from geometry element
  void SetExit(const TVector3 & p3exit) { fExit = p3exit;  }
  void SetExit(const Double_t * p3exit) 
     { fExit.SetXYZ(p3exit[0],p3exit[1],p3exit[2]); }

  /// info about the geometry element
  void SetGeo(const TGeoVolume * gvol, const TGeoMedium * gmed, 
              const TGeoMaterial * gmat)
  { fVolume = gvol; fMedium = gmed; fMaterial = gmat; }
#ifdef PATHSEG_KEEP_PATH
  void SetPath(const char* path) { fPathString = path; }
#endif

  /// step taken in the geometry element
  void SetStep(Double_t step, bool setlimits = true );

  bool IsTrimmedEmpty() const { return fStepRangeSet.empty(); }

  /// get the sum of all the step range (in case step has been trimmed or split)
  Double_t GetSummedStepRange() const;

  /// calculate position within allowed ranges passed on fraction of total
  TVector3 GetPosition(Double_t frac) const;

  /// perform cross check on segment, return differences
  void DoCrossCheck(const TVector3& startpos, double& ddist, double& dstep) const;

  //void Copy (const PathSegment & ps);
  //PathSegment& operator = (const PathSegment & ps);

  void Print (ostream & stream) const;
  friend ostream & operator << (ostream & stream, const PathSegment & list);
  friend bool      operator <  (const PathSegment &lhs, const PathSegment &rhs);
                                

  Double_t               fRayDist;      ///< distance from start of ray
  Double_t               fStepLength;   ///< total step size in volume
  const TGeoVolume *     fVolume;       ///< ref only ptr to TGeoVolume
  const TGeoMedium *     fMedium;       ///< ref only ptr to TGeoMedium
  const TGeoMaterial *   fMaterial;     ///< ref only ptr to TGeoMaterial
  TVector3               fEnter;        ///< top vol coordinates and units
  TVector3               fExit;         ///< top vol coordinates and units
#ifdef PATHSEG_KEEP_PATH
  std::string            fPathString;   ///< full path names
#endif
  StepRangeSet           fStepRangeSet; ///< collection of {steplo,stephi} pairs
};

inline bool operator < (const PathSegment &lhs, const PathSegment &rhs)
  { return ( lhs.fRayDist < rhs.fRayDist ); }


class PathSegmentList;
ostream & operator << (ostream & stream, const PathSegmentList & list);

class PathSegmentList {

public :
  PathSegmentList();
  PathSegmentList(const PathSegmentList & plist);
 ~PathSegmentList();

  void    SetDoCrossCheck (bool doit = true) { fDoCrossCheck = doit; }
  void    SetPrintVerbose (bool doit = true) { fPrintVerbose = doit; }
  void    SetAllToZero    (void);
  void    SetStartInfo    (const TVector3& pos = TVector3(0,0,1e37), 
                           const TVector3& dir = TVector3(0,0,0)     );
  bool    IsSameStart     (const TVector3& pos, const TVector3& dir) const;
  void    AddSegment      (const PathSegment& ps) { fSegmentList.push_back(ps); }

  const TVector3& GetDirection() const { return fDirection; }
  const TVector3& GetStartPos() const  { return fStartPos; }

  typedef std::list<PathSegment> PathSegmentV_t;
  typedef PathSegmentV_t::const_iterator PathSegVCItr_t;

  const   PathSegmentV_t&   GetPathSegmentV (void) const { return fSegmentList; }  
  size_t                    size(void) const { return fSegmentList.size(); }

  typedef std::map<const TGeoMaterial*,Double_t> MaterialMap_t;
  typedef MaterialMap_t::const_iterator MaterialMapCItr_t;

  void                      FillMatStepSum   (void);
  const   MaterialMap_t&    GetMatStepSumMap (void) const { return fMatStepSum; };

  void                      CrossCheck(double& mxddist, double& mxdstep) const;

#ifdef UNNEEDED_SEGFUNCS
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

  bool             fDoCrossCheck;
  bool             fPrintVerbose;

};

}      // geometry namespace
}      // genie    namespace

#endif // _PATH_SEGMENT_LIST_H_
