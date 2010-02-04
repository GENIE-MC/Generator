//____________________________________________________________________________
/*
 Copyright (c) 2003-2008, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

//#include "EVGDrivers/GeomVolSelectorI.h"
//rwh -- for now during test phase
#include "Geo/GeomVolSelectorI.h"
#include "Geo/PathSegmentList.h"

using namespace genie;
using namespace genie::geometry;

//____________________________________________________________________________
GeomVolSelectorI::GeomVolSelectorI() 
  : fRemoveEntries(false), fNeedPath(false), fName("no-name")
{
}

//____________________________________________________________________________
GeomVolSelectorI::GeomVolSelectorI(std::string name) 
  : fRemoveEntries(false), fNeedPath(false), fName(name)
{
}

//___________________________________________________________________________
GeomVolSelectorI::~GeomVolSelectorI()
{

}
//___________________________________________________________________________

PathSegmentList* 
GeomVolSelectorI::GenerateTrimmedList(const PathSegmentList* untrimmed) const
{
  PathSegmentList* trimmed = new PathSegmentList();

  const genie::geometry::PathSegmentList::PathSegmentV_t& segments = 
    untrimmed->GetPathSegmentV();
  genie::geometry::PathSegmentList::PathSegVCItr_t sitr = segments.begin();
  genie::geometry::PathSegmentList::PathSegVCItr_t sitr_end = segments.end();

  for ( ; sitr != sitr_end ; ++sitr ) {
    PathSegment ps = *sitr;  // new PathSegment is a copy of old
    this->TrimSegment(ps);
    if ( fRemoveEntries && ps.GetSummedStepRange() == 0 ) continue; // remove null segments
    // now put (adjusted) entry on trimmed list
    trimmed->AddSegment(ps);
  }

  return trimmed;
}
//___________________________________________________________________________
