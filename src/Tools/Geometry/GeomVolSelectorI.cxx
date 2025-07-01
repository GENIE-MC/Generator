//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Robert Hatcher <rhatcher@fnal.gov>
*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Tools/Geometry/GeomVolSelectorI.h"
#include "Tools/Geometry/PathSegmentList.h"

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
  this->BeginPSList(untrimmed);

  PathSegmentList* trimmed = new PathSegmentList();
  trimmed->SetStartInfo(untrimmed->GetStartPos(),untrimmed->GetDirection());

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

  this->EndPSList();

  return trimmed;
}
//___________________________________________________________________________
