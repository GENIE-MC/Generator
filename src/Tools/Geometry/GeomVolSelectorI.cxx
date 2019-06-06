//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ July 28, 2009 - RWH
   Was first added in v2.5.1.
   GeomVolSelectorI is the interface for a user defined volume selector.  
   Given a PathSegmentList* it returns a trimmed version (releasing ownership). 
   The selector can use info found in individual PathSegment element to 
   determine whether to keep, reject or trim elements.
 @ February 4, 2010 - RWH
   Make GenerateTrimmedList() a private method; users are now required to
   supply the TrimSegment() method for all derived classes. This way they 
   don't have to worry about maintaining the list only how to trim any 
   particular PathSegment.  Also, class now is configurable about whether to 
   retain a null segment; and also serves as a repository for the swimmer on 
   whether to fetch the geometry hierachy "path" (which turns out to be a 
   non-trivial overhead so we don't want to fetch it if we don't need to.

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
