//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ February 04, 2010 - RWH
   Was first added in v2.5.1.
   Introduce GeomVolSelectorBasic as concrete example of GeomVolSelectorI that 
   accepts/rejects path segments based on volume name, material, medium, or 
   path in geometry hierarchy.

*/
//____________________________________________________________________________

#include "Framework/Messenger/Messenger.h"
#include "Tools/Geometry/GeomVolSelectorBasic.h"
#include "Tools/Geometry/PathSegmentList.h"
#include "Framework/Utils/StringUtils.h"

using namespace genie;
using namespace genie::geometry;

#include <TGeoVolume.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>

//____________________________________________________________________________
GeomVolSelectorBasic::GeomVolSelectorBasic() 
  : GeomVolSelectorI("Basic")
{

}

//___________________________________________________________________________
GeomVolSelectorBasic::~GeomVolSelectorBasic()
{

}

//___________________________________________________________________________
void GeomVolSelectorBasic::SetVolumeSelection(string volstr)
{
  ParseSelection(volstr,fRequiredVol,fForbiddenVol);
}
void GeomVolSelectorBasic::SetMediumSelection(string medstr)
{
  ParseSelection(medstr,fRequiredMed,fForbiddenMed);
}
void GeomVolSelectorBasic::SetMaterialSelection(string matstr)
{
  ParseSelection(matstr,fRequiredMat,fForbiddenMat);
}
void GeomVolSelectorBasic::SetPathSelection(string pathstr)
{
  ParseSelection(pathstr,fRequiredPath,fForbiddenPath);
  if ( fRequiredPath.size()  > 0 || fForbiddenPath.size() > 0 ) {
#ifdef PATHSEG_KEEP_PATH
    SetNeedPath();
#else
    LOG("GeomVolSelectorBasic", pFATAL)
      << "PathSegment is not defined to hold fPathString -- selectors can not cut on it";
#endif
  }
}

//___________________________________________________________________________
// nothing special needs to be done at the beginning or end of a new PathSegmentList
void GeomVolSelectorBasic::BeginPSList(const PathSegmentList* /* untrimmed */ ) const 
{ ; }

void GeomVolSelectorBasic::EndPSList() const 
{ ; }

//___________________________________________________________________________
void GeomVolSelectorBasic::TrimSegment(PathSegment& ps) const
{
  bool reject = false;
  
  // not splitting PathSegment into 2 or more PathSegment elements
  // so either 
  //    - keep "as is"
  //    - adjust the low/high endpoints
  //    - not copy to output list (
  // be careful about steps outside all the geometry that might not
  // have an associated volume/media/material
  
  
  if ( ! reject ) {
    std::string volname = ( ps.fVolume) ? ps.fVolume->GetName() : "no-volume";
    reject = RejectString(volname,fRequiredVol,fForbiddenVol);
  }
  
  if ( ! reject ) {
    std::string medname = ( ps.fMedium) ? ps.fMedium->GetName() : "no-medium";
    reject = RejectString(medname,fRequiredMed,fForbiddenMed);
  }
  
  if ( ! reject ) {
    std::string matname = ( ps.fMaterial) ? ps.fMaterial->GetName() : "no-material";
    reject = RejectString(matname,fRequiredMat,fForbiddenMat);
  }
  
#ifdef PATHSEG_KEEP_PATH
  if ( ! reject ) {
    reject = RejectString(ps.fPathString,fRequiredPath,fForbiddenPath);
  }
#endif

  if ( reject ) ps.fStepRangeSet.clear();

}

//___________________________________________________________________________
void GeomVolSelectorBasic::ParseSelection(const string& strall, 
                                          vector<string>& required, 
                                          vector<string>& forbidden)
{
  required.clear();
  forbidden.clear();
  vector<string> pieces = genie::utils::str::Split(strall,":;,");
  size_t n = pieces.size();
  for ( size_t i = 0; i < n; ++i ) {
    string& strone = pieces[i];
    if ( strone == "" ) continue;  // reject null strings
    if      ( strone.find("-") == 0 ) forbidden.push_back(strone.substr(1,std::string::npos));
    else if ( strone.find("+") == 0 ) required.push_back(strone.substr(1,std::string::npos));
    else                              required.push_back(strone);
  }
}
//___________________________________________________________________________
bool GeomVolSelectorBasic::RejectString(const string& str, 
                                        const vector<string>& required, 
                                        const vector<string>& forbidden) const
{
  bool reject = false;

  // must have at least one of the required elements (if there are any)
  size_t nrequired = required.size();
  if ( nrequired > 0 ) {
    bool found = false;
    for (size_t jr = 0; jr < nrequired; ++jr) {
      if ( str.find(required[jr]) != std::string::npos ) {
        found = true;
        break; // found at least one case, so we're good
      }
    }
    if ( ! found ) reject = true;
  }

  // can not have any of the forbidden elements
  size_t nforbidden = forbidden.size();
  if ( nforbidden > 0 ) {
    for (size_t jf = 0; jf < nforbidden; ++jf) {
      if ( str.find(forbidden[jf]) != std::string::npos ) {
        reject = true;
        break; // found at least one case, so we can reject
      }
    }
  }

  return reject;
}
//___________________________________________________________________________

