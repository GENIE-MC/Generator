//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>
         FNAL - May 26, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ May 26, 2009 - RWH
   First added in v2.5.1
 @ July 16, 2009 - RWH
   Reworked ROOTGeomAnalyzer code so that each neutrino from the flux is only
   swum through the geometry once (not once per tgt PDG + 2 more to pick vertex 
   if interaction is selected). The PathSegmentList class holds all the individual 
   steps taken and is amenable to the next step which is trimming the list to 
   restrict the interactions to a non-physically realized volume.
 @ July 27, 2009 - RWH
   Print() also shows start pos and dir
 @ August 3, 2009 - RWH
   Individual path segments can be further trimmed to be less than the full   
   step (which corresponds to a volume).       
   The list is now internally a STL list (no longer a STL vector).
   This is in anticipation of a possible need for segment splitting if 
   trimming causes disjoint segments.
   Both classes have provisions for cross checking ray distance and step size 
   values after completion of the list. Generally I've found ray distance 
   errors of no more than ~8e-9 cm (80000fm) in a 42000cm swim through a 
   geometry. This is on order an atomic size so one shouldn't worry about it.  
   This precision will allow trimming to be based on the distance travelled 
   along the ray.  
 @ August 8, 2009 - RWH
   Push PathSegmentList down a namespace so it's in genie::geometry.  
   Fix const'ness of some functions. Add SetPath() function and fPathString 
   data member for recording geometry path info in case users want to cut on 
   something there.
 @ February 4, 2010 - RWH
   Correct statement about how overhead much fetching geometry path adds
   (significant, not negligable).     
   Generalize from (lo,hi) pair to vector of (lo,hi) pairs to allow the 
   geometry step to be split and not just squeezed.  This might not be 
   necessary and how much this adds to overhead isn't quite clear 
   (needs more testing).          
   New methods IsTrimmedEmpty() and GetSummedStepRange().  
   Also GetPosition() to pick position within vector of (lo,hi) pairs based on 
   fraction of total.  Needs more testing for case of split segments.

*/
//____________________________________________________________________________

#include <fstream>
#include <iomanip>
#include <sstream>
#include <cassert>

//#include "libxml/parser.h"
//#include "libxml/xmlmemory.h"

#include <TLorentzVector.h>
#include <TGeoVolume.h>
#include <TGeoMaterial.h>

#include "Tools/Geometry/PathSegmentList.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodeList.h"
#include "Framework/ParticleData/PDGLibrary.h"
#include "Framework/Utils/PrintUtils.h"
#include "Framework/Numerical/MathUtils.h"
//#include "Framework/Utils/XmlParserUtils.h"

using std::ofstream;
using std::setw;
using std::setfill;
using std::endl;

using namespace genie;
using namespace genie::geometry;

//#define RWH_DEBUG

//____________________________________________________________________________
namespace genie {
namespace geometry {

ostream & operator << (ostream & stream, const geometry::PathSegment & ps)
 {
   ps.Print(stream);
   return stream;
 }

ostream & operator << (ostream & stream, const geometry::PathSegmentList & list)
 {
   list.Print(stream);
   return stream;
 }

} // namespace geometry
} // namespace genie

//____________________________________________________________________________
namespace genie {
  namespace pathsegutils {
    string Vec3AsString (const TVector3 * vec)
    {
      int w=17, p=10;  // precision setting only affects ostringstream
      std::ostringstream fmt;
      fmt << "(" << std::setw(w)   << std::setprecision(p) << vec->x()
          << "," << std::setw(w)   << std::setprecision(p) << vec->y()
          << "," << std::setw(w+1) << std::setprecision(p) << vec->z() << ")";
      return fmt.str();
    }
  }
}

//___________________________________________________________________________
PathSegment::PathSegment(void) :
  fRayDist(0), fStepLength(0), 
  fVolume(0), fMedium(0), fMaterial(0),
  fEnter(), fExit()
{
}

//___________________________________________________________________________
void PathSegment::DoCrossCheck(const TVector3& startpos, 
                               double& ddist, double& dstep) const
{
  double dist_recalc = (fEnter-startpos).Mag();
  ddist = dist_recalc - fRayDist;

  double step_recalc = (fExit-fEnter).Mag();
  dstep = step_recalc - fStepLength;
}

//___________________________________________________________________________
void PathSegment::Print(ostream & stream) const
{
  const char* vname = (fVolume)   ? fVolume->GetName()   : "no volume";
  const char* mname = (fMaterial) ? fMaterial->GetName() : "no material";
  stream << genie::pathsegutils::Vec3AsString(&fEnter) << " " 
    //<< genie::pathsegutils::Vec3AsString(&fExit) 
         << " " // "raydist " 
         << std::setw(12) << fRayDist 
         << " " // "step " 
         << std::setw(12) << fStepLength << " "
         << std::left
         << std::setw(16) << vname << " '"
         << std::setw(18) << mname << "' ";
  size_t n = fStepRangeSet.size();
  const int rngw = 24;
  if ( n == 0 ) {
    stream << std::setw(rngw) << "[ ]";
  } else {
    std::ostringstream rngset;
    for ( size_t i = 0 ; i < n; ++i ) {
      const StepRange& sr = fStepRangeSet[i];
      rngset << "[" << sr.first << ":" << sr.second << "]";
    }
    stream << std::setw(rngw) << rngset.str();
  }
  stream << std::right;
#ifdef PATHSEG_KEEP_PATH
  stream << " " << fPathString;
#endif

}

//___________________________________________________________________________
void PathSegment::SetStep(Double_t step, bool setlimits)
{ 
  fStepLength = step; 
  if (setlimits) {
    fStepRangeSet.clear();
    fStepRangeSet.push_back(StepRange(0,step));
  }
}

//___________________________________________________________________________
Double_t PathSegment::GetSummedStepRange() const
{
  Double_t sum = 0;
  for ( size_t i = 0; i < fStepRangeSet.size(); ++i ) {
    const StepRange& sr = fStepRangeSet[i];
    sum += ( sr.second - sr.first );
  }
  return sum;
}

TVector3 PathSegment::GetPosition(Double_t fractrim) const
{
  /// calculate position within allowed ranges passed as
  /// fraction of trimmed segment
  ///      seg.fEnter + fractotal * ( seg.fExit - seg.fEnter );
  Double_t sumrange = GetSummedStepRange();
  if ( sumrange < 0.0 ) {
    LOG("PathS", pFATAL) << "GetPosition failed fractrim=" << fractrim
                         << " because sumrange = " << sumrange;
    return TVector3(0,0,0);
  }
  Double_t target   = fractrim * sumrange;
  Double_t sum      = 0;
  for ( size_t i = 0; i < fStepRangeSet.size(); ++i ) {
    const StepRange& sr = fStepRangeSet[i];
    Double_t ds = ( sr.second - sr.first );
    sum += ds;
#ifdef RWH_DEBUG
    LOG("PathS", pINFO) << "GetPosition fractrim=" << fractrim
                        << " target " << target << " [" << i << "] "
                        << " ds " << ds << " sum " << sum;
#endif
    if ( sum >= target ) {
      Double_t overstep  = sum - target;
      Double_t fractotal = (sr.second - overstep)/fStepLength;
#ifdef RWH_DEBUG
    LOG("PathS", pINFO) << "GetPosition fractrim=" << fractrim
                        << " overstep " << overstep
                        << " fractotal " << fractotal;
#endif
      return  fEnter + fractotal * ( fExit - fEnter );
    }
  }
  LOG("PathS", pFATAL) << "GetPosition failed fractrim=" << fractrim;
  return TVector3(0,0,0);
}


//===========================================================================
//___________________________________________________________________________
PathSegmentList::PathSegmentList(void)
  : fDoCrossCheck(false), fPrintVerbose(false)
{

}
//___________________________________________________________________________
PathSegmentList::PathSegmentList(const PathSegmentList & plist)
{
  this->Copy(plist);
}
//__________________________________________________________________________
PathSegmentList::~PathSegmentList()
{

}

//___________________________________________________________________________
void PathSegmentList::SetAllToZero(void)
{
  LOG("PathS", pDEBUG) << "SetAllToZero called";

  this->fStartPos.SetXYZ(0,0,1e37); // clear cache of position/direction
  this->fDirection.SetXYZ(0,0,0);   //
  this->fSegmentList.clear();       // clear the vector
  this->fMatStepSum.clear();        // clear the re-factorized info
}

//___________________________________________________________________________
void PathSegmentList::SetStartInfo(const TVector3& pos, const TVector3& dir)
{
  this->fStartPos  = pos;
  this->fDirection = dir;
}

//___________________________________________________________________________
bool PathSegmentList::IsSameStart(const TVector3& pos, const TVector3& dir) const
{
  return ( this->fStartPos == pos && this->fDirection == dir );
}

//___________________________________________________________________________
void PathSegmentList::FillMatStepSum(void) 
{
  fMatStepSum.clear();

  PathSegmentList::PathSegVCItr_t sitr = fSegmentList.begin();
  PathSegmentList::PathSegVCItr_t sitr_end = fSegmentList.end();
  for ( ; sitr != sitr_end ; ++sitr ) {
    const PathSegment& ps = *sitr;
    const TGeoMaterial* mat  = ps.fMaterial;
    // use the post-trim limits on how much material is stepped through
    fMatStepSum[mat] += ps.GetSummedStepRange();
  }

}

//___________________________________________________________________________
void PathSegmentList::Copy(const PathSegmentList & plist)
{
  fSegmentList.clear();
  fMatStepSum.clear();

  // copy the segments
  //vector<PathSegment>::const_iterator pl_iter;
  //for (pl_iter = plist.fSegmentList.begin(); pl_iter != plist.fSegmentList.end(); ++pl_iter) {
  //  this->fSegmentList.push_back( *pl_iter );
  //}

  // other elements
  fStartPos     = plist.fStartPos;
  fDirection    = plist.fDirection;
  fSegmentList  = plist.fSegmentList;
  fMatStepSum   = plist.fMatStepSum;
  fDoCrossCheck = plist.fDoCrossCheck;
  fPrintVerbose = plist.fPrintVerbose;
}

//___________________________________________________________________________
void PathSegmentList::CrossCheck(double& mxddist, double& mxdstep) const
{

  double dstep, ddist;
  mxdstep = 0;
  mxddist = 0;
  PathSegmentList::PathSegVCItr_t sitr = fSegmentList.begin();
  PathSegmentList::PathSegVCItr_t sitr_end = fSegmentList.end();
  for ( ; sitr != sitr_end ; ++sitr ) {
    const PathSegment& ps = *sitr;
    ps.DoCrossCheck(fStartPos,ddist,dstep);
    double addist = TMath::Abs(ddist);
    double adstep = TMath::Abs(dstep);
    if ( addist > mxddist ) mxddist = addist;
    if ( adstep > mxdstep ) mxdstep = adstep;
  }

}

//___________________________________________________________________________
void PathSegmentList::Print(ostream & stream) const
{
  stream << "\nPathSegmentList [-]" << endl;
  stream << " start " << pathsegutils::Vec3AsString(&fStartPos)
         << " dir " << pathsegutils::Vec3AsString(&fDirection) << endl;

  double dstep, ddist, mxdstep = 0, mxddist = 0;
  int k = 0, nseg = 0;
  PathSegmentList::PathSegVCItr_t sitr = fSegmentList.begin();
  PathSegmentList::PathSegVCItr_t sitr_end = fSegmentList.end();
  for ( ; sitr != sitr_end ; ++sitr, ++k ) {
    const PathSegment& ps = *sitr;
    ++nseg;
    stream << " [" << setw(4) << k << "] " << ps;
    if ( fDoCrossCheck ) {
      ps.DoCrossCheck(fStartPos,ddist,dstep);
      double addist = TMath::Abs(ddist);
      double adstep = TMath::Abs(dstep);
      if ( addist > mxddist ) mxddist = addist;
      if ( adstep > mxdstep ) mxdstep = adstep;
      stream << " recalc diff" 
             << " dist " << std::setw(12) << ddist
             << " step " << std::setw(12) << dstep;
        }
    stream << std::endl;
  }
  if ( nseg == 0 ) stream << " holds no segments." << std::endl;

  if ( fDoCrossCheck )
    stream << "PathSegmentList " 
           << " mxddist " << mxddist 
           << " mxdstep " << mxdstep 
           << std::endl;

  if ( fPrintVerbose ) {
    PathSegmentList::MaterialMapCItr_t mitr     = GetMatStepSumMap().begin();
    PathSegmentList::MaterialMapCItr_t mitr_end = GetMatStepSumMap().end();
    // loop over map to get tgt weight for each material (once)
    // steps outside the geometry may have no assigned material
    for ( ; mitr != mitr_end; ++mitr ) {
      const TGeoMaterial* mat = mitr->first;
      double sumsteps         = mitr->second;
      stream << " fMatStepSum[" << mat->GetName() << "] = " << sumsteps << std::endl;
    }
  }

}
//___________________________________________________________________________
#ifdef PATH_SEGMENT_SUPPORT_XML
XmlParserStatus_t PathSegmentList::LoadFromXml(string filename)
{
  this->clear();
  PDGLibrary * pdglib = PDGLibrary::Instance();

  LOG("PathS", pINFO)
               << "Loading PathSegmentList from XML file: " << filename;

  xmlDocPtr xml_doc = xmlParseFile(filename.c_str() );

  if(xml_doc==NULL) {
    LOG("PathS", pERROR)
           << "XML file could not be parsed! [filename: " << filename << "]";
    return kXmlNotParsed;
  }

  xmlNodePtr xmlCur = xmlDocGetRootElement(xml_doc);

  if(xmlCur==NULL) {
    LOG("PathS", pERROR)
        << "XML doc. has null root element! [filename: " << filename << "]";
    return kXmlEmpty;
  }

  if( xmlStrcmp(xmlCur->name, (const xmlChar *) "path_length_list") ) {
    LOG("PathS", pERROR)
     << "XML doc. has invalid root element! [filename: " << filename << "]";
    return kXmlInvalidRoot;
  }

  LOG("PathS", pINFO) << "XML file was successfully parsed";

  xmlCur = xmlCur->xmlChildrenNode; // <path_length>'s

  // loop over all xml tree nodes that are children of the root node
  while (xmlCur != NULL) {

    // enter everytime you find a <path_length> tag
    if( (!xmlStrcmp(xmlCur->name, (const xmlChar *) "path_length")) ) {

       xmlNodePtr xmlPlVal = xmlCur->xmlChildrenNode;

       string spdgc = utils::str::TrimSpaces(
                         XmlParserUtils::GetAttribute(xmlCur, "pdgc"));

       string spl = XmlParserUtils::TrimSpaces(
                           xmlNodeListGetString(xml_doc, xmlPlVal, 1));

       LOG("PathS", pDEBUG) << "pdgc = " << spdgc << " --> pl = " << spl;

       int    pdgc = atoi( spdgc.c_str() );
       double pl   = atof( spl.c_str()   );

       TParticlePDG * p = pdglib->Find(pdgc);
       if(!p) {
         LOG("PathS", pERROR)
            << "No particle with pdgc " << pdgc
                             << " found. Will not load its path length";
       } else
            this->insert( map<int, double>::value_type(pdgc, pl) );

       xmlFree(xmlPlVal);
    }
       xmlCur = xmlCur->next;
  } // [end of] loop over tags within root elements

  xmlFree(xmlCur);
  return kXmlOK;
}
//___________________________________________________________________________
void PathSegmentList::SaveAsXml(string filename) const
{
//! Save path length list to XML file

  LOG("PathS", pINFO)
          << "Saving PathSegmentList as XML in file: " << filename;

  ofstream outxml(filename.c_str());
  if(!outxml.is_open()) {
    LOG("PathS", pERROR) << "Couldn't create file = " << filename;
    return;
  }
  outxml << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>";
  outxml << endl << endl;
  outxml << "<!-- generated by PathSegmentList::SaveAsXml() -->";
  outxml << endl << endl;

  outxml << "<path_length_list>" << endl << endl;

  PathSegmentList::const_iterator pl_iter;

  for(pl_iter = this->begin(); pl_iter != this->end(); ++pl_iter) {

    int    pdgc = pl_iter->first;
    double pl   = pl_iter->second; // path length

    outxml << "   <path_length pdgc=\"" << pdgc << "\">"
                                 << pl << "</path_length>" << endl;
  }
  outxml << "</path_length_list>";
  outxml << endl;

  outxml.close();

}
#endif
//___________________________________________________________________________
PathSegmentList & PathSegmentList::operator = (const PathSegmentList & list)
{
  this->Copy(list);
  return (*this);
}
//___________________________________________________________________________
