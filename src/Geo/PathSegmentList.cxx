//____________________________________________________________________________
/*
 Copyright (c) 2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Robert Hatcher <rhatcher@fnal.gov>
         FNAL - May 26, 2009

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

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

#include "Geo/PathSegmentList.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"
#include "Utils/MathUtils.h"
//#include "Utils/XmlParserUtils.h"

using std::ofstream;
using std::setw;
using std::setfill;
using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const PathSegment & ps)
 {
   ps.Print(stream);
   return stream;
 }
}
//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const PathSegmentList & list)
 {
   list.Print(stream);
   return stream;
 }
}

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
  fEnter(), fExit(), fStepLength(0), fRayDist(0),
  fVolume(0), fMedium(0), fMaterial(0)
{
}

void PathSegment::Print(ostream & stream) const
{
  stream << genie::pathsegutils::Vec3AsString(&fEnter) << " " 
         << genie::pathsegutils::Vec3AsString(&fExit) 
         << " raydist " << std::setw(12) << fRayDist 
         << " step " << std::setw(12) << fStepLength << " "
         << std::setw(16) << fVolume->GetName() << " '"
         << std::setw(18) << fMaterial->GetName() << "'";
}

//===========================================================================
//___________________________________________________________________________
PathSegmentList::PathSegmentList(void)
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
  for (unsigned int indx = 0; indx < fSegmentList.size(); ++indx) {
    PathSegment & ps = fSegmentList[indx];
    const TGeoMaterial* mat  = ps.fMaterial;
    double              step = ps.fStepLength;
    fMatStepSum[mat] += step;
  }

}



#ifdef UNNEEDED_SEGFUNCS
//___________________________________________________________________________
void PathSegmentList::AddPathLength(int pdgc, double pl)
{
// Adds pl to the total path length for material with code = pdgc

LOG("PathS", pFATAL)
  << "AddPathLength not implemented";
  //if (this->count(pdgc) == 1) { (*this)[pdgc] += pl; }
  //else {
  //   LOG("PathS", pWARN)
  //       << "No material with PDG code = " << pdgc << " in path length list";
  //}
}
//___________________________________________________________________________
void PathSegmentList::SetPathLength(int pdgc, double pl)
{
// Sets the total path length for material with code = pdgc to be pl

  LOG("PathS", pFATAL)
    << "SetPathLength not implemented";
  //if (this->count(pdgc) == 1) { (*this)[pdgc] = pl; }
  //else {
  //   LOG("PathS", pWARN)
  //       << "No material with PDG code = " << pdgc << " in path length list";
  //}
}
//___________________________________________________________________________
void PathSegmentList::ScalePathLength(int pdgc, double scale)
{
// Scales pl for material with code = pdgc with the input scale factor

  LOG("PathS", pFATAL)
    << "ScalePathLength not implemented";
  //if (this->count(pdgc) == 1) {
  //   double pl = (*this)[pdgc];
  //   pl *= scale;
  //   (*this)[pdgc] = pl;
  //} else {
  //   LOG("PathS", pWARN)
  //       << "No material with PDG code = " << pdgc << " in path length list";
  //}
}
//___________________________________________________________________________
double PathSegmentList::PathLength(int pdgc) const
{
// Gets the total path length for material with code = pdgc to be pl

  LOG("PathS", pFATAL)
    << "PathLength not implemented";
  //if ( this->count(pdgc) == 1 ) {
  //   map<int, double>::const_iterator pl_iter = this->find(pdgc);
  //   return pl_iter->second;
  //} else {
  //   LOG("PathS", pWARN)
  //       << "No material with PDG code = " << pdgc << " in path length list";
  //}
  return 0;
}
//___________________________________________________________________________
bool PathSegmentList::AreAllZero(void) const
{
  bool allzero = true;

  LOG("PathS", pFATAL)
    << "AddPathLength not implemented";
  //PathSegmentList::const_iterator pl_iter;
  //
  //for(pl_iter = this->begin(); pl_iter != this->end(); ++pl_iter) {
  //  double pl = pl_iter->second;
  //  allzero = allzero && (utils::math::AreEqual(pl,0.));
  //}
  return allzero;
}
#endif
//___________________________________________________________________________
void PathSegmentList::Copy(const PathSegmentList & plist)
{
  fSegmentList.clear();
  fMatStepSum.clear();

  // copy the segments
  vector<PathSegment>::const_iterator pl_iter;
  for (pl_iter = plist.fSegmentList.begin(); pl_iter != plist.fSegmentList.end(); ++pl_iter) {
    this->fSegmentList.push_back( *pl_iter );
  }

  // other elements
  fStartPos   = plist.fStartPos;
  fDirection  = plist.fDirection;
  fMatStepSum = plist.fMatStepSum;

}
//___________________________________________________________________________
void PathSegmentList::Print(ostream & stream) const
{
  stream << "\nPathSegmentList [-]" << endl;

  //vector<PathSegment>::const_iterator ps_iter;
  //for(ps_iter = fSegmentList.begin(); ps_iter != fSegmentList.end(); ++ps_iter) {
  //  stream << *ps_iter;
  //}

  double mxdstep = 0;
  double mxddist = 0;
  size_t nps = fSegmentList.size();
  for(unsigned int k = 0; k < nps ; ++k ) {
    double raydist_recalc = (fSegmentList[k].fEnter - this->fStartPos).Mag();
    double step_recalc = (fSegmentList[k].fExit - fSegmentList[k].fEnter).Mag();
    double dstep = step_recalc-fSegmentList[k].fStepLength;
    double ddist = raydist_recalc-fSegmentList[k].fRayDist;
    if ( dstep > mxdstep ) mxdstep = dstep;
    if ( ddist > mxddist ) mxddist = ddist;
    stream << "[" << setw(4) << k << "] " << fSegmentList[k]
           << " raydist_recalc " << std::setw(12) << raydist_recalc 
           << " d " << std::setw(12) << ddist
           << " step " << std::setw(12) << step_recalc 
           << " " << std::setw(12) << dstep
           << std::endl;

  }
  //rwh
  stream << "PathSegmentList mxddist " << mxddist << " mxdstep " << mxdstep << endl;
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
