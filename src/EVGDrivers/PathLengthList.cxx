//____________________________________________________________________________
/*!

\class   genie::PathLengthList

\brief   Object to be filled with the neutrino path-length, for all detector
         geometry materials, when starting from a position x and travelling
         along the direction of the neutrino 4-momentum.

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 24, 2005

*/
//____________________________________________________________________________

#include <fstream>
#include <iomanip>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

#include <TLorentzVector.h>

#include "EVGDrivers/PathLengthList.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodeList.h"
#include "PDG/PDGLibrary.h"
#include "Utils/PrintUtils.h"
#include "Utils/MathUtils.h"
#include "Utils/XmlParserUtils.h"

using std::ofstream;
using std::setw;
using std::setfill;
using std::endl;

using namespace genie;

//____________________________________________________________________________
namespace genie {
 ostream & operator << (ostream & stream, const PathLengthList & list)
 {
   list.Print(stream);
   return stream;
 }
}
//___________________________________________________________________________
PathLengthList::PathLengthList(const PDGCodeList & pdg_list) :
map<int, double>()
{
  PDGCodeList::const_iterator pdg_iter;

  for(pdg_iter = pdg_list.begin(); pdg_iter != pdg_list.end(); ++pdg_iter) {

     int pdgc = *pdg_iter;
     this->insert( map<int, double>::value_type(pdgc, 0.) );
  }
}
//___________________________________________________________________________
PathLengthList::~PathLengthList()
{

}
//___________________________________________________________________________
void PathLengthList::AddPathLength(int pdgc, double pl)
{
// Adds pl to the total path length for material with code = pdgc

  if (this->count(pdgc) == 1) { (*this)[pdgc] += pl; }
  else {
     LOG("PathL", pWARN)
         << "No material with PDG code = " << pdgc << " in path length list";
  }
}
//___________________________________________________________________________
void PathLengthList::SetPathLength(int pdgc, double pl)
{
// Sets the total path length for material with code = pdgc to be pl

  if (this->count(pdgc) == 1) { (*this)[pdgc] = pl; }
  else {
     LOG("PathL", pWARN)
         << "No material with PDG code = " << pdgc << " in path length list";
  }
}
//___________________________________________________________________________
void PathLengthList::ScalePathLength(int pdgc, double scale)
{
// Scales pl for material with code = pdgc with the input scale factor

  if (this->count(pdgc) == 1) {
     double pl = (*this)[pdgc];
     pl *= scale;
     (*this)[pdgc] = pl;
  } else {
     LOG("PathL", pWARN)
         << "No material with PDG code = " << pdgc << " in path length list";
  }
}
//___________________________________________________________________________
double PathLengthList::PathLength(int pdgc) const
{
// Gets the total path length for material with code = pdgc to be pl

  if ( this->count(pdgc) == 1 ) {
     map<int, double>::const_iterator pl_iter = this->find(pdgc);
     return pl_iter->second;
  } else {
     LOG("PathL", pWARN)
         << "No material with PDG code = " << pdgc << " in path length list";
  }
  return 0;
}
//___________________________________________________________________________
void PathLengthList::SetAllToZero(void)
{
  PathLengthList::const_iterator pl_iter;

  for(pl_iter = this->begin(); pl_iter != this->end(); ++pl_iter) {
    int pdgc = pl_iter->first;
    (*this)[pdgc] = 0.;
  }
}
//___________________________________________________________________________
bool PathLengthList::AreAllZero(void) const
{
  bool allzero = true;

  PathLengthList::const_iterator pl_iter;

  for(pl_iter = this->begin(); pl_iter != this->end(); ++pl_iter) {
    double pl = pl_iter->second;
    allzero = allzero && (utils::math::AreEqual(pl,0.));
  }
  return allzero;
}
//___________________________________________________________________________
void PathLengthList::Print(ostream & stream) const
{
  stream << "\n[-]" << endl;

  PDGLibrary * pdglib = PDGLibrary::Instance();

  PathLengthList::const_iterator pl_iter;
  size_t nc = this->size();

  for(pl_iter = this->begin(); pl_iter != this->end(); ++pl_iter) {

    int    pdgc = pl_iter->first;
    double pl   = pl_iter->second; // path length

    TParticlePDG * p = pdglib->Find(pdgc);

    if(!p) {
      stream << " |---o ** ERR: no particle with PDG code: " << pdgc;
    } else {
      string name = p->GetName();
      stream << " |---o code: " << pdgc << " [" << setfill(' ')
         << setw(5) << name << "] " << "-----> path-length = " << pl;
    }
    if( (--nc) > 0) stream << endl;
  }
}
//___________________________________________________________________________
XmlParserStatus_t PathLengthList::LoadAsXml(string filename)
{
  LOG("PathL", pINFO)
                    << "Loading PathLengthList from XML file: " << filename;

  xmlDocPtr xml_doc = xmlParseFile(filename.c_str() );

  if(xml_doc==NULL) {
    LOG("PathL", pERROR)
           << "XML file could not be parsed! [filename: " << filename << "]";
    return kXmlNotParsed;
  }

  xmlNodePtr xmlCur = xmlDocGetRootElement(xml_doc);

  if(xmlCur==NULL) {
    LOG("PathL", pERROR)
        << "XML doc. has null root element! [filename: " << filename << "]";
    return kXmlEmpty;
  }

  if( xmlStrcmp(xmlCur->name, (const xmlChar *) "path_length_list") ) {
    LOG("PathL", pERROR)
     << "XML doc. has invalid root element! [filename: " << filename << "]";
    return kXmlInvalidRoot;
  }

  LOG("PathL", pINFO) << "XML file was successfully parsed";

  xmlCur = xmlCur->xmlChildrenNode; // <path_length>'s

  // loop over all xml tree nodes that are children of the root node
  while (xmlCur != NULL) {
    // enter everytime you find a <path_length> tag
    if( (!xmlStrcmp(xmlCur->name, (const xmlChar *) "path_length")) ) {

       string spdgc = utils::str::TrimSpaces(
                         XmlParserUtils::GetAttribute(xmlCur, "pdgc"));

       string spl = XmlParserUtils::TrimSpaces(
                             xmlNodeListGetString(xml_doc, xmlCur, 1));

       int    pdgc = atoi( spdgc.c_str() );
       double pl   = atof( spl.c_str()   );

       this->AddPathLength(pdgc, pl);

       xmlCur = xmlCur->next;
    }
  } // [end of] loop over tags within root elements

  xmlFree(xmlCur);
  return kXmlOK;
}
//___________________________________________________________________________
void PathLengthList::SaveAsXml(string filename) const
{
//! Save path length list to XML file

  LOG("PathL", pINFO)
          << "Saving PathLengthList as XML in file: " << filename;

  ofstream outxml(filename.c_str());
  if(!outxml.is_open()) {
    LOG("PathL", pERROR) << "Couldn't create file = " << filename;
    return;
  }
  outxml << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>";
  outxml << endl << endl;
  outxml << "<!-- generated by PathLengthList::SaveAsXml() -->";
  outxml << endl << endl;

  outxml << "<path_length_list>" << endl << endl;

  PathLengthList::const_iterator pl_iter;

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
//___________________________________________________________________________


