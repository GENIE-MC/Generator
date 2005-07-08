//____________________________________________________________________________
/*!

\class    genie::XSecSplineList

\brief    List of cross section vs energy splines

\author   Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
          CCLRC, Rutherford Appleton Laboratory

\created  May 12, 2005
*/
//____________________________________________________________________________

#include <fstream>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"

#include <TMath.h>
#include <TLorentzVector.h>

#include "Base/XSecAlgorithmI.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "Utils/StringUtils.h"
#include "Utils/XSecSplineList.h"
#include "XML/XmlParserUtils.h"

using std::ofstream;
using std::endl;

namespace genie {

//____________________________________________________________________________
ostream & operator<< (ostream& stream, const XSecSplineList & list)
{
  list.Print(stream);

  return stream;
}
//____________________________________________________________________________
XSecSplineList * XSecSplineList::fInstance = 0;
//____________________________________________________________________________
XSecSplineList::XSecSplineList()
{
  fInstance    =  0;
  fUseLogE     = true;
  fNKnots      = 100;
  fEmin        =   0.01; // GeV
  fEmax        = 100.00; // GeV
  fEExtrap     = -1.;
}
//____________________________________________________________________________
XSecSplineList::~XSecSplineList()
{
  fInstance = 0;
}
//____________________________________________________________________________
XSecSplineList * XSecSplineList::Instance()
{
  if(fInstance == 0) {

    static XSecSplineList::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();

    fInstance = new XSecSplineList;
  }
  return fInstance;
}
//____________________________________________________________________________
bool XSecSplineList::SplineExists(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string key = this->BuildSplineKey(alg,interaction);

  if(fSplineMap.count(key) == 1) return true;

  return false;
}
//____________________________________________________________________________
const Spline * XSecSplineList::GetSpline(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string key = this->BuildSplineKey(alg,interaction);

  if ( this->SplineExists(alg,interaction) ) {

     map<string, Spline *>::const_iterator iter = fSplineMap.find(key);
     return iter->second;

  } else {
    SLOG("XSecSplineList", pWARN) << "Couldn't find spline for key = " << key;
    return 0;
  }
  return 0;
}
//____________________________________________________________________________
void XSecSplineList::CreateSpline(const XSecAlgorithmI * alg,
        const Interaction * interaction, int nknots, double Emin, double Emax)
{
// Build a cross section spline for the input interaction using the input
// cross section algorithm and store in the list.
// For building this specific entry of the spline list, the user is allowed
// to override the list-wide nknots,Emin,Emax

  SLOG("XSecSplineList", pINFO)
             << "Creating cross section spline using the algorithm: " << *alg;

  string key = this->BuildSplineKey(alg,interaction);

  // if any of the nknots,Emin,Emax was not set or its value is not acceptable
  // use the list values
  if (Emin   < 0.) Emin   = this->Emin();
  if (Emax   < 0.) Emin   = this->Emax();
  if (nknots <= 2) nknots = this->NKnots();

  assert(Emin < Emax);

  double xsec[nknots];
  double E   [nknots];

  double dE = 0;
  if( this->UseLogE() )
       dE = (TMath::Log10(Emax) - TMath::Log10(Emin)) /(nknots-1);
  else dE = (Emax-Emin) /(nknots-1);

  for (int i = 0; i < nknots; i++) {

    if( this->UseLogE() )
          E[i] = TMath::Power(10., TMath::Log10(Emin) + i * dE);
    else  E[i] = Emin + i * dE;

    TLorentzVector p4(0,0,E[i],E[i]);
    interaction->GetInitialStatePtr()->SetProbeP4(p4);

    xsec[i] = alg->XSec(interaction);

    SLOG("XSecSplineList", pINFO)<< "xsec(E = " << E[i] << ") = " << xsec[i];
  }
  Spline * spline = new Spline(nknots, E, xsec);

  fSplineMap.insert( map<string, Spline *>::value_type(key, spline) );
}
//____________________________________________________________________________
void XSecSplineList::SetLogE(bool on)
{
  fUseLogE = on;
}
//____________________________________________________________________________
void XSecSplineList::SetNKnots(int nk)
{
  fNKnots = nk;

  if(fNKnots<2) fNKnots = 2; // minimum acceptable number of knots
}
//____________________________________________________________________________
void XSecSplineList::SetMinE(double Ev)
{
  fEmin = Ev;

  if(fEmin<0) fEmin = 0.;
}
//____________________________________________________________________________
void XSecSplineList::SetMaxE(double Ev)
{
  fEmax = Ev;

  if(fEmax<0) fEmax = 0.;
}
//____________________________________________________________________________
void XSecSplineList::SaveAsXml(string filename) const
{
//! Save XSecSplineList to XML file

  SLOG("XSecSplineList", pINFO)
                     << "Saving XSecSplineList as XML in file: " << filename;

  ofstream outxml(filename.c_str());
  if(!outxml.is_open()) {
    SLOG("XSecSplineList", pERROR) << "Couldn't create file = " << filename;
    return;
  }
  outxml << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>";
  outxml << endl << endl;
  outxml << "<!-- generated by genie::XSecSplineList::SaveSplineList() -->";
  outxml << endl << endl;

  int uselog = (fUseLogE ? 1 : 0);
  outxml << "<genie_xsec_spline_list "
                          << "version=\"1.00\" uselog=\"" << uselog << "\">";
  outxml << endl << endl;

  map<string, Spline *>::const_iterator mapiter;

  for(mapiter = fSplineMap.begin(); mapiter != fSplineMap.end(); ++mapiter) {

    string     key     = mapiter->first;
    Spline *   spline  = mapiter->second;

    spline->SaveAsXml(outxml,"E","xsec", key, true);
  }
  outxml << "</genie_xsec_spline_list>";
  outxml << endl;

  outxml.close();
}
//____________________________________________________________________________
XmlParserStatus_t XSecSplineList::LoadFromXml(string filename, bool keep)
{
//! Load XSecSplineList from ROOT file. If keep = true, then the loaded splines
//! are added to the existing list. If false, then the existing list is reseted
//! before loading the splines.

  SLOG("XSecSplineList", pINFO)
       << "Loading XSecSplineList from XML file: " << filename
                   << " with flag to keep existing splines switched "
                                                 << ( (keep) ? "ON" : "OFF" );
  if(!keep) fSplineMap.clear();

  xmlDocPtr xml_doc = xmlParseFile(filename.c_str() );

  if(xml_doc==NULL) {
    LOG("XSecSplineList", pERROR)
          << "\nXML file could not be parsed! [filename: " << filename << "]";
    return kXmlNotParsed;
  }

  xmlNodePtr xmlCur = xmlDocGetRootElement(xml_doc);

  if(xmlCur==NULL) {
    LOG("XSecSplineList", pERROR)
        << "\nXML doc. has null root element! [filename: " << filename << "]";
    return kXmlEmpty;
  }
  if( xmlStrcmp(xmlCur->name,
                 (const xmlChar *) "genie_xsec_spline_list") ) {
    LOG("XSecSplineList", pERROR)
     << "\nXML doc. has invalid root element! [filename: " << filename << "]";
    return kXmlInvalidRoot;
  }

  SLOG("XSecSplineList", pINFO) << "XML file was successfully parsed";

  // read/set the uselog attribute
  string uselog = string_utils::TrimSpaces(
                             XmlParserUtils::GetAttribute(xmlCur, "uselog"));
  if (atoi(uselog.c_str()) == 1) this->SetLogE(true);
  else this->SetLogE(false);

  xmlCur = xmlCur->xmlChildrenNode; // <spline>'s

  // loop over all xml tree nodes that are children of the root node
  while (xmlCur != NULL) {

    // enter everytime you find a <spline> tag
    if( (!xmlStrcmp(xmlCur->name, (const xmlChar *) "spline")) ) {

       string name = string_utils::TrimSpaces(
                             XmlParserUtils::GetAttribute(xmlCur, "name"));
       string snkn = string_utils::TrimSpaces(
                           XmlParserUtils::GetAttribute(xmlCur, "nknots"));
       int nknots = atoi( snkn.c_str() );

       LOG("XSecSplineList", pINFO)
             << "Parsing XML spline: " << name << ", nknots = " << nknots;

       int      iknot = 0;
       double * E     = new double[nknots];
       double * xsec  = new double[nknots];

       xmlNodePtr xmlSplChild = xmlCur->xmlChildrenNode; // <knots>'s

       // loop over all xml tree nodes that are children of the <spline> node
       while (xmlSplChild != NULL) {
         LOG("XSecSplineList", pDEBUG)
                      << "Got <spline> children node: " << xmlSplChild->name;

         // enter everytime you find a <knot> tag
         if( (!xmlStrcmp(xmlSplChild->name, (const xmlChar *) "knot")) ) {

           xmlNodePtr xmlKnotChild = xmlSplChild->xmlChildrenNode;

            // loop over all xml tree nodes that are children of this <knot>
            while (xmlKnotChild != NULL) {
               LOG("XSecSplineList", pDEBUG)
                       << "Got <knot> children node: "  << xmlKnotChild->name;

              // enter everytime you find a <E> or a <xsec> tag
              const xmlChar * tag = xmlKnotChild->name;
              bool is_E    = ! xmlStrcmp(tag,(const xmlChar *) "E");
              bool is_xsec = ! xmlStrcmp(tag,(const xmlChar *) "xsec");
              if (is_E || is_xsec) {
                 xmlNodePtr xmlValTagChild = xmlKnotChild->xmlChildrenNode;
                 string val = XmlParserUtils::TrimSpaces(
                             xmlNodeListGetString(xml_doc, xmlValTagChild, 1));

                 if (is_E   ) E   [iknot] = atof(val.c_str());
                 if (is_xsec) xsec[iknot] = atof(val.c_str());

                xmlFree(xmlValTagChild);
                LOG("XSecSplineList", pDEBUG)
                                       << "tag: " << tag << ", value: " << val;
              }//if current tag is <E>,<xsec>

              xmlKnotChild = xmlKnotChild->next;
             }//[end of] loop over tags within <knot>...</knot> tags
             xmlFree(xmlKnotChild);
             iknot++;
         } // if <knot>
         xmlSplChild = xmlSplChild->next;
       } //[end of] loop over tags within <spline>...</spline> tags
       xmlFree(xmlSplChild);

       // done looping over knots - build the spline
       Spline * spline = new Spline(nknots, E, xsec);

       // insert the spline to the list
       fSplineMap.insert( map<string, Spline *>::value_type(name,spline) );
      } // [end of] found & parsing a <spline> section

      xmlCur = xmlCur->next;
  } // [end of] loop over tags within root elements

  xmlFree(xmlCur);
  return kXmlOK;
}
//____________________________________________________________________________
string XSecSplineList::BuildSplineKey(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string alg_name  = alg->Name();
  string param_set = alg->ParamSet();
  string intkey    = interaction->AsString();

  string key = alg_name + "/" + param_set + "/" + intkey;

  return key;
}
//____________________________________________________________________________
void XSecSplineList::Print(ostream & stream) const
{
  stream << "\n ******************* XSecSplineList *************************";
  stream << "\n [-] Options:";
  stream << "\n  |";
  stream << "\n  |-----o  UseLogE..................." << fUseLogE;
  stream << "\n  |-----o  Spline Emin..............." << fNKnots;
  stream << "\n  |-----o  Spline Emax..............." << fEmin;
  stream << "\n  |-----o  Spline NKnots............." << fEmax;
  stream << "\n  |";
  stream << "\n [-] Available Splines:";
  stream << "\n  |";

  map<string, Spline *>::const_iterator mapiter;
  for(mapiter = fSplineMap.begin(); mapiter != fSplineMap.end(); ++mapiter) {
    string key = mapiter->first;
    stream << "\n  |-----o  " << key;
  }
  stream << "\n";
}
//___________________________________________________________________________

} // genie namespace
