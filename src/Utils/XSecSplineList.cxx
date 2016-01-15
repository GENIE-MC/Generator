//____________________________________________________________________________
/*
 Copyright (c) 2003-2016, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         University of Liverpool & STFC Rutherford Appleton Lab - May 12, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Jan 17, 2008 - CA
   Re-wrote LoadFromXml() and switched from the tree-based libxml2 API to
   the XmlTextReader API. That eliminates the need to load the whole XML file
   in memory and is faster. Works much better with the large XML spline files
   (~ 100 MB) used by MINOS.
 @ Jan 18, 2008 - CA
   Change the way spline knots are distributed in the given energy range so 
   that spline energy thresholds are handled more accurately. 
   One of the knots is sitting on the energy threshold for each interaction, 
   few knots are linearly spaced below threshold so that the spline behaves 
   ok for E<Ethr, and the bulk of the knots is spaced linearly or 
   logarithmically for E>Ethr.
 @ Jun 20, 2008 - CA
   Fix a memory leak in LoadFromXml(). Arrays were not deleted after splines
   instantiation. Also xmlChar* buffers got via xmlTextReaderGetAttribute()
   were not passed to xmlFree().
 @ Dec 06, 2008 - CA
   Tweak dtor so as not to clutter the output if GENIE exits in err so as to
   spot the fatal mesg immediately.
 @ Feb 25, 2010 - CA
   Exit immediately if the file pointed to by GSPLOAD isn't accessible.
 @ Sep 26, 2010 - CA
   Demote a few messages.
 @ Jan 24, 2013 - CA
   Use of variables $GSPLOAD and $GSPSAVE is no longer supported.

*/
//____________________________________________________________________________

#include <fstream>
#include <cstdlib>

#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#include "libxml/xmlreader.h"

#include <TSystem.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Units.h"
#include "Conventions/GBuild.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "Utils/StringUtils.h"
#include "Utils/PrintUtils.h"
#include "Utils/XSecSplineList.h"
#include "Utils/XmlParserUtils.h"

using std::ofstream;
using std::endl;

namespace genie {

//____________________________________________________________________________
ostream & operator << (ostream & stream, const XSecSplineList & list)
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
}
//____________________________________________________________________________
XSecSplineList::~XSecSplineList()
{
// Clean up. Don't clutter output if exiting in err.

  this->AutoSave();

  map<string, Spline *>::const_iterator spliter;
  for(spliter = fSplineMap.begin(); spliter != fSplineMap.end(); ++spliter) {
    Spline * spline = spliter->second;
    if(spline) {
      delete spline;
      spline = 0;
    }
  }
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
  return this->SplineExists(key);
}
//____________________________________________________________________________
bool XSecSplineList::SplineExists(string key) const
{
  SLOG("XSecSplLst", pDEBUG) << "Checking for spline with key = " << key;

  bool exists = (fSplineMap.count(key) == 1);
  SLOG("XSecSplLst", pDEBUG)
    << "Spline found?...." << utils::print::BoolAsYNString(exists);
  return exists;
}
//____________________________________________________________________________
const Spline * XSecSplineList::GetSpline(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  string key = this->BuildSplineKey(alg,interaction);
  return this->GetSpline(key);
}
//____________________________________________________________________________
const Spline * XSecSplineList::GetSpline(string key) const
{
  if ( this->SplineExists(key) ) {
     map<string, Spline *>::const_iterator iter = fSplineMap.find(key);
     return iter->second;
  } else {
    SLOG("XSecSplLst", pWARN) << "Couldn't find spline for key = " << key;
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

  double xsec[nknots];
  double E   [nknots];

  SLOG("XSecSplLst", pNOTICE)
     << "Creating cross section spline using the algorithm: " << *alg;

  string key = this->BuildSplineKey(alg,interaction);

  // If any of the nknots,Emin,Emax was not set or its value is not acceptable
  // use the list values
  //
  if (Emin   < 0.) Emin   = this->Emin();
  if (Emax   < 0.) Emax   = this->Emax();
  if (nknots <= 2) nknots = this->NKnots();
  assert(Emin < Emax);

  // Distribute the knots in the energy range (Emin,Emax) :
  // - Will use 5 knots linearly spaced below the energy thresholds so that the
  //   spline behaves correctly in (Emin,Ethr)
  // - Place 1 knot exactly on the input interaction threshold
  // - Place the remaining n-6 knots spaced either linearly or logarithmically 
  //   above the input interaction threshold
  // The above scheme schanges appropriately if Ethr<Emin (i.e. no knots
  // are computed below threshold)
  //
  double Ethr = interaction->PhaseSpace().Threshold();
  SLOG("XSecSplLst", pNOTICE)
    << "Energy threshold for current interaction = " << Ethr << " GeV";

  int nkb = (Ethr>Emin) ? 5 : 0; // number of knots <  threshold
  int nka = nknots-nkb;          // number of knots >= threshold

  // knots < energy threshold
  double dEb =  (Ethr>Emin) ? (Ethr - Emin) / nkb : 0;
  for(int i=0; i<nkb; i++) {     
     E[i] = Emin + i*dEb;
  }
  // knots >= energy threshold
  double E0  = TMath::Max(Ethr,Emin);
  double dEa = 0;
  if(this->UseLogE())
    dEa = (TMath::Log10(Emax) - TMath::Log10(E0)) /(nka-1);
  else 
    dEa = (Emax-E0) /(nka-1);

  for(int i=0; i<nka; i++) {
     if(this->UseLogE())
       E[i+nkb] = TMath::Power(10., TMath::Log10(E0) + i * dEa);
     else  
       E[i+nkb] = E0 + i * dEa;
  }

  // Compute cross sections for the input interaction at the selected
  // set of energies
  //
  for (int i = 0; i < nknots; i++) {
    TLorentzVector p4(0,0,E[i],E[i]);
    interaction->InitStatePtr()->SetProbeP4(p4);
    xsec[i] = alg->Integral(interaction);
    SLOG("XSecSplLst", pNOTICE)
            << "xsec(E = " << E[i] << ") = " 
                       << (1E+38/units::cm2)*xsec[i] << " x 1E-38 cm^2";
  }

  // Build & save the spline
  //
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
  if(fNKnots<10) fNKnots = 10; // minimum acceptable number of knots
}
//____________________________________________________________________________
void XSecSplineList::SetMinE(double Ev)
{
  if(Ev>0) fEmin = Ev;
}
//____________________________________________________________________________
void XSecSplineList::SetMaxE(double Ev)
{
  if(Ev>0) fEmax = Ev;
}
//____________________________________________________________________________
void XSecSplineList::SaveAsXml(string filename) const
{
//! Save XSecSplineList to XML file

  SLOG("XSecSplLst", pNOTICE)
       << "Saving XSecSplineList as XML in file: " << filename;

  ofstream outxml(filename.c_str());
  if(!outxml.is_open()) {
    SLOG("XSecSplLst", pERROR) << "Couldn't create file = " << filename;
    return;
  }
  outxml << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>";
  outxml << endl << endl;
  outxml << "<!-- generated by genie::XSecSplineList::SaveSplineList() -->";
  outxml << endl << endl;

  int uselog = (fUseLogE ? 1 : 0);
  outxml << "<genie_xsec_spline_list "
                          << "version=\"2.00\" uselog=\"" << uselog << "\">";
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

  SLOG("XSecSplLst", pNOTICE) << "Loading splines from: " << filename;
  SLOG("XSecSplLst", pINFO)
        << "Option to keep pre-existing splines is switched "
        << ( (keep) ? "ON" : "OFF" );

  if(!keep) fSplineMap.clear();

  const int kNodeTypeStartElement = 1;
  const int kNodeTypeEndElement   = 15;
  const int kKnotX                = 0;
  const int kKnotY                = 1;

  xmlTextReaderPtr reader;

  int ret = 0, val_type = -1, iknot = 0, nknots = 0;
  double * E = 0, * xsec = 0;
  string spline_name = "";

  reader = xmlNewTextReaderFilename(filename.c_str());
  if (reader != NULL) {
        ret = xmlTextReaderRead(reader);
        while (ret == 1) {
            xmlChar * name  = xmlTextReaderName     (reader);
            xmlChar * value = xmlTextReaderValue    (reader);
            int       type  = xmlTextReaderNodeType (reader);
            int       depth = xmlTextReaderDepth    (reader);

            if(depth==0 && type==kNodeTypeStartElement) {
               LOG("XSecSplLst", pDEBUG) << "Root element = " << name;
               if(xmlStrcmp(name, (const xmlChar *) "genie_xsec_spline_list")) {
                   LOG("XSecSplLst", pERROR)
                     << "\nXML doc. has invalid root element! [filename: " << filename << "]";
                   return kXmlInvalidRoot;
               }

               xmlChar * xvrs   = xmlTextReaderGetAttribute(reader,(const xmlChar*)"version");
               xmlChar * xinlog = xmlTextReaderGetAttribute(reader,(const xmlChar*)"uselog");
               string svrs      = utils::str::TrimSpaces((const char *)xvrs);
               string sinlog    = utils::str::TrimSpaces((const char *)xinlog);

               if (atoi(sinlog.c_str()) == 1) this->SetLogE(true);
               else this->SetLogE(false);

#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
               LOG("XSecSplLst", pDEBUG) << "Vrs   = " << svrs;
               LOG("XSecSplLst", pDEBUG) << "InLog = " << sinlog;
#endif
               xmlFree(xvrs);
               xmlFree(xinlog);
            }  

            if( (!xmlStrcmp(name, (const xmlChar *) "spline")) && type==kNodeTypeStartElement) {

               xmlChar * xname = xmlTextReaderGetAttribute(reader,(const xmlChar*)"name");
               xmlChar * xnkn  = xmlTextReaderGetAttribute(reader,(const xmlChar*)"nknots");
               string sname    = utils::str::TrimSpaces((const char *)xname);
               string snkn     = utils::str::TrimSpaces((const char *)xnkn);

               spline_name = sname;
               SLOG("XSecSplLst", pINFO) << "Loading spline: " << spline_name;

               nknots = atoi( snkn.c_str() );
               iknot=0;
               E     = new double[nknots];
               xsec  = new double[nknots];
  
               xmlFree(xname);
               xmlFree(xnkn);
          }
            if( (!xmlStrcmp(name, (const xmlChar *) "E"))    && type==kNodeTypeStartElement) { val_type = kKnotX; }
            if( (!xmlStrcmp(name, (const xmlChar *) "xsec")) && type==kNodeTypeStartElement) { val_type = kKnotY; }

            if( (!xmlStrcmp(name, (const xmlChar *) "#text")) && depth==4) {
                if      (val_type==kKnotX) E   [iknot] = atof((const char *)value);
                else if (val_type==kKnotY) xsec[iknot] = atof((const char *)value);
            }
            if( (!xmlStrcmp(name, (const xmlChar *) "knot")) && type==kNodeTypeEndElement) {
               iknot++;
            }
            if( (!xmlStrcmp(name, (const xmlChar *) "spline")) && type==kNodeTypeEndElement) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
               LOG("XSecSplLst", pINFO) << "Done with current spline";
               for(int i=0; i<nknots; i++) {
                  LOG("XSecSplLst", pINFO) << "xsec[E = " << E[i] << "] = " << xsec[i];
               }
#endif
               // done looping over knots - build the spline
               Spline * spline = new Spline(nknots, E, xsec);
               delete [] E;
               delete [] xsec;
               // insert the spline to the list
               fSplineMap.insert( map<string, Spline *>::value_type(spline_name,spline) );
            }
 
            xmlFree(name);
            xmlFree(value);
            ret = xmlTextReaderRead(reader);
        }
        xmlFreeTextReader(reader);
        if (ret != 0) {
          LOG("XSecSplLst", pERROR)
            << "\nXML file could not be parsed! [filename: " << filename << "]";
          return kXmlNotParsed;
        }
  } else {
    LOG("XSecSplLst", pERROR)
          << "\nXML file could not be found! [filename: " << filename << "]";
  }

  return kXmlOK;
}
//____________________________________________________________________________
/*
Below I am keeping a commented-out version the old LoadFromXml() method using 
the tree-based libxml2 API. That has been replaced by the above method using
the much faster XmlTextReader API.

XmlParserStatus_t XSecSplineList::LoadFromXml(string filename, bool keep)
{
//! Load XSecSplineList from ROOT file. If keep = true, then the loaded splines
//! are added to the existing list. If false, then the existing list is reseted
//! before loading the splines.

  SLOG("XSecSplLst", pNOTICE) << "Loading splines from: " << filename;
  SLOG("XSecSplLst", pINFO)
     << "Option to keep pre-existing splines is switched "
     << ( (keep) ? "ON" : "OFF" );

  if(!keep) fSplineMap.clear();

  xmlDocPtr xml_doc = xmlParseFile(filename.c_str() );

  if(xml_doc==NULL) {
    LOG("XSecSplLst", pERROR)
          << "\nXML file could not be parsed! [filename: " << filename << "]";
    return kXmlNotParsed;
  }

  xmlNodePtr xmlCur = xmlDocGetRootElement(xml_doc);

  if(xmlCur==NULL) {
    LOG("XSecSplLst", pERROR)
        << "\nXML doc. has null root element! [filename: " << filename << "]";
    return kXmlEmpty;
  }
  if( xmlStrcmp(xmlCur->name,
                 (const xmlChar *) "genie_xsec_spline_list") ) {
    LOG("XSecSplLst", pERROR)
     << "\nXML doc. has invalid root element! [filename: " << filename << "]";
    return kXmlInvalidRoot;
  }

  SLOG("XSecSplLst", pNOTICE) << "XML file was successfully parsed";

  // read/set the uselog attribute
  string uselog = utils::str::TrimSpaces(
                             XmlParserUtils::GetAttribute(xmlCur, "uselog"));
  if (atoi(uselog.c_str()) == 1) this->SetLogE(true);
  else this->SetLogE(false);

  xmlCur = xmlCur->xmlChildrenNode; // <spline>'s

  // loop over all xml tree nodes that are children of the root node
  while (xmlCur != NULL) {

    // enter everytime you find a <spline> tag
    if( (!xmlStrcmp(xmlCur->name, (const xmlChar *) "spline")) ) {

       string name = utils::str::TrimSpaces(
                             XmlParserUtils::GetAttribute(xmlCur, "name"));
       string snkn = utils::str::TrimSpaces(
                           XmlParserUtils::GetAttribute(xmlCur, "nknots"));
       int nknots = atoi( snkn.c_str() );

       BLOG("XSecSplLst", pNOTICE) << "Loading spline: " << name;
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
       SLOG("XSecSplLst", pDEBUG)  << "Number ok knots = " << nknots;
#endif
       int      iknot = 0;
       double * E     = new double[nknots];
       double * xsec  = new double[nknots];

       xmlNodePtr xmlSplChild = xmlCur->xmlChildrenNode; // <knots>'s

       // loop over all xml tree nodes that are children of the <spline> node
       while (xmlSplChild != NULL) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
         LOG("XSecSplLst", pDEBUG)
              << "Got <spline> children node: " << xmlSplChild->name;
#endif

         // enter everytime you find a <knot> tag
         if( (!xmlStrcmp(xmlSplChild->name, (const xmlChar *) "knot")) ) {

           xmlNodePtr xmlKnotChild = xmlSplChild->xmlChildrenNode;

            // loop over all xml tree nodes that are children of this <knot>
            while (xmlKnotChild != NULL) {
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
               LOG("XSecSplLst", pDEBUG)
                       << "Got <knot> children node: "  << xmlKnotChild->name;
#endif
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
#ifdef __GENIE_LOW_LEVEL_MESG_ENABLED__
                LOG("XSecSplLst", pDEBUG) << "tag: " << tag << ", value: " << val;
#endif
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
}*/
//____________________________________________________________________________
string XSecSplineList::BuildSplineKey(
            const XSecAlgorithmI * alg, const Interaction * interaction) const
{
  if(!alg) {
    LOG("XSecSplLst", pWARN)
            << "Null XSecAlgorithmI - Returning empty spline key";
    return "";
  }

  if(!interaction) {
    LOG("XSecSplLst", pWARN)
            << "Null Interaction - Returning empty spline key";
    return "";
  }

  string alg_name  = alg->Id().Name();
  string param_set = alg->Id().Config();
  string intkey    = interaction->AsString();

  string key = alg_name + "/" + param_set + "/" + intkey;

  return key;
}
//____________________________________________________________________________
bool XSecSplineList::AutoLoad(void)
{
/*
<disabled 24/01/13>

// Checks the $GSPLOAD env. variable and if found set reads the cross splines
// from the XML file it points to.
// Returns true if $GSPLOAD was set and splines were read / false otherwise.

  if( gSystem->Getenv("GSPLOAD") ) {

     string xmlfile = gSystem->Getenv("GSPLOAD");
     SLOG("XSecSplLst", pINFO) << "$GSPLOAD env.var = " << xmlfile;

     bool is_accessible = ! (gSystem->AccessPathName(xmlfile.c_str()));

     if(is_accessible) {
       SLOG("XSecSplLst", pINFO) << "Loading cross section splines";
       XmlParserStatus_t status = this->LoadFromXml(xmlfile);
       assert(status==kXmlOK);
       return true;
     } else {
       LOG("XSecSplLst", pFATAL)
         << "Specified XML file [" << xmlfile << "] is not accessible!";
       gAbortingInErr = true;
       exit(1);       
       return false;
     }
  }
  SLOG("XSecSplLst", pNOTICE)
      << "$GSPLOAD was not defined! No cross section splines will be loaded";
  return false;
*/

  if( gSystem->Getenv("GSPLOAD") ) {
     LOG("XSecSplLst", pFATAL)
      << "\n\n"
      << "********************************************************************************************** \n"
      << "The pre-computed cross-section data file can no longer be specified via the $GSPLOAD variable. \n"
      << "Please use the command-line option (typically --cross-sections) implemented in all GENIE apps \n"
      << "or, if in your user code you access XSecSplineList directly, use the \n"
      << "`XmlParserStatus_t XSecSplineList::LoadFromXml(string filename, bool keep)' method. \n"
      << "Unset $GSPLOAD to continue running GENIE. \n"
      << "********************************************************************************************** \n";
    gAbortingInErr = true;
    exit(1);
  }

  return true;
}
//____________________________________________________________________________
const vector<string> * XSecSplineList::GetSplineKeys(void) const
{
  vector<string> * keyv = new vector<string>(fSplineMap.size());

  unsigned int i=0;
  map<string, Spline *>::const_iterator mapiter;
  for(mapiter = fSplineMap.begin(); mapiter != fSplineMap.end(); ++mapiter) {
    string key = mapiter->first;
    (*keyv)[i++]=key;
  }
  return keyv;
}
//____________________________________________________________________________
void XSecSplineList::AutoSave(void)
{
/*
<disabled 24/01/13>

// Checks whether the $GSPSAVE env. variable and if found set it saves the
// cross section splines at the XML file this variable points to.
// You do not need to invoke this method (just set the env. variable). The
// method would be automatically called before the singleton destroys itself.

  if( gSystem->Getenv("GSPSAVE") ) {
     string xmlfile = gSystem->Getenv("GSPSAVE");
     cout << "Saving cross section splines to file: " << xmlfile << endl;
     this->SaveAsXml(xmlfile);
  }
*/

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

