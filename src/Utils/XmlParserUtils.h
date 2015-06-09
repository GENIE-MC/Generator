//____________________________________________________________________________
/*!

\  genie::utils::xml

\brief      XML utilities

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    May 04, 2004
 
\cpright    Copyright (c) 2003-2015, GENIE Neutrino MC Generator Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _XML_UTILS_H_
#define _XML_UTILS_H_

#include <string>
#include <vector>

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#endif

#include <TSystem.h>
#include <TVectorT.h>

#include "Utils/StringUtils.h"

class TFile;
class TH1F;
class TH1D;
class TH2D;
class TSpline3;

using std::string;
using std::vector;

namespace genie {
namespace utils {
namespace xml   {

#if !defined(__CINT__) && !defined(__MAKECINT__)
  inline string TrimSpaces(xmlChar * xmls)
  {
   // trim the leading/trailing spaces from an parsed xml string like in:
   //
   // "      I am a string with lots of spaces      " ---->
   //                                  "I am a string with lots of spaces"
   //
   // In this method, "\n" is treated as 'empty space' so as to trim not only
   // empty spaces in the line that contains the string but also all leading
   // and trailing empty lines
  
    string str = string( (const char *) xmls );
    return utils::str::TrimSpaces(str);
  }

  //_________________________________________________________________________
  inline string GetAttribute(xmlNodePtr xml_cur, string attr_name)
  {
    xmlChar * xmls = xmlGetProp(xml_cur, (const xmlChar *) attr_name.c_str());    
    string str = TrimSpaces(xmls);
    xmlFree(xmls);
    return str;
  }
#endif

  //_________________________________________________________________________
  inline string GetXMLPathList()
  {
    // Get a colon separated list of potential locations for xml files
    // e.g. ".:$MYSITEXML:/path/to/exp/version:$GALGCONF:$GENIE/config"
    // user additions should be in $GXMLPATH
    string pathlist; 
    const char* p1 = gSystem->Getenv("GXMLPATH");
    if ( p1 ) { pathlist = std::string(p1) + ":"; }
    const char* p2 = gSystem->Getenv("GXMLPATHS");  // handle extra 's'
    if ( p2 ) { pathlist = std::string(p2) + ":"; }
    // add originally supported alternative path
    const char* p3 = gSystem->Getenv("GALGCONF");
    if ( p3 ) { pathlist = std::string(p3) + ":"; }
    pathlist += "$GENIE/config";  // standard path in case no env
    pathlist += ":$GENIE/src/FluxDrivers/GNuMINtuple";  // special case
    return pathlist;
  }

  //_________________________________________________________________________
  inline string GetXMLFilePath(string basename)
  {
    // return a full path to a real XML file
    // e.g. passing in "GNuMIFlux.xml"
    //   will return   "/blah/GENIE/HEAD/config/GNuMIFlux.xml"
    // allow ::colon:: ::semicolon:: and ::comma:: as path item separators
    std::string pathlist = genie::utils::xml::GetXMLPathList();
    std::vector<std::string> paths = genie::utils::str::Split(pathlist,":;,");
    // expand any wildcards, etc.
    size_t np = paths.size();
    for ( size_t i=0; i< np; ++i ) {
      const char* tmppath = paths[i].c_str();
      std::string onepath = gSystem->ExpandPathName(tmppath);
      onepath += "/";
      onepath += basename;
      bool noAccess = gSystem->AccessPathName(onepath.c_str());
      if ( ! noAccess ) return onepath;  // found one
    }
    // didn't find any, return basename in case it is in "." and that
    // wasn't listed in the XML path list.   If you want "." to take
    // precedence then it needs to be explicitly listed in GXMLPATHS.
    return basename;  
  }

  //_________________________________________________________________________

#if !defined(__CINT__) && !defined(__MAKECINT__)

  // Find a particular node witin the input XML document.
  // The node is specified using the input path.
  // For example, to retrieve node <superk_energy_scale_err>
  // in XML doc below
  // <t2k>
  //   <systematics>
  //      <superk_energy_scale_err>
  //          0.015
  //      </superk_energy_scale_err>
  //   </systematics>
  // </t2k>
  // specify the path "t2k/systematics/superk_energy_scale_err"
  //
  xmlNodePtr FindNode(xmlDocPtr xml_doc, string node_path);

  //
  // Retrieve XML file data in various formats.
  // To retrieve a ROOT object from within a ROOT file, the following XML scheme is used
  // <some_node>
  //      <another_node>
  //            <filename> blah </filename>
  //            <objname>  blah </objname>
  //            <objtype>  blah </objtype>
  //      </another_node>
  // </some_node>
  //
  bool           GetBool        (xmlDocPtr xml_doc, string node_path);       
  int            GetInt         (xmlDocPtr xml_doc, string node_path);
  vector<int>    GetIntArray    (xmlDocPtr xml_doc, string node_path); // comma-separated values in XML file
  double         GetDouble      (xmlDocPtr xml_doc, string node_path);
  vector<double> GetDoubleArray (xmlDocPtr xml_doc, string node_path); // comma-separated values in XML file
  string         GetString      (xmlDocPtr xml_doc, string node_path);
  string         GetROOTFileName(xmlDocPtr xml_doc, string node_path);
  string         GetROOTObjName (xmlDocPtr xml_doc, string node_path);
  string         GetROOTObjType (xmlDocPtr xml_doc, string node_path);
  TFile *        GetTFile       (xmlDocPtr xml_doc, string node_path, string base_dir = "<env>");
  TH1F *         GetTH1F        (xmlDocPtr xml_doc, string node_path, string base_dir = "<env>");  
  TH1D *         GetTH1D        (xmlDocPtr xml_doc, string node_path, string base_dir = "<env>");  
  TH2D *         GetTH2D        (xmlDocPtr xml_doc, string node_path, string base_dir = "<env>");  
  TVectorD *     GetTVectorD    (xmlDocPtr xml_doc, string node_path, string base_dir = "<env>");  
/*
  TMatrixDSym *  GetTMatrixDSym (xmlDocPtr xml_doc, string node_path, string base_dir = "<env>");  
  TMatrixD *     GetTMatrixD    (xmlDocPtr xml_doc, string node_path, string base_dir = "<env>");  
  TSpline3 *     GetTSpline3    (xmlDocPtr xml_doc, string node_path, string base_dir = "<env>");
*/
#endif

}         // xml   namespace
}         // utils namespace
}         // genie namespace

#endif    // _XML_UTILS_H_

