//____________________________________________________________________________
/*!

\namespace  genie::utils::xml

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

#include "Utils/StringUtils.h"

#include <TSystem.h>

using std::string;

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

}         // xml   namespace
}         // utils namespace
}         // genie namespace

#endif    // _XML_UTILS_H_

