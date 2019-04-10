//____________________________________________________________________________
/*!

\  genie::utils::xml

\brief      XML utilities

\author     Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
            University of Liverpool & STFC Rutherford Appleton Lab

\created    May 04, 2004
 
\cpright    Copyright (c) 2003-2019, The GENIE Collaboration
            For the full text of the license visit http://copyright.genie-mc.org
            or see $GENIE/LICENSE
*/
//____________________________________________________________________________

#ifndef _XML_UTILS_H_
#define _XML_UTILS_H_

#include <string>
#include <vector>
#include <cstdlib>


#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "libxml/parser.h"
#include "libxml/xmlmemory.h"
#endif

#include <TSystem.h>
#include <TVectorT.h>

#include "Framework/Utils/StringUtils.h"
#include "Framework/Utils/RunOpt.h"


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

  string TrimSpaces(xmlChar * xmls) ;
  // trim the leading/trailing spaces from an parsed xml string like in:
  //
  // "      I am a string with lots of spaces      " ---->
  //                                  "I am a string with lots of spaces"
  //
  // In this method, "\n" is treated as 'empty space' so as to trim not only
  // empty spaces in the line that contains the string but also all leading
  // and trailing empty lines
  
  //_________________________________________________________________________

  string GetAttribute(xmlNodePtr xml_cur, string attr_name) ;
#endif

  //_________________________________________________________________________
  string GetXMLPathList( bool add_tune = true ) ;
  // Get a colon separated list of potential locations for xml files
  // e.g. ".:$MYSITEXML:/path/to/exp/version:$GALGCONF:$GENIE/config"
  // user additions should be in $GXMLPATH

  //_________________________________________________________________________
  inline string GetXMLDefaultPath() { return "$GENIE/config" ; }
  //standard path in case no env variable are set

  //_________________________________________________________________________
  string GetXMLFilePath(string basename) ;
  // return a full path to a real XML file
  // e.g. passing in "GNuMIFlux.xml"
  //   will return   "/blah/GENIE/HEAD/config/GNuMIFlux.xml"
  // allow ::colon:: ::semicolon:: and ::comma:: as path item separators
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

