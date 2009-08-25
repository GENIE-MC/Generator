//____________________________________________________________________________
/*!

\program   gnuvld_xmlread

\brief     Parses a NuValidator XML data file

\synopsis  gnuvld_xmlread -f filename.xml

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           STFC, Rutherford Appleton Laboratory

\created   October 05, 2004
*/
//____________________________________________________________________________

#include <iostream>
#include <string>

#include "Conventions/XmlParserStatus.h"
#include "Messenger/Messenger.h"
#include "ValidationTools/NuVld/NuVldXmlParser.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

using namespace genie;
using namespace genie::nuvld;

const char * get_filename(int argc, char ** argv);

//__________________________________________________________________________________________
int main(int argc, char ** argv)
{
 const char * filename = get_filename(argc, argv);

 LOG("NuVld", pNOTICE)  << "Parsing document.....: " << filename;

 NuVldXmlParser xml_parser;
 xml_parser.ParseXmlDocument(filename);

 if(xml_parser.GetXmlParsingStatus() == kXmlOK) {
    LOG("NuVld", pNOTICE)  << "DONE!";
    return 0;
 } 

 LOG("NuVld", pERROR) 
   << "XML file parsing failed!";
 LOG("NuVld", pERROR) 
   << "Parser status: " << XmlParserStatus::AsString( xml_parser.GetXmlParsingStatus());

 return 1;
}
//__________________________________________________________________________________________
const char * get_filename(int argc, char ** argv)
{
  for(int iarg = 0; iarg < argc-1; iarg++) {
     string argument(argv[iarg]);
     if( argument.compare("-f") == 0 ) return argv[++iarg];
  }
  return "";
}
//__________________________________________________________________________________________
