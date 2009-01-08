//____________________________________________________________________________
/*!

\program   parsexml

\brief     Parses a NuValidator XML data file

\synopsis  parsexml -f filename.xml

\author    Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
           STFC, Rutherford Appleton Laboratory

\created   October 05, 2004
*/
//____________________________________________________________________________

#include <iostream>
#include <string>

#include "XmlParser/NuVldXmlParser.h"
#include "XmlParser/ParserStatus.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using namespace genie::nuvld;

const char * get_filename(int argc, char ** argv);

//__________________________________________________________________________________________
int main(int argc, char ** argv)
{
 const char * _filename = get_filename(argc, argv);
 
 cout << "parsing document.................: " << _filename << endl;

 NuVldXmlParser _xml_parser;

 _xml_parser.ParseXmlDocument( _filename );

 if(_xml_parser.GetXmlParsingStatus() == eXml_OK) {

    cout << "DONE!" << endl;
    return 0;

 } else {

    cout << "problems :-( " << endl;
    cout << ParserStatus::AsString( _xml_parser.GetXmlParsingStatus() ) << endl;
    return 1;
 }

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
