//_____________________________________________________________________________
/*!

\class    genie::nuvld::RecordStructure

\brief

\author   Costas Andreopoulos (Rutherford Lab.)  <costas.andreopoulos \at stfc.ac.uk>

\created  August 2003          
*/
//_____________________________________________________________________________

#include <fstream>
#include <iostream>

#include "RecordStructure.h"

using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::ios;

namespace genie {
namespace nuvld {
  
//_____________________________________________________________________________
ostream & operator << (ostream & stream, const RecordStructure & recbase)
{
  map<string, vector<string> >::const_iterator element_iter;
  vector<string>::const_iterator               attribute_iter;

  stream << "printing record structure" << endl;
  for(element_iter = recbase._elements.begin(); 
                    element_iter != recbase._elements.end(); ++element_iter) {

      stream << "[-] element: " << (*element_iter).first << endl;;
      vector<string> attributes = (*element_iter).second;
      for(attribute_iter = attributes.begin(); 
                       attribute_iter != attributes.end(); ++attribute_iter) {

         stream << " |----> with attribute: " << (*attribute_iter) << endl;
      }
  }
  return stream;
}
//_____________________________________________________________________________
RecordStructure::RecordStructure()
{

}
//_____________________________________________________________________________
RecordStructure::RecordStructure(const char * filename)
{
  ReadElements(filename);
}
//_____________________________________________________________________________
RecordStructure::~RecordStructure()
{

}
//_____________________________________________________________________________
int RecordStructure::ReadElements(const char * filename)
{
  ifstream inp(filename,ios::in);

  if( inp.fail() ) {
    cerr << "Failed to open filename: " << filename << endl;
    exit(1);
  }

  string         name;
  vector<string> names;
  while ( inp >> name ) names.push_back(name);

  vector<string>::reverse_iterator name_iter;
  vector<string> attributes;

  for(name_iter = names.rbegin(); name_iter != names.rend(); ++name_iter) {
     if( (*name_iter)[0] == '+') {
         string element = (*name_iter).erase(0,1);
         _elements.insert( 
               map<string, vector<string> >::value_type(element, attributes));
         attributes.clear();                    
     } else {
        attributes.push_back( (*name_iter).erase(0,1) );
     }
  }

  inp.close();

  return 0;
}
//_____________________________________________________________________________

} // nuvld namespace
} // genie namespace
